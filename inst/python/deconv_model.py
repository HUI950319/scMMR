"""
DeconvModel: Deep neural network for bulk RNA-seq deconvolution.

Estimates cell type proportions from bulk RNA-seq expression using a
scRNA-seq reference. The model shares the same ResNet backbone architecture
as MultiTaskModel but uses:
- Continuous log-CPM input (NOT binary)
- Single regression task: proportions (softmax, sum=1)
- Training data: simulated pseudo-bulk from scRNA-seq

References:
- Scaden: Menden et al., 2020 (Science Advances)
- TAPE: Chen et al., 2022 (Nature Communications)
- DeSide: Li et al., 2024 (PNAS)
"""

import torch
import torch.nn as nn
import torch.nn.functional as F
import numpy as np
from torch.utils.data import Dataset, DataLoader
from sklearn.model_selection import train_test_split
import scipy.sparse as sp
from typing import List, Optional, Dict
import logging
import copy

from multi_task_model import InputBlock, ResNetBlock, ResNetLastBlock, TaskHead

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)


# =============================================================================
# Dataset
# =============================================================================

class DeconvDataset(Dataset):
    """Dataset for deconvolution: X = log-CPM, y = proportion vector."""

    def __init__(self, X: np.ndarray, y: np.ndarray):
        self.X = torch.FloatTensor(X)   # (n_samples, n_genes)
        self.y = torch.FloatTensor(y)   # (n_samples, n_cell_types)

    def __len__(self):
        return self.X.shape[0]

    def __getitem__(self, idx):
        return {'X': self.X[idx], 'y': self.y[idx]}


# =============================================================================
# Pseudo-bulk generation
# =============================================================================

def generate_pseudobulk(adata, label_column: str, n_samples: int = 5000,
                        n_cells_per_sample: int = 500,
                        var_genes: Optional[List[str]] = None,
                        target_sum: float = 1e6,
                        random_state: int = 42):
    """Generate simulated pseudo-bulk training data from scRNA-seq reference.

    For each pseudo-bulk sample:
    1. Sample proportions from Dirichlet(alpha=1)
    2. Allocate n_cells_per_sample cells across types
    3. Randomly sample cells (with replacement) and sum expression
    4. Normalize to log1p(CPM)

    Args:
        adata: AnnData with raw counts in .X
        label_column: obs column with cell type labels
        n_samples: number of pseudo-bulk samples to generate
        n_cells_per_sample: cells to mix per sample
        var_genes: subset of genes (if None, use all)
        target_sum: CPM normalization target (default 1e6)
        random_state: random seed

    Returns:
        dict with:
            X_pb: np.ndarray (n_samples, n_genes) float32 log-CPM
            y_pb: np.ndarray (n_samples, n_cell_types) float32 proportions
            cell_type_names: list of cell type names
    """
    import time
    t0 = time.time()
    rng = np.random.RandomState(random_state)

    # Get cell type info
    labels = adata.obs[label_column].values.astype(str)
    cell_types = sorted(list(set(labels)))
    n_types = len(cell_types)

    # Pre-index cells by type
    type_indices = {}
    for ct in cell_types:
        type_indices[ct] = np.where(labels == ct)[0]

    # ── Pre-convert expression matrix to dense + subset genes (ONCE) ──
    # This avoids repeated sparse→dense conversion inside the loop
    if var_genes is not None:
        var_genes = list(var_genes)
        # Use dict for O(1) lookup instead of O(n) list.index()
        gene_name_to_idx = {g: i for i, g in enumerate(adata.var_names)}
        shared = [g for g in var_genes if g in gene_name_to_idx]
        if len(shared) < len(var_genes):
            logger.info(f"  {len(var_genes) - len(shared)} var_genes missing "
                        f"from adata, using {len(shared)} shared genes")
        gene_idx = np.array([gene_name_to_idx[g] for g in shared])
        n_genes = len(shared)
        actual_var_genes = shared
    else:
        gene_idx = None
        n_genes = adata.n_vars
        actual_var_genes = list(adata.var_names)

    # One-time dense conversion: sparse→dense + gene subset
    # Following TAPE's approach: convert once, use C-contiguous array
    logger.info(f"  Pre-converting expression matrix to dense ...")
    if sp.issparse(adata.X):
        if gene_idx is not None:
            # Subset columns first (sparse), then convert to dense
            X_dense = adata.X[:, gene_idx].toarray()
        else:
            X_dense = adata.X.toarray()
    else:
        if gene_idx is not None:
            X_dense = np.array(adata.X[:, gene_idx])
        else:
            X_dense = np.array(adata.X)
    # C-contiguous float64 for fast row indexing (TAPE uses float32,
    # we use float64 for summation precision, cast to float32 after CPM)
    X_dense = np.ascontiguousarray(X_dense, dtype=np.float64)
    logger.info(f"  Dense matrix: {X_dense.shape} ({X_dense.nbytes / 1e6:.1f} MB)")

    logger.info(f"  Generating {n_samples} pseudo-bulk samples "
                f"({n_cells_per_sample} cells/sample, {n_types} types, "
                f"{n_genes} genes) ...")

    # ── Batch Dirichlet sampling ──
    all_proportions = rng.dirichlet(np.ones(n_types), size=n_samples)
    all_cell_counts = np.round(all_proportions * n_cells_per_sample).astype(int)
    # Fix rounding to ensure each row sums to n_cells_per_sample
    diffs = n_cells_per_sample - all_cell_counts.sum(axis=1)
    max_idx = np.argmax(all_cell_counts, axis=1)
    all_cell_counts[np.arange(n_samples), max_idx] += diffs

    # Actual proportions
    y_pb = all_cell_counts.astype(np.float32) / n_cells_per_sample

    # ── Generate pseudo-bulk samples ──
    X_pb = np.zeros((n_samples, n_genes), dtype=np.float32)

    for i in range(n_samples):
        # Collect all cell indices for this sample (no inner cell-type loop)
        all_idx = []
        for j, ct in enumerate(cell_types):
            n_pick = all_cell_counts[i, j]
            if n_pick > 0:
                all_idx.append(rng.choice(type_indices[ct], size=n_pick,
                                          replace=True))

        # Single dense indexing + sum (instead of per-type sparse→dense)
        all_idx = np.concatenate(all_idx)
        bulk_expr = X_dense[all_idx].sum(axis=0)

        # Normalize to log1p(CPM)
        total = bulk_expr.sum()
        if total > 0:
            bulk_cpm = bulk_expr / total * target_sum
        else:
            bulk_cpm = bulk_expr
        X_pb[i] = np.log1p(bulk_cpm).astype(np.float32)

        if (i + 1) % 1000 == 0:
            logger.info(f"    {i + 1}/{n_samples} samples generated")

    elapsed = time.time() - t0
    logger.info(f"  Done: X_pb {X_pb.shape}, y_pb {y_pb.shape} "
                f"({elapsed:.1f}s)")
    return {
        'X_pb': X_pb,
        'y_pb': y_pb,
        'cell_type_names': cell_types,
        'var_genes': actual_var_genes
    }


# =============================================================================
# DeconvModel
# =============================================================================

class DeconvModel(nn.Module):
    """Deep neural network for bulk RNA-seq deconvolution.

    Architecture (reuses building blocks from MultiTaskModel):

        Input: log-CPM expression (n_samples x n_genes)
          |
        InputBlock:     BatchNorm -> Dropout -> Linear -> SELU
          |
        ResNetBlock x4: Dropout -> Linear -> SELU + skip connection
          |
        ResNetLastBlock: Dropout -> Linear -> SELU (no skip)
          |
        prop_head (TaskHead): Linear -> SELU -> Dropout -> Linear
          |
        Softmax -> proportions (sum = 1)
    """

    def __init__(self,
                 input_size: int,
                 num_cell_types: int,
                 hidden_size: int = 512,
                 head_hidden_size: int = 256,
                 n_resnet_blocks: int = 4,
                 input_dropout_rate: float = 0.25,
                 resnet_dropout_rate: float = 0.1,
                 head_dropout_rate: float = 0.1,
                 batch_size: int = 256,
                 test_size: float = 0.1,
                 num_epochs: int = 50,
                 early_stopping_patience: int = 10,
                 learning_rate: float = 1e-3,
                 loss_function: str = 'mse',
                 use_lr_scheduler: bool = True,
                 random_state: int = 42,
                 recon_weight: float = 0.1):
        super().__init__()

        # Store config
        self.input_size = input_size
        self.num_cell_types = num_cell_types
        self.hidden_size = hidden_size
        self.head_hidden_size = head_hidden_size
        self.n_resnet_blocks = n_resnet_blocks
        self.input_dropout_rate = input_dropout_rate
        self.resnet_dropout_rate = resnet_dropout_rate
        self.head_dropout_rate = head_dropout_rate
        self.batch_size = batch_size
        self.test_size = test_size
        self.num_epochs = num_epochs
        self.early_stopping_patience = early_stopping_patience
        self.learning_rate = learning_rate
        self.loss_function = loss_function
        self.use_lr_scheduler = use_lr_scheduler
        self.random_state = random_state
        self.recon_weight = recon_weight

        # Build network
        self.input_block = InputBlock(input_size, hidden_size, input_dropout_rate)
        self.resnet_blocks = nn.Sequential(
            *[ResNetBlock(hidden_size, resnet_dropout_rate)
              for _ in range(n_resnet_blocks)]
        )
        self.resnet_last = ResNetLastBlock(hidden_size, resnet_dropout_rate)
        self.prop_head = TaskHead(hidden_size, head_hidden_size,
                                  num_cell_types, head_dropout_rate)

        # TAPE-style decoder: 5-layer Linear (no bias, no activation)
        # Weight product = cell-type-specific GEP (K × n_genes)
        self.decoder = nn.Sequential(
            nn.Linear(num_cell_types, 64, bias=False),
            nn.Linear(64, 128, bias=False),
            nn.Linear(128, 256, bias=False),
            nn.Linear(256, hidden_size, bias=False),
            nn.Linear(hidden_size, input_size, bias=False)
        )

        # Device
        self.device = torch.device('cuda' if torch.cuda.is_available()
                                   else 'cpu')
        self.to(self.device)

        # State
        self.is_trained = False
        self.var_genes = None
        self.cell_type_names = None
        self.target_sum = 1e6
        self.history = {}
        self.best_val_loss = float('inf')
        self.best_val_corr = 0.0
        self.metadata = {}

    def sigmatrix(self) -> torch.Tensor:
        """Compute cell-type-specific GEP via decoder weight product.

        Returns:
            Tensor of shape (num_cell_types, input_size), non-negative.
            Each row = gene expression profile for one cell type.
        """
        W = self.decoder[0].weight.T  # (K, 64)
        for layer in list(self.decoder)[1:]:
            W = W @ layer.weight.T
        return F.relu(W)  # (K, n_genes)

    def forward(self, x: torch.Tensor,
                mode: str = 'train') -> Dict[str, torch.Tensor]:
        """Forward pass with train/test mode support.

        Args:
            x: (batch, n_genes) log-CPM expression.
            mode: 'train' uses softmax; 'test' uses ReLU + L1 normalization
                  (TAPE style) for biologically interpretable proportions.
        """
        x_enc = self.input_block(x)
        x_enc = self.resnet_blocks(x_enc)
        x_enc = self.resnet_last(x_enc)
        embedding = x_enc
        logits = self.prop_head(x_enc)

        # Proportion constraint
        if mode == 'train':
            proportions = F.softmax(logits, dim=-1)
        else:
            proportions = F.relu(logits)
            proportions = proportions / (proportions.sum(dim=1, keepdim=True) + 1e-8)

        # Decoder reconstruction
        sigm = self.sigmatrix()
        x_recon = proportions @ sigm

        return {
            'logits': logits,
            'proportions': proportions,
            'embedding': embedding,
            'x_recon': x_recon,
            'sigmatrix': sigm
        }

    # -----------------------------------------------------------------
    # Loss computation
    # -----------------------------------------------------------------

    def _compute_loss(self, outputs, x_input, true_proportions):
        """Compute combined loss: proportion_loss + recon_weight * recon_loss.

        Proportion loss uses the configured loss_function (mse/kl/l1).
        Reconstruction loss always uses L1.

        Args:
            outputs: dict from forward() with 'proportions' and 'x_recon'.
            x_input: original input expression (for reconstruction target).
            true_proportions: ground-truth proportions.

        Returns:
            (total_loss, prop_loss, recon_loss) tuple.
        """
        # Proportion loss (configurable)
        if self.loss_function == 'l1':
            prop_loss = F.l1_loss(outputs['proportions'], true_proportions)
        elif self.loss_function == 'kl':
            pred_log = torch.log(outputs['proportions'].clamp(min=1e-8))
            true_clamped = true_proportions.clamp(min=1e-8)
            prop_loss = F.kl_div(pred_log, true_clamped,
                                 reduction='batchmean', log_target=False)
        else:  # 'mse' (default)
            prop_loss = F.mse_loss(outputs['proportions'], true_proportions)

        # Reconstruction loss (always L1)
        recon_loss = F.l1_loss(outputs['x_recon'], x_input)

        # Combined with configurable weight
        total_loss = prop_loss + self.recon_weight * recon_loss
        return total_loss, prop_loss, recon_loss

    # -----------------------------------------------------------------
    # Training
    # -----------------------------------------------------------------

    def fit(self, X_pb: np.ndarray, y_pb: np.ndarray,
            cell_type_names: List[str]):
        """Train the deconvolution model on pseudo-bulk data.

        Args:
            X_pb: (n_samples, n_genes) log-CPM pseudo-bulk expression
            y_pb: (n_samples, n_cell_types) true proportions
            cell_type_names: ordered list of cell type names
        """
        self.cell_type_names = list(cell_type_names)
        torch.manual_seed(self.random_state)
        np.random.seed(self.random_state)

        # Train / val split
        train_idx, val_idx = train_test_split(
            np.arange(X_pb.shape[0]),
            test_size=self.test_size,
            random_state=self.random_state
        )

        train_ds = DeconvDataset(X_pb[train_idx], y_pb[train_idx])
        val_ds = DeconvDataset(X_pb[val_idx], y_pb[val_idx])
        train_loader = DataLoader(train_ds, batch_size=self.batch_size,
                                   shuffle=True, drop_last=False)
        val_loader = DataLoader(val_ds, batch_size=self.batch_size,
                                 shuffle=False, drop_last=False)

        # Optimizer + scheduler
        optimizer = torch.optim.Adam(self.parameters(), lr=self.learning_rate)
        scheduler = None
        if self.use_lr_scheduler:
            from torch.optim.lr_scheduler import CosineAnnealingWarmRestarts
            scheduler = CosineAnnealingWarmRestarts(
                optimizer, T_0=10, T_mult=2, eta_min=1e-6)

        # History tracking
        history = {
            'train_loss': [], 'val_loss': [],
            'train_corr': [], 'val_corr': [],
            'train_recon_loss': [], 'val_recon_loss': []
        }

        best_val_loss = float('inf')
        best_state = None
        patience_counter = 0

        logger.info(f"  Training: {X_pb.shape[0]} samples "
                    f"(train={len(train_idx)}, val={len(val_idx)}), "
                    f"{self.num_epochs} epochs, loss={self.loss_function}")

        for epoch in range(self.num_epochs):
            # ── Train ──
            self.train()
            train_losses = []
            all_pred_train = []
            all_true_train = []

            train_recon_losses = []
            for batch in train_loader:
                X_batch = batch['X'].to(self.device)
                y_batch = batch['y'].to(self.device)

                outputs = self(X_batch, mode='train')
                loss, prop_loss, recon_loss = self._compute_loss(
                    outputs, X_batch, y_batch)

                optimizer.zero_grad()
                loss.backward()
                torch.nn.utils.clip_grad_norm_(self.parameters(), max_norm=1.0)
                optimizer.step()

                train_losses.append(prop_loss.item())
                train_recon_losses.append(recon_loss.item())
                all_pred_train.append(outputs['proportions'].detach().cpu().numpy())
                all_true_train.append(y_batch.cpu().numpy())

            if scheduler is not None:
                scheduler.step(epoch)

            # ── Validate ──
            self.eval()
            val_losses = []
            all_pred_val = []
            all_true_val = []

            val_recon_losses = []
            with torch.no_grad():
                for batch in val_loader:
                    X_batch = batch['X'].to(self.device)
                    y_batch = batch['y'].to(self.device)

                    outputs = self(X_batch, mode='train')
                    loss, prop_loss, recon_loss = self._compute_loss(
                        outputs, X_batch, y_batch)

                    val_losses.append(prop_loss.item())
                    val_recon_losses.append(recon_loss.item())
                    all_pred_val.append(outputs['proportions'].cpu().numpy())
                    all_true_val.append(y_batch.cpu().numpy())

            # Compute metrics
            train_loss = np.mean(train_losses)
            val_loss = np.mean(val_losses)
            train_recon = np.mean(train_recon_losses)
            val_recon = np.mean(val_recon_losses)

            pred_train = np.vstack(all_pred_train)
            true_train = np.vstack(all_true_train)
            pred_val = np.vstack(all_pred_val)
            true_val = np.vstack(all_true_val)

            train_corr = self._mean_pearson(pred_train, true_train)
            val_corr = self._mean_pearson(pred_val, true_val)

            history['train_loss'].append(float(train_loss))
            history['val_loss'].append(float(val_loss))
            history['train_corr'].append(float(train_corr))
            history['val_corr'].append(float(val_corr))
            history['train_recon_loss'].append(float(train_recon))
            history['val_recon_loss'].append(float(val_recon))

            # Logging
            lr_str = f"{optimizer.param_groups[0]['lr']:.2e}"
            logger.info(
                f"  Epoch {epoch + 1:3d}/{self.num_epochs} | "
                f"prop={train_loss:.5f}/{val_loss:.5f} "
                f"recon={train_recon:.5f}/{val_recon:.5f} | "
                f"r={train_corr:.4f}/{val_corr:.4f} | "
                f"lr={lr_str}"
            )

            # Early stopping (monitor val_loss)
            if val_loss < best_val_loss:
                best_val_loss = val_loss
                best_state = copy.deepcopy(self.state_dict())
                patience_counter = 0
                self.best_val_loss = float(val_loss)
                self.best_val_corr = float(val_corr)
            else:
                patience_counter += 1
                if patience_counter >= self.early_stopping_patience:
                    logger.info(f"  Early stopping at epoch {epoch + 1} "
                                f"(patience={self.early_stopping_patience})")
                    break

        # Restore best weights
        if best_state is not None:
            self.load_state_dict(best_state)

        self.history = history
        self.is_trained = True
        logger.info(f"  Training complete: best_val_loss={self.best_val_loss:.5f}, "
                    f"best_val_corr={self.best_val_corr:.4f}")

    # -----------------------------------------------------------------
    # Adaptive Stage (TAPE-style)
    # -----------------------------------------------------------------

    def adaptive_stage(self, X_bulk: np.ndarray,
                       max_iter: int = 3, step: int = 300,
                       lr: float = 1e-4) -> Dict[str, np.ndarray]:
        """TAPE adaptive stage: alternating decoder/encoder fine-tuning.

        Adapts the model to the target tissue by alternately optimizing:
        1. Decoder: minimize recon_loss + sigmatrix_drift
        2. Encoder: minimize proportion_drift + recon_loss

        Args:
            X_bulk: (n_samples, n_genes) log-CPM bulk expression.
            max_iter: Number of alternating rounds (default 3).
            step: Gradient steps per phase per round (default 300).
            lr: Learning rate for adaptive optimization (default 1e-4).

        Returns:
            dict with 'sigmatrix', 'proportions', 'embedding' numpy arrays.
        """
        X_t = torch.FloatTensor(X_bulk).to(self.device)

        # Save original predictions and sigmatrix
        self.eval()
        with torch.no_grad():
            ori = self(X_t, mode='test')
            ori_pred = ori['proportions'].detach()
            ori_sigm = ori['sigmatrix'].detach()

        # Separate decoder / encoder optimizers
        dec_params = [p for n, p in self.named_parameters() if 'decoder' in n]
        enc_params = [p for n, p in self.named_parameters() if 'decoder' not in n]
        optD = torch.optim.Adam(dec_params, lr=lr)
        optE = torch.optim.Adam(enc_params, lr=lr)

        for k in range(max_iter):
            self.train()
            # ── Decoder fine-tuning ──
            for _ in range(step):
                optD.zero_grad()
                out = self(X_t, mode='train')
                loss = (F.l1_loss(out['x_recon'], X_t)
                        + F.l1_loss(out['sigmatrix'], ori_sigm))
                loss.backward()
                optD.step()

            # ── Encoder fine-tuning ──
            for _ in range(step):
                optE.zero_grad()
                out = self(X_t, mode='train')
                loss = (F.l1_loss(ori_pred, out['proportions'])
                        + F.l1_loss(out['x_recon'], X_t))
                loss.backward()
                optE.step()

            logger.info(f"  Adaptive round {k+1}/{max_iter} done")

        # Final prediction in test mode
        self.eval()
        with torch.no_grad():
            final = self(X_t, mode='test')
        return {
            'sigmatrix': final['sigmatrix'].cpu().numpy(),
            'proportions': final['proportions'].cpu().numpy(),
            'embedding': final['embedding'].cpu().numpy()
        }

    # -----------------------------------------------------------------
    # Prediction
    # -----------------------------------------------------------------

    def predict(self, X: np.ndarray, adaptive: bool = False,
                mode: str = 'overall', max_iter: int = 3,
                step: int = 300, lr: float = 1e-4) -> Dict[str, np.ndarray]:
        """Predict cell type proportions from log-CPM expression.

        Args:
            X: (n_samples, n_genes) log-CPM expression matrix.
            adaptive: If True, run TAPE-style adaptive stage to refine
                      proportions and extract cell-type-specific GEP.
            mode: 'overall' (single GEP for all samples) or
                  'high-resolution' (per-sample GEP, slower).
            max_iter: Adaptive alternating rounds (default 3).
            step: Gradient steps per phase per round (default 300).
            lr: Adaptive learning rate (default 1e-4).

        Returns:
            dict with 'proportions', 'embedding', and optionally 'sigmatrix'.
        """
        if not adaptive:
            # ── Standard forward prediction (no adaptive) ──
            self.eval()
            X_tensor = torch.FloatTensor(X).to(self.device)
            with torch.no_grad():
                outputs = self(X_tensor, mode='test')
            return {
                'proportions': outputs['proportions'].cpu().numpy(),
                'embedding': outputs['embedding'].cpu().numpy()
            }
        else:
            # ── Adaptive prediction ──
            if mode == 'overall':
                logger.info("  Adaptive stage: overall mode")
                return self.adaptive_stage(X, max_iter, step, lr)
            elif mode == 'high-resolution':
                logger.info(f"  Adaptive stage: high-resolution "
                            f"({X.shape[0]} samples)")
                results_list = []
                original_state = copy.deepcopy(self.state_dict())
                for i in range(X.shape[0]):
                    self.load_state_dict(copy.deepcopy(original_state))
                    res = self.adaptive_stage(X[i:i+1], max_iter, step, lr)
                    results_list.append(res)
                    if (i + 1) % 10 == 0:
                        logger.info(f"    {i+1}/{X.shape[0]} samples done")
                self.load_state_dict(original_state)
                # Merge results
                return {
                    'sigmatrix': np.stack(
                        [r['sigmatrix'] for r in results_list]),  # (N, K, G)
                    'proportions': np.vstack(
                        [r['proportions'] for r in results_list]),
                    'embedding': np.vstack(
                        [r['embedding'] for r in results_list])
                }
            else:
                raise ValueError(f"Unknown adaptive mode: {mode}")

    # -----------------------------------------------------------------
    # Save / Load
    # -----------------------------------------------------------------

    def save(self, save_path: str, metadata: dict = None) -> None:
        """Save model and metadata to .pt file."""
        if not self.is_trained:
            raise ValueError("Model must be trained before saving")
        if not save_path.endswith('.pt'):
            save_path += '.pt'

        save_dict = {
            'model_type': 'DeconvModel',
            'model_state_dict': self.state_dict(),
            'model_config': {
                'input_size': self.input_size,
                'num_cell_types': self.num_cell_types,
                'hidden_size': self.hidden_size,
                'head_hidden_size': self.head_hidden_size,
                'n_resnet_blocks': self.n_resnet_blocks,
                'input_dropout_rate': self.input_dropout_rate,
                'resnet_dropout_rate': self.resnet_dropout_rate,
                'head_dropout_rate': self.head_dropout_rate,
                'batch_size': self.batch_size,
                'test_size': self.test_size,
                'num_epochs': self.num_epochs,
                'early_stopping_patience': self.early_stopping_patience,
                'learning_rate': self.learning_rate,
                'loss_function': self.loss_function,
                'use_lr_scheduler': self.use_lr_scheduler,
                'random_state': self.random_state,
                'recon_weight': self.recon_weight,
                'use_decoder': True,
            },
            'var_genes': self.var_genes,
            'cell_type_names': self.cell_type_names,
            'target_sum': self.target_sum,
            'history': self.history,
            'best_val_loss': self.best_val_loss,
            'best_val_corr': self.best_val_corr,
        }
        if metadata is not None:
            self.metadata.update(metadata)
        save_dict['metadata'] = self.metadata

        torch.save(save_dict, save_path)
        logger.info(f"DeconvModel saved to {save_path}")

    @classmethod
    def load(cls, load_path: str) -> 'DeconvModel':
        """Load a saved DeconvModel."""
        if not load_path.endswith('.pt'):
            load_path += '.pt'

        save_dict = torch.load(load_path, weights_only=False,
                               map_location=torch.device('cpu'))
        config = save_dict['model_config']

        model = cls(
            input_size=config['input_size'],
            num_cell_types=config['num_cell_types'],
            hidden_size=config.get('hidden_size', 512),
            head_hidden_size=config.get('head_hidden_size', 256),
            n_resnet_blocks=config.get('n_resnet_blocks', 4),
            input_dropout_rate=config.get('input_dropout_rate', 0.25),
            resnet_dropout_rate=config.get('resnet_dropout_rate', 0.1),
            head_dropout_rate=config.get('head_dropout_rate', 0.1),
            batch_size=config.get('batch_size', 256),
            test_size=config.get('test_size', 0.1),
            num_epochs=config.get('num_epochs', 50),
            early_stopping_patience=config.get('early_stopping_patience', 10),
            learning_rate=config.get('learning_rate', 1e-3),
            loss_function=config.get('loss_function', 'mse'),
            use_lr_scheduler=config.get('use_lr_scheduler', True),
            random_state=config.get('random_state', 42),
            recon_weight=config.get('recon_weight', 0.1),
        )

        # Handle backward compatibility: old models without decoder
        has_decoder = config.get('use_decoder', False)
        state_dict = save_dict['model_state_dict']

        if not has_decoder:
            # Old model without decoder weights → load encoder only
            logger.warning("  Loading old model without decoder. "
                           "Decoder initialized with random weights. "
                           "Re-train for TAPE-style features.")
            # Filter out decoder keys from expected state
            model_keys = set(model.state_dict().keys())
            saved_keys = set(state_dict.keys())
            # Only load matching keys
            compatible_state = {k: v for k, v in state_dict.items()
                                if k in model_keys}
            model.load_state_dict(compatible_state, strict=False)
        else:
            model.load_state_dict(state_dict)

        model.to(model.device)

        model.var_genes = save_dict.get('var_genes', None)
        model.cell_type_names = save_dict.get('cell_type_names', None)
        model.target_sum = save_dict.get('target_sum', 1e6)
        model.history = save_dict.get('history', {})
        model.best_val_loss = save_dict.get('best_val_loss', float('inf'))
        model.best_val_corr = save_dict.get('best_val_corr', 0.0)
        model.metadata = save_dict.get('metadata', {})
        model.is_trained = True

        logger.info(f"DeconvModel loaded from {load_path} "
                    f"(device: {model.device}, "
                    f"{model.num_cell_types} cell types)")
        return model

    # -----------------------------------------------------------------
    # Utilities
    # -----------------------------------------------------------------

    @staticmethod
    def _mean_pearson(pred: np.ndarray, true: np.ndarray) -> float:
        """Mean Pearson correlation across cell types (column-wise)."""
        n_types = pred.shape[1]
        corrs = []
        for j in range(n_types):
            p = pred[:, j]
            t = true[:, j]
            if np.std(p) < 1e-10 or np.std(t) < 1e-10:
                corrs.append(0.0)
            else:
                corrs.append(float(np.corrcoef(p, t)[0, 1]))
        return float(np.mean(corrs))

    def get_feature_names(self) -> List[str]:
        """Return the gene names used by the model."""
        return self.var_genes if self.var_genes is not None else []

    def print_summary(self) -> None:
        """Print model summary."""
        print("=" * 60)
        print("DeconvModel Summary")
        print("=" * 60)
        print(f"  Input size     : {self.input_size}")
        print(f"  Cell types     : {self.num_cell_types}")
        if self.cell_type_names:
            print(f"  Types          : {', '.join(self.cell_type_names)}")
        print(f"  Hidden size    : {self.hidden_size}")
        print(f"  ResNet blocks  : {self.n_resnet_blocks}")
        print(f"  Loss function  : {self.loss_function}")
        print(f"  Decoder        : TAPE-style (5 layers, bias=False)")
        print(f"  Trained        : {self.is_trained}")
        if self.is_trained:
            print(f"  Best val loss  : {self.best_val_loss:.5f}")
            print(f"  Best val corr  : {self.best_val_corr:.4f}")
        n_params = sum(p.numel() for p in self.parameters())
        n_dec = sum(p.numel() for p in self.decoder.parameters())
        print(f"  Parameters     : {n_params:,} (decoder: {n_dec:,})")
        print(f"  Device         : {self.device}")
        print("=" * 60)
