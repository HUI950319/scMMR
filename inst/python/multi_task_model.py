"""
MultiTaskModel: Unified cell type classification + embedding regression.

Combines both tasks into a single model with:
- ResNet backbone with skip connections (adapted from scMMT)
- GradNorm for dynamic multi-task loss balancing
- Shared embedding space (512-dim)
- Label smoothing for classification (from scMMT)
- Task-specific intermediate layers before output heads
- Configurable early stopping monitor (cls / reg / combined)
- Cosine annealing learning rate scheduler

v3 improvements:
- Temperature scaling for calibrated classification probabilities
- OOD (out-of-distribution) detection via confidence & is_ood flags
- Aleatoric uncertainty for regression (NLL loss, outputs mean + variance)
- Transfer learning: freeze_backbone / fine_tune for new atlases
- Integrated Gradients for gene-level interpretability (SHAP-like)

References:
- scMMT: Multi-modal transfer learning for single-cell data
- GradNorm: Chen et al., 2018 (arXiv:1711.02257)
- SuperCT: Xie et al., 2019 (NAR)
- Temperature scaling: Guo et al., 2017 (ICML)
- Aleatoric uncertainty: Kendall & Gal, 2017 (NeurIPS)
"""

import scanpy as sc
import torch
import torch.nn as nn
import torch.nn.functional as F
import numpy as np
from torch.utils.data import Dataset, DataLoader
from sklearn.preprocessing import LabelEncoder
from sklearn.model_selection import train_test_split
import scipy.sparse as sp
from typing import Tuple, List, Optional, Union, Dict
import logging
import copy

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)


# =============================================================================
# Layer definitions (adapted from scMMT Network/Layers.py)
# =============================================================================

class InputBlock(nn.Module):
    """Input block: BatchNorm -> Dropout -> Linear -> SELU."""

    def __init__(self, in_features: int, out_features: int, dropout_rate: float = 0.25):
        super().__init__()
        self.bnorm = nn.BatchNorm1d(in_features)
        self.dropout = nn.Dropout(dropout_rate)
        self.linear = nn.Linear(in_features, out_features)
        self.act = nn.SELU()

    def forward(self, x: torch.Tensor) -> torch.Tensor:
        x = self.bnorm(x)
        x = self.dropout(x)
        x = self.linear(x)
        return self.act(x)


class ResNetBlock(nn.Module):
    """Residual block: Dropout -> Linear -> SELU + skip connection."""

    def __init__(self, hidden_size: int, dropout_rate: float = 0.1):
        super().__init__()
        self.dropout = nn.Dropout(dropout_rate)
        self.linear = nn.Linear(hidden_size, hidden_size)
        self.act = nn.SELU()

    def forward(self, x: torch.Tensor) -> torch.Tensor:
        h = self.dropout(x)
        h = self.linear(h)
        h = self.act(h)
        return x + h  # skip connection


class ResNetLastBlock(nn.Module):
    """Final residual block: Dropout -> Linear -> SELU (no skip connection)."""

    def __init__(self, hidden_size: int, dropout_rate: float = 0.1):
        super().__init__()
        self.dropout = nn.Dropout(dropout_rate)
        self.linear = nn.Linear(hidden_size, hidden_size)
        self.act = nn.SELU()

    def forward(self, x: torch.Tensor) -> torch.Tensor:
        x = self.dropout(x)
        x = self.linear(x)
        return self.act(x)


class TaskHead(nn.Module):
    """Task-specific head: Linear -> SELU -> Dropout -> Linear.

    Adds a task-specific intermediate layer before the final output,
    giving each task its own representation space to reduce task interference.
    """

    def __init__(self, in_features: int, hidden_features: int, out_features: int,
                 dropout_rate: float = 0.1):
        super().__init__()
        self.fc1 = nn.Linear(in_features, hidden_features)
        self.act = nn.SELU()
        self.dropout = nn.Dropout(dropout_rate)
        self.fc2 = nn.Linear(hidden_features, out_features)

    def forward(self, x: torch.Tensor) -> torch.Tensor:
        x = self.fc1(x)
        x = self.act(x)
        x = self.dropout(x)
        return self.fc2(x)


# =============================================================================
# Dataset
# =============================================================================

class MultiTaskDataset(Dataset):
    """Dataset for multi-task learning with optional classification and regression targets."""

    def __init__(self, X: Union[np.ndarray, sp.spmatrix],
                 y_cls: Optional[np.ndarray] = None,
                 y_reg: Optional[np.ndarray] = None):
        if sp.issparse(X):
            self.X = torch.FloatTensor(X.toarray())
        else:
            self.X = torch.FloatTensor(X)
        self.y_cls = torch.LongTensor(y_cls) if y_cls is not None else None
        self.y_reg = torch.FloatTensor(y_reg) if y_reg is not None else None

    def __len__(self) -> int:
        return self.X.shape[0]

    def __getitem__(self, idx: int) -> Dict[str, torch.Tensor]:
        item = {'X': self.X[idx]}
        if self.y_cls is not None:
            item['y_cls'] = self.y_cls[idx]
        if self.y_reg is not None:
            item['y_reg'] = self.y_reg[idx]
        return item


# =============================================================================
# MultiTaskModel
# =============================================================================

class MultiTaskModel(nn.Module):
    """Unified multi-task model for cell type classification and embedding regression.

    Architecture:
        Input (binary matrix)
          -> InputBlock (BatchNorm -> Dropout -> Linear -> SELU)
          -> ResNetBlock x N (with skip connections)
          -> ResNetLastBlock (no skip)
          -> shared embedding (hidden_size dim)
          -> cls_head: TaskHead(512, 256, num_classes) [optional]
          -> reg_head: TaskHead(512, 256, embedding_dim) [optional]

    Key improvements over v1:
        1. Label smoothing (from scMMT) for better generalization
        2. Task-specific intermediate layers reduce task interference
        3. Configurable early_stopping_monitor: 'cls', 'reg', or 'combined'
        4. Cosine annealing LR scheduler for better convergence
        5. Fixed GradNorm L1Loss shape warning
    """

    def __init__(self, input_size: int,
                 num_classes: Optional[int] = None,
                 embedding_dim: int = 2,
                 hidden_size: int = 512,
                 head_hidden_size: int = 256,
                 n_resnet_blocks: int = 4,
                 input_dropout_rate: float = 0.25,
                 resnet_dropout_rate: float = 0.1,
                 head_dropout_rate: float = 0.1,
                 batch_size: int = 128,
                 test_size: float = 0.2,
                 num_epochs: int = 100,
                 early_stopping_patience: int = 10,
                 early_stopping_monitor: str = 'cls',
                 learning_rate: float = 1e-3,
                 use_gradnorm: bool = True,
                 gradnorm_alpha: float = 0.15,
                 label_smoothing: float = 0.1,
                 use_lr_scheduler: bool = True,
                 use_uncertainty: bool = False,
                 random_state: int = 42):
        """
        Args:
            input_size: Number of input genes.
            num_classes: Number of cell type classes. None = no classification head.
            embedding_dim: Dimension of embedding output (default 2 for UMAP).
            hidden_size: Width of shared backbone (default 512).
            head_hidden_size: Width of task-specific intermediate layer (default 256).
            n_resnet_blocks: Number of ResNet blocks (default 4).
            input_dropout_rate: Dropout for input block (default 0.25).
            resnet_dropout_rate: Dropout for ResNet blocks (default 0.1).
            head_dropout_rate: Dropout for task heads (default 0.1).
            batch_size: Training batch size.
            test_size: Validation split ratio.
            num_epochs: Maximum training epochs.
            early_stopping_patience: Epochs to wait before stopping.
            early_stopping_monitor: Which metric to monitor:
                'cls' = val accuracy only (recommended when classification is primary)
                'reg' = val RMSE only
                'combined' = weighted combination of both
            learning_rate: Initial learning rate for Adam.
            use_gradnorm: Enable GradNorm for multi-task loss balancing.
            gradnorm_alpha: GradNorm asymmetry parameter (0.15 = moderate).
            label_smoothing: Label smoothing for CrossEntropyLoss (0.1 from scMMT).
            use_lr_scheduler: Enable cosine annealing LR scheduler.
            use_uncertainty: Enable aleatoric uncertainty for regression head.
                When True, regression outputs mean + log_variance (NLL loss).
            random_state: Random seed.
        """
        super().__init__()

        # -- Shared backbone --
        self.input_block = InputBlock(input_size, hidden_size, input_dropout_rate)
        self.resnet_blocks = nn.Sequential(
            *[ResNetBlock(hidden_size, resnet_dropout_rate) for _ in range(n_resnet_blocks)]
        )
        self.resnet_last = ResNetLastBlock(hidden_size, resnet_dropout_rate)

        # -- Task heads --
        self.has_cls = num_classes is not None
        self.has_reg = True
        self.use_uncertainty = use_uncertainty

        # Regression output dim: 2x if uncertainty (mean + log_var)
        reg_out_dim = embedding_dim * 2 if use_uncertainty else embedding_dim

        if head_hidden_size > 0:
            # Task-specific intermediate layers (v2)
            if self.has_cls:
                self.cls_head = TaskHead(hidden_size, head_hidden_size, num_classes, head_dropout_rate)
            self.reg_head = TaskHead(hidden_size, head_hidden_size, reg_out_dim, head_dropout_rate)
        else:
            # Direct linear heads (v1 behavior)
            if self.has_cls:
                self.cls_head = nn.Linear(hidden_size, num_classes)
            self.reg_head = nn.Linear(hidden_size, reg_out_dim)

        # -- Architecture config --
        self.input_size = input_size
        self.num_classes = num_classes
        self.embedding_dim = embedding_dim
        self.hidden_size = hidden_size
        self.head_hidden_size = head_hidden_size
        self.n_resnet_blocks = n_resnet_blocks
        self.input_dropout_rate = input_dropout_rate
        self.resnet_dropout_rate = resnet_dropout_rate
        self.head_dropout_rate = head_dropout_rate

        # -- Training config --
        self.batch_size = batch_size
        self.test_size = test_size
        self.num_epochs = num_epochs
        self.early_stopping_patience = early_stopping_patience
        self.early_stopping_monitor = early_stopping_monitor
        self.learning_rate = learning_rate
        self.use_gradnorm = use_gradnorm
        self.gradnorm_alpha = gradnorm_alpha
        self.label_smoothing = label_smoothing
        self.use_lr_scheduler = use_lr_scheduler
        self.random_state = random_state

        # -- v3: Temperature scaling & OOD --
        self.temperature = nn.Parameter(torch.ones(1), requires_grad=False)  # calibrated post-training
        self.ood_threshold = 0.5  # default; user can tune after calibration

        # -- State --
        self.label_encoder = None
        self.var_genes = None
        self.embedding_key = None
        self._reg_mean = None   # UMAP z-score mean (saved for inverse-transform)
        self._reg_std = None    # UMAP z-score std
        self._ref_emb = None    # Reference shared embeddings for kNN UMAP projection
        self._ref_umap = None   # Reference original UMAP coordinates
        self._y_reg_raw = None  # Temporary: raw UMAP coords during training
        self.metadata = {}
        self.device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
        self.is_trained = False

        # -- Training history --
        self.history = {
            'train_loss': [], 'val_loss': [],
            'train_cls_loss': [], 'val_cls_loss': [],
            'train_reg_loss': [], 'val_reg_loss': [],
            'train_acc': [], 'val_acc': [],
            'train_rmse': [], 'val_rmse': [],
            'w_cls': [], 'w_reg': [],
            'lr': [],
        }
        self.best_val_acc = 0.0
        self.best_val_rmse = float('inf')

    def forward(self, x: torch.Tensor) -> Dict[str, torch.Tensor]:
        """Forward pass returning shared embedding and task-specific outputs.

        When use_uncertainty=True, regression output is split into:
            'reg_mean': predicted coordinates (embedding_dim)
            'reg_log_var': log variance for each coordinate (embedding_dim)
            'reg': same as 'reg_mean' (for backward compatibility)
        """
        x = self.input_block(x)
        x = self.resnet_blocks(x)
        x = self.resnet_last(x)
        result = {'embedding': x}
        if self.has_cls:
            result['cls'] = self.cls_head(x)
        if self.has_reg:
            reg_out = self.reg_head(x)
            if self.use_uncertainty:
                mu, log_var = reg_out.chunk(2, dim=-1)
                result['reg'] = mu
                result['reg_mean'] = mu
                result['reg_log_var'] = log_var
            else:
                result['reg'] = reg_out
        return result

    # =========================================================================
    # Data preprocessing
    # =========================================================================

    def preprocess_data(self, adata: sc.AnnData,
                        label_column: Optional[str] = None,
                        embedding_key: Optional[str] = None,
                        var_genes: Optional[List[str]] = None
                        ) -> Tuple[np.ndarray, Optional[np.ndarray], Optional[np.ndarray]]:
        """Preprocess AnnData for multi-task training."""
        if var_genes is not None:
            adata = adata[:, list(var_genes)].copy()
            self.var_genes = list(var_genes)
        else:
            self.var_genes = list(adata.var_names)

        logger.info("Converting counts to binary matrix...")
        X = (adata.X > 0).astype(np.float32)

        y_cls = None
        if label_column is not None:
            self.label_encoder = LabelEncoder()
            y_cls = self.label_encoder.fit_transform(adata.obs[label_column].values)

        y_reg = None
        if embedding_key is not None:
            if embedding_key not in adata.obsm:
                raise ValueError(
                    f"'{embedding_key}' not in adata.obsm. Available: {list(adata.obsm.keys())}"
                )
            y_reg_raw = adata.obsm[embedding_key].astype(np.float32)
            self.embedding_key = embedding_key
            self._y_reg_raw = y_reg_raw   # keep raw UMAP for kNN reference index

            # Z-score normalize UMAP coordinates for stable regression
            self._reg_mean = y_reg_raw.mean(axis=0)
            self._reg_std  = y_reg_raw.std(axis=0).clip(min=1e-6)
            y_reg = (y_reg_raw - self._reg_mean) / self._reg_std
            logger.info(f"UMAP z-score: mean={self._reg_mean}, std={self._reg_std}")

        return X, y_cls, y_reg

    # =========================================================================
    # Data loaders
    # =========================================================================

    def _prepare_data_loaders(self, X, y_cls, y_reg) -> Tuple[DataLoader, DataLoader]:
        """Split data and create train/val DataLoaders."""
        indices = np.arange(X.shape[0] if not sp.issparse(X) else X.shape[0])
        stratify = y_cls if y_cls is not None else None

        train_idx, val_idx = train_test_split(
            indices, test_size=self.test_size,
            random_state=self.random_state, stratify=stratify
        )

        def _index(arr, idx):
            return arr[idx] if arr is not None else None

        train_ds = MultiTaskDataset(_index(X, train_idx), _index(y_cls, train_idx), _index(y_reg, train_idx))
        val_ds = MultiTaskDataset(_index(X, val_idx), _index(y_cls, val_idx), _index(y_reg, val_idx))

        train_loader = DataLoader(train_ds, batch_size=self.batch_size, shuffle=True)
        val_loader = DataLoader(val_ds, batch_size=self.batch_size, shuffle=False)
        return train_loader, val_loader

    # =========================================================================
    # Training: multi-task with GradNorm
    # =========================================================================

    def _train_epoch_multitask(self, train_loader, cls_criterion, reg_criterion,
                                model_optimizer, gradnorm_optimizer,
                                w_cls, w_reg, l0_cls, l0_reg, alpha, epoch):
        """One epoch of GradNorm multi-task training."""
        self.train()
        total_loss = 0
        total_cls_loss = 0
        total_reg_loss = 0
        correct = 0
        total_samples = 0
        total_se = 0

        for batch in train_loader:
            batch_X = batch['X'].to(self.device)
            batch_y_cls = batch['y_cls'].to(self.device)
            batch_y_reg = batch['y_reg'].to(self.device)

            outputs = self(batch_X)

            raw_cls_loss = cls_criterion(outputs['cls'], batch_y_cls)
            # For uncertainty: NLL needs raw output (before split); for MSE: use mean
            if self.use_uncertainty:
                raw_reg_out = torch.cat([outputs['reg_mean'], outputs['reg_log_var']], dim=-1)
                raw_reg_loss = reg_criterion(raw_reg_out, batch_y_reg)
            else:
                raw_reg_loss = reg_criterion(outputs['reg'], batch_y_reg)

            weighted_cls_loss = w_cls * raw_cls_loss
            weighted_reg_loss = w_reg * raw_reg_loss
            loss = (weighted_cls_loss + weighted_reg_loss) / 2

            # Store initial losses at first batch of first epoch
            if epoch == 0 and l0_cls is None:
                l0_cls = weighted_cls_loss.data.clone().clamp(min=1e-8)
                l0_reg = weighted_reg_loss.data.clone().clamp(min=1e-8)

            # Backprop total loss
            model_optimizer.zero_grad()
            loss.backward(retain_graph=True)

            # GradNorm: compute gradient norms w.r.t. first shared parameter
            shared_param = next(self.input_block.parameters())
            G1R = torch.autograd.grad(
                weighted_cls_loss, shared_param,
                retain_graph=True, create_graph=True
            )
            G1 = torch.norm(G1R[0], 2)

            G2R = torch.autograd.grad(
                weighted_reg_loss, shared_param,
                retain_graph=True, create_graph=True
            )
            G2 = torch.norm(G2R[0], 2)
            G_avg = (G1 + G2) / 2

            # Relative loss ratios
            lhat_cls = weighted_cls_loss / l0_cls
            lhat_reg = weighted_reg_loss / l0_reg
            lhat_avg = (lhat_cls + lhat_reg) / 2

            inv_rate_cls = lhat_cls / lhat_avg
            inv_rate_reg = lhat_reg / lhat_avg

            # Target gradient norms
            C_cls = (G_avg * inv_rate_cls ** alpha).detach()
            C_reg = (G_avg * inv_rate_reg ** alpha).detach()

            # Update loss weights (fixed: ensure matching shapes for L1Loss)
            gradnorm_optimizer.zero_grad()
            grad_loss = torch.abs(G1 - C_cls) + torch.abs(G2 - C_reg)
            grad_loss.backward()
            gradnorm_optimizer.step()

            # Apply model gradients
            model_optimizer.step()

            # Renormalize weights to sum to 2
            with torch.no_grad():
                coef = 2.0 / (w_cls + w_reg)
                w_cls.data *= coef
                w_reg.data *= coef

            # Track metrics
            total_loss += loss.item() * batch_y_cls.size(0)
            total_cls_loss += raw_cls_loss.item() * batch_y_cls.size(0)
            total_reg_loss += raw_reg_loss.item() * batch_y_cls.size(0)

            _, predicted = outputs['cls'].max(1)
            correct += predicted.eq(batch_y_cls).sum().item()
            total_samples += batch_y_cls.size(0)

            se = ((outputs['reg'].detach() - batch_y_reg) ** 2).sum().item()
            total_se += se

        n = total_samples
        rmse = np.sqrt(total_se / (n * self.embedding_dim))

        return {
            'loss': total_loss / n,
            'cls_loss': total_cls_loss / n,
            'reg_loss': total_reg_loss / n,
            'accuracy': 100.0 * correct / n,
            'rmse': rmse,
            'l0_cls': l0_cls,
            'l0_reg': l0_reg,
            'w_cls': w_cls.item(),
            'w_reg': w_reg.item(),
        }

    def _train_epoch_single(self, train_loader, criterion, model_optimizer, task):
        """One epoch of single-task training (cls-only or reg-only)."""
        self.train()
        total_loss = 0
        correct = 0
        total_samples = 0
        total_se = 0

        for batch in train_loader:
            batch_X = batch['X'].to(self.device)
            outputs = self(batch_X)

            if task == 'cls':
                batch_y = batch['y_cls'].to(self.device)
                loss = criterion(outputs['cls'], batch_y)
                _, predicted = outputs['cls'].max(1)
                correct += predicted.eq(batch_y).sum().item()
            else:
                batch_y = batch['y_reg'].to(self.device)
                if self.use_uncertainty:
                    raw_reg_out = torch.cat([outputs['reg_mean'], outputs['reg_log_var']], dim=-1)
                    loss = criterion(raw_reg_out, batch_y)
                else:
                    loss = criterion(outputs['reg'], batch_y)
                se = ((outputs['reg'].detach() - batch_y) ** 2).sum().item()
                total_se += se

            model_optimizer.zero_grad()
            loss.backward()
            model_optimizer.step()

            total_loss += loss.item() * batch_X.size(0)
            total_samples += batch_X.size(0)

        n = total_samples
        result = {'loss': total_loss / n}
        if task == 'cls':
            result['accuracy'] = 100.0 * correct / n
            result['cls_loss'] = total_loss / n
        else:
            result['rmse'] = np.sqrt(total_se / (n * self.embedding_dim))
            result['reg_loss'] = total_loss / n
        return result

    # =========================================================================
    # Validation
    # =========================================================================

    def _validate(self, val_loader, cls_criterion, reg_criterion, do_cls, do_reg):
        """Validate on held-out data."""
        self.eval()
        total_loss = 0
        total_cls_loss = 0
        total_reg_loss = 0
        correct = 0
        total_samples = 0
        total_se = 0

        with torch.no_grad():
            for batch in val_loader:
                batch_X = batch['X'].to(self.device)
                outputs = self(batch_X)
                batch_loss = 0
                n = batch_X.size(0)

                if do_cls:
                    batch_y_cls = batch['y_cls'].to(self.device)
                    cls_loss = cls_criterion(outputs['cls'], batch_y_cls)
                    total_cls_loss += cls_loss.item() * n
                    batch_loss += cls_loss.item() * n
                    _, predicted = outputs['cls'].max(1)
                    correct += predicted.eq(batch_y_cls).sum().item()

                if do_reg:
                    batch_y_reg = batch['y_reg'].to(self.device)
                    if self.use_uncertainty:
                        raw_reg_out = torch.cat([outputs['reg_mean'], outputs['reg_log_var']], dim=-1)
                        reg_loss = reg_criterion(raw_reg_out, batch_y_reg)
                    else:
                        reg_loss = reg_criterion(outputs['reg'], batch_y_reg)
                    total_reg_loss += reg_loss.item() * n
                    batch_loss += reg_loss.item() * n
                    se = ((outputs['reg'] - batch_y_reg) ** 2).sum().item()
                    total_se += se

                total_loss += batch_loss
                total_samples += n

        n = total_samples
        result = {'loss': total_loss / n}
        if do_cls:
            result['cls_loss'] = total_cls_loss / n
            result['accuracy'] = 100.0 * correct / n
        if do_reg:
            result['reg_loss'] = total_reg_loss / n
            result['rmse'] = np.sqrt(total_se / (n * self.embedding_dim))
        return result

    # =========================================================================
    # Early stopping metric
    # =========================================================================

    def _compute_early_stopping_metric(self, val_metrics, do_cls, do_reg):
        """Compute early stopping metric based on monitor setting.

        Returns a scalar where HIGHER = BETTER.
        """
        monitor = self.early_stopping_monitor

        if monitor == 'cls' and do_cls:
            return val_metrics.get('accuracy', 0)
        elif monitor == 'reg' and do_reg:
            return -val_metrics.get('rmse', float('inf'))
        elif monitor == 'combined' and do_cls and do_reg:
            val_acc = val_metrics.get('accuracy', 0)
            val_rmse = val_metrics.get('rmse', float('inf'))
            return val_acc / 100.0 + 1.0 / (1.0 + val_rmse)
        else:
            # Fallback: use whatever is available
            if do_cls:
                return val_metrics.get('accuracy', 0)
            else:
                return -val_metrics.get('rmse', float('inf'))

    # =========================================================================
    # Main training entry point
    # =========================================================================

    def fit(self, adata: sc.AnnData,
            label_column: Optional[str] = None,
            embedding_key: Optional[str] = None,
            var_genes: Optional[List[str]] = None) -> None:
        """Train the multi-task model.

        Args:
            adata: AnnData with gene expression data.
            label_column: obs column for cell type labels (enables classification).
            embedding_key: obsm key for embedding coords (enables regression).
            var_genes: optional gene subset.
        """
        do_cls = label_column is not None
        do_reg = embedding_key is not None
        if not do_cls and not do_reg:
            raise ValueError("Must provide at least label_column or embedding_key")
        is_multitask = do_cls and do_reg

        # Preprocess
        X, y_cls, y_reg = self.preprocess_data(adata, label_column, embedding_key, var_genes)

        # Rebuild cls head if num_classes was not set at init
        if do_cls and not self.has_cls:
            n_classes = len(self.label_encoder.classes_)
            if self.head_hidden_size > 0:
                self.cls_head = TaskHead(self.hidden_size, self.head_hidden_size,
                                         n_classes, self.head_dropout_rate)
            else:
                self.cls_head = nn.Linear(self.hidden_size, n_classes)
            self.has_cls = True
            self.num_classes = n_classes

        # Rebuild reg head if uncertainty changed
        reg_out_dim = self.embedding_dim * 2 if self.use_uncertainty else self.embedding_dim
        if do_reg and hasattr(self, 'reg_head'):
            current_out = list(self.reg_head.parameters())[-1].shape[0]
            if current_out != reg_out_dim:
                if self.head_hidden_size > 0:
                    self.reg_head = TaskHead(self.hidden_size, self.head_hidden_size,
                                             reg_out_dim, self.head_dropout_rate)
                else:
                    self.reg_head = nn.Linear(self.hidden_size, reg_out_dim)

        # Prepare data
        train_loader, val_loader = self._prepare_data_loaders(X, y_cls, y_reg)

        # Setup criteria
        self.to(self.device)
        cls_criterion = nn.CrossEntropyLoss(
            label_smoothing=self.label_smoothing
        ) if do_cls else None

        # Regression loss: NLL if uncertainty, else MSE
        if do_reg and self.use_uncertainty:
            def _nll_loss(pred, target):
                """Negative log-likelihood with learned aleatoric uncertainty."""
                mu, log_var = pred.chunk(2, dim=-1)
                # NLL = 0.5 * (log_var + (y - mu)^2 / exp(log_var))
                precision = torch.exp(-log_var)
                return torch.mean(0.5 * (log_var + (target - mu) ** 2 * precision))
            reg_criterion = _nll_loss
        else:
            reg_criterion = nn.MSELoss() if do_reg else None

        model_optimizer = torch.optim.Adam(self.parameters(), lr=self.learning_rate)

        # LR scheduler
        scheduler = None
        if self.use_lr_scheduler:
            scheduler = torch.optim.lr_scheduler.CosineAnnealingWarmRestarts(
                model_optimizer, T_0=10, T_mult=2, eta_min=1e-6
            )

        # GradNorm weights
        w_cls = w_reg = None
        gradnorm_optimizer = None
        l0_cls = l0_reg = None
        if is_multitask and self.use_gradnorm:
            w_cls = torch.FloatTensor([1.0]).to(self.device).requires_grad_(True)
            w_reg = torch.FloatTensor([1.0]).to(self.device).requires_grad_(True)
            gradnorm_optimizer = torch.optim.Adam([w_cls, w_reg], lr=self.learning_rate)

        logger.info(f"Starting training on {self.device}")
        logger.info(f"Mode: {'multi-task (GradNorm)' if is_multitask and self.use_gradnorm else 'multi-task (fixed)' if is_multitask else 'single-task'}")
        logger.info(f"Input: {X.shape[1]} genes, Samples: {X.shape[0]}")
        logger.info(f"Early stopping monitor: {self.early_stopping_monitor}")
        logger.info(f"Label smoothing: {self.label_smoothing}")
        logger.info(f"LR scheduler: {'CosineAnnealingWarmRestarts' if self.use_lr_scheduler else 'None'}")
        logger.info(f"Uncertainty:  {'NLL (aleatoric)' if self.use_uncertainty else 'off (MSE)'}")
        if do_cls:
            logger.info(f"Classification: {len(self.label_encoder.classes_)} classes")
        if do_reg:
            logger.info(f"Regression: {y_reg.shape[1]}d embedding ({embedding_key})")

        # Training loop
        best_state = None
        best_metric = -float('inf')
        patience_counter = 0

        # Reset history
        for k in self.history:
            self.history[k] = []

        for epoch in range(self.num_epochs):
            # Train
            if is_multitask and self.use_gradnorm:
                metrics = self._train_epoch_multitask(
                    train_loader, cls_criterion, reg_criterion,
                    model_optimizer, gradnorm_optimizer,
                    w_cls, w_reg, l0_cls, l0_reg,
                    self.gradnorm_alpha, epoch
                )
                l0_cls = metrics['l0_cls']
                l0_reg = metrics['l0_reg']
                self.history['w_cls'].append(metrics['w_cls'])
                self.history['w_reg'].append(metrics['w_reg'])
            elif is_multitask:
                metrics = self._train_epoch_multitask_fixed(
                    train_loader, cls_criterion, reg_criterion, model_optimizer
                )
            else:
                task = 'cls' if do_cls else 'reg'
                criterion = cls_criterion if do_cls else reg_criterion
                metrics = self._train_epoch_single(train_loader, criterion, model_optimizer, task)

            # Step LR scheduler
            current_lr = model_optimizer.param_groups[0]['lr']
            self.history['lr'].append(current_lr)
            if scheduler is not None:
                scheduler.step()

            # Validate
            val_metrics = self._validate(val_loader, cls_criterion, reg_criterion, do_cls, do_reg)

            # Record history
            self.history['train_loss'].append(metrics.get('loss', 0))
            self.history['val_loss'].append(val_metrics.get('loss', 0))
            self.history['train_cls_loss'].append(metrics.get('cls_loss', 0))
            self.history['val_cls_loss'].append(val_metrics.get('cls_loss', 0))
            self.history['train_reg_loss'].append(metrics.get('reg_loss', 0))
            self.history['val_reg_loss'].append(val_metrics.get('reg_loss', 0))
            self.history['train_acc'].append(metrics.get('accuracy', 0))
            self.history['val_acc'].append(val_metrics.get('accuracy', 0))
            self.history['train_rmse'].append(metrics.get('rmse', 0))
            self.history['val_rmse'].append(val_metrics.get('rmse', 0))

            # Log
            msg = f"Epoch {epoch + 1}/{self.num_epochs} (lr={current_lr:.2e}): "
            if do_cls:
                msg += f"Train Acc: {metrics.get('accuracy', 0):.2f}%, Val Acc: {val_metrics.get('accuracy', 0):.2f}% | "
            if do_reg:
                msg += f"Train RMSE: {metrics.get('rmse', 0):.4f}, Val RMSE: {val_metrics.get('rmse', 0):.4f} | "
            if is_multitask and self.use_gradnorm:
                msg += f"w_cls: {metrics.get('w_cls', 1):.3f}, w_reg: {metrics.get('w_reg', 1):.3f}"
            logger.info(msg)

            # Early stopping
            current_metric = self._compute_early_stopping_metric(val_metrics, do_cls, do_reg)

            if current_metric > best_metric:
                best_metric = current_metric
                patience_counter = 0
                best_state = {k: v.clone() for k, v in self.state_dict().items()}
                if do_cls:
                    self.best_val_acc = val_metrics.get('accuracy', 0)
                if do_reg:
                    self.best_val_rmse = val_metrics.get('rmse', float('inf'))
            else:
                patience_counter += 1
                if patience_counter >= self.early_stopping_patience:
                    logger.info(f"Early stopping at epoch {epoch + 1}")
                    break

        # Restore best weights
        if best_state is not None:
            self.load_state_dict(best_state)
            logger.info(f"Restored best model weights")

        self.is_trained = True
        if do_cls:
            logger.info(f"Best val accuracy: {self.best_val_acc:.2f}%")
        if do_reg:
            logger.info(f"Best val RMSE: {self.best_val_rmse:.4f}")

        # --- Build reference index for kNN UMAP projection ---
        if do_reg and hasattr(self, '_y_reg_raw') and self._y_reg_raw is not None:
            logger.info("Building reference index for kNN UMAP projection...")
            self._build_ref_index(X, self._y_reg_raw, y_cls=y_cls)
            # Clean up temporary raw UMAP data
            del self._y_reg_raw
            self._y_reg_raw = None

    def _train_epoch_multitask_fixed(self, train_loader, cls_criterion, reg_criterion,
                                      model_optimizer):
        """Multi-task training with fixed equal weights (no GradNorm)."""
        self.train()
        total_loss = 0
        total_cls_loss = 0
        total_reg_loss = 0
        correct = 0
        total_samples = 0
        total_se = 0

        for batch in train_loader:
            batch_X = batch['X'].to(self.device)
            batch_y_cls = batch['y_cls'].to(self.device)
            batch_y_reg = batch['y_reg'].to(self.device)

            outputs = self(batch_X)
            cls_loss = cls_criterion(outputs['cls'], batch_y_cls)
            if self.use_uncertainty:
                raw_reg_out = torch.cat([outputs['reg_mean'], outputs['reg_log_var']], dim=-1)
                reg_loss = reg_criterion(raw_reg_out, batch_y_reg)
            else:
                reg_loss = reg_criterion(outputs['reg'], batch_y_reg)
            loss = (cls_loss + reg_loss) / 2

            model_optimizer.zero_grad()
            loss.backward()
            model_optimizer.step()

            n = batch_X.size(0)
            total_loss += loss.item() * n
            total_cls_loss += cls_loss.item() * n
            total_reg_loss += reg_loss.item() * n
            _, predicted = outputs['cls'].max(1)
            correct += predicted.eq(batch_y_cls).sum().item()
            total_samples += n
            total_se += ((outputs['reg'].detach() - batch_y_reg) ** 2).sum().item()

        n = total_samples
        return {
            'loss': total_loss / n,
            'cls_loss': total_cls_loss / n,
            'reg_loss': total_reg_loss / n,
            'accuracy': 100.0 * correct / n,
            'rmse': np.sqrt(total_se / (n * self.embedding_dim)),
        }

    # =========================================================================
    # Reference index for kNN UMAP projection (Seurat/Symphony-style)
    # =========================================================================

    def build_reference_index(self, adata: sc.AnnData,
                               embedding_key: str = 'umap',
                               label_column: Optional[str] = None,
                               max_ref_cells: int = 20000,
                               n_per_type: int = 1000) -> None:
        """Build reference index for kNN UMAP projection (public API).

        Use this to add kNN UMAP projection to an EXISTING trained model
        without retraining. After calling this, save the model to persist.

        Example:
            model = MultiTaskModel.load("model.pt")
            model.build_reference_index(ref_adata, "X_umap", "cell_type")
            model.save("model_with_ref.pt")

        Args:
            adata: Reference AnnData with raw counts and UMAP coordinates.
            embedding_key: obsm key for UMAP coordinates (default 'umap').
            label_column: obs column for stratified subsampling (optional).
        """
        if not self.is_trained:
            raise ValueError("Model must be trained first.")
        if embedding_key not in adata.obsm:
            raise ValueError(f"'{embedding_key}' not in adata.obsm. "
                             f"Available: {list(adata.obsm.keys())}")

        # Prepare binary feature matrix aligned to model genes
        var_genes = self.get_feature_names()
        query_genes = set(adata.var_names)
        shared = [g for g in var_genes if g in query_genes]
        logger.info(f"Building reference index: {len(shared)}/{len(var_genes)} genes shared")

        import scipy.sparse as _sp
        X_sub = adata[:, shared].X
        if _sp.issparse(X_sub):
            X_sub = X_sub.toarray()
        X_sub = (X_sub > 0).astype(np.float32)

        if len(shared) < len(var_genes):
            shared_pos = np.array([i for i, g in enumerate(var_genes) if g in query_genes])
            X = np.zeros((adata.n_obs, len(var_genes)), dtype=np.float32)
            X[:, shared_pos] = X_sub
        else:
            X = X_sub

        y_umap = adata.obsm[embedding_key].astype(np.float32)

        y_cls = None
        if label_column is not None and label_column in adata.obs:
            from sklearn.preprocessing import LabelEncoder as _LE
            le = _LE()
            y_cls = le.fit_transform(adata.obs[label_column].values)

        self._build_ref_index(X, y_umap, y_cls=y_cls,
                              max_ref_cells=max_ref_cells,
                              n_per_type=n_per_type)

    def _build_ref_index(self, X: np.ndarray, y_umap_raw: np.ndarray,
                         y_cls: Optional[np.ndarray] = None,
                         max_ref_cells: int = 20000,
                         n_per_type: int = 1000) -> None:
        """Build reference index for kNN-based UMAP projection.

        After training, computes shared embeddings for reference cells and
        stores a subsampled set for fast kNN projection at prediction time.
        Inspired by Seurat ProjectUMAP / Symphony mapQuery.

        Args:
            X: Binary expression matrix of reference cells.
            y_umap_raw: Original (un-normalized) UMAP coordinates.
            y_cls: Cell type labels (encoded), for stratified subsampling.
            max_ref_cells: Max reference cells to store.
            n_per_type: Max cells per cell type (stratified sampling).
        """
        self.eval()
        n_total = X.shape[0] if not sp.issparse(X) else X.shape[0]

        # --- Compute shared embeddings in batches ---
        emb_list = []
        batch_size = 2048
        with torch.no_grad():
            for i in range(0, n_total, batch_size):
                end = min(i + batch_size, n_total)
                if sp.issparse(X):
                    batch_X = torch.FloatTensor(X[i:end].toarray()).to(self.device)
                else:
                    batch_X = torch.FloatTensor(X[i:end]).to(self.device)
                outputs = self(batch_X)
                emb_list.append(outputs['embedding'].cpu().numpy())
        emb = np.concatenate(emb_list, axis=0)  # (n_total, hidden_size)

        # --- Stratified subsampling ---
        if n_total > max_ref_cells:
            rng = np.random.RandomState(self.random_state)
            if y_cls is not None:
                indices = []
                for cls_id in np.unique(y_cls):
                    cls_idx = np.where(y_cls == cls_id)[0]
                    n_sample = min(len(cls_idx), n_per_type)
                    indices.extend(rng.choice(cls_idx, n_sample, replace=False))
                indices = np.array(indices)
                # If still too many, randomly subsample
                if len(indices) > max_ref_cells:
                    indices = rng.choice(indices, max_ref_cells, replace=False)
            else:
                indices = rng.choice(n_total, max_ref_cells, replace=False)
            emb = emb[indices]
            y_umap_raw = y_umap_raw[indices]
            logger.info(f"Reference index subsampled: {n_total} -> {len(indices)} cells")
        else:
            logger.info(f"Reference index: using all {n_total} cells")

        # L2-normalize embeddings for cosine similarity (pre-compute)
        norms = np.linalg.norm(emb, axis=1, keepdims=True).clip(min=1e-8)
        emb_normed = emb / norms

        # Store as float16 for memory efficiency
        self._ref_emb = emb_normed.astype(np.float16)
        self._ref_umap = y_umap_raw.astype(np.float32)
        logger.info(f"Reference index built: {self._ref_emb.shape[0]} cells × "
                     f"{self._ref_emb.shape[1]}-dim embedding, "
                     f"storage ≈ {self._ref_emb.nbytes / 1024 / 1024:.1f} MB")

    def _knn_umap_project(self, query_emb: np.ndarray, k: int = 30,
                          temperature: float = 10.0) -> np.ndarray:
        """Project query cells to reference UMAP via kNN in shared embedding space.

        Uses cosine similarity with softmax-weighted averaging of reference
        UMAP coordinates. GPU-accelerated for speed.

        Args:
            query_emb: Query shared embeddings (n_query, hidden_size).
            k: Number of nearest neighbors.
            temperature: Softmax temperature (higher = sharper weights).

        Returns:
            Projected UMAP coordinates (n_query, embedding_dim).
        """
        device = self.device

        # Load reference data to device
        ref_emb = torch.from_numpy(self._ref_emb.astype(np.float32)).to(device)
        ref_umap = torch.from_numpy(self._ref_umap).to(device)

        # Normalize query embeddings
        q = torch.from_numpy(query_emb.astype(np.float32)).to(device)
        q = F.normalize(q, p=2, dim=1)

        projected = []
        batch_size = 4096  # large batches for GPU efficiency
        n_query = q.shape[0]

        for i in range(0, n_query, batch_size):
            batch_q = q[i:min(i + batch_size, n_query)]

            # Cosine similarity (ref is already L2-normalized)
            sim = batch_q @ ref_emb.T  # (batch, n_ref)

            # Top-k neighbors
            top_k_sim, top_k_idx = sim.topk(k, dim=1)  # (batch, k)

            # Softmax weights (temperature controls sharpness)
            weights = F.softmax(top_k_sim * temperature, dim=1)  # (batch, k)

            # Weighted average of reference UMAP coordinates
            top_k_umap = ref_umap[top_k_idx]  # (batch, k, 2)
            proj = (weights.unsqueeze(-1) * top_k_umap).sum(dim=1)  # (batch, 2)
            projected.append(proj)

        result = torch.cat(projected, dim=0).cpu().numpy()
        logger.info(f"kNN UMAP projection: {n_query} query cells -> "
                     f"{self._ref_emb.shape[0]} reference cells (k={k})")
        return result

    # =========================================================================
    # Prediction
    # =========================================================================

    def predict(self, X: Union[np.ndarray, sp.spmatrix],
                ood_threshold: Optional[float] = None,
                umap_method: str = 'auto',
                knn_k: int = 30,
                knn_temperature: float = 10.0) -> Dict[str, np.ndarray]:
        """Predict cell types, embeddings, shared representation, and uncertainty.

        Cell type predictions always retain the model's best-guess label.
        OOD status is provided as separate 'confidence' and 'is_ood' fields
        so that users can decide their own filtering strategy downstream.

        Args:
            X: Input binary gene expression matrix.
            ood_threshold: Confidence threshold for computing is_ood flag.
                Cells with max probability < threshold are flagged is_ood=True.
                If None, uses self.ood_threshold (default 0.5).
            umap_method: UMAP projection method:
                'auto' (default): use kNN if reference index available, else regression
                'knn': force kNN projection (error if no reference index)
                'regression': force regression head only
            knn_k: Number of nearest neighbors for kNN projection.
            knn_temperature: Softmax temperature for kNN weights.

        Returns:
            Dict with keys:
                'shared_embedding': backbone representation
                'cell_types': best-guess predicted labels (always concrete)
                'probabilities': temperature-scaled softmax probabilities
                'confidence': max probability per cell
                'is_ood': boolean mask (confidence < threshold)
                'embeddings': predicted coordinates (mean)
                'uncertainty': per-coordinate std (if use_uncertainty=True)
        """
        if not self.is_trained:
            raise ValueError("Model must be trained before making predictions")

        self.eval()
        self.to(self.device)

        if sp.issparse(X):
            X = X.toarray()
        X_tensor = torch.FloatTensor(X)

        result = {}
        with torch.no_grad():
            outputs = self(X_tensor.to(self.device))

        result['shared_embedding'] = outputs['embedding'].cpu().numpy()

        if 'cls' in outputs and self.has_cls:
            # Temperature-scaled softmax
            logits = outputs['cls'] / self.temperature
            probs = F.softmax(logits, dim=1)
            confidence, predictions = probs.max(1)

            confidence_np = confidence.cpu().numpy()
            predictions_np = predictions.cpu().numpy()

            # OOD detection
            threshold = ood_threshold if ood_threshold is not None else self.ood_threshold
            is_ood = confidence_np < threshold

            if self.label_encoder is not None:
                cell_types = self.label_encoder.inverse_transform(predictions_np)
                cell_types = cell_types.astype(object)
                # Keep the best-guess label for all cells;
                # users can filter via 'is_ood' or 'confidence' downstream.
            else:
                cell_types = predictions_np.copy()

            result['cell_types'] = cell_types
            result['probabilities'] = probs.cpu().numpy()
            result['confidence'] = confidence_np
            result['is_ood'] = is_ood

        if 'reg' in outputs and self.has_reg:
            # --- Determine UMAP projection method ---
            has_ref = (hasattr(self, '_ref_emb') and self._ref_emb is not None
                       and hasattr(self, '_ref_umap') and self._ref_umap is not None)

            use_knn = False
            if umap_method == 'auto':
                use_knn = has_ref
            elif umap_method == 'knn':
                if not has_ref:
                    raise ValueError("kNN UMAP projection requested but no reference "
                                     "index available. Retrain model or use 'auto'.")
                use_knn = True
            # umap_method == 'regression' -> use_knn = False

            if use_knn:
                # kNN-based projection (Seurat/Symphony-style)
                shared_emb = result['shared_embedding']
                emb = self._knn_umap_project(shared_emb, k=knn_k,
                                              temperature=knn_temperature)
                result['embeddings'] = emb
                result['umap_method'] = 'knn'
            else:
                # Regression head fallback
                emb = outputs['reg'].cpu().numpy()
                # Inverse z-score: map predicted normalized coords back to original UMAP space
                if hasattr(self, '_reg_mean') and self._reg_mean is not None:
                    emb = emb * self._reg_std + self._reg_mean
                result['embeddings'] = emb
                result['umap_method'] = 'regression'

            # Aleatoric uncertainty (std per coordinate, regression head only)
            if self.use_uncertainty and 'reg_log_var' in outputs:
                log_var = outputs['reg_log_var'].cpu().numpy()
                uncertainty = np.sqrt(np.exp(log_var))
                # Scale uncertainty back to original space
                if hasattr(self, '_reg_std') and self._reg_std is not None:
                    uncertainty = uncertainty * self._reg_std
                result['uncertainty'] = uncertainty

        return result

    def get_shared_embedding(self, X: Union[np.ndarray, sp.spmatrix]) -> np.ndarray:
        """Extract the shared embedding (hidden_size dim) for downstream analysis."""
        result = self.predict(X)
        return result['shared_embedding']

    # =========================================================================
    # v3: Temperature scaling (Guo et al., 2017)
    # =========================================================================

    def calibrate_temperature(self, adata: sc.AnnData,
                              label_column: str,
                              var_genes: Optional[List[str]] = None,
                              max_iter: int = 50,
                              lr: float = 0.01) -> float:
        """Post-hoc temperature scaling on a validation set.

        Optimizes a single scalar T to minimize NLL on held-out data,
        producing calibrated probabilities where confidence ≈ accuracy.

        Args:
            adata: Validation AnnData (should NOT be training data).
            label_column: obs column with ground-truth labels.
            var_genes: gene subset (must match training genes).
            max_iter: Optimization iterations.
            lr: Learning rate for temperature optimization.

        Returns:
            Optimized temperature value.
        """
        if not self.is_trained:
            raise ValueError("Model must be trained before calibration")

        genes = var_genes if var_genes is not None else self.var_genes
        adata_sub = adata[:, genes].copy()
        X = (adata_sub.X > 0).astype(np.float32)
        if sp.issparse(X):
            X = X.toarray()

        y = self.label_encoder.transform(adata_sub.obs[label_column].values)

        X_tensor = torch.FloatTensor(X).to(self.device)
        y_tensor = torch.LongTensor(y).to(self.device)

        self.eval()
        with torch.no_grad():
            outputs = self(X_tensor)
            logits = outputs['cls']

        # Optimize temperature
        temperature = nn.Parameter(torch.ones(1, device=self.device) * 1.5)
        optimizer = torch.optim.LBFGS([temperature], lr=lr, max_iter=max_iter)
        nll_criterion = nn.CrossEntropyLoss()

        def _eval():
            optimizer.zero_grad()
            loss = nll_criterion(logits / temperature, y_tensor)
            loss.backward()
            return loss

        optimizer.step(_eval)

        optimal_T = temperature.item()
        self.temperature.data.fill_(optimal_T)

        # Compute calibration stats
        with torch.no_grad():
            calibrated_probs = F.softmax(logits / optimal_T, dim=1)
            confidence, preds = calibrated_probs.max(1)
            acc = preds.eq(y_tensor).float().mean().item()
            mean_conf = confidence.mean().item()

        logger.info(f"Temperature calibration: T={optimal_T:.4f}")
        logger.info(f"  Calibrated accuracy: {acc * 100:.2f}%, mean confidence: {mean_conf * 100:.2f}%")
        logger.info(f"  Gap (|acc - conf|): {abs(acc - mean_conf) * 100:.2f}%")

        return optimal_T

    def set_ood_threshold(self, threshold: float) -> None:
        """Set the OOD detection confidence threshold.

        Args:
            threshold: Cells with max_prob < threshold are flagged
                is_ood=True. Typical range: 0.3 ~ 0.7. Use
                calibrate_temperature first for meaningful probability values.
                Note: cell_type predictions always retain the best-guess
                label; is_ood is an advisory flag for downstream filtering.
        """
        self.ood_threshold = threshold
        logger.info(f"OOD threshold set to {threshold:.4f}")

    # =========================================================================
    # v3: Transfer learning
    # =========================================================================

    def freeze_backbone(self) -> None:
        """Freeze shared backbone parameters (InputBlock + ResNet).

        After freezing, only task heads are trainable, enabling fast
        fine-tuning on a new atlas with the same feature space.
        """
        for param in self.input_block.parameters():
            param.requires_grad = False
        for param in self.resnet_blocks.parameters():
            param.requires_grad = False
        for param in self.resnet_last.parameters():
            param.requires_grad = False

        trainable = sum(p.numel() for p in self.parameters() if p.requires_grad)
        total = sum(p.numel() for p in self.parameters())
        logger.info(f"Backbone frozen: {trainable}/{total} params trainable ({100*trainable/total:.1f}%)")

    def unfreeze_backbone(self) -> None:
        """Unfreeze all backbone parameters for full fine-tuning."""
        for param in self.parameters():
            param.requires_grad = True
        logger.info("All parameters unfrozen")

    def fine_tune(self, adata: sc.AnnData,
                  label_column: Optional[str] = None,
                  embedding_key: Optional[str] = None,
                  var_genes: Optional[List[str]] = None,
                  freeze: bool = True,
                  new_num_epochs: Optional[int] = None,
                  new_learning_rate: Optional[float] = None) -> None:
        """Fine-tune model on a new atlas (transfer learning).

        Typical workflow:
            1. Train on reference atlas A  (model.fit)
            2. Fine-tune on atlas B        (model.fine_tune, freeze backbone)
            3. Optionally unfreeze for full fine-tuning

        Args:
            adata: New AnnData to fine-tune on.
            label_column: obs column for cell type labels.
            embedding_key: obsm key for embedding coords.
            var_genes: Gene subset (should match original training genes).
            freeze: If True, freeze backbone and only train heads.
            new_num_epochs: Override num_epochs for fine-tuning (default: original/2).
            new_learning_rate: Override learning rate (default: original/10).
        """
        if not self.is_trained:
            raise ValueError("Model must be trained before fine-tuning")

        # Save original config
        orig_epochs = self.num_epochs
        orig_lr = self.learning_rate
        orig_gradnorm = self.use_gradnorm

        # Fine-tuning defaults
        self.num_epochs = new_num_epochs or max(orig_epochs // 2, 10)
        self.learning_rate = new_learning_rate or orig_lr / 10.0

        # Use training genes if not specified
        if var_genes is None:
            var_genes = self.var_genes

        # Rebuild heads if new atlas has different classes
        if label_column is not None:
            new_classes = adata.obs[label_column].nunique()
            if self.num_classes is not None and new_classes != self.num_classes:
                logger.info(f"Rebuilding cls_head: {self.num_classes} -> {new_classes} classes")
                self.label_encoder = None  # Reset for new classes
                self.has_cls = False       # Will be rebuilt in fit()

        if freeze:
            self.freeze_backbone()
            # Disable GradNorm when backbone is frozen (no shared gradients)
            self.use_gradnorm = False

        logger.info(f"Fine-tuning: epochs={self.num_epochs}, lr={self.learning_rate:.1e}, freeze={freeze}")

        # Run training
        self.is_trained = False  # Allow fit() to run
        self.fit(adata, label_column=label_column,
                 embedding_key=embedding_key, var_genes=var_genes)

        # Restore original config
        self.num_epochs = orig_epochs
        self.learning_rate = orig_lr
        self.use_gradnorm = orig_gradnorm

        if freeze:
            self.unfreeze_backbone()

    # =========================================================================
    # v3: Interpretability (Integrated Gradients)
    # =========================================================================

    def explain(self, X: Union[np.ndarray, sp.spmatrix],
                task: str = 'cls',
                n_steps: int = 50,
                target_class: Optional[int] = None) -> np.ndarray:
        """Compute gene importance scores via Integrated Gradients.

        Measures the contribution of each input gene to the model's prediction
        by integrating gradients along the path from a zero baseline to the
        actual input. This is equivalent to a SHAP-like attribution.

        Args:
            X: Input binary gene expression matrix (n_cells x n_genes).
            task: 'cls' for classification attributions, 'reg' for regression.
            n_steps: Number of interpolation steps (higher = more precise).
            target_class: For classification, which class to explain.
                If None, uses the predicted class for each cell.

        Returns:
            attributions: np.ndarray (n_cells x n_genes) with per-gene importance.
                Positive = gene increases prediction; negative = decreases.
        """
        if not self.is_trained:
            raise ValueError("Model must be trained before explanation")

        self.eval()
        self.to(self.device)

        if sp.issparse(X):
            X = X.toarray()
        X_tensor = torch.FloatTensor(X).to(self.device)

        # Baseline: zero vector (no genes expressed)
        baseline = torch.zeros_like(X_tensor)

        # Get target classes if not specified
        if task == 'cls' and target_class is None:
            with torch.no_grad():
                outputs = self(X_tensor)
                target_class_tensor = outputs['cls'].argmax(dim=1)
        elif task == 'cls':
            target_class_tensor = torch.full((X_tensor.shape[0],), target_class,
                                              dtype=torch.long, device=self.device)

        # Integrated Gradients
        all_grads = torch.zeros_like(X_tensor)

        for step in range(n_steps + 1):
            alpha = step / n_steps
            interpolated = baseline + alpha * (X_tensor - baseline)
            interpolated.requires_grad_(True)

            outputs = self(interpolated)

            if task == 'cls':
                # Gather logits for target classes
                logits = outputs['cls']
                target_logits = logits.gather(1, target_class_tensor.unsqueeze(1)).squeeze(1)
                score = target_logits.sum()
            else:
                # Sum of all regression outputs
                score = outputs['reg'].sum()

            score.backward()
            all_grads += interpolated.grad.detach()
            interpolated.grad = None

        # Average gradients and multiply by input difference
        avg_grads = all_grads / (n_steps + 1)
        attributions = (X_tensor - baseline).detach() * avg_grads

        return attributions.cpu().numpy()

    def explain_top_genes(self, X: Union[np.ndarray, sp.spmatrix],
                          task: str = 'cls',
                          top_k: int = 20,
                          n_steps: int = 50,
                          target_class: Optional[int] = None) -> Dict[str, np.ndarray]:
        """Get top contributing genes with importance scores.

        Args:
            X: Input binary gene expression matrix.
            task: 'cls' or 'reg'.
            top_k: Number of top genes to return.
            n_steps: Integrated gradient steps.
            target_class: For classification, compute attributions w.r.t. a
                specific class index. If None, uses each cell's predicted class.

        Returns:
            Dict with:
                'gene_names': top gene names (top_k,)
                'importance': mean absolute importance scores (top_k,)
                'attributions': full attribution matrix (n_cells x n_genes)
        """
        attributions = self.explain(X, task=task, n_steps=n_steps,
                                    target_class=target_class)

        # Mean absolute attribution per gene (across all cells)
        mean_importance = np.abs(attributions).mean(axis=0)
        top_idx = np.argsort(mean_importance)[::-1][:top_k]

        result = {
            'importance': mean_importance[top_idx],
            'attributions': attributions,
        }

        if self.var_genes is not None:
            result['gene_names'] = np.array(self.var_genes)[top_idx]
        else:
            result['gene_names'] = np.array([f'gene_{i}' for i in top_idx])

        return result

    def explain_per_class(self, X: Union[np.ndarray, sp.spmatrix],
                          labels: np.ndarray,
                          top_k: int = 20,
                          n_steps: int = 50,
                          max_cells_per_class: int = 50
                          ) -> Dict[str, Dict[str, np.ndarray]]:
        """Compute per-class gene importance rankings.

        For each cell type, selects cells belonging to that class and computes
        Integrated Gradients w.r.t. that class's logit, then ranks genes.

        Args:
            X: Input binary gene expression matrix (n_cells x n_genes).
            labels: Array of cell type labels (n_cells,), matching label_encoder classes.
            top_k: Number of top genes to return per class.
            n_steps: Integrated gradient interpolation steps.
            max_cells_per_class: Max cells to use per class (random sample for speed).

        Returns:
            Dict[class_name -> Dict] with per-class results:
                'gene_names': np.ndarray (top_k,) top gene names
                'importance': np.ndarray (top_k,) importance scores
                'n_cells': int, number of cells used
        """
        if not self.is_trained or not hasattr(self, 'label_encoder'):
            raise ValueError("Model must be trained with classification task")

        if sp.issparse(X):
            X = X.toarray()
        labels = np.asarray(labels)

        class_names = self.label_encoder.classes_
        results = {}

        for cls_idx, cls_name in enumerate(class_names):
            # Select cells of this class
            mask = labels == cls_name
            if mask.sum() == 0:
                continue

            X_cls = X[mask]

            # Subsample for speed
            if X_cls.shape[0] > max_cells_per_class:
                rng = np.random.RandomState(42)
                idx = rng.choice(X_cls.shape[0], max_cells_per_class, replace=False)
                X_cls = X_cls[idx]

            # Compute IG for this specific class logit
            attributions = self.explain(X_cls, task='cls', n_steps=n_steps,
                                        target_class=cls_idx)

            # Rank genes by mean absolute attribution
            mean_imp = np.abs(attributions).mean(axis=0)
            top_idx = np.argsort(mean_imp)[::-1][:top_k]

            gene_names = (np.array(self.var_genes)[top_idx]
                          if self.var_genes else
                          np.array([f'gene_{i}' for i in top_idx]))

            results[cls_name] = {
                'gene_names': gene_names,
                'importance': mean_imp[top_idx],
                'n_cells': X_cls.shape[0],
            }

            logger.info(f"  {cls_name}: top gene={gene_names[0]} "
                         f"(imp={mean_imp[top_idx[0]]:.4f}), "
                         f"n_cells={X_cls.shape[0]}")

        return results

    # =========================================================================
    # Save / Load
    # =========================================================================

    def save(self, save_path: str, metadata: dict = None) -> None:
        """Save model and metadata to .pt file."""
        if not self.is_trained:
            raise ValueError("Model must be trained before saving")
        if not save_path.endswith('.pt'):
            save_path += '.pt'

        save_dict = {
            'model_state_dict': self.state_dict(),
            'model_config': {
                'input_size': self.input_size,
                'num_classes': self.num_classes,
                'embedding_dim': self.embedding_dim,
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
                'early_stopping_monitor': self.early_stopping_monitor,
                'learning_rate': self.learning_rate,
                'use_gradnorm': self.use_gradnorm,
                'gradnorm_alpha': self.gradnorm_alpha,
                'label_smoothing': self.label_smoothing,
                'use_lr_scheduler': self.use_lr_scheduler,
                'use_uncertainty': self.use_uncertainty,
                'random_state': self.random_state,
            },
            'has_cls': self.has_cls,
            'has_reg': self.has_reg,
            'history': self.history,
            'best_val_acc': self.best_val_acc,
            'best_val_rmse': self.best_val_rmse,
            'ood_threshold': self.ood_threshold,
        }
        if self.label_encoder is not None:
            save_dict['label_encoder'] = {'classes': self.label_encoder.classes_}
        if self.var_genes is not None:
            save_dict['var_genes'] = self.var_genes
        if self.embedding_key is not None:
            save_dict['embedding_key'] = self.embedding_key
        # Save UMAP z-score stats for inverse-transform during prediction
        if hasattr(self, '_reg_mean') and self._reg_mean is not None:
            save_dict['reg_mean'] = self._reg_mean
            save_dict['reg_std']  = self._reg_std
        # Save reference index for kNN UMAP projection
        if hasattr(self, '_ref_emb') and self._ref_emb is not None:
            save_dict['ref_emb']  = self._ref_emb   # float16, pre-normalized
            save_dict['ref_umap'] = self._ref_umap  # float32, original coords
            logger.info(f"Saving reference index: {self._ref_emb.shape[0]} cells, "
                         f"≈ {(self._ref_emb.nbytes + self._ref_umap.nbytes) / 1024 / 1024:.1f} MB")
        if metadata is not None:
            self.metadata.update(metadata)
        save_dict['metadata'] = self.metadata

        torch.save(save_dict, save_path)
        logger.info(f"Model saved to {save_path}")

    @classmethod
    def load(cls, load_path: str) -> 'MultiTaskModel':
        """Load a saved model (backward-compatible with legacy checkpoints)."""
        if not load_path.endswith('.pt'):
            load_path += '.pt'

        save_dict = torch.load(load_path, weights_only=False,
                               map_location=torch.device('cpu'))
        config = save_dict['model_config']

        # --- Backward compatibility: legacy key names ---
        # Old checkpoints used 'num_cell_types'; new code uses 'num_classes'
        num_classes = config.get('num_classes',
                                 config.get('num_cell_types', None))
        if num_classes is None and 'cell_type_names' in save_dict:
            num_classes = len(save_dict['cell_type_names'])

        model = cls(
            input_size=config['input_size'],
            num_classes=num_classes,
            embedding_dim=config.get('embedding_dim', 2),
            hidden_size=config['hidden_size'],
            head_hidden_size=config.get('head_hidden_size', 256),
            n_resnet_blocks=config['n_resnet_blocks'],
            input_dropout_rate=config['input_dropout_rate'],
            resnet_dropout_rate=config['resnet_dropout_rate'],
            head_dropout_rate=config.get('head_dropout_rate', 0.1),
            batch_size=config['batch_size'],
            test_size=config['test_size'],
            num_epochs=config['num_epochs'],
            early_stopping_patience=config['early_stopping_patience'],
            early_stopping_monitor=config.get('early_stopping_monitor', 'cls'),
            learning_rate=config['learning_rate'],
            use_gradnorm=config.get('use_gradnorm', False),
            gradnorm_alpha=config.get('gradnorm_alpha', 0.15),
            label_smoothing=config.get('label_smoothing', 0.0),
            use_lr_scheduler=config.get('use_lr_scheduler', False),
            use_uncertainty=config.get('use_uncertainty', False),
            random_state=config['random_state'],
        )

        model.has_cls = save_dict.get('has_cls', num_classes is not None)
        model.has_reg = save_dict.get('has_reg', False)
        model.ood_threshold = save_dict.get('ood_threshold', 0.5)

        # --- Backward compat: inject missing tensors ---
        state_dict = save_dict['model_state_dict']
        if 'temperature' not in state_dict:
            state_dict['temperature'] = torch.ones(1)

        # Legacy checkpoints may have extra/missing keys; load with strict=False
        missing, unexpected = model.load_state_dict(state_dict, strict=False)
        if missing:
            logger.warning(f"Missing keys in checkpoint (initialized to default): "
                           f"{missing}")
        if unexpected:
            logger.warning(f"Unexpected keys in checkpoint (ignored): {unexpected}")
        model.to(model.device)

        # --- Label encoder: new format ('label_encoder') or legacy ('cell_type_names') ---
        if 'label_encoder' in save_dict:
            model.label_encoder = LabelEncoder()
            model.label_encoder.classes_ = save_dict['label_encoder']['classes']
        elif 'cell_type_names' in save_dict:
            model.label_encoder = LabelEncoder()
            model.label_encoder.classes_ = np.array(save_dict['cell_type_names'])

        model.var_genes = save_dict.get('var_genes', None)
        model.embedding_key = save_dict.get('embedding_key', None)
        # Restore UMAP z-score stats for inverse-transform during prediction
        model._reg_mean = save_dict.get('reg_mean', None)
        model._reg_std  = save_dict.get('reg_std', None)
        # Restore reference index for kNN UMAP projection
        model._ref_emb  = save_dict.get('ref_emb', None)
        model._ref_umap = save_dict.get('ref_umap', None)
        if model._ref_emb is not None:
            logger.info(f"Reference index loaded: {model._ref_emb.shape[0]} cells "
                         f"for kNN UMAP projection")
        model.metadata = save_dict.get('metadata', {})
        model.history = save_dict.get('history', {})
        model.best_val_acc = save_dict.get('best_val_acc',
                                            save_dict.get('best_val_loss', 0.0))
        model.best_val_rmse = save_dict.get('best_val_rmse', float('inf'))
        model.is_trained = True

        logger.info(f"Model loaded from {load_path} (device: {model.device})")
        return model

    # =========================================================================
    # Utilities
    # =========================================================================

    def print_training_parameters(self) -> None:
        """Print model configuration and training status."""
        print("=" * 65)
        print("MultiTaskModel Training Parameters")
        print("=" * 65)

        print("\nModel Architecture:")
        print(f"  Input size:          {self.input_size}")
        print(f"  Hidden size:         {self.hidden_size}")
        print(f"  Head hidden size:    {self.head_hidden_size}")
        print(f"  ResNet blocks:       {self.n_resnet_blocks}")
        print(f"  Input dropout:       {self.input_dropout_rate}")
        print(f"  ResNet dropout:      {self.resnet_dropout_rate}")
        print(f"  Head dropout:        {self.head_dropout_rate}")
        if self.has_cls:
            if self.head_hidden_size > 0:
                print(f"  Classification head: {self.hidden_size} -> {self.head_hidden_size} -> {self.num_classes}")
            else:
                print(f"  Classification head: {self.hidden_size} -> {self.num_classes}")
        if self.has_reg:
            if self.head_hidden_size > 0:
                print(f"  Regression head:     {self.hidden_size} -> {self.head_hidden_size} -> {self.embedding_dim}")
            else:
                print(f"  Regression head:     {self.hidden_size} -> {self.embedding_dim}")

        print("\nTraining Hyperparameters:")
        print(f"  Batch size:          {self.batch_size}")
        print(f"  Test size:           {self.test_size}")
        print(f"  Num epochs:          {self.num_epochs}")
        print(f"  Early stopping:      patience={self.early_stopping_patience}, monitor='{self.early_stopping_monitor}'")
        print(f"  Learning rate:       {self.learning_rate}")
        print(f"  LR scheduler:        {'CosineAnnealingWarmRestarts' if self.use_lr_scheduler else 'None'}")
        print(f"  GradNorm:            {self.use_gradnorm} (alpha={self.gradnorm_alpha})")
        print(f"  Label smoothing:     {self.label_smoothing}")

        print("\nv3 Features:")
        print(f"  Temperature:         {self.temperature.item():.4f}")
        print(f"  OOD threshold:       {self.ood_threshold}")
        print(f"  Uncertainty:         {'NLL (aleatoric)' if self.use_uncertainty else 'off'}")

        print(f"\nDevice: {self.device}")
        print(f"Trained: {self.is_trained}")

        if self.is_trained:
            if self.best_val_acc > 0:
                print(f"Best val accuracy:     {self.best_val_acc:.2f}%")
            if self.best_val_rmse < float('inf'):
                print(f"Best val RMSE:         {self.best_val_rmse:.4f}")
            if self.label_encoder is not None:
                print(f"Classes ({len(self.label_encoder.classes_)}): {list(self.label_encoder.classes_)}")
            if self.var_genes is not None:
                print(f"Variable genes:        {len(self.var_genes)}")
            if len(self.history.get('w_cls', [])) > 0:
                print(f"Final GradNorm weights: w_cls={self.history['w_cls'][-1]:.3f}, w_reg={self.history['w_reg'][-1]:.3f}")

        print("=" * 65)

    def get_class_names(self) -> Optional[List[str]]:
        """Get class names from label encoder."""
        if self.label_encoder is not None:
            return list(self.label_encoder.classes_)
        return None

    def get_feature_names(self) -> Optional[List[str]]:
        """Get variable gene names used for training."""
        return self.var_genes
