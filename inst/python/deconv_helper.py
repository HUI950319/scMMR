"""
deconv_helper.py
─────────────────
Python functions exposed to R via reticulate::source_python().
Provides the R-Python bridge for bulk RNA-seq deconvolution.

All functions return plain Python dicts / numpy arrays so that
reticulate::py_to_r() can convert them to R lists / matrices directly.

Usage from R:
    library(reticulate)
    source_python("deconv_helper.py")
    adata  <- mt_load_query("reference.h5ad")
    hvg    <- mt_select_hvg(adata, 6000L)
    pb     <- deconv_generate_pseudobulk(adata, "cell_type", 5000L, 500L, hvg)
    model  <- deconv_create_model(length(hvg), length(pb$cell_type_names))
    deconv_train(model, pb$X_pb, pb$y_pb, pb$cell_type_names)
    deconv_save_model(model, "deconv_model.pt")

    # Prediction
    model  <- deconv_load_model("deconv_model.pt")
    X_bulk <- deconv_prepare_bulk(bulk_matrix, model, gene_names, sample_names)
    result <- deconv_predict(model, X_bulk)
"""

import os
import sys
import numpy as np
import scipy.sparse as sp

# ── Ensure deconv_model.py is importable ────────────────────────────────────
try:
    _dir = os.path.dirname(os.path.abspath(__file__))
except NameError:
    _dir = os.getcwd()
if _dir not in sys.path:
    sys.path.insert(0, _dir)

from deconv_model import DeconvModel, generate_pseudobulk

import torch


# ═════════════════════════════════════════════════════════════════════════════
# 1. Device Management
# ═════════════════════════════════════════════════════════════════════════════

_deconv_device = 'auto'


def deconv_set_device(device_str: str = 'auto'):
    """Set compute device for deconvolution models.

    Args:
        device_str: 'auto' (use GPU if available), 'cpu', or 'cuda'.
    """
    global _deconv_device
    _deconv_device = device_str

    if device_str == 'auto':
        dev = 'cuda' if torch.cuda.is_available() else 'cpu'
    else:
        dev = device_str

    print(f"  Deconv device: {dev}"
          f"{' (CUDA)' if dev == 'cuda' else ' (CPU)'}")


# ═════════════════════════════════════════════════════════════════════════════
# 2. Model Creation
# ═════════════════════════════════════════════════════════════════════════════

def deconv_create_model(input_size: int,
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
    """Create an untrained DeconvModel.

    Args:
        input_size: Number of genes (features).
        num_cell_types: Number of cell types (output dimension).
        **kwargs: Additional model hyperparameters.

    Returns:
        DeconvModel instance (untrained).
    """
    model = DeconvModel(
        input_size=int(input_size),
        num_cell_types=int(num_cell_types),
        hidden_size=int(hidden_size),
        head_hidden_size=int(head_hidden_size),
        n_resnet_blocks=int(n_resnet_blocks),
        input_dropout_rate=float(input_dropout_rate),
        resnet_dropout_rate=float(resnet_dropout_rate),
        head_dropout_rate=float(head_dropout_rate),
        batch_size=int(batch_size),
        test_size=float(test_size),
        num_epochs=int(num_epochs),
        early_stopping_patience=int(early_stopping_patience),
        learning_rate=float(learning_rate),
        loss_function=str(loss_function),
        use_lr_scheduler=bool(use_lr_scheduler),
        random_state=int(random_state),
        recon_weight=float(recon_weight),
    )

    # Set device
    if _deconv_device == 'auto':
        dev = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
    else:
        dev = torch.device(_deconv_device)
    model.device = dev
    model.to(dev)

    n_params = sum(p.numel() for p in model.parameters())
    print(f"  DeconvModel: {input_size} genes -> {num_cell_types} cell types "
          f"({n_params:,} params, device={dev})")
    return model


# ═════════════════════════════════════════════════════════════════════════════
# 3. Pseudo-bulk Generation
# ═════════════════════════════════════════════════════════════════════════════

def deconv_generate_pseudobulk(adata,
                                label_column: str,
                                n_samples: int = 5000,
                                n_cells_per_sample: int = 500,
                                var_genes=None,
                                target_sum: float = 1e6,
                                random_state: int = 42):
    """Generate pseudo-bulk training data from scRNA-seq reference.

    Args:
        adata: AnnData with raw counts.
        label_column: obs column for cell type labels.
        n_samples: Number of pseudo-bulk samples.
        n_cells_per_sample: Cells to mix per sample.
        var_genes: Gene list to use (from HVG selection).
        target_sum: CPM target (default 1e6).
        random_state: Random seed.

    Returns:
        dict with X_pb, y_pb, cell_type_names, var_genes.
    """
    if var_genes is not None:
        var_genes = list(np.array(var_genes, dtype=str))

    result = generate_pseudobulk(
        adata=adata,
        label_column=str(label_column),
        n_samples=int(n_samples),
        n_cells_per_sample=int(n_cells_per_sample),
        var_genes=var_genes,
        target_sum=float(target_sum),
        random_state=int(random_state)
    )
    return result


# ═════════════════════════════════════════════════════════════════════════════
# 4. Training
# ═════════════════════════════════════════════════════════════════════════════

def deconv_train(model, X_pb, y_pb, cell_type_names):
    """Train the deconvolution model on pseudo-bulk data.

    Args:
        model: DeconvModel instance.
        X_pb: (n_samples, n_genes) log-CPM pseudo-bulk matrix.
        y_pb: (n_samples, n_cell_types) proportion matrix.
        cell_type_names: Ordered list of cell type names.

    Returns:
        dict with best_val_loss, best_val_corr, n_epochs.
    """
    X_pb = np.array(X_pb, dtype=np.float32)
    y_pb = np.array(y_pb, dtype=np.float32)
    cell_type_names = list(np.array(cell_type_names, dtype=str))

    model.fit(X_pb, y_pb, cell_type_names)

    return {
        'best_val_loss': float(model.best_val_loss),
        'best_val_corr': float(model.best_val_corr),
        'n_epochs': len(model.history.get('train_loss', []))
    }


# ═════════════════════════════════════════════════════════════════════════════
# 5. Save / Load
# ═════════════════════════════════════════════════════════════════════════════

def deconv_save_model(model, path: str):
    """Save trained DeconvModel to .pt file.

    Args:
        model: Trained DeconvModel.
        path: Output file path.
    """
    # Create parent directory if needed
    parent = os.path.dirname(os.path.abspath(path))
    os.makedirs(parent, exist_ok=True)
    model.save(str(path))


def deconv_load_model(path: str):
    """Load a trained DeconvModel from .pt file.

    Args:
        path: Path to the saved model.

    Returns:
        DeconvModel instance.
    """
    model = DeconvModel.load(str(path))

    # Apply device setting
    if _deconv_device == 'auto':
        dev = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
    else:
        dev = torch.device(_deconv_device)
    model.device = dev
    model.to(dev)

    return model


# ═════════════════════════════════════════════════════════════════════════════
# 6. Bulk Expression Preparation
# ═════════════════════════════════════════════════════════════════════════════

def deconv_prepare_bulk(bulk_matrix, model, gene_names, sample_names=None,
                         target_sum: float = 1e6):
    """Prepare bulk expression matrix for prediction.

    Aligns genes to model's var_genes, normalizes to log1p(CPM).

    Args:
        bulk_matrix: numpy array (genes x samples) or (samples x genes).
                     Raw counts or pre-normalized.
        model: Trained DeconvModel (provides var_genes).
        gene_names: Gene names corresponding to bulk_matrix rows.
        sample_names: Sample names (optional).
        target_sum: CPM normalization target.

    Returns:
        numpy array (n_samples, n_genes) log-CPM, aligned to model genes.
    """
    bulk_matrix = np.array(bulk_matrix, dtype=np.float64)
    gene_names = list(np.array(gene_names, dtype=str))
    model_genes = model.get_feature_names()

    if len(model_genes) == 0:
        raise ValueError("Model has no var_genes. Was it trained properly?")

    # Detect orientation: if nrows matches n_model_genes, assume genes x samples
    # Otherwise if ncols matches, assume samples x genes
    n_model_genes = len(model_genes)

    if bulk_matrix.shape[0] == len(gene_names):
        # genes x samples -> transpose to samples x genes
        X = bulk_matrix.T
    elif bulk_matrix.shape[1] == len(gene_names):
        # already samples x genes
        X = bulk_matrix
    else:
        # Default: assume genes x samples (most common in R)
        X = bulk_matrix.T

    n_samples = X.shape[0]

    # Build gene name lookup
    gene_to_idx = {g: i for i, g in enumerate(gene_names)}

    # Align to model genes
    shared_genes = [g for g in model_genes if g in gene_to_idx]
    n_shared = len(shared_genes)
    n_missing = n_model_genes - n_shared
    overlap_pct = n_shared / n_model_genes * 100

    print(f"  Gene alignment: {n_shared}/{n_model_genes} genes matched "
          f"({overlap_pct:.1f}%), {n_missing} missing (zero-filled)")

    if overlap_pct < 50:
        print(f"  WARNING: Low gene overlap ({overlap_pct:.1f}%). "
              f"Results may be unreliable.")

    # Extract shared genes and zero-fill missing
    X_aligned = np.zeros((n_samples, n_model_genes), dtype=np.float64)
    for j, g in enumerate(model_genes):
        if g in gene_to_idx:
            X_aligned[:, j] = X[:, gene_to_idx[g]]

    # Normalize to log1p(CPM) — same as pseudo-bulk training
    row_sums = X_aligned.sum(axis=1, keepdims=True)
    row_sums = np.maximum(row_sums, 1.0)  # avoid division by zero
    X_cpm = X_aligned / row_sums * target_sum
    X_log_cpm = np.log1p(X_cpm).astype(np.float32)

    print(f"  Prepared: {n_samples} samples x {n_model_genes} genes (log-CPM)")
    return X_log_cpm


# ═════════════════════════════════════════════════════════════════════════════
# 7. Prediction
# ═════════════════════════════════════════════════════════════════════════════

def deconv_predict(model, X_bulk, adaptive=False, mode='overall',
                   max_iter=3, step=300, lr=1e-4):
    """Predict cell type proportions from prepared bulk expression.

    Args:
        model: Trained DeconvModel.
        X_bulk: (n_samples, n_genes) log-CPM expression matrix.
        adaptive: If True, run TAPE-style adaptive stage.
        mode: 'overall' or 'high-resolution' (per-sample GEP).
        max_iter: Adaptive alternating rounds.
        step: Gradient steps per phase per round.
        lr: Adaptive learning rate.

    Returns:
        dict with:
            proportions: (n_samples, n_cell_types) float32
            cell_type_names: list of cell type names
            embedding: (n_samples, hidden_size) float32
            sigmatrix: (optional) GEP matrix when adaptive=True
            var_genes: (optional) gene names when adaptive=True
    """
    X_bulk = np.array(X_bulk, dtype=np.float32)
    results = model.predict(
        X_bulk,
        adaptive=bool(adaptive),
        mode=str(mode),
        max_iter=int(max_iter),
        step=int(step),
        lr=float(lr)
    )

    ret = {
        'proportions': results['proportions'].astype(np.float32),
        'cell_type_names': model.cell_type_names,
    }
    if 'embedding' in results:
        ret['embedding'] = results['embedding'].astype(np.float32)
    if 'sigmatrix' in results:
        ret['sigmatrix'] = results['sigmatrix'].astype(np.float32)
        ret['var_genes'] = list(model.var_genes) if model.var_genes else []
    return ret


# ═════════════════════════════════════════════════════════════════════════════
# 8. Utilities
# ═════════════════════════════════════════════════════════════════════════════

def deconv_get_history(model):
    """Return training history as dict (for R data.frame conversion).

    Returns:
        dict with train_loss, val_loss, train_corr, val_corr arrays.
    """
    h = model.history
    return {
        'train_loss': np.array(h.get('train_loss', []), dtype=np.float32),
        'val_loss': np.array(h.get('val_loss', []), dtype=np.float32),
        'train_corr': np.array(h.get('train_corr', []), dtype=np.float32),
        'val_corr': np.array(h.get('val_corr', []), dtype=np.float32),
    }
