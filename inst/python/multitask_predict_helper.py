"""
multitask_predict_helper.py
────────────────────────────
Python functions exposed to R via reticulate::source_python().

All functions return plain Python dicts / numpy arrays so that
reticulate::py_to_r() can convert them to R lists / matrices directly.

Prediction usage (R side):
    library(reticulate)
    use_condaenv("/home/oyh/miniforge3/envs/scanpy-env", required = TRUE)
    source_python("multitask_predict_helper.py")
    model  <- mt_load_model("models/multi_task_model.pt")
    query  <- mt_load_query("../data/toy_test.h5ad")
    X      <- mt_prepare_query(query, model)
    out    <- mt_predict(model, query, X)

Training usage (R side):
    adata  <- mt_load_query("reference.h5ad")         # or mt_srt_to_adata(...)
    hvg    <- mt_select_hvg(adata, 6000L, "sample")
    model  <- mt_create_model(length(hvg), 18L)
    mt_train(model, adata, "cell_type", "umap", hvg)
    mt_set_ood(model, 0.5)   # optional: sets is_ood advisory flag threshold
    mt_save_model(model, "models/model.pt")
"""

import os
import sys
import numpy as np

# ── Ensure multi_task_model.py is importable ──────────────────────────────────
try:
    _dir = os.path.dirname(os.path.abspath(__file__))
except NameError:
    _dir = os.getcwd()
if _dir not in sys.path:
    sys.path.insert(0, _dir)

import scanpy as sc
import h5py
from multi_task_model import MultiTaskModel


# ═══════════════════════════════════════════════════════════════════════════════
# 1. Model & Data Loading
# ═══════════════════════════════════════════════════════════════════════════════

def mt_load_model(model_path: str):
    """Load a trained MultiTaskModel from a .pt file.

    Args:
        model_path: Path to the saved model (.pt).

    Returns:
        Loaded MultiTaskModel instance.
    """
    model = MultiTaskModel.load(model_path)
    return model


def mt_load_query(h5ad_path: str):
    """Load AnnData query from an h5ad file.

    Args:
        h5ad_path: Path to the .h5ad query file.

    Returns:
        AnnData object.
    """
    query = sc.read_h5ad(h5ad_path)
    return query


# ═══════════════════════════════════════════════════════════════════════════════
# 2. Feature Matrix Preparation
# ═══════════════════════════════════════════════════════════════════════════════

def mt_prepare_query(query, model) -> np.ndarray:
    """Prepare binary feature matrix aligned to model's variable genes.

    Handles missing genes by zero-filling (cells with no expression for those
    genes are treated as unexpressed). Uses O(n) numpy fancy indexing.

    Args:
        query: AnnData object.
        model: Loaded MultiTaskModel.

    Returns:
        X: np.ndarray of shape (n_cells, n_genes), dtype float32, values 0/1.
    """
    var_genes     = model.get_feature_names()
    query_gene_set = set(query.var_names)
    shared_vars   = [g for g in var_genes if g in query_gene_set]
    n_missing     = len(var_genes) - len(shared_vars)

    if n_missing > 0:
        print(f"  {n_missing}/{len(var_genes)} variable genes missing "
              f"→ zero-filled ({len(shared_vars)} shared genes used).")
    else:
        print(f"  All {len(var_genes)} variable genes found in query.")

    if not shared_vars:
        raise ValueError("No shared genes found between query and model. "
                         "Check the input .h5ad file.")

    X_sub = query[:, shared_vars].X
    if hasattr(X_sub, 'toarray'):
        X_sub = X_sub.toarray()
    X_sub = (X_sub > 0).astype(np.float32)

    if n_missing == 0:
        return X_sub  # fast path

    # O(n) fancy indexing: place shared columns in correct positions
    shared_pos = np.array([i for i, g in enumerate(var_genes) if g in query_gene_set])
    X_full = np.zeros((query.n_obs, len(var_genes)), dtype=np.float32)
    X_full[:, shared_pos] = X_sub
    return X_full


# ═══════════════════════════════════════════════════════════════════════════════
# 3. Prediction
# ═══════════════════════════════════════════════════════════════════════════════

def mt_predict(model, query, X: np.ndarray) -> dict:
    """Run MultiTaskModel prediction and return R-friendly output dict.

    Args:
        model: Loaded MultiTaskModel.
        query: AnnData object (provides cell IDs and optional true UMAP).
        X:     Feature matrix from mt_prepare_query().

    Returns:
        dict with keys:
          'cell_level'       - dict of 1-D arrays (cell_id, cell_type_pred,
                               confidence, is_ood, umap_1_pred, umap_2_pred)
          'shared_embedding' - np.ndarray (n_cells × 512)
    """
    results  = model.predict(X)
    cell_ids = np.array(list(query.obs_names.astype(str)))

    cell_level = {
        'cell_id':        cell_ids,
        'cell_type_pred': np.array(results['cell_types']),
        'confidence':     results['confidence'].astype(np.float32),
        'is_ood':         results['is_ood'],          # bool array → R logical
        'umap_1_pred':    results['embeddings'][:, 0].astype(np.float32),
        'umap_2_pred':    results['embeddings'][:, 1].astype(np.float32),
    }

    return {
        'cell_level':       cell_level,
        'shared_embedding': results['shared_embedding'].astype(np.float32),
    }


# ═══════════════════════════════════════════════════════════════════════════════
# 4. Gene Importance (Integrated Gradients)
# ═══════════════════════════════════════════════════════════════════════════════

def mt_explain_global(model, X: np.ndarray,
                      task: str = 'cls',
                      top_k: int = 15,
                      n_cells: int = 50) -> dict:
    """Compute global gene importance via Integrated Gradients.

    Args:
        model:   Loaded MultiTaskModel.
        X:       Feature matrix (full dataset).
        task:    'cls' (classification) or 'reg' (UMAP regression).
        top_k:   Number of top genes to return.
        n_cells: Number of cells to use as attribution baseline (speed tradeoff).

    Returns:
        dict with 'gene' (str array) and 'importance' (float32 array).
    """
    X_sub  = X[:int(n_cells)]
    result = model.explain_top_genes(X_sub, task=task,
                                     top_k=int(top_k), n_steps=50)
    return {
        'gene':       np.array(result['gene_names']),
        'importance': result['importance'].astype(np.float32),
    }


def mt_explain_per_class(model, X: np.ndarray,
                          labels,
                          top_k: int = 10) -> dict:
    """Compute per-class gene importance via Integrated Gradients.

    Args:
        model:  Loaded MultiTaskModel.
        X:      Feature matrix.
        labels: Cell type labels (R character vector → Python list/array).
        top_k:  Top genes per class.

    Returns:
        Flat dict with arrays: cell_type, n_cells, rank, gene, importance.
        (Flat structure is easier to convert to R data.frame.)
    """
    labels_arr = np.array(labels, dtype=str)
    per_class  = model.explain_per_class(
        X, labels=labels_arr,
        top_k=int(top_k), n_steps=50, max_cells_per_class=100
    )

    # Flatten Python dict-of-dicts to parallel arrays
    cls_names   = []
    n_cells_v   = []
    ranks       = []
    genes       = []
    importances = []

    for cls_name, info in per_class.items():
        for rank, (g, s) in enumerate(
                zip(info['gene_names'], info['importance']), start=1):
            cls_names.append(cls_name)
            n_cells_v.append(int(info['n_cells']))
            ranks.append(rank)
            genes.append(g)
            importances.append(float(s))

    return {
        'cell_type':  np.array(cls_names),
        'n_cells':    np.array(n_cells_v,   dtype=np.int32),
        'rank':       np.array(ranks,       dtype=np.int32),
        'gene':       np.array(genes),
        'importance': np.array(importances, dtype=np.float32),
    }


# ═══════════════════════════════════════════════════════════════════════════════
# 5. Save Results (HDF5 + CSV)
# ═══════════════════════════════════════════════════════════════════════════════

def mt_save_h5(cell_level_dict: dict, per_class_dict: dict, path: str) -> None:
    """Save results to R-readable HDF5 file (no pytables required).

    h5py 3.x requires string_dtype + Python list for variable-length strings.
    Numeric numpy arrays are written as-is.

    R read example:
        library(rhdf5)
        cell_df   <- as.data.frame(lapply(h5read(path, "cell_level"), identity))
        pclass_df <- as.data.frame(lapply(h5read(path, "per_class_df"), identity))

    Args:
        cell_level_dict: dict of 1-D arrays (columns of cell-level table).
        per_class_dict:  dict of 1-D arrays (columns of per-class table).
        path:            Output .h5 file path.
    """
    str_dt = h5py.string_dtype(encoding='utf-8')

    def _write_dataset(group, key, arr):
        arr = np.array(arr)
        if arr.dtype.kind in ('U', 'O', 'S'):           # string/object → list[str]
            group.create_dataset(key,
                                 data=list(arr.astype(str)),
                                 dtype=str_dt)
        elif arr.dtype.kind == 'b':                      # bool → int8
            group.create_dataset(key, data=arr.astype(np.int8))
        else:
            group.create_dataset(key, data=arr)

    out_dir = os.path.dirname(path)
    if out_dir:
        os.makedirs(out_dir, exist_ok=True)

    with h5py.File(path, 'w') as f:
        g = f.create_group('cell_level')
        for k, v in cell_level_dict.items():
            _write_dataset(g, k, v)

        p = f.create_group('per_class_df')
        for k, v in per_class_dict.items():
            _write_dataset(p, k, v)


# ═══════════════════════════════════════════════════════════════════════════════
# 6. Seurat → AnnData Conversion  (R Seurat object → Python AnnData)
# ═══════════════════════════════════════════════════════════════════════════════

def mt_srt_to_adata(counts_matrix, obs_df, var_names, obs_names,
                     obsm_dict=None):
    """Build an AnnData object from R Seurat components.

    Called when the user passes a Seurat object instead of an h5ad path.
    R side extracts counts matrix, metadata, and embeddings, then passes
    them to this function via reticulate.

    Args:
        counts_matrix: scipy.sparse.csr_matrix or np.ndarray (genes × cells).
                       NOTE: R Matrix is genes×cells; this function transposes.
        obs_df:        pandas DataFrame of cell metadata (from srt@meta.data).
        var_names:     list/array of gene names (rownames of counts matrix).
        obs_names:     list/array of cell barcodes (colnames of counts matrix).
        obsm_dict:     Optional dict of embeddings, e.g. {'umap': np.ndarray}.
                       Each value shape: (n_cells, n_dims).

    Returns:
        AnnData object ready for mt_select_hvg / mt_train.
    """
    import scipy.sparse as sp
    import pandas as pd

    # R stores genes×cells; AnnData wants cells×genes
    if sp.issparse(counts_matrix):
        X = counts_matrix.T.tocsr()
    else:
        X = np.array(counts_matrix).T

    var_names = list(np.array(var_names, dtype=str))
    obs_names = list(np.array(obs_names, dtype=str))

    adata = sc.AnnData(
        X=X,
        obs=pd.DataFrame(obs_df, index=obs_names) if obs_df is not None
            else pd.DataFrame(index=obs_names),
    )
    adata.var_names = var_names

    if obsm_dict is not None:
        for key, val in obsm_dict.items():
            adata.obsm[str(key)] = np.array(val, dtype=np.float32)

    print(f"  AnnData built: {adata.n_obs} cells × {adata.n_vars} genes")
    return adata


# ═══════════════════════════════════════════════════════════════════════════════
# 7. Training API
# ═══════════════════════════════════════════════════════════════════════════════

def mt_select_hvg(adata, n_top_genes: int = 6000,
                   batch_key: str = None) -> list:
    """Select highly variable genes using scanpy (seurat_v3 method).

    Modifies adata in-place (adds 'highly_variable' to var).

    Args:
        adata:       AnnData object with raw counts in .X.
        n_top_genes: Number of top variable genes to select.
        batch_key:   obs column for batch-aware HVG selection (optional).

    Returns:
        List of HVG gene names.
    """
    import copy as _copy
    adata_work = adata.copy()

    # Ensure counts are stored in a layer for seurat_v3
    if 'counts' not in adata_work.layers:
        adata_work.layers['counts'] = adata_work.X.copy()

    sc.pp.highly_variable_genes(
        adata_work,
        flavor='seurat_v3',
        n_top_genes=int(n_top_genes),
        layer='counts',
        batch_key=batch_key,
        subset=False,
    )
    hvg = list(adata_work.var_names[adata_work.var['highly_variable']])
    print(f"  HVG selected: {len(hvg)} genes"
          + (f" (batch_key={batch_key})" if batch_key else ""))
    return hvg


def mt_create_model(input_size: int,
                     num_classes: int,
                     embedding_dim: int = 2,
                     hidden_size: int = 512,
                     head_hidden_size: int = 256,
                     n_resnet_blocks: int = 4,
                     input_dropout_rate: float = 0.25,
                     resnet_dropout_rate: float = 0.1,
                     head_dropout_rate: float = 0.1,
                     batch_size: int = 256,
                     test_size: float = 0.2,
                     num_epochs: int = 50,
                     early_stopping_patience: int = 10,
                     early_stopping_monitor: str = 'cls',
                     learning_rate: float = 1e-3,
                     use_gradnorm: bool = True,
                     gradnorm_alpha: float = 0.15,
                     label_smoothing: float = 0.1,
                     use_lr_scheduler: bool = True,
                     use_uncertainty: bool = False,
                     random_state: int = 42):
    """Create a new (untrained) MultiTaskModel.

    Args mirror MultiTaskModel.__init__; see multi_task_model.py for details.

    Returns:
        MultiTaskModel instance (untrained).
    """
    model = MultiTaskModel(
        input_size=int(input_size),
        num_classes=int(num_classes),
        embedding_dim=int(embedding_dim),
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
        early_stopping_monitor=str(early_stopping_monitor),
        learning_rate=float(learning_rate),
        use_gradnorm=bool(use_gradnorm),
        gradnorm_alpha=float(gradnorm_alpha),
        label_smoothing=float(label_smoothing),
        use_lr_scheduler=bool(use_lr_scheduler),
        use_uncertainty=bool(use_uncertainty),
        random_state=int(random_state),
    )
    return model


def mt_train(model, adata, label_column: str, embedding_key: str,
              var_genes: list) -> dict:
    """Train the model on an AnnData reference dataset.

    Args:
        model:         MultiTaskModel instance (from mt_create_model).
        adata:         AnnData object with raw counts, obs metadata, and obsm.
        label_column:  obs column for cell type labels.
        embedding_key: obsm key for embedding coordinates (e.g. 'umap').
        var_genes:     List of variable gene names to use.

    Returns:
        dict with training summary: best_val_acc, best_val_rmse, n_epochs.
    """
    var_genes = list(np.array(var_genes, dtype=str))

    model.fit(
        adata=adata,
        label_column=label_column,
        embedding_key=embedding_key,
        var_genes=var_genes,
    )
    return {
        'best_val_acc':  float(model.best_val_acc) if model.best_val_acc else None,
        'best_val_rmse': float(model.best_val_rmse) if model.best_val_rmse else None,
        'n_epochs':      len(model.history.get('train_loss', [])),
    }


def mt_build_ref_index(model, adata, embedding_key: str = 'umap',
                        label_column: str = None,
                        max_ref_cells: int = 20000,
                        n_per_type: int = 1000) -> None:
    """Build reference index for kNN UMAP projection on an existing model.

    After calling this, use mt_save_model() to persist the reference index.

    Args:
        model:          Trained MultiTaskModel.
        adata:          Reference AnnData with UMAP in obsm.
        embedding_key:  obsm key for UMAP coordinates.
        label_column:   obs column for stratified subsampling.
        max_ref_cells:  Max reference cells to store.
        n_per_type:     Max cells per cell type.
    """
    model.build_reference_index(
        adata=adata,
        embedding_key=str(embedding_key),
        label_column=label_column,
        max_ref_cells=int(max_ref_cells),
        n_per_type=int(n_per_type),
    )
    n_ref = model._ref_emb.shape[0] if model._ref_emb is not None else 0
    print(f"  Reference index built: {n_ref} cells stored")


def mt_set_ood(model, threshold: float) -> None:
    """Set OOD detection threshold for the is_ood advisory flag.

    Cell type predictions always retain the best-guess label.
    The threshold only affects the is_ood boolean flag, allowing
    users to filter low-confidence cells downstream.

    Args:
        model:     Trained MultiTaskModel.
        threshold: Confidence cutoff (cells below flagged is_ood=True).
    """
    model.set_ood_threshold(float(threshold))


def mt_save_model(model, path: str) -> None:
    """Save trained model to a .pt file.

    Args:
        model: Trained MultiTaskModel.
        path:  Output file path (e.g. 'models/model.pt').
    """
    out_dir = os.path.dirname(path)
    if out_dir:
        os.makedirs(out_dir, exist_ok=True)
    model.save(path)


def mt_get_history(model) -> dict:
    """Return training history as dict of lists (→ R data.frame).

    Returns:
        dict with keys: epoch, train_loss, val_loss, train_acc, val_acc,
        train_rmse, val_rmse, w_cls, w_reg, lr.
    """
    h = model.history
    if not h:
        return {}
    n = len(h.get('train_loss', []))
    result = {'epoch': list(range(1, n + 1))}
    for key in ['train_loss', 'val_loss',
                'train_cls_loss', 'val_cls_loss',
                'train_reg_loss', 'val_reg_loss',
                'train_acc', 'val_acc',
                'train_rmse', 'val_rmse',
                'w_cls', 'w_reg', 'lr']:
        if key in h and len(h[key]) > 0:
            result[key] = [float(v) for v in h[key]]
    return result


def mt_fine_tune(model, adata, label_column: str, embedding_key: str,
                  var_genes: list,
                  freeze: bool = True,
                  new_num_epochs: int = None,
                  new_learning_rate: float = None) -> dict:
    """Fine-tune a trained model on new data (transfer learning).

    Args:
        model:             Trained MultiTaskModel.
        adata:             New AnnData dataset.
        label_column:      obs column for cell type labels.
        embedding_key:     obsm key for embedding coordinates.
        var_genes:         List of variable gene names.
        freeze:            If True, freeze backbone (train heads only).
        new_num_epochs:    Override number of epochs (default: half original).
        new_learning_rate: Override learning rate (default: 1/10 original).

    Returns:
        dict with fine-tuning summary.
    """
    var_genes = list(np.array(var_genes, dtype=str))

    kwargs = dict(
        adata=adata,
        label_column=label_column,
        embedding_key=embedding_key,
        var_genes=var_genes,
        freeze=bool(freeze),
    )
    if new_num_epochs is not None:
        kwargs['new_num_epochs'] = int(new_num_epochs)
    if new_learning_rate is not None:
        kwargs['new_learning_rate'] = float(new_learning_rate)

    model.fine_tune(**kwargs)
    return {
        'best_val_acc':  float(model.best_val_acc) if model.best_val_acc else None,
        'best_val_rmse': float(model.best_val_rmse) if model.best_val_rmse else None,
        'n_epochs':      len(model.history.get('train_loss', [])),
    }


# ═══════════════════════════════════════════════════════════════════════════════
# 8. Pathway Scoring (per-cell IG → GMT aggregation)
# ═══════════════════════════════════════════════════════════════════════════════

def mt_parse_gmt(gmt_path: str) -> dict:
    """Parse a GMT (Gene Matrix Transposed) file.

    GMT format: each line = pathway_name <TAB> description <TAB> gene1 <TAB> gene2 ...

    Args:
        gmt_path: Path to the GMT file.

    Returns:
        dict mapping pathway name → list of gene symbols.
    """
    pathways = {}
    with open(gmt_path, 'r', encoding='utf-8') as f:
        for line in f:
            parts = line.strip().split('\t')
            if len(parts) < 3:
                continue
            name = parts[0]
            genes = [g for g in parts[2:] if g.strip()]
            if genes:
                pathways[name] = genes
    print(f"  GMT parsed: {len(pathways)} pathways from {os.path.basename(gmt_path)}")
    return pathways


def mt_pathway_score(model, X: np.ndarray, gmt_path: str,
                     task: str = 'cls',
                     min_genes: int = 5,
                     n_steps: int = 50,
                     batch_size: int = 200) -> dict:
    """Compute per-cell pathway importance via IG + GMT aggregation.

    Pipeline:
      1. Compute per-cell Integrated Gradients attributions [n_cells × n_genes]
      2. Parse GMT to get gene-pathway membership
      3. Build binary mask matrix M [n_genes × n_pathways]
      4. scores = |attributions| @ M / pathway_sizes

    The resulting scores represent each pathway's discriminative contribution
    to the model's prediction for each individual cell — fundamentally different
    from expression-based scores (AUCell/GSVA/ssGSEA).

    Args:
        model:      Trained MultiTaskModel.
        X:          Binary feature matrix (n_cells × n_genes).
        gmt_path:   Path to GMT file.
        task:       'cls' (classification) or 'reg' (regression).
        min_genes:  Minimum gene overlap for a pathway to be retained.
        n_steps:    IG interpolation steps.
        batch_size: Cells per IG batch (memory tradeoff).

    Returns:
        dict with:
          'scores':        np.ndarray (n_cells × n_pathways), float32
          'pathway_names': list of pathway names
          'pathway_sizes': np.ndarray (n_pathways,), int — genes per pathway
    """
    # 1. Parse GMT
    pathways = mt_parse_gmt(gmt_path)

    # 2. Get model gene names and build index
    gene_names = model.get_feature_names()
    gene_set = set(gene_names)
    gene_to_idx = {g: i for i, g in enumerate(gene_names)}

    # 3. Filter pathways by overlap with model genes
    filtered = {}
    for pw_name, pw_genes in pathways.items():
        overlap = [g for g in pw_genes if g in gene_set]
        if len(overlap) >= min_genes:
            filtered[pw_name] = overlap

    if not filtered:
        raise ValueError(
            f"No pathways with >= {min_genes} overlapping genes found. "
            f"Model has {len(gene_names)} genes; GMT has {len(pathways)} pathways. "
            f"Try lowering min_genes or using a different GMT file."
        )

    print(f"  Pathways retained: {len(filtered)}/{len(pathways)} "
          f"(min_genes={min_genes})")

    # 4. Build binary mask matrix [n_genes × n_pathways]
    pw_names = sorted(filtered.keys())
    n_genes = len(gene_names)
    n_pw = len(pw_names)
    mask = np.zeros((n_genes, n_pw), dtype=np.float32)
    pw_sizes = np.zeros(n_pw, dtype=np.float32)

    for j, pw_name in enumerate(pw_names):
        for g in filtered[pw_name]:
            mask[gene_to_idx[g], j] = 1.0
        pw_sizes[j] = mask[:, j].sum()

    # 5. Compute IG in batches and aggregate to pathway scores
    n_cells = X.shape[0]
    scores = np.zeros((n_cells, n_pw), dtype=np.float32)
    n_batches = (n_cells + batch_size - 1) // batch_size

    print(f"  Computing IG for {n_cells} cells ({n_batches} batches) ...")
    for batch_idx, start in enumerate(range(0, n_cells, batch_size)):
        end = min(start + batch_size, n_cells)
        X_batch = X[start:end]

        # Per-cell IG attributions [batch_size × n_genes]
        attr_batch = model.explain(X_batch, task=task, n_steps=n_steps)

        # |attr| @ mask / sizes → [batch_size × n_pathways]
        scores[start:end] = np.abs(attr_batch) @ mask / pw_sizes[np.newaxis, :]

        print(f"    Batch {batch_idx + 1}/{n_batches} "
              f"({end}/{n_cells} cells)")

    print(f"  Pathway scoring complete: {n_cells} cells × {n_pw} pathways")

    return {
        'scores':        scores,
        'pathway_names': pw_names,
        'pathway_sizes': pw_sizes.astype(np.int32),
    }


# ═══════════════════════════════════════════════════════════════════════════════
# 9. Device Management (CPU / GPU)
# ═══════════════════════════════════════════════════════════════════════════════

def mt_set_device(device_str: str = "auto") -> str:
    """Set the compute device for model training and prediction.

    Must be called BEFORE mt_create_model() or mt_load_model().

    Args:
        device_str: "auto" (default, let PyTorch detect),
                    "cpu"  (force CPU even if GPU available),
                    "cuda" (require GPU, raises error if unavailable).

    Returns:
        str: The actual device that will be used ("cpu" or "cuda").
    """
    import torch

    if device_str == "cpu":
        os.environ["CUDA_VISIBLE_DEVICES"] = ""
        print("  Device forced to CPU (CUDA_VISIBLE_DEVICES='')")
        return "cpu"
    elif device_str == "cuda":
        # Remove any CPU-forcing override
        os.environ.pop("CUDA_VISIBLE_DEVICES", None)
        if not torch.cuda.is_available():
            raise RuntimeError(
                "device='cuda' requested but torch.cuda.is_available() is False. "
                "Install CUDA-enabled PyTorch or use device='auto'/'cpu'."
            )
        print(f"  Device: cuda ({torch.cuda.get_device_name(0)})")
        return "cuda"
    else:  # "auto"
        os.environ.pop("CUDA_VISIBLE_DEVICES", None)
        actual = "cuda" if torch.cuda.is_available() else "cpu"
        print(f"  Device (auto-detected): {actual}")
        return actual
