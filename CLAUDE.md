# scMMR — Single-Cell Multi-Task Model in R

## 项目概述

scMMR 是一个基于 PyTorch 深度学习的 R 包，通过 `reticulate` 桥接 R 与 Python。核心功能包括：
- **多任务 DNN**：联合学习细胞类型分类 + UMAP 嵌入回归（共享 ResNet 骨干）
- **Bulk 反卷积**：从 scRNA-seq 参考训练 DNN，估计 bulk RNA-seq 中的细胞类型比例
- **可解释性**：Integrated Gradients 基因/通路重要性分析
- **下游分析**：差异丰度检验、扰动排名、通路分析、基因集评分、CytoTRACE2

**版本**: 0.2.0
**作者**: Yuhao Ouyang
**License**: MIT

---

## 环境配置

| 项目 | 值 |
|------|-----|
| OS | WSL Ubuntu-22.04 on Windows |
| R | 4.4.3 |
| Conda 环境 | `scanpy-env` (`/home/oyh/miniforge3/envs/scanpy-env`) |
| Python 要求 | >= 3.8, torch (CPU/CUDA), scanpy, anndata, numpy, pandas, scipy, h5py |
| 包路径 | `/home/oyh/project/para/python_pacakge/scMMR/` |

**R 中设置 Python**：
```r
use_scMMR_python(condaenv = "/home/oyh/miniforge3/envs/scanpy-env")
```

**WSL 中执行 R 命令**：
```bash
wsl -d Ubuntu-22.04 -- Rscript -e '...'
wsl -d Ubuntu-22.04 -- bash -c '...'
```

---

## 包构建与安装

```r
# 每次修改 R 文件后必须执行：
devtools::document()          # 从 roxygen 生成 man/ 和 NAMESPACE
devtools::install(quick = TRUE)  # 安装到 R library
```

**注意事项**：
- 所有导出函数必须有 `@export` roxygen 标签
- 修改 R 文件后必须重新 `document()` + `install()`
- 新增 R 文件需在 `DESCRIPTION` 的 `Collate` 字段中添加（按依赖顺序）
- Python 文件修改后无需重新安装，但需重启 R session 或重新 source

---

## 目录结构

```
scMMR/
├── DESCRIPTION              # 包元数据、依赖、Collate 顺序
├── NAMESPACE                # 导出函数列表 (roxygen2 自动生成)
├── CLAUDE.md                # 本文件
├── R/                       # R 源代码 (25+ 文件)
│   ├── zzz.R                # 包初始化 (.onLoad, .scMMR_env)
│   ├── python_env.R         # Python 环境管理
│   ├── utils.R              # 内部辅助函数
│   ├── DNN_train.R          # 多任务 DNN 训练
│   ├── DNN_predict.R        # 多任务 DNN 预测
│   ├── DNN_deconv_train.R   # 反卷积 DNN 训练
│   ├── DNN_deconv_predict.R # 反卷积 DNN 预测
│   ├── plot_utils.R         # 绘图辅助函数 (.percentage_stat, .rank_scatter_parse, .extract_cellmeta)
│   ├── PlotAlluvia.R        # 冲积图
│   ├── PlotRankScatter.R    # 排名散点图
│   ├── PlotImportance.R     # 基因重要性柱状图
│   ├── PlotRoe.R            # O/E 比热力图 (含 p 值 + 分级符号)
│   ├── PlotMAP.R            # UMAP 对齐投影图
│   ├── PlotPerturbation.R   # 扰动检验结果图
│   ├── PlotPercent.R        # 差异丰度结果图
│   ├── PlotCorrelation.R    # 相关性热力图
│   ├── PlotPropCorrelation.R # 比例相关散点图
│   ├── PlotScatter.R        # 通用散点图
│   ├── PlotAnnotation.R     # 注释图
│   ├── PlotSankey.R         # 桑基图
│   ├── rank.R               # 扰动排名 + 差异丰度 (950 行)
│   ├── evaluate.R           # 嵌入质量评估 (660 行)
│   ├── compute_module_score.R # 基因集评分 (AUCell/Seurat/UCell)
│   ├── pathway_analysis.R   # GSEA/SCPA 通路分析 (738 行)
│   └── cytotrace.R          # CytoTRACE2 发育轨迹
├── inst/
│   ├── python/              # PyTorch 模型 + 桥接函数 (4 文件, ~3,143 行)
│   │   ├── multi_task_model.py       # MultiTaskModel 核心架构
│   │   ├── multitask_predict_helper.py # 多任务 R-Python 桥接
│   │   ├── deconv_model.py           # DeconvModel 反卷积架构
│   │   └── deconv_helper.py          # 反卷积 R-Python 桥接
│   ├── demo/                # 示例脚本 (9 个)
│   └── extdata/             # 预训练模型 + 测试数据 + GMT 文件
│       ├── model.pt         # 预训练多任务模型 (19 MB)
│       ├── toy_test.qs      # 测试 Seurat 对象 (19 MB)
│       ├── ref_umap.qs      # 参考 UMAP 坐标 (2.6 MB)
│       └── gmt/             # 16 个通路数据库文件
├── man/                     # roxygen2 自动生成的帮助文档
└── tests/                   # (暂无正式测试)
```

---

## 导出函数一览 (24 个)

### 核心训练与预测

| 函数 | 用途 | 输入 | 输出 |
|------|------|------|------|
| `DNN_train()` | 多任务 DNN 训练 | h5ad / Seurat | model, var_genes, history |
| `DNN_predict()` | 细胞类型预测 + 可解释性 | h5ad / Seurat + model.pt | predictions, embedding, importance |
| `DNN_deconv_train()` | 反卷积模型训练 | h5ad / Seurat (scRNA-seq 参考) | model, var_genes, cell_types |
| `DNN_deconv_predict()` | Bulk 反卷积预测 | matrix / CSV / h5ad / Seurat | proportions data.frame |

### 环境管理

| 函数 | 用途 |
|------|------|
| `install_scMMR_python()` | 创建 conda/virtualenv 并安装依赖 |
| `use_scMMR_python()` | 指定使用的 Python 环境 |

### 基因集与通路分析

| 函数 | 用途 |
|------|------|
| `read_gmt()` | 解析 GMT 文件 → named list |
| `parse_gene_sets()` | 统一转换基因集格式 (list / data.frame / GMT) |
| `ComputeModuleScore()` | 基因集评分 (S3 泛型: matrix / Seurat; AUCell / Seurat / UCell) |
| `RunPathwayAnalysis()` | GSEA 或 SCPA 通路分析 |
| `PlotPathwayBubble()` | 通路分析气泡图 |

### 排名与差异检验

| 函数 | 用途 |
|------|------|
| `RankPerturbation()` | 扰动排名 (Wasserstein / MMD / Energy / AUC) |
| `RankPercent()` | 差异丰度检验 (置换检验) |

### 可视化 (7 个 Plot 函数)

| 函数 | 用途 |
|------|------|
| `PlotAlluvia()` | 冲积图 (细胞组成变化) |
| `PlotRankScatter()` | 排名散点图 (基因/通路/regulon) |
| `PlotImportance()` | 基因重要性柱状图 |
| `PlotRoe()` | O/E 比热力图 (含 p 值 + 分级符号) |
| `PlotMAP()` | UMAP 对齐投影图 |
| `PlotPerturbation()` | 扰动检验结果图 |
| `PlotPercent()` | 差异丰度结果图 (beeswarm) |

### 嵌入评估

| 函数 | 用途 |
|------|------|
| `EvaluateEmbedding()` | DNN embedding vs PCA 质量评估 |
| `PlotEmbeddingEval()` | 评估结果 5 图 (elbow/KNN/距离相关) |

### CytoTRACE2

| 函数 | 用途 |
|------|------|
| `RunCytoTRACE2()` | 发育轨迹推断 (S3 泛型) |
| `PlotCytoTRACE2()` | CytoTRACE2 小提琴图 |

---

## 内部函数 (`.xxx` 前缀，不导出)

| 函数 | 文件 | 用途 |
|------|------|------|
| `.ensure_python()` | python_env.R | 惰性加载 Python helper (仅首次调用时 source) |
| `.resolve_input()` | utils.R | h5ad 路径或 Seurat → AnnData |
| `.seurat_to_adata()` | utils.R | Seurat → AnnData (counts + meta + reductions) |
| `.select_hvg()` | utils.R | HVG 选择 (batch-aware, tryCatch 回退) |
| `.prepare_bulk_input()` | utils.R | 多格式 bulk 输入 → matrix (genes×samples) |
| `.percentage_stat()` | plot_utils.R | 计算细胞群体百分比 |
| `.rank_scatter_parse()` | plot_utils.R | 解析多种数据格式用于散点图 |
| `.extract_cellmeta()` | plot_utils.R | 从 Seurat/data.frame 提取 metadata |
| `.validate_embedding()` | rank.R | 验证并对齐 embedding 与 metadata |
| `.wasserstein_1d()` | rank.R | 1D Wasserstein 距离 (排序法) |
| `.sliced_wasserstein()` | rank.R | 高维 Sliced Wasserstein (随机投影) |
| `.mmd_rbf()` | rank.R | MMD (多尺度 RBF 核) |
| `.energy_distance()` | rank.R | Energy distance (无参数) |
| `.auc_lda()` | rank.R | AUC via LDA (带 PCA 回退) |
| `.permutation_test()` | rank.R | 置换检验生成 p 值 |
| `.find_knee()` | evaluate.R | 方差解释率拐点检测 |
| `.knn_jaccard()` | evaluate.R | KNN Jaccard 重叠度 |
| `.dist_correlation()` | evaluate.R | 成对距离 Spearman 秩相关 |
| `.rv_coefficient_fast()` | evaluate.R | RV 系数 (trace trick) |
| `.score_aucell()` | compute_module_score.R | AUCell 批量评分 |
| `.score_seurat()` | compute_module_score.R | Seurat AddModuleScore 方法 |
| `.score_ucell()` | compute_module_score.R | UCell 签名评分 |
| `.clean_pathway_names()` | pathway_analysis.R | 去除通路前缀，Title Case |
| `.gs_to_term2gene()` | pathway_analysis.R | 转 clusterProfiler 格式 |
| `.run_gsea_across_celltypes()` | pathway_analysis.R | 逐细胞类型 GSEA |
| `.run_scpa_across_celltypes()` | pathway_analysis.R | 逐细胞类型 SCPA |

---

## Python 模块架构

### multi_task_model.py — 多任务 DNN (1,490 行)

```
输入 (binary: X > 0)
  ↓
InputBlock: BatchNorm → Dropout(0.25) → Linear(n_genes → 512) → SELU
  ↓
ResNetBlock ×4: Dropout(0.1) → Linear(512 → 512) → SELU + skip connection
  ↓
ResNetLastBlock: Dropout(0.1) → Linear(512 → 512) → SELU (无 skip)
  ↓ (512-dim shared embedding)
  ├─→ ClassificationHead: Linear(512→256) → SELU → Dropout → Linear(256→K)
  └─→ RegressionHead: Linear(512→256) → SELU → Dropout → Linear(256→D)
```

**关键特性**：
- GradNorm 多任务损失平衡
- Label smoothing 正则化
- OOD 检测 (置信度阈值)
- Integrated Gradients 可解释性
- CosineAnnealingWarmRestarts 学习率调度
- 支持 fine-tune 和 transfer learning

### deconv_model.py — 反卷积 DNN (587 行)

```
输入 (continuous: log1p(CPM))
  ↓
InputBlock → ResNetBlock ×4 → ResNetLastBlock
  ↓ (512-dim embedding)
  ↓
TaskHead: Linear(512→256) → SELU → Dropout → Linear(256→K)
  ↓
Softmax → proportions (sum = 1)
```

**关键差异** (vs 多任务 DNN)：
- 输入：连续 log-CPM（**不是**二值化 `X > 0`）
- 输出：细胞类型比例向量（softmax 归一化）
- 损失：MSE 或 KL 散度
- 训练数据：Dirichlet(α=1) 采样生成的 pseudo-bulk

### R-Python 桥接模式

```
R 函数 → .ensure_python() → .scMMR_env$py_func() → Python 执行 → 返回 dict/numpy → reticulate::py_to_r()
```

- `multitask_predict_helper.py`：`mt_*` 系列函数 (18 个)
- `deconv_helper.py`：`deconv_*` 系列函数 (9 个)
- 所有桥接函数返回 plain dict/numpy，便于 reticulate 自动转换

---

## 数据文件

### 内置数据 (`inst/extdata/`)

| 文件 | 大小 | 说明 |
|------|------|------|
| `model.pt` | 19 MB | 预训练多任务模型 (甲状旁腺 scRNA-seq) |
| `toy_test.qs` | 19 MB | 测试 Seurat 对象 (含 cell_type, group 列) |
| `ref_umap.qs` | 2.6 MB | 参考 UMAP 坐标 (umap_1, umap_2, cell_type) |

### GMT 通路数据库 (`inst/extdata/gmt/`)

| 文件 | 物种 | 来源 |
|------|------|------|
| `h.all.v2022.1.Hs.symbols.gmt` | Human | MSigDB Hallmark |
| `c2.cp.kegg.v2022.1.Hs.symbols.gmt` | Human | KEGG Legacy |
| `reactome.gmt` | Human | Reactome |
| `GO_bp.gmt` | Human | GO Biological Process |
| `TF.gmt` | Human | 转录因子 target |
| `immune.gmt` | Human | 免疫相关通路 |
| `collectri.human.gmt` | Human | CollecTRI regulon |
| `collectri.human.directional.gmt` | Human | CollecTRI (方向性) |
| `progeny.human.top100.gmt` | Human | PROGENy top 100 |
| `progeny.human.top500.gmt` | Human | PROGENy top 500 |
| `proliferation.combined.gmt` | Human | 增殖相关组合 |
| `m_reactome.gmt` | Mouse | Reactome |
| `m_GO_bp.gmt` | Mouse | GO Biological Process |
| `m_TF.gmt` | Mouse | 转录因子 target |

### Demo/测试常用外部数据

| 数据集 | 来源 | 用途 | 获取方式 |
|--------|------|------|----------|
| Baron Pancreas | scRNAseq Bioconductor | 反卷积 demo、分类 benchmark | `scRNAseq::BaronPancreasData("human")` |
| Muraro Pancreas | scRNAseq Bioconductor | 跨数据集 benchmark | `scRNAseq::MuraroPancreasData()` |
| Segerstolpe Pancreas | scRNAseq Bioconductor | 跨数据集 benchmark | `scRNAseq::SegerstolpePancreasData()` |
| Xin Pancreas | scRNAseq Bioconductor | 跨数据集 benchmark | `scRNAseq::XinPancreasData()` |
| 甲状旁腺 (seu_para) | 自有数据 | 项目核心数据 | `data/` 目录 |

---

## Demo 脚本 (`inst/demo/`, 9 个)

| 文件 | 功能 | 依赖数据 |
|------|------|----------|
| `demo_predict_and_plot.R` | 完整预测+可视化工作流 | toy_test.qs, model.pt, ref_umap.qs |
| `demo_deconvolution.R` | 反卷积训练→预测→评估→可视化 | Baron Pancreas (scRNAseq) |
| `demo_benchmark_classification.R` | DNN vs SVM/RF/KNN 基准测试 | Baron Pancreas |
| `demo_cross_dataset_benchmark.R` | 跨数据集泛化评估 | 4 个胰腺数据集 |
| `demo_pathway_analysis.R` | GSEA/SCPA 通路分析 | toy_test.qs, GMT files |
| `demo_compute_module_score.R` | 基因集评分 (3 种方法) | toy_test.qs, GMT files |
| `demo_evaluate_embedding.R` | DNN embedding 质量评估 | toy_test.qs, model.pt |
| `demo_cytotrace.R` | CytoTRACE2 发育轨迹 | toy_test.qs |
| `test_gpu_speed.R` | GPU vs CPU 性能对比 | toy_test.qs, model.pt |

---

## 关键编码约定

### R 侧

1. **输入灵活性**：所有核心函数同时接受 h5ad 路径和 Seurat 对象
2. **S3 泛型**：`ComputeModuleScore()` 和 `RunCytoTRACE2()` 支持 matrix/Seurat 分发
3. **Seurat v5 兼容**：写入用 `AddMetaData()`，读取用 `@meta.data[[col]]`
4. **msigdbr API**：使用 `collection`（非 `category`），`subcollection`（非 `subcategory`），`CP:KEGG_LEGACY`（非 `CP:KEGG`）
5. **reticulate 返回值**：Python 对象用 `reticulate::py_to_r()` 转换，list 元素用 `[["key"]]` 双括号提取（避免 `$` 与 Python 对象冲突）
6. **整数参数**：传给 Python 的整数必须用 `as.integer()` 或 `L` 后缀
7. **Roxygen 标签**：导出函数用 `@export`，不用 `@keywords internal`
8. **Collate 顺序**：`DESCRIPTION` 中按依赖顺序排列（zzz → python_env → utils → ... ）

### Python 侧

1. **模块复用**：`DeconvModel` 复用 `multi_task_model.py` 的 `InputBlock`、`ResNetBlock`、`TaskHead`
2. **桥接函数命名**：多任务用 `mt_*`，反卷积用 `deconv_*`
3. **返回值**：所有桥接函数返回 plain dict/numpy（不返回 PyTorch tensor 或 Python class）
4. **设备管理**：全局变量 `_DEVICE` 通过 `set_device()` 设置，所有函数共享
5. **稀疏矩阵**：自动检测 `scipy.sparse`，必要时 `.toarray()`

### 常见陷阱

| 问题 | 原因 | 解决 |
|------|------|------|
| `cat()` 报错 `type 'list'` | `result$key` 返回 Python 对象 | 用 `result[["key"]]` 双括号提取 |
| conda 环境找不到 | 路径不对 | 使用完整路径 `/home/oyh/miniforge3/envs/scanpy-env` |
| `colData(sce)` 返回 DFrame | Bioconductor 类不兼容 | 用 `colData(sce)[["col"]]` + `which()` 索引 |
| 反卷积预测偏差 | 输入已是 log-CPM 被重复归一化 | 传入 CPM（`expm1(log_cpm)`）或 raw counts |
| HVG 选择 batch 失败 | 数据中无 batch 信息 | `tryCatch` 自动回退到无 batch 模式 |
| `subscript out-of-bounds` | 逻辑向量长度不匹配 | 用 `which()` 转索引，或 `[["col"]]` 提取向量 |

---

## 修改指南

### 新增导出函数

1. 在 `R/` 下创建或编辑 `.R` 文件
2. 添加 roxygen 注释（`@export`, `@param`, `@return`, `@examples`）
3. 在 `DESCRIPTION` 的 `Collate` 字段中按依赖顺序添加文件名
4. 运行 `devtools::document()` → 检查 `NAMESPACE` 是否更新
5. 运行 `devtools::install(quick = TRUE)`

### 新增 Python 模块

1. 在 `inst/python/` 下创建 `.py` 文件
2. 桥接函数命名：`模块前缀_功能名()`，返回 dict/numpy
3. 在 `R/python_env.R` 的 `.ensure_python()` 中添加 `source_python()`：
   ```r
   new_file <- file.path(py_path, "new_helper.py")
   if (file.exists(new_file)) {
     reticulate::source_python(new_file, envir = .scMMR_env)
   }
   ```
4. R 侧通过 `.scMMR_env$new_func()` 调用

### 修改 Python 模型架构

- `multi_task_model.py`：修改 `MultiTaskModel` 类（分类+回归任务）
- `deconv_model.py`：修改 `DeconvModel` 类（反卷积任务）
- 共享模块：`InputBlock`, `ResNetBlock`, `ResNetLastBlock`, `TaskHead` 在 `multi_task_model.py` 中定义，`deconv_model.py` 通过 `from multi_task_model import ...` 复用
- **注意**：修改共享模块会同时影响两个模型

### 运行测试

```bash
# 运行单个 demo
wsl -d Ubuntu-22.04 -- bash -c 'cd /home/oyh/project/para/python_pacakge/scMMR && Rscript inst/demo/demo_deconvolution.R'

# 快速功能测试
wsl -d Ubuntu-22.04 -- Rscript -e '
library(scMMR)
use_scMMR_python(condaenv = "/home/oyh/miniforge3/envs/scanpy-env")
# ... 测试代码
'
```

---

## R 依赖

### Imports (必需)

```
reticulate (>= 1.24), Matrix, methods, parallel, ggplot2, dplyr, rlang
```

### Suggests (按功能)

| 功能 | 包 |
|------|-----|
| Seurat 支持 | Seurat (>= 4.1.0), SeuratObject |
| HDF5 读写 | rhdf5 |
| 热力图 | ComplexHeatmap |
| 颜色 | scales, RColorBrewer, colorspace |
| KNN | FNN |
| LDA 分类器 | MASS |
| Beeswarm 图 | ggbeeswarm |
| 基因集评分 | AUCell, UCell |
| 通路分析 | clusterProfiler, SCPA |
| 发育轨迹 | CytoTRACE2 |
| 进度条 | pbapply |
| 字符串处理 | stringr |
| 聚类评估 | cluster |

### Python 依赖

```
torch (CPU 或 CUDA), scanpy, anndata, numpy, pandas, scipy, h5py
```

---

## 版本历史 / 重要修改记录

### v0.1.0 (当前)

- **核心功能**：DNN_train, DNN_predict (多任务分类+嵌入)
- **反卷积**：DNN_deconv_train, DNN_deconv_predict (pseudo-bulk 训练, log-CPM 输入, softmax 输出)
- **可视化**：12 个 Plot 函数 (Alluvia, RankScatter, Importance, Roe, MAP, Perturbation, Percent, Correlation, PropCorrelation, Scatter, Annotation, Sankey)
- **下游分析**：RankPerturbation (5 种距离), RankPercent (置换检验), RunPathwayAnalysis (GSEA/SCPA), ComputeModuleScore (AUCell/Seurat/UCell), EvaluateEmbedding, RunCytoTRACE2
- **Demo**：9 个完整示例脚本
