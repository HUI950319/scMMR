# ============================================================================
# Demo: RunCytoTRACE2 & PlotCytoTRACE2
# ============================================================================

library(scMMR)
library(Seurat)
library(ggplot2)

# ── 1. 加载数据 ──────────────────────────────────────────────────────────────
# 替换为你自己的 Seurat 对象路径
seu <- readRDS("/home/oyh/project/para/GRN/scenic_pipe/input/your_seurat.rds")

cat("Cells:", ncol(seu), "\n")
cat("Genes:", nrow(seu), "\n")
cat("Cell types:", paste(unique(seu$cell_type), collapse = ", "), "\n")


# ── 2. 运行 CytoTRACE2 ──────────────────────────────────────────────────────
seu <- RunCytoTRACE2(seu, species = "human")

# 查看添加的 metadata 列
head(seu@meta.data[, grep("CytoTRACE", colnames(seu@meta.data))])


# ── 3. 可视化：PlotCytoTRACE2 (violin + box) ────────────────────────────────
p1 <- PlotCytoTRACE2(seu, group.by = "cell_type")
print(p1)

# 自定义颜色
p2 <- PlotCytoTRACE2(seu, group.by = "cell_type",
                      colors = c("#E64B35", "#4DBBD5", "#00A087",
                                 "#3C5488", "#F39B7F", "#8491B4"))
print(p2)


# ── 4. Seurat 原生可视化 ────────────────────────────────────────────────────
# UMAP 上显示 potency score
p3 <- FeaturePlot(seu, features = "CytoTRACE2_Score",
                  cols = c("lightgrey", "red")) +
  ggtitle("CytoTRACE2 Potency Score")
print(p3)

# Potency 类别分布
p4 <- DimPlot(seu, group.by = "CytoTRACE2_Potency") +
  ggtitle("CytoTRACE2 Potency Category")
print(p4)

# 分组 violin
p5 <- VlnPlot(seu, features = "CytoTRACE2_Score",
              group.by = "cell_type", pt.size = 0) +
  geom_boxplot(width = 0.1, outlier.shape = NA) +
  NoLegend()
print(p5)


# ── 5. 组间比较（如 PH vs SH） ──────────────────────────────────────────────
# 按条件拆分后比较 potency
if ("condition" %in% colnames(seu@meta.data)) {
  p6 <- VlnPlot(seu, features = "CytoTRACE2_Score",
                group.by = "cell_type", split.by = "condition",
                pt.size = 0) +
    ggtitle("CytoTRACE2 Score: PH vs SH")
  print(p6)
}


# ── 6. 用矩阵直接运行（非 Seurat） ─────────────────────────────────────────
mat <- GetAssayData(seu, layer = "counts")
result_df <- RunCytoTRACE2(mat, species = "human")
head(result_df)
