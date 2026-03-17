# ============================================================================
# scMMR Demo: RunPathwayAnalysis & PlotPathwayBubble
# ============================================================================

library(scMMR)
library(Seurat)
library(qs)

# ── 1. Load toy test data ──────────────────────────────────────────────────
cat("== 1. Loading toy test data ==\n")

test_path <- system.file("extdata", "toy_test.qs", package = "scMMR")
seu <- qread(test_path)
cat("Loaded:", ncol(seu), "cells,", nrow(seu), "genes\n")

DefaultAssay(seu) <- "RNA"
seu <- NormalizeData(seu, verbose = FALSE)

cat("Groups:", paste(sort(unique(seu$group)), collapse = ", "), "\n")
cat("Cell types:", length(unique(seu$cell_type)), "\n\n")
print(table(seu$cell_type, seu$group))

ct_tab <- table(seu$cell_type, seu$group)
valid_cts <- rownames(ct_tab)[ct_tab[, "PH"] >= 5 & ct_tab[, "PT"] >= 5]
cat("\nValid cell types:", paste(valid_cts, collapse = ", "), "\n")


# ── 2. Prepare gene sets ──────────────────────────────────────────────────
cat("\n== 2. Preparing gene sets ==\n")

gmt_path <- system.file("extdata", "gmt",  "proliferation.combined.gmt", # "h.all.v2022.1.Hs.symbols.gmt",
                         package = "scMMR")
gmt_sets <- parse_gene_sets(gmt_path)
cat("Total gene sets:", length(gmt_sets), "\n")

overlap <- vapply(gmt_sets, function(g) length(intersect(g, rownames(seu))), integer(1))
cat("Sets with >= 5 genes in data:", sum(overlap >= 5), "\n")


# ── 3. GSEA method ────────────────────────────────────────────────────────
cat("\n== 3. GSEA method ==\n")

test_cts <- head(valid_cts, 3)
test_cts <- "Parathyroid cells"

cat("Cell types:", paste(test_cts, collapse = ", "), "\n\n")

gsea_res <- RunPathwayAnalysis(
  seu,
  gene.sets    = gmt_path,
  method       = "GSEA",
  celltype.col = "cell_type",
  group.col    = "group",
  ident.1      = "PH",
  ident.2      = "PT",
  celltypes    = test_cts,
  pvalue.cutoff = 1,
  logfc.threshold = 0,
  minGSSize    = 5,
  clean.names  = TRUE,
  verbose      = TRUE
)

cat("\nGSEA:", nrow(gsea_res), "rows x", ncol(gsea_res), "cols\n")
if (nrow(gsea_res) > 0) {
  print(head(gsea_res[, c("Celltype", "Pathway", "Score", "AdjPValue", "Sign")], 5))
}


# ── 4. SCPA method ────────────────────────────────────────────────────────
cat("\n== 4. SCPA method ==\n")

scpa_res <- RunPathwayAnalysis(
  seu,
  gene.sets    = gmt_path,
  method       = "SCPA",
  celltype.col = "cell_type",
  group.col    = "group",
  ident.1      = "SH",
  ident.2      = "PT",
  celltypes    = test_cts,
  clean.names  = TRUE,
  verbose      = TRUE)


cat("\nSCPA:", nrow(scpa_res), "rows x", ncol(scpa_res), "cols\n")
if (nrow(scpa_res) > 0) {
  print(head(scpa_res[, c("Celltype", "Pathway", "Score", "AdjPValue", "Sign")], 5))
}


# ── 5. PlotPathwayBubble ──────────────────────────────────────────────────
cat("\n== 5. PlotPathwayBubble ==\n")

PlotPathwayBubble(gsea_res, top.n = 20, pvalue.cutoff = 1, facet.by.sign = T)
PlotPathwayBubble(gsea_res, top.n = 20, pvalue.cutoff = 0.1, size.by = "score")





PlotPathwayBubble(scpa_res, top.n = 5, pvalue.cutoff = 1, facet.by.sign = T)
PlotPathwayBubble(scpa_res, top.n = 5, pvalue.cutoff = 1, size.by = "qvalue")



# ── 6. Gene sets as named list ────────────────────────────────────────────
cat("\n== 6. Named list input ==\n")

gene_list <- list(
  GeneSet_A = head(rownames(seu), 50),
  GeneSet_B = rownames(seu)[51:100],
  GeneSet_C = rownames(seu)[101:150]
)

res_list <- RunPathwayAnalysis(
  seu,
  gene.sets    = gene_list,
  method       = "GSEA",
  celltype.col = "cell_type",
  group.col    = "group",
  ident.1      = "PH",
  ident.2      = "PT",
  celltypes    = test_cts[1],
  pvalue.cutoff = 1,
  logfc.threshold = 0,
  clean.names  = FALSE,
  verbose      = TRUE
)
cat("Named list input:", nrow(res_list), "results\n")


# ── Done ──────────────────────────────────────────────────────────────────
cat("\n== All tests completed! ==\n")






gmt_dir <- system.file("extdata", "gmt", package = "scMMR")
tf_sets <- parse_gene_sets(file.path(gmt_dir, "collectri.human.gmt"))

# PROGENy pathways
pw_sets <- parse_gene_sets(file.path(gmt_dir, "progeny.human.top100.gmt"))

# 可直接用于 RunPathwayAnalysis
res <- RunPathwayAnalysis(seu, 
                          gene.sets = file.path(gmt_dir, "progeny.human.top100.gmt"),
                          method       = "GSEA",
                          celltype.col = "cell_type",
                          group.col    = "group",
                          ident.1      = "PH",
                          ident.2      = "PT",
                          celltypes    = "Parathyroid cells",
                          pvalue.cutoff = 1,
                          logfc.threshold = 0)

res <- RunPathwayAnalysis(seu, 
                          gene.sets = file.path(gmt_dir, "progeny.human.top100.gmt"),
                          method       = "GSEA",
                          celltype.col = "cell_type",
                          group.col    = "group",
                          ident.1      = "SH",
                          ident.2      = "PT",
                          celltypes    = "Parathyroid cells",
                          pvalue.cutoff = 1,
                          logfc.threshold = 0)

PlotPathwayBubble(res,pvalue.cutoff = 1)



pw_sets <- parse_gene_sets(file.path(gmt_dir, "proliferation.combined.gmt"))



















