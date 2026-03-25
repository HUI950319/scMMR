# scMMR 0.2.0

## Breaking Changes

* Renamed `PlotGroupPreference()` to `PlotRoe()` with simplified parameters (`by`, `fill`).

## New Features

* `PlotRoe()` gains `display` parameter: `"value"` (numeric O/E), `"symbol"` (grade: +++/++/+/+−/−), or `"both"`.
* `PlotRoe()` gains `show.pvalue` and `p.method` parameters for per-cell significance testing via chi-squared residuals or Fisher's exact test.
* `PlotRoe()` gains `return.data` to export the O/E matrix and p-value matrix.
* `PlotSankey()` for Sankey/alluvial flow diagrams with aligned parameters to `PlotAlluvia()`.
* Added pkgdown website with 8 tutorial articles covering all major workflows.
* Added GitHub Actions CI for automatic website deployment.

## Internal

* Refactored monolithic `plot.R` (~1243 lines) into 12 individual `Plot*.R` files + `plot_utils.R`.

# scMMR 0.1.0

* Initial release.
* **Core DNN**: `DNN_train()`, `DNN_predict()` for multi-task cell type classification + embedding regression.
* **Deconvolution**: `DNN_deconv_train()`, `DNN_deconv_predict()` for bulk RNA-seq deconvolution.
* **Gene Set Scoring**: `ComputeModuleScore()` with AUCell, Seurat, and UCell methods.
* **Pathway Analysis**: `RunPathwayAnalysis()` (GSEA/SCPA), `PlotPathwayBubble()`.
* **Perturbation Ranking**: `RankPerturbation()` (5 distance metrics), `RankPercent()` (neighborhood-based DA).
* **Trajectory**: `RunCytoTRACE2()`, `PlotCytoTRACE2()`.
* **Embedding Evaluation**: `EvaluateEmbedding()`, `PlotEmbeddingEval()`.
* **Visualization**: 12 Plot functions (Alluvia, RankScatter, Importance, MAep, Perturbation, Percent, Correlation, PropCorrelation, Scatter, Annotation, DE, Gsea).
* **Quality Control**: `ComputeDoublets()`, `ComputeAmbientRNA()`.
* 14 bundled GMT databases (Hallmark, KEGG, Reactome, GO, CollecTRI, PROGENy, etc.).
* 9 demo scripts in `inst/demo/`.
