# Changelog

## scMMR 0.2.0

### Breaking Changes

- Renamed `PlotGroupPreference()` to
  [`PlotRoe()`](https://hui950319.github.io/scMMR/reference/PlotRoe.md)
  with simplified parameters (`by`, `fill`).

### New Features

- [`PlotRoe()`](https://hui950319.github.io/scMMR/reference/PlotRoe.md)
  gains `display` parameter: `"value"` (numeric O/E), `"symbol"` (grade:
  +++/++/+/+−/−), or `"both"`.
- [`PlotRoe()`](https://hui950319.github.io/scMMR/reference/PlotRoe.md)
  gains `show.pvalue` and `p.method` parameters for per-cell
  significance testing via chi-squared residuals or Fisher’s exact test.
- [`PlotRoe()`](https://hui950319.github.io/scMMR/reference/PlotRoe.md)
  gains `return.data` to export the O/E matrix and p-value matrix.
- [`PlotSankey()`](https://hui950319.github.io/scMMR/reference/PlotSankey.md)
  for Sankey/alluvial flow diagrams with aligned parameters to
  [`PlotAlluvia()`](https://hui950319.github.io/scMMR/reference/PlotAlluvia.md).
- Added pkgdown website with 8 tutorial articles covering all major
  workflows.
- Added GitHub Actions CI for automatic website deployment.

### Internal

- Refactored monolithic `plot.R` (~1243 lines) into 12 individual
  `Plot*.R` files + `plot_utils.R`.

## scMMR 0.1.0

- Initial release.
- **Core DNN**:
  [`DNN_train()`](https://hui950319.github.io/scMMR/reference/DNN_train.md),
  [`DNN_predict()`](https://hui950319.github.io/scMMR/reference/DNN_predict.md)
  for multi-task cell type classification + embedding regression.
- **Deconvolution**:
  [`DNN_deconv_train()`](https://hui950319.github.io/scMMR/reference/DNN_deconv_train.md),
  [`DNN_deconv_predict()`](https://hui950319.github.io/scMMR/reference/DNN_deconv_predict.md)
  for bulk RNA-seq deconvolution.
- **Gene Set Scoring**:
  [`ComputeModuleScore()`](https://hui950319.github.io/scMMR/reference/ComputeModuleScore.md)
  with AUCell, Seurat, and UCell methods.
- **Pathway Analysis**:
  [`RunPathwayAnalysis()`](https://hui950319.github.io/scMMR/reference/RunPathwayAnalysis.md)
  (GSEA/SCPA),
  [`PlotPathwayBubble()`](https://hui950319.github.io/scMMR/reference/PlotPathwayBubble.md).
- **Perturbation Ranking**:
  [`RankPerturbation()`](https://hui950319.github.io/scMMR/reference/RankPerturbation.md)
  (5 distance metrics),
  [`RankPercent()`](https://hui950319.github.io/scMMR/reference/RankPercent.md)
  (neighborhood-based DA).
- **Trajectory**:
  [`RunCytoTRACE2()`](https://hui950319.github.io/scMMR/reference/RunCytoTRACE2.md),
  [`PlotCytoTRACE2()`](https://hui950319.github.io/scMMR/reference/PlotCytoTRACE2.md).
- **Embedding Evaluation**:
  [`EvaluateEmbedding()`](https://hui950319.github.io/scMMR/reference/EvaluateEmbedding.md),
  [`PlotEmbeddingEval()`](https://hui950319.github.io/scMMR/reference/PlotEmbeddingEval.md).
- **Visualization**: 12 Plot functions (Alluvia, RankScatter,
  Importance, MAep, Perturbation, Percent, Correlation, PropCorrelation,
  Scatter, Annotation, DE, Gsea).
- **Quality Control**:
  [`ComputeDoublets()`](https://hui950319.github.io/scMMR/reference/ComputeDoublets.md),
  [`ComputeAmbientRNA()`](https://hui950319.github.io/scMMR/reference/ComputeAmbientRNA.md).
- 14 bundled GMT databases (Hallmark, KEGG, Reactome, GO, CollecTRI,
  PROGENy, etc.).
- 9 demo scripts in `inst/demo/`.
