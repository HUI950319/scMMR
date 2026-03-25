# Package index

## Python Environment

Set up and manage the Python environment required by scMMR.

- [`use_scMMR_python()`](https://hui950319.github.io/scMMR/reference/use_scMMR_python.md)
  : Set the Python environment for scMMR
- [`install_scMMR_python()`](https://hui950319.github.io/scMMR/reference/install_scMMR_python.md)
  : Install Python dependencies for scMMR

## Multi-Task DNN Training & Prediction

Train and predict with PyTorch multi-task deep neural networks for joint
cell type classification and embedding regression.

- [`DNN_train()`](https://hui950319.github.io/scMMR/reference/DNN_train.md)
  : Train a Multi-Task DNN Model
- [`DNN_predict()`](https://hui950319.github.io/scMMR/reference/DNN_predict.md)
  : Predict Cell Types Using a Trained Multi-Task DNN Model

## Deconvolution

Train deconvolution models from scRNA-seq references and estimate cell
type proportions in bulk RNA-seq.

- [`DNN_deconv_train()`](https://hui950319.github.io/scMMR/reference/DNN_deconv_train.md)
  : Train a Deconvolution DNN Model
- [`DNN_deconv_predict()`](https://hui950319.github.io/scMMR/reference/DNN_deconv_predict.md)
  : Predict Cell Type Proportions from Bulk RNA-seq

## Gene Set Scoring

Compute per-cell gene set activity scores using multiple methods.

- [`ComputeModuleScore()`](https://hui950319.github.io/scMMR/reference/ComputeModuleScore.md)
  : Calculate gene module scores
- [`read_gmt()`](https://hui950319.github.io/scMMR/reference/read_gmt.md)
  : Read a GMT (Gene Matrix Transposed) file
- [`parse_gene_sets()`](https://hui950319.github.io/scMMR/reference/parse_gene_sets.md)
  : Parse gene sets from various input formats

## Pathway Analysis

Run GSEA or SCPA pathway enrichment analysis across cell types.

- [`RunPathwayAnalysis()`](https://hui950319.github.io/scMMR/reference/RunPathwayAnalysis.md)
  : Run Pathway Analysis across Cell Types
- [`RunGsea()`](https://hui950319.github.io/scMMR/reference/RunGsea.md)
  : Run Gene Set Enrichment Analysis on DE Results
- [`RunGseaEnrich()`](https://hui950319.github.io/scMMR/reference/RunGseaEnrich.md)
  : Run compareCluster-based Enricher and GSEA Analysis
- [`CrossEnrichOverlap()`](https://hui950319.github.io/scMMR/reference/CrossEnrichOverlap.md)
  : Cross-comparison of Core Enrichment Genes between Two GSEA Results

## Differential Expression

Differential expression analysis between conditions.

- [`RunDE()`](https://hui950319.github.io/scMMR/reference/RunDE.md) :
  Run Differential Expression Test
- [`RunCorrelation()`](https://hui950319.github.io/scMMR/reference/RunCorrelation.md)
  : Correlation Ranking of One Variable Against Many

## Ranking & Perturbation

Rank genes, pathways, and regulons by perturbation strength or
differential abundance.

- [`RankPerturbation()`](https://hui950319.github.io/scMMR/reference/RankPerturbation.md)
  : Rank Cell Types by Transcriptional Perturbation
- [`RankPercent()`](https://hui950319.github.io/scMMR/reference/RankPercent.md)
  : Differential Abundance Testing via KNN Neighborhoods
- [`FindRegulonDrivers()`](https://hui950319.github.io/scMMR/reference/FindRegulonDrivers.md)
  : Find Upstream Regulon Drivers of Pathways

## Trajectory & Potency

Pseudotime trajectory and cellular potency inference.

- [`RunCytoTRACE2()`](https://hui950319.github.io/scMMR/reference/RunCytoTRACE2.md)
  : Run CytoTRACE 2 Cellular Potency Prediction
- [`RunTraceGene()`](https://hui950319.github.io/scMMR/reference/RunTraceGene.md)
  : Find Pseudotime-Associated Genes Along Trajectories
- [`RunTraceGSEA()`](https://hui950319.github.io/scMMR/reference/RunTraceGSEA.md)
  : Run Trajectory GSEA Along Pseudotime

## Embedding Evaluation

Evaluate DNN embedding quality vs PCA.

- [`EvaluateEmbedding()`](https://hui950319.github.io/scMMR/reference/EvaluateEmbedding.md)
  : Evaluate DNN Embedding Quality

## Quality Control

Doublet detection and ambient RNA estimation.

- [`ComputeDoublets()`](https://hui950319.github.io/scMMR/reference/ComputeDoublets.md)
  : Mark doublets via DoubletFinder
- [`ComputeAmbientRNA()`](https://hui950319.github.io/scMMR/reference/ComputeAmbientRNA.md)
  : Remove ambient RNA contamination via decontX

## Standard Pipeline

End-to-end analysis pipeline.

- [`StandardPipeline()`](https://hui950319.github.io/scMMR/reference/StandardPipeline.md)
  : Standard single-cell analysis pipeline

## Visualization

Publication-ready plots for single-cell analysis results.

- [`PlotRoe()`](https://hui950319.github.io/scMMR/reference/PlotRoe.md)
  : O/E Ratio Heatmap for Cell Type Composition
- [`PlotSankey()`](https://hui950319.github.io/scMMR/reference/PlotSankey.md)
  : Sankey Plot for Cell Population Composition
- [`PlotAlluvia()`](https://hui950319.github.io/scMMR/reference/PlotAlluvia.md)
  : Alluvial Plot for Cell Population Composition
- [`PlotAnnotation()`](https://hui950319.github.io/scMMR/reference/PlotAnnotation.md)
  : Annotation Heatmap for Single-Cell Metadata
- [`PlotScatter()`](https://hui950319.github.io/scMMR/reference/PlotScatter.md)
  : Scatter Correlation Plot
- [`PlotCorrelation()`](https://hui950319.github.io/scMMR/reference/PlotCorrelation.md)
  : Volcano-Style Correlation Plot
- [`PlotPropCorrelation()`](https://hui950319.github.io/scMMR/reference/PlotPropCorrelation.md)
  : Cell-Type Proportion vs Gene Expression Correlation
- [`PlotRankScatter()`](https://hui950319.github.io/scMMR/reference/PlotRankScatter.md)
  : Rank Scatter Plot for Feature Scoring
- [`PlotImportance()`](https://hui950319.github.io/scMMR/reference/PlotImportance.md)
  : Importance Plot (Lollipop / Bar)
- [`PlotMAP()`](https://hui950319.github.io/scMMR/reference/PlotMAP.md)
  : UMAP Projection Plot
- [`PlotPercent()`](https://hui950319.github.io/scMMR/reference/PlotPercent.md)
  : Differential Abundance Beeswarm Plot
- [`PlotPerturbation()`](https://hui950319.github.io/scMMR/reference/PlotPerturbation.md)
  : Perturbation Ranking Plot
- [`PlotDE()`](https://hui950319.github.io/scMMR/reference/PlotDE.md) :
  Plot Differential Expression Results
- [`PlotGsea()`](https://hui950319.github.io/scMMR/reference/PlotGsea.md)
  : Plot GSEA Results
- [`PlotPathwayBubble()`](https://hui950319.github.io/scMMR/reference/PlotPathwayBubble.md)
  : Pathway Bubble Plot
- [`PlotDynamicFeatures()`](https://hui950319.github.io/scMMR/reference/PlotDynamicFeatures.md)
  : Plot Dynamic Features Along Pseudotime
- [`PlotEmbeddingEval()`](https://hui950319.github.io/scMMR/reference/PlotEmbeddingEval.md)
  : Plot Embedding Evaluation Results
- [`PlotCytoTRACE2()`](https://hui950319.github.io/scMMR/reference/PlotCytoTRACE2.md)
  : CytoTRACE 2 Potency Plot

## Color Palettes

Built-in color palettes for consistent styling.

- [`palette_colors()`](https://hui950319.github.io/scMMR/reference/palette_colors.md)
  : Generate Named Color Vector from Palette
- [`show_palettes()`](https://hui950319.github.io/scMMR/reference/show_palettes.md)
  : Display Available Color Palettes

## Data

Built-in datasets.

- [`ribo.genes`](https://hui950319.github.io/scMMR/reference/ribo.genes.md)
  : Ribosomal gene symbols (human)
