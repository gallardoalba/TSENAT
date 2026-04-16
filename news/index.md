# Changelog

## TSENAT 0.99.0

- Initial Bioconductor release. TSENAT provides scale-dependent analysis
  of transcript isoform diversity using Tsallis entropy, enabling
  detection of splicing-driven regulatory changes orthogonal to
  count-based differential expression methods.

- Core analysis via
  [`calculate_diversity()`](https://gallardoalba.github.io/TSENAT/reference/calculate_diversity.md),
  [`calculate_divergence()`](https://gallardoalba.github.io/TSENAT/reference/calculate_divergence.md),
  [`calculate_srh()`](https://gallardoalba.github.io/TSENAT/reference/calculate_srh.md)
  for Q×Condition interaction testing,
  [`calculate_lm()`](https://gallardoalba.github.io/TSENAT/reference/calculate_lm.md)
  for linear models/GAM interaction analysis, and
  [`calculate_concordance()`](https://gallardoalba.github.io/TSENAT/reference/calculate_concordance.md)
  for method comparison. Supports Scheirer-Ray-Hare rank-based tests,
  permutation tests, and jackknife bootstrap.

- Parameter `q` tunes Tsallis entropy sensitivity: q \< 1 emphasizes
  rare isoforms, q ≈ 1 recovers Shannon entropy, q \> 1 emphasizes
  dominant isoforms.

- Unified S4 class `TSENATAnalysis` integrating SummarizedExperiment,
  configuration, results, and cached plots with full subsetting and
  accessor support.

- Data input:
  [`build_analysis()`](https://gallardoalba.github.io/TSENAT/reference/build_analysis.md)
  for RNA-seq count matrices or Salmon quantification (via internal
  utilities); supports GFF3 annotation files for transcript-to-gene
  mapping via `tx2gene` parameter.

- Quality control:
  [`filter_analysis()`](https://gallardoalba.github.io/TSENAT/reference/filter_analysis.md)
  for multi-criteria filtering; `jackknife_entropy_outliers_s4()` for
  outlier detection.

- Visualization suite (8 plot types): diversity q-curves, volcano/MA
  plots, violin/density plots, transcript composition heatmaps,
  divergence distance matrices, interaction surfaces, concordance plots,
  multi-gene q-spectrum plots. All ggplot2-based with publication-ready
  styling.

- High-level workflow orchestration via
  [`TSENAT()`](https://gallardoalba.github.io/TSENAT/reference/TSENAT.md)
  and flexible configuration management with
  [`TSENAT_config()`](https://gallardoalba.github.io/TSENAT/reference/TSENAT_config.md).

- Performance: Rcpp/RcppArmadillo for entropy calculations, BiocParallel
  support, lazy-loading of visualization libraries (~30% faster non-plot
  workflows).

- Comprehensive documentation: main vignette + 2 appendices
  (SplicingFactory validation, advanced workflows), 25 exported
  functions with complete roxygen2 documentation.
