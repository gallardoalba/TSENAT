# Changelog

## TSENAT 0.99.0

- Initial Bioconductor release. TSENAT provides scale-dependent analysis
  of transcript isoform diversity using Tsallis entropy, enabling
  detection of splicing-driven regulatory changes orthogonal to
  count-based differential expression methods.

- Core analysis via
  [`calculate_diversity()`](https://gallardoalba.github.io/TSENAT/reference/calculate_diversity.md)
  /
  [`calculate_diversity()`](https://gallardoalba.github.io/TSENAT/reference/calculate_diversity.md),
  [`calculate_difference()`](https://gallardoalba.github.io/TSENAT/reference/calculate_difference.md)
  /
  [`calculate_difference()`](https://gallardoalba.github.io/TSENAT/reference/calculate_difference.md),
  [`calculate_divergence()`](https://gallardoalba.github.io/TSENAT/reference/calculate_divergence.md)
  /
  [`calculate_divergence()`](https://gallardoalba.github.io/TSENAT/reference/calculate_divergence.md),
  and `calculate_lm_interaction()` / `calculate_lm_interaction_s4()`
  supporting diverse statistical tests (Wilcoxon, Friedman, permutation,
  jackknife bootstrap).

- Parameter `q` tunes Tsallis entropy sensitivity: q \< 1 emphasizes
  rare isoforms, q ≈ 1 recovers Shannon entropy, q \> 1 emphasizes
  dominant isoforms.

- Unified S4 class `TSENATAnalysis` integrating SummarizedExperiment,
  configuration, results, and cached plots with full subsetting and
  accessor support.

- Data input: `read_salmon_samples()` for direct Salmon/Sailfish
  quantification import with GFF3 annotation parsing.

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
  (SplicingFactory validation, advanced workflows), 46 exported
  functions with complete roxygen2 documentation.
