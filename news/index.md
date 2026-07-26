# Changelog

## TSENAT 0.99.33

- **Statistical implementation audit (July 2026)**: Fixed 22 bugs from
  systematic review of `sait_*`, `diversity_*`, `divergence_*`, and
  `entropy_*` modules.

  - **Critical fixes (4)**:
    - `.hochberg_stepup()` now uses correct Hochberg step-up
      (`rev(cummin(rev(...)))`) instead of Holm step-down (`cummax`).
      Matches `p.adjust(..., "hochberg")` exactly.
    - GEE bias-correction thresholds unified from mixed `<20`/`<30` to
      consistent `<30`.
    - Kauermann-Carroll HC1 multiplier now actually applied —
      `vcov(fit_alt)` and `coef_value` passed to `.kc_bias_correct()`.
    - Bootstrap CIs no longer invalidated by cross-gene normalization —
      remain on raw divergence scale.
  - **High-impact fixes (7)**:
    - `(method, regularization)` validation prevents `match.arg` errors
      with incompatible combos.
    - AR(1) design-effect replaced asymptotic formula with correct
      finite-m form; dead duplicate removed.
    - Westfall-Young permutation preserves paired structure (permutes
      within subjects).
    - ARIMA differencing uses `group[-1]` to preserve both condition
      levels.
    - `log_base` threaded through diversity computation chain (was
      silently nats).
    - `method="bca"` warns and reports `"percentile"` instead of
      silently substituting.
    - `log_odds_ratio` uses data-driven column-maximum normalization.
  - **Other fixes (11)**: Removed ~280 lines dead stationarity code;
    standardized `min_obs`; exposed `corstr="auto"`; FPCA MANOVA returns
    `NA` on failure; AR(1) ρ pooled within-cluster; q=0 divergence is
    support-difference; paired bootstrap warns on fallback; concordance
    handles empty results; plus 4 low-severity robustness fixes.

## TSENAT 0.99.31

- **Aligned Rank Transform (ART)**:
  [`calculate_rank_transform()`](https://gallardoalba.github.io/TSENAT/reference/calculate_rank_transform.md)
  now defaults to the Aligned Rank Transform via the ARTool package (Kay
  et al. 2021) for proper non-parametric interaction testing. ART strips
  main effects before ranking (“alignment”), preserving interaction
  structure — a known limitation of classical rank-transform methods.
  The Conover-Iman Rank Transform remains available via `method='rt'`.
  ARTool added to Imports.

- **Bug fixes (July 2026 audit)**: Comprehensive fixes from systematic
  code audit:

  - **Entropy & divergence core**: Fixed `.entropy_core()` q=0 species
    richness inflated by zero-count isoforms (B1); fixed
    `.entropy_max()` q=0 theoretical maximum off-by-one (B2); threaded
    `log_base` through `.normalize_log_odds_ratio()` pipeline (B3);
    restricted `log_base` to KL limit (q≈1) in Tsallis divergence (B10);
    eliminated pseudocount + min_prob double-correction (B11); corrected
    max divergence normalization for q≠1 using q-dependent formula
    (B12).

  - **C++ resampling kernels**: Fixed `entropy_cpp()` q=0 returning
    `log(n)` instead of `n-1`; removed `>1e-15` threshold so zeros
    contribute zero entropy; added `RNGseed` for reproducible parallel
    execution.

  - **Parallel reproducibility (R-layer)**: Fixed `.bplapply()` not
    propagating [`set.seed()`](https://rdrr.io/r/base/Random.html) to
    BiocParallel workers, causing non-deterministic results across
    parallel runs (Westfall-Young, bootstrap, etc.); now derives
    `RNGseed` from current RNG state so all parallel operations are
    reproducible with
    [`set.seed()`](https://rdrr.io/r/base/Random.html).

  - **Bootstrap & CI infrastructure**: Fixed BCa acceleration computed
    from bootstrap distribution instead of true jackknife (audit
    [\#3](https://github.com/gallardoalba/TSENAT/issues/3)); fixed BCa
    z0 using wrong point estimate with effective_length mismatch (audit
    [\#4](https://github.com/gallardoalba/TSENAT/issues/4)); added
    defensive copies in `bootstrap_compute_cpp_wrapper()` and
    `block_bootstrap_compute_cpp_wrapper()` (B6, I9); fixed BCa
    degenerate distribution handling and strict comparison with 0.5
    padding (audit
    [\#12](https://github.com/gallardoalba/TSENAT/issues/12));
    implemented replicate-level bootstrap with proper C++ path (audit
    [\#11](https://github.com/gallardoalba/TSENAT/issues/11)).

  - **Rank-based testing**: Renamed all Scheirer-Ray-Hare labels to
    Conover-Iman Rank Transform (B4, I3); fixed η² computed on ranks
    instead of original entropy scale via separate
    [`lm()`](https://rdrr.io/r/stats/lm.html) (B7, I5).

  - **Shrinkage & edge cases**: Eliminated silent zero-filling of NAs in
    `.apply_shrinkage()` by initializing `means_vector` to `NA_real_`
    (B5, I8).

  - **Vectorized divergence**: Fixed `.tsallis_divergence_vector()` q=0
    returning scalar 0 instead of per-pair divergences (B13); fixed
    `.divergence_normalize_log_odds_ratio()` hardcoded maximum
    divergence ignoring q-dependence (B17).

  - **GAM interaction models**: Identified GAM LRT statistically invalid
    due to REML + non-nested models + different k (audit
    [\#1](https://github.com/gallardoalba/TSENAT/issues/1)); identified
    missing main-effect smooth s(q) in interaction model making
    comparison non-nested (audit
    [\#2](https://github.com/gallardoalba/TSENAT/issues/2)); removed
    undeclared p-value bias correction multiplier (B14, audit
    [\#8](https://github.com/gallardoalba/TSENAT/issues/8)); scoped
    [`suppressWarnings()`](https://rdrr.io/r/base/warning.html) to only
    suppress rank-deficient warnings instead of all warnings globally
    (B15).

  - **Heteroscedasticity & weights**: Fixed
    `.handle_arima_and_weights()` silently dropping heteroscedasticity
    weights by directly assigning to correlation structure instead of
    calling [`update()`](https://rdrr.io/r/stats/update.html) (B16).

  - **Normalization bias**: Fixed `.normalize_range_matrix()` distorting
    CI bounds due to per-column normalization applied globally; changed
    to matrix-wide \[0,1\] normalization preserving rank order (B18).

  - **Convergence diagnostics**: Fixed `.mest_irls_location()`
    convergence check skipping first iteration (B8); removed allocated
    matrix immediately overwritten in `.estimate_storey_pi0()` (B9).

  - **FPCA gene-level testing**: Deprecated
    `min(BH-adjusted PC p-values)` as gene-level aggregate; replaced
    with sum-of-χ² pooling following Crainiceanu et al. (2009) guidance
    (audit [\#9](https://github.com/gallardoalba/TSENAT/issues/9)).

  - **JIS bootstrap p-values**: Fixed inverted JIS bootstrap p-values
    where `p_boot = 1 - mean(null ≥ obs)` never rejected; corrected to
    `mean(null ≥ obs)` (audit
    [\#7](https://github.com/gallardoalba/TSENAT/issues/7)).

- **Code quality**: Split monolithic `plots_helpers.R` (4,385 lines)
  into 5 thematic files; split `s4_functions.R` (2,561 lines, 11
  exports) into 11 per-function files; split bootstrap test file (7,752
  lines) into 6 focused test files.

- **Testing**: Added 85 golden reference tests for entropy/divergence
  formulas; added 39 ART-specific tests; test suite now exceeds 4,800
  test blocks.

- **Documentation**: Updated all vignettes and man pages to reflect ART
  as default method; corrected Hochberg vs. Benjamini-Hochberg labeling;
  fixed heteroscedasticity interpretation in main vignette; added ART
  citations (S271–S273) to TSENAT.bib.

## TSENAT 0.99.0

- Initial Bioconductor release. TSENAT provides scale-dependent analysis
  of transcript isoform diversity using Tsallis entropy, enabling
  detection of splicing-driven regulatory changes orthogonal to
  count-based differential expression methods.

- Core analysis via
  [`calculate_diversity()`](https://gallardoalba.github.io/TSENAT/reference/calculate_diversity.md),
  [`calculate_divergence()`](https://gallardoalba.github.io/TSENAT/reference/calculate_divergence.md),
  [`calculate_rank_transform()`](https://gallardoalba.github.io/TSENAT/reference/calculate_rank_transform.md)
  for Q×Condition interaction testing,
  [`calculate_sait()`](https://gallardoalba.github.io/TSENAT/reference/calculate_sait.md)
  for scale-adaptive interaction testing (GAM, LMM, GEE, FPCA) with
  AR(1) support for repeated measures, and
  [`calculate_concordance()`](https://gallardoalba.github.io/TSENAT/reference/calculate_concordance.md)
  for method comparison. Supports Conover-Iman Rank Transform tests,
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
