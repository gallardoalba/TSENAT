# Changelog

## TSENAT 0.99.35

- **Statistical inference hardening** (August 2026). Main changes:

  - **Functional reframing of q**: q is a deterministic functional
    argument of the Tsallis statistic, not a time index. ARIMA(1,1,0)
    differencing was removed from all confirmatory paths (GEE, GAM/GAMM,
    LMM, FPCA) — the interaction is now tested on the original H(q)
    curve (H0: beta(q) = 0 for all q); the legacy differencing helpers
    (`.compute_arima_differences`, `.apply_arima_differencing`,
    `.apply_arima_differencing_fpca`) and their tests were deleted, and
    the LMM fallback `slope_diff` is now extracted from fixed effects
    ([`nlme::fixef`](https://rdrr.io/pkg/nlme/man/fixed.effects.html))
    instead of per-subject coefficients..
  - **Correlation structure**: AR(1)-type within subject × condition
    modelled over **actual q distances**
    ([`nlme::corCAR1`](https://rdrr.io/pkg/nlme/man/corCAR1.html),
    Corr(e_i, e_j) = exp(-φ\|q_i - q_j\|)) as the primary structure in
    GAMM and LMM, so irregular q grids are handled correctly (the
    q-grid-index `rho^|Δgrid|` form, which is valid only on equally
    spaced grids, is kept as a documented fallback and for the legacy
    mgcv paths); GEE fits H(q) with a joint Wald test and a
    small-cluster F reference (df = n_clusters − p); the AR(1) design
    effect is descriptive-only.
  - **Salmon input integrity** (audit hardening): Salmon quantification
    files are now matched by **transcript ID** rather than positional
    order; transcript sets and sample identifiers must be unique and
    consistent across files (a mismatch is a hard error, not a warning),
    and negative or non-finite quantification values are rejected at
    input. Paired-design bootstrap also validates the 1-control +
    1-treatment-per-pair invariant before resampling. \*
    **TPM/effective-length contract**: diversity is computed from raw
    counts with effective-length correction; `tpm = TRUE` together with
    an `effective_length` (parameter or SummarizedExperiment metadata)
    is now a hard error, since TPM already incorporates effective-length
    normalization (double normalization rejected). TPM remains available
    for abundance-based filtering/QC.
  - **Pseudocount ordering and `'auto'` resolution**: effective-length
    normalization is now applied BEFORE the pseudocount, so
    regularization is constant on the effective-abundance scale
    (previously the count-space pseudocount was implicitly divided by
    transcript length, systematically boosting short isoforms);
    `pseudocount = 'auto'` is now estimated from the resolved raw-count
    matrix after input resolution and is rejected for TPM input. \*
    **Paired GAMM**:
    [`nlme::lme`](https://rdrr.io/pkg/nlme/man/lme.html) with
    `ns(q, df = 3) × condition` and a marginal F-test (mgcv gamm is
    singular on paired designs); no p-value underflow for strong signals
    (log-space recomputation) and pseudo-R² `effect_size`; `slope_diff`
    from population-level predictions; fit metadata records
    `model_used`/`fallback_level`/`correlation_structure`/`test_type`.
  - **Bootstrap resampling invariant**: the read-level bootstrap now
    resamples from exactly the point-estimate proportions
    `(x/l + c)/sum(x/l + c)` — the pseudocount is embedded on the
    effective-abundance scale BEFORE the depth rescale (previously it
    was added after the rescale, breaking the factorization for `c > 0`
    and shifting the resampling probabilities away from the assay
    estimate); the CI point estimate is computed by the estimator itself
    on the raw input, and q = 0 bootstrap replicates now carry the
    support distribution of the multinomial draws (`entropy_cpp` q = 0
    counts positive entries — zero-proportion bins no longer count as
    species). Locked by `test-bootstrap-invariant.R`.
  - **Provenance rename**: the primary paired GAMM metadata
    `fit_method`/`model_used` is now `lme_ns_car1` (continuous CAR(1)
    correlation over ACTUAL q distances), reserving
    `ar1_grid_within_subject_condition` for the grid-index fallback;
    previously the primary path was mislabelled `lme_ns_ar1`.
  - **Westfall–Young schemes**: `block_col`, `strata_col` and
    `permutation_scheme` with exchangeability validation (confounded
    blocks/strata rejected).
  - **Post-selection inference**: LASSO/ElasticNet selection is
    exploratory-only and never modifies the confirmatory model (testthat
    lock).
  - **Robust M-estimation**: sandwich variance with robust
    SE/p-values/95% CIs.
  - **Performance (ART)**: the Aligned Rank Transform no longer
    recomputes ANOVAs for all effects per gene — only the q × condition
    interaction row is extracted via
    [`ARTool::artlm()`](https://rdrr.io/pkg/ARTool/man/artlm.html) +
    `flat.anova()` (identical F and p-values, ~3× faster on the ANOVA
    step). Parallel per-gene workers now pin BLAS/OpenMP to a single
    thread to avoid core oversubscription (ART vignette step measured 42
    s → 3.3 s at nthreads = 3).
  - **Documented limitations**: unpaired ART slightly anti-conservative
    at small n (Conover-Iman `method='rt'` recommended for confirmatory
    unpaired inference); GEE mildly anti-conservative under strong
    heteroscedasticity + outliers; read-level bootstrap/divergence CIs
    under-cover at low depth; GEE QIC `corstr='auto'` validated under H0
    but pre-specification still recommended.
  - **Validation suite**: 22 Monte Carlo tests in `tests/testthat/`
    (type I, FWER/FDR, power, q-grid invariance, missing-q, filtering,
    effect sizes, bootstrap coverage, robustness), skipped on
    Bioconductor builds via `skip_on_bioc()` and runnable locally via
    `Rscript tests/testthat/run-validation.R`.
  - **Divergence reimplementation**: Tsallis divergence now compares
    per-condition **isoform** distributions (reads summed across samples
    within isoforms — the pooled condition-level composition) instead of
    sample-level aggregates; formula corrected to D_q(P\|\|Q) = (sum
    P_i^q Q_i^(1-q) - 1)/(q - 1) without
    [`abs()`](https://rdrr.io/r/base/MathFun.html) (KL only at q = 1, no
    ±0.01 band); q = 0 is evaluated on the RAW (pre-pseudocount) support
    — D_0 = 1 - sum\_{i: P_i\>0} Q_i — so the low-q end of the spectrum
    keeps its support-difference meaning under the default pseudocount
    instead of collapsing to 0 (review Option A); bootstrap CIs resample
    biological replicates (paired pairs as units; for paired designs the
    pairing is preserved when estimating uncertainty around the pooled
    estimate, not a subject-level divergence); `method='bca'` falls back
    to percentile with a warning and the effective method is recorded in
    metadata; default `norm = 'none'`; mathematical property test
    battery added (identity, non-negativity, KL limit, asymmetry,
    permutation/count-scale invariance, support, q=0, construct
    validity) plus coarse-graining, biological-replicate replication
    invariance, bootstrap consistency, directed support, and biological
    construct-validity scenarios (abundance-only, isoform switch, rare
    vs dominant remodeling).
  - **S4 layer hardening**: `paired`/`bootstrap`/ `norm`/`stringency`
    now default to `NULL` so `config` values are actually resolved
    (previously the literal defaults short-circuited the resolver);
    `calculate_assumptions(q=)` resolves q against available diversity
    keys with numeric tolerance (`q_1.000`/`q_0.5`/`q_1_00` conventions)
    and errors instead of silently analyzing another q;
    `TSENATAnalysis[i, j]` subsets divergence results by genes only
    (columns are q-values, not samples), subsets gene-level SAIT rows,
    and records `metadata$subset_applied`/`stale_results` so inferential
    results computed on the full dataset are explicitly flagged;
    [`calculate_jis()`](https://gallardoalba.github.io/TSENAT/reference/calculate_jis.md)
    restricts q to available diversity results (error if none); plot
    wrappers validate the object before touching `@config`; the
    constructor patches `colData` in place instead of rebuilding the
    SummarizedExperiment (preserves rowRanges/altExps); class validity
    no longer requires `sample_id`/ `gene_id`/`transcript_id` columns
    (module-level contracts instead); effect-size docs corrected to
    config-based resolution; 5 end-to-end S4 integration tests added
    (config `paired`/`bootstrap` reaching the core, exact q lookup, q=0
    preservation, subset q-column invariance, staleness flag).
  - **Plotting layer hardening**: the global divergence spectrum no
    longer presents averaged gene-wise CI bounds as a “Bootstrap 95% CI”
    — it computes a valid global bootstrap CI of the across-gene
    mean/median (shared resampling plan, quantiles of the aggregated
    statistic) and labels it explicitly; `metric="median"` is now
    honored when CIs exist;
    [`plot_sait()`](https://gallardoalba.github.io/TSENAT/reference/plot_sait.md)
    plots the STORED diversity results (via
    `.combine_diversity_results_for_sait()`) instead of recomputing with
    `norm=TRUE`; gene identity is never reconstructed by position
    ([`rep()`](https://rdrr.io/r/base/rep.html)) in q-curve plots
    (explicit error instead); IQR ribbons are labeled as descriptive
    spread (“Median ± IQR/2”) and the per-sample bootstrap CI ribbon is
    labeled as a descriptive aggregation, not a CI of the median;
    [`plot_diversity_violin_density()`](https://gallardoalba.github.io/TSENAT/reference/plot_diversity_violin_density.md)
    accepts `q=` and errors instead of silently using the first stored
    q;
    [`plot_expression()`](https://gallardoalba.github.io/TSENAT/reference/plot_expression.md)
    gained `quantity="usage"` (within-gene isoform fractions) vs
    `"abundance"`; heatmaps warn when requested genes cannot be plotted;
    [`plot_divergence_spectrum()`](https://gallardoalba.github.io/TSENAT/reference/plot_divergence_spectrum.md)
    returns the file path invisibly when saving (documented contract);
    SE dimension validation messages corrected (rows=genes,
    columns=samples).
  - **Performance**: removed the O(G×T) transcript scans —
    transcript→gene indices (`split(seq_along(genes), genes)`) are now
    built once and reused across diversity (`.tsallis_row`), divergence
    (`.compute_group_isoform_counts`/`.process_single_gene_div`), the
    global-divergence bootstrap, diversity bootstrap
    (`.bootstrap_diversity_ci`), gene aggregation
    (`.aggregate_counts_to_genes`), shrinkage and JIS summaries, and
    isoform filtering — each gene lookup is now O(1);
    `.calculate_tsallis_entropy()` computes only the requested quantity
    (`what="S"` no longer also evaluates Hill numbers and vice versa,
    ~2× on that path); `.tsallis_row()` extracts each gene block and
    applies pseudocount/effective-length once instead of per sample;
    `.estimate_shrinkage_params()` no longer does O(G²) rowname scans
    (vectorized row variances). Numerically identical outputs (locked by
    the existing test batteries: 2,313 expectations green in the
    affected suites). Benchmark (sequential): 2,000 genes × 20 samples ×
    5 isoforms × 5 q ≈ 2 s for diversity; 300 genes × 4 q divergence ≈ 1
    s.
  - **Performance**: new `tsallis_divergence_vector_cpp()` computes the
    whole q-spectrum from one normalization pass (exact mirror of the R
    semantics: pseudocount normalization, q=0 limit on unclamped
    probabilities, min_prob clamp only for pseudocount=0, KL limit with
    log_base correction, log-space fallback, roundoff clamp);
    `.tsallis_divergence_vector()` now calls it (2.3× on the kernel).
    New `bootstrap_compute_multi_q_cpp()` resamples ONCE per iteration
    and evaluates ALL q on the same resample (buffers preallocated
    outside the loop); the multi-q diversity bootstrap uses it as a fast
    path (`percentile`, unpaired, read-level), preserving the legacy
    per-q pipeline as an exact fallback for QC regeneration and
    degenerate cases (3.1× with 5 q; grows with the number of q). Point
    estimates bit-identical; bootstrap CIs now share one resampling plan
    across q, preserving the joint correlation structure. All affected
    suites green (2,917 expectations).
  - **Tables and data**: robust p/effect-size formatting in the vignette
    tables (no `0.00e+00`/`NA%`), top-10 concordance table,
    method-estimand table in README, and regenerated
    `inst/extdata/analysis_sait.rds`.

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
