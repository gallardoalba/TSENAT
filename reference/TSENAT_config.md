# Create and return TSENAT configuration

Builds a configuration list for use with
[`TSENAT`](https://gallardoalba.github.io/TSENAT/reference/TSENAT.md)().
Allows specifying analysis parameters once and reusing across multiple
analyses.

## Usage

``` r
TSENAT_config(
  q = 1,
  condition_col = "condition",
  subject_col = NULL,
  sample_col = "sample",
  paired = FALSE,
  control = NULL,
  p_threshold = 0.05,
  fdr_threshold = 0.05,
  significance_threshold = 0.05,
  bootstrap = FALSE,
  nboot = 1000,
  bootstrap_method = "percentile",
  stringency = "medium",
  nthreads = 1,
  norm = TRUE,
  bootstrap_ci = 0.95,
  bootstrap_include_diagnostics = TRUE,
  min_valid_frac = 0.75,
  norm_method = NULL,
  pseudocount = 0,
  shrinkage = "none",
  lm_method = "gam",
  lm_pcorr = "BH",
  jis_use_lm_fdr = TRUE,
  divergence_ci = 0.95,
  assumptions_checks = "all",
  ...
)
```

## Arguments

- q:

  `numeric`. Q-value(s) for Tsallis entropy (single value or vector).
  Default: 1.0 (Shannon entropy). Usage:
  `calculate_diversity/divergence` use this for spectrum computation (if
  vector) or as default fallback (if single).

- condition_col:

  `character`. Column name in colData containing conditions. Default:
  'condition'.

- subject_col:

  `character`. Column name in colData containing subject IDs (for paired
  designs). Default: NULL.

- sample_col:

  `character`. Column name in colData containing sample IDs. Default:
  'sample'.

- paired:

  `logical`. Whether samples are paired/repeated measures. Default:
  FALSE.

- control:

  `character`. Reference/control group label. Default: NULL.

- p_threshold:

  `numeric`. Raw p-value threshold. Default: 0.05.

- fdr_threshold:

  `numeric`. FDR-adjusted p-value threshold. Default: 0.05.

- significance_threshold:

  `numeric`. Significance cutoff. Default: 0.05.

- bootstrap:

  `logical`. Enable bootstrap CIs. Default: FALSE.

- nboot:

  `integer`. Bootstrap resamples for CIs. Default: 1000.

- bootstrap_method:

  `character`. Bootstrap method: 'percentile' or 'bca'. Default:
  'percentile'.

- stringency:

  `character`. Filtering stringency: 'lenient', 'medium', 'severe'.
  Default: 'medium'.

- nthreads:

  `integer`. Parallel threads. Default: 1.

- norm:

  `logical`. Enable normalization. Default: TRUE.

- bootstrap_ci:

  `numeric`. Confidence level (0-1). Default: 0.95.

- bootstrap_include_diagnostics:

  `logical`. Include diagnostics. Default: TRUE.

- min_valid_frac:

  `numeric`. Min valid replicate fraction. Default: 0.75.

- norm_method:

  `character`. Normalization: NULL, 'zscore', 'log_odds_ratio',
  'relative_reference'. Default: NULL.

- pseudocount:

  `numeric`. Pseudocount for sparse data. Default: 0.

- shrinkage:

  `character`. Variance reduction: 'none' or 'empirical_bayes'. Default:
  'none'.

- lm_method:

  `character`. LM method: 'gam', 'lmm', 'fpca', 'gee'. Default: 'gam'.

- lm_pcorr:

  `character`. P-value correction: 'BH', 'bonferroni', 'hochberg',
  'holm'. Default: 'BH'.

- jis_use_lm_fdr:

  `logical`. Filter jackknife genes using LM p-values. Default: TRUE.

- divergence_ci:

  `numeric`. Confidence level for divergence CIs. Default: 0.95.

- assumptions_checks:

  `character`. Which assumptions to test (default: 'all'). Presets: -
  'rank': core assumption checks (exchangeability, monotonicity,
  consistency) - 'all': all checks including method-specific diagnostics
  (GAM, GEE, LMM, FPCA) Explicit: character vector like
  `c('exchangeability', 'monotonicity')`.

- ...:

  Additional configuration parameters (stored as-is).

## Value

`list` with class `TSENATConfig` containing all specified parameters.

## Details

Configuration is stored in the TSENATAnalysis@config slot and used by
wrapper functions to configure analysis behavior. Note: Statistical
tests (Wilcoxon, shuffle) work on a single q-value, so only one q-value
is specified in config.

## Examples

``` r
# Default config with standard parameters (point estimates only)
cfg <- TSENAT_config()

# For Wilcoxon/shuffle tests (single q-value required in config)
cfg <- TSENAT_config(
  q = 1.0,                          # Shannon entropy - for rank tests
  condition_col = 'treatment',
  control = 'untreated'
)

# For Scheirer-Ray-Hare rank tests (multiple q-values)
cfg <- TSENAT_config(
  q = seq(0, 2, by = 0.5),          # Multiple q-values for spectrum or advanced testing
  condition_col = 'treatment',
  control = 'untreated'
)

# With bootstrap CIs for uncertainty quantification (recommended)
cfg <- TSENAT_config(
  bootstrap = TRUE,                # Enable bootstrap confidence intervals
  bootstrap_method = 'bca',         # Bias-corrected (better for skewed entropy)
  nboot = 1000,                     # 1000 resamples
  bootstrap_ci = 0.95               # 95% CI
)

# Custom with paired analysis, strict filtering, and normalization
cfg <- TSENAT_config(
  q = 1.0,                          # Shannon entropy
  condition_col = 'treatment',
  subject_col = 'subject_id',
  paired = TRUE,
  control = 'untreated',
  stringency = 'severe',            # High-confidence transcripts only
  norm_method = 'zscore',           # Cross-study standardization
  shrinkage = 'none',               # Empirical estimates
  bootstrap = TRUE,
  bootstrap_method = 'bca',
  nboot = 5000,                     # Higher precision
  pseudocount = 0,                  # Disabled by default; set > 0 to add pseudocount
  significance_threshold = 0.01     # Stricter significance level
)
```
