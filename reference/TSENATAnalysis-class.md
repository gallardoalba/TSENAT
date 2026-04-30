# TSENATAnalysis S4 Class Definition

Central container for unified analysis workflows in TSENAT.

## Usage

``` r
# S4 method for class 'TSENATAnalysis,ANY,ANY,ANY'
x[i, j, drop = TRUE]
```

## Arguments

- x:

  A `TSENATAnalysis` object to subset.

- i:

  Gene indices (numeric, logical, or character vector). Defaults to all
  genes.

- j:

  Sample indices (numeric, logical, or character vector). Defaults to
  all samples.

- drop:

  Ignored for TSENATAnalysis objects; included for S4 method signature
  compatibility.

## Value

A new `TSENATAnalysis` object containing only the specified genes and
samples, with all associated results and metadata preserved.

## Details

The `TSENATAnalysis` class encapsulates all components of a complete
TSENAT analysis: raw data, configuration metadata, and results from each
analytical step. This unified object ensures metadata is never lost
through the analysis pipeline and provides consistent accessor methods
for result retrieval.

Access results via the unified `results(obj, type = ...)` accessor
method: - `type='diversity'` for Tsallis entropy across q-values -
`type = "sait"` for regularized/penalized regression (GAM, LMM, GEE,
FPCA) interaction results - `type='rank_test'` for Scheirer-Ray-Hare
rank-based test results - `type='divergence'` for divergence metrics -
`type='jackknife'` for jackknife resampling results Use `metadata(obj)`
to access reproducibility metadata.

## Slots

- `se`:

  `SummarizedExperiment`. The base expression data object (genes x
  samples) with assays and colData.

- `config`:

  `list`. Configuration metadata specifying analysis parameters that
  persist through the workflow (q-values, sample grouping columns,
  etc.). Set once via
  [`TSENAT_config()`](https://gallardoalba.github.io/TSENAT/reference/TSENAT_config.md)
  and used by all downstream wrapper functions.

- `diversity_results`:

  `list`. Named list of diversity calculation results. Each name
  corresponds to a q-value (e.g., 'q_0.5', 'q_1.0'). Values are
  SummarizedExperiment objects or data.frames containing entropy values
  for each gene at that q-value.

- `sait_results`:

  `list`. Complex results from scale-adaptive interaction models (GAM
  with smoothing, LMM with random effects, GEE with working
  correlations, FPCA for functional data) (GAM, LMM, GEE, FPCA) and
  statistical testing. Top-level names identify analysis type:

  `sait_interaction`

  :   Regularized regression (GAM/LMM/GEE/FPCA) model results (list with
      `$results` data.frame, `$models` list, etc.)

  `rank_test`

  :   Scheirer-Ray-Hare rank-based test results

  `divergence_difference`

  :   Differential divergence comparison

- `pairwise_results`:

  `list`. Pairwise group comparison results. Results from running
  statistical tests for group comparisons.

- `rank_test_results`:

  `list`. Scheirer-Ray-Hare rank-based statistical test results. Names
  correspond to q-values (e.g., 'q_0.5', 'q_1.0'). Computed by
  [`calculate_srh()`](https://gallardoalba.github.io/TSENAT/reference/calculate_srh.md)
  as a non-parametric alternative to linear mixed model testing.

- `jackknife_results`:

  `list`. Resampling-based confidence intervals. Names correspond to
  q-values (e.g., 'q_0.5', 'q_1.0'). Values are jackknife result objects
  containing resamples, CI bounds, and diagnostics.

- `divergence_results`:

  `list`. Divergence metric calculations. Typically contains:

  `tsallis_divergence`

  :   SummarizedExperiment with divergence values

  `effect_sizes`

  :   data.frame with Cohen's d, etc.

- `plots`:

  `list`. Cached visualization objects (ggplot). Names identify plot
  type (e.g., 'q_curve', 'sait_interaction', 'influence'). Populated by
  [`TSENAT()`](https://gallardoalba.github.io/TSENAT/reference/TSENAT.md)
  if `generate_plots=TRUE`.

- `metadata`:

  `list`. Reproducibility and tracking metadata. Automatically
  maintained by wrapper functions. Includes:

  `created_at`

  :   Timestamp of object creation

  `function_calls`

  :   Vector of wrapper functions called

  `function_timestamps`

  :   Timestamps for each function call

  `package_version`

  :   TSENAT version at creation
