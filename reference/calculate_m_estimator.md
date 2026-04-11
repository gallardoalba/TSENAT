# M-Estimation for Sample Quality (S4 Wrapper)

S4 wrapper for `m_estimate` that performs robust M-estimation on
diversity results stored in a TSENATAnalysis object and stores results
back into the object.

## Usage

``` r
calculate_m_estimator(
  analysis,
  condition_col = NULL,
  loss_type = "huber",
  scale = NULL,
  max_iter = 50,
  tol = 1e-06,
  paired = NULL,
  pcorr = "BH",
  q_combine_method = "mean",
  influence_threshold = 0.75,
  scale_method = "mad",
  output_file = NULL,
  verbose = FALSE
)
```

## Arguments

- analysis:

  `TSENATAnalysis` object with diversity results (typically via
  [`calculate_diversity`](https://gallardoalba.github.io/TSENAT/reference/calculate_diversity.md)).

- condition_col:

  `character`. Column name in sample metadata indicating
  condition/sample grouping. Auto-detected from `@config$condition_col`
  if available.

- loss_type:

  `character`. Type of loss function: 'huber' (default), 'tukey', or
  'lsq'. Determines robustness vs efficiency trade-off.

- scale:

  `numeric`. Manual scale parameter. If NULL, estimated from data.

- max_iter:

  `integer`. Maximum iterations for M-estimation. Default: 50.

- tol:

  `numeric`. Convergence tolerance. Default: 1e-6.

- paired:

  `logical`. If TRUE, adjusts degrees of freedom for paired designs.
  Auto-detected from `@config$paired` if available. Default: FALSE.

- pcorr:

  `character`. P-value correction method. Default: 'BH'
  (Benjamini-Hochberg).

- q_combine_method:

  `character`. How to collapse multi-q results: 'mean' (default) or
  'median'.

- influence_threshold:

  `numeric`. Threshold for classifying samples as high-influence.
  Default: 0.75.

- scale_method:

  `character`. Scale estimation method: 'mad' (default), 'proposal2', or
  's-estimator'.

- output_file:

  `character` or `NULL`. Optional file path to save results. Supported
  formats: .rds (for S4 objects), .tsv, .csv, .txt (for tables).
  Default: NULL (no file output).

- verbose:

  `logical`. Print status messages. Default: TRUE.

## Value

Modified TSENATAnalysis object with M-estimation results stored in
`analysis@metadata$m_estimate_results`. Contains data frame with
influence scores, robustness weights, entropy statistics, and QC
classifications. Returns visibly to support method chaining and piping.

## Details

This wrapper extracts diversity results from
`analysis@diversity_results`, performs robust M-estimation on diversity
values (entropy), and stores results in the analysis object metadata.

\*\*M-Estimation:\*\* Robust regression technique that down-weights
outliers based on their residuals. Useful for detecting low-quality
samples that show unusual diversity patterns.

\*\*M-Estimation Results include:\*\*

- `sample_influence`: How much each sample affects the overall fit

- `robustness_weight`: Down-weighting factor (lower = more outlying)

- `entropy_mean`: Average entropy for the sample

- `entropy_sd`: Entropy variability within the sample

- `Status`: QC Classification ('OK' or 'Flag for QC' based on
  influence_threshold)

\*\*Parameter resolution priority\*\* (explicit \> @config \> error):

- `samples`: Uses explicit arg, else `@config$condition_col`, else
  error. Note: despite parameter name 'samples', maps to condition
  grouping column

\*\*Data Requirements:\*\*

- Diversity results must be computed via
  [`calculate_diversity()`](https://gallardoalba.github.io/TSENAT/reference/calculate_diversity.md)

- Sample grouping column required in colData (auto-detected from
  @config\$condition_col or via 'samples' parameter)

## See also

[`calculate_diversity`](https://gallardoalba.github.io/TSENAT/reference/calculate_diversity.md)
for computing diversity

## Examples

``` r
# Create test analysis and compute M-estimation
data(readcounts)
readcounts <- as.matrix(readcounts)
mode(readcounts) <- 'numeric'
metadata_df <- read.table(
  system.file('extdata', 'metadata.tsv', package = 'TSENAT'),
  header = TRUE, sep = '\t'
)
gff3_dataset <- system.file('extdata', 'annotation.gff3.gz', package =
'TSENAT')

# Create config (metadata passed as explicit parameter to build_analysis)
config <- TSENAT_config(
  sample_col = 'sample',
  condition_col = 'condition',
  q_values = seq(0, 2, by = 0.05),
  paired = FALSE
)

# Build analysis from vignette data
analysis <- build_analysis(
  readcounts = readcounts,
  tx2gene = gff3_dataset,
  metadata = metadata_df,
  config = config,
  tpm = tpm,
  effective_length = effective_length
)
analysis <- filter_analysis(analysis, min_samples = 1, subset_n_genes = 200)
analysis <- calculate_diversity(analysis, q = c(0.5, 1.0, 1.5))
#> Note: 37 genes excluded (< 75% valid values).
analysis <- calculate_m_estimator(
  analysis,
  condition_col = 'condition',
  loss_type = 'huber'
)
```
