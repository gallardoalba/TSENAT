# Compute concordance between two analysis methods in TSENATAnalysis

Compute concordance between two analysis methods in TSENATAnalysis

## Usage

``` r
compute_method_concordance_s4(analysis, ...)

# S4 method for class 'TSENATAnalysis'
compute_method_concordance_s4(
  analysis,
  gam_method = "rank_test",
  friedman_method = "rankbased",
  gam_results = NULL,
  verbose = FALSE,
  output_file = NULL
)
```

## Arguments

- analysis:

  `TSENATAnalysis` object with LM results (e.g., GAM).

- ...:

  Additional arguments for future extensibility.

- gam_method:

  `character`. Key for GAM/interaction results in `@lm_results`.
  Default: 'rank_test' (results from `rank_test_q_condition_s4`)

- friedman_method:

  `character`. Key for Friedman/rank-based results in `@lm_results`.
  Default: 'rankbased' (results from `test_rankbased_assumptions_s4`)

- gam_results:

  `data. frame` or `NULL`. Optional GAM results data frame to store in
  the analysis object. If provided, automatically stored in
  `@lm_results` under the key specified by `gam_method`. Useful for
  importing external results or results computed outside the S4 wrapper.
  Default: NULL (use existing results in analysis).

- verbose:

  `logical`. Print progress messages (default: FALSE).

- output_file:

  `character` or `NULL`. Optional file path to save results. Supported
  formats: .rds (for S4 objects). Default: NULL (no file output).

## Value

Modified TSENATAnalysis object with concordance results stored in:
`@metadata$method_concordance`:

- comparison_df:

  Data frame comparing results from both methods

- spearman_rho:

  Spearman correlation between adjusted p-values

- high_confidence:

  Genes with strong agreement

- agreement_table:

  Contingency table of significant/non-significant calls

- gam_method:

  Method name used for GAM analysis

- friedman_method:

  Method name used for Friedman analysis

- timestamp:

  When concordance was computed

## Details

Compares results from two different statistical methods (typically GAM
for continuous and Friedman/Kruskal-Wallis for rank-based analysis) on
the same data. Identifies: - Genes significant in both methods (high
confidence) - Genes detected by one method only (potential false
positives or method-specific signal) - Spearman correlation of p-values
(overall agreement trends)

## Examples

``` r
# Load example data (matching TSENAT.Rmd workflow)
data(readcounts)
readcounts <- as.matrix(readcounts)
mode(readcounts) <- 'numeric'
metadata_df <- read.table(
  system.file('extdata', 'metadata.tsv', package = 'TSENAT'),
  header = TRUE, sep = '\t'
)
gff3_dataset <- system.file('extdata', 'annotation.gff3.gz', package =
'TSENAT')

# Configure analysis parameters first (fail-fast principle)
config <- tsenat_config(
  sample_col = 'sample',
  condition_col = 'condition',
  subject_col = 'paired_samples',
  paired = TRUE,
  control = 'normal',
  q_values = seq(0, 2, by = 0.1)
)

# Build analysis with configured parameters and metadata as explicit parameter
analysis <- build_analysis_s4(
  readcounts = readcounts,
  tx2gene = gff3_dataset,
  metadata = metadata_df,
  config = config,
  tpm = tpm,
  effective_length = effective_length
)

analysis <- filter_analysis_s4(analysis, stringency = 'severe')
analysis <- calculate_diversity_s4(analysis, q = c(0.5, 1.0, 1.5, 2.0, 2.5))
#> Note: 7 genes excluded (< 75% valid values).
analysis <- calculate_divergence_s4(analysis, q = c(0.5, 1.0, 1.5, 2.0, 2.5))
analysis <- calculate_lm_interaction_s4(analysis, method = 'gam')
#> Warning: nlminb problem, convergence error code = 1
#>   message = singular convergence (7)
#> Warning: nlminb problem, convergence error code = 1
#>   message = iteration limit reached without convergence (10)
#> Warning: nlminb problem, convergence error code = 1
#>   message = iteration limit reached without convergence (10)
# Note: compute_method_concordance_s4 requires results from both
# rank_test_q_condition_s4 and test_rankbased_assumptions_s4
```
