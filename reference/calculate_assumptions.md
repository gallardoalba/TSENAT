# Test statistical assumptions on diversity data in TSENATAnalysis

Test statistical assumptions on diversity data in TSENATAnalysis

## Usage

``` r
calculate_assumptions(
  analysis,
  q = NULL,
  checks = "rank",
  alpha = 0.05,
  format = "text",
  ...
)

# S4 method for class 'TSENATAnalysis'
calculate_assumptions(
  analysis,
  q = NULL,
  checks = "rank",
  alpha = 0.05,
  format = "text",
  ...
)
```

## Arguments

- analysis:

  `TSENATAnalysis` object with diversity results stored in
  `@diversity_results`.

- q:

  `numeric`. Q-value(s) to extract from diversity results. If NULL, uses
  the first available diversity result or q=1.0.

- checks:

  `character`. Which assumptions to test (default: 'rank'). Presets: -
  'rank': core assumption checks (exchangeability, monotonicity,
  consistency) - 'all': all checks including GAM diagnostics Explicit:
  character vector like `c('exchangeability', 'monotonicity')`.

- alpha:

  `numeric`. Significance level for tests (default: 0.05).

- format:

  `character`. Output format when used with
  [`results()`](https://gallardoalba.github.io/TSENAT/reference/results.md).
  'text' (default): formatted text output for display 'list': returns
  structured list for programmatic access.

- ...:

  Additional arguments (output_file, verbose for file output).

## Value

Modified TSENATAnalysis object with assumption test results stored in
`@metadata$rankbased_assumptions`.

## Details

This wrapper calls `.calculate_assumptions()` on diversity data
extracted from the analysis object. Evaluates data stability and
consistency across dimensions. Results include:

- exchangeability:

  Permutation test for independence and temporal/spatial structure

- monotonicity:

  Spearman correlation consistency across rows

- consistency:

  Kendall's W concordance and ICC across samples

\*\*Data Extraction Priority:\*\* 1. If q specified: uses diversity
result for that q-value 2. If q NULL: uses first available diversity
result 3. If no diversity results: extracts from cached combined result
(`@metadata$diversity_combined`)

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

# Build analysis from vignette data and create small subset
config <- TSENAT_config(
  sample_col = 'sample',
  condition_col = 'condition',
  q_values = seq(0, 2, by = 0.05),
  paired = FALSE
)
analysis <- build_analysis(
  readcounts = readcounts,
  metadata = metadata_df,
  tx2gene = gff3_dataset,
  config = config,
  tpm = tpm,
  effective_length = effective_length
)
analysis <- filter_analysis(analysis, min_samples = 1, subset_n_genes = 200)
analysis <- calculate_diversity(analysis, q = c(0.5, 1.0, 1.5))
analysis <- calculate_assumptions(analysis, q = 1.0)
# Check results using rank_test accessor
results_df <- results(analysis, type = 'rank_test')
```
