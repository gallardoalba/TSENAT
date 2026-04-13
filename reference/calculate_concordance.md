# Compare method concordance for differential analysis results

Compares statistical results from two different methods (typically
LM/GAM for continuous data and Scheirer-Ray-Hare rank tests) to assess
agreement and identify genes detected by one method but not the other.

## Usage

``` r
calculate_concordance(analysis_lm, analysis_rank = NULL, ...)

# S4 method for class 'TSENATAnalysis'
calculate_concordance(
  analysis_lm,
  analysis_rank = NULL,
  verbose = FALSE,
  output_file = NULL
)
```

## Arguments

- analysis_lm:

  `TSENATAnalysis` object containing LM/GAM analysis results (from
  [`calculate_lm()`](https://gallardoalba.github.io/TSENAT/reference/calculate_lm.md)).

- analysis_rank:

  `TSENATAnalysis` object or NULL. If NULL, uses legacy single-object
  API with analysis_lm containing both results. If provided, compares LM
  results from analysis_lm with rank-test results from analysis_rank.

- ...:

  Additional arguments for future extensibility.

- verbose:

  `logical`. Print progress messages (default: FALSE).

- output_file:

  `character` or NULL. Optional file path to save results. Supported
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

- lm_method:

  Method name used for LM/GAM analysis

- rank_method:

  Method name used for rank-based analysis

- timestamp:

  When concordance was computed

## Details

Compares results from two different statistical methods (typically GAM
for continuous and Scheirer-Ray-Hare for rank-based analysis) on the
same data. Identifies: - Genes significant in both methods (high
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
config <- TSENAT_config(
  sample_col = 'sample',
  condition_col = 'condition',
  subject_col = 'paired_samples',
  paired = TRUE,
  control = 'normal',
  q = seq(0, 2, by = 0.1)
)

# Build analysis with configured parameters and metadata as explicit parameter
analysis <- build_analysis(
  readcounts = readcounts,
  tx2gene = gff3_dataset,
  metadata = metadata_df,
  config = config,
  tpm = tpm,
  effective_length = effective_length
)

analysis <- filter_analysis(analysis, stringency = 'severe')
analysis <- calculate_diversity(analysis, q = c(0.5, 1.0, 1.5, 2.0, 2.5))
analysis <- calculate_divergence(analysis, q = c(0.5, 1.0, 1.5, 2.0, 2.5))
analysis <- suppressWarnings(calculate_lm(analysis, method = 'gam'))
# Note: calculate_concordance requires results from both
# calculate_srh and calculate_assumptions
```
