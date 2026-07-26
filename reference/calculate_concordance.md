# Compare method concordance for differential analysis results

Compares statistical results from two different methods (typically
SAIT/GAM for continuous data and Conover-Iman Rank Transform tests) to
assess agreement and identify genes detected by one method but not the
other.

## Usage

``` r
calculate_concordance(analysis_sait, analysis_rank = NULL, ...)

# S4 method for class 'TSENATAnalysis'
calculate_concordance(
  analysis_sait,
  analysis_rank = NULL,
  verbose = FALSE,
  output_file = NULL,
  ...
)
```

## Arguments

- analysis_sait:

  `TSENATAnalysis` object containing SAIT/GAM analysis results (from
  [`calculate_sait()`](https://gallardoalba.github.io/TSENAT/reference/calculate_sait.md)).

- analysis_rank:

  `TSENATAnalysis` object or NULL. If NULL, uses legacy single-object
  API with analysis_sait containing both results. If provided, compares
  SAIT results from analysis_sait with rank-test results from
  analysis_rank.

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

- sait_method:

  Method name used for SAIT/GAM analysis

- rank_method:

  Method name used for rank-based analysis

- timestamp:

  When concordance was computed

## Details

Compares results from two different statistical methods (typically GAM
for continuous and Conover-Iman Rank Transform for rank-based analysis)
on the same data. Identifies: - Genes significant in both methods (high
confidence) - Genes detected by one method only (potential false
positives or method-specific signal) - Spearman correlation of p-values
(overall agreement trends)

## Examples

``` r
# Compare results from SAIT and rank-based testing
# (Requires pre-computed analysis objects from calculate_sait and calculate_rank_transform)
# results_df <- results(calculate_concordance(analysis_sait, analysis_rank))
```
