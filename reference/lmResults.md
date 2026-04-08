# Extract linear model interaction results

Extract linear model interaction results

## Usage

``` r
lmResults(object, component = NULL)

# S4 method for class 'TSENATAnalysis'
lmResults(object, component = NULL)

lmResults(object) <- value

# S4 method for class 'TSENATAnalysis'
lmResults(object) <- value
```

## Arguments

- object:

  `TSENATAnalysis` object.

- component:

  `character`. Which result component to extract. Options: NULL (all LM
  interaction results), 'lm_interaction', 'lm_interaction_model_data',
  'results', 'p_value', 'effect_size', etc.

- value:

  A list of LM interaction results to assign to the object.

## Value

List or data.frame of LM interaction results depending on component
requested.

## Details

Returns only LM interaction results stored in the `@lm_results` slot.
Note: Rank test q-value interaction results are retrieved separately via
[`rankResults()`](https://gallardoalba.github.io/TSENAT/reference/rankResults.md).
Use `lmResults(analysis)` to get all LM interaction components, or
specify component type for targeted extraction.

## Examples

``` r
# Load example data and run LM interaction analysis
data(readcounts)
readcounts <- as.matrix(readcounts)
mode(readcounts) <- 'numeric'
metadata_df <- read.table(
  system.file('extdata', 'metadata.tsv', package = 'TSENAT'),
  header = TRUE, sep = '\t'
)
gff3_dataset <- system.file('extdata', 'annotation.gff3.gz', package =
'TSENAT')

# Build analysis from vignette data
config <- tsenat_config(sample_col = 'sample', condition_col = 'condition')
analysis <- build_analysis_s4(readcounts = readcounts, tx2gene =
gff3_dataset, metadata = metadata_df, config = config,
  tpm = tpm, effective_length = effective_length)
analysis <- filter_analysis_s4(analysis, min_samples = 1, subset_n_genes
= 200)

# Compute diversity first (required for LM interaction)
analysis <- calculate_diversity_s4(analysis, q = c(0.5, 1.0, 1.5, 2.0, 2.5), verbose =
FALSE)

# Calculate LM interaction
analysis <- calculate_lm_interaction_s4(analysis, 
  condition_col = 'condition', method = 'gam')

# Extract and view LM interaction results
res <- lmResults(analysis, 'lm_interaction')
if (!is.null(res)) head(res, 3)
#>                 gene p_interaction        p_raw n_observations n_subjects
#> 1 ENSG00000107562.18  1.627847e-14 1.627847e-14             64         16
#> 2 ENSG00000140274.15  4.253465e-09 4.253465e-09             56         16
#> 3 ENSG00000057663.17  9.241246e-05 9.241246e-05             64         16
#>   n_effective    rho_ar1 test_statistic effect_size df_residual model_converged
#> 1   14.834530 0.05062463       63.49787   0.8651259    58.05485            TRUE
#> 2   12.754919 0.16228978       38.55106   0.7441343    51.03222            TRUE
#> 3    9.586938 0.35000000       18.57850   0.8466450    58.04111            TRUE
#>    slope_diff arima_transformation bounded_support_model family_used
#> 1 -0.05825786                 TRUE                 FALSE    Gaussian
#> 2  0.07862503                 TRUE                 FALSE    Gaussian
#> 3 -0.01899569                 TRUE                 FALSE    Gaussian
#>   heteroscedasticity_detected variance_ratio_q              fit_method
#> 1                       FALSE         2.092749 mgcv::gamm_arima(1,1,0)
#> 2                       FALSE         1.000000 mgcv::gamm_arima(1,1,0)
#> 3                       FALSE         2.272165 mgcv::gamm_arima(1,1,0)
#>   shapiro_p_value residuals_normal n_residuals_tested ci_weighted
#> 1    2.565713e-06            FALSE                 64       FALSE
#> 2    1.270953e-05            FALSE                 56       FALSE
#> 3    1.560005e-01             TRUE                 64       FALSE
#>   adj_p_interaction          gene_name            gene_id
#> 1      1.237163e-12 ENSG00000107562.18 ENSG00000107562.18
#> 2      3.190099e-07 ENSG00000140274.15 ENSG00000140274.15
#> 3      6.838522e-03 ENSG00000057663.17 ENSG00000057663.17
```
