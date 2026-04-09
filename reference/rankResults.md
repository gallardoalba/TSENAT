# Extract rank test q-value interaction results

Extract rank test q-value interaction results

Setter for rank test q-value interaction results stored in
@lm_results\$rank_test

## Usage

``` r
rankResults(object, component = NULL)

# S4 method for class 'TSENATAnalysis'
rankResults(object, component = NULL)

rankResults(object) <- value

# S4 method for class 'TSENATAnalysis'
rankResults(object) <- value
```

## Arguments

- object:

  TSENATAnalysis object

- component:

  `character`. Component to extract (default: all results).

- value:

  list or data.frame. Rank test q-value interaction results.

## Value

List or data.frame of rank test q-value interaction results, or NULL if
not computed.

## Details

Rank test results from
[`rank_test_q_condition_s4()`](https://gallardoalba.github.io/TSENAT/reference/rank_test_q_condition_s4.md)
are retrieved via this method. LM interaction results are retrieved
separately with
[`lmResults()`](https://gallardoalba.github.io/TSENAT/reference/lmResults.md).

## Examples

``` r
# Extract rank test results from TSENATAnalysis object
data(readcounts, package = 'TSENAT')
metadata <- read.table(
  system.file('extdata', 'metadata.tsv', package = 'TSENAT'),
  header = TRUE, sep = '\t'
)
gff3_file <- system.file('extdata', 'annotation.gff3.gz', package = 'TSENAT')

# TPM and effective_length REQUIRED for filter_analysis_s4()
tpm <- matrix(runif(nrow(readcounts) * ncol(readcounts), 0.1, 10),
              nrow = nrow(readcounts), ncol = ncol(readcounts),
              dimnames = dimnames(readcounts))
effective_length <- matrix(100, nrow = nrow(readcounts), ncol = ncol(readcounts))

config <- tsenat_config(q_values = c(0.5, 1.0), generate_plots = FALSE)
analysis <- build_analysis_s4(readcounts, tx2gene = gff3_file,
    metadata = metadata, tpm = tpm, effective_length = effective_length,
    config = config)
analysis <- filter_analysis_s4(analysis, stringency = 'severe')
analysis <- calculate_diversity_s4(analysis, norm = TRUE)
#> Note: 4 genes excluded (< 75% valid values).
analysis <- rank_test_q_condition_s4(analysis, q = 1.0)

# Retrieve all rank test results
results <- rankResults(analysis)
head(results)
#>              gene n_q_values_tested  f_statistic   p_value adj_p_value
#> 1            VRK2                 2 1.359619e-03 0.9708479           1
#> 2        FAM114A2                 2 6.051437e-03 0.9385476           1
#> 3         POLR2J4                 2 9.277872e-04 0.9759524           1
#> 4          TM7SF2                 2 2.980875e-03 0.9570009           1
#> 5 ENSG00000291072                 2 1.137271e-30 1.0000000           1
#> 6           MEF2A                 2 9.421371e-05 0.9923297           1
#>   ss_interaction ss_residual df_interaction df_residual effect_size_eta2
#> 1   2.279276e-02   0.4228992              1          30     0.0511401571
#> 2   1.412449e-03   1.1332383              1          30     0.0012448319
#> 3   8.311925e-04   0.7673871              1          30     0.0010819744
#> 4   6.505697e-04   1.2394128              1          30     0.0005246262
#> 5   5.851360e-05   0.2569950              1          30     0.0002276320
#> 6   3.742339e-06   0.2760756              1          30     0.0000135553
#>   interaction_class  test_method heteroscedastic boundary_clustered
#> 1   Robust across q srh_unpaired           FALSE              FALSE
#> 2   Robust across q srh_unpaired           FALSE              FALSE
#> 3   Robust across q srh_unpaired           FALSE              FALSE
#> 4   Robust across q srh_unpaired           FALSE              FALSE
#> 5   Robust across q srh_unpaired           FALSE              FALSE
#> 6   Robust across q srh_unpaired           FALSE              FALSE
#>   highly_skewed
#> 1         FALSE
#> 2         FALSE
#> 3         FALSE
#> 4         FALSE
#> 5         FALSE
#> 6         FALSE
```
