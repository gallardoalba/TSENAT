# Extract jackknife isoform switching results

Extract jackknife isoform switching results

## Usage

``` r
jisResults(object, q = NULL)

# S4 method for class 'TSENATAnalysis'
jisResults(object, q = NULL)
```

## Arguments

- object:

  `TSENATAnalysis` object.

- q:

  `numeric` or NULL. Q-value for specific results. If NULL, returns all
  isoform switching jackknife results.

## Value

List of isoform switching jackknife results, or NULL if not computed.

## Details

Jackknife isoform switching results are computed separately from entropy
outliers. Use this to access isoform switching analysis results.

Results include leave-one-out diagnostics for detecting genes with
condition-specific isoform switching patterns across q-values.

## Examples

``` r
# Load data and build analysis 
data(readcounts, package = 'TSENAT')
metadata_df <- read.table(system.file('extdata', 'metadata.tsv', package = 'TSENAT'),
  header = TRUE, sep = '\t')
gff3_file <- system.file('extdata', 'annotation.gff3.gz', package = 'TSENAT')

# TPM and effective_length REQUIRED for filter_analysis_s4()
tpm <- matrix(runif(nrow(readcounts) * ncol(readcounts), 0.1, 10),
              nrow = nrow(readcounts), ncol = ncol(readcounts),
              dimnames = dimnames(readcounts))
effective_length <- matrix(100, nrow = nrow(readcounts), ncol = ncol(readcounts))

# Build analysis object
config <- tsenat_config(sample_col = 'sample', condition_col = 'condition')
analysis <- build_analysis_s4(readcounts = readcounts,
                             tx2gene = gff3_file,
                             metadata = metadata_df,
                             tpm = tpm,
                             effective_length = effective_length,
                             config = config)

# Filter low-abundance transcripts
analysis <- filter_analysis_s4(analysis, stringency = 'medium')

# Calculate diversity
analysis <- calculate_diversity_s4(analysis, q = c(1.0, 2.0), norm = TRUE)
#> Note: 117 genes excluded (< 75% valid values).

# Run jackknife isoform switching analysis
analysis <- jackknife_isoform_switching_s4(analysis, q = 1.0)
#> Warning: no non-missing arguments to max; returning -Inf
#> Warning: no non-missing arguments to max; returning -Inf
#> Warning: no non-missing arguments to max; returning -Inf
#> Warning: no non-missing arguments to max; returning -Inf
#> Warning: no non-missing arguments to max; returning -Inf
#> Warning: no non-missing arguments to max; returning -Inf
#> Warning: no non-missing arguments to max; returning -Inf
#> Warning: no non-missing arguments to max; returning -Inf
#> Warning: no non-missing arguments to max; returning -Inf
#> Warning: no non-missing arguments to max; returning -Inf
#> Warning: no non-missing arguments to max; returning -Inf
#> Warning: no non-missing arguments to max; returning -Inf
#> Warning: no non-missing arguments to max; returning -Inf
#> Warning: no non-missing arguments to max; returning -Inf
#> Warning: no non-missing arguments to max; returning -Inf
#> Warning: no non-missing arguments to max; returning -Inf
#> Warning: no non-missing arguments to max; returning -Inf
#> Warning: no non-missing arguments to max; returning -Inf
#> Warning: no non-missing arguments to max; returning -Inf
#> Warning: no non-missing arguments to max; returning -Inf
#> Warning: no non-missing arguments to max; returning -Inf
#> Warning: no non-missing arguments to max; returning -Inf
#> Warning: no non-missing arguments to max; returning -Inf
#> Warning: no non-missing arguments to max; returning -Inf
#> Warning: no non-missing arguments to max; returning -Inf

# Extract isoform switching results for q = 1.0
jis_results <- jisResults(analysis, q = 1.0)
#> Warning: No jackknife isoform switching results found (multi_q key missing). Run jackknife_isoform_switching_s4() first.

# Extract all isoform switching results
all_jis <- jisResults(analysis)
#> Warning: No jackknife isoform switching results found (multi_q key missing). Run jackknife_isoform_switching_s4() first.
```
