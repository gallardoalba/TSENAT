# Extract jackknife resampling results

Extract jackknife resampling results

## Usage

``` r
jeoResults(object, q = NULL)

# S4 method for class 'TSENATAnalysis'
jeoResults(object, q = NULL)
```

## Arguments

- object:

  `TSENATAnalysis` object.

- q:

  `numeric`. Q-value for jackknife results (e.g., 1.0). If NULL
  (default), returns all q-values.

## Value

Jackknife result object (confidence intervals, resamples, etc.).

## Details

Jackknife results are stored per q-value. Use this to access confidence
intervals and diagnostic information from resampling.

## See also

Other TSENATAnalysis accessors:
[`diversity`](https://gallardoalba.github.io/TSENAT/reference/diversity.md),
[`divergence`](https://gallardoalba.github.io/TSENAT/reference/divergence.md),
[`lmResults`](https://gallardoalba.github.io/TSENAT/reference/lmResults.md),
[`se`](https://gallardoalba.github.io/TSENAT/reference/se.md),
[`metadata`](https://rdrr.io/pkg/S4Vectors/man/Annotated-class.html)

## Examples

``` r
# Load real TSENAT data and run jackknife analysis
data(readcounts)
metadata_df <- read.table(system.file('extdata', 'metadata.tsv', package
= 'TSENAT'),
  header = TRUE, sep = '\t')
gff3_file <- system.file('extdata', 'annotation.gff3.gz', package = 'TSENAT')
config <- tsenat_config(sample_col = 'sample', condition_col = 'condition')
analysis <- build_analysis_s4(readcounts = readcounts, tx2gene =
gff3_file, metadata = metadata_df, config = config,
  tpm = tpm, effective_length = effective_length)
analysis <- filter_analysis_s4(analysis, min_samples = 1, subset_n_genes
= 200)
analysis <- calculate_diversity_s4(analysis, q = c(0.5, 1.0, 1.5))
#> Note: 12 genes excluded (< 75% valid values).
analysis <- jackknife_entropy_outliers_s4(analysis, q = c(0.5, 1.0, 1.5))
jk_results <- jeoResults(analysis)
```
