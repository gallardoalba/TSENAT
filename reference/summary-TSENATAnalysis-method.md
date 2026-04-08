# Summary method for TSENATAnalysis

Summary method for TSENATAnalysis

Detailed summary of TSENATAnalysis

## Usage

``` r
# S4 method for class 'TSENATAnalysis'
summary(object)

# S4 method for class 'TSENATAnalysis'
summary(object)
```

## Arguments

- object:

  `TSENATAnalysis` object.

## Value

Prints detailed summary (invisibly returns object).

## Details

Provides comprehensive analysis summary including dimensions, results
counts, validation status, and metadata tracking.

## Examples

``` r
# Load real TSENAT data
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
summary(analysis)
#> === TSENAT Analysis Summary ===
#> DATA:
#>   Genes:       200
#>   Samples:      16
#>   Assays:        1 (counts)
#> 
#> CONFIGURATION:
#>   q_values: 0, 0.5, 1, 1.5, 2
#>   condition_col: condition
#>   subject_col: <NULL>
#>   sample_col: sample
#>   paired: <logical>
#>   control: <NULL>
#>   p_threshold: 0.05
#>   fdr_threshold: 0.05
#>   significance_threshold: 0.05
#>   nboot: 1000
#>   bootstrap_method: percentile
#>   stringency: medium
#>   nthreads: 1
#>   norm: <logical>
#>   bootstrap: <logical>
#>   bootstrap_ci: 0.95
#>   norm_method: <NULL>
#>   pseudocount: <NULL>
#>   metadata: <data.frame>
#> 
#> RESULTS:
#> 
#> METADATA:
#>   Created: 2026-04-08 04:30:24
#>   Package: TSENAT 0.99.0
#> 
```
