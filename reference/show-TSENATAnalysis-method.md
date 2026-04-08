# Show method for TSENATAnalysis

Show method for TSENATAnalysis

Display TSENATAnalysis object

## Usage

``` r
# S4 method for class 'TSENATAnalysis'
show(object)

# S4 method for class 'TSENATAnalysis'
show(object)
```

## Arguments

- object:

  `TSENATAnalysis` object to display.

## Value

Invisibly returns the `TSENATAnalysis` object (called for its side
effect of printing a formatted summary to the console).

## Details

Provides a concise summary of: data dimensions, completed analyses,
number of results, and metadata tracking.

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
show(analysis)
#> TSENATAnalysis object
#> =====================
#> SummarizedExperiment:
#>   Genes:   200
#>   Samples: 16
#> 
#> Configuration:
#>   q_values: 0, 0.5, 1, 1.5, 2
#>   condition_col: condition
#>   subject_col: NULL
#>   sample_col: sample
#>   paired: FALSE
#>   control: NULL
#>   p_threshold: 0.05
#>   fdr_threshold: 0.05
#>   significance_threshold: 0.05
#>   nboot: 1000
#>   bootstrap_method: percentile
#>   stringency: medium
#>   nthreads: 1
#>   norm: TRUE
#>   bootstrap: FALSE
#>   bootstrap_ci: 0.95
#>   norm_method: NULL
#>   pseudocount: NULL
#>   metadata: <data.frame>
#> 
#> Analysis Status:
#> 
```
