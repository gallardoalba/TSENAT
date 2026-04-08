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
#> TSENATAnalysis Object
#> ====================
#> Genes:       200
#> Samples:     16
#> Configuration: 26 parameters
#> Analysis status: EMPTY
#> Created: 2026-04-09 00:07:47.491723
#> 
```
