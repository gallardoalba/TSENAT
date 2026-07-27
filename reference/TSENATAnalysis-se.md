# Extract SummarizedExperiment from TSENATAnalysis

Extract SummarizedExperiment from TSENATAnalysis

## Usage

``` r
se(object)

# S4 method for class 'TSENATAnalysis'
se(object)
```

## Arguments

- object:

  `TSENATAnalysis` object.

## Value

SummarizedExperiment object

The `SummarizedExperiment` object stored in `@se` slot, containing raw
count data and sample metadata.

## Details

Provides read-only access to the underlying SummarizedExperiment. To
modify the SummarizedExperiment, use the configuration methods or
directly access `object@se`.

## Examples

``` r
# Create a simple TSENATAnalysis object
library(SummarizedExperiment)
counts <- matrix(rpois(100, lambda = 10), nrow = 10, ncol = 10)
colnames(counts) <- paste0('sample_', 1:10)
rownames(counts) <- paste0('gene_', 1:10)

se <- SummarizedExperiment(assays = list(counts = counts))
analysis <- new('TSENATAnalysis', se = se, config = list())

# Extract the SummarizedExperiment
extracted_se <- se(analysis)
dim(extracted_se)
#> [1] 10 10
```
