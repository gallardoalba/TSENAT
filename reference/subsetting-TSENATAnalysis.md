# Subset TSENATAnalysis Objects

Extract a subset of genes and/or samples from a TSENATAnalysis object.
Maintains consistency across the underlying SummarizedExperiment and all
computed results (diversity, LM, jackknife, divergence).

## Usage

``` r
# S4 method for class 'TSENATAnalysis,ANY,ANY,ANY'
x[i, j, drop = TRUE]
```

## Arguments

- x:

  TSENATAnalysis object

- i:

  Integer, logical, or character vector of rows (genes/transcripts) to
  retain. If missing, all rows are retained.

- j:

  Integer, logical, or character vector of columns (samples) to retain.
  If missing, all columns are retained.

- drop:

  Logical. Currently ignored (included for S4 compatibility). Always
  returns TSENATAnalysis (never drops to SE or vector).

## Value

A new TSENATAnalysis object containing only the specified genes and
samples

## Details

Subsetting preserves all analysis metadata and results while maintaining
consistency: - The underlying SummarizedExperiment is subset to the
specified genes/samples - Diversity and jackknife results are subset to
match sample selection - LM results are recalculated or removed if
sample structure changes - Divergence results are subset accordingly -
Analysis configuration is preserved

## Examples

``` r
# Create a minimal TSENATAnalysis object
library(SummarizedExperiment)
counts <- matrix(rpois(200, 10), nrow = 20, ncol = 10)
rownames(counts) <- paste0('TX_', 1:20)
colnames(counts) <- paste0('S', 1:10)
se <- SummarizedExperiment(
  assays = list(counts = counts),
  rowData = data.frame(gene_id = rep(paste0('G', 1:2), each = 10),
                       row.names = rownames(counts)),
  colData = data.frame(sample_id = colnames(counts),
                       condition = rep(c('A', 'B'), 5),
                       row.names = colnames(counts))
)
analysis <- TSENATAnalysis(se)

# Subset to first 10 genes and first 5 samples
analysis_subset <- analysis[1:10, 1:5]

# Subset by gene name
analysis_subset2 <- analysis[paste0('TX_', 1:5), ]

# Subset by sample condition (logical indexing)
keep_samples <- colData(se(analysis))$condition == 'A'
analysis_a <- analysis[, keep_samples]
```
