# Build a SummarizedExperiment from transcript readcounts and tx-\>gene map

This helper creates a \`SummarizedExperiment\` with an assay named
\`counts\`, stores the \`tx2gene\` table and the raw \`readcounts\` in
\`metadata()\`, and places the provided \`genes\` vector into
\`rowData(se)\$genes\`.

## Usage

``` r
build_se(tx2gene_tsv, readcounts, genes, assay_name = "counts")
```

## Arguments

- tx2gene_tsv:

  Path to a tab-separated \`tx2gene\` file or a data.frame with
  transcript-to-gene mapping. If a path is provided, it must contain at
  least two columns (transcript, gene).

- readcounts:

  Numeric matrix or data.frame of transcript-level counts (rows =
  transcripts, columns = samples).

- genes:

  Character vector of gene IDs assigning each transcript to a gene;
  length must equal nrow(readcounts).

- assay_name:

  Name for the assay to store readcounts (default: 'counts').

## Value

A \`SummarizedExperiment\` with assay, \`metadata()\$tx2gene\`,
\`metadata()\$readcounts\` and \`rowData(se)\$genes\` populated.

## Examples

``` r
tx2gene <- data.frame(Transcript = c('tx1', 'tx2'), Gene = c('g1', 'g2'), stringsAsFactors = FALSE)
readcounts <- matrix(c(10, 5, 2, 3), nrow = 2, dimnames = list(c('tx1', 'tx2'), c('s1', 's2')))
genes <- c('g1', 'g2')
se <- build_se(tx2gene, readcounts, genes)
SummarizedExperiment::assay(se, 'counts')
#>     s1 s2
#> tx1 10  2
#> tx2  5  3
```
