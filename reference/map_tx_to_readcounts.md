# Map transcript IDs from a tx2gene table to a readcounts matrix

Assign transcript identifiers as row names of a transcript-level read
counts matrix so downstream plotting functions can identify transcripts.
The tx2gene mapping can be provided as a file path (TSV) or a data.frame
with a \`Transcript\` column.

## Usage

``` r
map_tx_to_readcounts(
  readcounts,
  tx2gene,
  tx_col = "Transcript",
  verbose = FALSE
)
```

## Arguments

- readcounts:

  A numeric matrix or data.frame of read counts (rows = transcripts).

- tx2gene:

  Either a path to a tab-delimited file or a data.frame with at least a
  \`Transcript\` column.

- tx_col:

  Name of the transcript ID column in \`tx2gene\` (default:
  'Transcript').

- verbose:

  Logical; print informative messages (default: FALSE).

## Value

The input \`readcounts\` with rownames set to the transcript IDs.

## Examples

``` r
readcounts <- matrix(1:4, nrow = 2, dimnames = list(NULL, c('s1', 's2')))
tx2gene <- data.frame(Transcript = c('tx1', 'tx2'), Gene = c('g1', 'g2'), stringsAsFactors = FALSE)
rc <- map_tx_to_readcounts(readcounts, tx2gene)
rownames(rc)
#> [1] "tx1" "tx2"
```
