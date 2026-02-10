# Calculate Tsallis diversity per gene across samples

Calculate Tsallis diversity per gene across samples

## Usage

``` r
calculate_diversity(
  x,
  genes = NULL,
  norm = TRUE,
  tpm = FALSE,
  assayno = 1,
  verbose = FALSE,
  q = 2,
  what = c("S", "D")
)
```

## Arguments

- x:

  A numeric matrix or data.frame of transcript-level expression values
  (rows = transcripts, columns = samples), or a SummarizedExperiment-
  like object.

- genes:

  Character vector assigning each transcript (row) to a gene. Must have
  length equal to nrow(x) or the number of transcripts in \`x\`.

- norm:

  Logical; if TRUE, normalize Tsallis entropy to \[0,1\] per gene.

- tpm:

  Logical. If TRUE and \`x\` is a tximport-style list, use the
  \`\$abundance\` matrix instead of \`\$counts\`.

- assayno:

  Integer assay index to use when \`x\` is a SummarizedExperiment.

- verbose:

  Logical; print diagnostic messages when TRUE.

- q:

  Numeric scalar or vector of Tsallis q values to evaluate (q \> 0). If
  length(q) \> 1, the result will contain separate columns per sample
  and q.

- what:

  Which quantity to return: 'S' for Tsallis entropy or 'D' for Hill
  numbers.

## Value

A
[SummarizedExperiment](https://rdrr.io/pkg/SummarizedExperiment/man/SummarizedExperiment-class.html)
with assay \`diversity\` containing per-gene diversity values.

## Examples

``` r
data("tcga_brca_luma_dataset", package = "TSENAT")
rc <- as.matrix(tcga_brca_luma_dataset[1:20, -1, drop = FALSE])
gs <- tcga_brca_luma_dataset$genes[1:20]
se <- calculate_diversity(rc, gs, q = 0.1, norm = TRUE)
SummarizedExperiment::assay(se)[1:3, 1:3]
#>         TCGA-A7-A0CH_N TCGA-A7-A0CH_T TCGA-A7-A0D9_N
#> MXRA8        0.8409985      0.8149525      0.7861211
#> C1orf86      0.0000000      0.0000000      0.0000000
#> PDPN         0.3486679      0.3395884      0.0000000
```
