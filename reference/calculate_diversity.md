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
  (rows = transcripts, columns = samples), or a
  SummarizedExperiment-like object.

- genes:

  Character vector assigning each transcript (row) to a gene. Must have
  length equal to nrow(x) or the number of transcripts in \`x\`.

- norm:

  Logical; if TRUE, normalize Tsallis entropy to \[0,1\] per gene.

- tpm:

  Logical; if TRUE and \`x\` is a tximport-style list, use the
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

  Which quantity to return: "S" for Tsallis entropy, "D" for Hill
  numbers, or "both".

## Value

A
[SummarizedExperiment](https://rdrr.io/pkg/SummarizedExperiment/man/SummarizedExperiment-class.html)
with assay \`diversity\` containing per-gene diversity values.
