# Calculate Tsallis diversity values for transcripts grouped by gene

This helper computes per-gene Tsallis entropy across samples. The
trimmed package only supports the Tsallis method; other diversity
metrics were removed.

## Usage

``` r
calculate_method(
  x,
  genes,
  method = "tsallis",
  norm = TRUE,
  verbose = FALSE,
  q = 2,
  what = c("S", "D", "both")
)
```

## Arguments

- x:

  Numeric matrix or data.frame of transcript-level expression values
  (rows = transcripts, columns = samples).

- genes:

  Character vector with length equal to nrow(x) assigning each
  transcript to a gene.

- method:

  Only "tsallis" is supported (default).

- norm:

  Logical; if TRUE normalize Tsallis entropy values per gene.

- verbose:

  Logical; show diagnostic messages when TRUE.

- q:

  Numeric scalar or vector of q values to evaluate.

- what:

  Which quantity to return from \`calculate_tsallis_entropy\`: "S", "D"
  or "both" (default: "S").

## Value

A data.frame with genes in the first column and per-sample (and per-q)
Tsallis entropy values in subsequent columns.
