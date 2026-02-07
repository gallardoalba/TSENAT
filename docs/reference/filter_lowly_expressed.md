# Filter Lowly-Expressed Transcripts

Removes transcripts with low expression levels based on a count
threshold and minimum number of samples meeting that threshold.

Removes transcripts with low expression levels based on a count
threshold and minimum number of samples meeting that threshold.

## Usage

``` r
filter_lowly_expressed(
  counts,
  genes = NULL,
  count_threshold = 5,
  min_samples = 5,
  verbose = TRUE
)

filter_lowly_expressed(
  counts,
  genes = NULL,
  count_threshold = 5,
  min_samples = 5,
  verbose = TRUE
)
```

## Arguments

- counts:

  A matrix or data.frame with transcripts as rows and samples as
  columns. Typically raw read counts.

- genes:

  Optional character vector of gene names corresponding to rows of
  `counts`. If provided, will be filtered to match the kept transcripts.

- count_threshold:

  Minimum count value for a transcript in a sample to be considered
  expressed (default: 5).

- min_samples:

  Minimum number of samples in which a transcript must exceed
  `count_threshold` to be retained (default: 5).

- verbose:

  Logical; if TRUE, print summary statistics before and after filtering
  (default: TRUE).

## Value

A list with elements:

- counts:

  Filtered count matrix.

- genes:

  Filtered gene names (if provided), or NULL.

- n_before:

  Number of transcripts before filtering.

- n_after:

  Number of transcripts after filtering.

- keep_idx:

  Logical vector indicating which transcripts were retained.

A list with elements:

- counts:

  Filtered count matrix.

- genes:

  Filtered gene names (if provided), or NULL.

- n_before:

  Number of transcripts before filtering.

- n_after:

  Number of transcripts after filtering.

- keep_idx:

  Logical vector indicating which transcripts were retained.

## Details

This function is commonly used as a preprocessing step to reduce noise
and improve the stability of downstream analyses like diversity
estimation. Transcripts must have counts exceeding `count_threshold` in
at least `min_samples` samples to be retained.

This function is commonly used as a preprocessing step to reduce noise
and improve the stability of downstream analyses like diversity
estimation. Transcripts must have counts exceeding `count_threshold` in
at least `min_samples` samples to be retained.

## Examples

``` r
# Create example count matrix
counts <- matrix(c(0, 1, 10, 15, 2, 3, 20, 25, 1, 0, 5, 8), nrow = 3)
genes <- c("gene_A", "gene_B", "gene_C")

# Filter transcripts
filtered <- filter_lowly_expressed(counts, genes)
#> Transcripts: before = 3, after = 0
filtered$counts
#>      [,1] [,2] [,3] [,4]
filtered$genes
#> character(0)
message(sprintf("Kept %d of %d transcripts", filtered$n_after, filtered$n_before))
#> Kept 0 of 3 transcripts

# Create example count matrix
counts <- matrix(c(0, 1, 10, 15, 2, 3, 20, 25, 1, 0, 5, 8), nrow = 3)
genes <- c("gene_A", "gene_B", "gene_C")

# Filter transcripts
filtered <- filter_lowly_expressed(counts, genes)
#> Transcripts: before = 3, after = 0
filtered$counts
#>      [,1] [,2] [,3] [,4]
filtered$genes
#> character(0)
message(sprintf("Kept %d of %d transcripts", filtered$n_after, filtered$n_before))
#> Kept 0 of 3 transcripts
```
