# Plot method concordance results from TSENATAnalysis

Plot method concordance results from TSENATAnalysis

## Usage

``` r
plot_concordance(analysis, verbose = FALSE)

# S4 method for class 'TSENATAnalysis'
plot_concordance(analysis, verbose = FALSE)
```

## Arguments

- analysis:

  `TSENATAnalysis` object with computed method concordance (from
  [`calculate_concordance()`](https://gallardoalba.github.io/TSENAT/reference/calculate_concordance.md)).

- verbose:

  `logical`. Print progress messages. Default: FALSE

## Value

A ggplot/cowplot object showing:

- Panel 1:

  Scatter plot of -log10(p-values) comparing methods

- Panel 2:

  Histogram of p-value distributions by method

## Details

Creates visualization of method concordance including: - Comparison of
significance across two methods (with color-coded agreement) - P-value
distribution histograms for both methods - Significance threshold lines
at p \< 0.05

Requires that
[`calculate_concordance()`](https://gallardoalba.github.io/TSENAT/reference/calculate_concordance.md)
has already been run to populate `@metadata$method_concordance`.

## Examples

``` r
# Load example data (matching TSENAT.Rmd workflow)
data(readcounts)
readcounts <- as.matrix(readcounts)
mode(readcounts) <- 'numeric'
metadata_df <- read.table(
  system.file('extdata', 'metadata.tsv', package = 'TSENAT'),
  header = TRUE, sep = '\t'
)
gff3_dataset <- system.file('extdata', 'annotation.gff3.gz', package =
'TSENAT')

# Build analysis from vignette data and create small subset
config <- TSENAT_config(sample_col = 'sample', condition_col = 'condition')
analysis <- build_analysis(readcounts = readcounts, tx2gene =
gff3_dataset, metadata = metadata_df, config = config,
  tpm = tpm, effective_length = effective_length)
analysis <- filter_analysis(analysis, min_samples = 1, subset_n_genes
= 200)

# Note: calculate_concordance requires additional LM and Conover-Iman Rank Transform
# results computed. For demo purposes, we show that
# plot_concordance needs pre-computed concordance in @metadata
```
