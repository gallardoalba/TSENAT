# Plot method concordance results from TSENATAnalysis

Plot method concordance results from TSENATAnalysis

## Usage

``` r
plot_method_concordance_s4(analysis, verbose = FALSE)

# S4 method for class 'TSENATAnalysis'
plot_method_concordance_s4(analysis, verbose = FALSE)
```

## Arguments

- analysis:

  `TSENATAnalysis` object with computed method concordance (from
  [`compute_method_concordance_s4()`](https://gallardoalba.github.io/TSENAT/reference/compute_method_concordance_s4.md)).

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
[`compute_method_concordance_s4()`](https://gallardoalba.github.io/TSENAT/reference/compute_method_concordance_s4.md)
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
config <- tsenat_config(sample_col = 'sample', condition_col = 'condition')
analysis <- build_analysis_s4(readcounts = readcounts, tx2gene =
gff3_dataset, metadata = metadata_df, config = config,
  tpm = tpm, effective_length = effective_length)
analysis <- filter_analysis_s4(analysis, min_samples = 1, subset_n_genes
= 200)

# Note: compute_method_concordance_s4 requires additional LM and Friedman
# results computed. For demo purposes, we show that
# plot_method_concordance_s4 needs pre-computed concordance in @metadata
```
