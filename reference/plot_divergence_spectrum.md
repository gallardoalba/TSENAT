# Plot Global Divergence q-Curve Across All Genes (S4 Wrapper)

S4 wrapper that extracts divergence results from a TSENATAnalysis object
and visualizes the average Tsallis divergence across all genes (or
specified genes) as a function of q-value.

## Usage

``` r
plot_divergence_spectrum(
  analysis,
  gene = NULL,
  n_genes = 4,
  ncol = 2,
  metric = c("median", "mean"),
  variability_metric = c("iqr", "sd"),
  use_pvalue_ranking = FALSE,
  output_file = NULL,
  width = 12,
  height = NULL,
  verbose = FALSE,
  ...
)
```

## Arguments

- analysis:

  `TSENATAnalysis` object with divergence results (typically via
  [`calculate_divergence`](https://gallardoalba.github.io/TSENAT/reference/calculate_divergence.md)).

- gene:

  `character`. Optional specific gene name to plot. If NULL, plots
  global divergence curve (aggregated across all genes).

- n_genes:

  `integer`. Number of top genes to plot when showing multi-gene
  spectra. Default is 4. Genes are sorted by p-value significance.

- ncol:

  `integer`. Number of columns in grid layout for multi-gene plots.
  Default is 2. Number of rows is automatically calculated.

- metric:

  `character`. Summary statistic for global curve: 'median' (default) or
  'mean'. Only used when gene = NULL.

- variability_metric:

  `character`. Error bar type for global curve: 'iqr' (default) or 'sd'.
  Only used when gene = NULL.

- use_pvalue_ranking:

  `logical`. If TRUE, uses LM results to rank and display top n_genes by
  p-value significance. If FALSE (default), plots global divergence
  curve when gene = NULL. Default is FALSE.

- output_file:

  `character`. Optional file path to save the plot. If NULL, plot is
  returned but not saved.

- width:

  `numeric`. Plot width in inches. Default is 10.

- height:

  `numeric`. Plot height in inches. Default is 6.

- verbose:

  `logical`. Print status messages. Default is TRUE.

- ...:

  Additional arguments passed to the underlying plotting function.

## Value

Invisibly returns the file path if saved, otherwise the ggplot object.
If the plot cannot be created (missing data, ggplot2 not available),
returns NULL invisibly with an informative message.

## Details

This wrapper extracts the divergence SummarizedExperiment from
`analysis@divergence_results` and optionally the LM results from
`analysis@lm_results$lm_interaction` to pass to the base function.

\*\*Data Requirements:\*\*

- Divergence must be computed via
  [`calculate_divergence()`](https://gallardoalba.github.io/TSENAT/reference/calculate_divergence.md)

- `@divergence_results$divergence_se` or direct divergence SE

\*\*Modes:\*\*

- **Global mode** (gene = NULL): Shows median/mean divergence across all
  genes with variability bands

- **Gene-specific mode** (gene specified): Shows divergence spectrum for
  a single named gene

- **Top genes mode** (gene = NULL, lm_res provided): Shows top n_genes
  by significance

## See also

[`calculate_divergence`](https://gallardoalba.github.io/TSENAT/reference/calculate_divergence.md)
for computing divergence.

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
analysis <- filter_analysis(
  analysis,
  min_samples = 1,
  subset_n_genes = 200
)
analysis <- calculate_diversity(
  analysis,
  q = c(0.5, 1, 1.5),
  verbose = FALSE
)
analysis <- calculate_divergence(
  analysis,
  q = c(0.5, 1, 1.5),
  verbose = FALSE
)
p_global <- plot_divergence_spectrum(analysis)
# print(p_global)
```
