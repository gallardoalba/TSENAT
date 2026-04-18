# Plot Multi-Q Delta Influence Heatmaps from TSENATAnalysis Object

S4 wrapper for `. plot_jis_delta()` that extracts results directly from
a TSENATAnalysis object. Automatically retrieves jackknife switching
results from the analysis object slots.

## Usage

``` r
plot_jis_delta(
  analysis,
  n_genes = 4,
  rrm_results = NULL,
  verbose = FALSE,
  output_file = NULL,
  ...
)
```

## Arguments

- analysis:

  `TSENATAnalysis`. An S4 object containing completed jackknife isoform
  switching analysis across multiple q-values.

- n_genes:

  `numeric`. Number of top genes to display in heatmaps (default: 4).
  Genes are ranked by LM p-values if available, otherwise by order of
  appearance in results.

- rrm_results:

  `data.frame` or `NULL`. Optional RRM interaction results for ranking
  genes (default: NULL). If NULL, attempts to extract from
  `analysis@rrm_results$rrm_interaction`.

- verbose:

  `logical`. If `TRUE`, print diagnostic messages during plot generation
  (default: FALSE).

- output_file:

  `character` or `NULL`. Optional file path to save the plot. Default:
  NULL (no file output).

- ...:

  Additional arguments passed to the base function.

## Value

A file path (character) to the saved heatmap PNG file, invisibly.

## Details

This function extracts the following from `analysis`:

- Jackknife results:

  From `analysis@jackknife_results`, which should contain multi-q
  switching results keyed by q-value (e.g., 'q_1.00')

- RRM results:

  From `analysis@rrm_results$rrm_interaction` if not explicitly
  provided, for ranking genes by significance

The wrapper automatically handles parameter extraction and provides a
simplified interface compared to the base function.

## See also

[`calculate_jis`](https://gallardoalba.github.io/TSENAT/reference/calculate_jis.md)
for computing switching results

## Examples

``` r
# Plot 5: Multi-q delta influence (isoform switching) heatmaps
data(readcounts)
readcounts <- as.matrix(readcounts)
mode(readcounts) <- 'numeric'
metadata_df <- read.table(
  system.file('extdata', 'metadata.tsv', package = 'TSENAT'),
  header = TRUE, sep = '\t'
)
gff3_dataset <- system.file('extdata', 'annotation.gff3.gz', package =
'TSENAT')

# Configure analysis parameters first
config <- TSENAT_config(
  sample_col = 'sample',
  condition_col = 'condition',
  subject_col = 'paired_samples',
  paired = TRUE,
  control = 'normal',
  q = seq(0, 2, by = 0.1)
)

# Build analysis with configured parameters
analysis <- build_analysis(
  readcounts = readcounts,
  tx2gene = gff3_dataset,
  metadata = metadata_df,
  config = config,
  tpm = tpm,
  effective_length = effective_length
)

analysis <- filter_analysis(analysis, stringency = 'severe')
analysis <- calculate_diversity(
  analysis,
  q = c(0.5, 1.0, 1.5, 2.0, 2.5),
  verbose = FALSE
)
analysis <- calculate_divergence(
  analysis,
  q = c(0.5, 1.0, 1.5, 2.0, 2.5)
)
analysis <- suppressWarnings(calculate_rrm(analysis, method = 'gam'))
analysis <- calculate_jis(
  analysis,
  q = c(0.5, 1, 1.5),
  n_bootstrap = 50
)
heatmap_file <- plot_jis_delta(analysis, n_genes
= 2)

```
