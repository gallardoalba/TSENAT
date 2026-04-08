# Plot Multi-Q Delta Influence Heatmaps from TSENATAnalysis Object

S4 wrapper for `. plot_multiq_delta_influence_heatmaps()` that extracts
results directly from a TSENATAnalysis object. Automatically retrieves
jackknife switching results from the analysis object slots.

## Usage

``` r
plot_multiq_delta_influence_heatmaps_s4(
  analysis,
  n_genes = 4,
  lm_results = NULL,
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

- lm_results:

  `data.frame` or `NULL`. Optional LM interaction results for ranking
  genes (default: NULL). If NULL, attempts to extract from
  `analysis@lm_results$lm_interaction`.

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

- LM results:

  From `analysis@lm_results$lm_interaction` if not explicitly provided,
  for ranking genes by significance

The wrapper automatically handles parameter extraction and provides a
simplified interface compared to the base function.

## See also

[`jackknife_isoform_switching_s4`](https://gallardoalba.github.io/TSENAT/reference/jackknife_isoform_switching_s4.md)
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
config <- tsenat_config(
  sample_col = 'sample',
  condition_col = 'condition',
  subject_col = 'paired_samples',
  paired = TRUE,
  q_values = seq(0, 2, by = 0.1)
)
#> Warning: [tsenat_config] Paired design (paired=TRUE) requires complete configuration.
#>   Missing or incomplete parameters: control
#>   This will cause downstream analysis failure or empty results (LM interaction, plotting).
#>   Provide all parameters: 
#>     config <- tsenat_config(
#>       q_values = seq(0, 2, by = 0.05),       # At least 5 q-values (41 recommended)
#>       condition_col = 'condition',
#>       subject_col = 'paired_samples',        # Required for paired analysis
#>       paired = TRUE,
#>       control = 'normal'                     # Reference group for comparisons
#>     )

# Build analysis with configured parameters
analysis <- build_analysis_s4(
  readcounts = readcounts,
  tx2gene = gff3_dataset,
  metadata = metadata_df,
  config = config,
  tpm = tpm,
  effective_length = effective_length
)

analysis <- filter_analysis_s4(analysis, stringency = 'severe')
analysis <- calculate_diversity_s4(analysis, q = c(0.5, 1.0, 1.5, 2.0, 2.5), verbose
= FALSE)
analysis <- calculate_divergence_s4(analysis, q = c(0.5, 1.0, 1.5, 2.0, 2.5))
analysis <- calculate_lm_interaction_s4(analysis, method = 'gam')
#> Warning: nlminb problem, convergence error code = 1
#>   message = singular convergence (7)
#> Warning: nlminb problem, convergence error code = 1
#>   message = iteration limit reached without convergence (10)
#> Warning: nlminb problem, convergence error code = 1
#>   message = iteration limit reached without convergence (10)
analysis <- jackknife_isoform_switching_s4(analysis, q = c(0.5, 1, 1.5),
  n_bootstrap = 50)
heatmap_file <- plot_multiq_delta_influence_heatmaps_s4(analysis, n_genes
= 2)

```
