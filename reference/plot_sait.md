# Plot GAM q-curves from TSENATAnalysis object

S4 wrapper that accepts a TSENATAnalysis object and generates GAM
q-curve plots

## Usage

``` r
plot_sait(
  analysis,
  n_top = 6,
  genes = NULL,
  condition_col = NULL,
  sig_alpha = 0.05,
  assay_name = "diversity",
  output_file = NULL,
  width = 12,
  height = NULL,
  verbose = FALSE,
  ...
)
```

## Arguments

- analysis:

  `TSENATAnalysis` object with diversity and SAIT interaction results.

- n_top:

  `integer`. Number of top genes (by adjusted p-value) to plot (default:
  6). Only used if genes = NULL.

- genes:

  `character` vector. Optional specific gene names to plot. If provided,
  these genes are plotted directly regardless of significance. If NULL
  (default), top n_top significant genes are selected.

- condition_col:

  `character`. Column name in colData(se) specifying group assignments
  for samples. If NULL, attempts to auto-detect from `@config` (looks
  for 'condition_col' or 'condition'). If still NULL, defaults to
  'sample_type'.

- sig_alpha:

  `numeric`. Significance threshold for adjusted p-values (default:
  0.05). Only used if genes = NULL; filters sait_res to significant
  genes before selecting top n.

- assay_name:

  `character`. Name of the assay in se to extract (default:
  'diversity').

- output_file:

  `character` or `NULL`. Optional file path to save the plot. Default:
  NULL (no file output).

- width:

  `numeric`. Width of the plot in inches (default: 12). Only used if
  output_file is not NULL.

- height:

  `numeric` or `NULL`. Height of the plot in inches. Default: NULL
  (automatically calculated based on width and aspect ratio). Only used
  if output_file is not NULL.

- verbose:

  `logical`. If TRUE, print diagnostic messages during processing.
  (default: FALSE).

- ...:

  Additional arguments passed to the base function.

## Value

A single `ggplot` object with all selected genes arranged in a grid
layout. Can be saved with
[`ggplot2::ggsave()`](https://ggplot2.tidyverse.org/reference/ggsave.html).

## Details

This wrapper automatically: 1. Extracts SummarizedExperiment from `@se`
slot 2. Extracts SAIT results from `@sait_results$sait_interaction` slot
3. Detects condition_col from `@config` or uses default 4. Calls
`.plot_sait()` with extracted parameters

\*\*Parameter Resolution (condition_col):\*\*

1.  Explicit `condition_col` parameter (highest priority)

2.  `@config$condition_col` if available

3.  `@config$condition` if available

4.  Default: 'sample_type'

## See also

[`calculate_sait`](https://gallardoalba.github.io/TSENAT/reference/calculate_sait.md)
for running LM analysis on TSENATAnalysis.

## Examples

``` r
if (FALSE) { # \dontrun{
# Plot 3: GAM q-curves for genes with q-by-condition interactions
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
  q = seq(0.2, 2, by = 0.4)  # 5 unique q-values: 0.2, 0.6, 1.0, 1.4, 1.8
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
analysis <- calculate_diversity(analysis, q = seq(0.2, 2, by = 0.4))
analysis <- suppressWarnings(calculate_sait(analysis, method = 'gam'))

p_gam <- plot_sait(analysis, n_top = 2, sig_alpha = 0.15)
# print(p_gam)
} # }
```
