# Plot Tsallis Entropy q-Curve

Visualize Tsallis entropy (S_q) as a function of the diversity parameter
q across sample groups. Supports three modes: aggregate q-curves
(default), gene-specific q-curves (when \`gene\` provided), or bootstrap
confidence interval bands.

## Usage

``` r
plot_diversity_spectrum(
  se,
  assay_name = "diversity",
  condition_col = NULL,
  gene = NULL,
  sait_res = NULL,
  n_top = NULL,
  metric = "iqr",
  output_file = NULL,
  dev_width = NULL,
  dev_height = NULL
)
```

## Arguments

- se:

  A \`SummarizedExperiment\` returned by \`.calculate_diversity()\` with
  diversity assay. If pre-computed bootstrap confidence intervals are
  available (ci_lower/ci_upper assays), they will be displayed
  automatically. Otherwise, falls back to IQR visualization.

- assay_name:

  Character; name of the assay to plot (default: 'diversity').

- condition_col:

  Character or NULL; column name in colData indicating group/sample
  type. If NULL (default), reads from \`@config\$condition_col\` when
  input is TSENATAnalysis, otherwise defaults to 'sample_type'. Only
  used in aggregate and CI modes.

- gene:

  Character vector (optional); if provided, plot q-curves for specified
  gene(s). Overrides default aggregate behavior. When provided, uses
  median +/- SD for each gene.

- sait_res:

  Data frame (optional); gene interaction test results with \`gene\`
  column and p-value column. Accepts either: - Results from
  \`.calculate_sait()\` (has \`adj_p_interaction\` or \`p_interaction\`
  columns) - Results from \`.calculate_rank_transform()\` (has
  \`adj_p_value\` or \`p_value\` columns from Conover-Iman Rank
  Transform tests) If provided (and \`gene\` is NULL), plots top
  \`n_top\` genes ranked by p-value. Useful for plotting significant
  genes from any interaction analysis.

- n_top:

  Integer or NULL; number of top genes to select from \`sait_res\` when
  \`gene\` is NULL (default: NULL). When NULL, defaults to showing the
  single most significant gene (n_top=1), providing a conservative view
  of the strongest effect. Set to a numeric value to show that many top
  genes.

- metric:

  Character; when bootstrap CIs are NOT available, specifies the spread
  metric to display. Options: 'iqr' (default, Interquartile Range - more
  robust) or 'sd' (Standard Deviation). This parameter is ignored when
  bootstrap confidence intervals are available. Default: 'iqr'.

- output_file:

  `character` or `NULL`. Optional file path to save the plot. Default:
  NULL (no file output).

- dev_width:

  Numeric or NULL; width in inches for the graphics device when
  displaying the plot interactively. If specified (along with
  dev_height), creates a new device with this width. Default: NULL (uses
  current device).

- dev_height:

  Numeric or NULL; height in inches for the graphics device when
  displaying the plot interactively. If specified (along with
  dev_width), creates a new device with this height. Default: NULL (uses
  current device).

## Value

\*\*Aggregate mode (gene=NULL, sait_res=NULL)\*\*: - If bootstrap CI
assays available: A ggplot object showing median entropy with bootstrap
confidence interval bands. - If no CI assays: A ggplot object showing
median entropy with IQR ribbons (automatic fallback).

\*\*Gene-specific mode (gene or sait_res provided)\*\*: - Single gene: A
ggplot object showing median entropy +/- SD for that gene. - Multiple
genes: A grid plot object arranged in 2 rows x 2 columns with a shared
legend at the bottom. The legend appears once beneath the grid, avoiding
repetition across subplots.

## Details

\*\*Aggregate mode (default, gene=NULL, sait_res=NULL)\*\*: - Plots
median Tsallis entropy +/- IQR across all genes for each group - Works
with any SummarizedExperiment from .calculate_diversity() - Supports
single or multiple q values and any number of groups - No CI data
required for basic plots; bootstrap CIs optional

\*\*Gene-specific mode (gene or sait_res provided)\*\*: - Plots q-curve
separately for each selected gene - Shows median entropy +/- SD
(variance) for each gene across q-values and groups - When \`sait_res\`
provided: automatically ranks genes and selects top \`n_top\` by
p-value - Single gene: returns a ggplot object; multiple genes: returns
a grid plot (2 rows x 2 columns) with shared legend - For multiple
genes: legend appears once at the bottom of the grid to avoid repetition
and save space - Useful for highlighting specific genes of interest or
significant discoveries - Bootstrap mode not supported in gene-specific
mode

\*\*Automatic Bootstrap CI Detection\*\*: - When computing divergence or
diversity with bootstrap enabled, ci_lower and ci_upper assays are added
to the SummarizedExperiment. - This function automatically detects these
assays and displays bootstrap confidence interval bands instead of IQR.
No additional parameter needed. - For confidence bands to appear, use
\`calculate_diversity(..., bootstrap=TRUE, nboot=1000)\` or appropriate
divergence function with bootstrap enabled.

## Examples

``` r
# Plot 7: Tsallis entropy q-curve (combined across all sample diversity)
data(readcounts)
metadata_df <- read.table(
  system.file('extdata', 'metadata.tsv', package = 'TSENAT'),
  header = TRUE, sep = '\t'
)
gff3_dataset <- system.file('extdata', 'annotation.gff3.gz', package =
'TSENAT')
readcounts <- as.matrix(readcounts)
mode(readcounts) <- 'numeric'

# Create configuration (required when metadata is provided)
config <- TSENAT_config(sample_col = 'sample', condition_col = 'condition')
analysis <- build_analysis(readcounts = readcounts, tx2gene =
gff3_dataset, metadata = metadata_df, config = config,
  tpm = tpm, effective_length = effective_length)
analysis <- filter_analysis(analysis, min_samples = 1, subset_n_genes
= 200)
analysis <- calculate_diversity(analysis, q = seq(0, 2, by = 0.5),
)
p <- plot_diversity_spectrum(analysis)
# print(p)
```
