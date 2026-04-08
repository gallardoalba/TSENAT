# Plot Q-Spectrum Curves for Multiple Top Genes

Wrapper around .plot_tsallis_q_curve() that creates multi-gene
q-spectrum visualization from Tsallis divergence results. Shows how
divergence changes across q-values for top genes, revealing q-dependent
isoform switching.

## Usage

``` r
plot_multi_gene_q_spectrum_s4(
  eff_res = NULL,
  lm_res = NULL,
  divergence_results_se = NULL,
  n_genes = 9,
  ncol = 3,
  verbose = FALSE,
  output_file = NULL
)
```

## Arguments

- eff_res:

  Output from effect size computation OR `NULL`. If provided, must
  contain \`\$interaction_results\` with columns: gene,
  adj_p_interaction, per_q_pattern. If `NULL`, uses fallback with
  lm_res + divergence_results_se.

- lm_res:

  (Optional) Data frame from LMM analysis with columns: gene,
  adj_p_interaction. Only used if eff_res is NULL. Must be provided for
  fallback mode.

- divergence_results_se:

  (Optional) SummarizedExperiment from divergence calculation. Only used
  if eff_res is NULL. Must be provided for fallback mode.

- n_genes:

  Integer; number of top genes to plot (default: 9). Genes are sorted by
  increasing adjusted p-value (most significant first).

- ncol:

  Integer; number of columns in grid layout (default: 3). Number of rows
  is automatically calculated as ceiling(n_genes / ncol).

- verbose:

  Logical; if TRUE, print diagnostic messages (default: TRUE).

- output_file:

  Character or NULL. Optional file path to save the plot as an image. If
  provided, the plot will be saved with appropriate dimensions. Default:
  NULL (no file output, only return object).

## Value

A ggplot2 object created via `patchwork` combining all gene panels, or
NULL if gene data is unavailable. The function automatically handles
ggplot2 grid creation and returns a print-ready object.

## Details

Key Features:

- Multi-gene grid: Side-by-side comparison of top N genes (default: 9)

- q-spectrum curves: Per-gene divergence across q = 0.01 to 2.00

- Bootstrap CI bands: 95

- Statistical significance: Gene titles include adjusted p-values

- Automatic layout: Grid dimensions auto-calculated from n_genes

- q=1 reference line: KL divergence (information theory benchmark)

- Region labels: Rare (q\<1) vs Abundant (q\>1) isoform emphasis

Creates a multi-panel grid comparing per-q divergence profiles across
the top N genes identified by LMM interaction analysis. Each panel shows
the full q-spectrum divergence curve with the gene name and adjusted
p-value in the title.

\*\*Mathematical Background:\*\* Tsallis divergence as function of q:

- q \< 1: Emphasizes rare isoforms (sensitive to outliers)

- q = 1: Kullback-Leibler divergence (classical information theory)

- q \> 1: Emphasizes abundant isoforms (robust to rare variants)

D_q varies across q-spectrum, showing complexity of isoform differences.

\*\*Example Interpretation:\*\*

- Flat curve: Divergence stable across q (robust isoform difference)

- Curved pattern: q-dependent divergence (rare vs abundant isoforms
  differ)

- Peaks at high q: Main isoforms drive the difference, rare ones
  immaterial

\*\*Input Modes:\*\* - \*\*Mode 1 (Primary)\*\*: Pass eff_res directly
(from effect_sizes_divergence) - \*\*Mode 2 (Fallback)\*\*: Pass
lm_res + divergence_results_se instead

\*\*Gene Filtering:\*\* Genes are ranked by decreasing statistical
significance (increasing adj_p_interaction). Only genes with complete
per-q divergence data are included. If fewer than n_genes have valid
data, the function returns all available genes.

\*\*Plot Features:\*\* - Title shows: gene name and adjusted p-value
(q-value format) - Per-q divergence curve with point estimates and 95 -
Vertical reference line at q=1 (Kullback-Leibler divergence point) -
Region labels: 'Rare Isoforms' (q\<1), 'Balanced' (q~=1), 'Abundant
Isoforms' (q\>1) - All plots use consistent ggplot2 styling matching
plot_q_spectrum

## See also

[`calculate_divergence_s4`](https://gallardoalba.github.io/TSENAT/reference/calculate_divergence_s4.md)
for computing divergence values.

## Examples

``` r
# Plot 4: Multi-gene q-spectrum profiles
set.seed(42)
# Create robust test dataset with clear signal-to-noise ratio
n_genes <- 16
n_isoforms_per_gene <- 4
n_isoforms <- n_genes * n_isoforms_per_gene
n_samples_per_group <- 30  # Increased for statistical power
n_samples <- n_samples_per_group * 2

# Generate control and treatment with very strong separation
# This ensures sufficient statistical power for divergence tests
control_counts <- matrix(rpois(n_isoforms * n_samples_per_group, lambda =
100),
  nrow = n_isoforms, ncol = n_samples_per_group)
treatment_counts <- matrix(rpois(n_isoforms * n_samples_per_group, lambda
= 300),
  nrow = n_isoforms, ncol = n_samples_per_group)
counts <- cbind(control_counts, treatment_counts)
rownames(counts) <- paste0('TX_', 1:n_isoforms)
colnames(counts) <- paste0('Sample_', 1:n_samples)

se <- SummarizedExperiment::SummarizedExperiment(assays = list(counts =
counts))
tx2gene_df <- data.frame(Transcript = rownames(counts),
  Gene = rep(paste0('GENE_', 1:n_genes), each = n_isoforms_per_gene))
S4Vectors::metadata(se)$tx2gene <- tx2gene_df
SummarizedExperiment::colData(se) <- S4Vectors::DataFrame(
  sample_id = paste0('Sample_', 1:n_samples),
  condition = rep(c('control', 'treatment'), each = n_samples_per_group),
  row.names = colnames(se))
SummarizedExperiment::rowData(se)$transcript_id <- rownames(se)
SummarizedExperiment::rowData(se)$gene_id <-
tx2gene_df$Gene[match(rownames(se),
  tx2gene_df$Transcript)]

# Run complete analysis pipeline (skipped for speed in documentation)
# Uncomment to run actual analysis:
# analysis <- TSENATAnalysis(se)
# analysis <- calculate_diversity_s4(analysis, q = c(0.5, 1, 1.5), 
#   nboot = 50)
# analysis <- calculate_divergence_s4(analysis, q = c(0.5, 1, 1.5), 
#   nboot = 50)
# analysis <- calculate_lm_interaction_s4(analysis,
#   condition_col = 'condition')
# analysis <- effect_sizes_divergence_s4(analysis)
# p <- plot_multi_gene_q_spectrum_s4(analysis, n_genes = 4)
# if (!is.null(p)) print(p)
```
