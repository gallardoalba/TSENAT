# Package index

## All functions

- [`calculate_difference()`](https://gallardoalba.github.io/TSENAT/reference/calculate_difference.md)
  : Calculate splicing diversity changes between two conditions.
- [`calculate_diversity()`](https://gallardoalba.github.io/TSENAT/reference/calculate_diversity.md)
  : Calculate Tsallis diversity per gene across samples
- [`calculate_fc()`](https://gallardoalba.github.io/TSENAT/reference/calculate_fc.md)
  : Calculate splicing diversity changes between two conditions.
- [`calculate_lm_interaction()`](https://gallardoalba.github.io/TSENAT/reference/calculate_lm_interaction.md)
  : Linear-model interaction test for Tsallis entropy For each gene, fit
  a linear model of the form \`entropy ~ q \* group\` and extract the
  p-value for the interaction term (whether the effect of \`q\` differs
  between groups). The function expects a \`SummarizedExperiment\`
  produced by \`calculate_diversity()\` when multiple \`q\` values have
  been computed (column names contain \`\_q=\`).
- [`calculate_method()`](https://gallardoalba.github.io/TSENAT/reference/calculate_method.md)
  : Calculate Tsallis diversity values for transcripts grouped by gene
- [`calculate_tsallis_entropy()`](https://gallardoalba.github.io/TSENAT/reference/calculate_tsallis_entropy.md)
  : Calculate Tsallis entropy for a vector of transcript-level
  expression values of one gene.
- [`infer_sample_group()`](https://gallardoalba.github.io/TSENAT/reference/infer_sample_group.md)
  : Infer sample group from sample names
- [`label_shuffling()`](https://gallardoalba.github.io/TSENAT/reference/label_shuffling.md)
  : Calculate p-values using label shuffling.
- [`map_coldata_to_se()`](https://gallardoalba.github.io/TSENAT/reference/map_coldata_to_se.md)
  : Map external coldata into a SummarizedExperiment This helper maps an
  external \`coldata\` table (with sample IDs and a condition/label
  column) into a \`SummarizedExperiment\` produced by
  \`calculate_diversity()\`. It assigns a \`sample_type\` column to
  \`colData(ts_se)\` and falls back to \`infer_sample_group()\` for
  unmapped samples.
- [`plot_diversity_density()`](https://gallardoalba.github.io/TSENAT/reference/plot_diversity_density.md)
  : Plot diversity distributions (density) by sample type
- [`plot_ma()`](https://gallardoalba.github.io/TSENAT/reference/plot_ma.md)
  : Plot MA plot for difference results
- [`plot_mean_violin()`](https://gallardoalba.github.io/TSENAT/reference/plot_mean_violin.md)
  : Plot violin of per-gene mean diversity by sample type
- [`plot_top_transcripts`](https://gallardoalba.github.io/TSENAT/reference/plot_top_transcripts.md)
  : Plot top transcripts for a gene \#' For a given gene, find
  transcripts using a tx-\>gene mapping, compute per- Plot top
  transcripts for a gene For a given gene, find transcripts using a
  tx-\>gene mapping, compute per- transcript statistics between two
  sample groups, select the top N transcripts by p-value and plot their
  expression across groups.
- [`plot_tsallis_density_multq()`](https://gallardoalba.github.io/TSENAT/reference/plot_tsallis_density_multq.md)
  : Density plot of Tsallis entropy for multiple q values
- [`plot_tsallis_q_curve()`](https://gallardoalba.github.io/TSENAT/reference/plot_tsallis_q_curve.md)
  : Plot median +- IQR of Tsallis entropy across q values by group This
  reproduces the \`tsallis-q-curve-mean-sd\` plot from the vignette: for
  each q value, compute per-gene Tsallis entropy per sample, then
  summarize across genes by group (median and IQR) and plot median with
  a ribbon spanning median +- IQR/2.
- [`plot_tsallis_violin_multq()`](https://gallardoalba.github.io/TSENAT/reference/plot_tsallis_violin_multq.md)
  : Violin plot of Tsallis entropy for multiple q values
- [`plot_volcano()`](https://gallardoalba.github.io/TSENAT/reference/plot_volcano.md)
  : Volcano plot for differential results
- [`tcga_brca_luma_dataset`](https://gallardoalba.github.io/TSENAT/reference/tcga_brca_luma_dataset.md)
  : TCGA Luminal A breast cancer dataset
- [`wilcoxon()`](https://gallardoalba.github.io/TSENAT/reference/wilcoxon.md)
  : Calculate p-values using Wilcoxon rank sum test.
