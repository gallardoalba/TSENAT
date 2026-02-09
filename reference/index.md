# Package index

## All functions

- [`build_se()`](https://gallardoalba.github.io/TSENAT/reference/build_se.md)
  : Build a SummarizedExperiment from transcript readcounts and
  tx-\>gene map
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
- [`filter_se()`](https://gallardoalba.github.io/TSENAT/reference/filter_se.md)
  : Filter transcripts in a \`SummarizedExperiment\` by minimum
  count/sample
- [`label_shuffling()`](https://gallardoalba.github.io/TSENAT/reference/label_shuffling.md)
  : Calculate p-values using label shuffling.
- [`map_metadata()`](https://gallardoalba.github.io/TSENAT/reference/map_metadata.md)
  : Map external coldata into a SummarizedExperiment
- [`map_tx_to_readcounts()`](https://gallardoalba.github.io/TSENAT/reference/map_tx_to_readcounts.md)
  : Map transcript IDs from a tx2gene table to a readcounts matrix
- [`plot_diversity_density()`](https://gallardoalba.github.io/TSENAT/reference/plot_diversity_density.md)
  : Plot diversity distributions (density) by sample type
- [`plot_ma_expression()`](https://gallardoalba.github.io/TSENAT/reference/plot_ma_expression.md)
  : Plot MA using expression/readcount-based fold changes
- [`plot_ma_tsallis()`](https://gallardoalba.github.io/TSENAT/reference/plot_ma_tsallis.md)
  : Plot MA using Tsallis-based fold changes
- [`plot_mean_violin()`](https://gallardoalba.github.io/TSENAT/reference/plot_mean_violin.md)
  : Plot violin of per-gene mean diversity by sample type
- [`plot_top_transcripts()`](https://gallardoalba.github.io/TSENAT/reference/plot_top_transcripts.md)
  : Plot top transcripts for a gene
- [`plot_tsallis_density_multq()`](https://gallardoalba.github.io/TSENAT/reference/plot_tsallis_density_multq.md)
  : Density plot of Tsallis entropy for multiple q values
- [`plot_tsallis_gene_profile()`](https://gallardoalba.github.io/TSENAT/reference/plot_tsallis_gene_profile.md)
  : Plot q-curve profile for a single gene comparing groups
- [`plot_tsallis_q_curve()`](https://gallardoalba.github.io/TSENAT/reference/plot_tsallis_q_curve.md)
  : Plot median +- IQR of Tsallis entropy across q values by group
- [`plot_tsallis_violin_multq()`](https://gallardoalba.github.io/TSENAT/reference/plot_tsallis_violin_multq.md)
  : Violin plot of Tsallis entropy for multiple q values
- [`plot_volcano()`](https://gallardoalba.github.io/TSENAT/reference/plot_volcano.md)
  : Volcano plot for differential results
- [`tcga_brca_luma_dataset`](https://gallardoalba.github.io/TSENAT/reference/tcga_brca_luma_dataset.md)
  : TCGA Luminal A breast cancer dataset
- [`test_differential()`](https://gallardoalba.github.io/TSENAT/reference/test_differential.md)
  : Run a differential test by name: Wilcoxon or label-shuffle
- [`wilcoxon()`](https://gallardoalba.github.io/TSENAT/reference/wilcoxon.md)
  : Calculate p-values using Wilcoxon rank sum test.
