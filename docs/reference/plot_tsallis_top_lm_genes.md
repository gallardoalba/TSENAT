# Plot q-curves for top linear-model genes

Select the top genes from an lm interaction result (\`lm_res\`) by
adjusted p-value (robust to common column names). Compute median ± IQR
across samples per gene and plot the q-curve for each selected gene.
Returns a \`ggplot\` object.

## Usage

``` r
plot_tsallis_top_lm_genes(
  se,
  lm_res = NULL,
  lm_res2 = NULL,
  top_n = 6,
  assay_name = "diversity",
  sample_type_col = "sample_type"
)
```

## Arguments

- se:

  A \`SummarizedExperiment\` produced by \`calculate_diversity()\`.

- lm_res:

  Optional data.frame with linear-model results. If provided, the
  function will try to find a p-value column and select top genes by
  smallest adjusted p-value.

- top_n:

  Number of top genes to plot when \`lm_res\` is used (default: 6).

- assay_name:

  Assay name in \`se\` to use if \`hill\` is not available (default:
  "diversity").

- sample_type_col:

  Column name in \`colData(se)\` containing group labels (default:
  "sample_type").

## Value

A \`ggplot\` object plotting median ± IQR q-curves for selected genes.
