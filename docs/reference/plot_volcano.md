# Volcano plot for differential results

Volcano plot for differential results

## Usage

``` r
plot_volcano(
  diff_df,
  x_col = "mean_difference",
  padj_col = "adjusted_p_values",
  label_thresh = 0.1,
  padj_thresh = 0.05,
  top_n = 5
)
```

## Arguments

- diff_df:

  Data.frame from \`calculate_difference()\` containing at least
  \`mean_difference\` and an adjusted p-value column (default
  \`adjusted_p_values\`).

- x_col:

  Column name for x-axis (default \`mean_difference\`).

- padj_col:

  Column name for adjusted p-values (default \`adjusted_p_values\`).

- label_thresh:

  Absolute x threshold to mark significance (default 0.1).

- padj_thresh:

  Adjusted p-value cutoff (default 0.05).

- top_n:

  Integer; number of top genes to annotate by smallest adjusted p-value
  (default: 5).

## Value

ggplot volcano plot.
