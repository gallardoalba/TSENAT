# Volcano plot for differential results

Volcano plot for differential results

## Usage

``` r
plot_volcano(
  diff_df,
  x_col = NULL,
  padj_col = "adjusted_p_values",
  label_thresh = 0.1,
  padj_thresh = 0.05,
  top_n = 5,
  title = NULL
)
```

## Arguments

- diff_df:

  Data.frame from \`test_differential()\` or similar differential
  analysis results. Should contain p-value and optionally fold-change
  columns.

- x_col:

  Column name for x-axis values. If not specified, will auto-detect from
  available columns (e.g., "mean_difference", "median_difference").

- padj_col:

  Column name for adjusted p-values (default: "adjusted_p_values").

- label_thresh:

  Threshold for fold-change labeling (default: 0.1).

- padj_thresh:

  P-value threshold for significance (default: 0.05).

- top_n:

  Number of top significant genes to label (default: 5).

- title:

  Custom plot title (optional). If NULL, generates default title.

## Value

ggplot volcano plot.

A \`ggplot2\` object.

## Examples
