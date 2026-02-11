# Volcano plot for differential results

Create a volcano plot showing fold-change (x-axis) versus adjusted
p-value significance (y-axis). The function auto-detects a suitable
x-axis column if one is not provided and expects an adjusted p-value
column for significance coloring.

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

  Data.frame from \`test_differential()\` or similar results. Should
  contain p-values and optionally fold-change columns.

- x_col:

  Optional column name for the x-axis. If \`NULL\`, the function will
  try to auto-detect a suitable numeric column (excluding p-values).

- padj_col:

  Adjusted p-value column name (default: "adjusted_p_values").

- label_thresh:

  Fold-change threshold used to annotate points (default: 0.1).

- padj_thresh:

  Adjusted p-value cutoff for significance (default: 0.05).

- top_n:

  Number of top significant genes to label (default: 5).

- title:

  Optional plot title; if \`NULL\` a default title is used.

## Value

A \`ggplot2\` object.

## Examples

``` r
df <- data.frame(
    gene = paste0("g", seq_len(10)),
    mean_difference = runif(10),
    adjusted_p_values = runif(10)
)
# plot_volcano(df, x_col = "mean_difference", padj_col = "adjusted_p_values")
```
