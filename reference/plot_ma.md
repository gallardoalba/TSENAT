# Plot MA plot for difference results

Plot MA plot for difference results

## Usage

``` r
plot_ma(
  diff_df,
  mean_cols = NULL,
  fold_col = "log2_fold_change",
  padj_col = "adjusted_p_values",
  sig_alpha = 0.05
)
```

## Arguments

- diff_df:

  Data.frame from \`calculate_difference()\` containing mean columns, a
  \`log2_fold_change\` column, and \`adjusted_p_values\`.

- mean_cols:

  Optional character vector of length 2 naming the mean columns.
  Defaults to the first two columns that end with \`\_mean\`.

- fold_col:

  Name of the fold-change column (default: \`log2_fold_change\`).

- padj_col:

  Name of the adjusted p-value column (default: \`adjusted_p_values\`).

- sig_alpha:

  Threshold for significance (default: 0.05).

## Value

A \`ggplot\` MA-plot object.

## Examples

``` r
df <- data.frame(
    gene = paste0("g", seq_len(10)),
    sampleA_mean = runif(10),
    sampleB_mean = runif(10),
    log2_fold_change = rnorm(10),
    adjusted_p_values = runif(10)
)
plot_ma(df)
```
