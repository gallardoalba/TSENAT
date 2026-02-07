# Plot MA plot for difference or linear model results

Creates an MA plot (Mean vs log-fold-change) that works for both
differential analysis and linear model interaction results.

## Usage

``` r
plot_ma(
  x,
  fc_df = NULL,
  sig_alpha = 0.05,
  x_label = NULL,
  y_label = NULL,
  title = NULL
)
```

## Arguments

- x:

  Data.frame from \`test_differential()\` containing adjusted p-values
  and fold changes.

- fc_df:

  Optional data.frame with fold changes from alternative methods (e.g.,
  read count fold changes). Should have 'genes' and 'log2_fold_change'
  columns. If provided, will override fold changes in \`x\`.

- sig_alpha:

  Significance threshold for p-values (default: 0.05).

- x_label:

  Custom x-axis label (optional).

- y_label:

  Custom y-axis label (optional).

- title:

  Custom plot title (optional). If NULL, generates default title.

- mean_cols:

  Optional character vector of length 2 naming the mean columns.
  Defaults to the first two columns that end with \`\_mean\` or
  \`\_median\`. Ignored for linear model input if x-axis mapping is not
  applicable.

- fold_col:

  Name of the fold-change column (default: \`log2_fold_change\`). Also
  searches for \`logFC\`, \`estimate_interaction\`, etc.

- padj_col:

  Name of the adjusted p-value column (default: \`adjusted_p_values\`).
  Also searches for \`adj_p_interaction\`, etc.

- diff_res:

  Optional data.frame from \`calculate_difference()\`. Used with linear
  model results (\`x\`) to obtain fold changes and mean values for
  plotting. Allows plotting linear model p-values with differential fold
  changes.

## Value

A \`ggplot\` MA-plot object.

A \`ggplot2\` object.

## Details

Automatically detects the type of input data (differential vs linear
model) and creates appropriate MA plot. Can also combine results from
different workflows, e.g., plotting linear model p-values with
differential fold changes, or read count fold changes with
diversity-based p-values.

## Examples
