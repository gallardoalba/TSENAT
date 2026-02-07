# Plot MA-like and Volcano plots for lm_interaction results

These functions create diagnostic plots for linear model interaction
analysis results, showing interaction effect size vs adjusted p-value.

## Usage

``` r
plot_ma_lm_interaction(lm_res, diff_res = NULL, se = NULL)

plot_volcano_lm_interaction(lm_res)
```

## Arguments

- lm_res:

  A data frame of results from \`calculate_lm_interaction()\`, with
  columns including \`gene\`, \`estimate_interaction\`, and
  \`adj_p_interaction\`.

## Value

A ggplot2 object.
