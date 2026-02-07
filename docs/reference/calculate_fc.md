# Calculate fold changes independently of statistical testing.

Computes mean or median fold changes and differences between two sample
conditions. This function is independent of statistical testing and can
be used to calculate fold changes for any comparison of interest.

## Usage

``` r
calculate_fc(x, samples, control, method = "mean", pseudocount = 0)
```

## Arguments

- x:

  A `matrix` with the splicing diversity values (rows = features,
  columns = samples).

- samples:

  Character vector with length equal to the number of columns in `x`,
  specifying the category of each sample.

- control:

  Name of the control sample category, as defined in the `samples`
  vector (e.g., `'Normal'` or `'WT'`). The fold change is computed as
  treatment / control.

- method:

  Method to use for aggregating values within each condition. Can be
  `'mean'` or `'median'`. Default: `'mean'`.

- pseudocount:

  Numeric scalar. Small value added to non-positive observed group
  summaries to avoid zeros when computing differences and log2
  fold-changes. If `pseudocount <= 0` (default), the function
  automatically chooses a scale-aware value equal to half the smallest
  positive observed group summary across both groups. If no positive
  values are present, the fallback is `1e-6`. Rows with all `NA` values
  remain `NA` and are not imputed.

## Value

A `data.frame` with columns:

- \<treatment\>\_\<method\>:

  Aggregated value (mean or median) for the treatment condition.

- \<control\>\_\<method\>:

  Aggregated value (mean or median) for the control condition.

- \<method\>\_difference:

  Simple difference between treatment and control.

- log2_fold_change:

  log2(treatment / control) fold change.

## Details

The function computes group-wise summaries (mean or median) for each
condition, then calculates the difference and log2 fold change.
Pseudocounts are applied before log-transformation to handle
non-positive values. Rows with insufficient non-NA observations in
either condition remain `NA`.

## Examples

``` r
# Create example diversity matrix
mat <- matrix(c(2, 3, 1, 0.5, 4, 2.5, 1.5, 3.5), nrow = 2, ncol = 4)
samples <- c("control", "control", "treatment", "treatment")
calculate_fc(mat, samples, control = "control")
#>    treatment_mean control_mean mean_difference log2_fold_change
#> V1           2.75         1.50            1.25        0.8744691
#> V2           3.00         1.75            1.25        0.7776076
```
