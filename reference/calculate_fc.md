# Calculate splicing diversity changes between two conditions.

Calculate splicing diversity changes between two conditions.

## Usage

``` r
calculate_fc(x, samples, control, method = "mean", pseudocount = 0)
```

## Arguments

- x:

  A `matrix` with the splicing diversity values.

- samples:

  Character vector with an equal length to the number of columns in the
  input dataset, specifying the category of each sample.

- control:

  Name of the control sample category, defined in the `samples` vector,
  e.g. `control = 'Normal'` or `control = 'WT'`.

- method:

  Method to use for calculating the average splicing diversity value in
  a condition. Can be `'mean'` or `'median'`.

- pseudocount:

  Numeric scalar. Small value added to non-positive observed group
  summaries to avoid zeros when computing differences and log2
  fold-changes. If `pseudocount <= 0` the function will automatically
  choose a scale-aware value equal to half the smallest positive
  observed group summary (i.e. half the smallest observed mean/median
  across groups); if no positive values are present the fallback is
  `1e-6`. Rows with insufficient observations remain `NA` and are not
  imputed.

## Value

A `data.frame` with mean or median value of splicing diversity across
sample categories, the difference between these values and the log2 fold
change values.

## Details

The function uses a matrix of splicing diversity values in order to
calculate mean or median differences and log2 fold changes between two
conditions.

## Examples

``` r
# Simulate splicing diversity matrix (4 genes x 4 samples)
mat <- matrix(c(
  0.5, 0.6, 0.8, 0.9,  # gene1: low control, high treatment
  0.7, 0.75, 0.6, 0.5  # gene2: high control, low treatment
), nrow = 2, byrow = TRUE)
samples <- c('Normal', 'Normal', 'Tumor', 'Tumor')
result <- calculate_fc(mat, samples, control = 'Normal', method = 'mean')
head(result)
#>    Tumor_mean Normal_mean mean_difference log2_fold_change
#> V1       0.85       0.550           0.300        0.6280312
#> V2       0.55       0.725          -0.175       -0.3985494
```
