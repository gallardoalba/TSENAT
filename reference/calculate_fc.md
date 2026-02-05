# Calculate splicing diversity changes between two conditions.

Calculate splicing diversity changes between two conditions.

## Usage

``` r
calculate_fc(x, samples, control, method = "mean", pseudocount = 1e-06)
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
  fold-changes. If `pseudocount = 0` (default) the function will choose
  a small value equal to half the smallest positive observed value, or
  `1e-6` if no positive values are present. Rows with insufficient
  observations remain `NA` and are not imputed.

## Value

A `data.frame` with mean or median value of splicing diversity across
sample categories, the difference between these values and the log2 fold
change values.

## Details

The function uses a matrix of splicing diversity values in order to
calculate mean or median differences and log2 fold changes between two
conditions.
