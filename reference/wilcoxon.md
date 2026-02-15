# Calculate p-values using Wilcoxon rank sum test.

Calculate p-values using Wilcoxon rank sum test.

## Usage

``` r
wilcoxon(x, samples, pcorr = "BH", paired = FALSE, exact = FALSE, nthreads = 1)
```

## Arguments

- x:

  A `matrix` with the splicing diversity values.

- samples:

  Character vector with an equal length to the number of columns in the
  input dataset, specifying the category of each sample.

- pcorr:

  P-value correction method applied to the results, as defined in the
  `p.adjust` function.

- paired:

  If `TRUE`, the Wilcox-test will be paired, and therefore it will be a
  signed rank test instead of the rank sum test.

- exact:

  If `TRUE`, an exact p-value will be computed.

- nthreads:

  Number of threads for parallel processing (default: 1). Set to \> 1 to
  parallelize per-feature Wilcoxon tests.

## Value

Raw and corrected p-values in a matrix.

## Examples

``` r
# Create a matrix of splicing diversity values (3 genes x 6 samples)
mat <- matrix(rnorm(18), nrow = 3)
samples <- rep(c('Control', 'Treatment'), each = 3)

# Run Wilcoxon test
result <- wilcoxon(mat, samples, pcorr = 'BH')
head(result)
#>      raw_p_values adjusted_p_values
#> [1,]    0.0808556         0.2425668
#> [2,]    0.6625206         0.6625206
#> [3,]    0.3827331         0.5740996
```
