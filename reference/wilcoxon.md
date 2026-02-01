# Calculate p-values using Wilcoxon rank sum test.

Calculate p-values using Wilcoxon rank sum test.

## Usage

``` r
wilcoxon(x, samples, pcorr = "BH", paired = FALSE, exact = FALSE)
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

## Value

Raw and corrected p-values in a matrix.
