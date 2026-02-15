# Run a differential test by name: Wilcoxon or label-shuffle

Thin wrapper that selects between \`wilcoxon()\` and
\`label_shuffling()\`.

## Usage

``` r
test_differential(
  x,
  samples,
  control = NULL,
  method = c("wilcoxon", "shuffle"),
  fc_method = "mean",
  paired = FALSE,
  exact = FALSE,
  randomizations = 100,
  pcorr = "BH",
  seed = 123L,
  paired_method = c("swap", "signflip"),
  nthreads = 1
)
```

## Arguments

- x:

  A matrix with splicing diversity values (rows = features).

- samples:

  Character vector of sample group labels (length = ncol(x)).

- control:

  Name of the control group (required for label-shuffle).

- method:

  Character; one of \`'wilcoxon'\` or \`'shuffle'\`.

- fc_method:

  Character; aggregation method used by the permutation test when
  \`method = 'shuffle'\` ('mean' or 'median').

- paired:

  Logical passed to \`wilcoxon()\` when using the Wilcoxon test.

- exact:

  Logical passed to \`wilcoxon()\` to request exact p-values.

- randomizations:

  Integer number of permutations for \`label_shuffling()\`.

- pcorr:

  P-value adjustment method (passed to \`p.adjust\`).

- seed:

  Integer seed used to make permutations reproducible (default 123). The
  function sets a temporary RNG seed via \`withr::local_seed(seed)\`
  before running \`label_shuffling()\` when \`method = 'shuffle'\`.

- paired_method:

  Character; forwarded to \`label_shuffling()\` when \`method =
  'shuffle'\`. See \`label_shuffling()\` for details.

- nthreads:

  Number of threads for parallel processing (default: 1). Set to \> 1 to
  parallelize per-feature statistical tests in \`wilcoxon()\` or
  \`label_shuffling()\`.

## Value

A two-column matrix with raw and adjusted p-values (as returned by the
underlying functions).

## Examples

``` r
mat <- matrix(rnorm(20), nrow = 5)
samples <- rep(c('A', 'B'), length.out = ncol(mat))
test_differential(mat, samples, control = 'A', method = 'wilcoxon')
#>      raw_p_values adjusted_p_values
#> [1,]    0.6985354         0.8731692
#> [2,]    0.2452781         0.6131953
#> [3,]    0.2452781         0.6131953
#> [4,]    0.6985354         0.8731692
#> [5,]    1.0000000         1.0000000
```
