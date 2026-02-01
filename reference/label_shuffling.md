# Calculate p-values using label shuffling.

Calculate p-values using label shuffling.

## Usage

``` r
label_shuffling(
  x,
  samples,
  control,
  method,
  randomizations = 100,
  pcorr = "BH"
)
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

- randomizations:

  The number of random shuffles.

- pcorr:

  P-value correction method applied to the results, as defined in the
  `p.adjust` function.

## Value

Raw and corrected p-values.

## Details

The permutation p-values are computed two-sided as the proportion of
permuted log2 fold-changes at least as extreme as the observed value,
with a pseudocount added: (count + 1) / (n_perm + 1).

## Note

The permutation test returns two-sided empirical p-values using a
pseudocount to avoid zero p-values for small numbers of permutations.
See the function documentation for details.
