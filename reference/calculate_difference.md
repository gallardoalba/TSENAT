# Calculate splicing diversity changes between two conditions.

Calculate splicing diversity changes between two conditions.

## Usage

``` r
calculate_difference(
  x,
  samples = NULL,
  control,
  method = "mean",
  test = "wilcoxon",
  randomizations = 100,
  pcorr = "BH",
  assayno = 1,
  verbose = TRUE,
  paired = FALSE,
  exact = FALSE,
  pseudocount = 0
)
```

## Arguments

- x:

  A `SummarizedExperiment` with splicing diversity values for each gene
  in each sample or a `data.frame` with gene names in the first column
  and splicing diversity values for each sample in additional columns.

- samples:

  A vector of length one, specifying the column name of the `colData`
  annotation column from the `SummarizedExperiment` object, that should
  be used as the category column or a character vector with an equal
  length to the number of columns in the input dataset, specifying the
  category of each sample in the case of a `data.frame` input.

- control:

  Name of the control sample category, defined in the `samples` vector,
  e.g. `control = 'Normal'` or `control = 'WT'`.

- method:

  Method to use for calculating the average splicing diversity value in
  a condition. Can be `'mean'` or `'median'`.

- test:

  Method to use for p-value calculation: use `'wilcoxon'` for Wilcoxon
  rank sum test or `'shuffle'` for a label shuffling test.

- randomizations:

  Number of random shuffles, used for the label shuffling test (default
  = 100).

- pcorr:

  P-value correction method applied to the Wilcoxon rank sum test or
  label shuffling test results, as defined in the `p.adjust` function.

- assayno:

  An integer value. In case of multiple assays in a
  `SummarizedExperiment` input, the argument specifies the assay number
  to use for difference calculations.

- verbose:

  If `TRUE`, the function will print additional diagnostic messages.

- paired:

  Logical; if \`TRUE\`, run paired versions of tests when supported
  (default: \`FALSE\`).

- exact:

  Logical; passed to the Wilcoxon test to request exact p-values when
  supported (default: \`FALSE\`).

- pseudocount:

  Numeric scalar. Passed to `calculate_fc` and used to add a small value
  to non-positive group summaries before computing differences and log2
  fold-changes. Default `1e-6`. Rows excluded for low sample counts
  remain `NA`.

## Value

A `data.frame` with the mean or median values of splicing diversity
across sample categories and all samples, log2(fold change) of the two
different conditions, raw and corrected p-values.

## Details

The function calculates diversity changes between two sample conditions.
It uses the output of the diversity calculation function, which is a
`SummarizedExperiment` object of splicing diversity values.
Additionally, it can use a `data.frame` as input, where the first column
contains gene names, and all additional columns contain splicing
diversity values for each sample. A vector of sample conditions also
serves as input, used for aggregating the samples by condition. It
calculates the mean or median of the splicing diversity data per sample
condition, the difference of these values and the log2 fold change of
the two conditions. Furthermore, the user can select a statistical
method to calculate the significance of the changes. The p-values and
adjusted p-values are calculated using a Wilcoxon sum rank test or label
shuffling test. The function will exclude genes of low sample size from
the significance calculation, depending on which statistical test is
applied.

## Examples

``` r
x <- data.frame(Genes = letters[seq_len(10)], matrix(runif(80), ncol = 8))
samples <- c(rep("Healthy", 4), rep("Pathogenic", 4))
calculate_difference(x, samples,
    control = "Healthy", method = "mean", test =
        "wilcoxon"
)
#>    genes Pathogenic_mean Healthy_mean mean_difference log2_fold_change
#> 1      a       0.3939888    0.4673230    -0.073334138      -0.24626521
#> 2      b       0.4391230    0.4589390    -0.019816020      -0.06367737
#> 3      c       0.5746470    0.4147989     0.159848098       0.47026392
#> 4      d       0.7016257    0.4580544     0.243571275       0.61518255
#> 5      e       0.6093560    0.4994365     0.109919512       0.28698410
#> 6      f       0.3516991    0.6633813    -0.311682200      -0.91549669
#> 7      g       0.5652513    0.4689715     0.096279798       0.26939212
#> 8      h       0.3903022    0.3944323    -0.004130123      -0.01518619
#> 9      i       0.2585383    0.4280178    -0.169479477      -0.72729271
#> 10     j       0.4058325    0.4532731    -0.047440614      -0.15949619
#>    raw_p_values adjusted_p_values
#> 1     1.0000000                 1
#> 2     1.0000000                 1
#> 3     0.6650055                 1
#> 4     0.1939309                 1
#> 5     0.6650055                 1
#> 6     0.3123214                 1
#> 7     0.8852339                 1
#> 8     0.8852339                 1
#> 9     0.3123214                 1
#> 10    0.8852339                 1
```
