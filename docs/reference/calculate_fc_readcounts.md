# Compute fold changes from read counts with transcript collapsing.

This function aggregates read counts at the gene level (collapsing
multiple transcripts per gene) and then computes fold changes between
two conditions. Useful for plotting read count changes in MA plots
alongside diversity-based p-values.

## Usage

``` r
calculate_fc_readcounts(
  readcounts,
  genes = NULL,
  samples,
  control,
  agg_method = "mean",
  fc_method = "mean",
  pseudocount = 0
)
```

## Arguments

- readcounts:

  A data.frame or matrix with read counts. If a data.frame, the first
  column should contain gene identifiers. If a matrix, rownames should
  contain gene identifiers.

- genes:

  Character vector of gene names corresponding to each row of
  readcounts. If readcounts is a data.frame with gene names in the first
  column, this can be omitted.

- samples:

  Character vector with length equal to the number of count columns,
  specifying the category of each sample.

- control:

  Name of the control sample category, as defined in the \`samples\`
  vector (e.g., \`'Normal'\` or \`'Control'\`).

- agg_method:

  Method to use for aggregating transcripts to genes. Can be \`'mean'\`
  (default) or \`'median'\`.

- fc_method:

  Method for computing fold changes (default: \`'mean'\`). Can be
  \`'mean'\` or \`'median'\`.

- pseudocount:

  Numeric scalar passed to \`calculate_fc()\`. Default is 0 (scale-aware
  pseudocount).

## Value

A data.frame with columns:

- genes:

  Gene identifiers.

- \<treatment\>\_\<fc_method\>:

  Aggregated value for treatment condition.

- \<control\>\_\<fc_method\>:

  Aggregated value for control condition.

- \<fc_method\>\_difference:

  Simple difference.

- log2_fold_change:

  log2(treatment / control) fold change.

## Details

If readcounts contains multiple rows per gene (e.g., different
transcripts), they are first aggregated to the gene level using the
specified aggregation method (mean or median). Then fold changes are
computed using \`calculate_fc()\`.

## Examples

``` r
# Create example read count data with transcripts
rc_data <- data.frame(
    gene = c("GENE1", "GENE1", "GENE2", "GENE2"),
    sample1 = c(100, 50, 200, 100),
    sample2 = c(120, 60, 180, 90),
    sample3 = c(80, 40, 250, 120),
    sample4 = c(90, 45, 270, 130)
)
samples <- c("control", "control", "treatment", "treatment")
calculate_fc_readcounts(rc_data, samples = samples,
                        control = "control", agg_method = "mean")
#>       treatment_mean control_mean mean_difference log2_fold_change genes
#> GENE1          63.75         82.5          -18.75       -0.3719688 GENE1
#> GENE2         192.50        142.5           50.00        0.4338965 GENE2
```
