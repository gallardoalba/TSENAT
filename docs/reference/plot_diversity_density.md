# Plot diversity distributions (density) by sample type

Plot diversity distributions (density) by sample type

## Usage

``` r
plot_diversity_density(se, assay_name = "diversity", sample_type_col = NULL)
```

## Arguments

- se:

  A \`SummarizedExperiment\` returned by \`calculate_diversity\`.

- assay_name:

  Name of the assay to use (default: "diversity").

- sample_type_col:

  Optional column name in \`colData(se)\` that contains sample types. If
  missing, sample type will be inferred from column names (suffix after
  last underscore) or classified as 'Group'.

## Value

A \`ggplot\` object with layered density plots.
