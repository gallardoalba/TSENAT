# Plot violin of per-gene mean diversity by sample type

Plot violin of per-gene mean diversity by sample type

## Usage

``` r
plot_mean_violin(se, assay_name = "diversity", sample_type_col = NULL)
```

## Arguments

- se:

  A \`SummarizedExperiment\` returned by \`calculate_diversity\`.

- assay_name:

  Name of the assay to use (default: "diversity").

- sample_type_col:

  Optional column name in \`colData(se)\` containing sample types.

## Value

A \`ggplot\` violin plot object.
