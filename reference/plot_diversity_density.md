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

  Optional column name in \`colData(se)\` with sample types. If missing,
  sample types are inferred from column names (suffix after the last
  underscore) or set to 'Group'.

## Value

A \`ggplot\` object with layered density plots.

## Examples

``` r
data("tcga_brca_luma_dataset", package = "TSENAT")
rc <- as.matrix(tcga_brca_luma_dataset[1:20, -1, drop = FALSE])
gs <- tcga_brca_luma_dataset$genes[1:20]
se <- calculate_diversity(rc, gs, q = 0.1, norm = TRUE)
plot_diversity_density(se)
```
