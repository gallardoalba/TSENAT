# Violin plot of Tsallis entropy for multiple q values

Violin plot of Tsallis entropy for multiple q values

## Usage

``` r
plot_tsallis_violin_multq(se, assay_name = "diversity")
```

## Arguments

- se:

  A \`SummarizedExperiment\` returned by \`calculate_diversity\` with
  multiple q values (column names contain \`\_q=\`).

- assay_name:

  Name of the assay to use (default: "diversity").

## Value

A \`ggplot\` violin plot object faceted/colored by group and q.

## Examples

``` r
data("tcga_brca_luma_dataset", package = "TSENAT")
rc <- as.matrix(tcga_brca_luma_dataset[1:20, -1, drop = FALSE])
gs <- tcga_brca_luma_dataset$genes[1:20]
se <- calculate_diversity(rc, gs, q = c(0.1, 1), norm = TRUE)
plot_tsallis_violin_multq(se)
```
