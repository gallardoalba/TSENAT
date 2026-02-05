# Plot median +- IQR of Tsallis entropy across q values by group

This reproduces the \`tsallis-q-curve-mean-sd\` plot from the vignette:
for each q value, compute per-gene Tsallis entropy per sample, summarize
across genes by group (median and IQR) and plot median with a ribbon
spanning median +- IQR/2.

## Usage

``` r
plot_tsallis_q_curve(
  se,
  assay_name = "diversity",
  sample_type_col = "sample_type",
  y = c("S", "D")
)
```

## Arguments

- se:

  A \`SummarizedExperiment\` returned by \`calculate_diversity\`.

- assay_name:

  Name of the assay to use (default: "diversity").

- sample_type_col:

  Column name in \`colData(se)\` containing sample type labels (default:
  "sample_type").

## Value

A \`ggplot\` object showing median +- IQR across q values by group.

## Examples

``` r
data("tcga_brca_luma_dataset", package = "TSENAT")
rc <- as.matrix(tcga_brca_luma_dataset[1:40, -1, drop = FALSE])
gs <- tcga_brca_luma_dataset$genes[1:40]
se <- calculate_diversity(rc, gs,
   q = seq(0.01, 0.1, by = 0.03), norm = TRUE)
p <- plot_tsallis_q_curve(se)
p
```
