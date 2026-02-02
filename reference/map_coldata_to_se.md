# Map external coldata into a SummarizedExperiment

This helper maps an external \`coldata\` table (with sample IDs and a
condition/label column) into a \`SummarizedExperiment\` produced by
\`calculate_diversity()\`. It assigns a \`sample_type\` column to
\`colData(ts_se)\` and falls back to \`infer_sample_group()\` for
unmapped samples.

## Usage

``` r
map_coldata_to_se(
  ts_se,
  coldata,
  coldata_sample_col = "Sample",
  coldata_condition_col = "Condition"
)
```

## Arguments

- ts_se:

  A \`SummarizedExperiment\` object with assay columns named possibly
  including \`\_q=\` suffixes.

- coldata:

  A data.frame with sample metadata or \`NULL\`.

- coldata_sample_col:

  Name of the column in \`coldata\` with sample IDs.

- coldata_condition_col:

  Name of the column in \`coldata\` with condition/labels.

## Value

The input \`ts_se\` with \`colData(ts_se)\$sample_type\` set when
possible.

## Examples

``` r
data("tcga_brca_luma_dataset", package = "TSENAT")
rc <- as.matrix(tcga_brca_luma_dataset[1:20, -1, drop = FALSE])
gs <- tcga_brca_luma_dataset$genes[1:20]
se <- calculate_diversity(rc, gs, q = 0.1, norm = TRUE)
sample_names <- sub("_q=.*", "", colnames(SummarizedExperiment::assay(se)))
coldata_df <- data.frame(
  Sample = sample_names,
  Condition = rep(c("A", "B"), length.out = ncol(se))
)
map_coldata_to_se(se, coldata_df)
#> class: SummarizedExperiment 
#> dim: 6 40 
#> metadata(4): method norm q what
#> assays(1): diversity
#> rownames(6): MXRA8 C1orf86 ... HNRNPR C1orf213
#> rowData names(1): genes
#> colnames(40): TCGA-A7-A0CH_N TCGA-A7-A0CH_T ... TCGA-BH-A0BV_T
#>   TCGA-BH-A0BV_N
#> colData names(2): samples sample_type
```
