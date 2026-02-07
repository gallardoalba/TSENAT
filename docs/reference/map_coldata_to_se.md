# Map external coldata into a SummarizedExperiment

Map an external \`coldata\` table (sample IDs and a condition/label
column) into a \`SummarizedExperiment\` produced by
\`calculate_diversity()\`. The helper sets
\`colData(ts_se)\$sample_type\` when possible and records a
\`sample_base\` identifier for each column. Unmatched entries remain NA
and no automatic inference is attempted.

## Usage

``` r
map_coldata_to_se(
  ts_se,
  coldata,
  coldata_sample_col = "Sample",
  coldata_condition_col = "Condition",
  paired = FALSE
)
```

## Arguments

- ts_se:

  A \`SummarizedExperiment\` whose assay columns are named with sample
  identifiers. Names may include per-sample annotations; mapping is
  exact and case-sensitive unless you normalize identifiers beforehand.

- coldata:

  A data.frame with sample metadata (or \`NULL\`).

- coldata_sample_col:

  Name of the column in \`coldata\` containing sample identifiers
  (default: "Sample").

- coldata_condition_col:

  Name of the column in \`coldata\` with condition/labels (default:
  "Condition").

- paired:

  Logical; if \`TRUE\`, validate pairing and reorder columns so that
  matched samples for each base are adjacent (default: \`FALSE\`). Use
  \`paired = TRUE\` for paired analyses so downstream paired tests see
  aligned samples regardless of the original \`coldata\` order.

## Value

The input \`ts_se\` with \`colData(ts_se)\$sample_type\` populated when
possible. \`colData(ts_se)\$sample_base\` is also added containing the
base sample identifier; rownames of \`colData(ts_se)\` are aligned to
the assay column names.

## Examples

``` r
data("tcga_brca_luma_dataset", package = "TSENAT")
rc <- as.matrix(tcga_brca_luma_dataset[1:20, -1, drop = FALSE])
gs <- tcga_brca_luma_dataset$genes[1:20]
se <- calculate_diversity(rc, gs, q = 0.1, norm = TRUE)
sample_names <- sub("_q=.*", "", colnames(SummarizedExperiment::assay(se)))
coldata_df <- data.frame(Sample = sample_names, Condition = rep(c("A", "B"),
    length.out = ncol(se)
))
map_coldata_to_se(se, coldata_df)
#> class: SummarizedExperiment 
#> dim: 6 40 
#> metadata(4): method norm q what
#> assays(1): diversity
#> rownames(6): MXRA8 C1orf86 ... HNRNPR C1orf213
#> rowData names(1): genes
#> colnames(40): TCGA-A7-A0CH_N TCGA-A7-A0CH_T ... TCGA-BH-A0BV_T
#>   TCGA-BH-A0BV_N
#> colData names(3): samples sample_type sample_base
# Optionally validate pairs when appropriate
map_coldata_to_se(se, coldata_df, paired = TRUE)
#> class: SummarizedExperiment 
#> dim: 6 40 
#> metadata(4): method norm q what
#> assays(1): diversity
#> rownames(6): MXRA8 C1orf86 ... HNRNPR C1orf213
#> rowData names(1): genes
#> colnames(40): TCGA-A7-A0CH_N TCGA-A7-A0CH_T ... TCGA-BH-A0BV_T
#>   TCGA-BH-A0BV_N
#> colData names(3): samples sample_type sample_base
```
