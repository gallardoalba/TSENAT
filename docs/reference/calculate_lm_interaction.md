# Linear-model interaction test for Tsallis entropy For each gene, fit a linear model of the form \`entropy ~ q \* group\` and extract the p-value for the interaction term (whether the effect of \`q\` differs between groups). The function expects a \`SummarizedExperiment\` produced by \`calculate_diversity()\` when multiple \`q\` values have been computed (column names contain \`\_q=\`).

Linear-model interaction test for Tsallis entropy For each gene, fit a
linear model of the form \`entropy ~ q \* group\` and extract the
p-value for the interaction term (whether the effect of \`q\` differs
between groups). The function expects a \`SummarizedExperiment\`
produced by \`calculate_diversity()\` when multiple \`q\` values have
been computed (column names contain \`\_q=\`).

## Usage

``` r
calculate_lm_interaction(
  se,
  sample_type_col = NULL,
  min_obs = 10,
  method = c("linear", "gam", "fpca"),
  nthreads = 1,
  assay_name = "diversity"
)
```

## Arguments

- se:

  A \`SummarizedExperiment\` containing a \`diversity\` assay produced
  by \`calculate_diversity(..., q = \<vector\>)\`.

- sample_type_col:

  Optional column name in \`colData(se)\` that contains a grouping
  factor for samples (character). If \`NULL\`, the function will attempt
  to infer group from column names (suffix \`\_N\` interpreted as
  "Normal").

- min_obs:

  Minimum number of non-NA observations required to fit a model for a
  gene (default: 10).

- method:

  Modeling method to use for interaction testing: one of
  `c("linear", "gam", "fpca")`.

- nthreads:

  Number of threads (mc.cores) to use for parallel processing (default:
  1).

- assay_name:

  Name of the assay in the SummarizedExperiment to use (default:
  "diversity").

## Value

A data.frame with columns \`gene\`, \`p_interaction\`, and
\`adj_p_interaction\`, ordered by ascending \`p_interaction\`.

## Examples

``` r
data("tcga_brca_luma_dataset", package = "TSENAT")
rc <- as.matrix(tcga_brca_luma_dataset[1:20, -1, drop = FALSE])
gs <- tcga_brca_luma_dataset$genes[1:20]
se <- calculate_diversity(rc, gs, q = c(0.1, 1), norm = TRUE)
calculate_lm_interaction(se)
#> [calculate_lm_interaction] method=linear
#> [calculate_lm_interaction] parsed samples and groups
#>       gene p_interaction adj_p_interaction
#> 1  C1orf86     0.2415774         0.6590893
#> 2     PDPN     0.4229195         0.6590893
#> 3   HNRNPR     0.4658805         0.6590893
#> 4 C1orf213     0.5487154         0.6590893
#> 5    MXRA8     0.5492411         0.6590893
#> 6  ALDH4A1     0.6606980         0.6606980
```
