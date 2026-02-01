# TSENAT: Tsallis Entropy Analysis Toolbox

TSENAT analyze expression/transcript differences and compute diversity
metrics.

## Tsallis theory

Tsallis entropy generalizes Shannon entropy. For a probability vector
$`p = (p_1, \dots, p_n)`$ (with $`p_i \ge 0`$ and $`\sum_i p_i = 1`$)
the Tsallis entropy of order $`q`$ is defined for $`q \ne 1`$ as

``` math
S_q(p) = \frac{1 - \sum_{i} p_i^q}{q - 1}.
```

In the limit $`q \to 1`$ this recovers the Shannon entropy

``` math
\lim_{q \to 1} S_q(p) = -\sum_i p_i \log p_i.
```

In the transcript-expression context, $`p_i`$ are the normalized
abundances of isoforms for a gene (i.e. nonnegative and summing to 1).
The parameter `q` controls sensitivity to isoform abundance: $`q<1`$
emphasizes rare isoforms, $`q>1`$ emphasizes dominant isoforms. Common
interpretations are that $`q=0`$ corresponds to isoform richness (the
count of nonzero isoforms), and $`q=2`$ relates to the inverse Simpson
index. In TSENAT we compute the Tsallis entropy $`S_q`$ per gene to
provide an interpretable diversity measure; users can compute $`S_q`$
for multiple $`q`$ values and choose whether to normalize raw counts to
proportions before analysis.

## Integrated features

- Tsallis entropy and diversity calculations:
  - `calculate_tsallis_entropy`: computes S_q and/or D_q for a numeric
    vector of expression values (supports normalization, multiple `q`
    values and the q→1 limit). Returns numeric vectors or a list
    depending on `what`.
  - `calculate_diversity`: applies the calculation across
    transcripts/genes for matrices, `tximport`-style lists or
    `SummarizedExperiment` objects and returns a `SummarizedExperiment`
    with assay `diversity` (S_q) or `hill` (D_q).
- Differential and statistical analyses:
  - `calculate_difference` and helpers in `difference_functions` compute
    group means, differences (or log2 fold-changes), p-values and
    adjusted p-values. These functions are designed to work with
    diversity summaries as well as expression matrices.
- Plotting and visualization:
  - `plot_tsallis_q_curve`: median ± IQR of Tsallis entropy across
    q-values by group.

## Installation and documentation

Install from GitHub during development:

``` r

if (!requireNamespace("remotes", quietly = TRUE)) install.packages("remotes")
remotes::install_github("gallardoalba/TSENAT")

## Recommended (reproducible): use `renv`.

install.packages("renv")
renv::init()
# when dependencies are set, record lock
renv::snapshot()
```

## Quick start

Compute Tsallis diversity for a single `q` and plot a q-curve across
multiple `q` values (small, focused example taken from the vignette):

``` r

library(TSENAT)
data("tcga_brca_luma_dataset", package = "TSENAT")

# compute for q = 0.1 and normalize
readcounts <- tcga_brca_luma_dataset$counts
genes <- tcga_brca_luma_dataset$gene
ts_se <- calculate_diversity(readcounts, genes, method = "tsallis", q = 0.1, norm = TRUE)

# compute for multiple q values and plot the q-curve
qvec <- seq(0.01, 2, by = 0.1)
ts_multi <- calculate_diversity(readcounts, genes, method = "tsallis", q = qvec, norm = TRUE)
plot_tsallis_q_curve(ts_multi, group = colData(ts_multi)$sample_type)
```

![q-curve example](articles/TSENAT_files/figure-html/plot-q-curve-1.png)

q-curve example

## Licence, citation and atribution

This project is licensed under the GNU General Public License v3.0
(GPL-3). A copy of the full license text is included in the repository
at `LICENSE` and is installed with the package under `inst/LICENSE`.

To cite TSENAT in publications, run the following in R:

``` r

citation("TSENAT")
```

This command shows the recommended bibliographic entry; a
machine-readable `inst/CITATION` file is included with the package for
easy citation export.

TSENAT builds upon and adapts substantial portions of code from the
[SplicingFactory package](https://github.com/esebesty/SplicingFactory).
The codebase has been extended with additional utilities, bug fixes and
plotting helpers focused on Tsallis based transcript-level diversity
analysis.
