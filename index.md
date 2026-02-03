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

In transcript-expression data we compute Tsallis entropy per gene from
the isoform-level relative abundances. For a gene with isoform counts or
expression values x_i, convert to proportions

``` math
p_i = x_i / \sum_j x_j
```

so that $`p_i \ge 0`$ and $`\sum_i p_i = 1`$.

The parameter `q` controls how the entropy weights isoforms by
abundance: values $`q < 1`$ emphasize low-abundance (rare) isoforms and
therefore capture richness, while values $`q > 1`$ emphasize
high-abundance (dominant) isoforms and therefore capture dominance.

Practical choices: `q = 0` counts expressed isoforms (richness), `q = 1`
recovers Shannon entropy (overall uncertainty of isoform usage), and
`q = 2` is closely related to the inverse Simpson index (sensitive to
dominance). In practice, choose smaller `q` to highlight rare isoforms,
`q` near 1 for general diversity, and larger `q` to focus on dominant
isoforms.

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
readcounts <- as.matrix(tcga_brca_luma_dataset[, -1, drop = FALSE])
genes <- tcga_brca_luma_dataset$genes
ts_se <- calculate_diversity(readcounts, genes, q = 0.1, norm = TRUE)

# compute for multiple q values and plot the q-curve
qvec <- seq(0.01, 2, by = 0.1)
ts_multi <- calculate_diversity(readcounts, genes, q = qvec, norm = TRUE)
plot_tsallis_q_curve(ts_multi, group = colData(ts_multi)$sample_type)
```

For a detailed, reproducible example see the package vignette:
[HTML](https://gallardoalba.github.io/TSENAT/articles/TSENAT.html) \|
[PDF](https://gallardoalba.github.io/TSENAT/articles/TSENAT.pdf).

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
