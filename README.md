[![R-CMD-check](https://github.com/gallardoalba/TSENAT/actions/workflows/check-standard.yaml/badge.svg)](https://github.com/gallardoalba/TSENAT/actions) [![pkgdown](https://github.com/<OWNER>/<REPO>/actions/workflows/pkgdown.yaml/badge.svg)](https://<OWNER>.github.io/<REPO>/)  [![License: GPL-3](https://img.shields.io/badge/License-GPLv3-blue.svg)](LICENSE)

TSENAT:  Tsallis Entropy Analysis Toolbox
=========================================

TSENAT analyze expression/transcript differences and compute diversity metrics.


Origin and attribution
-----------------------

TSENAT builds upon and adapts substantial portions of code from the
[SplicingFactory package](https://github.com/esebesty/SplicingFactory). The
codebase has been extended with additional utilities, bug fixes and plotting
helpers focused on Tsallis based transcript-level diversity analysis.

Tsallis theory
---------------

Tsallis entropy generalizes Shannon entropy and is defined as
S_q = (1 - sum p^q) / (q - 1) for a probability vector p. In the limit
q -> 1 it recovers Shannon entropy. In TSENAT we compute Tsallis entropy
and the related Hill numbers (D_q) per gene to measure isoform diversity.
The parameter `q` tunes sensitivity to rare vs abundant isoforms (q < 1
emphasizes rare isoforms; q > 1 emphasizes abundant ones).

Integrated features
-------------------

- Tsallis entropy and diversity calculations:
    - `calculate_tsallis_entropy`: computes S_q and/or D_q for a numeric vector of expression values (supports normalization, multiple `q` values and the q→1 limit). Returns numeric vectors or a list depending on `what`.
    - `calculate_diversity`: applies the calculation across transcripts/genes for matrices, `tximport`-style lists or `SummarizedExperiment` objects and returns a `SummarizedExperiment` with assay `diversity` (S_q) or `hill` (D_q).

- Differential and statistical analyses:
    - `calculate_difference` and helpers in `difference_functions` compute group means, differences (or log2 fold-changes), p-values and adjusted p-values. These functions are designed to work with diversity summaries as well as expression matrices.

- Method wrappers and utilities:
    - `calculate_method` provides a wrapper to run the chosen diversity method per gene, format outputs and evaluate multiple `q` values in one pass.

- Plotting and visualization:
    - `plot_tsallis_q_curve`: median ± IQR of Tsallis entropy across q-values by group.
    - `plot_tsallis_violin_multq`: violin plots of Tsallis entropy for multiple q-values and groups.
    - `plot_diversity_density`: density plots of diversity by sample type.
    - `plot_mean_violin`: violin plot of per-gene mean diversity by sample type.
    - `plot_ma`: MA-plot for differential results.
    - `plot_volcano`: volcano plot with labeling of top genes.
    - `plot_top_transcripts`: summary/visualization for top transcripts.

- Internal helpers:
    - Utilities such as `prepare_tsallis_long` format results for `ggplot2` and work with `tidyr`/`dplyr` to produce ready-to-plot data frames.

- Example dataset: `tcga_brca_luma_dataset`, included for vignette examples and tests.

Installation and documentation
------------

Install from GitHub during development:

```r
if (!requireNamespace("remotes", quietly = TRUE)) install.packages("remotes")
remotes::install_github("gallardoalba/TSENAT")

## Development / Reproducible setup
## Recommended (reproducible): use `renv`.

install.packages("renv")
renv::init()
# when dependencies are set, record lock
renv::snapshot()
```

Quick start
-----------

```r
library(TSENAT)
data("tcga_brca_luma_dataset", package = "TSENAT")

res <- calculate_difference(tcga_brca_luma_dataset$counts,
                            group = tcga_brca_luma_dataset$group)
head(res)
```

The complete documentation can be found at [inst/doc/TSENAT.html](inst/doc/TSENAT.html).

Citation and license
--------------------

This project is licensed under the GNU General Public License v3.0 (GPL-3).
A copy of the full license text is included in the repository at `LICENSE`
and is installed with the package under `inst/LICENSE`.

See `citation("TSENAT")` or [inst/CITATION](inst/CITATION). License in [LICENSE](LICENSE).