TSENAT:  Tsallis Entropy Analysis Toolbox
=========================================

TSENAT analyze expression/transcript differences and compute diversity metrics.

Badges
------

- GitHub Actions (CI):

        [![R-CMD-check](https://github.com/gallardoalba/TSENAT/actions/workflows/check-standard.yaml/badge.svg)](https://github.com/gallardoalba/TSENAT/actions)

- Español README: [README.es.md](README.es.md)
- Deutsch README: [README.de.md](README.de.md)

Origin and attribution
-----------------------

TSENAT builds upon and adapts substantial portions of code from the
[SplicingFactory package](https://github.com/esebesty/SplicingFactory). The 
codebase has been extended with additional utilities, bug fixes and plotting 
helpers focused on Tsallis based transcript-level diversity analysis.

Tsallis theory
----------------------

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
                        ````markdown
                        TSENAT
                        =========

                        R package to detect and analyze expression/transcript differences and compute diversity metrics (includes Tsallis functions, statistical tests, and visualization helpers).

                        Origin and attribution
                        -----------------------

                        TSENAT builds upon and adapts substantial portions of code from the SplicingFactory project (credit to the original authors). The codebase has been extended with additional utilities, bug fixes and plotting helpers focused on transcript-level diversity analysis.

                        Tsallis theory (brief)
                        ----------------------

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

                        - Sample group inference and metadata helpers:
                                - `infer_sample_group` attempts to infer sample classes (e.g. Normal/Tumor or TCGA barcodes) from column names or from provided metadata when explicit group labels are absent.

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


                        Installation
                        ------------

                        Install from GitHub during development:

                        ```r
                        if (!requireNamespace("remotes", quietly = TRUE)) install.packages("remotes")
                        remotes::install_github("gallardoalba/TSENAT")
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

                        Vignette and documentation
                        --------------------------

                        View the package vignette with:

                        ```r
                        browseVignettes("TSENAT")
                        vignette("TSENAT")
                        ```

                        Tests and development
                        ---------------------

                        Run tests with `devtools::test()` or:

                        ```sh
                        R -e 'devtools::test()'
                        ```

                        Check the package:

                        ```sh
                        R CMD check --as-cran .
                        ```

                        Citation and license
                        --------------------

                        See `citation("TSENAT")` or [inst/CITATION](inst/CITATION). License in [LICENSE](LICENSE).

                        Contributing
                        ------------

                        Open an issue or submit a pull request on GitHub. Please run `devtools::check()` before submitting changes.

                        Publishing to GitHub (quick steps)
                        --------------------------------

                        1. Create remote repo and push:

                        ```sh
                        git init
                        git add .
                        git commit -m "Initial commit: add package source and README"
                        gh repo create YOUR_GITHUB_USERNAME/TSENAT --public --source=. --remote=origin
                        git push -u origin main
                        ```

                        2. Add CI using the included GitHub Actions workflow or call `usethis::use_github_action_check_standard()`.

                        Support
                        -------

                        Report issues on the GitHub repository.

                        — The TSENAT team
                        ````