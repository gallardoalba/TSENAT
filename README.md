[![CircleCI](https://circleci.com/gh/gallardoalba/TSENAT.svg?style=svg)](https://app.circleci.com/pipelines/github/gallardoalba/TSENAT) [![pkgdown](https://img.shields.io/badge/docs-pkgdown-blue.svg)](https://gallardoalba.github.io/TSENAT/) [![License: GPL-3](https://img.shields.io/badge/License-GPLv3-blue.svg)](LICENSE) ![GitHub last commit](https://img.shields.io/github/last-commit/gallardoalba/TSENAT) ![GitHub R package version](https://img.shields.io/github/r-package/v/gallardoalba/TSENAT)

# TSENAT: Tsallis Entropy Analysis Toolbox

Comprehensive R package for analyzing isoform diversity in transcript-level expression data using Tsallis entropy. TSENAT provides a unified framework for computing, testing, and visualizing scale-dependent isoform heterogeneity across biological samples.

## Overview

TSENAT analyzes expression and transcript differences to compute diversity metrics. Key capabilities:

- **Scale-dependent diversity analysis**: Evaluate isoform heterogeneity at different sensitivity levels using the parameter `q`
- **Statistical testing**: Compare diversity measures between sample groups using Wilcoxon tests or permutation-based approaches
- **Reproducible workflows**: From raw counts to publication-ready visualizations with paired sample support

## Tsallis Theory

Tsallis entropy generalizes Shannon entropy. For a probability vector $p = (p_1, \dots, p_n)$ (with $p_i \ge 0$ and $\sum_i p_i = 1$) the Tsallis entropy of order $q$ is defined for $q \ne 1$ as

$$
S_q(p) = \frac{1 - \sum_{i} p_i^q}{q - 1}.
$$

In the limit $q \to 1$ this recovers the Shannon entropy

$$
\lim_{q \to 1} S_q(p) = -\sum_i p_i \log p_i.
$$

### Application to Transcript Expression

In transcript-expression data we compute Tsallis entropy per gene from the isoform-level relative abundances. For a gene with isoform counts or expression values $x_i$, convert to proportions

$$
p_i = x_i / \sum_j x_j
$$

so that $p_i \ge 0$ and $\sum_i p_i = 1$. 

The parameter `q` controls how the entropy weights isoforms by abundance: 
- **$q < 1$**: emphasizes low-abundance (rare) isoforms; captures richness
- **$q \approx 1$**: Shannon entropy; overall uncertainty of isoform usage
- **$q > 1$**: emphasizes high-abundance (dominant) isoforms; captures dominance

## Features

### Tsallis Entropy and Diversity Calculations

- `calculate_tsallis_entropy()`: Computes $S_q$ and/or $D_q$ (Hill numbers) for numeric vectors
  - Supports normalization, multiple `q` values, and the $q \to 1$ limit
  - Returns numeric vectors or structured lists
  
- `calculate_diversity()`: Applies calculations across transcripts/genes
  - Works with matrices, `tximport`-style lists, or `SummarizedExperiment` objects
  - Returns `SummarizedExperiment` with assay `diversity` ($S_q$) or `hill` ($D_q$)

### Differential and Statistical Analyses

- `calculate_difference()`: Computes group means, differences, log₂-fold-changes, and p-values
- Support for paired and unpaired designs
- Wilcoxon rank-sum or permutation-based hypothesis tests
- Adjusted p-values and effect size reporting

### Plotting and Visualization

- `plot_tsallis_q_curve()`: q-curves showing median ± IQR across sample groups
- `plot_ma()`, `plot_volcano()`: Effect size and significance summaries
- `plot_top_transcripts()`: Transcript-level patterns for candidate genes
- `plot_tsallis_gene_profile()`: Per-gene q-curve profiles

## Installation and Documentation

Install from GitHub during development:

```r
if (!requireNamespace("remotes", quietly = TRUE)) install.packages("remotes")
remotes::install_github("gallardoalba/TSENAT")

## Recommended (reproducible environment): use renv
install.packages("renv")
renv::init()
renv::snapshot()
```

**Documentation**: [pkgdown site](https://gallardoalba.github.io/TSENAT/) | [Vignette](https://gallardoalba.github.io/TSENAT/articles/TSENAT.html)

## Quick Start

Compute Tsallis diversity for a single `q` and plot a q-curve across multiple `q` values:

```r
library(TSENAT)
data("tcga_brca_luma_dataset", package = "TSENAT")

# Compute normalized diversity for q = 0.1
readcounts <- as.matrix(tcga_brca_luma_dataset[, -1, drop = FALSE])
genes <- tcga_brca_luma_dataset$genes
ts_se <- calculate_diversity(readcounts, genes, q = 0.1, norm = TRUE)

# Compute across multiple q values
qvec <- seq(0.01, 2, by = 0.1)
ts_multi <- calculate_diversity(readcounts, genes, q = qvec, norm = TRUE)

# Visualize
p_qcurve <- plot_tsallis_q_curve(ts_multi)
print(p_qcurve)
```

For a detailed, reproducible workflow see the [package vignette](https://gallardoalba.github.io/TSENAT/articles/TSENAT.html).

## Citation

If you use TSENAT in your research, please cite:

```r
citation("TSENAT")
```

This command displays the recommended bibliographic entry. A machine-readable `CITATION` file is included with the package for easy export to reference managers.

**BibTeX entry:**
```bibtex
@software{gallardo2026tsenat,
  title={TSENAT: Tsallis Entropy Analysis Toolbox},
  author={Gallardo Alba, Cristóbal},
  url={https://github.com/gallardoalba/TSENAT},
  year={2026}
}
```

## Contributing

We welcome contributions! Please follow these guidelines:

1. **Fork** the repository
2. **Create a feature branch** (`git checkout -b feature/your-feature-name`)
3. **Make your changes** and test locally with `R CMD check`
4. **Commit with clear messages** (`git commit -m 'Add feature: description'`)
5. **Push to your fork** (`git push origin feature/your-feature-name`)
6. **Open a Pull Request** describing your changes

### Local Testing

Ensure all checks pass before submitting:

```r
devtools::check()
devtools::test()
```

## CI and Local Checks

Continuous integration is configured with CircleCI to install all suggested packages for comprehensive testing. To reproduce a CI-like environment locally:

```r
# Install all suggested dependencies
remotes::install_deps(dependencies = c("Suggests"))

# Run checks
R CMD check --as-cran
```

## License and Attribution

This project is licensed under the GNU General Public License v3.0 (GPL-3). See [LICENSE](LICENSE) for details.

**Attribution**: TSENAT builds upon and adapts code from the [SplicingFactory package](https://github.com/esebesty/SplicingFactory), extended with additional functionality, bug fixes, and visualization tools.

---

**Maintainer**: Cristóbal Gallardo (gallardoalba@pm.me)  
**Repository**: [github.com/gallardoalba/TSENAT](https://github.com/gallardoalba/TSENAT)
