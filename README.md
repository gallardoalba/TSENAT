[![CircleCI](https://circleci.com/gh/gallardoalba/TSENAT.svg?style=svg)](https://app.circleci.com/pipelines/github/gallardoalba/TSENAT) [![pkgdown](https://img.shields.io/badge/docs-pkgdown-blue.svg)](https://gallardoalba.github.io/TSENAT/) [![License: GPL-3](https://img.shields.io/badge/License-GPLv3-blue.svg)](LICENSE) ![GitHub last commit](https://img.shields.io/github/last-commit/gallardoalba/TSENAT) ![GitHub R package version](https://img.shields.io/github/r-package/v/gallardoalba/TSENAT) [![coverage](https://codecov.io/gh/gallardoalba/TSENAT/branch/stable/graph/badge.svg)](https://codecov.io/gh/gallardoalba/TSENAT/branch/stable)

# TSENAT: Tsallis Entropy Analysis Toolbox

TSENAT is a Bioconductor package for quantifying and modeling **isoform-usage diversity** across RNA-seq samples using **Tsallis entropy**—a scale-dependent information-theoretic measure of transcript heterogeneity. 

## The Problem

Standard differential expression tools (DESeq2, edgeR) detect changes in total transcript abundance. However, genes often reorganize their isoform diversity *without* changing total abundance: they may shift from a balanced isoform distribution to dominance by a single isoform, or vice versa. This **isoform switching and splicing-driven regulation** is biologically important for cell state and function but invisible to abundance-focused methods.

## The Solution

TSENAT captures **isoform complexity** independently of which specific isoforms are abundant. The method uses **Tsallis entropy** with a sensitivity parameter `q` that acts like a lens:
- **Low q** (e.g., 0.5): Focuses on rare isoforms—detects if diversity is maintained or collapsed
- **Mid q** (e.g., 1.0): Balanced view (Shannon entropy)—overall isoform complexity
- **High q** (e.g., 2.0): Focuses on dominant isoforms—detects dominance shifts

By examining diversity across multiple q-values, you identify **scale-dependent** diversity changes—the hallmark of coordinate isoform switching.

## Installation

**Requirements:** R >= 4.5.0

Install from [Bioconductor](https://bioconductor.org/packages/TSENAT) (recommended):

```r
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("TSENAT")
```

Or the development version from GitHub:

```r
remotes::install_github("gallardoalba/TSENAT")
```

## Workflow

TSENAT follows a streamlined pipeline:

```
Transcript Counts → Build Analysis → Filter → Configure → 
Compute Diversity → Test Differences → Visualize
```

### Quick Start: Orchestration Function

For a complete analysis with default parameters, use the `tsenat()` orchestration function:

```r
library(TSENAT)

# Build SummarizedExperiment from raw counts
se <- build_analysis_s4(
  readcounts = readcounts,
  tx2gene = gff3_file,
  metadata = metadata_df
)

# Configure and run complete pipeline
cfg <- tsenat_config(
  q_values = seq(0, 2, by = 0.1),
  condition_col = "treatment"
)

analysis <- tsenat(se, config = cfg)
# Returns: Fully configured TSENATAnalysis object with diversity, testing, and plots
```

### Detailed Step-by-Step Workflow

For customization at each stage, use individual functions:

### 1. Load Data & Configure

```r
library(TSENAT)

# Load transcript counts, annotation (GFF3), and sample metadata
analysis <- build_analysis_s4(
  readcounts = readcounts,
  tx2gene = gff3_file,
  metadata = metadata_df
)

# Configure analysis parameters once (used throughout pipeline)
analysis <- tsenat_config(
  analysis,
  q_values = seq(0, 2, by = 0.1),       # q-spectrum for scale analysis
  condition_col = "treatment",          # experimental groups
  subject_col = "patient",              # for paired designs
  control = "control"
)
```

### 2. Filter & Compute Diversity

```r
# Remove low-abundance transcripts
analysis <- filter_analysis_s4(analysis, stringency = "medium")

# Compute Tsallis entropy across q-spectrum
analysis <- calculate_diversity_s4(analysis, norm = TRUE)
```

### 3. Statistical Testing

```r
# Test for diversity differences between groups
analysis <- calculate_difference_s4(analysis, test = "wilcox")

# Fit linear models to detect q×condition interactions
analysis <- calculate_lm_interaction_s4(analysis, method = "gam")

# Identify isoform switching via jackknife diagnostics
analysis <- jackknife_isoform_switching_s4(analysis)
```

### 4. Visualize Results

```r
# Q-curve profile (how diversity changes across q-spectrum)
plot_tsallis_q_curve(analysis, gene = "your_gene")

# Volcano plot (significance vs. effect size)
plot_volcano_ma_grid_s4(analysis)

# Isoform switching heatmaps
plot_multiq_delta_influence_heatmaps_s4(analysis, n_genes = 4)
```

## Core Features

### Diversity Analysis
- **Multi-scale q-curves**: Examine isoform heterogeneity from rare to dominant isoforms
- **Tsallis entropy**: Scale-dependent complexity measure capturing beyond Shannon entropy
- **Normalization options**: Account for gene-specific isoform potential

### Statistical Inference  
- **Paired designs**: Account for repeated measures, subject random effects (via LMM/GEE)
- **Multiple testing methods**: Wilcoxon, permutation, linear models, GAM, robust M-estimation
- **Confidence intervals**: Bootstrap (percentile, BCA) and jackknife resampling

### Advanced Analysis
- **Divergence metrics**: Pairwise information-theoretic distance with effect sizes
- **Q×condition interactions**: Detect scale-dependent group differences via GAM/GEE
- **Isoform switching**: Jackknife-based diagnostics identifying transcript shifts
- **Robust methods**: M-estimation (Huber, Tukey) for outlier-resistant analysis

### Data Integration
- **Unified object**: `TSENATAnalysis` encapsulates data, config, and all results
- `SummarizedExperiment` foundation: Full Bioconductor ecosystem compatibility
- Accessor functions: `diversity()`, `divergence()`, `lmResults()`, etc.

## Related Packages

TSENAT complements other Bioconductor RNA-seq analysis tools:

- **DESeq2, edgeR, limma**: These tools detect *abundance* differences. TSENAT detects **isoform diversity** changes independent of total abundance.
- **SplicingFactory, DRIMSeq**: These tools test for shifts in *individual transcript proportions* (differential transcript usage). TSENAT quantifies **overall isoform heterogeneity/complexity** as a unified measure.
- **Bioconductor standard**: All tools use `SummarizedExperiment`, enabling seamless integration in multi-tool workflows.

## Example Data

TSENAT includes a reproducible example dataset with transcript counts and sample metadata. Load it with:

```r
data("readcounts", package = "TSENAT")
meta_file <- system.file("extdata", "metadata.tsv", package = "TSENAT")
gff3_file <- system.file("extdata", "annotation.gff3.gz", package = "TSENAT")
```

## Loading Salmon Quantification Data

TSENAT supports raw **Salmon quantification results** for length-normalized diversity analysis. This is particularly important for accurate entropy calculations across transcripts of different lengths.

### From Salmon Output Directory

If you have Salmon quantification output (one `quant.sf` file per sample):

```r
library(TSENAT)

# 1. Aggregate Salmon quant.sf files into a counts matrix
readcounts <- tximport::tximport(
  files = list.files("salmon_output", pattern = "quant.sf$", full.names = TRUE),
  type = "salmon",
  txOut = TRUE,  # Keep transcript-level (not gene-aggregated)
  ignoreTxVersion = TRUE
)

# Extract count data (NumReads column)
counts_matrix <- readcounts$counts

# 2. Get TPM and effective length from Salmon
tpm_matrix <- readcounts$abundance  # TPM values
eff_length <- rowMeans(readcounts$length)  # Median effective length per transcript

# 3. Load transcript annotation (GFF3) and sample metadata
gff3_file <- "path/to/annotation.gff3.gz"  # or .gff3
metadata_df <- read.table("metadata.tsv", header = TRUE, sep = "\t")

# 4. Build analysis with Salmon data
analysis <- build_analysis_s4(
  readcounts = counts_matrix,
  tx2gene = gff3_file,
  metadata = metadata_df,
  tpm = tpm_matrix,                # Salmon TPM for filtering
  effective_length = eff_length     # Salmon effective length for normalization
)

# Configure and run analysis
cfg <- tsenat_config(
  q_values = seq(0, 2, by = 0.1),
  condition_col = "treatment"
)
analysis <- setConfig(analysis, cfg)
analysis <- filter_analysis_s4(analysis, stringency = "medium")
analysis <- calculate_diversity_s4(analysis)  # TPM-informed filtering + length-normalized entropy
```

### Key Parameters for Salmon Data:

- **`tpm`**: Transcripts Per Million from Salmon (`abundance` column from `tximport`)
  - Used for TPM-based filtering in `filter_analysis_s4()` to remove very low abundance transcripts
  - Optional but recommended for improved filtering accuracy

- **`effective_length`**: Effective transcript length from Salmon (`length` column from `tximport`)
  - Used for length-normalization in entropy calculations
  - Accounts for transcript GC-bias and variable sequencing depth
  - Optional but recommended for cross-study comparability

See `build_analysis_s4()` documentation for details on data format requirements.

## Tests coverage

Testing is vital in research as it ensures the validity and reliability of results, which is essential for accurately interpreting findings. The report about the current testing coverage can be found [here](https://app.codecov.io/gh/gallardoalba/TSENAT).

## Learn More

### Comprehensive Workflow
See the package vignette for detailed examples, theory background, and typical workflows:

```r
vignette("TSENAT")
```

### Function Reference
Interactive help for functions and classes:

```r
?build_analysis_s4
?calculate_diversity_s4
?TSENATAnalysis-class
```

### Online Documentation
Full documentation and examples: [gallardoalba.github.io/TSENAT](https://gallardoalba.github.io/TSENAT)

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

1. Fork the repository
2. Create a feature branch (`git checkout -b feature/your-feature-name`)
3. Make your changes and test locally with `R CMD check`
4. Commit with clear messages (`git commit -m 'Add feature: description'`)
5. Push to your fork (`git push origin feature/your-feature-name`)
6. Open a Pull Request describing your changes

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
## Learn More

For methodology details and a comprehensive bibliography with 50+ peer-reviewed citations, see the [TSENAT vignette](vignettes/TSENAT.Rmd):

```r
vignette("TSENAT")
```

## License and Attribution


This project is licensed under the GNU General Public License v3.0 (GPL-3). See [LICENSE](LICENSE) for details.

Attribution: TSENAT builds upon the [SplicingFactory package](https://github.com/esebesty/SplicingFactory), extending it with specialized focus on Tsallis entropy analysis.

> **“If I ever come back from the past, it's to create a cyclone.”**
>
> — Juan José Lozano