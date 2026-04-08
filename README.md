[![CircleCI](https://circleci.com/gh/gallardoalba/TSENAT.svg?style=svg)](https://app.circleci.com/pipelines/github/gallardoalba/TSENAT) [![pkgdown](https://img.shields.io/badge/docs-pkgdown-blue.svg)](https://gallardoalba.github.io/TSENAT/) [![License: GPL-3](https://img.shields.io/badge/License-GPLv3-blue.svg)](LICENSE) ![GitHub last commit](https://img.shields.io/github/last-commit/gallardoalba/TSENAT) ![GitHub R package version](https://img.shields.io/github/r-package/v/gallardoalba/TSENAT) [![coverage](https://codecov.io/gh/gallardoalba/TSENAT/branch/stable/graph/badge.svg)](https://codecov.io/gh/gallardoalba/TSENAT/branch/stable)

# TSENAT: Tsallis Entropy Analysis Toolbox

TSENAT is a R package for quantifying and modeling **isoform-usage diversity** across RNA-seq samples using **Tsallis entropy** - a scale-dependent information-theoretic measure of transcript heterogeneity. 

## The Problem

Standard differential expression tools (DESeq2, edgeR) detect changes in total transcript abundance. However, genes often reorganize their isoform diversity *without* changing total abundance: they may shift from a balanced isoform distribution to dominance by a single isoform, or vice versa. This **isoform switching and splicing-driven regulation** is biologically important for cell state and function but invisible to abundance-focused methods.

## The Solution

TSENAT captures **isoform complexity** independently of which specific isoforms are abundant. The method uses **Tsallis entropy** with a sensitivity parameter `q` that acts like a lens:

- **Low q** (e.g., 0.5): Focuses on rare isoforms - detects if diversity is maintained or collapsed

- **Mid q** (e.g., 1.0): Balanced view (Shannon entropy) - overall isoform complexity

- **High q** (e.g., 2.0): Focuses on dominant isoforms - detects dominance shifts

By examining diversity across multiple q-values, you identify **scale-dependent** diversity changes - the hallmark of coordinate isoform switching.

## The Mathematics Behind Tsallis Entropy

Tsallis entropy is a parametric family of diversity measures that generalizes Shannon entropy and enables tuning sensitivity to different scales of isoform organization.

**Mathematical Definition**: For a discrete probability vector $p = (p_1, \ldots, p_n)$ representing isoform proportions within a gene, Tsallis entropy is:

$$S_q = \frac{1 - \sum_{i=1}^{n} p_i^q}{q - 1}$$

This elegant formula unifies diverse diversity concepts at specific q-values:

- **q = 0**: Richness — Simple count of expressed isoforms; emphasizes rare variants most strongly.
- **q = 1**: Shannon entropy — Standard information-theoretic measure; balanced weighting across scales.
- **q = 2**: Gini-Simpson index — Probability that two randomly-drawn transcripts are different; robust to rare variants.

By examining diversity across multiple q-values, you identify **scale-dependent** diversity changes—the hallmark of coordinate isoform switching. See **vignette("TSENAT")** for the complete mathematical treatment and information-theoretic interpretation.

### Divergence Analysis: Measuring Information-Theoretic Distance Between Conditions

While Tsallis entropy quantifies diversity *within* a single distribution, **Tsallis divergence** $D_q$ measures the information-theoretic distance *between* two distributions. 

**Mathematical Definition**: For two probability distributions $P$ and $Q$ representing isoform proportions in control and treatment conditions, Tsallis divergence is:

$$D_q(P||Q) = \frac{\sum_i p_i^q - \sum_i p_i \cdot q_i^{q-1}}{(q-1) \sum_i p_i}$$

Tsallis divergence enables the quantification of how fundamentally different the isoform complexity patterns are between experimental conditions.


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

# Quick Start

### Load Example Data

Start by loading the built-in example dataset from TSENAT, which includes transcript-level read counts, TPM values, and effective lengths from Salmon quantification. Then load the sample metadata and annotation file that describe your experimental design.

```r
suppressMessages({
  library(TSENAT)
  library(SummarizedExperiment)})

# Load example dataset (includes readcounts, tpm, and effective_length)
data(readcounts)
readcounts <- as.matrix(readcounts)

# Load sample metadata and annotation
metadata_df <- read.table(
  system.file("extdata", "metadata.tsv", package = "TSENAT"),
  header = TRUE,
  sep = "\t")

gff3_file <- system.file("extdata", "annotation.gff3.gz", package = "TSENAT")
```

### Create configuration and build analysis

Create a configuration object specifying your experimental design parameters (sample/condition columns from metadata) before building the analysis object. This fail-fast pattern ensures invalid parameters are caught immediately before processing begins.

```r
## Create configuration file
config <- tsenat_config(
  sample_col = "sample",
  condition_col = "condition",
  subject_col = "paired_samples",
  q_values = seq(0, 2, by = 0.05),
  nthreads = 2,
  paired = TRUE,
  control = "normal")

## Build TSENATAnalysis object
analysis <- build_analysis_s4(
  readcounts = readcounts, 
  tx2gene = gff3_file,
  metadata = metadata_df,
  config = config,
  tpm = tpm,
  effective_length = effective_length)
```

### Orchestration Function

The `tsenat()` function provides a complete, automated analysis pipeline in a single call. It takes your configured `TSENATAnalysis` object and executes all downstream analysis steps: entropy computation, statistical testing for q×condition interactions, and rich visualization. This is the recommended entry point for most users—it orchestrates the full workflow while respecting your configuration parameters (q-values, design, bootstrap settings, etc.) and handles output management seamlessly. For advanced customization, use individual functions directly as shown in the step-by-step workflow below.

```r
# Returns: Fully configured TSENATAnalysis object
result <- tsenat(analysis)
```

## Detailed Step-by-Step Workflow

For fine-grained control over your analysis, TSENAT also provides individual functions for each major step. This modular approach allows you to apply custom parameters at each stage, inspect intermediate results, or skip certain components entirely. The workflow below demonstrates the core analysis pipeline when using individual function calls—useful for exploratory analysis, parameter optimization, or integrating TSENAT results into larger custom workflows.

### 1. Filter & Compute Diversity

Remove low-abundance transcripts that may contribute noise to entropy calculations, then compute Tsallis entropy across your specified q-spectrum. This produces normalized diversity scores for each gene across all samples and q-values.

```r
# Remove low-abundance transcripts
analysis <- filter_analysis_s4(analysis, stringency = "medium")

# Compute Tsallis entropy across q-spectrum
analysis <- calculate_diversity_s4(analysis, norm = TRUE)
# Compute Jackknife
analysis <- jackknife_isoform_switching_s4(analysis,
    nboot = 100,
    lm_p_threshold = 0.05,
    threshold = 90)

# Compute divergence
analysis <- calcualte_divergence_s4(analysis)
```

### 2. Statistical Testing

Perform statistical testing to identify significant differences in entropy between experimental groups across your q-spectrum. TSENAT supports diverse statistical methods (linear models, rank-based tests, and more) to accommodate different study designs and data characteristics, ensuring you detect robust biological signals at the appropriate diversity scales.

```r
# Fit linear models to detect qxcondition interactions
analysis <- calculate_lm_interaction_s4(
  analysis,
  method = "gam")
```

### 3. Visualize Results

Create diverse visualizations to explore and communicate your analysis results. TSENAT provides multiple plotting functions including q-curves, volcano plots, heatmaps, and interaction plots to reveal scale-dependent diversity patterns and statistical findings across different aspects of your data.

```r
# Plot overall Tsallis q-spectrum for all genes 
p_qcurve <- plot_tsallis_q_curve_s4(analysis)
print(p_qcurve)

# Plot q-curve profiles for the top 4 genes
combined_plot <- plot_lm_interaction_gam_s4(
    analysis,
    n_top = 4)
print(combined_plot)
```


## Statistical Inference Methods

TSENAT provides a flexible statistical framework optimized for entropy-based diversity analysis.

### Parametric methods:

- **Generalized Additive Models (GAM) with ARIMA Differencing**: GAM is a semi-parametric approach using flexible smooth functions of *q* and group with ARIMA(1,1,0) first-differencing to remove monotone trend and achieve stationarity. Estimates fixed effect parameters via smooth basis functions; automatically detects heteroscedasticity and applies optimal variance weighting.

- **Linear Mixed Models (LMM) with AR(1)**: Fits mixed-effects models using `nlme::lme()` with random intercept by subject and AR(1) correlation structure for observation-level errors with explicit time-ordering by *q*-value.

- **Generalized Estimating Equations (GEE) with Multiple Correlation Structures**: Fits marginal models with three selectable correlation structures: AR(1) for *q*-ordered data, exchangeable for unordered measurements, or independence. Includes Kauermann-Carroll bias correction for small number of clusters.

- **Functional Principal Component Analysis (FPCA) with Regularization**: Treats entropy values as ordered curves across *q*-values; applies ARIMA differencing then performs PCA on the curve matrix to extract orthogonal smooth principal components.

- **Huber M-Estimation with Leave-One-Out Diagnostics**: Iteratively re-weighted least squares (IRLS) using Huber loss function.

### Non-parametric methods:

- **Friedman Rank-Based Test with Hochberg Correction**: Non-parametric alternative to repeated-measures ANOVA that operates entirely on ranks, requiring no distributional assumptions.

- **Jackknife via Delta Influence**: Leave-one-out resampling to identify transcript-level contributors to entropy changes via delta influence.

- **Wilcoxon Rank Sum Test with Multiple Hypothesis Correction**: Unpaired or paired (Wilcoxon signed-rank) non-parametric alternative to t-tests.

- **Permutation/Label Shuffling Test with FDR Control**: Non-parametric exact test via label shuffling.

## Data Integration

- **Unified object**: `TSENATAnalysis` encapsulates sequencing data, configuration, and all statistical results
- `SummarizedExperiment` foundation: Full Bioconductor ecosystem compatibility
- Accessor functions: `diversity()`, `divergence()`, `lmResults()`, `jisResults()`, etc.

## Related Packages

TSENAT addresses a fundamental but underappreciated question in transcriptomic analysis: **How do genes reorganize their isoform usage patterns, independent of changes in total abundance?** This question is distinct from standard differential expression analysis and reveals a layer of biological complexity—coordinated isoform switching—that conventional methods overlook. TSENAT fills a specific niche in the Bioconductor ecosystem by measuring scale-dependent isoform diversity rather than abundance or individual transcript shifts. Below is how TSENAT complements other Bioconductor tools:

| Tool | Answers | TSENAT Difference |
|------|---------|-------------------|
| **DESeq2, edgeR, limma** | Which genes change in *total abundance*? | TSENAT detects isoform diversity changes **independent of total abundance** |
| **DRIMSeq** | Which *individual transcripts* shift usage? | TSENAT measures overall isoform diversity, not individual transcript shifts |
| **IsoformSwitchAnalyzeR** | Which *individual isoforms* switch; what are the *functional consequences*? | TSENAT measures overall isoform diversity and diversity **shifts** rather than cataloging individual transcript switches or predicting functional consequences; complements switch identification with diversity patterns |
| **SplicingFactory** | What is the overall isoform diversity? | TSENAT extends with **scale-dependent diversity** (q-spectrum) vs fixed measures |
| **Kallisto, Salmon** | How many reads per transcript? | TSENAT uses their quantification as input; adds diversity analysis layer |

## Native Salmon Integration

TSENAT is specifically engineered to work seamlessly with Salmon quantification output. Rather than requiring manual parsing or format conversion, TSENAT automatically discovers transcript-level quantification files across your Salmon output directory and integrates them directly into the analysis pipeline. This tight integration means you can move from Salmon quantification to entropy analysis without intermediate data manipulation—the raw `quant.sf` files are all you need. TSENAT discovers these files automatically, validates their compatibility with your experimental design, and handles length-correction and normalization as part of the diversity computation workflow.

To get started with Salmon-quantified data:

```r
suppressMessages(library(TSENAT))

# Prepare configuration FIRST
config <- tsenat_config(
  q_values = seq(0, 2, by = 0.1),
  condition_col = "condition"
)

# Build analysis directly from Salmon output directory
# IMPORTANT: Use named parameters (salmon_dir=, tx2gene=)
analysis <- build_analysis_s4(
  salmon_dir = "path/to/salmon_output",  # Auto-discovers all quant.sf files
  tx2gene = "path/to/annotation.gff3.gz",  # Transcript-to-gene mapping
  metadata = metadata_df,  # Sample metadata (see structure below)
  config = config
)

# Run analysis pipeline
analysis <- filter_analysis_s4(analysis, stringency = "severe")
analysis <- calculate_diversity_s4(analysis)  # Salmon-informed length-normalized entropy
```

### Salmon Directory Structure

Your Salmon output must be organized with one folder per sample, each containing a quant.sf file with transcript-level quantification. Sample folder names must exactly match the row names in your metadata file (case-sensitive) to ensure correct sample attribution.

TSENAT expects Salmon output organized with one subdirectory per sample:

```
salmon_output/
├── Sample_1/
│   └── quant.sf
├── Sample_2/
│   └── quant.sf
├── Sample_3/
│   └── quant.sf
└── Sample_4/
    └── quant.sf
```

**Critical requirement**: Folder names (e.g., `Sample_1`, `Sample_2`) must **exactly match** the row names in your metadata file (case-sensitive).

## Metadata File Structure

Example metadata structure:

```r
# Load metadata from TSV file
metadata_df <- read.table("metadata.tsv", header = TRUE, sep = "\t")
```

Expected TSV file format (`metadata.tsv`):

```
sample        condition    paired_samples
SRR14800481   normal       A
SRR14800480   normal       B
SRR14800479   tumor        A
SRR14800478   tumor        B
```

**Key requirements:**
- First column: `sample` (must match Salmon folder names exactly)
- Second column: `condition` (experimental groups: normal, tumor, treated, control, etc.)
- Third column: `paired_samples` (required if using paired designs; identifier for matched samples)


## Tests coverage

Testing is vital in research as it ensures the validity and reliability of results, which is essential for accurately interpreting findings. The report about the current testing coverage can be found [here](https://app.codecov.io/gh/gallardoalba/TSENAT).

## Learn More

### Comprehensive Workflow

For a complete walkthrough of the analysis pipeline with real biological examples, see the main package vignette. This includes theory background, step-by-step explanations of each analysis function, and interpretation guidance for understanding your results.

See the package vignette for detailed examples, theory background, and typical workflows:

```r
vignette("TSENAT")
```

### Function Reference

Use R's built-in help system to explore detailed documentation for individual TSENAT functions and S4 classes. Each help page includes function arguments, return values, and practical examples of usage.

Interactive help for functions and classes:

```r
?build_analysis_s4
?calculate_diversity_s4
?TSENATAnalysis-class
```

For methodology details and a comprehensive bibliography, see the [TSENAT vignette](vignettes/TSENAT.Rmd):

```r
vignette("TSENAT")
```

## Citation

If you use TSENAT in your research, please cite:

```r
citation("TSENAT")
```

BibTeX entry:

```bibtex
@software{gallardo2026tsenat,
  title={TSENAT: Tsallis Entropy Analysis Toolbox},
  author={Gallardo Alba, Cristóbal},
  url={https://github.com/gallardoalba/TSENAT},
  year={2026}
}
```

## License and Attribution

This project is licensed under the GNU General Public License v3.0 (GPL-3). See [LICENSE](LICENSE) for details.

Attribution: TSENAT builds upon the [SplicingFactory package](https://github.com/esebesty/SplicingFactory), extending it with specialized focus on Tsallis entropy analysis.

> **“If I ever come back from the past, it's to create a cyclone.”**
>
> - Juan José Lozano
