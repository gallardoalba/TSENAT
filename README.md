[![CircleCI](https://circleci.com/gh/gallardoalba/TSENAT.svg?style=svg)](https://app.circleci.com/pipelines/github/gallardoalba/TSENAT) [![pkgdown](https://img.shields.io/badge/docs-pkgdown-blue.svg)](https://gallardoalba.github.io/TSENAT/) [![License: GPL-3](https://img.shields.io/badge/License-GPLv3-blue.svg)](LICENSE) ![GitHub last commit](https://img.shields.io/github/last-commit/gallardoalba/TSENAT) ![GitHub R package version](https://img.shields.io/github/r-package/v/gallardoalba/TSENAT) [![coverage](https://codecov.io/gh/gallardoalba/TSENAT/branch/stable/graph/badge.svg)](https://codecov.io/gh/gallardoalba/TSENAT/branch/stable)

# TSENAT: Tsallis Entropy Analysis Toolbox

TSENAT is a Bioconductor package for quantifying and modeling **isoform-usage diversity** across RNA-seq samples using **Tsallis entropy** - a scale-dependent information-theoretic measure of transcript heterogeneity. 

## The Problem

Standard differential expression tools (DESeq2, edgeR) detect changes in total transcript abundance. However, genes often reorganize their isoform diversity *without* changing total abundance: they may shift from a balanced isoform distribution to dominance by a single isoform, or vice versa. This **isoform switching and splicing-driven regulation** is biologically important for cell state and function but invisible to abundance-focused methods.

## The Solution

TSENAT captures **isoform complexity** independently of which specific isoforms are abundant. The method uses **Tsallis entropy** with a sensitivity parameter `q` that acts like a lens:

- **Low q** (e.g., 0.5): Focuses on rare isoforms - detects if diversity is maintained or collapsed

- **Mid q** (e.g., 1.0): Balanced view (Shannon entropy) - overall isoform complexity

- **High q** (e.g., 2.0): Focuses on dominant isoforms - detects dominance shifts

By examining diversity across multiple q-values, you identify **scale-dependent** diversity changes - the hallmark of coordinate isoform switching.

## The Mathematics Behind Tsallis Entropy

Tsallis entropy is defined as: $S_q = (1 - \sum p_i^q) / (q-1)$, where $p_i$ represents isoform proportions within a gene. This elegant equation generalizes Shannon entropy (which is recovered when $q \to 1$) and enables tuning sensitivity to different scales of isoform organization:

- **q = 0**: Richness (count of expressed isoforms)
- **q = 1**: Shannon entropy (balanced view across all abundance scales)
- **q = 2**: Gini-Simpson index (robust to rare variants, focuses on dominant isoforms)

This parametric family is the key innovation: by sliding q across scales, you zoom from rare isoform variants to dominant transcript patterns, capturing biological signal invisible to fixed-scale methods. See **vignette("TSENAT")** for the complete mathematical treatment and information-theoretic interpretation.

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
## Quick Start


# Quick Start

## Load Data & Configure

```r
library(TSENAT)
library(SummarizedExperiment)

# Load example dataset (includes readcounts, tpm, and effective_length)
data(readcounts)
readcounts <- as.matrix(readcounts)

# Load sample metadata and annotation
metadata_df <- read.table(
  system.file("extdata", "metadata.tsv",
  package = "TSENAT"), header = TRUE, row.names = 1,
  sep = "\t")

gff3_file <- system.file("extdata",
  "annotation.gff3.gz", package = "TSENAT")

## Configure analysis parameters first (best practice: fail-fast principle)
## This validates all parameters against metadata before object creation
config <- tsenat_config(
  condition_col = "condition",
  subject_col = "paired_samples",
  q_values = seq(0, 2, by = 0.05),
  nthreads = 2,
  paired = TRUE,
  control = "normal",
  metadata = metadata_df
)

## Build a complete `TSENATAnalysis` object from readcounts + GFF3.gz annotation
## Pass config at construction (Bioconductor pattern): immutable object creation
analysis <- build_analysis_s4(
  config = config,
  readcounts = readcounts, 
  tx2gene = gff3_file, 
  tpm = tpm,
  effective_length = effective_length
)
```

### Orchestration Function

For a complete analysis with default parameters, use the `tsenat()` orchestration function:

```r
# Returns: Fully configured TSENATAnalysis object with diversity, testing, and plots
analysis <- tsenat(se, config = cfg)
```

### Detailed Step-by-Step Workflow

For customization at each stage, use individual functions:

### 1. Filter & Compute Diversity

```r
# Remove low-abundance transcripts
analysis <- filter_analysis_s4(analysis, stringency = "medium")

# Compute Tsallis entropy across q-spectrum
analysis <- calculate_diversity_s4(analysis, norm = TRUE)
```

### 2. Statistical Testing

```r
# Fit linear models to detect qxcondition interactions
analysis <- calculate_lm_interaction_s4(
  analysis,
  method = "lmm")
```

### 3. Visualize Results

```r
# Plot overall q-curve
p_qcurve <- plot_tsallis_q_curve_s4(analysis)

print(p_qcurve)
```


## Statistical Inference Methods

- **Multiple testing methods**:
    - *Wilcoxon/Permutation*: Distribution-free testing for pairwise comparisons
    - *Linear Mixed Models (LMM)*: Parametric testing with AR(1) correlation structure for repeated measures; ideal when residuals are approximately normal
    - *GAM*: Non-parametric for detecting nonlinear scale-dependent q×condition interactions with adaptive smooth splines
    - *GEE*: Semi-parametric for clustered data; robust to variance misspecification with sandwich standard errors
    - *Friedman rank tests*: Maximal robustness for paired designs; ideal for bounded distributions like entropy
    - *M-estimation*: Outlier-resistant effect size calculations (Huber, Tukey weights)
- **Confidence intervals**:
    - BCA (bias-corrected and accelerated) bootstrap correction for skewed distributions like Tsallis entropy
    - Percentile bootstrap for symmetric distributions
    - Jackknife leave-one-out for identifying outlier-influential samples


### Data Integration
- **Unified object**: `TSENATAnalysis` encapsulates data, config, and all results
- `SummarizedExperiment` foundation: Full Bioconductor ecosystem compatibility
- Accessor functions: `diversity()`, `divergence()`, `lmResults()`, etc.

## Related Packages

TSENAT answers a unique question: **How do isoforms reorganize, independent of abundance changes?** It complements other Bioconductor tools:

| Tool | Answers | TSENAT Difference |
|------|---------|-------------------|
| **DESeq2, edgeR, limma** | Which genes change in *total abundance*? | TSENAT detects isoform diversity changes **independent of total abundance** |
| **DRIMSeq** | Which *individual transcripts* shift usage? | TSENAT measures overall isoform diversity, not individual transcript shifts |
| **SplicingFactory** | What is the overall isoform diversity? | TSENAT extends with **scale-dependent diversity** (q-spectrum) vs fixed measures |
| **Kallisto, Salmon** | How many reads per transcript? | TSENAT uses their quantification as input; adds diversity analysis layer |

## Loading Salmon Quantification Data

TSENAT automatically discovers and reads Salmon output when you provide a directory:

```r
library(TSENAT)

# Prepare configuration FIRST
config <- tsenat_config(
  q_values = seq(0, 2, by = 0.1),
  condition_col = "treatment"
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
metadata_df <- read.table("metadata.tsv", header = TRUE, sep = "\t", row.names = 1)
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

For methodology details and a comprehensive bibliography, see the [TSENAT vignette](vignettes/TSENAT.Rmd):

```r
vignette("TSENAT")
```

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

## License and Attribution

This project is licensed under the GNU General Public License v3.0 (GPL-3). See [LICENSE](LICENSE) for details.

Attribution: TSENAT builds upon the [SplicingFactory package](https://github.com/esebesty/SplicingFactory), extending it with specialized focus on Tsallis entropy analysis.

> **“If I ever come back from the past, it's to create a cyclone.”**
>
> - Juan José Lozano