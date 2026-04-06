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

## Workflow

TSENAT follows a streamlined pipeline:

```
Transcript Counts -> Build Analysis -> Filter -> Configure -> Compute Diversity -> Test Differences -> Visualize
```

### Quick Start: Orchestration Function

For a complete analysis with default parameters, use the `tsenat()` orchestration function:

```r
library(TSENAT)

# Build SummarizedExperiment from raw counts
se <- build_analysis_s4(
  readcounts = readcounts,
  tx2gene = gff3_file,
  metadata = metadata_df,
  tpm = tpm_matrix,
  effective_length = tx_lengths
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

**⚠️ Important**: Always use **named parameters** when calling `build_analysis_s4()`. The optional `salmon_dir` parameter comes before the required `tx2gene` parameter, so positional arguments may be misinterpreted. Use `tx2gene = ` and `salmon_dir = ` explicitly.

### 1. Load Data & Configure

```r
library(TSENAT)

# Load transcript counts, annotation (GFF3), and sample metadata
# CORRECT: Use named parameters
analysis <- build_analysis_s4(
  readcounts = readcounts,
  tx2gene = gff3_file,
  metadata = metadata_df,
  tpm = tpm_matrix,
  effective_length = tx_lengths
)

# WRONG: Do not use positional arguments
# analysis <- build_analysis_s4(readcounts, gff3_file, metadata = metadata_df)
# The above would fail because gff3_file is interpreted as salmon_dir

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

# Fit linear models to detect qxcondition interactions
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

### Statistical Inference  
- **Paired designs**: Account for repeated measures with subject random effects (via LMM)
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

### Confidence Intervals & Effect Sizes
- **Bootstrap confidence intervals**: Automatic BCA correction for asymmetric entropy distributions
- **Effect size interpretation**: Standardized measures enabling cross-study comparison
- **Jackknife diagnostics**: Identify which samples drive isoform switching signals

### Advanced Analysis
- **Divergence metrics**: Pairwise Kullback-Leibler and Jensen-Shannon divergence with effect sizes
- **Qxcondition interactions**: Detect scale-dependent group differences via GAM (smooth nonlinear patterns) or Friedman rank tests (maximal robustness)
- **Isoform switching**: Jackknife-based diagnostics identifying transcript shifts and influence plots
- **Robust methods**: M-estimation (Huber, Tukey) for outlier-resistant analysis

### Data Integration
- **Unified object**: `TSENATAnalysis` encapsulates data, config, and all results
- `SummarizedExperiment` foundation: Full Bioconductor ecosystem compatibility
- Accessor functions: `diversity()`, `divergence()`, `lmResults()`, etc.

## Related Packages

TSENAT answers a unique question: **How do isoforms reorganize, independent of abundance changes?** It complements other Bioconductor tools:

| Tool | Answers | TSENAT Difference |
|------|---------|-------------------|
| **DESeq2, edgeR, limma** | Which genes change in *total abundance*? | TSENAT detects isoform diversity changes **independent of total abundance** |
| **DRIMSeq, SplicingFactory** | Which *individual transcripts* shift usage? | TSENAT quantifies **overall heterogeneity** as a unified complexity measure |
| **Kallisto, Salmon** | How many reads per transcript? | TSENAT uses their quantification as input; adds diversity analysis layer |
| **All Bioconductor tools** | Various RNA-seq questions | TSENAT integrates via `SummarizedExperiment` for seamless multi-tool workflows |

## Loading Salmon Quantification Data

TSENAT automatically discovers and reads Salmon output when you provide a directory:

**Note**: When using `salmon_dir`, the parameter must be explicitly named (not positional) because `salmon_dir` is optional and comes before the required `tx2gene` parameter.

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
analysis <- filter_analysis_s4(analysis, stringency = "medium")
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

### Metadata File Structure

The metadata file must be a data frame with:
- **Row names**: Sample identifiers that exactly match Salmon folder names
- **Columns**: Experimental factors and sample information

Example metadata structure:

```r
# Load metadata from TSV file
metadata_df <- read.table("metadata.tsv", header = TRUE, sep = "\t", row.names = 1)

# Or create manually:
metadata_df <- data.frame(
  treatment = c("control", "control", "treated", "treated"),
  batch = c("batch1", "batch1", "batch2", "batch2"),
  patient_id = c("P001", "P002", "P001", "P002"),
  row.names = c("Sample_1", "Sample_2", "Sample_3", "Sample_4")  # Must match folder names!
)
```

Expected TSV file format (`metadata.tsv`):

```
sample_id    treatment    batch       patient_id
Sample_1     control      batch1      P001
Sample_2     control      batch1      P002
Sample_3     treated      batch2      P001
Sample_4     treated      batch2      P002
```

When reading from TSV:
```r
metadata_df <- read.table("metadata.tsv", header = TRUE, sep = "\t", row.names = 1)
# row.names = 1 uses first column (sample IDs) as row names
```

### Key Design Principles

- **Configure first**: Use `tsenat_config()` BEFORE `build_analysis_s4()` for fail-fast validation
- **Folder-metadata matching**: Sample folder names MUST exactly match metadata row names (case-sensitive); no auto-mapping
- **Row names in metadata**: Use `read.table(..., row.names = 1)` when reading TSV to set sample identifiers as row names
- **Salmon auto-discovery**: Function discovers all `quant.sf` files recursively; ensure one per sample subdirectory
- **Condition column**: Specify the metadata column containing experimental groups in `tsenat_config(condition_col = "...")`

See `build_analysis_s4()` documentation for complete parameter details.

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

### Complementary Validation Methods

For users interested in validating results across statistical frameworks:

```r
vignette("TSENAT_appendix_B")  # Compares linear models vs GAM vs Friedman rank-based tests
```

Appendix B demonstrates that discoveries generalize across non-parametric alternatives, providing critical validation that findings are robust to modeling assumptions.

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
> - Juan José Lozano