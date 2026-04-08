# TSENAT: Tsallis Entropy Analysis Toolbox

TSENAT is a R package for quantifying and modeling **isoform-usage
diversity** across RNA-seq samples using **Tsallis entropy** - a
scale-dependent information-theoretic measure of transcript
heterogeneity.

## The Problem

Standard differential expression tools (DESeq2, edgeR) detect changes in
total transcript abundance. However, genes often reorganize their
isoform diversity *without* changing total abundance: they may shift
from a balanced isoform distribution to dominance by a single isoform,
or vice versa. This **isoform switching and splicing-driven regulation**
is biologically important for cell state and function but invisible to
abundance-focused methods.

## The Solution

TSENAT captures **isoform complexity** independently of which specific
isoforms are abundant. The method uses **Tsallis entropy** with a
sensitivity parameter `q` that acts like a lens:

- **Low q** (e.g., 0.5): Focuses on rare isoforms - detects if diversity
  is maintained or collapsed

- **Mid q** (e.g., 1.0): Balanced view (Shannon entropy) - overall
  isoform complexity

- **High q** (e.g., 2.0): Focuses on dominant isoforms - detects
  dominance shifts

By examining diversity across multiple q-values, you identify
**scale-dependent** diversity changes - the hallmark of coordinate
isoform switching.

## The Mathematics Behind Tsallis Entropy

Tsallis entropy is defined as: $`S_q = (1 - \sum p_i^q) / (q-1)`$, where
$`p_i`$ represents isoform proportions within a gene. This elegant
equation generalizes Shannon entropy (which is recovered when
$`q \to 1`$) and enables tuning sensitivity to different scales of
isoform organization:

- **q = 0**: Richness
- **q = 1**: Shannon entropy
- **q = 2**: Gini-Simpson index

This parametric family is the key innovation: by sliding q across
scales, you zoom from rare isoform variants to dominant transcript
patterns, capturing biological signal invisible to fixed-scale methods.
See **vignette(“TSENAT”)** for the complete mathematical treatment and
information-theoretic interpretation.

### Divergence Analysis: Measuring Information-Theoretic Distance Between Conditions

While Tsallis entropy quantifies diversity *within* a single
distribution, **Tsallis divergence** $`D_q`$ measures the
information-theoretic distance *between* two distributions. This enables
quantification of how fundamentally different the isoform complexity
patterns are between experimental conditions, automatically accounting
for scale-dependent effects.

**Mathematical Definition**: For two probability distributions $`P`$ and
$`Q`$ representing isoform proportions in control and treatment
conditions, Tsallis divergence is:

``` math
D_q(P||Q) = \frac{\sum_i p_i^q - \sum_i p_i \cdot q_i^{q-1}}{(q-1) \sum_i p_i}
```

This measure unifies several well-known divergence concepts: - **q =
1**: Recovers Kullback-Leibler divergence (relative entropy) - **q =
0.5**: Emphasizes rare isoforms, sensitive to minority variants - **q =
2**: Emphasizes dominant isoforms, robust to rare variants

**Key Properties**: - **Scale-dependent sensitivity**: Different q
values reveal whether diversity shifts occur in rare (low q) or abundant
(high q) fractions of the transcriptome - **Non-symmetry**:
$`D_q(P||Q) \ne D_q(Q||P)`$, reflecting the directional nature of
information comparison (important for paired designs) - **Effect size
interpretation**: Values $`D > 0.1`$ indicate meaningful biological
separation between conditions; values near 0 suggest similar isoform
complexity patterns

**Application in TSENAT**: For paired study designs, TSENAT computes
divergence separately for each pair, then averages to create a robust,
paired-design-aware effect size. This respects within-pair correlation
while accounting for between-pair variation. The multi-q divergence
profile reveals whether group differences are concentrated at specific
diversity scales (indicating mechanism-specific isoform shifts) or
distributed uniformly (indicating broad-spectrum reorganization). See
[`calculate_divergence_s4()`](https://gallardoalba.github.io/TSENAT/reference/calculate_divergence_s4.md)
and
[`effect_sizes_divergence_s4()`](https://gallardoalba.github.io/TSENAT/reference/effect_sizes_divergence_s4.md)
for implementation details.

### Jackknife Isoform Switching: Identifying Robust Transcript-Level Contributors

Beyond group-level diversity statistics, researchers often need to
identify *which individual transcripts* drive observed isoform
complexity changes. TSENAT uses **jackknife leave-one-out resampling**
combined with **delta influence** weighting to robustly identify
transcript-level switching with stability assessment.

**Jackknife Delta Influence**: For each transcript $`i`$ in a gene, the
delta influence quantifies how that transcript’s relative contribution
to isoform diversity changes between conditions:

``` math
\Delta I_i = \text{Mean influence in condition 1} - \text{Mean influence in condition 2}
```

where influence is computed across bootstrap replicates. Values are: -
**Positive**: Transcript more influential (abundant, stable) in first
condition - **Negative**: Transcript more influential in second
condition - **Magnitude**: Larger absolute values indicate more robust,
consistent switching across replicates

**Robustness Weighting**: To distinguish signal from noise, TSENAT
weights delta influence by the **support frequency** across bootstrap
iterations. A transcript switching identified in 95% of bootstrap
samples receives higher confidence than one identified in only 60%, even
if both have similar magnitude. This operationalizes the principle that
robust signals persist across resampling, while artifacts disappear.

**Scale-Dependent Switching Patterns**: By computing jackknife delta
influence at each q-value, TSENAT reveals whether: - **Consistent
switching**: The same transcripts dominate (high delta influence) across
all q-values, indicating scale-independent isoform shifts -
**Scale-dependent switching**: Different transcripts show high delta
influence at different q-values (e.g., different transcripts drive
changes in rare vs. abundant fractions), indicating complex regulatory
mechanisms

**Key Parameters**: - `nboot`: Number of bootstrap replicates (default
100-1000; higher values increase stability assessment precision) -
`threshold`: Minimum support frequency to classify a transcript as
“robustly switching” (default 90%; e.g., must appear in ≥90% of
bootstrap samples) - `lm_p_threshold`: Pre-filter genes before jackknife
analysis (only test genes with significant q×condition interaction, p \<
0.05)

**Implementation**: See
[`jackknife_isoform_switching_s4()`](https://gallardoalba.github.io/TSENAT/reference/jackknife_isoform_switching_s4.md)
for details on compute;
[`plot_multiq_delta_influence_heatmaps_s4()`](https://gallardoalba.github.io/TSENAT/reference/plot_multiq_delta_influence_heatmaps_s4.md)
for visualization; and
[`prepare_gene_switching_tables_s4()`](https://gallardoalba.github.io/TSENAT/reference/prepare_gene_switching_tables_s4.md)
for extracting results as publication-ready tables.

**Biological Interpretation**: Genes with robust jackknife switching
signals (high support, large magnitude) represent high-confidence
isoform reorganization events. These are candidate targets for
functional validation, because they represent coordinate, reproducible
transcript usage changes that are unlikely to be driven by experimental
noise.

## Installation

**Requirements:** R \>= 4.5.0

Install from [Bioconductor](https://bioconductor.org/packages/TSENAT)
(recommended):

``` r

if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("TSENAT")
```

Or the development version from GitHub:

``` r

remotes::install_github("gallardoalba/TSENAT")
```

# Quick Start

### Load Example Data

Start by loading the built-in example dataset from TSENAT, which
includes transcript-level read counts, TPM values, and effective lengths
from Salmon quantification. Then load the sample metadata and annotation
file that describe your experimental design.

``` r

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

Create a configuration object specifying your experimental design
parameters (sample/condition columns from metadata) before building the
analysis object. This fail-fast pattern ensures invalid parameters are
caught immediately before processing begins.

``` r

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

The
[`tsenat()`](https://gallardoalba.github.io/TSENAT/reference/tsenat.md)
function provides a complete, automated analysis pipeline in a single
call. It takes your configured `TSENATAnalysis` object and executes all
downstream analysis steps: entropy computation, statistical testing for
q×condition interactions, and rich visualization. This is the
recommended entry point for most users—it orchestrates the full workflow
while respecting your configuration parameters (q-values, design,
bootstrap settings, etc.) and handles output management seamlessly. For
advanced customization, use individual functions directly as shown in
the step-by-step workflow below.

``` r

# Returns: Fully configured TSENATAnalysis object
result <- tsenat(analysis)
```

## Detailed Step-by-Step Workflow

For fine-grained control over your analysis, TSENAT also provides
individual functions for each major step. This modular approach allows
you to apply custom parameters at each stage, inspect intermediate
results, or skip certain components entirely. The workflow below
demonstrates the core analysis pipeline when using individual function
calls—useful for exploratory analysis, parameter optimization, or
integrating TSENAT results into larger custom workflows.

### 1. Filter & Compute Diversity

Remove low-abundance transcripts that may contribute noise to entropy
calculations, then compute Tsallis entropy across your specified
q-spectrum. This produces normalized diversity scores for each gene
across all samples and q-values.

``` r

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

Perform statistical testing to identify significant differences in
entropy between experimental groups across your q-spectrum. TSENAT
supports diverse statistical methods (linear models, rank-based tests,
and more) to accommodate different study designs and data
characteristics, ensuring you detect robust biological signals at the
appropriate diversity scales.

``` r

# Fit linear models to detect qxcondition interactions
analysis <- calculate_lm_interaction_s4(
  analysis,
  method = "gam")
```

### 3. Visualize Results

Create diverse visualizations to explore and communicate your analysis
results. TSENAT provides multiple plotting functions including q-curves,
volcano plots, heatmaps, and interaction plots to reveal scale-dependent
diversity patterns and statistical findings across different aspects of
your data.

``` r

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

TSENAT provides a flexible statistical framework optimized for
entropy-based diversity analysis. The default configuration uses
Generalized Additive Models (GAM) with Benjamini-Hochberg multiple
testing correction and Friedman rank tests for paired designs. Users can
adjust the statistical methodology via the `lm_method`, `lm_pcorr`, and
`bootstrap_method` configuration parameters to suit their study design
and data characteristics.

### Linear Modeling Approaches

TSENAT supports four distinct parametric/semi-parametric linear modeling
frameworks, selectable via `lm_method`:

- **Generalized Additive Models (GAM/GAMM)** \[DEFAULT\]: Flexible
  smoothing with locally-weighted basis functions. For unpaired designs,
  uses standard GAM with F-tests. For paired/repeated measures,
  automatically switches to GAMM with AR(1) autocorrelation structure.
  Automatically selects appropriate family (Beta for bounded \[0,1\]
  entropy values, Gamma for heteroscedastic data, Gaussian otherwise)
  with smoothing bias correction for small sample sizes.

- **Linear Mixed Models (LMM)**: Parametric framework with AR(1)
  correlation structure for repeated measures designs. Appropriate when
  study design assumptions align with mixed model requirements and
  parametric inference is preferred. Supports multiple p-value
  computation methods (Satterthwaite, Likelihood Ratio Test, or both).

- **Functional Principal Component Analysis (FPCA)**: Treats entropy
  curves across the q-spectrum as functional objects, extracting
  orthogonal functional principal components. Performs ANOVA on
  component scores. Ideal for investigations emphasizing scale-dependent
  (multi-q) diversity patterns.

- **Generalized Estimating Equations (GEE)**: Semi-parametric approach
  using working correlation structures for correlated non-normal data.
  Provides robust alternative to mixed models when parametric
  distributional assumptions are uncertain.

**Method Selection Guidance**: GAM is recommended for entropy analysis
due to its flexibility with bounded, non-normal distributions. Choose
LMM if parametric assumptions are well-justified. Choose FPCA for
multi-scale comparative analysis across the q-spectrum. Choose GEE for
robustness to model misspecification. See
[`vignette("TSENAT")`](https://gallardoalba.github.io/TSENAT/articles/TSENAT.md)
for detailed method comparisons and examples.

### Rank-Based Tests for Paired Designs

For paired study designs with non-normal data distributions:

- **Friedman Rank Test** \[DEFAULT for paired\]: Non-parametric test
  comparing entropy distributions across conditions within matched
  pairs. Provides maximal robustness for bounded, non-normal
  distributions characteristic of entropy metrics.

- **Kendall Test**: Non-parametric alternative with different power
  characteristics; useful as sensitivity analysis.

- **Conditioned Rank Tests**: Stratified rank testing within
  experimental strata for complex designs.

### Distribution-Free and Permutation Methods

- **Wilcoxon Signed-Rank Test**: For paired pairwise comparisons of
  entropy between samples or conditions.

- **Permutation Testing**: Distribution-free testing via sample
  permutation, implemented for pairwise comparisons.

- **Westfall-Young Permutation Correction**: Advanced permutation-based
  multiple testing procedure controlling family-wise error rate with
  stepdown, more powerful than traditional Bonferroni correction. Useful
  for comprehensive hypothesis testing scenarios.

### Multiple Testing Corrections

- **Benjamini-Hochberg (BH)** \[DEFAULT\]: Controls False Discovery Rate
  (FDR) across multiple hypothesis tests. Appropriate for exploratory
  analysis where false positive control is balanced against power.

- **Bonferroni**: Conservative family-wise error rate (FWER) control;
  use when strict false positive protection is required.

- **Holm Step-Down**: Less conservative FWER control; compromise between
  Bonferroni and BH.

- **Westfall-Young Permutation**: Permutation-based stepdown FWER
  procedure (see Permutation Methods).

### Bootstrap Confidence Intervals

TSENAT incorporates bootstrap resampling for robust uncertainty
quantification:

- **Percentile Bootstrap** \[DEFAULT\]: Computes confidence intervals
  from empirical quantiles of bootstrap distribution. Fast and
  assumption-free but assumes symmetric sampling distribution.

- **Bias-Corrected and Accelerated (BCA) Bootstrap** \[RECOMMENDED for
  entropy\]: Adjusts for distribution skewness and bias through
  acceleration factor computation. Ideal for bounded, skewed entropy
  distributions. Provides diagnostic output (skewness, acceleration) to
  assess bootstrap assumption satisfaction.

**Bootstrap Recommendation for Entropy Data**: Use
`bootstrap_method = "bca"` (set via
[`tsenat_config()`](https://gallardoalba.github.io/TSENAT/reference/tsenat_config.md))
for entropy analysis, as entropy metrics are inherently bounded \[0, log
N\] and often exhibit positive skewness. The BCA correction is
specifically designed to handle this scenario.

### Transcript-Level Switching and Influence Assessment

- **Jackknife Isoform Switching Analysis**: Bootstrap-based
  leave-one-transcript-out resampling identifying which individual
  transcripts drive observed isoform complexity changes. Computes delta
  influence (transcript contribution shift between conditions) weighted
  by support frequency across bootstrap replicates. Identifies “robust
  switching”—transcripts showing consistent switching across ≥90% of
  bootstrap samples (configurable threshold). Complements group-level
  diversity tests by pinpointing mechanism.

- **M-Estimation Robustness Weighting** (Tukey biweight, Huber weights):
  Provides outlier-resistant influence assessment and effect size
  calculations for divergence analysis. Identifies samples with
  disproportionate effect on results during quality control stage
  (incorporated as step 4/14 of automated workflow).

### Data Integration

- **Unified object**: `TSENATAnalysis` encapsulates sequencing data,
  configuration, and all statistical results
- `SummarizedExperiment` foundation: Full Bioconductor ecosystem
  compatibility
- Accessor functions:
  [`diversity()`](https://gallardoalba.github.io/TSENAT/reference/diversity.md),
  [`divergence()`](https://gallardoalba.github.io/TSENAT/reference/divergence.md),
  [`lmResults()`](https://gallardoalba.github.io/TSENAT/reference/lmResults.md),
  [`jisResults()`](https://gallardoalba.github.io/TSENAT/reference/jisResults.md),
  etc.

## Related Packages

TSENAT addresses a fundamental but underappreciated question in
transcriptomic analysis: **How do genes reorganize their isoform usage
patterns, independent of changes in total abundance?** This question is
distinct from standard differential expression analysis and reveals a
layer of biological complexity—coordinated isoform switching—that
conventional methods overlook. TSENAT fills a specific niche in the
Bioconductor ecosystem by measuring scale-dependent isoform diversity
rather than abundance or individual transcript shifts. Below is how
TSENAT complements other Bioconductor tools:

| Tool | Answers | TSENAT Difference |
|----|----|----|
| **DESeq2, edgeR, limma** | Which genes change in *total abundance*? | TSENAT detects isoform diversity changes **independent of total abundance** |
| **DRIMSeq** | Which *individual transcripts* shift usage? | TSENAT measures overall isoform diversity, not individual transcript shifts |
| **IsoformSwitchAnalyzeR** | Which *individual isoforms* switch; what are the *functional consequences*? | TSENAT measures overall isoform diversity and diversity **shifts** rather than cataloging individual transcript switches or predicting functional consequences; complements switch identification with diversity patterns |
| **SplicingFactory** | What is the overall isoform diversity? | TSENAT extends with **scale-dependent diversity** (q-spectrum) vs fixed measures |
| **Kallisto, Salmon** | How many reads per transcript? | TSENAT uses their quantification as input; adds diversity analysis layer |

**Note on SplicingFactory:** TSENAT shares Shannon and Simpson diversity
metrics (identical mathematical definitions) with SplicingFactory but
extends with **scale-dependent analysis** (q-spectrum 0≤q≤2) and
advanced statistical inference. **Important limitation:** TSENAT cannot
compute Gini index. Choose TSENAT for comprehensive q-spectrum analysis
and robust inference; keep SplicingFactory if Gini index analysis is
required. See [Appendix
A](#appendix-a-detailed-tsenat-vs-splicingfactory-comparison) for
detailed comparison and migration guidance.

## Native Salmon Integration

TSENAT is specifically engineered to work seamlessly with Salmon
quantification output. Rather than requiring manual parsing or format
conversion, TSENAT automatically discovers transcript-level
quantification files across your Salmon output directory and integrates
them directly into the analysis pipeline. This tight integration means
you can move from Salmon quantification to entropy analysis without
intermediate data manipulation—the raw `quant.sf` files are all you
need. TSENAT discovers these files automatically, validates their
compatibility with your experimental design, and handles
length-correction and normalization as part of the diversity computation
workflow.

To get started with Salmon-quantified data:

``` r

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

Your Salmon output must be organized with one folder per sample, each
containing a quant.sf file with transcript-level quantification. Sample
folder names must exactly match the row names in your metadata file
(case-sensitive) to ensure correct sample attribution.

TSENAT expects Salmon output organized with one subdirectory per sample:

    salmon_output/
    ├── Sample_1/
    │   └── quant.sf
    ├── Sample_2/
    │   └── quant.sf
    ├── Sample_3/
    │   └── quant.sf
    └── Sample_4/
        └── quant.sf

**Critical requirement**: Folder names (e.g., `Sample_1`, `Sample_2`)
must **exactly match** the row names in your metadata file
(case-sensitive).

## Metadata File Structure

Example metadata structure:

``` r

# Load metadata from TSV file
metadata_df <- read.table("metadata.tsv", header = TRUE, sep = "\t")
```

Expected TSV file format (`metadata.tsv`):

    sample        condition    paired_samples
    SRR14800481   normal       A
    SRR14800480   normal       B
    SRR14800479   tumor        A
    SRR14800478   tumor        B

**Key requirements:** - First column: `sample` (must match Salmon folder
names exactly) - Second column: `condition` (experimental groups:
normal, tumor, treated, control, etc.) - Third column: `paired_samples`
(required if using paired designs; identifier for matched samples)

## Tests coverage

Testing is vital in research as it ensures the validity and reliability
of results, which is essential for accurately interpreting findings. The
report about the current testing coverage can be found
[here](https://app.codecov.io/gh/gallardoalba/TSENAT).

## Learn More

### Comprehensive Workflow

For a complete walkthrough of the analysis pipeline with real biological
examples, see the main package vignette. This includes theory
background, step-by-step explanations of each analysis function, and
interpretation guidance for understanding your results.

See the package vignette for detailed examples, theory background, and
typical workflows:

``` r

vignette("TSENAT")
```

### Function Reference

Use R’s built-in help system to explore detailed documentation for
individual TSENAT functions and S4 classes. Each help page includes
function arguments, return values, and practical examples of usage.

Interactive help for functions and classes:

``` r

?build_analysis_s4
?calculate_diversity_s4
?TSENATAnalysis-class
```

For methodology details and a comprehensive bibliography, see the
[TSENAT
vignette](https://gallardoalba.github.io/TSENAT/vignettes/TSENAT.Rmd):

``` r

vignette("TSENAT")
```

## Citation

If you use TSENAT in your research, please cite:

``` r

citation("TSENAT")
```

BibTeX entry:

``` bibtex
@software{gallardo2026tsenat,
  title={TSENAT: Tsallis Entropy Analysis Toolbox},
  author={Gallardo Alba, Cristóbal},
  url={https://github.com/gallardoalba/TSENAT},
  year={2026}
}
```

## License and Attribution

This project is licensed under the GNU General Public License v3.0
(GPL-3). See [LICENSE](https://gallardoalba.github.io/TSENAT/LICENSE)
for details.

Attribution: TSENAT builds upon the [SplicingFactory
package](https://github.com/esebesty/SplicingFactory), extending it with
specialized focus on Tsallis entropy analysis.

> **“If I ever come back from the past, it’s to create a cyclone.”**
>
> - Juan José Lozano
