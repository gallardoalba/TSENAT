# TSENAT: Tsallis Entropy Analysis Toolbox

## Overview

TSENAT (Tsallis Entropy Analysis Toolbox) quantifies **isoform
heterogeneity** in RNA-seq data using Tsallis entropy-a scale-dependent
diversity measure distinct from abundance-focused tools (DESeq2, Salmon)
and differential transcript usage packages (*SplicingFactory*). By
tuning a sensitivity parameter `q`, users examine diversity at different
scales: rare variants (low `q`) or dominant isoforms (high `q`). This
package enables computing Tsallis entropy from transcript-level
abundance estimates, comparing measures between groups, and visualizing
scale-dependent differences via q-curves. Results integrate seamlessly
with Bioconductor’s `SummarizedExperiment`, capturing scale-dependent
isoform complexity that often signals cellular state changes independent
of total abundance shifts.

### Motivation and Bioconductor Contribution

Common RNA-seq tools focus on either total gene abundance changes or
individual transcript usage shifts. Yet genes often remodel their
isoform landscapes without changing overall expression-a phenomenon
missed by these approaches. TSENAT quantifies this *isoform complexity*
directly via Tsallis entropy, a scale-dependent diversity framework
tuned by parameter `q`. By sliding `q` across scales, researchers zoom
between rare variants (low `q`) and dominant isoforms (high `q`),
capturing biological signal invisible to abundance- or proportion-based
summaries.

In this guide, we demonstrate the complete workflow: preprocessing
transcript counts, computing entropy across q-values, testing for
between-group differences, and visualizing scale-dependent complexity.
All results integrate seamlessly with Bioconductor’s
`SummarizedExperiment`, making TSENAT a natural complement to existing
DTU and abundance-focused tools in the ecosystem.

### High-level workflow

1.  **Load and preprocess** transcript counts using
    [`build_analysis_s4()`](https://gallardoalba.github.io/TSENAT/reference/build_analysis_s4.md),
    then filter low-abundance transcripts with
    [`filter_analysis_s4()`](https://gallardoalba.github.io/TSENAT/reference/filter_analysis_s4.md).
2.  **Select q-values** tuned to your research question: low q (rare
    isoforms) vs. high q (dominant isoforms). Compute Tsallis entropy
    across these scales with
    [`calculate_diversity_s4()`](https://gallardoalba.github.io/TSENAT/reference/calculate_diversity_s4.md).
3.  **Test for differences** in entropy between groups using
    [`rank_test_q_condition_s4()`](https://gallardoalba.github.io/TSENAT/reference/rank_test_q_condition_s4.md)
    (Wilcoxon or permutation tests) (Kerby 2014; Saulsbury 2020).
4.  **Visualize results** as q-curves and inspect transcript counts for
    genes with the strongest entropy shifts using
    [`plot_tsallis_q_curve_s4()`](https://gallardoalba.github.io/TSENAT/reference/plot_tsallis_q_curve_s4.md)
    and related functions.

**Design assumptions**: This guide assumes paired or longitudinal
designs with \>=6-8 samples per group for adequate power to detect
entropy shifts while controlling false discovery rate. Smaller sample
sizes may be underpowered to detect subtle isoform complexity changes.

### Installation

To install `TSENAT` from Bioconductor:

``` r

if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager", repos = "http://cran.us.r-project.org")
BiocManager::install("TSENAT")
```

Then load the package:

``` r

library(TSENAT)
```

### Quick Start

For users who want to see results quickly, here is a minimal workflow:

``` r

library(TSENAT)
library(SummarizedExperiment)

# Load example data
data("readcounts", package = "TSENAT")
readcounts <- as.matrix(readcounts)

# Load sample metadata and annotation
metadata_df <- read.table(
  system.file("extdata", "metadata.tsv", package = "TSENAT"),
  header = TRUE, sep = "\t"
)
gff3_file <- system.file("extdata", "annotation.gff3.gz", package = "TSENAT")

# Create analysis object with transcript counts, annotation, and metadata
config <- tsenat_config(
  sample_col = "sample",
  condition_col = "condition",
  paired = TRUE,
  subject_col = "paired_samples",
  control = "normal",
  q_values = seq(0, 2, by = 0.05)
)

## Build TSENATAnalysis object
analysis <- build_analysis_s4(
  readcounts = readcounts, 
  tx2gene = gff3_file,
  metadata = metadata_df,
  config = config,
  tpm = tpm,
  effective_length = effective_length)

# Filter low-abundance transcripts
analysis <- filter_analysis_s4(analysis)

# Compute Tsallis entropy using S4 wrapper (using single q value for quick start)
analysis <- calculate_diversity_s4(analysis)

# Plot overall q-curve
p_qcurve <- plot_tsallis_q_curve_s4(analysis)
print(p_qcurve)
```

### Data Structure Overview

**Input Data**: TSENAT expects transcript-level read counts (rows =
transcripts, columns = samples) from quantification tools like SALMON or
kallisto. A GFF3 annotation file maps transcripts to genes.

**Processing**: The
[`build_analysis_s4()`](https://gallardoalba.github.io/TSENAT/reference/build_analysis_s4.md)
function accepts these inputs and creates a **TSENATAnalysis** S4
object, which is the central data container throughout your analysis.

**Organization**: The TSENATAnalysis object encapsulates:

- `@se`: `SummarizedExperiment` storing counts and metadata.
- `@diversity_results`: Entropy values across q-values.
- `@jackknife_results`: Jackknife confidence intervals and isoform
  switching results.
- `@lm_results`: LM interaction statistics and rank test q-value
  effects.
- `@divergence_results`: Pairwise divergence metrics.
- `@plots`: Generated visualizations.

For more details on the SummarizedExperiment class, see
[SummarizedExperiment
documentation](https://bioconductor.org/packages/release/bioc/html/SummarizedExperiment.html).

## What is Entropy and Tsallis Entropy?

### The Problem: Measuring Isoform Complexity

Standard RNA-seq analysis measures *whether* transcript abundance
changes between conditions. But a critical complementary question
remains underexplored: *how* does the diversity of isoforms change? A
gene may show little change in total abundance while dramatically
reshuffling its isoform repertoire-a phenomenon that current methods
largely miss.

**Entropy** quantifies precisely this: the complexity, richness, and
balance of isoform heterogeneity. By measuring entropy across different
biological scales, researchers can detect whether changes are driven by
shifts in rare variants or reorganization of dominant isoforms.

### Tsallis Entropy: A Scale-Dependent Diversity Measure

**Tsallis entropy** (Tsallis 2006) is a one-parameter family of
diversity measures that generalizes Shannon entropy (Shannon 1948).
Unlike Shannon entropy (which treats all isoforms equally), Tsallis
entropy lets you tune a sensitivity parameter `q` to zoom into different
aspects of isoform complexity:

- **q \< 1**: Emphasizes rare, low-abundance isoforms (discovery mode).
- **q = 1**: Recovers Shannon entropy (balanced across all abundance
  scales).
- **q = 2**: Emphasizes dominant isoforms (robustness mode).
- **q \> 2**: Focuses almost exclusively on the most abundant species.

This is the key innovation: by computing entropy across a range of
q-values (a “q-curve”), you obtain a complete picture of isoform
heterogeneity (Rényi 1961). In TSENAT, you’ll explore multiple q-values
using
[`calculate_diversity_s4()`](https://gallardoalba.github.io/TSENAT/reference/calculate_diversity_s4.md)
and visualize results with
[`plot_tsallis_q_curve_s4()`](https://gallardoalba.github.io/TSENAT/reference/plot_tsallis_q_curve_s4.md).

### Mathematical Definition

For a discrete probability vector $`p = (p_1, \ldots, p_n)`$
representing isoform proportions within a gene, Tsallis entropy is:

``` math
S_q(p) = \frac{1-\sum_{i=1}^n p_i^q}{q-1}.
```

This parametric family unites diverse entropy concepts under a single
framework:

- **Generalization**: Extends beyond Shannon entropy to capture
  scale-dependent phenomena
- **Mathematical elegance**: Reduces to well-known diversity indices at
  specific q-values
- **Practical flexibility**: Enables data-driven exploration across the
  full diversity spectrum

### Information-Theoretic Interpretation

From an information theory perspective (Shannon 1948; Tsallis 2006),
entropy measures the uncertainty or surprise when drawing a single
transcript from an isoform distribution. Higher entropy means the draw
is less predictable (many similarly abundant isoforms), while lower
entropy means one or a few isoforms dominate. This distinction is
crucial: two genes with identical total abundance may have dramatically
different isoform complexity.

The q parameter acts as a **sensitivity dial**. Consider a simple
example: - A gene with 5 equally abundant isoforms: all q-values give
high entropy - A gene where 1 isoform dominates: low-q entropy stays
high (rare variants matter), while high-q entropy drops (rare variants
ignored)

This scale-dependent nature reveals biological signal invisible to
abundance- or proportion-based summaries alone.

#### Special Cases and Limiting Behavior

Tsallis entropy exhibits important special cases (Masi 2005) that appear
when you set specific q-values:

- **q = 0 (Richness)**: $`S_0 = m - 1`$ (number of expressed isoforms).
  Pure species count, most minimal assumption.
- **q = 1 (Shannon)**: $`S_1(p) = -\sum_i p_i \log p_i`$. Standard
  information entropy (Shannon 1948); optimal estimator of information
  content.
- **q = 2 (Gini-Simpson)**: $`S_2 = 1 - \sum_i p_i^2`$. Probability that
  two randomly drawn transcripts differ (Simpson 1949); robust to rare
  variants.

**Normalized entropy**: For genes with different isoform counts,
normalize by maximum possible entropy:
$`\tilde{S}_q = \frac{S_q(p)}{S_{q,\max}(m)}`$ where
$`S_{q,\max}(m) = \frac{1-m^{1-q}}{q-1}`$. This enables fair comparison
across genes.

### Biological Potentiality: Why TSENAT Matters

Tsallis entropy is grounded in Shannon’s foundational information theory
(Shannon 1948), generalized by Renyi’s family of entropies (Rényi 1961),
and extended by Tsallis (Tsallis 2006). The key innovation is a tunable
parameter `q` that controls sensitivity to different aspects of the
diversity distribution-which aspects of isoform heterogeneity become
visible depends on the q-value chosen (Anastasiadis 2012; Ramírez-Reyes
et al. 2016; Alomani and Kayid 2023).

As stated explicitly in the mathematical literature (Masi 2005): “This
introduces the formal possibility not to set rare and common events on
the same footing, as in BG or Shannon statistics, but it enhances or
depresses them according to the parameter chosen.” This principle is
formalized in the Hill numbers framework (Chao et al. 2010), which
unifies diverse entropy-based measures (richness, Shannon, Simpson) as
different manifestations of the same parametric family with varying
sensitivity to abundance scales.

The flexibility is crucial for isoform analysis. The concept of “true
diversity” (Chao et al. 2010) emphasizes that diversity can be
decomposed into two independent components: species richness (how many
distinct isoforms exist) and evenness (how evenly distributed they are
across the population). Two genes can have identical Shannon entropy yet
differ dramatically in isoform structure: one might be dominated by a
single abundant isoform (maximizing apparent heterogeneity at low `q`
where rare variants matter), while another distributes transcripts
across many isoforms equally (maximizing heterogeneity across all `q`
values). By tuning the q-parameter, we operationalize what Masi (2005)
describes: depressing the weight on rare events (high q) reveals
dominant isoform patterns, while enhancing rare-event weight (low q)
reveals cryptic isoform complexity. TSENAT’s multi-*q* approach captures
this full spectrum of isoform diversity, enabling researchers to detect
both rare isoform innovations and robust isoform usage patterns (Drost
2018) by exploring the q-curve across multiple scales.

#### Potential Applications in Genomics

The multi-scale nature of Tsallis entropy makes it suited for exploring
isoform complexity across diverse biological contexts. By measuring
information content at different q-values, researchers can detect
patterns invisible to traditional transcript abundance measures alone:

**Isoform complexity as a biological signal**: Isoform
switching-reorganization of the isoform landscape without necessarily
changing total gene abundance-reflects strategic shifts in protein
function driven by splicing regulation. Evidence from single-cell
transcriptomics demonstrates that transcript-level complexity varies
systematically across cell types and developmental states (Cao et al.
2017), validating that isoform heterogeneity is a genuine biological
phenomenon rather than noise. Increased entropy in gene regulatory
networks drives phenotypic heterogeneity and cellular plasticity (Nijman
2020), suggesting that transcript-level entropy captures similar
organizational principles. TSENAT enables detection of these changes
through entropy-based approaches, which capture whether complexity is
increasing (diversity spreading across isoforms) or decreasing
(consolidation onto dominant isoforms).

**Scale-dependent organization**: Different biological processes may
prioritize different scales of organization (Tarabichi et al. 2013): -
Changes in rare isoform usage (revealed through low `q` sensitivity)
might reflect exploratory or error-correction mechanisms - Shifts in
dominant isoform selection (revealed through high `q` sensitivity) might
reflect functional specialization or robustness demands - The full
q-curve reveals whether cellular transitions involve wholesale
reorganization or targeted adjustments, paralleling how systems biology
approaches dissect biological disorder and organization at multiple
scales

**Beyond abundance measures**: Traditional analysis focuses on
fold-changes and differential abundance. Since isoform reorganization
can occur independently of total abundance changes, entropy-based
approaches complement classical methods by detecting complexity shifts
that transcript-level statistics alone cannot reveal.

#### Evidence from Literature

The recent emphasis on information-theoretic approaches in computational
biology Bajić (2024) reflects broader recognition that complex
biological systems encode information across multiple organizational
scales. Peer-reviewed literature provides strong empirical support for
these applications:

**Entropy and cancer biology**: Cancer cells accumulate genetic and
epigenetic perturbations that systematically increase disorder in gene
regulatory networks. As Tarabichi and colleagues demonstrate, “Increased
entropy of signaling (or gene interaction networks) has been well
studied as a cancer characteristic: Network entropy increases along with
cancer progresses” (Tarabichi et al. 2013). Systems biology approaches
reveal that this entropy increase, rather than being an incidental
feature, actively drives cancer progression through selection of cells
with greater network flexibility and adaptability.

**Perturbation-driven mechanisms of heterogeneity**: Nijman’s analysis
of perturbation-driven entropy proposes a complementary mechanism:
“cancer-associated perturbations collectively disrupt normal gene
regulatory networks by increasing their entropy. Importantly, in this
model both somatic driver and passenger alterations contribute to
‘perturbation-driven entropy’, thereby increasing phenotypic
heterogeneity and evolvability” (Nijman 2020). This framework elegantly
explains observed cancer heterogeneity without requiring that every
genetic change confers a selective advantage-some mutations contribute
entropy directly through network disruption.

**Single-cell validation of transcript diversity**: Cao and colleagues’
landmark single-cell transcriptomics study provided empirical validation
that transcript-level organization varies systematically across cell
types: “expression levels of mRNA species are linked to cellular
function and therefore can be used to classify cell types” (Cao et al.
2017). Their comprehensive profiling of C. elegans demonstrates that
individual cells maintain specific, consistent isoform compositions
reflecting cellular identity-establishing that isoform complexity is not
noise but a fundamental aspect of cellular differentiation and function.

These three perspectives-entropy as driver of cancer evolution,
perturbations as network disruption mechanisms, and single-cell evidence
for systematic isoform organization-converge on a framework where
measuring entropy at multiple scales (via TSENAT’s multi-*q* approach)
captures biologically meaningful variation in cellular organization and
adaptation. In cancer biology specifically, entropy-driven mechanisms
explain how perturbations increase heterogeneity and plasticity (Nijman
2020), demonstrating that entropy concepts have mechanistic relevance
beyond abstract information theory. Tsallis entropy, through the
parametric q-spectrum, provides a systematic framework for exploring
this multi-scale organization at the transcript and isoform level.

## Isoform Switching Workflow

We now demonstrate a complete workflow for detecting **isoform
switching**-when cells reorganize their isoform landscape without
necessarily changing total gene abundance. This often reflects strategic
shifts in protein function driven by splicing regulation.

In this workflow, you will:

1.  Load transcript counts and sample metadata
2.  Identify genes with significant scale-dependent isoform complexity
    changes (using linear model interaction testing)
3.  Assess which individual transcripts drive those changes (using
    jackknife resampling (Efron and Tibshirani 1993))
4.  Interpret results in the context of paired experimental designs

We’ll work with a paired experimental design where treated and control
samples are linked, enabling detection of robust biological signals
while controlling for subject-level variability. By the end, you’ll know
not just *which genes* undergo isoform switching, but *which
transcripts* are responsible and *at which diversity scales* the
switching occurs.

### Load data and metadata

An example dataset is included for demonstration.

``` r

# Load packages
suppressPackageStartupMessages({
    library(TSENAT)
    library(ggplot2)
    library(SummarizedExperiment)
})

# Setup random seed for reproducibility
set.seed(42)
```

Now we will load the example dataset and associated metadata:

``` r

# Load example dataset (lazy-loaded as 'readcounts' by default)
# This includes SALMON preprocessing outputs:
# - readcounts: transcript-level NumReads (raw fragment counts)
# - tpm: TPM matrix (length- and library-normalized estimates)
# - effective_length: EffectiveLength vector (read-length corrected)
data(readcounts)

readcounts <- as.matrix(readcounts)

metadata_df <- read.table(
    system.file("extdata", "metadata.tsv", package = "TSENAT"),
    header = TRUE, sep = "\t"
)

gff3_file <- system.file("extdata", "annotation.gff3.gz", package = "TSENAT")
```

### Data preprocessing and filtering

Next, we will create a `TSENATAnalysis` object containing a
`SummarizedExperiment` data container (a standard Bioconductor class for
storing high-throughput assay data along with associated metadata). You
can read more about it in the [Bioconductor
documentation](https://bioconductor.org/packages/release/bioc/html/SummarizedExperiment.html).

We will use the GFF3.gz annotation file for extracting the
transcript-to-gene mapping and building the analysis object. We’ll also
include sample metadata for improved analysis. Following Bioconductor
best practices, we create the configuration FIRST, then pass it to the
builder function (fail-fast principle):

``` r

## Configure analysis parameters first (best practice: fail-fast principle)
## This validates all parameters before object creation
config <- tsenat_config(
  sample_col = "sample",
  condition_col = "condition",
  subject_col = "paired_samples",
  q_values = seq(0, 2, by = 0.05),
  nthreads = 2,
  paired = TRUE,
  control = "normal"
)

## Build a complete `TSENATAnalysis` object from readcounts + GFF3.gz annotation
## Pass config at construction (Bioconductor pattern): immutable object creation
## Note: metadata is passed as explicit parameter (Bioconductor best practice)
analysis <- build_analysis_s4(
  config = config,
  readcounts = readcounts,
  metadata = metadata_df,
  tx2gene = gff3_file, 
  tpm = tpm,
  effective_length = effective_length
)
```

The `TSENATAnalysis` object contains a `SummarizedExperiment` with
transcript-level counts and gene annotations extracted from the GFF3
file.

To reduce noise and improve statistical power, we filter out
lowly-expressed transcripts. Expression estimates of transcript isoforms
with zero or low expression might be highly variable (Jose and Lal 2013;
Chakraborty 2019). For more details on the effect of transcript isoform
prefiltering on differential transcript usage, see [this
paper](https://doi.org/10.1186/s13059-015-0862-3).

``` r

analysis <- filter_analysis_s4(analysis, stringency = "medium")
cat("Retained", nrow(se(analysis)), "transcripts after filtering\n")
#> Retained 341 transcripts after filtering
```

The `stringency = "medium"` parameter keeps transcripts present in
\>=50% of samples with TPM values above the median of mean gene TPM.
This balances noise reduction with preservation of isoform diversity
needed for meaningful entropy calculations.

### Compute Tsallis Entropy

We now compute Tsallis entropy across your configured q-spectrum.
Diversity measures provide a comprehensive framework for assessing
isoform heterogeneity at multiple scales (Ramírez-Reyes et al. 2016; Gao
et al. 2019). We normalize entropy to \[0,1\] range (`norm = TRUE`) for
comparable cross-q assessment.

``` r

# Set random seed for reproducible bootstrap CI calculations
set.seed(12345)

# Compute diversity using S4 wrapper with bootstrap confidence intervals [@S232; @S115]
analysis <- calculate_diversity_s4(
    analysis,
    norm = TRUE,
    output_file = "vignette_diversity_results.tsv"
)
#> Note: 12 genes excluded (< 75% valid values).
#> [calculate_diversity_s4] Saved diversity spectrum to: vignette_diversity_results_spectrum.tsv
#> [calculate_diversity_s4] Saved table to: vignette_diversity_results.tsv
#> [calculate_diversity_s4] Saved diversity results to: vignette_diversity_results.tsv
```

| Gene    |   q=0.0 |   q=0.5 |   q=1.0 |   q=1.5 |   q=2.0 |
|:--------|--------:|--------:|--------:|--------:|--------:|
| FOXJ2   | 0.66667 | 0.53796 | 0.48150 | 0.46975 | 0.48518 |
| TMEM38A | 1.00000 | 0.82704 | 0.72526 | 0.66965 | 0.64401 |
| CCNE1   | 0.33333 | 0.32256 | 0.32749 | 0.34571 | 0.37418 |
| MCOLN1  | 1.00000 | 0.39297 | 0.21567 | 0.16305 | 0.15195 |

**Table 1:** Tsallis entropy for first 4 genes across three diversity
scales (sample: SRR14800481) {.table .table .table-striped .table-hover
.table-condensed style="margin-left: auto; margin-right: auto;"}

With the q-spectrum we can produce a q-curve per sample and gene. These
curves show how diversity emphasis shifts from rare to dominant isoforms
as `q` increases and form the basis for interaction tests. The q-curve
shows entropy as a function of `q`. Diverging curves between groups
indicate scale-dependent diversity differences: separation at low `q`
implies differences in rare isoforms, while separation at high `q`
signals differences in dominant isoforms.

``` r

# Plot overall q-curve
p_qcurve <- plot_tsallis_q_curve_s4(analysis)

print(p_qcurve)
```

![\*\*Figure 1:\*\* Isoform diversity profiles across q-values. Lines
show normalized Tsallis entropy (0-1) for each sample, blue control and
red
treatment.](TSENAT_files/figure-html/fig-1-isoform-diversity-profiles-1.png)

**Figure 1:** Isoform diversity profiles across q-values. Lines show
normalized Tsallis entropy (0-1) for each sample, blue control and red
treatment.

### Quality Control: Sample Influence Assessment

Before proceeding with group-level comparisons, we should assess whether
individual samples exert disproportionate influence on our entropy
estimates (Efron, Bradley and Tibshirani, Robert J. 1993). To address
this, we employ a leave-one-out influence assessment combined with
robust M-estimation (Phipson and Smyth 2010) (iteratively re-weighted
least squares with Huber loss). This approach quantifies how much each
sample’s removal affects the estimated location differences across the
q-spectrum, providing a sample-level quality control metric independent
of group assignment.

Samples with high influence scores exert disproportionate leverage on
entropy estimates and may warrant deeper investigation for technical
quality issues, subject-specific outliers, or genuine biological signals
that require careful interpretation.

We configure three key parameters: `loss_type = "huber"` employs Huber
M-estimation (robust to outliers with bounded influence),
`q_combine_method = "mean"` averages influence metrics across the entire
q-spectrum for a single score per sample, and
`influence_threshold = 0.75` designates the percentile cutoff above
which samples are flagged (such that the top 25% most influential
samples are reviewed).

``` r

# Perform multi-q sample influence analysis via S4 wrapper
analysis <- m_estimate_s4(
    analysis,
    loss_type = "huber",
    q_combine_method = "mean",
    influence_threshold = 0.75
)

# Extract results from metadata using accessor function
sample_qc <- metadata(analysis, "m_estimate_results")
```

| Sample | Condition | Proportion_Affected | Genes_Affected | Entropy_Mean | Entropy_SD | Distance_from_Centroid | Status |
|:---|:---|---:|---:|---:|---:|---:|:---|
| SRR14800481 | normal | 0.847 | 64.4 | 0.6936 | 0.2364 | 1.7930 | OK |
| SRR14800480 | normal | 0.806 | 61.2 | 0.7345 | 0.1972 | 1.2384 | OK |
| SRR14800477 | normal | 0.833 | 63.3 | 0.6896 | 0.2852 | 1.5809 | OK |
| SRR14800476 | normal | 0.889 | 67.6 | 0.7310 | 0.2329 | 0.9118 | OK |
| SRR14800475 | normal | 0.861 | 65.4 | 0.6774 | 0.2435 | 1.6796 | OK |
| SRR14800490 | normal | 0.849 | 64.5 | 0.6843 | 0.2357 | 1.3522 | OK |
| SRR14800489 | normal | 0.875 | 66.5 | 0.7339 | 0.2426 | 1.0969 | OK |
| SRR14800488 | normal | 0.861 | 65.4 | 0.6738 | 0.2649 | 1.6536 | OK |
| SRR14800479 | tumor | 0.861 | 65.4 | 0.6842 | 0.2721 | 2.2114 | OK |
| SRR14800478 | tumor | 0.903 | 68.6 | 0.7345 | 0.2440 | 1.5721 | Flag for QC |
| SRR14800487 | tumor | 0.861 | 65.4 | 0.7257 | 0.2193 | 1.2201 | OK |
| SRR14800486 | tumor | 0.903 | 68.6 | 0.7207 | 0.2497 | 1.4798 | Flag for QC |
| SRR14800485 | tumor | 0.861 | 65.4 | 0.7334 | 0.2131 | 1.4140 | OK |
| SRR14800484 | tumor | 0.889 | 67.6 | 0.7382 | 0.2095 | 1.1499 | OK |
| SRR14800483 | tumor | 0.861 | 65.4 | 0.7425 | 0.2283 | 1.6020 | OK |
| SRR14800482 | tumor | 0.917 | 69.7 | 0.7186 | 0.2446 | 1.6711 | Flag for QC |

**Table 2 \| Sample Influence Assessment via M-estimation.** Ranked by
magnitude of resampling-based influence values. Columns: Sample
identifier; mean and standard deviation of Tsallis entropy; distance
from centroid; outlier status. M-estimation identifies influential
samples for sensitivity analysis. {.table .table .table-striped
.table-hover .table-condensed
style="margin-left: auto; margin-right: auto;"}

**Interpretation:** Samples with distance from centroid \>1.5 or
proportion of affected genes \>0.85 are flagged for quality control
review. Flagged samples may reflect genuine biological heterogeneity or
technical artifacts requiring careful examination before proceeding with
downstream analysis.

### Linear-Model interaction

We can test for interactions between `q` and sample groups across the
q-sequence using one of five methods: linear, lmm (linear mixed model
(Le Roux and Rouanet 2011; Phipson and Smyth 2010)), gam (generalized
additive model (Chapman 2017; Goude 2024)), gee (generalized estimating
equations for paired data), or fpca (functional principal component
analysis). The gam method flexibly captures nonlinear q-response
patterns and is particularly useful for complex interaction structures.

``` r

# Linear-model interaction test across q values using S4 wrapper
if (requireNamespace("mgcv", quietly = TRUE)) {
    analysis <- calculate_lm_interaction_s4(
        analysis,
        method = "gam",
        multicorr = "hochberg"
    )
}
```

| Gene | Gene Name | P-value | Adj. P-value | Effect Size | Test Statistic | Model Converged | Heteroscedasticity |
|:---|:---|:---|:---|---:|---:|:---|:---|
| CXCL12 | CXCL12 | 3.93e-164 | 2.99e-162 | 0.8127 | 752.5086 | TRUE | TRUE |
| THY1 | THY1 | 9.80e-97 | 7.35e-95 | 0.6008 | 442.1369 | TRUE | TRUE |
| ING3 | ING3 | 2.72e-77 | 2.02e-75 | 0.3533 | 352.5945 | TRUE | TRUE |
| SNHG10 | SNHG10 | 9.81e-75 | 7.16e-73 | 0.0515 | 340.8200 | TRUE | TRUE |
| LINC03040 | LINC03040 | 1.87e-57 | 1.35e-55 | 0.7380 | 261.2380 | TRUE | TRUE |
| HDAC2 | HDAC2 | 2.61e-51 | 1.85e-49 | 0.2844 | 232.9431 | TRUE | TRUE |

**Table 3 \| Leading genes with scale-dependent condition effects from
generalized additive models.** We display the 6 genes with lowest
adjusted p-values (*q* \< 0.05). Benjamini-Hochberg correction applied;
columns show gene identifier, effect size, test statistic, convergence
status, and heteroscedasticity detection. Methods: GAMs with *q* and
condition as smooth predictors (Benjamini and Hochberg 1995). {.table
.table .table-striped .table-hover .table-condensed
style="margin-left: auto; margin-right: auto;"}

**Interpretation:** TRUE in the Heteroscedasticity column indicates that
the model satisfies homogeneity-of-variance assumptions across the
q-spectrum; FALSE values suggest variance heterogeneity and warrant
caution in result interpretation.

Now we will plot the q-curve profile for the top genes identified by the
linear-model interaction test.

``` r

# Plot q-curve profiles for the top 4 genes using the S4 wrapper
combined_plot <- plot_lm_interaction_gam_s4(
    analysis,
    n_top = 4
)
print(combined_plot)
```

![\*\*Figure 2:\*\* Scale-dependent interaction analysis. GAM-identified
genes showing significant q\$\times\$condition effects
(Benjamini-Hochberg q \<
0.05).](TSENAT_files/figure-html/fig-2-scale-dependent-genes-1.png)

**Figure 2:** Scale-dependent interaction analysis. GAM-identified genes
showing significant q$`\times`$condition effects (Benjamini-Hochberg q
\< 0.05).

### Transcript Switching Across Diversity Scales

#### Why Scale-Dependent Analysis Matters

Biological systems operate through processes at fundamentally different
organizational scales. A single gene’s expression can be viewed through
different “lenses”-from the perspective of dominant isoforms (high *q*,
emphasizing abundant variants) to rare isoforms that may carry
functional significance in specific contexts (low *q*, emphasizing rare
variants). The LM interaction test revealed that condition effects vary
across these diversity scales (q parameter), reflecting how biological
reorganization can manifest differently depending on which aspect of
isoform complexity we examine.

The theoretical foundation for scale-dependent analysis comes from
generalized entropy measures. Renyi’s entropy family (Rényi 1961)
introduced a continuous spectrum of diversity measures parameterized by
a sensitivity index, which Tsallis later extended and popularized
(Alomani and Kayid 2023). Unlike Shannon entropy (which represents a
single point on this spectrum), Tsallis entropy provides a family of
measures where the parameter *q* acts as a “magnifying glass”-adjusting
what organizational scales become visible. Mathematically, different
*q*-values weight rare versus abundant isoforms differently:

- Low *q* (0-0.5): Emphasized rare isoforms, revealing exploratory or
  conditional splicing programs
- Mid *q* (0.5-1.5): Balanced sensitivity across the full isoform
  distribution  
- High *q* (1.5-2.0): Emphasized abundant isoforms, revealing core
  splicing decisions

In biological systems, these different scales often correspond to
different regulatory mechanisms. A transcript may be consistently rare
across conditions (no switching at low *q*) while being stably abundant
at high *q*. Detecting such scale-specific switching patterns reveals
the mechanistic structure of splicing regulation (Tarabichi et al.
2013).

#### Bridging Theory: Multi-*q* Analysis to Fine-Scale Jackknife Testing

The current analysis bridges two complementary approaches:

1.  **Multi-*q* LM interaction test**: Tests whether the condition
    effect on divergence depends on *q* (overall pattern across
    diversity scales)

2.  **Single-*q* jackknife resampling**: Tests which individual
    transcripts drive those patterns at each *q*-value, with robustness
    assessment

This two-stage approach combines hypothesis testing with robustness
validation (Efron, Bradley and Tibshirani, Robert J. 1993). The LM test
identifies *where* (at which *q*-values) significant switching occurs.
The jackknife-based delta influence then identifies *which transcripts*
contribute to those patterns, weighting confidence by operational
stability rather than simple frequency.

#### Interpreting Scale-Dependent Switching Patterns

Here we assess whether these switching patterns are consistent across
different q values or vary by scale. The key question is: **Do the same
transcripts switch across all diversity scales, or does the identity of
switching transcripts depend on which scale we examine?**

- **Consistent switching**: The same transcripts show significant delta
  influence across *q*-values suggests a robust, scale-independent
  splicing shift
- **Mixed/scale-dependent switching**: Different transcripts show
  importance at different *q*-values suggests regulatory complexity
  where the mechanism differs by scale

This classification reveals whether isoform reorganization is a
coordinated, global phenomenon (same transcripts across scales) or a
layered process where different regulatory inputs dominate at different
organizational scales.

``` r

# Multi-Q analysis: Test switching at different q values
analysis <- jackknife_isoform_switching_s4(
    analysis,
    q = c(0, 0.5, 1, 1.5, 2),
    nboot = 100,
    lm_p_threshold = 0.05,
    threshold = 90
)
```

We configure three key parameters for the switching analysis:
`nboot = 100` specifies 100 bootstrap resamples to estimate robust
confidence intervals for delta influence values (Efron, Bradley and
Tibshirani, Robert J. 1993); `lm_p_threshold = 0.05` filters genes to
those with significant q$`\times`$condition interaction effects (p \<
0.05) before assessing transcript switching; and `threshold = 90`
designates transcripts with support $`\geq`$ 90% across bootstrap
samples as robust switches. These parameters balance sensitivity
(detecting switching) with specificity (avoiding false positives from
noise).

The tables below examine transcript switching patterns for the top genes
with strongest q$`\times`$condition interactions. Delta influence,
derived from jackknife resampling (Efron, Bradley and Tibshirani, Robert
J. 1993), quantifies how each transcript’s relative importance changes
between conditions, weighted by stability across bootstrap iterations
(positive = more influential in first condition; negative = second
condition). The `Direction Consistency` column classifies each
transcript: “Consistent” indicates stable switching across q-values,
while “Mixed” indicates scale-dependent switching. This reveals whether
isoform shifts involve the same transcripts across scales or whether
different q-values emphasize different transcripts.

``` r

# Prepare comparison tables for top genes across q-values using S4 accessor
tables_result <- prepare_gene_switching_tables_s4(
    analysis, 
    n_top_genes = 2
)
```

##### Gene: CXCL12 (ENSG00000107562.18)

| Transcript         | q=0.00 | q=0.50 | q=1.00 | q=1.50 | q=2.00 |     | Direction Consistency |
|:-------------------|-------:|-------:|-------:|-------:|-------:|:----|:----------------------|
| ENST00000343575.11 | -0.033 |  0.000 |  0.041 |  0.062 |  0.073 |     | Mixed directions      |
| ENST00000374426.6  |  0.046 | -0.063 | -0.051 | -0.040 | -0.035 |     | Mixed directions      |
| ENST00000374429.6  | -0.033 |  0.164 |  0.156 |  0.123 |  0.112 |     | Mixed directions      |

##### Gene: THY1 (ENSG00000154096.15)

| Transcript        | q=0.00 | q=0.50 | q=1.00 | q=1.50 | q=2.00 |     | Direction Consistency |
|:------------------|-------:|-------:|-------:|-------:|-------:|:----|:----------------------|
| ENST00000524970.5 | -0.198 | -0.046 |  0.002 |  0.014 |  0.020 |     | Mixed directions      |
| ENST00000900758.1 |  0.196 |  0.045 |  0.031 |  0.019 |  0.016 |     | Consistent positive   |
| ENST00000956364.1 |  0.118 | -0.075 | -0.092 | -0.102 | -0.108 |     | Mixed directions      |

#### Delta Influence Across Diversity Scales

Visualize switching patterns for the top genes identified by the LM
interaction test. This shows which transcripts are switching in genes
with significant q \* condition interaction effects. The heatmaps below
display **jackknife delta influence** (Efron and Tibshirani 1993)-a
resampling-based measure of how robustly each transcript’s relative
contribution changes between conditions-across different q-values (rows:
0 to 2.0, representing increasing diversity scales from rare to abundant
isoforms) and transcripts (columns) for each gene. Delta influence
values are computed from bootstrap resampling iterations, providing
robust, non-parametric estimates of transcript switching significance
weighted by consistency across replicates.

``` r

# Generate multi-Q heatmaps for top 4 genes using S4 wrapper
plot_multiq_delta_influence_heatmaps_s4(
    analysis, 
    n_genes = 4
)
```

![\*\*Figure 3:\*\* Transcript-level switching patterns via jackknife
resampling. Rows q-values (0-2.0), columns transcripts. Red/blue
indicates condition
influence.](TSENAT_files/figure-html/fig-3-transcript-jackknife-influence-1.png)

**Figure 3:** Transcript-level switching patterns via jackknife
resampling. Rows q-values (0-2.0), columns transcripts. Red/blue
indicates condition influence.

**Interpretation:** Rows represent q-values across the diversity
spectrum (0 to 2.0, with low q emphasizing rare isoforms and high q
emphasizing abundant isoforms). Red cells indicate transcripts with
strong positive delta influence in the first condition (meaning their
relative importance increased); blue indicates negative influence in the
second condition (meaning their importance decreased). **Patterns to
look for:**

- **Strong color intensity (bright red or blue)**: Robust switching
  signal, consistently observed across bootstrap resamples
- **Weak color intensity (pale/white)**: Uncertain or noisy switching
  signal, high variability across resamples
- **Consistent colors down a column**: That transcript switches stably
  across all q-values (scale-independent switching)
- **Varied colors down a column**: That transcript shows scale-dependent
  switching (different roles at rare vs. abundant scales)
- **Color differences across rows**: Different q-values emphasize
  different transcripts, revealing how regulatory mechanisms vary by
  topological scale

This robustness-weighted visualization reveals not just *which*
transcripts switch, but *how reliably* they switch, allowing distinction
between robust biological signals and artifacts of sampling variation.

``` r

# Generate transcript abundance heatmap for top 4 transcripts using hierarchical clustering
# metric = "median" aggregates expression values across paired conditions
plot_top_transcripts_s4(
    analysis, 
    top_n = 4, 
    metric = "median"
)
```

![\*\*Figure 4:\*\* Transcript abundance heatmap with hierarchical
clustering. Rows are transcripts, columns are conditions. Blue low, red
high
expression.](TSENAT_files/figure-html/fig-4-transcript-abundance-heatmap-1.png)

**Figure 4:** Transcript abundance heatmap with hierarchical clustering.
Rows are transcripts, columns are conditions. Blue low, red high
expression.

**Interpretation:** Transcripts (rows) are hierarchically clustered by
expression similarity. Red indicates high expression, blue indicates low
expression. Distinct horizontal bands reveal condition-specific
isoforms; uniform coloring suggests constitutive expression. Clustering
reveals which transcripts co-vary in expression patterns across
conditions.

#### Effect size analysis

To understand which diversity scales show the most pronounced biological
differences between groups, we compute effect size metrics (Tsallis
divergence $`D_q`$) across the q-spectrum while respecting the paired
sample design. This reveals whether group differences are driven by rare
isoforms (low q), dominant isoforms (high q), or uniformly across
scales.

The approach uses linear mixed-effects regression (LMM) with q as
continuous predictor to assess slope differences in entropy change
across q-spectrum between treatment groups, with effect size computed as
Tsallis divergence $`D_q`$(control\|\|treatment) (Kullback and Leibler
1951; Tsallis 2006; Rényi 1961; Sason 2022b, 2022a). This
information-theoretic measure automatically respects Tsallis entropy
properties and adapts to each q value (Shiner et al. 2002).

**Theoretical Justification**

Tsallis entropy theory Sason (2022b) demonstrates that different q
values reveal fundamentally different aspects of the isoform
distribution. As documented in the foundational divergence literature
(Sason 2022b, 2022a; Ré and Azad 2014; R0̆0e9 and Azad 2014),
$`D_{0.5}(P||Q)`$ often differs substantially from $`D_{2}(P||Q)`$, with
q=0.5 emphasizing rare variants and q=2 emphasizing abundant isoforms.
Computing effect sizes independently at each q value thus captures which
diversity scales (rare vs. dominant isoforms) show the strongest
biological differentiation between conditions (Sfetcu et al. 2022; R0̆0e9
and Azad 2014). Re and Azad (2014) demonstrated that Tsallis divergence
generalizations improve discrimination of genomic sequences, validating
the multi-scale divergence approach for biological applications (R0̆0e9
and Azad 2014).

**Effect Size Interpretation**

**Tsallis divergence** $`D_q`$ quantifies the information-theoretic
distance between two distributions’ entropy patterns across the
q-spectrum (Yulmetyev et al. 2004a; Shiner et al. 2007; Wang et al.
2021; Chernyshov 2009; Sason 2022a). For paired designs, divergence is
computed separately for each pair and then averaged, accounting for
within-pair correlation and reducing noise from between-pair variation
(Yulmetyev et al. 2004b). The mean divergence across q values serves as
the effect size, automatically capturing how entropy distributions
differ between control and treatment groups at all scales (rare to
abundant isoforms). Values D \> 0.1 indicate meaningful divergence,
following effect size classification thresholds (Kerby 2014; Chao et al.
2010; Sason 2022a). Genes with D \> 0.1 demonstrate significant
information-theoretic separation between conditions, revealing that
isoform complexity patterns (captured by multi-q Tsallis entropy)
fundamentally differ between treatment groups across all q-dependent
scales (R0̆0e9 and Azad 2014).

``` r

# Calculate divergence: computes pairwise information-theoretic distance (Tsallis divergence) between conditions
# using bootstrap confidence intervals across all configured q-values from getConfig(analysis)$q_values
analysis <- calculate_divergence_s4(
    analysis)
```

``` r

# Compute effect sizes: significance_threshold = 0.05 filters genes to those with q*condition interaction p < 0.05;
# enrich_per_q_pattern = TRUE classifies each gene by its divergence pattern (Rare-driven/Balanced/Abundant-driven)
analysis <- effect_sizes_divergence_s4(
    analysis,
    significance_threshold = 0.05,
    enrich_per_q_pattern = TRUE
)
```

| Gene | Mean Divergence | Q-Pattern | D_rare | D_abundant | Ratio | LM adj. p-value |
|:---|---:|:---|---:|---:|---:|:---|
| CXCL12 | 0.2846 | Balanced | 0.3063 | 0.2712 | 1.13 | 3.0e-162 |
| THY1 | 0.1279 | Balanced | 0.1455 | 0.1142 | 1.27 | 7.3e-95 |
| ING3 | 0.0101 | Rare driven | 0.0120 | 0.0085 | 1.42 | 2.0e-75 |
| SNHG10 | 0.0714 | Rare driven | 0.0861 | 0.0591 | 1.46 | 7.2e-73 |
| LINC03040 | 0.2398 | Rare driven | 0.2813 | 0.2057 | 1.37 | 1.3e-55 |
| HDAC2 | 0.0525 | Rare driven | 0.0642 | 0.0425 | 1.51 | 1.9e-49 |
| ENSG00000274322 | 0.1506 | Rare driven | 0.1852 | 0.1206 | 1.54 | 3.3e-46 |
| MEF2A | 0.0298 | Rare driven | 0.0382 | 0.0225 | 1.70 | 9.5e-42 |
| RAP1GDS1 | 0.0139 | Rare driven | 0.0168 | 0.0114 | 1.47 | 5.2e-36 |
| PDE7A | 0.0419 | Rare driven | 0.0490 | 0.0360 | 1.36 | 1.1e-35 |

**Table 5 \| Top genes by linear model significance with effect sizes
and q-spectrum patterns.** Ranked by statistical significance (ascending
*P*-values, Benjamini-Hochberg *q*-value \< 0.05). Columns: gene
identifier; effect size (mean divergence across q-spectrum); pattern
classification (Rare-driven/Balanced/Abundant-driven); pattern
metrics-D_rare: median divergence for rare isoforms (q\<1), D_abundant:
median divergence for abundant isoforms (q\>1), Ratio: D_rare/D_abundant
(\>1.3 indicates rare-driven, \<0.77 indicates abundant-driven);
adjusted P-value for linear model interaction test. {.table .table
.table-striped .table-hover .table-condensed
style="margin-left: auto; margin-right: auto;"}

Interpretation of Pattern Types:

| Pattern Type | Signature | Biological Meaning |
|----|----|----|
| Rare driven | D(q=0.5) \>\> D(q=2.0) | Low-abundance isoforms are condition-specific; treatment group preferentially expresses rare transcripts not seen in control |
| Abundant driver | D(q=0.5) \<\< D(q=2.0) | High-abundance isoforms shift between conditions; treatment remodels the dominant transcript landscape without rare variants changing |
| Balanced | D(q=0.5) ~= D(q=2.0) | All isoforms shift proportionally; no preferential weighting to rare or abundant forms; suggests systematic rebalancing |

The effect sizes reveal which genes show the most information-theoretic
separation across the q-spectrum **between paired treatment groups**.
Genes with Tsallis divergence D \> 0.1 show substantial divergence,
indicating that entropy distributions differ fundamentally between
conditions across the full q-spectrum.

**Biological Interpretation:** Following Tsallis entropy theory and
information-theoretic principles (Tsallis 2006; Rényi 1961; Sason 2022b;
Hyndman and Athanasopoulos 2018), effect size quantification (Unknown
2020) genes with large D values are those where isoform complexity
distributions differ qualitatively between conditions when examined at
all sensitivity scales (q = rare -\> abundant isoforms). For example, a
gene might show high entropy at low q (emphasizing rare isoforms) in one
condition but low entropy at high q (emphasizing abundant isoforms) in
another, revealing scale-dependent regulatory mechanisms. These genes
are candidates for investigation of condition-specific splicing
architecture and dynamic isoform switching.

``` r

# Visualize the distribution of Tsallis divergence effect sizes across genes
# The red dashed line marks D = 0.1 (information-theoretic significance threshold)
# Generate divergence distribution plot (returns plot object)
# Optionally save to file for external use
plot_obj <- plot_divergence_distribution_s4(
    analysis, 
    threshold = 0.05
)

print(plot_obj)
```

![](TSENAT_files/figure-html/effect-size-plot-multiq-1.png)

``` r

# Visualize q-spectrum curves for top 4 genes by significance
# Plot is rendered directly by knitr (no file dependency)
p_multi <- plot_divergence_spectrum_s4(
    analysis, 
    n_genes = 4, 
    use_pvalue_ranking = TRUE
)

print(p_multi)
```

![\*\*Figure 5:\*\* Q-spectrum curves for top genes by significance.
Each line represents divergence as a function of q-value; shape reveals
biological pattern (rare-driven, balanced, or
abundant-driven).](TSENAT_files/figure-html/compare-q-spectra-plot-1.png)

**Figure 5:** Q-spectrum curves for top genes by significance. Each line
represents divergence as a function of q-value; shape reveals biological
pattern (rare-driven, balanced, or abundant-driven).

Interpreting the Q-Spectrum Curve:

- Shape matters more than single value: The **q-spectrum** curve is the
  complete biological story. A flat curve (balanced) vs. declining curve
  (rare driven) vs. rising curve (abundant driven) encode fundamentally
  different mechanisms.

- q=1 (KL divergence): The middle point corresponds to ordinary
  Kullback-Leibler divergence (Kullback and Leibler 1951; R0̆0e9 and Azad
  2014), where in the limit q-\>1 the Tsallis divergence reduces to
  standard KL divergence (R0̆0e9 and Azad 2014). This weights all
  isoforms equally. This is the “average” effect.

- q \< 1 (rare isoforms): Emphasizes how much low-abundance transcripts
  differ between groups. If D(q=0.5) \>\> D(q=1), the divergence is
  driven by changes in rare variants.

- q \> 1 (abundant isoforms): Emphasizes dominant transcripts. If D(q=2)
  \>\> D(q=1), major isoforms rebalance while rare ones stay similar.

- Bootstrap CI bands: Shaded region shows 95% confidence limits on per-q
  estimates. Narrow bands indicate robust estimates; wide bands suggest
  noisy data or small sample sizes.

We can also plot the global divergence spectrum across all genes. This
curve is the full biological signature of **isoform switching**-the
pattern encodes which abundance scales (rare vs. abundant) are affected.

``` r

# Visualize global divergence q-spectrum across all genes (aggregated)
# Plot is rendered directly by knitr (no file dependency)
p <- plot_divergence_spectrum_s4(
    analysis
)

print(p)
```

![\*\*Figure 6:\*\* Global divergence spectrum across all genes.
Aggregated q-spectrum pattern shows which abundance scales (rare vs.
abundant isoforms) are affected
genome-wide.](TSENAT_files/figure-html/visualize-q-spectrum-1.png)

**Figure 6:** Global divergence spectrum across all genes. Aggregated
q-spectrum pattern shows which abundance scales (rare vs. abundant
isoforms) are affected genome-wide.

Statistical Validation During Interpretation:

1.  Bootstrap CI validity: Per-q CIs should narrow as q increases
    (Efron, Bradley and Tibshirani, Robert J. 1993) (higher q =
    aggregation effect, less variance). If CIs grow with q, check for
    data quality issues.

2.  Monotonicity check: By **Tsallis entropy** theory (Tsallis 2006),
    entropy is monotone decreasing in q. Divergence should NOT show
    erratic increases with q. Small fluctuations are normal, but large
    spikes indicate numerical instability.

3.  **Comparison with genome-wide patterns**: Compute q-spectra for
    housekeeping genes (GAPDH, ACTB, etc.). These should show BALANCED
    patterns. If not, revisit normalization parameters.

### S4 Workflow: Unified Analysis via TSENATAnalysis Objects

TSENAT provides an integrated S4-based workflow for coordinated analysis
following three principles:

1.  **Configure once**: Use
    [`tsenat_config()`](https://gallardoalba.github.io/TSENAT/reference/tsenat_config.md)
    to specify analysis parameters (q-values, metadata columns)
2.  **Run pipeline**: Use
    [`tsenat()`](https://gallardoalba.github.io/TSENAT/reference/tsenat.md)
    to orchestrate all analysis steps in sequence
3.  **Access results**: Use S4 accessor methods for type-safe result
    retrieval

The
[`tsenat()`](https://gallardoalba.github.io/TSENAT/reference/tsenat.md)
function accepts a configured `TSENATConfig` object and orchestrates a
complete pipeline: diversity computation -\> quality control -\>
statistical testing -\> effect sizes (as specified by the `methods`
parameter). Parameter `verbose = TRUE` outputs progress messages showing
data filtering steps and method completion. This single-function
approach ensures consistency between preprocessing and analysis steps,
eliminating boilerplate and reducing risk of parameter mismatches.

**Accessor methods return types:** Most accessors return
`SummarizedExperiment` objects (following Bioconductor conventions)
containing results as assays and metadata:

- `diversity(analysis)` -\> SummarizedExperiment (genes x samples,
  multiple assays for each q-value)

- `lmResults(analysis)` -\> list containing data.frame of LM interaction
  results

- `jeoResults(analysis, q = 1.0)` -\> SummarizedExperiment of jackknife
  entropy outlier confidence intervals

- `jisResults(analysis, q = 1.0)` -\> list of jackknife isoform
  switching results

- `rankResults(analysis)` -\> data.frame of rank test results

- `divergence(analysis)` -\> SummarizedExperiment of pairwise divergence
  metrics

If a requested method was not run during
[`tsenat()`](https://gallardoalba.github.io/TSENAT/reference/tsenat.md)
orchestration, the accessor returns NULL-users should check results
before downstream use.

``` r

# Configure analysis parameters once (best practice for reproducibility)
tsenat_config <- tsenat_config(
  q_values = seq(0, 2, by = 0.05),
  sample_col = "sample",
  condition_col = "condition",
  subject_col = "paired_samples",
  nthreads = 2
)

## Build a complete `TSENATAnalysis` object
analysis <- build_analysis_s4(
  config = config,
  readcounts = readcounts,
  metadata = metadata_df,
  tx2gene = gff3_file, 
  tpm = tpm,
  effective_length = effective_length
)

# Run complete pipeline with one function call orchestrating all specified methods
# methods parameter: vector of analysis steps to execute ("diversity", "lm_interaction", "jackknife", "divergence")
# verbose = TRUE: prints progress messages showing gene filtering, step completion, and result summaries
analysis <- tsenat(analysis)

# Access results via consistent S4 methods (type-safe, mutually-exclusive access instead of nested list indexing)
diversity_results <- diversity(analysis)
lm_table <- lmResults(analysis)
jk_entropy <- jeoResults(analysis, q = 1.0)
jk_iso <- jisResults(analysis, q = 1.0)
rank_test <- rankResults(analysis)
divergence_table <- divergence(analysis)
```

For detailed step-by-step implementation, see the sections above
demonstrating each analysis component and its accessor patterns.

## Appendices

This main vignette is complemented by two comprehensive appendices:

#### Appendix A: Equivalence Validation - TSENAT vs SplicingFactory

Validates that TSENAT’s Shannon and Simpson entropy implementations are
mathematically equivalent to SplicingFactory. Includes:

- Mathematical proof of equivalence for q=1 (Shannon) and q=2 (Simpson)
- Benchmarking on TCGA BRCA RNA-seq data
- Detailed comparison tables and visualizations
- When to use TSENAT vs SplicingFactory

**[View Appendix
A](https://gallardoalba.github.io/TSENAT/articles/TSENAT_appendix_A.md)**

#### TSENAT Appendix B: Non-Parametric Validation of Linear Model Results via GAM and Rank-Based Methods

**Objective**: Validate q-value \* group interaction detection results
from linear models using complementary non-parametric statistical
approaches.

TSENAT’s default linear modeling approach (via
[`calculate_difference_s4()`](https://gallardoalba.github.io/TSENAT/reference/calculate_difference_s4.md))
provides parametric tests for detecting scale-dependent differences in
Tsallis entropy across q-values and experimental groups. To ensure
robustness and generalization, Appendix B presents two independent,
non-parametric alternatives that serve as validation methods:

**Approach 1: Generalized Additive Model (GAM) with ARIMA-Ordered Time
Series:**

- Treats q-values as explicit sequential measurements (ordered from low
  to high)
- Applies ARIMA(1,1,0) differencing to account for autocorrelation
  between adjacent q-points
- Fits smooth additive functions of q and group effects (non-parametric
  basis functions)
- Automatically detects heteroscedasticity and applies optimal weighting

**Approach 2: Rank-Based Friedman Test with Hochberg Correction:**

- Non-parametric alternative requiring no distributional assumptions
- Friedman test for detecting differences across multiple q-values
  within each group
- Accounts for paired/blocked structure in q-ordered measurements
- Hochberg step-up procedure for controlling family-wise error rate
  (FWER)
- References: Efron & Tibshirani (1993); Benjamini & Hochberg (1995)

Both methods incorporate principled advances to handle real-world
RNA-seq data: - **Heteroscedasticity detection and weighting**: Entropy
variance often increases with q-value; both approaches adaptively weight
measurements - **Boundary condition handling**: Low and high q-values
exhibit different variance properties; explicit boundary adaptation
improves reliability - **Automatic test selection**: Model diagnostics
trigger appropriate method switches based on data structure

**[View Appendix
B](https://gallardoalba.github.io/TSENAT/articles/TSENAT_appendix_B.md)**

### References

1.  Adami, C. (2004). “Information Theory in Molecular Biology.” *arXiv
    Preprint q-Bio/0405004*.

2.  Bajic, D. (2024). “Information Theory, Living Systems, and
    Communication Engineering.” *Entropy*, 26(5), 430.
    <https://doi.org/10.3390/e26050430>.

3.  Bartal, A., & Jagodnik, K. M. (2022). “Progress in and Opportunities
    for Applying Information Theory to Computational Biology and
    Bioinformatics.” *Entropy*, 24(7), 925.
    <https://doi.org/10.3390/e24070925>.

4.  Chanda, P., Costa, E., Hu, J., Sukumar, S., Van Hemert, J., &
    Walia, R. (2020). “Information Theory in Computational Biology:
    Where We Stand Today.” *Entropy*, 22(6), 627.
    <https://doi.org/10.3390/e22060627>.

5.  Cover, T. M., & Thomas, J. A. (2006). *Elements of Information
    Theory* (2nd ed.). Wiley-Interscience.

6.  Derian, N., Pham, H.-P., Nehar-Belaid, D., et al. (2022). “The
    Tsallis Generalized Entropy Enhances the Interpretation of
    Transcriptomics Datasets.” *PLOS ONE*, 17(4), e0266618.
    <https://doi.org/10.1371/journal.pone.0266618>.

7.  Furuichi, S. (2006). “Information Theoretical Properties of Tsallis
    Entropies.” *Journal of Mathematical Physics*, 47, 023302.
    <https://doi.org/10.1063/1.2165744>.

8.  Gandrillon, O., Gaillard, M., Espinasse, T., et al. (2021). “Entropy
    as a Measure of Variability and Stemness in Single-Cell
    Transcriptomics.” *Entropy*, 21(5), 450.

9.  Golomb, R., Yoles, M., Fishilevich, S., et al. (2026). “An
    Information Content Principle Explains Regulatory Patterns of Gene
    Expression Across Human Tissues.” *bioRxiv*, ahead of print.
    <https://doi.org/10.64898/2026.02.19.706555>.

&nbsp;

11. Jost, L. (2006). “Entropy and Diversity.” *Oikos*, 113(2), 363-375.
    <https://doi.org/10.1111/j.2006.0030-1299.14714.x>.

12. Renyi, A. (1961). “On Measures of Entropy and Information.”
    *Proceedings of the Fourth Berkeley Symposium on Mathematical
    Statistics and Probability*, 1, 547-561.

13. Seweryn, M. T., Pietrzak, M., & Ma, Q. (2020). “Application of
    Information Theoretical Approaches to Assess Diversity and
    Similarity in Single-Cell Transcriptomics.” *Computational and
    Structural Biotechnology Journal*, 18, 1830-1837.
    <https://doi.org/10.1016/j.csbj.2020.05.006>.

14. Shannon, C. E. (1948). “A Mathematical Theory of Communication.”
    *The Bell System Technical Journal*, 27(3-4), 379-423.

15. Simpson, E. H. (1949). “Measurement of Diversity.” *Nature*,
    163, 688. <https://doi.org/10.1038/163688a0>.

16. Tsallis, C. (1988). “Possible Generalization of Boltzmann-Gibbs
    Statistics.” *Journal of Statistical Physics*, 52(1), 479-487.
    <https://doi.org/10.1007/BF01016429>.

### Session Information

``` r

sessionInfo()
#> R version 4.5.2 (2025-10-31)
#> Platform: x86_64-conda-linux-gnu
#> Running under: Ubuntu 22.04.5 LTS
#> 
#> Matrix products: default
#> BLAS/LAPACK: /home/nouser/miniconda3/lib/libopenblasp-r0.3.30.so;  LAPACK version 3.12.0
#> 
#> locale:
#>  [1] LC_CTYPE=es_ES.UTF-8       LC_NUMERIC=C              
#>  [3] LC_TIME=de_DE.UTF-8        LC_COLLATE=es_ES.UTF-8    
#>  [5] LC_MONETARY=de_DE.UTF-8    LC_MESSAGES=es_ES.UTF-8   
#>  [7] LC_PAPER=de_DE.UTF-8       LC_NAME=C                 
#>  [9] LC_ADDRESS=C               LC_TELEPHONE=C            
#> [11] LC_MEASUREMENT=de_DE.UTF-8 LC_IDENTIFICATION=C       
#> 
#> time zone: Europe/Berlin
#> tzcode source: system (glibc)
#> 
#> attached base packages:
#> [1] stats4    stats     graphics  grDevices utils     datasets  methods  
#> [8] base     
#> 
#> other attached packages:
#>  [1] SummarizedExperiment_1.40.0 Biobase_2.70.0             
#>  [3] GenomicRanges_1.62.1        Seqinfo_1.0.0              
#>  [5] IRanges_2.44.0              S4Vectors_0.48.0           
#>  [7] BiocGenerics_0.56.0         generics_0.1.4             
#>  [9] MatrixGenerics_1.22.0       matrixStats_1.5.0          
#> [11] ggplot2_4.0.2               TSENAT_0.99.0              
#> [13] kableExtra_1.4.0           
#> 
#> loaded via a namespace (and not attached):
#>  [1] gtable_0.3.6        xfun_0.56           bslib_0.10.0       
#>  [4] htmlwidgets_1.6.4   lattice_0.22-9      vctrs_0.7.2        
#>  [7] tools_4.5.2         parallel_4.5.2      tibble_3.3.1       
#> [10] pkgconfig_2.0.3     pheatmap_1.0.13     Matrix_1.7-4       
#> [13] RColorBrewer_1.1-3  S7_0.2.1            desc_1.4.3         
#> [16] lifecycle_1.0.5     compiler_4.5.2      farver_2.1.2       
#> [19] stringr_1.6.0       textshaping_1.0.5   codetools_0.2-20   
#> [22] htmltools_0.5.9     sass_0.4.10         yaml_2.3.12        
#> [25] pkgdown_2.2.0       pillar_1.11.1       jquerylib_0.1.4    
#> [28] tidyr_1.3.2         BiocParallel_1.44.0 cachem_1.1.0       
#> [31] DelayedArray_0.36.0 abind_1.4-8         nlme_3.1-168       
#> [34] tidyselect_1.2.1    digest_0.6.39       stringi_1.8.7      
#> [37] dplyr_1.2.0         purrr_1.2.1         splines_4.5.2      
#> [40] labeling_0.4.3      cowplot_1.2.0       fastmap_1.2.0      
#> [43] grid_4.5.2          cli_3.6.5           SparseArray_1.10.9 
#> [46] magrittr_2.0.4      S4Arrays_1.10.1     withr_3.0.2        
#> [49] scales_1.4.0        rmarkdown_2.30      XVector_0.50.0     
#> [52] otel_0.2.0          ragg_1.5.0          memoise_2.0.1      
#> [55] evaluate_1.0.5      knitr_1.51          viridisLite_0.4.3  
#> [58] mgcv_1.9-4          rlang_1.1.7         Rcpp_1.1.1         
#> [61] glue_1.8.0          xml2_1.5.2          svglite_2.2.2      
#> [64] rstudioapi_0.18.0   jsonlite_2.0.0      R6_2.6.1           
#> [67] systemfonts_1.3.2   fs_1.6.7
```
