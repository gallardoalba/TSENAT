# TSENAT Feature Roadmap: Opportunities from the Literature Database

*Generated from `database/tsenat_papers.db` — 486 papers across 14
categories*

------------------------------------------------------------------------

## Tier 1: High-Impact, Directly Aligned (Near-Term)

### 1.1 Jensen-Tsallis Divergence (Symmetrized)

**Database support**: Papers I041 (Quantum Tsallis-Jensen-Shannon
divergence), I045, I060, I050

**What**: The current
[`calculate_divergence()`](https://gallardoalba.github.io/TSENAT/reference/calculate_divergence.md)
uses the asymmetric Tsallis divergence $`D_q(p\|r)`$. A symmetrized
Jensen-Tsallis divergence:

``` math
JT_q(p, r) = \frac{1}{2}\left[D_q\left(p \,\middle\|\, \frac{p+r}{2}\right) + D_q\left(r \,\middle\|\, \frac{p+r}{2}\right)\right]
```

would provide a proper distance metric (symmetric, non-negative,
satisfies triangle inequality for certain q). This is important
because: - Asymmetric divergence requires designating a “reference”
group — Jensen-Tsallis removes this requirement - Enables distance-based
clustering/ordination of samples without arbitrary reference choice -
Papers I041 and I050 contain the mathematical framework

**Effort**: Low — reuses existing divergence machinery; add a
`method='jensen'` option.

------------------------------------------------------------------------

### 1.2 Renyi Entropy as Alternative Parameterization

**Database support**: Papers I019, I059, I064

**What**: Rényi entropy is the monotonic transformation of Tsallis
entropy:

``` math
R_q = \frac{1}{1-q}\log\sum_i p_i^q
```

while Hill numbers $`D_q = \exp(R_q)`$ are already available via
`what='D'`. Rényi entropy has better additivity properties for
independent systems and is the standard in ecology (Jost 2006). Adding a
`what='R'` option alongside `'S'` (Tsallis) and `'D'` (Hill) would
increase interoperability with the ecology literature.

**Effort**: Very low — one-line transformation:
$`R_q = \log(1 + (1-q)S_q) / (1-q)`$

------------------------------------------------------------------------

### 1.3 Power Analysis and Sample Size Planning

**Database support**: 7 dedicated power-analysis papers (S063-S067,
C106-C109) plus 50+ method references

**What**: A `calculate_power()` function that simulates power for
detecting Q×Condition interactions at specified effect sizes. The
database contains: - `RnaSeqSampleSize` (C108): Sample size for RNA-seq
with NB distribution - `PROPER` (C107): Prospective power evaluation for
RNA-seq - `ssizeRNA` (C106): Power and sample size for RNA-seq
differential expression - `ERSSA` (C109): Empirical RNA-seq sample size
analysis

An entropy-specific power analysis would: 1. Simulate isoform count data
with known switching patterns 2. Apply the TSENAT pipeline (diversity →
rank test) 3. Estimate power as a function of sample size, effect size,
and q-range

**Effort**: Medium — requires simulation framework but is
self-contained.

------------------------------------------------------------------------

### 1.4 Optimal q-Value Selection

**Database support**: Papers on model selection (Cross-validation: 60
refs, AIC/BIC)

**What**: Currently users must manually specify q-values. An automated
procedure that identifies the most informative q-values for a given
dataset would: 1. Compute entropy profiles across a dense q-grid 2. Use
functional PCA or mutual information to identify q-values where entropy
curves maximally discriminate between conditions 3. Return recommended
q-values for subsequent analysis

This addresses a real user pain point: “what q-values should I use?”

**Effort**: Medium — statistical core is straightforward; needs good
heuristics.

------------------------------------------------------------------------

## Tier 2: Methodological Extensions (Medium-Term)

### 2.1 Bayesian Entropy Inference

**Database support**: 9 Bayesian Statistics papers + 2,159 Bayesian
method references

**What**: Replace point estimates of entropy with full posterior
distributions: 1. Model isoform counts as Dirichlet-multinomial with a
prior 2. Compute posterior distribution of $`S_q`$ (not just point
estimate) 3. Provides credible intervals, posterior probability of
$`S_q^{(A)} > S_q^{(B)}`$, and Bayes factors for differential entropy

The database contains: - BY007: Bayesian nonparametric isoform
discovery - Papers on Bayesian posterior predictive checks (4 papers) -
Bayesian credible intervals and hypothesis testing

This is a genuine methodological advance — no existing tool provides
Bayesian uncertainty quantification for multi-q entropy analysis.

**Effort**: High — requires MCMC or variational inference; significant
implementation effort.

------------------------------------------------------------------------

### 2.2 Jensen-Shannon Divergence for Multi-Group Comparisons

**Database support**: I041 (Quantum Tsallis-Jensen-Shannon), general JS
divergence papers

**What**: Extend divergence to handle **more than 2 groups**. The
generalized Jensen-Shannon divergence for $`K`$ groups with mixture
weights $`\pi_k`$:

``` math
JS_q(\{p_k\}) = \sum_k \pi_k D_q(p_k \| \bar{p}), \quad \bar{p} = \sum_k \pi_k p_k
```

This enables ANOVA-like divergence analysis for multi-condition
experiments (e.g., drug dose-response, time series, multiple tissue
types).

**Effort**: Medium — generalization of existing machinery.

------------------------------------------------------------------------

### 2.3 Gene Regulatory Network Entropy

**Database support**: C117 (Gene regulatory network inference), papers
on Network Entropy, network modularity

**What**: Tsallis entropy applied to **gene regulatory networks** rather
than isoform proportions. The “entropy of signaling” framework (papers
in database) computes entropy of information flow through regulatory
networks. This would complement the isoform-level analysis with a
systems-level perspective.

**Effort**: High — requires network inference + new entropy formulation.

------------------------------------------------------------------------

### 2.4 Effect Size Standardization and Meta-Analysis Support

**Database support**: Papers on Effect Size Calculation (7), Effect Size
Conversion, Standardized Effect Sizes, Rank-Biserial Effect Size

**What**: Standardize TSENAT effect sizes ($`\eta^2`$, slope
differences, divergence magnitudes) into common metrics (Cohen’s d,
Hedges’ g, log-odds ratios) to enable: 1. Cross-study meta-analysis of
isoform diversity 2. Power analysis with standardized effect sizes 3.
Comparison with DGE/DTE effect sizes from DESeq2/edgeR

**Effort**: Low-Medium — conversion formulas are well-established.

------------------------------------------------------------------------

## Tier 3: Data Type Expansions (Longer-Term)

### 3.1 Single-Cell RNA-seq Support

**Database support**: S066 (scRNA-seq power analysis), C117 (scRNA-seq
gene regulatory networks), ISO012 (IsoformSwitchAnalyzeR for scRNA-seq),
BY018

**What**: Adapt TSENAT for single-cell data. Key challenges: 1.
**Sparsity**: Many isoforms have zero counts in most cells — the
pseudocount regularization needs cell-type-aware tuning 2.
**Zero-inflation**: Standard entropy formulas behave poorly with excess
zeros 3. **Cell-type heterogeneity**: Need to account for cell-type
composition

The database contains methods for scRNA-seq power analysis (S066) and
existing isoform tools that handle single-cell data (ISO012).

**Effort**: Very High — requires new normalization, zero-inflation
models, and cell-type deconvolution.

------------------------------------------------------------------------

### 3.2 Spatial Transcriptomics Integration

**Database support**: C114 (Information-Theoretic Methods in Spatial
Transcriptomics, 2026)

**What**: Extend Tsallis entropy to spatial contexts — compute **local
entropy** in spatial neighborhoods to detect regions of isoform
diversity change. Paper C114 provides the theoretical framework.

**Effort**: Very High — new data structure, spatial statistics,
visualization.

------------------------------------------------------------------------

### 3.3 Long-Read Sequencing (Iso-Seq, Nanopore)

**Database support**: ISO012 (IsoformSwitchAnalyzeR v2 for long-read and
scRNA-seq)

**What**: Long-read technologies provide **full-length isoform**
resolution, eliminating the transcript assembly problem. TSENAT’s
entropy framework would be directly applicable to full-length isoform
counts, with higher accuracy since isoform identification is
unambiguous. Key additions: 1. Direct import of Iso-Seq/Nanopore isoform
counts 2. Full-length isoform-aware entropy normalization 3. Detection
of novel isoform switching events invisible to short-read data

**Effort**: Medium — mostly import/format work; entropy machinery is
ready.

------------------------------------------------------------------------

## Tier 4: Quality-of-Life and Ecosystem

### 4.1 Interactive Visualization (Shiny/plotly)

Add an interactive `explore()` Shiny app for: - Interactive q-spectrum
exploration - Gene-level drill-down from multi-gene views - Dynamic
p-value threshold adjustment

### 4.2 Cross-Tool Benchmarking Framework

Add
[`calculate_concordance()`](https://gallardoalba.github.io/TSENAT/reference/calculate_concordance.md)
enhancements to benchmark TSENAT against: - DGE tools: DESeq2, edgeR,
limma-voom - DTU tools: DEXSeq, DRIMSeq, SUPPA2, IsoformSwitchAnalyzeR -
Splicing diversity: SplicingFactory

### 4.3 Containerized Reproducibility

Docker/Singularity images with all dependencies pre-installed, plus a
Nextflow pipeline for automated TSENAT workflows from FASTQ → results.

### 4.4 q-Profile Database

A curated database of reference entropy q-profiles across tissues,
species, and conditions (like recount3 but for isoform diversity),
enabling: - “Is this gene’s q-profile abnormal for this tissue?” -
Reference-based normalization - Cross-study meta-analysis

------------------------------------------------------------------------

## Prioritization Matrix

| Feature | Impact | Effort | Database Support | Priority |
|----|----|----|----|----|
| Jensen-Tsallis divergence | High | Low | Strong (4 papers) | ⭐⭐⭐ |
| Rényi entropy option | Medium | Very Low | Moderate (3 papers) | ⭐⭐⭐ |
| Power analysis | High | Medium | Strong (7+ papers) | ⭐⭐⭐ |
| Optimal q-value selection | High | Medium | Moderate (60 CV refs) | ⭐⭐ |
| Effect size standardization | Medium | Low-Medium | Strong (7+ papers) | ⭐⭐ |
| Bayesian entropy inference | Very High | High | Strong (9+ papers) | ⭐⭐ |
| Multi-group JS divergence | Medium | Medium | Moderate | ⭐⭐ |
| Single-cell support | High | Very High | Moderate (4 papers) | ⭐ |
| Spatial transcriptomics | Medium | Very High | Limited (1 paper) | ⭐ |
| Long-read support | Medium | Medium | Limited (1 paper) | ⭐ |
| Network entropy | Medium | High | Moderate | ⭐ |
| Interactive viz (Shiny) | Medium | Medium | N/A | ⭐ |
| Cross-tool benchmarking | Medium | Low | N/A | ⭐ |

------------------------------------------------------------------------

*Generated: 2026-07-25 from tsenat_papers.db (486 papers)*
