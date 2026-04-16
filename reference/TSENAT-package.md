# TSENAT: Tsallis Entropy for Systems Evolutionary Analysis and Network Analysis Tool

TSENAT is a Bioconductor package that quantifies transcriptomic
complexity at the isoform level using Tsallis entropy, a generalization
of Shannon entropy that captures rare-event dynamics and scale-dependent
information structure. The package integrates diversity quantification,
statistical resampling, linear mixed models, and information-theoretic
divergence metrics into a unified analytical framework for comparative
transcriptomics.

## Scientific Rationale

Classical Shannon entropy (a special case of Tsallis with q=1) treats
all probability states equally. Tsallis entropy with order parameter
\\q\\ weights rare and common transcripts differently:

- \\q \< 1\\: Emphasizes rare isoforms, capturing tail behavior

- \\q = 1\\: Shannon entropy (equipartition weighting)

- \\q = 2\\: Simpson entropy (dominance weighting, less sensitive to
  rare species)

- \\q \to \infty\\: Reciprocal of the maximum transcript fraction

This multiscale perspective reveals condition-specific dependencies in
splicing architecture and detects genes with asymmetric isoform
distributions that traditional metrics miss.

## Workflow Overview

TSENAT implements a complete analytical pipeline:

1.  **Diversity Quantification**: Compute Tsallis entropy H_q across a
    range of q-values (0 to \\\infty\\) for each gene-condition group.
    Returns Hill numbers (\\D_q = \exp(H_q)\\) for interpretable
    effective isoform counts.

2.  **Jackknife Resampling**: Generate confidence intervals via
    leave-one-out resampling. Identifies genes with statistically
    significant entropy outliers and condition-biased isoform usage
    (isoform switching).

3.  **Linear Mixed Models**: Tests gene-by-condition interaction effects
    on diversity landscapes using nlme, MASS, or geepack backends.
    Supports paired designs and random effect structures.

4.  **Divergence Analysis**: Computes Tsallis divergence (generalized
    Kullback-Leibler) between condition groups. Detects differential
    transcript complexity via effect size estimation.

5.  **Visualization**: Generates publication-ready plots: q-curves,
    volcano plots, heatmaps, violin distributions, and dimension-reduced
    embeddings.

## Core Functions

- [`TSENAT`](https://gallardoalba.github.io/TSENAT/reference/TSENAT.md):

  Main orchestration function. Executes complete pipeline with
  configurable method backends and parallelization.

- [`TSENAT_config`](https://gallardoalba.github.io/TSENAT/reference/TSENAT_config.md):

  Configure analysis parameters: q-values, condition grouping, sample
  pairing, and method selection.

- [`calculate_diversity`](https://gallardoalba.github.io/TSENAT/reference/calculate_diversity.md):

  Compute Tsallis entropy and Hill numbers. Vectorized over genes and
  q-values with numerical stability safeguards.

- [`calculate_jeo`](https://gallardoalba.github.io/TSENAT/reference/calculate_jeo.md):

  Resampling-based confidence intervals for diversity. Identifies
  outlier genes and validates q-value signal.

- [`calculate_lm`](https://gallardoalba.github.io/TSENAT/reference/calculate_lm.md):

  Statistical testing for gene-by-condition interactions. LMM, GAM, and
  rank-based approaches.

- [`calculate_divergence`](https://gallardoalba.github.io/TSENAT/reference/calculate_divergence.md):

  Tsallis divergence between groups with effect size (Cohen's d,
  rank-biserial) and hypothesis testing.

- [`filter_analysis`](https://gallardoalba.github.io/TSENAT/reference/filter_analysis.md):

  Filter genes by diversity, effect size, or interaction significance
  (adjustable FDR, Benjamini-Hochberg).

## Statistical Methods

- **Diversity Metrics**:

  Tsallis Entropy: \\H_q = \frac{1}{q-1} \left(1 - \sum\_{i}
  p_i^q\right)\\ Hill Numbers: \\D_q = \left(\sum\_{i}
  p_i^q\right)^{\frac{1}{1-q}}\\

- **Divergence**:

  Tsallis Divergence: \\D_q(p\|\|r) = \frac{1}{q-1}\left(1 - \sum\_{i}
  p_i^q r_i^{1-q}\right)\\ Kullback-Leibler distance (q=1 limit):
  \\D\_{KL}(p\|\|r) = \sum_i p_i \log(p_i/r_i)\\

- **Resampling**:

  Jackknife leave-one-out, block bootstrap (paired data), and
  percentile/BCA confidence intervals with accuracy enhancement for
  skewed distributions.

- **Interaction Testing**:

  LMM: \\H\_{g,c} = \beta_0 + \beta_g(\text{gene}) +
  \beta_c(\text{condition}) + \beta\_{gc}(\text{gene} \times
  \text{condition}) + \epsilon\\ With random effects for paired/nested
  designs.

## Key Features

- **Multiscale q-parameter sweep**: Detects biological signals at
  different probability scales without multiple testing correction for
  q-values.

- **Robust resampling**: Leave-one-out jackknife + block bootstrap
  handles paired designs and maintains correlation structure.

- **Multiple method backends**: nlme, MASS, geepack, or rank-based
  (Scheirer-Ray-Hare, sign test) for model flexibility.

- **FDR-corrected testing**: Benjamini-Hochberg adjustment across genes
  and q-values.

- **Production-grade visualization**: ggplot2-based plots optimized for
  large gene sets and publication quality.

- **Reproducibility tracking**: Automatic logging of function calls,
  parameters, timestamps, and package version for audit trails.

- **S4 object-oriented design**: Unified
  [`TSENATAnalysis`](https://gallardoalba.github.io/TSENAT/reference/TSENATAnalysis.md)
  container ensures metadata integrity throughout the pipeline.

## Data Input & Output

- **Input**:

  SummarizedExperiment object with gene-by-sample count matrix (raw or
  normalized) and sample metadata (condition, pairing, batch).

- **Output**:

  TSENATAnalysis S4 object containing: (1) Diversity values across
  q-parameter range, (2) Jackknife CIs and switching indicators, (3)
  LMM/interaction statistics with p-values and effect sizes, (4)
  Divergence metrics with hypothesis test results, (5) Cached
  publication-ready plots, (6) Reproducibility metadata.

## Interpreter Notes

- **q-value sweeps reveal biological complexity**: Genes showing entropy
  peaks at intermediate q-values (0.5-1.5) suggest balanced isoform
  distributions; peaks at q \to \infty indicate dominant-isoform
  architectures.

- **Interaction q-profiles**: q-dependent gene-by-condition interactions
  signal condition-specific isoform switching; flat profiles indicate
  constitutive splicing.

- **Divergence as biological distance**: Tsallis divergence quantifies
  transcript composition distance; high divergence (D_q \> 0.5)
  indicates distinct splicing programs.

- **Statistical interpretation**: Jackknife CIs that exclude the point
  estimate (valid for skewed distributions) occur with bounded
  entropies; not a software error.

## References

Tsallis entropy and information theory:

- Tsallis, C. (1988). Possible generalization of Boltzmann-Gibbs
  statistics. *Journal of Statistical Physics*, 52(1-2), 479-487.

- Furuichi, S. (2006). Information theoretical properties of Tsallis
  entropies. *Journal of Mathematical Physics*, 47(2), 023302.

Hill numbers and diversity indices:

- Hill, M. O. (1973). Diversity and evenness: A unifying notation and
  its consequences. *Ecology*, 54(2), 427-432.

- Chao, A., Chiu, C. H., & Jost, L. (2014). Unifying species diversity,
  phylogenetic diversity, functional diversity, and related similarities
  into a single framework. *Annual Review of Ecology, Evolution, and
  Systematics*, 45, 297-324.

Isoform switching and transcriptomics:

- Vitting-Seerup, K., et al. (2017). IsoformSwitchAnalyzeR enables
  robust detection of isoform switches and novel isoforms in the human
  transcriptome from long-read cDNA sequencing data. *Genome Biology*,
  18(1), 122.

## S4 Container

[`TSENATAnalysis-class`](https://gallardoalba.github.io/TSENAT/reference/TSENATAnalysis-class.md)
— Central container unifying all analysis components. Access results via
slots: `@diversity_results`, `@lm_results`, `@jackknife_results`,
`@divergence_results`, or accessor methods
[`se`](https://gallardoalba.github.io/TSENAT/reference/TSENATAnalysis-se.md),
`metadata<-`, and
[`results`](https://gallardoalba.github.io/TSENAT/reference/results.md).

## See also

[`TSENAT`](https://gallardoalba.github.io/TSENAT/reference/TSENAT.md)
for running the complete analysis pipeline.
[`se`](https://gallardoalba.github.io/TSENAT/reference/TSENATAnalysis-se.md),
`metadata<-`, and
[`results`](https://gallardoalba.github.io/TSENAT/reference/results.md)
for accessor functions.
[`TSENATAnalysis-class`](https://gallardoalba.github.io/TSENAT/reference/TSENATAnalysis-class.md)
for the S4 object structure.
[`build_analysis`](https://gallardoalba.github.io/TSENAT/reference/build_analysis.md)
for creating TSENATAnalysis objects.

## Author

**Maintainer**: Cristóbal Gallardo <gallardoalba@pm.me>
([ORCID](https://orcid.org/0000-0002-5752-2155))
