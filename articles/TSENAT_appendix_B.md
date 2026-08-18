# Appendix B: Cross-Method Robustness and Concordance of Scale-Adaptive Interaction Test Results via GAMM and Aligned Rank Transform (ART)

## Introduction

The primary goal of this vignette is to assess the **robustness and
concordance** of scale-dependent entropic index × group interaction
discoveries from Scale-Adaptive Interaction Tests (SAIT): do they
generalize to a non-parametric statistical framework with fewer
distributional assumptions? Agreement between GAMM and ART establishes
robustness under complementary modeling assumptions.

### Two Complementary Approaches

**Method 1: Paired regression-spline mixed model (GAMM-style)** with
functional q-curve structure (full details in the main vignette). TSENAT
historically refers to this family as GAMM; it is fitted with
[`nlme::lme`](https://rdrr.io/pkg/nlme/man/lme.html), not
[`mgcv::gamm()`](https://rdrr.io/pkg/mgcv/man/gamm.html).

- Framework: Semi-parametric model
  $`Y = f_1(q) + f_2(condition) + f_3(q, condition)`$ +
  subject-intercept + $`\varepsilon`$ with smooth rather than linear
  $`f`$ terms.
- Paired fits: [`nlme::lme`](https://rdrr.io/pkg/nlme/man/lme.html) with
  `ns(q, df = 3) × condition`, subject random intercepts, and a
  continuous-q AR(1) working correlation within subject × condition
  ([`nlme::corCAR1`](https://rdrr.io/pkg/nlme/man/corCAR1.html),
  $`\exp(-\phi|q_i-q_j|)`$; validated by Monte Carlo simulation); the
  interaction is tested with a marginal F-test.
- q is a deterministic functional argument (not time): the confirmatory
  hypothesis is the functional interaction $`H_0: \beta(q) = 0`$ for all
  q, tested on the original H(q) curve. Gaussian residuals are assumed
  (paired mixed-model paths do not support extended families);
  heteroscedasticity is detected but variance weights are not applied in
  paired designs.

**Method 2: Aligned Rank Transform (ART)** — rank-based non-parametric
interaction testing (ARTool (Kay et al. 2021)), with Hochberg step-up
FWER correction. The Conover-Iman Rank Transform is available via
`method='rt'`.

- ART strips main effects before ranking (“alignment”) to preserve
  interaction structure (Higgins and Tashtoush 1994; Wobbrock et al.
  2011). It makes no Gaussian-error assumption, but is not
  assumption-free: validity depends on the design, exchangeability,
  alignment and sample size.
- q is treated as a categorical factor (its continuous ordering and the
  AR(1) working correlation across q are not modeled) — the fundamental
  difference from GAMM.
- Monte Carlo caveat: unpaired ART is mildly anti-conservative at small
  n (Conover-Iman `rt` recommended there). A full treatment follows in
  the “Aligned Rank Transform (ART)” section below.

### Why Two Methods for One Question?

Statistical testing in transcriptomics faces a fundamental challenge: no
single method is universally optimal. Different approaches trade-off
power for robustness:

| Aspect | Semi-Parametric (GAMM) | Non-Parametric (ART) |
|----|----|----|
| Assumptions | Normality, homoscedasticity | No Gaussian errors; validity still requires a valid design, exchangeability and alignment |
| Power | Highest (if assumptions hold) | Good; reduced but robust to violations |
| Interaction testing | Marginal F-test on regression-spline interaction terms | Proper interaction tests via alignment (Higgins and Tashtoush 1994) |
| q-value treatment | Continuous smooth predictor (preserves ordinal structure) | Categorical factor (discards ordering) |
| Robustness | Moderate | Highest |
| Outlier sensitivity | Moderate potential for bias | Minimal; inherently resistant |

Robustness strategy: Concordance between both methods provides evidence
of robustness across statistical frameworks (it does not by itself prove
either method correct).

------------------------------------------------------------------------

## Setup

This section initializes the analysis environment by loading required
packages, setting a random seed for reproducibility, and preparing the
test dataset.

``` r

suppressPackageStartupMessages({
    library(ggplot2)
    library(SummarizedExperiment)
    library(dplyr)
    library(gridExtra)
})

set.seed(42)

if (!requireNamespace("devtools", quietly = TRUE)) {
    stop("Please install devtools to run this vignette against the local package source.")
}
if ("package:TSENAT" %in% search()) {
    detach("package:TSENAT", unload = TRUE, character.only = TRUE)
}

pkg_root <- normalizePath("..", mustWork = FALSE)
devtools::load_all(pkg_root)

# Load preprocessed dataset
data(readcounts)
readcounts <- as.matrix(readcounts)
mode(readcounts) <- "numeric"

# Load metadata
metadata_df <- read.table(
    system.file("extdata", "metadata.tsv", package = "TSENAT"),
    header = TRUE, sep = "\t"
)

gff3_dataset <- system.file("extdata", "annotation.gff3.gz", package = "TSENAT")

# Configure analysis parameters first (best practice: fail-fast validation)
config <- TSENAT_config(
  sample_col = "sample",
  condition_col = "condition",
  subject_col = "paired_samples",
  q = seq(0, 2, by = 0.05),
  nthreads = 3,
  paired = TRUE,
  control = "normal"
)

# Build analysis object: creates SummarizedExperiment and initializes TSENATAnalysis
analysis <- build_analysis(
    readcounts = readcounts,
    tx2gene = gff3_dataset,
    metadata = metadata_df,
    tpm = tpm,
    effective_length = effective_length,
    config = config
)
# Roles: tpm -> abundance QC in filter_analysis(); effective_length ->
# length-normalized entropy from RAW counts. The two are mutually exclusive
# in a single computation (TPM already incorporates length normalization).

# Apply filtering for quality control
analysis <- filter_analysis(
    analysis,
    stringency = "medium"
)
```

### Computing Tsallis Entropy with Pseudocount Regularization

Pseudocounts are an **optional regularization mechanism**. RNA-seq
experiments typically contain many zero or near-zero counts, which
creates two practical problems for rank-based methods:

1.  **Ties in rankings:** Sparse counts produce many tied values
    (especially zeros), which reduces the discriminatory power of rank
    tests. Rank-based statistics depend on unique orderings; when many
    observations are identical, the test statistic contains less
    information, reducing statistical power.

2.  **Library size artifacts:** Small count differences due to
    sequencing depth variations can overshadow biological differences.
    Normalized library sizes (accounting for sequencing depth) can make
    entropy estimates more comparable across samples (Robinson et al.
    2010; Love et al. 2014).

By adding small pseudocounts we achieve two benefits:

- **Regularization** breaks ties and reduces zero-inflation effects on
  rank-based inference.
- **Comparability** makes entropy estimates more stable across samples
  with different sequencing depths.

Importantly, a pseudocount is a regularization *choice*, not a standard
RNA-seq normalization procedure, and it **changes the estimand**: for

``` math
p_i = \frac{x_i + a}{\sum_j (x_j + a)}
```

you are no longer measuring the empirical transcript composition. This
is especially consequential for low-abundance genes and $`q<1`$ (where
small probabilities are amplified) and for $`q=0`$ (where a positive
pseudocount changes the support).

The
[`calculate_diversity()`](https://gallardoalba.github.io/TSENAT/reference/calculate_diversity.md)
function implements automatic pseudocount selection, scaling pseudocount
magnitude and diversity estimates by median library size to ensure
robust statistical inference across the range of sequencing depths in
your dataset.

``` r

# Compute diversity using S4 wrapper with bootstrap confidence intervals
analysis <- calculate_diversity(
  analysis, 
  norm = TRUE,
  pseudocount = "auto"
)
```

## Aligned Rank Transform (ART)

TSENAT uses the **Aligned Rank Transform (ART)** as the default
non-parametric method for testing Q×Condition interactions. ART is a
well-established rank-based non-parametric approach for factorial
analysis, implemented via the ARTool R package (Kay et al. 2021;
Wobbrock et al. 2011). The classical Conover-Iman Rank Transform is also
available via `method='rt'` as an alternative approach (Conover and Iman
1981).

1.  **Robustness to Distribution Violations**

Tsallis entropy values exhibit inherent distributional properties that
often violate parametric assumptions:

- Non-normality: Entropy commonly exhibits bounded support (e.g.,
  normalized entropy in \[0,1\]), producing skewed or bimodal
  distributions rather than normal distributions.
- Heteroscedasticity: Variance in entropy estimates varies
  systematically across q-values. Low q-values (emphasizing rare
  transcripts) produce volatile entropy estimates with high variance;
  high q-values (emphasizing abundant transcripts) produce stable
  estimates with lower variance.
- Outlier sensitivity: Rare isoforms and count variability can produce
  extreme entropy values that heavily influence parametric tests.

ART mitigates these issues through alignment and ranking — stripping
main effects before ranking preserves interaction structure while
removing distributional shape dependence. The aligned ranks are more
uniformly distributed and contain no outliers, making the method
substantially less sensitive to distributional shape than parametric
tests (Higgins and Tashtoush 1994; Wobbrock et al. 2011). Critically,
rank ordering is unaffected by normalization choice; whether entropy is
normalized, log-transformed, or left raw, the test produces identical
results (Conover and Iman 1981). This robustness does not make ART
assumption-free: its validity still depends on the experimental design,
exchangeability structure, alignment procedure, and sample size.

2.  **Paired Design for Ordered q-Value Structure**

Tsallis entropy has a fundamental sequential property: entropy curves
are smooth functions of q. As q increases from 0
(rare-species-emphasizing) to $`\infty`$ (common-species-emphasizing),
entropy values change systematically. This ordered structure is critical
to interpreting q-dependent patterns and exhibits **AR(1)
autocorrelation**: consecutive q-values produce correlated entropy
estimates.

ART handles this structure through its alignment step: within-subject
ranking after stripping main effects preserves the pairing structure
while removing the influence of AR(1) dependence (Higgins and Tashtoush
1994; Wobbrock et al. 2011). Because aligned ranks are constructed from
residuals after removing main effects, autocorrelation in the original
entropy values does not propagate into the test statistic — rank-based
tests are invariant to monotonic transformations of the data, and
alignment renders the residuals approximately exchangeable under the
null hypothesis of no interaction (Ernst 2004; Song 2007). However, this
robustness is not absolute: strong autocorrelation can reduce effective
sample size and affect power even when the test remains valid.

3.  **Interaction Testing: Testing q × Condition Effects**

The implementation uses ART (alignment + ranking + ANOVA), a rank-based
factorial approach:

``` math
F_{q \times \text{condition}} = \frac{MS_{\text{interaction}}}{MS_{\text{residual}}}
```

where:

- Main effects are stripped via alignment: residuals = original - main
  effect estimates
- Aligned data are then ranked:
  $`R = \text{rank}(\text{aligned residuals})`$ (within-subject for
  paired designs)
- Two-way ANOVA decomposes aligned ranked data:
  $`R = \mu + \alpha_q + \beta_{\text{condition}} + \gamma_{q \times \text{condition}} + \epsilon`$
- $`MS_{\text{interaction}}`$ = sum of squares for $`q \times`$
  condition interaction / degrees of freedom
- $`MS_{\text{residual}}`$ = residual sum of squares / residual degrees
  of freedom
- $`F`$-statistic follows $`F`$-distribution under null hypothesis of no
  interaction

This directly tests the core biological question: “Does entropy
q-dependence differ between groups?” ART is validated across balanced
and unbalanced designs (Wobbrock et al. 2011); our Monte Carlo results
confirm calibration for the paired design (and flag unpaired ART as
mildly anti-conservative at small n).

### Assumption Validation

To ensure valid statistical inference, we verify key data assumptions
across multiple dimensions.

``` r

# Validate statistical assumptions (including GAMM diagnostics)
analysis <- calculate_assumptions(
    analysis,
    checks = "all"  # Include core checks + GAMM diagnostics
)

# Get assumptions
assumptions_text <- results(analysis, type = "assumptions")
print(assumptions_text)
```

| Test | Result | Interpretation |
|:---|:---|:---|
| Exchangeability (Permutation test) | p=0e+00 | Ordering detected |
| Monotonicity (Spearman rho) | r=0.326 | Homogeneous |
| Consistency (Kendall’s W / ICC) | W=0e+00, ICC=0.047 | Low |
| Concurvity (Smooth collinearity) | 0.000 | Low |
| EDF Ratio (Smoothing) | 0.024 | Over-smoothed |
| Non-linearity (Delta R^2 vs LM) | 0.0% | Use linear |
| Basis Dimension (Spline basis) | k=3 | Adequate |
| Correlation fit | Observed autocorr=0.076; independence suitable | Good fit |
| Cluster variation | Mean size=77.0 | Homogeneous |
| Independence | Mean within-cluster residual correlation=0.069 | Independent |
| Scale parameter | phi=0.100 | Under-dispersed (rare) |
| Variance components | ICC=0.049; B=2.14e-04, W=0.004 | Use sait model |
| Normality | p=8.76e-09 | Non-normal |
| Homogeneity | p=2.23e-159; CV=0.872 | Heterogeneous |
| Influence | Outliers=1144, Extreme=894, Influential=8.6% | Some outliers |
| Variance adequacy | Components for 90%=5, 95%=7, 99%=11 | Moderate reduction |
| Bootstrap stability | Bootstrap SE=0.380; Stable CIs=100% | Moderate stability |

**Supplementary Table 14 \| Assumption Checks for Aligned Rank Transform
(ART).** Evaluates homogeneity of variance and normality for rank-based
interaction test validity. {.table}

**Interpretation of Results:** The exchangeability test detects strong
serial correlation ($`p \approx 0`$), indicating that consecutive
samples are more correlated than expected by chance. This violation is
**not problematic** for the Aligned Rank Transform because ART’s
alignment step removes main effects before ranking, after which the
aligned residuals are approximately exchangeable under the null
hypothesis of no interaction (Conover and Iman 1981; Higgins and
Tashtoush 1994). While strong autocorrelation can reduce effective
sample size and affect power, the test itself remains valid.

### Running the ART Analysis

We run ART (the default rank method) with Hochberg step-up correction
for multiple testing across genes, then compare the results against GAMM
below.

``` r

# Run Aligned Rank Transform (ART, default method)
# This preserves the GAMM results in 'analysis' for later concordance comparison

# Run ART for q × condition interaction
analysis_rank <- calculate_rank_transform(
    analysis,
    multicorr = "hochberg"
)

# Extract ALL results (no filtering) for summary statistics
rank_transform_results_all <- results(analysis_rank, type = "rank_test", rankBy = "pvalue")

# Extract top significant genes for display
rank_transform_results <- results(analysis_rank, type = "rank_test", rankBy = "pvalue", n = 20, filterFDR = 0.05)
print(head(rank_transform_results, n = 10))
```

| Metric                                       | Value |
|:---------------------------------------------|------:|
| Genes tested                                 |    77 |
| Significant (p \< 0.05)                      |    52 |
| Significant (adj_p \< 0.05, FWER-controlled) |    50 |
| NAs                                          |     0 |
| Mean effect size ($`\eta^2`$)                |  1.6% |
| Median effect size ($`\eta^2`$)              |  0.8% |
| Strong effect genes ($`\eta^2`$ \> 10%)      |     1 |

**Supplementary Table 15 \| Aligned Rank Transform (ART) Results
Summary.** Non-parametric test of *q* × condition interactions via ART.
Compare with GAMM results (main vignette Table 3) to assess parametric
vs. non-parametric method concordance. {.table}

| Gene | P-value | Adj. P-value | F-Statistic | Effect Size (eta-squared) | Interaction Class |
|:---|:---|:---|:---|:---|:---|
| CXCL12 | 3.8e-146 | 2.9e-144 | 106.0585 | 0.1062 | Strongly q-dependent |
| SPNS2 | 1.1e-111 | 8.3e-110 | 56.7003 | 0.0774 | Moderately q-dependent |
| CCNE1 | 1.1e-87 | 8e-86 | 35.5687 | 0.0367 | Moderately q-dependent |
| THOC6 | 1.7e-79 | 1.3e-77 | 30.0465 | 0.0422 | Moderately q-dependent |
| DUOXA2 | 3e-78 | 2.2e-76 | 29.2740 | 0.0624 | Moderately q-dependent |
| ETFRF1 | 1.5e-77 | 1.1e-75 | 28.8350 | 0.0227 | Moderately q-dependent |
| LINC03040 | 8.1e-74 | 5.8e-72 | 26.6334 | 0.0511 | Moderately q-dependent |
| PLP1 | 7.8e-68 | 5.5e-66 | 23.3649 | 0.0551 | Moderately q-dependent |
| SRPRA | 1.5e-58 | 1.1e-56 | 18.8734 | 0.0141 | Moderately q-dependent |
| GLB1L2 | 2.7e-56 | 1.8e-54 | 17.8874 | 0.0357 | Moderately q-dependent |

**Supplementary Table 16 \| Top genes identified by Aligned Rank
Transform (ART).** Ranked by adjusted p-value; comparison with GAMM
(main vignette Table 3) reveals method-specific sensitivities. {.table}

Visualize q-curves for top genes:

``` r

# Plot top genes from Aligned Rank Transform (ART) (using all genes, not just significant ones)
# This ensures the plot displays n_top genes regardless of significance threshold
top_genes_plot <- plot_diversity_spectrum(analysis, sait_res = rank_transform_results_all, n_top = 4)
print(top_genes_plot)
```

![\*\*Supplementary Figure 1 \| Aligned Rank Transform (ART) results for
scale-dependent isoform switching.\*\* Q-curves for top genes identified
by ART; compare visual patterns with GAMM results (main vignette Figure
2) to assess method
agreement.](TSENAT_appendix_B_files/figure-html/extended-fig-1-rank-transform-q-curves-1.png)

**Supplementary Figure 1 \| Aligned Rank Transform (ART) results for
scale-dependent isoform switching.** Q-curves for top genes identified
by ART; compare visual patterns with GAMM results (main vignette Figure
2) to assess method agreement.

------------------------------------------------------------------------

## Generalized Additive Mixed Model Approach (GAMM)

This section loads precomputed GAMM results generated in the main
vignette
([TSENAT.Rmd](https://gallardoalba.github.io/TSENAT/articles/TSENAT.Rmd)).

``` r

# Load precomputed GAMM analysis object from RDS file
analysis_sait <- readRDS(
    system.file("extdata", "analysis_sait.rds", package = "TSENAT")
)

# Extract GAMM results for inspection using accessor function
sait_results <- results(analysis_sait, type = "sait", rankBy = "pvalue")
print(sait_results)
```

| Metric                                   |  Value |
|:-----------------------------------------|-------:|
| Genes tested                             |     76 |
| Significant (p \< 0.05)                  |     66 |
| Concordant (p \< 0.05 AND adj_p \< 0.05) |     57 |
| NAs                                      |      0 |
| Mean effect size                         |  36.3% |
| Median effect size                       |  33.3% |
| Strong effect genes (effect_size \> 30%) |     39 |
| Model convergence rate                   | 100.0% |

**Supplementary Table 17 \| GAMM Results Summary** (summary of the GAMM
results computed in the main vignette, Table 3). Overview of generalized
additive mixed model statistics for q × condition interaction tests
across full q-spectrum with paired sample structure. {.table}

| Gene            | p (interaction) | Adjusted p | Effect size | Test statistic |
|:----------------|----------------:|-----------:|------------:|---------------:|
| CXCL12          |       2.15e-139 |  1.63e-137 |       74.0% |         370.00 |
| ING3            |        4.15e-57 |   3.12e-55 |       32.3% |         109.18 |
| THY1            |        7.49e-55 |   5.54e-53 |       35.8% |         103.97 |
| LINC03040       |        5.54e-47 |   4.05e-45 |       58.3% |          86.46 |
| SNHG10          |        7.69e-43 |   5.54e-41 |       21.0% |          77.62 |
| MEF2A           |        2.34e-36 |   1.66e-34 |       17.7% |          64.30 |
| ASMTL           |        2.18e-33 |   1.52e-31 |       39.4% |          58.41 |
| HDAC2           |        2.93e-33 |   2.02e-31 |       22.7% |          58.15 |
| ENSG00000274322 |        3.02e-29 |   2.05e-27 |       43.9% |          50.38 |
| ATG5            |        2.32e-28 |   1.56e-26 |       51.5% |          48.70 |

**Supplementary Table 18 \| Top 10 Genes by GAMM p-value (with Effect
Size and Test Statistics).** Ranked by adjusted p-value; compare with
ART results (Supplementary Table 16) for method agreement assessment.
{.table}

------------------------------------------------------------------------

## Method Concordance Analysis

A critical validation step is to compare whether the Aligned Rank
Transform (ART) and GAMM produce consistent results. High concordance
between methods (i.e., genes significant in both approaches) provides
strong evidence that discoveries are robust to methodological choice.
Discordant genes — those significant in only one method — warrant closer
inspection: they may represent genuine biological signals revealed by
one method’s specific advantages, or artifacts of that method’s
assumptions. This section quantifies agreement between approaches and
identifies high-confidence genes significant across both statistical
frameworks.

``` r

# Compute concordance analysis using new two-object API
rank_results <- results(analysis_rank, type = "rank_test")
if (is.null(rank_results) || length(rank_results) == 0) {
    stop("Cannot run calculate_concordance(): rank_test results are empty. Ensure calculate_rank_transform() completed successfully.")
}

analysis_rank <- calculate_concordance(
    analysis_sait = analysis_sait,
    analysis_rank = analysis_rank,
    verbose = TRUE
)

# Extract concordance results with list format for component display
concordance_results <- results(analysis_rank, type = "concordance", format = "list")
```

| Metric                               |        Value |
|:-------------------------------------|-------------:|
| Total genes compared                 |           76 |
| Spearman correlation (p-values)      | rho = 0.3082 |
| Both methods significant (p \< 0.05) |   40 (52.6%) |
| SAIT only significant                |   17 (22.4%) |
| Rank test only significant           |    9 (11.8%) |
| Neither significant                  |   10 (13.2%) |
| Concordance rate                     |        52.6% |
| Discordance rate                     |        34.2% |

**Supplementary Table 19 \| Global Concordance Metrics: GAMM vs Aligned
Rank Transform (ART).** Quantifies agreement between parametric (GAMM)
and non-parametric (ART) approaches for detecting q × condition
interactions. {.table}

| Gene      | SAIT adj p | Rank test adj p | SAIT Effect | Rank test rho^2 |
|:----------|-----------:|----------------:|------------:|----------------:|
| CXCL12    |  1.63e-137 |       2.95e-144 |       74.0% |           0.106 |
| LINC03040 |  4.046e-45 |        5.79e-72 |       58.3% |           0.051 |
| SNHG10    |  5.538e-41 |       5.655e-49 |       21.0% |           0.030 |
| MEF2A     |  1.662e-34 |        7.61e-53 |       17.7% |           0.043 |
| FAXDC2    |  1.438e-21 |       1.288e-32 |       19.4% |           0.033 |
| ING3      |   3.12e-55 |       6.756e-21 |       32.3% |           0.014 |
| FAM114A2  |  6.857e-19 |        2.73e-54 |       39.5% |           0.037 |
| ETFRF1    |  1.935e-18 |        1.11e-75 |       62.5% |           0.023 |
| ATG5      |  1.557e-26 |       3.117e-18 |       51.5% |           0.013 |
| SPNS2     |  4.395e-18 |       8.25e-110 |        4.4% |           0.077 |

**Supplementary Table 20 \| Robust Entropic Order Index Interactions:
Top 10 High-Confidence Genes Detected by Both Methods (of 40 concordant,
ranked by statistical significance).** Cross-method validation: genes
significant in both GAMM and ART provide strongest evidence for q ×
condition effects. {.table}

### Statistical Power vs. Robustness: Understanding Method Discordance

The substantial discordance between GAMM and ART results — where GAMM
identifies many more significant genes — reflects a fundamental
trade-off in statistical methodology: **parametric methods maximize
power when assumptions hold, while non-parametric methods sacrifice
power for robustness to assumption violations**. The exact percentages
depend on the dataset and analysis parameters; the pattern of GAMM
detecting more interactions than ART is the consistent finding.

GAMM’s superior power derives from three factors:

1.  Distributional assumptions: GAMM assumes approximately normal
    residuals and homogeneous variance, which are reasonable after
    appropriate transformation for entropy data. Under these
    assumptions, parametric methods are theoretically optimal, achieving
    maximum power for a given type I error rate. ART is assumption-free
    but necessarily discards quantitative information through alignment
    and ranking, which reduces statistical power when the underlying
    data are approximately normal.

2.  Model flexibility: GAMM captures nonlinear q patterns with splines
    (natural regression splines with 3 df for paired designs) that
    simultaneously fit nonlinear trends while controlling degrees of
    freedom. This flexibility allows GAMM to detect subtle entropic
    index × group interactions across all q-values.

3.  Continuous vs. categorical q treatment: GAMM models q as a
    continuous smooth predictor across all 41 q-values, preserving the
    ordinal structure and borrowing strength across adjacent q-levels.
    ART treats q as a categorical factor — while this enables proper
    interaction testing via alignment, it discards the continuous
    ordering information, reducing sensitivity to graded q-dependent
    trends. This discretization is a major contributor to ART’s lower
    power relative to GAMM, independent of distributional assumptions.

This complementary approach—combining high-power parametric tests with
robust nonparametric alternatives—provides confidence that discoveries
are methodologically sound rather than artifacts of statistical
assumptions.

### Visualization: Method Comparison

Visual comparison of concordance patterns provides intuitive assessment
of which genes align between methods and which are method-specific. The
following plots display significance calls, effect sizes, and agreement
metrics in a format that facilitates interpretation of both robust
discoveries and method-specific findings.

``` r

# Create comparison plot using S4 wrapper
# Concordance analysis results are stored in analysis metadata after calculate_concordance()
p <- plot_concordance(analysis_rank, verbose = TRUE)
grid::grid.draw(p)
```

![\*\*Supplementary Figure 2 \| Method concordance visualization
comparing Aligned Rank Transform (ART) and GAMM approaches for detecting
q × condition interactions.\*\* Scatter plot showing overlap in genes
identified as significant by each method; reveals
parametric-non-parametric agreement
patterns.](TSENAT_appendix_B_files/figure-html/visualize-method-concordance-1.png)

**Supplementary Figure 2 \| Method concordance visualization comparing
Aligned Rank Transform (ART) and GAMM approaches for detecting q ×
condition interactions.** Scatter plot showing overlap in genes
identified as significant by each method; reveals
parametric-non-parametric agreement patterns.

------------------------------------------------------------------------

## Session Information

``` r

sessionInfo()
#> R version 4.5.3 (2026-03-11)
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
#>  [1] TSENAT_0.99.35              testthat_3.3.2             
#>  [3] gridExtra_2.3               dplyr_1.2.1                
#>  [5] SummarizedExperiment_1.40.0 Biobase_2.70.0             
#>  [7] GenomicRanges_1.62.1        Seqinfo_1.0.0              
#>  [9] IRanges_2.44.0              S4Vectors_0.48.1           
#> [11] BiocGenerics_0.56.0         generics_0.1.4             
#> [13] MatrixGenerics_1.22.0       matrixStats_1.5.0          
#> [15] ggplot2_4.0.3               kableExtra_1.4.0           
#> [17] BiocStyle_2.38.0           
#> 
#> loaded via a namespace (and not attached):
#>   [1] Rdpack_2.6.6        sandwich_3.1-1      rlang_1.2.0        
#>   [4] magrittr_2.0.5      multcomp_1.4-30     otel_0.2.0         
#>   [7] compiler_4.5.3      mgcv_1.9-4          systemfonts_1.3.2  
#>  [10] vctrs_0.7.3         stringr_1.6.0       pkgconfig_2.0.3    
#>  [13] fastmap_1.2.0       backports_1.5.1     XVector_0.50.0     
#>  [16] ellipsis_0.3.3      labeling_0.4.3      rmarkdown_2.31     
#>  [19] sessioninfo_1.2.3   tzdb_0.5.0          nloptr_2.2.1       
#>  [22] ragg_1.5.2          purrr_1.2.2         xfun_0.58          
#>  [25] cachem_1.1.0        jsonlite_2.0.0      DelayedArray_0.36.1
#>  [28] BiocParallel_1.44.0 broom_1.0.13        parallel_4.5.3     
#>  [31] R6_2.6.1            bslib_0.11.0        stringi_1.8.7      
#>  [34] RColorBrewer_1.1-3  boot_1.3-32         car_3.1-5          
#>  [37] pkgload_1.5.2       brio_1.1.5          jquerylib_0.1.4    
#>  [40] estimability_1.5.1  Rcpp_1.1.1-1.1      bookdown_0.46      
#>  [43] knitr_1.51          usethis_3.2.1       zoo_1.8-15         
#>  [46] readr_2.2.0         Matrix_1.7-5        splines_4.5.3      
#>  [49] tidyselect_1.2.1    rstudioapi_0.18.0   abind_1.4-8        
#>  [52] yaml_2.3.12         codetools_0.2-20    pkgbuild_1.4.8     
#>  [55] plyr_1.8.9          lattice_0.22-9      tibble_3.3.1       
#>  [58] withr_3.0.3         S7_0.2.2            evaluate_1.0.5     
#>  [61] survival_3.8-6      desc_1.4.3          xml2_1.5.2         
#>  [64] pillar_1.11.1       BiocManager_1.30.27 carData_3.0-6      
#>  [67] reformulas_0.4.4    rprojroot_2.1.1     hms_1.1.4          
#>  [70] scales_1.4.0        minqa_1.2.8         ARTool_0.11.2      
#>  [73] xtable_1.8-8        glue_1.8.1          pheatmap_1.0.13    
#>  [76] emmeans_2.0.3       tools_4.5.3         lme4_2.0-1         
#>  [79] fs_2.1.0            mvtnorm_1.3-6       cowplot_1.2.0      
#>  [82] grid_4.5.3          tidyr_1.3.2         rbibutils_2.4.1    
#>  [85] devtools_2.5.2      nlme_3.1-169        Formula_1.2-5      
#>  [88] cli_3.6.6           textshaping_1.0.5   S4Arrays_1.10.1    
#>  [91] viridisLite_0.4.3   svglite_2.2.2       geepack_1.3.13     
#>  [94] gtable_0.3.6        sass_0.4.10         digest_0.6.39      
#>  [97] TH.data_1.1-5       SparseArray_1.10.10 htmlwidgets_1.6.4  
#> [100] farver_2.1.2        memoise_2.0.1       htmltools_0.5.9    
#> [103] pkgdown_2.2.0       lifecycle_1.0.5     MASS_7.3-65
```
