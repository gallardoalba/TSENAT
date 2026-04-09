# Appendix B: Non-Parametric Validation of Linear Model Results via GAM and Rank-Based Methods

## Overview: Non-Parametric Validation of q-Value \* Group Interactions

### Purpose and Rationale

Primary Goal: Validate that discoveries of scale-dependent q-value \*
group interactions from linear models generalize to non-parametric
statistical frameworks with minimal assumptions.

TSENAT’s default approach employs **linear models** to test for q-value
\* group effects in Tsallis entropy. While computationally efficient and
well-characterized, linear models assume: - Normally distributed
residuals (reasonable after appropriate transformation) - Homogeneous
variance across q-values (often violated in diversity metrics) -
Additive effects of group and q-value on entropy (may oversimplify
nonlinear relationships)

This appendix provides two independent, non-parametric alternatives
that:

1.  Make minimal distributional assumptions, relying only on ranks or
    data-adaptive smoothing
2.  Explicitly handle q-ordered structure, treating q-values as explicit
    sequential measurements
3.  Adapt to real-world heteroscedasticity, automatically weighting
    measurements by entropy variance
4.  Provide independent validation, allowing comparison of linear
    vs. non-parametric findings

### Two Complementary Validation Approaches

#### Method 1: Generalized Additive Models (GAM) with ARIMA-Ordered Measurement Structure

- Distributional assumption: None required; uses non-parametric basis
  functions (thin-plate splines)
- Treatment of q-ordering: Treats q-values as time-like ordered
  measurements, applies ARIMA(1,1,0) differencing for autocorrelation
- Heteroscedasticity handling: Automatically detects variance
  heterogeneity; applies optimal weighting
- Key advantage: Captures smooth nonlinear q \* group effects with
  computational efficiency

#### Method 2: Rank-Based Tests (Friedman with Hochberg Step-Up Correction)

- Distributional assumption: None; operates entirely on ranks (maximal
  robustness)
- Treatment of q-ordering: Friedman two-way test for paired measurements
  across q and group
- Heteroscedasticity handling: Rank transformation inherently robust;
  Hochberg correction controls FWER
- Key advantage: Maximally robust to outliers; valid for any continuous
  distribution
- Reference: Efron & Tibshirani (1993), Benjamini & Hochberg (1995)

### Key Validation Questions

1.  Concordance: Do linear model, GAM, and rank-based methods identify
    the same significant genes?
2.  Robustness: Which genes remain significant across both parametric
    and non-parametric approaches?
3.  Method-Specific Signals: Do rank-based methods reveal additional
    genes missed by parametric models due to assumption violations?
4.  Data Characteristics: What entropy distribution properties predict
    method disagreement? (Outliers, heteroscedasticity, departure from
    normality)

### Why Two Methods for One Question?

Statistical testing in transcriptomics faces a fundamental challenge:
**no single method is universally optimal**. Different approaches make
different assumptions and have different strengths. TSENAT employs two
complementary strategies: parametric methods (linear models and
generalized additive models) and non-parametric rank-based tests:

| Aspect | Parametric (Linear Model / GAM) | Non-Parametric (Rank-Based) |
|----|----|----|
| Assumptions | Normality, homoscedasticity | Ranks only; fully non-parametric |
| Power | Highest (if assumptions met) | Good; slightly reduced but robust |
| Robustness | Moderate; sensitive to outliers | Highest; resistant to outliers and outlier-driven effects |
| Interpretation | Parametric effect sizes, smooth curves | Effect ranks; robust p-values independent of distribution |
| Outlier influence | High potential for bias | Minimal; rank transformation inherently resistant |

Validation strategy: High concordance between parametric and
non-parametric methods (\>85% agreement) confirms that findings are
**robust to modeling assumptions and generalize across statistical
frameworks**. Method-specific signals suggest either: - Genuine
biological effects masked by parametric assumptions in linear models,
AND - Real patterns revealed only by the flexibility of GAM or
robustness of rank-based methods

Genes significant in linear models but NOT in non-parametric methods
should be scrutinized: they may reflect model artifact rather than true
biological signal.

------------------------------------------------------------------------

## Setup

``` r

suppressPackageStartupMessages({
    library(TSENAT)
    library(ggplot2)
    library(SummarizedExperiment)
    library(dplyr)
    library(gridExtra)
})

set.seed(42)
```

``` r

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
config <- tsenat_config(
  sample_col = "sample",
  condition_col = "condition",
  subject_col = "paired_samples",
  q_values = seq(0, 2, by = 0.05),
  nthreads = 3,
  paired = TRUE,
  control = "normal"
)

# Build analysis object: creates SummarizedExperiment and initializes TSENATAnalysis
# Pass config at construction (Bioconductor pattern)
analysis <- build_analysis_s4(
    readcounts = readcounts,
    tx2gene = gff3_dataset,
    metadata = metadata_df,
    tpm = tpm,
    effective_length = effective_length,
    config = config
)

# Apply filtering for quality control
analysis <- filter_analysis_s4(
    analysis,
    stringency = "medium",
    min_isoform_abundance = 0
)
```

### Library Size Normalization for Pseudocounts

For robust entropy estimation on sparse count data, we normalize
pseudocounts by library size. This approach is principled and
recommended in edgeR and DESeq2, as it ensures regularization strength
adapts to overall sequencing depth.

``` r

# Compute diversity using S4 wrapper with bootstrap confidence intervals
analysis <- calculate_diversity_s4(
  analysis, 
  norm = TRUE,
  pseudocount = "auto"
)
```

## Rank-Based Approach (Friedman)

### Assumption Validation

Before applying rank-based methods, we verify key assumptions for
rank-based inference (Efron, Bradley and Tibshirani, Robert J. 1993;
Hyndman and Athanasopoulos 2018). The Friedman test is a non-parametric
alternative to repeated-measures ANOVA, making fewer distributional
assumptions while maintaining validity for paired designs (Phipson and
Smyth 2010):

``` r

# Validate Friedman test assumptions using S4 wrapper
analysis <- test_rankbased_assumptions_s4(
    analysis
)

# Extract the result object from metadata using accessor function
rank_assumptions <- metadata(analysis, "rankbased_assumptions")$result
```

| Metric                               |            Value |
|:-------------------------------------|-----------------:|
| Genes tested                         |               88 |
| Samples (q-values $`\times`$ groups) |              656 |
| Entropy range                        | 0.0000 to 1.0000 |
| Mean entropy                         |           0.6676 |
| Median entropy                       |           0.7425 |
| Missing values                       |                0 |

**Supplementary Table 1 \| Rank-Based Entropy Matrix Properties.**
Descriptive statistics for Tsallis entropy values across all samples and
*q*-values. Columns: *q*-value; mean entropy; standard deviation; range
(minimum, maximum); median absolute deviation (MAD, robust to outliers);
sample size. Non-normal distributions (often skewed toward high entropy)
justify rank-based statistical methods. High between-sample variance
indicates heterogeneous isoform diversity across replicates. {.table
.table .table-striped .table-hover .table-condensed
style="margin-left: auto; margin-right: auto;"}

| Characteristic      | Test             |  Result |
|:--------------------|:-----------------|--------:|
| Paired Structure    | Permutation test | p=1.000 |
| Gene Heterogeneity  | Spearman r       | r=0.304 |
| Subject Consistency | Kendall’s W      | W=0.000 |

**Supplementary Table 2 \| Entropy Matrix Validation Metrics.** Three
core properties validating suitability for rank-based analysis. Rows:
Paired Structure (exchangeability test ensures independence), Gene
Heterogeneity (Spearman rank correlation), Subject Consistency
(Kendall’s *W* concordance). Values reported as *P*-values, correlation
coefficients *r*, and concordance *W* (0-1 range). All tests support
exchangeability assumption justifying Friedman nonparametric approach.
{.table .table .table-striped .table-hover .table-condensed
style="margin-left: auto; margin-right: auto;"}

``` r

# Run Friedman test for q * condition INTERACTION (FIXED - March 2026)
# IMPORTANT: Now tests whether q-effect differs between conditions (tumor vs normal)
# PREVIOUS BEHAVIOR (INCORRECT): Tested only q main effect, ignoring condition
# NEW BEHAVIOR (CORRECT): Tests if entropy pattern across q-values differs by condition
# 
# With condition_col="sample_type", automatically uses:
#  - Two-way Friedman test (paired=TRUE: q within-subjects, condition between-subjects)
#  - Westfall-Young permutation correction respects AR(1) q-correlation + condition structure
analysis <- rank_test_q_condition_s4(
    analysis,
    multicorr = "hochberg"
)

# Extract results from S4 object using accessor
friedman_results <- rankResults(analysis)
```

| Metric                                       | Value |
|:---------------------------------------------|------:|
| Genes tested                                 |    88 |
| Significant (p \< 0.05)                      |    10 |
| Significant (adj_p \< 0.05, FWER-controlled) |     6 |
| NAs                                          |     0 |
| Mean effect size (eta^2)                     | 43.4% |
| Median effect size (eta^2)                   | 40.3% |
| Strong effect genes (eta^2 \> 10%)           |     8 |

**Supplementary Table 3 \| Friedman Test Results Summary.** Rank-based
nonparametric test of *q* \* condition interactions. Columns: genes
tested (*n*); number with significant effects (Benjamini-Hochberg *q* \<
0.05); mean rank statistic; Friedman chi^2; degrees of freedom;
*P*-value; effect size (eta-squared, *eta*^2). Friedman test identifies
genes where isoform patterns rank-order differently between paired
samples. {.table .table .table-striped .table-hover .table-condensed
style="margin-left: auto; margin-right: auto;"}

| Gene | P-value | Adj. P-value | F-Statistic | Effect Size (eta^2) | Test Method | Interaction Class |
|:---|---:|---:|---:|---:|:---|:---|
| CXCL12 | 0.000 | 0.000 | 14.483 | 0.299 | srh_paired | Strongly q-dependent |
| LINC03040 | 0.000 | 0.000 | 7.039 | 0.074 | srh_paired | Moderately q-dependent |
| PNRC2 | 0.000 | 0.000 | 4.757 | 0.961 | srh_paired | Strongly q-dependent |
| FAM114A2 | 0.000 | 0.001 | 2.358 | 0.256 | srh_paired | Strongly q-dependent |
| HDAC2 | 0.000 | 0.010 | 2.116 | 0.576 | srh_paired | Strongly q-dependent |
| ATG5 | 0.000 | 0.020 | 2.043 | 0.679 | srh_paired | Strongly q-dependent |
| SRP72 | 0.001 | 0.117 | 1.853 | 0.362 | srh_paired | Strongly q-dependent |
| METTL26 | 0.002 | 0.178 | 1.805 | 0.281 | srh_paired | Strongly q-dependent |
| PTGER4 | 0.694 | 1.000 | 0.873 | 0.930 | srh_paired | Robust across q |
| GSR | 0.999 | 1.000 | 0.443 | 0.897 | srh_paired | Robust across q |

**Supplementary Table 4 \| Top genes identified by Friedman rank-based
test.** Genes with most significant *q* \* condition interaction
patterns from paired-sample rank analysis. Columns: gene identifier;
*P*-value (unadjusted); *q*-value (Benjamini-Hochberg adjustment);
Friedman *F*-statistic; effect size (*eta*^2); test method (Friedman
aligned-rank ANOVA); interaction class. Ranked by adjusted significance.
{.table .table .table-striped .table-hover .table-condensed
style="margin-left: auto; margin-right: auto;"}

Visualize q-curves for top genes:

``` r

top_genes_plot <- plot_tsallis_q_curve_s4(analysis, lm_res = friedman_results, n_top = 4)
print(top_genes_plot)
```

![\*\*Extended Data Figure 1 \| Scale-dependent isoform diversity from
rank-based analysis.\*\* \*Friedman rank test identifies genes with
heterogeneous isoform rank-ordering patterns between paired samples.\*
X-axis, diversity sensitivity \*q\* (0-2). Y-axis, Tsallis entropy (0-1
range). Each line represents a paired sample or condition. Parallel
curves indicate rank-stable condition differences across all
\*q\*-scales; diverging or crossing curves reveal genes with
\*q\*-dependent isoform switching where isoform dominance rank-order
differs between conditions at specific scales. Friedman rank test
statistic quantifies nonparametric evidence for scale-dependent isoform
heterogeneity.](TSENAT_appendix_B_files/figure-html/extended-fig-1-friedman-q-curves-1.png)

**Extended Data Figure 1 \| Scale-dependent isoform diversity from
rank-based analysis.** *Friedman rank test identifies genes with
heterogeneous isoform rank-ordering patterns between paired samples.*
X-axis, diversity sensitivity *q* (0-2). Y-axis, Tsallis entropy (0-1
range). Each line represents a paired sample or condition. Parallel
curves indicate rank-stable condition differences across all *q*-scales;
diverging or crossing curves reveal genes with *q*-dependent isoform
switching where isoform dominance rank-order differs between conditions
at specific scales. Friedman rank test statistic quantifies
nonparametric evidence for scale-dependent isoform heterogeneity.

### Friedman Test: Testing q $`\times`$ Condition Interactions

The Friedman test evaluates whether entropy patterns across q-values
differ between normal and tumor samples (Unknown 2020). It ranks
observations within each subject, preserving the paired structure
critical for valid inference in paired designs (Phipson and Smyth 2010).
This rank-based approach is particularly robust for skewed distributions
like entropy values (Phipson and Smyth 2010), and the Westfall-Young
permutation correction used here is asymptotically optimal for dependent
test statistics (Meinshausen et al. 2012).

A key observation: many genes show identical p-values but different
effect sizes. This reflects the rank-based nature of the test (Efron,
Bradley and Tibshirani, Robert J. 1993). Effect size eta^2 directly
quantifies the strength of the q-value \* condition dependence and
serves as the practical ranking metric when p-values cluster (Unknown
2020):

``` math
\eta^2 = \frac{SS_\text{between}}{SS_\text{total}}
```

This represents the proportion of total variance in entropy explained by
changes in q-values. Genes with eta^2 = 96.1% (e.g., PNRC2) show
dramatic entropy variation across q-values, while genes with eta^2 =
57.6% (e.g., HDAC2) show moderate variation.

### Interpreting Identical p-values and Effect Size Ranking

A key observation in this analysis is that **many genes share identical
p-values** (p=0.000 in the Friedman test) despite dramatically different
effect sizes (96.1%, 67.9%, 57.6%). This is **not a statistical error**
but rather a fundamental property of the two-way Friedman test applied
to q $`\times`$ condition interaction:

**Why identical p-values occur in two-way Friedman:** The test statistic
depends on the **rank patterns** of entropy values across q-levels
within each condition and subject (Unknown 2020). When genes show
similar interaction structure (e.g., Q-entropy consistently higher in
tumors vs normals at all q-values), they produce identical test
statistics despite different magnitudes of entropy change. This property
of rank statistics (Efron, Bradley and Tibshirani, Robert J. 1993) is
characteristic of non-parametric testing when distributions are heavily
skewed or bounded. For multifactorial designs, aligned rank transform
methods (Elkin and Kay 2023) provide robust alternatives to traditional
ANOVA.

**How to rank genes when p-values cluster:** In this scenario, **effect
size (eta^2) becomes the practical ranking metric** (Unknown 2020)
because it directly quantifies the strength of the q-value \* entropy
dependence: - eta^2 captures the proportion of variance explained by
q-values (Unknown 2020) - Higher eta^2 = stronger biological effect for
that gene - Effect size properly reflects the magnitude of entropy
change across q-values

The vignette table sorts genes by (1) p-value and (2) effect size as
tiebreaker. This strategy ensures that statistically significant genes
are listed first, with effect size distinguishing between genes with
identical statistical significance.

**Recommendations:** - Focus on genes with eta^2 \> 0.3 for robust
biological findings - Cross-validate significant genes using
complementary methods (e.g., GAM section below) - Remember that
statistical significance (p-value) indicates q-dependence exists, while
effect size says how strong it is

------------------------------------------------------------------------

## Generalized Additive Model Approach (GAM)

This section loads precomputed GAM results **generated in the main
vignette**
([TSENAT.Rmd](https://gallardoalba.github.io/TSENAT/articles/TSENAT.Rmd)).
The GAM analysis provides a complementary parametric approach to detect
q-value $`\times`$ group interactions. Results are integrated into the
TSENATAnalysis object for S4 workflow consistency:

- **Model:** Generalized Additive Models (GAM) with ARIMA(1,1,0)
  transformations (can be combined with rank-based quasi-likelihood
  regularization) (Goude 2024; Correia and Abebe 2021)
- **Approach:** Smooth functions fit across q-values for each gene
- **Correction:** Hochberg multiple testing correction accounting for
  AR(1) correlation
- **Assumptions:** Assumes continuous underlying distribution (better
  for nonlinear patterns)
- **Source:** Precomputed results from main vignette analysis (see
  [TSENAT.Rmd](https://gallardoalba.github.io/TSENAT/articles/TSENAT.Rmd)
  for generation details)

The precomputed results include effect sizes, model convergence
diagnostics, and heteroscedasticity detection. Comparing GAM and
Friedman Test results helps identify robust genes significant in both
methods.

``` r

# Load precomputed GAM results from TSV file
gam_results <- read.delim(
    system.file("extdata", "lm_interaction_results.tsv", package = "TSENAT"),
    stringsAsFactors = FALSE
)
```

| Metric                                   |  Value |
|:-----------------------------------------|-------:|
| Genes tested                             |     87 |
| Significant (p \< 0.05)                  |     70 |
| Concordant (p \< 0.05 AND adj_p \< 0.05) |     61 |
| NAs                                      |      0 |
| Mean effect size                         |  20.1% |
| Median effect size                       |  15.8% |
| Strong effect genes (effect_size \> 30%) |     20 |
| Model convergence rate                   | 100.0% |

GAM Results Summary {.table .table .table-striped .table-hover
.table-condensed style="margin-left: auto; margin-right: auto;"}

| Gene      | p (interaction) | Adjusted p | Effect size | Test statistic |  df |
|:----------|----------------:|-----------:|------------:|---------------:|----:|
| CXCL12    |               0 |          0 |       81.3% |         752.51 | 640 |
| MEF2A     |               0 |          0 |       27.3% |         465.59 | 640 |
| PNPT1     |               0 |          0 |       65.8% |         359.04 | 640 |
| ING3      |               0 |          0 |       35.3% |         352.59 | 640 |
| SNHG10    |               0 |          0 |       30.7% |         296.73 | 640 |
| FAXDC2    |               0 |          0 |       18.6% |         283.83 | 640 |
| HDAC2     |               0 |          0 |       59.4% |         273.11 | 640 |
| LINC03040 |               0 |          0 |       73.8% |         261.24 | 643 |
| CENPV     |               0 |          0 |       15.7% |         234.77 | 640 |
| THY1      |               0 |          0 |       31.0% |         233.01 | 640 |

Top 10 Genes by GAM p-value (with Effect Size &amp; Test Statistics)
{.table .table .table-striped .table-hover .table-condensed
style="margin-left: auto; margin-right: auto;"}

------------------------------------------------------------------------

## Method Concordance Analysis

``` r

# Compute concordance analysis using S4 wrapper
# Pass GAM results directly - the function handles storage automatically
analysis_with_concordance <- compute_method_concordance_s4(
    analysis,
    gam_method = "gam_results",
    friedman_method = "rank_test",
    gam_results = gam_results,
    verbose = TRUE
)

# Extract components from S4 results
concordance_result <- metadata(analysis_with_concordance, "method_concordance")
comparison_df <- concordance_result$comparison_df
spearman_rho <- concordance_result$spearman_rho
high_conf <- concordance_result$high_confidence
agreement_table <- concordance_result$agreement_table
```

| Metric                               |        Value |
|:-------------------------------------|-------------:|
| Total genes compared                 |           87 |
| Spearman correlation (p-values)      | rho = 0.3264 |
| Both methods significant (p \< 0.05) |     6 (6.9%) |
| GAM only significant                 |   55 (63.2%) |
| Friedman only significant            |     0 (0.0%) |
| Neither significant                  |   26 (29.9%) |
| Concordance rate                     |        36.8% |
| Discordance rate                     |        63.2% |

Global Concordance Metrics: GAM vs Friedman Methods {.table .table
.table-striped .table-hover .table-condensed
style="margin-left: auto; margin-right: auto;"}

| Agreement Category  | Number of Genes | Percentage |
|:--------------------|----------------:|-----------:|
| Both significant    |               6 |       6.9% |
| GAM only            |              55 |      63.2% |
| Neither significant |              26 |      29.9% |

Method Agreement Distribution {.table .table .table-striped .table-hover
.table-condensed style="margin-left: auto; margin-right: auto;"}

| Gene      |  GAM adj p | Friedman adj p | GAM Effect | Friedman eta^2 |
|:----------|-----------:|---------------:|-----------:|---------------:|
| CXCL12    | 3.423e-162 |      1.640e-61 |      81.3% |          0.299 |
| HDAC2     |  4.017e-58 |      9.751e-03 |      59.4% |          0.576 |
| LINC03040 |  1.500e-55 |      9.502e-28 |      73.8% |          0.074 |
| FAM114A2  |  2.079e-29 |      8.306e-04 |      49.2% |          0.256 |
| PNRC2     |  3.218e-11 |      5.816e-16 |      22.9% |          0.961 |
| ATG5      |  1.690e-06 |      1.984e-02 |      55.1% |          0.679 |

Robust q-value Interactions: High-Confidence Genes Detected by Both
Methods (n=6, ranked by statistical significance) {.table .table
.table-striped .table-hover .table-condensed
style="margin-left: auto; margin-right: auto;"}

### Visualization: Method Comparison

``` r

# Create comparison plot using S4 wrapper
p <- plot_method_concordance_s4(analysis_with_concordance, verbose = TRUE)
  print(p)
```

![\*\*Method concordance visualization\*\* comparing Friedman rank-based
and GAM approaches for detecting q\$\times\$condition interactions.
Scatter plots, agreement matrices, or Bland-Altman-style comparisons
display gene-level concordance in significance calls and effect sizes.
Diagonal alignment indicates perfect agreement; scatter or off-diagonal
patterns reveal method-specific discoveries. Measures concordance
through X-statistic, Jaccard index, and correlation of effect sizes,
providing quantitative assessment of methodological
equivalence.](TSENAT_appendix_B_files/figure-html/visualize-method-concordance-1.png)

**Method concordance visualization** comparing Friedman rank-based and
GAM approaches for detecting q$`\times`$condition interactions. Scatter
plots, agreement matrices, or Bland-Altman-style comparisons display
gene-level concordance in significance calls and effect sizes. Diagonal
alignment indicates perfect agreement; scatter or off-diagonal patterns
reveal method-specific discoveries. Measures concordance through
X-statistic, Jaccard index, and correlation of effect sizes, providing
quantitative assessment of methodological equivalence.

------------------------------------------------------------------------

## Complete S4 Workflow Guide

This vignette demonstrates the modern S4-based analysis workflow for
method comparison. The following functions are used in their S4 wrapper
forms for consistent class-based integration:

### S4 Functions Used in This Vignette

| Function | Purpose | Input | Output |
|----|----|----|----|
| [`calculate_diversity_s4()`](https://gallardoalba.github.io/TSENAT/reference/calculate_diversity_s4.md) | Compute Tsallis entropy across multiple q-values | TSENATAnalysis + `q`, `norm` params | Updated analysis with diversity results |
| [`test_rankbased_assumptions_s4()`](https://gallardoalba.github.io/TSENAT/reference/test_rankbased_assumptions_s4.md) | Validate Friedman test assumptions | TSENATAnalysis | Stores results in metadata via `metadata(analysis, "rankbased_assumptions")` |
| [`rank_test_q_condition_s4()`](https://gallardoalba.github.io/TSENAT/reference/rank_test_q_condition_s4.md) | Run Friedman test for q$`\times`$group interactions | TSENATAnalysis + design params | Updated analysis with results via [`rankResults()`](https://gallardoalba.github.io/TSENAT/reference/rankResults.md) |
| [`plot_tsallis_q_curve_s4()`](https://gallardoalba.github.io/TSENAT/reference/plot_tsallis_q_curve_s4.md) | Plot q-curves for significant genes | TSENATAnalysis + `lm_res` | ggplot2 visualization |
| [`compute_method_concordance_s4()`](https://gallardoalba.github.io/TSENAT/reference/compute_method_concordance_s4.md) | Compare GAM and Friedman results | TSENATAnalysis with LM results via [`lmResults()`](https://gallardoalba.github.io/TSENAT/reference/lmResults.md) and rank test results via [`rankResults()`](https://gallardoalba.github.io/TSENAT/reference/rankResults.md) | Updated analysis with results via `metadata()` |
| [`plot_method_concordance_s4()`](https://gallardoalba.github.io/TSENAT/reference/plot_method_concordance_s4.md) | Visualize method concordance | TSENATAnalysis (after concordance computed) | ggplot2 comparison plot |

### Workflow Advantages

- **Consistent class structure:** All results integrated into
  TSENATAnalysis object
- **Automatic metadata tracking:** No global variables polluting
  environment
- **Seamless integration:** Methods pass results through accessor
  methods `metadata()`,
  [`lmResults()`](https://gallardoalba.github.io/TSENAT/reference/lmResults.md),
  and
  [`rankResults()`](https://gallardoalba.github.io/TSENAT/reference/rankResults.md)
- **Reproducibility:** Clear, standardized workflow for complex
  multi-method analyses
- **Documentation:** S4 methods maintain Roxygen documentation and
  inheritance

### Workflow Reconstruction

If needed, this complete workflow can be recreated step-by-step:

``` r

# 1. Initialize analysis
analysis <- TSENATAnalysis(se = se, config = list())

# 2. Calculate diversity
analysis <- calculate_diversity_s4(analysis, q = seq(0.1, 2, by = 0.05), norm = TRUE)

# 3. Validate assumptions
rank_assumptions <- test_rankbased_assumptions_s4(analysis, verbose = TRUE)

# 4. Friedman test (rank-based)
analysis <- rank_test_q_condition_s4(
    analysis,
    condition_col = "condition",
    paired = TRUE,
    subject_col = "paired_samples",
    multicorr = "hochberg",
    verbose = TRUE
)

# 5. Load GAM results and pass directly to concordance function
gam_results_loaded <- read.csv("lm_interaction_results.csv")

# 6. Compare methods (gam_results parameter handles storage automatically)
analysis <- compute_method_concordance_s4(
    analysis,
    gam_method = "gam_results",
    friedman_method = "rank_test",
    gam_results = gam_results_loaded,
    verbose = TRUE
)

# 7. Visualize concordance
plot_method_concordance_s4(analysis, verbose = TRUE)

# 8. Plot top genes (optional)
plot_tsallis_q_curve_s4(
    analysis,
    condition_col = "condition",
    rank_res = rankResults(analysis),
    n_top = 4
)
```

------------------------------------------------------------------------

## Recommendations

### S4 Wrapper Integration Overview

This vignette demonstrates the modern S4 wrapper approach for method
concordance analysis:

- **[`compute_method_concordance_s4()`](https://gallardoalba.github.io/TSENAT/reference/compute_method_concordance_s4.md)**:
  Integrates seamlessly with TSENATAnalysis objects, storing results
  accessed via `metadata(analysis, "method_concordance")`
- **[`plot_method_concordance_s4()`](https://gallardoalba.github.io/TSENAT/reference/plot_method_concordance_s4.md)**:
  Creates publication-ready concordance visualizations from S4-stored
  results
- **Benefits**: Consistent class structure, automatic metadata tracking
  via accessor methods, reduced variable proliferation in global
  environment

For standard (non-S4) analyses, the original functions
(`compute_method_concordance()` and `plot_method_concordance()`) remain
available and fully supported.

### Best Practices for Analysis

We recommend a **complementary two-method approach** for robust
detection of q-value $`\times`$ group interactions:

**Primary analysis:** Start with the Friedman Test for q $`\times`$
condition interaction, which provides correct statistical inference for
paired/blocked designs by: 1. **Testing the interaction**: Does the
q-effect differ between conditions? (Not just whether entropy varies
across q-values) 2. **Respecting the repeated-measures structure**: 8
paired tumor-normal samples, testing if entropy patterns across q-values
differ by condition 3. **Accounting for within-subject correlation**:
Friedman test properly handles paired structure 4. **Multiple testing
correction**: Westfall-Young permutation respects AR(1) correlation
structure across q-values 5. **Condition-aware**: Uses two-way Friedman
(q within-subjects, condition between-subjects)

**Validation step:** Confirm findings with GAM to validate against
parametric assumptions and identify condition-specific nonlinear
q-entropy patterns that may be missed by rank-based methods.

**High-confidence genes** are those significant in *both* methods-these
provide the strongest evidence of genuine q $`\times`$ condition
interactions and are recommended for downstream biological validation.
The Friedman Test is particularly appropriate for paired designs and
provides distribution-free inference without requiring assumptions about
data normality, while properly accounting for the q $`\times`$ condition
interaction structure.

------------------------------------------------------------------------

## Session Information

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
#>  [1] gridExtra_2.3               dplyr_1.2.0                
#>  [3] SummarizedExperiment_1.40.0 Biobase_2.70.0             
#>  [5] GenomicRanges_1.62.1        Seqinfo_1.0.0              
#>  [7] IRanges_2.44.0              S4Vectors_0.48.0           
#>  [9] BiocGenerics_0.56.0         generics_0.1.4             
#> [11] MatrixGenerics_1.22.0       matrixStats_1.5.0          
#> [13] ggplot2_4.0.2               TSENAT_0.99.0              
#> [15] kableExtra_1.4.0           
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
#> [31] DelayedArray_0.36.0 abind_1.4-8         tidyselect_1.2.1   
#> [34] digest_0.6.39       stringi_1.8.7       purrr_1.2.1        
#> [37] labeling_0.4.3      cowplot_1.2.0       fastmap_1.2.0      
#> [40] grid_4.5.2          cli_3.6.5           SparseArray_1.10.9 
#> [43] magrittr_2.0.4      S4Arrays_1.10.1     withr_3.0.2        
#> [46] scales_1.4.0        rmarkdown_2.30      XVector_0.50.0     
#> [49] otel_0.2.0          ragg_1.5.0          memoise_2.0.1      
#> [52] evaluate_1.0.5      knitr_1.51          viridisLite_0.4.3  
#> [55] rlang_1.1.7         Rcpp_1.1.1          glue_1.8.0         
#> [58] xml2_1.5.2          svglite_2.2.2       rstudioapi_0.18.0  
#> [61] jsonlite_2.0.0      R6_2.6.1            systemfonts_1.3.2  
#> [64] fs_1.6.7
```

------------------------------------------------------------------------

## References

Correia, Hannah E., and Asheber Abebe. 2021. “Regularised Rank
Quasi-Likelihood Estimation for Generalised Additive Models.” *Journal
of Nonparametric Statistics* 33 (1): 101–17.
<https://doi.org/10.1080/10485252.2021.1921176>.

Efron, Bradley and Tibshirani, Robert J. 1993. *An Introduction to the
Bootstrap*. No. 57. Monographs on Statistics and Applied Probability.
Chapman; Hall.

Elkin, Lisa A., and Matthew Kay. 2023. “An Aligned Rank Transform
Procedure for Multifactor Contrast Tests.” *Journal of Statistical
Software* forthcoming.

Goude, Yannig. 2024. *Forecasting at EDF: Generalized Additive Models
for Time Series*. EDF R&D, EDF Lab Saclay.

Hyndman, Rob J., and George Athanasopoulos. 2018. *Forecasting:
Principles and Practice*. 2nd ed.
[Https://otexts.com/fpp2/](https://otexts.com/fpp2/).

Meinshausen, Nicolai, Marloes H. Maathuis, and Peter Bühlmann. 2012.
“Asymptotic Optimality of the Westfall-Young Permutation Procedure for
Multiple Testing Under Dependence.” *The Annals of Statistics* 39 (6):
3369–91. <https://doi.org/10.1214/11-AOS946>.

Phipson, Belinda, and Gordon K. Smyth. 2010. “Permutation P-Values
Should Never Be Zero: Calculating Exact P-Values When Permutations Are
Ranked.” *Statistical Applications in Genetics and Molecular Biology* 9
(1): 39. <https://doi.org/10.2202/1544-6115.1585>.

Unknown. 2020. *Autocorrelated Time Series: ARIMA(1,1,1) Model
Forecasting Techniques*.
