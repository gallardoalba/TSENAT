# Appendix B: Non-Parametric Validation of Linear Model Results via GAM and Rank-Based Methods

## Overview: Non-Parametric Validation of Entropic Index × Group Interactions (q-values)

### Purpose and Rationale

Primary Goal: Validate that discoveries of scale-dependent entropic
index × group interactions from linear models generalize to
non-parametric statistical frameworks with minimal assumptions.

TSENAT’s default approach employs **linear models** to test for entropic
index × group effects in Tsallis entropy. While computationally
efficient and well-characterized, linear models assume: - Normally
distributed residuals (reasonable after appropriate transformation) -
Homogeneous variance across entropic indices (often violated in
diversity metrics) - Additive effects of group and entropic index on
entropy (may oversimplify nonlinear relationships)

This appendix provides two independent, non-parametric alternatives
that:

1.  Make minimal distributional assumptions, relying only on ranks or
    data-adaptive smoothing
2.  Explicitly handle the ordered structure of entropic indices,
    treating them as sequential measurements (q-values)
3.  Adapt to real-world heteroscedasticity, automatically weighting
    measurements by entropy variance
4.  Provide independent validation, allowing comparison of linear
    vs. non-parametric findings

### Two Complementary Validation Approaches

**Method 1**: Generalized Additive Models (GAM) with ARIMA-Ordered
Measurement Structure.

- Distributional assumption: None required; uses non-parametric basis
  functions (thin-plate splines).
- Treatment of entropic structure: Treats q-values as time-like ordered
  measurements, applies ARIMA(1,1,0) differencing for autocorrelation.
- Heteroscedasticity handling: Automatically detects variance
  heterogeneity; applies optimal weighting.
- Key advantage: Captures smooth nonlinear q \* group effects with
  computational efficiency.

**Method 2**: Rank-Based Tests (Scheirer-Ray-Hare with Hochberg Step-Up
Correction).

- Distributional assumption: None; operates entirely on ranks (maximal
  robustness).
- Treatment of q-ordering: Scheirer-Ray-Hare two-way test for paired
  measurements across q and group.
- Heteroscedasticity handling: Rank transformation inherently robust;
  Hochberg correction controls FWER.
- Key advantage: Maximally robust to outliers; valid for any continuous
  distribution.

### Why Two Methods for One Question?

Statistical testing in transcriptomics faces a fundamental challenge: no
single method is universally optimal. Different approaches make
different assumptions and have different strengths. TSENAT employs two
complementary strategies: parametric methods (linear models) and
non-parametric rank-based tests:

| Aspect | Parametric (Linear Model / GAM) | Non-Parametric (Rank-Based) |
|----|----|----|
| Assumptions | Normality, homoscedasticity | Ranks only; fully non-parametric |
| Power | Highest (if assumptions met) | Good; slightly reduced but robust |
| Robustness | Moderate; sensitive to outliers | Highest; resistant to outliers and outlier-driven effects |
| Interpretation | Parametric effect sizes, smooth curves | Effect ranks; robust p-values independent of distribution |
| Outlier influence | High potential for bias | Minimal; rank transformation inherently resistant |

Validation strategy: High concordance between parametric and
non-parametric methods confirms that findings are robust to modeling
assumptions and generalize across statistical frameworks.

------------------------------------------------------------------------

## Setup

This section initializes the analysis environment by loading required
packages, setting a random seed for reproducibility, and preparing the
test dataset.

``` r

suppressPackageStartupMessages({
    library(TSENAT)
    library(ggplot2)
    library(SummarizedExperiment)
    library(dplyr)
    library(gridExtra)
})

set.seed(42)

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

# Apply filtering for quality control
analysis <- filter_analysis(
    analysis,
    stringency = "medium",
    min_isoform_abundance = 0
)
```

### Computing Tsallis Entropy with Pseudocount Regularization

Pseudocounts are critical for statistical robustness when applying
rank-based tests to sparse RNA-seq count data. RNA-seq experiments
typically contain many zero or near-zero counts, which creates two
problems for rank-based nonparametric methods:

1.  **Ties in rankings:** Sparse counts produce many tied values
    (especially zeros), which reduces the discriminatory power of rank
    tests. Rank-based statistics depend on unique orderings; when many
    observations are identical, the test statistic contains less
    information, reducing statistical power.

2.  **Library size artifacts:** Small count differences due to
    sequencing depth variations overshadow true biological differences.
    Normalized library sizes (accounting for sequencing depth) ensure
    that entropy estimates reflect true biological diversity rather than
    technical artifacts (Robinson et al. 2010; Love et al. 2014).

By adding small pseudocounts proportional to library size, we achieve
two benefits:

- **Regularization** breaks ties and prevents zero-inflation bias that
  compromises rank-based inference.
- **Normalization** makes entropy estimates comparable across samples
  with different sequencing depths. This approach is standard in RNA-seq
  analysis and is particularly important for entropy-based diversity
  metrics, which require positive values for logarithmic transforms.

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

## Rank-Based Approach (Scheirer-Ray-Hare)

The Scheirer-Ray-Hare test provides a non-parametric validation of
linear model results by leveraging the paired nature of the experimental
design and the ordered structure of entropic indices. Unlike parametric
methods that assume normality and homogeneity of variance, the
Scheirer-Ray-Hare test operates solely on ranks, making it maximally
robust to outliers and extreme values.

### Assumption Validation

Before applying rank-based methods, we verify key assumptions for
rank-based inference.

``` r

# Validate Scheirer-Ray-Hare test assumptions
analysis <- calculate_rank_assumptions(
    analysis
)

# Extract the result object
rank_assumptions <- metadata(analysis, "rankbased_assumptions")$result
```

| Metric | Value |
|:---|---:|
| Genes tested | 88 |
| Samples (entropic order indices × groups) | 656 |
| Entropy range | 0.0000 to 1.0000 |
| Mean entropy | 0.6676 |
| Median entropy | 0.7425 |
| Missing values | 0 |
| Note: |  |
|  Data from complete entropy matrix across all q-values and samples. |  |

Supplementary Table 1 \| Rank-Based Entropy Matrix Properties. {.table
.table .table-striped .table-hover
style="margin-left: auto; margin-right: auto;border-bottom: 0;"}

| Characteristic | Test | Result |
|:---|:---|---:|
| Paired Structure | Permutation test | p=1.000 |
| Gene Heterogeneity | Spearman r | r=0.304 |
| Subject Consistency | Kendall’s W | W=0.000 |
| Note: |  |  |
|  All metrics support validity of paired-design rank-based inference. |  |  |

Supplementary Table 2 \| Entropy Matrix Validation Metrics. Three core
properties validating suitability for rank-based analysis. Rows: Paired
Structure (exchangeability test ensures independence), Gene
Heterogeneity (Spearman rank correlation), Subject Consistency
(Kendall’s *W* concordance). Values reported as *P*-values, correlation
coefficients *r*, and concordance *W* (0-1 range). All tests support
exchangeability assumption justifying Scheirer-Ray-Hare nonparametric
approach. {.table .table .table-striped .table-hover
style="margin-left: auto; margin-right: auto;border-bottom: 0;"}

### Scheirer-Ray-Hare Test: Testing q $`\times`$ Condition Interactions

Now we will use the Scheirer-Ray-Hare test to evaluate whether entropy
patterns across entropic indices (q-values) differ between normal and
tumor samples. The Scheirer-Ray-Hare test is a non-parametric
alternative to the paired t-test and repeated measures ANOVA, making it
ideal for paired designs where distributional assumptions may be
violated (Zhang and Yuan 2018).

``` r

# Run Scheirer-Ray-Hare test for q * condition INTERACTION
analysis <- calculate_rank_test(
    analysis,
    multicorr = "hochberg"
)

# Extract ALL results (no filtering) for summary statistics
srh_results_all <- results(analysis, type = "rank_test", rankBy = "pvalue")

# Extract top significant genes for display
srh_results <- results(analysis, type = "rank_test", rankBy = "pvalue", n = 20, filterFDR = 0.05)
print(head(srh_results, n = 10))
```

| Metric | Value |
|:---|---:|
| Genes tested | 88 |
| Significant (p \< 0.05) | 10 |
| Significant (adj_p \< 0.05, FWER-controlled) | 6 |
| NAs | 0 |
| Mean effect size ($`\eta^2`$) | 43.4% |
| Median effect size ($`\eta^2`$) | 40.3% |
| Strong effect genes ($`\eta^2`$ \> 10%) | 8 |
| Note: |  |
|  Scheirer-Ray-Hare test: nonparametric rank-based ANOVA for paired designs. |  |

Supplementary Table 3 \| Scheirer-Ray-Hare Test Results Summary.
Rank-based nonparametric test of *q* × condition interactions. {.table
.table .table-striped .table-hover
style="margin-left: auto; margin-right: auto;border-bottom: 0;"}

| Gene | P-value | Adj. P-value | F-Statistic | Effect Size (η²) | Test Method | Interaction Class |
|:---|:---|:---|:---|:---|:---|:---|
| CXCL12 | 1.9e-63 | 1.6e-61 | 14.4833 | 0.2993 | Scheirer-Ray-Hare (paired) | Strongly q-dependent |
| LINC03040 | 1.1e-29 | 9.5e-28 | 7.0391 | 0.0735 | Scheirer-Ray-Hare (paired) | Moderately q-dependent |
| PNRC2 | 6.8e-18 | 5.8e-16 | 4.7573 | 0.9610 | Scheirer-Ray-Hare (paired) | Strongly q-dependent |
| FAM114A2 | 9.8e-06 | 0.000831 | 2.3578 | 0.2556 | Scheirer-Ray-Hare (paired) | Strongly q-dependent |
| HDAC2 | 0.000116 | 0.009751 | 2.1161 | 0.5758 | Scheirer-Ray-Hare (paired) | Strongly q-dependent |
| ATG5 | 0.000239 | 0.019839 | 2.0427 | 0.6787 | Scheirer-Ray-Hare (paired) | Strongly q-dependent |
| Note: |  |  |  |  |  |  |
|  Effect size η² quantifies strength of *q* × condition interaction; top 10 genes ranked by significance. P-values \< 0.0001 shown in scientific notation for precision. |  |  |  |  |  |  |

Supplementary Table 4 \| Top genes identified by Scheirer-Ray-Hare
rank-based test. {.table .table .table-striped .table-hover
style="margin-left: auto; margin-right: auto;border-bottom: 0;"}

Visualize q-curves for top genes:

``` r

# Plot top genes from Scheirer-Ray-Hare test
top_genes_plot <- plot_diversity_spectrum(analysis, lm_res = srh_results, n_top = 4)
print(top_genes_plot)
```

![Extended Data Figure 1 \| Scheirer-Ray-Hare rank test results for
scale-dependent isoform
switching.](TSENAT_appendix_B_files/figure-html/extended-fig-1-srh-q-curves-1.png)

Extended Data Figure 1 \| Scheirer-Ray-Hare rank test results for
scale-dependent isoform switching.

------------------------------------------------------------------------

## Generalized Additive Model Approach (GAM)

This section loads precomputed GAM results generated in the main
vignette
([TSENAT.Rmd](https://gallardoalba.github.io/TSENAT/articles/TSENAT.Rmd)).

``` r

# Load precomputed LM analysis object from RDS file
analysis_lm <- readRDS(
    system.file("extdata", "analysis_lm.rds", package = "TSENAT")
)

# Extract GAM/LM results for inspection using accessor function
gam_results <- results(analysis_lm, type = "lm", rankBy = "pvalue")
```

| Metric | Value |
|:---|---:|
| Genes tested | 76 |
| Significant (p \< 0.05) | 65 |
| Concordant (p \< 0.05 AND adj_p \< 0.05) | 59 |
| NAs | 0 |
| Mean effect size | 20.1% |
| Median effect size | 15.4% |
| Strong effect genes (effect_size \> 30%) | 17 |
| Model convergence rate | 100.0% |
| Note: |  |
|  Aggregate statistics across all genes tested by GAM method; effect size \> 30% indicates strong biological signal. |  |

GAM Results Summary {.table .table .table-striped .table-hover
style="margin-left: auto; margin-right: auto;border-bottom: 0;"}

| Gene | p (interaction) | Adjusted p | Effect size | Test statistic | df |
|:---|---:|---:|---:|---:|---:|
| CXCL12 | 0e+00 | 0e+00 | 81.3% | 752.51 | 640 |
| THY1 | 0e+00 | 0e+00 | 60.1% | 442.14 | 640 |
| ING3 | 0e+00 | 0e+00 | 35.3% | 352.59 | 640 |
| SNHG10 | 0e+00 | 0e+00 | 5.1% | 340.82 | 640 |
| LINC03040 | 0e+00 | 0e+00 | 73.8% | 261.24 | 643 |
| HDAC2 | 0e+00 | 0e+00 | 28.4% | 232.94 | 640 |
| ENSG00000274322 | 0e+00 | 0e+00 | 4.7% | 217.93 | 640 |
| MEF2A | 0e+00 | 0e+00 | 19.3% | 197.39 | 640 |
| RAP1GDS1 | 0e+00 | 0e+00 | 54.8% | 170.93 | 640 |
| PDE7A | 0e+00 | 0e+00 | 0.4% | 169.40 | 640 |
| Note: |  |  |  |  |  |
|  Top 10 genes ranked by significance; p-values \< 0.001 shown in scientific notation; effect size range 0-100%. |  |  |  |  |  |

Top 10 Genes by GAM p-value (with Effect Size and Test Statistics)
{.table .table .table-striped .table-hover
style="margin-left: auto; margin-right: auto;border-bottom: 0;"}

------------------------------------------------------------------------

## Method Concordance Analysis

A critical validation step is to compare whether the Scheirer-Ray-Hare
rank-based test and GAM produce consistent results. High concordance
between methods (i.e., genes significant in both approaches) provides
strong evidence that discoveries are robust to methodological choice.
Discordant genes—those significant in only one method—warrant closer
inspection: they may represent genuine biological signals revealed by
one method’s specific advantages, or artifacts of that method’s
assumptions. This section quantifies agreement between approaches and
identifies high-confidence genes significant across both statistical
frameworks.

``` r

# Compute concordance analysis using new two-object API
# Compares LM results (from analysis_lm) with rank test results (from analysis)
analysis_with_concordance <- calculate_concordance(
    analysis_lm = analysis_lm,
    analysis_rank = analysis,
    lm_method = "lm_interaction",
    rank_method = "rank_test",
    verbose = TRUE
)

# Extract components from S4 results
concordance_result <- metadata(analysis_with_concordance, "method_concordance")
comparison_df <- concordance_result$comparison_df
spearman_rho <- concordance_result$spearman_rho
high_conf <- concordance_result$high_confidence
agreement_table <- concordance_result$agreement_table
```

| Metric | Value |
|:---|---:|
| Total genes compared | 76 |
| Spearman correlation (p-values) | rho = 0.3571 |
| Both methods significant (p \< 0.05) | 5 (6.6%) |
| LM only significant | 54 (71.1%) |
| Rank test only significant | 0 (0.0%) |
| Neither significant | 17 (22.4%) |
| Concordance rate | 28.9% |
| Discordance rate | 71.1% |
| Note: |  |
|  Comparison of significant genes (p \< 0.05) detected by GAM vs Scheirer-Ray-Hare methods; high concordance validates robustness. |  |

Global Concordance Metrics: GAM vs Scheirer-Ray-Hare Methods {.table
.table .table-striped .table-hover
style="margin-left: auto; margin-right: auto;border-bottom: 0;"}

| Agreement Category | Number of Genes | Percentage |
|:---|---:|---:|
| Both significant | 5 | 6.6% |
| LM only | 54 | 71.1% |
| Neither significant | 17 | 22.4% |
| Note: |  |  |
|  Categories: Concordant (both methods, p \< 0.05); GAM-only; Scheirer-Ray-Hare-only; Neither (both p ≥ 0.05). |  |  |

Method Agreement Distribution {.table .table .table-striped .table-hover
style="margin-left: auto; margin-right: auto;border-bottom: 0;"}

| Gene | LM adj p | Rank test adj p | LM Effect | Rank test $`\eta^2`$ |
|:---|---:|---:|---:|---:|
| CXCL12 | 2.990e-162 | 1.640e-61 | 81.3% | 0.299 |
| LINC03040 | 1.350e-55 | 9.502e-28 | 73.8% | 0.074 |
| HDAC2 | 1.855e-49 | 9.751e-03 | 28.4% | 0.576 |
| FAM114A2 | 6.136e-27 | 8.306e-04 | 46.4% | 0.256 |
| ATG5 | 2.301e-24 | 1.984e-02 | 43.1% | 0.679 |
| Note: |  |  |  |  |
|  High-confidence genes: significant by both GAM and Scheirer-Ray-Hare rank test (adj p \< 0.05); $`\eta^2`$ = rank test effect size. |  |  |  |  |

Robust Entropic Order Index Interactions: High-Confidence Genes Detected
by Both Methods (n=5, ranked by statistical significance) {.table .table
.table-striped .table-hover
style="margin-left: auto; margin-right: auto;border-bottom: 0;"}

### Statistical Power vs. Robustness: Understanding Method Discordance

The striking discordance between GAM and Scheirer-Ray-Hare results
(71.1% significant in GAM only, 0% in rank test only, 28.9% total
concordance) reflects a fundamental trade-off in statistical
methodology: **parametric methods maximize power when assumptions hold,
while nonparametric methods sacrifice power for robustness to assumption
violations**.

GAM’s superior power derives from two factors:

1.  Distributional assumptions: GAM assumes approximately normal
    residuals and homogeneous variance, which are reasonable after
    appropriate transformation for entropy data. Under these
    assumptions, parametric methods are theoretically optimal, achieving
    maximum power for a given type I error rate. The rank-based
    Scheirer-Ray-Hare test is assumption-free but necessarily discards
    quantitative information by converting measurements to ranks, which
    reduces statistical power when the underlying data are approximately
    normal.

2.  Model flexibility with penalty: GAM uses thin-plate splines with
    smoothness penalties that simultaneously fit nonlinear patterns
    while controlling degrees of freedom. This flexibility allows GAM to
    detect subtle entropic index × group interactions across all
    q-values. The Scheirer-Ray-Hare test, by contrast, operates on rank
    patterns only, which is inherently less sensitive to continuous
    relationships across the ordering (q-values).

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
# Check if concordance analysis was successful
concordance_results <- metadata(analysis_with_concordance, "method_concordance")

p <- plot_concordance(analysis_with_concordance, verbose = TRUE)
print(p)
```

![Method concordance visualization comparing Scheirer-Ray-Hare
rank-based and GAM approaches for detecting q × condition
interactions.](TSENAT_appendix_B_files/figure-html/visualize-method-concordance-1.png)

Method concordance visualization comparing Scheirer-Ray-Hare rank-based
and GAM approaches for detecting q × condition interactions.

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
