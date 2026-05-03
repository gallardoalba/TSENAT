# Appendix A: TSENAT vs SplicingFactory Validation

## Introduction

This appendix validates that **TSENAT’s transcript diversity measures
(Shannon entropy, Simpson index) are mathematically equivalent to
SplicingFactory’s implementations**. We benchmark both packages on TCGA
BRCA RNA-seq data to demonstrate:

- Shannon equivalence: TSENAT Tsallis q=1 ~ SplicingFactory Shannon
- Simpson equivalence: TSENAT Tsallis q=2 ~ SplicingFactory Simpson  
- Identical statistical results: Same gene rankings, p-values, fold
  changes

**Who needs this appendix?** Anyone evaluating TSENAT for transcript
diversity analysis or validating cross-package comparisons.

------------------------------------------------------------------------

### Document Structure

- Background: Conceptual overview of transcript diversity and the
  Tsallis framework
- Benchmarking (main section): Hands-on comparison using both packages
- Results: Side-by-side analysis tables (Shannon and Simpson)
- Analysis: Technical explanation of equivalence, normalization
  differences, and when to use each method
- Conclusion

------------------------------------------------------------------------

### Background: Transcript Diversity Concepts

**Transcript diversity** measures how variable isoform expression is
within genes. Two standard approaches:

- Shannon entropy: Overall diversity (number and evenness of isoforms)
- Simpson index: Dominance (whether one isoform monopolizes expression)

TSENAT’s approach: Generalized Tsallis entropy framework where q=1 gives
Shannon and q=2 gives Simpson. This allows sensitivity analysis across
the q continuum while maintaining equivalence with established methods.

**Key question**: Does TSENAT’s mathematical generalization produce
identical statistical results?

------------------------------------------------------------------------

### Dataset and Methods

We use TCGA BRCA RNA-seq data as our validation dataset, which is
included in the **SplicingFactory package**:

- Public, reproducible resource.
- Multiple transcript isoforms per gene.
- Real biological signal.
- Variable sequencing depth and isoform abundance.

Analysis pipeline:

1.  Load transcript-level read counts.
2.  Filter genes.
3.  Compute entropy values.
4.  Statistical testing (Wilcoxon rank-sum) for group differences.
5.  Compare results (concordance, effect sizes, rankings).

------------------------------------------------------------------------

## Benchmarking

#### Importing example data

We begin by loading the required libraries and the TCGA BRCA dataset
that will serve as our validation benchmark.

``` r

suppressPackageStartupMessages({
    library("SplicingFactory")
    library(TSENAT)
    library("SummarizedExperiment")
})
set.seed(12345)

# Load dataset
data(tcga_brca_luma_dataset)

# Extract gene names
genes <- tcga_brca_luma_dataset[, 1]

# Extract read count data without gene names
readcounts <- tcga_brca_luma_dataset[, -1]
```

#### Data filtering and preprocessing

We filter genes to retain only those with sufficient coverage across
samples, ensuring robust diversity estimates.

``` r

tokeep <- rowSums(readcounts > 5) > 5
readcounts <- readcounts[tokeep, ]
genes      <- genes[tokeep]

# Create SummarizedExperiment object for cleaner data handling
# Important: Set gene names as rownames for compatibility with SplicingFactory
# TSENATAnalysis requires 'gene_id' in rowData
se_data <- SummarizedExperiment(
  assays = list(counts = as.matrix(readcounts)),
  rowData = DataFrame(gene_id = genes)
)
rownames(se_data) <- genes
```

#### Transcript diversity calculation

To calculate transcript diversity measures (Shannon entropy, Simpson
index, and Tsallis entropy at q=1 and q=2), use:

``` r

# Create sample group classification based on sample names
# (Normal samples end with _N, Tumor samples don't)
sample_group <- ifelse(grepl("_N$", colnames(se_data)), "Normal", "Tumor")

# Add group metadata to se_data for all downstream analyses
colData(se_data)$group <- sample_group

# SplicingFactory: Calculate diversity for both Shannon (naive) and Simpson methods
sf_methods <- c("naive", "simpson")
splicing_results <- setNames(vector("list", length(sf_methods)), sf_methods)

for (method in sf_methods) {
  splicing_results[[method]] <- SplicingFactory::calculate_diversity(
    x = se_data,
    method = method,
    norm = TRUE,
    verbose = FALSE
  )
  # Add sample_type and group metadata to SplicingFactory objects
  # (sample_type matches SplicingFactory naming; group matches TSENAT naming for consistency)
  colData(splicing_results[[method]])$sample_type <- sample_group
  colData(splicing_results[[method]])$group <- sample_group
}

# TSENAT: Single config for q-value sweep
# (Configuration is identical for all q-values)
config <- TSENAT_config(condition_col = "group")

# TSENAT: Create analysis objects and compute diversity for q=1 and q=2
q_values <- c(1, 2)
tsenat_results <- setNames(vector("list", length(q_values)), paste0("q", q_values))

for (q in q_values) {
  # Create analysis object (following Bioconductor pattern: config first, then constructor)
  tsenat_results[[paste0("q", q)]] <- TSENATAnalysis(se_data, config = config)
  
  # Compute diversity with dynamic output filename using sprintf
  # Use explicit namespace to avoid SplicingFactory::calculate_diversity masking
  tsenat_results[[paste0("q", q)]] <- TSENAT::calculate_diversity(
    analysis = tsenat_results[[paste0("q", q)]],
    q = q,
    norm = TRUE,
    verbose = FALSE
  )
}

# Extract for downstream use (aliases maintain backward compatibility with existing code)
shannon_entropy <- splicing_results$naive
simpson_index <- splicing_results$simpson
tsenat_analysis_q1 <- tsenat_results$q1
tsenat_analysis_q2 <- tsenat_results$q2
```

Both packages return a `SummarizedExperiment` object, that you can
investigate further with the `assay` function.

#### Differential analysis

We apply statistical testing to identify genes with significantly
different transcript diversity between Normal and Tumor samples using
both SplicingFactory and TSENAT frameworks.

``` r

# SplicingFactory: Differential analysis for Shannon entropy
entropy_significance <- SplicingFactory::calculate_difference(
  x = shannon_entropy,
  samples = "sample_type",
  control = "Normal",
  method = "mean",
  test = "wilcoxon",
  verbose = FALSE
)

# SplicingFactory: Differential analysis for Simpson index
simpson_significance <- SplicingFactory::calculate_difference(
  x = simpson_index,
  samples = "sample_type",
  control = "Normal",
  method = "mean",
  test = "wilcoxon",
  verbose = FALSE
)


# TSENAT: Differential analysis for Tsallis q=1 (Shannon equivalent)
# Extract the diversity results for q=1 using the results() accessor
# Returns SummarizedExperiment with diversity values (genes x samples)
div_result_q1 <- results(tsenat_analysis_q1, type = "diversity", q = 1.0, format = "se")

# Apply SplicingFactory's calculate_difference to TSENAT q=1 diversity
# Use "group" column already in colData (from se_data)
tsenat_shannon_diff <- SplicingFactory::calculate_difference(
  x = div_result_q1,
  samples = "group",
  control = "Normal",
  method = "mean",
  test = "wilcoxon",
  verbose = FALSE
)


# TSENAT: Differential analysis for Tsallis q=2 (Simpson equivalent)
# Extract the diversity results for q=2 using the results() accessor
# Returns SummarizedExperiment with diversity values (genes x samples)
div_result_q2 <- results(tsenat_analysis_q2, type = "diversity", q = 2.0, format = "se")

# Apply SplicingFactory's calculate_difference to TSENAT q=2 diversity
# Use "group" column already in colData (from se_data)
tsenat_simpson_diff <- SplicingFactory::calculate_difference(
  x = div_result_q2,
  samples = "group",
  control = "Normal",
  method = "mean",
  test = "wilcoxon",
  verbose = FALSE
)
```

For the benchmark comparison with TSENAT, we will compute Tsallis
entropy at q=1 (equivalent to Shannon) and q=2 (equivalent to Simpson).
This allows direct comparison with the SplicingFactory results shown
above.

## Benchmark Results Summary

#### Shannon Entropy (SplicingFactory Laplace vs TSENAT Tsallis q=1)

| Method | Genes Tested | Significant (padj\<0.05) | Mean log2FC | SD log2FC |
|:---|---:|---:|---:|---:|
| SplicingFactory (Shannon/naive) | 214 | 26 | -0.0009354 | 1.582967 |
| TSENAT (Tsallis q=1) | 213 | 26 | 0.0003618 | 1.586675 |

**Supplementary Table 5 \| Shannon Entropy Analysis Summary
Statistics.** Differential analysis comparing Normal (N=8) vs Tumor
(N=8) samples using SplicingFactory’s
[`calculate_difference()`](https://rdrr.io/pkg/SplicingFactory/man/calculate_difference.html)
method. {.table}

**Interpretation**: Both methods similarly identify significantly
different transcript diversity between Normal and Tumor samples. The
comparable number of significant genes (padj\<0.05) and mean log2FC
values demonstrate equivalent statistical power and effect size
detection.

##### Top 10 Significant Genes - SplicingFactory (Laplace)

| Gene     | Normal Mean | Tumor Mean | Mean Diff |  log2FC |  p-value | adj p-value |
|:---------|------------:|-----------:|----------:|--------:|---------:|------------:|
| C1orf213 |      0.8122 |     0.1400 |   -0.6722 | -2.5365 | 3.97e-07 |    8.50e-05 |
| HAPLN3   |      0.7030 |     0.4299 |   -0.2731 | -0.7096 | 4.17e-05 |    4.46e-03 |
| COL1A2   |      0.0000 |     0.1643 |    0.1643 |     Inf | 6.68e-05 |    4.77e-03 |
| F10      |      0.1311 |     0.3720 |    0.2409 |  1.5050 | 1.37e-04 |    7.31e-03 |
| DUSP14   |      0.2989 |     0.1302 |   -0.1688 | -1.1995 | 1.89e-04 |    8.10e-03 |
| MBD2     |      0.0027 |     0.0371 |    0.0344 |  3.7876 | 2.39e-04 |    8.54e-03 |
| C1orf86  |      0.1479 |     0.0029 |   -0.1450 | -5.6777 | 4.59e-04 |    9.85e-03 |
| GFPT1    |      0.0231 |     0.0025 |   -0.0207 | -3.2253 | 3.80e-04 |    9.85e-03 |
| HNRNPR   |      0.5162 |     0.5907 |    0.0746 |  0.1947 | 4.60e-04 |    9.85e-03 |
| OSR1     |      0.0497 |     0.2891 |    0.2394 |  2.5411 | 4.10e-04 |    9.85e-03 |

**Supplementary Table 6a \| SplicingFactory method - Top 10 genes
identified by Shannon entropy (naive mode).** Ranked by adjusted
*P*-value. {.table}

##### Top 10 Significant Genes - TSENAT (Tsallis q=1)

| Gene     | Normal Mean | Tumor Mean | Mean Diff |  log2FC |  p-value | adj p-value |
|:---------|------------:|-----------:|----------:|--------:|---------:|------------:|
| C1orf213 |      0.8122 |     0.1400 |   -0.6722 | -2.5365 | 3.97e-07 |    8.46e-05 |
| HAPLN3   |      0.7030 |     0.4299 |   -0.2731 | -0.7096 | 4.17e-05 |    4.44e-03 |
| COL1A2   |      0.0000 |     0.1643 |    0.1643 |     Inf | 6.68e-05 |    4.74e-03 |
| F10      |      0.1311 |     0.3720 |    0.2409 |  1.5050 | 1.37e-04 |    7.27e-03 |
| DUSP14   |      0.2989 |     0.1302 |   -0.1688 | -1.1995 | 1.89e-04 |    8.06e-03 |
| MBD2     |      0.0027 |     0.0371 |    0.0344 |  3.7876 | 2.39e-04 |    8.50e-03 |
| C1orf86  |      0.1479 |     0.0029 |   -0.1450 | -5.6777 | 4.59e-04 |    9.80e-03 |
| HNRNPR   |      0.5162 |     0.5907 |    0.0746 |  0.1947 | 4.60e-04 |    9.80e-03 |
| OSR1     |      0.0497 |     0.2891 |    0.2394 |  2.5411 | 4.10e-04 |    9.80e-03 |
| GFPT1    |      0.0231 |     0.0025 |   -0.0207 | -3.2253 | 3.80e-04 |    9.80e-03 |

**Supplementary Table 6b \| TSENAT method - Top 10 genes identified by
Tsallis entropy at q=1 (Shannon equivalence).** Results obtained using
SplicingFactory’s
[`calculate_difference()`](https://rdrr.io/pkg/SplicingFactory/man/calculate_difference.html)
method applied to TSENAT Tsallis q=1 diversity values. Ranked by
adjusted *P*-value. {.table}

------------------------------------------------------------------------

#### Simpson Index (SplicingFactory vs TSENAT Tsallis q=2)

| Method | Genes Tested | Significant (padj\<0.05) | Mean log2FC | SD log2FC |
|:---|---:|---:|---:|---:|
| SplicingFactory (Simpson) | 214 | 26 | -0.0012687 | 1.812772 |
| TSENAT (Tsallis q=2) | 213 | 26 | -0.0000521 | 1.817060 |

**Supplementary Table 8 \| Simpson Index Analysis Summary Statistics.**
Differential analysis using Simpson diversity metric (Gini-Simpson
index) comparing Normal (N=8) vs Tumor (N=8) samples using
SplicingFactory’s
[`calculate_difference()`](https://rdrr.io/pkg/SplicingFactory/man/calculate_difference.html)
method. {.table}

**Interpretation**: Simpson index results parallel Shannon entropy
findings, with both methods identifying similar numbers of significantly
different isoform dominance patterns. The agreement across both
diversity measures provides strong evidence for TSENAT’s mathematical
and statistical equivalence to SplicingFactory.

##### Top 10 Significant Genes - SplicingFactory (Simpson)

| Gene     | Normal Mean | Tumor Mean | Mean Diff |  log2FC |  p-value | adj p-value |
|:---------|------------:|-----------:|----------:|--------:|---------:|------------:|
| C1orf213 |      0.3852 |     0.0515 |   -0.3337 | -2.9037 | 3.97e-07 |    8.50e-05 |
| HAPLN3   |      0.4823 |     0.2558 |   -0.2266 | -0.9152 | 5.17e-06 |    5.53e-04 |
| COL1A2   |      0.0000 |     0.0626 |    0.0626 |     Inf | 6.68e-05 |    4.77e-03 |
| F10      |      0.0697 |     0.2494 |    0.1797 |  1.8400 | 1.10e-04 |    5.86e-03 |
| DUSP14   |      0.1686 |     0.0617 |   -0.1069 | -1.4499 | 2.11e-04 |    7.32e-03 |
| HNRNPR   |      0.5376 |     0.6134 |    0.0758 |  0.1903 | 1.79e-04 |    7.32e-03 |
| MBD2     |      0.0005 |     0.0107 |    0.0102 |  4.5028 | 2.39e-04 |    7.32e-03 |
| CXorf40A |      0.7455 |     0.7939 |    0.0484 |  0.0907 | 2.75e-04 |    7.34e-03 |
| GFPT1    |      0.0051 |     0.0004 |   -0.0047 | -3.5271 | 3.80e-04 |    8.78e-03 |
| OSR1     |      0.0124 |     0.1156 |    0.1032 |  3.2168 | 4.10e-04 |    8.78e-03 |

**Supplementary Table 9a \| SplicingFactory method - Top 10 genes by
Simpson index (Gini-Simpson, emphasizes dominant isoforms).** Ranked by
adjusted *P*-value. {.table}

##### Top 10 Significant Genes - TSENAT (Tsallis q=2)

| Gene     | Normal Mean | Tumor Mean | Mean Diff |  log2FC |  p-value | adj p-value |
|:---------|------------:|-----------:|----------:|--------:|---------:|:------------|
| C1orf213 |      0.7703 |     0.1029 |   -0.6674 | -2.9037 | 3.97e-07 | 8.46e-05    |
| HAPLN3   |      0.7235 |     0.3836 |   -0.3398 | -0.9152 | 5.17e-06 | 5.50e-04    |
| COL1A2   |      0.0000 |     0.1251 |    0.1251 |     Inf | 6.68e-05 | 4.74e-03    |
| F10      |      0.1045 |     0.3741 |    0.2696 |  1.8400 | 1.10e-04 | 5.83e-03    |
| HNRNPR   |      0.6144 |     0.7010 |    0.0866 |  0.1903 | 1.79e-04 | 7.28e-03    |
| DUSP14   |      0.2529 |     0.0926 |   -0.1603 | -1.4499 | 2.11e-04 | 7.28e-03    |
| MBD2     |      0.0009 |     0.0213 |    0.0204 |  4.5028 | 2.39e-04 | 7.28e-03    |
| CXorf40A |      0.8387 |     0.8931 |    0.0544 |  0.0907 | 2.75e-04 | 7.31e-03    |
| OSR1     |      0.0249 |     0.2312 |    0.2063 |  3.2168 | 4.10e-04 | 8.74e-03    |
| GFPT1    |      0.0103 |     0.0009 |   -0.0094 | -3.5271 | 3.80e-04 | 8.74e-03    |

**Supplementary Table 9b \| TSENAT method - Top 10 genes by Tsallis
entropy at q=2 (Simpson equivalence).** Results obtained using
SplicingFactory’s
[`calculate_difference()`](https://rdrr.io/pkg/SplicingFactory/man/calculate_difference.html)
method applied to TSENAT Tsallis q=2 diversity values. Ranked by
adjusted *P*-value. {.table}

------------------------------------------------------------------------

### Analysis: Simpson Index and Tsallis Entropy at q=2

Both methods share the same mathematical formula:

``` math
Simpson = Tsallis_{q=2} = 1 - \sum_{i=1}^n p_i^2
```

where $`p_i = \frac{x_i}{\sum_i x_i}`$ (transcript proportions).

Howerver, SplicingFactory returns raw values, while TSENAT normalizes by
the theoretical maximum $`(1-1/n)`$:

|                  | SplicingFactory    | TSENAT                             |
|------------------|--------------------|------------------------------------|
| Formula          | $`1 - \sum p_i^2`$ | $`\frac{1 - \sum p_i^2}{1 - 1/n}`$ |
| Range            | 0 to ~1            | 0 to 1                             |
| Scale difference | Raw                | ~2x for 2-isoform genes            |

This explains observed differences like *C1orf213* Normal: 0.3852
(SplicingFactory) vs 0.7703 (TSENAT). Since C1orf213 has 2 expressed
isoforms, the normalization factor is $`1 - 1/2 = 0.5`$. The
calculation: $`0.3852 / 0.5 = 0.7704`$ confirms this scaling factor.
Despite raw value differences, statistical analysis yields identical
results:

1.  log2 fold changes: Normalization factor cancels in ratios:
    $`\log_2\left(\frac{S_{norm,t}}{S_{norm,n}}\right) = \log_2\left(\frac{S_{raw,t}/K}{S_{raw,n}/K}\right) = \log_2\left(\frac{S_{raw,t}}{S_{raw,n}}\right)`$

2.  P-values and rankings: Wilcoxon rank-sum tests depend only on gene
    ranks, not absolute values, so both methods identify the same
    significant genes

3.  Statistical validity: Both approaches are correct-SplicingFactory
    favors unbounded diversity, TSENAT favors interpretable \[0,1\]
    scaling

------------------------------------------------------------------------

## Conclusion: Equivalence Validation Summary

This appendix has demonstrated that TSENAT’s implementation of
transcript diversity measures is mathematically and statistically
equivalent to SplicingFactory’s established methods for both Shannon
entropy and Simpson index calculations.

------------------------------------------------------------------------

### Session Information

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
#>  [1] knitr_1.51                  dplyr_1.2.1                
#>  [3] SummarizedExperiment_1.40.0 Biobase_2.70.0             
#>  [5] GenomicRanges_1.62.1        Seqinfo_1.0.0              
#>  [7] IRanges_2.44.0              S4Vectors_0.48.1           
#>  [9] BiocGenerics_0.56.0         generics_0.1.4             
#> [11] MatrixGenerics_1.22.0       matrixStats_1.5.0          
#> [13] TSENAT_0.99.0               SplicingFactory_1.18.0     
#> [15] kableExtra_1.4.0           
#> 
#> loaded via a namespace (and not attached):
#>  [1] gtable_0.3.6        xfun_0.57           bslib_0.10.0       
#>  [4] ggplot2_4.0.3       htmlwidgets_1.6.4   lattice_0.22-9     
#>  [7] vctrs_0.7.3         tools_4.5.3         tibble_3.3.1       
#> [10] pkgconfig_2.0.3     pheatmap_1.0.13     Matrix_1.7-5       
#> [13] RColorBrewer_1.1-3  S7_0.2.2            desc_1.4.3         
#> [16] lifecycle_1.0.5     compiler_4.5.3      farver_2.1.2       
#> [19] stringr_1.6.0       textshaping_1.0.5   htmltools_0.5.9    
#> [22] sass_0.4.10         yaml_2.3.12         pkgdown_2.2.0      
#> [25] pillar_1.11.1       jquerylib_0.1.4     tidyr_1.3.2        
#> [28] DelayedArray_0.36.1 cachem_1.1.0        abind_1.4-8        
#> [31] tidyselect_1.2.1    digest_0.6.39       stringi_1.8.7      
#> [34] purrr_1.2.2         cowplot_1.2.0       fastmap_1.2.0      
#> [37] grid_4.5.3          cli_3.6.6           SparseArray_1.10.10
#> [40] magrittr_2.0.5      S4Arrays_1.10.1     withr_3.0.2        
#> [43] scales_1.4.0        rmarkdown_2.31      XVector_0.50.0     
#> [46] otel_0.2.0          ragg_1.5.0          memoise_2.0.1      
#> [49] evaluate_1.0.5      viridisLite_0.4.3   rlang_1.2.0        
#> [52] Rcpp_1.1.1-1.1      glue_1.8.1          xml2_1.5.2         
#> [55] svglite_2.2.2       rstudioapi_0.18.0   jsonlite_2.0.0     
#> [58] R6_2.6.1            systemfonts_1.3.2   fs_2.1.0
```
