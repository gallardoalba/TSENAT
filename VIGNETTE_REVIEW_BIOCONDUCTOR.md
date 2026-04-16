# TSENAT Vignette Review for Bioconductor Submission

**Date:** April 16, 2026  
**Scope:** Review of 3 vignette files against R Markdown best practices
and Bioconductor standards  
**Status:** Comprehensive analysis with actionable recommendations

------------------------------------------------------------------------

## Executive Summary

Your vignettes are **high quality and well-structured**, with excellent
mathematical rigor, comprehensive appendices, and strong biological
motivation. However, **8 specific improvements** would strengthen the
package for Bioconductor review:

| Priority | Category | Issue | Impact | Effort |
|----|----|----|----|----|
| 🔴 HIGH | Best Practice | Vignette metadata incomplete | Package build may fail | 5 min |
| 🔴 HIGH | Bioconductor | Missing session info in main vignette | Reproducibility concern | 2 min |
| 🟡 MEDIUM | Style | Code output suppression inconsistent | Reader experience | 30 min |
| 🟡 MEDIUM | Documentation | Key functions not explicitly documented | Missing API reference | 1 hour |
| 🟡 MEDIUM | Best Practice | Figure captions missing alt-text | Accessibility requirement | 20 min |
| 🟢 LOW | Polish | LaTeX rendering edge cases | PDF output quality | 1 hour |
| 🟢 LOW | Clarity | Some code blocks need comments | Code readability | 30 min |
| 🟢 LOW | Enhancement | Troubleshooting section missing | User support reduction | 1 hour |

------------------------------------------------------------------------

## Detailed Findings

### 1. 🔴 HIGH PRIORITY: Vignette Metadata Completeness

**Files affected:** All three vignettes  
**Bioconductor requirement:** All vignette metadata must be complete and
valid

#### Current YAML Header (TSENAT.Rmd):

``` yaml
---
title: "TSENAT: Tsallis Entropy Analysis Toolbox"
author: "Cristobal Gallardo <gallardoalba@pm.me>"
date: "`r Sys.Date()`"
output:
  rmarkdown::html_document:
    toc: true
  rmarkdown::pdf_document:
    toc: true
vignette: |
  %\VignetteIndexEntry{TSENAT}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
---
```

#### ✅ Recommendations:

**Add to main vignette (TSENAT.Rmd):**

``` yaml
vignette: |
  %\VignetteIndexEntry{TSENAT: Complete Workflow}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
  %\VignetteKeyword{RNASeq}
  %\VignetteKeyword{Transcriptomics}
  %\VignetteKeyword{EntropyAnalysis}
```

**Why:** - Keywords help Bioconductor indexing/discovery - Full title
disambiguates from package name - Bioconductor check requires these
fields

**Time to fix:** 2 minutes per file

------------------------------------------------------------------------

### 2. 🔴 HIGH PRIORITY: Missing Session Info in Main Vignette

**File:** `vignettes/TSENAT.Rmd`  
**Location:** Should be at end of document (currently missing)

#### Current state:

The main vignette ends with “## References” section but **NO session
information block**.

#### ✅ Recommendation:

Add this section before final References:

``` rmarkdown
## Session Information

```{r session-info, results='markup'}
sessionInfo()
```


    #### Why:
    - **Bioconductor requirement:** Reproducibility tracking
    - Shows R version, package versions, platform
    - Helps users debug version-specific issues
    - Both appendices already have this correctly implemented

    **Time to fix:** 2 minutes

    ---

    ### 3. 🟡 MEDIUM PRIORITY: Code Output Suppression Inconsistency

    **Files affected:** All vignettes
    **Current situation:**
    - Some chunks use `results = 'hide'`
    - Others use `eval = FALSE` inappropriately
    - Some use `echo = FALSE` when `results = 'asis'` would be clearer

    #### Current patterns:

    **TSENAT.Rmd line ~35:**
    ```rmarkdown
    ```{r filter-validate-main}
    # Filter lowly-expressed transcripts
    analysis <- filter_analysis(analysis, stringency = "medium")

    ❌ **Problem:** No output control → prints console output in rendered vignette

    **TSENAT_appendix_B.Rmd line ~220:**
    ```rmarkdown
    ```{r srh-rank-test-analysis, message=FALSE, results='hide'}
    # IMPORTANT: Create a copy of analysis before running SRH
    analysis <- calculate_srh(analysis, multicorr = "hochberg")

    ✅ **Good:** Explicitly hides output

    #### ✅ Recommendation - Standardize to this pattern:

    ```rmarkdown
    # For code that produces console output you want hidden:
    ```{r chunk-name, message=FALSE, warning=FALSE, results='hide'}
    # Code here

# For code requiring manual inspection:

`{r chunk-name, message=FALSE, warning=FALSE} # Code here result <- some_function() print(result)`

# For plots (already mostly correct):

`{r fig-chunk, fig.width=10, fig.height=6, out.width="98%"} # Plotting code`


    #### Specific locations needing adjustment:

    | Line Range | Chunk Name | Current | Recommend |
    |-----------|-----------|---------|-----------|
    | ~150-155 | `tsallis-entropy-computation` | Missing output control | Add `results='hide'` |
    | ~280-285 | `test-q-condition-interaction` | No control | Add `message=FALSE, results='hide'` |
    | ~710-715 | `effect-size-merge` | Missing | Add `results='hide'` |

    **Time to fix:** 30 minutes

    ---

    ### 4. 🟡 MEDIUM PRIORITY: Key Functions Not Explicitly Documented in Vignette

    **Issue:** Core functions introduced but not linked to their documentation

    #### Functions introduced without explicit reference to `?function`:

    | Function | First mention | Recommendation |
    |----------|---------------|-----------------|
    | `build_analysis()` | Line ~80 | Add: "See `?build_analysis` for full parameter details" |
    | `filter_analysis()` | Line ~150 | Add: "See `?filter_analysis` for stringency options" |
    | `calculate_diversity()` | Line ~160 | Add: "See `?calculate_diversity` for bootstrap and normalization parameters" |
    | `calculate_lm()` | Line ~280 | Add: "See `?calculate_lm` for available methods (GAM, LMM, GEE, FPCA)" |
    | `calculate_jis()` | Line ~370 | Add: "See `?calculate_jis` for jackknife options" |

    #### ✅ Recommendation:

    After first mention of each major function, add inline reference with backticks:

    **Current (line ~160):**

We now compute Tsallis entropy across your configured q-spectrum.
Diversity measures provide a comprehensive framework for assessing
isoform heterogeneity at multiple scales \[@I007; @I010\].


    **Improved:**

We now compute Tsallis entropy across your configured q-spectrum using
[`calculate_diversity()`](https://gallardoalba.github.io/TSENAT/reference/calculate_diversity.md).
Diversity measures provide a comprehensive framework for assessing
isoform heterogeneity at multiple scales \[@I007; @I010\]. See
[`?calculate_diversity`](https://gallardoalba.github.io/TSENAT/reference/calculate_diversity.md)
for normalization and bootstrap options.


    #### Functions to document:
    1. `build_analysis()` - Line ~80
    2. `filter_analysis()` - Line ~150
    3. `calculate_diversity()` - Line ~160
    4. `calculate_m_estimator()` - Line ~250
    5. `calculate_lm()` - Line ~280
    6. `plot_diversity_spectrum()` - Line ~125
    7. `plot_lm()` - Line ~310
    8. `calculate_jis()` - Line ~370
    9. `calculate_divergence()` - Line ~600
    10. `calculate_effect_sizes()` - Line ~650

    **Time to fix:** 1-2 hours

    ---

    ### 5. 🟡 MEDIUM PRIORITY: Figure Captions Accessibility

    **Files affected:** All vignettes
    **Issue:** Figure captions are descriptive but lack accessibility metadata

    #### Current figure captions (example):
    ```rmarkdown
    ```{r fig-1-isoform-diversity-profiles, fig.width=10, fig.height=5.2,
         out.width="98%", fig.pos="H",
         fig.cap = "**Figure 1:** Isoform diversity profiles..."}

#### ✅ Recommendation - Add alt-text:

``` rmarkdown
```{r fig-1-isoform-diversity-profiles, fig.width=10, fig.height=5.2,
     out.width="98%", fig.pos="H",
     fig.cap = "**Figure 1:** Isoform diversity profiles across entropic indices.
                Normalized Tsallis entropy (0-1) plotted as q-curves for control (blue)
                and treatment (red) samples, showing scale-dependent diversity patterns.",
     fig.alt = "Line plot showing entropy values from q=0 to q=2. Two distinct curves
                per sample, one blue (control) and one red (treatment), demonstrating
                entropy separation across diversity scales."}
```

#### Figures requiring updates:

1.  Figure 1 (isoform diversity profiles) - Line ~445
2.  Figure 2 (scale-dependent genes) - Line ~510
3.  Figure 3 (transcript jackknife influence) - Line ~590
4.  Figure 4 (transcript abundance heatmap) - Line ~620
5.  Figure 5 (divergence distribution) - Line ~730
6.  Figure 6 (global divergence spectrum) - Line ~750

**Why:** - Bioconductor values accessibility compliance (WCAG 2.1 AA) -
Screen readers benefit from alt-text - Reader skimming benefits from
descriptive captions

**Time to fix:** 20 minutes

------------------------------------------------------------------------

### 6. 🟢 LOW PRIORITY: LaTeX Rendering Edge Cases

**Files affected:** TSENAT_appendix_B.Rmd (PDF output)

#### Issue 1: Math delimiters in YAML metadata

**Location:** Line 21

``` yaml
fig.cap = "**Figure 2:** Scale-dependent interaction analysis. 
           GAM-identified genes showing significant q$\\times$condition effects..."
```

**Problem:** LaTeX interpretations of `\times` may vary between PDF
engines  
**Fix:** Use HTML entity in markdown, native LaTeX in math mode:

``` rmarkdown
fig.cap = "**Figure 2:** Scale-dependent interaction analysis.
           GAM-identified genes showing significant q × condition effects
           (p < 0.05 after Benjamini-Hochberg correction)."
```

#### Issue 2: Special characters in mathematical notation

**Location:** Multiple math blocks  
**Current:** Some use `\ne`, `\approx`, etc.  
**Improvement:** Use `≠`, `≈` (unicode) or ensure LaTeX packages loaded

In YAML header, verify:

``` yaml
header-includes:
  - \usepackage{amssymb}
  - \usepackage{amsmath}
```

**Time to fix:** 1 hour (thorough PDF testing)

------------------------------------------------------------------------

### 7. 🟢 LOW PRIORITY: Code Comment Clarity

**Files affected:** All vignettes

#### Example of code that needs inline comments:

**Current (line ~490, TSENAT.Rmd):**

``` r

# Extract ALL results (no filtering) for summary statistics
analysis <- calculate_srh(
    analysis,
    multicorr = "hochberg"
)
srh_results_all <- results(analysis, type = "rank_test", rankBy = "pvalue")
```

**Improved:**

``` r

# Extract ALL results (no filtering) for summary statistics
analysis <- calculate_srh(
    analysis,
    multicorr = "hochberg"  # Step-up correction controls family-wise error rate
)
# Retrieve all genes tested, ranked by raw p-value for summary statistics
srh_results_all <- results(analysis, type = "rank_test", rankBy = "pvalue")
```

#### Locations needing comments:

1.  **TSENAT.Rmd:**
    - Lines ~100-105: config setup rationale
    - Lines ~250-260: M-estimator parameters explanation
    - Lines ~600-610: divergence computation parameters
2.  **TSENAT_appendix_A.Rmd:**
    - Lines ~95-110: dataset filtering rationale
    - Lines ~200-220: cross-package comparison approach
3.  **TSENAT_appendix_B.Rmd:**
    - Lines ~70-90: pseudocount rationale
    - Lines ~300-320: GAM assumption explanation

**Time to fix:** 30 minutes

------------------------------------------------------------------------

### 8. 🟢 LOW PRIORITY: Add Troubleshooting Section

**Files affected:** Main vignette (TSENAT.Rmd)

#### Suggested addition - New section before Appendices:

``` rmarkdown
# Troubleshooting Common Issues

## Problem: "Object of class TSENATAnalysis not found"
**Cause:** Analysis object not properly created with `build_analysis()`
**Solution:**
1. Verify `config <- TSENAT_config(...)` completed without errors
2. Ensure all input files (readcounts, metadata, tx2gene) are valid
3. Check that `metadata_df` columns match `sample_col` and `condition_col`

## Problem: "Filter removed all genes"
**Cause:** `filter_analysis(..., stringency="medium")` too strict for low-coverage data
**Solution:**
```r
# Use lenient filtering for low-coverage datasets
analysis <- filter_analysis(analysis, stringency = "lenient")
```

## Problem: “Bootstrap confidence intervals are very wide”

**Cause:** Small sample size or high entropy variance **Solution:** 1.
Increase bootstrap replicates: `calculate_diversity(..., nboot=1000)` 2.
Use BCA bootstrap for skewed distributions: `bootstrap_method="bca"` 3.
Check sample quality with
[`calculate_m_estimator()`](https://gallardoalba.github.io/TSENAT/reference/calculate_m_estimator.md)

## Problem: “NAs in LinearModel output”

**Cause:** Model convergence issues (common with sparse entropy data)
**Solution:** 1. Verify sufficient samples per group (≥6-8 recommended)
2. Try alternative method: `calculate_lm(..., method="lmm")` 3. Inspect
model diagnostics: `results(analysis, type="assumptions")`

## Problem: “Plots are hard to read / overlapping labels”

**Cause:** Default sizing for large number of genes **Solution:**

``` r

# Increase figure dimensions
plot_lm(analysis, n_top = 4)  # Reduce to top 4 genes
```

\`\`\`

**Why:** Reduces user support burden, improves experience  
**Time to fix:** 1 hour (requires testing each scenario)

------------------------------------------------------------------------

## Implementation Priority & Timeline

### Phase 1: Critical (Do Now - 10 minutes)

Add session info to main vignette

Complete vignette metadata keywords

### Phase 2: Important (Before submission - 2 hours)

Standardize code output suppression

Add function documentation cross-references

Add figure alt-text

### Phase 3: Polish (Recommended - 3 hours)

Test LaTeX rendering thoroughly

Add inline code comments

Create troubleshooting section

**Total estimated time:** ~5-6 hours for all improvements

------------------------------------------------------------------------

## Bioconductor-Specific Checklist

All `VignetteKeyword{}` tags added to each vignette

Session information included (all vignettes)

No external image files (all figures generated in-code) ✅

Bibliography in bibtex format ✅

Vignette compilation time \<5 minutes (verify during check)

No hard-coded file paths (all use
[`system.file()`](https://rdrr.io/r/base/system.file.html)) ✅

No interactive content requiring manual input (all examples
self-contained) ✅

All functions called are exported/documented (verify with `NAMESPACE`)
✅

Figures have meaningful captions (✅ mostly done, add alt-text)

Code chunks are reproducible (✅ seed set, data included)

------------------------------------------------------------------------

## Quality Assessment

| Aspect | Rating | Notes |
|----|----|----|
| **Mathematical rigor** | ⭐⭐⭐⭐⭐ | Excellent theoretical foundation |
| **Code clarity** | ⭐⭐⭐⭐ | Good, needs minor inline comments |
| **Documentation** | ⭐⭐⭐⭐ | Strong biological motivation, needs API refs |
| **Reproducibility** | ⭐⭐⭐⭐ | Excellent data management, seed control |
| **Accessibility** | ⭐⭐⭐ | Good figures, needs alt-text |
| **Organization** | ⭐⭐⭐⭐⭐ | Excellent structure with appendices |
| **User support** | ⭐⭐⭐ | Good examples, needs troubleshoot section |
| **Bioconductor compliance** | ⭐⭐⭐⭐ | Nearly ready, minor metadata fixes needed |

------------------------------------------------------------------------

## Key Strengths to Maintain

1.  ✅ **Dual validation appendices** (vs SplicingFactory,
    non-parametric methods) - unique and valuable
2.  ✅ **Comprehensive bibliography** (35+ citations) - strong
    foundation reference
3.  ✅ **Progressive workflow** (simple to advanced) - good pedagogical
    structure  
4.  ✅ **Real biological data** - reproducible with vignette data
5.  ✅ **Multiple statistical frameworks** (GAM, LMM, GEE, FPCA, SRH) -
    robust methodology

------------------------------------------------------------------------

## R Markdown Best Practices Applied

Based on official guides in your `/rmarkdown/` folder:

| Best Practice | Current Status | Recommendation |
|----|----|----|
| YAML specification complete | ⚠️ Partial | Complete with keywords |
| Code output suppression deliberate | ⚠️ Inconsistent | Standardize patterns |
| Figure sizing explicit | ✅ Yes | Keep current approach |
| Cross-references present | ⚠️ Limited | Add function [`?name`](https://rdrr.io/r/base/name.html) refs |
| Table formatting consistent | ✅ Yes | Keep `kableExtra` styling |
| Math expressions validated | ⚠️ Minor issues | Test PDF rendering |
| Chunk naming descriptive | ✅ Yes | Good convention |
| Document structure hierarchical | ✅ Yes | Excellent TOC depth |

------------------------------------------------------------------------

## Conclusion

Your vignettes are **publication-ready with minor refinements**. The 8
identified improvements are **achievable in \<6 hours** and will
significantly strengthen your Bioconductor submission.

**Recommended next steps:** 1. Run Phase 1 improvements first (10 min)
2. Test vignette compilation: `R CMD build --no-build-vignettes` then
`R CMD Rd2pdf` 3. Implement Phase 2 before submission 4. Phase 3
improvements are optional but recommended

The package will be stronger for Bioconductor review with these
enhancements.
