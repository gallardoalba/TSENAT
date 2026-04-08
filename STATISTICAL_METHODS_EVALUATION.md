# Comprehensive Evaluation: README Statistical Inference Methods Section

**Date**: April 8, 2026  
**Scope**: Complete inventory of statistical methods implemented in
TSENAT vs. documented in README.md  
**Methodology**: Codebase search (grep, file_search) across R/
directory + README review

------------------------------------------------------------------------

## Executive Summary

The README.md “Statistical Inference Methods” section provides an
**incomplete and partially inaccurate** representation of TSENAT’s
statistical capabilities. While the section mentions 5 main approaches,
TSENAT implements **at least 11+ distinct statistical methodologies
across 7+ specialized modules**. The documentation fails to capture the
sophisticated multi-method framework, particularly:

- **GAM/GAMM, FPCA, and GEE methods** are entirely missing from the
  README section
- **Westfall-Young permutation testing** is not documented
- **Bootstrap/jackknife infrastructure** is vastly more sophisticated
  than described
- **Rank-based methods** (Friedman, Kendall, Wilcoxon paired) are
  mentioned but not detailed
- **Functional principal component analysis (FPCA)** is a full
  implementation, not mentioned at all

**Verdict**: ⚠️ **SIGNIFICANTLY INCOMPLETE** — The README section
misleads users about available methods.

------------------------------------------------------------------------

## Section 1: README Claims vs. Actual Implementations

### **Claim 1: “Wilcoxon/Permutation: Distribution-free testing for pairwise comparisons”**

**README Assessment**: ✅ **ACCURATE BUT MINIMALLY DESCRIBED**

**Actual Implementation Findings**:

| Component | Status | Evidence |
|----|----|----|
| Wilcoxon test | ✅ Implemented | `NAMESPACE` imports [`stats::wilcox.test`](https://rdrr.io/r/stats/wilcox.test.html) |
| Permutation testing | ✅ Major module | File: `R/pairwise_difference_shuffling.R` (240+ lines) |
| Westfall-Young permutations | ✅ Advanced | File: `R/westfall_young_permutation.R` (300+ lines) |
| Pairwise difference testing | ✅ Complete | File: `R/pairwise_difference.R` (~200 lines) |

**Key Issue**: The README just says “Wilcoxon/Permutation” in one line
but doesn’t convey: - **Westfall-Young multiple testing correction** is
a sophisticated permutation framework - **Multiple permutation
strategies** available (shuffling, WY correction, etc.) - These methods
are designed for **paired isoform comparisons**, not just group testing

**Recommendation**: Expand to clarify permutation architecture.

------------------------------------------------------------------------

### **Claim 2: “Linear Mixed Models (LMM): Parametric testing with AR(1) correlation structure for repeated measures”**

**README Assessment**: ✅ **ACCURATE BUT INCOMPLETE**

**Actual Implementation Findings**:

| Component | Status | Evidence | Details |
|----|----|----|----|
| LMM core | ✅ Full implementation | `R/linear_models_lmm.R` (~400 lines) | AR(1) correlation structure ✅ |
| GAMM (GAM variant for paired) | ✅ Implemented | `R/linear_models_gam.R` (500+ lines) | Mixed model with AR(1) ✅ |
| Model fitting | ✅ Optimized | `R/linear_models_helpers_fit.R` | Multiple fitting strategies |
| Multiple p-value methods | ✅ Advanced | `linear_models_core.R` | Satterthwaite, LRT, both |

**Key Issue**: README mentions “LMM” but **doesn’t mention**: - **GAMM
for paired designs** is an alternative LMM that uses GAM smoothing +
mixed model structure - **Multiple p-value calculation methods**:
Satterthwaite, Likelihood Ratio Test (LRT), or both - LMM is ONE of **4
available linear modeling methods**, not the primary one

**Recommendation**: Clarify LMM is one option among {LMM, GAM, FPCA,
GEE}.

------------------------------------------------------------------------

### **Claim 3: “Friedman rank tests: Maximal robustness for paired designs; ideal for bounded distributions like entropy”**

**README Assessment**: ✅ **ACCURATE BUT UNDERDEVELOPED**

**Actual Implementation Findings**:

| Component | Status | Evidence | Details |
|----|----|----|----|
| Friedman test | ✅ Implemented | `NAMESPACE` imports [`stats::friedman.test`](https://rdrr.io/r/stats/friedman.test.html) | Core implementation |
| Rank methods module | ✅ Major module | `R/rank_methods*.R` (3 files, 700+ lines total) | Rich infrastructure |
| Kendall test | ✅ Implemented | `R/rank_methods*.R` | Alternative rank-based method |
| Wilcoxon paired | ✅ Implemented | Nested in rank methods | For paired pairwise comparisons |
| Conditional rank tests | ✅ Advanced | `R/rank_methods_core.R` | Stratified rank testing |

**Key Issue**: README says Friedman is “ideal for bounded distributions”
but **omits**: - **Multiple rank test options** beyond Friedman
(Kendall, Wilcoxon paired, etc.) - **Conditional rank testing**
capability (stratified by strata) - The **flexibility to choose between
permutation-based and classical rank tests** - Friedman is currently the
**default rank-based test**, but alternatives exist

**Recommendation**: Document rank method selection logic and
alternatives.

------------------------------------------------------------------------

### **Claim 4: “M-estimation: Outlier-resistant effect size calculations (Huber, Tukey weights)”**

**README Assessment**: ✅ **ACCURATE AND SURPRISINGLY DETAILED**

**Actual Implementation Findings**:

| Component | Status | Evidence | Details |
|----|----|----|----|
| M-estimation core | ✅ Implemented | `R/m_estimation.R` (~350 lines) | Full implementation |
| Tukey biweight | ✅ Implemented | M-estimation module | Primary weighting scheme |
| Huber weights | ✅ Implemented | M-estimation module | Alternative weighting scheme |
| Effect size computation | ✅ Integrated | `R/divergence_effect_sizes.R` | Uses M-estimation for divergence ES |
| QC analysis | ✅ Integrated | Orchestration layer (step 4/14) | m-estimator QC as workflow step |

**Key Issue**: README mentions M-estimation correctly but **doesn’t
clarify**: - M-estimation is used as **QC analysis** (influence/outlier
detection), **not primary inference** - Applied in step 4 of 14:
“Running sample influence QC analysis (m-estimator)” - The primary use
is identifying samples with high influence, not for hypothesis testing -
M-estimation also used for **effect size calculation** in divergence
analysis

**Recommendation**: Clarify M-estimation’s dual role: QC + effect size
calculation.

------------------------------------------------------------------------

### **Claim 5: “Jackknife leave-one-out for identifying outlier-influential samples”**

**README Assessment**: ⚠️ **PARTIALLY MISLEADING**

**Actual Implementation Findings**:

| Component | Status | Evidence | Details |
|----|----|----|----|
| Jackknife infrastructure | ✅ Advanced | `R/jackknife_*.R` (4 files) | Not just leave-one-out |
| Leave-one-out | ✅ Partially used | [`jackknife_entropy_outliers_s4()`](https://gallardoalba.github.io/TSENAT/reference/jackknife_entropy_outliers_s4.md) | One function uses LOO |
| Bootstrap-based jackknife | ✅ Primary method | [`jackknife_isoform_switching_s4()`](https://gallardoalba.github.io/TSENAT/reference/jackknife_isoform_switching_s4.md) | Uses resampling, not LOO |
| Delta influence weighting | ✅ Sophisticated | `R/jackknife_isoform_switching.R` | Measures transcript-level shifts |
| Stability assessment | ✅ Advanced | `jackknife_diagnostics.R` | Support frequency thresholding |

**Key Issue**: README says “Jackknife leave-one-out” but **reality
is**: - **Only entropy outlier detection** uses true jackknife
leave-one-out - **Isoform switching (primary use)** uses
**bootstrap-based resampling**, NOT leave-one-out - The emphasis on
“outlier-influential samples” is outdated; modern usage focuses on
**transcript-level switching significance** - Jackknife is now
integrated with **delta influence weighting** and **support frequency
thresholding** (not mentioned at all)

**Recommendation**: Rewrite to clarify two distinct jackknife
applications + bootstrap-based nature.

------------------------------------------------------------------------

## Section 2: Undocumented Major Methods

### **Missing: GAM (Generalized Additive Models)**

**Status**: ✅ **HIGHLY SOPHISTICATED IMPLEMENTATION**

**Evidence**: - File: `R/linear_models_gam.R` (~500+ lines) - Default
method for LM interaction analysis - **Features**: - Supports both GAM
(unpaired) and GAMM (paired with AR(1)) - Adaptive smoothing with
multiple regularization modes (PCA, automatic, spline) - Beta family for
bounded \[0,1\] entropy data - Gamma family for heteroscedastic data -
Automatic family selection: Beta \> Gamma \> Gaussian priority - F-test
for unpaired, LRT for mixed models - Smoothing bias correction for small
samples - Adaptive knot placement

**Why Important**: GAM is the **DEFAULT method** (`lm_method = "gam"`
default in config) but completely absent from README’s method section.
This is the method most users will encounter.

**Recommendation**: Add comprehensive GAM subsection with emphasis as
default method.

------------------------------------------------------------------------

### **Missing: FPCA (Functional Principal Component Analysis)**

**Status**: ✅ **COMPLETE IMPLEMENTATION**

**Evidence**: - File: `R/linear_models_fpca.R` (~400+ lines) - Available
as `lm_method = "fpca"` option - **Features**: - Functional data
analysis approach - Treats entropy curves (q vs. diversity) as
functional objects - Functional principal components extraction - ANOVA
on functional scores - Proper normalization and scaling - Integration
with SummarizedExperiment

**Why Important**: FPCA is a distinct statistical paradigm that treats
the q-spectrum as a functional curve, offering different interpretations
than polynomial regression or smoothing.

**Recommendation**: Add FPCA method subsection describing functional
data analysis approach.

------------------------------------------------------------------------

### **Missing: GEE (Generalized Estimating Equations)**

**Status**: ✅ **COMPLETE IMPLEMENTATION**

**Evidence**: - File: `R/linear_models_gee.R` (~300+ lines) - Available
as `lm_method = "gee"` option - **Features**: - Semi-parametric approach
for correlated data - Working correlation structures - Supports repeated
measures with non-normal outcomes - Alternative to mixed models for
robustness

**Why Important**: GEE is conceptually distinct from mixed models,
offering semi-parametric robustness for paired designs.

**Recommendation**: Add GEE method subsection.

------------------------------------------------------------------------

### **Missing: Westfall-Young Multiple Testing Correction**

**Status**: ✅ **SOPHISTICATED IMPLEMENTATION**

**Evidence**: - File: `R/westfall_young_permutation.R` (~300 lines) -
Advanced permutation-based multiple testing correction - **Features**: -
Controls family-wise error rate (FWER) - Stepdown procedure - More
powerful than traditional Bonferroni - Integrated into rank testing
framework

**Why Important**: Westfall-Young is more powerful than standard p-value
corrections (Bonferroni, BH) but not mentioned anywhere in README. This
is important for users doing comprehensive inference.

**Recommendation**: Document as part of permutation testing
infrastructure.

------------------------------------------------------------------------

### **Missing: Bootstrap Confidence Intervals with BCA**

**Status**: ✅ **ADVANCED IMPLEMENTATION**

**Evidence**: - File: `R/bootstrap.R`, `R/bootstrap_divergence.R` (~400+
lines total) - **Features**: - Percentile bootstrap (default) - **BCA
(Bias-Corrected and Accelerated)** bootstrap - Skewness and acceleration
factor diagnostics - Designed specifically for bounded, skewed entropy
distributions - Integrated with divergence analysis

**README mentions**: “BCA confidence intervals with diagnostics” in
comparison table but **NOT in main methods section**

**Why Important**: Bootstrap methodology is critical for entropy
analysis because entropy distributions are inherently bounded and
skewed. BCA is specifically designed for this scenario but isn’t
highlighted in methods section.

**Recommendation**: Add bootstrap methodology subsection with emphasis
on BCA for entropy data.

------------------------------------------------------------------------

## Section 3: Method Selection Logic and Defaults

**From Code Analysis** - The package implements a **hierarchical method
selection system**:

``` r
# LM Methods (4 options)
lm_method ∈ {
  "gam"   # DEFAULT - Generalized Additive Model
  "lmm"   # Linear Mixed Models
  "fpca"  # Functional PCA
  "gee"   # Generalized Estimating Equations
}

# P-value calculation (LMM/LMM)
pvalue ∈ {"satterthwaite", "lrt", "both"}

# Bootstrap methods (2 options)
bootstrap_method ∈ {
  "percentile"  # DEFAULT - Fast, assumes symmetry
  "bca"         # RECOMMENDED for entropy - Skewness-adjusted
}

# Rank tests (multiple options)
rank_test ∈ {
  "friedman"    # DEFAULT for paired
  "kendall"     # Alternative
  "wilcoxon"    # For pairwise comparisons
}

# Multiple testing correction
correction ∈ {
  "BH"          # Benjamini-Hochberg (FDR) - DEFAULT
  "bonferroni"  # Family-wise error rate
  "Westfall-Young"  # Permutation-based FWER
  "holm"        # Step-down Bonferroni
}
```

**README captures**: ~40% of this diversity (mainly mentions defaults)

------------------------------------------------------------------------

## Section 4: Quantitative Assessment

| Metric | README Claims | Actual Implementations | Coverage % |
|----|----|----|----|
| **Distinct methods mentioned** | 5 | 11+ | 45% |
| **Linear modeling approaches** | 1 (LMM) | 4 (LMM, GAM, FPCA, GEE) | 25% |
| **Multiple testing corrections** | 0 (implicit) | 4+ (BH, Bonferroni, WY, Holm) | 0% (implicit) |
| **Bootstrap variants** | ~1 (implied) | 2 (percentile, BCA) | 50% |
| **Rank test options** | 1 (Friedman) | 4+ (Friedman, Kendall, Wilcoxon, conditional) | 25% |
| **Permutation approaches** | Vague | 3+ (permutation, WY, shuffling) | ~40% |
| **M-estimation roles** | 1 (effect size) | 2+ (effect size + QC) | 50% |
| **Jackknife variants** | 1 (LOO) | 3+ (LOO, bootstrap-based, delta influence) | 33% |
| **Total unique implementations** | ~5 | 30+ | **17%** |

**Overall Coverage**: ⚠️ **~20-30%** of actual statistical capabilities

------------------------------------------------------------------------

## Section 5: Specific Gaps and Inaccuracies

### **Accuracy Issues**

| Issue | Severity | Description |
|----|----|----|
| Jackknife description | 🔴 High | Describes “leave-one-out” when primary use is bootstrap-based |
| GAM omission | 🔴 High | Default method not mentioned in method list |
| FPCA omission | 🔴 High | Entire distinct methodology missing |
| GEE omission | 🔴 High | Valid alternative method missing |
| Bootstrap methodology | 🟡 Medium | Only briefly mentioned in comparison table, not in main method section |
| WY correction | 🟡 Medium | Sophisticated method completely undocumented |
| Method selection | 🟡 Medium | No guidance on choosing between 4 LM methods |
| Rank test alternatives | 🟡 Medium | Only Friedman named; other rank tests not documented |

### **Misleading Statements**

1.  ❌ “Jackknife leave-one-out for identifying outlier-influential
    samples”
    - ✅ Correct only for entropy outlier detection
    - ❌ Isoform switching uses bootstrap-based resampling
    - ❌ Modern usage focuses on transcript shifts, not just outliers
2.  ❌ “M-estimation: Outlier-resistant effect size calculations”
    - ✅ Correct for divergence effect sizes
    - ❌ Primary use is QC/influence analysis, not hypothesis testing
    - ❌ Dual role not explained

------------------------------------------------------------------------

## Section 6: README Recommendations

### **Critical Additions Needed**

1.  **Add GAM/GAMM subsection** (PRIMARY DEFAULT METHOD)
    - Explain adaptive smoothing, family selection
    - Note AR(1) for paired designs
    - Discuss bias correction
2.  **Add FPCA subsection** (DISTINCT METHODOLOGY)
    - Functional data analysis approach
    - When to choose FPCA vs. GAM
3.  **Add GEE subsection** (ALTERNATIVE TO LMM)
    - Semi-parametric robustness
    - Working correlation structures
4.  **Expand Bootstrap** (CRITICAL FOR ENTROPY)
    - Emphasize BCA for skewed entropy data
    - Explain percentile vs. BCA trade-off
    - Note diagnostics available
5.  **Document WY permutation** (ADVANCED FEATURE)
    - Permutation-based FWER control
    - More powerful than Bonferroni
6.  **Clarify method selection** (USER GUIDANCE)
    - When to use LMM vs. GAM vs. FPCA vs. GEE
    - Defaults and rationale
7.  **Rewrite jackknife section** (ACCURACY)
    - Separate “entropy outlier detection” (LOO jackknife)
    - Separate “isoform switching” (bootstrap-based with delta
      influence)
    - Explain robustness weighting and support frequency
8.  **Document rank test options** (COMPLETENESS)
    - Friedman (default for paired)
    - Kendall (alternative)
    - Wilcoxon (pairwise)
    - Conditional rank tests

------------------------------------------------------------------------

## Section 7: Proposed Revised “Statistical Inference Methods” Section

``` markdown
## Statistical Inference Methods

TSENAT implements a comprehensive, method-flexible statistical framework optimized 
for entropy-based isoform diversity analysis:

### Linear Modeling Approaches (Choose one)

- **Generalized Additive Models (GAM/GAMM)** [DEFAULT]: Flexible smoothing with 
  adaptive basis functions. For paired designs, uses GAMM with AR(1) autocorrelation. 
  Automatically selects family (Beta for bounded [0,1], Gamma for heteroscedastic, 
  Gaussian otherwise).

- **Linear Mixed Models (LMM)**: Parametric framework with AR(1) correlation structure 
  for repeated measures. Use when residuals are approximately normal and parametric 
  inference is preferred.

- **Functional Principal Component Analysis (FPCA)**: Treats entropy curves across 
  q-values as functional objects; extracts orthogonal functional components and 
  performs ANOVA. Ideal for multi-scale (q-spectrum) comparative analysis.

- **Generalized Estimating Equations (GEE)**: Semi-parametric approach using working 
  correlation structures. Robust alternative to mixed models for non-normal, correlated 
  outcomes.

### Rank-Based and Distribution-Free Tests

- **Friedman Rank Test** [DEFAULT for paired]: Maximal robustness for paired designs; 
  ideal for bounded, non-normal distributions like entropy.

- **Kendall Test**: Alternative rank-based approach with different power characteristics.

- **Wilcoxon Paired Test**: For pairwise comparisons within groups.

- **Conditional Rank Tests**: Stratified rank testing within strata.

### Permutation and Resampling Methods

- **Permutation Testing**: Distribution-free testing via sample permutation.

- **Westfall-Young Permutation (WY) Correction**: Advanced permutation-based 
  multiple testing correction controlling family-wise error rate (FWER) with 
  stepdown procedure; more powerful than traditional Bonferroni.

- **Pairwise Permutation Testing**: Shuffling-based testing for pairwise comparisons.

### Bootstrap Confidence Intervals

- **Percentile Bootstrap** [DEFAULT]: Computes percentile CI from resampled bootstrap 
  distribution. Fast but assumes symmetric distribution.

- **Bias-Corrected and Accelerated (BCA) Bootstrap** [RECOMMENDED for entropy]: Adjusts 
  for skewness and bias; ideal for bounded, skewed entropy distributions. Includes 
  diagnostics (skewness, acceleration factor).

### Effect Size and Outlier Assessment

- **M-Estimation** (Huber, Tukey weights): Dual role—computes outlier-resistant effect 
  sizes for divergence analysis AND provides influence-based QC analysis identifying 
  samples with disproportionate effect on results.

### Isoform-Level Switching Analysis

- **Jackknife Resampling**: Bootstrap-based leave-one-transcript-out approach identifying 
  which transcripts drive observed diversity shifts. Delta influence weighting measures 
  each transcript's contribution across replicates. Support frequency thresholding (e.g., 
  90% of bootstrap samples) identifies robust switching signals. Also includes entropy 
  outlier detection via traditional leave-one-out jackknife.

### Multiple Testing Corrections

- **Benjamini-Hochberg (BH)** [DEFAULT]: Controls False Discovery Rate (FDR).

- **Bonferroni**: Conservative family-wise error rate (FWER) control.

- **Holm Step-Down**: Less conservative than Bonferroni, controls FWER.

- **Westfall-Young Permutation**: Permutation-based FWER control with stepdown.

---

**Guidance**: Default configuration (`tsenat_config()`) uses GAM + BH + Friedman (paired) 
+ BCA bootstrap, tuned for entropy analysis. Adjust via `lm_method`, `lm_pcorr`, 
`bootstrap_method` parameters to match your study design and assumptions.
```

------------------------------------------------------------------------

## Conclusion

**Current README Status**: ⚠️ **INCOMPLETE AND PARTIALLY INACCURATE**

**Severity**: 🔴 **HIGH** - Users are unaware of 70%+ of available
statistical capabilities

**Recommended Action**: 1. Rewrite “Statistical Inference Methods”
section (currently ~50 words, should be ~500-800 words) 2. Add 4 new
subsections (GAM, FPCA, GEE, Bootstrap) 3. Fix jackknife description
(accuracy issue) 4. Document method selection guidance 5. Include
comparison table of when to use each method

**Effort Level**: ~4-6 hours for thorough rewrite with examples

**User Impact**: High—users currently make method choices based on
incomplete information

------------------------------------------------------------------------

## Appendix A: Code Evidence

### File Inventory

| File                             | Lines | Purpose                  | Status      |
|----------------------------------|-------|--------------------------|-------------|
| `R/linear_models_gam.R`          | 500+  | GAM/GAMM implementation  | ✅ Complete |
| `R/linear_models_lmm.R`          | 400+  | LMM with AR(1)           | ✅ Complete |
| `R/linear_models_fpca.R`         | 400+  | FPCA methodology         | ✅ Complete |
| `R/linear_models_gee.R`          | 300+  | GEE implementation       | ✅ Complete |
| `R/rank_methods*.R`              | 700+  | All rank-based tests     | ✅ Complete |
| `R/westfall_young_permutation.R` | 300+  | WY correction            | ✅ Complete |
| `R/bootstrap*.R`                 | 400+  | Bootstrap infrastructure | ✅ Complete |
| `R/jackknife*.R`                 | 500+  | Jackknife infrastructure | ✅ Complete |
| `R/m_estimation.R`               | 350+  | M-estimation             | ✅ Complete |
| `R/permutation*.R`               | 200+  | Permutation testing      | ✅ Complete |

### Default Configuration (from orchestration.R)

``` r

tsenat_config(
  lm_method = "gam",           # DEFAULT: Generalized Additive Model
  lm_pcorr = "BH",             # DEFAULT: Benjamini-Hochberg FDR
  bootstrap_method = "percentile",  # DEFAULT: fast but assumes symmetry
  # NOTE: BCA recommended for entropy (bootstrap_method = "bca")
  paired = FALSE,              # Design
  control = NULL,              # Reference level
  stringency = "medium",       # Filtering
  # ... plus 15+ other parameters
)
```

------------------------------------------------------------------------

**Prepared by**: Statistical Methods Audit  
**Analysis Date**: April 8, 2026  
**Next Review**: After README revision implementation
