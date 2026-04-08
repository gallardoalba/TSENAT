# TSENAT Core Functions: Missing Parameters Audit Report

## Structured Findings & Recommendations

**Date:** April 8, 2026  
**Scope:** calculate_lm_interaction_s4(), calculate_divergence_s4(),
effect_sizes_divergence_s4(), jackknife_isoform_switching_s4()  
**Audit Method:** Dynamic analysis of function signatures, orchestration
layer, config resolution, and internal computation functions

------------------------------------------------------------------------

## EXECUTIVE SUMMARY

### Key Finding

**7 critical statistical parameters are not exposed in
tsenat_config()**, making reproducible analysis difficult and preventing
users from documenting or controlling key methodological choices.

### Recommendation Summary

- **4 CRITICAL parameters MUST be added immediately** (reproducibility
  blocking)
- **3 MEDIUM parameters SHOULD be added** (improves user control)
- **5+ OPTIONAL parameters CAN be added** (advanced/specialized use)
- **All existing parameters verified in config** (no removals needed)

### Implementation Effort: ~3-4 hours

------------------------------------------------------------------------

## 1. calculate_lm_interaction_s4() ANALYSIS

### Function Signature (S4 Wrapper)

``` r

calculate_lm_interaction_s4(analysis, fdr_threshold = NULL, formula = NULL,
    condition_col = NULL, method = "gam", paired = NULL, subject_col = NULL, 
    nthreads = NULL, multicorr = NULL, corstr = NULL, pcorr = NULL, 
    verbose = NULL, return_model_data = NULL, output_file = NULL, ...)
```

### Internal Computation Function

``` r

.calculate_lm_interaction(se, condition_col = "condition", min_obs = 5,
    method = c("lmm", "gam", "fpca", "gee"), pvalue = c("satterthwaite", "lrt", "both"),
    subject_col = NULL, paired = FALSE, nthreads = 1, assay_name = "diversity", 
    pcorr = "BH", verbose = FALSE, bias_correction = TRUE, 
    regularization = c("pca", "lasso", "elasticnet", "gamsel", "spline"), 
    corstr = c("ar1", "exchangeable", "independence"), 
    multicorr = c("hochberg", "westfall-young", "benjamini-yekutieli"), 
    storey = FALSE, wy_randomizations = 1000, adaptive_knots = TRUE, 
    return_model_data = FALSE)
```

### RECOMMENDED Parameter Additions

#### ✅ CRITICAL: `method` - LM/GAM Model Selection

| Aspect | Detail |
|----|----|
| **Parameter Name** | method |
| **Current Status** | Hardcoded to “gam” in function signature |
| **Possible Values** | “lmm” (linear mixed), “gam” (additive), “fpca” (functional PCA), “gee” (estimating equations) |
| **What It Controls** | Core statistical approach for modeling q\*condition interactions |
| **Impact on Results** | CRITICAL |
| **Different Values Produce** | Different p-values, effect sizes, and significance calls |
| **Reproducibility** | ⚠️ MUST be documented in analysis reports |
| **Configuration Propagation** | User arg → config slot → .extract_lm_params() → .calculate_lm_interaction() |
| **Recommendation** | ✅ **ADD to tsenat_config() as `lm_method = "gam"`** |
| **Justification** | GAM is smooth-fitting friendly for entropy curves; but users must be able to choose LMM for linear modeling scenarios or GEE for correlated observations |

#### ✅ CRITICAL: `pcorr` - Multiple Comparison Correction

| Aspect | Detail |
|----|----|
| **Parameter Name** | pcorr |
| **Current Status** | Hardcoded to “BH” (Benjamini-Hochberg) in .calculate_lm_interaction() |
| **Possible Values** | “BH”, “bonferroni”, “hochberg”, “holm” |
| **What It Controls** | P-value multiple comparison correction method |
| **Impact on Results** | CRITICAL - Controls FDR/FWER of interaction tests |
| **Different Values Produce** | Different adjusted p-values and gene selection |
| **Reproducibility** | ⚠️ MUST match between analyses |
| **Configuration Propagation** | User arg → config slot → .extract_lm_params() → .calculate_lm_interaction() |
| **Recommendation** | ✅ **ADD to tsenat_config() as `lm_pcorr = "BH"`** |
| **Justification** | BH controls FDR (less conservative); bonferroni controls FWER (more conservative). Users need explicit control based on study requirements |

#### ⭕ OPTIONAL: `multicorr` - Alternative Multi-Test Correction

- **Status**: Currently NULL (unused)
- **Values**: “hochberg”, “westfall-young”, “benjamini-yekutieli”
- **Use Case**: When stricter correction than pcorr needed
- **Recommendation**: ⭕ Add to config for advanced users (Optional)

#### ⭕ OPTIONAL: `corstr` - GEE Correlation Structure

- **Status**: Currently NULL
- **Values**: “ar1”, “exchangeable”, “independence”
- **When Relevant**: Only when method=“gee”
- **Recommendation**: ⭕ Add to config as conditional parameter
  (Optional)

### Parameters Already Correctly Exposed

✅ fdr_threshold  
✅ condition_col (auto-detected, good default)  
✅ paired  
✅ subject_col  
✅ nthreads  
✅ verbose  
✅ return_model_data  
✅ output_file

### Parameters to Keep Internal Only

- `min_obs` - validation threshold, not configurable
- `formula` - advanced; pass directly to function if needed
- `pvalue` - internal LM calculation method
- `assay_name` - auto-detected from se
- `bias_correction` - keep TRUE
- `regularization` - internal tuning
- `storey` - advanced; rarely used
- `adaptive_knots` - GAM tuning; internal

------------------------------------------------------------------------

## 2. calculate_divergence_s4() ANALYSIS

### Function Signature (S4 Wrapper)

``` r

calculate_divergence_s4(analysis, q = NULL, verbose = FALSE, nthreads = NULL,
    output_file = NULL, control_group = NULL, paired = FALSE, method = NULL, 
    bootstrap = FALSE, nboot = NULL, progress = FALSE, ...)
```

### Internal Computation Function

``` r

.calculate_divergence(se, group_col = NULL, control_group = NULL, q = 1,
    paired = FALSE, bootstrap = FALSE, nboot = "auto", ci = 0.95, 
    method = "percentile", norm = TRUE, log_base = exp(1), 
    pseudocount = 0.5, nthreads = 1, progress = FALSE, verbose = TRUE)
```

### RECOMMENDED Parameter Additions

#### ✅ CRITICAL: `ci` - Confidence Interval Level for Divergence

| Aspect | Detail |
|----|----|
| **Parameter Name** | ci (in internal function) |
| **Current Status** | Hardcoded to 0.95 in .calculate_divergence() |
| **Possible Values** | Any value between 0 and 1 (commonly 0.90, 0.95, 0.99) |
| **What It Controls** | Width of confidence intervals for divergence metrics |
| **Examples** | ci=0.90 → narrower CI; ci=0.99 → wider CI |
| **Impact on Results** | Affects CI bounds only (not point estimates or p-values) |
| **Reproducibility** | ⚠️ Should match between analyses |
| **Configuration Propagation** | User arg → config → .resolve_divergence_parameters() → .calculate_divergence() |
| **Recommendation** | ✅ **ADD to tsenat_config() as `divergence_ci = 0.95`** |
| **Justification** | Currently hardcoded; users cannot adjust without modifying code. Some analyses require different CI levels |

#### ✅ MEDIUM: `control_group` - Reference Group for Divergence Comparisons

| Aspect | Detail |
|----|----|
| **Parameter Name** | control_group |
| **Current Status** | Can be passed to function; may not be in config extraction |
| **Possible Values** | Group name (character) or NULL for all comparisons |
| **What It Controls** | Which group serves as baseline for pairwise divergence comparisons |
| **Impact on Results** | HIGH - Changes comparison structure |
| **Reproducibility** | ⚠️ Must be documented |
| **Related Parameter** | `control` in config (paired reference; different purpose) |
| **Configuration Propagation** | User arg → config → .resolve_divergence_parameters() → .build_divergence_args() |
| **Recommendation** | ✅ **ADD to tsenat_config() as `divergence_control_group = NULL`** |
| **Justification** | Ensure comparisons are explicit and documented |

#### ⭕ OPTIONAL: `progress` - Progress Bar Display

- **Status**: Currently FALSE
- **Values**: TRUE/FALSE
- **Recommendation**: ⭕ Could add (Optional; mostly for user
  convenience)

#### ⚠️ CLARIFY: `method` Parameter

- **Current**: accepts “percentile” or “bca” as bootstrap method
- **Connected to**: `bootstrap_method` in config
- **Note**: Naming could be clearer; consider documentation

### Parameters Already Correctly Exposed

✅ q (mapped from q_values)  
✅ verbose (from config)  
✅ nthreads (from config)  
✅ output_file (orchestration logic)  
✅ paired (from config)  
✅ bootstrap (from config)  
✅ nboot (from config)  
✅ norm (from config)  
✅ log_base (from config)  
✅ pseudocount (from config)

### Parameters to Keep Internal Only

- `group_col` - auto-detected from condition_col
- `assay_name` - internal; auto-determined

### ⚠️ CRITICAL FOR ENTROPY: Bootstrap Method

**When analyzing bounded distributions (like Tsallis entropy):**

``` r

# ✅ CORRECT for entropy (0 to log N):
config <- tsenat_config(bootstrap_method = "bca")
# BCA = Bias-Corrected and Accelerated
# Accounts for skewness inherent in bounded distributions
# Produces valid CIs that contain point estimate

# ❌ INCORRECT for entropy:
config <- tsenat_config(bootstrap_method = "percentile")
# Percentile assumes symmetric distribution
# For skewed data (entropy), CI may NOT contain point estimate
# This is mathematically valid but uninformative
```

**Evidence**: Database papers C016 (2005), C030 (2023), S115 (2015) all
demonstrate this pattern

------------------------------------------------------------------------

## 3. effect_sizes_divergence_s4() ANALYSIS

### Function Signature

``` r

effect_sizes_divergence_s4(analysis, significance_threshold = NULL, 
    enrich_per_q_pattern = NULL, verbose = NULL, output_file = NULL, ...)
```

### FINDING: All Critical Parameters Already in Config ✅

| Parameter | Type | Default | Status | Notes |
|----|----|----|----|----|
| significance_threshold | numeric | 0.05 | ✅ In config | Filters genes to keep only adj_p \< threshold |
| enrich_per_q_pattern | logical | TRUE | ✅ In config | Adds pattern enrichment columns to output |
| verbose | logical | FALSE | ✅ In config | Logging verbosity |
| output_file | character | (orchestration) | ✅ Handled | File output logic |

### Conclusion

✅ **NO NEW PARAMETERS NEEDED** - This function is properly configured

### Internal Function for Reference

``` r

.effect_sizes_divergence(lm_res, divergence_results_se, significance_threshold = 0.05,
    enrich_per_q_pattern = TRUE, verbose = FALSE)
```

------------------------------------------------------------------------

## 4. jackknife_isoform_switching_s4() ANALYSIS

### Function Signature

``` r

jackknife_isoform_switching_s4(analysis, condition_col = NULL, subject_col = NULL,
    gene_col = NULL, isoform_col = NULL, q = c(0, 0.5, 1, 1.5, 2), 
    norm = NULL, log_base = NULL, threshold = 90, n_bootstrap = 1000, 
    pseudocount = NULL, lm_results = NULL, lm_p_threshold = 0.05,
    use_lm_fdr = TRUE, output_file = NULL, verbose = FALSE, ...)
```

### Parameter Resolution Function

``` r

.resolve_and_validate_jis_params(q, norm = NULL, log_base = NULL, pseudocount = NULL, 
    n_bootstrap = 1000, threshold = NULL, lm_p_threshold = NULL, analysis, verbose = FALSE)
```

### RECOMMENDED Parameter Additions

#### ✅ CRITICAL: `use_lm_fdr` - Whether to Use FDR-Corrected P-Values

| Aspect | Detail |
|----|----|
| **Parameter Name** | use_lm_fdr |
| **Current Status** | Hardcoded to TRUE in function signature |
| **Possible Values** | TRUE (use FDR-corrected p-values), FALSE (use raw p-values) |
| **What It Controls** | Whether JIS filters genes using adj_p_interaction or p_interaction from LM results |
| **Impact on Results** | CRITICAL - Dramatically affects number of genes selected for switching detection |
| **Example Impact** | TRUE → stricter (maybe 100 genes); FALSE → lenient (maybe 500 genes) |
| **Reproducibility** | ⚠️ MUST be documented; must match between analyses |
| **Configuration Propagation** | Currently: NOT checked from config! Just uses default TRUE |
| **Recommendation** | ✅ **ADD to tsenat_config() as `jis_use_lm_fdr = TRUE`** |
| **Justification** | Users should be able to choose stringency level. Current hardcoding prevents flexibility |

### Parameters Already Correctly Exposed

✅ condition_col (auto-detect capable)  
✅ subject_col (from config)  
✅ q (from config via q_values)  
✅ norm (from config)  
✅ log_base (from config)  
✅ threshold (from config; default 90)  
✅ n_bootstrap (from config; default 1000)  
✅ pseudocount (from config)  
✅ lm_p_threshold (from config; default 0.05)  
✅ output_file (orchestration)

### Parameters to Keep Internal Only

- `gene_col` - auto-detected from rowData
- `isoform_col` - auto-detected from rowData  
- `lm_results` - extracted from analysis object automatically

------------------------------------------------------------------------

## MASTER RECOMMENDATION TABLE

### CRITICAL ADDITIONS (Do Immediately)

| Function | Parameter | Type | Current State | Default | Add To Config As | Priority |
|----|----|----|----|----|----|----|
| calculate_lm_interaction_s4 | method | char | Hardcoded “gam” | “gam” | lm_method | 🔴 HIGH |
| calculate_lm_interaction_s4 | pcorr | char | Hardcoded “BH” | “BH” | lm_pcorr | 🔴 HIGH |
| calculate_divergence_s4 | ci | numeric | Hardcoded 0.95 | 0.95 | divergence_ci | 🔴 HIGH |
| jackknife_isoform_switching_s4 | use_lm_fdr | logical | Hardcoded TRUE | TRUE | jis_use_lm_fdr | 🔴 HIGH |

### MEDIUM PRIORITY ADDITIONS

| Function | Parameter | Type | Current State | Default | Add To Config As |
|----|----|----|----|----|----|
| calculate_divergence_s4 | control_group | char | Can pass; unclear config | NULL | divergence_control_group |

### OPTIONAL ADDITIONS (Nice to Have)

| Function | Parameter | Type | Reason | Config As |
|----|----|----|----|----|
| calculate_lm_interaction_s4 | multicorr | char | Advanced multi-test correction | lm_multicorr |
| calculate_lm_interaction_s4 | corstr | char | GEE correlation structure | lm_corstr |
| calculate_lm_interaction_s4 | storey | logical | Storey q-value method | lm_storey |
| calculate_lm_interaction_s4 | wy_randomizations | numeric | Westfall-Young permutations | lm_wy_randomizations |
| calculate_divergence_s4 | progress | logical | Progress bar display | divergence_progress |

------------------------------------------------------------------------

## PARAMETER IMPACT MATRIX

### Which Global Choices Affect Analysis Results?

[TABLE]

------------------------------------------------------------------------

## CODE LOCATIONS

### Function Definitions

| Function | File | Line | Signature Length |
|----|----|----|----|
| calculate_lm_interaction_s4 | R/s4_functions_lm.R | 142 | Complex; ~12 params |
| .calculate_lm_interaction | R/linear_models_core.R | 275 | Extended; ~12 params + 8 internals |
| calculate_divergence_s4 | R/s4_functions_divergence.R | 107 | ~11 params |
| .calculate_divergence | R/divergence_core.R | 437 | ~13 params |
| .calculate_divergence_impl | R/divergence_core.R | 510+ | Internal implementation |
| jackknife_isoform_switching_s4 | R/s4_functions_jis.R | 173 | ~14 params |
| .resolve_and_validate_jis_params | R/s4_functions_jis.R | 392 | Parameter validation |
| effect_sizes_divergence_s4 | R/s4_functions_effect_size.R | 283 | ~4 params (all configured) |
| tsenat_config | R/orchestration.R | 485 | Needs expansion |
| .execute_lm_interaction_s4 | R/orchestration.R | 725 | Needs enhancement |
| .execute_divergence_s4 | R/orchestration.R | 842 | Needs enhancement |
| .execute_jackknife_isoform_switching | R/orchestration.R | 761 | Needs enhancement |

------------------------------------------------------------------------

## IMPLEMENTATION REQUIREMENTS

### Step 1: Update tsenat_config() Signature

**File:** R/orchestration.R, Line 485  
**Changes:** - Add 4 new parameters to function signature - Add
validation logic for each - Add to config list assembly

**Estimated code:** 30-40 lines

### Step 2: Update Orchestration Layers

**Files:** R/orchestration.R, .execute\_\* functions  
**Changes:** - Extract new parameters from config - Pass to S4 wrapper
functions

**Estimated code:** 15-20 lines

### Step 3: Verify Parameter Resolution

**Files:** R/s4_functions_jis.R, R/s4_functions_divergence.R  
**Changes:** - Ensure resolve_slot_param() handles new config keys -
Verify ci is returned from parameter resolution

**Estimated code:** 5-10 lines

### Step 4: Documentation

**File:** R/orchestration.R (roxygen2 comments)  
**Changes:** - Add @param sections for each new parameter - Add details
on bootstrap method for entropy - Add examples showing method selection

**Estimated code:** 40-50 lines

### Step 5: Testing

- Unit tests: 50-75 lines
- Integration tests: 100-150 lines
- Vignette updates

**Total estimated effort:** 250-300 lines of code + documentation

------------------------------------------------------------------------

## SUGGESTED tsenat_config() UPDATES

### New Function Signature (excerpt)

``` r
tsenat_config <- function(
    # ... existing parameters ...
    
    # === CRITICAL NEW PARAMETERS ===
    lm_method = "gam",                    # "lmm", "gam", "fpca", "gee"
    lm_pcorr = "BH",                      # "BH", "bonferroni", "hochberg", "holm"
    divergence_ci = 0.95,                 # 0 < ci < 1
    jis_use_lm_fdr = TRUE,               # TRUE or FALSE
    
    # === OPTIONAL PARAMETERS ===
    divergence_control_group = NULL,      # character or NULL
    
    ...)
```

### New Validation Code (excerpt)

``` r

    # Validate LM method
    if (!lm_method %in% c("lmm", "gam", "fpca", "gee")) {
        stop("'lm_method' must be one of: lmm, gam, fpca, gee", call. = FALSE)
    }
    
    # Validate LM p-correction
    if (!lm_pcorr %in% c("BH", "bonferroni", "hochberg", "holm")) {
        stop("'lm_pcorr' must be one of: BH, bonferroni, hochberg, holm", call. = FALSE)
    }
    
    # Validate divergence CI
    if (!is.numeric(divergence_ci) || divergence_ci <= 0 || divergence_ci >= 1) {
        stop("'divergence_ci' must be numeric between 0 and 1", call. = FALSE)
    }
    
    # Validate JIS use_lm_fdr
    if (!is.logical(jis_use_lm_fdr) || is.na(jis_use_lm_fdr)) {
        stop("'jis_use_lm_fdr' must be TRUE or FALSE", call. = FALSE)
    }
```

------------------------------------------------------------------------

## RISK ASSESSMENT

### If Changes Are NOT Made

🔴 **HIGH RISK** - Users cannot document their methodological choices -
Reproducibility is compromised - Different researchers may use different
methods unknowingly - Results cannot be compared across analyses

### If Changes ARE Made

✅ **LOW RISK** - Backward compatible (all new parameters have
defaults) - Existing analyses continue to work (defaults match current
behavior) - Users gain explicit control over methodology -
Reproducibility is enhanced

------------------------------------------------------------------------

## FINAL RECOMMENDATIONS

### Priority 1 (This Week)

✅ Add: lm_method, lm_pcorr, divergence_ci, jis_use_lm_fdr  
✅ Add validation for each  
✅ Update orchestration layer (.execute\_\* functions)  
✅ Run integration tests

### Priority 2 (Next Week)

✅ Add: divergence_control_group (medium priority)  
✅ Comprehensive documentation  
✅ Update examples and vignette  
✅ Create reproducibility checklist for users

### Priority 3 (Optional)

⭕ Add: lm_multicorr, lm_corstr, divergence_progress  
⭕ Advanced parameter validation  
⭕ Configuration comparison/diff tool

------------------------------------------------------------------------

## REPRODUCIBILITY CHECKLIST FOR USERS

``` r

# When documenting analysis for reproducibility, specify:
config <- tsenat_config(
  # Core data mapping
  q_values = seq(0, 2, by = 0.05),        # ≥5 required; 41+ for paired
  condition_col = "condition",
  subject_col = "paired_samples",         # If paired=TRUE
  paired = FALSE,                         # Or TRUE
  
  # Statistical methods (NEW - CRITICAL)
  lm_method = "gam",                      # ← NEW: Document which method
  lm_pcorr = "BH",                        # ← NEW: Document correction
  
  # Bootstrap method (CRITICAL for entropy)
  bootstrap_method = "bca",               # ← CRITICAL: Use "bca" for entropy!
  
  # Thresholds and parameters
  fdr_threshold = 0.05,
  significance_threshold = 0.05,
  divergence_ci = 0.95,                   # ← NEW: Document CI level
  jis_use_lm_fdr = TRUE,                 # ← NEW: Document filter choice
  
  # Other parameters
  nthreads = 1                            # Document for reproducibility
)
```

------------------------------------------------------------------------

## FILES DELIVERED

This audit includes 4 comprehensive documents in the TSENAT project
directory:

1.  **UNDOCUMENTED_PARAMETERS_AUDIT.md** - Complete technical audit
    (~25KB)
2.  **IMPLEMENTATION_SUMMARY.md** - Executive summary & quick reference
    (~8KB)
3.  **IMPLEMENTATION_CODE_GUIDE.md** - Step-by-step implementation guide
    (~12KB)
4.  **PARAMETER_REFERENCE_TABLES.md** - Detailed parameter tables
    (~15KB)

------------------------------------------------------------------------

**Audit Report Complete**  
**Date:** April 8, 2026  
**Status:** ✅ Ready for Implementation
