# Parameter Reference Tables

## Complete listing of parameters by function

------------------------------------------------------------------------

## TABLE 1: calculate_lm_interaction_s4()

### S4 Wrapper Function Signature

``` r

calculate_lm_interaction_s4(analysis, fdr_threshold = NULL, formula = NULL,
    condition_col = NULL, method = "gam", paired = NULL, subject_col = NULL, 
    nthreads = NULL, multicorr = NULL, corstr = NULL, pcorr = NULL, 
    verbose = NULL, return_model_data = NULL, output_file = NULL, ...)
```

| Parameter | Type | Current Source | Config Available | Recommendation |
|----|----|----|----|----|
| fdr_threshold | numeric | config (fdr_threshold) | ✅ | Keep in config |
| formula | character | direct call | ❌ | INTERNAL ONLY |
| condition_col | character | auto-detect / config | ✅ | Keep auto-detect |
| method | character | **hardcoded “gam”** | ❌ | ✅ **ADD to config** as `lm_method` |
| paired | logical | config(paired) | ✅ | Keep in config |
| subject_col | character | config(subject_col) | ✅ | Keep in config |
| nthreads | numeric | config(nthreads) | ✅ | Keep in config |
| multicorr | character | hardcoded NULL | ❌ | ⭕ OPTIONAL add to config |
| corstr | character | hardcoded NULL | ❌ | ⭕ OPTIONAL add to config |
| pcorr | character | **hardcoded “BH”** | ❌ | ✅ **ADD to config** as `lm_pcorr` |
| verbose | logical | config(verbose) or NULL | ✅ | Keep in config |
| return_model_data | logical | config(return_model_data) | ✅ | Keep as-is |
| output_file | character | config(output_file) | ✅ | Keep in config |

### Internal Function Signature (.calculate_lm_interaction)

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

| Parameter | Internal Only | Values | Default | Notes |
|----|:--:|----|----|----|
| se | ✅ | SummarizedExperiment | \- | Input data |
| min_obs | ✅ | numeric | 5 | Validation threshold |
| method | ❌ | “lmm”, “gam”, “fpca”, “gee” | “gam” | **EXPOSE as config** |
| pvalue | ✅ | “satterthwaite”, “lrt”, “both” | \- | P-value calculation method |
| subject_col | ❌ | character or NULL | NULL | From config |
| paired | ❌ | logical | FALSE | From config |
| nthreads | ❌ | numeric | 1 | From config |
| assay_name | ✅ | character | “diversity” | Internal; validation only |
| pcorr | ❌ | “BH”, “bonferroni”, “hochberg”, “holm” | “BH” | **EXPOSE as config** |
| verbose | ❌ | logical | FALSE | From config |
| bias_correction | ✅ | logical | TRUE | INTERNAL; keep as-is |
| regularization | ✅ | “pca”, “lasso”, “elasticnet”, “gamsel”, “spline” | \- | INTERNAL; advanced tuning |
| corstr | ✅ | “ar1”, “exchangeable”, “independence” | \- | INTERNAL; GEE-specific |
| multicorr | ✅ | “hochberg”, “westfall-young”, “benjamini-yekutieli” | \- | INTERNAL; alternative to pcorr |
| storey | ✅ | logical | FALSE | INTERNAL; advanced |
| wy_randomizations | ✅ | numeric | 1000 | INTERNAL; Westfall-Young specific |
| adaptive_knots | ✅ | logical | TRUE | INTERNAL; GAM-specific |

### Configuration Propagation Path

    user calls: calculate_lm_interaction_s4(analysis, method="lmm", pcorr="bonferroni")
        ↓
    config slot: analysis@config$lm_method = "lmm"
                analysis@config$lm_pcorr = "bonferroni"
        ↓
    .extract_lm_params(): resolves from config
        ↓
    .calculate_lm_interaction(): receives method, pcorr as arguments
        ↓
    LM test executed with specified method and correction

------------------------------------------------------------------------

## TABLE 2: calculate_divergence_s4()

### S4 Wrapper Function Signature

``` r

calculate_divergence_s4(analysis, q = NULL, verbose = FALSE, nthreads = NULL,
    output_file = NULL, control_group = NULL, paired = FALSE, method = NULL, 
    bootstrap = FALSE, nboot = NULL, progress = FALSE, ...)
```

| Parameter | Type | Current Source | Config Available | Recommendation |
|----|----|----|----|----|
| q | numeric | config(q_values) | ✅ | Keep in config |
| verbose | logical | hardcoded FALSE | ❌ | Config has verbose slot |
| nthreads | numeric | config(nthreads) | ✅ | Keep in config |
| output_file | character | orchestration logic | Partial | Keep as-is |
| control_group | character | **hardcoded NULL** | ❌ | ✅ **ADD to config** as `divergence_control_group` |
| paired | logical | config(paired) | ✅ | Keep in config |
| method | character | config(bootstrap_method) | ✅ | Clarify: use bootstrap_method for percentile/bca |
| bootstrap | logical | config(bootstrap) | ✅ | Keep in config |
| nboot | numeric | config(nboot) | ✅ | Keep in config |
| progress | logical | **hardcoded FALSE** | ❌ | ✅ **ADD to config** as `divergence_progress` (OPTIONAL) |

### Internal Function Signature (.calculate_divergence)

``` r

.calculate_divergence(se, group_col = NULL, control_group = NULL, q = 1,
    paired = FALSE, bootstrap = FALSE, nboot = "auto", ci = 0.95, 
    method = "percentile", norm = TRUE, log_base = exp(1), 
    pseudocount = 0.5, nthreads = 1, progress = FALSE, verbose = TRUE)
```

| Parameter | Internal Only | Values | Default | Notes |
|----|:--:|----|----|----|
| se | ✅ | SummarizedExperiment | \- | Input data |
| group_col | ✅ | character or NULL | NULL | Auto-detected |
| control_group | ❌ | character or NULL | NULL | **EXPOSE as config** |
| q | ❌ | numeric vector | 1 | From config(q_values) |
| paired | ❌ | logical | FALSE | From config |
| bootstrap | ❌ | logical | FALSE | From config |
| nboot | ❌ | numeric or “auto” | “auto” | From config |
| ci | ❌ | numeric (0-1) | 0.95 | **EXPOSE as config** as `divergence_ci` |
| method | ❌ | “percentile”, “bca” | “percentile” | From config(bootstrap_method) |
| norm | ❌ | logical | TRUE | From config(norm) |
| log_base | ❌ | numeric | exp(1) | From config(log_base) |
| pseudocount | ❌ | numeric | 0.5 | From config(pseudocount) |
| nthreads | ❌ | numeric | 1 | From config(nthreads) |
| progress | ❌ | logical | FALSE | From config or hardcoded |
| verbose | ❌ | logical | TRUE | From config(verbose) |

### Critical: Bootstrap Method Recommendation

    For BOUNDED DISTRIBUTIONS (like entropy: 0 to log N):
      config <- tsenat_config(bootstrap_method = "bca")
      # NOT "percentile" - percentile CI may not contain point estimate

    For UNBOUNDED DISTRIBUTIONS (like counts, continuous):
      config <- tsenat_config(bootstrap_method = "percentile")
      # Faster and sufficient for symmetric distributions

------------------------------------------------------------------------

## TABLE 3: jackknife_isoform_switching_s4()

### Function Signature

``` r

jackknife_isoform_switching_s4(analysis, condition_col = NULL, subject_col = NULL,
    gene_col = NULL, isoform_col = NULL, q = c(0, 0.5, 1, 1.5, 2), 
    norm = NULL, log_base = NULL, threshold = 90, n_bootstrap = 1000, 
    pseudocount = NULL, lm_results = NULL, lm_p_threshold = 0.05,
    use_lm_fdr = TRUE, output_file = NULL, verbose = FALSE, ...)
```

| Parameter | Type | Current Source | Config Available | Recommendation |
|----|----|----|----|----|
| condition_col | character | auto-detect / config | ✅ | Keep in config |
| subject_col | character | config(subject_col) | ✅ | Keep in config |
| gene_col | character | auto-detect | ❌ | Auto-detect is sufficient |
| isoform_col | character | auto-detect | ❌ | Auto-detect is sufficient |
| q | numeric | function default; can use config | ✅ | Can use config(q_values) |
| norm | logical | config(norm) | ✅ | Keep in config |
| log_base | numeric | config(log_base) | ✅ | Keep in config |
| threshold | numeric | config(threshold) | ✅ | Keep in config |
| n_bootstrap | numeric | config(n_bootstrap) | ✅ | Keep in config |
| pseudocount | numeric | config(pseudocount) | ✅ | Keep in config |
| lm_results | data.frame | extracted from analysis | N/A | Internal; extracted automatically |
| lm_p_threshold | numeric | config(lm_p_threshold) | ✅ | Keep in config |
| use_lm_fdr | logical | **hardcoded TRUE** | ❌ | ✅ **ADD to config** as `jis_use_lm_fdr` |
| output_file | character | orchestration logic | Partial | Keep as-is |
| verbose | logical | hardcoded FALSE | ❌ | Config has verbose slot |

### Parameter Resolution Priority

    .resolve_and_validate_jis_params():
      q           → config$q_values or function default
      norm        → config$norm or TRUE
      log_base    → config$log_base or exp(1)
      pseudocount → config$pseudocount or 0
      n_bootstrap → config$n_bootstrap or 1000
      threshold   → config$threshold or 90
      lm_p_threshold → config$lm_p_threshold or 0.05
      use_lm_fdr  → (currently not checked!) should check config$jis_use_lm_fdr

------------------------------------------------------------------------

## TABLE 4: effect_sizes_divergence_s4()

### Function Signature

``` r

effect_sizes_divergence_s4(analysis, significance_threshold = NULL, 
    enrich_per_q_pattern = NULL, verbose = NULL, output_file = NULL, ...)
```

| Parameter | Type | Current Source | Config Available | Recommendation |
|----|----|----|----|----|
| significance_threshold | numeric | config(significance_threshold) | ✅ | Keep in config |
| enrich_per_q_pattern | logical | config(enrich_per_q_pattern) | ✅ | Keep in config |
| verbose | logical | config(verbose) or NULL | ✅ | Keep in config |
| output_file | character | orchestration logic | Partial | Keep as-is |

**Status:** ✅ All critical parameters already in config!

------------------------------------------------------------------------

## TABLE 5: Current tsenat_config() Parameters

### Already Exposed (No Changes Needed)

| Parameter | Type | Default | Controls | Category |
|----|----|----|----|----|
| q_values | numeric | seq(0,2,0.5) | Tsallis q-spectrum | Analysis |
| condition_col | character | “condition” | Treatment/condition variable | Metadata |
| subject_col | character | NULL | Sample/subject ID for paired | Metadata |
| sample_col | character | “sample” | Sample identifier | Metadata |
| paired | logical | FALSE | Paired design flag | Design |
| control | character | NULL | Reference group (paired) | Design |
| p_threshold | numeric | 0.05 | P-value threshold | Inference |
| fdr_threshold | numeric | 0.05 | FDR threshold for LM | Inference |
| significance_threshold | numeric | 0.05 | Effect size filter | Inference |
| bootstrap | logical | FALSE | Enable bootstrap CIs | Inference |
| nboot | numeric | 1000 | Bootstrap resamples | Inference |
| bootstrap_method | string | “percentile” | “percentile” or “bca” | Inference |
| bootstrap_ci | numeric | 0.95 | CI level (e.g., 0.95) | Inference |
| bootstrap_include_diagnostics | logical | TRUE | Include diagnostic info | Inference |
| stringency | string | “medium” | Filter stringency level | Analysis |
| nthreads | numeric | 1 | Parallel threads | Computation |
| norm | logical | TRUE | Normalization flag | Data |
| norm_method | character | NULL | Normalization method | Data |
| pseudocount | numeric | NULL | Small value for log(0) | Data |
| shrinkage | string | “none” | Shrinkage method | Data |
| min_valid_frac | numeric | 0.75 | Min valid fraction | Filtering |
| n_bootstrap | numeric | ? | JIS bootstrap resamples | Inference |
| threshold | numeric | 90 | JIS switching threshold | Inference |
| lm_p_threshold | numeric | 0.05 | JIS LM p-value filter | Inference |
| log_base | numeric | exp(1) | Log base for JIS | Data |
| enrich_per_q_pattern | logical | TRUE | Effect size enrichment | Analysis |

------------------------------------------------------------------------

## TABLE 6: New Parameters to Add

### CRITICAL (Must Add)

| Config Key | Type | Default | Function | Maps To | Validation |
|----|----|----|----|----|----|
| lm_method | character | “gam” | calculate_lm_interaction_s4 | method | ∈ {“lmm”, “gam”, “fpca”, “gee”} |
| lm_pcorr | character | “BH” | calculate_lm_interaction_s4 | pcorr | ∈ {“BH”, “bonferroni”, “hochberg”, “holm”} |
| divergence_ci | numeric | 0.95 | calculate_divergence_s4 | ci | 0 \< x \< 1 |
| jis_use_lm_fdr | logical | TRUE | jackknife_isoform_switching_s4 | use_lm_fdr | ∈ {TRUE, FALSE} |

### OPTIONAL (Nice to Have)

| Config Key | Type | Default | Function | Maps To | Validation |
|----|----|----|----|----|----|
| divergence_control_group | character | NULL | calculate_divergence_s4 | control_group | character or NULL |
| lm_multicorr | character | NULL | calculate_lm_interaction_s4 | multicorr | ∈ {NULL, “hochberg”, “westfall-young”, “benjamini-yekutieli”} |
| lm_corstr | character | NULL | calculate_lm_interaction_s4 | corstr | ∈ {NULL, “ar1”, “exchangeable”, “independence”} |
| lm_storey | logical | FALSE | calculate_lm_interaction_s4 | storey | ∈ {TRUE, FALSE} |
| divergence_progress | logical | FALSE | calculate_divergence_s4 | progress | ∈ {TRUE, FALSE} |

------------------------------------------------------------------------

## TABLE 7: Implementation Checklist

### File Changes Required

| File | Function | Changes | Lines | Status |
|----|----|----|----|----|
| R/orchestration.R | tsenat_config | Add 4 params + validation | ~30 | ⏳ |
| R/orchestration.R | .execute_lm_interaction_s4 | Extract & pass method, pcorr | ~5 | ⏳ |
| R/orchestration.R | .execute_divergence_s4 | Extract & pass ci, control_group | ~5 | ⏳ |
| R/orchestration.R | .execute_jackknife_isoform_switching | Extract & pass use_lm_fdr | ~5 | ⏳ |
| R/s4_functions_divergence.R | .resolve_divergence_parameters | Add ci to return list | ~3 | ⏳ |
| R/s4_functions_jis.R | jackknife_isoform_switching_s4 | Verify use_lm_fdr resolution | ~1 | ✅ |
| R/orchestration.R | tsenat_config (roxygen) | Add parameter documentation | ~40 | ⏳ |
| tests/ | test_tsenat_config.R | Add unit tests for new params | ~50 | ⏳ |

### Testing Checklist

| Test | Function | Scenarios | Status |
|----|----|----|----|
| Parameter acceptance | tsenat_config | All 4 new params accepted | ⏳ |
| Parameter validation | tsenat_config | Invalid values rejected | ⏳ |
| Configuration propagation | .execute\_\* functions | Config → function parameter | ⏳ |
| Result consistency | Full pipeline | Same config → same results | ⏳ |
| Result differences | Full pipeline | Different method → different results | ⏳ |
| Bootstrap methods | Full pipeline | “percentile” vs “bca” behavior | ⏳ |
| LM methods | Full pipeline | “lmm” vs “gam” vs “fpca” | ⏳ |

------------------------------------------------------------------------

## TABLE 8: Impact Matrix

### Which Parameters Affect Which Outputs?

| Parameter | P-values | Effect Sizes | CI Width | Gene Selection | False Discovery | Switching Calls |
|----|:--:|:--:|:--:|:--:|:--:|:--:|
| lm_method | ✅ | ✅ | ❌ | ✅ | ✅ | ✅ |
| lm_pcorr | ✅ | ❌ | ❌ | ✅ | ✅ | ✅ |
| divergence_ci | ❌ | ❌ | ✅ | ❌ | ❌ | ❌ |
| jis_use_lm_fdr | ❌ | ⭐ | ❌ | ✅ | ⭐ | ✅ |
| divergence_control_group | ✅ | ✅ | ✅ | ✅ | ✅ | ⭐ |
| bootstrap_method | ✅ | ✅ | ✅ | ✅ | ✅ | ⭐ |
| threshold | ❌ | ❌ | ❌ | ✅ | ❌ | ✅ |

**Legend:** ✅ = Major effect \| ⭐ = Critical effect \| ❌ = No effect

------------------------------------------------------------------------

**Document Updated:** April 8, 2026 **Version:** 1.0
