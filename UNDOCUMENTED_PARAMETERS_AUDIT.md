# TSENAT Undocumented Parameters Audit

## Statistical Parameters NOT Currently Exposed in tsenat_config()

**Date:** April 8, 2026  
**Scope:** Core functions: calculate_lm_interaction_s4(),
calculate_divergence_s4(), effect_sizes_divergence_s4(),
jackknife_isoform_switching_s4()

------------------------------------------------------------------------

## EXECUTIVE SUMMARY

### Current tsenat_config() Parameters (Exposed)

``` r
tsenat_config <- function(
  q_values = NULL,              # Tsallis q-spectrum
  condition_col = "condition",  # Treatment/condition column
  subject_col = NULL,           # Sample/subject ID for paired designs
  sample_col = "sample",        # Sample identifier
  paired = FALSE,               # Paired design flag
  control = NULL,               # Reference control group
  p_threshold = 0.05,           # P-value threshold
  fdr_threshold = 0.05,         # FDR threshold
  significance_threshold = 0.05,# Effect size significance filter
  bootstrap = FALSE,            # Bootstrap flag
  nboot = 1000,                 # Bootstrap resamples
  bootstrap_method = "percentile", # percentile or bca
  stringency = "medium",        # Filter stringency (low/medium/high)
  nthreads = 1,                 # Parallel threads
  norm = TRUE,                  # Normalization flag
  bootstrap_ci = 0.95,          # Bootstrap confidence interval
  bootstrap_include_diagnostics = TRUE,
  min_valid_frac = 0.75,        # Minimum valid fraction
  norm_method = NULL,           # Normalization method
  pseudocount = NULL,           # Small value to add to pseudocounts
  shrinkage = "none",           # Shrinkage method
  ...                           # Extra args accepted
)
```

------------------------------------------------------------------------

## 1. RECOMMENDED ADDITIONS TO tsenat_config()

### (Critical for reproducibility and statistical inference)

### 1.1 LM INTERACTION ANALYSIS PARAMETERS

#### **method** (LM method selection)

| Aspect | Value |
|----|----|
| **Function** | calculate_lm_interaction_s4() |
| **Parameter** | method |
| **Type** | character |
| **Current Value** | “gam” (default) |
| **Possible Values** | “lmm” (linear mixed model), “gam” (generalized additive model), “fpca” (functional PCA), “gee” (generalized estimating equations) |
| **What it Controls** | Core statistical approach for modeling q\*condition interactions across q-spectrum |
| **Impact on Results** | CRITICAL - Different methods yield different effect estimates and p-values |
| **Reproducibility** | ⚠️ MUST be configurable for reproducibility |
| **Current Implementation** | Hardcoded in orchestration.R, passed only via direct function call |
| **Recommendation** | ✅ **RECOMMENDED** - Add to tsenat_config() |
| **Default Rationale** | GAM is smooth-fitting friendly for entropy curves; LMM for traditional linear modeling |

#### **pcorr** (P-value correction method)

| Aspect | Value |
|----|----|
| **Function** | calculate_lm_interaction_s4() |
| **Parameter** | pcorr |
| **Type** | character |
| **Current Value** | “BH” (Benjamini-Hochberg, default) |
| **Possible Values** | “BH”, “bonferroni”, “hochberg”, “holm” |
| **What it Controls** | Multiple comparison correction for interaction p-values across genes |
| **Impact on Results** | CRITICAL - Controls FDR/family-wise error rate |
| **Reproducibility** | ⚠️ MUST be configurable |
| **Current Implementation** | Hardcoded default; not extracted from config |
| **Recommendation** | ✅ **RECOMMENDED** - Add to tsenat_config() |
| **Related Parameters** | fdr_threshold (already in config) |
| **Note** | Different methods have different stringency levels |

#### **paired** (Already partially in config)

| Aspect | Value |
|----|----|
| **Status** | Already in tsenat_config() ✅ |
| **Enhancement Needed** | Schema validation - when paired=TRUE, enforce q_values \>= 5 |
| **Issue** | Paired designs are mentioned but validation is weakly enforced |

#### **subject_col** (Already in config)

| Aspect            | Value                                               |
|-------------------|-----------------------------------------------------|
| **Status**        | Already in tsenat_config() ✅                       |
| **Current Usage** | Mapped through resolve_slot_param()                 |
| **Enhancement**   | Ensure it’s documented as REQUIRED when paired=TRUE |

#### **condition_col** (Already in config)

| Aspect            | Value                                |
|-------------------|--------------------------------------|
| **Status**        | Already in tsenat_config() ✅        |
| **Current Usage** | Auto-detected if not provided        |
| **Enhancement**   | Should be explicit in paired designs |

------------------------------------------------------------------------

### 1.2 DIVERGENCE COMPUTATION PARAMETERS

#### **bootstrap_method** (Already partially in config)

| Aspect | Value |
|----|----|
| **Function** | calculate_divergence_s4() |
| **Parameter** | method (divergence calculation method) |
| **Type** | character |
| **Current Value** | “percentile” (default) |
| **Possible Values** | “percentile”, “bca” (bias-corrected and accelerated) |
| **What it Controls** | Bootstrap confidence interval calculation method |
| **Critical For** | Skewed distributions (e.g., Tsallis entropy bounded 0 to log N) |
| **Impact on Results** | ⚠️ BCA adjusts for skewness; percentile assumes symmetry |
| **Reproducibility** | CRITICAL for bounded distributions like entropy |
| **Current Implementation** | bootstrap_method in tsenat_config() but “method” parameter in divergence_s4 is different |
| **Recommendation** | ✅ **CLARIFY AND RECOMMEND** - Rename to “divergence_method” or separate bootstrap specifics |
| **Evidence** | Database papers C016 (2005), C030 (2023), S115 (2015) confirm BCA needed for entropy |

#### **ci** (Confidence interval level)

| Aspect | Value |
|----|----|
| **Function** | calculate_divergence_s4() via .calculate_divergence() |
| **Parameter** | ci |
| **Type** | numeric |
| **Current Value** | 0.95 (hardcoded in internal function) |
| **Possible Values** | 0.90, 0.95, 0.99, etc. |
| **What it Controls** | Confidence interval width (e.g., 95% vs 99%) |
| **Impact on Results** | Changes width of divergence CIs (wider at 0.99, narrower at 0.90) |
| **Reproducibility** | ⚠️ Should be configurable, currently hardcoded |
| **Current Implementation** | Not exposed; hardcoded as 0.95 |
| **Recommendation** | ✅ **RECOMMENDED** - Add to tsenat_config() |
| **Default Rationale** | 0.95 is standard, but should be user-controllable |

#### **control_group** (Treatment reference level)

| Aspect | Value |
|----|----|
| **Function** | calculate_divergence_s4() |
| **Parameter** | control_group |
| **Type** | character |
| **Current Value** | NULL (all groups compared) |
| **What it Controls** | Which group is used as reference for divergence comparisons |
| **Impact on Results** | Determines comparison baseline (e.g., normal vs disease vs disease vs normal) |
| **Reproducibility** | CRITICAL - Must be explicit for pairwise comparisons |
| **Current Implementation** | Accepted via direct function call, can be in config |
| **Recommendation** | ✅ **RECOMMENDED** - Document clearly in tsenat_config() |
| **Related** | control parameter in tsenat_config() has different purpose (paired control group) |
| **Note** | Distinct from “control” (paired reference); may need separate slot |

------------------------------------------------------------------------

### 1.3 JACKKNIFE ISOFORM SWITCHING PARAMETERS

#### **n_bootstrap** (Jackknife bootstrap resamples)

| Aspect | Value |
|----|----|
| **Function** | jackknife_isoform_switching_s4() |
| **Parameter** | n_bootstrap |
| **Type** | numeric |
| **Current Value** | 1000 (default) |
| **Possible Values** | 50, 100, 500, 1000, 2000, … (≥ 50 recommended) |
| **What it Controls** | Number of bootstrap replicates for isoform switching confidence estimation |
| **Impact on Results** | Affects CI width and stability of switching detection |
| **Reproducibility** | CRITICAL - Must match between runs |
| **Current Implementation** | Resolved from config via resolve_slot_param() |
| **Config Status** | n_bootstrap in tsenat_config() ✅ ALREADY EXPOSED |
| **Recommendation** | ✅ Already correct; Ensure validation for ≥ 50 |

#### **threshold** (Isoform switching detection threshold)

| Aspect | Value |
|----|----|
| **Function** | jackknife_isoform_switching_s4() |
| **Parameter** | threshold |
| **Type** | numeric (percentage) |
| **Current Value** | 90 (function default) |
| **Possible Values** | 50-95 (common: 70, 80, 90) |
| **What it Controls** | Percentile threshold for isoform switching confidence |
| **Interpretation** | E.g., 90 = “90% of bootstrap replicates show consistent switching” |
| **Impact on Results** | CRITICAL - Stringency of switching detection |
| **Reproducibility** | ⚠️ Must be configurable, currently hardcoded |
| **Current Implementation** | Resolved from config via resolve_slot_param() |
| **Config Status** | threshold in tsenat_config() ✅ ALREADY EXPOSED (but may need explicit documentation) |
| **Recommendation** | ✅ Already in config; Clarify in documentation |

#### **lm_p_threshold** (LM interaction p-value filter)

| Aspect | Value |
|----|----|
| **Function** | jackknife_isoform_switching_s4() |
| **Parameter** | lm_p_threshold |
| **Type** | numeric |
| **Current Value** | 0.05 (function default) |
| **Possible Values** | 0.01, 0.05, 0.10 |
| **What it Controls** | P-value threshold for including LM interaction results in switching analysis |
| **Impact on Results** | Filters which genes’ interactions are used for switching detection |
| **Reproducibility** | CRITICAL - Must match between runs |
| **Current Implementation** | Resolved from config via resolve_slot_param() |
| **Config Status** | lm_p_threshold in tsenat_config() ✅ ALREADY EXPOSED (but may be undocumented) |
| **Recommendation** | ✅ Already in config; Ensure prominent documentation |

#### **log_base** (Log transformation base)

| Aspect | Value |
|----|----|
| **Function** | jackknife_isoform_switching_s4() |
| **Parameter** | log_base |
| **Type** | numeric |
| **Current Value** | exp(1) (natural log, ~2.718) |
| **Possible Values** | 2, exp(1), 10 |
| **What it Controls** | Logarithm base for isoform count normalization |
| **Impact on Results** | Changes scale of isoform abundance calculations |
| **Reproducibility** | CRITICAL for reproducibility |
| **Current Implementation** | Resolved from config via resolve_slot_param() |
| **Config Status** | log_base in tsenat_config() ✅ ALREADY EXPOSED |
| **Recommendation** | ✅ Already in config; Clarify documentation |

#### **pseudocount** (Zero handling for log transformation)

| Aspect | Value |
|----|----|
| **Function** | jackknife_isoform_switching_s4() |
| **Parameter** | pseudocount |
| **Type** | numeric |
| **Current Value** | 0 (function default) or NULL (config default) |
| **Possible Values** | 0, 0.5, 1.0 |
| **What it Controls** | Small constant added before log transformation to avoid log(0) |
| **Impact on Results** | Affects stability of low-abundance isoforms |
| **Reproducibility** | CRITICAL - Must match between runs |
| **Current Implementation** | Resolved from config via resolve_slot_param() |
| **Config Status** | pseudocount in tsenat_config() ✅ ALREADY EXPOSED |
| **Recommendation** | ✅ Already in config; Ensure default is consistent |

#### **use_lm_fdr** (LM FDR vs raw p-value)

| Aspect | Value |
|----|----|
| **Function** | jackknife_isoform_switching_s4() |
| **Parameter** | use_lm_fdr |
| **Type** | logical |
| **Current Value** | TRUE (function default) |
| **What it Controls** | Whether to use FDR-corrected (TRUE) or raw p-values (FALSE) from LM results |
| **Impact on Results** | CRITICAL - FDR vs raw determines stringency of filtering |
| **Reproducibility** | ⚠️ Should be configurable |
| **Current Implementation** | Hardcoded in function signature |
| **Config Status** | NOT in tsenat_config() ❌ |
| **Recommendation** | ✅ **RECOMMENDED** - Add to tsenat_config() |

------------------------------------------------------------------------

### 1.4 EFFECT SIZE PARAMETERS

#### **enrich_per_q_pattern** (Q-pattern enrichment)

| Aspect | Value |
|----|----|
| **Function** | effect_sizes_divergence_s4() |
| **Parameter** | enrich_per_q_pattern |
| **Type** | logical |
| **Current Value** | TRUE (default) |
| **What it Controls** | Whether to compute per-q divergence patterns as enrichment features |
| **Impact on Results** | Adds pattern columns (div_increase_with_q, etc.) to effect size output |
| **Reproducibility** | Should be configurable for consistency |
| **Current Implementation** | Resolved from config via resolve_slot_param() |
| **Config Status** | enrich_per_q_pattern in tsenat_config() ✅ ALREADY EXPOSED |
| **Recommendation** | ✅ Already in config |

------------------------------------------------------------------------

## 2. OPTIONAL ADDITIONS TO tsenat_config()

### (Improve control; less critical for general use)

### 2.1 LM ANALYSIS - ADVANCED OPTIONS

#### **multicorr** (Multi-test correction method)

| Aspect | Value |
|----|----|
| **Parameter** | multicorr |
| **Type** | character |
| **Current Value** | NULL (uses pcorr instead) |
| **Possible Values** | “hochberg”, “westfall-young”, “benjamini-yekutieli” |
| **What it Controls** | Alternative multiple comparison correction (more conservative than pcorr) |
| **Use Case** | When stricter control than FDR is needed |
| **Recommendation** | ⭕ **OPTIONAL** - Advanced users only |
| **Current Status** | Can be passed via direct call or … args |

#### **corstr** (Correlation structure for GEE)

| Aspect | Value |
|----|----|
| **Parameter** | corstr |
| **Type** | character |
| **Current Value** | NULL (not used unless method=“gee”) |
| **Possible Values** | “ar1” (autoregressive), “exchangeable”, “independence” |
| **What it Controls** | Assumed correlation structure in GEE models |
| **Use Case** | When method=“gee” is selected |
| **Recommendation** | ⭕ **OPTIONAL** - Only relevant for GEE method |
| **Current Status** | Hardcoded; could be configurable |

#### **storey** (Storey’s q-value method)

| Aspect | Value |
|----|----|
| **Parameter** | storey |
| **Type** | logical |
| **Current Value** | FALSE (default) |
| **What it Controls** | Whether to use Storey’s q-value method for p-value adjustment |
| **Use Case** | For large-scale multiple testing scenarios |
| **Recommendation** | ⭕ **OPTIONAL** - Specialized use case |

#### **wy_randomizations** (Westfall-Young permutations)

| Aspect | Value |
|----|----|
| **Parameter** | wy_randomizations |
| **Type** | numeric |
| **Current Value** | 1000 (default) |
| **What it Controls** | Number of permutations for Westfall-Young correction |
| **Use Case** | When multicorr=“westfall-young” |
| **Recommendation** | ⭕ **OPTIONAL** - Only used with specific multicorr method |

#### **adaptive_knots** (GAM adaptive knot placement)

| Aspect | Value |
|----|----|
| **Parameter** | adaptive_knots |
| **Type** | logical |
| **Current Value** | TRUE (default) |
| **What it Controls** | Whether GAM uses adaptive knot placement vs fixed spacing |
| **Use Case** | Fine-tuning GAM smoothness when method=“gam” |
| **Recommendation** | ⭕ **OPTIONAL** - Advanced GAM tuning |

#### **bias_correction** (LM bias correction)

| Aspect | Value |
|----|----|
| **Parameter** | bias_correction |
| **Type** | logical |
| **Current Value** | TRUE (default) |
| **What it Controls** | Whether to apply bias correction in LM coefficient estimates |
| **Recommendation** | ⭕ **OPTIONAL** - Keep TRUE for most analyses |

#### **regularization** (LM regularization method)

| Aspect               | Value                                            |
|----------------------|--------------------------------------------------|
| **Parameter**        | regularization                                   |
| **Type**             | character                                        |
| **Possible Values**  | “pca”, “lasso”, “elasticnet”, “gamsel”, “spline” |
| **Current Value**    | Depends on implementation                        |
| **What it Controls** | Regularization approach to prevent overfitting   |
| **Recommendation**   | ⭕ **OPTIONAL** - Advanced model tuning          |

------------------------------------------------------------------------

### 2.2 DIVERGENCE - ADVANCED OPTIONS

#### **pvalue** (P-value calculation method in LM)

| Aspect               | Value                                         |
|----------------------|-----------------------------------------------|
| **Parameter**        | pvalue                                        |
| **Type**             | character                                     |
| **Possible Values**  | “satterthwaite”, “lrt”, “both”                |
| **Current Value**    | Depends on context                            |
| **What it Controls** | Method for p-value calculation from LM        |
| **Recommendation**   | ⭕ **OPTIONAL** - LM-specific advanced option |

#### **norm** (Normalization flag)

| Aspect               | Value                                              |
|----------------------|----------------------------------------------------|
| **Function**         | calculate_divergence_s4() (and others)             |
| **Parameter**        | norm                                               |
| **Type**             | logical                                            |
| **Current Value**    | TRUE (default)                                     |
| **What it Controls** | Whether to normalize before divergence calculation |
| **Status**           | Already in tsenat_config() ✅                      |
| **Recommendation**   | ✅ Keep normalized                                 |

------------------------------------------------------------------------

## 3. PARAMETERS THAT SHOULD STAY IN DIRECT FUNCTION CALLS ONLY

### (Internal/advanced; not suitable for config)

### 3.1 IMPLEMENTATION DETAILS

| Parameter | Function | Reason | Alternative |
|----|----|----|----|
| min_obs | .calculate_lm_interaction | Internal validation (minimum observations per group) | Not configurable; hardcoded for quality |
| pvalue | .calculate_lm_interaction | Internal p-value calculation method | Advanced LM option; leave as internal |
| assay_name | .calculate_lm_interaction | Which assay matrix to use | Auto-detected; should not be user-configurable |
| return_model_data | calculate_lm_interaction_s4 | Whether to return raw model objects (memory overhead) | Advanced users can pass directly |
| formula | calculate_lm_interaction_s4 | Custom formula specification | Advanced; pass directly to function |
| group_col | .calculate_divergence | Auto-detected from config | Use condition_col from config |
| verbose | All functions | Logging verbosity | Handled per-function call |
| progress | .calculate_divergence | Progress bar display | Handled per-function call |
| gene_col, isoform_col | jackknife_isoform_switching_s4 | Auto-detected from rowData | Auto-detection sufficient; no config needed |
| subject_col (in JIS) | jackknife_isoform_switching_s4 | Auto-detected if not paired | General subject_col from config |

------------------------------------------------------------------------

## SUMMARY TABLE: Parameter Categorization

### RECOMMENDED Additions (Must Add)

| Function | Parameter | Type | Default | Priority |
|----|----|----|----|----|
| calculate_lm_interaction_s4 | method | character | “gam” | HIGH |
| calculate_lm_interaction_s4 | pcorr | character | “BH” | HIGH |
| calculate_divergence_s4 | ci | numeric | 0.95 | HIGH |
| effect_sizes_divergence_s4 | enrich_per_q_pattern | logical | TRUE | ✅ Already in |
| jackknife_isoform_switching_s4 | use_lm_fdr | logical | TRUE | HIGH |
| jackknife_isoform_switching_s4 | threshold | numeric | 90 | ✅ Already in |
| jackknife_isoform_switching_s4 | lm_p_threshold | numeric | 0.05 | ✅ Already in |

### OPTIONAL Additions (Nice to Have)

| Function                    | Parameter         | Type      | Default | Priority |
|-----------------------------|-------------------|-----------|---------|----------|
| calculate_lm_interaction_s4 | multicorr         | character | NULL    | LOW      |
| calculate_lm_interaction_s4 | corstr            | character | NULL    | LOW      |
| calculate_lm_interaction_s4 | storey            | logical   | FALSE   | LOW      |
| calculate_lm_interaction_s4 | wy_randomizations | numeric   | 1000    | LOW      |
| calculate_lm_interaction_s4 | adaptive_knots    | logical   | TRUE    | LOW      |

### ALREADY EXPOSED (No Change Needed)

``` r

# In tsenat_config() - confirmed already exposed:
q_values
condition_col
subject_col
paired
nthreads
nboot
bootstrap
bootstrap_method
norm
pseudocount
bootstrap_ci
significance_threshold
fdr_threshold
n_bootstrap
threshold (jackknife)
lm_p_threshold (jackknife)
log_base
enrich_per_q_pattern
```

------------------------------------------------------------------------

## IMPLEMENTATION RECOMMENDATIONS

### Step 1: Add Critical Parameters to tsenat_config()

``` r
tsenat_config <- function(
    # ... existing parameters ...
    
    # LM INTERACTION PARAMETERS
    lm_method = "gam",                    # NEW: "lmm", "gam", "fpca", "gee"
    lm_pcorr = "BH",                      # NEW: "BH", "bonferroni", "hochberg", "holm"
    
    # DIVERGENCE PARAMETERS
    divergence_ci = 0.95,                 # NEW: Confidence interval level
    divergence_control = NULL,            # NEW: Reference group for divergence
    divergence_bootstrap_method = NULL,   # CLARIFY: Use bootstrap_method instead
    
    # JACKKNIFE PARAMETERS
    jis_use_lm_fdr = TRUE,               # NEW: Use FDR vs raw p-values
    
    # ... other parameters ...
)
```

### Step 2: Update Documentation

- Clarify which parameters are **reproducibility-critical** (MUST match
  between runs)
- Add examples showing paired vs unpaired configurations
- Document q-value spectrum requirements (≥5 or ≥41 for paired designs)

### Step 3: Enhance Validation

- Enforce q_values \>= 5 when paired=TRUE
- Validate lm_method against installed packages (nlme, mgcv, geepack)
- Warn if divergence_bootstrap_method=“percentile” with bounded data
  (entropy)

### Step 4: Orchestration Updates

- Extract new parameters from config in .execute\_\* functions
- Pass through to S4 wrappers instead of hardcoding

------------------------------------------------------------------------

## REPRODUCIBILITY CHECKLIST

When setting up TSENAT for reproducible analysis, users MUST specify: -
\[ \] q_values (≥5, ≥41 for paired) - \[ \] condition_col - \[ \]
subject_col (if paired=TRUE) - \[ \] paired design flag - \[ \] control
group reference - \[ \] lm_method (gam vs lmm vs others) - \[ \] pcorr
(p-value correction) - \[ \] fdr_threshold - \[ \] divergence_ci
(confidence interval) - \[ \] bootstrap_method (percentile vs bca, **BCA
for entropy!**) - \[ \] n_bootstrap replicates - \[ \] threshold for
isoform switching - \[ \] lm_p_threshold for filtering - \[ \] nthreads
(document for reproducibility)

------------------------------------------------------------------------

## REFERENCES

### Database Papers Confirming BCA Need for Entropy

- **C016 (2005)**: Bootstrap analysis of Tsallis entropy
- **C030 (2023)**: Confidence interval construction for bounded
  statistics
- **S115 (2015)**: Skewness-adjusted bootstrap methods

### Internal Code References

- LM parameters: `R/linear_models_core.R` L275-340
- Divergence parameters: `R/divergence_core.R` L437-500
- Effect size parameters: `R/divergence_effect_sizes.R` L114-180
- Jackknife parameters: `R/s4_functions_jis.R` L173-450
- Orchestration calls: `R/orchestration.R` L725-875
- Parameter resolution: `R/s4_functions_lm.R` L255-340,
  `R/s4_functions_divergence.R` L181-250

------------------------------------------------------------------------

**End of Audit Report**
