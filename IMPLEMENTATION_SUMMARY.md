# QUICK REFERENCE: Missing Parameters by Function

## What Needs to be Added to tsenat_config()

------------------------------------------------------------------------

### 🔴 CRITICAL ADDITIONS (Affect reproducibility & results significantly)

#### calculate_lm_interaction_s4()

``` r

# Add to tsenat_config():
lm_method = "gam"              # Currently hardcoded; should be configurable
                               # Choices: "lmm", "gam", "fpca", "gee"
                               # Controls: Statistical approach for q*condition interactions
                               # HIGH IMPACT: Different methods give different p-values

lm_pcorr = "BH"                # Currently hardcoded; should be configurable  
                               # Choices: "BH", "bonferroni", "hochberg", "holm"
                               # Controls: Multiple comparison correction
                               # HIGH IMPACT: Changes FDR/FWER control
```

#### calculate_divergence_s4()

``` r

# Add to tsenat_config():
divergence_ci = 0.95           # Currently hardcoded; should be configurable
                               # Controls: Confidence interval width (e.g., 95% vs 99%)
                               # HIGH IMPACT: Changes CI bounds
                               
divergence_control_group = NULL # Currently can be passed, should be in config
                               # Controls: Reference group for pairwise comparisons
                               # HIGH IMPACT: Changes comparison baseline
```

#### jackknife_isoform_switching_s4()

``` r

# Add to tsenat_config():
jis_use_lm_fdr = TRUE          # Currently hardcoded; should be configurable
                               # Controls: Use FDR-corrected (TRUE) vs raw p-values (FALSE)
                               # HIGH IMPACT: Changes filtering stringency
```

#### effect_sizes_divergence_s4()

``` r

# Already in config ✅:
enrich_per_q_pattern = TRUE
significance_threshold = 0.05
```

------------------------------------------------------------------------

### 🟡 OPTIONAL ADDITIONS (Improve control; specialized use)

#### calculate_lm_interaction_s4() - Advanced options

``` r

lm_multicorr = NULL            # Multi-test correction: "hochberg", "westfall-young", "benjamini-yekutieli"
lm_corstr = NULL               # Correlation structure (GEE): "ar1", "exchangeable", "independence"
lm_storey = FALSE              # Use Storey's q-value method
lm_wy_randomizations = 1000    # Westfall-Young permutation count
lm_adaptive_knots = TRUE       # GAM adaptive knot placement
```

------------------------------------------------------------------------

### ✅ ALREADY IN CONFIG (No changes needed)

``` r

# These parameters are already exposed in tsenat_config():
q_values                        # Tsallis q-spectrum
condition_col                   # Treatment/condition column
subject_col                     # Sample ID for paired designs
paired                          # Paired design flag
control                         # Reference group (paired)
p_threshold                     # P-value threshold
fdr_threshold                   # FDR threshold  
significance_threshold          # Effect size threshold
bootstrap                       # Bootstrap flag
nboot                          # Bootstrap resamples
bootstrap_method               # "percentile" or "bca"
bootstrap_ci                   # Confidence interval (0.95)
nthreads                        # Parallel threads
norm                           # Normalization flag
pseudocount                    # Small value for log(0) prevention
log_base                       # Log transformation base
n_bootstrap                    # JIS bootstrap resamples
threshold                      # JIS switching detection threshold (default: 90)
lm_p_threshold                 # JIS LM p-value filter (default: 0.05)
enrich_per_q_pattern           # Effect size enrichment flag
```

------------------------------------------------------------------------

### 🔧 IMPLEMENTATION PRIORITY

| Priority | Parameter | Function | Add To Config | Lines of Code |
|----|----|----|----|----|
| 1 (HIGH) | lm_method | calculate_lm_interaction_s4 | ✅ YES | ~10 |
| 1 (HIGH) | lm_pcorr | calculate_lm_interaction_s4 | ✅ YES | ~10 |
| 1 (HIGH) | divergence_ci | calculate_divergence_s4 | ✅ YES | ~10 |
| 1 (HIGH) | jis_use_lm_fdr | jackknife_isoform_switching_s4 | ✅ YES | ~10 |
| 2 (MEDIUM) | divergence_control_group | calculate_divergence_s4 | ✅ YES | ~5 |
| 3 (LOW) | lm_multicorr | calculate_lm_interaction_s4 | ⭕ OPT | ~5 |
| 3 (LOW) | lm_corstr | calculate_lm_interaction_s4 | ⭕ OPT | ~5 |

**Total implementation effort:** ~30 lines of code + documentation
updates

------------------------------------------------------------------------

### ⚠️ CRITICAL NOTES FOR USERS

#### When using Bounded Distributions (Entropy)

- **MUST set** `bootstrap_method = "bca"` instead of “percentile”
- **Reason:** Percentile bootstrap assumes symmetric distribution; BCA
  corrects for skewness
- **Impact:** With percentile, CI may not contain point estimate for
  skewed distributions
- **Evidence:** Database papers C016, C030, S115 confirm this pattern

#### When using Paired Designs

- **MUST set** `q_values = seq(0, 2, by = 0.05)` (41+ values) instead of
  default seq(0, 2, by = 0.5)
- **MUST set** `subject_col = "paired_samples"` or appropriate column
  name
- **MUST set** `control = "normal"` or reference group name
- **MUST set** `paired = TRUE`
- **Reason:** Interaction testing requires sufficient q-spectrum
  granularity

#### When Comparing Methods

- **Document** the `lm_method` used (“gam” vs “lmm” produce different
  results)
- **Document** the `pcorr` used (“BH” vs “bonferroni” have different
  stringency)
- **Document** the `bootstrap_method` used (“percentile” vs “bca” for
  skewed data)
- **Use** `git` or version control to track config variations

------------------------------------------------------------------------

### 📋 CHECKLIST FOR NEW ANALYSES

``` r

# Recommended minimal config for reproducible analysis:
config <- tsenat_config(
  # Metadata mapping
  q_values = seq(0, 2, by = 0.05),     # ≥5 required; 41+ for paired
  condition_col = "condition",          
  subject_col = "paired_samples",       # If paired=TRUE
  sample_col = "sample",
  paired = FALSE,                       # Or TRUE if applicable
  control = NULL,                       # Reference group if applicable
  
  # Statistical thresholds
  p_threshold = 0.05,
  fdr_threshold = 0.05,
  significance_threshold = 0.05,
  lm_p_threshold = 0.05,                # For JIS filtering
  
  # Methods (NEW - MUST SPECIFY)
  lm_method = "gam",                    # ← NEW: Add this
  lm_pcorr = "BH",                      # ← NEW: Add this
  bootstrap_method = "bca",             # ← CRITICAL: Use "bca" for entropy!
  
  # Bootstrap/resampling
  bootstrap = FALSE,
  nboot = 1000,
  bootstrap_ci = 0.95,
  n_bootstrap = 1000,
  
  # Divergence (NEW - OPTIONAL but recommended)
  divergence_ci = 0.95,                 # ← NEW: Add this
  divergence_control_group = NULL,      # ← NEW: Add this if needed
  
  # Isoform switching
  threshold = 90,                       # Switching detection threshold
  jis_use_lm_fdr = TRUE,                # ← NEW: Add this
  
  # Data processing
  norm = TRUE,
  norm_method = NULL,
  pseudocount = 0,
  log_base = exp(1),
  shrinkage = "none",
  stringency = "medium",
  
  # Computation
  nthreads = 1
)

# Then create analysis:
analysis <- build_analysis_s4(
  readcounts = readcounts,
  metadata = metadata,
  tx2gene = annotation,
  config = config
)

# And run pipeline:
result <- tsenat(analysis)
```

------------------------------------------------------------------------

### 📊 PARAMETER IMPACT MATRIX

[TABLE]

------------------------------------------------------------------------

**Document Generated:** April 8, 2026  
**Audit Scope:** Core statistical parameters in TSENAT v0.1  
**Status:** Ready for implementation
