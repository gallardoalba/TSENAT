# TSENAT Statistical Parameters Audit

## Complete Documentation of Accessible vs. Documented Parameters

**Scope**: Analyzes statistically important parameters used in: -
[`calculate_diversity_s4()`](https://gallardoalba.github.io/TSENAT/reference/calculate_diversity_s4.md) -
[`calculate_divergence_s4()`](https://gallardoalba.github.io/TSENAT/reference/calculate_divergence_s4.md) -
[`tsenat_config()`](https://gallardoalba.github.io/TSENAT/reference/tsenat_config.md)
(documentation source) - Core diversity/divergence computation functions

**Methodology**: Source code search + parameter resolution tracking from
config objects.

------------------------------------------------------------------------

## SECTION 1: DOCUMENTED PARAMETERS (Already in tsenat_config() Help)

These parameters are fully documented in the
[`tsenat_config()`](https://gallardoalba.github.io/TSENAT/reference/tsenat_config.md)
function documentation (man/tsenat_config.Rd).

### 1. q_values

- **Default value**: `seq(0, 2, by = 0.5)` (i.e., c(0, 0.5, 1.0, 1.5,
  2.0))
- **Type**: numeric vector
- **What it controls**: The Tsallis entropy q-spectrum range - computes
  entropy at each value of q
- **Statistical importance**:
  - q=0: Species richness (presence/absence only)
  - q=1: Shannon entropy (balanced weighting)
  - q=2: Gini-Simpson index (emphasizes dominant isoforms)
  - q\>1: Increasingly weighted toward dominant isoforms
- **Impact on results**: Different q-values capture different aspects of
  isoform complexity. Multiple q-values provide complete picture of
  diversity spectrum
- **Citation**: Papers I001-I004 (Tsallis entropy theory), S063-S067
  (power analysis showing multi-q informativeness)

### 2. condition_col

- **Default value**: `"condition"`
- **Type**: character
- **What it controls**: Column name in `colData()` containing
  experimental group labels (e.g., “treated”, “control”)
- **Statistical importance**: Defines grouping for divergence/difference
  calculations; errors here cause analysis to fail
- **Impact on results**: All downstream analyses use this for group
  comparisons; wrong column → wrong results
- **Required**: TRUE for most workflows

### 3. subject_col

- **Default value**: `NULL`
- **Type**: character or NULL
- **What it controls**: Column name for paired/repeated design subject
  identifiers
- **Statistical importance**: Enables paired t-test, paired bootstrap,
  repeated measures analysis
- **Impact on results**: Paired designs have higher statistical power
  than unpaired; wrong pairing → invalid p-values
- **Required if**: `paired = TRUE`

### 4. sample_col

- **Default value**: `"sample"`
- **Type**: character
- **What it controls**: Column name for sample identifiers in metadata
- **Statistical importance**: Used for sample-level filtering and
  tracking
- **Impact on results**: Incorrect mapping → samples mismatched or lost

### 5. paired

- **Default value**: `FALSE`
- **Type**: logical
- **What it controls**: Whether to use paired/repeated measures design
- **Statistical importance**: Changes statistical test (paired t-test
  vs. Welch’s t-test); affects power, validity of p-values
- **Impact on results**: Paired=TRUE requires paired samples with
  subject_col; improves power when samples are truly paired
- **Warning in code**: If `paired=TRUE` without `subject_col` or
  `control`, function warns and analysis may fail

### 6. control

- **Default value**: `NULL`
- **Type**: character or NULL
- **What it controls**: Reference group label for difference/divergence
  calculations
- **Statistical importance**: Defines the baseline for comparisons
- **Impact on results**: Choosing different control group changes all
  fold-changes/divergence directions (relative effect)
- **Example**: If condition_col has “control”, “treated1”, “treated2”,
  set `control="control"`

### 7. p_threshold

- **Default value**: `0.05`
- **Type**: numeric
- **What it controls**: Raw p-value threshold for significance in LM
  interaction testing
- **Statistical importance**: Standard alpha level before multiple
  testing correction
- **Impact on results**: Lower threshold → fewer false positives but
  more false negatives; typically 0.05 for FDR pre-filtering

### 8. fdr_threshold

- **Default value**: `0.05`
- **Type**: numeric
- **What it controls**: Adjusted p-value (FDR via Benjamini-Hochberg)
  threshold after multiple testing correction
- **Statistical importance**: Controls family-wise error rate; standard
  for genomic studies is 0.05 or 0.01
- **Impact on results**: Results passing both p_threshold AND
  fdr_threshold are “significant”

### 9. significance_threshold

- **Default value**: `0.05`
- **Type**: numeric
- **What it controls**: Cutoff for effect sizes, assumptions testing
  results, and final filtering
- **Statistical importance**: Used when p_threshold/fdr_threshold not
  applicable; flexible significance standard
- **Impact on results**: Filters final results; 0.05 is conservative,
  0.10 is lenient

### 10. nboot

- **Default value**: `1000`
- **Type**: integer
- **What it controls**: Number of bootstrap resamples for confidence
  interval computation
- **Statistical importance**: More bootstrap samples → wider, more
  accurate CI bounds (diminishing returns after 500)
- **Impact on results**:
  - 100: Fast, CI may be unreliable (minimum)
  - 500-1000: Standard, good balance of speed/accuracy
  - 5000+: Conservative (slower), useful for publication-quality results
- **Mathematical basis**: BCa method convergence requires ~500+ samples
  (papers S111, S114); percentile method slightly less demanding

### 11. bootstrap_method

- **Default value**: `"percentile"`
- **Type**: character
- **What it controls**: Bootstrap confidence interval construction
  algorithm
- **Options**:
  - `"percentile"`: Simple percentile bootstrap, assumes symmetric
    distribution
  - `"bca"`: Bias-corrected and accelerated bootstrap, adjusts for
    skewness and bias
- **Statistical importance**:
  - Distribution shape matters: Tsallis entropy is naturally bounded
    \[0, log(m)\], often skewed
  - Percentile CI may not contain point estimate when distribution is
    skewed (invalid coverage)
  - BCA corrects this but costs ~20-50% more computation
- **Impact on results**:
  - percentile: Fast but potentially invalid CIs for skewed data (common
    with entropy)
  - bca: More valid CIs especially for bounded distributions
    (RECOMMENDED for Tsallis entropy)
- **Citation**: Papers C016 (2005), C030 (2023), S115 (2015) show
  bounded Tsallis entropy exhibits skewness requiring BCA
- **Decision rule**: Use `"bca"` for entropy analysis; use
  `"percentile"` only if you’ve verified symmetry

### 12. stringency

- **Default value**: `"medium"`
- **Type**: character
- **Options**: `"lenient"`, `"medium"` (default), `"severe"`
- **What it controls**: Transcript filtering level during
  [`filter_analysis_s4()`](https://gallardoalba.github.io/TSENAT/reference/filter_analysis_s4.md)
  execution
- **Statistical importance**: Filters low-abundance transcripts before
  analysis; affects downstream gene list, statistical power
- **Impact on results**:
  - lenient: Keeps more transcripts, higher noise, may include spurious
    signals
  - medium: Balance of noise/signal (recommended baseline)
  - severe: Strict filtering, high confidence only, may lose weak
    signals
- **Example**: “severe” might require ≥5 reads per transcript in ≥75% of
  samples; “lenient” might require ≥1 read in ≥1 sample

### 13. nthreads

- **Default value**: `1`
- **Type**: integer
- **What it controls**: Number of CPU threads for parallel computation
- **Statistical importance**: No statistical effect (only speed)
- **Impact on results**: Speed improvement ~linear up to number of
  available cores; no accuracy loss
- **Practical guidance**: Set to `detectCores() - 1` on multi-core
  systems

------------------------------------------------------------------------

## SECTION 2: UNDOCUMENTED BUT ACCESSIBLE PARAMETERS (Should Be in tsenat_config() Help)

These parameters **ARE** used by core functions and **ARE**
configurable, but **NOT** documented in
[`tsenat_config()`](https://gallardoalba.github.io/TSENAT/reference/tsenat_config.md)
function documentation. Users must discover them by reading individual
function help pages or source code.

### TIER 1: CRITICAL - MUST DOCUMENT (HIGH IMPACT ON RESULTS)

#### 1. bootstrap

- **Current location**: Only documented in
  [`calculate_diversity_s4()`](https://gallardoalba.github.io/TSENAT/reference/calculate_diversity_s4.md)
  @param section  
- **Default value**: `FALSE`
- **Type**: logical
- **What it controls**: Whether to compute bootstrap confidence
  intervals around entropy/divergence estimates
- **Statistical importance**: WITHOUT bootstrap: only point estimates
  (no uncertainty quantification). WITH bootstrap: point + CI bounds
- **Impact on results**:
  - `FALSE`: Fast (1x computation), output has only diversity/divergence
    assay, NO ci_lower/ci_upper assays
  - `TRUE`: Slower (5-10x depending on nboot), adds ci_lower/ci_upper
    assays with bounds
- **Accessibility**:
  - Via function:
    `calculate_diversity_s4(analysis, bootstrap=TRUE, nboot=1000)`
  - Via config:
    `seq_config <- tsenat_config(); seq_config$bootstrap <- TRUE`
  - PROBLEM: NOT in tsenat_config() signature, can’t pass directly
    during config creation
- **When to use**: Always TRUE for publication results; FALSE for
  exploratory/quick checks
- **Recommendation**: ADD to tsenat_config() with default=FALSE for
  backward compatibility

#### 2. pseudocount

- **Current location**: Only documented in
  [`calculate_diversity_s4()`](https://gallardoalba.github.io/TSENAT/reference/calculate_diversity_s4.md)
  and `.calculate_diversity()` help
- **Default value**: `0` or `NULL` (both interpreted as no pseudocount)
- **Type**: numeric or character (“auto”)
- **What it controls**: Pseudocount added to transcript counts before
  entropy calculation
- **Statistical importance**: CRITICAL for numerical stability
  - With zero counts: log(0) = -∞ causing NaN/Inf entropy values
  - Pseudocount prevents zero-division: (count + pseudocount) / (total +
    pseudo\*n_transcripts)
  - Acts as regularization: lower values = more variance, higher values
    = more bias
- **Options**:
  - `0` or `NULL`: No pseudocount (raw entropy, risky for sparse genes
    with zero-count samples)
  - `0.5`: Half-count pseudocount (edgeR default, standard practice)
  - `1.0`: Full pseudocount (DESeq2 default, more conservative)
  - `"auto"`: Automatically estimated based on library size via
    `.estimate_pseudocount()` (RECOMMENDED)
- **Impact on results**:
  - `0`: Entropy values for genes with any zero-count sample will be NaN
  - `0.5-1.0`: Prevents NaN, adds small bias, reduces variance (Bayesian
    regularization)
  - `"auto"`: Adapts pseudocount to sequencing depth (principled
    approach, used by edgeR/DESeq2)
- **Accessibility**:
  - Via function: `calculate_diversity_s4(analysis, pseudocount=0.5)`
  - Via config: Would need to expose in tsenat_config()
- **Citation**: Robinson et al. 2010 (edgeR paper S001), Li et al. 2023
  (bootstrapped entropy validation)
- **Recommendation**: ADD to tsenat_config() with three options:
  numeric, “auto”, NULL
  - CRITICAL: Document that “auto” is RECOMMENDED for real data

#### 3. bootstrap_ci

- **Current location**: Only documented in
  [`calculate_diversity_s4()`](https://gallardoalba.github.io/TSENAT/reference/calculate_diversity_s4.md)
  @param section
- **Default value**: `0.95`
- **Type**: numeric (0 \< bootstrap_ci \< 1)
- **What it controls**: Confidence interval level (0.90 = 90% CI, 0.95 =
  95% CI, 0.99 = 99% CI)
- **Statistical importance**: Users may require different CI levels
  depending on publication standards or risk tolerance
- **Impact on results**:
  - `0.90`: Narrower CI (faster computation), used when resources
    limited or preliminary results
  - `0.95`: Standard/default (widest common CI), best balance for
    publication
  - `0.99`: Conservative, very wide CI, used when false positive risk
    must be minimized
- **Accessibility**:
  - Via function: `calculate_diversity_s4(analysis, bootstrap_ci=0.99)`
  - Would need tsenat_config() field
- **Recommendation**: ADD to tsenat_config() with default=0.95

#### 4. shrinkage

- **Current location**: Only documented in
  [`calculate_diversity_s4()`](https://gallardoalba.github.io/TSENAT/reference/calculate_diversity_s4.md)
  and `.calculate_diversity()` help
- **Default value**: `"none"`
- **Type**: character
- **Options**: `"none"`, `"empirical_bayes"`
- **What it controls**: Entropy estimation variance reduction via
  Bayesian shrinkage toward global mean
- **Statistical importance**: Improves stability for genes with few
  expressed isoforms (1-2 isoforms)
  - Gene with 1 isoform: entropy = 0 (deterministic), high variance
    estimate
  - Shrinkage: borrows information from all other genes, reduces
    estimate variance
- **Impact on results**:
  - `"none"`: Raw entropy per gene, no information borrowing (higher
    variance on sparse genes)
  - `"empirical_bayes"`: Shrunk estimates toward global mean (lower
    variance, reduced false positives on noise)
- **When to use**:
  - Use `"none"` for small datasets (\<50 genes) where global mean
    unreliable
  - Use `"empirical_bayes"` for moderate/large datasets (\>100 genes)
    where borrowing strength valid
- **Mathematical basis**: Papers S004-S006 validate empirical Bayes
  strengthening for entropy estimation
- **Accessibility**: Via function:
  `calculate_diversity_s4(analysis, shrinkage="empirical_bayes")`
- **Recommendation**: ADD to tsenat_config() with default=“none”

#### 5. norm_method

- **Current location**: PARTIALLY documented in calculate_diversity_s4()
  help ONLY, NOT in tsenat_config()
- **Default value**: `NULL` (no post-hoc normalization)
- **Type**: character or NULL
- **Options**:
  - `NULL`: No post-hoc normalization (after initial norm applied)
  - `"default"`: Simple range normalization \[0,1\] by theoretical
    maximum
  - `"zscore"`: Z-score normalization per q-value: (diversity - mean) /
    sd
  - `"log_odds_ratio"`: Log-odds ratio relative to max entropy
  - `"relative_reference"`: Divide by reference group mean (requires
    reference_group param)
- **What it controls**: POST-HOC normalization applied AFTER primary
  diversity calculation
- **Statistical importance**: Different methods enable different
  downstream comparisons
  - No normalization: Entropy scale depends on isoform count, not
    comparable across genes
  - Range \[0,1\]: Standardizes within gene, comparable within-gene
    across q-values
  - Z-score: Enables cross-q comparison, removes q-specific scale
    effects
  - Log-odds-ratio: Interpretation relative to random expectation
    (isoform-aware)
  - Relative-reference: Normalization by reference group (control group)
    for comparative analysis
- **Impact on results**: Dramatically changes scale and interpretability
- **Accessibility**: Via function:
  `calculate_diversity_s4(analysis, norm_method="zscore")`
- **Recommendation**: ADD to tsenat_config() with clear documentation of
  when each is appropriate
- **Citation**: Papers B002-B007 validate different normalization
  methods; I023 (Hill numbers) discusses standardization

#### 6. bootstrap_include_diagnostics

- **Current location**: Only documented in
  [`calculate_diversity_s4()`](https://gallardoalba.github.io/TSENAT/reference/calculate_diversity_s4.md)
  and core functions
- **Default value**: `TRUE`
- **Type**: logical
- **What it controls**: Whether to include diagnostic metadata for
  bootstrap CI quality assessment
- **What it enables**: Computation of effective sample size, skewness,
  bias, acceleration factor (BCA term)
- **Statistical importance**: Diagnostics reveal CI reliability
  - High skewness → CI may be asymmetric/invalid
  - High acceleration factor → CI highly sensitive to distribution tails
  - Low effective sample size → too many bootstrap failures, CI
    unreliable
- **Impact on results**:
  - `TRUE`: Adds 5-10% computation cost, includes diagnostic columns in
    output
  - `FALSE`: Faster, cleaner output, diagnostic information lost
- **Use cases**:
  - Research/publication: Use TRUE for quality control
  - Bulk screening: Use FALSE for speed
- **Accessibility**: Via function:
  `calculate_diversity_s4(analysis, bootstrap_include_diagnostics=FALSE)`
- **Recommendation**: ADD to tsenat_config() with clear guidance

#### 7. min_valid_frac

- **Current location**: Only documented in `.calculate_diversity()` core
  function help
- **Default value**: `0.75`
- **Type**: numeric (0-1 range)
- **What it controls**: Minimum fraction of valid bootstrap resamples
  required to flag a CI as reportable
- **Statistical importance**: If \>25% of bootstrap resamples fail (due
  to insufficient data), CI is unreliable
- **Impact on results**: Filtering threshold for confidence bounds that
  are output vs. flagged as NA/invalid
- **Example**: If gene has only 2 isoforms with rare counts, 30% of
  resamples might fail (\>25% threshold) → CI flagged as invalid
- **Accessibility**: Via function:
  `calculate_diversity_s4(analysis, min_valid_frac=0.80)`
- **Recommendation**: ADD to tsenat_config() but mark as “advanced
  parameter”

### TIER 2: IMPORTANT - SHOULD DOCUMENT (MEDIUM IMPACT)

#### 8. norm (already semi-documented)

- **Current location**: Documented in calculate_diversity_s4() but NOT
  in tsenat_config()  
- **Default value**: `TRUE`
- **Type**: logical or character
- **Options**:
  - `TRUE` or `"range"`: Range normalization \[0,1\] by theoretical
    maximum log(m)
  - `FALSE` or `"none"`: No normalization (raw entropy)
  - `"zscore"`: Z-score normalization (backward compat: send to
    norm_method instead)
  - `"log_odds_ratio"`: Log-odds ratio (backward compat: send to
    norm_method instead)
  - `"relative_reference"`: Relative to reference group (backward
    compat: send to norm_method instead)
- **What it controls**: Type of entropy standardization
- **Statistical importance**: Affects scale and interpretability
- **Impact on results**:
  - `TRUE`: Entropy \[0,1\], interpretable (0=mono-, 1=equi-isoform
    expression)
  - `FALSE`: Raw entropy, scale depends on isoform count, not
    interpretable cross-genes
- **Recommendation**: Already exposed functionally but should ENHANCE
  tsenat_config() documentation

#### 9. tpm

- **Current location**: Only documented in
  [`calculate_diversity_s4()`](https://gallardoalba.github.io/TSENAT/reference/calculate_diversity_s4.md)
  help
- **Default value**: `FALSE`
- **Type**: logical
- **What it controls**: Use TPM-normalized abundance data vs. raw counts
  for entropy calculation
- **Statistical importance**: Removes transcript length bias
  - Raw counts: Long transcripts artificially inflate counts/entropy due
    to sequencing depth allocation
  - TPM: Normalized by transcript length per SALMON method, removes bias
- **Impact on results**:
  - `FALSE`: Uses raw counts (typical if already handling length bias
    separately)
  - `TRUE`: Uses assay=“tpm” if available, else raw counts with warning
- **When to use**: Set to TRUE if:
  - Using raw RNA-seq read counts (not pseudo-aligned)
  - Want to remove transcript-length bias
- **Accessibility**: Via function:
  `calculate_diversity_s4(analysis, tpm=TRUE)`
- **Recommendation**: ADD to tsenat_config() with documentation

#### 10. what

- **Current location**: Only documented in `.calculate_diversity()` core
  function
- **Default value**: `"S"`
- **Type**: character
- **Options**:
  - `"S"`: Tsallis entropy H_q
  - `"D"`: Hill numbers (true diversity, numbers equivalent)
- **What it controls**: Output type - which diversity index to compute
- **Statistical importance**: Different indices have different
  interpretations
  - Entropy: Information measure, on log scale
  - Hill numbers: Interpretable as “equivalent number of
    equally-abundant isoforms”
- **Impact on results**: Scale/interpretation changes; Hill=N^(H/log(N))
  approximately
- **Status**: Mostly tested for “S”; “D” may have edge-case issues
- **Recommendation**: Do NOT expose yet; too experimental. Add to
  advanced/internal parameters section

#### 11. effective_length

- **Current location**: Only documented in `.calculate_diversity()` core
  function
- **Default value**: `NULL`
- **Type**: numeric vector or matrix or NULL
- **What it controls**: Effective transcript lengths for length-bias
  correction
- **Statistical importance**: Same motivation as `tpm` parameter -
  removes length bias
- **How it works**: Proportions calculated as (normalized_counts /
  effective_length) to remove length effects
- **Typical source**: SALMON/Kallisto quantification EffectiveLength
  column
- **Impact on results**: With proper effective_length, more balanced
  entropy estimates across genes of different isoform compositions
- **Accessibility**: Via function:
  `calculate_diversity_s4(analysis, effective_length=my_lengths)`
- **Recommendation**: Document in separate “Advanced Length Bias
  Correction” guide, NOT in tsenat_config() (too specialized)

### TIER 3: INTERNAL/ADVANCED - REFERENCE ONLY

#### 12. reference_group

- **Current location**: Only documented in
  [`calculate_diversity_s4()`](https://gallardoalba.github.io/TSENAT/reference/calculate_diversity_s4.md)
  help
- **Default value**: `NULL`
- **Type**: character or NULL
- **What it controls**: Reference group column name when using
  `norm_method="relative_reference"`
- **Impact**: Enables relative normalization analysis (e.g., normalize
  all samples to control group mean)
- **Status**: Internal helper, not user-facing; should NOT expose in
  tsenat_config()
- **Recommendation**: Document as “internal/advanced, for use only with
  norm_method=‘relative_reference’”

#### 13. assayno / genes / metadata

- **Default values**: assayno=1, genes=NULL, metadata=NULL
- **Status**: Internal parameters for advanced/power-user workflows
- **Recommendation**: Document separately as “Advanced Internal
  Parameters”

------------------------------------------------------------------------

## SECTION 3: SUMMARY TABLE - WHAT’S DOCUMENTED VS. WHAT SHOULD BE

| Parameter | Tier | Doc Status | Default | Type | Action |
|----|----|----|----|----|----|
| q_values | Core | ✅ DOCUMENTED | seq(0,2,0.5) | numeric | No change |
| condition_col | Core | ✅ DOCUMENTED | “condition” | character | No change |
| subject_col | Core | ✅ DOCUMENTED | NULL | character | No change |
| sample_col | Core | ✅ DOCUMENTED | “sample” | character | No change |
| paired | Core | ✅ DOCUMENTED | FALSE | logical | No change |
| control | Core | ✅ DOCUMENTED | NULL | character | No change |
| p_threshold | Core | ✅ DOCUMENTED | 0.05 | numeric | No change |
| fdr_threshold | Core | ✅ DOCUMENTED | 0.05 | numeric | No change |
| significance_threshold | Core | ✅ DOCUMENTED | 0.05 | numeric | No change |
| nboot | Core | ✅ DOCUMENTED | 1000 | integer | No change |
| bootstrap_method | Core | ✅ DOCUMENTED | “percentile” | character | **IMPROVE DOCS** (add BCA guidance) |
| stringency | Core | ✅ DOCUMENTED | “medium” | character | No change |
| nthreads | Core | ✅ DOCUMENTED | 1 | integer | No change |
| **bootstrap** | **Tier 1** | ❌ MISSING | FALSE | logical | **ADD TO tsenat_config()** |
| **pseudocount** | **Tier 1** | ❌ MISSING | 0/“auto” | numeric/char | **ADD TO tsenat_config()** |
| **bootstrap_ci** | **Tier 1** | ❌ MISSING | 0.95 | numeric | **ADD TO tsenat_config()** |
| **shrinkage** | **Tier 1** | ❌ MISSING | “none” | character | **ADD TO tsenat_config()** |
| **norm_method** | **Tier 1** | ❌ MISSING | NULL | character | **ADD TO tsenat_config()** |
| **bootstrap_include_diagnostics** | **Tier 1** | ❌ MISSING | TRUE | logical | **ADD TO tsenat_config()** |
| **min_valid_frac** | **Tier 1** | ❌ MISSING | 0.75 | numeric | **ADD TO tsenat_config()** |
| norm | Core | ⚠️ PARTIAL | TRUE | logical/char | **ENHANCE DOCS in config** |
| tpm | Tier 2 | ⚠️ PARTIAL | FALSE | logical | **ADD TO tsenat_config()** |
| effective_length | Tier 2 | ⚠️ PARTIAL | NULL | numeric/matrix | Document separately (advanced) |
| what | Tier 2 | ⚠️ EXPERIMENTAL | “S” | character | Do NOT expose yet |
| reference_group | Tier 3 | ⚠️ INTERNAL | NULL | character | Internal only |
| genes / metadata / assayno | Tier 3 | ⚠️ INTERNAL | NULL/1 | various | Internal only |
| verbose / show_messages | Tier 3 | ⚠️ INTERNAL | TRUE/FALSE | logical | Internal only |

------------------------------------------------------------------------

## SECTION 4: RECOMMENDED DOCUMENTATION UPDATES

### Priority 1: Add Missing Parameters to tsenat_config()

**File to edit**: `R/orchestration.R` function
[`tsenat_config()`](https://gallardoalba.github.io/TSENAT/reference/tsenat_config.md)

Add these parameters to function signature:

``` r
tsenat_config <- function(
  # ... existing parameters ...
  bootstrap = FALSE,              # NEW
  bootstrap_ci = 0.95,            # NEW
  pseudocount = NULL,             # NEW - Can be numeric, "auto", or NULL
  shrinkage = "none",             # NEW
  norm_method = NULL,             # NEW
  bootstrap_include_diagnostics = TRUE,  # NEW
  min_valid_frac = 0.75,          # NEW
  ...
)
```

### Priority 2: Improve Existing Documentation

- **bootstrap_method**: Add section “When to use ‘percentile’ vs ‘bca’”:
  - Percentile: Fast, assumes symmetric distribution, typical use
  - BCA: Recommended for bounded/skewed distributions (like Tsallis
    entropy)
  - Guidance: Use BCA for entropy; use percentile for
    normally-distributed measures
- **norm**: Clarify relationship to norm_method in documentation:
  - `norm`: Primary normalization on/off or method (backward compat)
  - `norm_method`: Post-hoc normalization specifics (more flexible)

### Priority 3: Add Separate Documentation

Create new documentation section in vignette or separate guide:

**“Advanced Parameters Reference”**: - What each parameter does - When
to change defaults - Example: “For sparse GTEx data, use
pseudocount=‘auto’ and shrinkage=‘empirical_bayes’” - Example: “For
publication-quality results, use bootstrap=TRUE, bootstrap_method=‘bca’,
nboot=5000”

------------------------------------------------------------------------

## SECTION 5: CODE LOCATIONS & EVIDENCE

### Parameter Usage Evidence

**bootstrap** - Used in: - R/s4_functions_diversity.R line 206 (param in
calculate_diversity_s4) - R/diversity_core.R line 206 (param in
.calculate_diversity) - R/bootstrap_diversity.R (all bootstrap CI
functions) - Tests: test-core-functions-diversity.R ~20+ tests for
bootstrap=TRUE/FALSE

**pseudocount** - Used in: - R/s4_functions_diversity.R line 51 (@param
documented) - R/diversity_core.R line 206 (param in
.calculate_diversity) - R/bootstrap_diversity.R (all CI computations) -
Tests: test-core-functions-diversity.R bootstrap tests verify
pseudocount behavior

**bootstrap_ci** - Used in: - R/s4_functions_diversity.R line 55 (@param
documented) - R/diversity_core.R line 206 (param in
.calculate_diversity) - All bootstrapping code defaults to 0.95

**norm_method** - Used in: - R/s4_functions_diversity.R line 29-32
(@param documented but only there) - R/bootstrap_diversity.R
post-processing - Tests: test-core-functions-diversity.R ~10+ tests for
different norm_methods

### Configuration Resolution Pattern

All parameters follow this pattern in **R/s4_functions_diversity.R**
lines 706-784:

``` r

# Parameter priority: explicit > @config > default
param_value <- resolve_slot_param(param, analysis@config, "param_name", default)
```

This means adding parameter to tsenat_config() automatically enables: 1.
Setting default in tsenat_config() 2. Overriding in function call 3.
Reading from <analysis@config> slot

------------------------------------------------------------------------

## SECTION 6: RISKS OF CURRENT UNDOCUMENTED STATE

### For Users:

1.  **Bootstrap not enabled by default**: Users don’t get confidence
    intervals unless they discover `bootstrap=TRUE` parameter
2.  **Pseudocount=0 causes NaN for sparse genes**: Silent failures on
    sparse data without explicit pseudocount setting
3.  **Percentile bootstrap used for skewed data**: Invalid CI coverage
    for Tsallis entropy (bounded distribution)
4.  **No guidance on shrinkage**: Users with sparse genes have inflated
    noise, unaware of shrinkage option

### For Maintainers:

1.  Parameter documentation scattered across 3+ files (tsenat_config,
    calculate_diversity_s4, .calculate_diversity)
2.  Inconsistent parameter documentation across functions
3.  Hard to maintain consistency when updating parameter defaults

### For reproducibility:

1.  Users unsure which parameters were used in previous analysis (not
    recorded in config documentation)
2.  Different users make different parameter choices → non-reproducible
    results without explicit config recording

------------------------------------------------------------------------

## FINAL RECOMMENDATION

**Immediate Actions** (before next release): 1. Add Tier 1 parameters to
tsenat_config() signature and help (7 parameters) 2. Add @examples
section to tsenat_config() showing these parameters in use 3. Update
DESCRIPTION/NEWS.md noting parameter additions

**Medium-term** (next release cycle): 1. Create “Advanced Parameters”
vignette explaining all parameters 2. Add parameter audit table to
documentation 3. Refactor parameter resolution to single unified
function

**Long-term** (future releases): 1. Consider YAML config file support
(for reproducibility) 2. Add parameter validation schema 3. Create
web-based configuration builder tool
