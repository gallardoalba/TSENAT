# Code Implementation Guide

## Adding Missing Parameters to tsenat_config() and Orchestration

------------------------------------------------------------------------

## FILE 1: R/orchestration.R - tsenat_config() function

### Current State (Lines 485-550)

``` r
tsenat_config <- function(q_values = NULL, condition_col = "condition", subject_col = NULL,
    sample_col = "sample", paired = FALSE, control = NULL, p_threshold = 0.05, fdr_threshold = 0.05,
    significance_threshold = 0.05, bootstrap = FALSE, nboot = 1000, bootstrap_method = "percentile",
    stringency = "medium", nthreads = 1, norm = TRUE, 
    bootstrap_ci = 0.95, bootstrap_include_diagnostics = TRUE, min_valid_frac = 0.75,
    norm_method = NULL, pseudocount = NULL, shrinkage = "none", ...)
```

### Recommended Changes

#### Step 1: Add new parameters to function signature

``` r
tsenat_config <- function(
    q_values = NULL, 
    condition_col = "condition", 
    subject_col = NULL,
    sample_col = "sample", 
    paired = FALSE, 
    control = NULL, 
    p_threshold = 0.05, 
    fdr_threshold = 0.05,
    significance_threshold = 0.05, 
    bootstrap = FALSE, 
    nboot = 1000, 
    bootstrap_method = "percentile",
    stringency = "medium", 
    nthreads = 1, 
    norm = TRUE, 
    bootstrap_ci = 0.95, 
    bootstrap_include_diagnostics = TRUE, 
    min_valid_frac = 0.75,
    norm_method = NULL, 
    pseudocount = NULL, 
    shrinkage = "none",
    
    # === NEW PARAMETERS (CRITICAL) ===
    lm_method = "gam",                    # NEW: LM interaction method
    lm_pcorr = "BH",                      # NEW: P-value correction
    divergence_ci = 0.95,                 # NEW: Divergence CI width
    divergence_control_group = NULL,      # NEW: Divergence reference group
    jis_use_lm_fdr = TRUE,               # NEW: JIS FDR usage
    
    ...)
```

#### Step 2: Add validation for new parameters

``` r

    # Validate bootstrap_method
    valid_bootstrap_methods <- c("percentile", "bca")
    if (!bootstrap_method %in% valid_bootstrap_methods) {
        stop("'bootstrap_method' must be 'percentile' or 'bca'", call. = FALSE)
    }
    
    # === NEW VALIDATION ===
    # Validate LM method
    valid_lm_methods <- c("lmm", "gam", "fpca", "gee")
    if (!lm_method %in% valid_lm_methods) {
        stop("'lm_method' must be one of: ", paste(valid_lm_methods, collapse = ", "), 
             call. = FALSE)
    }
    
    # Validate LM p-correction method
    valid_pcorr_methods <- c("BH", "bonferroni", "hochberg", "holm")
    if (!lm_pcorr %in% valid_pcorr_methods) {
        stop("'lm_pcorr' must be one of: ", paste(valid_pcorr_methods, collapse = ", "), 
             call. = FALSE)
    }
    
    # Validate divergence CI
    if (!is.numeric(divergence_ci) || divergence_ci <= 0 || divergence_ci >= 1) {
        stop("'divergence_ci' must be a numeric value between 0 and 1", call. = FALSE)
    }
    
    # Validate JIS use_lm_fdr
    if (!is.logical(jis_use_lm_fdr) || is.na(jis_use_lm_fdr)) {
        stop("'jis_use_lm_fdr' must be a logical value (TRUE or FALSE)", call. = FALSE)
    }
```

#### Step 3: Add parameters to config list

``` r

    # Build config list with all parameters
    config <- list(
        q_values = q_values,
        condition_col = condition_col,
        subject_col = subject_col,
        sample_col = sample_col,
        paired = paired,
        control = control,
        p_threshold = p_threshold,
        fdr_threshold = fdr_threshold,
        significance_threshold = significance_threshold,
        nboot = nboot,
        bootstrap_method = bootstrap_method,
        bootstrap = bootstrap,
        bootstrap_ci = bootstrap_ci,
        bootstrap_include_diagnostics = bootstrap_include_diagnostics,
        min_valid_frac = min_valid_frac,
        stringency = stringency,
        nthreads = nthreads,
        norm = norm,
        norm_method = norm_method,
        pseudocount = pseudocount,
        shrinkage = shrinkage,
        
        # === NEW PARAMETERS ===
        lm_method = lm_method,
        lm_pcorr = lm_pcorr,
        divergence_ci = divergence_ci,
        divergence_control_group = divergence_control_group,
        jis_use_lm_fdr = jis_use_lm_fdr
    )
```

------------------------------------------------------------------------

## FILE 2: R/orchestration.R - .execute_lm_interaction_s4() function

### Current State (Lines 725-740)

``` r

.execute_lm_interaction_s4 <- function(analysis, verbose, output_dir, output_format) {
    if (verbose)
        message(sprintf("[>] [%2d/14] Testing LM interactions with GAM smoother", 5))
    tryCatch({
        cfg <- getConfig(analysis)
        fdr <- cfg$fdr_threshold %||% 0.05
        output_file <- .build_output_file("lm_interaction_results", output_dir, output_format)
        analysis <- calculate_lm_interaction_s4(analysis, fdr_threshold = fdr, output_file = output_file)
        if (verbose)
            message("          [OK] LM interaction analysis complete")
    }, error = function(e) warning("LM interaction analysis failed:\n", e$message, call. = FALSE))
    analysis
}
```

### Recommended Changes

``` r

.execute_lm_interaction_s4 <- function(analysis, verbose, output_dir, output_format) {
    if (verbose)
        message(sprintf("[>] [%2d/14] Testing LM interactions", 5))
    tryCatch({
        cfg <- getConfig(analysis)
        fdr <- cfg$fdr_threshold %||% 0.05
        lm_method <- cfg$lm_method %||% "gam"              # NEW: Extract method
        lm_pcorr <- cfg$lm_pcorr %||% "BH"                 # NEW: Extract pcorr
        output_file <- .build_output_file("lm_interaction_results", output_dir, output_format)
        
        # NEW: Pass method and pcorr to function
        analysis <- calculate_lm_interaction_s4(
            analysis, 
            fdr_threshold = fdr, 
            method = lm_method,                            # NEW
            pcorr = lm_pcorr,                              # NEW
            output_file = output_file
        )
        if (verbose)
            message(sprintf("          [OK] LM interaction analysis complete (method=%s, pcorr=%s)", 
                           lm_method, lm_pcorr))           # NEW: Enhanced message
    }, error = function(e) warning("LM interaction analysis failed:\n", e$message, call. = FALSE))
    analysis
}
```

------------------------------------------------------------------------

## FILE 3: R/orchestration.R - .execute_divergence_s4() function

### Current State (Lines 842-855)

``` r

.execute_divergence_s4 <- function(analysis, q_vals, verbose, output_dir, output_format) {
    if (verbose)
        message(sprintf("[>] [%2d/14] Computing divergence metrics", 11))
    tryCatch({
        output_file <- .build_output_file("divergence_results", output_dir, output_format)
        analysis <- calculate_divergence_s4(analysis, q = q_vals, output_file = output_file, verbose = FALSE)
        if (verbose)
            message("          [OK] Divergence computed")
    }, error = function(e) {
        if (verbose)
            warning("Divergence failed: ", e$message, call. = FALSE)
    })
    analysis
}
```

### Recommended Changes

``` r

.execute_divergence_s4 <- function(analysis, q_vals, verbose, output_dir, output_format) {
    if (verbose)
        message(sprintf("[>] [%2d/14] Computing divergence metrics", 11))
    tryCatch({
        cfg <- getConfig(analysis)
        output_file <- .build_output_file("divergence_results", output_dir, output_format)
        
        # NEW: Extract divergence parameters from config
        div_ci <- cfg$divergence_ci %||% 0.95
        div_method <- cfg$bootstrap_method %||% "percentile"
        div_control <- cfg$divergence_control_group %||% NULL
        
        analysis <- calculate_divergence_s4(
            analysis, 
            q = q_vals, 
            output_file = output_file, 
            verbose = FALSE,
            method = div_method,                           # NEW: Pass CI method
            ci = div_ci,                                   # NEW: Pass CI level
            control_group = div_control                    # NEW: Pass control group
        )
        if (verbose)
            message(sprintf("          [OK] Divergence computed (ci=%s, method=%s)", 
                           div_ci, div_method))            # NEW: Enhanced message
    }, error = function(e) {
        if (verbose)
            warning("Divergence failed: ", e$message, call. = FALSE)
    })
    analysis
}
```

------------------------------------------------------------------------

## FILE 4: R/orchestration.R - .execute_jackknife_isoform_switching() function

### Current State (Lines 761-780)

``` r

.execute_jackknife_isoform_switching <- function(analysis, q_vals, condition_col,
    verbose, output_dir, output_format) {
    if (verbose)
        message(sprintf("[>] [%2d/14] Computing jackknife isoform switching analysis", 7))
    tryCatch({
        cfg <- getConfig(analysis)
        output_file <- .build_output_file("jackknife_isoform_switching", output_dir, output_format)
        analysis <- jackknife_isoform_switching_s4(analysis, condition_col = condition_col,
            output_file = output_file, verbose = FALSE)
        if (verbose)
            message("          [OK] Jackknife isoform switching complete")
    }, error = function(e) {
        if (verbose)
            warning("Jackknife isoform switching failed: ", e$message, call. = FALSE)
    })
    analysis
}
```

### Recommended Changes

``` r

.execute_jackknife_isoform_switching <- function(analysis, q_vals, condition_col,
    verbose, output_dir, output_format) {
    if (verbose)
        message(sprintf("[>] [%2d/14] Computing jackknife isoform switching analysis", 7))
    tryCatch({
        cfg <- getConfig(analysis)
        output_file <- .build_output_file("jackknife_isoform_switching", output_dir, output_format)
        
        # NEW: Extract JIS parameters from config
        jis_use_fdr <- cfg$jis_use_lm_fdr %||% TRUE
        
        analysis <- jackknife_isoform_switching_s4(
            analysis, 
            condition_col = condition_col,
            use_lm_fdr = jis_use_fdr,                      # NEW: Pass FDR usage flag
            output_file = output_file, 
            verbose = FALSE
        )
        if (verbose)
            message(sprintf("          [OK] Jackknife isoform switching complete (use_lm_fdr=%s)", 
                           jis_use_fdr))                   # NEW: Enhanced message
    }, error = function(e) {
        if (verbose)
            warning("Jackknife isoform switching failed: ", e$message, call. = FALSE)
    })
    analysis
}
```

------------------------------------------------------------------------

## FILE 5: R/s4_functions_divergence.R - Update function signature

### Current State (Line 107)

``` r
calculate_divergence_s4 <- function(analysis, q = NULL, verbose = FALSE, nthreads = NULL,
    output_file = NULL, control_group = NULL, paired = FALSE, method = NULL, bootstrap = FALSE,
    nboot = NULL, progress = FALSE, ...)
```

### Updates Needed

The function already accepts `ci`, `method` parameters but we need to
verify they are passed through properly:

**Verify that `ci` parameter is passed to .calculate_divergence():**

``` r

# In .build_divergence_args() function, ensure:
if (!is.null(params$ci)) {
    args$ci <- params$ci  # Already handles boot CI?
}
```

**Verify parameter resolution includes ci:**

``` r

# In .resolve_divergence_parameters(), add:
ci <- resolve_slot_param(ci, analysis@config, "divergence_ci", 0.95)

# Then return it
list(q = q, control_group = control_group, method = method, nthreads = nthreads,
    nboot = nboot, paired = paired, bootstrap = bootstrap, ci = ci)  # ADD ci
```

------------------------------------------------------------------------

## FILE 6: R/s4_functions_jis.R - Update function signature

### Current State (Line 173)

``` r
jackknife_isoform_switching_s4 <- function(analysis, condition_col = NULL, subject_col = NULL,
    gene_col = NULL, isoform_col = NULL, q = c(0, 0.5, 1, 1.5, 2), norm = NULL, log_base = NULL, threshold = 90,
    n_bootstrap = 1000, pseudocount = NULL, lm_results = NULL, lm_p_threshold = 0.05,
    use_lm_fdr = TRUE, output_file = NULL, verbose = FALSE, ...)
```

### Status

✅ **Good news:** The `use_lm_fdr` parameter is already in the
signature!

### Verification Needed

Ensure parameter resolution handles it:

``` r

# In parameter resolution section, confirm:
use_lm_fdr <- resolve_slot_param(use_lm_fdr, analysis@config, "jis_use_lm_fdr", TRUE)
```

------------------------------------------------------------------------

## FILE 7: Documentation Update - R/orchestration.R tsenat_config() roxygen

### Add to roxygen2 documentation

``` r

#' @param lm_method character. Linear model method for interaction testing.
#'   Default: "gam" (generalized additive model).
#'   Choices: "lmm" (linear mixed model), "gam", "fpca" (functional PCA), "gee" (generalized estimating equations).
#'   
#'   **Important:** Different methods produce different p-values and effect sizes.
#'   Reproducible analyses must document which method was used.
#'
#' @param lm_pcorr character. P-value correction method for LM interaction tests.
#'   Default: "BH" (Benjamini-Hochberg, controls FDR).
#'   Choices: "BH", "bonferroni" (controls FWER, most conservative), 
#'            "hochberg", "holm".
#'   
#'   **Important:** Different methods have different stringency. Must match
#'   between analyses for reproducibility.
#'
#' @param divergence_ci numeric. Confidence interval level for divergence calculations.
#'   Default: 0.95 (95% CI).
#'   Range: 0 to 1.
#'   
#'   **Important:** Currently hardcoded to 0.95 in core function; this parameter
#'   allows override without modifying code.
#'
#' @param divergence_control_group character or NULL. Reference group for pairwise
#'   divergence comparisons. Default: NULL (compares all groups).
#'   
#'   **Important:** Specifying a control group ensures all comparisons are
#'   relative to a baseline.
#'
#' @param jis_use_lm_fdr logical. Whether to use FDR-corrected (TRUE) or raw
#'   p-values (FALSE) from LM results when filtering in isoform switching analysis.
#'   Default: TRUE (use FDR-corrected).
#'   
#'   **Important:** Different choices dramatically affect number of genes selected
#'   for switching detection. Must match between analyses.
```

### Add to Details section

    **Bootstrap Method Recommendation for Bounded Data:**

    When analyzing Tsallis entropy (or other bounded distributions like proportions or probabilities):
    \itemize{
      \item \strong{Use} \code{bootstrap_method = "bca"} (bias-corrected and accelerated)
      \item \strong{NOT} \code{bootstrap_method = "percentile"}
      \item \strong{Reason:} Percentile bootstrap assumes symmetric distribution;
        entropy is inherently bounded (0 to log N) and skewed. BCA adjusts for this skewness
        and produces valid confidence intervals that contain the point estimate.
      \item \strong{Citation:} References C016 (2005), C030 (2023), S115 (2015)
    }

    **Method Selection Impact:**

    The choice of \code{lm_method} significantly affects results:
    \itemize{
      \item \code{lm_method = "gam"}: Smooth fit over q-spectrum; good for curved relationships
      \item \code{lm_method = "lmm"}: Mixed effects; better for complex hierarchical structure
      \item \code{lm_method = "fpca"}: Functional principal component analysis; good for functional data
      \item \code{lm_method = "gee"}: Generalized estimating equations; good for correlated observations
    }

    Different methods may produce different p-values and effect sizes. Document your choice
    in analysis reports and use the same method across comparable analyses.

------------------------------------------------------------------------

## TESTING CHECKLIST

When implementing these changes, test:

### ✓ Unit Tests

`test_tsenat_config_new_params.R` - Validates new parameters are
accepted

`test_tsenat_config_validation.R` - Validates parameter constraints

`test_lm_method_selection.R` - Ensures lm_method is passed through
correctly

`test_divergence_ci_levels.R` - Ensures divergence CI is respected

`test_jis_fdr_usage.R` - Ensures JIS respects use_lm_fdr flag

### ✓ Integration Tests

Run full pipeline with each `lm_method` value (“gam”, “lmm”, “fpca”,
“gee”)

Run full pipeline with each `lm_pcorr` value

Run full pipeline with `bootstrap_method = "bca"` (entropy case)

Run full pipeline with `bootstrap_method = "percentile"` (non-entropy
case)

Verify divergence CI changes with `divergence_ci` values (0.90, 0.95,
0.99)

Verify JIS filtering changes with `jis_use_lm_fdr = TRUE/FALSE`

### ✓ Reproducibility Tests

Save config from run A

Repeat analysis with identical config → results should match exactly

Change one parameter → results should differ predictably

Document which parameter changes affect which output

### ✓ Documentation Tests

roxygen2 build succeeds without errors

[`?tsenat_config`](https://gallardoalba.github.io/TSENAT/reference/tsenat_config.md)
shows all new parameters

Examples in documentation run without errors

Vignette examples work with new parameters

------------------------------------------------------------------------

## ESTIMATED IMPLEMENTATION TIME

| Task                                        | Time           |
|---------------------------------------------|----------------|
| Update tsenat_config() signature and list   | 15 min         |
| Add parameter validation                    | 20 min         |
| Update .execute\_\* functions (4 functions) | 30 min         |
| Add/update roxygen2 documentation           | 30 min         |
| Create unit tests                           | 45 min         |
| Integration testing                         | 60 min         |
| Documentation/examples                      | 30 min         |
| **TOTAL**                                   | **~3.5 hours** |

------------------------------------------------------------------------

**End of Implementation Guide**
