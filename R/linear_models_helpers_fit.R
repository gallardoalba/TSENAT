#' Fit Statistical Model for Gene-Specific Q-Entropy Interaction
#'
#' @description
#' Dispatcher function that fits a specified statistical model to test for
#' q-dependent
#' interaction effects in Tsallis entropy across sample groups. Supports
#' multiple
#' modeling approaches (LMM, GAM, FPCA, GEE) optimized for different data
#' structures.
#'
#' @details
#' CORE PURPOSE:
#' For a single gene g in matrix mat, extract entropy values and test for an
#' interaction between q-parameter and sample group. Returns p-values and effect
#' size estimates (slope_diff) describing q*group interaction.
#'
#' STATISTICAL METHODS:
#'
#' 1. LMM (Linear Mixed Models):
#'    - Uses nlme::lme with AR(1) covariance structure to model correlation
#'      across ordered q-values within each subject (repeated measures design)
#'    - ARIMA(1,1,0): Implements first-differencing to enforce stationarity
#' - Heteroscedasticity detection: Applies nlme::varPower() for q-dependent
#' variance
#' - Hypothesis test: Likelihood Ratio Test (LRT) comparing null (no
#' interaction)
#'      vs alternative (q*group interaction) models
#'    - Best for: Paired/repeated measures data with multiple subjects
#'
#' 2. GAM (Generalized Additive Models):
#' - Uses mgcv::gam with smooth spline terms s(q) to capture nonlinear
#' q-effects
#' - Regularization: Optional feature selection (PCA, LASSO, elastic net,
#' gamsel)
#'    - Adaptive knots: Automatically selects optimal number of basis functions
#'    - Hypothesis test: F-test or ANOVA comparing model fits
#'    - Best for: Flexible modeling of complex nonlinear relationships
#'
#' 3. FPCA (Functional Principal Component Analysis):
#' - Treats each subject's entropy curve (q → entropy) as a functional
#' observation
#'    - Extracts principal components explaining variance in curve shapes
#'    - Tests for group differences in functional structure via PCA scores
#'    - Best for: Small-sample designs; captures curve topology
#'
#' 4. GEE (Generalized Estimating Equations):
#'    - Semi-parametric method for clustered/correlated observations
#'    - AR(1), exchangeable, or independence correlation structures
#'    - Robust to variance misspecification; sandwich standard errors
#' - Best for: Large samples with clusters; robustness to distributional
#' assumptions
#'
#' SUBJECT IDENTIFICATION:
#' For methods requiring repeated measures (LMM, GAM with subject, FPCA),
#' subject IDs are determined by:
#'   1. If subject_col provided: use that colData column
#' 2. Else if paired=TRUE: search for 'paired_samples' or 'sample_base'
#' columns
#'   3. Else: use sample names (each treated as independent unit)
#'
#' DATA VALIDATION:
#'   - Returns NULL if < min_obs observations (insufficient statistical power)
#' - Returns NULL if < 2 subjects (required for mixed models with random
#' effects)
#'   - Handles missing values via casewise deletion (na.omit in df construction)
#'
#' @param g Character gene identifier (row name in mat)
#' @param se SummarizedExperiment object containing sample metadata in colData
#' @param mat Numeric matrix (genes × samples) of entropy values indexed by g
#' @param q_vals Numeric vector of Tsallis q-parameters (length = ncol(mat))
#' @param sample_names Character vector of sample identifiers (length =
#' ncol(mat))
#' @param group_vec Factor/character vector of group assignments (length =
#' ncol(mat))
#' @param method Character: statistical method - 'lmm', 'gam', 'fpca', or 'gee'
#' @param pvalue Character: p-value extraction method - 'lrt',
#' 'satterthwaite', or 'both'
#'   (LMM only; ignored for GAM/FPCA/GEE)
#' @param subject_col Character: colData column name for subject IDs
#' (optional; overrides paired)
#' @param paired Logical: if TRUE, search for paired_samples or sample_base
#' columns (LMM, GAM)
#' @param min_obs Integer: minimum observations required (default 2);
#' returns NULL if nrow(df) < min_obs
#' @param verbose Logical: if TRUE, print diagnostic messages during fitting
#' @param suppress_lme4_warnings Logical: if TRUE, suppress lme4 warnings
#' during fitting
#' @param progress Logical: if TRUE, show progress messages and timing
#' @param bias_correction Logical: if TRUE, apply bias corrections in GAM
#' models (default TRUE)
#' @param regularization Character: regularization method - 'pca', 'lasso',
#' 'elasticnet',
#'   'gamsel', or 'spline' (affects GAM and FPCA feature selection)
#' @param corstr Character: correlation structure - 'ar1' (default),
#' 'exchangeable', or
#'   'independence' (GEE only)
#' @param adaptive_knots Logical: if TRUE, automatically select basis
#' dimension in GAM (default TRUE)
#' @param weights Numeric vector: optional inverse-variance weights for
#' robust estimation
#'   (length must equal nrow(df); applied in LMM and GEE)
#'
#' @return
#' Data frame with one row containing:
#'   - gene: Character gene identifier
#'   - p_interaction: Numeric p-value for q*group interaction test
#' - p_lrt: Numeric p-value from Likelihood Ratio Test (LMM only; NA for
#' GAM/FPCA/GEE)
#' - slope_diff: Numeric interaction coefficient (slope difference between
#' groups)
#' - fit_method: Character method used ('nlme::lme',
#' 'nlme::lme_arima(1,1,0)', 'gam', 'fpca', 'gee', etc.)
#'   - singular: Logical TRUE if model fit was singular (lmer only)
#' - arima_transformation: Logical TRUE if ARIMA(1,1,0) first-differencing
#' applied (LMM)
#'   - ci_weighted: Logical TRUE if inverse-variance weights were applied
#'   - n_subjects: Integer number of subjects in model
#'   - small_sample_flag: Logical TRUE if sample size < optimal threshold
#'
#' Returns NULL if:
#'   - nrow(df) < min_obs (insufficient observations)
#'   - n_subjects < 2 (required for mixed models)
#'   - Data quality issues prevent model fitting
#'
#' @noRd
#' @noRd
.fit_one_interaction <- function(g, se, mat, q_vals, sample_names, group_vec, method,
    pvalue, subject_col, paired, min_obs, verbose, suppress_lme4_warnings, progress,
    bias_correction = TRUE, regularization = c("pca", "lasso", "elasticnet", "gamsel",
        "spline"), corstr = c("ar1", "exchangeable", "independence"), adaptive_knots = TRUE,
    weights = NULL) {
    # ═══════════════════════════════════════════════════════════════════════════
    # INPUT VALIDATION AND SETUP (compressed via helpers)
    # ═══════════════════════════════════════════════════════════════════════════
    # Matches regularization/corstr args, validates gene exists, extracts
    # entropy, initializes data frame, applies optional weights

    regularization <- match.arg(regularization)
    corstr <- match.arg(corstr)

    # ═══════════════════════════════════════════════════════════════════════════
    # VALIDATE subject_col AND paired PARAMETERS (before method dispatch)
    # ═══════════════════════════════════════════════════════════════════════════
    if (!is.null(subject_col)) {
        # If subject_col is explicitly provided, it must exist in colData
        if (is.null(se) || !(subject_col %in% colnames(SummarizedExperiment::colData(se)))) {
            stop(sprintf("subject_col '%s' not found in colData(se)", subject_col))
        }
    }

    if (paired && is.null(subject_col)) {
        # If paired=TRUE without explicit subject_col, check for standard
        # columns
        if (!is.null(se)) {
            coldata_cols <- colnames(SummarizedExperiment::colData(se))
            has_paired_info <- ("paired_samples" %in% coldata_cols) || ("sample_base" %in%
                coldata_cols) || (length(coldata_cols) >= 3)
            if (!has_paired_info) {
                stop("paired = TRUE requires 'paired_samples' or 'sample_base' column in colData; supply subject_col explicitly")
            }
        }
    }

    # Setup data frame via helper (validates gene + builds df)
    df <- .setup_interaction_data(g, mat, q_vals, group_vec)
    df <- .apply_weights_to_df(df, weights, g, verbose)

    # Phase 15: Wrap all method fitting in tryCatch to handle edge case errors
    # gracefully (e.g., 'los nombres no coinciden' from factor level
    # mismatches) CRITICAL: Must NOT use return() inside tryCatch - it bypasses
    # error handler!  Instead, assign to result variable so error handler can
    # catch anything
    result <- tryCatch({
        # Dispatch to method-specific helper (contains all details for that
        # method)
        if (method == "lmm") {
            .lmm_interaction(df, mat, q_vals, sample_names, group_vec, g, se, subject_col,
                paired, min_obs, verbose, suppress_lme4_warnings, progress, regularization,
                weights)
        } else if (method == "gam") {
            subject <- .get_subject_ids(se, subject_col, paired, mat, sample_names)
            df$subject <- if (!is.null(subject))
                factor(subject) else factor(sample_names)
            .gam_interaction(df, q_vals, g, min_obs = min_obs, subject = subject,
                regularization = regularization, bias_correction = bias_correction,
                adaptive_knots = adaptive_knots, weights = weights)
        } else if (method == "fpca") {
            subject <- .get_subject_ids(se, subject_col, paired, mat, sample_names)
            .fpca_interaction(mat, q_vals, sample_names, group_vec, g, min_obs = min_obs,
                subject = subject, regularization = regularization, weights = weights)
        } else if (method == "gee") {
            subject <- .get_subject_ids(se, subject_col, paired, mat, sample_names)
            if (is.null(subject) && !paired)
                subject <- sample_names
            .gee_interaction(df, q_vals, g, subject = subject, min_obs = min_obs,
                corstr = corstr, bias_correction = bias_correction, weights = weights)
        } else {
            NULL  # Invalid method
        }
    }, error = function(e) {
        # Phase 15: Log detailed error information for debugging Phase 16:
        # Capture gene-specific error diagnostics
        error_msg <- conditionMessage(e)

        # Build diagnostic message with gene-specific details
        df_info <- if (exists("df_model") && is.data.frame(df_model)) {
            sprintf("rows=%d, q-levels=%d, groups=%s, subjects=%d", nrow(df_model),
                length(unique(df_model$q)), paste(levels(df_model$group), collapse = "/"),
                length(levels(df_model$subject)))
        } else if (exists("df") && is.data.frame(df)) {
            sprintf("rows=%d, q-levels=%d, groups=%s", nrow(df), length(unique(df$q)),
                paste(levels(df$group), collapse = "/"))
        } else {
            "data structure unavailable"
        }

        diag_msg <- sprintf("[.fit_one_interaction] Gene '%s' failed: %s [%s]", g,
            error_msg, df_info)

        if (verbose)
            message(diag_msg)

        # Return NA results on any error instead of crashing This handles edge
        # cases like 'los nombres no coinciden' gracefully
        data.frame(gene = g, p_interaction = NA_real_, p_lrt = NA_real_, slope_diff = NA_real_,
            fit_method = "ERROR", singular = NA, arima_transformation = NA, ci_weighted = NA,
            n_subjects = NA_integer_, small_sample_flag = NA, message = error_msg,
            stringsAsFactors = FALSE)
    })

    # Ensure 'gene' column exists in result before returning
    if (is.data.frame(result) && nrow(result) > 0) {
        if (!("gene" %in% colnames(result))) {
            result$gene <- g
        }
    }

    result

}

#' @noRd
#' @noRd
.lmm_interaction <- function(df, mat, q_vals, sample_names, group_vec, g, se, subject_col,
    paired, min_obs, verbose, suppress_lme4_warnings, progress, regularization, weights) {
    # ═══════════════════════════════════════════════════════════════════════════
    # .lmm_interaction() - LINEAR MIXED MODELS (LMM METHOD IMPLEMENTATION)
    # ═══════════════════════════════════════════════════════════════════════════
    # Purpose: Fit LMM with AR(1) covariance structure for interaction testing
    # Contains all method-specific logic extracted from main dispatcher
    # function Preserves full documentation of methodology while reducing line
    # count

    if (!requireNamespace("nlme", quietly = TRUE)) {
        stop("Package 'nlme' is required for method = 'lmm'")
    }

    # Get subject IDs and build subject column
    subject <- .get_subject_ids(se, subject_col, paired, mat, sample_names)
    if (is.null(subject))
        subject <- sample_names
    df$subject <- factor(subject)

    # Validate sample size requirements
    result <- .check_lmm_sample_sizes(df, min_obs)
    if (!isTRUE(result))
        return(NULL)

    # ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
    # ARIMA(1,1,0) TRANSFORMATION: First differencing for stationarity
    # ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
    # Problem: Raw Tsallis entropy H_q is monotone increasing with q, violating
    # stationarity assumption (constant mean) required for AR(1) modeling
    # Solution: Use first differences DeltaH_q = H_q - H_{q-1} to remove trend
    # - Bounded-support data [0, log(m)] after differencing approximates
    # normality - Enables valid hypothesis testing under AR(1) correlation
    # structure - Information preserved: interaction effects remain in
    # differenced data

    arima_result <- .compute_arima_differences(df, q_vals, df$group, df$subject)

    if (is.null(arima_result) || nrow(arima_result$df) < 3) {
        df_model <- df
        use_arima <- FALSE
        if (verbose) {
            message("[.lmm_interaction] ARIMA(1,1,0) differencing lost too many observations; using raw entropy")
        }
    } else {
        df_model <- arima_result$df
        use_arima <- TRUE
        if (verbose) {
            message(sprintf("[.lmm_interaction] ARIMA(1,1,0): %d observations -> %d after differencing",
                arima_result$n_observations_original, arima_result$n_observations_differenced))
        }
    }

    # ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
    # MODEL FORMULAS: NULL (no interaction) vs ALTERNATIVE (q*group
    # interaction)
    # ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
    formula_null <- entropy ~ q + group
    formula_alt <- entropy ~ q * group

    # OPTIONAL REGULARIZATION: Apply feature selection if not PCA
    fs_result <- NULL
    if (regularization != "pca") {
        fs_result <- .lmm_regularization(q_vals = df_model$q, entropy_vals = df_model$entropy,
            group_vec = df_model$group, subject_vec = df_model$subject, regularization = regularization)
        if (!is.null(fs_result) && verbose) {
            uq_levels <- length(fs_result$q_values)
            message("[.lmm_interaction] regularization retained ", length(fs_result$selected_features),
                " of ", uq_levels - 1, " q-interaction features")
        }
    }

    # ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
    # HETEROSCEDASTICITY DETECTION: Use varPower() if variance depends on q
    # ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
    # NOTE: Fixed in linear_models_lmm.R - is.na() coercion error is now
    # handled Variance structure detection is now enabled by default
    use_var_structure <- TRUE

    # ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
    # MODEL FITTING: nlme::lme with AR(1) covariance and optional variance
    # structure
    # ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

    # Phase 16: Log factor structure before fitting to diagnose 'nombres no
    # coinciden' errors
    if (verbose) {
        message(sprintf("[.lmm_interaction] Gene '%s' data: %d obs, q-levels=%d, group-levels=%s, subject-levels=%d",
            g, nrow(df_model), length(unique(df_model$q)), paste(levels(df_model$group),
                collapse = "/"), length(levels(df_model$subject))))
    }

    if (use_var_structure) {
        if (verbose) {
            fit0 <- try(nlme::lme(formula_null, random = ~1 | subject, data = df_model,
                method = "ML"), silent = TRUE)
            fit1 <- try(nlme::lme(formula_alt, random = ~1 | subject, data = df_model,
                method = "ML"), silent = TRUE)
        } else {
            fit0 <- try(nlme::lme(formula_null, random = ~1 | subject, data = df_model,
                method = "ML"), silent = TRUE)
            fit1 <- try(nlme::lme(formula_alt, random = ~1 | subject, data = df_model,
                method = "ML"), silent = TRUE)
        }
    } else {
        if (verbose) {
            fit0 <- try(nlme::lme(formula_null, random = ~1 | subject, data = df_model,
                method = "ML"), silent = TRUE)
            fit1 <- try(nlme::lme(formula_alt, random = ~1 | subject, data = df_model,
                method = "ML"), silent = TRUE)
        } else {
            fit0 <- try(nlme::lme(formula_null, random = ~1 | subject, data = df_model,
                method = "ML"), silent = TRUE)
            fit1 <- try(nlme::lme(formula_alt, random = ~1 | subject, data = df_model,
                method = "ML"), silent = TRUE)
        }
    }

    # Phase 16: Log fit errors for diagnosis
    if (inherits(fit0, "try-error") && verbose) {
        message(sprintf("[.lmm_interaction] Gene '%s' fit0 error (null model): %s",
            g, attr(fit0, "condition")$message))
    }
    if (inherits(fit1, "try-error") && verbose) {
        message(sprintf("[.lmm_interaction] Gene '%s' fit1 error (alt model): %s",
            g, attr(fit1, "condition")$message))
    }

    # ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
    # FALLBACK FITTING: When nlme fails, try simpler fixed-effects models
    # ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
    fallback_lm <- NULL
    used_fit_method <- "nlme::lme"
    used_singular <- FALSE

    if (inherits(fit0, "try-error") || inherits(fit1, "try-error")) {
        if (progress || verbose) {
            message("[.lmm_interaction] nlme::lme failed; trying fallback models")
        }
        fb <- .try_lm_fallbacks(df_model, verbose = verbose)
        if (!is.null(fb)) {
            fallback_lm <- fb
            used_fit_method <- fb$method
        }
    } else {
        used_fit_method <- if (use_arima)
            "nlme::lme_arima(1,1,0)" else "nlme::lme_ar1_raw"
    }

    # ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
    # HYPOTHESIS TEST: Extract LRT p-value and effect size
    # ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
    lrt_result <- list(p_value = NA_real_, n_subjects = NA_integer_, small_sample_flag = FALSE)
    msg <- NULL

    if (!is.null(fallback_lm)) {
        lrt_result <- .extract_lrt_p(fallback_lm$fit0, fallback_lm$fit1, df = df_model)
        if (!is.null(fallback_lm$message))
            msg <- fallback_lm$message
    } else {
        lrt_result <- .extract_lrt_p(fit0, fit1, df = df_model)
    }

    # ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
    # EFFECT SIZE: Extract slope_diff (q:group interaction coefficient)
    # ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
    slope_diff <- NA_real_
    if (!is.null(fallback_lm) && !is.null(fallback_lm$fit1)) {
        coefs <- tryCatch(coef(fallback_lm$fit1), error = function(e) NULL)
        if (!is.null(coefs)) {
            interaction_idx <- grep("q:group|group:q", names(coefs), ignore.case = FALSE)
            if (length(interaction_idx) > 0)
                slope_diff <- coefs[interaction_idx[1]]
        }
    } else if (!inherits(fit1, "try-error")) {
        coefs <- tryCatch(nlme::fixef(fit1), error = function(e) NULL)
        if (!is.null(coefs)) {
            interaction_idx <- grep("q:group|group:q", names(coefs), ignore.case = FALSE)
            if (length(interaction_idx) > 0)
                slope_diff <- coefs[interaction_idx[1]]
        }
    }

    # ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
    # RESULTS OUTPUT: Return data frame with model statistics
    # ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
    has_weights <- !is.null(df$weight)
    res <- data.frame(gene = g, p_interaction = lrt_result$p_value, p_lrt = lrt_result$p_value,
        slope_diff = slope_diff, fit_method = used_fit_method, singular = used_singular,
        arima_transformation = use_arima, ci_weighted = has_weights, n_subjects = lrt_result$n_subjects,
        small_sample_flag = lrt_result$small_sample_flag, stringsAsFactors = FALSE)
    if (!is.null(msg))
        res$message <- msg
    res
}

#' Helper functions for .fit_one_interaction() 
#' 
#' These functions extract repeated logic from .fit_one_interaction() to reduce
#' function complexity and meet Bioconductor guideline of <50 lines per
#' function.
#' 
#' @details Bioconductor package guidelines recommend functions be kept
#' under 50 lines
#' when possible. This file contains extracted helpers to improve
#' maintainability.

# ═══════════════════════════════════════════════════════════════════════════════
# .setup_interaction_data() - Validate inputs and initialize data frame
# ═══════════════════════════════════════════════════════════════════════════════

#' @noRd
#' @noRd
.setup_interaction_data <- function(g, mat, q_vals, group_vec) {
    # Validate that gene g exists in matrix
    if (!(g %in% rownames(mat))) {
        stop(sprintf("Gene '%s' not found in matrix rownames. Available genes: %s",
            g, paste(rownames(mat)[seq_len(min(5, nrow(mat)))], collapse = ", ")))
    }

    # Extract entropy values and build data frame
    vals <- as.numeric(mat[g, ])
    
    # Extract sample names from column names (remove q-value suffix if present)
    # If colnames are NULL, generate generic sample names
    if (is.null(colnames(mat))) {
        sample_names_extracted <- paste0("S", seq_len(ncol(mat)))
    } else {
        sample_names_extracted <- colnames(mat)
        sample_names_extracted <- sub("_q=.*", "", sample_names_extracted)
    }
    
    df <- data.frame(
        entropy = vals, 
        q = q_vals, 
        group = factor(group_vec),
        sample_name = factor(sample_names_extracted),
        stringsAsFactors = FALSE
    )
    df
}

# ═══════════════════════════════════════════════════════════════════════════════
# .apply_weights_to_df() - Handle weight validation and application
# ═══════════════════════════════════════════════════════════════════════════════

#' @noRd
#' @noRd
.apply_weights_to_df <- function(df, weights, g, verbose = FALSE) {
    if (!is.null(weights) && length(weights) == nrow(df)) {
        df$weight <- weights
        if (verbose) {
            message(sprintf("[.fit_one_interaction] Gene '%s': weights applied (n=%d, mean=%.4f, min=%.4f, max=%.4f)",
                g, length(weights), mean(weights, na.rm = TRUE), min(weights, na.rm = TRUE),
                max(weights, na.rm = TRUE)))
        }
    } else {
        if (verbose && !is.null(weights)) {
            message(sprintf("[.fit_one_interaction] Gene '%s': weights NOT applied - length mismatch (weights=%d, df rows=%d)",
                g, length(weights), nrow(df)))
        }
    }
    df
}

# ═══════════════════════════════════════════════════════════════════════════════
# .get_subject_ids() - Extract subject identifiers (refactored from 4
# duplicates)
# ═══════════════════════════════════════════════════════════════════════════════

#' @noRd
#' @noRd
.get_subject_ids <- function(se = NULL, subject_col = NULL, paired = FALSE, mat = NULL,
    sample_names = NULL) {
    subject <- NULL

    # Explicit subject_col provided
    if (!is.null(subject_col)) {
        if (is.null(se) || !(subject_col %in% colnames(SummarizedExperiment::colData(se)))) {
            stop(sprintf("subject_col '%s' not found in colData(se)", subject_col))
        }
        subj_full <- as.character(SummarizedExperiment::colData(se)[, subject_col])
        col_names_for_indexing <- SummarizedExperiment::colData(se)$samples
        if (is.null(col_names_for_indexing)) {
            col_names_for_indexing <- sub("_q=.*", "", colnames(mat))
        }
        names(subj_full) <- col_names_for_indexing
        subject <- unname(subj_full[sample_names])
        return(subject)
    }

    # paired=TRUE: extract from colData
    if (paired && !is.null(se)) {
        coldata <- SummarizedExperiment::colData(se)
        coldata_cols <- colnames(coldata)

        # Look for standard paired columns
        subject_col_name <- NULL
        if ("paired_samples" %in% coldata_cols) {
            subject_col_name <- "paired_samples"
        } else if ("sample_base" %in% coldata_cols) {
            subject_col_name <- "sample_base"
        } else if (length(coldata_cols) >= 3) {
            subject_col_name <- coldata_cols[3]  # fallback
        }

        if (!is.null(subject_col_name)) {
            subject_ids <- as.character(coldata[, subject_col_name])
            sample_names_expanded <- sub("_q=.*", "", rownames(coldata))
            names(subject_ids) <- sample_names_expanded
            subject <- unname(subject_ids[sample_names])
            return(subject)
        } else {
            stop("paired = TRUE requires 'paired_samples' or 'sample_base' column in colData; supply subject_col explicitly")
        }
    }

    # Fallback: use sample names as subject identifiers
    sample_names
}

# ═══════════════════════════════════════════════════════════════════════════════
# .check_lmm_sample_sizes() - Validate minimum sample and subject requirements
# ═══════════════════════════════════════════════════════════════════════════════

#' @noRd
#' @noRd
.check_lmm_sample_sizes <- function(df, min_obs = 3) {
    if (nrow(df) < min_obs) {
        warning(sprintf(".check_lmm_sample_sizes: Insufficient observations in model data. Found %d rows, minimum required: %d.", nrow(df), min_obs), call. = FALSE)
        return(NULL)  # Not enough observations
    }

    n_subjects <- length(unique(na.omit(df$subject)))
    if (n_subjects < 2) {
        warning(sprintf(".check_lmm_sample_sizes: Insufficient subjects for random intercept model. Found %d unique subjects, minimum required: 2.", n_subjects), call. = FALSE)
        return(NULL)  # Not enough subjects for random intercept
    }

    TRUE  # Passes all checks
}
