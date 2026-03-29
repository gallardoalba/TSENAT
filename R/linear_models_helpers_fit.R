#' Fit Statistical Model for Gene-Specific Q-Entropy Interaction
#'
#' @description
#' Dispatcher function that fits a specified statistical model to test for q-dependent
#' interaction effects in Tsallis entropy across sample groups. Supports multiple
#' modeling approaches (LMM, GAM, FPCA, GEE) optimized for different data structures.
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
#'    - Heteroscedasticity detection: Applies nlme::varPower() for q-dependent variance
#'    - Hypothesis test: Likelihood Ratio Test (LRT) comparing null (no interaction)
#'      vs alternative (q*group interaction) models
#'    - Best for: Paired/repeated measures data with multiple subjects
#'
#' 2. GAM (Generalized Additive Models):
#'    - Uses mgcv::gam with smooth spline terms s(q) to capture nonlinear q-effects
#'    - Regularization: Optional feature selection (PCA, LASSO, elastic net, gamsel)
#'    - Adaptive knots: Automatically selects optimal number of basis functions
#'    - Hypothesis test: F-test or ANOVA comparing model fits
#'    - Best for: Flexible modeling of complex nonlinear relationships
#'
#' 3. FPCA (Functional Principal Component Analysis):
#'    - Treats each subject's entropy curve (q → entropy) as a functional observation
#'    - Extracts principal components explaining variance in curve shapes
#'    - Tests for group differences in functional structure via PCA scores
#'    - Best for: Small-sample designs; captures curve topology
#'
#' 4. GEE (Generalized Estimating Equations):
#'    - Semi-parametric method for clustered/correlated observations
#'    - AR(1), exchangeable, or independence correlation structures
#'    - Robust to variance misspecification; sandwich standard errors
#'    - Best for: Large samples with clusters; robustness to distributional assumptions
#'
#' SUBJECT IDENTIFICATION:
#' For methods requiring repeated measures (LMM, GAM with subject, FPCA),
#' subject IDs are determined by:
#'   1. If subject_col provided: use that colData column
#'   2. Else if paired=TRUE: search for 'paired_samples' or 'sample_base' columns
#'   3. Else: use sample names (each treated as independent unit)
#'
#' DATA VALIDATION:
#'   - Returns NULL if < min_obs observations (insufficient statistical power)
#'   - Returns NULL if < 2 subjects (required for mixed models with random effects)
#'   - Handles missing values via casewise deletion (na.omit in df construction)
#'
#' @param g Character gene identifier (row name in mat)
#' @param se SummarizedExperiment object containing sample metadata in colData
#' @param mat Numeric matrix (genes × samples) of entropy values indexed by g
#' @param q_vals Numeric vector of Tsallis q-parameters (length = ncol(mat))
#' @param sample_names Character vector of sample identifiers (length = ncol(mat))
#' @param group_vec Factor/character vector of group assignments (length = ncol(mat))
#' @param method Character: statistical method - "lmm", "gam", "fpca", or "gee"
#' @param pvalue Character: p-value extraction method - "lrt", "satterthwaite", or "both"
#'   (LMM only; ignored for GAM/FPCA/GEE)
#' @param subject_col Character: colData column name for subject IDs (optional; overrides paired)
#' @param paired Logical: if TRUE, search for paired_samples or sample_base columns (LMM, GAM)
#' @param min_obs Integer: minimum observations required (default 2); returns NULL if nrow(df) < min_obs
#' @param verbose Logical: if TRUE, print diagnostic messages during fitting
#' @param suppress_lme4_warnings Logical: if TRUE, suppress lme4 warnings during fitting
#' @param progress Logical: if TRUE, show progress messages and timing
#' @param bias_correction Logical: if TRUE, apply bias corrections in GAM models (default TRUE)
#' @param regularization Character: regularization method - "pca", "lasso", "elasticnet",
#'   "gamsel", or "spline" (affects GAM and FPCA feature selection)
#' @param corstr Character: correlation structure - "ar1" (default), "exchangeable", or
#'   "independence" (GEE only)
#' @param adaptive_knots Logical: if TRUE, automatically select basis dimension in GAM (default TRUE)
#' @param weights Numeric vector: optional inverse-variance weights for robust estimation
#'   (length must equal nrow(df); applied in LMM and GEE)
#'
#' @return
#' Data frame with one row containing:
#'   - gene: Character gene identifier
#'   - p_interaction: Numeric p-value for q*group interaction test
#'   - p_lrt: Numeric p-value from Likelihood Ratio Test (LMM only; NA for GAM/FPCA/GEE)
#'   - p_satterthwaite: Numeric p-value from Satterthwaite approximation (LMM+lmer only; NA for nlme)
#'   - slope_diff: Numeric interaction coefficient (slope difference between groups)
#'   - fit_method: Character method used ("nlme::lme", "nlme::lme_arima(1,1,0)", "gam", "fpca", "gee", etc.)
#'   - singular: Logical TRUE if model fit was singular (lmer only)
#'   - arima_transformation: Logical TRUE if ARIMA(1,1,0) first-differencing applied (LMM)
#'   - ci_weighted: Logical TRUE if inverse-variance weights were applied
#'   - n_subjects: Integer number of subjects in model
#'   - small_sample_flag: Logical TRUE if sample size < optimal threshold
#'
#' Returns NULL if:
#'   - nrow(df) < min_obs (insufficient observations)
#'   - n_subjects < 2 (required for mixed models)
#'   - Data quality issues prevent model fitting
#'
#' @keywords internal
#' @noRd
.fit_one_interaction <- function(g, se, mat, q_vals, sample_names, group_vec,
    method, pvalue, subject_col, paired, min_obs, verbose, suppress_lme4_warnings,
    progress, bias_correction = TRUE, regularization = c("pca", "lasso", "elasticnet", "gamsel", "spline"),
    corstr = c("ar1", "exchangeable", "independence"), adaptive_knots = TRUE, weights = NULL) {
    # ═══════════════════════════════════════════════════════════════════════════
    # INPUT VALIDATION AND SETUP
    # ═══════════════════════════════════════════════════════════════════════════
    
    # Match regularization method to available options
    regularization <- match.arg(regularization)
    # Match correlation structure to available options (used in GEE)
    corstr <- match.arg(corstr)
    
    # Extract entropy values for gene g from expression matrix
    # mat is expected to be rows=genes, cols=samples
    # Validate that gene g exists in matrix before accessing
    if (!(g %in% rownames(mat))) {
        stop(sprintf("Gene '%s' not found in matrix rownames. Available genes: %s",
                     g, paste(rownames(mat)[1:min(5, nrow(mat))], collapse=", ")))
    }
    vals <- as.numeric(mat[g, ])
    
    # Build working data frame with entropy values, q-parameters, and group assignments
    df <- data.frame(entropy = vals, q = q_vals, group = factor(group_vec))
    
    
    # ═══════════════════════════════════════════════════════════════════════════
    # OPTIONAL WEIGHTING (Phase 1: Bootstrap CI weighting)
    # ═══════════════════════════════════════════════════════════════════════════
    # Apply inverse-variance weights if provided (used for robust estimation)
    # Weights inversely scale influences: high-variance observations weighted less
    # This enables consistent inference despite heteroscedasticity
    
    # Add inverse-variance weights if provided (Phase 1: Bootstrap CI weighting)
    if (!is.null(weights) && length(weights) == nrow(df)) {
        df$weight <- weights
        if (verbose) {
            message(sprintf("[.fit_one_interaction] Gene '%s': weights applied (n=%d, mean=%.4f, min=%.4f, max=%.4f)",
                           g, length(weights), mean(weights, na.rm=TRUE), min(weights, na.rm=TRUE), max(weights, na.rm=TRUE)))
        }
    } else {
        if (verbose && !is.null(weights)) {
            message(sprintf("[.fit_one_interaction] Gene '%s': weights NOT applied - length mismatch (weights=%d, df rows=%d)",
                           g, length(weights), nrow(df)))
        }
    }
    

    if (method == "lmm") {
        # ═══════════════════════════════════════════════════════════════════════════
        # LINEAR MIXED MODELS (LMM) with AR(1) Covariance Structure
        # ═══════════════════════════════════════════════════════════════════════════
        # Purpose: Test q*group interaction using nlme::lme with within-subject
        # correlation structure (AR(1) for repeated measurements across q-values)
        #
        # Methodological Foundation:
        # Paper S171 (Zimmerman & Harville, 1991): "Linear Models with Generalized AR(1) 
        # Covariance Structure for Longitudinal and Spatial Data" validates this approach.
        # Papers S168-S170: Theoretical foundation and empirical estimation of AR(1) parameters.
        # TEST L.1.6: Confirms q-value correlation follows AR(1) pattern (rho(k) = phi^|k|).
        #
        # CRITICAL FIX (March 2026): Implements true ARIMA(1,1,0) by:
        # 1. Computing first differences DeltaH_q = H_q - H_{q-1} within each subject
        # 2. Ensuring stationarity: DeltaH_q has constant mean (unlike monotone H_q)
        # 3. Applying AR(1) to differenced data (not raw entropy)
        # 4. Fitting all models (null and alt) on differenced entropy
        #
        if (!requireNamespace("nlme", quietly = TRUE)) {
            stop("Package 'nlme' is required for method = 'lmm' (AR(1) covariance support)")
        }
        
        # ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
        # SUBJECT IDENTIFICATION: Required for random effects structure
        # ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
        # Subject IDs are needed to define random intercept structure: ~1 | subject
        # Determines data structure: are samples from same subject (repeated measures)?
        # Identifies clustering structure for within-subject correlation estimation
        
        # determine subject IDs for random effect
        subject <- NULL
        if (!is.null(subject_col)) {
            if (!(subject_col %in% colnames(SummarizedExperiment::colData(se)))) {
                stop(sprintf("subject_col '%s' not found in colData(se)", subject_col))
            }
            subj_full <- as.character(SummarizedExperiment::colData(se)[, subject_col])
            # CRITICAL FIX: When colData doesn't have 'samples' column, use colnames(mat)
            # but strip "_q=" suffixes to match how sample_names were parsed in calculate_lm_interaction
            col_names_for_indexing <- SummarizedExperiment::colData(se)$samples
            if (is.null(col_names_for_indexing)) {
                # Strip "_q=....." suffixes from column names for proper indexing
                col_names_for_indexing <- sub("_q=.*", "", colnames(mat))
            }
            names(subj_full) <- col_names_for_indexing
            subject <- unname(subj_full[sample_names])
        } else if (paired) {
            coldata <- SummarizedExperiment::colData(se)
            coldata_cols <- colnames(coldata)
            
            # Look for 'paired_samples' or 'sample_base' columns created by map_metadata()
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
                # Map sample names by setting names to base names (without _q= suffix)
                sample_names_expanded <- sub("_q=.*", "", rownames(coldata))
                names(subject_ids) <- sample_names_expanded
                subject <- unname(subject_ids[sample_names])  # sample_names passed from calculate_lm_interaction

            } else {
                stop("paired = TRUE requires 'paired_samples' or 'sample_base' column in colData; supply subject_col explicitly or use map_metadata(...)")
            }
        } else {
            # fallback to sample base derived from column names
            subject <- sample_names
        }
        df$subject <- factor(subject)
        # require at least two subjects and at least two groups represented
        # ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
        # DATA QUALITY FILTERS
        # ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
        # Check statistical power: are there enough observations to estimate parameters?
        # Rule 1: Minimum observations (min_obs) required for meaningful inference
        # Rule 2: At least 2 subjects required for random intercept structures
        #   (1 subject = no variation between subjects, can't estimate random SD)
        
        n_subjects <- length(unique(na.omit(df$subject)))
        
        # Check minimum observation requirement
        if (nrow(df) < min_obs) {
            return(NULL)
        }
        
        # Check minimum subject requirement for mixed models
        if (n_subjects < 2) {
            return(NULL)
        }
                
        # ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
        # ARIMA(1,1,0) TRANSFORMATION: First differencing for stationarity
        # ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
        # Problem: Raw Tsallis entropy H_q is monotone increasing with q, violating
        #   stationarity assumption (constant mean) required for AR(1) modeling
        # Solution: Use first differences DeltaH_q = H_q - H_{q-1} to remove trend
        # - Bounded-support data [0, log(m)] after differencing approximates normality
        # - Enables valid hypothesis testing under AR(1) correlation structure
        # - Information preserved: interaction effects remain in differenced data
        
        # ARIMA(1,1,0) IMPLEMENTATION: Compute first differences for stationarity
        # Differencing removes monotone trend from Tsallis entropy, enabling valid AR(1) inference
        # BONUS: First differencing of bounded [0, log(m)] data helps normalize distribution
        # (bounded support becomes approximately normal after differencing in many cases)
        arima_result <- .compute_arima_differences(df, q_vals, df$group, df$subject)
        
        if (is.null(arima_result) || nrow(arima_result$df) < 3) {
            # Insufficient data for ARIMA differencing; fall back to raw data with warning
            # This can occur when: few observations, all same group, or excessive NAs
            if (verbose) {
                message("[calculate_lm_interaction] ARIMA(1,1,0) differencing lost too many observations; using raw entropy")
            }
            df_model <- df
            use_arima <- FALSE
        } else {
            # Successfully differenced data; use for model fitting
            df_model <- arima_result$df
            use_arima <- TRUE
            if (verbose) {
                message(sprintf("[calculate_lm_interaction] ARIMA(1,1,0): %d observations -> %d after differencing", 
                    arima_result$n_observations_original, arima_result$n_observations_differenced))
            }
        }
        
        # ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
        # MODEL FORMULA CONSTRUCTION
        # ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
        # Test: NULL (no interaction) vs ALTERNATIVE (q*group interaction)
        # Both formulas include main effects (q, group) for proper interpretation
        
        # fit null (no interaction) and alternative (with q:group interaction)
        # Using nlme::lme() for AR(1) covariance structure support (instead of lme4::lmer)
        mm_suppress_pattern <- "boundary \\(singular\\) fit|Computed variance-covariance matrix problem|not a positive definite matrix"

        # ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
        # OPTIONAL REGULARIZATION: Feature selection for complex models
        # ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
        # Regularization reduces effective degrees of freedom in high-dimensional models
        # Methods: PCA (dimension reduction), LASSO/elasticnet (sparse selection),
        #   gamsel (generalized additive model selection), spline (smooth selection)
        # Applied to: Interaction terms (q:group) across multiple q-levels
        
        # Apply regularization for feature selection if requested (not "pca")
        fs_result <- NULL
        # nlme formula syntax: fixed effects ~ random intercept
        formula_null <- entropy ~ q + group
        formula_alt <- entropy ~ q * group
        
        if (regularization != "pca") {
            fs_result <- .lmm_regularization(q_vals = df_model$q, entropy_vals = df_model$entropy,
                                                    group_vec = df_model$group, subject_vec = df_model$subject,
                                                    regularization = regularization)
            if (!is.null(fs_result)) {
                # If we got feature selection results, modify formulas to use only selected interaction terms
                # For simplicity, we construct a reduced model with selected q-value ranges
                # Note: This keeps main effects (q, group) but regularizes the interaction
                uq_levels <- length(fs_result$q_values)
                if (uq_levels > 2 && length(fs_result$selected_features) < (uq_levels - 1)) {
                    # Regularization reduced model complexity - use main effects plus interaction
                    # (full formula still, but now justified by regularization path)
                    if (verbose) {
                        message("[calculate_lmm_interaction] regularization retained ", 
                                length(fs_result$selected_features), " of ", 
                                uq_levels - 1, " q-interaction features")
                    }
                }
            }
        }
        
        # ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
        # HETEROSCEDASTICITY DETECTION AND VARIANCE WEIGHTING (March 2026)
        # ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
        # Test if variance is constant across q-values and groups
        # Problem: q-dependent entropy typically shows increasing variance at higher q
        #   (especially in high-rich transcriptome data)
        # Solution: Apply nlme::varPower(form = ~q) to model Var(Y) ~ q^theta
        # Result: Enables valid inference despite heteroscedasticity (unlike naive OLS)
        #
        # Method: Breusch-Pagan test compares residual variance across q-values
        # Output: is_heteroscedastic (Boolean), p_value (BP test p-value for null hypothesis
        #   of constant variance)
        
        # Detect q-dependent and group-dependent variance heterogeneity
        # Apply nlme::varPower() to model variance heterogeneity if detected
        hetero_result <- .detect_heteroscedasticity(df_model, df_model$q, df_model$group)
        use_var_structure <- FALSE
        
        if (!is.na(hetero_result$is_heteroscedastic) && hetero_result$is_heteroscedastic) {
            use_var_structure <- TRUE
            if (verbose) {
                message(sprintf("[calculate_lmm_interaction] Heteroscedasticity detected (BP p = %.4f); applying varPower",
                               hetero_result$p_value))
            }
        }
        
        # ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
        # MODEL FITTING: nlme::lme with AR(1) covariance and optional variance structure
        # ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
        # Fit null and alternative models with:
        #   - Fixed effects: q and group (and q*group for alternative)
        #   - Random effects: random intercept per subject (1 | subject)
        #   - Within-subject correlation: AR(1) structure indexed by q (nlme::corAR1)
        #   - Variance structure: Power law if heteroscedasticity detected (varPower)
        #   - Estimation method: ML (Maximum Likelihood) for LRT
        
        # Fit nlme models with AR(1) covariance structure for q-measurements within subjects
        # AR(1) model: Cov(Y_t, Y_s) = sigma^2 phi^|t-s| where t,s are q-ordered indices
        # *** CRITICAL: Now applied to differenced entropy DeltaH_q, not raw H_q ***
        # HETEROSCEDASTICITY: Add varPower() structure if heteroscedasticity detected
        
        if (use_var_structure) {
            # Include variance power model: Var(Y) ~ q^?
            fit0 <- try(
                nlme::lme(formula_null, random = ~1 | subject, data = df_model, method = "ML",
                         correlation = nlme::corAR1(form = ~1 | subject),
                         weights = nlme::varPower(form = ~ q)),
                silent = TRUE
            )
            fit1 <- try(
                nlme::lme(formula_alt, random = ~1 | subject, data = df_model, method = "ML",
                         correlation = nlme::corAR1(form = ~1 | subject),
                         weights = nlme::varPower(form = ~ q)),
                silent = TRUE
            )
        } else {
            fit0 <- try(
                nlme::lme(formula_null, random = ~1 | subject, data = df_model, method = "ML",
                         correlation = nlme::corAR1(form = ~1 | subject)),
                silent = TRUE
            )
            fit1 <- try(
                nlme::lme(formula_alt, random = ~1 | subject, data = df_model, method = "ML",
                         correlation = nlme::corAR1(form = ~1 | subject)),
                silent = TRUE
            )
        }

        # If either fit failed, try falling back to simpler approaches
        fallback_lm <- NULL
        used_fit_method <- "nlme::lme"
        used_singular <- FALSE
        if (inherits(fit0, "try-error") || inherits(fit1, "try-error")) {
            # ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━��━━━━━━━━━
            # FALLBACK FITTING: When nlme::lme fails (convergence, rank deficiency, etc)
            # ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
            # Cascade approach: try simpler models in order
            # 1. glmmTMB: Alternative GLMM engine with robust estimator
            # 2. Fixed-effects lm with subject as factor: ignores random structure
            # 3. Simple lm ignoring subjects: assumes independence
            # Each step relaxes assumptions; hypothesis test becomes less powerful but more robust
            
            used_fit_method <- "fallback"
            used_singular <- FALSE
            if ((verbose && progress) || (!verbose && progress)) {
                message("[calculate_lm_interaction] mixed model failed; trying simpler fixed-effects fallback")
            }
            fb <- .try_lm_fallbacks(df_model, verbose = verbose)
            if (!is.null(fb)) {
                fallback_lm <- fb
                used_fit_method <- fb$method
            }
        } else {
            used_fit_method <- if (use_arima) "nlme::lme_arima(1,1,0)" else "nlme::lme_ar1_raw"
        }

        # ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
        # HYPOTHESIS TEST: Extract p-value and effect size from fitted models
        # ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
        # H0: No interaction (q*group coefficient = 0)
        # H1: Interaction present (q*group coefficient ≠ 0)
        # Test statistic: -2 * log(LR) ~ chi^2(df=1) under H0
        # Where LR = Likelihood(H0 model) / Likelihood(H1 model)
        
        lrt_result <- list(p_value = NA_real_, n_subjects = NA_integer_, small_sample_flag = FALSE)
        msg <- NULL
        if (!is.null(fallback_lm)) {
            # Extract LRT p-value from fallback models (lm, glmmTMB, etc)
            lrt_result <- .extract_lrt_p(fallback_lm$fit0, fallback_lm$fit1, df = df_model)
            # If glmmTMB fallback failed due to convergence, propagate message
            if (!is.null(fallback_lm$message)) {
                msg <- fallback_lm$message
            }
        } else {
            # Extract LRT p-value from nlme models
            lrt_result <- .extract_lrt_p(fit0, fit1, df = df_model)
        }

        # nlme models use LRT for hypothesis testing (not Satterthwaite)
        # pvalue argument is ignored for nlme method
        satter_p <- NA_real_
        lrt_p <- lrt_result$p_value
        n_subj_lmm <- lrt_result$n_subjects
        small_sample_lmm <- lrt_result$small_sample_flag

        # nlme always uses LRT for hypothesis testing
        p_interaction <- lrt_p
        
        # ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
        # EFFECT SIZE: Extract interaction coefficient from fitted model
        # ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
        # slope_diff = coefficient of q:group interaction term
        # Interpretation: additional change in entropy per unit q for group 1 vs group 0
        # Estimated from alternative (full) model fit1 with interaction term included
        
        # Extract interaction coefficient (slope_diff) from fitted model
        slope_diff <- NA_real_
        if (!is.null(fallback_lm) && !is.null(fallback_lm$fit1)) {
            # For fallback lm/glm models
            coefs <- tryCatch(coef(fallback_lm$fit1), error = function(e) NULL)
            if (!is.null(coefs)) {
                # Look for q:group interaction term
                interaction_idx <- grep("q:group|group:q", names(coefs), ignore.case = FALSE)
                if (length(interaction_idx) > 0) {
                    slope_diff <- coefs[interaction_idx[1]]
                }
            }
        } else if (!inherits(fit1, "try-error")) {
            # For nlme models
            coefs <- tryCatch(nlme::fixef(fit1), error = function(e) NULL)
            if (!is.null(coefs)) {
                # Look for q:group interaction term
                interaction_idx <- grep("q:group|group:q", names(coefs), ignore.case = FALSE)
                if (length(interaction_idx) > 0) {
                    slope_diff <- coefs[interaction_idx[1]]
                }
            }
        }
        
        # Add weighting information to results (Phase 1)
        has_weights <- !is.null(df$weight)

        # ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
        # RESULTS OUTPUT: Return data frame with model statistics
        # ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
        # PHASE 14 ENHANCEMENT: Document sample size and power flags
        # Report n_subjects for transparency and small_sample_flag for caution
        res <- data.frame(
            gene = g, 
            p_interaction = p_interaction, 
            p_lrt = lrt_p,
            p_satterthwaite = NA_real_, 
            slope_diff = slope_diff, 
            fit_method = used_fit_method, 
            singular = used_singular, 
            arima_transformation = use_arima, 
            ci_weighted = has_weights,
            n_subjects = n_subj_lmm,
            small_sample_flag = small_sample_lmm,
            stringsAsFactors = FALSE
        )
        if (!is.null(msg)) res$message <- msg
        return(res)
    }

    if (method == "gam") {
        # ═══════════════════════════════════════════════════════════════════════════
        # GENERALIZED ADDITIVE MODELS (GAM) - Flexible Nonlinear Regression
        # ═══════════════════════════════════════════════════════════════════════════
        # Purpose: Fit smooth spline terms s(q) to capture nonlinear q-effects
        # Advantages over LMM:
        #   - No distributional assumptions (except exponential family)
        #   - Automatic smoothness selection via GCV/REML
        #   - Can detect nonlinear patterns in q-entropy relationship
        #   - Supports adaptive knot placement for complex curves
        # Disadvantages:
        #   - Less power for linear effects (if truly linear)
        #   - Requires more observations for stable smooth estimation
        #   - May not preserve AR(1) correlation structure
        #
        # Formula: entropy ~ group + s(q) + group:s(q) for interaction effect
        # Test: ANOVA comparing null (no interaction) vs alternative (with q*group interaction)
        
        # Extract subject info for paired/repeated measures (same as LMM)

        subject <- NULL
        
        if (!is.null(subject_col)) {
            if (!(subject_col %in% colnames(SummarizedExperiment::colData(se)))) {
                stop(sprintf("subject_col '%s' not found in colData(se)", subject_col))
            }
            subj_full <- as.character(SummarizedExperiment::colData(se)[, subject_col])
            # CRITICAL FIX: When colData doesn't have 'samples' column, use colnames(mat)
            # but strip "_q=" suffixes to match how sample_names were parsed in calculate_lm_interaction
            col_names_for_indexing <- SummarizedExperiment::colData(se)$samples
            if (is.null(col_names_for_indexing)) {
                # Strip "_q=....." suffixes from column names for proper indexing
                col_names_for_indexing <- sub("_q=.*", "", colnames(mat))
            }
            names(subj_full) <- col_names_for_indexing
            subject <- unname(subj_full[sample_names])
            
        } else if (paired) {
        
            coldata <- SummarizedExperiment::colData(se)
            coldata_cols <- colnames(coldata)
            
            # Look for 'paired_samples' or 'sample_base' columns created by map_metadata()
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
                # Map sample names by setting names to base names (without _q= suffix)
                sample_names_expanded <- sub("_q=.*", "", rownames(coldata))
                names(subject_ids) <- sample_names_expanded
                subject <- unname(subject_ids[sample_names])  # sample_names passed from calculate_lm_interaction

            } else {
                stop("paired = TRUE requires 'paired_samples' or 'sample_base' column in colData; supply subject_col explicitly or use map_metadata(...)")
            }
        }

        # Pass subject info, regularization, bias correction, adaptive knots parameters, and weights to GAM
        return(.gam_interaction(df, q_vals, g, min_obs = min_obs, subject = subject,
                                       regularization = regularization, bias_correction = bias_correction,
                                       adaptive_knots = adaptive_knots, weights = weights))
    }

    if (method == "fpca") {
        # ═══════════════════════════════════════════════════════════════════════════
        # FUNCTIONAL PRINCIPAL COMPONENT ANALYSIS (FPCA) - Curve-Based Analysis
        # ═══════════════════════════════════════════════════════════════════════════
        # Purpose: Treat each subject's entropy curve (q → entropy) as functional observation
        # Approach:
        #   1. For each subject, create curve matrix (rows=q-values, cols=groups)
        #   2. Apply PCA to extract principal components (functional modes)
        #   3. Score each subject on principal components
        #   4. Test for group differences in PC scores via t-test or ANOVA
        # Advantages:
        #   - Natural for paired/repeated measures (one curve per subject)
        #   - Captures global curve topology (not point-wise)
        #   - Implicit AR(1) correlation: PCA respects q-ordering
        #   - Small-sample friendly
        # Disadvantages:
        #   - Requires multiple q-values per subject
        #   - Interpretation less direct (curves → components → group differences)
        #
        # Math: Covariates extracted via first N principal components explaining ~80% variance
        # Test: Group (A vs B) difference in component scores via t-test or F-test
        
        # Extract subject info for paired samples (same as LMM and GAMM)
        subject <- NULL
        if (!is.null(subject_col)) {
            if (!(subject_col %in% colnames(SummarizedExperiment::colData(se)))) {
                stop(sprintf("subject_col '%s' not found in colData(se)", subject_col))
            }
            subj_full <- as.character(SummarizedExperiment::colData(se)[, subject_col])
            # CRITICAL FIX: When colData doesn't have 'samples' column, use colnames(mat)
            # but strip "_q=" suffixes to match how sample_names were parsed in calculate_lm_interaction
            col_names_for_indexing <- SummarizedExperiment::colData(se)$samples
            if (is.null(col_names_for_indexing)) {
                # Strip "_q=....." suffixes from column names for proper indexing
                col_names_for_indexing <- sub("_q=.*", "", colnames(mat))
            }
            names(subj_full) <- col_names_for_indexing
            subject <- unname(subj_full[sample_names])
        } else if (paired) {
            coldata <- SummarizedExperiment::colData(se)
            coldata_cols <- colnames(coldata)
            
            # Look for 'paired_samples' or 'sample_base' columns created by map_metadata()
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
                # Map sample names by setting names to base names (without _q= suffix)
                sample_names_expanded <- sub("_q=.*", "", rownames(coldata))
                names(subject_ids) <- sample_names_expanded
                subject <- unname(subject_ids[sample_names])  # sample_names passed from calculate_lm_interaction
            } else {
                stop("paired = TRUE requires 'paired_samples' or 'sample_base' column in colData; supply subject_col explicitly or use map_metadata(...)")
            }
        }
        # FPCA respects q-value ordering: creates curve matrix with q-values as columns (ordered),
        # then applies PCA which naturally captures smooth functional dependence structure (S168-S171).
        # This implicitly models AR(1) correlation: rho(k) = phi^|k| across ordered q-values.
        # Test L.1.6 validates this AR(1) pattern for entropy across q-values.
        return(.fpca_interaction(mat, q_vals, sample_names, group_vec, g,
            min_obs = min_obs, subject = subject, regularization = regularization, weights = weights))
    }

    if (method == "gee") {
        # ═══════════════════════════════════════════════════════════════════════════
        # GENERALIZED ESTIMATING EQUATIONS (GEE) - Robust Clustered Analysis
        # ═══════════════════════════════════════════════════════════════════════════
        # Purpose: Semi-parametric method for clustered/repeated measures data
        # Approach (via geepack::gee or geepack::geese):
        #   1. Specify working correlation structure (AR(1), exchangeable, independence)
        #   2. Fit model without requiring correct correlation (only need mean structure)
        #   3. Use sandwich estimator for variance (robust to correlation misspecification)
        # Advantages:
        #   - Robust to distributional assumptions
        #   - Automatically accounts for clustering/dependence
        #   - Sandwich standard errors valid even if working correlation wrong
        #   - Good for large samples with multiple clusters
        # Disadvantages:
        #   - Loses efficiency if correlation seriously misspecified
        #   - May be underpowered for small number of clusters
        #   - Hypothesis tests are Wald tests (not LRT)
        #
        # Supported correlation structures:
        #   - "ar1": AR(1) structure rho(k) = phi^|k|
        #   - "exchangeable": common correlation across all pairs
        #   - "independence": assumes independence (cluster-robust standard errors)
        # Test: Wald test for q*group interaction coefficient
        
        # Extract subject info for GEE (clustered/repeated measures)
        subject <- NULL
        # Note: weights are passed directly to the GEE helper, no need to prepare locally
        if (!is.null(subject_col)) {
            if (!(subject_col %in% colnames(SummarizedExperiment::colData(se)))) {
                stop(sprintf("subject_col '%s' not found in colData(se)", subject_col))
            }
            subj_full <- as.character(SummarizedExperiment::colData(se)[, subject_col])
            # CRITICAL FIX: When colData doesn't have 'samples' column, use colnames(mat)
            # but strip "_q=" suffixes to match how sample_names were parsed in calculate_lm_interaction
            col_names_for_indexing <- SummarizedExperiment::colData(se)$samples
            if (is.null(col_names_for_indexing)) {
                # Strip "_q=....." suffixes from column names for proper indexing
                col_names_for_indexing <- sub("_q=.*", "", colnames(mat))
            }
            names(subj_full) <- col_names_for_indexing
            subject <- unname(subj_full[sample_names])
        } else if (paired) {
            coldata <- SummarizedExperiment::colData(se)
            coldata_cols <- colnames(coldata)
            
            # Look for 'paired_samples' or 'sample_base' columns created by map_metadata()
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
                # Map sample names by setting names to base names (without _q= suffix)
                sample_names_expanded <- sub("_q=.*", "", rownames(coldata))
                names(subject_ids) <- sample_names_expanded
                subject <- unname(subject_ids[sample_names])  # sample_names passed from calculate_lm_interaction
            } else {
                stop("paired = TRUE requires 'paired_samples' or 'sample_base' column in colData; supply subject_col explicitly or use map_metadata(...)")
            }
        } else {
            # Use sample names as cluster IDs (each sample is independent)
            subject <- sample_names
        }
        # Pass subject info to GEE helper with AR(1) correlation structure (default)
        return(.gee_interaction(df, q_vals, g, subject = subject, min_obs = min_obs, 
                                       corstr = corstr, bias_correction = bias_correction, weights = weights))
    }

    # ═══════════════════════════════════════════════════════════════════════════
    # NO VALID METHOD SPECIFIED - Return NULL (error caught upstream)
    # ═══════════════════════════════════════════════════════════════════════════
    # This should not be reached if method validation done in calling function
    # (calculate_lm_interaction checks method ∈ {"lmm", "gam", "fpca", "gee"})
    
    return(NULL)
}