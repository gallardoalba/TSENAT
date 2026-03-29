# LMM regularization helper: performs feature selection on q-value interactions
# before fitting mixed model. Reduces overfitting with high-dimensional q-interaction terms.
.lmm_regularization <- function(q_vals, entropy_vals, group_vec, subject_vec = NULL,
                                      regularization = c("pca", "lasso", "elasticnet")) {
    regularization <- match.arg(regularization)
    
    if (regularization == "pca") {
        # PCA mode: no regularization, return NULL to skip feature selection
        return(NULL)
    }
    
    # Build feature matrix: create q-by-group interactions
    uq <- sort(unique(q_vals))
    group_levels <- unique(as.character(na.omit(group_vec)))
    
    if (length(group_levels) < 2) {
        return(NULL)
    }
    
    # Create design matrix with q and group:q interactions
    X <- cbind(q_vals)  # Include q as baseline
    
    # Add group indicator variable (for first group comparison)
    group_indicator <- as.numeric(factor(group_vec)) - 1
    X <- cbind(X, group_indicator)
    
    # Add q:group interaction terms
    for (qv in head(uq, -1)) {  # Avoid perfect collinearity with all q values
        X <- cbind(X, q_vals * group_indicator * (q_vals == qv))
    }
    
    # Remove columns with zero variance
    col_vars <- apply(X, 2, var, na.rm = TRUE)
    X <- X[, col_vars > 1e-10, drop = FALSE]
    
    if (ncol(X) < 2) {
        return(NULL)  # Not enough features for regularization
    }
    
    if (!requireNamespace("glmnet", quietly = TRUE)) {
        stop("Package 'glmnet' is required for LMM regularization")
    }
    
    # Fit regularized regression to identify important q:group interactions
    alpha_val <- if (regularization == "lasso") 1 else 0.5  # 1 for LASSO, 0.5 for Elastic Net
    
    cv_fit <- try(
        glmnet::cv.glmnet(
            x = X,
            y = entropy_vals,
            family = "gaussian",
            alpha = alpha_val,
            nfolds = min(5, length(entropy_vals) - 1),
            standardize = TRUE
        ),
        silent = TRUE
    )
    
    if (inherits(cv_fit, "try-error")) {
        return(NULL)
    }
    
    # Extract feature selection: which columns have non-zero coefficients at lambda.min
    coef_lambda_min <- stats::coef(cv_fit, s = "lambda.min")
    selected_features <- which(as.numeric(coef_lambda_min[-1]) != 0)  # Exclude intercept
    
    # If all features selected or none selected, return NULL to use full model
    if (length(selected_features) == 0 | length(selected_features) >= ncol(X) - 1) {
        return(NULL)
    }
    
    # Return the selected feature indices (these correspond to q-value interaction terms)
    list(
        selected_features = selected_features,
        feature_names = colnames(X)[selected_features],
        q_values = uq
    )
}

# Try lme4::lmer with multiple optimizers and controlled warnings.
.try_lmer <- function(formula, data, suppress_lme4_warnings = TRUE, verbose = FALSE,
    mm_suppress_pattern = "boundary \\(singular\\) fit|Computed variance-covariance matrix problem|not a positive definite matrix") {
    if (!requireNamespace("lme4", quietly = TRUE)) {
        stop("Package 'lme4' is required for mixed-model fitting")
    }
    opts <- list(list(optimizer = "bobyqa", optCtrl = list(maxfun = 2e+05)), list(optimizer = "nloptwrap",
        optCtrl = list(maxfun = 5e+05)))
    for (o in opts) {
        ctrl <- lme4::lmerControl(optimizer = o$optimizer, optCtrl = o$optCtrl)
        muffle_cond <- suppress_lme4_warnings || (!verbose)
        fit_try <- withCallingHandlers(try(lme4::lmer(formula, data = data, REML = FALSE,
            control = ctrl), silent = TRUE), warning = function(w) {
            if (muffle_cond && grepl(mm_suppress_pattern, conditionMessage(w), ignore.case = TRUE)) {
                invokeRestart("muffleWarning")
            }
        }, message = function(m) {
            if (muffle_cond && grepl(mm_suppress_pattern, conditionMessage(m), ignore.case = TRUE)) {
                invokeRestart("muffleMessage")
            }
        })
        if (!inherits(fit_try, "try-error")) {
            # check singularity if function available
            is_sing <- FALSE
            if (exists("isSingular", where = asNamespace("lme4"), inherits = FALSE)) {
                is_sing <- tryCatch(lme4::isSingular(fit_try, tol = 1e-04), error = function(e) FALSE)
            }
            attr(fit_try, "singular") <- is_sing
            return(fit_try)
        }
    }
    # all attempts failed
    return(structure("error", class = "try-error"))
}

## AR(1) correlation structure helper for ordered q-values (Phase 14 enhancement)
## Implements autocorrelated errors for entropy curves respecting q-order
.try_lmm_ar1 <- function(df, verbose = FALSE) {
    if (!requireNamespace("nlme", quietly = TRUE)) {
        return(NULL)
    }
    
    # AR(1) structure: Cov(ε_i,j, ε_i,k) = σ² * φ^|j-k|
    # Appropriate for entropy curves where H(q) is ordered and autocorrelated
    # Based on validation: papers confirm AR(1) decreasing covariance structure
    tryCatch({
        # Create time index for AR(1) ordering by q within each subject
        df <- df[order(df$subject, df$q), ]
        df$time_idx <- sequence(rle(as.character(df$subject))$lengths)
        
        fit0_ar1 <- nlme::lme(
            entropy ~ q + group,
            random = ~1 | subject,
            correlation = nlme::corAR1(form = ~time_idx | subject),
            data = df,
            method = "ML"
        )
        fit1_ar1 <- nlme::lme(
            entropy ~ q * group,
            random = ~1 | subject,
            correlation = nlme::corAR1(form = ~time_idx | subject),
            data = df,
            method = "ML"
        )
        
        if (!inherits(fit0_ar1, "try-error") && !inherits(fit1_ar1, "try-error")) {
            if (verbose) message("[.try_lmm_ar1] AR(1) correlation structure fitted successfully")
            return(list(fit0 = fit0_ar1, fit1 = fit1_ar1, method = "nlme_ar1"))
        }
        NULL
    }, error = function(e) {
        if (verbose) message("[.try_lmm_ar1] AR(1) fitting failed: ", conditionMessage(e))
        NULL
    })
}

## Consolidated helpers for calculate_lm_interaction fallbacks, LRT and Satterthwaite
## Improved mixed model handling with multiple fallback strategies (Phase 14 enhanced)
.try_lm_fallbacks <- function(df, verbose = FALSE) {
    # Strategy 0: Try AR(1) correlation for ordered q-values (NEW - Phase 14)
    ar1_result <- .try_lmm_ar1(df, verbose = verbose)
    if (!is.null(ar1_result)) {
        if (verbose) message("[.try_lm_fallbacks] Strategy 0 SUCCESS: AR(1) correlated random intercept")
        return(ar1_result)
    }
    
    # Strategy 1: Try nlme::lme() - more stable than lme4 for some datasets
    if (requireNamespace("nlme", quietly = TRUE)) {
        fit0_nlme <- try(nlme::lme(entropy ~ q + group, random = ~1 | subject, data = df,
            method = "ML"), silent = TRUE)
        fit1_nlme <- try(nlme::lme(entropy ~ q * group, random = ~1 | subject, data = df,
            method = "ML"), silent = TRUE)
        if (!inherits(fit0_nlme, "try-error") && !inherits(fit1_nlme, "try-error")) {
            if (verbose) message("[.try_lm_fallbacks] Strategy 1 SUCCESS: nlme random intercept")
            return(list(fit0 = fit0_nlme, fit1 = fit1_nlme, method = "nlme"))
        }
    }

    # Strategy 2: Try glmmTMB::glmmTMB() - newer, often more robust
    if (requireNamespace("glmmTMB", quietly = TRUE)) {
        fit0_tmb <- try(glmmTMB::glmmTMB(entropy ~ q + group + (1 | subject), data = df,
            REML = FALSE, verbose = FALSE), silent = TRUE)
        fit1_tmb <- try(glmmTMB::glmmTMB(entropy ~ q * group + (1 | subject), data = df,
            REML = FALSE, verbose = FALSE), silent = TRUE)
        if (!inherits(fit0_tmb, "try-error") && !inherits(fit1_tmb, "try-error")) {
            # Check for model convergence for both fits
            conv0 <- tryCatch({
                c0 <- fit0_tmb$fit$converged
                if (is.null(c0)) FALSE else isTRUE(c0)
            }, error = function(e) FALSE)
            conv1 <- tryCatch({
                c1 <- fit1_tmb$fit$converged
                if (is.null(c1)) FALSE else isTRUE(c1)
            }, error = function(e) FALSE)
            if (conv0 && conv1) {
                return(list(fit0 = fit0_tmb, fit1 = fit1_tmb, method = "glmmTMB"))
            } else {
                msg <- paste0("glmmTMB model did not converge: ",
                              "fit0 converged=", conv0, ", fit1 converged=", conv1)
                if (verbose) message("[.try_lm_fallbacks] ", msg)
                return(list(fit0 = NA, fit1 = NA, method = "glmmTMB", message = msg))
            }
        }
    }

    # Strategy 3: Linear model with subject as fixed effect (treated as factor)
    # Use factor() to ensure proper dummy variable coding, not raw numeric
    # Apply inverse-variance weights if available (Phase 1 weighting)
    fit0_lm <- try(stats::lm(entropy ~ q + group + factor(subject), data = df,
                             weights = if (!is.null(df$weight)) df$weight else NULL),
                   silent = TRUE)
    fit1_lm <- try(stats::lm(entropy ~ q * group + factor(subject), data = df,
                             weights = if (!is.null(df$weight)) df$weight else NULL),
                   silent = TRUE)
    if (!inherits(fit0_lm, "try-error") && !inherits(fit1_lm, "try-error")) {
        if (verbose) {
            message("[.try_lm_fallbacks] Using fixed-effect lm with factor(subject)")
        }
        return(list(fit0 = fit0_lm, fit1 = fit1_lm, method = "lm_subject_fixed"))
    }

    # Strategy 4: Last resort - drop subject entirely
    fit0_lm2 <- try(stats::lm(entropy ~ q + group, data = df,
                              weights = if (!is.null(df$weight)) df$weight else NULL),
                    silent = TRUE)
    fit1_lm2 <- try(stats::lm(entropy ~ q * group, data = df,
                              weights = if (!is.null(df$weight)) df$weight else NULL),
                    silent = TRUE)
    if (!inherits(fit0_lm2, "try-error") && !inherits(fit1_lm2, "try-error")) {
        if (verbose) {
            message("[.try_lm_fallbacks] Strategy 4: Subject removed - reduced power expected")
        }
        return(list(fit0 = fit0_lm2, fit1 = fit1_lm2, method = "lm_nosubject"))
    }

    if (verbose) message("[.try_lm_fallbacks] ALL STRATEGIES FAILED - no model fitted")
    return(NULL)
}


.extract_lrt_p <- function(fit0, fit1, df = NULL) {
    # Extract LRT p-value and related statistics from nested model comparison
    # NEW (Phase 14): Also return sample size info and power flags
    an <- try(stats::anova(fit0, fit1), silent = TRUE)
    if (!inherits(an, "try-error") && nrow(an) >= 2) {
        pcol <- grep("Pr\\(>F\\)|Pr\\(>Chisq\\)|Pr\\(>Chi\\)", colnames(an), value = TRUE)
        pval <- if (length(pcol) == 0) {
            as.numeric(an[2, ncol(an)])
        } else {
            as.numeric(an[2, pcol[1]])
        }
        
        # Compute sample size info for flagging
        n_subjects <- NA_integer_
        if (!is.null(df) && "subject" %in% colnames(df)) {
            n_subjects <- length(unique(df$subject))
        }
        
        # Flag: Type I error may be inflated with n_subjects < 5 (LOW POWER)
        small_sample_flag <- if (!is.na(n_subjects) && n_subjects < 5) TRUE else FALSE
        
        return(list(
            p_value = pval,
            n_subjects = n_subjects,
            small_sample_flag = small_sample_flag
        ))
    }
    
    # Fallback if anova fails
    n_subjects <- NA_integer_
    if (!is.null(df) && "subject" %in% colnames(df)) {
        n_subjects <- length(unique(df$subject))
    }
    
    return(list(
        p_value = NA_real_,
        n_subjects = n_subjects,
        small_sample_flag = if (!is.na(n_subjects) && n_subjects < 5) TRUE else FALSE
    ))
}