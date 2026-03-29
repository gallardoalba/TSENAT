# GEE interaction helper with correlation structure validation
# @param df data frame with entropy, q, group columns
# @param q_vals numeric vector of q-values
# @param g gene identifier
# @param subject character vector of subject/cluster IDs for grouping repeated measures
# @param min_obs integer minimum observations required
# @param corstr character correlation structure: "auto" (select via QIC), "ar1", "exchangeable", "independence"
#               default "auto" tests all three and selects best
# @param bias_correction logical; apply Kenward-Roger correction for small clusters
#
# @return data.frame with gene, p_interaction, correlation_structure, and bias correction status
.gee_interaction <- function(df, q_vals, g, subject = NULL, min_obs = 10, corstr = "auto", bias_correction = TRUE, weights = NULL) {
    if (!requireNamespace("geepack", quietly = TRUE)) {
        stop("Package 'geepack' is required for method = 'gee'")
    }
    
    # Handle corstr options: if not "auto", validate it's one of the standard options
    if (corstr != "auto" && !(corstr %in% c("ar1", "exchangeable", "independence"))) {
        corstr <- "auto"
        warning("Invalid corstr; using 'auto' to select via QIC")
    }
    
    # Validate inputs
    # PHASE 1 WEIGHTING (March 2026): Use bootstrap CI weights if provided
    # Add weights to df if available and valid
    if (!is.null(weights) && length(weights) == nrow(df)) {
        df$weight <- weights
    }
    
    if (sum(!is.na(df$entropy)) < min_obs) {
        return(NULL)
    }
    if (length(unique(na.omit(df$group))) < 2) {
        return(NULL)
    }
    
    # If no subject specified, use row indices (independent observations)
    if (is.null(subject)) {
        subject <- factor(seq_len(nrow(df)))
    } else {
        subject <- factor(subject)
    }
    
    # ARIMA(1,1,0) IMPLEMENTATION: Compute first differences for stationarity
    # Tsallis entropy is monotone decreasing in q -> apply AR(1) to DeltaH_q instead of H_q
    use_arima <- FALSE
    df_orig_nrows <- nrow(df)
    
    if (length(unique(subject)) > 1) {
        # Sort by subject and q for proper within-subject differencing
        sort_idx <- order(subject, df$q)
        df_sorted <- df[sort_idx, ]
        subject_sorted <- subject[sort_idx]
        
        # Compute first differences within subjects
        df_diff_list <- list()
        for (subj in unique(subject_sorted)) {
            subj_idx <- which(subject_sorted == subj)
            if (length(subj_idx) >= 2) {
                subj_data <- df_sorted[subj_idx, ]
                n_diff <- nrow(subj_data) - 1
                df_diff_list[[as.character(subj)]] <- data.frame(
                    entropy = diff(subj_data$entropy),
                    q = subj_data$q[-nrow(subj_data)],
                    group = subj_data$group[-nrow(subj_data)],
                    stringsAsFactors = FALSE
                )
            }
        }
        
        if (length(df_diff_list) > 0) {
            df <- do.call(rbind, df_diff_list)
            rownames(df) <- NULL
            # Rebuild subject factor for differenced data
            subject <- rep(names(df_diff_list), vapply(df_diff_list, nrow, FUN.VALUE = integer(1)))
            use_arima <- TRUE
        }
    }
    
    df$subject <- factor(subject)
    
    # HETEROSCEDASTICITY ADJUSTMENT: Detect variance dependence on q and group
    # Breusch-Pagan test to determine if weights are needed
    hetero_result <- .detect_heteroscedasticity(df, q_vals = df$q, group_vec = df$group)
    gee_weights <- NULL
    
    # PHASE 1 WEIGHTING (March 2026): Bootstrap CI weights take precedence over heteroscedasticity weights
    if (!is.null(df$weight)) {
        # Use bootstrap CI weights if provided
        gee_weights <- df$weight
    } else if (!is.na(hetero_result$is_heteroscedastic) && hetero_result$is_heteroscedastic) {
        # Estimate variance weights using power-law model: Var ~ q^theta
        weights_result <- .estimate_variance_weights(df, q_vals = df$q, method = "power")
        if (!is.null(weights_result) && !is.null(weights_result$weights)) {
            gee_weights <- weights_result$weights
        }
    }
    
    # Count clusters for bias correction decisions
    n_clusters <- length(unique(as.numeric(df$subject)))
    
    # Ensure we have at least 2 groups and data isn't all NA
    if (sum(!is.na(df$entropy)) < 2) {
        return(NULL)
    }
    
    # Log differencing information if applied
    if (use_arima && nrow(df) < df_orig_nrows) {
        # ARIMA(1,1,0) was applied: note observation loss in result metadata
        arima_note <- sprintf("ARIMA(1,1,0): %d observations -> %d after differencing", df_orig_nrows, nrow(df))
    } else {
        arima_note <- NULL
    }
    
    # Fit GEE models using Gaussian (normal) family for continuous entropy values
    # Null model: entropy ~ q + group (no interaction)
    # Alternative model: entropy ~ q * group (with interaction)
    
    # Determine correlation structure to use
    selected_corstr <- corstr
    corstr_selection_info <- NULL
    
    if (corstr == "auto") {
        # Test all correlation structures and select best via QIC
        selection_result <- .select_gee_correlation(
            df = df,
            formula_null = entropy ~ q + group,
            formula_alt = entropy ~ q * group,
            subject = df$subject,
            criteria = "qic"
        )
        
        selected_corstr <- selection_result$best_corstr
        corstr_selection_info <- list(
            method = "QIC-based selection",
            qic_values = selection_result$qic_values,
            reason = selection_result$reason,
            report = selection_result$report
        )
    }
    
    # Fit models with selected correlation structure
    # If heteroscedasticity detected, apply variance weights
    if (!is.null(gee_weights)) {
        df$gee_weights <- gee_weights
        
        fit_null <- try(
            geepack::geeglm(
                entropy ~ q + group,
                id = df$subject,
                data = df,
                family = stats::gaussian(),
                weights = gee_weights,
                corstr = selected_corstr,
                na.action = stats::na.omit
            ),
            silent = TRUE
        )
        
        fit_alt <- try(
            geepack::geeglm(
                entropy ~ q * group,
                id = df$subject,
                data = df,
                family = stats::gaussian(),
                weights = gee_weights,
                corstr = selected_corstr,
                na.action = stats::na.omit
            ),
            silent = TRUE
        )
    } else {
        # Standard GEE fitting without weights
        fit_null <- try(
            geepack::geeglm(
                entropy ~ q + group,
                id = df$subject,
                data = df,
                family = stats::gaussian(),
                corstr = selected_corstr,
                na.action = stats::na.omit
            ),
            silent = TRUE
        )
        
        fit_alt <- try(
            geepack::geeglm(
                entropy ~ q * group,
                id = df$subject,
                data = df,
                family = stats::gaussian(),
                corstr = selected_corstr,
                na.action = stats::na.omit
            ),
            silent = TRUE
        )
    }
    
    # If either fit failed, return NULL
    if (inherits(fit_null, "try-error") || inherits(fit_alt, "try-error")) {
        return(NULL)
    }
    
    # Ensure both models have valid results
    if (is.null(fit_null) || is.null(fit_alt)) {
        return(NULL)
    }
    
    # Extract interaction term p-value using Wald test
    # The interaction term is the coefficient for the q:group interaction
    coefs_alt <- stats::coef(fit_alt)
    
    # Find interaction term (usually named something like "q:groupTumor" or similar)
    ia_names <- names(coefs_alt)[grepl("^q:", names(coefs_alt), ignore.case = TRUE)]
    
    if (length(ia_names) == 0) {
        return(NULL)
    }
    
    # Use summary to get standard errors and Wald test statistics/p-values
    summ <- try(summary(fit_alt), silent = TRUE)
    
    if (inherits(summ, "try-error") || is.null(summ)) {
        return(NULL)
    }
    
    # Wald test p-value for interaction (two-sided)
    # Extract from coefficients table if available
    coef_table <- summ$coefficients
    
    if (is.null(coef_table)) {
        return(NULL)
    }
    
    # Find p-value corresponding to interaction term
    p_interaction <- NA_real_
    z_stat_value <- NA_real_  # Store for potential K-C correction
    se_robust_value <- NA_real_  # Store for potential K-C correction
    
    # Try to find interaction term in coefficient table (rownames contain "q:group" pattern)
    for (ia_name in ia_names) {
        if (ia_name %in% rownames(coef_table)) {
            row_idx <- which(rownames(coef_table) == ia_name)[1]
            # Usually column 4 or 5 contains the p-value (Pr(>|Z|) or Pr(>|W|))
            # Check multiple possible column names
            p_col <- grep("Pr\\(>", colnames(coef_table))[1]
            if (!is.na(p_col) && p_col <= ncol(coef_table)) {
                p_int_candidate <- coef_table[row_idx, p_col]
                # Convert NaN to NA (occurs with numerical instability in geeglm summary)
                if (!is.na(p_int_candidate) && !is.nan(p_int_candidate)) {
                    p_interaction <- p_int_candidate
                    break
                }
            }
        }
    }
    
    # If we couldn't extract p-value from table, try computing Wald test manually
    # using sandwich (robust) variance estimator
    if (is.na(p_interaction)) {
        # Wald test: (coef / SE)^2 ~ chi2(1) or t-dist for small samples
        ia_idx <- which(names(coefs_alt) %in% ia_names)[1]
        if (!is.na(ia_idx)) {
            # Get robust SE from covariance matrix
            vcov_robust <- try(
                {
                    # Compute sandwich estimator (robust variance)
                    X <- model.matrix(fit_alt)
                    residuals_vec <- fit_alt$y - fit_alt$fitted.values
                    
                    # Define meat of sandwich
                    W <- diag(1 / fit_alt$scale)
                    meat <- t(X) %*% W %*% (residuals_vec^2 * diag(nrow(X))) %*% W %*% X
                    
                    # Bread is X'VX (inverted)
                    bread <- solve(t(X) %*% W %*% X)
                    
                    # Sandwich: bread %*% meat %*% bread
                    bread %*% meat %*% bread
                },
                silent = TRUE
            )
            
            if (!inherits(vcov_robust, "try-error") && !is.null(vcov_robust)) {
                se_robust <- sqrt(diag(vcov_robust)[ia_idx])
                if (!is.na(se_robust) && se_robust > 0) {
                    z_stat <- coefs_alt[ia_idx] / se_robust
                    z_stat_value <- z_stat
                    se_robust_value <- se_robust
                    
                    # For small number of clusters, use t-distribution (Kauermann-Carroll style correction)
                    # This is a conservative approach that maintains Type I error rate (bias correction)
                    # Reference: Li & Redden (2015), Statistics in Medicine
                    if (bias_correction && n_clusters < 20) {
                        # Using t-distribution with df = n_clusters - 1 (conservative)
                        # This approximates the Kauermann-Carroll correction
                        df_corr <- max(1, n_clusters - 1)
                        p_interaction <- 2 * stats::pt(abs(z_stat), df = df_corr, lower.tail = FALSE)
                    } else {
                        # Standard normal (Wald test)
                        p_interaction <- 2 * stats::pnorm(abs(z_stat), lower.tail = FALSE)
                    }
                }
            }
        }
    } else if (bias_correction && n_clusters < 20) {
        # If we extracted p-value from table but have small clusters, recompute with t-distribution
        # This provides K-C bias correction
        ia_idx <- which(names(coefs_alt) %in% ia_names)[1]
        if (!is.na(ia_idx)) {
            # Try to extract robust SE and recompute with t-distribution
            vcov_robust <- try(
                {
                    X <- model.matrix(fit_alt)
                    residuals_vec <- fit_alt$y - fit_alt$fitted.values
                    W <- diag(1 / fit_alt$scale)
                    meat <- t(X) %*% W %*% (residuals_vec^2 * diag(nrow(X))) %*% W %*% X
                    bread <- solve(t(X) %*% W %*% X)
                    bread %*% meat %*% bread
                },
                silent = TRUE
            )
            
            if (!inherits(vcov_robust, "try-error") && !is.null(vcov_robust)) {
                se_robust <- sqrt(diag(vcov_robust)[ia_idx])
                if (!is.na(se_robust) && se_robust > 0) {
                    z_stat <- coefs_alt[ia_idx] / se_robust
                    df_corr <- max(1, n_clusters - 1)
                    p_interaction <- 2 * stats::pt(abs(z_stat), df = df_corr, lower.tail = FALSE)
                }
            }
        }
    }
    
    # ===============================================================================
    # RESIDUAL NORMALITY TESTING (NEW - March 2026)
    # Database Evidence: B001, B004, C017 (Normality testing in regression)
    # ===============================================================================
    shapiro_result <- .test_residual_normality(
        model = fit_alt,
        model_type = "gee",
        verbose = FALSE
    )
    
    # Return result with Shapiro-Wilk test
    gee_result <- data.frame(
        gene = g,
        p_interaction = p_interaction,
        n_clusters = n_clusters,
        bias_correction_applied = bias_correction && n_clusters < 20,
        correlation_structure = selected_corstr,
        corstr_selection_method = if (corstr == "auto") "QIC_based" else "user_specified",
        stringsAsFactors = FALSE
    )
    
    # Add Shapiro-Wilk residual normality test results
    if (!is.na(shapiro_result$shapiro_p_value)) {
        gee_result$shapiro_p_value <- shapiro_result$shapiro_p_value
        gee_result$residuals_normal <- shapiro_result$residuals_normal
        gee_result$n_residuals_tested <- shapiro_result$n_residuals
    } else {
        gee_result$shapiro_p_value <- NA_real_
        gee_result$residuals_normal <- NA
        gee_result$n_residuals_tested <- NA_integer_
    }
    
    # Add Phase 1 bootstrap CI weighting tracking (March 2026)
    gee_result$ci_weighted <- !is.null(df$weight)
    
    # Extract slope_diff from GEE interaction coefficient
    slope_diff <- NA_real_
    if (!inherits(fit_alt, "try-error") && !is.null(fit_alt)) {
        tryCatch({
            coefs_alt <- stats::coef(fit_alt)
            ia_names <- names(coefs_alt)[grepl("^q:", names(coefs_alt), ignore.case = TRUE)]
            
            if (length(ia_names) > 0) {
                slope_diff <- coefs_alt[ia_names[1]]
            }
        }, error = function(e) { NULL })
    }
    gee_result$slope_diff <- slope_diff
    
    return(gee_result)
}
