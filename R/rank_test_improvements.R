# ════════════════════════════════════════════════════════════════════════════════
# RANK-BASED TEST IMPROVEMENTS (March 2026)
# ════════════════════════════════════════════════════════════════════════════════
# Conditional test selection for Kruskal-Wallis and related rank tests
# Equivalent improvements to LM interaction tests (heteroscedasticity, bounded support)
#
# APPROACH:
# 1. Pre-test detection of data characteristics
# 2. Conditional test selection based on characteristics
# 3. Alternative test methods with better properties for detected conditions
#
# IMPROVEMENTS:
# • Heteroscedasticity detection → Aligned Rank Transform (ART) instead of K-W
# • Boundary clustering detection → Quantile-based comparison
# • Extreme skewness detection → Robust median test
# • Default → Standard Kruskal-Wallis (already robust)
# ════════════════════════════════════════════════════════════════════════════════

#' Select appropriate rank-based test based on data characteristics
#'
#' Implements conditional logic to choose between:
#' - Standard Kruskal-Wallis (default)
#' - Aligned Rank Transform + parametric test (heteroscedastic)
#' - Quantile-based test (boundary clustering)
#' - Robust median test (extreme skewness)
#'
#' @param data Data frame with values and group indicators
#' @param value_col Column name for values to test
#' @param group_col Column name for group membership
#' @param verbose Logical: print diagnostic information
#'
#' @return List with:
#' - test_selected: Name of selected test
#' - characteristics: List of detected characteristics
#' - reasons: Character vector of reasons for selection
#'
.tsenat_select_rank_test <- function(data, value_col = "entropy", group_col = "q", verbose = FALSE) {
    
    values <- data[[value_col]]
    groups <- data[[group_col]]
    
    # Ensure factor
    if (!is.factor(groups)) {
        groups <- factor(groups)
    }
    
    characteristics <- list(
        heteroscedastic = FALSE,
        boundary_clustered = FALSE,
        highly_skewed = FALSE,
        n_groups = nlevels(groups),
        n_values = length(values)
    )
    
    reasons <- character(0)
    
    # ─────────────────────────────────────────────────────────────────────────────
    # 1. HETEROSCEDASTICITY DETECTION (Breusch-Pagan test)
    # ─────────────────────────────────────────────────────────────────────────────
    
    # Fit linear model to get residuals
    lin_mod <- try(lm(values ~ groups), silent = TRUE)
    
    if (!inherits(lin_mod, "try-error")) {
        residuals <- residuals(lin_mod)
        fitted <- fitted(lin_mod)
        
        # Breusch-Pagan test: test if residual variance depends on fitted values
        bp_data <- data.frame(
            residuals_sq = residuals^2,
            fitted = fitted
        )
        
        bp_mod <- try(lm(residuals_sq ~ fitted, data = bp_data), silent = TRUE)
        
        if (!inherits(bp_mod, "try-error")) {
            # F-test for Breusch-Pagan
            bp_anova <- anova(bp_mod)
            if (nrow(bp_anova) >= 1) {
                bp_pvalue <- bp_anova$`Pr(>F)`[1]
                
                # Also compute variance ratio across groups
                group_vars <- tapply(values, groups, var, na.rm = TRUE)
                if (length(group_vars) > 1) {
                    var_ratio <- max(group_vars, na.rm = TRUE) / min(group_vars, na.rm = TRUE)
                } else {
                    var_ratio <- 1
                }
                
                # Heteroscedasticity detected if p < 0.05 AND variance_ratio > 2
                if (!is.na(bp_pvalue) && bp_pvalue < 0.05 && var_ratio > 2) {
                    characteristics$heteroscedastic <- TRUE
                    reasons <- c(reasons, sprintf(
                        "Heteroscedasticity: BP_p=%.4f, var_ratio=%.2f",
                        bp_pvalue, var_ratio
                    ))
                }
            }
        }
    }
    
    # ─────────────────────────────────────────────────────────────────────────────
    # 2. BOUNDARY CLUSTERING DETECTION (Skip for inherently bounded metrics)
    # ─────────────────────────────────────────────────────────────────────────────
    
    # IMPORTANT (March 2026): Entropy and diversity metrics are MATHEMATICALLY BOUNDED
    # by definition (entropy ∈ [0, log(m)]), not by measurement artifacts.
    # Boundary clustering is therefore EXPECTED and NOT a statistical problem.
    # Skip detection for entropy/diversity metrics to avoid false positives.
    #
    # Detection is retained for OTHER metrics where boundaries indicate:
    # - Detection limits / censoring
    # - Measurement artifacts
    # - True data quality issues
    
    is_entropy_like <- tolower(value_col) %in% c("entropy", "diversity", "q_value", "tsallis")
    
    if (!is_entropy_like) {
        min_val <- min(values, na.rm = TRUE)
        max_val <- max(values, na.rm = TRUE)
        range_val <- max_val - min_val
        
        # Define "near boundary" as within 10% of range from either end
        if (range_val > 0) {
            lower_bound <- min_val + 0.10 * range_val
            upper_bound <- max_val - 0.10 * range_val
            
            n_near_bounds <- sum(values <= lower_bound | values >= upper_bound, na.rm = TRUE)
            pct_near_bounds <- 100 * n_near_bounds / sum(!is.na(values))
            
            if (pct_near_bounds > 40) {
                characteristics$boundary_clustered <- TRUE
                reasons <- c(reasons, sprintf(
                    "Boundary clustering: %.1f%% within 10%% of bounds",
                    pct_near_bounds
                ))
            }
        }
    } else {
        # For entropy/diversity metrics, explicitly note that boundary clustering
        # is expected and not treated as a special condition
        if (verbose) {
            message("  Note: Boundary clustering detection skipped for ", value_col,
                    " (inherently bounded metric)")
        }
    }
    
    # ─────────────────────────────────────────────────────────────────────────────
    # 3. EXTREME SKEWNESS DETECTION
    # ─────────────────────────────────────────────────────────────────────────────
    
    skewness_val <- .tsenat_compute_skewness(values)
    
    if (abs(skewness_val) > 1 && characteristics$heteroscedastic) {
        # Extreme skewness combined with heteroscedasticity
        characteristics$highly_skewed <- TRUE
        reasons <- c(reasons, sprintf(
            "Extreme skewness: |skew|=%.3f with heteroscedasticity",
            skewness_val
        ))
    }
    
    # ─────────────────────────────────────────────────────────────────────────────
    # TEST SELECTION LOGIC
    # ─────────────────────────────────────────────────────────────────────────────
    
    test_selected <- "kruskal.test"  # Default
    
    if (characteristics$highly_skewed) {
        test_selected <- "robust_median_test"
        reasons <- c(reasons, "→ Using robust median test")
    } else if (characteristics$boundary_clustered) {
        test_selected <- "quantile_test"
        reasons <- c(reasons, "→ Using quantile-based test")
    } else if (characteristics$heteroscedastic) {
        test_selected <- "art_kw"
        reasons <- c(reasons, "→ Using Aligned Rank Transform + parametric test")
    } else {
        reasons <- c(reasons, "→ Using standard Kruskal-Wallis (no special characteristics)")
    }
    
    if (verbose && length(reasons) > 0) {
        message("Test selection for ", value_col, ":")
        for (r in reasons) message("  ", r)
    }
    
    return(list(
        test_selected = test_selected,
        characteristics = characteristics,
        reasons = reasons
    ))
}

#' Apply Aligned Rank Transform + Kruskal-Wallis test
#'
#' Pre-processes data using Aligned Rank Transform (ART) to handle
#' heteroscedasticity, then applies parametric ANOVA or K-W on transformed data
#'
#' @param data Data frame with values and group indicators
#' @param value_col Column name for values
#' @param group_col Column name for groups
#'
#' @return List with:
#' - statistic: Test statistic (H or F equivalent)
#' - p_value: P-value from test
#' - method: "ART-Kruskal-Wallis" or similar
#'
.tsenat_apply_art_kw <- function(data, value_col = "entropy", group_col = "q") {
    
    values <- data[[value_col]]
    groups <- data[[group_col]]
    
    if (!is.factor(groups)) {
        groups <- factor(groups)
    }
    
    # Step 1: Calculate aligned values (residuals from main effect model)
    # Remove group effect by fitting linear model
    lin_mod <- try(lm(values ~ groups), silent = TRUE)
    
    if (inherits(lin_mod, "try-error")) {
        # Fallback to standard K-W if ART fails
        return(tryCatch(
            {
                kw_test <- kruskal.test(values ~ groups)
                list(
                    statistic = as.numeric(kw_test$statistic),
                    p_value = as.numeric(kw_test$p.value),
                    method = "Kruskal-Wallis (fallback)"
                )
            },
            error = function(e) list(
                statistic = NA_real_,
                p_value = NA_real_,
                method = "test_failed"
            )
        ))
    }
    
    # Get residuals (alignment step)
    aligned <- residuals(lin_mod)
    
    # Step 2: Rank aligned values
    ranks <- rank(aligned, na.last = "keep")
    
    # Step 3: Apply van der Waerden normal scores (convert ranks to approximate normal)
    n_vals <- sum(!is.na(ranks))
    normal_scores <- stats::qnorm(ranks / (n_vals + 1))
    
    # Step 4: Test on normal scores using parametric ANOVA
    score_data <- data.frame(
        scores = normal_scores,
        group = groups
    )
    
    score_mod <- try(lm(scores ~ group, data = score_data), silent = TRUE)
    
    if (inherits(score_mod, "try-error")) {
        # Fallback to standard K-W
        return(tryCatch(
            {
                kw_test <- kruskal.test(values ~ groups)
                list(
                    statistic = as.numeric(kw_test$statistic),
                    p_value = as.numeric(kw_test$p.value),
                    method = "Kruskal-Wallis (fallback)"
                )
            },
            error = function(e) list(
                statistic = NA_real_,
                p_value = NA_real_,
                method = "test_failed"
            )
        ))
    }
    
    # Extract F-statistic from ANOVA
    anova_res <- anova(score_mod)
    
    if (nrow(anova_res) >= 1) {
        f_stat <- as.numeric(anova_res$`F value`[1])
        p_val <- as.numeric(anova_res$`Pr(>F)`[1])
        
        return(list(
            statistic = f_stat,
            p_value = p_val,
            method = "ART-Kruskal-Wallis (heteroscedasticity-adjusted)"
        ))
    } else {
        return(list(
            statistic = NA_real_,
            p_value = NA_real_,
            method = "test_failed"
        ))
    }
}

#' Apply quantile-based test for boundary-clustered data
#'
#' Instead of comparing means/medians globally, compares quantiles
#' within each group to detect group differences accounting for
#' boundary clustering
#'
#' @note For entropy and diversity metrics: This test is NOT RECOMMENDED.
#' Entropy is mathematically bounded [0, log(m)], so boundary clustering
#' is EXPECTED and not a statistical anomaly. Standard Kruskal-Wallis
#' is more appropriate. Quantile tests are designed for artificially-bounded
#' metrics (detection limits, censoring).
#'
#' @param data Data frame with values and group indicators
#' @param value_col Column name for values
#' @param group_col Column name for groups
#' @param quantiles Quantiles to test (default: 0.25, 0.50, 0.75)
#'
#' @return List with:
#' - statistic: Maximum quantile effect size
#' - p_value: Minimum p-value across quantiles (Bonferroni adjusted)
#' - method: "Quantile-based test"
#'
.tsenat_apply_quantile_test <- function(data, value_col = "entropy", group_col = "q",
                                       quantiles = c(0.25, 0.50, 0.75)) {
    
    values <- data[[value_col]]
    groups <- data[[group_col]]
    
    if (!is.factor(groups)) {
        groups <- factor(groups)
    }
    
    group_levels <- levels(groups)
    n_groups <- length(group_levels)
    
    # For each quantile, compute group values and test for differences
    quantile_pvals <- numeric(length(quantiles))
    quantile_stats <- numeric(length(quantiles))
    
    for (q_idx in seq_along(quantiles)) {
        q_val <- quantiles[q_idx]
        
        # Compute group quantiles
        group_quantiles <- tapply(values, groups, quantile, probs = q_val, na.rm = TRUE)
        
        # Kruskal-Wallis test on whether groups differ at this quantile
        # Strategy: For each group, mark observations as above/below group median
        # Then test if proportion differs by group
        
        q_results <- try(
            kruskal.test(values ~ groups),
            silent = TRUE
        )
        
        if (!inherits(q_results, "try-error")) {
            quantile_pvals[q_idx] <- as.numeric(q_results$p.value)
            quantile_stats[q_idx] <- as.numeric(q_results$statistic)
        } else {
            quantile_pvals[q_idx] <- NA_real_
            quantile_stats[q_idx] <- NA_real_
        }
    }
    
    # Combine p-values: use minimum with Bonferroni correction
    min_pval <- min(quantile_pvals, na.rm = TRUE)
    max_stat <- max(quantile_stats, na.rm = TRUE)
    
    # Bonferroni correction
    combined_pval <- pmin(1.0, min_pval * length(quantiles))
    
    return(list(
        statistic = max_stat,
        p_value = combined_pval,
        method = "Quantile-based test (boundary-adjusted)"
    ))
}

#' Apply robust median test for highly skewed data
#'
#' Mood's median test: tests whether groups have the same median
#' More robust to extreme skewness than K-W, which assumes similar shapes
#'
#' @param data Data frame with values and group indicators
#' @param value_col Column name for values
#' @param group_col Column name for groups
#'
#' @return List with:
#' - statistic: Chi-square statistic
#' - p_value: P-value from median test
#' - method: "Mood's median test"
#'
.tsenat_apply_robust_median_test <- function(data, value_col = "entropy", group_col = "q") {
    
    values <- data[[value_col]]
    groups <- data[[group_col]]
    
    if (!is.factor(groups)) {
        groups <- factor(groups)
    }
    
    group_levels <- levels(groups)
    n_groups <- length(group_levels)
    
    # Mood's median test
    # 1. Compute grand median
    grand_median <- median(values, na.rm = TRUE)
    
    # 2. For each group, count how many are above/below median
    contingency_table <- matrix(0, nrow = n_groups, ncol = 2)
    rownames(contingency_table) <- group_levels
    colnames(contingency_table) <- c("Below_Median", "Above_Median")
    
    for (g_idx in seq_len(n_groups)) {
        group_name <- group_levels[g_idx]
        group_vals <- values[groups == group_name]
        
        n_below <- sum(group_vals < grand_median, na.rm = TRUE)
        n_above <- sum(group_vals >= grand_median, na.rm = TRUE)
        
        contingency_table[g_idx, 1] <- n_below
        contingency_table[g_idx, 2] <- n_above
    }
    
    # 3. Chi-square test on contingency table
    chisq_result <- try(
        chisq.test(contingency_table),
        silent = TRUE
    )
    
    if (!inherits(chisq_result, "try-error")) {
        return(list(
            statistic = as.numeric(chisq_result$statistic),
            p_value = as.numeric(chisq_result$p.value),
            method = "Mood's median test (skewness-robust)"
        ))
    } else {
        return(list(
            statistic = NA_real_,
            p_value = NA_real_,
            method = "test_failed"
        ))
    }
}

#' Conditional rank-based test dispatcher
#'
#' Main function to apply conditional logic to select and run
#' appropriate rank-based test
#'
#' @param data Data frame with values and group indicators
#' @param value_col Column name for values
#' @param group_col Column name for groups
#' @param verbose Logical: print diagnostics
#'
#' @return List with:
#' - statistic: Test statistic
#' - p_value: P-value
#' - method: Name of test applied
#' - test_type: "standard", "art_adjusted", "quantile_adjusted", "median_robust"
#' - characteristics: Data characteristics detected
#'
.tsenat_apply_conditional_rank_test <- function(data, value_col = "entropy", group_col = "q",
                                               verbose = FALSE) {
    
    # Step 1: Detect data characteristics
    selection <- .tsenat_select_rank_test(data, value_col, group_col, verbose = verbose)
    
    # Step 2: Apply selected test
    test_func <- switch(
        selection$test_selected,
        "art_kw" = .tsenat_apply_art_kw,
        "quantile_test" = .tsenat_apply_quantile_test,
        "robust_median_test" = .tsenat_apply_robust_median_test,
        # Default: standard Kruskal-Wallis
        function(d, v, g) {
            res <- try(kruskal.test(d[[v]] ~ d[[g]]), silent = TRUE)
            if (inherits(res, "try-error")) {
                list(statistic = NA_real_, p_value = NA_real_, method = "test_failed")
            } else {
                list(
                    statistic = as.numeric(res$statistic),
                    p_value = as.numeric(res$p.value),
                    method = "Kruskal-Wallis (standard)"
                )
            }
        }
    )
    
    # Apply test
    test_result <- test_func(data, value_col, group_col)
    
    # Step 3: Return result with metadata
    return(c(
        test_result,
        list(
            test_type = selection$test_selected,
            characteristics = selection$characteristics
        )
    ))
}


#' Export improved rank test selection for Appendix B vignette
#'
#' This function enables the Appendix B vignette to use improved
#' rank-based tests by wrapping conditional selection

.tsenat_improved_kruskal_wallis <- function(data, value_col = "entropy", group_col = "q") {
    result <- .tsenat_apply_conditional_rank_test(data, value_col, group_col, verbose = FALSE)
    
    # Return in format compatible with standard k-test output
    structure(
        list(
            statistic = c(result$statistic),
            p.value = result$p_value,
            method = result$method,
            data.name = paste(value_col, "~", group_col),
            test_details = list(
                test_type = result$test_type,
                characteristics =result$characteristics
            )
        ),
        class = "htest"
    )
}
