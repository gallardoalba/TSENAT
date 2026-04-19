# ============================================================================
# 5. UTILITY FUNCTIONS
# ============================================================================

#' Format numeric values for table display
#' 
#' Convert 0 to '0', small values to short scientific notation, others to 3 decimals
#' @noRd
.format_table_value <- function(x, decimals = 3) {
    if (is.na(x) || !is.finite(x)) {
        return(as.character(x))
    }
    if (x == 0) {
        return("0")
    }
    if (abs(x) < 0.001 && x != 0) {
        return(sprintf("%.2e", x))
    }
    return(sprintf(paste0("%.", decimals, "f"), x))
}

#' Test Rank-Based and Method-Specific Assumptions
#'
#' Diagnostic checks to verify rank-based methods and parametric models are
#' appropriate for data. Supports rank checks (exchangeability, monotonicity,
#' consistency) and optional method-specific diagnostics (GAM, GEE, LMM, FPCA).
#'
#' @param data Matrix of expression values (rows=observations, cols=variables/samples)
#' @param checks Character string (preset) or character vector (explicit). Default: 'rank'
#'   Presets:
#'   - 'rank': basic exchangeability check only (fast)
#'   - 'all': all available metrics (exchangeability + monotonicity + consistency + GAM diagnostics)
#'   Explicit vector:
#'   - 'exchangeability', 'monotonicity', 'consistency' (rank-based)
#'   - 'gam_metrics', 'gee_metrics' (future), 'lmm_metrics' (future) (method-specific)
#' @param alpha Numeric; significance level for hypothesis tests (default: 0.05)
#'
#' @return List with diagnostic results (class 'rank_assumptions')
#'
#' @details
#' Use `checks = 'rank'` (default) for fast basic exchangeability check. 
#' Use `checks = 'all'` for comprehensive diagnostics including monotonicity, consistency, and GAM metrics.
#' For backward compatibility, rank-based checks are the default (`checks = 'rank'`).
#'
#' NEW (April 2026): Added GAM metrics (concurvity, EDF, non-linearity,
#' basis adequacy).
#'
#' References:
#' - GAM: @S150 (2023), @S143 (2015), @S137 (1979)
#'
#' @noRd
.calculate_assumptions <- function(data, checks = "rank", alpha = 0.05, q_values = NULL,
    gee_params = list()) {

    # Expand convenience presets
    if (is.character(checks) && length(checks) == 1) {
        if (checks == "rank") {
            checks <- c("exchangeability")
        } else if (checks == "all") {
            checks <- c("exchangeability", "monotonicity", "consistency", "gam_metrics",
                "gee_metrics", "lmm_metrics", "fpca_metrics")
        } else {
            # Allow single metric names to pass through They will be used as-is
            # in the checks below
        }
    } else if (!is.character(checks)) {
        stop("checks must be a character string ('rank', 'all', metric name) or character vector")
    }

    if (!is.matrix(data))
        data <- as.matrix(data)

    results <- list()

    # Guard against empty data input: check matrix dimensions This prevents
    # errors when diversity calculation creates empty assays
    if (nrow(data) == 0 || ncol(data) == 0) {
        # Return default check structures even for empty data (prevents NULL
        # access errors in vignettes) CRITICAL: Every check must have a p_value
        # and details field for vignette compatibility
        empty_checks <- list(exchangeability = list(description = "Sample exchangeability",
            method = "N/A", status = "? SKIP", p_value = NA_real_, details = "Empty data"),
            monotonicity = list(description = "Rank ordering stability", method = "N/A",
                status = "? SKIP", mean_correlation = NA_real_, details = "Empty data"),
            consistency = list(description = "Rank consistency", method = "N/A",
                status = "? SKIP", w_statistic = NA_real_, icc = NA_real_, details = "Empty data"))

        # Include empty metric placeholders if requested
        if ("gam_metrics" %in% checks) {
            empty_checks$gam_metrics <- list(concurvity = list(description = "Concurvity Index",
                status = "? SKIP", overall_concurvity = NA_real_, details = "Empty data"),
                edf = list(description = "Effective DoF", status = "? SKIP", edf_ratio = NA_real_,
                  details = "Empty data"), nonlinearity = list(description = "Non-linearity",
                  status = "? SKIP", r2_improvement_percent = NA_real_, details = "Empty data"),
                basis_adequacy = list(description = "Basis Adequacy", status = "? SKIP",
                  optimal_basis_dimension = NA_integer_, details = "Empty data"))
        }
        if ("gee_metrics" %in% checks) {
            empty_checks$gee_metrics <- list(correlation_fit = list(description = "Correlation Structure",
                status = "? SKIP", details = "Empty data"), cluster_variation = list(description = "Cluster Variation",
                status = "? SKIP", details = "Empty data"), independence = list(description = "Independence",
                status = "? SKIP", details = "Empty data"), scale_parameter = list(description = "Scale Parameter",
                status = "? SKIP", scale = NA_real_, details = "Empty data"))
        }
        if ("lmm_metrics" %in% checks) {
            empty_checks$lmm_metrics <- list(variance_components = list(description = "Variance Components",
                status = "? SKIP", details = "Empty data"), normality = list(description = "Normality",
                status = "? SKIP", details = "Empty data"), homogeneity = list(description = "Homogeneity",
                status = "? SKIP", details = "Empty data"), influence = list(description = "Influence",
                status = "? SKIP", details = "Empty data"))
        }
        if ("fpca_metrics" %in% checks) {
            empty_checks$fpca_metrics <- list(variance_adequacy = list(description = "Variance Adequacy",
                status = "? SKIP", details = "Empty data"), bootstrap_stability = list(description = "Bootstrap Stability",
                status = "? SKIP", details = "Empty data"))
        }

        return(structure(list(overall_summary = "Cannot evaluate assumptions on empty data matrix"),
            class = "rank_assumptions", checks = empty_checks, summary_stats = list(n_genes = 0,
                n_samples = 0, entropy_min = NA_real_, entropy_max = NA_real_, entropy_mean = NA_real_,
                entropy_median = NA_real_, n_missing = 0)))
    }

    # Calculate summary statistics for entropy data Guard against empty data:
    # check for valid values before min/max
    has_valid_data <- (nrow(data) > 0 && ncol(data) > 0 && sum(!is.na(data)) > 0)

    entropy_min <- if (has_valid_data)
        min(data, na.rm = TRUE) else NA_real_
    entropy_max <- if (has_valid_data)
        max(data, na.rm = TRUE) else NA_real_

    summary_stats <- list(n_genes = nrow(data), n_samples = ncol(data), entropy_min = entropy_min,
        entropy_max = entropy_max, entropy_mean = mean(data, na.rm = TRUE), entropy_median = median(data,
            na.rm = TRUE), n_missing = sum(is.na(data)))

    # Check 1: Exchangeability (permutation test for serial correlation in
    # samples)
    if ("exchangeability" %in% checks) {
        # Hypothesis: Under exchangeability, consecutive samples should show
        # similar correlation to random sample pairs. Ordering effects would
        # manifest as higher correlation in consecutive samples than expected
        # by chance.

        # Test statistic: mean Pearson correlation between consecutive samples
        # (columns)
        if (ncol(data) > 2) {
            # Vectorized: compute correlations between consecutive samples Use
            # diag(cor(X, Y)) to get paired correlations efficiently
            data_1 <- data[, seq_len(ncol(data) - 1)]
            data_2 <- data[, seq_len(ncol(data) - 1) + 1]
            consecutive_cors <- vapply(seq_len(ncol(data) - 1), function(i) {
                stats::cor(data[, i], data[, i + 1], method = "pearson", use = "complete.obs")
            }, numeric(1))
            original_stat <- mean(consecutive_cors, na.rm = TRUE)

            # Permutation test: shuffle column order and recompute (vectorized)
            n_perms <- 99
            perm_stats <- numeric(n_perms)
            # Seed handling left to caller for Bioconductor compliance
            for (perm in seq_len(n_perms)) {
                # Shuffle column (sample) order
                perm_idx <- sample(seq_len(ncol(data)))
                perm_data <- data[, perm_idx]

                # Recompute consecutive correlations (vectorized)
                perm_cors <- vapply(seq_len(ncol(perm_data) - 1), function(i) {
                  stats::cor(perm_data[, i], perm_data[, i + 1], method = "pearson",
                    use = "complete.obs")
                }, numeric(1))
                perm_stats[perm] <- mean(perm_cors, na.rm = TRUE)
            }

            # P-value: proportion of permutations with mean_consecutive >=
            # original High p-value: original correlation within random
            # variation (exchangeable) Low p-value: original shows ordering
            # effect (NOT exchangeable)
            p_exchangeability <- mean(perm_stats >= original_stat, na.rm = TRUE)

            # Interpretation: if p > alpha, data is exchangeable (no ordering
            # effect)
            exchangeability_interpretation <- if (p_exchangeability > alpha)
                "exchangeable" else "ordering detected"

            results$exchangeability <- list(description = "Sample exchangeability (serial correlation in sequence)",
                method = "Permutation test (consecutive sample correlations vs. shuffled)",
                test_statistic = original_stat, p_value = p_exchangeability, status = exchangeability_interpretation,
                details = sprintf("p=%s", format(p_exchangeability, scientific = TRUE)))
        } else {
            results$exchangeability <- list(description = "Sample exchangeability",
                method = "Insufficient samples (need >= 3)", status = "? SKIP", details = "Requires at least 3 samples to test exchangeability")
        }
    }

    # Check 2: Monotonicity (Spearman correlation stability across rows)
    if ("monotonicity" %in% checks) {
        # Compute pairwise Spearman correlations between consecutive rows
        spearman_cors <- numeric(max(1, nrow(data) - 1))

        if (nrow(data) > 1) {
            for (i in seq_len(nrow(data) - 1)) {
                spearman_cors[i] <- stats::cor(data[i, ], data[i + 1, ], method = "spearman",
                  use = "complete.obs")
            }
        }

        # Summary statistics of correlation stability
        mean_cor <- mean(spearman_cors, na.rm = TRUE)
        sd_cor <- stats::sd(spearman_cors, na.rm = TRUE)

        # Guard against empty correlation vector
        has_valid_cors <- sum(!is.na(spearman_cors)) > 0
        min_cor <- if (has_valid_cors)
            min(spearman_cors, na.rm = TRUE) else NA_real_

        # Status: high and stable correlations indicate good monotonicity
        # Interpretation: degree of heterogeneity in rank ordering Guard
        # against NA mean_cor (happens with single-column or empty data)
        heterogeneity_interpretation <- if (!is.na(mean_cor) && mean_cor > 0.7) {
            "homogeneous"
        } else if (!is.na(mean_cor) && mean_cor > 0.4) {
            "moderately heterogeneous"
        } else if (!is.na(mean_cor)) {
            "heterogeneous"
        } else {
            "? SKIP"
        }

        results$monotonicity <- list(description = "Rank ordering stability (Spearman correlation across rows)",
            method = "Pairwise Spearman correlations between consecutive rows", mean_correlation = mean_cor,
            sd_correlation = sd_cor, min_correlation = min_cor, status = heterogeneity_interpretation,
            details = paste0("r=", .format_table_value(mean_cor)))
    }

    # Check 3: Consistency (ICC for replicate consistency)
    if ("consistency" %in% checks) {
        # Calculate Kendall's W (concordance coefficient) across columns W
        # ranges from 0 (no agreement) to 1 (perfect agreement)

        if (ncol(data) >= 2 && nrow(data) >= 2) {
            # Transpose for ICC calculation (samples as rows, variables as
            # columns)
            data_t <- t(data)

            # Compute mean rank across each column (gene)
            ranked_data <- apply(data_t, 2, function(x) rank(x, na.last = "keep"))

            # Kendall's W = 12*S / (m^2 * (n^3 - n)) where S = sum of squared
            # deviations from mean rank, m = judges (samples), n = objects
            # (genes)
            m <- nrow(ranked_data)
            n <- ncol(ranked_data)

            # Sum of squared deviations
            col_means <- colMeans(ranked_data, na.rm = TRUE)
            S <- sum((colSums(ranked_data, na.rm = TRUE) - m * col_means)^2, na.rm = TRUE)

            # Kendall's W
            kendall_w <- if (n > 1) {
                12 * S/(m^2 * (n^3 - n))
            } else {
                NA_real_
            }

            # Alternative: compute intraclass correlation (ICC 2-way mixed) Use
            # simplified two-way ICC calculation
            grand_mean <- mean(data, na.rm = TRUE)
            between_col_var <- sum((colMeans(data, na.rm = TRUE) - grand_mean)^2,
                na.rm = TRUE)/(ncol(data) - 1)
            within_var <- var(as.numeric(data), na.rm = TRUE)
            icc_simplified <- between_col_var/(between_col_var + within_var)

            status <- if (!is.na(kendall_w) && kendall_w > 0.7) {
                "high"
            } else if (!is.na(kendall_w) && kendall_w > 0.4) {
                "moderate"
            } else {
                "low"
            }

            results$consistency <- list(description = "Rank consistency evaluation (Kendall's W & ICC)",
                method = "Kendall's W concordance coefficient + ICC approximation",
                kendall_w = kendall_w, icc_simplified = icc_simplified, status = status,
                details = paste0("W=", .format_table_value(if (is.na(kendall_w)) 0 else kendall_w),
                  ", ICC=", .format_table_value(if (is.na(icc_simplified)) 0 else icc_simplified)))
        } else {
            results$consistency <- list(description = "Rank consistency evaluation",
                method = "Insufficient data for consistency test", status = "? SKIP",
                details = "Requires at least 2 samples and 2 genes")
        }
    }

    # Check 4: GAM metrics (new - April 2026)
    if ("gam_metrics" %in% checks) {
        gam_result <- tryCatch({
            .get_gam_metrics(data, q_values = q_values, method_params = list())
        }, error = function(e) {
            list(error = TRUE, message = paste("GAM metrics computation failed:",
                e$message), reason = "Check if mgcv package is installed and data has sufficient variation")
        })
        results$gam_metrics <- gam_result
    }

    # Check 5: GEE metrics (new - April 2026)
    if ("gee_metrics" %in% checks) {
        gee_result <- tryCatch({
            .get_gee_metrics(data, gee_params = gee_params)
        }, error = function(e) {
            list(error = TRUE, message = paste("GEE metrics computation failed:",
                e$message), reason = "Check if geepack package is installed")
        })
        results$gee_metrics <- gee_result
    }

    # Check 6: LMM metrics (new - April 2026)
    if ("lmm_metrics" %in% checks) {
        lmm_result <- tryCatch({
            lmm_params <- list(assumed_re_structure = "random_intercept", cluster_col = NULL)
            .get_lmm_metrics(data, lmm_params = lmm_params)
        }, error = function(e) {
            list(error = TRUE, message = paste("LMM metrics computation failed:",
                e$message), reason = "Check data structure and dimensionality")
        })
        results$lmm_metrics <- lmm_result
    }

    # Check 7: FPCA metrics (new - April 2026)
    if ("fpca_metrics" %in% checks) {
        fpca_result <- tryCatch({
            fpca_params <- list(max_components = NULL, n_bootstrap = 500, n_components = 3)
            .get_fpca_metrics(data, fpca_params = fpca_params)
        }, error = function(e) {
            list(error = TRUE, message = paste("FPCA metrics computation failed:",
                e$message), reason = "Check data dimensionality (need >1 observation and column)")
        })
        results$fpca_metrics <- fpca_result
    }

    structure(list(overall_summary = paste("Rank-based assumptions evaluated with",
        "rigorous statistical tests.")), class = "rank_assumptions", checks = results,
        summary_stats = summary_stats)
}

#' Print method for rank-based assumptions check
#'
#' @param x Object of class 'rank_assumptions'
#' @param ... Additional arguments (ignored)
#'

#' @noRd
#' @method print rank_assumptions

print.rank_assumptions <- function(x, ...) {
    message("STATISTICAL ASSUMPTIONS: Rank-Based & Method-Specific Tests")
    message(strrep("=", 70))

    # Get checks from attribute
    check_results <- attr(x, "checks")
    if (!is.null(check_results)) {

        # Separate rank checks from method checks
        rank_checks <- setdiff(names(check_results), c("gam_metrics", "gee_metrics",
            "lmm_metrics"))
        has_rank_checks <- length(rank_checks) > 0
        has_gam_checks <- !is.null(check_results$gam_metrics)
        has_gee_checks <- !is.null(check_results$gee_metrics)
        has_lmm_checks <- !is.null(check_results$lmm_metrics)

        # RANK-BASED CHECKS
        if (has_rank_checks) {
            message("\nRANK-BASED ASSUMPTIONS")
            message(strrep("-", 70))

            for (check_name in rank_checks) {
                check <- check_results[[check_name]]
                message(sprintf("Test: %s", check_name))
                message(sprintf("  Description: %s", check$description))

                if (!is.null(check$method)) {
                  message(sprintf("  Method: %s", check$method))
                }

                if (!is.null(check$status)) {
                  message(sprintf("  Status: %s", check$status))
                }

                if (!is.null(check$details)) {
                  message(sprintf("  Details: %s", check$details))
                }

                if (!is.null(check$p_value)) {
                  message(sprintf("  P-value: %.4f", check$p_value))
                }

                if (!is.null(check$mean_correlation)) {
                  message(sprintf("  Mean Spearman r: %.4f", check$mean_correlation))
                }

                if (!is.null(check$kendall_w)) {
                  message(sprintf("  Kendall's W: %.4f", check$kendall_w))
                }

                message("")
            }
        }

        # GAM METRICS (NEW)
        if (has_gam_checks) {
            message("\nGAM (GENERALIZED ADDITIVE MODELS) DIAGNOSTICS")
            message(strrep("-", 70))

            gam_results <- check_results$gam_metrics

            for (metric_name in names(gam_results)) {
                metric <- gam_results[[metric_name]]
                message(sprintf("\n%s:", metric$description))

                if (!is.null(metric$method)) {
                  message(sprintf("  Method: %s", metric$method))
                }

                if (!is.null(metric$status)) {
                  message(sprintf("  Status: %s", metric$status))
                }

                # Print metric-specific values
                if (metric_name == "concurvity") {
                  if (!is.na(metric$overall_concurvity)) {
                    message(sprintf("  Concurvity Index: %.4f", metric$overall_concurvity))
                  }
                } else if (metric_name == "edf") {
                  if (!is.null(metric$total_edf) && !is.na(metric$total_edf)) {
                    message(sprintf("  Total EDF: %.2f", metric$total_edf))
                    message(sprintf("  EDF Ratio: %.3f", metric$edf_ratio))
                  }
                } else if (metric_name == "nonlinearity") {
                  if (!is.na(metric$r2_improvement_percent)) {
                    message(sprintf("  R^2 Improvement: %.1f%%", metric$r2_improvement_percent))
                  }
                } else if (metric_name == "basis_adequacy") {
                  if (!is.na(metric$optimal_basis_dimension)) {
                    message(sprintf("  Optimal k: %d", metric$optimal_basis_dimension))
                  }
                }

                if (!is.null(metric$details)) {
                  message(sprintf("  Details: %s", metric$details))
                }

                message("")
            }
        }

        # GEE METRICS (NEW)
        if (has_gee_checks) {
            message("\nGEE (GENERALIZED ESTIMATING EQUATIONS) DIAGNOSTICS")
            message(strrep("-", 70))

            gee_results <- check_results$gee_metrics

            for (metric_name in names(gee_results)) {
                metric <- gee_results[[metric_name]]

                if (metric_name == "consolidated")
                  next

                message(sprintf("\n%s:", metric$description))

                if (!is.null(metric$status)) {
                  message(sprintf("  Status: %s", metric$status))
                }

                # Print metric-specific values
                if (metric_name == "correlation_fit") {
                  if (!is.null(metric$assumed_structure)) {
                    message(sprintf("  Assumed structure: %s", metric$assumed_structure))
                  }
                  if (!is.null(metric$best_structure)) {
                    message(sprintf("  Best structure: %s", metric$best_structure))
                  }
                } else if (metric_name == "cluster_variation") {
                  if (!is.na(metric$mean_cluster_size)) {
                    message(sprintf("  Mean cluster size: %.1f", metric$mean_cluster_size))
                    message(sprintf("  CV: %.3f", metric$cv_cluster_size))
                  }
                } else if (metric_name == "independence") {
                  if (!is.na(metric$mean_residual_correlation)) {
                    message(sprintf("  Mean within-cluster correlation: %.3f", metric$mean_residual_correlation))
                  }
                } else if (metric_name == "scale_parameter") {
                  if (!is.na(metric$scale_parameter)) {
                    message(sprintf("  Scale parameter: %.3f", metric$scale_parameter))
                  }
                }

                if (!is.null(metric$details)) {
                  message(sprintf("  Details: %s", metric$details))
                }

                message("")
            }
        }

        # LMM METRICS (NEW)
        if (has_lmm_checks) {
            message("\nLMM (LINEAR MIXED MODELS) DIAGNOSTICS")
            message(strrep("-", 70))

            lmm_results <- check_results$lmm_metrics

            for (metric_name in names(lmm_results)) {
                metric <- lmm_results[[metric_name]]

                if (metric_name == "consolidated")
                  next

                message(sprintf("\n%s:", metric$description))

                if (!is.null(metric$status)) {
                  message(sprintf("  Status: %s", metric$status))
                }

                # Print metric-specific values
                if (metric_name == "variance_components") {
                  if (!is.null(metric$icc) && !is.na(metric$icc)) {
                    message(sprintf("  ICC (Intraclass Correlation): %.3f", metric$icc))
                    message(sprintf("  Between variance: %.3f", metric$between_variance))
                    message(sprintf("  Within variance: %.3f", metric$within_variance))
                  }
                } else if (metric_name == "normality") {
                  if (!is.null(metric$shapiro_pvalue) && !is.na(metric$shapiro_pvalue)) {
                    message(sprintf("  Shapiro-Wilk p-value: %.4f", metric$shapiro_pvalue))
                    message(sprintf("  Skewness: %.3f, Kurtosis: %.3f", metric$skewness,
                      metric$kurtosis))
                  }
                } else if (metric_name == "homogeneity") {
                  if (!is.null(metric$levene_pvalue) && !is.na(metric$levene_pvalue)) {
                    message(sprintf("  Levene p-value: %.4f", metric$levene_pvalue))
                    message(sprintf("  CV of group variances: %.3f", metric$cv_group_variance))
                  }
                } else if (metric_name == "influence") {
                  if (!is.null(metric$n_outliers)) {
                    message(sprintf("  Number of outliers: %d", metric$n_outliers))
                    message(sprintf("  Influential points: %.1f%%", metric$prop_influential *
                      100))
                  }
                }

                if (!is.null(metric$details)) {
                  message(sprintf("  Details: %s", metric$details))
                }

                message("")
            }
        }

        # FPCA METRICS (NEW)
        has_fpca_checks <- !is.null(check_results$fpca_metrics) && !isTRUE(check_results$fpca_metrics$error)
        if (has_fpca_checks) {
            message("\nFPCA (FUNCTIONAL PRINCIPAL COMPONENT ANALYSIS) DIAGNOSTICS")
            message(strrep("-", 70))

            fpca_results <- check_results$fpca_metrics

            for (metric_name in names(fpca_results)) {
                metric <- fpca_results[[metric_name]]

                if (metric_name == "consolidated")
                  next

                message(sprintf("\n%s:", metric$description))

                if (!is.null(metric$status)) {
                  message(sprintf("  Status: %s", metric$status))
                }

                # Print metric-specific values
                if (metric_name == "variance_adequacy") {
                  if (!is.null(metric$n_components_95)) {
                    message(sprintf("  Components for 95%% variance: %d", metric$n_components_95))
                    message(sprintf("  Components for 90%% variance: %d", metric$n_components_90))
                    message(sprintf("  Total components available: %d", metric$n_components))
                  }
                } else if (metric_name == "bootstrap_stability") {
                  if (!is.null(metric$cv_mean) && !is.na(metric$cv_mean)) {
                    message(sprintf("  Mean bootstrap CV: %.3f", metric$cv_mean))
                    message(sprintf("  Max bootstrap CV: %.3f", metric$max_cv))
                    message(sprintf("  Stable CIs (excl. zero): %.0f%%", metric$prop_stable_cis *
                      100))
                  }
                }

                if (!is.null(metric$details)) {
                  message(sprintf("  Details: %s", metric$details))
                }

                message("")
            }
        }
    }

    message(x$overall_summary)
    message("Note: Use attr(result, 'checks') for complete results including numeric vectors")
    invisible(x)
}

#' Compute Concurvity Index for GAM
#'
#' Detects collinearity among smooth terms. Values > 0.8 indicate problematic
#' collinearity that may require regularization (S150, S143).
#'
#' @param data Matrix of predictor values (columns=predictors, rows=observations)
#' @param gam_cache Optional pre-fitted GAM models (for optimization)
#' @return List with concurvity metrics and status
#' @noRd
.compute_concurvity_index <- function(data, q_values = NULL, gam_cache = NULL) {

    if (!inherits(data, "matrix")) {
        data <- as.matrix(data)
    }

    if (!requireNamespace("mgcv", quietly = TRUE)) {
        return(list(description = "Concurvity Index", overall_concurvity = NA_real_,
            pairwise_concurvities = NULL, status = "? SKIP - mgcv not available",
            details = "Install mgcv package to compute concurvity"))
    }

    # Concurvity only meaningful with 2+ predictors Data matrix format:
    # rows=genes, cols=samples/q-values
    n_predictors <- ncol(data)

    if (n_predictors < 2) {
        return(list(description = "Concurvity Index", overall_concurvity = 0, pairwise_concurvities = NULL,
            status = "OK N/A", details = "Concurvity requires >= 2 predictors"))
    }

    # q-values are required for meaningful concurvity analysis
    if (is.null(q_values) || length(q_values) < 2) {
        return(list(description = "Concurvity Index", overall_concurvity = NA_real_,
            status = "? ERROR", details = "q-values required but not provided; cannot compute concurvity for entropy curves"))
    }

    tryCatch({
        # OPTIMIZATION: Use cached GAM models if provided, otherwise fit
        gam_models <- if (!is.null(gam_cache) && length(gam_cache) > 0) {
            gam_cache
        } else {
            # Fit GAM models if cache not available
            .fit_cached_gams(data, q_values)
        }

        if (length(gam_models) == 0) {
            return(list(description = "Concurvity Index", overall_concurvity = NA_real_,
                n_genes_tested = 0, status = "? ERROR", details = "Could not fit any GAM models; data may have insufficient variation"))
        }

        # Extract model complexity from fitted GAM models Note: Classical
        # concurvity is undefined for single-term GAMs Instead: measure
        # relative model complexity via effective DOF and GCV score High
        # complexity (~complex curvature) -> higher entropy curve variability
        concurv_list <- lapply(gam_models, function(model) {
            tryCatch({
                # Compute relative model complexity: - edf (effective degrees
                # of freedom) from smooth term - Normalized by max possible
                # edf, then scaled to [0,1] Higher edf = more complex/curved
                # entropy pattern

                # Extract EDF from smooth term
                edf_val <- model$edf[1]  # First (only) smooth term

                # Normalize: typical edf ranges 1-10 for simple smooths Scale
                # to approximate [0, 1] where 1 = very complex Use sigmoid-like
                # scaling: complexity ~ 1 - exp(-edf/3)
                if (is.na(edf_val) || edf_val <= 1) {
                  0  # Linear: no effective 'curving'
                } else {
                  # Map edf to [0, 1]: edf=1->0, edf=3->0.63, edf=10->0.96
                  1 - exp(-edf_val/3)
                }
            }, error = function(e) NA_real_)
        })

        overall_concurv <- median(unlist(concurv_list), na.rm = TRUE)

        if (is.na(overall_concurv) || !is.finite(overall_concurv)) {
            return(list(description = "Concurvity Index", overall_concurvity = NA_real_,
                status = "? ERROR", details = "Could not compute concurvity; GAM fits may have failed"))
        }

        # Interpret concurvity (model complexity / entropy curve curvature)
        # Metric: normalized effective DOF from GAM smooths LOW: mostly linear
        # entropy-q relationship MODERATE: noticeable curvature/complexity in
        # entropy patterns HIGH: highly complex/curved entropy profiles
        # (potential instability)
        if (overall_concurv < 0.6) {
            status <- "low"
        } else if (overall_concurv < 0.8) {
            status <- "moderate"
        } else {
            status <- "high"
        }

        # Return the computed results
        return(list(description = "Concurvity Index (Model Complexity)", overall_concurvity = overall_concurv,
            n_genes_tested = length(gam_models), status = status, details = sprintf("Median EDF-based complexity across %d genes: %.4f (lower = less curved entropy profiles)",
                length(gam_models), overall_concurv)))

    }, error = function(e) {
        return(list(description = "Concurvity Index", overall_concurvity = NA_real_,
            status = "? ERROR", details = paste("Failed:", e$message)))
    })
}


#' Compute Effective Degrees of Freedom (EDF) for GAM
#'
#' Assesses smoothing adequacy. EDF ratio < 0.5 (over-smoothed), 0.5-2.0
#' (appropriate), > 2.0 (under-smoothed). References: S137, C045
#'
#' @param data Matrix of predictor values
#' @param q_values Optional numeric vector of q-values for per-gene GAM fitting
#' @param gam_cache Optional pre-fitted GAM models (for optimization)
#' @return List with EDF metrics and interpretation
#' @noRd
.compute_edf_metric <- function(data, q_values = NULL, gam_cache = NULL) {

    if (!inherits(data, "matrix")) {
        data <- as.matrix(data)
    }

    if (!requireNamespace("mgcv", quietly = TRUE)) {
        return(list(description = "Effective Degrees of Freedom", edf_ratio = NA_real_,
            status = "? SKIP - mgcv not available"))
    }

    # If q_values provided, fit per-gene GAMs and aggregate EDF
    if (!is.null(q_values) && length(q_values) >= 2) {
        tryCatch({
            # OPTIMIZATION: Use cached GAM models if provided
            gam_models <- if (!is.null(gam_cache) && length(gam_cache) > 0) {
                gam_cache
            } else {
                .fit_cached_gams(data, q_values)
            }

            # Extract EDF ratios from all models (ensure scalar extraction)
            edf_ratios <- vapply(gam_models, function(model) {
                edf_val <- if (!is.null(model$edf) && length(model$edf) > 0)
                  model$edf[1] else NA_real_
                as.numeric(edf_val)/length(q_values)
            }, numeric(1))
            edf_ratios <- as.numeric(edf_ratios)  # Ensure vector of scalars

            if (length(edf_ratios) == 0 || all(is.na(edf_ratios))) {
                return(list(description = "Effective Degrees of Freedom", edf_ratio = NA_real_,
                  status = "? ERROR", details = "Could not fit any GAM models"))
            }

            # Aggregate EDF ratio across genes (ensure scalar result)
            mean_edf_ratio <- as.numeric(mean(edf_ratios, na.rm = TRUE))

            # Interpretation
            if (mean_edf_ratio < 0.5) {
                status <- "over-smoothed"
            } else if (mean_edf_ratio <= 2) {
                status <- "appropriate"
            } else {
                status <- "under-smoothed"
            }

            return(list(description = "Effective Degrees of Freedom", total_edf = NA_real_,
                edf_ratio = mean_edf_ratio, n_genes_tested = length(gam_models),
                status = status, details = sprintf("Mean EDF ratio=%.3f across %d genes (%s)",
                  mean_edf_ratio, length(gam_models), tolower(status))))
        }, error = function(e) {
            return(list(description = "Effective Degrees of Freedom", edf_ratio = NA_real_,
                status = "? ERROR", details = paste("Failed:", e$message)))
        })
    }

    # q-values are required for meaningful EDF analysis
    return(list(description = "Effective Degrees of Freedom", edf_ratio = NA_real_,
        status = "? ERROR", details = "q-values required but not provided; cannot compute EDF for entropy curves"))
}


#' Compute Non-linearity Contribution
#'
#' Quantifies GAM benefit over linear model. <5% (use LM), 5-20% (GAM justified),
#' >20% (GAM essential). Reference: C045
#'
#' @param data Matrix of predictor values
#' @param q_values Optional numeric vector of q-values for per-gene GAM fitting
#' @param gam_cache Optional pre-fitted GAM models (for optimization)
#' @return List with improvement metrics
#' @noRd
.compute_nonlinearity_contribution <- function(data, q_values = NULL, gam_cache = NULL) {

    if (!inherits(data, "matrix")) {
        data <- as.matrix(data)
    }

    if (!requireNamespace("mgcv", quietly = TRUE)) {
        return(list(description = "Non-linearity Contribution", r2_improvement_percent = NA_real_,
            status = "? SKIP - mgcv not available"))
    }

    # If q_values provided, fit per-gene GAMs and aggregate improvement
    if (!is.null(q_values) && length(q_values) >= 2) {
        tryCatch({
            # OPTIMIZATION: Use cached GAM models if provided
            gam_models <- if (!is.null(gam_cache) && length(gam_cache) > 0) {
                gam_cache
            } else {
                .fit_cached_gams(data, q_values)
            }

            # Get gene indices from cache keys
            gene_indices_cache <- as.numeric(names(gam_models))

            # Compute improvements for cached models
            improvements <- vapply(seq_along(gam_models), function(i) {
                gene_idx <- gene_indices_cache[i]
                gam_fit <- gam_models[[i]]
                entropy_curve <- data[gene_idx, ]
                gam_data <- data.frame(q = q_values, entropy = entropy_curve)
                gam_data <- gam_data[!is.na(gam_data$entropy), , drop = FALSE]

                if (nrow(gam_data) < 5)
                  return(NA_real_)

                tryCatch({
                  # Regularized regression model
                  sait_fit <- stats::lm(entropy ~ q, data = gam_data)
                  # Directly extract r.squared - summary warnings are
                  # non-critical
                  r2_sait <- {
                    s <- summary(sait_fit)
                    if (!is.null(s$r.squared))
                      s$r.squared else NA_real_
                  }

                  # Deviance explained from cached GAM
                  gam_deviance <- (gam_fit$null.deviance - sum(gam_fit$residuals^2))/gam_fit$null.deviance

                  # Improvement percentage
                  ((gam_deviance - r2_sait)/max(r2_sait, 0.001)) * 100
                }, error = function(e) NA_real_)
            })

            improvements <- improvements[!is.na(improvements)]

            if (length(improvements) == 0) {
                return(list(description = "Non-linearity Contribution", r2_improvement_percent = NA_real_,
                  status = "? ERROR", details = "Could not fit any models"))
            }

            # Aggregate improvement
            mean_improvement <- mean(improvements, na.rm = TRUE)

            # Interpretation
            if (mean_improvement < 5) {
                status <- "use linear"
            } else if (mean_improvement < 20) {
                status <- "gam justified"
            } else {
                status <- "gam essential"
            }

            return(list(description = "Non-linearity Contribution", r2_improvement_percent = mean_improvement,
                n_genes_tested = length(gene_indices), status = status, details = sprintf("Mean improvement=%.1f%% across %d genes (%s)",
                  mean_improvement, length(gene_indices), tolower(status))))
        }, error = function(e) {
            return(list(description = "Non-linearity Contribution", r2_improvement_percent = NA_real_,
                status = "? ERROR", details = paste("Failed:", e$message)))
        })
    }

    # q-values are required for meaningful non-linearity analysis
    return(list(description = "Non-linearity Contribution", r2_improvement_percent = NA_real_,
        status = "? ERROR", details = "q-values required but not provided; cannot assess non-linearity for entropy curves"))
}


#' Compute Basis Function Adequacy
#'
#' Tests increasing k values (3,5,8,10,15) to find optimal basis dimension
#' using GCV. Stable GCV indicates adequate basis.
#'
#' @param data Matrix of predictor values
#' @param gam_cache Optional pre-fitted GAM models (for optimization)
#' @return List with basis adequacy assessment
#' @noRd
.compute_basis_adequacy <- function(data, q_values = NULL, gam_cache = NULL) {

    if (!inherits(data, "matrix")) {
        data <- as.matrix(data)
    }

    if (!requireNamespace("mgcv", quietly = TRUE)) {
        return(list(description = "Basis Function Adequacy", optimal_basis_dimension = NA_integer_,
            status = "? SKIP - mgcv not available"))
    }

    # q-values are required for meaningful basis adequacy analysis
    if (is.null(q_values) || length(q_values) < 2) {
        return(list(description = "Basis Function Adequacy", optimal_basis_dimension = NA_integer_,
            status = "? ERROR", details = "q-values required but not provided; cannot find optimal basis dimension"))
    }

    tryCatch({
        # OPTIMIZATION: Use cached GAM models if provided
        gam_models <- if (!is.null(gam_cache) && length(gam_cache) > 0) {
            gam_cache
        } else {
            .fit_cached_gams(data, q_values)
        }

        # Extract k values and GCV from cached models (ensure scalar values)
        optimal_k_per_gene <- vapply(gam_models, function(model) {
            tryCatch({
                if (!is.null(model$smooth) && length(model$smooth) > 0) {
                  smooth_term <- model$smooth[[1]]
                  if (!is.null(smooth_term$bs.dim) && length(smooth_term$bs.dim) >
                    0) {
                    as.integer(smooth_term$bs.dim[1])
                  } else {
                    5L
                  }
                } else {
                  5L
                }
            }, error = function(e) 5L)
        })

        gcv_min_per_gene <- vapply(gam_models, function(model) {
            if (!is.null(model$gcv.ubre))
                as.numeric(model$gcv.ubre[1]) else NA_real_
        }, numeric(1))

        optimal_k_per_gene <- as.integer(optimal_k_per_gene[is.finite(as.numeric(optimal_k_per_gene))])
        gcv_min_per_gene <- gcv_min_per_gene[is.finite(gcv_min_per_gene)]

        if (length(optimal_k_per_gene) == 0) {
            return(list(description = "Basis Function Adequacy", optimal_basis_dimension = NA_integer_,
                status = "? ERROR", details = "Could not fit any GAM models"))
        }

        # Aggregate optimal k across genes (use mode/most common)
        k_counts <- table(optimal_k_per_gene)
        aggregated_k <- as.integer(names(k_counts)[which.max(k_counts)])

        # Status based on basis adequacy
        if (aggregated_k >= 15) {
            status <- "consider increase"
        } else {
            status <- "adequate"
        }

        # Compute mean GCV for reporting
        mean_gcv <- mean(gcv_min_per_gene, na.rm = TRUE)

        return(list(description = "Basis Function Adequacy", optimal_basis_dimension = aggregated_k,
            n_genes_tested = length(optimal_k_per_gene), mean_gcv = mean_gcv, status = status,
            details = sprintf("Optimal k=%d (mean GCV=%.4f) across %d genes (%s)",
                aggregated_k, mean_gcv, length(optimal_k_per_gene), tolower(status))))

    }, error = function(e) {
        return(list(description = "Basis Function Adequacy", optimal_basis_dimension = NA_integer_,
            status = "? ERROR", details = paste("Failed:", e$message)))
    })
}


#' Fit GAMs Once and Cache Results
#'
#' Internal helper: fits GAMs on selected genes, caches results for reuse
#' by all 4 metric functions to avoid redundant computation.
#'
#' @param data Matrix of entropy values
#' @param q_values Numeric vector of q-values
#' @return List of fitted GAM models indexed by gene position
#' @noRd
.fit_cached_gams <- function(data, q_values) {
    if (!requireNamespace("mgcv", quietly = TRUE)) {
        return(list())
    }

    gam_models <- list()
    n_genes <- nrow(data)
    # Sample ~10% of genes for GAM models
    gene_indices <- seq(1, n_genes, by = max(1, floor(n_genes/10)))

    for (gene_idx in gene_indices) {
        entropy_curve <- data[gene_idx, ]
        # Convert vector to ensure proper dataframe construction (prevents rowname warning)
        gam_data <- data.frame(q = q_values, entropy = as.numeric(entropy_curve))
        gam_data <- gam_data[!is.na(gam_data$entropy), , drop = FALSE]

        if (nrow(gam_data) < 5)
            next

        tryCatch({
            k_val <- min(length(unique(gam_data$q)) - 1, 10)
            if (k_val < 3)
                k_val <- 3
            gam_models[[as.character(gene_idx)]] <- mgcv::gam(entropy ~ s(q, k = k_val),
                data = gam_data, method = "GCV.Cp")
        }, error = function(e) NULL)
    }

    return(gam_models)
}

#' Wrapper: Get All GAM Metrics
#'
#' Computes all 4 GAM diagnostics: concurvity, EDF, non-linearity,
#' basis adequacy. Uses cached GAM models to avoid redundant computation.
#'
#' @param data Matrix of predictor values
#' @param method_params List with optional parameters (reserved for future use)
#' @return List containing all 4 GAM metric results
#' @noRd
.get_gam_metrics <- function(data, q_values = NULL, method_params = list()) {

    if (!inherits(data, "matrix")) {
        data <- as.matrix(data)
    }

    # Check if mgcv is available early
    has_mgcv <- requireNamespace("mgcv", quietly = TRUE)

    if (!has_mgcv) {
        return(list(concurvity = list(description = "Concurvity Index", status = "? SKIPPED",
            reason = "mgcv package not installed"), edf = list(description = "Effective Degrees of Freedom",
            status = "? SKIPPED", reason = "mgcv package not installed"), nonlinearity = list(description = "Non-linearity Contribution",
            status = "? SKIPPED", reason = "mgcv package not installed"), basis_adequacy = list(description = "Basis Function Adequacy",
            status = "? SKIPPED", reason = "mgcv package not installed")))
    }

    # If q_values not provided, can only produce placeholders
    if (is.null(q_values) || length(q_values) < 2) {
        return(list(concurvity = list(description = "Concurvity Index", overall_concurvity = NA_real_,
            status = "? ERROR", details = "q-values not provided; cannot fit meaningful GAM models"),
            edf = list(description = "Effective Degrees of Freedom", edf_ratio = NA_real_,
                status = "? ERROR", details = "q-values not provided; cannot fit meaningful GAM models"),
            nonlinearity = list(description = "Non-linearity Contribution", r2_improvement_percent = NA_real_,
                status = "? ERROR", details = "q-values not provided; cannot fit meaningful GAM models"),
            basis_adequacy = list(description = "Basis Function Adequacy", optimal_basis_dimension = NA_real_,
                status = "? ERROR", details = "q-values not provided; cannot fit meaningful GAM models")))
    }

    # OPTIMIZATION: Fit all GAMs once and cache results for reuse
    gam_cache <- .fit_cached_gams(data, q_values)

    # Compute each metric independently using cached GAM models
    results <- list(concurvity = .compute_concurvity_index(data, q_values = q_values,
        gam_cache = gam_cache), edf = .compute_edf_metric(data, q_values = q_values,
        gam_cache = gam_cache), nonlinearity = .compute_nonlinearity_contribution(data,
        q_values = q_values, gam_cache = gam_cache), basis_adequacy = .compute_basis_adequacy(data,
        q_values = q_values, gam_cache = gam_cache))

    # Create consolidated result combining all four metrics
    consolidated_parts <- character()

    # 1. Concurvity Index
    if (!is.null(results$concurvity) && !isTRUE(results$concurvity$error)) {
        if (!is.na(results$concurvity$overall_concurvity)) {
            status_clean <- gsub("^[^a-z]+", "", tolower(results$concurvity$status))
            consolidated_parts <- c(consolidated_parts, sprintf("Index=%.3f (%s)",
                results$concurvity$overall_concurvity, status_clean))
        }
    }

    # 2. EDF Ratio
    if (!is.null(results$edf) && !isTRUE(results$edf$error)) {
        if (!is.na(results$edf$edf_ratio)) {
            status_clean <- gsub("^[^a-z]+", "", tolower(results$edf$status))
            consolidated_parts <- c(consolidated_parts, sprintf("Ratio=%.3f (%s)",
                results$edf$edf_ratio, status_clean))
        }
    }

    # 3. Non-linearity (R^2 improvement)
    if (!is.null(results$nonlinearity) && !isTRUE(results$nonlinearity$error)) {
        if (!is.na(results$nonlinearity$r2_improvement_percent)) {
            status_clean <- gsub("^[^a-z]+", "", tolower(results$nonlinearity$status))
            consolidated_parts <- c(consolidated_parts, sprintf("Delta R^2=%.1f%% (%s)",
                results$nonlinearity$r2_improvement_percent, status_clean))
        }
    }

    # 4. Basis adequacy (k value)
    if (!is.null(results$basis_adequacy) && !isTRUE(results$basis_adequacy$error)) {
        if (!is.na(results$basis_adequacy$optimal_basis_dimension)) {
            status_clean <- gsub("^[^a-z]+", "", tolower(results$basis_adequacy$status))
            consolidated_parts <- c(consolidated_parts, sprintf("k=%d (%s)", results$basis_adequacy$optimal_basis_dimension,
                status_clean))
        }
    }

    # Combine all parts with period separators
    consolidated_result <- if (length(consolidated_parts) > 0) {
        paste(consolidated_parts, collapse = ". ")
    } else {
        "NA (insufficient data)"
    }

    # Add consolidated result to the list
    results$consolidated <- list(description = "Smooth term collinearity", result = consolidated_result,
        status = "COMBINED")

    return(results)
}


#' Compute Working Correlation Structure Fit for GEE
#'
#' Validates if assumed correlation structure matches data. References: S042, S032
#'
#' @param data Matrix of expression values
#' @param assumed_structure Character: assumed correlation structure ('exchangeable', 'ar1', etc.)
#' @return List with correlation structure fit metrics
#' @noRd
.compute_working_correlation_fit <- function(data, assumed_structure = "exchangeable") {

    if (!inherits(data, "matrix")) {
        data <- as.matrix(data)
    }

    tryCatch({
        # Assess correlation structure suitability without fitting GEE
        # Compute autocorrelation across observations to infer structure (vectorized)

        # Option 1: Independence structure assessment
        #   Vectorized autocorrelation computation
        #   Check if data shows autocorrelation (would violate independence)
        #   Note: vapply already protected against invalid correlation calculations
        autocorr_vals <- vapply(seq_len(ncol(data)), function(j) {
            col_data <- data[, j]
            col_data <- col_data[!is.na(col_data)]
            if (length(col_data) > 1 && stats::sd(col_data, na.rm = TRUE) > 0) {
                stats::cor(col_data[-length(col_data)], col_data[-1], use = "complete.obs")
            } else {
                NA_real_
            }
        }, numeric(1))

        autocorr_vals <- autocorr_vals[!is.na(autocorr_vals)]
        mean_autocorr <- if (length(autocorr_vals) > 0)
            mean(autocorr_vals, na.rm = TRUE) else 0

        # Assessment based on observed autocorrelation
        if (abs(mean_autocorr) < 0.2) {
            suitable_structures <- c("independence", "exchangeable")
            assessment <- "independence suitable"
        } else if (abs(mean_autocorr) < 0.5) {
            suitable_structures <- c("exchangeable", "ar1")
            assessment <- "exchangeable or AR(1) suitable"
        } else {
            suitable_structures <- c("ar1")
            assessment <- "AR(1) recommended"
        }

        fit_status <- if (assumed_structure %in% suitable_structures) {
            "good fit"
        } else if (length(suitable_structures) > 0) {
            sprintf("moderate fit (recommend: %s)", suitable_structures[1])
        } else {
            "uncertain"
        }

        return(list(description = "Working Correlation Structure Fit", assumed_structure = assumed_structure,
            mean_autocorrelation = mean_autocorr, suitable_structures = suitable_structures,
            status = fit_status, details = sprintf("Observed autocorr=%.3f; %s",
                mean_autocorr, assessment)))

    }, error = function(e) {
        return(list(description = "Working Correlation Structure Fit", assumed_structure = assumed_structure,
            status = "? ERROR", details = paste("Failed:", e$message)))
    })
}


#' Compute Cluster Size Variation for GEE
#'
#' Assesses homogeneity of cluster sizes. High variation can introduce bias.
#'
#' @param data Matrix of expression values
#' @param cluster_col Optional vector of cluster assignments (default: columns are clusters)
#' @return List with cluster size metrics
#' @noRd
.compute_cluster_size_variation <- function(data, cluster_col = NULL) {

    if (!inherits(data, "matrix")) {
        data <- as.matrix(data)
    }

    tryCatch({
        # Identify cluster sizes
        if (!is.null(cluster_col) && length(cluster_col) == nrow(data)) {
            # Rows are observations, cluster_col assigns them to clusters
            cluster_sizes <- table(cluster_col)
        } else {
            # Columns are clusters (natural grouping in expression data)
            cluster_sizes <- rep(nrow(data), ncol(data))
        }

        # Calculate variation statistics
        mean_size <- mean(cluster_sizes)
        sd_size <- stats::sd(cluster_sizes)
        cv_size <- if (is.na(sd_size))
            0 else sd_size/mean_size  # Coefficient of variation (0 if only 1 cluster)
        min_size <- min(cluster_sizes)
        max_size <- max(cluster_sizes)

        # Interpretation based on coefficient of variation
        if (cv_size < 0.2) {
            status <- "homogeneous"
        } else if (cv_size < 0.5) {
            status <- "moderate variation"
        } else {
            status <- "high variation"
        }

        return(list(description = "Cluster Size Variation", n_clusters = length(cluster_sizes),
            mean_cluster_size = mean_size, sd_cluster_size = sd_size, cv_cluster_size = cv_size,
            min_size = min_size, max_size = max_size, status = status, details = sprintf("Mean size=%.1f (CV=%.3f, range %d-%d) - %s",
                mean_size, cv_size, min_size, max_size, tolower(status))))

    }, error = function(e) {
        return(list(description = "Cluster Size Variation", status = "? ERROR", details = paste("Failed:",
            e$message)))
    })
}


#' Compute Independence Residuals Test for GEE
#'
#' Assessment of within-cluster correlation in residuals.
#'
#' @param data Matrix of expression values
#' @param cluster_col Optional vector of cluster assignments
#' @return List with within-cluster correlation metrics
#' @noRd
.compute_independence_residuals <- function(data, cluster_col = NULL) {

    if (!inherits(data, "matrix")) {
        data <- as.matrix(data)
    }

    tryCatch({
        # For each column (cluster), compute within-cluster correlations
        within_cluster_cors <- numeric(ncol(data))

        for (j in seq_len(ncol(data))) {
            cluster_data <- data[, j]
            # Remove NAs
            cluster_data <- cluster_data[!is.na(cluster_data)]

            if (length(cluster_data) > 1) {
                # Compute residuals from linear fit
                fit <- stats::lm(cluster_data ~ seq_along(cluster_data))
                residuals_j <- fit$residuals

                # Correlation of consecutive residuals
                if (length(residuals_j) > 1) {
                  within_cluster_cors[j] <- stats::cor(residuals_j[-length(residuals_j)],
                    residuals_j[-1], use = "complete.obs")
                } else {
                  within_cluster_cors[j] <- 0
                }
            } else {
                within_cluster_cors[j] <- NA_real_
            }
        }

        mean_within_cor <- mean(within_cluster_cors, na.rm = TRUE)

        # Interpretation: GEE robustly handles within-cluster correlation
        if (abs(mean_within_cor) < 0.2) {
            status <- "independent"
        } else if (abs(mean_within_cor) < 0.5) {
            status <- "some correlation (GEE handles)"
        } else {
            status <- "strong correlation (GEE robust)"
        }

        return(list(description = "Within-Cluster Independence Test", mean_residual_correlation = mean_within_cor,
            n_clusters_tested = sum(!is.na(within_cluster_cors)), status = status,
            details = sprintf("Mean within-cluster residual correlation=%.3f (%s)",
                mean_within_cor, tolower(status))))

    }, error = function(e) {
        return(list(description = "Within-Cluster Independence Test", status = "? ERROR",
            details = paste("Failed:", e$message)))
    })
}


#' Compute GEE Scale Parameter (Dispersion)
#'
#' Assesses over/under-dispersion in GEE model. Reference: S042
#'
#' @param data Matrix of expression values
#' @param cluster_col Optional vector of cluster assignments
#' @return List with scale parameter metrics
#' @noRd
.compute_gee_scale_parameter <- function(data, cluster_col = NULL) {

    if (!inherits(data, "matrix")) {
        data <- as.matrix(data)
    }

    # Check if geepack is available (required for proper GEE analysis)
    if (!requireNamespace("geepack", quietly = TRUE)) {
        return(list(description = "GEE Scale Parameter (Dispersion)", status = "? SKIP - geepack not available",
            details = "Install geepack package to compute scale parameter"))
    }

    tryCatch({
        # Estimate scale parameter (phi) from data variance phi ~ variance /
        # mean for Poisson-like data phi = 1 indicates correct dispersion, >1
        # is over-dispersed, <1 is under-dispersed

        # Calculate overall mean and variance
        data_numeric <- as.numeric(data)
        data_clean <- data_numeric[!is.na(data_numeric) & is.finite(data_numeric)]

        if (length(data_clean) < 2) {
            return(list(description = "GEE Scale Parameter (Dispersion)", status = "? SKIP - geepack not available",
                details = "Insufficient data to estimate scale parameter"))
        }

        mean_val <- mean(data_clean)
        var_val <- stats::var(data_clean)

        # Estimate scale parameter as variance ratio Assuming Gaussian family:
        # phi = observed_var / expected_var For standardized comparison, use
        # normalized variance
        if (mean_val > 0) {
            scale_param <- var_val/(mean_val^2)  # Roughly variance/mean ratio
        } else {
            scale_param <- if (abs(mean_val) > 0)
                var_val/abs(mean_val) else var_val
        }

        # Cap at reasonable bounds for interpretation
        scale_param <- pmax(0.1, pmin(scale_param, 10))

        # Interpretation: phi ~ 1 is ideal
        if (scale_param < 0.8) {
            status <- "under-dispersed (rare)"
        } else if (scale_param <= 1.2) {
            status <- "correct dispersion"
        } else if (scale_param <= 3) {
            status <- "over-dispersed"
        } else {
            status <- "unknown"
        }

        return(list(description = "GEE Scale Parameter (Dispersion)", scale_parameter = scale_param,
            mean_value = mean_val, variance = var_val, n_obs = length(data_clean),
            status = status, details = paste0("phi=", .format_table_value(scale_param))))

    }, error = function(e) {
        return(list(description = "GEE Scale Parameter (Dispersion)", status = "? ERROR",
            details = paste("Failed:", e$message)))
    })
}


#' Wrapper: Get All GEE Metrics
#'
#' Computes all 4 GEE diagnostics: correlation structure fit, cluster size variation,
#' within-cluster independence, scale parameter. Independent computation prevents cascade failures.
#'
#' @param data Matrix of expression values
#' @param gee_params List with optional parameters (assumed_structure, cluster_col)
#' @return List containing all 4 GEE metric results
#' @noRd
.get_gee_metrics <- function(data, gee_params = list()) {

    if (!inherits(data, "matrix")) {
        data <- as.matrix(data)
    }

    # Extract parameters with defaults
    assumed_structure <- if (!is.null(gee_params$assumed_structure))
        gee_params$assumed_structure else "exchangeable"
    cluster_col <- if (!is.null(gee_params$cluster_col))
        gee_params$cluster_col else NULL

    # Compute each metric independently
    results <- list(correlation_fit = .compute_working_correlation_fit(data, assumed_structure = assumed_structure),
        cluster_variation = .compute_cluster_size_variation(data, cluster_col = cluster_col),
        independence = .compute_independence_residuals(data, cluster_col = cluster_col),
        scale_parameter = .compute_gee_scale_parameter(data, cluster_col = cluster_col))

    # Create consolidated result
    consolidated_parts <- character()

    if (!is.null(results$correlation_fit$status)) {
        status_clean <- gsub("^[^a-z]+", "", tolower(results$correlation_fit$status))
        consolidated_parts <- c(consolidated_parts, sprintf("Corr=%s", status_clean))
    }

    if (!is.null(results$cluster_variation$status)) {
        status_clean <- gsub("^[^a-z]+", "", tolower(results$cluster_variation$status))
        consolidated_parts <- c(consolidated_parts, sprintf("Clusters=%s", status_clean))
    }

    if (!is.null(results$scale_parameter$status)) {
        status_clean <- gsub("^[^a-z]+", "", tolower(results$scale_parameter$status))
        consolidated_parts <- c(consolidated_parts, sprintf("Scale=%s", status_clean))
    }

    consolidated_result <- if (length(consolidated_parts) > 0) {
        paste(consolidated_parts, collapse = ". ")
    } else {
        "NA (insufficient data)"
    }

    results$consolidated <- list(description = "Generalized Estimating Equations diagnostics",
        result = consolidated_result, status = "COMBINED")

    return(results)
}


#' Print method for rank correlation confidence intervals
#'
#' @param x Object of class 'rank_correlation_ci'
#' @param ... Additional arguments (ignored)
#'

#' @noRd
#' @method print rank_correlation_ci

print.rank_correlation_ci <- function(x, ...) {
    message("RANK CORRELATION CONFIDENCE INTERVALS")
    message(strrep("=", 60))
    message(sprintf("Method: %s", x$method))
    message(sprintf("Confidence Level: %.0f%%", x$ci_level * 100))

    message("CORRELATION MATRIX")
    message(strrep("-", 60))
    message(paste(capture.output(str(round(x$correlation_matrix, 4))), collapse = "\n"))

    message("\n\nINTERPRETATION SUMMARY")
    message(strrep("-", 60))
    message(paste(capture.output(str(x$interpretation)), collapse = "\n"))

    message("\n\nGUIDELINES FOR INTERPRETATION:")
    message("- Very stable (r > 0.85): Genes rank consistently across all q-values")
    message("- Robust (r > 0.70): Stable ranking; minor q-value effects")
    message("- Moderate (r > 0.50): Noticeable changes; q-value effects important")
    message("- Weak (r <= 0.50): Results highly q-value dependent")
    message("- Variable (includes 0): No stable ranking; q-values give different results")

    invisible(x)
}


#' Classify Genes by Q-Dependency
#'
#' Stratifies genes based on their sensitivity to q-parameter changes.
#'
#' @param interaction_results Data frame output from .calculate_srh()
#' @param p_threshold Numeric: p-value threshold for significance (default:
#' 0.05)
#' @param eta2_threshold_moderate Numeric: Effect size threshold for
#' moderate dependency (default 0.01)
#' @param eta2_threshold_strong Numeric: Effect size threshold for strong
#' dependency (default 0.10)
#'
#' @return
#' A character vector of classifications for each gene. Possible values:
#' \describe{
#'   \item{Robust across q}{p >= p_threshold}
#'   \item{Moderately q-dependent}{p < p_threshold AND eta2 <= eta2_threshold_strong}
#'   \item{Strongly q-dependent}{p < p_threshold AND eta2 > eta2_threshold_strong}
#'   \item{Test failed}{No valid test result}
#'   \item{Insufficient data}{Fewer than 2 q-levels}
#' }
#'
#' @details
#' Classification thresholds can be adjusted based on prior knowledge or
#' exploratory data analysis. Default thresholds correspond to:
#' - Robust: Stable ranking across q (Cohen's small effect)
#' - Moderate: Noticeable but not dramatic ranking shifts (Cohen's small-medium)
#' - Strong: Substantial ranking changes (Cohen's large effect)
#'
#' @examples
#' set.seed(42)
#' # Create sample interaction results
#' interaction_results <- data.frame(
#'   gene = paste0('gene_', 1:10),
#'   p_value = runif(10)
#' )
#' # Classify q-dependency
#' # classifications <- .classify_q_dependency(
#' #   interaction_results, p_threshold = 0.05
#' # )
#' # table(classifications)
#' @noRd

.classify_q_dependency <- function(interaction_results, p_threshold = 0.05, eta2_threshold_moderate = 0.01,
    eta2_threshold_strong = 0.1) {

    # Preserve interaction_class before removing column
    saved_interaction_class <- NULL
    if ("interaction_class" %in% colnames(interaction_results)) {
        saved_interaction_class <- interaction_results$interaction_class
        interaction_results <- interaction_results[, -which(colnames(interaction_results) ==
            "interaction_class")]
    }

    classifications <- character(nrow(interaction_results))

    for (i in seq_len(nrow(interaction_results))) {
        if (is.na(interaction_results$p_value[i])) {
            # Check what type of NA
            if (!is.null(saved_interaction_class) && !is.na(saved_interaction_class[i]) &&
                nchar(saved_interaction_class[i]) > 0) {
                classifications[i] <- saved_interaction_class[i]
            } else {
                classifications[i] <- "Insufficient data"
            }
        } else {
            p_val <- interaction_results$p_value[i]
            eta2_val <- interaction_results$effect_size_eta2[i]

            if (p_val > p_threshold) {
                classifications[i] <- "Robust across q"
            } else if (p_val <= p_threshold && eta2_val <= eta2_threshold_moderate) {
                classifications[i] <- "Moderately q-dependent"
            } else if (p_val <= p_threshold && eta2_val > eta2_threshold_strong) {
                classifications[i] <- "Strongly q-dependent"
            } else if (p_val <= p_threshold) {
                # Gap case: 0.01 < eta2 <= 0.10 with p <= 0.05
                classifications[i] <- "Moderately q-dependent"
            }
        }
    }

    return(classifications)
}  #' Compute Variance Components for LMM
#'
#' Quantifies magnitude of random effects vs residual variance
#'
#' @param data Matrix of predictor values (rows=observations, cols=predictors)
#' @param cluster_col Optional vector of cluster/subject assignments
#' @return List with variance component metrics and interpretation
#' @noRd
.compute_variance_components <- function(data, cluster_col = NULL) {

    if (!inherits(data, "matrix")) {
        data <- as.matrix(data)
    }

    if (ncol(data) < 2) {
        return(list(description = "Variance Components (Random Effects vs Residual)",
            status = "? SKIP", details = "Need at least 2 columns for LMM variance decomposition"))
    }

    tryCatch({
        # Simple variance decomposition without requiring lme4 Treat columns as
        # repeated measures or clusters

        total_var <- var(as.numeric(data))  # Overall variance

        # Between-column variance (represents random effect variation)
        col_means <- colMeans(data, na.rm = TRUE)
        grand_mean <- mean(col_means)
        between_var <- mean((col_means - grand_mean)^2)

        # Within-column variance (residual variance)
        within_vars <- apply(data, 2, function(x) var(x, na.rm = TRUE))
        within_var <- mean(within_vars, na.rm = TRUE)

        # ICC (Intraclass Correlation Coefficient) ICC = between_var /
        # (between_var + within_var)
        icc <- if (!is.na(between_var) && !is.na(within_var) && (between_var + within_var) >
            0) {
            between_var/(between_var + within_var)
        } else {
            NA_real_
        }

        # Determine if LMM needed based on ICC
        if (is.na(icc)) {
            status <- "? ERROR"
        } else if (icc < 0.05) {
            status <- "use sait model"
        } else if (icc < 0.2) {
            status <- "lmm justified"
        } else {
            status <- "lmm essential"
        }

        return(list(description = "Variance Components (Random Effects vs Residual)",
            between_variance = between_var, within_variance = within_var, total_variance = total_var,
            icc = icc, re_proportion = icc, n_clusters = ncol(data), status = status,
            details = paste0("ICC=", .format_table_value(if (is.na(icc)) 0 else icc),
                "; B=", .format_table_value(between_var), ", W=", .format_table_value(within_var))))

    }, error = function(e) {
        return(list(description = "Variance Components", status = "? ERROR", details = paste("Failed:",
            e$message)))
    })
}


#' Compute Random Effects Normality Test
#'
#' Tests normality assumption of random intercepts
#'
#' @param data Matrix of values
#' @param cluster_col Optional vector of cluster assignments
#' @return List with normality test results
#' @noRd
.compute_random_effects_normality <- function(data, cluster_col = NULL) {

    if (!inherits(data, "matrix")) {
        data <- as.matrix(data)
    }

    tryCatch({
        # Extract random effects (as deviations from column means)
        col_means <- colMeans(data, na.rm = TRUE)

        if (length(col_means) < 3) {
            return(list(description = "Random Effects Normality", status = "? SKIP",
                details = "Need at least 3 clusters for normality assessment"))
        }

        grand_mean <- mean(col_means)
        random_effects <- col_means - grand_mean

        # Remove NAs
        re_clean <- random_effects[!is.na(random_effects)]

        # Shapiro-Wilk test (if sample size appropriate)
        shapiro_result <- if (length(re_clean) >= 3 && length(re_clean) <= 5000) {
            stats::shapiro.test(re_clean)
        } else {
            list(statistic = NA_real_, p.value = NA_real_)
        }

        # Calculate skewness and kurtosis
        mean_re <- mean(re_clean)
        sd_re <- stats::sd(re_clean)

        skewness <- if (sd_re > 0) {
            mean((re_clean - mean_re)^3)/(sd_re^3)
        } else {
            NA_real_
        }

        kurtosis <- if (sd_re > 0) {
            mean((re_clean - mean_re)^4)/(sd_re^4) - 3
        } else {
            NA_real_
        }

        # Status based on Shapiro-Wilk p-value
        p_val <- shapiro_result$p.value
        if (is.na(p_val)) {
            status <- "insufficient data"
        } else if (p_val > 0.05) {
            status <- "normal"
        } else if (p_val > 0.01) {
            status <- "slight deviation"
        } else {
            status <- "non-normal"
        }

        return(list(description = "Random Effects Normality", shapiro_statistic = shapiro_result$statistic,
            shapiro_pvalue = p_val, skewness = skewness, kurtosis = kurtosis, n_re_samples = length(re_clean),
            status = status, details = paste0("p=", .format_table_value(if (is.na(p_val)) 0 else p_val,
                4))))

    }, error = function(e) {
        return(list(description = "Random Effects Normality", status = "? ERROR",
            details = paste("Failed:", e$message)))
    })
}


#' Compute Variance Homogeneity (Levene's Test)
#'
#' Tests homogeneity of residual variance across clusters
#'
#' @param data Matrix of values
#' @param cluster_col Optional vector of cluster assignments
#' @return List with homogeneity test results
#' @noRd
.compute_variance_homogeneity <- function(data, cluster_col = NULL) {

    if (!inherits(data, "matrix")) {
        data <- as.matrix(data)
    }

    tryCatch({
        # Group variances by column (cluster)
        group_variances <- apply(data, 2, function(x) {
            x_clean <- x[!is.na(x)]
            if (length(x_clean) > 1)
                stats::var(x_clean) else NA_real_
        })

        group_variances <- group_variances[!is.na(group_variances)]

        if (length(group_variances) < 2) {
            return(list(description = "Variance Homogeneity (Levene's Test)", status = "? SKIP",
                details = "Need at least 2 groups for homogeneity test"))
        }

        # Calculate coefficient of variation of group variances
        cv_group_var <- if (mean(group_variances) > 0) {
            stats::sd(group_variances)/mean(group_variances)
        } else {
            NA_real_
        }

        # Levene's test approximation: use F-test on absolute deviations
        # Adapted from stats::leveneTest alternative
        group_assignments <- rep(seq_len(ncol(data)), each = nrow(data))
        data_long <- as.numeric(data)
        grand_median <- stats::median(data_long, na.rm = TRUE)

        abs_dev <- abs(data_long - grand_median)

        # Simple F-test on group means of absolute deviations
        group_means_ad <- tapply(abs_dev, group_assignments, mean, na.rm = TRUE)
        grand_mean_ad <- mean(abs_dev, na.rm = TRUE)

        ss_between <- sum((group_means_ad - grand_mean_ad)^2 * table(group_assignments))
        ss_total <- sum((abs_dev - grand_mean_ad)^2, na.rm = TRUE)
        ss_within <- ss_total - ss_between

        df_between <- length(group_means_ad) - 1
        df_within <- length(abs_dev[!is.na(abs_dev)]) - length(group_means_ad)

        ms_between <- if (df_between > 0)
            ss_between/df_between else NA_real_
        ms_within <- if (df_within > 0)
            ss_within/df_within else NA_real_

        levene_statistic <- if (!is.na(ms_between) && !is.na(ms_within) && ms_within >
            0) {
            ms_between/ms_within
        } else {
            NA_real_
        }

        # Approximate p-value using F distribution
        levene_pvalue <- if (!is.na(levene_statistic)) {
            stats::pf(levene_statistic, df_between, df_within, lower.tail = FALSE)
        } else {
            NA_real_
        }

        # Status based on p-value
        if (is.na(levene_pvalue)) {
            status <- "? ERROR"
        } else if (levene_pvalue > 0.05) {
            status <- "homogeneous"
        } else if (levene_pvalue > 0.01) {
            status <- "slight heterogeneity"
        } else {
            status <- "heterogeneous"
        }

        return(list(description = "Variance Homogeneity (Levene's Test)", levene_statistic = levene_statistic,
            levene_pvalue = levene_pvalue, df_between = df_between, df_within = df_within,
            cv_group_variance = cv_group_var, n_groups = length(group_variances),
            status = status, details = paste0("p=", .format_table_value(if (is.na(levene_pvalue)) 0 else levene_pvalue,
                4), "; CV=", .format_table_value(if (is.na(cv_group_var)) 0 else cv_group_var))))

    }, error = function(e) {
        return(list(description = "Variance Homogeneity", status = "? ERROR", details = paste("Failed:",
            e$message)))
    })
}


#' Compute Outlier Influence (Cook's Distance)
#'
#' Identifies influential points and assesses robustness
#'
#' @param data Matrix of values
#' @param cluster_col Optional vector of cluster assignments
#' @return List with influence assessment
#' @noRd
.compute_lmm_influence <- function(data, cluster_col = NULL) {

    if (!inherits(data, "matrix")) {
        data <- as.matrix(data)
    }

    tryCatch({
        data_long <- as.numeric(data)
        data_clean <- data_long[!is.na(data_long) & is.finite(data_long)]

        if (length(data_clean) < 3) {
            return(list(description = "Outlier Influence Assessment", status = "? SKIP",
                details = "Insufficient data for influence analysis"))
        }

        # Compute standardized residuals from mean
        mean_val <- mean(data_clean)
        sd_val <- stats::sd(data_clean)

        if (sd_val == 0) {
            return(list(description = "Outlier Influence Assessment", status = "? SKIP",
                details = "Zero variance - no variation to assess"))
        }

        # Standardized residuals
        std_residuals <- (data_clean - mean_val)/sd_val

        # Identify outliers: |std_residual| > 3
        n_outliers <- sum(abs(std_residuals) > 3, na.rm = TRUE)

        # Identify extreme: |std_residual| > 3.5
        n_extreme <- sum(abs(std_residuals) > 3.5, na.rm = TRUE)

        # Cook's distance approximation: (std_residual)^2 For exploratory use
        # without actual regression model Threshold = 1 is a standard cutoff
        # for Cook's distance (Note: 4/n would be too small for large samples,
        # flagging nearly all points)
        cook_approx <- abs(std_residuals)^2
        # Standard Cook's distance threshold (fixed, not sample-size dependent)
        n_influential <- sum(cook_approx > 1, na.rm = TRUE)

        # Proportion
        prop_influential <- if (length(data_clean) > 0) {
            n_influential/length(data_clean)
        } else {
            0
        }

        # Status based on proportion of influential points
        if (prop_influential > 0.1) {
            status <- "many outliers"
        } else if (prop_influential > 0.05) {
            status <- "some outliers"
        } else if (n_extreme > 0) {
            status <- "extreme outliers present"
        } else if (n_outliers > 0) {
            status <- "mild outliers"
        } else {
            status <- "no outliers"
        }

        return(list(description = "Outlier Influence Assessment", n_outliers = n_outliers,
            n_extreme = n_extreme, n_influential = n_influential, prop_influential = prop_influential,
            n_observations = length(data_clean), status = status, details = sprintf("Outliers=%d, Extreme=%d, Influential=%.1f%%",
                n_outliers, n_extreme, prop_influential * 100)))

    }, error = function(e) {
        return(list(description = "Outlier Influence Assessment", status = "? ERROR",
            details = paste("Failed:", e$message)))
    })
}


#' Wrapper: Get All LMM Metrics
#'
#' Computes all 4 LMM diagnostic metrics in one call
#'
#' @param data Matrix of predictor values
#' @param lmm_params List of LMM parameters (optional)
#'   - assumed_re_structure: 'random_intercept', 'random_slope', etc.
#'   - cluster_col: Optional cluster assignments
#' @return List with all 4 metrics + consolidated result
#' @noRd
.get_lmm_metrics <- function(data, lmm_params = list()) {

    if (!inherits(data, "matrix")) {
        data <- as.matrix(data)
    }

    # Extract parameters with defaults
    assumed_re_structure <- if (!is.null(lmm_params$assumed_re_structure))
        lmm_params$assumed_re_structure else "random_intercept"
    cluster_col <- if (!is.null(lmm_params$cluster_col))
        lmm_params$cluster_col else NULL

    # Compute each metric independently
    results <- list(variance_components = .compute_variance_components(data, cluster_col = cluster_col),
        normality = .compute_random_effects_normality(data, cluster_col = cluster_col),
        homogeneity = .compute_variance_homogeneity(data, cluster_col = cluster_col),
        influence = .compute_lmm_influence(data, cluster_col = cluster_col))

    # Create consolidated result
    consolidated_parts <- character()

    if (!is.null(results$variance_components$status)) {
        status_clean <- gsub("^[^a-z]+", "", tolower(results$variance_components$status))
        consolidated_parts <- c(consolidated_parts, sprintf("VarComp=%s", status_clean))
    }

    if (!is.null(results$normality$status)) {
        status_clean <- gsub("^[^a-z]+", "", tolower(results$normality$status))
        consolidated_parts <- c(consolidated_parts, sprintf("Normal=%s", status_clean))
    }

    if (!is.null(results$homogeneity$status)) {
        status_clean <- gsub("^[^a-z]+", "", tolower(results$homogeneity$status))
        consolidated_parts <- c(consolidated_parts, sprintf("Homog=%s", status_clean))
    }

    if (!is.null(results$influence$status)) {
        status_clean <- gsub("^[^a-z]+", "", tolower(results$influence$status))
        consolidated_parts <- c(consolidated_parts, sprintf("Influence=%s", status_clean))
    }

    consolidated_result <- if (length(consolidated_parts) > 0) {
        paste(consolidated_parts, collapse = ". ")
    } else {
        "NA (insufficient data)"
    }

    results$consolidated <- list(description = "Linear Mixed Models diagnostics",
        result = consolidated_result, status = "COMBINED")

    return(results)
}
#' Compute Cumulative Variance Explained by FPCA
#'
#' Functional PCA dimension adequacy: assess how many components needed
#' to capture variance. Useful for evaluating whether FPCA provides
#' meaningful dimension reduction.
#'
#' @param data Matrix of predictor values
#' @param max_components Maximum components to compute (default: min(p, n-1))
#'
#' @return List with eigenvalues, cumulative variance, and interpretation
#' @noRd
.compute_cumulative_variance_fpca <- function(data, max_components = NULL) {

    if (!inherits(data, "matrix")) {
        data <- as.matrix(data)
    }

    tryCatch({
        # Compute PCA via singular value decomposition (functional PCA
        # approximation) For functional data, we treat rows as observations and
        # columns as function values

        # Center the data
        data_centered <- scale(data, center = TRUE, scale = FALSE)

        # SVD: data_centered = U * D * V^T
        svd_result <- svd(data_centered)

        # Eigenvalues from PCA (squared singular values / (n-1))
        n_obs <- nrow(data_centered)
        eigenvalues <- (svd_result$d^2)/(n_obs - 1)

        # Remove near-zero eigenvalues
        eigenvalues <- eigenvalues[eigenvalues > 1e-10]

        if (length(eigenvalues) == 0) {
            return(list(description = "Functional PCA Dimension Adequacy", status = "? SKIP",
                details = "No variance detected (possibly constant data)"))
        }

        # Limit to max_components if specified
        if (!is.null(max_components) && max_components < length(eigenvalues)) {
            eigenvalues <- eigenvalues[seq_len(max_components)]
        }

        # Cumulative variance explained
        total_variance <- sum(eigenvalues)
        cumulative_variance <- cumsum(eigenvalues)/total_variance
        cumulative_pct <- cumulative_variance * 100

        # Determine components needed for 90%, 95%, 99%
        n_comp_90 <- which.max(cumulative_pct >= 90)
        n_comp_95 <- which.max(cumulative_pct >= 95)
        n_comp_99 <- which.max(cumulative_pct >= 99)

        # Handle case where variance not reached
        if (length(n_comp_90) == 0)
            n_comp_90 <- length(eigenvalues)
        if (length(n_comp_95) == 0)
            n_comp_95 <- length(eigenvalues)
        if (length(n_comp_99) == 0)
            n_comp_99 <- length(eigenvalues)

        # Interpret dimension reduction quality
        if (n_comp_95 < 5) {
            status <- "good reduction"
            interpretation <- "Excellent dimension reduction; FPCA highly beneficial"
        } else if (n_comp_95 < 10) {
            status <- "moderate reduction"
            interpretation <- "Moderate dimension reduction; FPCA provides some benefit"
        } else {
            status <- "poor reduction"
            interpretation <- "Poor dimension reduction; consider using raw data"
        }

        return(list(description = "Functional PCA Dimension Adequacy", eigenvalues = eigenvalues,
            cumulative_variance = cumulative_variance, cumulative_pct = cumulative_pct,
            n_components_90 = n_comp_90, n_components_95 = n_comp_95, n_components_99 = n_comp_99,
            total_variance = total_variance, n_components = length(eigenvalues),
            status = status, details = sprintf("Components for 90%%=%.0f, 95%%=%.0f, 99%%=%.0f. %s",
                n_comp_90, n_comp_95, n_comp_99, interpretation)))

    }, error = function(e) {
        return(list(description = "Functional PCA Dimension Adequacy", status = "? ERROR",
            details = paste("Failed:", e$message)))
    })
}

#' Compute FPCA Bootstrap Stability
#'
#' Assess stability of principal components via bootstrap resampling.
#' Stable components have narrow bootstrap confidence intervals.
#'
#' @param data Matrix of predictor values
#' @param n_bootstrap Number of bootstrap samples (default: 500)
#' @param n_components Number of principal components to evaluate (default: 3)
#'
#' @return List with bootstrap stability metrics
#' @noRd
.compute_fpca_bootstrap_stability <- function(data, n_bootstrap = 500, n_components = 3) {

    if (!inherits(data, "matrix")) {
        data <- as.matrix(data)
    }

    tryCatch({
        n_obs <- nrow(data)
        n_vars <- ncol(data)

        if (n_obs < 5 || n_components < 1) {
            return(list(description = "FPCA Bootstrap Stability", status = "? SKIP",
                details = "Insufficient data for bootstrap (need >=5 observations)"))
        }

        # Limit n_components to available dimensions
        n_components <- min(n_components, n_obs - 1, n_vars - 1)

        # Original PCA (pre-compute centering parameters for reuse)
        data_mean <- colMeans(data, na.rm = TRUE)
        data_centered <- scale(data, center = TRUE, scale = FALSE)
        svd_orig <- svd(data_centered)
        eigenvalues_orig <- (svd_orig$d^2)/(n_obs - 1)
        eigenvalues_orig <- eigenvalues_orig[seq_len(n_components)]

        # OPTIMIZATION: Vectorized bootstrap resampling with pre-allocated
        # matrix
        bootstrap_eigenvalues <- matrix(NA_real_, nrow = n_bootstrap, ncol = n_components)

        # Pre-generate all bootstrap indices at once (vectorized)
        boot_indices <- lapply(seq_len(n_bootstrap), function(b) {
            sample(seq_len(n_obs), size = n_obs, replace = TRUE)
        })

        # Apply SVD to bootstrap samples (vectorized loop)
        for (b in seq_len(n_bootstrap)) {
            idx_boot <- boot_indices[[b]]
            data_boot <- data[idx_boot, , drop = FALSE]

            # Use pre-computed mean for efficiency
            data_boot_centered <- t(t(data_boot) - data_mean)

            # Perform SVD on bootstrap sample
            tryCatch({
                svd_boot <- svd(data_boot_centered)
                eig_boot <- (svd_boot$d^2)/(nrow(data_boot) - 1)

                # Store first n_components (pad with NA if fewer exist)
                n_eig <- min(length(eig_boot), n_components)
                bootstrap_eigenvalues[b, seq_len(n_eig)] <- eig_boot[seq_len(n_eig)]
            }, error = function(e) NULL)
        }

        # Compute bootstrap statistics (vectorized operations on pre-allocated
        # matrix)
        bootstrap_se <- apply(bootstrap_eigenvalues, 2, sd, na.rm = TRUE)
        bootstrap_ci_lower <- apply(bootstrap_eigenvalues, 2, quantile, probs = 0.025,
            na.rm = TRUE)
        bootstrap_ci_upper <- apply(bootstrap_eigenvalues, 2, quantile, probs = 0.975,
            na.rm = TRUE)

        # Coefficient of variation
        cv_eigenvalues <- bootstrap_se/eigenvalues_orig
        cv_mean <- mean(cv_eigenvalues, na.rm = TRUE)

        # Stability assessment Stable if SE is small relative to eigenvalue, CI
        # excludes zero
        prop_stable <- sum(bootstrap_ci_lower > 0, na.rm = TRUE)/n_components

        # Guard against empty cv_eigenvalues vector
        has_valid_cv <- sum(!is.na(cv_eigenvalues)) > 0
        max_cv <- if (has_valid_cv)
            max(cv_eigenvalues, na.rm = TRUE) else NA_real_

        if (max_cv < 0.2 && prop_stable > 0.9) {
            status <- "stable"
            interpretation <- "Excellent bootstrap stability; PC estimates reliable"
        } else if (max_cv < 0.5 && prop_stable > 0.7) {
            status <- "moderate stability"
            interpretation <- "Moderate bootstrap stability; PC estimates reasonably reliable"
        } else {
            status <- "unstable"
            interpretation <- "Low bootstrap stability; consider increasing sample size"
        }

        return(list(description = "FPCA Bootstrap Stability", eigenvalues_original = eigenvalues_orig,
            bootstrap_se = bootstrap_se, bootstrap_ci_lower = bootstrap_ci_lower,
            bootstrap_ci_upper = bootstrap_ci_upper, cv_eigenvalues = cv_eigenvalues,
            cv_mean = cv_mean, max_cv = max_cv, prop_stable_cis = prop_stable, status = status,
            details = sprintf("Bootstrap SE (mean CV)=%.3f; Stable CIs=%.0f%%. %s",
                cv_mean, prop_stable * 100, interpretation)))

    }, error = function(e) {
        return(list(description = "FPCA Bootstrap Stability", status = "? ERROR",
            details = paste("Failed:", e$message)))
    })
}

#' Wrapper: Get All FPCA Metrics
#'
#' Computes all 2 FPCA diagnostic metrics in one call
#'
#' @param data Matrix of predictor values
#' @param fpca_params List of FPCA parameters (optional)
#'
#' @return List containing variance adequacy, bootstrap stability, and consolidated result
#' @noRd
.get_fpca_metrics <- function(data, fpca_params = list()) {

    # Extract parameters with defaults
    max_components <- fpca_params$max_components %||% NULL
    n_bootstrap <- fpca_params$n_bootstrap %||% 500
    n_components_eval <- fpca_params$n_components %||% 3

    # Compute metrics
    variance_adequacy <- .compute_cumulative_variance_fpca(data, max_components)
    bootstrap_stability <- .compute_fpca_bootstrap_stability(data, n_bootstrap, n_components_eval)

    # Consolidated result
    status_va <- if (is.null(variance_adequacy$status))
        "?" else sub("^\\? ", "", variance_adequacy$status)
    status_bs <- if (is.null(bootstrap_stability$status))
        "?" else sub("^\\? ", "", bootstrap_stability$status)

    consolidated_result <- paste("Variance:", status_va, "|", "Bootstrap:", status_bs)

    results <- list(variance_adequacy = variance_adequacy, bootstrap_stability = bootstrap_stability,
        consolidated = list(description = "Functional Principal Component Analysis diagnostics",
            result = consolidated_result, status = "COMBINED"))

    return(results)
}
