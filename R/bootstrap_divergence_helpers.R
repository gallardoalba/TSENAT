# ============================================================================
# HELPER FUNCTIONS
# ============================================================================

#' Compute Tsallis Divergence Between Two Distributions
#'
#' @param p Numeric vector; first probability distribution (must sum to 1)
#' @param r Numeric vector; second probability distribution (must sum to 1)
#' @param q Numeric; Tsallis parameter (q > 0)
#' @param log_base Base of logarithm (default: exp(1))
#' @param norm Logical; normalize by theoretical maximum if TRUE
#'
#' @return Numeric scalar; Tsallis divergence value
#'

#' @noRd
.compute_tsallis_divergence <- function(p, r, q, log_base = exp(1), norm = FALSE) {

    # Validate lengths
    if (length(p) != length(r)) {
        return(NA_real_)
    }

    # Handle edge cases
    if (abs(q - 1) < 1e-10) {
        # KL divergence (q -> 1 limit), natural log, then log_base conversion
        # Single definition, no double log_base division).
        idx <- p > 0
        if (sum(idx) == 0)
            return(NA_real_)
        divergence <- sum(p[idx] * log(p[idx]/r[idx]))
        divergence <- divergence/log(log_base)
    } else if (q > 0) {
        # General Tsallis divergence (Furuichi 2006)
        # D_q(p||r) = (sum(p^q * r^(1-q)) - 1) / (q - 1)
        # Canonical SIGNED form matching .tsallis_divergence_scalar/
        # .tsallis_divergence_vector and the C++ kernel: NO abs(). Only tiny
        # negative numerical roundoff is clamped to zero.

        p_power <- p^q
        r_power <- r^(1 - q)

        # Check for numerical issues (inf, nan, underflow)
        if (any(is.nan(p_power)) || any(is.infinite(p_power)) || any(is.nan(r_power)) ||
            any(is.infinite(r_power))) {
            # Log-space computation for numerical stability when q is far from
            # 1
            log_p_power <- q * log(pmax(p, 1e-10))
            log_r_power <- (1 - q) * log(pmax(r, 1e-10))
            sum_term <- sum(exp(log_p_power + log_r_power), na.rm = TRUE)
        } else {
            sum_term <- sum(p_power * r_power, na.rm = TRUE)
        }

        # Apply Furuichi formula
        divergence <- (sum_term - 1)/(q - 1)
        if (divergence < 0 && divergence > -1e-12)
            divergence <- 0
    } else {
        # Invalid q value
        return(NA_real_)
    }

    # Handle invalid results
    if (is.nan(divergence) || !is.finite(divergence)) {
        return(NA_real_)
    }

    # Normalize if requested
    # Max divergence normalization depends on q.
    # - For q≈1 (KL limit): max = log(n), using Shannon-style max. Correct.
    # - For q>1: Tsallis divergence is bounded by 1/(q-1) when p and r are
    #   maximally different (one element concentrates in p, another in r).
    # - For 0<q<1: Tsallis divergence is unbounded — normalization skipped.
    # Previously used log(n) for all q, which is only valid for the KL limit.
    if (norm && divergence > 0) {
        if (abs(q - 1) < 1e-10) {
            max_div <- log(length(p), base = log_base)
        } else if (q > 1) {
            max_div <- 1 / (q - 1)
        } else {
            max_div <- NA_real_  # q<1: unbounded, skip normalization
        }
        if (!is.na(max_div) && is.finite(max_div) && max_div > 0) {
            divergence <- divergence / max_div
        }
    }

    return(as.numeric(divergence))
}


#' Bias-Corrected and Accelerated (BCa) Confidence Intervals
#'
#' Compute BCa CIs using jackknife for acceleration and bias correction
#'
#' @param boot_dist Numeric vector; bootstrap distribution
#' @param theta_hat Numeric; point estimate (from original data)
#' @param alpha Numeric; significance level (1 - ci)
#'
#' @return List with \code{lower} and \code{upper} CI bounds
#'

#' @noRd
.bca_ci <- function(boot_dist, theta_hat, alpha, jackknife_estimates = NULL) {

    n <- length(boot_dist)

    # If theta_hat is infinite or the bootstrap distribution is degenerate,
    # fall back to percentile method. Use scale-relative threshold (sd/|mean|)
    # instead of absolute sd < 1e-10 to handle both small and large entropy scales.
    theta_abs <- abs(theta_hat)
    sd_boot <- sd(boot_dist, na.rm = TRUE)
    is_degenerate <- if (theta_abs > 1e-10) {
        sd_boot / theta_abs < 1e-10
    } else {
        sd_boot < 1e-10
    }
    if (!is.finite(theta_hat) || is_degenerate) {
        lower <- stats::quantile(boot_dist, alpha/2, na.rm = TRUE)
        upper <- stats::quantile(boot_dist, 1 - alpha/2, na.rm = TRUE)
        return(list(lower = as.numeric(lower), upper = as.numeric(upper)))
    }

    # Use strict < comparison with +0.5/B padding.
    # Clamp to avoid qnorm(0) = -Inf and qnorm(1) = Inf/NaN.
    prop_less <- (sum(boot_dist < theta_hat, na.rm = TRUE) + 0.5) / n
    prop_less <- pmax(0.001, pmin(0.999, prop_less))
    z0 <- stats::qnorm(prop_less)

    # Handle case where z0 is infinite
    if (!is.finite(z0)) {
        z0 <- 0
    }

    # BCa acceleration MUST be computed from true leave-one-out
    # jackknife on the ORIGINAL data, not from the bootstrap distribution.
    # When jackknife_estimates is provided (e.g., from divergence jackknife),
    # use those. Otherwise fall back to bootstrap-based acceleration with a
    # warning for backwards compatibility.
    if (!is.null(jackknife_estimates) && length(jackknife_estimates) >= 3) {
        # True jackknife: use leave-one-out estimates from original data
        theta_jack <- jackknife_estimates[is.finite(jackknife_estimates)]
        if (length(theta_jack) < 3) {
            acceleration <- 0
        } else {
            theta_bar <- mean(theta_jack)
            deviations <- theta_bar - theta_jack
            numerator <- sum(deviations^3)
            denom_base <- sum(deviations^2)
            denominator <- 6 * (denom_base)^(3/2)
            acceleration <- if (denominator > 1e-10 && is.finite(denominator))
                numerator / denominator else 0
        }
    } else {
        # Legacy: bootstrap-based acceleration (biased toward 0 for large B)
        # Reference: Efron & Tibshirani (1993), "An Introduction to the Bootstrap", Ch. 14
        theta_bar <- mean(boot_dist, na.rm = TRUE)
        total_sum <- sum(boot_dist, na.rm = TRUE)
        n_valid <- sum(!is.na(boot_dist))

        if (n_valid > 1) {
            theta_jack <- (total_sum - boot_dist)/(n_valid - 1)
        } else {
            theta_jack <- rep(boot_dist[!is.na(boot_dist)][1], length(boot_dist))
        }

        deviations <- theta_bar - theta_jack
        numerator <- sum(deviations^3, na.rm = TRUE)
        denom_base <- sum(deviations^2, na.rm = TRUE)
        denominator <- 6 * (denom_base)^(3/2)

        acceleration <- if (denominator > 1e-10 && is.finite(denominator))
            numerator/denominator else 0
    }

    # Handle invalid acceleration
    if (!is.finite(acceleration)) {
        acceleration <- 0
    }

    # Adjusted quantiles
    z_alpha_lower <- stats::qnorm(alpha/2)
    z_alpha_upper <- stats::qnorm(1 - alpha/2)

    # Calculate adjusted probabilities, handling division by zero
    denom_lower <- 1 - acceleration * (z0 + z_alpha_lower)
    denom_upper <- 1 - acceleration * (z0 + z_alpha_upper)

    if (abs(denom_lower) < 1e-10) {
        # Fall back to unadjusted if denominator near zero
        p_lower <- alpha/2
    } else {
        p_lower <- stats::pnorm(z0 + (z0 + z_alpha_lower)/denom_lower)
    }

    if (abs(denom_upper) < 1e-10) {
        # Fall back to unadjusted if denominator near zero
        p_upper <- 1 - alpha/2
    } else {
        p_upper <- stats::pnorm(z0 + (z0 + z_alpha_upper)/denom_upper)
    }

    # Clamp to valid quantile range
    p_lower <- pmax(0.001, pmin(0.999, p_lower))
    p_upper <- pmax(0.001, pmin(0.999, p_upper))

    lower <- stats::quantile(boot_dist, p_lower, na.rm = TRUE)
    upper <- stats::quantile(boot_dist, p_upper, na.rm = TRUE)

    list(lower = as.numeric(lower), upper = as.numeric(upper))
}


# ============================================================================
# S3 METHODS FOR PRINT AND SUMMARY
# ============================================================================

#' @noRd
#' @exportS3Method
print.tsenat_divergence_bootstrap_ci <- function(x, ...) {
    ci_pct <- if (!is.null(x$ci_level)) round(x$ci_level * 100, 1) else 95
    message("Tsallis Divergence Bootstrap Confidence Interval")
    if (!is.null(x$gene_name)) message("Gene: ", x$gene_name)
    if (!is.null(x$q)) message("q-parameter: ", x$q)
    message("Point estimate: ", sprintf("%.6f", x$estimate))
    message(ci_pct, "% CI: [", sprintf("%.6f", x$lower_ci), ", ", sprintf("%.6f", x$upper_ci), "]")
    message("Method: ", x$method %||% "N/A", " | Replicates: ", x$nboot)
    invisible(x)
}

#' Summary method for divergence bootstrap CI results
#'
#' Provides summary statistics and diagnostic information for divergence
#' bootstrap confidence interval objects.
#'
#' @return
#' Invisibly returns the object itself.
#' Prints to console: q-values, sample information, and CI statistics.
#'
# ============================================================================
# CONSOLIDATED BOOTSTRAP UTILITIES (moved from other files, March 2026)
# ============================================================================

#' Internal: Aggregate Bootstrap Confidence Intervals
#'
#' Consolidates bootstrap CI results from SummarizedExperiment format into
#' a long-format data.frame suitable for visualization.
#'
#' @param se SummarizedExperiment with ci_lower and ci_upper assays
#' @param long data.frame with columns: q, group, and other metadata
#'
#' @return data.frame with aggregated CI values across samples and q-values
#'
#' @noRd
# From compute_stats.R: Aggregate bootstrap CIs
.bootstrap_aggregate_ci <- function(se, long) {

    ci_lower_mat <- SummarizedExperiment::assay(se, "ci_lower")
    ci_upper_mat <- SummarizedExperiment::assay(se, "ci_upper")

    sample_names <- colnames(ci_lower_mat)
    if (is.null(sample_names)) {
        sample_names <- paste0("Sample", seq_len(ncol(ci_lower_mat)))
    }

    groups <- unique(sort(long$group))
    # BUGFIX: long$q is a factor - must convert to numeric!
    unique_q <- sort(as.numeric(as.character(unique(long$q))))

    plot_df <- data.frame(q = numeric(), median = numeric(), ci_lower = numeric(),
        ci_upper = numeric(), group = character(), stringsAsFactors = FALSE)

    # BUGFIX (April 2026): Pair each gene's median with its corresponding CI
    # bounds Previously: computed median from ALL genes, but CI bounds from ALL
    # genes separately Result: CI displacement when genes have different
    # variability Fix: aggregate PER GENE first, then combine across genes
    # using median of medians

    for (group_val in groups) {
        for (q_val in unique_q) {
            # Convert q_val to numeric for comparison with long$q
            matching_q_val <- as.numeric(as.character(q_val))

            group_q_data <- long %>%
                dplyr::filter(group == group_val, as.numeric(as.character(q)) ==
                  matching_q_val)

            if (nrow(group_q_data) > 0) {
                # CORRECTED: Group by Gene first to keep each gene's median
                # with its CI bounds
                if ("Gene" %in% colnames(group_q_data)) {
                  genes_in_group <- unique(as.character(group_q_data$Gene))

                  all_gene_medians <- c()
                  all_gene_ci_lower <- c()
                  all_gene_ci_upper <- c()

                  # Process each gene separately to maintain pairing
                  for (gene_val in genes_in_group) {
                    gene_data <- group_q_data[as.character(group_q_data$Gene) ==
                      gene_val, ]

                    # Get median tsallis for this SPECIFIC gene
                    gene_median <- median(gene_data$tsallis, na.rm = TRUE)
                    all_gene_medians <- c(all_gene_medians, gene_median)

                    # Get samples for this gene
                    gene_samples <- unique(gene_data$sample)
                    gene_ci_lower <- c()
                    gene_ci_upper <- c()

                    for (samp in gene_samples) {
                      # Construct the expected column name Format q-value with
                      # 3 decimal places to match colname format
                      q_formatted <- formatC(matching_q_val, format = "f", digits = 3)
                      expected_col_name <- paste0(samp, "_q=", q_formatted)

                      samp_idx <- which(sample_names == expected_col_name)

                      if (length(samp_idx) > 0) {
                        # Get CI for this SPECIFIC gene
                        gene_idx <- match(gene_val, rownames(se))
                        if (!is.na(gene_idx)) {
                          gene_ci_lower <- c(gene_ci_lower, ci_lower_mat[gene_idx,
                            samp_idx[1]])
                          gene_ci_upper <- c(gene_ci_upper, ci_upper_mat[gene_idx,
                            samp_idx[1]])
                        }
                      }
                    }

                    # Aggregate CI for this gene across samples
                    if (length(gene_ci_lower) > 0) {
                      all_gene_ci_lower <- c(all_gene_ci_lower, median(gene_ci_lower,
                        na.rm = TRUE))
                      all_gene_ci_upper <- c(all_gene_ci_upper, median(gene_ci_upper,
                        na.rm = TRUE))
                    }
                  }

                  # Final aggregation: median of per-gene medians and CIs This
                  # ensures alignment: each point estimate has CI bounds from
                  # the same gene
                  if (length(all_gene_medians) > 0) {
                    median_val <- median(all_gene_medians, na.rm = TRUE)
                    ci_lower_final <- median(all_gene_ci_lower, na.rm = TRUE)
                    ci_upper_final <- median(all_gene_ci_upper, na.rm = TRUE)
                  } else {
                    median_val <- NA_real_
                    ci_lower_final <- NA_real_
                    ci_upper_final <- NA_real_
                  }
                } else {
                  # Fallback for data without Gene column (original behavior)
                  median_val <- median(group_q_data$tsallis, na.rm = TRUE)

                  group_samples <- unique(group_q_data$sample)
                  all_ci_lower <- c()
                  all_ci_upper <- c()

                  for (samp in group_samples) {
                    q_formatted <- formatC(matching_q_val, format = "f", digits = 3)
                    expected_col_name <- paste0(samp, "_q=", q_formatted)
                    samp_idx <- which(sample_names == expected_col_name)

                    if (length(samp_idx) > 0) {
                      all_ci_lower <- c(all_ci_lower, ci_lower_mat[, samp_idx[1]])
                      all_ci_upper <- c(all_ci_upper, ci_upper_mat[, samp_idx[1]])
                    }
                  }

                  if (length(all_ci_lower) > 0) {
                    ci_lower_final <- median(all_ci_lower, na.rm = TRUE)
                    ci_upper_final <- median(all_ci_upper, na.rm = TRUE)
                  } else {
                    ci_lower_final <- median(ci_lower_mat, na.rm = TRUE)
                    ci_upper_final <- median(ci_upper_mat, na.rm = TRUE)
                  }
                }

                if (!is.na(median_val)) {
                  plot_df <- rbind(plot_df, data.frame(q = matching_q_val, median = median_val,
                    ci_lower = ci_lower_final, ci_upper = ci_upper_final, group = group_val,
                    stringsAsFactors = FALSE))
                }
            }
        }
    }

    plot_df
}

# From diversity_core.R: Compute bootstrap CI for diversity measures
.bootstrap_diversity_ci <- function(bootstrap, result, genes, se_assay_mat, bootstrap_method = c("percentile", "bca"),
    bootstrap_ci, bootstrap_nboot, q, pseudocount, nthreads, bootstrap_include_diagnostics,
    verbose, effective_length = NULL, show_messages = FALSE, min_valid_frac = 0.75) {

    # Validate bootstrap_method parameter per Bioconductor code syntax standards
    bootstrap_method <- match.arg(bootstrap_method)

    bootstrap_ci_results <- NULL

    if (!bootstrap)
        return(NULL)

    if (verbose && show_messages)
        message("Computing bootstrap confidence intervals...")

    if (!is.numeric(bootstrap_ci) || bootstrap_ci <= 0 || bootstrap_ci >= 1) {
        stop("bootstrap_ci must be a probability in (0, 1)", call. = FALSE)
    }

    # Auto-suggest nboot if needed
    if (is.null(bootstrap_nboot)) {
        n_genes_filtered <- nrow(result) - 1
        if (n_genes_filtered < 1) {
            stop("After filtering, no genes remain. Try relaxing filter parameters.",
                call. = FALSE)
        }
        bootstrap_nboot <- .suggest_nboot(n_genes_filtered, use_bca = (bootstrap_method ==
            "bca"))
        if (verbose)
            message(sprintf("  -> Auto-suggested nboot = %d for %d genes", bootstrap_nboot,
                n_genes_filtered))
    }

    # Prepare data and compute bootstrap CIs For each (gene x sample) pair, we
    # compute one CI from bootstrap resampling of transcripts


    filtered_genes <- as.character(result[, 1])

    # Precompute transcript->gene index once (O(T)); each gene lookup is O(1)
    # instead of re-scanning the full genes vector (O(G*T)).
    gene_index <- split(seq_along(genes), genes)

    # Create a list where each element is bootstrap results for one (gene,
    # sample) pair
    bootstrap_results_list <- list()
    pair_metadata <- data.frame(gene = character(), sample_idx = integer())



    for (g_idx in seq_along(filtered_genes)) {
        g <- filtered_genes[g_idx]
        tx_mask <- gene_index[[as.character(g)]]

        if (is.null(tx_mask) || length(tx_mask) == 0)
            next

        # Get effective_length normalization for this gene's transcripts if
        # available
        el_for_gene_txs <- NULL
        if (!is.null(effective_length)) {
            # effective_length is indexed by transcript position
            el_for_gene_txs <- effective_length[tx_mask]
        }

        # For each sample, compute bootstrap CI on this gene's transcripts in
        # that sample
        for (s in seq_len(ncol(se_assay_mat))) {
            # Get transcript counts for this gene in this sample
            counts_vec <- se_assay_mat[tx_mask, s]



            # Compute bootstrap CI for this (gene, sample) pair Pass raw counts
            # AND effective_length separately so bootstrap handles both
            # correctly
            tryCatch({
                boot_result <- .calculate_tsallis_entropy_bootstrap(x = counts_vec,
                  q = q, norm = TRUE, nboot = bootstrap_nboot, ci = bootstrap_ci,
                  method = bootstrap_method, pseudocount = pseudocount, nthreads = nthreads,
                  verbose = FALSE, include_diagnostics = bootstrap_include_diagnostics,
                  effective_length = el_for_gene_txs, show_messages = show_messages,
                  min_valid_frac = min_valid_frac)

                bootstrap_results_list[[length(bootstrap_results_list) + 1]] <- boot_result
                pair_metadata <- rbind(pair_metadata, data.frame(gene = g, sample_idx = s))
            }, error = function(e) {
                if (verbose && show_messages)
                  warning("Bootstrap failed for ", g, " sample ", s, ": ", conditionMessage(e))
            })
        }
    }

    # Convert results list to named list for easier mapping
    if (nrow(pair_metadata) > 0) {
        result_names <- paste0(pair_metadata$gene, "_sample_", pair_metadata$sample_idx)
    } else {
        result_names <- character(0)
    }



    if (length(result_names) > 0) {
        names(bootstrap_results_list) <- result_names
    }

    if (verbose) {
        message("  Bootstrap data prepared:")
        message("    filtered_genes: ", length(filtered_genes))
        message("    Pairs analyzed: ", nrow(pair_metadata))
    }

    # Set the bootstrap_ci_results to our pre-computed list
    bootstrap_ci_results <- bootstrap_results_list

    list(bootstrap_ci_results = bootstrap_ci_results, bootstrap_nboot = bootstrap_nboot,
        bootstrap_method = bootstrap_method, bootstrap_ci = bootstrap_ci)
}

# From divergence_core.R: Configure parallel bootstrap execution
.bootstrap_configure_parallel <- function(bootstrap, nboot, method, num_genes, nthreads,
    progress) {
    # Validate bootstrap flag - use isTRUE to safely handle NA
    if (!isTRUE(bootstrap) && !isFALSE(bootstrap)) {
        bootstrap <- FALSE  # Default to no bootstrap if invalid
    }

    if (!isTRUE(bootstrap)) {
        nboot <- 0
    }

    # AUTO-SELECT NBOOT WHEN 'auto'
    if (isTRUE(bootstrap) && identical(nboot, "auto")) {
        use_bca <- !is.null(method) && identical(method, "bca")
        nboot <- .suggest_nboot(num_genes, use_bca = use_bca, nthreads = nthreads)
        if (isTRUE(progress)) {
            message("Auto-selected nboot = ", nboot, " for ", num_genes, " genes")
        }
    }

    # Configure parallel execution (fixes line 257 bug by using num_genes
    # parameter)
    parallel_config <- .configure_parallel(nthreads, num_genes)

    list(nboot = nboot, nthreads = parallel_config$nthreads, use_parallel = parallel_config$use_parallel)
}

# From divergence_core.R: Build bootstrap arguments for divergence
.bootstrap_build_args <- function(x, y, q_val, nboot, ci, method, log_base, pseudocount,
    gene_name, pair_ids = NULL) {
    args <- list(x = x, y = y, q = q_val, nboot = nboot, ci = ci, method = method,
        log_base = log_base, pseudocount = pseudocount, gene_name = gene_name, verbose = FALSE,
        paired = !is.null(pair_ids))

    if (!is.null(pair_ids)) {
        args$pair_ids <- pair_ids
    }

    args
}

# ============================================================================

#' @param object An object of class \code{tsenat_divergence_bootstrap_ci}.
#' @param ... Additional arguments (ignored).
#'
#' @return Invisibly returns the object after printing summary statistics.
#' Displays bootstrap distribution summary (mean, median, SD, min, max),
#' diagnostics (skewness, effective sample size), and stability metrics.
#'

#' @noRd
#' @method summary tsenat_divergence_bootstrap_ci

summary.tsenat_divergence_bootstrap_ci <- function(object, ...) {
    message("")
    message("=== Summary of Divergence Bootstrap ===")
    message("Bootstrap distribution:")
    message("  Mean:", round(mean(object$bootstrap_dist), 4))
    message("  Median:", round(stats::median(object$bootstrap_dist), 4))
    message("  SD:", round(stats::sd(object$bootstrap_dist), 4))
    message("  Min:", round(min(object$bootstrap_dist, na.rm = TRUE), 4))
    message("  Max:", round(max(object$bootstrap_dist, na.rm = TRUE), 4))

    # Diagnostics section
    message("")
    message("Diagnostics:")

    # Simple skewness calculation
    m <- mean(object$bootstrap_dist)
    s <- stats::sd(object$bootstrap_dist)
    if (s > 0) {
        n <- length(object$bootstrap_dist)
        skew <- (sum((object$bootstrap_dist - m)^3)/n)/s^3
        message("  Skewness:", round(skew, 4))
    } else {
        message("  Skewness: N/A (no variation)")
    }

    # Effective sample size (ESS) - simplified as ratio of bootstrap replicates
    # with unique values
    n_unique <- length(unique(round(object$bootstrap_dist, 6)))
    n_total <- length(object$bootstrap_dist)
    ess <- (n_unique/n_total) * 100
    message("  Effective sample size:", round(ess, 1), "%")

    message("")
    message("Stability metrics:")
    message(sprintf("  CI width to estimate ratio: %.2f", (object$upper_ci - object$lower_ci)/pmax(object$estimate,
        0.01)))

    # Check for multimodality (simple approximation)
    modes <- length(unique(round(object$bootstrap_dist, 3)))
    message(sprintf("  Unique rounded values: %d", modes))

    invisible(object)
}

