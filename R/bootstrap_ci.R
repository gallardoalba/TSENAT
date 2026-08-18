.calculate_tsallis_entropy_bootstrap <- function(x = NULL, se = NULL, res = NULL,
    top_n = 1, q = 2, norm = TRUE, nboot = "auto", ci = 0.95, method = c("percentile",
        "bca"), log_base = exp(1), pseudocount = 0, what = c("S", "D"), gene_name = NULL,
    verbose = TRUE, include_diagnostics = TRUE, use_job = FALSE, nthreads = 1, paired = FALSE,
    effective_length = NULL, show_messages = FALSE, min_valid_frac = 0.75,
    resample_by = c("read", "replicate"), counts_matrix = NULL) {

    method <- match.arg(method)
    what <- match.arg(what)
    resample_by <- match.arg(resample_by)

    nboot <- ._bootstrap_resolve_nboot(nboot, x, se, method, nthreads)

    # Dispatch by input type: matrix → bulk, SE+res → single-SE, vector → single
    result <- ._bootstrap_dispatch_input(x, se, res, top_n, q, norm, nboot, ci,
        method, log_base, pseudocount, what, gene_name, verbose, include_diagnostics,
        use_job, nthreads, paired, resample_by, counts_matrix)
    if (!is.null(result)) return(invisible(result))

    # If SE/result path returned NULL (e.g., no genes with sufficient counts),
    # the warning was already issued; do not fall through to single-vector path
    if (is.null(x)) return(invisible(NULL))

    # Single vector path: validate, compute, print
    ._bootstrap_single_vector(x, q, nboot, ci, paired, show_messages,
        effective_length, pseudocount, norm, method, log_base, what,
        min_valid_frac, resample_by, counts_matrix,
        include_diagnostics, use_job, gene_name, verbose)
}

#' Resolve nboot: auto-select or use provided value
#' @noRd
._bootstrap_resolve_nboot <- function(nboot, x, se, method, nthreads) {
    if (!identical(nboot, "auto")) return(nboot)
    n_genes <- if (!is.null(x) && is.matrix(x)) nrow(x)
               else if (!is.null(se)) nrow(se) else 1
    .bootstrap_auto_select_nboot(n_genes, method == "bca", nthreads)
}

#' Dispatch: matrix (bulk), SE+res, or return NULL for single-vector path
#' @noRd
._bootstrap_dispatch_input <- function(x, se, res, top_n, q, norm, nboot, ci,
    method, log_base, pseudocount, what, gene_name, verbose, include_diagnostics,
    use_job, nthreads, paired, resample_by, counts_matrix) {
    if (!is.null(x) && is.matrix(x)) {
        return(.bootstrap_process_matrix(x, q, norm, nboot, ci, method,
            log_base, pseudocount, what, gene_name, verbose, include_diagnostics,
            use_job, nthreads, paired, resample_by = resample_by,
            counts_matrix = counts_matrix))
    }
    if (!is.null(se) && !is.null(res)) {
        return(.bootstrap_process_se(se, res, top_n, q, norm, nboot, ci, method,
            log_base, pseudocount, what, gene_name, verbose, include_diagnostics,
            use_job, paired))
    }
    if (is.null(x)) stop("Either 'x' or both 'se' and 'res' must be provided")
    NULL
}

#' Single-vector bootstrap pipeline: validate → compute → print
#' @noRd
._bootstrap_single_vector <- function(x, q, nboot, ci, paired, show_messages,
    effective_length, pseudocount, norm, method, log_base, what,
    min_valid_frac, resample_by, counts_matrix,
    include_diagnostics, use_job, gene_name, verbose) {
    .bootstrap_validate_inputs(x, q, nboot, ci, paired, show_messages)
    .validate_bootstrap_data(x, effective_length = effective_length, pseudocount = pseudocount)

    if (length(q) > 1) {
        return(._bootstrap_multi_q(x, q, norm, nboot, ci, method, log_base,
            pseudocount, what, gene_name, verbose, include_diagnostics, use_job,
            paired, effective_length, min_valid_frac, resample_by, counts_matrix))
    }

    ci_data <- .bootstrap_compute_ci(x, q, norm, nboot, ci, method, log_base, pseudocount,
        what, paired, effective_length, min_valid_frac, resample_by = resample_by,
        counts_matrix = counts_matrix)

    diag_list <- .bootstrap_compute_diag(ci_data$point_est, ci_data$bootstrap_dist,
        use_job, paired, x, q, norm, nboot, ci, method, log_base, pseudocount, what,
        ci_data$accel_factor)

    result <- .bootstrap_assemble_result(ci_data$point_est, ci_data$ci_result,
        ci_data$bootstrap_dist, ci, method, nboot, diag_list, include_diagnostics, use_job)

    .bootstrap_print_results(result, gene_name, ci, verbose)
    invisible(result)
}

#' Multi-q bootstrap pipeline
#' @noRd
._bootstrap_multi_q <- function(x, q, norm, nboot, ci, method, log_base,
    pseudocount, what, gene_name, verbose, include_diagnostics, use_job,
    paired, effective_length, min_valid_frac, resample_by, counts_matrix) {
    result <- .bootstrap_process_multiple_q(x, q, norm, nboot, ci, method, log_base,
        pseudocount, what, gene_name, verbose, include_diagnostics, use_job,
        paired, effective_length, min_valid_frac, resample_by = resample_by,
        counts_matrix = counts_matrix)
    if (verbose && !is.null(gene_name)) {
        message("Bootstrap Confidence Intervals for ", gene_name, " (multiple q values)")
        for (i in seq_along(result)) {
            res <- result[[i]]
            message("q=", q[i], ": [", sprintf("%.6f", res$lower_ci), ", ",
                sprintf("%.6f", res$upper_ci), "]")
        }
    }
    invisible(result)
}

#' Summary and Printing for Bootstrap CI Results
#'
#' @param object An object of class \code{tsenat_bootstrap_ci} from
#'   \code{calculate_tsallis_entropy_bootstrap}.
#' @param x An object of class \code{tsenat_bootstrap_ci}.
#' @param \ldots Additional arguments (unused).
#'
#' @return Invisibly returns the object.
#'

#' @noRd
#' @method summary tsenat_bootstrap_ci

summary.tsenat_bootstrap_ci <- function(object, ...) {
    message("=== Tsallis Entropy Bootstrap Confidence Interval ===")
    message("Method: ", object$method)
    message("Bootstrap replicates: ", object$nboot)
    message("Confidence level: ", object$ci_level * 100, "%")
    message("Point estimate (S_q): ", sprintf("%.6f", object$estimate))
    message("Lower CI: ", sprintf("%.6f", object$lower_ci))
    message("Upper CI: ", sprintf("%.6f", object$upper_ci))
    message("CI width: ", sprintf("%.6f", object$upper_ci - object$lower_ci))
    message("Bootstrap distribution summary:")
    stats <- summary(object$bootstrap_dist)
    message(paste(capture.output(str(stats)), collapse = "\n"))

    # Display diagnostics if available (Zhang & Cao 2023; Friedl & Stampfer 2002)
    if (!is.null(object$diagnostics)) {
        message("")
        message("=== CI Quality Diagnostics (jackknife resampling, e.g., Zhang & Cao 2023) ===")
        message("Effective sample size: ", sprintf("%.1f", object$diagnostics$effective_sample_size),
            " (>= n * 0.5 is good)")
        message("Skewness: ", sprintf("%.4f", object$diagnostics$skewness), " (|.| > 2 suggests unreliability)")
        message("Bias: ", sprintf("%.6f", object$diagnostics$bias), " (distance from median to estimate)")
        if (!is.na(object$diagnostics$acceleration_factor)) {
            message("Acceleration (BCa): ", sprintf("%.6f", object$diagnostics$acceleration_factor),
                " (skewness correction factor)")
        }
        message("")
        message("Interpretation: Check effective_sample_size and skewness to assess CI reliability.")
    }
    invisible(object)
}


#' @noRd
#' @exportS3Method
print.tsenat_bootstrap_ci <- function(x, ...) {
    ci_pct <- if (!is.null(x$ci_level)) round(x$ci_level * 100, 1) else 95
    message("Tsallis Entropy Bootstrap Confidence Interval")
    message("Point estimate: ", sprintf("%.6f", x$estimate))
    message(ci_pct, "% CI: [", sprintf("%.6f", x$lower_ci), ", ", sprintf("%.6f", x$upper_ci),
        "]")
    invisible(x)
}

#' @noRd
#' @exportS3Method
print.tsenat_bootstrap_ci_list <- function(x, ...) {
    message("Bootstrap Confidence Intervals for Multiple q Values")
    message("Number of q values: ", length(x))
    for (i in seq_along(x)) {
        ci_pct <- if (!is.null(x[[i]]$ci_level)) round(x[[i]]$ci_level * 100, 1) else 95
        message("\n  q = ", names(x)[i], ":")
        message("    Estimate: ", sprintf("%.6f", x[[i]]$estimate))
        message("    ", ci_pct, "% CI: [", sprintf("%.6f", x[[i]]$lower_ci), ", ", sprintf("%.6f",
            x[[i]]$upper_ci), "]")
    }
    invisible(x)
}

#' Jackknife-of-Bootstrap (JOB) CI Stability Assessment
#'
#' Compute bootstrap CIs leaving out each observation and assess stability.
#' JOB is a hybrid approach combining jackknife (leave-one-out) validation with
#' bootstrap confidence intervals, providing more robust estimates when data
#' is limited (Zhang & Cao 2023).
#'
#' @param x Numeric vector: original transcript counts
#' @param q Numeric: Tsallis entropy parameter
#' @param norm Logical: normalize entropy calculation
#' @param nboot Integer: bootstrap replicates per jackknife sample
#' @param ci Numeric: confidence level
#' @param method Character: 'percentile' or 'bca'
#' @param log_base Numeric: logarithm base
#' @param pseudocount Numeric: added to proportions
#' @param what Character: 'S' (entropy) or 'D' (Hill numbers)
#'
#' @return List with:
#'   \describe{
#'     \item{ci_lower_stable}{Lower CI bound (stability-adjusted)}
#'     \item{ci_upper_stable}{Upper CI bound (stability-adjusted)}
#'     \item{ci_width_variation}{Coefficient of variation of CI widths across jackknife samples}
#'     \item{bound_variability}{Max relative change in bounds across jackknife samples}
#'     \item{n_outlier_bounds}{Count of jackknife samples with 
#' outlier CI bounds}
#'   }
#'
#' @details
#' JOB Procedure (Zhang & Cao 2023):
#' 1. Compute bootstrap CI on full dataset
#' 2. For each observation i, remove it and compute bootstrap CI on
#' remaining data
#' 3. Track CI stability: are bounds consistent across leave-one-out replicates?
#' 4. Return stability metrics assessing robustness
#'
#' Conservative estimate: use maximum of lower bounds and minimum of upper
#' bounds
#' across all jackknife replicates to get widest CI (most conservative).
#'

#' @noRd
.compute_job <- function(x, q, norm, nboot, ci, method, log_base, pseudocount, what,
    paired = FALSE) {
    n <- length(x)
    if (n < 3) {
        warning("JOB requires n >= 3. Skipping JOB computation.")
        return(NULL)
    }

    # Store bootstrap CI from full dataset and each jackknife sample
    job_cis <- list()

    # Full dataset CI (index 0)
    full_ci <- .ci_from_bootstrap(x, q = q, norm = norm, nboot = nboot, ci = ci,
        method = method, log_base = log_base, pseudocount = pseudocount, what = what,
        paired = paired)
    job_cis[[1]] <- data.frame(lower = full_ci$lower, upper = full_ci$upper, width = full_ci$upper -
        full_ci$lower, label = "full")

    # Leave-one-out jackknife replicates
    for (i in seq_len(n)) {
        x_minus_i <- x[-i]

        loo_ci <- .ci_from_bootstrap(x_minus_i, q = q, norm = norm, nboot = nboot,
            ci = ci, method = method, log_base = log_base, pseudocount = pseudocount,
            what = what, paired = paired)
        job_cis[[i + 1]] <- data.frame(lower = loo_ci$lower, upper = loo_ci$upper,
            width = loo_ci$upper - loo_ci$lower, label = paste0("LOO_", i))
    }

    # Convert to data frame for easier analysis
    job_df <- do.call(rbind, job_cis)

    # Compute stability metrics Conservative bounds: widest interval from all
    # jackknife samples
    ci_lower_stable <- min(job_df$lower, na.rm = TRUE)
    ci_upper_stable <- max(job_df$upper, na.rm = TRUE)

    # Stability metrics
    widths <- job_df$width[-1]  # Exclude full dataset width
    ci_width_variation <- sd(widths, na.rm = TRUE)/mean(widths, na.rm = TRUE)

    # BUG FIX (March 2026): Numerical stability in bound variation calculation
    # Previous code divided by abs(value) + 1e-10 which is unstable for small
    # values New approach: Use ratio of range to mean absolute value (robust to
    # scale)
    lower_range <- max(job_df$lower[-1], na.rm = TRUE) - min(job_df$lower[-1], na.rm = TRUE)
    lower_mean_abs <- mean(abs(job_df$lower[-1]), na.rm = TRUE)
    lower_variation <- if (lower_mean_abs > 1e-08)
        lower_range/lower_mean_abs else 0

    upper_range <- max(job_df$upper[-1], na.rm = TRUE) - min(job_df$upper[-1], na.rm = TRUE)
    upper_mean_abs <- mean(abs(job_df$upper[-1]), na.rm = TRUE)
    upper_variation <- if (upper_mean_abs > 1e-08)
        upper_range/upper_mean_abs else 0

    bound_variability <- max(lower_variation, upper_variation, na.rm = TRUE)

    # Count outlier CI bounds (> 2 SD from jackknife mean)
    lower_mean <- mean(job_df$lower[-1], na.rm = TRUE)
    lower_sd <- sd(job_df$lower[-1], na.rm = TRUE)
    upper_mean <- mean(job_df$upper[-1], na.rm = TRUE)
    upper_sd <- sd(job_df$upper[-1], na.rm = TRUE)

    lower_outliers <- sum(abs(job_df$lower[-1] - lower_mean) > 2 * lower_sd, na.rm = TRUE)
    upper_outliers <- sum(abs(job_df$upper[-1] - upper_mean) > 2 * upper_sd, na.rm = TRUE)
    n_outlier_bounds <- lower_outliers + upper_outliers

    return(list(ci_lower_stable = ci_lower_stable, ci_upper_stable = ci_upper_stable,
        ci_width_variation = ci_width_variation, bound_variability = bound_variability,
        n_outlier_bounds = as.numeric(n_outlier_bounds), jackknife_cis = job_df))
}

#' Helper: Compute bootstrap CI from data (internal utility)
#'

#' @noRd
.ci_from_bootstrap <- function(x, q, norm, nboot, ci, method, log_base, pseudocount,
    what, paired = FALSE) {
    bootstrap_dist <- .bootstrap_resample_optimized(x, q = q, norm = norm, nboot = nboot,
        log_base = log_base, pseudocount = pseudocount, what = what, paired = paired)

    if (method == "percentile") {
        ci_result <- .ci_percentile(bootstrap_dist, ci = ci)
    } else {
        ci_result <- .ci_bca(x, bootstrap_dist, q = q, norm = norm, ci = ci, log_base = log_base,
            pseudocount = pseudocount, what = what)
    }

    return(ci_result)
}

#' Bootstrap Confidence Intervals for Q-curve Data
#'
#' Helper function to compute bootstrap confidence intervals for Tsallis entropy
#' across multiple q-values and groups. Resamples genes (not individual
#' transcripts)
#' with replacement and computes quantile-based confidence intervals for
#' medians.
#'
#' @param long Data frame in long format with columns: Gene, q, tsallis, group.
#'   Typically output from \code{.prepare_tsallis_long()}.
#' @param unique_q Numeric vector of unique q-values (sorted).
#' @param groups Character vector of group names (e.g., c('group1', 'group2')).
#' @param ci_level Numeric; confidence level (default: 0.95 for 95% CI).
#' @param n_bootstrap Integer; number of bootstrap replicates (default: 500).
#'
#' @return A nested list structure:
#'  \code{[[group]][[q_string]]} where each element
#'   contains a list with:
#'   \describe{
#'   \item{median}{Median entropy from original data.}
#'   \item{ci_lower}{Lower confidence bound.}
#'   \item{ci_upper}{Upper confidence bound.}
#'   \item{n}{Number of genes in group at that q-value.}
#'   }
#'
#' @details
#' **Bootstrap Procedure:**
#' For each group and q-value combination:
#' 1. Extract all gene-level entropy values
#' 2. Resample genes with replacement n_bootstrap times
#' 3. For each resample, compute the median entropy
#' 4. Calculate CI bounds from percentiles of bootstrap distribution
#'
#' Uses percentile method with quantile type 7 (recommended by Hyndman &
#' Fan, 1996).
#'

#' @noRd
#' @examples
#' # After .prepare_tsallis_long()
#' set.seed(42)
#' # Create sample long-format diversity data
#' long_data <- data.frame(
#'   entropy = runif(100, 0, 5),
#'   q = rep(c(0.5, 1.0, 1.5, 2.0, 2.5), 20),
#'   group = rep(c('control', 'treatment'), 50)
#' )
#' unique_q <- c(0.5, 1.0, 1.5, 2.0, 2.5)
#' ci_results <- .compute_bootstrap_qcurve_cis(
#'   long = long_data, unique_q = unique_q,
#'   groups = c('control', 'treatment'), ci_level = 0.95, n_bootstrap = 100
#' )

.compute_bootstrap_qcurve_cis <- function(long, unique_q, groups, ci_level = 0.95,
    n_bootstrap = 500) {

    bootstrap_results <- list()

    # For each group and q-value, compute bootstrap CI
    for (g in groups) {
        group_data <- long %>%
            dplyr::filter(group == g) %>%
            dplyr::select(Gene, q, tsallis)

        bootstrap_results[[g]] <- list()

        for (q_val in unique_q) {
            q_data <- group_data %>%
                dplyr::filter(q == q_val) %>%
                dplyr::pull(tsallis)

            if (length(q_data) < 2) {
                warning("Group '", g, "' has < 2 samples at q=", q_val)
                bootstrap_results[[g]][[as.character(q_val)]] <- list(median = NA_real_,
                  ci_lower = NA_real_, ci_upper = NA_real_, n = length(q_data))
                next
            }

            # Compute bootstrap replicates (resample genes with replacement)
            boot_replicates <- numeric(n_bootstrap)
            for (b in seq_len(n_bootstrap)) {
                boot_idx <- sample(seq_along(q_data), replace = TRUE)
                boot_replicates[b] <- median(q_data[boot_idx], na.rm = TRUE)
            }

            # Compute confidence interval from percentiles
            ci_lower <- as.numeric(quantile(boot_replicates, (1 - ci_level)/2, type = 7,
                na.rm = TRUE))
            ci_upper <- as.numeric(quantile(boot_replicates, 1 - (1 - ci_level)/2,
                type = 7, na.rm = TRUE))

            bootstrap_results[[g]][[as.character(q_val)]] <- list(median = median(q_data,
                na.rm = TRUE), ci_lower = ci_lower, ci_upper = ci_upper, n = length(q_data))
        }
    }

    return(bootstrap_results)
}

#' Suggest Adaptive Bootstrap Sample Size
#'
#' Recommends an appropriate number of bootstrap replicates based on the
#' number of genes
#' being analyzed and the method (percentile vs BCa). This helps balance
#' computational
#' efficiency with statistical accuracy.
#'
#' @param n_genes Integer: Number of genes to be analyzed simultaneously.
#'                If analyzing a single gene, use n_genes=1. For multiple genes,
#'                provide the total count.
#' @param use_bca Logical: If TRUE (default FALSE), recommends higher sample
#' sizes
#'                suitable for the more computationally intensive BCa method.
#'                If FALSE, recommends for the faster percentile method.
#'
#' @return Integer: Recommended number of bootstrap replicates.
#'
#' @details
#' **Rationale (from paper Springer Handbook (2006) - Bootstrap computational methods):**
#'
#' The BCa (bias-corrected and accelerated) method is more accurate but requires
#' higher computational cost due to jackknife calculations. For datasets
#' with many genes,
#' the percentile method offers a good accuracy-to-speed trade-off.
#'
#' **Recommendations by scenario:**
#' - Single gene analysis: 1000-2000 replicates (detailed inference)
#' - Small gene sets (2-5 genes): 500-1000 replicates (balanced)
#' - Large gene sets (>10 genes): 250-500 replicates (speed-prioritized)
#'
#' **Usage:**
#' ```
#' # For analyzing 3 genes with percentile method
#' nboot <- .suggest_nboot(n_genes = 3, use_bca = FALSE)
#' # Returns 500
#'
#' # For single gene with BCa method (more precise inference)
#' nboot <- .suggest_nboot(n_genes = 1, use_bca = TRUE)
#' # Returns 2000
#' ```
#'
#' @references
#' Paper Springer Handbook (2006): Bootstrap computational methods and efficiency trade-offs.
#' Discusses how sample size affects accuracy and speed of bootstrap inference.
#'
#' @examples
#' .suggest_nboot(1, use_bca = FALSE)   # Single gene, percentile: 1000
#' .suggest_nboot(1, use_bca = TRUE)    # Single gene, BCa: 1500
#' .suggest_nboot(3, use_bca = FALSE)   # 3 genes, percentile: 500
#' .suggest_nboot(15, use_bca = FALSE)  # 15 genes, percentile: 250
#'
#' @noRd

.suggest_nboot <- function(n_genes, use_bca = FALSE, nthreads = 1) {

    # Input validation
    if (!is.numeric(n_genes) || n_genes < 1 || n_genes != as.integer(n_genes)) {
        stop("'n_genes' must be a positive integer")
    }
    if (!is.logical(use_bca)) {
        stop("'use_bca' must be logical")
    }
    if (nthreads < 1 || nthreads != as.integer(nthreads)) {
        stop("'nthreads' must be a positive integer")
    }

    # Base recommendations with smooth scaling (avoids discontinuous jumps)
    # Recommendations follow Springer Handbook (2006) (Bootstrap computational
    # methods) efficiency guidelines
    base_nboot <- if (n_genes == 1) {
        # Single gene: detailed inference justified
        1000
    } else if (n_genes <= 5) {
        # Small gene set: balance accuracy and speed
        500
    } else if (n_genes <= 20) {
        # Medium gene set: smooth interpolation (500 → 250 as genes 5 → 20)
        round(500 - (n_genes - 5) * 16.67)
    } else {
        # Large gene set: prioritize speed
        250
    }

    # Adjust for BCa (bias-corrected and accelerated method) BCa requires
    # jackknife calculations, approximately 50% more replicates needed
    if (use_bca) {
        base_nboot <- round(base_nboot * 1.5)
    }

    # Account for parallelization (more threads = less need for huge sample
    # sizes) Diminishing returns after ~4 threads, floor at 0.75x
    parallel_factor <- max(0.75, 1 - log(nthreads)/12)
    base_nboot <- round(base_nboot * parallel_factor)

    # Enforce minimum (need >= 100 for meaningful percentile CIs)
    max(100, base_nboot)
}

# NOTE: .compute_skewness defined in calc_sait_helpers.R Bootstrap distributions
# are clean (generated from rmultinom + entropy calculations), so the version
# with na.rm parameter is safe to use with default na.rm=TRUE

#' Compute Effective Sample Size from Bootstrap Data
#'
#' Estimate effective sample size using the relationship between bootstrap
#' replicates autocorrelation and true sample size. Higher values indicate
#' more independent bootstrap samples (better CI reliability).
#'

#' @noRd
.compute_effective_n <- function(x) {
    n <- length(x)
    if (n < 2)
        return(n)

    # Compute lag-1 autocorrelation
    x_centered <- x - mean(x, na.rm = TRUE)
    acf_1 <- sum(x_centered[-n] * x_centered[-1], na.rm = TRUE)/sum(x_centered^2,
        na.rm = TRUE)
    acf_1 <- max(-0.999, min(0.999, acf_1))  # Bound to (-1, 1)

    # Effective sample size accounting for autocorrelation magnitude Uses
    # absolute value following GEE standard (Liang & Zeger, 2001; Pan 2001) and
    # modern variance estimation methodology (Lai et al. 2011; Demidenko 2013;
    # Bates et al. 2015; Li et al. 2021; Hardin & Hilbe 2013).
    # Correlation magnitude (not sign) affects variance structure
    # symmetrically.
    n_eff <- n/(1 + 2 * abs(acf_1))

    return(min(n, max(1, n_eff)))  # Ensure 1 <= n_eff <= n
}




