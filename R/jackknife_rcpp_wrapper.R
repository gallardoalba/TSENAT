# Rcpp wrapper and fallback infrastructure for jackknife resampling Purpose:
# Accelerate jackknife computation with optional C++ implementation Fallback:
# Pure R version available if Rcpp unavailable

# Check if Rcpp implementation is available Since the package compiles
# successfully, Rcpp is always available
.rcpp_available <- function() {
    TRUE
}

# Initialize Rcpp check on first use (cached via options)
.initialize_rcpp_check <- function() {
    # Check if already cached
    cached <- getOption("tsenat.rcpp_initialized", default = FALSE)

    if (!cached) {
        # Perform the check
        rcpp_available <- .rcpp_available()

        # Cache the result using options (not subject to namespace locking)
        options(tsenat.rcpp_initialized = TRUE, tsenat.rcpp_is_available = rcpp_available)

        if (rcpp_available) {
            packageStartupMessage("Rcpp acceleration enabled for jackknife resampling")
        }
    }

    # Return the cached result
    return(getOption("tsenat.rcpp_is_available", default = FALSE))
}


#' Jackknife Resampling Computation
#'
#' @description
#' C++ wrapper for leave-one-out jackknife resampling to estimate
#' influence and confidence intervals for diversity metrics.
#'
#' @param counts \code{matrix}. Count matrix (rows=features, columns=samples).
#' @param q \code{numeric}. Tsallis q parameter. Default: 1.0.
#' @param normalize \code{logical}. Normalize entropy? Default: TRUE.
#' @param log_base \code{numeric}. Logarithm base. Default: e (natural log).
#' @param pseudocount \code{numeric}.  Pseudocount for  abundance inflation.
#'  Default:  0. 0.
#'
#' @return \code{list}.  Results including jackknife estimates,
#'  confidence intervals,
#'   and influence diagnostics.
#'
#' @details
#' Performs leave-one-out jackknife resampling with C++ acceleration for
#' efficiency.
#' Estimates confidence intervals and diagnostic influence scores.
#'
#' @keywords internal
#' @noRd
jackknife_resampling_cpp <- function(counts, q = 1, normalize = TRUE, log_base = exp(1),
    pseudocount = 0) {
    # Internal wrapper - calls the C++ function directly via .Call() Not
    # exported to user (internal use only)
    .Call("_TSENAT_jackknife_resampling_cpp", as.matrix(counts), as.numeric(q), as.logical(normalize),
        as.numeric(log_base), as.numeric(pseudocount), PACKAGE = "TSENAT")
}

# Hybrid wrapper: Use C++ implementation (no fallback)
.jackknife_resampling_hybrid <- function(counts, entropy_fn = .entropy_core, q = 1,
    threshold = 90, norm = TRUE, log_base = exp(1), pseudocount = 0) {

    # Initialize Rcpp check once
    .initialize_rcpp_check()

    # Call internal C++ function via wrapper
    result_cpp <- jackknife_resampling_cpp(counts, q = q, normalize = norm, log_base = log_base,
        pseudocount = pseudocount)

    # Check for errors from C++
    if (!is.null(result_cpp$error)) {
        stop(result_cpp$error)
    }

    # Add missing fields to match original function output
    outlier_influence <- result_cpp$influence
    outlier_cutoff <- stats::quantile(outlier_influence, threshold/100, na.rm = TRUE)
    outlier_mask <- outlier_influence > outlier_cutoff & !is.na(outlier_influence)
    outlier_indices <- which(outlier_mask)

    result <- list(estimate = result_cpp$estimate, jackknife_estimates = result_cpp$jackknife_estimates,
        influence = result_cpp$influence, jackknife_se = result_cpp$jackknife_se,
        outlier_indices = outlier_indices, outlier_cutoff_value = unname(outlier_cutoff),
        outlier_threshold = threshold, n_observations = nrow(counts))

    class(result) <- c("tsenat_jackknife", "list")
    return(result)
}

# Set the hybridversion as the primary implementation Note: R implementation
# has been removed and moved to C++
.jackknife_resampling <- .jackknife_resampling_hybrid

# Provide user-facing function to check Rcpp status
viz_rcpp_status <- function() {
    list(rcpp_available = .initialize_rcpp_check(), description = "Rcpp acceleration for jackknife resampling computation")
}
