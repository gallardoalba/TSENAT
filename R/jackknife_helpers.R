################################################################################
#
#' Internal: Jackknife leave-one-out resampling
#'
#' Consolidated jackknife procedure used by multiple functions to eliminate duplication.
#' Performs leave-one-out resampling and computes influence values.
#'
#' @param counts Matrix/data.frame where rows are genes/observations, columns are samples
#' @param entropy_fn Function to compute entropy (default: .entropy_core)
#' @param q Numeric. Tsallis parameter. Default: 1.0
#' @param threshold Numeric. Percentile cutoff for outlier detection (0-100). Default: 90
#' @param norm Logical. Normalize entropy. Default: TRUE
#' @param log_base Numeric. Logarithm base. Default: exp(1)
#' @param pseudocount Numeric. Add to counts before normalization. Default: 0
#'
#' @return List with elements:
#'   - estimate: original entropy value
#'   - jackknife_estimates: leave-one-out entropy values
#'   - influence: influence of each observation
#'   - jackknife_se: standard error from jackknife procedure
#'   - outlier_indices: which observations are outliers
#'   - outlier_cutoff_value: threshold value used for outlier detection
#'   - outlier_threshold: percentile used
#'

#' @noRd
.jackknife_resampling <- function(counts, entropy_fn = .entropy_core, 
                                   q = 1, threshold = 90, norm = TRUE, 
                                   log_base = exp(1), pseudocount = 0) {
  # Ensure matrix format
  counts <- as.matrix(counts)
  
  # Input validation
  if (nrow(counts) == 0 || ncol(counts) == 0) {
    return(NULL)
  }
  
  if (nrow(counts) < 2) {
    warning("Insufficient observations for jackknife (need >= 2)")
    return(NULL)
  }
  
  # Compute original entropy estimate
  total_count <- sum(counts, na.rm = TRUE) + ncol(counts) * pseudocount
  p_full <- (colSums(counts, na.rm = TRUE) + pseudocount) / total_count
  
  estimate <- entropy_fn(p_full, q = q, norm = norm, log_base = log_base)
  
  if (is.na(estimate) || is.nan(estimate)) {
    return(NULL)
  }
  
  # Leave-one-out jackknife loop
  n_obs <- nrow(counts)
  jackknife_estimates <- numeric(n_obs)
  
  for (i in seq_len(n_obs)) {
    # Remove observation i
    counts_minus_i <- counts[-i, , drop = FALSE]
    total_minus_i <- sum(counts_minus_i, na.rm = TRUE) + ncol(counts_minus_i) * pseudocount
    
    if (total_minus_i <= 0) {
      jackknife_estimates[i] <- NA_real_
      next
    }
    
    p_minus_i <- (colSums(counts_minus_i, na.rm = TRUE) + pseudocount) / total_minus_i
    jackknife_estimates[i] <- entropy_fn(p_minus_i, q = q, norm = norm, log_base = log_base)
  }
  
  # Compute influence (absolute change when removing each observation)
  influence <- abs(jackknife_estimates - estimate)
  
  # Compute jackknife standard error (bias-corrected)
  theta_jack_mean <- mean(jackknife_estimates, na.rm = TRUE)
  jackknife_se <- sqrt(((n_obs - 1) / n_obs) * sum((jackknife_estimates - theta_jack_mean)^2, na.rm = TRUE))
  
  # Outlier detection: observations with influence > threshold percentile
  outlier_cutoff <- stats::quantile(influence, threshold / 100, na.rm = TRUE)
  outlier_mask <- influence > outlier_cutoff & !is.na(influence)
  outlier_indices <- which(outlier_mask)
  
  return(list(
    estimate = estimate,
    jackknife_estimates = jackknife_estimates,
    influence = influence,
    jackknife_se = jackknife_se,
    outlier_indices = outlier_indices,
    outlier_cutoff_value = as.numeric(outlier_cutoff),
    outlier_threshold = threshold,
    n_observations = n_obs,
    q = q,
    norm = norm
  ))
}

#' Internal: Batch jackknife for multiple observations (genes)
#'
#' Apply jackknife resampling to multiple genes/observations efficiently
#'
#' @param counts_matrix Matrix where rows = samples, columns = genes/species
#' @param entropy_fn Function to compute entropy. Default: .entropy_core
#' @param q Numeric. Tsallis parameter. Default: 1.0
#' @param threshold Numeric. Outlier detection percentile. Default: 90
#' @param norm Logical. Normalize entropy. Default: TRUE
#' @param log_base Numeric. Log base. Default: exp(1)
#' @param pseudocount Numeric. Pseudocount. Default: 0
#' @param verbose Logical. Print progress. Default: FALSE
#'
#' @return List of jackknife results (one per gene), with class "tsenat_jackknife_list"
#'

#' @noRd
.jackknife_batch <- function(counts_matrix, entropy_fn = .entropy_core,
                              q = 1, threshold = 90, norm = TRUE,
                              log_base = exp(1), pseudocount = 0, verbose = FALSE) {
  counts_matrix <- as.matrix(counts_matrix)
  
  if (ncol(counts_matrix) == 0) {
    return(list())
  }
  
  # Apply jackknife to each gene (column)
  results <- vector("list", ncol(counts_matrix))
  names(results) <- colnames(counts_matrix)
  
  for (gene_idx in seq_len(ncol(counts_matrix))) {
    gene_counts <- counts_matrix[, gene_idx, drop = FALSE]
    
    result <- .jackknife_resampling(
      gene_counts,
      entropy_fn = entropy_fn,
      q = q,
      threshold = threshold,
      norm = norm,
      log_base = log_base,
      pseudocount = pseudocount
    )
    
    # Handle cases where jackknife fails
    if (is.null(result)) {
      result <- list(
        estimate = NA_real_,
        jackknife_estimates = NA_real_,
        influence = NA_real_,
        jackknife_se = NA_real_,
        outlier_indices = integer(0),
        n_observations = nrow(gene_counts)
      )
    }
    
    class(result) <- c("tsenat_jackknife", "list")
    results[[gene_idx]] <- result
    
    if (verbose && gene_idx %% 10 == 0) {
      message(sprintf("Jackknife processed %d/%d genes", gene_idx, ncol(counts_matrix)))
    }
  }
  
  class(results) <- c("tsenat_jackknife_list", "list")
  return(results)
}
