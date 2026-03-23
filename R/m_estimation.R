#' Internal Helper Functions for M-Estimation
#' @keywords internal
#' @noRd
NULL

#' Robust Statistical Methods for Differential Analysis
#'
#' Alternative statistical approaches that are resistant to outliers and
#' extreme values, particularly useful for RNA-seq data with unusual
#' distributions or unexpected outliers.
#'
#' @details
#' **Trimmed Wilcoxon Test:**
#' Combines the robustness of trimming extreme values with the power of
#' non-parametric rank-based testing. Extreme values are excluded from
#' the analysis before computing ranks, making the test resistant to
#' influential outliers.
#'
#' **M-Estimation:**
#' Uses iteratively re-weighted least squares (IRLS) to estimate location
#' differences while down-weighting outliers. The Huber loss function
#' provides a compromise between least squares (sensitive to outliers)
#' and absolute deviations (less efficient).
#'
#' **Influence Diagnostics (DFBETA Standardization):**
#' Identifies samples that have disproportionate influence on the fitted model.
#' Uses leave-one-out (LOO) analysis with standardized DFBETA statistics:
#'   DFBETA_i = (coef_full - coef_{-i}) / SE(coef_full)
#' 
#' This approach standardizes influence by the precision of the estimate,
#' allowing meaningful comparison across genes with different levels of 
#' variability. Samples with |DFBETA| > 2/sqrt(n) are flagged as problematic.
#' 
#' This is analogous to classical regression diagnostics (Cook's distance,
#' DFBETA) but applied in the robust M-estimation context.
#'
#' @references
#' Wilkinson, L. (2005). The grammar of graphics. Springer.
#' Huber, P. J. (1981). Robust Statistics. John Wiley & Sons.
#' Maronna, R. A., Martin, R. D., & Yohai, V. J. (2006).
#' Robust Statistics: Theory and Methods. John Wiley & Sons.
#' Cook, R. D., & Weisberg, S. (1982). Residuals and influence in regression.
#' Chapman & Hall.
#' Fox, J. (2016). Applied regression analysis and generalized linear models
#' (3rd ed.). SAGE Publications.

#' Helper function: Huber's Proposal 2 scale
#' @keywords internal
#' @noRd
.huber_proposal2_scale <- function(y) {
  # Iteratively determines optimal scale for M-estimation
  # Based on finding scale k that balances efficiency and robustness
  
  # Step 1: Start with MAD estimate
  k0 <- median(abs(y - median(y)), na.rm = TRUE) * 1.345
  if (k0 == 0) k0 <- 1
  
  # Step 2: Huber's Proposal 2 iteration
  # Find k such that E[rho(y/k)] ≈ 0.5 (Huber's requirement)
  # Using 4-point iteration
  for (iter in seq_len(10)) {
    # Weight using Huber function
    weights <- ifelse(abs(y/k0) <= 1, 1, 1/pmax(abs(y/k0), 0.01))
    
    # Update scale using M-estimate of scale
    c2 <- sum(weights * pmin((y/k0)^2, 1), na.rm = TRUE) / sum(weights, na.rm = TRUE)
    
    # Update k to target b = 0.5 (Huber's choice)
    k_new <- k0 * sqrt(c2)
    
    # Check convergence
    if (abs(k_new - k0) / (k0 + 1e-6) < 1e-4) break
    k0 <- k_new
  }
  
  return(pmax(k0, 1e-6))
}

# Helper function: S-estimator scale
.s_estimator_scale <- function(y, b = 0.5, max_iter = 20) {
  # S-estimator: minimizes scale with high breakdown point
  # Target: E[rho(y/s)] = b (typically b=0.5)
  
  # Start with initial scale
  s0 <- median(abs(y - median(y)), na.rm = TRUE) * 1.345
  if (s0 == 0) s0 <- 1
  
  # Minimize s such that mean Huber rho = b
  for (iter in seq_len(max_iter)) {
    # Residuals standardized by current scale
    u <- y / s0
    
    # Huber loss: rho(u) = u^2/2 if |u|<=1, |u|-0.5 if |u|>1
    rho <- ifelse(abs(u) <= 1, u^2 / 2, abs(u) - 0.5)
    
    # Mean rho should equal b
    mean_rho <- mean(rho, na.rm = TRUE)
    
    # Check convergence
    if (abs(mean_rho - b) < 0.01) break
    
    # Update s (bisection-like approach)
    if (mean_rho > b) {
      # Scale too small, increase it
      s0 <- s0 * 1.1
    } else {
      # Scale too large, decrease it
      s0 <- s0 * 0.9
    }
  }
  
  return(pmax(s0, 1e-6))
}

# ============================================================================
# MODULAR IRLS CORE (Approach 3)
# ============================================================================
# Single vector robust location estimation using Iteratively Re-Weighted
# Least Squares. Extracted into standalone function for reuse in:
#   1. m_estimate() - full diagnostic output
#   2. calculate_difference() - efficient group summaries
#
# @param y Numeric vector of observations (may contain NA)
# @param loss_type Character: "huber" (default), "tukey", or "lsq"
# @param scale Numeric scale parameter (default: auto-computed)
# @param scale_method Character: "mad" (default), "proposal2", "s-estimator"
# @param max_iter Integer: max IRLS iterations (default: 50)
# @param tol Numeric: convergence tolerance (default: 1e-6)
# @param return_weights Logical: if TRUE, return final weights; else just estimate
#
# @return If return_weights=FALSE: location_diff (numeric scalar)
#         If return_weights=TRUE: list(location_diff, weights, scale_used)
#
# @keywords internal
.irls_estimate_location <- function(y, loss_type = "huber", scale = NULL,
                                      scale_method = "mad", max_iter = 50,
                                      tol = 1e-6, return_weights = FALSE) {
  # Input validation
  y <- as.numeric(y)
  n_obs <- length(y)
  
  if (n_obs < 1) {
    return(if (return_weights) 
      list(location_diff = NA, weights = numeric(0), scale_used = NA)
    else NA)
  }
  
  # Handle all-NA case
  if (all(is.na(y))) {
    return(if (return_weights)
      list(location_diff = NA, weights = rep(NA, n_obs), scale_used = NA)
    else NA)
  }
  
  # Determine scale if not provided
  if (is.null(scale)) {
    if (scale_method == "proposal2") {
      scale_local <- .huber_proposal2_scale(y)
    } else if (scale_method == "s-estimator") {
      scale_local <- .s_estimator_scale(y)
    } else {
      # Default to MAD
      mad_y <- median(abs(y - median(y, na.rm = TRUE)), na.rm = TRUE)
      scale_local <- 1.345 * mad_y
      if (scale_local == 0) scale_local <- 1
    }
  } else {
    scale_local <- scale
  }
  
  # Initialize IRLS
  weights <- rep(1, n_obs)
  location_prev <- median(y, na.rm = TRUE)
  
  # IRLS iterations for single location (intercept model)
  for (iter in seq_len(max_iter)) {
    # Weighted mean as location estimate
    valid_idx <- !is.na(y)
    if (sum(weights[valid_idx]) > 0) {
      location <- sum(weights[valid_idx] * y[valid_idx], na.rm = TRUE) / 
                   sum(weights[valid_idx], na.rm = TRUE)
    } else {
      location <- median(y, na.rm = TRUE)
    }
    
    # Compute residuals and standardize
    residuals <- y - location
    standardized_resid <- residuals / scale_local
    
    # Update weights based on loss function
    if (loss_type == "huber") {
      abs_resid <- abs(standardized_resid)
      abs_resid[is.na(abs_resid)] <- 0
      weights <- ifelse(abs_resid <= 1, 1, 1 / pmax(abs_resid, 0.01))
    } else if (loss_type == "tukey") {
      abs_resid <- abs(standardized_resid)
      abs_resid[is.na(abs_resid)] <- 0
      weights <- ifelse(abs_resid <= 1, (1 - abs_resid^2)^2, 0)
    } else if (loss_type == "lsq") {
      weights <- rep(1, n_obs)
    }
    
    # Check convergence
    if (iter > 1) {
      diff_current <- abs(location - location_prev)
      if (!is.na(diff_current) && !is.infinite(diff_current) && diff_current < tol) {
        break
      }
    }
    location_prev <- location
  }
  
  # Return based on request
  if (return_weights) {
    return(list(location_diff = location, weights = weights, scale_used = scale_local))
  } else {
    return(location)
  }
}

# ============================================================================
# M-ESTIMATION FOR ROBUST GROUP COMPARISON
# ============================================================================

#' M-Estimation for Robust Location Comparison
#'
#' Estimates location differences between groups using M-estimation
#' (iteratively re-weighted least squares), which is more robust to
#' outliers than standard least squares.
#'
#' @param x Matrix of values (rows = features, columns = samples), or a
#'   SummarizedExperiment object with multi-q entropy data
#' @param samples Character vector indicating group membership. If x is a
#'   SummarizedExperiment, this should be a column name in colData.
#'   For multi-q data, can also specify "multi_q_analysis" to automatically
#'   handle q-value collapsing and leave-one-out influence analysis.
#' @param loss_type Type of loss function: "huber" (default, robust),
#'        "tukey" (more aggressive), or "lsq" (least squares, for comparison)
#' @param scale Numeric. Scale parameter for Huber loss (default: 1.345*MAD).
#'        Controls how much weight is given to outliers.
#' @param max_iter Integer. Maximum iterations for IRLS. Default: 50
#' @param tol Numeric. Convergence tolerance. Default: 1e-6
#' @param paired Logical. If TRUE, use paired design. Default: FALSE
#' @param pcorr P-value correction method. Default: "BH"
#' @param q_combine_method Character. For multi-q data: "mean" (default) or 
#'   "median" for summarizing across q values
#' @param influence_threshold Numeric. Quantile threshold (0-1) for flagging high-influence
#'   samples in multi-q analysis. Default: 0.75 (75th percentile)
#' @param scale_method Character. Scale selection method: "mad" (default, Median Absolute Deviation),
#'   "proposal2" (Huber's Proposal 2 for automatic scale selection), or 
#'   "s-estimator" (S-estimator for high breakdown point). Default: "mad"
#'
#' @return Data frame with columns:
#'   - location_diff: Estimated location difference (from M-estimation)
#'   - se_diff: Standard error of difference
#'   - t_stat: t-statistic
#'   - pvalue: Two-tailed p-value
#'   - padj: Adjusted p-value
#'   - n_down_weighted: Number of observations down-weighted as outliers
#'   - max_weight: Maximum weight assigned (1 = no down-weighting)
#'   
#'   For multi-q analysis on SummarizedExperiment, returns sample-level
#'   influence scores (proportion of genes with >2% change when sample removed).
#'
#' @references
#' Huber, P. J. (1981). Robust Statistics. John Wiley & Sons.
#' Maronna, R. A., Martin, R. D., & Yohai, V. J. (2006).
#' Robust Statistics: Theory and Methods. John Wiley & Sons.
#' Lopuhaä, H. P., & Rousseeuw, P. J. (1991). Breakdown points of affine equivariant 
#' estimators of multivariate location and covariance matrices. Annals of Statistics, 19(1), 229-248.
#'
#' @keywords internal
#' @noRd
#' @details
#' M-estimation uses the Huber loss function by default:
#' L(u) = u^2/2 if |u| <= k (quadratic, like LSQ)
#' L(u) = k|u| - k^2/2 if |u| > k (linear, like absolute value)
#'
#' This provides a compromise: near the center, it's as efficient as LSQ,
#' but observations far from the center (outliers) have reduced influence.
#'
#' The default scale k = 1.345 * MAD detects outliers beyond 1.345 standard
#' deviations (scaled by the median absolute deviation).
#'
#' For multi-q SummarizedExperiment data, the function automatically:
#' 1. Extracts the multi-q entropy assay
#' 2. Collapses samples across q values (using mean or median)
#' 3. Performs leave-one-out influence analysis
#' 4. Returns sample influence scores
#'
#' **Scale Estimation Methods:**
#' - **mad (default):** Scale = 1.345 * MAD (Median Absolute Deviation).
#'   Fast, consistent for normal data. Detects outliers at ~1.345 sigma.
#'
#' - **proposal2:** Huber's Proposal 2. Iteratively selects optimal k 
#'   to balance efficiency and robustness. More adaptive but slower.
#'   Good for data with unknown error distribution.
#'
#' - **s-estimator:** S-estimator with high breakdown point (~50%).
#'   More robust to extreme contamination than M-estimation (~25%).
#'   Recommended when data contamination is suspected.
#'
m_estimate <- function(x, samples, loss_type = "huber", scale = NULL,
                       max_iter = 50, tol = 1e-6, paired = FALSE, pcorr = "BH",
                       q_combine_method = "mean", influence_threshold = 0.75,
                       scale_method = "mad", verbose = FALSE) {
  # Handle SummarizedExperiment input with multi-q analysis
  if (inherits(x, "SummarizedExperiment")) {
    entropy_matrix <- SummarizedExperiment::assay(x)
    sample_info <- SummarizedExperiment::colData(x)
    
    # Auto-detect paired from SE metadata if paired parameter is default FALSE
    if (!isTRUE(paired) && length(metadata(x)) > 0 && "paired" %in% names(metadata(x))) {
      paired_meta <- metadata(x)$paired
      if (is.logical(paired_meta) && length(paired_meta) == 1) {
        paired <- paired_meta
      }
    }
    
    # Get group assignment
    if (!(samples %in% colnames(sample_info))) {
      stop(sprintf("Column '%s' not found in colData", samples))
    }
    
    # Extract column - note: as.vector with mode can fail on certain S4 objects
    # Use fallback approach if direct conversion fails
    group_assignment <- tryCatch({
      as.vector(sample_info[[samples]])
    }, error = function(e) {
      col_val <- sample_info[[samples]]
      if (is.atomic(col_val)) {
        col_val
      } else {
        as.character(col_val)
      }
    })
    
    col_names <- colnames(entropy_matrix)
    sample_names_full <- sub("_q=.*$", "", col_names)
    unique_samples <- unique(sample_names_full)
    
    # Create sample-level entropy matrix (collapse across q values)
    entropy_by_sample <- matrix(0, nrow = nrow(entropy_matrix), ncol = length(unique_samples),
                                dimnames = list(rownames(entropy_matrix), unique_samples))
    
    for (j in seq_along(unique_samples)) {
      cols_for_sample <- which(sample_names_full == unique_samples[j])
      if (q_combine_method == "median") {
        entropy_by_sample[, j] <- apply(entropy_matrix[, cols_for_sample, drop = FALSE], 1, median)
      } else {
        entropy_by_sample[, j] <- rowMeans(entropy_matrix[, cols_for_sample, drop = FALSE])
      }
    }
    
    # Get group assignment for unique samples
    # Initialize directly without type inference to avoid issues with S4 objects
    group_assignment_unique <- rep(NA, length(unique_samples))
    for (j in seq_along(unique_samples)) {
      first_col_idx <- which(sample_names_full == unique_samples[j])[1]
      group_assignment_unique[j] <- group_assignment[first_col_idx]
    }
    
    # Perform leave-one-out influence analysis
    # First, calculate M-estimate with ALL samples as the baseline
    # NOTE: Use paired mode for full fit to leverage any pairing in the design
    m_est_full <- tryCatch({
      m_estimate(entropy_by_sample, samples = group_assignment_unique,
                 loss_type = loss_type, scale = scale,
                 max_iter = max_iter, tol = tol, paired = paired, pcorr = pcorr,
                 scale_method = scale_method)
    }, error = function(e) {
      stop(e)
    })

    sample_influence <- numeric(length(unique_samples))
    sample_robustness_weights <- numeric(length(unique_samples))
    entropy_means <- numeric(length(unique_samples))
    entropy_sds <- numeric(length(unique_samples))
    names(sample_influence) <- unique_samples
    names(sample_robustness_weights) <- unique_samples
    names(entropy_means) <- unique_samples
    names(entropy_sds) <- unique_samples
    
    # Store gene-level location_diff for debugging: list with one entry per sample
    gene_level_changes <- list()
    
    # Calculate DFBETA threshold based on sample size (standard regression diagnostics)
    # Threshold of 2/sqrt(n) is commonly used (see Fox 2016, Regression Diagnostics)
    n_samples <- length(unique_samples)
    dfbeta_threshold <- 2 / sqrt(max(n_samples, 2))
    
    for (i in seq_along(unique_samples)) {
      # Calculate entropy statistics for this sample FIRST (before any skips)
      sample_entropy_vals <- entropy_by_sample[, i]
      entropy_means[i] <- mean(sample_entropy_vals, na.rm = TRUE)
      entropy_sds[i] <- sd(sample_entropy_vals, na.rm = TRUE)
      
      entropy_subset <- entropy_by_sample[, -i, drop = FALSE]
      group_subset <- group_assignment_unique[-i]
      
      # Check if we still have both groups represented
      unique_groups_subset <- unique(group_subset)
      if (length(unique_groups_subset) < 2) {
        # If removing sample i leaves only one group, this sample has very HIGH influence
        # Set to maximum influence (all genes affected) 
        sample_influence[i] <- 1.0
        # For high-influence samples, use maximum observed robustness weight variation
        sample_robustness_weights[i] <- NA  # Mark as NA to handle separately later
        # Record that this sample's removal breaks group structure
        gene_level_changes[[unique_samples[i]]] <- data.frame(
          gene = rownames(entropy_by_sample),
          full_location_diff = NA_real_,
          loo_location_diff = NA_real_,
          dfbeta = NA_real_,
          reason = "Only one group remains after removing sample"
        )
        next
      }
      
      # Call m_estimate recursively on matrix data WITHOUT sample i
      # IMPORTANT: Use paired=FALSE for leave-one-out diagnostics
      # REASON: Removing one observation breaks paired structure. Influence diagnostics
      # must work on data as-is. The FULL fit uses paired=TRUE; comparison vs LOO is valid.
      # Bibliography reference: Standard statistical practice (Fox 2016, Cook & Weisberg 1982)
      m_est_subset <- m_estimate(entropy_subset, samples = group_subset,
                                 loss_type = loss_type, scale = scale,
                                 max_iter = max_iter, tol = tol, paired = FALSE, pcorr = pcorr,
                                 scale_method = scale_method)
      
      # Calculate standardized influence (DFBETA) for each gene
      # DFBETA = (coef_full - coef_loo) / SE(coef_full)
      # This scales the change by the precision/uncertainty of the estimate
      # See: Fox (2016), Regression Diagnostics; Cook & Weisberg (1982)
      dfbeta <- (m_est_full$location_diff - m_est_subset$location_diff) / 
                pmax(m_est_full$se_diff, 1e-6)  # Use full model's SE for standardization
      
      # Get gene names - use rownames if available, otherwise use rownames from m_est_full or generate
      gene_names <- rownames(entropy_by_sample)
      if (is.null(gene_names) || length(gene_names) == 0) {
        # Try to get from m_est_full result
        gene_names <- rownames(m_est_full)
      }
      if (is.null(gene_names) || length(gene_names) == 0) {
        # Generate default names if still missing
        gene_names <- paste0("Gene_", seq_len(nrow(entropy_by_sample)))
      }
      
      # Store gene-level changes for this sample
      gene_level_changes[[unique_samples[i]]] <- data.frame(
        gene = gene_names,
        full_location_diff = m_est_full$location_diff,
        loo_location_diff = m_est_subset$location_diff,
        full_se_diff = m_est_full$se_diff,
        dfbeta = dfbeta,
        exceeds_threshold = abs(dfbeta) > dfbeta_threshold,
        removed_sample_idx = i,
        removed_sample_name = unique_samples[i],
        full_n_normal = sum(group_assignment_unique == "normal"),
        full_n_tumor = sum(group_assignment_unique == "tumor"),
        loo_n_normal = sum(group_subset == "normal"),
        loo_n_tumor = sum(group_subset == "tumor"),
        stringsAsFactors = FALSE
      )
      
      # Influence = proportion of genes with problematic DFBETA values
      # DFBETA > threshold indicates substantive influence on that gene
      sample_influence[i] <- mean(abs(dfbeta) > dfbeta_threshold, na.rm = TRUE)
    }
    
    # Extract robustness weights from baseline M-estimate 
    # (how much was each sample downweighted in the full fit?)
    # For SummarizedExperiment input, m_est_full is a data.frame with sample-level results
    # We need to get per-gene weights from the baseline regression
    # As approximation, use the max_weight from the baseline fit
    if (is.data.frame(m_est_full) && "max_weight" %in% colnames(m_est_full)) {
      # In SummarizedExperiment mode, max_weight represents robustness per gene
      # Average across genes to get overall robustness
      mean_weight <- mean(m_est_full$max_weight, na.rm = TRUE)
      # Apply to all samples, then override NA (high-influence) samples with minimum weight
      sample_robustness_weights[is.na(sample_robustness_weights)] <- 
        min(m_est_full$max_weight, na.rm = TRUE)
      sample_robustness_weights[sample_robustness_weights == 0] <- mean_weight
    } else {
      # Default: all samples fully trusted, except high-influence samples get lower weight
      sample_robustness_weights[is.na(sample_robustness_weights)] <- 0.5
      sample_robustness_weights[sample_robustness_weights == 0] <- 1.0
    }
    
    # Calculate distance from each sample to its group centroid
    centroid_distances <- numeric(length(unique_samples))
    names(centroid_distances) <- unique_samples
    
    for (group in unique(group_assignment_unique)) {
      group_samples_idx <- which(group_assignment_unique == group)
      
      if (length(group_samples_idx) > 0) {
        # Compute group centroid as median entropy pattern
        if (length(group_samples_idx) == 1) {
          group_centroid <- entropy_by_sample[, group_samples_idx, drop = FALSE]
        } else {
          group_centroid <- apply(entropy_by_sample[, group_samples_idx, drop = FALSE], 1, median)
        }
        
        # Calculate Euclidean distance from each sample in group to centroid
        for (sample_idx in group_samples_idx) {
          sample_entropy <- entropy_by_sample[, sample_idx]
          if (is.matrix(group_centroid)) {
            dist <- sqrt(sum((sample_entropy - group_centroid[, 1])^2, na.rm = TRUE))
          } else {
            dist <- sqrt(sum((sample_entropy - group_centroid)^2, na.rm = TRUE))
          }
          centroid_distances[sample_idx] <- dist
        }
      }
    }
    
    # Identify problematic samples (using specified threshold)
    high_influence_threshold <- quantile(sample_influence, influence_threshold, na.rm = TRUE)
    
    # Extract condition/group information for each unique sample
    condition_for_samples <- character(length(unique_samples))
    for (j in seq_along(unique_samples)) {
      first_col_idx <- which(sample_names_full == unique_samples[j])[1]
      condition_for_samples[j] <- as.character(group_assignment[first_col_idx])
    }
    
    # Extract paired sample information if available
    paired_sample_info <- NULL
    # Check for various column name conventions
    pair_col_names <- c("pair_id", "paired_samples", "Pair_ID", "Pair")
    pair_col <- NULL
    for (col in pair_col_names) {
      if (col %in% colnames(sample_info)) {
        pair_col <- col
        break
      }
    }
    
    if (!is.null(pair_col)) {
      paired_sample_info <- character(length(unique_samples))
      for (j in seq_along(unique_samples)) {
        first_col_idx <- which(sample_names_full == unique_samples[j])[1]
        pair_id <- sample_info[first_col_idx, pair_col]
        paired_sample_info[j] <- as.character(pair_id)
      }
    }
    
    # Return influence scores as data frame with condition, robustness, entropy info, and centroid distance
    result_df <- data.frame(
      Sample = names(sample_influence),
      Condition = condition_for_samples,
      Proportion_Affected = sample_influence,
      Genes_Affected = sample_influence * nrow(entropy_by_sample),
      Robustness_Weight = sample_robustness_weights,
      Entropy_Mean = entropy_means,
      Entropy_SD = entropy_sds,
      Distance_from_Centroid = centroid_distances,
      Status = ifelse(sample_influence > high_influence_threshold, "Flag for QC", "OK"),
      stringsAsFactors = FALSE,
      row.names = NULL
    )
    
    # Add paired sample information if available
    if (!is.null(paired_sample_info)) {
      result_df$Pair_ID <- paired_sample_info
    }
    return(result_df)
  }
  
  # Original matrix-based m_estimate logic

  # Input validation
  if (!is.matrix(x) && !is.data.frame(x)) {
    x <- as.matrix(x)
  }

  if (length(samples) != ncol(x)) {
    stop("Length of samples must equal number of columns in x")
  }

  groups <- unique(samples)
  if (length(groups) != 2) {
    stop("Must have exactly 2 groups")
  }

  if (!(loss_type %in% c("huber", "tukey", "lsq"))) {
    stop("loss_type must be 'huber', 'tukey', or 'lsq'")
  }

  if (!(scale_method %in% c("mad", "proposal2", "s-estimator"))) {
    stop("scale_method must be 'mad', 'proposal2', or 's-estimator'")
  }

  n_features <- nrow(x)
  results_list <- list()

  # Process each feature
  for (i in seq_len(n_features)) {
    feature_vals <- as.numeric(x[i, ])
    group1_idx <- which(samples == groups[1])
    group2_idx <- which(samples == groups[2])

    group1 <- feature_vals[group1_idx]
    group2 <- feature_vals[group2_idx]

    # Prepare data for regression
    if (paired) {
      # Paired: use differences
      if (length(group1) != length(group2)) {
        min_n <- min(length(group1), length(group2))
        group1 <- group1[seq_len(min_n)]
        group2 <- group2[seq_len(min_n)]
      }
      y <- group1 - group2
      X <- rep(1, length(y))
      use_intercept <- FALSE
    } else {
      # Unpaired: use group indicator
      y <- c(group1, group2)
      X <- c(rep(0, length(group1)), rep(1, length(group2)))
      use_intercept <- TRUE
    }

    n_obs <- length(y)

    # IRLS regression to estimate location_diff (slope)
    # For unpaired: y ~ intercept + location_diff * X
    # For paired: y ~ location_diff (no intercept, since y = group1 - group2)
    
    # Determine scale
    if (is.null(scale)) {
      if (scale_method == "proposal2") {
        scale_local <- .huber_proposal2_scale(y)
      } else if (scale_method == "s-estimator") {
        scale_local <- .s_estimator_scale(y)
      } else {
        mad_y <- median(abs(y - median(y, na.rm = TRUE)), na.rm = TRUE)
        scale_local <- 1.345 * mad_y
        if (scale_local == 0) scale_local <- 1
      }
    } else {
      scale_local <- scale
    }
    
    # GUARD: If scale is NA or infinite (but allow scale <= 0 to fallback to 1)
    if (is.na(scale_local) || !is.finite(scale_local)) {
      results_list[[i]] <- data.frame(
        location_diff = NA_real_,
        se_diff = NA_real_,
        t_stat = NA_real_,
        pvalue = NA_real_,
        n_down_weighted = NA_integer_,
        max_weight = NA_real_,
        row.names = rownames(x)[i]
      )
      next
    }
    
    # Final safeguard: if scale is still 0, set to 1
    if (scale_local == 0) {
      scale_local <- 1
    }
    
    # Initialize regression coefficients
    if (use_intercept) {
      # Start with group medians
      intercept <- median(y[X == 0], na.rm = TRUE)
      location_diff <- median(y[X == 1], na.rm = TRUE) - intercept
      coef <- c(intercept, location_diff)
    } else {
      # Paired: just estimate the mean of differences
      location_diff <- median(y, na.rm = TRUE)
      coef <- location_diff
    }
    
    weights <- rep(1, n_obs)
    
    # IRLS iterations for regression
    for (iter in seq_len(max_iter)) {
      # Compute fitted values
      if (use_intercept) {
        fitted <- coef[1] + coef[2] * X
      } else {
        fitted <- coef[1] * X  # For paired, X = 1 for all
      }
      
      # Residuals
      residuals <- y - fitted
      
      # GUARD: If all residuals are NA, skip this gene
      if (all(is.na(residuals))) {
        results_list[[i]] <- data.frame(
          location_diff = NA_real_,
          se_diff = NA_real_,
          t_stat = NA_real_,
          pvalue = NA_real_,
          n_down_weighted = NA_integer_,
          max_weight = NA_real_,
          row.names = rownames(x)[i]
        )
        break
      }
      
      standardized_resid <- residuals / scale_local
      
      # Compute weights based on loss function
      if (loss_type == "huber") {
        weights <- ifelse(abs(standardized_resid) <= 1, 1, 1 / pmax(abs(standardized_resid), 0.01))
      } else if (loss_type == "tukey") {
        weights <- ifelse(abs(standardized_resid) <= 1, (1 - standardized_resid^2)^2, 0)
      } else if (loss_type == "lsq") {
        weights <- rep(1, n_obs)
      }
      
      # SAFEGUARD: If all weights are NA or non-finite, skip this gene
      if (all(is.na(weights)) || all(!is.finite(weights))) {
        # This gene has no valid residuals/scale - skip to next gene
        results_list[[i]] <- data.frame(
          location_diff = NA_real_,
          se_diff = NA_real_,
          t_stat = NA_real_,
          pvalue = NA_real_,
          n_down_weighted = NA_integer_,
          max_weight = NA_real_,
          row.names = rownames(x)[i]
        )
        next
      }
      
      # Update coefficients using weighted least squares
      if (use_intercept) {
        # Weighted regression: y ~ 1 + X
        X_design <- cbind(1, X)
        XtWX <- crossprod(X_design, weights * X_design)
        XtWy <- crossprod(X_design, weights * y)
        
        # Solve for coefficients
        tryCatch({
          coef_new <- solve(XtWX, XtWy)
          coef_prev <- coef
          coef <- as.numeric(coef_new)
          
          # Check convergence
          if (max(abs(coef - coef_prev)) < tol) break
        }, error = function(e) {
          # Singular matrix, keep current coefficients
        })
      } else {
        # Paired: just estimate mean of differences
        valid_idx <- !is.na(y) & is.finite(weights)
        if (sum(weights[valid_idx]) > 0) {
          coef_new <- sum(weights[valid_idx] * y[valid_idx]) / sum(weights[valid_idx])
          if (abs(coef_new - coef) < tol) break
          coef <- coef_new
        }
      }
    }
    
    # Extract final estimates
    if (use_intercept) {
      intercept <- coef[1]
      location_diff <- coef[2]
      fitted <- intercept + location_diff * X
    } else {
      location_diff <- coef[1]
      fitted <- location_diff * X
    }
    
    residuals <- y - fitted
    
    # Estimate of error variance
    rss <- sum(weights * residuals^2, na.rm = TRUE)
    sigma_sq <- rss / max(1, n_obs - 2)
    if (is.na(sigma_sq) || is.infinite(sigma_sq) || sigma_sq <= 0) {
      sigma_sq <- var(y, na.rm = TRUE)
      if (is.na(sigma_sq) || sigma_sq <= 0) sigma_sq <- 1
    }

    # Standard error depends on design
    if (use_intercept) {
      # SE for slope in grouped design
      w_mean_x <- weighted.mean(X, weights, na.rm = TRUE)
      sx_squared <- sum(weights * (X - w_mean_x)^2, na.rm = TRUE)
      if (sx_squared > 0) {
        se_diff <- sqrt(sigma_sq / sx_squared)
      } else {
        se_diff <- sqrt(sigma_sq / length(group1) + sigma_sq / length(group2))
      }
    } else {
      # SE for mean of differences
      se_diff <- sqrt(sigma_sq / n_obs)
    }

    # t-statistic and p-value
    if (se_diff > 0) {
      t_stat <- location_diff / se_diff
      df <- n_obs - 2
      pvalue <- 2 * pt(abs(t_stat), df = max(1, df), lower.tail = FALSE)
    } else {
      t_stat <- Inf
      pvalue <- if (abs(location_diff) > 0) 0 else 1
    }

    # Count outliers
    n_down_weighted <- sum(weights < 0.99, na.rm = TRUE)
    max_weight <- max(weights, na.rm = TRUE)

    results_list[[i]] <- data.frame(
      location_diff = location_diff,
      se_diff = se_diff,
      t_stat = t_stat,
      pvalue = pvalue,
      n_down_weighted = n_down_weighted,
      max_weight = max_weight,
      row.names = rownames(x)[i]
    )
  }

  # Combine and adjust p-values
  result <- do.call(rbind, results_list)
  result$padj <- p.adjust(result$pvalue, method = pcorr)

  return(result)
}

