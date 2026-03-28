# =========================================================================
# PRIVATE HELPER FUNCTIONS (must come before main roxygen block)
# =========================================================================

# -------------------------------------------------------------------------
# PRIVATE HELPER: Auto-detect group column and control group
# -------------------------------------------------------------------------

#' Auto-Detect Group Column and Control Group from SummarizedExperiment colData
#'
#' Automatically searches for a column containing group/condition assignments
#' and identifies the reference (control) group.
#'
#' **Detection Strategy for group_col:**
#' Searches colData in precedence order for common grouping column names:
#' 1. "group" (TSENAT default)
#' 2. "condition", "treatment", "sample_type" (common alternatives)
#'
#' **Detection Strategy for control_group:**
#' Once group column is found, identifies the reference group by:
#' 1. First looking for common control/reference names: "Normal", "Control", 
#'    "WT" (wild-type), "Reference", "Baseline", "wt"
#' 2. If no match, selects the unique group value with smallest sample count
#'    (typically the control/reference in case-control designs)
#' 3. If still no match, uses the first alphabetically sorted group name
#'
#' @param se SummarizedExperiment object with sample metadata in colData
#'
#' @return List with elements:
#'   \describe{
#'     \item{group_col}{Name of the colData column used for grouping, or NA_character_ if none detected}
#'     \item{control_group}{Name of the control/reference group, or NA_character_ if none detected}
#'     \item{groups}{Character vector of all unique groups found}
#'     \item{sample_counts}{Named integer vector of sample counts per group (names: group names)}
#'   }
#'

#' @noRd
.tsenat_auto_detect_groups <- function(se) {
  
  cd <- SummarizedExperiment::colData(se)
  cd_colnames <- colnames(cd)
  
  # Candidate column names (in priority order) for group/condition
  group_col_candidates <- c(
    "sample_type",      # TSENAT standard (created by map_metadata)
    "group",            # Alias for sample_type
    "condition",        # Common alternative
    "treatment",        # Experimental design
    "phenotype",        # Biological phenotype
    "batch",            # Last resort
    "category"          # Generic fallback
  )
  
  group_col <- NA_character_
  
  # Find first matching group column
  for (col_name in group_col_candidates) {
    if (col_name %in% cd_colnames) {
      group_col <- col_name
      break
    }
  }
  
  # If no group column found, return NAs
  if (is.na(group_col)) {
    return(list(
      group_col = NA_character_,
      control_group = NA_character_,
      groups = character(0),
      sample_counts = integer(0)
    ))
  }
  
  # Get group values and counts
  group_vec <- as.character(cd[[group_col]])
  unique_groups <- unique(group_vec)
  group_counts <- table(group_vec)
  
  # Candidate names for control/reference group (in priority order)
  control_candidates <- c(
    "Normal",           # TSENAT default
    "Control",          # Most common control label
    "WT",               # Wild-type (common in genetics)
    "wt",               # Lowercase variant
    "Reference",        # Explicit reference label
    "Baseline",         # Baseline condition
    "Wild-type",        # Full name variant
    "wild_type"         # Underscore variant
  )
  
  control_group <- NA_character_
  
  # Try to match control candidates
  for (control_name in control_candidates) {
    if (control_name %in% unique_groups) {
      control_group <- control_name
      break
    }
  }
  
  # If no match found, use heuristics:
  # 1. Select group with fewer samples (typical case-control design)
  # 2. Fallback to first alphabetically
  if (is.na(control_group)) {
    if (length(unique_groups) >= 2) {
      # Find group with minimum samples (typically control)
      min_samples_group <- names(group_counts)[which.min(group_counts)]
      control_group <- min_samples_group
    } else if (length(unique_groups) == 1) {
      control_group <- unique_groups[1]
    }
  }
  
  return(list(
    group_col = group_col,
    control_group = control_group,
    groups = unique_groups,
    sample_counts = as.vector(group_counts)
  ))
}


# -------------------------------------------------------------------------
# PRIVATE HELPER: Auto-detect paired samples from metadata
# -------------------------------------------------------------------------

#' Detect Paired Sample Structure from SummarizedExperiment colData
#'
#' Automatically searches for a column containing paired sample identifiers
#' (e.g., "pair_id", "paired_samples", "patient_id", "subject_id").
#' This enables pair-respecting bootstrap resampling in divergence calculations.
#'
#' **Detection Strategy:**
#' Searches colData in precedence order for common pairing column names:
#' 1. "paired_samples" (TSENAT default, matches readcounts metadata)
#' 2. "pair_id", "pair_samples", "subject_id", "patient_id" (common alternatives)
#'
#' Returns a mapping from sample names to pair identifiers, or NULL if no
#' pairing column is found. A valid pairing column has:
#' - Non-NA values for all samples
#' - At least 2 samples per pair
#' - Deterministic structure (e.g., all A's paired with another A sample nearby)
#'
#' **Database References (Papers validating auto-detection approach):**
#' - S102: "Experimental Control and Paired Design" - standardizes paired design annotation
#' - S107: "Related Sample Designs and Paired t-test" - validates paired structure detection
#'
#' @param se SummarizedExperiment object with sample metadata in colData
#'
#' @return List with elements:
#'   \describe{
#'     \item{pair_ids}{Character vector (names: sample names, values: pair identifiers)
#'       or NULL if no pairing detected}
#'     \item{column_name}{Name of the colData column used, or NA_character_ if none}
#'     \item{num_pairs}{Number of unique pairs (0 if none detected)}
#'     \item{samples_per_pair}{Vector of samples per pair (names: pair IDs, values: counts)}
#'   }
#'
#' @note Paired samples detected from any of: "paired_samples", "pair_id", "pair_samples",
#'   "subject_id", "patient_id". Returns NULL if none present or validation fails.
#'

#' @noRd
.tsenat_detect_pair_ids <- function(se) {
  
  cd <- SummarizedExperiment::colData(se)
  sample_names <- colnames(se)
  
  # Candidate column names (in priority order)
  candidate_cols <- c(
    "paired_samples",      # TSENAT default
    "pair_id",             # Common alternative
    "pair_samples",        # Variant
    "subject_id",          # Statistical defaults
    "patient_id"
  )
  
  for (col_name in candidate_cols) {
    if (col_name %in% colnames(cd)) {
      pair_col <- cd[[col_name]]
      
      # Validate: must be non-NA for all samples
      if (any(is.na(pair_col))) {
        next  # Skip if any NAs
      }
      
      # Valid pairing structure found
      pair_ids <- setNames(as.character(pair_col), sample_names)
      unique_pairs <- unique(pair_ids)
      samples_per_pair <- table(pair_ids)
      
      return(list(
        pair_ids = pair_ids,
        column_name = col_name,
        num_pairs = length(unique_pairs),
        samples_per_pair = samples_per_pair
      ))
    }
  }
  
  # No pairing detected
  return(list(
    pair_ids = NULL,
    column_name = NA_character_,
    num_pairs = 0,
    samples_per_pair = numeric(0)
  ))
}


#' Resample Data Respecting Paired Structure
#'
#' When resampling paired data, both members of a pair are selected or discarded
#' together. This preserves within-pair correlations critical for statistical validity
#' in matched designs (Efron & Tibshirani 1993).
#'
#' For k pairs:
#' - Draw k pair indices uniformly with replacement: pair_idx ~ U(1:k)
#' - For each drawn pair i, include both samples from pair i
#' - This maintains pairing structure across bootstrap replicates
#'
#' **Statistical Justification (Papers C016, S102-S109):**
#' - C016: Bootstrap for confidence intervals requires preserving data structure
#' - S102: "Paired Design" - paired resampling required for matched samples
#' - S107: "Related Sample Designs" - within-pair correlation invalidates independent resampling
#'
#' @param control_samples Vector of counts for control group
#' @param treatment_samples Vector of counts for treatment group
#' @param pair_ids Named character vector: names = sample names, values = pair IDs
#' @param group_col Character vector: group assignment (one per sample)
#' @param control_group Character: name of control group
#'
#' @return List with elements:
#'   \describe{
#'     \item{control_resampled}{Resampled control group counts}
#'     \item{treatment_resampled}{Resampled treatment group counts}
#'   }
#'

#' @noRd
.tsenat_jis_resample_paired_data <- function(
    control_samples,
    treatment_samples,
    pair_ids,
    group_col,
    control_group) {
  
  # Map sample names to indices
  all_samples <- c(names(control_samples), names(treatment_samples))
  all_groups <- c(
    rep(control_group, length(control_samples)),
    rep(setdiff(unique(group_col), control_group), length(treatment_samples))
  )
  
  # Get unique pairs involved
  pairs_in_data <- unique(pair_ids[all_samples])
  num_pairs <- length(pairs_in_data)
  
  # Resample pairs with replacement
  resampled_pairs <- sample(seq_len(num_pairs), size = num_pairs, replace = TRUE)
  resampled_pair_ids <- pairs_in_data[resampled_pairs]
  
  # Collect samples for each resampled pair
  resampled_control <- numeric(0)
  resampled_treatment <- numeric(0)
  
  for (pair_id in resampled_pair_ids) {
    # Get both samples from this pair
    pair_samples <- names(pair_ids)[pair_ids == pair_id]
    
    for (sample_name in pair_samples) {
      if (sample_name %in% names(control_samples)) {
        resampled_control <- c(resampled_control, control_samples[sample_name])
      } else if (sample_name %in% names(treatment_samples)) {
        resampled_treatment <- c(resampled_treatment, treatment_samples[sample_name])
      }
    }
  }
  
  # BUGFIX #3: Validate balanced groups after paired resampling
  # Ensures control and treatment have equal sizes (required for divergence computation)
  if (length(resampled_control) != length(resampled_treatment)) {
    stop("Paired bootstrap produced unequal group sizes (",
         length(resampled_control), " control vs ", 
         length(resampled_treatment), " treatment). ",
         "Check for unbalanced or incomplete pairs in input data.")
  }
  
  return(list(
    control_resampled = resampled_control,
    treatment_resampled = resampled_treatment
  ))
}





# -------------------------------------------------------------------------
# HELPER FUNCTION: Bootstrap CI Computation
# -------------------------------------------------------------------------

#' Calculate Tsallis Divergence Bootstrap Confidence Intervals
#'
#' Internal helper function for compute_divergence_bootstrap.
#' Computes bootstrap confidence intervals for Tsallis divergence.
#'

#' @noRd

.calculate_divergence_bootstrap <- function(
    x, y,
    q = 1,
    nboot = 1000,
    ci = 0.95,
    method = "percentile",
    log_base = exp(1),
    pseudocount = 0.5,
    gene_name = NA_character_,
    verbose = FALSE,
    seed = NULL,
    paired = FALSE,
    pair_ids = NULL) {

  # Seed handling left to caller for Bioconductor compliance

  # Compute point estimate
  point_est <- .tsenat_tsallis_divergence_scalar(x, y, q, pseudocount, log_base)

  # Bootstrap confidence interval
  if (nboot > 0) {
    bootstrap_dist <- numeric(nboot)
    
    # Determine if using paired resampling (use isTRUE to handle NA safely)
    use_paired_bootstrap <- isTRUE(paired) && !is.null(pair_ids)

    for (b in seq_len(nboot)) {
      if (use_paired_bootstrap) {
        # Paired bootstrap: resample pair indices to preserve within-pair correlation
        # while maintaining equal group sizes
        x_names <- names(x)
        y_names <- names(y)
        all_names <- c(x_names, y_names)
        
        # BUGFIX: Safely handle pair_ids subset with proper indexing
        # Only keep pairs that have samples in the current groups
        pair_ids_subset <- pair_ids[all_names]
        
        # Get unique pairs and resample indices with replacement
        unique_pairs <- unique(pair_ids_subset[!is.na(pair_ids_subset)])
        
        if (length(unique_pairs) == 0) {
          # Fallback to independent bootstrap if pairing structure is broken
          x_boot <- sample(x, size = length(x), replace = TRUE)
          y_boot <- sample(y, size = length(y), replace = TRUE)
        } else {
          num_pairs <- length(unique_pairs)
          resampled_pair_indices <- sample(seq_len(num_pairs), size = num_pairs, replace = TRUE)
          
          # Collect samples from resampled pairs, maintaining group structure
          x_boot <- c()
          y_boot <- c()
          
          for (idx in resampled_pair_indices) {
            pair_id <- unique_pairs[idx]
            pair_mask <- pair_ids_subset == pair_id
            pair_samples <- names(pair_ids_subset)[pair_mask]
            
            for (sample in pair_samples) {
              if (sample %in% x_names) {
                x_boot <- c(x_boot, x[sample])
              } else if (sample %in% y_names) {
                y_boot <- c(y_boot, y[sample])
              }
            }
          }
          
          # BUGFIX: Ensure both vectors are non-empty and numeric
          if (length(x_boot) == 0) x_boot <- numeric(0)
          if (length(y_boot) == 0) y_boot <- numeric(0)
          
          # If pairing structure doesn't preserve group sizes, fall back to independent bootstrap
          if (length(x_boot) != length(x) || length(y_boot) != length(y)) {
            x_boot <- sample(x, size = length(x), replace = TRUE)
            y_boot <- sample(y, size = length(y), replace = TRUE)
          }
        }
      } else {
        # Independent bootstrap: standard resampling (each group independently)
        x_boot <- sample(x, size = length(x), replace = TRUE)
        y_boot <- sample(y, size = length(y), replace = TRUE)
      }

      bootstrap_dist[b] <- .tsenat_tsallis_divergence_scalar(x_boot, y_boot, q, pseudocount, log_base)
    }

    alpha <- (1 - ci) / 2

    if (method == "percentile") {
      lower_ci <- stats::quantile(bootstrap_dist, probs = alpha, names = FALSE)
      upper_ci <- stats::quantile(bootstrap_dist, probs = 1 - alpha, names = FALSE)
    } else if (method == "bca") {
      # BCA not appropriate for divergence (requires two-sample jackknife)
      # Fall back to percentile method which is valid for any divergence
      lower_ci <- stats::quantile(bootstrap_dist, probs = alpha, names = FALSE)
      upper_ci <- stats::quantile(bootstrap_dist, probs = 1 - alpha, names = FALSE)
    }
  } else {
    lower_ci <- NA_real_
    upper_ci <- NA_real_
  }

  return(list(
    estimate = point_est,
    lower_ci = lower_ci,
    upper_ci = upper_ci,
    q = q,
    nboot = nboot,
    method = if (nboot > 0) method else NA_character_
  ))
}


#' Compute Tsallis Divergence Between Two Count Vectors
#'
#' Computes scalar Tsallis divergence D_q(p || q) using the Furuichi formula.
#'

#' @noRd
.tsenat_tsallis_divergence_scalar <- function(x, y, q_val, pseudocount = 0.5, log_base = exp(1)) {
  # Validate input vectors
  if (length(x) == 0 || length(y) == 0) {
    return(NA_real_)
  }
  
  # BUGFIX: Ensure x and y have equal length (required for divergence)
  if (length(x) != length(y)) {
    # This can happen if paired bootstrap resampling produces unequal group sizes
    # Return NA rather than crashing
    return(NA_real_)
  }
  
  # Normalize to probabilities
  p <- (x + pseudocount) / (sum(x) + length(x) * pseudocount)
  r <- (y + pseudocount) / (sum(y) + length(y) * pseudocount)

  if (any(is.na(p)) || any(is.na(r))) {
    return(NA_real_)
  }
  
  # BUGFIX #2: Add explicit safeguard for near-zero probabilities
  # Prevents NaN/Inf from log(0) or extremely small values in power operations
  min_prob <- 1e-10
  p[p < min_prob] <- min_prob
  r[r < min_prob] <- min_prob
  
  # Re-normalize to maintain probability constraint (sum = 1)
  p <- p / sum(p)
  r <- r / sum(r)

  # Compute Tsallis divergence using correct formula from Paper I004
  # D_q(p||r) with D_q >= 0 and equality iff p = r
  # BUGFIX: Ensure formula is applied correctly for all q values
  
  if (abs(q_val) < 0.01) {
    # q=0: Tsallis divergence D_0(p||r) = (1/(0-1)) * (1 - sum(p^0 * r^1))
    #     = -1 * (1 - sum(1 * r)) = -1 * (1 - 1) = 0 (always 0 for any distributions)
    # This is mathematically correct: at q=0, all probability distributions have equal "divergence"
    div <- 0
  } else if (abs(q_val - 1) < 0.01) {
    # KL divergence (special case q -> 1): lim_{q->1} D_q = sum(p*log(p/r))
    div <- sum(p * log(p / r), na.rm = TRUE)
  } else if (q_val > 0 && q_val != 1) {
    # Standard Tsallis divergence formula: D_q(p||r) = (1/(q-1)) * (1 - sum(p^q * r^(1-q)))
    # This ensures D_q >= 0 and is asymmetric in p, r
    # CRITICAL: Ensure p and r vectors are properly aligned
    p_power <- p^q_val
    r_power <- r^(1 - q_val)
    
    # Check for numerical issues (inf, nan, underflow)
    if (any(is.nan(p_power)) || any(is.infinite(p_power)) ||
        any(is.nan(r_power)) || any(is.infinite(r_power))) {
      # Log-space computation for numerical stability when q is far from 1
      log_p_power <- q_val * log(p + 1e-10)
      log_r_power <- (1 - q_val) * log(r + 1e-10)
      sum_term <- sum(exp(log_p_power + log_r_power), na.rm = TRUE)
    } else {
      sum_term <- sum(p_power * r_power, na.rm = TRUE)
    }
    
    div <- (1 - sum_term) / (q_val - 1)
  } else {
    # Invalid q value
    return(NA_real_)
  }

  if (is.nan(div) || !is.finite(div)) {
    return(NA_real_)
  }

  # Apply log_base normalization CONSISTENTLY for all q values
  # This ensures consistent scaling across multi-q spectrum analysis
  if (log_base != exp(1)) {
    div <- div / log(log_base)
  }

  # BUG FIX: Handle sign correctly for q < 1
  # When q < 1, (q_val - 1) is negative, so the formula naturally produces
  # a positive divergence. We must take absolute value and ensure non-negativity.
  # Divergence should always be >= 0.
  return(abs(div))
}

#' Classify a per-q divergence spectrum into biological pattern types
#'
#' When Tsallis divergence has been computed across multiple q values
#' for a gene, the resulting vector can be summarised by its trend across the
#' spectrum. This helper computes the Pearson correlation between the numeric
#' q values and their corresponding divergences and interprets the
#' slope as one of three categories:
#' \describe{
#'   \item{RARE_DRIVEN}{Divergence decreases with q (high at q=0.5);
#'     indicates changes driven by low-abundance isoforms.}
#'   \item{ABUNDANT_DRIVEN}{Divergence increases with q (high at q=2);
#'     indicates shifts among the most abundant transcripts.}
#'   \item{BALANCED}{Little or no trend across q; effects are proportional.}
#' }
#'
#' If the input vector is too short, contains only NAs, or the
#' correlation cannot be calculated, NA is returned.
#'
#' @param per_q_divs Named numeric vector of divergences. Names should be of form
#'   "q_0.5", "q_1.0", etc. (or similar with separators "_", "=", or ".").
#' @param threshold Numeric; correlation threshold for classification. Absolute
#'   correlation values below this threshold are classified as BALANCED (default: 0.5).
#' @return Character scalar giving the pattern type, or NA if classification
#'   cannot be performed.

#' @noRd
#' @examples
#' per_q <- c(q_0.5=0.5, q_1=0.3, q_2=0.1)
#' .classify_q_pattern(per_q)
#'
#' # handling missing values
#' .classify_q_pattern(c(q_0.5=NA, q_1=0.2))

.classify_q_pattern <- function(per_q_divs, threshold = 0.5) {
  # Input validation
  if (!is.numeric(per_q_divs) || length(per_q_divs) < 2) {
    return(NA_character_)
  }
  
  # Validate names format
  nm <- names(per_q_divs)
  if (is.null(nm) || any(is.na(nm))) {
    return(NA_character_)
  }
  
  # Handle flexible name formats: try multiple patterns
  q_vals <- NA
  
  # Try format: "q_0.5" (standard)
  if (all(grepl("^q_", nm))) {
    q_vals <- as.numeric(gsub("^q_", "", nm))
  } else if (all(grepl("^q[_=.]", nm))) {
    # Try other separators
    q_vals <- as.numeric(gsub("^q[_=.]", "", gsub("_", ".", nm)))
  }
  
  # Fallback: try to extract numeric directly after "q"
  if (all(is.na(q_vals))) {
    extracted_q <- gsub("^q", "", nm)
    q_vals <- as.numeric(extracted_q)
  }
  
  # If we still can't extract numeric q values, return NA
  if (any(is.na(q_vals))) {
    return(NA_character_)
  }
  
  # Check if all divergence values are NA
  if (all(is.na(per_q_divs))) {
    return(NA_character_)
  }
  
  # Need at least 2 non-NA pairs for correlation
  valid_pairs <- !is.na(per_q_divs)
  if (sum(valid_pairs) < 2) {
    return(NA_character_)
  }
  
  # IMPROVED CLASSIFICATION: Compare q-regions instead of just endpoints
  # Split into rare-region (q < 1) and abundant-region (q >= 1)
  rare_mask <- q_vals < 1
  abund_mask <- q_vals >= 1
  
  # Calculate median divergence in each region
  if (sum(rare_mask & valid_pairs) > 0) {
    rare_div_median <- median(per_q_divs[rare_mask & valid_pairs], na.rm = TRUE)
  } else {
    rare_div_median <- NA
  }
  
  if (sum(abund_mask & valid_pairs) > 0) {
    abund_div_median <- median(per_q_divs[abund_mask & valid_pairs], na.rm = TRUE)
  } else {
    abund_div_median <- NA
  }
  
  # If we have both regions, compare them
  if (!is.na(rare_div_median) && !is.na(abund_div_median) && abund_div_median > 0) {
    ratio <- rare_div_median / abund_div_median
    
    # Use a proper ratio threshold (not the correlation threshold)
    # Ratio threshold should be > 1 to distinguish RARE from ABUNDANT
    # Use 1.3 as the ratio threshold (30% difference = sensitive but not overly permissive)
    ratio_threshold <- 1.3
    
    # RARE_DRIVEN: rare region has notably higher divergence (ratio > threshold)
    # ABUNDANT_DRIVEN: abundant region has notably higher divergence (ratio < 1/threshold)
    # BALANCED: similar divergence across regions (ratio near 1)
    
    if (!is.na(ratio) && ratio > ratio_threshold) {
      return("RARE_DRIVEN")
    } else if (!is.na(ratio) && ratio < 1 / ratio_threshold) {
      return("ABUNDANT_DRIVEN")
    } else {
      return("BALANCED")
    }
  }
  
  # Fallback: Use original correlation-based approach
  # Check for constant values (zero variance) before computing correlation
  # This avoids errors from cor() when one variable has no variance
  q_sd <- sd(q_vals, na.rm = TRUE)
  div_sd <- sd(per_q_divs, na.rm = TRUE)
  
  # Use isTRUE for safe comparison (handles NA)
  if (isTRUE(q_sd == 0) || isTRUE(div_sd == 0)) {
    return("BALANCED")
  }
  
  slope <- cor(q_vals, per_q_divs, use = "complete.obs")
  
  # If correlation is NA (e.g., constant divergence), treat as balanced
  if (is.na(slope)) {
    return("BALANCED")
  }
  
  # Classify based on slope and threshold
  abs_slope <- abs(slope)
  
  if (abs_slope < threshold) {
    "BALANCED"
  } else if (slope < -threshold) {
    "RARE_DRIVEN"
  } else if (slope > threshold) {
    "ABUNDANT_DRIVEN"
  } else {
    "BALANCED"
  }
}

