################################################################################
#
#' Internal: Core Tsallis entropy calculation (consolidated)
#' 
#' Centralized entropy calculation used by all functions to eliminate duplication.
#' Supports Shannon (q=1), Tsallis (q≠1), and species richness (q=0).
#'
#' @param proportions Numeric vector of species proportions (must sum to ~1)
#' @param q Numeric. Generalization parameter. Default: 1.0 (Shannon entropy)
#' @param norm Logical. Normalize by maximum entropy. Default: FALSE
#' @param log_base Numeric. Logarithm base. Default: exp(1) (natural log)
#' @param q_tol Numeric. Tolerance for detecting q=1 case. Default: 1e-6
#'
#' @return Numeric. Entropy value
#'

#' @noRd
.entropy_core <- function(proportions, q = 1, norm = FALSE, 
                                   log_base = exp(1), q_tol = 1e-6) {
  # Input validation
  if (!is.numeric(proportions) || length(proportions) == 0) {
    return(NA_real_)
  }
  
  if (any(is.na(proportions)) || any(is.infinite(proportions))) {
    return(NA_real_)
  }
  
  if (!is.numeric(q) || q < 0) {
    stop("q must be non-negative")
  }
  
  # Filter out zeros (standard in entropy)
  p_nonzero <- proportions[proportions > 1e-15]
  
  if (length(p_nonzero) == 0) {
    return(NA_real_)
  }
  
  # Normalize to sum to 1 (handle numerical errors)
  p <- p_nonzero / sum(p_nonzero)
  
  # Species richness (q=0): just count species
  if (q < q_tol) {
    H <- (log(length(p)) - 1) / log(log_base)
    return(H)
  }
  
  # Shannon entropy (q=1): use -sum(p*log(p))
  if (abs(q - 1) < q_tol) {
    H <- -sum(p * log(p)) / log(log_base)
  } else {
    # Tsallis entropy: (1 - sum(p^q)) / (q-1)
    H <- (1.0 - sum(p^q)) / ((q - 1.0) * log(log_base))
  }
  
  # Normalize by maximum entropy if requested
  if (norm) {
    n <- length(p)
    if (q < q_tol) {
      H_max <- (log(n) - 1) / log(log_base)
    } else if (abs(q - 1) < q_tol) {
      H_max <- log(n) / log(log_base)
    } else {
      H_max <- (1.0 - n^(1.0 - q)) / ((q - 1.0) * log(log_base))
    }
    
    if (!is.na(H_max) && !is.nan(H_max) && H_max > 0 && is.finite(H_max)) {
      H <- H / H_max
    }
  }
  
  return(H)
}

#' Internal: Vectorized Tsallis entropy calculation
#'
#' Compute entropy for multiple observations (rows = observations, cols = species)
#'
#' @param counts Matrix/data.frame where rows are observations, columns are species
#' @param q Numeric. Generalization parameter. Default: 1.0
#' @param norm Logical. Normalize by maximum entropy. Default: FALSE
#' @param log_base Numeric. Logarithm base. Default: exp(1)
#' @param pseudocount Numeric. Add to counts before normalization. Default: 0
#'
#' @return Numeric vector of entropy values (one per observation)
#'

#' @noRd
.entropy_vectorized <- function(counts, q = 1, norm = FALSE, 
                                        log_base = exp(1), pseudocount = 0) {
  counts <- as.matrix(counts)
  
  if (nrow(counts) == 0 || ncol(counts) == 0) {
    return(numeric(0))
  }
  
  # Apply to each row
  entropy_vals <- apply(counts, 1, function(row) {
    # Add pseudocount and normalize
    total <- sum(row, na.rm = TRUE) + length(row) * pseudocount
    if (total <= 0) return(NA_real_)
    
    p <- (row + pseudocount) / total
    .entropy_core(p, q = q, norm = norm, log_base = log_base)
  })
  
  return(unname(entropy_vals))
}

#' Internal: Maximum Tsallis entropy for n species
#'
#' Compute the theoretical maximum entropy for uniform distribution of n species
#'
#' @param n_species Integer. Number of species
#' @param q Numeric. Generalization parameter. Default: 1.0
#' @param log_base Numeric. Logarithm base. Default: exp(1)
#' @param q_tol Numeric. Tolerance for q=1 detection. Default: 1e-6
#'
#' @return Numeric. Maximum entropy value
#'

#' @noRd
.entropy_max <- function(n_species, q = 1, log_base = exp(1), q_tol = 1e-6) {
  if (n_species < 1) return(NA_real_)
  
  if (q < q_tol) {
    # Species richness max: log(n)
    H_max <- (log(n_species) - 1) / log(log_base)
  } else if (abs(q - 1) < q_tol) {
    # Shannon max: log(n)
    H_max <- log(n_species) / log(log_base)
  } else {
    # Tsallis max: (1 - n^(1-q)) / (q-1)
    H_max <- (1.0 - n_species^(1.0 - q)) / ((q - 1.0) * log(log_base))
  }
  
  return(H_max)
}
