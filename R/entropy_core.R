################################################################################
#' Internal: Core Tsallis entropy calculation (consolidated)
#' 
#' Centralized entropy calculation used by all functions to eliminate
#' duplication.
#' Supports Shannon (q=1), Tsallis (q!=1), and species richness (q=0).
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
.entropy_core <- function(proportions, q = 1, norm = FALSE, log_base = exp(1), q_tol = 1e-06) {
    # Input validation
    if (!is.numeric(proportions) || length(proportions) == 0) {
        return(NA_real_)
    }

    if (any(is.na(proportions)) || any(is.infinite(proportions))) {
        return(NA_real_)
    }

    if (!is.numeric(q)) {
        stop("q must be numeric")
    }
    if (length(q) > 1) {
        stop("q must be a scalar (length 1), not a vector")
    }
    if (q < 0) {
        stop("q must be non-negative")
    }

        # Filter out only negative values (zeros contribute 0 to entropy)
    # Filter out only negative values (zeros contribute 0 to entropy)
    # AUDIT FIX R15: Removed >1e-15 threshold — zeros are valid, consistent with C++ fix #17
    p_nonzero <- proportions[proportions >= 0]

    if (length(p_nonzero) == 0) {
        return(NA_real_)
    }

    # Normalize to sum to 1 (handle numerical errors)
    p <- p_nonzero/sum(p_nonzero)

    # Species richness (q=0): S_0 = n_nonzero - 1
    # Mathematically: S_0 = (1 - Σp_i^0)/(-1) = n - 1 (Tsallis 1988).
    # NOTE: This is the Tsallis ENTROPY (n-1), NOT the Hill number/effective
    # richness D_0 = n. The entropy value n-1 is consistent with the C++
    # implementation (entropy_cpp, audit fix #5).
    # AUDIT FIX R13: Changed from length(p) to length(p)-1.
    if (q < q_tol) {
        H <- length(p) - 1
        return(H)
    }

    # Shannon entropy (q=1): -sum(p*log(p)) for p > 0 only
    if (abs(q - 1) < q_tol) {
        p_pos <- p[p > 0]
        H <- -sum(p_pos * log(p_pos))/log(log_base)
    } else {
        # Tsallis entropy: (1 - sum(p^q)) / (q-1)
        # NOTE: log_base is NOT applied to Tsallis entropy, consistent with
        # .entropy_single() and the standard mathematical definition.
        # The Tsallis formula (1 - Σp^q)/(q-1) is scale-invariant and does
        # not use a logarithm base. Only Shannon entropy (q→1 limit) uses
        # log_base for cross-study comparability.
        H <- (1 - sum(p^q))/(q - 1)
    }

    # Normalize by maximum entropy if requested
    if (norm) {
        n <- length(p)
        if (q < q_tol) {
            # AUDIT FIX R14: Tsallis q=0 max = n - 1, not log(n)
            H_max <- length(p) - 1
        } else if (abs(q - 1) < q_tol) {
            H_max <- log(n)/log(log_base)
        } else {
            # Tsallis: max entropy (no log_base applied, consistent with
            # .entropy_single() and standard Tsallis definition)
            H_max <- (1 - n^(1 - q))/(q - 1)
        }

        if (!is.na(H_max) && !is.nan(H_max) && H_max > 0 && is.finite(H_max)) {
            H <- H/H_max
        }
    }

    return(H)
}

#' Internal: Vectorized Tsallis entropy calculation
#'
#' Compute entropy for multiple observations (rows = observations, cols =
#' species)
#'
#' @param counts Matrix/data.frame where rows are observations, columns are
#' species
#' @param q Numeric. Generalization parameter. Default: 1.0
#' @param norm Logical. Normalize by maximum entropy. Default: FALSE
#' @param log_base Numeric. Logarithm base. Default: exp(1)
#' @param pseudocount Numeric. Add to counts before normalization. Default: 0
#'
#' @return Numeric vector of entropy values (one per observation)
#'

#' @noRd
.entropy_vectorized <- function(counts, q = 1, norm = FALSE, log_base = exp(1), pseudocount = 0) {
    counts <- as.matrix(counts)

    if (nrow(counts) == 0 || ncol(counts) == 0) {
        return(numeric(0))
    }

    # Apply to each row
    entropy_vals <- apply(counts, 1, function(row) {
        # Add pseudocount and normalize
        total <- sum(row, na.rm = TRUE) + length(row) * pseudocount
        if (total <= 0)
            return(NA_real_)

        p <- (row + pseudocount)/total
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
.entropy_max <- function(n_species, q = 1, log_base = exp(1), q_tol = 1e-06) {
    if (n_species < 1)
        return(NA_real_)

    if (q < q_tol) {
        # Species richness max: total number of species
        H_max <- n_species
    } else if (abs(q - 1) < q_tol) {
        # Shannon max: log(n)
        H_max <- log(n_species)/log(log_base)
    } else {
        # Tsallis max: (1 - n^(1-q)) / (q-1) [log_base NOT applied to Tsallis]
        H_max <- (1 - n_species^(1 - q))/(q - 1)
    }

    return(H_max)
}

#' Internal: Entropy calculation for a single vector with pseudocount support
#'
#' Compute Tsallis, Shannon, or species richness entropy for a single counts
#' vector,
#' with optional pseudocount and log base parameters. Used by jackknife and
#' other
#' internal calculations.
#'
#' @param counts Numeric vector of (non-negative) counts
#' @param q Numeric. Generalization parameter. Default: 1.0 (Shannon entropy)
#' @param norm Logical. Normalize by maximum entropy [0,1]. Default: TRUE
#' @param log_base Numeric. Logarithm base. Default: exp(1) (natural log)
#' @param pseudocount Numeric. Add to each count before normalization.
#' Default: 0
#' @param q_tol Numeric. Tolerance for detecting q=1 case. Default: 1e-6
#'
#' @return Numeric scalar: entropy value
#'
#' @references
#'   Originally from jackknife_diagnostics.R, consolidated into entropy_core.R
#'   for unified entropy calculations across the package.
#'
#' @noRd
.entropy_single <- function(counts, q = 1, norm = TRUE, log_base = exp(1), pseudocount = 0,
    q_tol = 1e-06) {
    # Normalize to proportions with pseudocount
    total <- sum(counts) + length(counts) * pseudocount
    if (total <= 0)
        return(NA_real_)

    p <- (counts + pseudocount)/total
    n <- length(p)

    # Calculate entropy using standardized core logic
    if (q < q_tol) {
        # Tsallis entropy at q=0: S_0 = n_nonzero - 1
        # Consistent with .entropy_core() (AUDIT FIX R13) and entropy_cpp
        # (AUDIT FIX #5).  Not the Hill number / effective richness D_0 = n.
        p_nonzero <- p[p > 0]
        entropy <- length(p_nonzero) - 1
    } else if (abs(q - 1) < q_tol) {
        # Shannon entropy as q -> 1
        p_nonzero <- p[p > 0]
        if (length(p_nonzero) > 0) {
            entropy <- -sum(p_nonzero * log(p_nonzero)/log(log_base))
        } else {
            entropy <- 0
        }
    } else {
        # Generalized Tsallis entropy: (1 - sum(p^q)) / (q-1) [log_base NOT
        # applied]
        entropy <- (1/(q - 1)) * (1 - sum(p^q))
    }

    # Normalize to [0, 1] if requested
    if (norm) {
        # Maximum entropy achieved with uniform distribution
        if (abs(q - 1) < q_tol) {
            max_entropy <- log(n)/log(log_base)
        } else {
            # Tsallis: no log_base applied to maximum
            max_entropy <- (1/(q - 1)) * (1 - n^(1 - q))
        }

        if (!is.na(max_entropy) && !is.nan(max_entropy) && max_entropy > 0 && is.finite(max_entropy)) {
            entropy <- entropy/max_entropy
        }
    }

    return(as.numeric(entropy))
}
