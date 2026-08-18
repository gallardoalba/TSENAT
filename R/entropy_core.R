
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
#' @param q_tol Numeric. Tolerance for detecting q=1 case. Default: TSENAT_Q_TOL (1e-6)
#' @param n_present Integer or NULL. Raw support count (number of categories
#' with positive RAW counts, BEFORE pseudocount regularization). When not NULL,
#' the q=0 branch returns n_present - 1 instead of counting positive
#' proportions, so pseudocounts never alter the q=0 support statistic (audit3).
#'
#' @return Numeric. Entropy value
#'

#' @noRd
.entropy_core <- function(proportions, q = 1, norm = FALSE, log_base = exp(1), q_tol = TSENAT_Q_TOL,
    n_present = NULL) {
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

    # Negative proportions are invalid input: silently discarding them would
    # compute entropy on a different abundance vector than the user supplied
    # (audit 2026-08-17: reject instead of filter).
    if (any(proportions < 0)) {
        stop("Proportions must be non-negative.", call. = FALSE)
    }

    # Keep only positive proportions (zeros contribute 0 to entropy).
    # Removed >1e-15 threshold — zeros are valid, consistent with C++ fix #17
    p_nonzero <- proportions[proportions >= 0]

    if (length(p_nonzero) == 0) {
        return(NA_real_)
    }

    # Normalize to sum to 1 (handle numerical errors)
    p <- p_nonzero/sum(p_nonzero)

    # Species richness (q=0): S_0 = n_present - 1
    # Mathematically: S_0 = (1 - Σp_i^0)/(-1) = n_present - 1 (Tsallis 1988).
    # NOTE: This is the Tsallis ENTROPY (n_present-1), NOT the Hill number/effective
    # richness D_0 = n. The entropy value n_present-1 is consistent with the C++
    # implementation (entropy_cpp).
    # Changed from length(p) to length(p)-1.
    # Count only species with p > 0 (zeros contribute nothing
    # to species richness). Using sum(p > 0) instead of length(p) to exclude
    # zero-proportion isoforms.
    if (q < q_tol) {
        # AUDIT3 RED 2: q=0 uses RAW support (pre-pseudocount) when known, so
        # pseudocount regularization never inflates the support statistic to
        # the annotated universe.
        H <- (if (!is.null(n_present)) n_present else sum(p > 0)) - 1
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
        # Numerically stable evaluation near q = 1. The raw
        # form (1 - Σp^q)/(q-1) becomes 0/0 as q→1. With t = q-1 and
        # Σ p exp(t log p) = 1 + Σ p expm1(t log p), the entropy equals
        # -Σ p expm1(t log p) / t, which is stable for arbitrarily small t.
        p_pos <- p[p > 0]
        t_q <- q - 1
        if (abs(t_q) < 1e-07) {
            H <- -sum(p_pos * log(p_pos))/log(log_base)
        } else {
            H <- -sum(p_pos * expm1(t_q * log(p_pos)))/t_q
        }
    }

    # Normalize by maximum entropy if requested
    if (norm) {
        n <- length(p)
        if (q < q_tol) {
            # Tsallis q=0 max = n - 1, not log(n)
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

################################################################################
# SINGLE SOURCE OF TRUTH: q-tolerance for the entropy layer
#
# Previously .entropy_core() used q_tol = 1e-6 while .calc_S()/
# .calc_D() used sqrt(.Machine$double.eps) (~1.5e-8), so the same conceptual
# quantity could be evaluated differently depending on the internal route
# (e.g. .entropy_core(x, q = 1e-7) returned S_0 but
# .calculate_tsallis_entropy(x, q = 1e-7) returned S(1e-7)).
#
# TSENAT_Q_TOL is the ONE tolerance used by the entropy layer to decide
#   q ~ 0  -> Tsallis S_0 = n_present - 1  (support convention)
#   q ~ 1  -> Shannon limit
# The divergence layer keeps its own documented DIVERGENCE_Q_TOL = 1e-10
# (validated against the C++ kernels in test-divergence-kernel-validation.R).
################################################################################
TSENAT_Q_TOL <- 1e-06


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
        # Raw support BEFORE pseudocount (q=0 policy)
        support_raw <- sum(row > 0)
        # Add pseudocount and normalize
        total <- sum(row, na.rm = TRUE) + length(row) * pseudocount
        if (total <= 0)
            return(NA_real_)

        p <- (row + pseudocount)/total
        .entropy_core(p, q = q, norm = norm, log_base = log_base, n_present = support_raw)
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
#' @param q_tol Numeric. Tolerance for q=1 detection. Default: TSENAT_Q_TOL (1e-6)
#'
#' @return Numeric. Maximum entropy value
#'

#' @noRd
.entropy_max <- function(n_species, q = 1, log_base = exp(1), q_tol = TSENAT_Q_TOL) {
    if (n_species < 1)
        return(NA_real_)

    if (q < q_tol) {
        # Tsallis q=0 max = n_species - 1 (not n_species).
        # The maximum Tsallis entropy at q=0 for n species is n-1, achieved
        # when all n species have equal (non-zero) proportions.  Using n_species
        # would cause inconsistent normalization with .entropy_core() which
        # correctly uses length(p)-1.
        H_max <- n_species - 1
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
#' @param q_tol Numeric. Tolerance for detecting q=1 case. Default: TSENAT_Q_TOL (1e-6)
#'
#' @return Numeric scalar: entropy value
#'
#' @references
#'   Originally from jackknife_diagnostics.R, consolidated into entropy_core.R
#'   for unified entropy calculations across the package.
#'
#' @noRd
.entropy_single <- function(counts, q = 1, norm = TRUE, log_base = exp(1), pseudocount = 0,
    q_tol = TSENAT_Q_TOL) {
    # Pseudocount contract: .entropy_single adds pseudocount and normalizes to
    # proportions before calling .entropy_core. .entropy_core expects already-
    # proportioned input with no pseudocount handling. .entropy_vectorized also
    # adds pseudocount at the row level. All three functions form a stack where
    # pseudocount is applied at the entry points (.entropy_single, .entropy_vectorized)
    # but not in the core computation (.entropy_core).
    # Normalize to proportions with pseudocount
    total <- sum(counts) + length(counts) * pseudocount
    if (total <= 0)
        return(NA_real_)

    # Raw support BEFORE pseudocount (q=0 policy, audit3)
    support_raw <- sum(counts > 0)

    p <- (counts + pseudocount)/total
    n <- length(p)

    # Calculate entropy using standardized core logic
    if (q < q_tol) {
        # Tsallis entropy at q=0: S_0 = n_nonzero - 1 on RAW support, so a
        # positive pseudocount never inflates richness to the annotated
        # universe.
        entropy <- support_raw - 1
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
        # Stable evaluation near q = 1 via
        # H = -Σ p expm1(t log p) / t with t = q - 1.
        p_nonzero <- p[p > 0]
        t_q <- q - 1
        if (abs(t_q) < 1e-07) {
            entropy <- if (length(p_nonzero) > 0) {
                -sum(p_nonzero * log(p_nonzero)/log(log_base))
            } else 0
        } else {
            entropy <- -sum(p_nonzero * expm1(t_q * log(p_nonzero)))/t_q
        }
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
