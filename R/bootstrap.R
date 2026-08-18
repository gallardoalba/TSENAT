# ============================================================================
# SHARED BOOTSTRAP INPUT VALIDATION (July 2026: metrics.json refactoring)
# ============================================================================

#' Validate bootstrap input data
#'
#' Centralized validation for all bootstrap C++ wrappers. Eliminates duplicated
#' validation logic across `bootstrap.R` and `bootstrap_divergence.R`.
#' Handles common checks: zero-length, NA/Inf, negative values, all-zero,
#' paired even-length, and vector pseudocount validation.
#'
#' @param x Numeric vector. Primary data for validation.
#' @param context Character. Label for error messages (e.g., "paired bootstrap").
#' @param paired Logical. If TRUE, validates even-length requirement for paired designs.
#' @param pseudocount Numeric vector or scalar. If vector, validated against x length.
#' @param allow_empty Logical. If TRUE, empty vectors pass validation (for edge cases).
#'
#' @return Invisibly returns the (possibly adjusted) x. Throws on validation failure.
#'
#' @noRd
.validate_bootstrap_input <- function(x, context = "bootstrap", paired = FALSE,
    pseudocount = 0, allow_empty = FALSE) {
    # NULL guard: treat as empty (as.numeric(NULL) → numeric(0))
    if (is.null(x)) {
        if (allow_empty) return(invisible(numeric(0)))
        stop("For ", context, ", input vector cannot be empty")
    }

    # Type check: reject non-numeric before silent coercion
    if (!is.numeric(x)) {
        stop("For ", context, ", input must be numeric")
    }

    # Defensive copy to prevent accidental modification of caller's data
    x <- as.numeric(x)

    # Zero-length check (numeric(0) from empty numeric input)
    if (length(x) == 0) {
        if (allow_empty) return(invisible(x))
        stop("For ", context, ", input vector cannot be empty")
    }

    # NA/Inf check
    na_count <- sum(is.na(x))
    inf_count <- sum(is.infinite(x))
    if (na_count > 0) {
        warning("Input vector contains ", na_count, " NA values. ",
                "These will affect bootstrap resampling. ",
                "Consider removing NA values before calling ", context, ".")
    }
    if (inf_count > 0) {
        stop("For ", context, ", input vector contains ", inf_count,
             " infinite values. Cannot compute meaningful estimates.")
    }

    # Negative values check
    neg_count <- sum(x < 0, na.rm = TRUE)
    if (neg_count > 0) {
        stop("For ", context, ", input contains ", neg_count,
             " negative values. Count data must be non-negative.")
    }

    # All-zero check
    if (all(x == 0, na.rm = TRUE)) {
        stop("All values in ", context, " data are zero. ",
             "Cannot compute meaningful entropy/divergence estimates.")
    }

    # Paired: even-length requirement
    if (paired && length(x) %% 2 != 0) {
        stop("For ", context, " with paired=TRUE, input vector must have ",
             "even length (pairs). Got length=", length(x),
             ". Please verify pairing structure.")
    }

    # Vector pseudocount validation
    if (length(pseudocount) > 1 && length(pseudocount) != length(x)) {
        stop("For ", context, ", pseudocount must have length 1 or equal to x length")
    }

    invisible(x)
}

# ============================================================================
# C++ BOOTSTRAP WRAPPERS (Optimized resampling)
# ============================================================================

#' Block Bootstrap Entropy Computation
#'
#' @description
#' C++ wrapper for block bootstrap computation on paired resampling data.
#' Performs multiple bootstrap iterations for entropy estimation.
#'
#' @param x \code{numeric}.  Data vector (must have even length for 
#' paired design).
#' @param q \code{numeric}. Tsallis q parameter. Default: 1.0.
#' @param normalize \code{logical}. Normalize entropy? Default: TRUE.
#' @param nboot \code{integer}. Number of bootstrap samples. Default: 1000.
#' @param log_base \code{numeric}. Logarithm base. Default: e (natural log).
#' @param pseudocount \code{numeric}.  Pseudocount for  abundance inflation.
#'  Default:  0.
#'
#' @return \code{numeric}. Vector of nboot bootstrap entropy estimates.
#'
#' @details
#' Performs block bootstrap for paired samples with C++ acceleration.
#' Input must have even length (pairs). Accelerated for speed.
#'
#' @noRd
block_bootstrap_compute_cpp_wrapper <- function(x, q = 1, normalize = TRUE, nboot = 1000L,
    log_base = exp(1), pseudocount = 0) {
    # Centralized validation via shared helper (July 2026 refactoring)
    .validate_bootstrap_input(x, context = "paired block bootstrap", paired = TRUE,
        pseudocount = pseudocount)
    x <- as.numeric(x)

    # Handle vector pseudocount
    if (length(pseudocount) > 1) {
        if (length(pseudocount) != length(x)) {
            stop("pseudocount must have length 1 or equal to x length")
        }
        # For bootstrap, apply vector pseudocount once upfront
        x_adj <- x + pseudocount
        pseudocount_scalar <- 0  # Already applied above
    } else {
        x_adj <- x
        pseudocount_scalar <- pseudocount
    }

    .Call("_TSENAT_block_bootstrap_compute_cpp", PACKAGE = "TSENAT", as.numeric(x_adj),
        as.integer(nboot), as.numeric(q), as.logical(normalize), as.numeric(log_base),
        as.numeric(pseudocount_scalar))
}

# ============================================================================
# C++ BOOTSTRAP WRAPPERS (Optimized resampling)
# ============================================================================

#' Standard Bootstrap Entropy Computation
#'
#' @description
#' C++ wrapper for standard bootstrap computation with independent resampling.
#' Performs multiple bootstrap iterations for entropy estimation.
#'
#' @param x \code{numeric}. Data vector.
#' @param q \code{numeric}. Tsallis q parameter. Default: 1.0.
#' @param normalize \code{logical}. Normalize entropy? Default: TRUE.
#' @param nboot \code{integer}. Number of bootstrap samples. Default: 1000.
#' @param log_base \code{numeric}. Logarithm base. Default: e (natural log).
#' @param pseudocount \code{numeric}.  Pseudocount for  abundance inflation.
#'  Default:  0.
#'
#' @return \code{numeric}. Vector of nboot bootstrap entropy estimates.
#'
#' @details
#' Performs standard (independent) bootstrap with C++ acceleration.
#' Handles vector pseudocounts by applying them upfront.
#'
#' @noRd
bootstrap_compute_cpp_wrapper <- function(x, q = 1, normalize = TRUE, nboot = 1000L,
    log_base = exp(1), pseudocount = 0) {
    # Centralized validation via shared helper (July 2026 refactoring)
    .validate_bootstrap_input(x, context = "standard bootstrap",
        pseudocount = pseudocount)
    x <- as.numeric(x)

    # Handle vector pseudocount by converting to scalar (sum per-element
    # effects)
    if (length(pseudocount) > 1) {
        if (length(pseudocount) != length(x)) {
            stop("pseudocount must have length 1 or equal to x length")
        }
        # For bootstrap, apply vector pseudocount once upfront
        x_adj <- x + pseudocount
        pseudocount_scalar <- 0  # Already applied above
    } else {
        x_adj <- x
        pseudocount_scalar <- pseudocount
    }

    .Call("_TSENAT_bootstrap_compute_cpp", PACKAGE = "TSENAT", as.numeric(x_adj),
        as.integer(nboot), as.numeric(q), as.logical(normalize), as.numeric(log_base),
        as.numeric(pseudocount_scalar))
}

#' Multi-q standard bootstrap (one resample -> all q)
#'
#' Returns an `nboot x length(q)` matrix. The resampling plan is drawn ONCE
#' per iteration and every q is evaluated on the same resample, removing the
#' O(Q*B) resampling cost and preserving the joint correlation across q.
#'
#' @noRd
bootstrap_compute_multi_q_cpp_wrapper <- function(x, q, normalize = TRUE, nboot = 1000L,
    log_base = exp(1), pseudocount = 0) {
    .validate_bootstrap_input(x, context = "standard multi-q bootstrap",
        pseudocount = pseudocount)
    x <- as.numeric(x)

    if (length(pseudocount) > 1) {
        if (length(pseudocount) != length(x)) {
            stop("pseudocount must have length 1 or equal to x length")
        }
        x_adj <- x + pseudocount
        pseudocount_scalar <- 0
    } else {
        x_adj <- x
        pseudocount_scalar <- pseudocount
    }

    .Call("_TSENAT_bootstrap_compute_multi_q_cpp", PACKAGE = "TSENAT", as.numeric(x_adj),
        as.numeric(q), as.integer(nboot), as.logical(normalize), as.numeric(log_base),
        as.numeric(pseudocount_scalar))
}

# ============================================================================
# REPLICATE-LEVEL BOOTSTRAP WRAPPER
# ============================================================================

#' Replicate-Level Bootstrap Entropy Computation
#'
#' @description
#' C++ wrapper for replicate-level bootstrap: resamples entire samples (columns)
#' with replacement instead of individual reads. Captures biological variability.
#'
#' @param counts \code{numeric matrix}. Counts matrix (transcripts × samples).
#' @param q \code{numeric}. Tsallis q parameter. Default: 1.0.
#' @param normalize \code{logical}. Normalize entropy? Default: TRUE.
#' @param nboot \code{integer}. Number of bootstrap samples. Default: 1000.
#' @param log_base \code{numeric}. Logarithm base. Default: e (natural log).
#' @param pseudocount \code{numeric}. Pseudocount. Default: 0.0.
#' @param block_ids \code{integer}. Block IDs for block-aware resampling. Default: none.
#'
#' @return \code{numeric}. Vector of nboot bootstrap entropy estimates.
#'
#' @noRd
bootstrap_replicate_cpp_wrapper <- function(counts, q = 1, normalize = TRUE, nboot = 1000L,
    log_base = exp(1), pseudocount = 0.0, block_ids = integer(0)) {
    if (!is.matrix(counts)) {
        stop("counts must be a matrix (transcripts x samples) for replicate bootstrap")
    }
    if (nrow(counts) < 2 || ncol(counts) < 2) {
        stop("Need at least 2 transcripts and 2 samples for replicate bootstrap. ",
             "Got ", nrow(counts), " x ", ncol(counts), ".")
    }
    bootstrap_replicate_cpp(counts = counts, nboot = nboot, q = q,
        normalize = normalize, log_base = log_base, pseudocount = pseudocount,
        block_ids = block_ids)
}

#' Bootstrap Divergence Computation
#'
#' @description
#' C++ wrapper for bootstrap resampling of Tsallis divergence between two
#' distributions. Supports independent and paired (block) bootstrap modes.
#'
#' @param x \code{numeric}. Count vector for first distribution.
#' @param y \code{numeric}. Count vector for second distribution.
#' @param q \code{numeric}. Tsallis q parameter. Default: 1.0 (KL divergence).
#' @param nboot \code{integer}. Number of bootstrap replicates. Default: 1000L.
#' @param paired \code{logical}. Use paired (block) bootstrap? Default: FALSE.
#' @param pseudocount \code{numeric}. Pseudocount to add. Default: 0.0.
#' @param log_base \code{numeric}. Logarithm base. Default: e (natural log).
#'
#' @return \code{numeric}. Vector of nboot bootstrap divergence estimates.
#'
#' @details
#' Performs bootstrap divergence computation with C++ acceleration.
#' For paired=TRUE, resamples pairs as units maintaining correlation structure.
#' For paired=FALSE (default), resamples x and y independently.
#'
#' @noRd
divergence_bootstrap_compute_cpp_wrapper <- function(x, y, q = 1, nboot = 1000L,
    paired = FALSE, pseudocount = 0, log_base = exp(1)) {
    # Divergence-specific validation (must precede shared validation)
    if (!is.logical(paired) || length(paired) != 1) {
        stop("paired must be a single logical value")
    }
    if (!is.numeric(q) || q < 0) {
        stop("q must be a non-negative numeric value")
    }
    if (!is.numeric(nboot) || nboot < 1 || nboot != as.integer(nboot)) {
        stop("nboot must be a positive integer")
    }

    # Centralized validation via shared helper (July 2026 refactoring)
    .validate_bootstrap_input(x, context = "divergence bootstrap (x)", paired = paired,
        pseudocount = pseudocount)
    .validate_bootstrap_input(y, context = "divergence bootstrap (y)", paired = paired,
        pseudocount = pseudocount)

    # Cross-validation between x and y
    if (length(x) != length(y)) {
        stop("x and y must have the same length. Got length(x)=", length(x),
             ", length(y)=", length(y), ". ",
             "For mixed paired/unpaired designs, use divergence_bootstrap_flexible_cpp_wrapper instead.")
    }

    # Handle vector pseudocount
    if (length(pseudocount) > 1) {
        if (length(pseudocount) != length(x)) {
            stop("pseudocount must have length 1 or equal to x length")
        }
        # Apply vector pseudocount upfront
        x_adj <- x + pseudocount
        y_adj <- y + pseudocount
        pseudocount_scalar <- 0
    } else {
        x_adj <- x
        y_adj <- y
        pseudocount_scalar <- pseudocount
    }

    # Call C++ function
    .Call("_TSENAT_divergence_bootstrap_compute_cpp", PACKAGE = "TSENAT", as.numeric(x_adj),
        as.numeric(y_adj), as.integer(nboot), as.numeric(q), as.logical(paired),
        as.numeric(pseudocount_scalar), as.numeric(log_base))
}

#' Paired Divergence Bootstrap (C++ Optimized)
#'
#' @description
#' C++ wrapper for paired divergence bootstrap with explicit pair structure
#' handling.
#' Used for paired/matched study designs where samples are linked across groups.
#'
#' @param x numeric. Control group counts.
#' @param y numeric. Treatment group counts.
#' @param pair_ids integer. Pair identifiers matching(length = length(x)).
#' @param nboot integer. Number of bootstrap replicates. Default: 1000.
#' @param q numeric. Tsallis q parameter. Default: 1.0.
#' @param pseudocount numeric. Pseudocount adjustment. Default: 0.0.
#' @param log_base numeric. Logarithm base. Default: e (natural log).
#'
#' @return numeric. Vector of nboot bootstrap divergence estimates.
#'
#' @details
#' Performs pair-respecting bootstrap by:
#' 1. Extracting pair structure from pair_ids
#' 2. Resampling pairs (not individual samples)
#' 3. Aggregating counts per resampled pair
#' 4. Computing divergence between resampled distributions
#'
#' @noRd
divergence_bootstrap_paired_cpp_wrapper <- function(x, y, pair_ids, nboot = 1000L,
    q = 1, pseudocount = 0, log_base = exp(1)) {

    # Input validation
    if (length(x) != length(y)) {
        stop("x and y must have equal length")
    }
    if (length(x) != length(pair_ids)) {
        stop("pair_ids must have same length as x and y")
    }

    # Paired designs require even-length vectors for proper pairing
    if (length(x) %% 2 != 0) {
        stop("x and y must have same length, with even number of elements for paired bootstrap")
    }

    # Handle vector pseudocount
    if (length(pseudocount) > 1) {
        if (length(pseudocount) != length(x)) {
            stop("pseudocount must have length 1 or equal to x length")
        }
        # Apply vector pseudocount upfront
        x_adj <- x + pseudocount
        y_adj <- y + pseudocount
        pseudocount_scalar <- 0
    } else {
        x_adj <- x
        y_adj <- y
        pseudocount_scalar <- pseudocount
    }

    # Ensure pair_ids is integer
    pair_ids_int <- as.integer(pair_ids)

    # Call C++ function
    .Call("_TSENAT_divergence_bootstrap_paired_cpp", PACKAGE = "TSENAT", as.numeric(x_adj),
        as.numeric(y_adj), pair_ids_int, as.integer(nboot), as.numeric(q), as.numeric(pseudocount_scalar),
        as.numeric(log_base))
}

# ============================================================================
# ENHANCED: Flexible Paired/Unpaired Bootstrap (NEW - MARCH 2026)
# ============================================================================
#' C++ Flexible Paired/Unpaired Divergence Bootstrap
#'
#' @description
#' Enhanced C++ implementation supporting complete pairs, incomplete pairs, 
#' and unpaired samples. Handles arbitrary mixing of paired and unpaired data
#' with efficient correlated resampling for pairs.
#'
#' @param x Numeric vector of control group counts
#' @param y Numeric vector of treatment group counts  
#' @param x_pair_ids Integer vector of pair IDs for x (0/NA = unpaired)
#' @param y_pair_ids Integer vector of pair IDs for y (0/NA = unpaired)
#' @param nboot Number of bootstrap iterations (default: 1000)
#' @param q Tsallis entropy order (default: 1.0 for Shannon entropy)
#' @param pseudocount Pseudocount to add (scalar or vector, default: 0.0)
#' @param log_base Logarithm base (default: exp(1) for natural log)
#'
#' @details
#' This function extends paired bootstrap to handle:
#' - Complete pairs: both pair_ids present in x AND y (resampled as units)
#' - Unpaired x: pair_id present in x but not y (resampled independently)
#' - Unpaired y: pair_id present in y but not x (resampled independently)
#' 
#' Allows different numbers of samples per group and arbitrary pairing patterns.
#' ~10x speedup compared to R implementation.
#'
#' @return Numeric vector of Tsallis divergence bootstrap estimates
#'
#' @noRd
divergence_bootstrap_flexible_cpp_wrapper <- function(x, y, x_pair_ids, y_pair_ids,
    nboot = 1000L, q = 1, pseudocount = 0, log_base = exp(1)) {

    # Input validation
    if (length(x) != length(x_pair_ids)) {
        stop("x and x_pair_ids must have same length")
    }
    if (length(y) != length(y_pair_ids)) {
        stop("y and y_pair_ids must have same length")
    }

    # Handle vector pseudocount for x and y (combined length)
    if (length(pseudocount) > 1) {
        expected_len <- length(x) + length(y)
        if (length(pseudocount) != expected_len) {
            stop("pseudocount must have length 1 or ", expected_len,
                 " (combined x+y length), got ", length(pseudocount))
        }
        # Split pseudocount using safe indexing (seq_along for efficiency)
        x_pseudo <- pseudocount[seq_along(x)]
        y_pseudo <- pseudocount[length(x) + seq_along(y)]
        x_adj <- x + x_pseudo
        y_adj <- y + y_pseudo
        pseudocount_scalar <- 0
    } else {
        x_adj <- x
        y_adj <- y
        pseudocount_scalar <- pseudocount
    }

    # Convert pair_ids to integer for C++ with explicit type conversion
    # Avoid silent truncation/corruption by using as.integer
    x_pair_ids_int <- as.integer(x_pair_ids)
    y_pair_ids_int <- as.integer(y_pair_ids)
    
    # Validate conversion succeeded (check for NAs from coercion)
    if (any(is.na(x_pair_ids_int) & !is.na(x_pair_ids))) {
        stop("Failed to convert x_pair_ids to integer. ",
             "Check that pair_ids contain only whole numbers.")
    }
    if (any(is.na(y_pair_ids_int) & !is.na(y_pair_ids))) {
        stop("Failed to convert y_pair_ids to integer. ",
             "Check that pair_ids contain only whole numbers.")
    }

    # Call C++ function
    .Call("_TSENAT_divergence_bootstrap_flexible_cpp", PACKAGE = "TSENAT", as.numeric(x_adj),
        as.numeric(y_adj), x_pair_ids_int, y_pair_ids_int, as.integer(nboot), as.numeric(q),
        as.numeric(pseudocount_scalar), as.numeric(log_base))
}

# ============================================================================
# BOOTSTRAP RESAMPLING INPUT PREPARATION (auditx follow-up, 2026-08)
# ============================================================================
# Prepares the resampling input so that the bootstrap's multinomial
# probabilities equal the point-estimate proportions EXACTLY.
#
# Point estimator: T(x, l, c) = (x/l + c) / sum(x/l + c)
# (effective-length division FIRST, pseudocount on the effective-abundance
# scale, then normalization to proportions).
#
# The read-level bootstrap draws a multinomial whose total must equal the
# original depth sum(x). Embedding the pseudocount BEFORE the depth rescale
#   x_eff = (x/l + c) * sum(x) / sum(x/l + c)
# preserves both properties: the probabilities equal T(x, l, c) exactly (the
# scale factor cancels in the ratio) and the total equals sum(x). The returned
# pseudocount is therefore 0 -- it is already embedded in x_eff.
#
# The previous implementation rescaled the normalized abundances FIRST and
# only then added the pseudocount: ((x/l)*k + c) does not factor with k, so
# the resampling probabilities disagreed with the point estimate whenever
# c > 0 and the effective lengths vary within the unit.
#
# Without effective lengths nothing is transformed: the kernels' own
# pseudocount addition already yields T(x, NULL, c) exactly.
#
# Returns list(x, pseudocount, scale_k); scale_k is the depth-rescaling factor
# (1 when no rescaling occurs).
.prepare_bootstrap_resample <- function(x, effective_length = NULL, pseudocount = 0) {
    scale_k <- 1
    x_eff <- x
    pseudocount_eff <- pseudocount
    if (!is.null(effective_length) && length(effective_length) == length(x)) {
        x_abund <- x/effective_length
        # Invalid effective lengths must FAIL, not silently become zeros
        # (audit 2026-08-17: a zero length previously turned C/0 -> Inf -> 0,
        # silently dropping the transcript from the resampling input).
        if (any(!is.finite(x_abund))) {
            stop("[.prepare_bootstrap_resample] effective_length normalization produced non-finite abundances. effective_length must contain finite positive values.",
                call. = FALSE)
        }
        x_adj <- x_abund + pseudocount
        sum_original <- sum(x)
        sum_adj <- sum(x_adj)
        if (sum_adj > 0) {
            scale_k <- sum_original/sum_adj
            x_eff <- x_adj * scale_k
            pseudocount_eff <- 0  # already embedded in x_eff
        } else {
            # Degenerate (all-zero after adjustment): keep the normalized
            # abundances; downstream guards handle the zero total.
            x_eff <- x_abund
            scale_k <- 1
        }
    }
    list(x = x_eff, pseudocount = pseudocount_eff, scale_k = scale_k)
}

# ============================================================================
# OPTIMIZED BOOTSTRAP RESAMPLE (C++ accelerated when available)
# ============================================================================

#' C++ Accelerated Bootstrap Resampling
#'
#' @description
#' Optimized bootstrap resampling using C++ via Rcpp. Handles both independent
#' and paired (block) bootstrap with entropy computation.
#'
#' @noRd
.bootstrap_resample_optimized <- function(x, q, norm, nboot, log_base, pseudocount,
    what, paired = FALSE, effective_length = NULL, resample_by = c("read", "replicate"),
    counts_matrix = NULL) {
    resample_by <- match.arg(resample_by)
    # ==================================================================
    # BOOTSTRAP RESAMPLING TRANSFORMATION (auditx follow-up, 2026-08)
    # ==================================================================
    # INVARIANT: the bootstrap must resample from EXACTLY the point-estimate
    # proportions T(x, l, c) = (x/l + c) / sum(x/l + c).
    #
    # .prepare_bootstrap_resample() embeds the pseudocount on the
    # effective-abundance scale and rescales to the original depth, so the
    # multinomial probabilities equal T(x, l, c) exactly and the total equals
    # sum(x). The previous implementation rescaled the normalized abundances
    # but added the pseudocount AFTER the rescale, which broke the
    # factorization when c > 0 (the pseudocount no longer cancelled with the
    # scale factor). With c = 0 the two formulations are identical.
    if (!is.null(effective_length) && length(effective_length) != length(x)) {
        warning("effective_length provided but length mismatch: length(effective_length)=",
            if (!is.null(effective_length))
                length(effective_length) else "NULL", " vs length(x)=", length(x), call. = FALSE)
    }
    pseudocount_original <- pseudocount
    prep <- .prepare_bootstrap_resample(x, effective_length, pseudocount)
    x_for_bootstrap <- prep$x
    pseudocount <- prep$pseudocount

    # Dispatch to C++ block bootstrap for paired samples
    if (paired) {
        if (length(x_for_bootstrap)%%2 != 0) {
            stop("For paired=TRUE, data must have even length (n_pairs * 2)")
        }

        # Block bootstrap for paired samples
        if (what == "S") {
            # For entropy
            bootstrap_dist <- block_bootstrap_compute_cpp_wrapper(x = x_for_bootstrap,
                q = q, normalize = norm, nboot = nboot, log_base = log_base, pseudocount = pseudocount)
        } else if (what == "D") {
            # For Hill numbers: compute entropy then convert
            bootstrap_dist <- block_bootstrap_compute_cpp_wrapper(x = x_for_bootstrap,
                q = q, normalize = FALSE, nboot = nboot, log_base = log_base, pseudocount = pseudocount)
            # Hill number conversion: D_q = (1 - (q-1) * H_q)^(1/(1-q))
            if (abs(q - 1) < 1e-06) {
                bootstrap_dist <- exp(bootstrap_dist)  # exp(H) for q=1
            } else {
                bootstrap_dist <- (1 - (q - 1) * bootstrap_dist)^(1/(1 - q))
            }
        } else {
            stop("Invalid 'what' parameter: must be 'S' (entropy) or 'D' (Hill numbers)")
        }

        return(bootstrap_dist)
    }

    # Replicate-level bootstrap — resample entire samples
    # (columns) with replacement to capture biological variability between
    # replicates. This is in contrast to the default "read" mode which
    # resamples individual reads (multinomial) and only captures sampling noise.
    if (resample_by == "replicate") {
        if (!is.null(counts_matrix) && is.matrix(counts_matrix)) {
            # Replicate-level resampling applies the point-estimator
            # transformation T(x, l, c) = (x/l + c)/sum(x/l + c) to the
            # resampled samples. The vector preparation above does not touch
            # counts_matrix, so apply the effective-length division here and
            # embed the pseudocount on the abundance scale (the C++ kernel
            # then adds 0).
            if (!is.null(effective_length) && length(effective_length) == nrow(counts_matrix)) {
                counts_matrix <- sweep(counts_matrix, 1, effective_length, "/") + pseudocount_original
                pseudocount_matrix <- 0
            } else {
                pseudocount_matrix <- pseudocount_original
            }
            # Use C++ accelerated replicate bootstrap on transcript × sample matrix
            if (what == "S") {
                bootstrap_dist <- bootstrap_replicate_cpp_wrapper(
                    counts = counts_matrix, q = q, normalize = norm,
                    nboot = nboot, log_base = log_base, pseudocount = pseudocount_matrix)
            } else if (what == "D") {
                bootstrap_dist <- bootstrap_replicate_cpp_wrapper(
                    counts = counts_matrix, q = q, normalize = FALSE,
                    nboot = nboot, log_base = log_base, pseudocount = pseudocount_matrix)
                if (abs(q - 1) < 1e-06) {
                    bootstrap_dist <- exp(bootstrap_dist)
                } else {
                    base <- 1 - (q - 1) * bootstrap_dist
                    base <- pmax(base, 1e-10)
                    bootstrap_dist <- base^(1/(1 - q))
                }
            } else {
                stop("Invalid 'what' parameter: must be 'S' or 'D'")
            }
            return(bootstrap_dist)
        } else {
            # R-level replicate resampling for vector input:
            # Resample indices with replacement and aggregate
            n_obs <- length(x_for_bootstrap)
            if (n_obs < 2) {
                stop("Need at least 2 observations for replicate bootstrap")
            }
            bootstrap_dist <- numeric(nboot)
            for (b in seq_len(nboot)) {
                idx <- sample(seq_len(n_obs), size = n_obs, replace = TRUE)
                x_boot <- x_for_bootstrap[idx]
                x_boot_sum <- sum(x_boot)
                if (x_boot_sum <= 1e-10) {
                    bootstrap_dist[b] <- NA_real_
                    next
                }
                p_boot <- x_boot / x_boot_sum
                if (what == "S") {
                    bootstrap_dist[b] <- entropy_cpp(p_boot, q, norm, log_base)
                } else {
                    h <- entropy_cpp(p_boot, q, FALSE, log_base)
                    if (!is.finite(h)) {
                        bootstrap_dist[b] <- NA_real_
                    } else if (abs(q - 1) < 1e-06) {
                        bootstrap_dist[b] <- exp(h)
                    } else {
                        base <- 1 - (q - 1) * h
                        bootstrap_dist[b] <- if (base > 0) base^(1/(1 - q)) else NA_real_
                    }
                }
            }
            return(bootstrap_dist)
        }
    }

    # Standard (independent) bootstrap resampling
    if (what == "S") {
        # For entropy
        bootstrap_dist <- bootstrap_compute_cpp_wrapper(x = x_for_bootstrap, q = q,
            normalize = norm, nboot = nboot, log_base = log_base, pseudocount = pseudocount)
    } else if (what == "D") {
        # For Hill numbers: compute entropy then convert
        bootstrap_dist <- bootstrap_compute_cpp_wrapper(x = x_for_bootstrap, q = q,
            normalize = FALSE, nboot = nboot, log_base = log_base, pseudocount = pseudocount)
        # Hill number conversion: D_q = (1 - (q-1) * H_q)^(1/(1-q))
        if (abs(q - 1) < 1e-06) {
            bootstrap_dist <- exp(bootstrap_dist)  # exp(H) for q=1
        } else {
            # Compute base: 1 - (q-1) * H
            base <- 1 - (q - 1) * bootstrap_dist
            
            # Defensive: check for negative bases (indicates entropy outside valid range)
            n_negative <- sum(base < 0, na.rm = TRUE)
            if (n_negative > 0) {
                warning("Hill number conversion produced ", n_negative, 
                        " negative base values for q=", q, ". ",
                        "This suggests entropy values exceed valid range. ",
                        "Clamping to small positive value (1e-10).")
                base <- pmax(base, 1e-10)
            }
            
            # Apply Hill number transformation
            bootstrap_dist <- base^(1/(1 - q))
        }
    } else {
        stop("Invalid 'what' parameter: must be 'S' (entropy) or 'D' (Hill numbers)")
    }

    # Validation: Check for NaN/Inf in bootstrap distribution
    n_nan <- sum(is.nan(bootstrap_dist))
    n_inf <- sum(is.infinite(bootstrap_dist))
    n_total <- length(bootstrap_dist)

    if (n_nan > 0 || n_inf > 0) {
        warning("Bootstrap resampling produced ", n_nan, " NaN and ", n_inf, " Inf values ",
            "out of ", n_total, " replicates. ", "This typically indicates all-zero counts or numerical instability. ",
            "Consider checking input data or adding pseudocount.")
    }

    return(bootstrap_dist)
}

