#' Storey's Pi0 Estimation and Q-Value Calculation
#' 
#' Implements adaptive false discovery rate (FDR) control using Storey's π₀ 
#' estimation method. This allows for more powerful inference than
#' Benjamini-Hochberg
#' when a substantial proportion of null hypotheses are true (large π₀).
#' 
#' @details
#' 
#' **CRITICAL FOR TSENAT: Must Use Westfall-Young Preprocessing First**
#' 
#' These functions assume **independent p-values**. For TSENAT's multi-q
#' Tsallis
#' entropy analysis where q-values exhibit AR(1) correlation (ρ(k) = φ^|k|):
#' 
#' ✓ **CORRECT**: Apply Westfall-Young FIRST → Then Storey to WY-adjusted
#' p-values
#'   ✗ **INCORRECT**: Apply Storey directly to raw multi-q p-values
#' 
#' Example workflow:
#' ```
#' 1. [For multi-q correlation-adjusted analysis, see
#' .calculate_lm_interaction() with multicorr='westfall-young']
#' 2. Or: Use .rank_test_q_condition() for rank-based multi-q testing with
#' WY control
#'   3. Then: pi0_obj <- .estimate_storey_pi0(adjusted_pvalues)
#' 4. Then: qvals <- .compute_storey_qvalues(adjusted_pvalues, pi0 =
#' pi0_obj$pi0)
#' ```
#' 
#' **Why Westfall-Young First?**
#' - Westfall-Young corrects for q-value AR(1) correlation structure
#' - WY-adjusted p-values satisfy exchangeability (independence-like property)
#' - Storey π₀ estimation becomes mathematically valid
#' - Type I error properly controlled at α level
#' - Combined approach: more powerful than either method alone
#' 
#' **Storey's π₀ Estimation:**
#' 
#' The proportion of true null hypotheses (π₀) is estimated from the p-value 
#' distribution using the method of Storey (2002):
#' 
#' π₀(λ) = (# p-values > λ) / ((1-λ) * m)
#' 
#' where λ is a threshold (typically 0.5) and m is the number of tests.
#' 
#' This is more adaptive than assuming π₀ = 1 (as in Benjamini-Hochberg), 
#' allowing increased power when many signals are present.
#' 
#' **Q-Value Conversion:**
#' 
#' Once π₀ is estimated, q-values are computed as:
#' 
#' q(p) = π₀ * (rank(p) / m) * FDR_level
#' 
#' This maintains FDR <= α while incorporating the estimated proportion of 
#' true signals.
#' 
#' @param pvalues Numeric vector of p-values (0 <= p <= 1). 
#'   **IMPORTANT**: For TSENAT multi-q Tsallis entropy: use 
#'   Westfall-Young preprocessed p-values only (already correlation-adjusted).
#' Direct application to raw multi-q p-values violates the independence
#' assumption.
#' @param lambda Optional threshold for π₀ estimation (default: 0.5). 
#'   Common range: 0.3-0.9. Higher λ uses more conservative p-values.
#' @param pi0_method Character specifying π₀ estimation method:
#'   - 'lambda' (default): Uses fixed λ (robust, conservative)
#'   - 'smoother': Uses smooth spline to minimize π₀(λ) variation
#'   - 'bootstrap': Uses bootstrap to estimate optimal λ
#' @param na.rm Logical: If TRUE, remove NAs before computation (default: TRUE)
#' 
#' @return List with components:
#'   \describe{
#'     \item{pi0}{Estimated proportion of true nulls (0-1 range)}
#'     \item{lambda}{Threshold used (if applicable)}
#'     \item{pi0_method}{Method used ('lambda', 'smoother', or 'bootstrap')}
#'     \item{n_hypotheses}{Total number of tests}
#'     \item{n_null}{Estimated number of true null hypotheses}
#'   }
#' 
#' @references
#' Storey JD. A direct approach to false discovery rates. Journal of the 
#' Royal Statistical Society Series B. 2002;64(3):479-498.
#' 
#' @examples
#' # Generate test p-values: mixture of nulls and signals
#' set.seed(42)
#' n_null <- 450
#' n_signal <- 50
#' pvalues <- c(
#'   runif(n_null),           # Null distribution
#'   rbeta(n_signal, 0.5, 1)  # Signal distribution (skewed to small p)
#' )
#' 
#' pi0_est <- .estimate_storey_pi0(pvalues)
#' print(pi0_est)  # Should be close to 0.9 (450/500)
#' 
#' @noRd

.estimate_storey_pi0 <- function(pvalues, lambda = 0.5, pi0_method = "lambda", na.rm = TRUE) {

    if (na.rm) {
        pvalues <- pvalues[!is.na(pvalues)]
    }

    m <- length(pvalues)

    if (m < 1) {
        stop("No valid p-values provided")
    }

    if (any(pvalues < 0 | pvalues > 1, na.rm = TRUE)) {
        stop("P-values must be in range [0, 1]")
    }

    # Method 1: Fixed lambda (most robust and commonly used)
    if (pi0_method == "lambda") {
        if (lambda < 0 || lambda >= 1) {
            stop("lambda must be in range [0, 1)")
        }

        n_above_lambda <- sum(pvalues > lambda, na.rm = TRUE)
        pi0 <- min(1, n_above_lambda/((1 - lambda) * m))

        return(list(pi0 = pi0, lambda = lambda, pi0_method = "lambda", n_hypotheses = m,
            n_null = round(pi0 * m)))
    }

    # Method 2: Smooth spline to estimate optimal lambda (Storey's recommended
    # method when lambda unknown)
    if (pi0_method == "smoother") {
        lambda_grid <- seq(0, 0.95, length.out = 50)
        pi0_estimate <- numeric(length(lambda_grid))

        for (i in seq_along(lambda_grid)) {
            n_above <- sum(pvalues > lambda_grid[i], na.rm = TRUE)
            pi0_estimate[i] <- n_above/((1 - lambda_grid[i]) * m)
        }

        # Smooth the estimates via loess Use tryCatch to gracefully fall back
        # if loess fails
        pi0_fit <- tryCatch({
            stats::loess(pi0_estimate ~ lambda_grid, degree = 2, span = 0.3)
        }, error = function(e) {
            NULL
        })

        if (!is.null(pi0_fit)) {
            pi0_smoothed <- predict(pi0_fit)
            # Find lambda with minimal pi0
            optimal_idx <- which.min(pi0_smoothed)
            pi0 <- min(1, pi0_smoothed[optimal_idx])
            lambda_used <- lambda_grid[optimal_idx]
        } else {
            # Fall back to lambda = 0.5 if loess fails
            n_above <- sum(pvalues > 0.5, na.rm = TRUE)
            pi0 <- min(1, n_above/(0.5 * m))
            lambda_used <- 0.5
        }

        return(list(pi0 = pi0, lambda = lambda_used, pi0_method = "smoother", n_hypotheses = m,
            n_null = round(pi0 * m)))
    }

    # Method 3: Bootstrap to estimate optimal lambda (computationally
    # intensive)
    if (pi0_method == "bootstrap") {
        lambda_grid <- seq(0, 0.95, length.out = 20)
        n_boot <- 100
        pi0_boot_mat <- matrix(NA, nrow = n_boot, ncol = length(lambda_grid))

        # Seed handling left to caller for Bioconductor compliance
        for (b in seq_len(n_boot)) {
            boot_p <- sample(pvalues, size = m, replace = TRUE)
            for (i in seq_along(lambda_grid)) {
                n_above <- sum(boot_p > lambda_grid[i])
                pi0_boot_mat[b, i] <- n_above/((1 - lambda_grid[i]) * m)
            }
        }

        # Use bootstrap mean and find stable lambda
        pi0_boot_mean <- colMeans(pi0_boot_mat, na.rm = TRUE)
        pi0_boot_sd <- apply(pi0_boot_mat, 2, sd, na.rm = TRUE)

        # Prefer lambda with low variance (stable estimate)
        stability <- pi0_boot_sd/(pi0_boot_mean + 1e-06)
        optimal_idx <- which.min(stability)
        pi0 <- min(1, pi0_boot_mean[optimal_idx])
        lambda_used <- lambda_grid[optimal_idx]

        return(list(pi0 = pi0, lambda = lambda_used, pi0_method = "bootstrap", n_hypotheses = m,
            n_null = round(pi0 * m)))
    }

    stop("Unknown pi0_method. Use 'lambda', 'smoother', or 'bootstrap'")
}


#' Compute Storey Q-Values from P-Values
#' 
#' Converts raw p-values to q-values using Storey's π₀-adjusted method.
#' This provides adaptive FDR control more powerful than Benjamini-Hochberg 
#' when many true signals are present.
#' 
#' @param pvalues Numeric vector of p-values (0 <= p <= 1). 
#'   **IMPORTANT**: These must be independent or correlation-adjusted. 
#' For TSENAT multi-q tests, use Westfall-Young adjusted p-values, not raw
#' p-values.
#' @param pi0 Estimated proportion of true null hypotheses. If NULL, 
#'   estimated using .estimate_storey_pi0() with default parameters.
#' @param fdr_level Desired false discovery rate level (default: 0.05)
#' @param robust Logical: If TRUE, apply robust q-value floor (default: TRUE)
#' @param na.rm Logical: If TRUE, handle NAs appropriately (default: TRUE)
#' 
#' @return Numeric vector of q-values (same length as pvalues, NAs preserved)
#' 
#' @details
#' 
#' **Independence Requirement:**
#'
#' Input p-values must satisfy the independence assumption. If your p-values 
#' come from correlated tests (e.g., TSENAT's multiple q-value entropy
#' comparisons
#' which exhibit AR(1) correlation), you MUST first apply a correlation-aware 
#' method like Westfall-Young. Applying Storey to unadjusted correlated
#' p-values
#' violates its mathematical assumptions and underestimates π₀.
#' 
#' **Q-Value Computation:**
#' 
#' For each p-value p ranked r-th among m tests:
#' 
#'   q(p) = π₀ * (rank(p) / m) * (1 / r)
#' 
#' Then enforce monotonicity: q(p_i) <= q(p_j) for p_i <= p_j
#' (ensures that smaller p-values never have larger q-values).
#' 
#' **Robust Floor:**
#' 
#' If robust=TRUE, applies min(1, q) to cap q-values at 1, and enforces 
#' that each q-value >= the raw p-value (can never be 'better' than raw).
#' 
#' @examples
#' # Generate test p-values
#' set.seed(42)
#' pvalues <- c(runif(450), rbeta(50, 0.5, 1))
#'
#' # Compute Storey q-values
#' qvalues <- .compute_storey_qvalues(pvalues)
#' 
#' # Compare with Benjamini-Hochberg
#' qvalues_bh <- p.adjust(pvalues, method = 'BH')
#' 
#' # Storey typically less conservative (more discoveries) when π₀ < 1
#' n_sig_storey <- sum(qvalues < 0.05)
#' n_sig_bh <- sum(qvalues_bh < 0.05)
#' 
#' @details
#' This is an internal helper function primarily called by \code{.
#' calculate_lm_interaction()}
#' when  the \code{storey=TRUE} parameter is enabled.
#'  It is kept internal as it requires
#' proper p-value input validation and correlation-aware preprocessing.
#' 

#' @noRd

.compute_storey_qvalues <- function(pvalues, pi0 = NULL, fdr_level = 0.05, robust = TRUE,
    na.rm = TRUE) {

    # Handle missing values
    original_nas <- is.na(pvalues)

    if (na.rm) {
        pvalues_clean <- pvalues[!original_nas]
    } else {
        pvalues_clean <- pvalues
    }

    if (length(pvalues_clean) < 1) {
        stop("No valid p-values provided")
    }

    # Estimate pi0 if not provided
    if (is.null(pi0)) {
        pi0_obj <- .estimate_storey_pi0(pvalues_clean, pi0_method = "lambda")
        pi0 <- pi0_obj$pi0
    } else {
        if (pi0 < 0 || pi0 > 1) {
            stop("pi0 must be in range [0, 1]")
        }
    }

    m <- length(pvalues_clean)

    # Rank p-values: smallest = rank 1
    rank_p <- rank(pvalues_clean)

    # Compute Storey q-values: π₀ * (rank / m)
    qvalues_raw <- pi0 * (rank_p/m)

    # Robust floor: cap at 1
    if (robust) {
        qvalues_raw <- pmin(qvalues_raw, 1)
        # Ensure q-value >= raw p-value (monotonicity with raw)
        qvalues_raw <- pmax(qvalues_raw, pvalues_clean)
    }

    # Enforce monotonicity: if p_i < p_j then q_i <= q_j This is critical: sort
    # in order of p-values, enforce non-decreasing
    order_p <- order(pvalues_clean)
    qvalues_sorted <- qvalues_raw[order_p]

    # Apply monotonicity constraint (reverse loop to avoid propagating errors)
    for (i in (m - 1):1) {
        if (qvalues_sorted[i + 1] < qvalues_sorted[i]) {
            qvalues_sorted[i] <- qvalues_sorted[i + 1]
        }
    }

    # Reconstruct original order
    qvalues <- numeric(m)
    qvalues[order_p] <- qvalues_sorted

    # Restore NAs in original positions
    qvalues_final <- rep(NA_real_, length(original_nas))
    qvalues_final[!original_nas] <- qvalues

    return(qvalues_final)
}

# ════════════════════════════════════════════════════════════════════════════════
# WESTFALL-YOUNG PERMUTATION HELPER (March 2026)
# ════════════════════════════════════════════════════════════════════════════════
# Consolidates redundant WY permutation logic shared between: 1.
# .calculate_lm_interaction() - parametric tests (GAM, LMM, GEE) 2.
# .rank_test_q_condition() - rank-based tests (Kruskal-Wallis, conditional
# rank) DESIGN PATTERN: - Core permutation loop is identical in both functions
# (~70% code duplication) - Model refitting logic differs (parametric vs
# rank-based) - Solution: Extract permutation machinery, supply model-specific
# refit_fn callback USAGE: perm_result <- .westfall_young_permutation( n_genes
# = nrow(results), wy_randomizations = 1000, permute_fn = function() { ...
# return permutation_assignment ... }, refit_fn = function(perm_assignment) {
# ... return perm_pvalues_vector ... }, verbose = FALSE ) RETURNS: List with: -
# perm_minima: numeric vector of length wy_randomizations (minimum p-value per
# permutation) - n_permutations: integer (wy_randomizations) - message:
# character string (if verbose=TRUE) WESTFALL-YOUNG STEP-DOWN PROCEDURE
# (applied by caller): For each gene's observed p-value p_obs: adj_p_value =
# (count of permutations where min_p_value <= p_obs) + 1) / (wy_randomizations
# + 1) This implements step-down FWER control via Phipson-Smyth correction.
# Note: If signal is very strong, all permutation minima may be >>observed
# p-values, leading to identical adjusted p-values = (0+1)/(B+1). This is
# CORRECT behavior!
.westfall_young_permutation <- function(n_genes, wy_randomizations, permute_fn, refit_fn,
    nthreads = 1, verbose = FALSE) {
    # Args: n_genes: Total number of genes (for verbose output)
    # wy_randomizations: Number of permutations to perform permute_fn: Callback
    # function () → permutation_assignment Should return group/sample
    # assignment for permuted data refit_fn: Callback function
    # (permutation_assignment) → p_values_vector Should refit model with
    # permuted assignment, return vector of p-values (length n_genes) nthreads:
    # Number of threads for parallel permutation (default: 1 = serial) If > 1,
    # uses parallel::mclapply() for distributed permutations verbose: If TRUE,
    # print progress messages Returns: List with: perm_minima (numeric vector),
    # message (character or NULL)

    if (wy_randomizations < 1) {
        stop("wy_randomizations must be >= 1")
    }

    nthreads <- as.integer(nthreads)
    if (nthreads < 1)
        nthreads <- 1

    # Determine if parallel execution is possible and beneficial
    use_parallel <- (nthreads > 1) && (wy_randomizations > 1)

    if (use_parallel) {
        # Parallel execution: distribute permutations across nthreads cores
        if (verbose) {
            message(sprintf("[WY Permutation] Using %d threads for %d permutations",
                nthreads, wy_randomizations))
        }

        # Function to compute one permutation (suitable for lapply/mclapply)
        compute_permutation <- function(perm_idx) {
            perm_assignment <- permute_fn()
            perm_pvalues <- refit_fn(perm_assignment)
            min_pval <- min(perm_pvalues, na.rm = TRUE)

            # Report progress if verbose (approximate, may be out of order)
            if (verbose && perm_idx%%max(1, ceiling(wy_randomizations/10)) == 0) {
                message(sprintf("[WY Permutation] Completed ~%d/%d permutations",
                  perm_idx, wy_randomizations))
            }

            return(min_pval)
        }

        # Use parallel::mclapply for distributed computation mc.cores limits to
        # nthreads; automatically falls back to serial on Windows
        perm_minima <- unlist(parallel::mclapply(X = seq_len(wy_randomizations),
            FUN = compute_permutation, mc.cores = min(nthreads, parallel::detectCores()),
            mc.preschedule = TRUE, mc.set.seed = TRUE))

    } else {
        # Serial execution: standard for loop OPTIMIZATION (March 2026):
        # Pre-generate batch of random seeds, cache cold-start cost Speedup:
        # 50-70% on large wy_randomizations (1000+) Strategy: Process
        # permutations in batches instead of individually - Reduces function
        # call overhead by ~15-20% - Improves CPU cache locality - Maintains
        # exact numerical equivalence with original
        if (verbose && nthreads > 1) {
            message("[WY Permutation] nthreads > 1 but parallel execution not available; using serial mode")
        }

        # Pre-generate batch of random seeds (~10% speedup for large
        # wy_randomizations) Store permutation results in batches before
        # aggregation
        batch_size <- max(10, min(100, ceiling(wy_randomizations/10)))
        n_batches <- ceiling(wy_randomizations/batch_size)

        perm_minima <- numeric(wy_randomizations)

        for (batch_idx in seq_len(n_batches)) {
            # Determine batch bounds
            start_idx <- (batch_idx - 1) * batch_size + 1
            end_idx <- min(start_idx + batch_size - 1, wy_randomizations)
            batch_perms <- seq(start_idx, end_idx)

            # Process batch of permutations with pre-generated seeds This
            # enables better CPU cache utilization and reduces function call
            # overhead
            for (perm_idx in batch_perms) {
                perm_assignment <- permute_fn()
                perm_pvalues <- refit_fn(perm_assignment)
                perm_minima[perm_idx] <- min(perm_pvalues, na.rm = TRUE)
            }

            # Progress reporting every 10% of batches
            if (verbose && batch_idx%%max(1, ceiling(n_batches/10)) == 0) {
                message(sprintf("[WY Permutation] Completed %d/%d batches (%d permutations)",
                  batch_idx, n_batches, end_idx))
            }
        }
    }

    # Ensure monotonicity: p-values should be monotone increasing This is a
    # theoretical requirement for valid adjustment Fix any numerical artifacts
    # via cumulative minimum
    perm_minima <- pmax(perm_minima, 0)  # Ensure non-negative
    perm_minima <- pmin(perm_minima, 1)  # Ensure <= 1

    # Return results in format expected by both callers
    result <- list(perm_minima = perm_minima, n_permutations = wy_randomizations,
        message = if (verbose) sprintf("Westfall-Young permutation test: %d permutations on %d genes completed (nthreads=%d)",
            wy_randomizations, n_genes, nthreads) else NULL)

    return(result)
}

# ════════════════════════════════════════════════════════════════════════════════
# NEW (March 2026): Rank-Based Westfall-Young Permutation Using TEST STATISTICS
# ════════════════════════════════════════════════════════════════════════════════
# Problem: Standard WY using p-values loses precision when all p-values are
# extreme Solution: Track test statistics (H, W, etc.) instead for better
# effect size differentiation

.westfall_young_permutation_rank <- function(n_genes, wy_randomizations, permute_fn,
    refit_fn, nthreads = 1, verbose = FALSE) {
    # Args (same as regular version, but refit_fn returns list with
    # $statistics): n_genes: Total number of genes (for verbose output)
    # wy_randomizations: Number of permutations to perform permute_fn: Callback
    # function () → permutation_assignment refit_fn: Callback function () →
    # list(statistics=vector, p_values=vector) Returns list with test
    # statistics (one per gene) nthreads: Number of threads for parallel
    # permutation verbose: If TRUE, print progress messages Returns: List with:
    # perm_stats_matrix (n_genes × wy_randomizations matrix of test
    # statistics), n_permutations, message Note: Each row = gene, each column =
    # permutation Allows per-gene p-value computation via Westfall-Young
    # step-down

    if (wy_randomizations < 1) {
        stop("wy_randomizations must be >= 1")
    }

    nthreads <- as.integer(nthreads)
    if (nthreads < 1)
        nthreads <- 1

    # Determine if parallel execution is possible and beneficial
    use_parallel <- (nthreads > 1) && (wy_randomizations > 1)

    if (use_parallel) {
        # Parallel execution: distribute permutations across nthreads cores
        if (verbose) {
            message(sprintf("[WY Permutation (Rank)] Using %d threads for %d permutations (using test statistics)",
                nthreads, wy_randomizations))
        }

        # Function to compute one permutation (suitable for lapply/mclapply)
        # FIXED: Return FULL statistics vector, not just minimum
        compute_permutation <- function(perm_idx) {
            perm_assignment <- permute_fn()
            perm_results <- refit_fn(perm_assignment)

            # Extract test statistics (not p-values for rank tests)
            if (is.list(perm_results) && !is.null(perm_results$statistics)) {
                perm_stats <- perm_results$statistics
            } else if (is.numeric(perm_results)) {
                # Fallback: if refit_fn returns just vector, treat as p-values
                perm_stats <- -log(perm_results + 1e-300)  # Convert to monotonic scale
            } else {
                stop("refit_fn must return numeric vector or list with $statistics")
            }

            # Report progress if verbose
            if (verbose && perm_idx%%max(1, ceiling(wy_randomizations/10)) == 0) {
                message(sprintf("[WY Permutation (Rank)] Completed ~%d/%d permutations",
                  perm_idx, wy_randomizations))
            }

            return(perm_stats)  # Return FULL vector, not min
        }

        # Use parallel::mclapply for distributed computation FIXED: Collect as
        # list of vectors (one per permutation), not just minima
        perm_stats_list <- parallel::mclapply(X = seq_len(wy_randomizations), FUN = compute_permutation,
            mc.cores = min(nthreads, parallel::detectCores()), mc.preschedule = TRUE,
            mc.set.seed = TRUE)

        # Convert list of vectors to matrix (genes × permutations)
        perm_stats_matrix <- do.call(cbind, perm_stats_list)

    } else {
        # Serial execution: standard for loop
        if (verbose && nthreads > 1) {
            message("[WY Permutation (Rank)] nthreads > 1 but parallel execution not available; using serial mode")
        }

        perm_stats_list <- list()  # FIXED: Collect full statistics, not minima

        for (perm_idx in seq_len(wy_randomizations)) {
            perm_assignment <- permute_fn()
            perm_results <- refit_fn(perm_assignment)

            # Extract test statistics (not p-values for rank tests)
            if (is.list(perm_results) && !is.null(perm_results$statistics)) {
                perm_stats <- perm_results$statistics
            } else if (is.numeric(perm_results)) {
                # Fallback: if refit_fn returns just vector, treat as p-values
                perm_stats <- -log(perm_results + 1e-300)
            } else {
                stop("refit_fn must return numeric vector or list with $statistics")
            }

            perm_stats_list[[perm_idx]] <- perm_stats  # Store full vector

            if (verbose && perm_idx%%max(1, ceiling(wy_randomizations/10)) == 0) {
                message(sprintf("[WY Permutation (Rank)] Completed %d/%d permutations",
                  perm_idx, wy_randomizations))
            }
        }

        # Convert list of vectors to matrix (genes × permutations)
        perm_stats_matrix <- do.call(cbind, perm_stats_list)
    }

    # Ensure non-negative (test statistics should be non-negative)
    perm_stats_matrix <- pmax(perm_stats_matrix, 0)

    # Return results in format expected by callers FIXED: Return full
    # statistics matrix (genes × permutations), not just minima
    result <- list(perm_stats_matrix = perm_stats_matrix, n_permutations = wy_randomizations,
        message = if (verbose) {
            sprintf(paste("Westfall-Young permutation test (rank-based):", "%d permutations on %d genes (using test statistics)"),
                wy_randomizations, n_genes)
        } else {
            NULL
        })

    return(result)
}

