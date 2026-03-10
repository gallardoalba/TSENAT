#' Storey's Pi0 Estimation and Q-Value Calculation
#' 
#' Implements adaptive false discovery rate (FDR) control using Storey's π₀ 
#' estimation method. This allows for more powerful inference than Benjamini-Hochberg 
#' when a substantial proportion of null hypotheses are true (large π₀).
#' 
#' @details
#' 
#' **CRITICAL FOR TSENAT: Must Use Westfall-Young Preprocessing First**
#' 
#' These functions assume **independent p-values**. For TSENAT's multi-q Tsallis 
#' entropy analysis where q-values exhibit AR(1) correlation (ρ(k) = φ^|k|):
#' 
#'   ✓ **CORRECT**: Apply Westfall-Young FIRST → Then Storey to WY-adjusted p-values
#'   ✗ **INCORRECT**: Apply Storey directly to raw multi-q p-values
#' 
#' Example workflow:
#' ```
#'   1. Run: wy_result <- label_shuffling_westfall_young(...)
#'   2. Get: wy_pvalues <- wy_result$pvalue_raw (correlation-adjusted)
#'   3. Then: pi0_obj <- estimate_storey_pi0(wy_pvalues)
#'   4. Then: qvals <- compute_storey_qvalues(wy_pvalues, pi0 = pi0_obj$pi0)
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
#'   Direct application to raw multi-q p-values violates the independence assumption.
#' @param lambda Optional threshold for π₀ estimation (default: 0.5). 
#'   Common range: 0.3-0.9. Higher λ uses more conservative p-values.
#' @param pi0_method Character specifying π₀ estimation method:
#'   - "lambda" (default): Uses fixed λ (robust, conservative)
#'   - "smoother": Uses smooth spline to minimize π₀(λ) variation
#'   - "bootstrap": Uses bootstrap to estimate optimal λ
#' @param na.rm Logical: If TRUE, remove NAs before computation (default: TRUE)
#' 
#' @return List with components:
#'   \describe{
#'     \item{pi0}{Estimated proportion of true nulls (0-1 range)}
#'     \item{lambda}{Threshold used (if applicable)}
#'     \item{pi0_method}{Method used ("lambda", "smoother", or "bootstrap")}
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
#' pi0_est <- estimate_storey_pi0(pvalues)
#' print(pi0_est)  # Should be close to 0.9 (450/500)
#' 
#' @export
estimate_storey_pi0 <- function(pvalues, lambda = 0.5, pi0_method = "lambda", 
                                 na.rm = TRUE) {
  
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
    pi0 <- min(1, n_above_lambda / ((1 - lambda) * m))
    
    return(list(
      pi0 = pi0,
      lambda = lambda,
      pi0_method = "lambda",
      n_hypotheses = m,
      n_null = round(pi0 * m)
    ))
  }
  
  # Method 2: Smooth spline to estimate optimal lambda
  # (Storey's recommended method when lambda unknown)
  if (pi0_method == "smoother") {
    lambda_grid <- seq(0, 0.95, length.out = 50)
    pi0_estimate <- numeric(length(lambda_grid))
    
    for (i in seq_along(lambda_grid)) {
      n_above <- sum(pvalues > lambda_grid[i], na.rm = TRUE)
      pi0_estimate[i] <- n_above / ((1 - lambda_grid[i]) * m)
    }
    
    # Smooth the estimates via loess
    # Use tryCatch to gracefully fall back if loess fails
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
      pi0 <- min(1, n_above / (0.5 * m))
      lambda_used <- 0.5
    }
    
    return(list(
      pi0 = pi0,
      lambda = lambda_used,
      pi0_method = "smoother",
      n_hypotheses = m,
      n_null = round(pi0 * m)
    ))
  }
  
  # Method 3: Bootstrap to estimate optimal lambda (computationally intensive)
  if (pi0_method == "bootstrap") {
    lambda_grid <- seq(0, 0.95, length.out = 20)
    n_boot <- 100
    pi0_boot_mat <- matrix(NA, nrow = n_boot, ncol = length(lambda_grid))
    
    set.seed(12345)  # For reproducibility
    for (b in seq_len(n_boot)) {
      boot_p <- sample(pvalues, size = m, replace = TRUE)
      for (i in seq_along(lambda_grid)) {
        n_above <- sum(boot_p > lambda_grid[i])
        pi0_boot_mat[b, i] <- n_above / ((1 - lambda_grid[i]) * m)
      }
    }
    
    # Use bootstrap mean and find stable lambda
    pi0_boot_mean <- colMeans(pi0_boot_mat, na.rm = TRUE)
    pi0_boot_sd <- apply(pi0_boot_mat, 2, sd, na.rm = TRUE)
    
    # Prefer lambda with low variance (stable estimate)
    stability <- pi0_boot_sd / (pi0_boot_mean + 1e-6)
    optimal_idx <- which.min(stability)
    pi0 <- min(1, pi0_boot_mean[optimal_idx])
    lambda_used <- lambda_grid[optimal_idx]
    
    return(list(
      pi0 = pi0,
      lambda = lambda_used,
      pi0_method = "bootstrap",
      n_hypotheses = m,
      n_null = round(pi0 * m)
    ))
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
#'   For TSENAT multi-q tests, use Westfall-Young adjusted p-values, not raw p-values.
#' @param pi0 Estimated proportion of true null hypotheses. If NULL, 
#'   estimated using estimate_storey_pi0() with default parameters.
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
#' come from correlated tests (e.g., TSENAT's multiple q-value entropy comparisons
#' which exhibit AR(1) correlation), you MUST first apply a correlation-aware 
#' method like Westfall-Young. Applying Storey to unadjusted correlated p-values 
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
#' that each q-value >= the raw p-value (can never be "better" than raw).
#' 
#' @examples
#' # Generate test p-values
#' set.seed(42)
#' pvalues <- c(runif(450), rbeta(50, 0.5, 1))
#'
#' # Compute Storey q-values
#' qvalues <- compute_storey_qvalues(pvalues)
#' 
#' # Compare with Benjamini-Hochberg
#' qvalues_bh <- p.adjust(pvalues, method = "BH")
#' 
#' # Storey typically less conservative (more discoveries) when π₀ < 1
#' n_sig_storey <- sum(qvalues < 0.05)
#' n_sig_bh <- sum(qvalues_bh < 0.05)
#' 
#' @export
compute_storey_qvalues <- function(pvalues, pi0 = NULL, fdr_level = 0.05, 
                                    robust = TRUE, na.rm = TRUE) {
  
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
    pi0_obj <- estimate_storey_pi0(pvalues_clean, pi0_method = "lambda")
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
  qvalues_raw <- pi0 * (rank_p / m)
  
  # Robust floor: cap at 1
  if (robust) {
    qvalues_raw <- pmin(qvalues_raw, 1)
    # Ensure q-value >= raw p-value (monotonicity with raw)
    qvalues_raw <- pmax(qvalues_raw, pvalues_clean)
  }
  
  # Enforce monotonicity: if p_i < p_j then q_i <= q_j
  # This is critical: sort in order of p-values, enforce non-decreasing
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


#' Integrated Storey Analysis: π₀ Estimation + Q-Value Computation
#' 
#' Convenient wrapper that estimates π₀ and computes Storey q-values,
#' returning all information needed for interpretation.
#' 
#' @param pvalues Numeric vector of p-values. **CRITICAL**: Must be independent 
#'   or correlation-adjusted. For TSENAT multi-q entropy tests, use 
#'   Westfall-Young adjusted p-values. See \code{\link{estimate_storey_pi0}} 
#'   for details on the independence requirement.
#' @param pi0_method Character: "lambda", "smoother", or "bootstrap"
#' @param lambda Threshold for π₀ estimation (if pi0_method = "lambda")
#' @param fdr_level Desired FDR level (default: 0.05)
#' @param na.rm Logical: Handle NAs (default: TRUE)
#' 
#' @return List with components:
#'   \describe{
#'     \item{pvalues}{Input p-values}
#'     \item{qvalues}{Computed Storey q-values}
#'     \item{pi0}{Estimated proportion of true nulls}
#'     \item{lambda}{Lambda parameter used for \eqn{\pi_0}{pi_0} estimation}
#'     \item{pi0_method}{Method used}
#'     \item{n_significant}{Number of significant tests at given FDR level}
#'     \item{fdr_level}{FDR level used}
#'     \item{power_gain}{Comparison with Benjamini-Hochberg ((n_sig_storey - n_sig_bh) / n_sig_bh)}
#'   }
#'   
#' @export
storey_analysis <- function(pvalues, pi0_method = "lambda", lambda = 0.5, 
                            fdr_level = 0.05, na.rm = TRUE) {
  
  # Estimate π₀
  pi0_result <- estimate_storey_pi0(pvalues, lambda = lambda, 
                                     pi0_method = pi0_method, na.rm = na.rm)
  
  # Compute Storey q-values
  qvalues <- compute_storey_qvalues(pvalues, pi0 = pi0_result$pi0, 
                                    fdr_level = fdr_level, na.rm = na.rm)
  
  # Compare with Benjamini-Hochberg
  qvalues_bh <- p.adjust(pvalues, method = "BH")
  
  n_sig_storey <- sum(qvalues < fdr_level, na.rm = TRUE)
  n_sig_bh <- sum(qvalues_bh < fdr_level, na.rm = TRUE)
  
  power_gain <- NA_real_
  if (n_sig_bh > 0) {
    power_gain <- (n_sig_storey - n_sig_bh) / n_sig_bh
  } else if (n_sig_storey > 0) {
    power_gain <- Inf
  }
  
  return(list(
    pvalues = pvalues,
    qvalues = qvalues,
    pi0 = pi0_result$pi0,
    lambda = pi0_result$lambda,
    pi0_method = pi0_result$pi0_method,
    n_hypotheses = length(pvalues),
    n_significant = n_sig_storey,
    n_significant_bh = n_sig_bh,
    fdr_level = fdr_level,
    power_gain = power_gain
  ))
}
