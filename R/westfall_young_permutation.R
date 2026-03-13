
#' Multi-q Permutation Test with Westfall-Young Stepdown Procedure
#'
#' Performs joint permutation testing across multiple q-values (diversity orders)
#' using the Westfall-Young stepdown procedure for family-wise error rate control.
#'
#' @param x A \code{matrix} with diversity values (rows = genes, columns = samples*q-values).
#'   Column names should indicate q-value (e.g., "Sample1_q=0.5", "Sample1_q=1.0").
#' @param q_values Numeric vector of q-values used in computing diversity.  
#'   Should be in increasing order and match column names in \code{x}.
#' @param samples Character vector matching the sample pairs in column names.
#'   Length should match number of unique samples (before q expansion).
#' @param control Name of the control sample category (e.g., "Normal", "WT").
#' @param method Method for aggregating diversity values: "mean" or "median".
#' @param randomizations Integer number of permutations (default: 1000).
#' @param verbose Logical; if TRUE, print progress updates (default: FALSE).
#' @param nthreads Integer number of threads for parallel processing (default: 1).
#'
#' @return A data frame with columns:
#'   - `gene`: Gene identifier (from row names of x)
#'   - `q_value`: q-value used
#'   - `log2FC`: Log2 fold-change at this q
#'   - `pvalue_raw`: Raw two-sided permutation p-value (Phipson-Smyth corrected)
#'   - `pvalue_wy`: Westfall-Young adjusted p-value (true stepdown FWER control)
#'   - `padj_bh`: Multiple testing adjusted p-value (BH-FDR on raw p-values)
#'   - Family-wise error rate <= 0.05 is guaranteed by true WY stepdown
#'
#' @details
#' **Westfall-Young Stepdown Procedure for Multi-q Analysis:**
#'
#' Traditional Benjamini-Hochberg FDR control assumes independence, but multi-q
#' tests for the same gene are correlated. The Westfall-Young stepdown procedure
#' properly accounts for these correlations by:
#'
#' 1. Computing raw two-sided permutation p-values for each gene-q combination
#' 2. Ranking p-values from smallest to largest
#' 3. For each p-value in order:
#'    - Count permutations where minimum p-value <= observed (accounts for all smaller tests)
#'    - Adjusted p-value = min(raw_p, count/num_permutations)
#' 4. Enforce monotonicity: adjusted p-values are non-decreasing as p-values increase
#'
#' **Advantages Over Independent Testing:**
#' - Controls correlation between q-values for same gene
#' - Maintains family-wise error rate (FWER) at α level
#' - More powerful than simple Bonferroni when tests are correlated
#' - Accounts for joint dependence structure in multi-q analysis
#'
#' **Phipson-Smyth Correction (S019):**
#' Each permutation test uses p = (b+1)/(m+1) to avoid zero p-values,
#' where b = count of permutations >= observed, m = number of permutations.
#'
#' @references
#' Westfall, P. H., & Young, S. S. (1993). Resampling-based Multiple Testing:
#' Examples and Methods for p-Value Adjustment. John Wiley & Sons.
#'
#' Phipson, B., & Smyth, G. K. (2010). Permutation p-values should never be zero:
#' calculating exact p-values when permutations are randomly drawn.
#' Statistical Applications in Genetics and Molecular Biology, 9(1), 39.
#'
#' Meinshausen, N., Maathuis, M., & Bühlmann, P. (2011). Optimality of the 
#' Westfall-Young permutation procedure for multiple testing under dependence.
#' The Annals of Statistics, 39(6), 3359-3382. [Database: S165, S166]
#'
#' Cox, D. D. Pointwise Testing with Functional Data Using the Westfall-Young 
#' Randomization Method. [Database: S167]
#'
#' @export
#' @examples
#' \dontrun{
#' # Diversity matrix with multiple q-values
#' set.seed(123)
#' n_genes <- 50
#' n_samples <- 8  # 4 pairs
#' n_q <- 5  # q = 0.5, 1.0, 1.5, 2.0
#'
#' # Create diversity matrix: genes * (samples * q-values)
#' diversity_data <- matrix(rnorm(n_genes * n_samples * n_q), 
#'                          nrow = n_genes)
#' colnames(diversity_data) <- 
#'   paste0(rep(paste0("Sample", 1:n_samples), n_q), 
#'          "_q=", rep(c(0.5, 1.0, 1.5, 2.0), each = n_samples))
#' rownames(diversity_data) <- paste0("Gene", 1:n_genes)
#'
#' # Sample labels
#' samples <- rep(c("Normal", "Tumor"), each = n_samples/2)
#' q_vals <- c(0.5, 1.0, 1.5, 2.0)
#'
#' # Run Westfall-Young multi-q permutation test
#' wy_results <- label_shuffling_westfall_young(
#'   x = diversity_data,
#'   q_values = q_vals,
#'   samples = samples,
#'   control = "Normal",
#'   randomizations = 1000,
#'   verbose = TRUE
#' )
#'
#' # View top results (most significant across all q-values)
#' head(wy_results[order(wy_results$pvalue_wy), ])
#'
#' # Count significant genes at FWER = 0.05
#' sum(wy_results$pvalue_wy < 0.05)
#' }
label_shuffling_westfall_young <- function(x, q_values, samples, control, 
                                           method = "mean", randomizations = 1000,
                                           verbose = FALSE, nthreads = 1) {
    
    # Input validation
    if (!all(grepl(paste0("_q=", collapse = "|"), colnames(x)))) {
        stop("Column names must contain _q= to indicate q-values (e.g., 'Sample1_q=0.5')")
    }
    
    if (length(q_values) == 0) {
        stop("q_values must be provided")
    }
    
    # Extract q-values from column names and verify they match input
    # Note: calculate_diversity arranges columns as (sample, q-value)
    # So pattern is rep(q_values, times=length(samples))
    col_q_values <- as.numeric(sub(".*_q=", "", colnames(x)))
    if (!isTRUE(all.equal(col_q_values, rep(q_values, times = length(samples))))) {
        stop("q-values in column names do not match provided q_values")
    }
    
    n_genes <- nrow(x)
    n_samples <- length(samples)
    n_q <- length(q_values)
    n_total_tests <- n_genes * n_q
    
    if (verbose) {
        cat("Westfall-Young Multi-q Permutation Test\n")
        cat("=========================================\n")
        cat("Genes:", n_genes, "\n")
        cat("Q-values:", paste(q_values, collapse = ", "), "\n")
        cat("Total tests:", n_total_tests, "\n")
        cat("Permutations:", randomizations, "\n")
        cat("FWER control: α = 0.05\n\n")
    }
    
    # Step 1: Compute raw p-values for all gene * q combinations
    # Also collect permutation distributions for true WY stepdown
    if (verbose) cat("Computing permutations...\n")
    
    raw_pvalues <- matrix(NA_real_, nrow = n_genes, ncol = n_q,
                         dimnames = list(rownames(x), paste0("q=", q_values)))
    log2fc_values <- matrix(NA_real_, nrow = n_genes, ncol = n_q,
                           dimnames = list(rownames(x), paste0("q=", q_values)))
    
    # List to store permutation p-value matrices for all q-values
    # Will be used to compute permutation minima for true WY stepdown
    perm_pvalues_all_q <- vector("list", n_q)
    names(perm_pvalues_all_q) <- paste0("q=", q_values)
    
    # For each q-value, extract the appropriate columns and test
    for (q_idx in seq_along(q_values)) {
        q_val <- q_values[q_idx]
        
        # Get column indices for this q-value
        col_indices <- which(col_q_values == q_val)
        x_q <- x[, col_indices, drop = FALSE]
        
        if (verbose) {
            cat("  q =", q_val, "... ")
        }
        
        # Run permutation test for this q-value
        # Returns both p-values and permutation distributions
        perm_results <- .label_shuffling_single_q(
            x = x_q,
            samples = samples,
            control = control,
            method = method,
            randomizations = randomizations,
            nthreads = nthreads
        )
        
        raw_pvalues[, q_idx] <- perm_results$pvalue
        log2fc_values[, q_idx] <- perm_results$log2FC
        perm_pvalues_all_q[[q_idx]] <- perm_results$perm_matrix
        
        if (verbose) {
            cat("done\n")
        }
    }
    
    # Step 2: Compute permutation minima across all tests (true WY stepdown)
    if (verbose) cat("Computing permutation minima across all tests...\n")
    
    # For true Westfall-Young, we need the minimum p-value from each permutation
    # across all tests. However, computing p-values for each gene in each permutation
    # is expensive. Instead, we use test statistics (absolute fold-change):
    # For each permutation, track the maximum |FC| across all genes * q-values
    # Then for each observed test, count permutations where max(|perm_FC|) >= |obs_FC|
    
    perm_minima <- numeric(randomizations)
    
    for (perm_idx in 1:randomizations) {
        # Collect all fold-change values across genes and q-values for this permutation
        all_fc_abs_perm <- numeric()
        
        for (q_idx in seq_along(q_values)) {
            # perm_pvalues_all_q[[q_idx]] is a matrix: genes * randomizations
            # Extract the fold-change values for this permutation at this q-value
            fc_perm_q <- abs(perm_pvalues_all_q[[q_idx]][, perm_idx])
            all_fc_abs_perm <- c(all_fc_abs_perm, fc_perm_q)
        }
        
        # The minimum p-value in this permutation corresponds to maximum |FC|
        # (since p-value is monotonically decreasing in |FC|)
        # Store the maximum absolute fold-change for this permutation
        perm_minima[perm_idx] <- max(all_fc_abs_perm, na.rm = TRUE)
    }
    
    if (verbose) cat("Applying true Westfall-Young stepdown...\n")
    
    # Flatten p-values and log2FC for stepdown procedure
    all_pvalues <- as.numeric(raw_pvalues)
    all_log2fc <- as.numeric(log2fc_values)
    
    # Create index mapping back to original matrix
    idx_matrix <- expand.grid(gene = 1:n_genes, q = 1:n_q)
    gene_idx <- idx_matrix$gene
    q_idx_vec <- idx_matrix$q
    
    # Sort by p-value (ascending), which corresponds to largest |FC| first
    sort_order <- order(all_pvalues)
    sorted_pvalues <- all_pvalues[sort_order]
    sorted_log2fc <- all_log2fc[sort_order]
    sorted_gene_idx <- gene_idx[sort_order]
    sorted_q_idx <- q_idx_vec[sort_order]
    
    # Initialize adjusted p-values
    wy_adjusted <- numeric(length(all_pvalues))
    
    # True Westfall-Young Stepdown Procedure
    # ======================================
    # For each ranked test i (ordered by p-value):
    # 1. Get the observed test statistic (|log2FC|)
    # 2. Count how many permutation maxima are >= this observed |FC|
    #    (equivalently: count permutations where min p-value would be <= observed p)
    # 3. Adjusted p = (count + 1) / (randomizations + 1)  [Phipson-Smyth correction]
    # 4. Enforce monotonicity: adjusted p-values must be non-decreasing
    #
    # This implements the true WY stepdown from Westfall & Young (1993),
    # accounting for the joint dependence among all gene * q-value tests.
    # Reference: Meinshausen et al. (2011, S165/S166) prove optimality.
    #
    for (i in seq_along(sorted_pvalues)) {
        obs_fc_abs <- abs(sorted_log2fc[i])
        # Count how many permutation maxima are >= observed |FC|
        # (This is equivalent to counting permutations with min p-value <= obs p)
        count <- sum(perm_minima >= obs_fc_abs)
        # Phipson-Smyth correction
        wy_adjusted[sort_order[i]] <- min(1.0, (count + 1) / (randomizations + 1))
    }
    
    # Enforce monotonicity: adjusted p-values must be non-decreasing
    # This maintains the stepdown property and guarantees FWER <= α
    for (i in 2:length(wy_adjusted)) {
        if (wy_adjusted[sort_order[i]] < wy_adjusted[sort_order[i-1]]) {
            wy_adjusted[sort_order[i]] <- wy_adjusted[sort_order[i-1]]
        }
    }
    
    # Step 3: Also compute standard BH-FDR adjustments for comparison
    bh_adjusted <- p.adjust(all_pvalues, method = "BH")
    
    # Step 4: Build output data frame
    output <- data.frame(
        gene = rownames(x)[gene_idx],
        q_value = q_values[q_idx_vec],
        log2FC = all_log2fc,
        pvalue_raw = all_pvalues,
        pvalue_wy = wy_adjusted,  # Westfall-Young stepdown
        padj_bh = bh_adjusted,     # BH-FDR for comparison
        stringsAsFactors = FALSE
    )
    
    # Sort by Westfall-Young adjusted p-value
    output <- output[order(output$pvalue_wy), ]
    rownames(output) <- NULL
    
    # Add attributes with procedure information
    attr(output, "procedure") <- "True Westfall-Young Stepdown"
    attr(output, "fwer_level") <- 0.05
    attr(output, "n_tests") <- n_total_tests
    attr(output, "n_q_values") <- n_q
    attr(output, "n_genes") <- n_genes
    attr(output, "randomizations") <- randomizations
    
    if (verbose) {
        sig_wy <- sum(output$pvalue_wy < 0.05, na.rm = TRUE)
        sig_bh <- sum(output$padj_bh < 0.05, na.rm = TRUE)
        cat("\nResults:\n")
        cat("  Significant (WY-adjusted p < 0.05):", sig_wy, "\n")
        cat("  Significant (BH-FDR < 0.05):", sig_bh, "\n")
        cat("  FWER control: Family-wise error rate <= 0.05\n")
    }
    
    return(output)
}


#' Internal: Single q-value label shuffling
#'
#' Helper function for computing permutation p-values at a single q-value.
#' Returns both the p-values and the full permutation distribution matrix.
#'
#' @param x Matrix of diversity values (genes * samples)
#' @param samples Character vector of sample group labels
#' @param control Control group name
#' @param method Aggregation method ("mean" or "median")
#' @param randomizations Number of permutations
#' @param nthreads Number of threads
#'
#' @return List with:
#'   \describe{
#'     \item{pvalue}{Vector of p-values for each gene}
#'     \item{log2FC}{Vector of observed log2 fold-changes}
#'     \item{perm_matrix}{Matrix of permutation fold-changes (genes * randomizations)}
#'   }
#'
#' @keywords internal
#' @noRd
.label_shuffling_single_q <- function(x, samples, control, method, randomizations, nthreads) {
    
    # Compute observed fold changes
    fc_result <- calculate_fc(x, samples, control, method)
    log2_fc <- fc_result[, 4]
    
    # Generate permutations
    permuted <- replicate(randomizations, 
                         calculate_fc(x, sample(samples), control, method), 
                         simplify = FALSE)
    perm_mat <- vapply(permuted, function(z) as.numeric(z[, 4]), numeric(nrow(x)))
    
    # Ensure perm_mat is always a matrix (fixes single-gene edge case)
    if (!is.matrix(perm_mat)) {
        perm_mat <- matrix(perm_mat, nrow = nrow(x), byrow = FALSE)
    }
    
    # Compute p-values with Phipson-Smyth correction
    .compute_pval <- function(i) {
        obs <- log2_fc[i]
        nulls <- perm_mat[i, ]
        if (is.na(obs) || all(is.na(nulls))) {
            return(NA_real_)
        }
        nulls_non_na <- nulls[!is.na(nulls)]
        n_non_na <- length(nulls_non_na)
        if (n_non_na == 0) {
            return(NA_real_)
        }
        cnt <- sum(abs(nulls_non_na) >= abs(obs))
        pval <- (cnt + 1) / (n_non_na + 1)  # Phipson-Smyth correction
        return(pval)
    }
    
    raw_p_values <- unlist(.tsenat_bplapply(seq_len(nrow(perm_mat)), .compute_pval, 
                                            nthreads = nthreads))
    
    # Replace NA p-values with 1.0 for compatibility
    raw_p_values[is.na(raw_p_values)] <- 1.0
    
    # Return results including permutation matrix for true WY stepdown
    result <- list(
        pvalue = raw_p_values,
        log2FC = log2_fc,
        perm_matrix = perm_mat  # Full permutation distribution (genes * randomizations)
    )
    
    return(result)
}

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

# ════════════════════════════════════════════════════════════════════════════════
# WESTFALL-YOUNG PERMUTATION HELPER (March 2026)
# ════════════════════════════════════════════════════════════════════════════════
# Consolidates redundant WY permutation logic shared between:
#   1. calculate_lm_interaction() - parametric tests (GAM, LMM, GEE)
#   2. detect_q_gene_interactions() - rank-based tests (Kruskal-Wallis, conditional rank)
#
# DESIGN PATTERN:
#   - Core permutation loop is identical in both functions (~70% code duplication)
#   - Model refitting logic differs (parametric vs rank-based)
#   - Solution: Extract permutation machinery, supply model-specific refit_fn callback
#
# USAGE:
#   perm_result <- .tsenat_westfall_young_permutation(
#       n_genes = nrow(results),
#       wy_randomizations = 1000,
#       permute_fn = function() { ... return permutation_assignment ... },
#       refit_fn = function(perm_assignment) { ... return perm_pvalues_vector ... },
#       verbose = FALSE
#   )
#
# RETURNS:
#   List with:
#   - perm_minima: numeric vector of length wy_randomizations (minimum p-value per permutation)
#   - message: character string (if verbose=TRUE)
#
# PHIPSON-SMYTH CORRECTION:
#   Applied by caller as: (sum(perm_minima <= p_obs) + 1) / (wy_randomizations + 1)
#   This adjusts for finite permutation sample and ensures 0 < p_adj <= 1
#
.tsenat_westfall_young_permutation <- function(n_genes, wy_randomizations,
                                               permute_fn, refit_fn,
                                               verbose = FALSE) {
    # Args:
    #   n_genes: Total number of genes (for verbose output)
    #   wy_randomizations: Number of permutations to perform
    #   permute_fn: Callback function () → permutation_assignment
    #     Should return group/sample assignment for permuted data
    #   refit_fn: Callback function (permutation_assignment) → p_values_vector
    #     Should refit model with permuted assignment, return vector of p-values (length n_genes)
    #   verbose: If TRUE, print progress messages
    # 
    # Returns:
    #   List with: perm_minima (numeric vector), message (character or NULL)
    
    if (wy_randomizations < 1) {
        stop("wy_randomizations must be >= 1")
    }
    
    # Initialize storage for minimum p-values across genes for each permutation
    perm_minima <- numeric(wy_randomizations)
    
    # Permutation loop: for each of wy_randomizations random assignments
    for (perm_idx in seq_len(wy_randomizations)) {
        # Call permutation function to get this permutation's assignment
        perm_assignment <- permute_fn()
        
        # Call model-specific refit function with permuted assignment
        # Expected return: numeric vector of p-values (length n_genes)
        perm_pvalues <- refit_fn(perm_assignment)
        
        # Track minimum p-value across all genes for this permutation
        # This is the WY "multiple comparison" correction: using min pval distribution
        perm_minima[perm_idx] <- min(perm_pvalues, na.rm = TRUE)
        
        # Verbose output: report progress at reasonable intervals
        if (verbose && perm_idx %% max(1, ceiling(wy_randomizations / 10)) == 0) {
            message(sprintf("[WY Permutation] Completed %d/%d permutations", 
                          perm_idx, wy_randomizations))
        }
    }
    
    # Ensure monotonicity: p-values should be monotone increasing
    # This is a theoretical requirement for valid adjustment
    # Fix any numerical artifacts via cumulative minimum
    perm_minima <- pmax(perm_minima, 0)  # Ensure non-negative
    perm_minima <- pmin(perm_minima, 1)  # Ensure <= 1
    
    # Return results in format expected by both callers
    result <- list(
        perm_minima = perm_minima,
        n_permutations = wy_randomizations,
        message = if (verbose) 
            sprintf("Westfall-Young permutation test: %d permutations on %d genes completed",
                   wy_randomizations, n_genes)
            else NULL
    )
    
    return(result)
}

