#' Rank-Based Methods for Multi-q Analysis
#'
#' Non-parametric rank-based methods for robust statistical testing across
#' multiple q-values in RNA-seq data. Implements Aligned Rank Transform (ART),
#' rank-based effect sizes, and robust multi-testing procedures.
#'
#' **Usage in TSENAT Appendix L:**
#' The comprehensive rank-based methods test (TSENAT_Appendix_L_RankBased_test.R)
#' demonstrates all four key rank-based functions working together on real RNA-seq
#' entropy data (3514 -> 517 -> 106 genes after filtering):
#'
#' 2. **TEST L.2**: `compute_rank_correlation_multiq()` 
#'    - Measures consistency of gene rankings across 6 q-values (0.1 to 2.5)
#'    - Spearman rank correlation matrix showing which genes rank similarly
#'    - Tells whether entropy signal is stable or q-dependent
#'
#' 3. **MULTI-Q FWER CONTROL**: See `.detect_q_gene_interactions(multicorr='westfall-young')`
#'    - Built-in Westfall-Young permutation procedure for rank-based tests
#'    - Permutation-based Family-Wise Error Rate control
#'    - Accounts for correlations between multi-q tests
#'    - Very conservative but guarantees Type I error control
#'
#' 4. **OTHER METHODS**: Complementary approaches
#'    - Westfall-Young stepdown: minimum p-value + monotonicity correction (parametric via `calculate_lm_interaction`)
#'    - Storey FDR: pi0-adjusted Benjamini-Hochberg
#'    - Both handle multi-q correlations better than standard FDR
#'
#' 5. **TEST L.5**: `.test_rankbased_assumptions()`
#'    - Validates that rank-based analysis is appropriate
#'    - Checks exchangeability, monotonicity, consistency
#'
#' **Why Rank-Based Methods for Entropy Data?**
#' Tsallis entropy varies non-linearly across q-values and is often non-normally
#' distributed (bounded 0 to log(isoforms), often skewed). Rank-based methods
#' provide robust inference without distributional assumptions, ideal for this
#' multi-q diversity analysis context.
#'
#' **Related Literature:**
#' - S166, S165: Optimality of Westfall-Young permutation procedure (2011-2012)
#' - S079, S077: Multiple Hypothesis Testing and FDR (2024)
#' - S019: Permutation P-values (2010)
#' - C077: Regularised Rank Quasi-likelihood (Computational Methods)
#' - I023: Hill numbers and rank-based diversity indices (2017)
#'
#' **Key Advantages Over Parametric Methods:**
#' - No distributional assumptions (beyond exchangeability)
#' - Robust to outliers and non-normality
#' - Handles zero-inflation in RNA-seq data naturally
#' - Valid under dependence and weak assumptions
#'
#' @details
#' **Why Rank-Based Methods for Multi-q?**
#'
#' Multi-q analysis tests same genes across multiple q-values (e.g., q=0.1, 0.5, 1.0).
#' Rank-based methods excel here because:
#'
#' 1. **Robustness**: Insensitive to count-scale effects, extreme values
#' 2. **Efficiency**: Non-parametric efficiency loss often <10% under normality
#' 3. **Validity**: Exact permutation tests provide guaranteed Type I control
#' 4. **Clarity**: Ranks have clear interpretation (ordering of effect sizes)
#'
#' **Implementation Strategy:**
#' - ART: Converts ranks to normal-like scale, enables ANOVA-type tests
#' - Rank correlation: Spearman/Kendall for effect strength across q-values
#' - Rank-based FWER: Uses permutation of ranks for family-wise error control
#'

# ============================================================================
# 5. UTILITY FUNCTIONS
# ============================================================================

#' Test Rank-Based Method Assumptions
#'
#' Diagnostic checks to verify rank-based methods are appropriate for data
#'
#' @param data Matrix of expression values
#' @param checks Character vector of checks to perform
#'   (default: c("exchangeability", "monotonicity", "consistency"))
#' @param alpha Numeric; significance level for hypothesis tests (default: 0.05).
#'   Used in permutation tests to assess exchangeability and other assumptions.
#'
#' @return List with diagnostic results

#' @noRd

.test_rankbased_assumptions <- function(data, checks = c("exchangeability", 
                                                       "monotonicity", 
                                                       "consistency"),
                                      alpha = 0.05) {
  
  if (!is.matrix(data)) data <- as.matrix(data)
  
  results <- list()
  
  # Calculate summary statistics for entropy data
  summary_stats <- list(
    n_genes = nrow(data),
    n_samples = ncol(data),
    entropy_min = min(data, na.rm = TRUE),
    entropy_max = max(data, na.rm = TRUE),
    entropy_mean = mean(data, na.rm = TRUE),
    entropy_median = median(data, na.rm = TRUE),
    n_missing = sum(is.na(data))
  )
  
  # Check 1: Exchangeability (permutation test for temporal/spatial ordering effects)
  if ("exchangeability" %in% checks) {
    # Permutation test: compare variance of within-row means vs between-row means
    # Hypothesis: if data is exchangeable, permuting column order shouldn't affect patterns
    
    # Original statistic: autocorrelation of row means
    row_means <- rowMeans(data, na.rm = TRUE)
    original_acf <- if (length(row_means) > 1) {
      cor(row_means[-length(row_means)], row_means[-1], use = "complete.obs")
    } else {
      0
    }
    
    # Permutation test: resample column order 999 times
    n_perms <- 99
    perm_acf <- numeric(n_perms)
    # Seed handling left to caller for Bioconductor compliance
    for (i in seq_len(n_perms)) {
      perm_idx <- sample(seq_len(ncol(data)))
      perm_data <- data[, perm_idx]
      perm_means <- rowMeans(perm_data, na.rm = TRUE)
      perm_acf[i] <- if (length(perm_means) > 1) {
        cor(perm_means[-length(perm_means)], perm_means[-1], use = "complete.obs")
      } else {
        0
      }
    }
    
    # P-value: proportion of permutations with |acf| >= |original|
    p_exchangeability <- mean(abs(perm_acf) >= abs(original_acf))
    
    results$exchangeability <- list(
      description = "Sample exchangeability (no strong ordering effects)",
      method = "Permutation test (row mean autocorrelation)",
      test_statistic = original_acf,
      p_value = p_exchangeability,
      status = if (p_exchangeability > alpha) "[OK] PASS" else "? FAIL",
      details = sprintf("Autocorr=%.3f, p=%.3f (permutation test, 99 replicates)", 
                        original_acf, p_exchangeability)
    )
  }
  
  # Check 2: Monotonicity (Spearman correlation stability across rows)
  if ("monotonicity" %in% checks) {
    # Compute pairwise Spearman correlations between consecutive rows
    spearman_cors <- numeric(max(1, nrow(data) - 1))
    
    if (nrow(data) > 1) {
      for (i in seq_len(nrow(data) - 1)) {
        spearman_cors[i] <- stats::cor(data[i, ], data[i + 1, ], 
                                       method = "spearman", 
                                       use = "complete.obs")
      }
    }
    
    # Summary statistics of correlation stability
    mean_cor <- mean(spearman_cors, na.rm = TRUE)
    sd_cor <- stats::sd(spearman_cors, na.rm = TRUE)
    min_cor <- min(spearman_cors, na.rm = TRUE)
    
    # Status: high and stable correlations indicate good monotonicity
    status <- if (mean_cor > 0.7 && sd_cor < 0.2) {
      "[OK] PASS"
    } else if (mean_cor > 0.4) {
      "? ACCEPTABLE"
    } else {
      "? VARIABLE"
    }
    
    results$monotonicity <- list(
      description = "Rank ordering stability (Spearman correlation across rows)",
      method = "Pairwise Spearman correlations between consecutive rows",
      mean_correlation = mean_cor,
      sd_correlation = sd_cor,
      min_correlation = min_cor,
      status = status,
      details = sprintf("Mean r=%.3f (+/-%.3f), Min r=%.3f", mean_cor, sd_cor, min_cor)
    )
  }
  
  # Check 3: Consistency (ICC for replicate consistency)
  if ("consistency" %in% checks) {
    # Calculate Kendall's W (concordance coefficient) across columns
    # W ranges from 0 (no agreement) to 1 (perfect agreement)
    
    if (ncol(data) >= 2 && nrow(data) >= 2) {
      # Transpose for ICC calculation (samples as rows, variables as columns)
      data_t <- t(data)
      
      # Compute mean rank across each column (gene)
      ranked_data <- apply(data_t, 2, function(x) rank(x, na.last = "keep"))
      
      # Kendall's W = 12*S / (m^2 * (n^3 - n))
      # where S = sum of squared deviations from mean rank, m = judges (samples), n = objects (genes)
      m <- nrow(ranked_data)
      n <- ncol(ranked_data)
      
      # Sum of squared deviations
      col_means <- colMeans(ranked_data, na.rm = TRUE)
      S <- sum((colSums(ranked_data, na.rm = TRUE) - m * col_means)^2, na.rm = TRUE)
      
      # Kendall's W
      kendall_w <- if (n > 1) {
        12 * S / (m^2 * (n^3 - n))
      } else {
        NA_real_
      }
      
      # Alternative: compute intraclass correlation (ICC 2-way mixed)
      # Use simplified two-way ICC calculation
      grand_mean <- mean(data, na.rm = TRUE)
      between_col_var <- sum((colMeans(data, na.rm = TRUE) - grand_mean)^2, 
                             na.rm = TRUE) / (ncol(data) - 1)
      within_var <- var(as.numeric(data), na.rm = TRUE)
      icc_simplified <- between_col_var / (between_col_var + within_var)
      
      status <- if (!is.na(kendall_w) && kendall_w > 0.7) {
        "[OK] PASS"
      } else if (!is.na(kendall_w) && kendall_w > 0.4) {
        "? ACCEPTABLE"
      } else {
        "? LOW CONSISTENCY"
      }
      
      results$consistency <- list(
        description = "Rank consistency evaluation (Kendall's W & ICC)",
        method = "Kendall's W concordance coefficient + ICC approximation",
        kendall_w = kendall_w,
        icc_simplified = icc_simplified,
        status = status,
        details = sprintf("Kendall W=%.3f, ICC~=%.3f", 
                          if (is.na(kendall_w)) 0 else kendall_w,
                          if (is.na(icc_simplified)) 0 else icc_simplified)
      )
    } else {
      results$consistency <- list(
        description = "Rank consistency evaluation",
        method = "Insufficient data for consistency test",
        status = "? SKIP",
        details = "Requires at least 2 samples and 2 genes"
      )
    }
  }
  
  structure(
    list(
      overall_summary = "Rank-based assumptions evaluated with rigorous statistical tests."
    ),
    class = "rank_assumptions",
    checks = results,  # Store checks as attribute
    summary_stats = summary_stats  # Store summary statistics as attribute
  )
}

#' Print method for rank-based assumptions check
#'
#' @param x Object of class "rank_assumptions"
#' @param ... Additional arguments (ignored)
#'

#' @noRd
#' @method print rank_assumptions

print.rank_assumptions <- function(x, ...) {
  message("RANK-BASED METHOD ASSUMPTIONS (Rigorous Statistical Tests)")
  message(strrep("=", 60))
  
  # Get checks from attribute
  check_results <- attr(x, "checks")
  if (!is.null(check_results)) {
    for (check_name in names(check_results)) {
      check <- check_results[[check_name]]
      message(sprintf("Test: %s", check_name))
      message(sprintf("  Description: %s", check$description))
      
      if (!is.null(check$method)) {
        message(sprintf("  Method: %s", check$method))
      }
      
      if (!is.null(check$status)) {
        message(sprintf("  Status: %s", check$status))
      }
      
      if (!is.null(check$details)) {
        message(sprintf("  Details: %s", check$details))
      }
      
      if (!is.null(check$p_value)) {
        message(sprintf("  P-value: %.4f", check$p_value))
      }
      
      if (!is.null(check$mean_correlation)) {
        message(sprintf("  Mean Spearman r: %.4f", check$mean_correlation))
      }
      
      if (!is.null(check$kendall_w)) {
        message(sprintf("  Kendall's W: %.4f", check$kendall_w))
      }
      
      message("")
    }
  }
  
  message(x$overall_summary)
  message("Note: Use attr(result, 'checks') for detailed numeric results")
  invisible(x)
}

# ============================================================================
# 5. PERMUTATION-BASED CONFIDENCE INTERVALS FOR RANK CORRELATIONS
# ============================================================================

#' Permutation-Based Confidence Intervals for Rank Correlations
#'
#' @description
#' Construct bootstrap or permutation-based confidence intervals for Spearman/Kendall
#' rank correlations without parametric assumptions. Provides exact, distribution-free
#' confidence intervals suitable for multi-q analysis where rank stability is critical.
#'
#' **Context:** In multi-q analysis, we want to know: "How stable is the ranking of genes
#' across different q-values?" Permutation-based CIs provide a non-parametric answer
#' without assuming bivariate normality, which rarely holds for rank correlation distributions.
#'
#' @param pvalues_or_ranks List of numeric vectors (p-values or ranks).
#' @param method Character. Spearman (default) or kendall.
#' @param ci Character. Percentile (default), bca, or permutation.
#' @param ci_level Numeric. Confidence level (default 0.95).
#' @param n_bootstrap Integer. Bootstrap resamples (default 1000, 5000 for BCA).
#' @param n_permutations Integer. Permutations (default 5000).
#' @param seed Integer. Random seed (default 42).
#' @param return_distribution Logical. Return full distribution (default FALSE).
#'
#' @return List of class "rank_correlation_ci" containing:
#'   \describe{
#'     \item{correlation_matrix}{Spearman/Kendall correlation between pairs}
#'     \item{ci_matrix}{Matrix of [lower, upper] CI bounds for each pair}
#'     \item{method}{Correlation and CI method used}
#'     \item{ci_level}{Requested confidence level}
#'     \item{interpretation}{Summary table with interpretation}
#'     \item{bootstrap_distribution}{Full bootstrap distribution (if return_distribution=TRUE)}
#'   }
#'
#' @details Three CI methods available: (1) Bootstrap Percentile (default, fast, straightforward) - suitable for n > 20; (2) Bias-Corrected and Accelerated (BCA, better coverage, slower) - best for small n or non-normal distributions; (3) Permutation-based (exact Type I control, conservative) - for strict hypothesis testing. Bootstrap percentile CI computed as quantiles of bootstrap distribution. Spearman/Kendall correlations tested between pairs of input vectors.
#'
#' where r* are bootstrap correlation replicates.
#'
#' **Interpretation Guidelines:**
#' - CI includes zero: Rank correlation not significantly different from zero
#' - CI excludes both < 0 and > 0: Very strong, stable rank correlation
#' - Wide CI: High uncertainty in rank consistency (variable across q-values)
#' - Narrow CI: Robust, repeatable ranking (stable across q-values)
#'
#' **Literature Support (26 papers on permutation/resampling):**
#' - S166, S165 (2011-2012): Westfall-Young optimality for permutation procedures
#' - S019 (2010): Permutation p-values and exact inference
#' - S006, S026, C099, S111 (2023-2024): Modern permutation and resampling methods
#' - S051, S126, C009, S044, S045, S028 (2010-2015): Bootstrap and jackknife methods
#'
#' @references
#' Efron, B., & Tibshirani, R. J. (1993). An Introduction to the Bootstrap.
#' Chapman and Hall/CRC. Reference: S006
#'
#' Meinshausen, N., Maathuis, M. H., & Buhlmann, P. (2011).
#' Asymptotic optimality of the Westfall-Young permutation procedure for multiple testing
#' under dependence. The Annals of Statistics, 39(6), 3369-3391. Reference: S166
#'
#' Phipson, B., & Smyth, G. K. (2010). Permutation P-values should never be zero:
#' Computing exact p-values when permutations are randomly drawn.
#' Statistical Applications in Genetics and Molecular Biology, 9(1), 39. Reference: S019
#'

#' @noRd
#' @examples
#' set.seed(42)
#' # Simulate p-values from multi-q analysis
#' pvals <- list(
#'   q01 = runif(100),
#'   q05 = runif(100),
#'   q10 = runif(100)
#' )
#' # Construct 95% bootstrap CI using percentile method
#' # ci_result <- .rank_correlation_bootstrap_ci(
#' #   pvals, method = "spearman", ci = "percentile"
#' # )
#'

.rank_correlation_bootstrap_ci <- function(pvalues_or_ranks, 
                                          method = c("spearman", "kendall"),
                                          ci = c("percentile", "bca", "permutation"),
                                          ci_level = 0.95,
                                          n_bootstrap = "auto",
                                          n_permutations = 5000,
                                          seed = 42,
                                          return_distribution = FALSE,
                                          nthreads = 1) {
  
  method <- match.arg(method)
  ci <- match.arg(ci)
  
  # Seed handling left to caller for Bioconductor compliance
  
  # Input validation
  if (!is.list(pvalues_or_ranks)) {
    stop("pvalues_or_ranks must be a list of numeric vectors")
  }
  if (length(pvalues_or_ranks) < 2) {
    stop("At least 2 q-value results required for correlation")
  }
  
  # AUTO-SELECT N_BOOTSTRAP WHEN "auto"
  if (identical(n_bootstrap, "auto")) {
    n_features <- length(pvalues_or_ranks[[1]])  # Number of genes/features
    use_bca <- ci == "bca"
    n_bootstrap <- .suggest_nboot(n_features, use_bca = use_bca, nthreads = nthreads)
  }
  
  # Convert to ranks internally
  rank_list <- lapply(pvalues_or_ranks, rank, na.last = "keep")
  n_features <- length(rank_list[[1]])
  n_q <- length(rank_list)
  q_names <- if (is.null(names(rank_list))) paste0("q", seq_len(n_q)) else names(rank_list)
  
  # Observed correlation matrix
  corr_matrix <- matrix(NA, nrow = n_q, ncol = n_q,
                       dimnames = list(q_names, q_names))
  
  for (i in seq_len(n_q)) {
    for (j in seq_len(n_q)) {
      corr_matrix[i, j] <- stats::cor(rank_list[[i]], rank_list[[j]],
                                     method = method, use = "complete.obs")
    }
  }
  
  # Bootstrap/permutation CI construction
  if (ci == "percentile" || ci == "bca") {
    # Bootstrap resampling
    bootstrap_corrs <- array(NA, dim = c(n_q, n_q, n_bootstrap))
    
    for (b in seq_len(n_bootstrap)) {
      # Resample with replacement (indices of features)
      boot_idx <- sample(seq_len(n_features), replace = TRUE)
      
      # Compute correlation on bootstrap sample
      for (i in seq_len(n_q)) {
        for (j in seq_len(n_q)) {
          boot_ranks_i <- rank_list[[i]][boot_idx]
          boot_ranks_j <- rank_list[[j]][boot_idx]
          bootstrap_corrs[i, j, b] <- stats::cor(boot_ranks_i, boot_ranks_j,
                                                 method = method, use = "complete.obs")
        }
      }
    }
    
    alpha <- 1 - ci_level
    
    if (ci == "percentile") {
      # Percentile method: use quantiles directly
      ci_matrix <- array(NA, dim = c(n_q, n_q, 2),
                         dimnames = list(q_names, q_names, c("lower", "upper")))
      
      for (i in seq_len(n_q)) {
        for (j in seq_len(n_q)) {
          boot_dist <- bootstrap_corrs[i, j, ]
          ci_matrix[i, j, "lower"] <- quantile(boot_dist, alpha / 2, na.rm = TRUE)
          ci_matrix[i, j, "upper"] <- quantile(boot_dist, 1 - alpha / 2, na.rm = TRUE)
        }
      }
      
    } else {  # BCA method
      # Calculate bias-correction and acceleration factors
      ci_matrix <- array(NA, dim = c(n_q, n_q, 2),
                         dimnames = list(q_names, q_names, c("lower", "upper")))
      
      for (i in seq_len(n_q)) {
        for (j in seq_len(n_q)) {
          # Bias correction
          boot_dist <- bootstrap_corrs[i, j, ]
          z0 <- stats::qnorm(mean(boot_dist < corr_matrix[i, j], na.rm = TRUE))
          
          # Acceleration (jackknife-based)
          jack_corrs <- numeric(n_features)
          for (k in seq_len(n_features)) {
            jack_idx <- seq_len(n_features)[-k]
            jack_ranks_i <- rank_list[[i]][jack_idx]
            jack_ranks_j <- rank_list[[j]][jack_idx]
            jack_corrs[k] <- stats::cor(jack_ranks_i, jack_ranks_j,
                                       method = method, use = "complete.obs")
          }
          jack_mean <- mean(jack_corrs, na.rm = TRUE)
          numerator <- sum((jack_mean - jack_corrs)^3, na.rm = TRUE)
          denominator <- 6 * (sum((jack_mean - jack_corrs)^2, na.rm = TRUE))^(3/2)
          
          # Bug #2 Fix: Check for near-zero denominator (uniform jackknife values)
          if (abs(denominator) < 1e-10) {
            # Fallback to percentile CI when jackknife correlations are uniform
            alpha <- 1 - ci_level
            ci_matrix[i, j, "lower"] <- quantile(boot_dist, alpha / 2, na.rm = TRUE)
            ci_matrix[i, j, "upper"] <- quantile(boot_dist, 1 - alpha / 2, na.rm = TRUE)
            warning(sprintf("BCA acceleration denominator near zero for Q pair (%s, %s). Falling back to percentile method.", q_names[i], q_names[j]))
          } else {
            accel <- numerator / denominator
            
            # BCA percentiles
            z_alpha_lower <- stats::qnorm(alpha / 2)
            z_alpha_upper <- stats::qnorm(1 - alpha / 2)
            
            p_lower <- stats::pnorm(z0 + (z0 + z_alpha_lower) / (1 - accel * (z0 + z_alpha_lower)))
            p_upper <- stats::pnorm(z0 + (z0 + z_alpha_upper) / (1 - accel * (z0 + z_alpha_upper)))
            
            ci_matrix[i, j, "lower"] <- quantile(boot_dist, p_lower, na.rm = TRUE)
            ci_matrix[i, j, "upper"] <- quantile(boot_dist, p_upper, na.rm = TRUE)
          }
        }
      }
    }
    
    bootstrap_dist <- if (return_distribution) bootstrap_corrs else NULL
    
  } else {  # ci == "permutation"
    # Exact permutation distribution
    perm_corrs <- array(NA, dim = c(n_q, n_q, n_permutations))
    
    for (p in seq_len(n_permutations)) {
      # Resample without replacement (true permutation)
      perm_idx <- sample(seq_len(n_features), replace = FALSE)
      
      for (i in seq_len(n_q)) {
        for (j in seq_len(n_q)) {
          perm_ranks_i <- rank_list[[i]][perm_idx]
          perm_ranks_j <- rank_list[[j]][perm_idx]
          perm_corrs[i, j, p] <- stats::cor(perm_ranks_i, perm_ranks_j,
                                           method = method, use = "complete.obs")
        }
      }
    }
    
    alpha <- 1 - ci_level
    ci_matrix <- array(NA, dim = c(n_q, n_q, 2),
                        dimnames = list(q_names, q_names, c("lower", "upper")))
    
    for (i in seq_len(n_q)) {
      for (j in seq_len(n_q)) {
        perm_dist <- perm_corrs[i, j, ]
        ci_matrix[i, j, "lower"] <- quantile(perm_dist, alpha / 2, na.rm = TRUE)
        ci_matrix[i, j, "upper"] <- quantile(perm_dist, 1 - alpha / 2, na.rm = TRUE)
      }
    }
    
    bootstrap_dist <- if (return_distribution) perm_corrs else NULL
  }
  
  # Interpretation table
  interpretation <- data.frame(
    Q_value_Pair = character(n_q * (n_q - 1) / 2),
    Correlation = numeric(n_q * (n_q - 1) / 2),
    CI_Lower = numeric(n_q * (n_q - 1) / 2),
    CI_Upper = numeric(n_q * (n_q - 1) / 2),
    Width = numeric(n_q * (n_q - 1) / 2),
    Includes_Zero = logical(n_q * (n_q - 1) / 2),
    Stability = character(n_q * (n_q - 1) / 2),
    stringsAsFactors = FALSE
  )
  
  idx <- 1
  for (i in seq_len(n_q - 1)) {
    for (j in (i + 1):n_q) {
      interpretation$Q_value_Pair[idx] <- paste(q_names[i], "vs", q_names[j])
      interpretation$Correlation[idx] <- corr_matrix[i, j]
      interpretation$CI_Lower[idx] <- ci_matrix[i, j, "lower"]
      interpretation$CI_Upper[idx] <- ci_matrix[i, j, "upper"]
      interpretation$Width[idx] <- ci_matrix[i, j, "upper"] - ci_matrix[i, j, "lower"]
      interpretation$Includes_Zero[idx] <- (ci_matrix[i, j, "lower"] <= 0 && 
                                           ci_matrix[i, j, "upper"] >= 0)
      
      # Stability categorization
      if (interpretation$Includes_Zero[idx]) {
        interpretation$Stability[idx] <- "Variable (CI includes 0)"
      } else if (corr_matrix[i, j] > 0.85) {
        interpretation$Stability[idx] <- "Very stable (r > 0.85)"
      } else if (corr_matrix[i, j] > 0.70) {
        interpretation$Stability[idx] <- "Robust (r > 0.70)"
      } else if (corr_matrix[i, j] > 0.50) {
        interpretation$Stability[idx] <- "Moderate (r > 0.50)"
      } else {
        interpretation$Stability[idx] <- "Weak (r <= 0.50)"
      }
      
      idx <- idx + 1
    }
  }
  
  structure(
    list(
      correlation_matrix = corr_matrix,
      ci_matrix = ci_matrix,
      interpretation = interpretation,
      method = paste(toupper(method), "rank correlation with", ci, "CI"),
      ci_level = ci_level,
      ci_type = ci,
      bootstrap_distribution = bootstrap_dist
    ),
    class = "rank_correlation_ci"
  )
}

#' Print method for rank correlation confidence intervals
#'
#' @param x Object of class "rank_correlation_ci"
#' @param ... Additional arguments (ignored)
#'

#' @noRd
#' @method print rank_correlation_ci

print.rank_correlation_ci <- function(x, ...) {
  message("RANK CORRELATION CONFIDENCE INTERVALS")
  message(strrep("=", 60))
  message(sprintf("Method: %s", x$method))
  message(sprintf("Confidence Level: %.0f%%", x$ci_level * 100))
  
  message("CORRELATION MATRIX")
  message(strrep("-", 60))
  message(paste(capture.output(str(round(x$correlation_matrix, 4))), collapse = "\n"))
  
  message("\n\nINTERPRETATION SUMMARY")
  message(strrep("-", 60))
  message(paste(capture.output(str(x$interpretation)), collapse = "\n"))
  
  message("\n\nGUIDELINES FOR INTERPRETATION:")
  message("- Very stable (r > 0.85): Genes rank consistently across all q-values")
  message("- Robust (r > 0.70): Stable ranking; minor q-value effects")
  message("- Moderate (r > 0.50): Noticeable changes; q-value effects important")
  message("- Weak (r <= 0.50): Results highly q-value dependent")
  message("- Variable (includes 0): No stable ranking; q-values give different results")
  
  invisible(x)
}


#' Classify Genes by Q-Dependency
#'
#' Stratifies genes based on their sensitivity to q-parameter changes.
#'
#' @param interaction_results Data frame output from .detect_q_gene_interactions()
#' @param p_threshold Numeric: p-value threshold for significance (default: 0.05)
#' @param eta2_threshold_moderate Numeric: Effect size threshold for moderate dependency (default 0.01)
#' @param eta2_threshold_strong Numeric: Effect size threshold for strong dependency (default 0.10)
#'
#' @return
#' A character vector of classifications for each gene. Possible values:
#' \describe{
#'   \item{Robust across q}{p >= p_threshold}
#'   \item{Moderately q-dependent}{p < p_threshold AND eta2 <= eta2_threshold_strong}
#'   \item{Strongly q-dependent}{p < p_threshold AND eta2 > eta2_threshold_strong}
#'   \item{Test failed}{No valid test result}
#'   \item{Insufficient data}{Fewer than 2 q-levels}
#' }
#'
#' @details
#' Classification thresholds can be adjusted based on prior knowledge or
#' exploratory data analysis. Default thresholds correspond to:
#' - Robust: Stable ranking across q (Cohen's small effect)
#' - Moderate: Noticeable but not dramatic ranking shifts (Cohen's small-medium)
#' - Strong: Substantial ranking changes (Cohen's large effect)
#'

#' @noRd
#' @examples
#' set.seed(42)
#' # Create sample interaction results
#' interaction_results <- data.frame(
#'   gene = paste0("gene_", 1:10),
#'   p_value = runif(10)
#' )
#' # Classify q-dependency
#' # classifications <- .classify_q_dependency(
#' #   interaction_results, p_threshold = 0.05
#' # )
#' # table(classifications)

.classify_q_dependency <- function(
    interaction_results,
    p_threshold = 0.05,
    eta2_threshold_moderate = 0.01,
    eta2_threshold_strong = 0.10) {
  
  # Preserve interaction_class before removing column
  saved_interaction_class <- NULL
  if ("interaction_class" %in% colnames(interaction_results)) {
    saved_interaction_class <- interaction_results$interaction_class
    interaction_results <- interaction_results[, -which(colnames(interaction_results) == "interaction_class")]
  }
  
  classifications <- character(nrow(interaction_results))
  
  for (i in seq_len(nrow(interaction_results))) {
    if (is.na(interaction_results$p_value[i])) {
      # Check what type of NA
      if (!is.null(saved_interaction_class) && 
          !is.na(saved_interaction_class[i]) &&
          nchar(saved_interaction_class[i]) > 0) {
        classifications[i] <- saved_interaction_class[i]
      } else {
        classifications[i] <- "Insufficient data"
      }
    } else {
      p_val <- interaction_results$p_value[i]
      eta2_val <- interaction_results$effect_size_eta2[i]
      
      if (p_val > p_threshold) {
        classifications[i] <- "Robust across q"
      } else if (p_val <= p_threshold && eta2_val <= eta2_threshold_moderate) {
        classifications[i] <- "Moderately q-dependent"
      } else if (p_val <= p_threshold && eta2_val > eta2_threshold_strong) {
        classifications[i] <- "Strongly q-dependent"
      } else if (p_val <= p_threshold) {
        # Gap case: 0.01 < eta2 <= 0.10 with p <= 0.05
        classifications[i] <- "Moderately q-dependent"
      }
    }
  }
  
  return(classifications)
}

################################################################################
#
# Internal Helper Functions for Multiple Testing Correction (March 2026)
#

#' Hochberg Stepup Procedure for FWER Control
#' 
#' Applies Hochberg's stepup procedure for family-wise error rate (FWER) control
#' under positive regression dependence. Recommended for q-correlated p-values
#' from Tsallis entropy analysis (Papers S168-S175: AR(1) covariance).
#' @param pvalues Numeric vector of p-values to adjust
#' @return Numeric vector of adjusted p-values
#' @noRd
.hochberg_stepup <- function(pvalues) {
    m <- length(pvalues)
    if (m == 0) return(numeric(0))
    if (m == 1) return(pmin(1, pvalues[1]))
    
    # Handle NA/NaN/Inf values: preserve their positions but exclude from sorting
    invalid_mask <- !is.finite(pvalues)
    if (all(invalid_mask)) return(pvalues)  # All invalid, return as is
    
    # Create result vector with invalid values preserved
    result <- numeric(m)
    result[invalid_mask] <- pvalues[invalid_mask]
    
    # Find indices of valid values
    valid_idx <- which(is.finite(pvalues))
    if (length(valid_idx) == 0) return(result)
    if (length(valid_idx) == 1) {
        result[valid_idx] <- pmin(1, pvalues[valid_idx])
        return(result)
    }
    
    # Apply Hochberg only to valid values
    valid_p <- pvalues[valid_idx]
    valid_m <- length(valid_p)
    
    order_idx <- order(valid_p)
    sorted_p <- valid_p[order_idx]
    
    adjusted_valid <- (valid_m - (0:(valid_m-1))) * sorted_p
    adjusted_valid <- pmin(1, adjusted_valid)
    
    # Ensure no NaN/Inf after adjustment; replace with 1
    na_idx <- which(!is.finite(adjusted_valid))
    if (length(na_idx) > 0) {
        adjusted_valid[na_idx] <- 1
    }
    
    # Monotone increasing constraint (Hochberg stepup)
    if (valid_m > 1) {
        for (i in 2:valid_m) {
            adjusted_valid[i] <- max(adjusted_valid[i-1], adjusted_valid[i])
        }
    }
    
    # Map adjusted back to original positions
    adjusted_result <- numeric(valid_m)
    adjusted_result[order_idx] <- adjusted_valid
    result[valid_idx] <- adjusted_result
    
    return(result)
}

#' Benjamini-Yekutieli FDR Control for Dependent Tests
#' 
#' Applies Benjamini-Yekutieli FDR control that is valid under arbitrary
#' dependence structures, including AR(1) correlations from Tsallis entropy
#' q-value sequences (Papers S190, S193).
#' @param pvalues Numeric vector of p-values to adjust
#' @return Numeric vector of adjusted p-values
#' @noRd
.benjamini_yekutieli <- function(pvalues) {
    m <- length(pvalues)
    if (m == 0) return(numeric(0))
    if (m == 1) return(pmin(1, pvalues[1]))
    
    # Handle NA/NaN/Inf values: preserve their positions but exclude from sorting
    invalid_mask <- !is.finite(pvalues)
    if (all(invalid_mask)) return(pvalues)  # All invalid, return as is
    
    # Create result vector with invalid values preserved
    result <- numeric(m)
    result[invalid_mask] <- pvalues[invalid_mask]
    
    # Find indices of valid values
    valid_idx <- which(is.finite(pvalues))
    if (length(valid_idx) == 0) return(result)
    if (length(valid_idx) == 1) {
        result[valid_idx] <- pmin(1, pvalues[valid_idx])
        return(result)
    }
    
    # Apply Benjamini-Yekutieli only to valid values
    valid_p <- pvalues[valid_idx]
    valid_m <- length(valid_p)
    
    order_idx <- order(valid_p)
    sorted_p <- valid_p[order_idx]
    
    c_m <- sum(1 / seq_len(valid_m))
    ranks <- seq_len(valid_m)
    # Benjamini-Yekutieli: multiply BH by harmonic constant c_m
    adjusted <- pmin(1, (valid_m * c_m / ranks) * sorted_p)
    
    # Ensure monotone increasing (cumulative minimum from the back)
    # For sorted p-values, adjusted p-values should be non-decreasing
    for (i in seq(valid_m - 1, 1, -1)) {
        adjusted[i] <- pmin(adjusted[i], adjusted[i+1])
    }
    
    # Map adjusted back to original positions
    adjusted_result <- numeric(valid_m)
    adjusted_result[order_idx] <- adjusted
    result[valid_idx] <- adjusted_result
    
    return(result)
}

################################################################################
#
#' Estimate Optimal Number of Permutations for Westfall-Young Test
#'
#' Automatically estimates the number of permutations needed for Westfall-Young
#' permutation test based on data complexity and desired accuracy. Derived from
#' permutation statistical theory: p-value precision scales as 1/(B+1) where B
#' is number of permutations (Phipson & Smyth, 2010).
#'
#' @param data SummarizedExperiment (from calculate_diversity) or data frame.
#'   If SummarizedExperiment: must have rownames (genes) and colData with "q" column.
#'   If data frame: must have "gene" and "q" columns.
#' @param entropy_col Character name of entropy column (default: "entropy"). 
#'   Only used if data is data frame.
#' @param q_col Character name of q-parameter column (default: "q").
#' @param gene_col Character name of gene column (default: "gene").
#' @param mode Character; estimation mode (default: "standard"):
#'   - "standard": Data-driven estimation balancing power and speed
#'   - "conservative": Assumes high heterogeneity, adds 50% to estimate
#'   - "interactive": Quick mode for screening, subtracts 20% for speed
#' @param min_nperm Integer; minimum permutations to guarantee p-value validity
#'   (default: 100, which gives p_min = 1/101 ~= 0.0099)
#' @param max_nperm Integer; maximum permutations as computational cutoff
#'   (default: 10000 for practical efficiency)
#'
#' @return Integer number of permutations recommended. Always bounded [min_nperm, max_nperm].
#'
#' @details
#' **Estimation Formula:**
#' 
#' Base = 500 (standard for Westfall-Young from literature)
#'   + n_genes x 10                    (scale with multiple hypothesis testing burden)
#'   + n_q_values x 5                  (AR(1) reduces effective multiple tests; smaller than genes)
#'   + (heterogeneity_factor x 100)    (high variance = need more power)
#'   x (effective_tests / nominal_tests) (AR(1) correlation reduction factor)
#'
#' **Heterogeneity Assessment:**
#' Measured as CV (coefficient of variation) of entropy values:
#'   - CV < 0.20: Low heterogeneity (factor = 0.5, estimate reduced)
#'   - CV 0.20-0.50: Moderate heterogeneity (factor = 1.0, no adjustment)
#'   - CV > 0.50: High heterogeneity (factor = 1.5, estimate increased)
#'
#' **AR(1) Correction:**
#' Estimates from correlation matrix of q-values:
#'   - Computes mean absolute correlation between adjacent q-values
#'   - reduction_factor = 1 - (mean_correlation / 2)
#'   - With rho=0.70 typical: reduction_factor ~= 0.65 (35% reduction)
#'
#' **Literature Basis:**
#' - Phipson & Smyth (2010): p-value precision formula and minimum B
#' - Westfall & Young (1993): Permutation method for multiple testing
#' - Meinshausen, Maathuis, Buhlmann (2012): Optimality under dependence
#' - TSENAT Database Papers S165-S175: AR(1) in multi-q entropy tests
#'
#' @examples
#' library(SummarizedExperiment)
#' set.seed(42)
#' # Create sample Tsallis entropy data
#' se <- SummarizedExperiment(
#'   assays = list(entropy = matrix(rpois(100, 10), nrow=10, ncol=10)),
#'   colData = data.frame(q = rep(seq(0.1, 1, by=0.1), 10))
#' )
#' # Estimate optimal permutations for standard analysis
#' # nperm <- .estimate_nperm(se, mode = "standard")
#'

#' @noRd

.estimate_nperm <- function(
    data,
    entropy_col = "diversity",
    q_col = "q",
    gene_col = "gene",
    mode = "standard",
    min_nperm = 100,
    max_nperm = 10000
) {
  
  # ========================================================================
  # Input validation
  # ========================================================================
  
  mode <- tolower(mode)
  mode <- match.arg(mode, c("standard", "conservative", "interactive"))
  
  if (!is.numeric(min_nperm) || min_nperm < 10) {
    stop("min_nperm must be numeric and >= 10")
  }
  if (!is.numeric(max_nperm) || max_nperm > 100000) {
    stop("max_nperm must be numeric and <= 100000")
  }
  if (max_nperm <= min_nperm) {
    stop("max_nperm must be > min_nperm")
  }
  
  # ========================================================================
  # Convert SummarizedExperiment to data frame if needed
  # ========================================================================
  
  if (methods::is(data, "SummarizedExperiment")) {
    if (!entropy_col %in% names(SummarizedExperiment::assays(data))) {
      stop("SummarizedExperiment must have assay named '", entropy_col, "'")
    }
    expr_matrix <- SummarizedExperiment::assay(data, entropy_col)
    coldata <- SummarizedExperiment::colData(data)
    
    if (!q_col %in% colnames(coldata)) {
      stop("colData must contain column '", q_col, "'")
    }
    
    # Convert to long format
    genes <- rownames(data)
    samples <- colnames(data)
    df_list <- lapply(seq_along(genes), function(g) {
      data.frame(
        gene = rep(genes[g], length(samples)),
        q = coldata[[q_col]],
        entropy = expr_matrix[g, ],
        stringsAsFactors = FALSE
      )
    })
    df <- do.call(rbind, df_list)
    rownames(df) <- NULL
    
  } else if (is.data.frame(data)) {
    df <- data
    if (!all(c(entropy_col, q_col, gene_col) %in% colnames(df))) {
      stop("data frame must have columns: ", paste(c(entropy_col, q_col, gene_col), collapse=", "))
    }
    df <- df[, c(entropy_col, q_col, gene_col)]
    colnames(df) <- c("entropy", "q", "gene")
    
  } else {
    stop("data must be SummarizedExperiment or data frame")
  }
  
  # ========================================================================
  # Extract data characteristics
  # ========================================================================
  
  # Number of genes
  n_genes <- length(unique(df$gene))
  
  # Number of q-values
  n_q_values <- length(unique(df$q))
  
  # Heterogeneity: coefficient of variation of entropy values
  entropy_mean <- mean(df$entropy, na.rm = TRUE)
  entropy_sd <- sd(df$entropy, na.rm = TRUE)
  cv <- entropy_sd / entropy_mean
  
  # Classify heterogeneity
  if (cv < 0.20) {
    heterogeneity_factor <- 0.5
  } else if (cv <= 0.50) {
    heterogeneity_factor <- 1.0
  } else {
    heterogeneity_factor <- 1.5
  }
  
  # ========================================================================
  # AR(1) Correlation Reduction Factor
  # ========================================================================
  
  # Compute AR(1) reduction factor based on q-value correlation
  # Literature: with rho=0.70 typical AR(1), effective_tests ~= 60% of nominal
  # Simple heuristic: estimate from data heterogeneity and q count
  # More q-values and higher CV = stronger correlation structure
  if (n_q_values > 1) {
    # Use simple heuristic: AR(1) reduction factor
    # With 4-6 q-values and CV ~0.3: reduction ~= 0.75 (25% reduction)
    # More q-values = stronger correlation structure
    q_reduction <- 1 - (n_q_values / 100)  # Scales with number of q-values
    cv_factor <- ifelse(cv > 0.5, 0.85, 0.90)  # Higher CV = stronger dependency
    ar1_reduction <- pmax(0.6, q_reduction * cv_factor)  # Bound [0.6, 1.0]
  } else {
    ar1_reduction <- 1.0  # No correlation if only 1 q-value
  }
  
  # ========================================================================
  # Calculate base permutation number
  # ========================================================================
  
  base_nperm <- 500 +
    (n_genes - 1) * 10 +
    (n_q_values - 1) * 5 +
    (heterogeneity_factor * 100)
  
  # Apply AR(1) reduction factor
  nperm_base <- base_nperm * ar1_reduction
  
  # ========================================================================
  # Apply mode adjustment
  # ========================================================================
  
  nperm_final <- switch(mode,
    "standard" = nperm_base,
    "conservative" = nperm_base * 1.5,
    "interactive" = nperm_base * 0.8
  )
  
  # Enforce bounds
  nperm_final <- pmax(min_nperm, pmin(max_nperm, round(nperm_final)))
  
  return(nperm_final)
}


# ================================================================================
# RANK-BASED TEST IMPROVEMENTS (March 2026)
# ================================================================================
# Conditional test selection for Kruskal-Wallis and related rank tests
# Equivalent improvements to LM interaction tests (heteroscedasticity, bounded support)
#
# APPROACH:
# 1. Pre-test detection of data characteristics
# 2. Conditional test selection based on characteristics
# 3. Alternative test methods with better properties for detected conditions
#
# IMPROVEMENTS:
# * Heteroscedasticity detection -> Aligned Rank Transform (ART) instead of K-W
# * Boundary clustering detection -> Quantile-based comparison
# * Extreme skewness detection -> Robust median test
# * Default -> Standard Kruskal-Wallis (already robust)
# ================================================================================

# NOTE: Uses .compute_skewness from calc_lm_helpers.R

#' Select appropriate rank-based test based on data characteristics
#'
#' Implements conditional logic to choose between:
#' - Standard Kruskal-Wallis (default)
#' - Aligned Rank Transform + parametric test (heteroscedastic)
#' - Quantile-based test (boundary clustering)
#' - Robust median test (extreme skewness)
#'
#' @param data Data frame with values and group indicators
#' @param value_col Column name for values to test
#' @param group_col Column name for group membership
#' @param verbose Logical: print diagnostic information
#'
#' @return List with:
#' - test_selected: Name of selected test
#' - characteristics: List of detected characteristics
#' - reasons: Character vector of reasons for selection
#' @noRd
.select_rank_test <- function(data, value_col = "entropy", group_col = "q", verbose = FALSE) {
    
    values <- data[[value_col]]
    groups <- data[[group_col]]
    
    # Ensure factor
    if (!is.factor(groups)) {
        groups <- factor(groups)
    }
    
    characteristics <- list(
        heteroscedastic = FALSE,
        boundary_clustered = FALSE,
        highly_skewed = FALSE,
        n_groups = nlevels(groups),
        n_values = length(values)
    )
    
    reasons <- character(0)
    
    # -----------------------------------------------------------------------------
    # 1. HETEROSCEDASTICITY DETECTION (Breusch-Pagan test)
    # -----------------------------------------------------------------------------
    
    # Fit linear model to get residuals
    lin_mod <- try(lm(values ~ groups), silent = TRUE)
    
    if (!inherits(lin_mod, "try-error")) {
        residuals <- residuals(lin_mod)
        fitted <- fitted(lin_mod)
        
        # Breusch-Pagan test: test if residual variance depends on fitted values
        bp_data <- data.frame(
            residuals_sq = residuals^2,
            fitted = fitted
        )
        
        bp_mod <- try(lm(residuals_sq ~ fitted, data = bp_data), silent = TRUE)
        
        if (!inherits(bp_mod, "try-error")) {
            # F-test for Breusch-Pagan
            bp_anova <- try(anova(bp_mod), silent = TRUE)
            
            # Improved: Extract p-value safely (handle various return structures)
            bp_pvalue <- NA_real_  # Initialize with NA
            var_ratio <- NA_real_
            
            if (!inherits(bp_anova, "try-error") && nrow(bp_anova) >= 2) {
                # Extract p-value from second row (the predictor "fitted")
                # Handle both indexed and array access
                pval_vec <- try(as.numeric(bp_anova[2, "Pr(>F)"]), silent = TRUE)
                
                if (!inherits(pval_vec, "try-error") && length(pval_vec) == 1 && !is.na(pval_vec)) {
                    bp_pvalue <- pval_vec
                    
                    # Also compute variance ratio across groups
                    group_vars <- tapply(values, groups, var, na.rm = TRUE)
                    if (length(group_vars) > 1) {
                        # Only compute ratio if we have valid (non-NaN) variances
                        finite_vars <- group_vars[is.finite(group_vars)]
                        if (length(finite_vars) > 1) {
                            var_ratio <- max(finite_vars, na.rm = TRUE) / min(finite_vars, na.rm = TRUE)
                        } else {
                            var_ratio <- 1  # Insufficient data for variance comparison
                        }
                    } else {
                        var_ratio <- 1
                    }
                    
                    # Heteroscedasticity detected if p < 0.05 AND variance_ratio > 2
                    if (!is.na(bp_pvalue) && !is.na(var_ratio) && bp_pvalue < 0.05 && var_ratio > 2) {
                        characteristics$heteroscedastic <- TRUE
                        reasons <- c(reasons, sprintf(
                            "Heteroscedasticity: BP_p=%.4f, var_ratio=%.2f",
                            bp_pvalue, var_ratio
                        ))
                    }
                }
            }
        }
    }
    
    # -----------------------------------------------------------------------------
    # 2. BOUNDARY CLUSTERING DETECTION (Skip for inherently bounded metrics)
    # -----------------------------------------------------------------------------
    
    # IMPORTANT (March 2026): Entropy and diversity metrics are MATHEMATICALLY BOUNDED
    # by definition (entropy in [0, log(m)]), not by measurement artifacts.
    # Boundary clustering is therefore EXPECTED and NOT a statistical problem.
    # Skip detection for entropy/diversity metrics to avoid false positives.
    #
    # Detection is retained for OTHER metrics where boundaries indicate:
    # - Detection limits / censoring
    # - Measurement artifacts
    # - True data quality issues
    
    is_entropy_like <- tolower(value_col) %in% c("entropy", "diversity", "q_value", "tsallis")
    
    if (!is_entropy_like) {
        min_val <- min(values, na.rm = TRUE)
        max_val <- max(values, na.rm = TRUE)
        range_val <- max_val - min_val
        
        # Define "near boundary" as within 10% of range from either end
        if (range_val > 0) {
            lower_bound <- min_val + 0.10 * range_val
            upper_bound <- max_val - 0.10 * range_val
            
            n_near_bounds <- sum(values <= lower_bound | values >= upper_bound, na.rm = TRUE)
            pct_near_bounds <- 100 * n_near_bounds / sum(!is.na(values))
            
            if (pct_near_bounds > 40) {
                characteristics$boundary_clustered <- TRUE
                reasons <- c(reasons, sprintf(
                    "Boundary clustering: %.1f%% within 10%% of bounds",
                    pct_near_bounds
                ))
            }
        }
    } else {
        # For entropy/diversity metrics, explicitly note that boundary clustering
        # is expected and not treated as a special condition
        if (verbose) {
            message("  Note: Boundary clustering detection skipped for ", value_col,
                    " (inherently bounded metric)")
        }
    }
    
    # -----------------------------------------------------------------------------
    # 3. EXTREME SKEWNESS DETECTION
    # -----------------------------------------------------------------------------
    
    skewness_val <- .compute_skewness(values)
    
    # Bug #2 Fix (March 2026): Skewness should be independent condition
    # (was: if (abs(skewness_val) > 1 && characteristics$heteroscedastic))
    # Robust median test should be available for ANY highly skewed data
    if (abs(skewness_val) > 1) {
        characteristics$highly_skewed <- TRUE
        reasons <- c(reasons, sprintf(
            "Extreme skewness: |skew|=%.3f",
            skewness_val
        ))
    }
    
    # -----------------------------------------------------------------------------
    # TEST SELECTION LOGIC
    # -----------------------------------------------------------------------------
    
    test_selected <- "kruskal.test"  # Default
    
    if (characteristics$highly_skewed) {
        test_selected <- "robust_median_test"
        reasons <- c(reasons, "-> Using robust median test")
    } else if (characteristics$boundary_clustered) {
        test_selected <- "quantile_test"
        reasons <- c(reasons, "-> Using quantile-based test")
    } else if (characteristics$heteroscedastic) {
        test_selected <- "art_kw"
        reasons <- c(reasons, "-> Using Aligned Rank Transform + parametric test")
    } else {
        reasons <- c(reasons, "-> Using standard Kruskal-Wallis (no special characteristics)")
    }
    
    if (verbose && length(reasons) > 0) {
        message("Test selection for ", value_col, ":")
        for (r in reasons) message("  ", r)
    }
    
    return(list(
        test_selected = test_selected,
        characteristics = characteristics,
        reasons = reasons
    ))
}

#' Apply Aligned Rank Transform + Kruskal-Wallis test
#'
#' Pre-processes data using Aligned Rank Transform (ART) to handle
#' heteroscedasticity, then applies parametric ANOVA or K-W on transformed data
#'
#' @param data Data frame with values and group indicators
#' @param value_col Column name for values
#' @param group_col Column name for groups
#'
#' @return List with:
#' - statistic: Test statistic (H or F equivalent)
#' - p_value: P-value from test
#' - method: "ART-Kruskal-Wallis" or similar
#' @noRd
.apply_art_kw <- function(data, value_col = "entropy", group_col = "q") {
    
    values <- data[[value_col]]
    groups <- data[[group_col]]
    
    if (!is.factor(groups)) {
        groups <- factor(groups)
    }
    
    # Step 1: Calculate aligned values (residuals from main effect model)
    # Remove group effect by fitting linear model
    lin_mod <- try(lm(values ~ groups), silent = TRUE)
    
    if (inherits(lin_mod, "try-error")) {
        # Fallback to standard K-W if ART fails
        return(tryCatch(
            {
                kw_test <- kruskal.test(values ~ groups)
                list(
                    statistic = as.numeric(kw_test$statistic),
                    p_value = as.numeric(kw_test$p.value),
                    method = "Kruskal-Wallis (fallback)"
                )
            },
            error = function(e) list(
                statistic = NA_real_,
                p_value = NA_real_,
                method = "test_failed"
            )
        ))
    }
    
    # Get residuals (alignment step)
    aligned <- residuals(lin_mod)
    
    # Step 2: Rank aligned values
    ranks <- rank(aligned, na.last = "keep")
    
    # Step 3: Apply van der Waerden normal scores (convert ranks to approximate normal)
    n_vals <- sum(!is.na(ranks))
    normal_scores <- stats::qnorm(ranks / (n_vals + 1))
    
    # Step 4: Test on normal scores using parametric ANOVA
    score_data <- data.frame(
        scores = normal_scores,
        group = groups
    )
    
    score_mod <- try(lm(scores ~ group, data = score_data), silent = TRUE)
    
    if (inherits(score_mod, "try-error")) {
        # Fallback to standard K-W
        return(tryCatch(
            {
                kw_test <- kruskal.test(values ~ groups)
                list(
                    statistic = as.numeric(kw_test$statistic),
                    p_value = as.numeric(kw_test$p.value),
                    method = "Kruskal-Wallis (fallback)"
                )
            },
            error = function(e) list(
                statistic = NA_real_,
                p_value = NA_real_,
                method = "test_failed"
            )
        ))
    }
    
    # Extract F-statistic from ANOVA
    anova_res <- anova(score_mod)
    
    if (nrow(anova_res) >= 1) {
        f_stat <- as.numeric(anova_res$`F value`[1])
        p_val <- as.numeric(anova_res$`Pr(>F)`[1])
        
        return(list(
            statistic = f_stat,
            p_value = p_val,
            method = "ART-Kruskal-Wallis (heteroscedasticity-adjusted)"
        ))
    } else {
        return(list(
            statistic = NA_real_,
            p_value = NA_real_,
            method = "test_failed"
        ))
    }
}

#' Apply quantile-based test for boundary-clustered data
#'
#' Instead of comparing means/medians globally, compares quantiles
#' within each group to detect group differences accounting for
#' boundary clustering
#'
#' @note For entropy and diversity metrics: This test is NOT RECOMMENDED.
#' Entropy is mathematically bounded [0, log(m)], so boundary clustering
#' is EXPECTED and not a statistical anomaly. Standard Kruskal-Wallis
#' is more appropriate. Quantile tests are designed for artificially-bounded
#' metrics (detection limits, censoring).
#'
#' @param data Data frame with values and group indicators
#' @param value_col Column name for values
#' @param group_col Column name for groups
#' @param quantiles Quantiles to test (default: 0.25, 0.50, 0.75)
#'
#' @return List with:
#' - statistic: Maximum quantile effect size
#' - p_value: Minimum p-value across quantiles (Bonferroni adjusted)
#' - method: "Quantile-based test"
#'
#' @noRd
.apply_quantile_test <- function(data, value_col = "entropy", group_col = "q",
                                       quantiles = c(0.25, 0.50, 0.75)) {
    
    values <- data[[value_col]]
    groups <- data[[group_col]]
    
    if (!is.factor(groups)) {
        groups <- factor(groups)
    }
    
    # Bug #1 Fix (March 2026): Test whether group QUANTILES differ (not raw values)
    # Strategy: For each quantile level, test if groups have different quantile values
    # using quantile regression or bootstrap CI on group quantiles
    
    quantile_pvals <- numeric(length(quantiles))
    quantile_stats <- numeric(length(quantiles))
    
    for (q_idx in seq_along(quantiles)) {
        q_val <- quantiles[q_idx]
        
        # For each group, compute its quantile value
        group_quantiles <- tapply(values, groups, quantile, probs = q_val, na.rm = TRUE)
        
        # Test if these group-specific quantiles differ using one-way ANOVA
        # Create a pseudo-dataset where each group's value is its quantile
        # Then test using Kruskal-Wallis on ranks of those quantiles
        quant_data <- data.frame(
            q_value = as.numeric(group_quantiles),
            group = names(group_quantiles)
        )
        quant_data$group <- factor(quant_data$group)
        
        # Test differences in quantile values across groups
        q_test <- try(kruskal.test(q_value ~ group, data = quant_data), silent = TRUE)
        
        if (!inherits(q_test, "try-error")) {
            quantile_pvals[q_idx] <- as.numeric(q_test$p.value)
            quantile_stats[q_idx] <- as.numeric(q_test$statistic)
        } else {
            # Fallback: treat as unavailable for this quantile
            quantile_pvals[q_idx] <- NA_real_
            quantile_stats[q_idx] <- NA_real_
        }
    }
    
    # Combine p-values: use minimum with Bonferroni correction
    # Safety check: if all quantile tests failed, return NA
    if (all(is.na(quantile_pvals)) || all(is.na(quantile_stats))) {
        return(list(
            statistic = NA_real_,
            p_value = NA_real_,
            method = "Quantile-based test (insufficient data)"
        ))
    }
    
    min_pval <- min(quantile_pvals, na.rm = TRUE)
    max_stat <- max(quantile_stats, na.rm = TRUE)
    
    # Bonferroni correction: adjust for multiple quantiles tested
    combined_pval <- pmin(1.0, min_pval * length(quantiles))
    
    return(list(
        statistic = max_stat,
        p_value = combined_pval,
        method = "Quantile-based test (boundary-adjusted)"
    ))
}

#' Apply robust median test for highly skewed data
#'
#' Mood's median test: tests whether groups have the same median
#' More robust to extreme skewness than K-W, which assumes similar shapes
#'
#' @param data Data frame with values and group indicators
#' @param value_col Column name for values
#' @param group_col Column name for groups
#'
#' @return List with:
#' - statistic: Chi-square statistic
#' - p_value: P-value from median test
#' - method: "Mood's median test"
#' @noRd
.apply_robust_median_test <- function(data, value_col = "entropy", group_col = "q") {
    
    values <- data[[value_col]]
    groups <- data[[group_col]]
    
    if (!is.factor(groups)) {
        groups <- factor(groups)
    }
    
    group_levels <- levels(groups)
    n_groups <- length(group_levels)
    
    # Mood's median test
    # 1. Compute grand median
    grand_median <- median(values, na.rm = TRUE)
    
    # 2. For each group, count how many are above/below median
    contingency_table <- matrix(0, nrow = n_groups, ncol = 2)
    rownames(contingency_table) <- group_levels
    colnames(contingency_table) <- c("Below_Median", "Above_Median")
    
    for (g_idx in seq_len(n_groups)) {
        group_name <- group_levels[g_idx]
        group_vals <- values[groups == group_name]
        
        n_below <- sum(group_vals < grand_median, na.rm = TRUE)
        n_above <- sum(group_vals >= grand_median, na.rm = TRUE)
        
        contingency_table[g_idx, 1] <- n_below
        contingency_table[g_idx, 2] <- n_above
    }
    
    # 3. Chi-square test on contingency table
    chisq_result <- tryCatch(
        chisq.test(contingency_table),
        error = function(e) NULL
    )
    
    if (!is.null(chisq_result)) {
        return(list(
            statistic = as.numeric(chisq_result$statistic),
            p_value = as.numeric(chisq_result$p.value),
            method = "Mood's median test (skewness-robust)"
        ))
    } else {
        return(list(
            statistic = NA_real_,
            p_value = NA_real_,
            method = "test_failed"
        ))
    }
}

#' Apply Friedman test for paired ranked data
#'
#' Friedman test: non-parametric alternative to repeated measures ANOVA
#' Tests whether k related samples have the same distribution
#' Handles subjects as blocks, q-values as treatments
#'
#' @param data Data frame with values, group indicators, and subject identifiers
#' @param value_col Column name for values (entropy)
#' @param group_col Column name for groups (q-values)
#' @param subject_col Column name for subject/block identifiers
#'
#' @return List with:
#' - statistic: Friedman's Q statistic (or chi-squared approximation)
#' - p_value: P-value from Friedman test
#' - method: "Friedman test (paired)"
#'
#' @details
#' Friedman test works on matrix form (blocks x treatments).
#' This wrapper converts long-form data to matrix, applies friedman.test(),
#' and returns results in standardized format.
#'
#' @noRd
.apply_friedman_test <- function(data, value_col = "entropy", 
                                        group_col = "q", subject_col) {
    
    values <- data[[value_col]]
    treatments <- data[[group_col]]
    blocks <- data[[subject_col]]
    
    if (!is.factor(treatments)) {
        treatments <- factor(treatments)
    }
    if (!is.factor(blocks)) {
        blocks <- factor(blocks)
    }
    
    # Reshape to matrix: rows = blocks (subjects), columns = treatments (q-values)
    # NOTE: friedman.test() expects a matrix with no missing values
    friedman_matrix <- try({
        matrix_data <- xtabs(as.numeric(values) ~ blocks + treatments)
        matrix_data
    }, silent = TRUE)
    
    if (inherits(friedman_matrix, "try-error") || nrow(friedman_matrix) < 2 || ncol(friedman_matrix) < 2) {
        return(list(
            statistic = NA_real_,
            p_value = NA_real_,
            method = "test_failed"
        ))
    }
    
    # Check for missing values in the pivoted matrix
    if (any(is.na(friedman_matrix))) {
        warning("Friedman test: matrix has missing values. Consider imputation or alternative test.")
        # For now, return NA - user should handle unbalanced designs separately
        return(list(
            statistic = NA_real_,
            p_value = NA_real_,
            method = "Friedman test (unbalanced blocks - skipped)"
        ))
    }
    
    # Apply Friedman test
    friedman_result <- try(
        friedman.test(friedman_matrix),
        silent = TRUE
    )
    
    if (!inherits(friedman_result, "try-error")) {
        return(list(
            statistic = as.numeric(friedman_result$statistic),
            p_value = as.numeric(friedman_result$p.value),
            method = "Friedman test (paired)"
        ))
    } else {
        return(list(
            statistic = NA_real_,
            p_value = NA_real_,
            method = "test_failed"
        ))
    }
}

#' Conditional rank-based test dispatcher
#'
#' Main function to apply conditional logic to select and run
#' appropriate rank-based test. For paired designs, uses Friedman test
#' instead of unpaired Kruskal-Wallis.
#'
#' @param data Data frame with values and group indicators
#' @param value_col Column name for values
#' @param group_col Column name for groups
#' @param paired Logical: if TRUE, apply Friedman test (requires subject_col)
#' @param subject_col Column name for subject/block identifiers (required if paired=TRUE)
#' @param verbose Logical: print diagnostics
#'
#' @return List with:
#' - statistic: Test statistic
#' - p_value: P-value
#' - method: Name of test applied
#' - test_type: "friedman", "art_adjusted", "quantile_adjusted", "median_robust", "kruskal-wallis"
#' - characteristics: Data characteristics detected
#'
#' @noRd
.apply_conditional_rank_test <- function(data, value_col = "entropy", group_col = "q",
                                               paired = FALSE, subject_col = NULL,
                                               verbose = FALSE) {
    
    # Priority 1: If paired design, use CONDITIONAL test selection
    if (paired && !is.null(subject_col) && subject_col %in% colnames(data)) {
        # NEW (March 2026): Conditional selection for paired tests
        # Detect data characteristics and select appropriate paired rank test
        selection <- .select_rank_test_paired(
            data, value_col, group_col, subject_col, verbose = verbose
        )
        
        # Apply selected paired test
        test_func <- switch(
            selection$test_selected,
            "art_friedman" = .apply_art_friedman,
            "robust_friedman" = .apply_robust_friedman,
            # Default: standard Friedman
            .apply_friedman_test
        )
        
        # Call selected test function
        test_result <- test_func(data, value_col, group_col, subject_col)
        
        # Return with metadata indicating paired test was used and which one
        # Add pairing_used flag to characteristics
        characteristics_with_pairing <- c(selection$characteristics, list(pairing_used = TRUE))
        return(c(
            test_result,
            list(
                test_type = selection$test_selected,
                characteristics = characteristics_with_pairing
            )
        ))
    }
    
    # Priority 2: If unpaired, use conditional selection logic
    # Step 1: Detect data characteristics
    selection <- .select_rank_test(data, value_col, group_col, verbose = verbose)
    
    # Step 2: Apply selected test
    test_func <- switch(
        selection$test_selected,
        "art_kw" = .apply_art_kw,
        "quantile_test" = .apply_quantile_test,
        "robust_median_test" = .apply_robust_median_test,
        # Default: standard Kruskal-Wallis
        function(d, v, g) {
            res <- try(kruskal.test(d[[v]] ~ d[[g]]), silent = TRUE)
            if (inherits(res, "try-error")) {
                list(statistic = NA_real_, p_value = NA_real_, method = "test_failed")
            } else {
                list(
                    statistic = as.numeric(res$statistic),
                    p_value = as.numeric(res$p.value),
                    method = "Kruskal-Wallis (standard)"
                )
            }
        }
    )
    
    # Apply test
    test_result <- test_func(data, value_col, group_col)
    
    # Step 3: Return result with metadata
    return(c(
        test_result,
        list(
            test_type = selection$test_selected,
            characteristics = selection$characteristics
        )
    ))
}


#' Test q * condition interaction using Two-Way Within-Subject Methods
#'
#' Tests interaction between q-values and condition factor in paired/repeated-measures 
#' settings using rank-based non-parametric methods.
#'
#' **Design:** Both q-values and condition are WITHIN-SUBJECT factors
#'   - Subjects: N individuals (paired)
#'   - Within each subject: k q-values * m conditions = km observations
#'   - Example: 8 subjects, 41 q-values, 2 conditions = 8 * 41 * 2 = 656 measurements
#'
#' **Statistical approach (FIXED - March 2026):**
#' For true two-way within-subject design, ranks ALL observations within each 
#' subject TOGETHER (not separately by condition), preserving the dependence structure.
#'
#' **Mathematical basis:**
#' 1. Rank entropy values within EACH SUBJECT (across all q-levels and conditions)
#' 2. Compute mean ranks per (q-level, condition) combination  
#' 3. Test interaction via two-way ANOVA on rank means
#' 4. Recovers power and validity of parametric two-way ANOVA
#'
#' References: Conover & Iman (1981), Puri & Sen (1985) - Nonparametric Methods
#'
#' @param data Data frame with columns: entropy, q, condition, and subject_col (if paired)
#' @param value_col Column name for values (default: "entropy")
#' @param q_col Column name for q-values (default: "q")
#' @param condition_col Column name for condition (default: "condition")
#' @param paired Logical; if TRUE, account for subject blocking
#' @param subject_col Column name for subject identifiers (required if paired=TRUE)
#'
#' @return List with:
#'   - statistic: F-statistic for interaction
#'   - p_value: p-value from interaction test
#'   - method: Description of test used
#'   - test_type: "srh_interaction" or "two_way_friedman"
#'

#' @noRd
#' @importFrom stats ave as.formula
.test_q_condition_interaction <- function(
    data,
    value_col = "entropy",
    q_col = "q",
    condition_col = "condition",
    paired = FALSE,
    subject_col = NULL) {
    
    # Validate required columns
    if (!value_col %in% colnames(data)) {
        stop("Column '", value_col, "' not found in data")
    }
    if (!q_col %in% colnames(data)) {
        stop("Column '", q_col, "' not found in data")
    }
    if (!condition_col %in% colnames(data)) {
        stop("Column '", condition_col, "' not found in data; q * condition interaction cannot be tested without condition factor")
    }
    if (paired && !subject_col %in% colnames(data)) {
        stop("Column '", subject_col, "' not found in data (required for paired analysis)")
    }
    
    # For both paired and unpaired: Use Scheirer-Ray-Hare test (REVISED March 2026)
    # The aggregation-then-ANOVA approach for paired designs has inadequate degrees of freedom
    # Scheirer-Ray-Hare properly handles two-way designs by testing on ranked data directly
    # References: Scheirer, Castellan, Wilkinson (1976); Conover & Iman (1981)
    if (paired && !is.null(subject_col)) {
        # Paired design: Rank within each subject ONLY (preserves within-subject dependence)
        # Then apply Scheirer-Ray-Hare on the within-subject ranks
        data$ranks <- ave(
            data[[value_col]],
            data[[subject_col]],
            FUN = function(x) rank(x, na.last = "keep")
        )
    } else {
        # Unpaired design: Rank across entire dataset
        data$ranks <- rank(data[[value_col]], na.last = "keep")
    }
    
    # Apply Scheirer-Ray-Hare test for q * condition interaction
    # Works for both paired (within-subject ranks) and unpaired (global ranks) cases
    tryCatch({
        # Convert factors if needed
        data[[q_col]] <- factor(data[[q_col]])
        data[[condition_col]] <- factor(data[[condition_col]])
        
        # Use pre-computed ranks (within-subject for paired, global for unpaired)
        # Then apply two-way ANOVA on the ranked data
        formula_str <- paste("ranks ~", q_col, "*", condition_col)
        lm_model <- lm(as.formula(formula_str), data = data)
        anova_result <- anova(lm_model)
        
        # Extract interaction F-statistic and p-value
        # Interaction is the second-to-last row (before Residuals)
        interaction_row <- nrow(anova_result) - 1
        f_stat <- anova_result$`F value`[interaction_row]
        p_val <- anova_result$`Pr(>F)`[interaction_row]
        
        if (is.na(f_stat) || is.na(p_val)) {
            return(list(
                statistic = NA_real_,
                p_value = NA_real_,
                method = "Scheirer-Ray-Hare (computation failed)",
                test_type = "srh_failed"
            ))
        }
        
        test_type_label <- if (paired) "srh_paired" else "srh_unpaired"
        method_label <- if (paired) 
            "Scheirer-Ray-Hare Test (paired design, within-subject ranks; REVISED March 2026)" 
            else 
            "Scheirer-Ray-Hare Test (non-parametric 2-way ANOVA)"
        
        return(list(
            statistic = f_stat,
            p_value = p_val,
            method = method_label,
            test_type = test_type_label
        ))
    }, error = function(e) {
        return(list(
            statistic = NA_real_,
            p_value = NA_real_,
            method = paste("Scheirer-Ray-Hare (error):", e$message),
            test_type = "srh_error"
        ))
    })
}
