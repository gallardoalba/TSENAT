// =============================================================================
// Core Resampling and Entropy Computation in Rcpp + Armadillo
// =============================================================================
// 
// Purpose: Accelerate core statistical computations with Rcpp/Armadillo
// Target speedup: 2-3x~ compared to pure R implementation
// 
// This file contains:
//   1. Entropy computations: entropy_cpp, hill_number_cpp
//   2. Jackknife resampling: jackknife_resampling_cpp, jis_jackknife_influences_cpp
//   3. Bootstrap resampling: bootstrap_compute_cpp, block_bootstrap_compute_cpp
//   4. Bootstrap entropy: bootstrap_entropy_vec_cpp, jis_bootstrap_delta_cpp
//   5. Divergence computations: tsallis_divergence_cpp, divergence_bootstrap_compute_cpp
//   6. Utility functions: check_rcpp_available
//
// Implementation uses Armadillo (via RcppArmadillo) for efficient matrix operations
// and vectorized computation to replace R's row-by-row loops.
// =============================================================================

#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;

// Internal helper: Compute Shannon or Tsallis entropy for proportions
// [[Rcpp::export(rng = false)]]
double entropy_cpp(NumericVector p, double q = 1.0, bool normalize = true, double log_base = 2.718281828) {
  // BUG FIX: Validate parameters before computation
  // log_base must be > 1 to avoid division by zero or negative log values
  if (log_base <= 0 || std::abs(log_base - 1.0) < 1e-10) {
    Rcpp::warning("Invalid log_base (must be > 1, not equal to 1)");
    return NA_REAL;
  }
  
  // q must be non-negative for Tsallis entropy to be defined
  if (q < 0) {
    Rcpp::warning("Invalid q parameter (must be non-negative)");
    return NA_REAL;
  }
  
  // Remove NA values
  LogicalVector not_na = !is_na(p);
  NumericVector p_clean = p[not_na];
  
  if (p_clean.size() == 0) return NA_REAL;
  
  // Filter out zeros and negative values
  LogicalVector valid = (p_clean > 1e-10);
  NumericVector p_valid = p_clean[valid];
  
  if (p_valid.size() == 0) return NA_REAL;
  
  // Convert size to int safely
  int n = static_cast<int>(p_valid.size());
  if (n <= 0) {
    Rcpp::warning("Input size invalid");
    return NA_REAL;
  }
  
  double entropy = 0.0;
  double q_tol = 1e-6;
  
  // Special case for q ≈ 0 (species richness)
  if (q < q_tol) {
    // D_0 = effective richness = number of species (entropy of richness is just species count)
    entropy = std::log(n);  // Raw richness in natural log
    if (normalize) {
      // BUG FIX: Handle n=1 case where log(n) = 0 to avoid 0/0 = NaN
      if (n == 1) {
        entropy = 0.0;  // Normalized richness of 1 species is 0
      } else {
        entropy = entropy / std::log(n);  // Normalized: already 1.0 since max richness = n
      }
    }
    return entropy;
  }
  
  if (std::abs(q - 1.0) < q_tol) {
    // Shannon entropy (q = 1)
    for (int i = 0; i < n; i++) {
      double pi = p_valid[i];
      if (pi > 1e-15) {  // BUG FIX: Use machine epsilon threshold not arbitrary 0
        entropy -= pi * std::log(pi) / std::log(log_base);
      }
    }
  } else {
    // Tsallis entropy (q != 1): (1 - sum(p^q)) / (q - 1)
    // NOTE: This formula does NOT include log_base
    double sum_pq = 0.0;
    for (int i = 0; i < n; i++) {
      sum_pq += std::pow(p_valid[i], q);
    }
    
    entropy = (1.0 - sum_pq) / (q - 1.0);  // BUG FIX: Removed * std::log(log_base)
  }
  
  // Normalize by maximum entropy if requested
  if (normalize) {
    // BUG FIX: Use size_t consistently, avoid int conversion issues
    int n = static_cast<int>(p_valid.size());
    
    // Edge case: only one species has entropy 0
    if (n == 1) {
      entropy = 0.0;  // Set to exactly 0 for single species
    } else {
      double max_entropy;
      
      if (q < q_tol) {
        // Species richness: max = log(n)
        max_entropy = std::log(n);
      } else if (std::abs(q - 1.0) < q_tol) {
        // Shannon: max = log(n)
        max_entropy = std::log(n) / std::log(log_base);
      } else {
        // Tsallis: max = (1 - n^(1-q)) / (q - 1)
        // BUG FIX: Removed * std::log(log_base) from denominator
        max_entropy = (1.0 - std::pow(n, 1.0 - q)) / (q - 1.0);
      }
      
      // BUG FIX: Threshold should be machine epsilon, not arbitrary 1e-10
      if (max_entropy > std::numeric_limits<double>::epsilon() && std::isfinite(max_entropy)) {
        entropy = entropy / max_entropy;
      }
      // If max_entropy <= 0 or non-finite, skip normalization (entropy stays unnormalized)
    }
  }
  
  return entropy;
}

// ============================================================================
// HILL NUMBERS (Diversity indices derived from entropy)
// ============================================================================
// Hill numbers (D_q) convert Tsallis entropy to effective number of species
// Formula: D_q = (Σp^q)^(1/(1-q))
// Special cases:
//   D_0 = count(p > 0)  [true species richness]
//   D_1 = exp(H_1)      [exponential of Shannon entropy]
//   D_q = general form for q ≠ 1
// [[Rcpp::export(rng = false)]]
double hill_number_cpp(NumericVector p, double q = 1.0, double log_base = 2.718281828) {
  // BUG FIX: Validate log_base parameter (see entropy_cpp for details)
  if (log_base <= 0 || std::abs(log_base - 1.0) < 1e-10) {
    Rcpp::warning("Invalid log_base (must be > 1, not equal to 1)");
    return NA_REAL;
  }
  
  // q must be non-negative
  if (q < 0) {
    Rcpp::warning("Invalid q parameter (must be non-negative)");
    return NA_REAL;
  }
  
  // Remove NA and filter invalid values
  LogicalVector not_na = !is_na(p);
  NumericVector p_clean = p[not_na];
  
  if (p_clean.size() == 0) return NA_REAL;
  
  LogicalVector valid = (p_clean > 1e-10);
  NumericVector p_valid = p_clean[valid];
  
  if (p_valid.size() == 0) return NA_REAL;
  
  double q_tol = 1e-6;
  int n = p_valid.size();
  
  if (q < q_tol) {
    // D_0: true richness = number of nonzero species
    return static_cast<double>(n);
  } else if (std::abs(q - 1.0) < q_tol) {
    // D_1 = exp(Shannon entropy)
    double shannon = 0.0;
    for (int i = 0; i < n; i++) {
      double pi = p_valid[i];
      if (pi > 1e-15) {  // BUG FIX: Use machine epsilon threshold not arbitrary 0
        shannon -= pi * std::log(pi) / std::log(log_base);
      }
    }
    return std::pow(log_base, shannon);
  } else {
    // D_q = (Σp^q)^(1/(1-q))
    double sum_pq = 0.0;
    for (int i = 0; i < n; i++) {
      sum_pq += std::pow(p_valid[i], q);
    }
    
    // BUG FIX: Check for numerical stability - avoid log domain if not needed
    // Also validate sum_pq is in valid range (0, infinity)
    if (sum_pq <= 0 || !std::isfinite(sum_pq)) {
      return NA_REAL;
    }
    
    double exponent = 1.0 / (1.0 - q);
    
    // For very large exponents, use log-domain computation to avoid overflow
    // log(result) = exponent * log(sum_pq)
    double log_result = exponent * std::log(sum_pq);
    if (std::abs(log_result) > 700) {  // 700 ~ ln(DBL_MAX), avoid overflow
      return std::exp(log_result);  // Let overflow be handled by exp()
    }
    return std::pow(sum_pq, exponent);
  }
}
// Internal C++ function (registered but not exported to R NAMESPACE)
// Leave-one-out jackknife with vectorized Armadillo implementation
// Callable via .Call("_TSENAT_jackknife_resampling_cpp") but not in user namespace
// [[Rcpp::export]]
List jackknife_resampling_cpp(NumericMatrix counts, double q = 1.0, 
                              bool normalize = true, double log_base = 2.718281828,
                              double pseudocount = 0.0) {
  // BUG FIX: Validate q and log_base parameters at entry
  if (log_base <= 0 || std::abs(log_base - 1.0) < 1e-10) {
    return List::create(Named("error") = "Invalid log_base");
  }
  if (q < 0) {
    return List::create(Named("error") = "Invalid q parameter");
  }
  
  int n_obs = counts.nrow();
  int n_samples = counts.ncol();
  
  if (n_obs < 2) {
    Rcpp::warning("Insufficient observations for jackknife (need >= 2)");
    return List::create(Named("error") = "Insufficient observations");
  }
  
  // BUG FIX: Validate n_samples > 0 (critical for column operations)
  if (n_samples < 1) {
    Rcpp::warning("No samples in count matrix");
    return List::create(Named("error") = "Invalid count matrix");
  }
  
  // Convert to Armadillo matrix for efficient operations
  arma::mat counts_arma(counts.begin(), n_obs, n_samples, false);
  
  // Compute original full estimate
  arma::vec col_sums = arma::sum(counts_arma, 0).t();
  // BUG FIX: Use double multiplication to avoid integer overflow with huge n_samples
  double total = arma::accu(counts_arma) + (double)n_samples * pseudocount;
  
  if (total <= 0) {
    Rcpp::warning("Total count is zero or negative");
    return List::create(Named("error") = "Invalid total count");
  }
  
  arma::vec p_full = (col_sums + pseudocount) / total;
  // BUG FIX: Avoid double conversion. Use implicit conversion or direct NumericVector constructor
  double estimate = entropy_cpp(NumericVector(p_full.begin(), p_full.end()), q, normalize, log_base);
  
  if (std::isnan(estimate) || std::isinf(estimate)) {
    Rcpp::warning("Invalid entropy estimate");
    return List::create(Named("error") = "Invalid estimate");
  }
  
  // TRUE VECTORIZATION: Compute row sums once
  // Then for each row i, compute adjusted sums without creating new matrices
  arma::vec row_sums = arma::sum(counts_arma, 1);
  
  // Leave-one-out jackknife computation (VECTORIZED - no matrix creation)
  NumericVector jackknife_estimates(n_obs);
  
  for (int i = 0; i < n_obs; i++) {
    // Compute adjusted sums WITHOUT creating a new matrix
    // col_sums_minus_i = col_sums - counts[i,]
    arma::vec col_sums_minus_i = col_sums - counts_arma.row(i).t();
    
    // total_minus_i = total - row_sum[i]
    double total_minus_i = total - row_sums[i];
    
    if (total_minus_i <= 0) {
      jackknife_estimates[i] = NA_REAL;
      continue;
    }
    
    // Compute proportions for this leave-one-out estimate
    arma::vec p_minus_i = (col_sums_minus_i + pseudocount) / total_minus_i;
    // BUG FIX: Avoid double conversion wrap() -> as<>(). Pass arma::vec directly via implicit conversion
    jackknife_estimates[i] = entropy_cpp(NumericVector(p_minus_i.begin(), p_minus_i.end()), q, normalize, log_base);
  }
  
  // Compute influence (fully vectorized)
  NumericVector influence(n_obs);
  for (int i = 0; i < n_obs; i++) {
    if (std::isnan(jackknife_estimates[i])) {
      influence[i] = NA_REAL;
    } else {
      influence[i] = std::abs(jackknife_estimates[i] - estimate);
    }
  }
  
  // Compute jackknife standard error (vectorized)
  double theta_mean = 0.0;
  int valid_count = 0;
  for (int i = 0; i < n_obs; i++) {
    if (!std::isnan(jackknife_estimates[i])) {
      theta_mean += jackknife_estimates[i];
      valid_count++;
    }
  }
  
  if (valid_count > 0) {
    theta_mean /= valid_count;
  }
  
  double se_sum = 0.0;
  for (int i = 0; i < n_obs; i++) {
    if (!std::isnan(jackknife_estimates[i])) {
      double diff = jackknife_estimates[i] - theta_mean;
      se_sum += diff * diff;
    }
  }
  
  double jackknife_se;
  if (valid_count == 0) {
    jackknife_se = NA_REAL;  // All estimates are NaN, so SE is undefined
  } else {
    jackknife_se = std::sqrt(((n_obs - 1.0) / n_obs) * se_sum);
  }
  
  return List::create(
    Named("estimate") = estimate,
    Named("jackknife_estimates") = jackknife_estimates,
    Named("influence") = influence,
    Named("jackknife_se") = jackknife_se,
    Named("q") = q,
    Named("normalize") = normalize
  );
}

// ============================================================================
// BOOTSTRAP RESAMPLING WITH ENTROPY COMPUTATION (C++ OPTIMIZATION)
// ============================================================================
// Purpose: Accelerate bootstrap resampling by computing entropy for each
//          replicate in C++ rather than via R apply() loop.
// Target speedup: 2-3x compared to pure R implementation
// 
// The bottleneck in R is: apply(bootstrap_samples, 2, entropy_calc)
// This eliminates R function call overhead and vectorizes the entropy loop.

// [[Rcpp::export]]
NumericVector bootstrap_compute_cpp(NumericVector x, int nboot = 1000, 
                                    double q = 1.0, bool normalize = true, 
                                    double log_base = 2.718281828,
                                    double pseudocount = 0.0) {
  // BUG FIX: Validate parameters at entry
  if (log_base <= 0 || std::abs(log_base - 1.0) < 1e-10) {
    Rcpp::stop("Invalid log_base (must be > 1, not equal to 1)");
  }
  if (q < 0) {
    Rcpp::stop("Invalid q parameter (must be non-negative)");
  }
  
  // Input validation
  int n = x.size();
  if (n < 1) {
    Rcpp::stop("Input vector x must have at least 1 element");
  }
  if (nboot < 1) {
    Rcpp::stop("nboot must be at least 1");
  }
  
  // Adjust counts with pseudocount (scalar)
  NumericVector x_adj = x + pseudocount;
  double total = sum(x_adj);
  
  if (total <= 0) {
    // Return vector of NAs for zero-count genes (consistent with non-bootstrap behavior)
    // This allows bootstrap to handle pseudocount=0 with zero-count genes gracefully
    return NumericVector(nboot, NA_REAL);
  }
  
  // Compute proportions from original data (for resampling)
  NumericVector p_hat = x_adj / total;
  
  // Pre-allocate result vector
  NumericVector boot_dist(nboot);
  
  // Get R's RNG state for reproducibility
  GetRNGstate();
  
  // Main bootstrap loop (tight C++ loop for speed)
  // Each iteration: resample counts -> normalize to proportions -> compute entropy
  for (int b = 0; b < nboot; b++) {
    // Generate multinomial bootstrap sample using R::rmultinom
    IntegerVector boot_sample_int(n);
    // BUG FIX: Check total against INT_MAX before conversion to avoid overflow
    if (total > static_cast<double>(INT_MAX)) {
      Rcpp::stop("Total count exceeds maximum integer value (%d)", INT_MAX);
    }
    int total_int = static_cast<int>(std::round(total));
    if (total_int <= 0) total_int = 1;  // Safety check
    R::rmultinom(total_int, p_hat.begin(), n, boot_sample_int.begin());
    
    // Convert to NumericVector for proportion calculation
    NumericVector boot_counts = as<NumericVector>(boot_sample_int);
    
    // CRITICAL: Normalize bootstrap sample from counts to proportions
    // entropy_cpp() expects proportions, not counts
    double boot_total = sum(boot_counts);
    
    // Safety check: avoid division by zero
    if (boot_total <= 0) {
      boot_dist[b] = NA_REAL;
      continue;
    }
    
    NumericVector boot_props = boot_counts / boot_total;
    
    // Compute Tsallis entropy for this bootstrap replicate
    boot_dist[b] = entropy_cpp(boot_props, q, normalize, log_base);
    
    // Check for invalid estimates
    if (!std::isfinite(boot_dist[b])) {
      boot_dist[b] = NA_REAL;
    }
  }
  
  // Restore R's RNG state
  PutRNGstate();
  
  return boot_dist;
}

// ============================================================================
// BOOTSTRAP ENTROPY VECTORIZATION (for pre-generated samples)
// ============================================================================
// Purpose: Fast entropy computation across multiple bootstrap samples
//          Takes pre-computed bootstrap samples matrix and returns entropy
//          for each column (replicate).
// Use case: When bootstrap samples are generated at R level.

// [[Rcpp::export]]
NumericVector bootstrap_entropy_vec_cpp(NumericMatrix boot_samples, 
                                        double q = 1.0, 
                                        bool normalize = true, 
                                        double log_base = 2.718281828) {
  int nboot = boot_samples.ncol();
  
  if (nboot < 1) {
    Rcpp::stop("bootstrap_samples matrix must have at least 1 column");
  }
  
  NumericVector boot_dist(nboot);
  
  // Vectorized entropy computation across all bootstrap replicates
  // Each column = one bootstrap sample
  for (size_t b = 0; b < static_cast<size_t>(nboot); b++) {  // BUG FIX: Type consistency
    NumericVector boot_sample = boot_samples(_, (int)b);  // Cast back for indexing
    boot_dist[b] = entropy_cpp(boot_sample, q, normalize, log_base);
    
    if (!std::isfinite(boot_dist[b])) {
      boot_dist[b] = NA_REAL;
    }
  }
  
  return boot_dist;
}

// ============================================================================
// BLOCK BOOTSTRAP FOR PAIRED SAMPLES (C++ OPTIMIZATION)
// ============================================================================
// Purpose: Accelerate block bootstrap (paired sample resampling) by generating
//          bootstrap replicates in C++ and computing entropy efficiently.
// 
// Block bootstrap maintains correlation structure by resampling pairs as units.
// Input x has pairs (x[0],x[1]), (x[2],x[3]), ..., (x[2n-2],x[2n-1])
// Each bootstrap replicate resamples n_pairs pairs with replacement.

// [[Rcpp::export]]
NumericVector block_bootstrap_compute_cpp(NumericVector x, int nboot = 1000, 
                                          double q = 1.0, bool normalize = true, 
                                          double log_base = 2.718281828,
                                          double pseudocount = 0.0) {
  // BUG FIX: Validate parameters before computation
  if (log_base <= 0 || std::abs(log_base - 1.0) < 1e-10) {
    Rcpp::stop("Invalid log_base (must be > 1, not equal to 1)");
  }
  if (q < 0) {
    Rcpp::stop("Invalid q parameter (must be non-negative)");
  }
  if (pseudocount < 0) {
    Rcpp::stop("pseudocount must be non-negative");
  }
  
  // Input validation
  int n = x.size();
  if (n < 2) {
    Rcpp::stop("Input vector x must have at least 2 elements (1 pair)");
  }
  if (n % 2 != 0) {
    Rcpp::stop("Input vector x must have even length for paired samples");
  }
  if (nboot < 1) {
    Rcpp::stop("nboot must be at least 1");
  }
  
  int n_pairs = n / 2;
  
  // Adjust counts with pseudocount (scalar) - but keep raw x for resampling
  // pseudocount will be applied during entropy calculation, not here
  NumericVector x_adj = x + pseudocount;
  double total = sum(x_adj);  // Total with pseudocount for computing proportions
  
  if (total <= 0) {
    // Return vector of NAs for zero-count genes (consistent with non-bootstrap behavior)
    // This allows block bootstrap to handle pseudocount=0 with zero-count genes gracefully
    return NumericVector(nboot, NA_REAL);
  }
  
  // Pre-allocate result vector
  NumericVector boot_dist(nboot);
  
  // Get R's RNG state for reproducibility
  GetRNGstate();
  
  // Main block bootstrap loop
  // For each bootstrap replicate, resample pairs with replacement
  for (int b = 0; b < nboot; b++) {
    // Initialize bootstrap sample for this replicate
    NumericVector boot_sample(n, 0.0);
    
    // BUG FIX: Use Rcpp::sample() directly for efficient pair index resampling
    // Much faster than rmultinom(1, ...) in a loop, and avoids precision issues
    IntegerVector pair_indices = Rcpp::sample(n_pairs, n_pairs, true) - 1;  // -1 for 0-based indexing
    
    // Reconstruct bootstrap sample from sampled pair indices
    for (int p = 0; p < n_pairs; p++) {
      int sampled_pair_idx = pair_indices[p];  // Already 0-based
      int original_idx1 = 2 * sampled_pair_idx;
      int original_idx2 = 2 * sampled_pair_idx + 1;
      
      boot_sample[2 * p] = x[original_idx1];  // Use raw x, not x_adj
      boot_sample[2 * p + 1] = x[original_idx2];
    }
    
    // BUG FIX: Apply pseudocount to bootstrap sample proportions for entropy calc
    // Add pseudocount for consistency with bootstrap_compute_cpp
    NumericVector boot_sample_adj = boot_sample + pseudocount;
    double boot_total = sum(boot_sample_adj);
    // Compute proportions from adjusted bootstrap sample
    NumericVector boot_props = boot_sample_adj / boot_total;
    
    // Compute Tsallis entropy for this block bootstrap replicate
    boot_dist[b] = entropy_cpp(boot_props, q, normalize, log_base);
    
    // Check for invalid estimates
    if (!std::isfinite(boot_dist[b])) {
      boot_dist[b] = NA_REAL;
    }
  }
  
  // Restore R's RNG state
  PutRNGstate();
  
  return boot_dist;
}

// Fallback: R wrapper for graceful degradation if Rcpp not available
//' @keywords internal
// [[Rcpp::export]]
bool check_rcpp_available() {
  return true;  // If this function exists, Rcpp compilation succeeded
}

// ============================================================================
// ISOFORM SWITCHING JIS OPTIMIZATION (C++ ACCELERATED)
// ============================================================================
// Purpose: Accelerate jackknife_isoform_switching computation by moving
//          intensive matrix operations and loops to C++
// Target speedup: 2-3x compared to pure R implementation
// Candidate functions: 
//   1. jis_tsallis_entropy_cpp() - Fast entropy for all samples in matrix
//   2. jis_jackknife_influences_cpp() - Leave-one-out for all transcripts
//   3. jis_bootstrap_delta_cpp() - Bootstrap delta statistics

// [[Rcpp::export]]
NumericVector jis_tsallis_entropy_cpp(NumericMatrix counts, 
                                      double q = 1.0, 
                                      bool normalize = true, 
                                      double log_base = 2.718281828,
                                      double pseudocount = 1e-8,
                                      int n_tx_fixed = -1) {
  // BUG FIX: Validate parameters before computation
  if (log_base <= 0 || std::abs(log_base - 1.0) < 1e-10) {
    Rcpp::warning("Invalid log_base parameter");
    return NumericVector(0);  // Return empty on invalid params
  }
  if (q < 0) {
    Rcpp::warning("Invalid q parameter (must be non-negative)");
    return NumericVector(0);
  }
  
  // Vectorized Tsallis entropy computation for each column of counts matrix
  // counts: n_transcripts (rows) × n_samples (columns) matrix
  // Returns: numeric vector of length n_samples with entropy for each sample
  
  int n_tx = counts.nrow();
  int n_samples = counts.ncol();
  
  // BUG FIX: Validate input dimensions before using them
  if (n_tx < 1 || n_samples < 1) {
    return NumericVector(0);  // Empty output for empty input
  }
  
  // Get fixed n_tx if provided, else use current dimensions
  int n_tx_use = (n_tx_fixed > 0) ? n_tx_fixed : n_tx;
  
  NumericVector entropy_result(n_samples);
  arma::mat counts_arma(counts.begin(), n_tx, n_samples, false);
  
  for (int sample = 0; sample < n_samples; sample++) {
    // Get column (sample)
    arma::vec col = counts_arma.col(sample);
    
    // BUG FIX: Apply pseudocount BEFORE computing sum (was missing!)
    // This is consistent with jackknife_resampling_cpp behavior
    col = col + pseudocount;
    double col_sum = arma::accu(col);
    
    // Safety check: ensure col_sum is positive
    if (col_sum <= 0) {
      entropy_result[sample] = NA_REAL;
      continue;
    }
    
    // Compute proportions (now includes pseudocount effect)
    arma::vec p = col / col_sum;
    
    // Compute Tsallis entropy
    double h = std::nan("");  // BUG FIX #17: Initialize to NaN to catch logic errors
    double q_tol = 1e-6;
    
    if (q < q_tol) {
      // BUG FIX: Add missing q≈0 case (richness) in JIS entropy
      // D_0 = log(n_tx_use) (richness entropy)
      int n_nonzero = 0;
      for (size_t i = 0; i < p.n_elem; i++) {
        if (p[i] > 1e-15) n_nonzero++;
      }
      h = std::log(std::max(1, n_nonzero));
    } else if (std::abs(q - 1.0) < q_tol) {
      // Shannon entropy
      h = 0.0;
      for (size_t i = 0; i < p.n_elem; i++) {  // BUG FIX: Use armadillo size, not separate variable
        if (p[i] > 1e-15) {  // BUG FIX: Machine epsilon threshold
          h -= p[i] * std::log(p[i]) / std::log(log_base);
        }
      }
    } else {
      // Tsallis entropy: H_q = (1 - sum(p_i^q)) / (q - 1)
      double p_q_sum = 0.0;
      for (size_t i = 0; i < p.n_elem; i++) {  // BUG FIX: Use armadillo size
        p_q_sum += std::pow(p[i], q);
      }
      // Use proper formula: (1 - sum(p^q)) / (q - 1)
      // Note: Unlike Shannon, log_base doesn't appear in the same way
      h = (1.0 - p_q_sum) / (q - 1.0);
    }
    
    // Set to NA if invalid
    if (!std::isfinite(h)) {
      entropy_result[sample] = NA_REAL;
      continue;
    }
    
    // Normalize if requested
    if (normalize) {
      double max_h;
      double q_tol_norm = 1e-6;
      
      if (std::abs(q - 1.0) < q_tol_norm) {
        // Shannon: max = log(n_tx_use)
        max_h = std::log(static_cast<double>(n_tx_use)) / std::log(log_base);
      } else if (q < q_tol_norm) {
        // BUG FIX: Add missing q≈0 case (richness)
        max_h = std::log(static_cast<double>(n_tx_use));
      } else {
        // Tsallis: max = (1 - n^(1-q)) / (q - 1)
        double numerator = 1.0 - std::pow(static_cast<double>(n_tx_use), 1.0 - q);
        double denominator = q - 1.0;
        max_h = numerator / denominator;  // Result is positive by Tsallis definition
      }
      
      // BUG FIX: Use machine epsilon instead of arbitrary 1e-10 threshold
      if (!std::isnan(max_h) && std::isfinite(max_h) && max_h > std::numeric_limits<double>::epsilon()) {
        h = h / max_h;
      }
    }
    
    entropy_result[sample] = h;
  }
  
  return entropy_result;
}

// [[Rcpp::export]]
NumericVector jis_jackknife_influences_cpp(NumericMatrix counts, 
                                           double q = 1.0, 
                                           bool normalize = true, 
                                           double log_base = 2.718281828,
                                           double pseudocount = 1e-8,
                                           int n_tx_fixed = -1) {
  // BUG FIX: Validate parameters before computation
  if (log_base <= 0 || std::abs(log_base - 1.0) < 1e-10) {
    Rcpp::warning("Invalid log_base parameter");
    return NumericVector(0);
  }
  if (q < 0) {
    Rcpp::warning("Invalid q parameter (must be non-negative)");
    return NumericVector(0);
  }
  
  // Vectorized leave-one-out jackknife influence computation
  // counts: n_transcripts (rows) × n_samples (columns) matrix
  // Returns: influence vector of length n_transcripts
  //  Each value = mean absolute change in entropy when that transcript is removed
  
  int n_tx = counts.nrow();
  int n_samples = counts.ncol();
  
  // BUG FIX: Return correctly-sized empty vector if empty input
  // (not max(1, n_tx) which corrupts empty case)
  if (n_tx < 1 || n_samples < 1) {
    return NumericVector(std::max(0, n_tx));  // Return vector of actual size, even if 0
  }
  
  // If only 1 transcript, return 0 influence
  if (n_tx == 1) {
    return NumericVector::create(0.0);
  }
  
  // Compute full entropy for all samples
  NumericVector h_full = jis_tsallis_entropy_cpp(counts, q, normalize, log_base, pseudocount, n_tx_fixed);
  
  // BUG FIX #11: MEDIUM - Add defensive validation for h_full size
  // If entropy computation fails, h_full could be empty or smaller than expected
  if ((int)h_full.size() != n_samples) {
    Rcpp::warning("Entropy vector size mismatch: expected %d, got %d", 
                  n_samples, (int)h_full.size());
    return NumericVector(std::max(0, n_tx));  // Return empty on error
  }
  
  NumericVector influences(n_tx);  // One per transcript (row)
  arma::mat counts_arma(counts.begin(), n_tx, n_samples, false);
  
  // Leave-one-transcript-out loop
  for (size_t i = 0; i < static_cast<size_t>(n_tx); i++) {  // BUG FIX: Use size_t or cast explicitly
    // Create leave-i-out matrix by removing row i
    arma::mat counts_leave_i(n_tx - 1, n_samples);
    int row_idx = 0;
    for (size_t j = 0; j < static_cast<size_t>(n_tx); j++) {  // BUG FIX: Type consistency
      if (j != i) {
        counts_leave_i.row(row_idx) = counts_arma.row(j);
        row_idx++;
      }
    }
    
    // BUG FIX: Avoid creating temporary NumericMatrix for conversion to arma::mat
    // Instead, directly create NumericMatrix from arma::mat pointer
    // NumericMatrix constructor can wrap the underlying data
    Rcpp::NumericMatrix counts_leave_i_R = Rcpp::wrap(counts_leave_i);
    
    // Compute entropy without this transcript (still per sample)
    NumericVector h_leave_i = jis_tsallis_entropy_cpp(counts_leave_i_R, q, normalize, log_base, 
                                                      pseudocount, n_tx_fixed);
    
    // Compute mean absolute difference in entropy across samples
    double diffs_sum = 0.0;
    size_t valid_count = 0;  // BUG FIX: Use size_t for counting to match loop type
    for (size_t s = 0; s < static_cast<size_t>(n_samples); s++) {
      if (!std::isnan(h_full[s]) && !std::isnan(h_leave_i[s]) && 
          std::isfinite(h_full[s]) && std::isfinite(h_leave_i[s])) {
        diffs_sum += std::abs(h_full[s] - h_leave_i[s]);
        valid_count++;
      }
    }
    
    if (valid_count > 0) {
      influences[i] = diffs_sum / valid_count;
    } else {
      influences[i] = NA_REAL;
    }
  }
  
  return influences;
}

// [[Rcpp::export]]
List jis_bootstrap_delta_cpp(NumericMatrix counts_A, 
                             NumericMatrix counts_B,
                             NumericVector delta_influence,
                             double q = 1.0, 
                             bool normalize = true, 
                             double log_base = 2.718281828,
                             double pseudocount = 1e-8,
                             int n_bootstrap = 1000,
                             double confidence = 0.95,
                             String method = "percentile",
                             int n_transcripts_fixed = -1) {
  // BUG FIX: Validate parameters before computation
  if (log_base <= 0 || std::abs(log_base - 1.0) < 1e-10) {
    return List::create(Named("error") = "Invalid log_base parameter");
  }
  if (q < 0) {
    return List::create(Named("error") = "Invalid q parameter");
  }
  if (confidence <= 0 || confidence >= 1) {
    return List::create(Named("error") = "Confidence level must be in (0, 1)");
  }
  if (n_bootstrap < 1) {
    return List::create(Named("error") = "n_bootstrap must be >= 1");
  }
  
  // Bootstrap computation for delta statistics (confidence intervals and p-values)
  // counts_A, counts_B: n_transcripts (rows) × n_samples (columns) matrices
  // delta_influence: per-transcript vector of jackknife influences
  // Returns: per-transcript bootstrap statistics
  
  int n_tx = counts_A.nrow();
  int n_samples_A = counts_A.ncol();
  int n_samples_B = counts_B.ncol();
  
  // CRITICAL VALIDATION: Check all dimension consistency
  if (n_tx == 0 || delta_influence.size() != n_tx) {
    return List::create(Named("error") = "Dimension mismatch: counts and delta_influence");
  }
  if (counts_B.nrow() != n_tx) {
    return List::create(Named("error") = "Dimension mismatch: counts_A and counts_B have different numbers of transcripts");
  }
  
  arma::mat counts_A_arma(counts_A.begin(), n_tx, n_samples_A, false);
  arma::mat counts_B_arma(counts_B.begin(), counts_B.nrow(), n_samples_B, false);
  
  // Pre-allocate bootstrap delta matrix (n_bootstrap × n_tx)
  // Initialize with NaN (not using fill::value() which doesn't work with NA_REAL)
  arma::mat bootstrap_deltas(n_bootstrap, n_tx);
  bootstrap_deltas.fill(std::nan(""));
  
  GetRNGstate();
  
  // Bootstrap loop - resample columns (samples) with replacement
  for (int b = 0; b < n_bootstrap; b++) {
    // BUG FIX: Validate sample sizes > 0 before resampling
    // Empty matrices cause undefined behavior with Rcpp::sample(0, 0, true)
    if (n_samples_A < 1 || n_samples_B < 1) {
      // Skip this bootstrap iteration (leave NaN values)
      continue;
    }
    
    // Create bootstrap samples by resampling columns
    arma::mat boot_A(n_tx, n_samples_A);
    arma::mat boot_B(counts_B.nrow(), n_samples_B);
    
    // Use Rcpp::sample for proper uniform resampling with replacement
    IntegerVector idx_A_vec = Rcpp::sample(n_samples_A, n_samples_A, true) - 1;  // -1 for 0-based indexing
    IntegerVector idx_B_vec = Rcpp::sample(n_samples_B, n_samples_B, true) - 1;  // -1 for 0-based indexing
    
    for (int i = 0; i < n_samples_A; i++) {
      boot_A.col(i) = counts_A_arma.col(idx_A_vec[i]);
    }
    
    for (int i = 0; i < n_samples_B; i++) {
      boot_B.col(i) = counts_B_arma.col(idx_B_vec[i]);
    }
    
    // Convert to NumericMatrix for R functions
    // BUG FIX: Use wrap() for efficient conversion instead of element-by-element copy
    NumericMatrix boot_A_R = Rcpp::wrap(boot_A);
    NumericMatrix boot_B_R = Rcpp::wrap(boot_B);
    
    // Compute jackknife influences for bootstrap samples
    NumericVector jack_A = jis_jackknife_influences_cpp(boot_A_R, q, normalize, log_base, 
                                                        pseudocount, n_transcripts_fixed);
    NumericVector jack_B = jis_jackknife_influences_cpp(boot_B_R, q, normalize, log_base, 
                                                        pseudocount, n_transcripts_fixed);
    
    // BUG FIX #8: CRITICAL - Add defensive size validation before array access
    // Jackknife can return empty vectors on error, which would cause out-of-bounds access
    if ((int)jack_A.size() != n_tx || (int)jack_B.size() != n_tx) {
      // Skip this bootstrap iteration if sizes don't match
      continue;
    }
    
    for (int i = 0; i < n_tx; i++) {
      if (!std::isnan(jack_A[i]) && !std::isnan(jack_B[i])) {
        bootstrap_deltas(b, i) = jack_A[i] - jack_B[i];
      }
    }
  }
  
  PutRNGstate();
  
  // Compute confidence intervals and p-values
  double alpha = 1.0 - confidence;  // Use confidence parameter
  NumericVector ci_lower(n_tx);
  NumericVector ci_upper(n_tx);
  NumericVector pvalues(n_tx);
  NumericVector variances(n_tx);
  NumericVector effect_sizes(n_tx);
  NumericVector ci_widths(n_tx);
  NumericVector rel_ci_widths(n_tx);
  
  for (int i = 0; i < n_tx; i++) {
    arma::vec col_deltas = bootstrap_deltas.col(i);
    arma::vec valid_deltas = col_deltas(arma::find_finite(col_deltas));
    
    if (valid_deltas.n_elem < 2) {
      ci_lower[i] = NA_REAL;
      ci_upper[i] = NA_REAL;
      pvalues[i] = NA_REAL;
      variances[i] = NA_REAL;
      effect_sizes[i] = NA_REAL;
      ci_widths[i] = NA_REAL;
      rel_ci_widths[i] = NA_REAL;
      continue;
    }
    
    // Sort for quantile computation
    arma::vec sorted_deltas = arma::sort(valid_deltas);
    // BUG FIX #13: Add explicit check for size truncation before int conversion
    size_t n_valid_size = sorted_deltas.n_elem;
    if (n_valid_size > static_cast<size_t>(INT_MAX)) {
      Rcpp::warning("Bootstrap sample size exceeds integer maximum, using capped value");
    }
    int n_valid = static_cast<int>(std::min(n_valid_size, static_cast<size_t>(2000000000)));  // Avoid int overflow
    
    // BUG FIX #9: CRITICAL - Safely compute quantile indices with overflow protection
    // Use nearest-rank method for quantile: ceil(p * n) - 1 (for 0-based indexing)
    // This is the standard method used by R's quantile() function with type=1
    // Prevent integer overflow when computing index positions
    double lower_pos = alpha / 2.0;  // As fraction instead of absolute position
    double upper_pos = 1.0 - alpha / 2.0;
    
    // Compute indices safely: rank = max(1, ceil(p * n))  [1-based], then -1 for 0-based
    int lower_rank = std::max(1, (int)std::ceil(lower_pos * n_valid));
    int upper_rank = std::max(1, (int)std::ceil(upper_pos * n_valid));
    
    int lower_idx = std::min(lower_rank - 1, n_valid - 1);  // Convert to 0-based and bound
    int upper_idx = std::min(upper_rank - 1, n_valid - 1);  // Convert to 0-based and bound
    lower_idx = std::max(0, lower_idx);
    upper_idx = std::max(0, upper_idx);
    
    ci_lower[i] = sorted_deltas[lower_idx];
    ci_upper[i] = sorted_deltas[upper_idx];
    
    // Compute mean and variance
    double mean_delta = arma::mean(valid_deltas);
    double var_delta = arma::var(valid_deltas);
    variances[i] = var_delta;
    
    // BUG FIX: CI width and relative CI width with correct numerical stability
    // Relative CI: CI_width / |mean_delta| (prevent division by zero with small epsilon)
    double width = ci_upper[i] - ci_lower[i];
    ci_widths[i] = width;
    double abs_mean_delta = std::abs(mean_delta);
    rel_ci_widths[i] = (width > 0) ? width / (abs_mean_delta + 1e-10) : width;
    
    // Effect size (normalized mean delta)
    effect_sizes[i] = std::abs(mean_delta);
    
    // Two-tailed p-value: proportion of bootstrap samples as/more extreme than observed
    double obs_abs = std::abs(delta_influence[i]);
    int extreme_count = 0;
    for (size_t j = 0; j < n_valid_size; j++) {  // BUG FIX: Use size_t loop to match array type
      if (std::abs(sorted_deltas[j]) >= obs_abs) {
        extreme_count++;
      }
    }
    // Divide by n_valid (number of valid bootstrap replicates actually computed)
    // BUG FIX #12: Use explicit float division to avoid integer division
    // Set minimum p-value as 1/n_valid to avoid spurious zero p-values
    pvalues[i] = std::max(1.0 / n_valid, (double)extreme_count / (double)n_valid);
  }
  
  return List::create(
    Named("delta_influence") = delta_influence,
    Named("variance") = variances,
    Named("ci_lower") = ci_lower,
    Named("ci_upper") = ci_upper,
    Named("p_value") = pvalues,
    Named("effect_size") = effect_sizes,
    Named("ci_width") = ci_widths,
    Named("relative_ci_width") = rel_ci_widths
  );
}

// ============================================================================
// DIVERGENCE BOOTSTRAP COMPUTATION (C++ OPTIMIZATION)
// ============================================================================
// Purpose: Accelerate divergence bootstrap by computing Tsallis divergence for each
//          replicate in C++ rather than via R loops.
// Target speedup: 10-15x compared to pure R implementation (nested loops eliminated)
//
// Supports:
//   1. Independent bootstrap: resample x and y independently
//   2. Paired/block bootstrap: resample pairs as units while maintaining correlation

// Internal helper: Compute Tsallis divergence between two probability distributions
// [[Rcpp::export(rng = false)]]
double tsallis_divergence_cpp(NumericVector p, NumericVector r, double q = 1.0, 
                              double log_base = 2.718281828) {
  // Validate inputs
  // BUG FIX: Allow different-sized vectors; compute over min(p.size(), r.size())
  if (log_base <= 0 || std::abs(log_base - 1.0) < 1e-10) {
    return NA_REAL;
  }
  
  int n = std::min(static_cast<int>(p.size()), static_cast<int>(r.size()));
  if (n == 0) return NA_REAL;
  
  double divergence = 0.0;
  double q_tol = 1e-6;
  
  // Special case: q ≈ 0 (Tsallis at q=0 always equals 0)
  if (q < q_tol) {
    // D_0(p||r) = (1/(0-1)) * (1 - sum(p^0 * r^1))
    //           = (-1) * (1 - 1) = 0
    return 0.0;
  } 
  // Special case: q ≈ 1 (KL divergence)
  else if (std::abs(q - 1.0) < q_tol) {
    // For KL: sum(p * log(p/r)) where we skip p=0 or r=0
    for (int i = 0; i < n; i++) {
      double p_i = p[i];
      double r_i = r[i];
      // Only compute where both p and r are positive
      if (p_i > 1e-15 && r_i > 1e-15) {
        divergence += p_i * std::log(p_i / r_i) / std::log(log_base);
      }
    }
  } else if (q > 0) {
    // General Tsallis divergence (Furuichi formula)
    // D_q(p||r) = (1/(q-1)) * (1 - sum(p^q * r^(1-q)))
    
    double sum_pq_r = 0.0;
    bool has_valid_term = false;
    
    for (int i = 0; i < n; i++) {
      double p_i = p[i];
      double r_i = r[i];
      
      // Skip zero/negative values
      if (p_i <= 1e-15 || r_i <= 1e-15) continue;
      
      has_valid_term = true;
      
      // Handle extreme values in log-space for numerical stability
      if (q > 2.0 || q < 0.5) {
        // Log-space computation for large |q|
        double log_term = q * std::log(std::max(p_i, 1e-10)) + 
                         (1.0 - q) * std::log(std::max(r_i, 1e-10));
        sum_pq_r += std::exp(log_term);
      } else {
        // Direct computation for q near 1
        sum_pq_r += std::pow(p_i, q) * std::pow(r_i, 1.0 - q);
      }
    }
    
    // If no valid terms found, divergence is 0 by convention
    // (distributions are orthogonal in the support)
    if (!has_valid_term) {
      divergence = 0.0;
    } else {
      divergence = (1.0 - sum_pq_r) / (q - 1.0);
    }
  } else {
    return NA_REAL;
  }
  
  // Ensure divergence is non-negative (mathematical property)
  divergence = std::abs(divergence);
  
  // Handle invalid results
  if (std::isnan(divergence) || std::isinf(divergence)) {
    return NA_REAL;
  }
  
  return divergence;
}

// Main divergence bootstrap function
// [[Rcpp::export]]
NumericVector divergence_bootstrap_compute_cpp(NumericVector x, NumericVector y, 
                                              int nboot = 1000, double q = 1.0,
                                              bool paired = false,
                                              double pseudocount = 0.0,
                                              double log_base = 2.718281828) {
  // Input validation
  if (log_base <= 0 || std::abs(log_base - 1.0) < 1e-10) {
    Rcpp::stop("Invalid log_base (must be > 1, not equal to 1)");
  }
  if (q < 0) {
    Rcpp::stop("Invalid q parameter (must be non-negative)");
  }
  if (nboot < 1) {
    Rcpp::stop("nboot must be at least 1");
  }
  
  int n_x = x.size();
  int n_y = y.size();
  
  if (n_x < 1 || n_y < 1) {
    Rcpp::stop("Input vectors x and y must have at least 1 element");
  }
  
  // Paired bootstrap requires equal-length vectors
  if (paired && n_x != n_y) {
    Rcpp::stop("For paired bootstrap, x and y must have equal length");
  }
  
  // Adjust counts with pseudocount (scalar)
  NumericVector x_adj = x + pseudocount;
  NumericVector y_adj = y + pseudocount;
  
  double total_x = sum(x_adj);
  double total_y = sum(y_adj);
  
  if (total_x <= 0 || total_y <= 0) {
    Rcpp::stop("Total count (after pseudocount) must be positive for both groups");
  }
  
  // Compute proportions from original data (for resampling)
  NumericVector p_hat = x_adj / total_x;
  NumericVector r_hat = y_adj / total_y;
  
  // Pre-allocate result vector
  NumericVector boot_divs(nboot);
  
  // Get R's RNG state for reproducibility
  GetRNGstate();
  
  // Main bootstrap loop
  for (int b = 0; b < nboot; b++) {
    NumericVector x_boot(n_x);
    NumericVector y_boot(n_y);
    
    if (paired) {
      // Paired/block bootstrap: resample pairs as units
      // For paired data, both x and y have n_x == n_y samples
      
      // Generate multinomial bootstrap samples from pair indices
      IntegerVector pair_indices(n_x);
      int total_int = n_x;  // Number of pairs
      
      // Resample from uniform distribution over pairs (equal weighting)
      NumericVector pair_probs(n_x, 1.0 / n_x);
      R::rmultinom(total_int, pair_probs.begin(), n_x, pair_indices.begin());
      
      // Now aggregate resampled pairs
      NumericVector x_boot_agg(n_x, 0.0);
      NumericVector y_boot_agg(n_y, 0.0);
      
      for (int i = 0; i < n_x; i++) {
        // pair_indices[i] tells us how many times to include pair i
        int pair_count = pair_indices[i];
        x_boot_agg[i] += x_adj[i] * pair_count;
        y_boot_agg[i] += y_adj[i] * pair_count;
      }
      
      x_boot = x_boot_agg;
      y_boot = y_boot_agg;
    } else {
      // Independent bootstrap: resample x and y independently
      IntegerVector x_boot_int(n_x);
      IntegerVector y_boot_int(n_y);
      
      // BUG FIX: Check totals against INT_MAX before conversion to avoid overflow
      if (total_x > static_cast<double>(INT_MAX) || total_y > static_cast<double>(INT_MAX)) {
        Rcpp::stop("Total count exceeds maximum integer value (%d)", INT_MAX);
      }
      int total_x_int = static_cast<int>(std::round(total_x));
      int total_y_int = static_cast<int>(std::round(total_y));
      if (total_x_int <= 0) total_x_int = 1;
      if (total_y_int <= 0) total_y_int = 1;
      
      // Resample from multinomial distributions
      R::rmultinom(total_x_int, p_hat.begin(), n_x, x_boot_int.begin());
      R::rmultinom(total_y_int, r_hat.begin(), n_y, y_boot_int.begin());
      
      x_boot = as<NumericVector>(x_boot_int);
      y_boot = as<NumericVector>(y_boot_int);
    }
    
    // Normalize bootstrap samples to probability distributions
    double x_boot_sum = sum(x_boot);
    double y_boot_sum = sum(y_boot);
    
    // Safety check: avoid division by zero
    if (x_boot_sum <= 0 || y_boot_sum <= 0) {
      boot_divs[b] = NA_REAL;
      continue;
    }
    
    NumericVector p_boot = x_boot / x_boot_sum;
    NumericVector r_boot = y_boot / y_boot_sum;
    
    // Compute Tsallis divergence for this bootstrap sample
    boot_divs[b] = tsallis_divergence_cpp(p_boot, r_boot, q, log_base);
  }
  
  PutRNGstate();
  
  return boot_divs;
}

// ============================================================================
// Paired Divergence Bootstrap with Explicit Pair Structure
// ============================================================================
// [[Rcpp::export]]
NumericVector divergence_bootstrap_paired_cpp(
    NumericVector x,          // Control group counts (length = n_pairs)
    NumericVector y,          // Treatment group counts (length = n_pairs)
    IntegerVector pair_ids,   // Pair identifiers (length = n_pairs) - just for validation
    int nboot = 1000,
    double q = 1.0,
    double pseudocount = 0.0,
    double log_base = 2.718281828) {
  
  // Input validation
  if (log_base <= 0 || std::abs(log_base - 1.0) < 1e-10) {
    Rcpp::stop("Invalid log_base (must be > 1, not equal to 1)");
  }
  if (q < 0) {
    Rcpp::stop("Invalid q parameter (must be non-negative)");
  }
  if (nboot < 1) {
    Rcpp::stop("nboot must be at least 1");
  }
  
  int n_pairs = x.size();
  if (n_pairs < 2) {
    Rcpp::stop("paired bootstrap requires at least 2 pairs");
  }
  if (y.size() != n_pairs || pair_ids.size() != n_pairs) {
    Rcpp::stop("x, y, and pair_ids must have the same length");
  }
  
  // Adjust counts with pseudocount
  NumericVector x_adj = x + pseudocount;
  NumericVector y_adj = y + pseudocount;
  
  // Pre-allocate result vector
  NumericVector boot_divs(nboot);
  
  // BUG FIX #10: Pre-allocate vectors OUTSIDE bootstrap loop to avoid redundant allocation
  // These will be reset on each iteration instead of recreated
  IntegerVector pair_counts(n_pairs);
  NumericVector pair_probs(n_pairs, 1.0 / n_pairs);
  NumericVector x_boot(n_pairs);
  NumericVector y_boot(n_pairs);
  
  // Get R's RNG state for reproducibility
  GetRNGstate();
  
  // Main bootstrap loop: resample pairs with replacement
  for (int b = 0; b < nboot; b++) {
    // Reset vectors instead of reallocating
    std::fill(x_boot.begin(), x_boot.end(), 0.0);
    std::fill(y_boot.begin(), y_boot.end(), 0.0);
    std::fill(pair_counts.begin(), pair_counts.end(), 0);
    
    // Multinomial resampling of pair indices
    // Each pair is selected with equal probability
    // rmultinom(n, prob, K, counts) samples n items from K categories with probabilities prob
    R::rmultinom(n_pairs, pair_probs.begin(), n_pairs, pair_counts.begin());
    
    // Aggregate counts from resampled pairs
    double x_boot_sum = 0.0;
    double y_boot_sum = 0.0;
    
    // For each original pair, add its contribution to the bootstrap sample
    for (int i = 0; i < n_pairs; i++) {
      int resample_count = pair_counts[i];
      
      if (resample_count > 0) {
        // This pair was selected resample_count times in the bootstrap
        // Add its adjusted counts to the bootstrap sample
        x_boot[i] = x_adj[i] * resample_count;
        y_boot[i] = y_adj[i] * resample_count;
        
        x_boot_sum += x_boot[i];
        y_boot_sum += y_boot[i];
      }
    }
    
    // Safety check: avoid division by zero
    if (x_boot_sum <= 1e-10 || y_boot_sum <= 1e-10) {
      boot_divs[b] = NA_REAL;
      continue;
    }
    
    // Normalize to probability distributions
    NumericVector p_boot = x_boot / x_boot_sum;
    NumericVector r_boot = y_boot / y_boot_sum;
    
    // Compute Tsallis divergence for this bootstrap sample
    boot_divs[b] = tsallis_divergence_cpp(p_boot, r_boot, q, log_base);
  }
  
  PutRNGstate();
  
  return boot_divs;
}

// =============================================================================
// ENHANCED: Flexible paired/unpaired divergence bootstrap (MARCH 2026)
// =============================================================================
// Handles complete pairs, incomplete pairs, and unpaired samples
// Supports arbitrary mixing of paired and unpaired data
// 
// Strategy:
//   1. Identify complete pairs (pair_id present in both x AND y)
//   2. Identify unpaired x samples (x_pair_id = 0 or not in y_pair_ids)
//   3. Identify unpaired y samples (y_pair_id = 0 or not in x_pair_ids)
//   4. For each bootstrap:
//      - Resample complete pairs as units (preserve correlation)
//      - Resample unpaired x independently
//      - Resample unpaired y independently
//   5. Compute divergence on aggregated counts

// [[Rcpp::export(rng = false)]]
NumericVector divergence_bootstrap_flexible_cpp(
    NumericVector x,           // All control samples (any length)
    NumericVector y,           // All treatment samples (any length)
    IntegerVector x_pair_ids,  // Pair IDs for x (0/NA = unpaired)
    IntegerVector y_pair_ids,  // Pair IDs for y (0/NA = unpaired)
    int nboot = 1000,
    double q = 1.0,
    double pseudocount = 0.0,
    double log_base = 2.718281828) {
  
  // Input validation
  if (log_base <= 0 || std::abs(log_base - 1.0) < 1e-10) {
    Rcpp::stop("Invalid log_base (must be > 1, not equal to 1)");
  }
  if (q < 0) {
    Rcpp::stop("Invalid q parameter (must be non-negative)");
  }
  if (nboot < 1) {
    Rcpp::stop("nboot must be at least 1");
  }
  if (x.size() != x_pair_ids.size() || y.size() != y_pair_ids.size()) {
    Rcpp::stop("x and x_pair_ids must have same length; y and y_pair_ids must have same length");
  }
  
  int nx = x.size();
  int ny = y.size();
  
  if (nx == 0 || ny == 0) {
    Rcpp::stop("x and y must have at least one sample each");
  }
  
  // Adjust counts with pseudocount
  NumericVector x_adj = x + pseudocount;
  NumericVector y_adj = y + pseudocount;
  
  // Step 2: Index samples by pair membership (use std::vector for reliable resizing)
  std::map<int, std::vector<int>> x_indices_by_pair;  // pair_id -> indices in x
  std::map<int, std::vector<int>> y_indices_by_pair;  // pair_id -> indices in y
  
  std::vector<int> unpaired_x_indices;  // indices of unpaired x samples
  std::vector<int> unpaired_y_indices;  // indices of unpaired y samples
  
  // Categorize x samples
  for (int i = 0; i < nx; i++) {
    int pid = x_pair_ids[i];
    if (pid <= 0 || pid == NA_INTEGER) {
      unpaired_x_indices.push_back(i);
    } else {
      // Store index for this pair
      if (x_indices_by_pair.find(pid) == x_indices_by_pair.end()) {
        x_indices_by_pair[pid] = std::vector<int>();
      }
      x_indices_by_pair[pid].push_back(i);
    }
  }
  
  // Categorize y samples
  for (int i = 0; i < ny; i++) {
    int pid = y_pair_ids[i];
    if (pid <= 0 || pid == NA_INTEGER) {
      unpaired_y_indices.push_back(i);
    } else {
      // Store index for this pair
      if (y_indices_by_pair.find(pid) == y_indices_by_pair.end()) {
        y_indices_by_pair[pid] = std::vector<int>();
      }
      y_indices_by_pair[pid].push_back(i);
    }
  }
  
  // Step 3: Move incomplete pairs to unpaired pools
  // Incomplete = pair exists in only one group
  for (auto& item : x_indices_by_pair) {
    int pid = item.first;
    if (y_indices_by_pair.find(pid) == y_indices_by_pair.end()) {
      // This pair only exists in x - move to unpaired
      for (int idx : item.second) {
        unpaired_x_indices.push_back(idx);
      }
    }
  }
  
  for (auto& item : y_indices_by_pair) {
    int pid = item.first;
    if (x_indices_by_pair.find(pid) == x_indices_by_pair.end()) {
      // This pair only exists in y - move to unpaired
      for (int idx : item.second) {
        unpaired_y_indices.push_back(idx);
      }
    }
  }
  
  // Remove incomplete pairs from the pair maps
  std::vector<int> incomplete_x_pairs;
  for (auto& item : x_indices_by_pair) {
    if (y_indices_by_pair.find(item.first) == y_indices_by_pair.end()) {
      incomplete_x_pairs.push_back(item.first);
    }
  }
  for (int pid : incomplete_x_pairs) {
    x_indices_by_pair.erase(pid);
  }
  
  std::vector<int> incomplete_y_pairs;
  for (auto& item : y_indices_by_pair) {
    if (x_indices_by_pair.find(item.first) == x_indices_by_pair.end()) {
      incomplete_y_pairs.push_back(item.first);
    }
  }
  for (int pid : incomplete_y_pairs) {
    y_indices_by_pair.erase(pid);
  }
  
  // Step 4: Build list of remaining complete pairs
  // (now that incomplete pairs have been removed and moved to unpaired pools)
  std::vector<int> complete_pairs;
  for (auto& item : x_indices_by_pair) {
    int pid = item.first;
    // Should exist in y since we removed incomplete pairs
    if (y_indices_by_pair.find(pid) != y_indices_by_pair.end()) {
      complete_pairs.push_back(pid);
    }
  }
  
  // Pre-allocate result
  NumericVector boot_divs(nboot);
  
  GetRNGstate();
  
  // Main bootstrap loop
  for (int b = 0; b < nboot; b++) {
    double x_boot_sum = 0.0;
    double y_boot_sum = 0.0;
    std::vector<double> x_boot_counts(nx, 0.0);
    std::vector<double> y_boot_counts(ny, 0.0);
    bool skip_iteration = false;  // BUG FIX #15, #19: Flag to skip on errors
    
    // Step 3a: Resample complete pairs as units
    if (complete_pairs.size() > 0 && !skip_iteration) {
      // Generate random weights for sampling pairs with replacement
      std::vector<double> pair_probs(complete_pairs.size(), 1.0 / complete_pairs.size());
      std::vector<int> pair_counts(complete_pairs.size(), 0);
      // BUG FIX: Cast size_t to int for R::rmultinom, remove unsafe (int*) cast
      int n_complete_pairs = static_cast<int>(complete_pairs.size());
      R::rmultinom(n_complete_pairs, pair_probs.data(), 
                   n_complete_pairs, pair_counts.data());
      
      for (size_t p = 0; p < complete_pairs.size(); p++) {
        int pid = complete_pairs[p];
        int count = pair_counts[p];
        
        if (count > 0) {
          // BUG FIX #15: CRITICAL - Add defensive validation before map access
          // Using [] operator creates empty vectors if key doesn't exist
          if (x_indices_by_pair.find(pid) == x_indices_by_pair.end()) {
            Rcpp::warning("Pair ID %d not found in x_indices_by_pair (logic error in pair matching)", pid);
            continue;
          }
          if (y_indices_by_pair.find(pid) == y_indices_by_pair.end()) {
            Rcpp::warning("Pair ID %d not found in y_indices_by_pair (logic error in pair matching)", pid);
            continue;
          }
          
          // Add this pair's contribution (all samples for this pair_id)
          for (int idx : x_indices_by_pair.at(pid)) {
            // BUG FIX #19: CRITICAL - Cast count to double early to avoid integer overflow
            double count_d = static_cast<double>(count);
            double contribution = x_adj[idx] * count_d;
            
            // Safety check for numeric overflow
            if (!std::isfinite(contribution)) {
              Rcpp::warning("Numeric overflow in x pair contribution calculation for pair %d", pid);
              boot_divs[b] = NA_REAL;
              skip_iteration = true;
              break;
            }
            
            x_boot_counts[idx] = contribution;
            x_boot_sum += contribution;
          }
          for (int idx : y_indices_by_pair.at(pid)) {
            double count_d = static_cast<double>(count);
            double contribution = y_adj[idx] * count_d;
            
            if (!std::isfinite(contribution)) {
              Rcpp::warning("Numeric overflow in y pair contribution calculation for pair %d", pid);
              boot_divs[b] = NA_REAL;
              skip_iteration = true;
              break;
            }
            
            y_boot_counts[idx] = contribution;
            y_boot_sum += contribution;
          }
        }
      }
    }
    
    // Step 3b: Resample unpaired x samples independently
    if (unpaired_x_indices.size() > 0 && !skip_iteration) {
      std::vector<double> x_unp_probs(unpaired_x_indices.size(), 1.0 / unpaired_x_indices.size());
      std::vector<int> x_unp_counts(unpaired_x_indices.size(), 0);
      // BUG FIX: Cast size_t to int for R::rmultinom, remove unsafe (int*) cast
      int n_unp_x = static_cast<int>(unpaired_x_indices.size());
      R::rmultinom(n_unp_x, x_unp_probs.data(), 
                   n_unp_x, x_unp_counts.data());
      
      for (size_t u = 0; u < unpaired_x_indices.size(); u++) {
        int idx = unpaired_x_indices[u];
        int count = x_unp_counts[u];
        if (count > 0) {
          double contribution = x_adj[idx] * count;
          x_boot_counts[idx] = contribution;
          x_boot_sum += contribution;
        }
      }
    }
    
    // Step 3c: Resample unpaired y samples independently
    if (unpaired_y_indices.size() > 0 && !skip_iteration) {
      std::vector<double> y_unp_probs(unpaired_y_indices.size(), 1.0 / unpaired_y_indices.size());
      std::vector<int> y_unp_counts(unpaired_y_indices.size(), 0);
      // BUG FIX: Cast size_t to int for R::rmultinom, remove unsafe (int*) cast
      int n_unp_y = static_cast<int>(unpaired_y_indices.size());
      R::rmultinom(n_unp_y, y_unp_probs.data(), 
                   n_unp_y, y_unp_counts.data());
      
      for (size_t u = 0; u < unpaired_y_indices.size(); u++) {
        int idx = unpaired_y_indices[u];
        int count = y_unp_counts[u];
        if (count > 0) {
          double contribution = y_adj[idx] * count;
          y_boot_counts[idx] = contribution;
          y_boot_sum += contribution;
        }
      }
    }
    
    // Safety check: avoid division by zero
    if (x_boot_sum <= 1e-10 || y_boot_sum <= 1e-10) {
      boot_divs[b] = NA_REAL;
      continue;
    }
    
    // BUG FIX #15, #19: Skip divergence computation if we hit an error in pair resampling
    if (skip_iteration) {
      continue;
    }
    
    // Step 4: Normalize and compute divergence
    // BUG FIX: Use actual vector sizes without zero-padding
    // Divergence is computed only over overlapping indices
    NumericVector p_boot(nx);
    NumericVector r_boot(ny);
    
    for (int i = 0; i < nx; i++) {
      p_boot[i] = x_boot_counts[i] / x_boot_sum;
    }
    for (int i = 0; i < ny; i++) {
      r_boot[i] = y_boot_counts[i] / y_boot_sum;
    }
    
    boot_divs[b] = tsallis_divergence_cpp(p_boot, r_boot, q, log_base);
  }
  
  PutRNGstate();
  
  return boot_divs;
}
