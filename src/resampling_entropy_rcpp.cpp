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
//   5. Utility functions: check_rcpp_available
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
  
  // Convert sizes from size_t to int safely (check against INT_MAX)
  if (p_valid.size() > static_cast<size_t>(INT_MAX)) {
    Rcpp::warning("Input size exceeds maximum integer value");
    return NA_REAL;
  }
  
  double entropy = 0.0;
  double q_tol = 1e-6;
  
  // Special case for q ≈ 0 (species richness)
  if (q < q_tol) {
    int n = p_valid.size();
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
    for (size_t i = 0; i < p_valid.size(); i++) {
      double pi = p_valid[i];
      if (pi > 1e-15) {  // BUG FIX: Use machine epsilon threshold not arbitrary 0
        entropy -= pi * std::log(pi) / std::log(log_base);
      }
    }
  } else {
    // Tsallis entropy (q != 1): (1 - sum(p^q)) / (q - 1)
    // NOTE: This formula does NOT include log_base
    double sum_pq = 0.0;
    for (size_t i = 0; i < p_valid.size(); i++) {
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
    for (size_t i = 0; i < p_valid.size(); i++) {
      double pi = p_valid[i];
      if (pi > 1e-15) {  // BUG FIX: Use machine epsilon threshold not arbitrary 0
        shannon -= pi * std::log(pi) / std::log(log_base);
      }
    }
    return std::pow(log_base, shannon);
  } else {
    // D_q = (Σp^q)^(1/(1-q))
    double sum_pq = 0.0;
    for (size_t i = 0; i < p_valid.size(); i++) {
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
    Rcpp::stop("Total count (after pseudocount) must be positive");
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
    int total_int = (int)std::round(total);  // BUG FIX: Use round() not cast
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
    Rcpp::stop("Total count (after pseudocount) must be positive");
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
    double h;
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
    
    // BUG FIX: Simplified size validation - jackknife always returns n_tx elements
    // Trust function contract: jack_A and jack_B have size == n_tx (from n_tx-row matrices)
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
    // BUG FIX: Use size_t to avoid integer cast overflow on huge arrays
    size_t n_valid_size = sorted_deltas.n_elem;
    int n_valid = static_cast<int>(std::min(n_valid_size, static_cast<size_t>(2000000000)));  // Avoid int overflow
    
    // Use nearest-rank method for quantile: ceil(p * n) - 1 (for 0-based indexing)
    // This is the standard method used by R's quantile() function with type=1
    // BUG FIX: Ensure indices stay in [0, n_valid-1] bounds
    double lower_pos = n_valid * alpha / 2.0;
    double upper_pos = n_valid * (1.0 - alpha / 2.0);
    int lower_idx = std::max(0, std::min((int)std::ceil(lower_pos) - 1, n_valid - 1));
    int upper_idx = std::max(0, std::min((int)std::ceil(upper_pos) - 1, n_valid - 1));
    
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
    // Set minimum p-value as 1/n_valid to avoid spurious zero p-values
    pvalues[i] = std::max(1.0 / n_valid, (double)extreme_count / n_valid);
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
