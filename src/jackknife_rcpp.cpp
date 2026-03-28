// Jackknife resampling implementation in Rcpp + Armadillo
// 
// Purpose: Accelerate leave-one-out jackknife computation
// Target speedup: 2-3x compared to pure R implementation
// 
// This implementation uses Armadillo (via RcppArmadillo) for efficient matrix operations
// and vectorized computation to replace the R's row-by-row loop.

#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;

// Internal helper: Compute Shannon or Tsallis entropy for proportions
// [[Rcpp::export(rng = false)]]
double entropy_cpp(NumericVector p, double q = 1.0, bool normalize = true, double log_base = 2.718281828) {
  // Remove NA values
  LogicalVector not_na = !is_na(p);
  NumericVector p_clean = p[not_na];
  
  if (p_clean.size() == 0) return NA_REAL;
  
  // Filter out zeros and negative values
  LogicalVector valid = (p_clean > 1e-10);
  NumericVector p_valid = p_clean[valid];
  
  if (p_valid.size() == 0) return NA_REAL;
  
  double entropy = 0.0;
  double q_tol = 1e-6;
  
  // Special case for q ≈ 0 (species richness)
  if (q < q_tol) {
    int n = p_valid.size();
    entropy = (std::log(n) - 1.0) / std::log(log_base);
    return entropy;
  }
  
  if (std::abs(q - 1.0) < q_tol) {
    // Shannon entropy (q = 1)
    for (int i = 0; i < p_valid.size(); i++) {
      double pi = p_valid[i];
      if (pi > 0) {
        entropy -= pi * std::log(pi) / std::log(log_base);
      }
    }
  } else {
    // Tsallis entropy (q != 1): (1 - sum(p^q)) / (q - 1)
    double sum_pq = 0.0;
    for (int i = 0; i < p_valid.size(); i++) {
      sum_pq += std::pow(p_valid[i], q);
    }
    
    entropy = (1.0 - sum_pq) / ((q - 1.0) * std::log(log_base));
  }
  
  // Normalize by maximum entropy if requested
  if (normalize) {
    int n = p_valid.size();
    double max_entropy;
    
    if (q < q_tol) {
      // Species richness: max = (log(n) - 1)
      max_entropy = (std::log(n) - 1.0) / std::log(log_base);
    } else if (std::abs(q - 1.0) < q_tol) {
      // Shannon: max = log(n)
      max_entropy = std::log(n) / std::log(log_base);
    } else {
      // Tsallis: max = (1 - n^(1-q)) / (q - 1)
      max_entropy = (1.0 - std::pow(n, 1.0 - q)) / ((q - 1.0) * std::log(log_base));
    }
    
    if (max_entropy > 0 && std::isfinite(max_entropy)) {
      entropy = entropy / max_entropy;
    }
    // If max_entropy <= 0 or non-finite, skip normalization (entropy stays unnormalized)
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
      if (pi > 0) {
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
    
    // Compute D_q = (sum_pq)^(1/(1-q))
    double exponent = 1.0 / (1.0 - q);
    return std::pow(sum_pq, exponent);
  }
}
// [[Rcpp::export]]
List jackknife_resampling_cpp(NumericMatrix counts, double q = 1.0, 
                              bool normalize = true, double log_base = 2.718281828,
                              double pseudocount = 0.0) {
  
  int n_obs = counts.nrow();
  int n_samples = counts.ncol();
  
  if (n_obs < 2) {
    Rcpp::warning("Insufficient observations for jackknife (need >= 2)");
    return List::create(Named("error") = "Insufficient observations");
  }
  
  // Convert to Armadillo matrix for efficient operations
  arma::mat counts_arma(counts.begin(), n_obs, n_samples, false);
  
  // Compute original full estimate
  arma::vec col_sums = arma::sum(counts_arma, 0).t();
  double total = arma::accu(counts_arma) + n_samples * pseudocount;
  
  if (total <= 0) {
    Rcpp::warning("Total count is zero or negative");
    return List::create(Named("error") = "Invalid total count");
  }
  
  arma::vec p_full = (col_sums + pseudocount) / total;
  double estimate = entropy_cpp(as<NumericVector>(wrap(p_full)), q, normalize, log_base);
  
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
    jackknife_estimates[i] = entropy_cpp(as<NumericVector>(wrap(p_minus_i)), q, normalize, log_base);
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
    R::rmultinom((int)(total), p_hat.begin(), n, boot_sample_int.begin());
    
    // Convert to NumericVector for proportion calculation
    NumericVector boot_counts = as<NumericVector>(boot_sample_int);
    
    // CRITICAL: Normalize bootstrap sample from counts to proportions
    // entropy_cpp() expects proportions, not counts
    double boot_total = sum(boot_counts);
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
  for (int b = 0; b < nboot; b++) {
    NumericVector boot_sample = boot_samples(_, b);
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
  
  // Adjust counts with pseudocount (scalar)
  NumericVector x_adj = x + pseudocount;
  double total = sum(x_adj);
  
  if (total <= 0) {
    Rcpp::stop("Total count (after pseudocount) must be positive");
  }
  
  // Pre-allocate result vector
  NumericVector boot_dist(nboot);
  
  // Create uniform probability vector for sampling pairs (each pair equally likely)
  NumericVector pair_probs(n_pairs);
  for (int i = 0; i < n_pairs; i++) {
    pair_probs[i] = 1.0 / n_pairs;
  }
  
  // Get R's RNG state for reproducibility
  GetRNGstate();
  
  // Main block bootstrap loop
  // For each bootstrap replicate, resample pairs with replacement
  for (int b = 0; b < nboot; b++) {
    // Initialize bootstrap sample for this replicate
    NumericVector boot_sample(n, 0.0);
    
    // Resample n_pairs pair indices with replacement using multinomial sampling
    // Create pairs dynamically by sampling pair indices one at a time
    IntegerVector pair_indices(n_pairs);
    for (int p = 0; p < n_pairs; p++) {
      // Sample from categorical distribution with uniform probs over pairs
      IntegerVector temp_count(n_pairs);
      R::rmultinom(1, pair_probs.begin(), n_pairs, temp_count.begin());
      // Find which category was sampled
      for (int k = 0; k < n_pairs; k++) {
        if (temp_count[k] == 1) {
          pair_indices[p] = k;
          break;
        }
      }
    }
    
    // Reconstruct bootstrap sample from sampled pair indices
    for (int p = 0; p < n_pairs; p++) {
      int sampled_pair_idx = pair_indices[p];
      int original_idx1 = 2 * sampled_pair_idx;
      int original_idx2 = 2 * sampled_pair_idx + 1;
      
      boot_sample[2 * p] = x[original_idx1];
      boot_sample[2 * p + 1] = x[original_idx2];
    }
    
    // Normalize bootstrap sample to proportions
    double boot_total = sum(boot_sample);
    NumericVector boot_props = boot_sample / boot_total;
    
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
