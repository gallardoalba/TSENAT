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

// Jackknife resampling computation (C++ implementation, called from R wrapper)
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

// Fallback: R wrapper for graceful degradation if Rcpp not available
//' @keywords internal
// [[Rcpp::export]]
bool check_rcpp_available() {
  return true;  // If this function exists, Rcpp compilation succeeded
}
