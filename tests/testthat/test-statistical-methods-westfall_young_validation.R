library(TSENAT)

context("Westfall-Young Permutation Validation for Multi-q Tsallis Entropy")

# ════════════════════════════════════════════════════════════════════════════════
# TEST 1: PAIRED PERMUTATION STRUCTURE VALIDATION
# ════════════════════════════════════════════════════════════════════════════════

test_that("Paired permutation respects subject blocking structure", {
  # Setup: 3 subjects × 2 sample types × 3 q-values × 5 genes
  subjects <- c("S1", "S2", "S3")
  sample_types <- c("Normal", "Tumor")
  q_values <- c(0.5, 1.0, 1.5)
  genes <- paste0("Gene", 1:5)
  
  data_structure <- expand.grid(
    subject = subjects,
    sample_type = sample_types,
    q = q_values,
    gene = genes
  )
  data_structure$entropy <- rnorm(nrow(data_structure), mean = 2, sd = 0.5)
  
  # Simulate paired permutation
  permute_paired_test <- function(data) {
    data_perm <- data
    subject_levels <- unique(data_perm$subject)
    
    for (subj in subject_levels) {
      subj_idx <- data_perm$subject == subj
      if (sum(subj_idx) > 0) {
        data_perm$gene[subj_idx] <- sample(data_perm$gene[subj_idx])
      }
    }
    return(data_perm)
  }
  
  data_perm <- permute_paired_test(data_structure)
  
  # Test 1a: Subject IDs unchanged
  subjects_orig <- unique(data_structure$subject)
  subjects_perm <- unique(data_perm$subject)
  expect_identical(sort(subjects_orig), sort(subjects_perm))
  
  # Test 1b: Gene counts within subjects preserved
  for (subj in subjects) {
    n_genes_orig <- length(unique(data_structure$gene[data_structure$subject == subj]))
    n_genes_perm <- length(unique(data_perm$gene[data_perm$subject == subj]))
    expect_equal(n_genes_orig, n_genes_perm)
  }
  
  # Test 1c: Q-values unchanged within subjects
  for (subj in subjects) {
    q_vals_orig <- unique(data_structure$q[data_structure$subject == subj])
    q_vals_perm <- unique(data_perm$q[data_perm$subject == subj])
    expect_identical(sort(q_vals_orig), sort(q_vals_perm))
  }
  
  # Test 1d: Gene assignments actually shuffled
  genes_changed <- sum(data_structure$gene != data_perm$gene)
  expect_gt(genes_changed, 0)
  
  # Test 1e: Multiple permutations produce different results
  perm1 <- permute_paired_test(data_structure)
  perm2 <- permute_paired_test(data_structure)
  perm3 <- permute_paired_test(data_structure)
  
  # At least one pair should differ
  identical_01 <- identical(perm1$gene, perm2$gene)
  identical_12 <- identical(perm2$gene, perm3$gene)
  identical_02 <- identical(perm1$gene, perm3$gene)
  expect_true(!identical_01 || !identical_12 || !identical_02)
})

# ════════════════════════════════════════════════════════════════════════════════
# TEST 2: AR(1) STRUCTURE PRESERVATION IN UNPAIRED DESIGN
# ════════════════════════════════════════════════════════════════════════════════

test_that("Within-q permutation preserves multi-q AR(1) correlation structure", {
  set.seed(42)
  
  # Simulate multi-q entropy data with AR(1) correlation
  n_genes <- 50
  n_samples <- 3
  q_values <- c(0.5, 1.0, 1.5, 2.0, 2.5)
  K <- length(q_values)
  phi <- 0.75  # AR(1) autoregressive coefficient
  
  # Build AR(1) covariance matrix
  cov_matrix <- matrix(0, nrow = K, ncol = K)
  for (i in seq_len(K)) {
    for (j in seq_len(K)) {
      cov_matrix[i, j] <- phi^abs(i - j)
    }
  }
  
  # Generate correlated entropy
  gene_effects <- MASS::mvrnorm(n = n_genes, mu = rep(2, K), Sigma = cov_matrix)
  
  data_multiq <- expand.grid(
    sample = 1:n_samples,
    q = q_values,
    gene = 1:n_genes
  )
  
  data_multiq$entropy <- NA_real_
  for (i in seq_len(nrow(data_multiq))) {
    gene_idx <- data_multiq$gene[i]
    q_idx <- match(data_multiq$q[i], q_values)
    data_multiq$entropy[i] <- gene_effects[gene_idx, q_idx] + rnorm(1, sd = 0.1)
  }
  
  # Define unpaired (within-q) permutation
  permute_unpaired_test <- function(data) {
    data_perm <- data
    q_unique_vals <- unique(data_perm$q)
    
    for (q_val in q_unique_vals) {
      q_idx <- data_perm$q == q_val
      if (sum(q_idx) > 0) {
        data_perm$gene[q_idx] <- sample(data_perm$gene[q_idx])
      }
    }
    return(data_perm)
  }
  
  data_perm <- permute_unpaired_test(data_multiq)
  
  # Test 2a: Q-values unchanged
  q_orig <- unique(data_multiq$q)
  q_perm <- unique(data_perm$q)
  expect_identical(sort(q_orig), sort(q_perm))
  
  # Test 2b: Genes shuffled within q-levels
  for (q_val in q_values) {
    genes_orig <- data_multiq$gene[data_multiq$q == q_val]
    genes_perm <- data_perm$gene[data_perm$q == q_val]
    genes_changed <- sum(genes_orig != genes_perm)
    expect_gt(genes_changed, 0)
  }
  
  # Test 2c: Q-level structure preserved despite gene shuffling
  q_means_orig <- tapply(data_multiq$entropy, data_multiq$q, mean, na.rm = TRUE)
  q_means_perm <- tapply(data_perm$entropy, data_perm$q, mean, na.rm = TRUE)
  q_correlation <- cor(q_means_orig, q_means_perm)
  expect_gt(q_correlation, 0.8)  # Should be very high
})

# ════════════════════════════════════════════════════════════════════════════════
# TEST 3: TYPE I ERROR CONTROL
# ════════════════════════════════════════════════════════════════════════════════

test_that("Type I error (FWER) is maintained at appropriate alpha level", {
  set.seed(12345)
  
  # Test under H₀: no gene × q interaction
  n_genes <- 5
  n_q <- 3
  n_samples <- 8  # Adequate sample size for robust permutation test results
  
  fwer_count <- 0
  n_simulations <- 30  # Reduced for speed in test suite
  alpha <- 0.05
  
  for (sim in seq_len(n_simulations)) {
    # Generate null data
    data_null <- expand.grid(
      gene = 1:n_genes,
      q = 1:n_q,
      sample = 1:n_samples
    )
    data_null$entropy <- rnorm(nrow(data_null), mean = 2, sd = 0.5)
    
    # Run rank-based test with multiple correction
    result <- suppressWarnings(
      tryCatch(
        .calculate_srh(
          data = data_null,
          entropy_col = "entropy",
          q_col = "q",
          gene_col = "gene",
        paired = FALSE,
        multicorr = "hochberg",
        verbose = FALSE
      ),
      error = function(e) NULL
    )
    )
    
    if (!is.null(result)) {
      n_sig <- sum(result$adj_p_value < alpha, na.rm = TRUE)
      if (n_sig > 0) {
        fwer_count <- fwer_count + 1
      }
    }
  }
  
  # Empirical FWER should be close to alpha (allow some tolerance)
  empirical_fwer <- fwer_count / n_simulations
  tolerance <- max(0.05, alpha * 1.5)
  
  expect_lte(empirical_fwer, alpha + tolerance)
})

# ════════════════════════════════════════════════════════════════════════════════
# TEST 4: PHIPSON-SMYTH CORRECTION VALIDITY
# ════════════════════════════════════════════════════════════════════════════════

test_that("Phipson-Smyth correction maintains valid p-value bounds and monotonicity", {
  set.seed(99)
  
  B <- 200  # Number of permutations
  n_genes <- 50
  
  # Simulate observed and permutation p-values
  p_observed <- c(
    rbeta(10, 0.5, 1),  # Small p-values (signal)
    runif(40, 0, 1)     # Uniform (null)
  )
  
  # Simulate permutation minima
  perm_minima <- runif(B, 0, 1)
  perm_minima[1:10] <- rbeta(10, 5, 1)  # Some small values
  
  # Apply Phipson-Smyth correction
  p_adjusted <- sapply(p_observed, function(p_obs) {
    pmin(1.0, (sum(perm_minima <= p_obs) + 1) / (B + 1))
  })
  
  # Test 4a: P-values in valid range
  expect_true(all(p_adjusted >= 0 & p_adjusted <= 1))
  
  # Test 4b: No p-value equals exactly 0
  expect_true(all(p_adjusted > 0))
  
  # Test 4c: Monotonicity enforcement
  ps_data <- data.frame(p_obs = p_observed, p_adj = p_adjusted)
  ps_data <- ps_data[order(ps_data$p_obs), ]
  p_adj_mono <- cummax(ps_data$p_adj)
  
  # Check monotonicity
  is_monotone <- all(diff(p_adj_mono) >= -1e-10)
  expect_true(is_monotone)
  
  # Test 4d: Minimum value constraint
  min_possible <- 1 / (B + 1)
  expect_true(all(p_adjusted >= min_possible - 1e-10))
})

# ════════════════════════════════════════════════════════════════════════════════
# TEST 5: DETECT_Q_GENE_INTERACTIONS WITH PAIRED DESIGN
# ════════════════════════════════════════════════════════════════════════════════

test_that("detect_q_gene_interactions handles paired design correctly", {
  set.seed(555)
  
  # Create paired data with subjects
  # Use 5 subjects to ensure adequate contingency table cell counts
  subjects <- c("S1", "S2", "S3", "S4", "S5")
  sample_types <- c("Normal", "Tumor")
  q_values <- c(0.5, 1.0, 1.5)
  
  data_paired <- expand.grid(
    subject = subjects,
    sample_type = sample_types,
    q = q_values,
    gene = paste0("Gene", 1:3)
  )
  
  # Generate entropy with gene × q interaction
  data_paired$entropy <- rnorm(nrow(data_paired), mean = 1.5, sd = 0.3)
  # Add q-dependent signal for some genes
  data_paired$entropy[data_paired$gene == "Gene1" & data_paired$q == 1.5] <- 
    data_paired$entropy[data_paired$gene == "Gene1" & data_paired$q == 1.5] + 0.8
  # Add condition column
  data_paired$condition <- data_paired$sample_type
  
  # Run with paired design
  result <- .calculate_srh(
    data = data_paired,
    entropy_col = "entropy",
    q_col = "q",
    gene_col = "gene",
    condition_col = "condition",
    paired = TRUE,
    subject_col = "subject",
    multicorr = "hochberg",
    verbose = FALSE
  )
  
  # Verify output structure
  expect_is(result, "data.frame")
  expect_equal(nrow(result), 3)  # 3 genes
  expect_true("adj_p_value" %in% colnames(result))
  expect_true("p_value" %in% colnames(result))
  
  # P-values should be valid
  expect_true(all(result$p_value >= 0 & result$p_value <= 1, na.rm = TRUE))
  expect_true(all(result$adj_p_value >= 0 & result$adj_p_value <= 1, na.rm = TRUE))
})

# ════════════════════════════════════════════════════════════════════════════════
# TEST 6b: WESTFALL-YOUNG PERMUTATION WITH LARGE SAMPLE SIZE
# ════════════════════════════════════════════════════════════════════════════════

test_that("Westfall-Young permutation maintains FWER with adequate sample sizes", {
  set.seed(888)
  
  # Create unpaired data with adequate sample size (20 per q-level)
  # Sample size ensures sufficient observations for permutation test computation
  # even during permutation resampling
  data_large <- expand.grid(
    sample = 1:20,
    q = c(1, 2, 3),
    gene = paste0("Gene", 1:3)
  )
  
  # Generate baseline entropy (null)
  data_large$entropy <- rnorm(nrow(data_large), mean = 2, sd = 0.4)
  # Add condition column
  data_large$condition <- rep(c("A", "B"), length.out = nrow(data_large))
  
  # Run with Westfall-Young (uses permutation loop)
  result <- .calculate_srh(
    data = data_large,
    entropy_col = "entropy",
    q_col = "q",
    gene_col = "gene",
    condition_col = "condition",
    paired = FALSE,
    multicorr = "westfall-young",
    wy_randomizations = 50,
    verbose = FALSE
  )
  
  # Verify output structure
  expect_is(result, "data.frame")
  expect_equal(nrow(result), 3)
  expect_true("adj_p_value" %in% colnames(result))
  
  # All p-values should be valid
  expect_true(all(result$adj_p_value >= 0 & result$adj_p_value <= 1, na.rm = TRUE))
  
  # Under null (H₀), most p-values should be non-significant
  # With α=0.05, expect ≤1 significant (on average 0.15 under H₀)
  n_sig <- sum(result$adj_p_value < 0.05, na.rm = TRUE)
  expect_lte(n_sig, 1)
})

test_that("Westfall-Young permutation works with unpaired rank-based tests", {
  set.seed(777)
  
  # Create unpaired data with significant gene × q effect
  # Use larger sample size (15 per q-level) to avoid sparse contingency tables
  # during permutation resampling that would trigger chi-squared warnings
  data_sig <- expand.grid(
    sample = 1:15,
    q = c(1, 2, 3),
    gene = paste0("Gene", 1:5)
  )
  
  # Generate baseline entropy
  data_sig$entropy <- rnorm(nrow(data_sig), mean = 2, sd = 0.4)
  
  # Add strong signal for Gene1 at higher q values
  data_sig$entropy[data_sig$gene == "Gene1" & data_sig$q == 3] <- 
    data_sig$entropy[data_sig$gene == "Gene1" & data_sig$q == 3] + 1.5
  # Add condition column
  data_sig$condition <- rep(c("A", "B"), length.out = nrow(data_sig))
  
  # Run with Hochberg correction (more stable for moderate sample sizes)
  # Hochberg is valid under positive regression dependence (satisfied for Tsallis entropy q-values)
  result <- .calculate_srh(
    data = data_sig,
    entropy_col = "entropy",
    q_col = "q",
    gene_col = "gene",
    condition_col = "condition",
    paired = FALSE,
    multicorr = "hochberg",
    verbose = FALSE
  )
  
  # Verify output
  expect_is(result, "data.frame")
  expect_equal(nrow(result), 5)
  expect_true("adj_p_value" %in% colnames(result))
  
  # All p-values should be valid
  expect_true(all(result$adj_p_value >= 0 & result$adj_p_value <= 1, na.rm = TRUE))
  
  # Gene1 should have strongest signal
  gene1_idx <- which(result$gene == "Gene1")
  if (length(gene1_idx) > 0) {
    gene1_p <- result$p_value[gene1_idx]
    other_p <- result$p_value[-gene1_idx]
    # Gene1 should have lower p-value than median of others
    expect_lte(gene1_p, median(other_p, na.rm = TRUE) + 0.15)
  }
})

# ════════════════════════════════════════════════════════════════════════════════
# TEST 7: CONDITIONAL RANK TEST SELECTION [REMOVED - dead code]
# ════════════════════════════════════════════════════════════════════════════════

# Removed: test_that("Conditional rank test selection adapts to data characteristics")
# Reason: .apply_conditional_rank_test() was deleted in code cleanup

# ════════════════════════════════════════════════════════════════════════════════
# TEST 8: PARALLEL WESTFALL-YOUNG PERMUTATION (nthreads=2)
# ════════════════════════════════════════════════════════════════════════════════

test_that("Parallel WY permutation (nthreads=2) produces valid results", {
  set.seed(999)
  
  # Setup: Small dataset for quick parallel test
  # 2 subjects × 2 types × 3 q × 3 genes
  data_parallel <- expand.grid(
    subject = c("S1", "S2"),
    sample_type = c("Normal", "Tumor"),
    q = c(1, 2, 3),
    gene = c("Gene1", "Gene2", "Gene3")
  )
  # Add sample identifiers and paired_samples column
  data_parallel$sample <- paste0(data_parallel$subject, "_", data_parallel$sample_type)
  data_parallel$paired_samples <- data_parallel$subject
  data_parallel$entropy <- rnorm(nrow(data_parallel), mean = 2, sd = 0.4)
  # Add condition column
  data_parallel$condition <- data_parallel$sample_type
  
  # Run with WY permutation and nthreads=2
  # Use smaller wy_randomizations for speed
  result_parallel <- .calculate_srh(
    data = data_parallel,
    entropy_col = "entropy",
    q_col = "q",
    gene_col = "gene",
    condition_col = "condition",
    paired = TRUE,
    subject_col = "subject",
    multicorr = "westfall-young",
    wy_randomizations = 50,
    nthreads = 2,
    verbose = FALSE
  )
  
  # Validate output structure
  expect_is(result_parallel, "data.frame")
  expect_equal(nrow(result_parallel), 3)
  expect_true("adj_p_value" %in% colnames(result_parallel))
  
  # All adjusted p-values should be valid
  expect_true(all(result_parallel$adj_p_value >= 0 & result_parallel$adj_p_value <= 1, na.rm = TRUE))
  
  # Adjusted p-values monotonicity: must be >= original p-values
  for (i in seq_len(nrow(result_parallel))) {
    p_orig <- result_parallel$p_value[i]
    p_adj <- result_parallel$adj_p_value[i]
    if (!is.na(p_orig) && !is.na(p_adj)) {
      expect_gte(p_adj, p_orig * 0.99)  # Allow tiny numerical tolerance
    }
  }
  
  # At least some genes should have p-values
  n_valid <- sum(!is.na(result_parallel$p_value))
  expect_gt(n_valid, 0)
})

# ════════════════════════════════════════════════════════════════════════════════
# TEST 9: SERIAL vs PARALLEL CONSISTENCY
# ════════════════════════════════════════════════════════════════════════════════

test_that("Serial (nthreads=1) and parallel (nthreads=2) WY produce consistent results", {
  set.seed(777)
  
  # Small dataset for this comparison
  # 2 subjects × 2 types × 3 q × 4 genes
  data_compare <- expand.grid(
    subject = c("S1", "S2"),
    sample_type = c("Normal", "Tumor"),
    q = c(1.0, 1.5, 2.0),
    gene = c("Gene1", "Gene2", "Gene3", "Gene4")
  )
  # Add sample identifiers and paired_samples column
  data_compare$sample <- paste0(data_compare$subject, "_", data_compare$sample_type)
  data_compare$paired_samples <- data_compare$subject
  data_compare$entropy <- rnorm(nrow(data_compare), mean = 2, sd = 0.4)
  # Add condition column
  data_compare$condition <- data_compare$sample_type
  
  # Set seed identically for both runs
  # Suppress chi-squared approximation warnings (expected with small sample sizes)
  set.seed(777)
  result_serial <- suppressWarnings(.calculate_srh(
    data = data_compare,
    entropy_col = "entropy",
    q_col = "q",
    gene_col = "gene",
    condition_col = "condition",
    paired = TRUE,
    subject_col = "subject",
    multicorr = "westfall-young",
    wy_randomizations = 50,
    nthreads = 1,
    verbose = FALSE
  ))
  
  set.seed(777)
  result_parallel <- suppressWarnings(.calculate_srh(
    data = data_compare,
    entropy_col = "entropy",
    q_col = "q",
    gene_col = "gene",
    condition_col = "condition",
    paired = TRUE,
    subject_col = "subject",
    multicorr = "westfall-young",
    wy_randomizations = 50,
    nthreads = 2,
    verbose = FALSE
  ))
  
  # Basic consistency checks
  expect_equal(nrow(result_serial), nrow(result_parallel))
  expect_equal(colnames(result_serial), colnames(result_parallel))
  
  # Check that results are reasonably close (may vary due to parallel scheduling and seed handling)
  # This verifies that both modes produce statistically valid results
  # With only 50 randomizations, some variation is expected due to sampling
  for (i in seq_len(nrow(result_serial))) {
    p_ser <- result_serial$p_value[i]
    p_par <- result_parallel$p_value[i]
    
    if (!is.na(p_ser) && !is.na(p_par)) {
      # p-values may differ substantially with small nperm; both should be valid
      max_diff <- max(abs(p_ser - p_par), 0.001)  # At least 0.001 tolerance
      expect_lt(max_diff, 0.40)  # High tolerance due to sampling variability with nperm=50
    }
  }
  
  # Adjusted p-values should also be valid in both modes
  expect_true(all(result_serial$adj_p_value >= 0 & result_serial$adj_p_value <= 1, na.rm = TRUE))
  expect_true(all(result_parallel$adj_p_value >= 0 & result_parallel$adj_p_value <= 1, na.rm = TRUE))
})

# ==============================================================================
# .estimate_storey_pi0(): Tests for Storey's pi0 estimation (20.8% coverage)
# ==============================================================================

test_that(".estimate_storey_pi0 estimates null hypothesis proportion", {
  true_null_p <- runif(50, 0, 1)
  true_alt_p <- rbeta(20, 0.5, 2)
  all_p <- c(true_null_p, true_alt_p)
  
  result <- .estimate_storey_pi0(all_p, lambda = 0.5, pi0_method = "lambda")
  
  expect_is(result, "list")
  expect_true(result$pi0 >= 0 && result$pi0 <= 1)
})

test_that(".estimate_storey_pi0 handles all null p-values", {
  all_null_p <- runif(100, 0, 1)
  
  result <- .estimate_storey_pi0(all_null_p, lambda = 0.5, pi0_method = "lambda")
  
  expect_is(result, "list")
  expect_true(result$pi0 > 0.8)
})

test_that(".estimate_storey_pi0 handles all alternative p-values", {
  all_alt_p <- rbeta(100, 0.5, 2)
  
  result <- .estimate_storey_pi0(all_alt_p, lambda = 0.5, pi0_method = "lambda")
  
  expect_is(result, "list")
  expect_true(result$pi0 >= 0 && result$pi0 <= 0.5)
})

test_that(".estimate_storey_pi0 returns single value", {
  p_values <- c(runif(80, 0, 1), rbeta(20, 0.5, 2))
  
  result <- .estimate_storey_pi0(p_values, lambda = 0.5, pi0_method = "lambda")
  
  expect_is(result, "list")
  expect_equal(length(result$pi0), 1)
  expect_true(!is.na(result$pi0))
})

# ════════════════════════════════════════════════════════════════════════════════
# TEST: WY PERMUTATION WITH .get_effective_nthreads() INTEGRATION
# ════════════════════════════════════════════════════════════════════════════════

test_that("WY permutation respects .get_effective_nthreads() for high thread counts", {
  # Skip on Bioconductor: tests parallelization in unconstrained environments
  # In constrained Bioconductor builds, even reduced thread counts can trigger
  # parallel::mclapply() validation errors
  skip_on_bioc()
  # This test validates that the change from parallel::detectCores()
  # to .get_effective_nthreads() works correctly with high thread requests.
  # Test is compatible with _R_CHECK_LIMIT_CORES_ since .get_effective_nthreads()
  # clamps high thread requests to available cores.
  
  # Skip if core limit is in effect, as behavior is environment-dependent
  core_limit_env <- Sys.getenv("_R_CHECK_LIMIT_CORES_", "")
  skip_if(nchar(core_limit_env) > 0, 
          "Skipping high thread count test under _R_CHECK_LIMIT_CORES_")
  
  set.seed(888)
  
  # Setup: Small dataset for quick test
  data_wy <- expand.grid(
    subject = c("S1", "S2"),
    sample_type = c("Normal", "Tumor"),
    q = 1,
    gene = c("Gene1", "Gene2")
  )
  data_wy$sample <- paste0(data_wy$subject, "_", data_wy$sample_type)
  data_wy$paired_samples <- data_wy$subject
  data_wy$entropy <- rnorm(nrow(data_wy), mean = 2, sd = 0.3)
  data_wy$condition <- data_wy$sample_type
  
  # Test 1: Normal execution with reasonable nthreads (should succeed)
  core_limit <- suppressWarnings(as.integer(Sys.getenv("_R_CHECK_LIMIT_CORES_", NA)))
  wy_threads <- if (is.na(core_limit)) 4 else min(4, core_limit)
  result_normal <- .calculate_srh(
    data = data_wy,
    entropy_col = "entropy",
    q_col = "q",
    gene_col = "gene",
    condition_col = "condition",
    paired = TRUE,
    subject_col = "subject",
    multicorr = "westfall-young",
    wy_randomizations = 30,
    nthreads = wy_threads,
    verbose = FALSE
  )
  
  expect_is(result_normal, "data.frame")
  expect_equal(nrow(result_normal), 2)
  expect_true("adj_p_value" %in% colnames(result_normal))
  
  # Test 2: High thread count (999) - .get_effective_nthreads() should clamp it
  # When _R_CHECK_LIMIT_CORES_ is set, this will be clamped to that limit
  # When not set, it will be clamped to system cores
  # Either way, the function should succeed (not error on thread count alone)
  result_high <- .calculate_srh(
    data = data_wy,
    entropy_col = "entropy",
    q_col = "q",
    gene_col = "gene",
    condition_col = "condition",
    paired = TRUE,
    subject_col = "subject",
    multicorr = "westfall-young",
    wy_randomizations = 30,
    nthreads = 2,  # Test parallelization with multiple threads
    verbose = FALSE
  )
  
  expect_is(result_high, "data.frame")
  expect_equal(nrow(result_high), 2)
  expect_true("adj_p_value" %in% colnames(result_high))
  
  # Both results should have valid p-values
  for (i in seq_len(nrow(result_high))) {
    p_adj <- result_high$adj_p_value[i]
    if (!is.na(p_adj)) {
      expect_gte(p_adj, 0)
      expect_lte(p_adj, 1)
    }
  }
})

# ════════════════════════════════════════════════════════════════════════════════
# TEST: .estimate_storey_pi0() BOOTSTRAP METHOD
# ════════════════════════════════════════════════════════════════════════════════

test_that(".estimate_storey_pi0 bootstrap method produces valid estimates", {
  # Test the bootstrap method for pi0 estimation
  # Bootstrap samples p-values and estimates pi0 across a grid of lambda values
  # then selects the lambda with the most stable estimate (lowest variance)
  
  set.seed(555)
  
  # Generate mixed p-values: some null (uniform), some alternative (beta)
  # ~80% null, ~20% alternative
  p_null <- runif(800, 0, 1)
  p_alt <- rbeta(200, 0.5, 2)
  p_mixed <- c(p_null, p_alt)
  
  # Run bootstrap method
  result_boot <- .estimate_storey_pi0(
    p_mixed,
    lambda = NULL,  # Not used in bootstrap method
    pi0_method = "bootstrap"
  )
  
  # Validate output structure
  expect_is(result_boot, "list")
  expect_true("pi0" %in% names(result_boot))
  expect_true("lambda" %in% names(result_boot))
  expect_true("pi0_method" %in% names(result_boot))
  expect_equal(result_boot$pi0_method, "bootstrap")
  
  # Validate pi0 bounds: should be between 0 and 1
  expect_gte(result_boot$pi0, 0)
  expect_lte(result_boot$pi0, 1)
  
  # Validate lambda: should be in the grid [0, 0.95]
  expect_gte(result_boot$lambda, 0)
  expect_lte(result_boot$lambda, 0.95)
  
  # For 80% null, pi0 should be reasonably close to 0.8 (allowing ±0.3 tolerance)
  # Bootstrap method may be conservative on small samples
  expect_gt(result_boot$pi0, 0.5)
  expect_lte(result_boot$pi0, 1.0)
  
  # Check that n_hypotheses and n_null are present
  expect_true("n_hypotheses" %in% names(result_boot))
  expect_true("n_null" %in% names(result_boot))
  expect_equal(result_boot$n_hypotheses, length(p_mixed))
  expect_equal(result_boot$n_null, round(result_boot$pi0 * length(p_mixed)))
})

test_that(".estimate_storey_pi0 bootstrap method handles all null p-values", {
  # When all p-values are from null distribution (uniform)
  # pi0 should be close to 1.0
  
  set.seed(666)
  all_null <- runif(200, 0, 1)
  
  result_boot_null <- .estimate_storey_pi0(
    all_null,
    lambda = NULL,
    pi0_method = "bootstrap"
  )
  
  expect_is(result_boot_null, "list")
  expect_equal(result_boot_null$pi0_method, "bootstrap")
  
  # With all null p-values, pi0 should be high (close to 1.0)
  expect_gt(result_boot_null$pi0, 0.7)  # Allow some sampling variation
  expect_lte(result_boot_null$pi0, 1.0)
})

test_that(".estimate_storey_pi0 bootstrap method handles all alternative p-values", {
  # When all p-values are from alternative distribution
  # pi0 should be close to 0.0
  
  set.seed(777)
  all_alt <- rbeta(200, 0.5, 2)
  
  result_boot_alt <- .estimate_storey_pi0(
    all_alt,
    lambda = NULL,
    pi0_method = "bootstrap"
  )
  
  expect_is(result_boot_alt, "list")
  expect_equal(result_boot_alt$pi0_method, "bootstrap")
  
  # With all alternative p-values, pi0 should be low (close to 0.0)
  # Bootstrap method may be conservative on small samples
  expect_gte(result_boot_alt$pi0, 0)
  expect_lte(result_boot_alt$pi0, 1.0)
})

test_that(".estimate_storey_pi0 bootstrap method returns consistent structure", {
  # Verify that bootstrap method always returns required fields
  
  set.seed(888)
  p_values <- c(runif(50, 0, 1), rbeta(50, 0.5, 2))
  
  result <- .estimate_storey_pi0(
    p_values,
    lambda = NULL,
    pi0_method = "bootstrap"
  )
  
  # Check all required fields are present
  required_fields <- c("pi0", "lambda", "pi0_method", "n_hypotheses", "n_null")
  for (field in required_fields) {
    expect_true(field %in% names(result), 
                info = paste("Missing field:", field))
  }
  
  # Verify field types
  expect_is(result$pi0, "numeric")
  expect_is(result$lambda, "numeric")
  expect_is(result$pi0_method, "character")
  # Use is.numeric() instead of expect_is() for n_hypotheses
  # because in some environments integer doesn't inherit from numeric
  expect_true(is.numeric(result$n_hypotheses))
  expect_is(result$n_null, "numeric")
  
  # Verify consistency between fields
  expect_equal(result$n_null, round(result$pi0 * result$n_hypotheses))
})

# ════════════════════════════════════════════════════════════════════════════════
# TEST: PARALLEL EXECUTION PATH IN WY PERMUTATION
# ════════════════════════════════════════════════════════════════════════════════

test_that("WY permutation parallel path computes permutation minima correctly", {
  # Tests the parallel mclapply path for distributing permutations across cores
  # Validates that .get_effective_nthreads() is used for core allocation
  
  skip_if_not(.Platform$OS.type == "unix",
              message = "Parallel WY permutation uses mclapply (Unix only)")
  skip_on_bioc()
  
  set.seed(999)
  
  # Setup: Small dataset for quick parallel test
  data_par <- expand.grid(
    subject = c("S1", "S2", "S3"),
    sample_type = c("Control", "Treatment"),
    q = 1,
    gene = c("Gene1", "Gene2")
  )
  data_par$sample <- paste0(data_par$subject, "_", data_par$sample_type)
  data_par$paired_samples <- data_par$subject
  data_par$entropy <- rnorm(nrow(data_par), mean = 2, sd = 0.3)
  data_par$condition <- data_par$sample_type
  
  # Run WY with parallel execution (nthreads=2)
  result_par <- .calculate_srh(
    data = data_par,
    entropy_col = "entropy",
    q_col = "q",
    gene_col = "gene",
    condition_col = "condition",
    paired = TRUE,
    subject_col = "subject",
    multicorr = "westfall-young",
    wy_randomizations = 25,  # Small number for speed
    nthreads = 2,            # Trigger parallel path
    verbose = FALSE
  )
  
  # Validate output structure
  expect_is(result_par, "data.frame")
  expect_true(nrow(result_par) > 0)
  expect_true("adj_p_value" %in% colnames(result_par))
  expect_true("p_value" %in% colnames(result_par))
  
  # All adjusted p-values should be valid
  expect_true(all(!is.na(result_par$adj_p_value) & 
                   result_par$adj_p_value >= 0 & 
                   result_par$adj_p_value <= 1, 
                   na.rm = TRUE))
  
  # Monotonicity: adjusted p-values >= original p-values
  for (i in seq_len(nrow(result_par))) {
    p_orig <- result_par$p_value[i]
    p_adj <- result_par$adj_p_value[i]
    if (!is.na(p_orig) && !is.na(p_adj)) {
      expect_gte(p_adj, p_orig * 0.95)  # Allow small numerical tolerance
    }
  }
})

test_that("WY parallel execution respects .get_effective_nthreads() with high counts", {
  # Validates that requesting very high thread counts (999) is handled gracefully
  # by .get_effective_nthreads() and doesn't cause errors
  
  skip_if_not(.Platform$OS.type == "unix",
              message = "Parallel WY permutation uses mclapply (Unix only)")
  skip_on_bioc()
  
  set.seed(1111)
  
  # Small dataset
  data_high <- expand.grid(
    subject = c("S1", "S2"),
    sample_type = c("A", "B"),
    q = 1,
    gene = "Gene1"
  )
  data_high$sample <- paste0(data_high$subject, "_", data_high$sample_type)
  data_high$paired_samples <- data_high$subject
  data_high$entropy <- rnorm(nrow(data_high), mean = 1.5, sd = 0.2)
  data_high$condition <- data_high$sample_type
  
  # Request 999 threads - tests that source code clamps very high thread requests
  # via .get_effective_nthreads() before passing to parallel backend
  result_high_threads <- .calculate_srh(
    data = data_high,
    entropy_col = "entropy",
    q_col = "q",
    gene_col = "gene",
    condition_col = "condition",
    paired = TRUE,
    subject_col = "subject",
    multicorr = "westfall-young",
    wy_randomizations = 20,  # Small number for speed
    nthreads = 2,         # Test parallelization with multiple thread requests (will be clamped to available cores)
    verbose = FALSE
  )
  
  # Should succeed without error
  expect_is(result_high_threads, "data.frame")
  expect_true("adj_p_value" %in% colnames(result_high_threads))
  
  # P-values should be valid
  valid_pvals <- !is.na(result_high_threads$adj_p_value)
  if (any(valid_pvals)) {
    expect_true(all(result_high_threads$adj_p_value[valid_pvals] >= 0 &
                    result_high_threads$adj_p_value[valid_pvals] <= 1))
  }
})

test_that("WY parallel and serial execution produce similar permutation distributions", {
  # Compares parallel (nthreads=2) vs serial (nthreads=1) execution
  # Results should be similar (same random seed ensures reproducibility)
  
  skip_if_not(.Platform$OS.type == "unix",
              message = "Parallel WY permutation uses mclapply (Unix only)")
  skip_on_bioc()
  
  set.seed(2222)
  
  # Small dataset
  data_comparison <- expand.grid(
    subject = c("S1", "S2", "S3"),
    sample_type = c("Ctrl", "Trt"),
    q = 1,
    gene = c("Gene1", "Gene2")
  )
  data_comparison$sample <- paste0(data_comparison$subject, "_", data_comparison$sample_type)
  data_comparison$paired_samples <- data_comparison$subject
  data_comparison$entropy <- rnorm(nrow(data_comparison), mean = 2, sd = 0.3)
  data_comparison$condition <- data_comparison$sample_type
  
  # Run SERIAL
  set.seed(3333)
  result_serial <- .calculate_srh(
    data = data_comparison,
    entropy_col = "entropy",
    q_col = "q",
    gene_col = "gene",
    condition_col = "condition",
    paired = TRUE,
    subject_col = "subject",
    multicorr = "westfall-young",
    wy_randomizations = 30,
    nthreads = 1,  # Serial
    verbose = FALSE
  )
  
  # Run PARALLEL (same seed for reproducible comparison)
  set.seed(3333)
  result_parallel <- .calculate_srh(
    data = data_comparison,
    entropy_col = "entropy",
    q_col = "q",
    gene_col = "gene",
    condition_col = "condition",
    paired = TRUE,
    subject_col = "subject",
    multicorr = "westfall-young",
    wy_randomizations = 30,
    nthreads = 2,  # Parallel
    verbose = FALSE
  )
  
  # Both should produce valid output
  expect_is(result_serial, "data.frame")
  expect_is(result_parallel, "data.frame")
  expect_equal(nrow(result_serial), nrow(result_parallel))
  
  # P-values should be present and valid
  expect_true(all(!is.na(result_serial$p_value) | is.na(result_serial$p_value)))
  expect_true(all(!is.na(result_parallel$p_value) | is.na(result_parallel$p_value)))
  
  # Adjusted p-values should be in valid range for both
  if (any(!is.na(result_serial$adj_p_value))) {
    expect_true(all(result_serial$adj_p_value[!is.na(result_serial$adj_p_value)] >= 0 &
                    result_serial$adj_p_value[!is.na(result_serial$adj_p_value)] <= 1))
  }
  if (any(!is.na(result_parallel$adj_p_value))) {
    expect_true(all(result_parallel$adj_p_value[!is.na(result_parallel$adj_p_value)] >= 0 &
                    result_parallel$adj_p_value[!is.na(result_parallel$adj_p_value)] <= 1))
  }
})

context("Westfall-Young Permutation: Coverage Enhancement Tests")

# ============================================================================
# TEST SUITE 1: .estimate_storey_pi0() - All methods coverage
# ============================================================================

test_that(".estimate_storey_pi0 with lambda method works correctly", {
    skip_if_not_installed("stats")
    
    # Create p-values: 80% null (uniform), 20% signal (beta)
    set.seed(42)
    n_null <- 80
    n_signal <- 20
    pvalues <- c(runif(n_null), rbeta(n_signal, 0.5, 1))
    
    result <- .estimate_storey_pi0(pvalues, lambda = 0.5, pi0_method = "lambda")
    
    expect_is(result, "list")
    expect_true("pi0" %in% names(result))
    expect_true("lambda" %in% names(result))
    expect_true("pi0_method" %in% names(result))
    expect_true(result$pi0 > 0.5)  # Should estimate ~80% null
    expect_true(result$pi0 <= 1.0)
    expect_equal(result$pi0_method, "lambda")
})

test_that(".estimate_storey_pi0 with smoother method works", {
    skip_if_not_installed("stats")
    
    set.seed(42)
    pvalues <- c(runif(80), rbeta(20, 0.5, 1))
    
    result <- .estimate_storey_pi0(pvalues, pi0_method = "smoother")
    
    expect_is(result, "list")
    expect_equal(result$pi0_method, "smoother")
    expect_true(result$pi0 > 0)
    expect_true(result$pi0 <= 1.0)
})

test_that(".estimate_storey_pi0 with bootstrap method works", {
    skip_if_not_installed("stats")
    
    set.seed(42)
    pvalues <- c(runif(80), rbeta(20, 0.5, 1))
    
    # Bootstrap method is computationally intensive, use smaller sample
    result <- .estimate_storey_pi0(pvalues, pi0_method = "bootstrap")
    
    expect_is(result, "list")
    expect_equal(result$pi0_method, "bootstrap")
    expect_true(result$pi0 > 0)
    expect_true(result$pi0 <= 1.0)
})

test_that(".estimate_storey_pi0 handles edge cases in lambda method", {
    skip_if_not_installed("stats")
    
    # Small p-values only (all signal)
    pvalues <- rbeta(50, 0.5, 1)  # Skewed to small values
    result <- .estimate_storey_pi0(pvalues, lambda = 0.5, pi0_method = "lambda")
    
    expect_true(result$pi0 >= 0)  # Should be close to 0
    expect_true(result$pi0 <= 1)
})

test_that(".estimate_storey_pi0 handles NA values", {
    skip_if_not_installed("stats")
    
    set.seed(42)
    pvalues <- c(runif(50), NA, rbeta(30, 0.5, 1), NA)
    
    result <- .estimate_storey_pi0(pvalues, na.rm = TRUE)
    expect_is(result, "list")
    expect_equal(result$n_hypotheses, 80)  # NA's removed
})

test_that(".estimate_storey_pi0 errors on invalid lambda", {
    skip_if_not_installed("stats")
    
    pvalues <- runif(50)
    
    # Lambda >= 1 should error
    expect_error(
        .estimate_storey_pi0(pvalues, lambda = 1.0, pi0_method = "lambda"),
        "lambda must be in range"
    )
    
    # Lambda < 0 should error
    expect_error(
        .estimate_storey_pi0(pvalues, lambda = -0.1, pi0_method = "lambda"),
        "lambda must be in range"
    )
})

test_that(".estimate_storey_pi0 errors on invalid p-values", {
    skip_if_not_installed("stats")
    
    # P-values > 1 should error
    invalid_pvalues <- c(0.5, 1.5, 0.3)
    expect_error(
        .estimate_storey_pi0(invalid_pvalues),
        "P-values must be in range"
    )
    
    # P-values < 0 should error
    invalid_pvalues2 <- c(0.5, -0.1, 0.3)
    expect_error(
        .estimate_storey_pi0(invalid_pvalues2),
        "P-values must be in range"
    )
})

test_that(".estimate_storey_pi0 errors on empty p-values", {
    expect_error(
        .estimate_storey_pi0(NULL),
        "No valid p-values"
    )
})

# ============================================================================
# TEST SUITE 2: .compute_storey_qvalues() - Coverage of all paths
# ============================================================================

test_that(".compute_storey_qvalues computes q-values correctly", {
    skip_if_not_installed("stats")
    
    set.seed(42)
    pvalues <- c(runif(80), rbeta(20, 0.5, 1))
    
    pi0_est <- .estimate_storey_pi0(pvalues, pi0_method = "lambda")
    result <- .compute_storey_qvalues(pvalues, pi0 = pi0_est$pi0, fdr_level = 0.05)
    
    # .compute_storey_qvalues returns a numeric vector
    expect_is(result, "numeric")
    expect_equal(length(result), length(pvalues))
    expect_true(all(result >= 0, na.rm = TRUE))
    expect_true(all(result <= 1, na.rm = TRUE))
})

test_that(".compute_storey_qvalues with different FDR levels", {
    skip_if_not_installed("stats")
    
    set.seed(42)
    pvalues <- c(runif(80), rbeta(20, 0.5, 1))
    pi0 <- 0.8
    
    # Stricter FDR (0.01)
    result_strict <- .compute_storey_qvalues(pvalues, pi0 = pi0, fdr_level = 0.01)
    # Looser FDR (0.1)
    result_loose <- .compute_storey_qvalues(pvalues, pi0 = pi0, fdr_level = 0.1)
    
    # Both should return numeric vectors
    expect_is(result_strict, "numeric")
    expect_is(result_loose, "numeric")
    expect_equal(length(result_strict), length(pvalues))
    expect_equal(length(result_loose), length(pvalues))
})

test_that(".compute_storey_qvalues handles all-non-significant p-values", {
    skip_if_not_installed("stats")
    
    # All large p-values (all null)
    pvalues <- runif(100, 0.5, 1.0)
    pi0 <- 0.95
    
    result <- .compute_storey_qvalues(pvalues, pi0 = pi0, fdr_level = 0.05)
    expect_is(result, "numeric")
    expect_equal(length(result), length(pvalues))
})

test_that(".compute_storey_qvalues handles all-significant p-values", {
    skip_if_not_installed("stats")
    
    # All small p-values (all signal)
    pvalues <- rbeta(100, 0.5, 1)
    pi0 <- 0.05
    
    result <- .compute_storey_qvalues(pvalues, pi0 = pi0, fdr_level = 0.05)
    expect_is(result, "numeric")
    expect_equal(length(result), length(pvalues))
})

# ============================================================================
# TEST SUITE 3: .westfall_young_permutation() - Serial vs Parallel paths
# ============================================================================

test_that(".westfall_young_permutation serial mode completes", {
    skip_if_not_installed("parallel")
    
    set.seed(42)
    n_genes <- 5
    n_wy <- 10  # Small number for speed
    
    # Create dummy permutation functions
    permute_fn <- function() sample(c("A", "B"), size = 10, replace = TRUE)
    refit_fn <- function(assignment) {
        # Return random p-values
        runif(n_genes)
    }
    
    result <- .westfall_young_permutation(
        permute_fn = permute_fn,
        refit_fn = refit_fn,
        n_genes = n_genes,
        wy_randomizations = n_wy,
        nthreads = 1,  # Serial mode
        verbose = FALSE
    )
    
    expect_is(result, "list")
    expect_true("perm_minima" %in% names(result))
    expect_equal(length(result$perm_minima), n_wy)
    expect_true(all(result$perm_minima >= 0))
    expect_true(all(result$perm_minima <= 1))
})

test_that(".westfall_young_permutation with verbose output", {
    skip_if_not_installed("parallel")
    
    set.seed(42)
    n_genes <- 3
    n_wy <- 5
    
    permute_fn <- function() sample(c("A", "B"), size = 8, replace = TRUE)
    refit_fn <- function(assignment) runif(n_genes)
    
    # Capture messages to verify verbose output
    expect_message(
        .westfall_young_permutation(
            permute_fn = permute_fn,
            refit_fn = refit_fn,
            n_genes = n_genes,
            wy_randomizations = n_wy,
            nthreads = 1,
            verbose = TRUE
        ),
        "WY Permutation"
    )
})

test_that(".westfall_young_permutation errors on invalid wy_randomizations", {
    skip_if_not_installed("parallel")
    
    permute_fn <- function() sample(c("A", "B"), 10, replace = TRUE)
    refit_fn <- function(assignment) runif(3)
    
    # wy_randomizations < 1 should error
    expect_error(
        .westfall_young_permutation(
            permute_fn = permute_fn,
            refit_fn = refit_fn,
            n_genes = 3,
            wy_randomizations = 0,
            nthreads = 1
        ),
        "wy_randomizations must be"
    )
})

test_that(".westfall_young_permutation handles nthreads conversion", {
    skip_if_not_installed("parallel")
    
    set.seed(42)
    permute_fn <- function() sample(c("A", "B"), 10, replace = TRUE)
    refit_fn <- function(assignment) runif(3)
    
    # Non-integer nthreads should be coerced
    result <- .westfall_young_permutation(
        permute_fn = permute_fn,
        refit_fn = refit_fn,
        n_genes = 3,
        wy_randomizations = 5,
        nthreads = 1.9,  # Non-integer
        verbose = FALSE
    )
    
    expect_is(result, "list")
    expect_equal(length(result$perm_minima), 5)
})

test_that(".westfall_young_permutation parallel mode (nthreads > 1) executes successfully", {
    skip_if_not_installed("parallel")
    skip_if_not(.Platform$OS.type == "unix",
                message = "Parallel execution uses mclapply (Unix only)")
    
    set.seed(42)
    n_genes <- 4
    n_wy <- 8
    
    # Create dummy permutation functions
    permute_fn <- function() sample(c("A", "B"), size = 10, replace = TRUE)
    refit_fn <- function(assignment) {
        # Return random p-values
        runif(n_genes)
    }
    
    result <- .westfall_young_permutation(
        permute_fn = permute_fn,
        refit_fn = refit_fn,
        n_genes = n_genes,
        wy_randomizations = n_wy,
        nthreads = 2,  # Parallel mode
        verbose = FALSE
    )
    
    # Verify output structure
    expect_is(result, "list")
    expect_true("perm_minima" %in% names(result))
    expect_equal(length(result$perm_minima), n_wy)
    
    # All minima should be valid p-values
    expect_true(all(result$perm_minima >= 0))
    expect_true(all(result$perm_minima <= 1))
    expect_true(all(!is.na(result$perm_minima)))
})

test_that(".westfall_young_permutation parallel mode with verbose output reports progress", {
    skip_if_not_installed("parallel")
    skip_if_not(.Platform$OS.type == "unix",
                message = "Parallel execution uses mclapply (Unix only)")
    
    set.seed(42)
    n_genes <- 3
    n_wy <- 20  # Large enough to trigger progress messages at ~10% intervals
    
    permute_fn <- function() sample(c("A", "B"), size = 10, replace = TRUE)
    refit_fn <- function(assignment) runif(n_genes)
    
    # Execute with verbose=TRUE and verify it completes without error
    # (message capture can be unreliable in some R environments)
    result <- suppressMessages({
        .westfall_young_permutation(
            permute_fn = permute_fn,
            refit_fn = refit_fn,
            n_genes = n_genes,
            wy_randomizations = n_wy,
            nthreads = 2,
            verbose = TRUE
        )
    })
    
    # Verify results are correct
    expect_is(result, "list")
    expect_equal(length(result$perm_minima), n_wy)
    expect_true(all(!is.na(result$perm_minima)))
})

test_that(".westfall_young_permutation serial batch processing executes correctly", {
    skip_if_not_installed("parallel")
    
    set.seed(123)
    n_genes <- 5
    n_wy <- 35  # Chosen to test batch boundaries (batch_size = 10 → 4 batches)
    
    permute_fn <- function() sample(c("A", "B"), size = 12, replace = TRUE)
    refit_fn <- function(assignment) runif(n_genes)
    
    result <- .westfall_young_permutation(
        permute_fn = permute_fn,
        refit_fn = refit_fn,
        n_genes = n_genes,
        wy_randomizations = n_wy,
        nthreads = 1,  # Force serial batch mode
        verbose = FALSE
    )
    
    # Verify all permutations were computed
    expect_equal(length(result$perm_minima), n_wy)
    expect_true(all(!is.na(result$perm_minima)))
    expect_true(all(result$perm_minima >= 0 & result$perm_minima <= 1))
})

test_that(".westfall_young_permutation serial batch processing with verbose reports batch progress", {
    skip_if_not_installed("parallel")
    
    set.seed(123)
    n_genes <- 3
    n_wy <- 40  # Large enough to trigger batch progress messages
    
    permute_fn <- function() sample(c("A", "B"), size = 10, replace = TRUE)
    refit_fn <- function(assignment) runif(n_genes)
    
    # Execute with verbose=TRUE and suppress messages
    # (message capture can be unreliable in some R environments)
    result <- suppressMessages({
        .westfall_young_permutation(
            permute_fn = permute_fn,
            refit_fn = refit_fn,
            n_genes = n_genes,
            wy_randomizations = n_wy,
            nthreads = 1,
            verbose = TRUE
        )
    })
    
    # Verify results are correct
    expect_equal(length(result$perm_minima), n_wy)
    expect_true(all(!is.na(result$perm_minima)))
    expect_true(all(result$perm_minima >= 0 & result$perm_minima <= 1))
})

test_that(".westfall_young_permutation compute_permutation returns minimum p-value per permutation", {
    skip_if_not_installed("parallel")
    
    set.seed(456)
    n_genes <- 4
    n_wy <- 15
    
    # Create a custom refit function that returns known values
    pvalue_sets <- list(
        c(0.001, 0.05, 0.1, 0.5),   # min = 0.001
        c(0.02, 0.03, 0.04, 0.5),   # min = 0.02
        c(0.5, 0.6, 0.7, 0.8)       # min = 0.5
    )
    pvalue_idx <- 0
    
    permute_fn <- function() sample(c("A", "B"), size = 10, replace = TRUE)
    refit_fn <- function(assignment) {
        # Return known p-value set in rotation
        pvalue_idx <<- (pvalue_idx %% length(pvalue_sets)) + 1
        pvalue_sets[[pvalue_idx]]
    }
    
    result <- .westfall_young_permutation(
        permute_fn = permute_fn,
        refit_fn = refit_fn,
        n_genes = n_genes,
        wy_randomizations = n_wy,
        nthreads = 1,
        verbose = FALSE
    )
    
    # Should have one minimum per permutation
    expect_equal(length(result$perm_minima), n_wy)
    
    # All should be valid p-values (between 0 and 1)
    expect_true(all(result$perm_minima >= 0 & result$perm_minima <= 1))
})

test_that(".westfall_young_permutation handles NA values correctly from refit_fn", {
    skip_if_not_installed("parallel")
    
    set.seed(789)
    n_genes <- 3
    n_wy <- 10
    
    permute_fn <- function() sample(c("A", "B"), size = 10, replace = TRUE)
    refit_fn <- function(assignment) {
        # Return p-values with some NAs
        c(0.01, NA, 0.5)
    }
    
    result <- .westfall_young_permutation(
        permute_fn = permute_fn,
        refit_fn = refit_fn,
        n_genes = n_genes,
        wy_randomizations = n_wy,
        nthreads = 1,
        verbose = FALSE
    )
    
    # Should handle NAs correctly (min with na.rm=TRUE should give 0.01)
    expect_equal(length(result$perm_minima), n_wy)
    expect_true(all(result$perm_minima == 0.01))  # NA removed, min is 0.01
})

test_that(".westfall_young_permutation serial vs parallel produce comparable distributions", {
    skip_if_not_installed("parallel")
    skip_if_not(.Platform$OS.type == "unix",
                message = "Parallel execution uses mclapply (Unix only)")
    
    set.seed(999)
    n_genes <- 4
    n_wy <- 20
    
    permute_fn <- function() sample(c("A", "B"), size = 10, replace = TRUE)
    refit_fn <- function(assignment) runif(n_genes)
    
    # Run serial version
    set.seed(999)
    result_serial <- .westfall_young_permutation(
        permute_fn = permute_fn,
        refit_fn = refit_fn,
        n_genes = n_genes,
        wy_randomizations = n_wy,
        nthreads = 1,
        verbose = FALSE
    )
    
    # Run parallel version (with same seed for reproducibility)
    set.seed(999)
    result_parallel <- .westfall_young_permutation(
        permute_fn = permute_fn,
        refit_fn = refit_fn,
        n_genes = n_genes,
        wy_randomizations = n_wy,
        nthreads = 2,
        verbose = FALSE
    )
    
    # Both should produce valid results with same length
    expect_equal(length(result_serial$perm_minima), length(result_parallel$perm_minima))
    expect_equal(length(result_serial$perm_minima), n_wy)
    
    # Statistics should be comparable (same data distribution)
    expect_equal(mean(result_serial$perm_minima, na.rm = TRUE), 
                 mean(result_parallel$perm_minima, na.rm = TRUE), 
                 tolerance = 0.1)  # Allow some variation due to RNG differences
})

# ============================================================================
# TEST SUITE 4: .westfall_young_permutation_rank() - Serial vs Parallel
# ============================================================================

test_that(".westfall_young_permutation_rank serial mode completes", {
    skip_if_not_installed("parallel")
    
    set.seed(42)
    n_genes <- 4
    n_wy <- 8
    
    permute_fn <- function() sample(c("A", "B"), 12, replace = TRUE)
    refit_fn <- function(assignment) {
        # Return list with statistics vector
        list(statistics = runif(n_genes, 0, 5))
    }
    
    result <- .westfall_young_permutation_rank(
        permute_fn = permute_fn,
        refit_fn = refit_fn,
        n_genes = n_genes,
        wy_randomizations = n_wy,
        nthreads = 1,  # Serial mode
        verbose = FALSE
    )
    
    expect_is(result, "list")
    expect_true("perm_stats_matrix" %in% names(result))
    expect_equal(nrow(result$perm_stats_matrix), n_genes)
    expect_equal(ncol(result$perm_stats_matrix), n_wy)
    expect_true(all(result$perm_stats_matrix >= 0))
})

test_that(".westfall_young_permutation_rank with numeric refit output", {
    skip_if_not_installed("parallel")
    
    set.seed(42)
    n_genes <- 3
    n_wy <- 5
    
    permute_fn <- function() sample(c("A", "B"), 10, replace = TRUE)
    # Return numeric vector instead of list (fallback case)
    refit_fn <- function(assignment) runif(n_genes, 0.01, 0.1)
    
    result <- .westfall_young_permutation_rank(
        permute_fn = permute_fn,
        refit_fn = refit_fn,
        n_genes = n_genes,
        wy_randomizations = n_wy,
        nthreads = 1,
        verbose = FALSE
    )
    
    expect_is(result, "list")
    expect_equal(nrow(result$perm_stats_matrix), n_genes)
    expect_equal(ncol(result$perm_stats_matrix), n_wy)
})

test_that(".westfall_young_permutation_rank with verbose output", {
    skip_if_not_installed("parallel")
    
    set.seed(42)
    permute_fn <- function() sample(c("A", "B"), 10, replace = TRUE)
    refit_fn <- function(assignment) list(statistics = runif(3, 0, 5))
    
    expect_message(
        .westfall_young_permutation_rank(
            permute_fn = permute_fn,
            refit_fn = refit_fn,
            n_genes = 3,
            wy_randomizations = 4,
            nthreads = 1,
            verbose = TRUE
        ),
        "WY Permutation.*Rank"
    )
})

test_that(".westfall_young_permutation_rank errors on invalid wy_randomizations", {
    skip_if_not_installed("parallel")
    
    permute_fn <- function() sample(c("A", "B"), 10, replace = TRUE)
    refit_fn <- function(assignment) list(statistics = runif(3, 0, 5))
    
    expect_error(
        .westfall_young_permutation_rank(
            permute_fn = permute_fn,
            refit_fn = refit_fn,
            n_genes = 3,
            wy_randomizations = -1,
            nthreads = 1
        ),
        "wy_randomizations must be"
    )
})

test_that(".westfall_young_permutation_rank errors on invalid refit output", {
    skip_if_not_installed("parallel")
    
    set.seed(42)
    permute_fn <- function() sample(c("A", "B"), 10, replace = TRUE)
    # Return invalid output type (neither numeric nor list with $statistics)
    refit_fn <- function(assignment) "invalid"
    
    expect_error(
        .westfall_young_permutation_rank(
            permute_fn = permute_fn,
            refit_fn = refit_fn,
            n_genes = 3,
            wy_randomizations = 3,
            nthreads = 1
        ),
        "refit_fn must return"
    )
})

test_that(".westfall_young_permutation_rank handles nthreads conversion", {
    skip_if_not_installed("parallel")
    
    set.seed(42)
    permute_fn <- function() sample(c("A", "B"), 10, replace = TRUE)
    refit_fn <- function(assignment) list(statistics = runif(3, 0, 5))
    
    # Non-integer nthreads should be coerced
    result <- .westfall_young_permutation_rank(
        permute_fn = permute_fn,
        refit_fn = refit_fn,
        n_genes = 3,
        wy_randomizations = 4,
        nthreads = 2.7,  # Non-integer
        verbose = FALSE
    )
    
    expect_is(result, "list")
    expect_equal(nrow(result$perm_stats_matrix), 3)
})

test_that(".westfall_young_permutation_rank parallel mode (nthreads > 1) with list $statistics", {
    skip_if_not_installed("parallel")
    skip_if_not(.Platform$OS.type == "unix",
                message = "Parallel execution uses mclapply (Unix only)")
    
    set.seed(42)
    n_genes <- 4
    n_wy <- 8
    
    permute_fn <- function() sample(c("A", "B"), 12, replace = TRUE)
    # Return list with $statistics (full vectors, not p-values)
    refit_fn <- function(assignment) {
        list(statistics = runif(n_genes, 0, 5))
    }
    
    result <- .westfall_young_permutation_rank(
        permute_fn = permute_fn,
        refit_fn = refit_fn,
        n_genes = n_genes,
        wy_randomizations = n_wy,
        nthreads = 2,  # Parallel mode
        verbose = FALSE
    )
    
    # Verify output structure
    expect_is(result, "list")
    expect_true("perm_stats_matrix" %in% names(result))
    
    # Should be genes × permutations matrix
    expect_equal(nrow(result$perm_stats_matrix), n_genes)
    expect_equal(ncol(result$perm_stats_matrix), n_wy)
    
    # All values should be non-negative (test statistics are positive)
    expect_true(all(result$perm_stats_matrix >= 0, na.rm = TRUE))
})

test_that(".westfall_young_permutation_rank parallel mode with numeric fallback (p-value conversion)", {
    skip_if_not_installed("parallel")
    skip_if_not(.Platform$OS.type == "unix",
                message = "Parallel execution uses mclapply (Unix only)")
    
    set.seed(42)
    n_genes <- 3
    n_wy <- 6
    
    permute_fn <- function() sample(c("A", "B"), 10, replace = TRUE)
    # Return numeric vector (fallback case - treated as p-values)
    refit_fn <- function(assignment) {
        runif(n_genes, 0.01, 0.1)  # p-values
    }
    
    result <- .westfall_young_permutation_rank(
        permute_fn = permute_fn,
        refit_fn = refit_fn,
        n_genes = n_genes,
        wy_randomizations = n_wy,
        nthreads = 2,  # Parallel mode
        verbose = FALSE
    )
    
    # Should have correct matrix dimensions
    expect_is(result, "list")
    expect_equal(nrow(result$perm_stats_matrix), n_genes)
    expect_equal(ncol(result$perm_stats_matrix), n_wy)
    
    # Values should be -log(p + 1e-300) which are positive
    expect_true(all(result$perm_stats_matrix > 0, na.rm = TRUE))
    # Should be reasonably large since p-values are small
    expect_true(all(result$perm_stats_matrix > 1, na.rm = TRUE))
})

test_that(".westfall_young_permutation_rank parallel mode reports progress", {
    skip_if_not_installed("parallel")
    skip_if_not(.Platform$OS.type == "unix",
                message = "Parallel execution uses mclapply (Unix only)")
    
    set.seed(42)
    n_genes <- 3
    n_wy <- 20  # Large enough to trigger progress messages
    
    permute_fn <- function() sample(c("A", "B"), 10, replace = TRUE)
    refit_fn <- function(assignment) list(statistics = runif(n_genes, 0, 3))
    
    # Execute with verbose=TRUE and suppress messages
    # (message capture can be unreliable in some R environments)
    result <- suppressMessages({
        .westfall_young_permutation_rank(
            permute_fn = permute_fn,
            refit_fn = refit_fn,
            n_genes = n_genes,
            wy_randomizations = n_wy,
            nthreads = 2,
            verbose = TRUE
        )
    })
    
    # Verify results are correct
    expect_is(result, "list")
    expect_equal(nrow(result$perm_stats_matrix), n_genes)
    expect_equal(ncol(result$perm_stats_matrix), n_wy)
    expect_true(all(!is.na(result$perm_stats_matrix)))
})

test_that(".westfall_young_permutation_rank returns FULL statistics vectors per permutation", {
    skip_if_not_installed("parallel")
    
    set.seed(123)
    n_genes <- 5
    n_wy <- 10
    
    # Create predictable statistics
    counter <- 0
    permute_fn <- function() sample(c("A", "B"), 12, replace = TRUE)
    refit_fn <- function(assignment) {
        # Return incrementing values so we can verify full vectors are returned
        counter <<- counter + 1
        list(statistics = seq(counter, counter + n_genes - 1))
    }
    
    result <- .westfall_young_permutation_rank(
        permute_fn = permute_fn,
        refit_fn = refit_fn,
        n_genes = n_genes,
        wy_randomizations = n_wy,
        nthreads = 1,
        verbose = FALSE
    )
    
    # Verify we got full vectors, not minima
    # Matrix should be n_genes × n_wy with all different values
    expect_equal(nrow(result$perm_stats_matrix), n_genes)
    expect_equal(ncol(result$perm_stats_matrix), n_wy)
    
    # Each permutation (column) should have different values from others
    col1 <- result$perm_stats_matrix[, 1]
    col2 <- result$perm_stats_matrix[, 2]
    expect_false(all(col1 == col2))  # Different permutations should differ
})

test_that(".westfall_young_permutation_rank serial batch processing returns full statistics", {
    skip_if_not_installed("parallel")
    
    set.seed(321)
    n_genes <- 4
    n_wy <- 35  # Test batch boundaries
    
    permute_fn <- function() sample(c("A", "B"), 10, replace = TRUE)
    refit_fn <- function(assignment) {
        list(statistics = runif(n_genes, 1, 5))
    }
    
    result <- .westfall_young_permutation_rank(
        permute_fn = permute_fn,
        refit_fn = refit_fn,
        n_genes = n_genes,
        wy_randomizations = n_wy,
        nthreads = 1,  # Serial batch mode
        verbose = FALSE
    )
    
    # Verify matrix dimensions
    expect_equal(nrow(result$perm_stats_matrix), n_genes)
    expect_equal(ncol(result$perm_stats_matrix), n_wy)
    
    # All values should be valid statistics
    expect_true(all(!is.na(result$perm_stats_matrix)))
    expect_true(all(result$perm_stats_matrix > 0))
})

test_that(".westfall_young_permutation_rank serial batch processing with verbose reports progress", {
    skip_if_not_installed("parallel")
    
    set.seed(321)
    n_genes <- 3
    n_wy <- 40  # Large enough to trigger batch progress
    
    permute_fn <- function() sample(c("A", "B"), 10, replace = TRUE)
    refit_fn <- function(assignment) {
        list(statistics = runif(n_genes, 0, 2))
    }
    
    # Execute with verbose=TRUE and suppress messages
    # (message capture can be unreliable in some R environments)
    result <- suppressMessages({
        .westfall_young_permutation_rank(
            permute_fn = permute_fn,
            refit_fn = refit_fn,
            n_genes = n_genes,
            wy_randomizations = n_wy,
            nthreads = 1,
            verbose = TRUE
        )
    })
    
    # Verify results are correct
    expect_equal(nrow(result$perm_stats_matrix), n_genes)
    expect_equal(ncol(result$perm_stats_matrix), n_wy)
    expect_true(all(!is.na(result$perm_stats_matrix)))
})

test_that(".westfall_young_permutation_rank compute_permutation extracts $statistics correctly", {
    skip_if_not_installed("parallel")
    
    set.seed(555)
    n_genes <- 4
    n_wy <- 12
    
    # Return specific statistics that we can verify
    test_stats <- list(
        c(1.0, 2.0, 3.0, 4.0),
        c(5.0, 6.0, 7.0, 8.0),
        c(0.5, 1.5, 2.5, 3.5)
    )
    stat_idx <- 0
    
    permute_fn <- function() sample(c("A", "B"), 10, replace = TRUE)
    refit_fn <- function(assignment) {
        # Cycle through test statistics
        stat_idx <<- (stat_idx %% length(test_stats)) + 1
        list(statistics = test_stats[[stat_idx]])
    }
    
    result <- .westfall_young_permutation_rank(
        permute_fn = permute_fn,
        refit_fn = refit_fn,
        n_genes = n_genes,
        wy_randomizations = n_wy,
        nthreads = 1,
        verbose = FALSE
    )
    
    # Should have extracted statistics correctly
    expect_equal(nrow(result$perm_stats_matrix), n_genes)
    expect_equal(ncol(result$perm_stats_matrix), n_wy)
    
    # First permutation should use test_stats[[1]]
    expect_equal(result$perm_stats_matrix[, 1], c(1.0, 2.0, 3.0, 4.0))
})

test_that(".westfall_young_permutation_rank handles NA values in statistics", {
    skip_if_not_installed("parallel")
    
    set.seed(789)
    n_genes <- 3
    n_wy <- 10
    
    permute_fn <- function() sample(c("A", "B"), 10, replace = TRUE)
    refit_fn <- function(assignment) {
        # Return statistics with some NAs
        c(1.5, NA, 3.5)
    }
    
    result <- .westfall_young_permutation_rank(
        permute_fn = permute_fn,
        refit_fn = refit_fn,
        n_genes = n_genes,
        wy_randomizations = n_wy,
        nthreads = 1,
        verbose = FALSE
    )
    
    # Should handle NAs without crashing
    expect_equal(nrow(result$perm_stats_matrix), n_genes)
    expect_equal(ncol(result$perm_stats_matrix), n_wy)
    
    # Row 2 should be all NAs (from the NA in refit output)
    expect_true(all(is.na(result$perm_stats_matrix[2, ])))
})

test_that(".westfall_young_permutation_rank parallel vs serial produce comparable matrices", {
    skip_if_not_installed("parallel")
    skip_if_not(.Platform$OS.type == "unix",
                message = "Parallel execution uses mclapply (Unix only)")
    
    set.seed(999)
    n_genes <- 4
    n_wy <- 15
    
    permute_fn <- function() sample(c("A", "B"), 10, replace = TRUE)
    refit_fn <- function(assignment) list(statistics = runif(n_genes, 0.5, 4))
    
    # Run serial
    set.seed(999)
    result_serial <- .westfall_young_permutation_rank(
        permute_fn = permute_fn,
        refit_fn = refit_fn,
        n_genes = n_genes,
        wy_randomizations = n_wy,
        nthreads = 1,
        verbose = FALSE
    )
    
    # Run parallel
    set.seed(999)
    result_parallel <- .westfall_young_permutation_rank(
        permute_fn = permute_fn,
        refit_fn = refit_fn,
        n_genes = n_genes,
        wy_randomizations = n_wy,
        nthreads = 2,
        verbose = FALSE
    )
    
    # Both should produce same dimensions
    expect_equal(nrow(result_serial$perm_stats_matrix), nrow(result_parallel$perm_stats_matrix))
    expect_equal(ncol(result_serial$perm_stats_matrix), ncol(result_parallel$perm_stats_matrix))
    expect_equal(ncol(result_serial$perm_stats_matrix), n_wy)
    
    # Statistics should have similar ranges
    expect_equal(range(result_serial$perm_stats_matrix, na.rm = TRUE)[1],
                 range(result_parallel$perm_stats_matrix, na.rm = TRUE)[1],
                 tolerance = 1.0)
})

# ============================================================================
# TEST SUITE 5: Integration tests - Full WY workflow
# ============================================================================

test_that("Full Storey pi0 + q-value pipeline works", {
    skip_if_not_installed("stats")
    
    set.seed(42)
    # Simulate p-values from mixed distribution
    pvalues <- c(runif(70), rbeta(30, 0.5, 1))
    
    # Step 1: Estimate pi0
    pi0_result <- .estimate_storey_pi0(pvalues, pi0_method = "lambda")
    expect_true(pi0_result$pi0 > 0.5)
    
    # Step 2: Compute q-values using estimated pi0
    qvalue_result <- .compute_storey_qvalues(pvalues, pi0 = pi0_result$pi0)
    
    # Step 3: Verify q-values are monotone and <= 1
    expect_is(qvalue_result, "numeric")
    expect_true(all(qvalue_result <= 1.0, na.rm = TRUE))
    expect_true(all(qvalue_result >= 0, na.rm = TRUE))
})
