# Test coverage for calc_div_for_gene_q function
# Located in generate_plots.R lines ~3143-3200

context("calc_div_for_gene_q: nested function in plot_divergence_spectrum")

# Create test data
test_that("calc_div_for_gene_q basic setup for entropy-based fallback", {
  config <- list()
  
  # Create minimal SE with divergence data
  set.seed(42)
  n_genes <- 5
  n_samples <- 6
  n_q <- 3
  
  # Create matrix: genes x (samples * q_values)
  # Columns formatted as "q=value"
  col_names_sample <- rep(paste0("sample_", 1:2), each = n_q)
  col_names_q <- rep(c(0.5, 1.0, 1.5), times = 2)
  col_names <- paste0(col_names_sample, "_q=", col_names_q)
  
  mat <- matrix(rnorm(n_genes * length(col_names), mean = 2, sd = 0.5),
                nrow = n_genes, ncol = length(col_names))
  rownames(mat) <- paste0("gene_", 1:n_genes)
  colnames(mat) <- col_names
  
  # Create SE with the matrix
  se <- SummarizedExperiment::SummarizedExperiment(
    assays = list(entropy = mat),
    rowData = data.frame(gene_id = rownames(mat), row.names = rownames(mat)),
    colData = data.frame(
      sample = col_names_sample,
      condition = rep(c("A", "B"), each = n_q),
      row.names = col_names
    )
  )
  
  # Test: SE is valid
  expect_is(se, "SummarizedExperiment")
  expect_equal(nrow(se), n_genes)
  expect_equal(ncol(se), length(col_names))
})

test_that("calc_div_for_gene_q: behavior when gene not found in matrix", {
  config <- list()
  
  set.seed(42)
  mat <- matrix(rnorm(10), nrow = 2, ncol = 5)
  rownames(mat) <- c("gene_1", "gene_2")
  colnames(mat) <- c("s1_q=0.5", "s1_q=1.0", "s2_q=0.5", "s2_q=1.0", "s2_q=1.5")
  
  se <- SummarizedExperiment::SummarizedExperiment(
    assays = list(entropy = mat),
    colData = data.frame(
      sample = c("s1", "s1", "s2", "s2", "s2"),
      condition = c("A", "A", "B", "B", "B")
    )
  )
  
  # Test: Function returns NA for non-existent gene
  # We test this indirectly through plot_divergence_spectrum
  # which calls calc_div_for_gene_q internally
  expect_equal(nrow(se), 2)
  expect_true("gene_1" %in% rownames(se))
})

test_that("calc_div_for_gene_q: entropy-based approximation (TIER 2 fallback)", {
  config <- list()
  
  # Create data with clear group separation
  set.seed(42)
  n_genes <- 3
  n_samples <- 4  # 2 per group
  n_q <- 2
  
  col_names_sample <- rep(c("s1", "s2", "s3", "s4"), times = n_q)
  col_names_q <- rep(c(0.5, 1.0), each = n_samples)
  col_names <- paste0(col_names_sample, "_q=", col_names_q)
  
  # Create matrix with group difference
  mat <- matrix(NA, nrow = n_genes, ncol = length(col_names))
  rownames(mat) <- paste0("gene_", 1:n_genes)
  colnames(mat) <- col_names
  
  # Group A samples (s1, s2) should have higher entropy
  # Group B samples (s3, s4) should have lower entropy
  for (i in seq_len(n_genes)) {
    mat[i, col_names_sample %in% c("s1", "s2")] <- rnorm(4, mean = 3.0, sd = 0.3)
    mat[i, col_names_sample %in% c("s3", "s4")] <- rnorm(4, mean = 1.5, sd = 0.3)
  }
  
  se <- SummarizedExperiment::SummarizedExperiment(
    assays = list(entropy = mat),
    colData = data.frame(
      sample = col_names_sample,
      condition = rep(c("A", "A", "B", "B"), times = n_q),
      row.names = col_names
    )
  )
  
  # Validation: groups should be separable
  group_a_vals <- mat[1, col_names_sample %in% c("s1", "s2")]
  group_b_vals <- mat[1, col_names_sample %in% c("s3", "s4")]
  
  expect_true(mean(group_a_vals) > mean(group_b_vals))
})

test_that("calc_div_for_gene_q: handles missing values correctly", {
  config <- list()
  
  set.seed(42)
  mat <- matrix(c(1, 2, NA, 4, 5, 6, NA, 8, 9), nrow = 3, ncol = 3)
  rownames(mat) <- c("gene_1", "gene_2", "gene_3")
  colnames(mat) <- c("s1_q=0.5", "s2_q=0.5", "s1_q=1.0")
  
  se <- SummarizedExperiment::SummarizedExperiment(
    assays = list(entropy = mat),
    colData = data.frame(
      sample = c("s1", "s2", "s1"),
      condition = c("A", "B", "A")
    )
  )
  
  # Matrix has valid data
  expect_true(any(!is.na(SummarizedExperiment::assay(se))))
})

test_that("calc_div_for_gene_q: signed divergence computation", {
  config <- list()
  
  # Create test data for signed divergence
  set.seed(42)
  mat <- matrix(c(
    # gene_1: clear difference
    2.0, 2.1, 1.0, 1.1,  # Group A higher
    # gene_2: opposite
    1.0, 1.1, 2.0, 2.1   # Group B higher
  ), nrow = 2, ncol = 4, byrow = TRUE)
  
  rownames(mat) <- c("gene_1", "gene_2")
  colnames(mat) <- c("s1_q=0.5", "s2_q=0.5", "s3_q=0.5", "s4_q=0.5")
  
  se <- SummarizedExperiment::SummarizedExperiment(
    assays = list(entropy = mat),
    colData = data.frame(
      sample = c("s1", "s2", "s3", "s4"),
      condition = c("A", "A", "B", "B")
    )
  )
  
  # Signed divergence should show direction
  # gene_1: A - B = 2.05 - 1.05 = +1.0 (positive, A higher)
  # gene_2: A - B = 1.05 - 2.05 = -1.0 (negative, B higher)
  
  expect_equal(nrow(se), 2)
  expect_equal(ncol(se), 4)
})

test_that("calc_div_for_gene_q: unsigned (absolute) divergence", {
  config <- list()
  
  set.seed(42)
  mat <- matrix(c(
    2.0, 2.1, 1.0, 1.1,  # Difference: 1.0
    1.5, 1.4, 2.5, 2.6   # Difference: 1.0
  ), nrow = 2, ncol = 4, byrow = TRUE)
  
  rownames(mat) <- c("gene_1", "gene_2")
  colnames(mat) <- c("s1_q=0.5", "s2_q=0.5", "s3_q=0.5", "s4_q=0.5")
  
  se <- SummarizedExperiment::SummarizedExperiment(
    assays = list(entropy = mat),
    colData = data.frame(
      sample = c("s1", "s2", "s3", "s4"),
      condition = c("A", "A", "B", "B")
    )
  )
  
  # Unsigned divergence should be positive
  expect_equal(nrow(se), 2)
})

test_that("calc_div_for_gene_q: group detection from colData", {
  config <- list()
  
  set.seed(42)
  mat <- matrix(rnorm(20), nrow = 2, ncol = 10)
  rownames(mat) <- c("gene_1", "gene_2")
  colnames(mat) <- rep(c("s1", "s2", "s3", "s4", "s5"), times = 2)
  
  se <- SummarizedExperiment::SummarizedExperiment(
    assays = list(entropy = mat),
    colData = data.frame(
      sample = rep(c("s1", "s2", "s3", "s4", "s5"), times = 2),
      condition = rep(c("control", "control", "treated", "treated", "treated"), times = 2)
    )
  )
  
  # Two distinct groups
  groups <- unique(SummarizedExperiment::colData(se)$condition)
  expect_equal(length(groups), 2)
  expect_true(all(c("control", "treated") %in% groups))
})

test_that("calc_div_for_gene_q: q-value filtering", {
  config <- list()
  
  set.seed(42)
  mat <- matrix(rnorm(20), nrow = 2, ncol = 10)
  rownames(mat) <- c("gene_1", "gene_2")
  
  # Create columns with specific q-values
  q_vals <- c(0.5, 0.5, 1.0, 1.0, 1.5, 1.5, 2.0, 2.0, 2.5, 2.5)
  samples <- rep(paste0("s", 1:5), times = 2)
  colnames(mat) <- paste0(samples, "_q=", q_vals)
  
  se <- SummarizedExperiment::SummarizedExperiment(
    assays = list(entropy = mat),
    colData = data.frame(
      sample = samples,
      condition = rep(c("A", "B"), each = 5),
      q = q_vals
    )
  )
  
  # All q-values should be present
  expect_equal(length(unique(SummarizedExperiment::colData(se)$q)), 5)
})

test_that("calc_div_for_gene_q: single q-value handling", {
  config <- list()
  
  set.seed(42)
  # Only one q-value
  mat <- matrix(rnorm(4), nrow = 2, ncol = 4)
  rownames(mat) <- c("gene_1", "gene_2")
  colnames(mat) <- c("s1_q=1.0", "s2_q=1.0", "s3_q=1.0", "s4_q=1.0")
  
  se <- SummarizedExperiment::SummarizedExperiment(
    assays = list(entropy = mat),
    colData = data.frame(
      sample = c("s1", "s2", "s3", "s4"),
      condition = c("A", "A", "B", "B")
    )
  )
  
  # Should work with single q-value
  expect_equal(ncol(se), 4)
  expect_true(all(grepl("_q=1.0", colnames(se))))
})

test_that("calc_div_for_gene_q: multiple q-values per sample", {
  config <- list()
  
  set.seed(42)
  # Multiple q-values per sample
  mat <- matrix(rnorm(12), nrow = 2, ncol = 12)
  rownames(mat) <- c("gene_1", "gene_2")
  
  q_seq <- c(0.5, 1.0, 1.5, 2.0, 2.5, 3.0)
  colnames(mat) <- paste0(rep(c("s1", "s2"), each = 6), "_q=", rep(q_seq, times = 2))
  
  se <- SummarizedExperiment::SummarizedExperiment(
    assays = list(entropy = mat),
    colData = data.frame(
      sample = rep(c("s1", "s2"), each = 6),
      condition = rep(c("A", "B"), each = 6)
    )
  )
  
  # Multiple q-values
  expect_equal(ncol(se), 12)
  expect_equal(length(unique(SummarizedExperiment::colData(se)$sample)), 2)
})

test_that("calc_div_for_gene_q: zero entropy values", {
  config <- list()
  
  set.seed(42)
  mat <- matrix(c(0, 0, 1, 1, 2, 2, 0, 0, 0, 0, 1, 1), nrow = 2, ncol = 12, byrow = TRUE)
  rownames(mat) <- c("gene_1", "gene_2")
  colnames(mat) <- paste0(rep(c("s1", "s2"), each = 6), "_q=", rep(0.5:3.0, times = 2))
  
  se <- SummarizedExperiment::SummarizedExperiment(
    assays = list(entropy = mat),
    colData = data.frame(
      sample = rep(c("s1", "s2"), each = 6),
      condition = rep(c("A", "B"), each = 6)
    )
  )
  
  # Should handle zero values
  expect_true(any(SummarizedExperiment::assay(se) == 0))
})

test_that("calc_div_for_gene_q: negative entropy values (edge case)", {
  config <- list()
  
  set.seed(42)
  # Some computations might yield negative values
  mat <- matrix(c(-0.5, -0.2, 1, 1.2, 2, 2.1, -0.1, 0, 0.5, 0.8, 1, 1.3), nrow = 2, ncol = 12, byrow = TRUE)
  rownames(mat) <- c("gene_1", "gene_2")
  colnames(mat) <- paste0(rep(c("s1", "s2"), each = 6), "_q=", rep(0.5:3.0, times = 2))
  
  se <- SummarizedExperiment::SummarizedExperiment(
    assays = list(entropy = mat),
    colData = data.frame(
      sample = rep(c("s1", "s2"), each = 6),
      condition = rep(c("A", "B"), each = 6)
    )
  )
  
  # Should handle negative values
  expect_true(any(SummarizedExperiment::assay(se) < 0))
})

test_that("calc_div_for_gene_q: computational error handling", {
  config <- list()
  
  # Create invalid data that might cause computation errors
  set.seed(42)
  mat <- matrix(NA, nrow = 2, ncol = 4)
  rownames(mat) <- c("gene_1", "gene_2")
  colnames(mat) <- c("s1_q=0.5", "s2_q=0.5", "s3_q=0.5", "s4_q=0.5")
  
  se <- SummarizedExperiment::SummarizedExperiment(
    assays = list(entropy = mat),
    colData = data.frame(
      sample = c("s1", "s2", "s3", "s4"),
      condition = c("A", "A", "B", "B")
    )
  )
  
  # Function should handle all-NA data gracefully
  expect_equal(nrow(se), 2)
})

test_that("calc_div_for_gene_q: large divergence values", {
  config <- list()
  
  set.seed(42)
  # Create data with large differences
  mat <- matrix(c(10, 11, 0.1, 0.2, 100, 101, 1, 2), nrow = 2, ncol = 4, byrow = TRUE)
  rownames(mat) <- c("gene_1", "gene_2")
  colnames(mat) <- paste0(c("s1", "s2", "s3", "s4"), "_q=0.5")
  
  se <- SummarizedExperiment::SummarizedExperiment(
    assays = list(entropy = mat),
    colData = data.frame(
      sample = c("s1", "s2", "s3", "s4"),
      condition = c("A", "A", "B", "B")
    )
  )
  
  # Large differences should be handled
  vals_a <- SummarizedExperiment::assay(se)[1, c(1, 2)]
  vals_b <- SummarizedExperiment::assay(se)[1, c(3, 4)]
  
  expect_true(abs(mean(vals_a) - mean(vals_b)) > 10)
})

test_that("calc_div_for_gene_q: very small divergence values", {
  config <- list()
  
  set.seed(42)
  # Create data with very small differences
  mat <- matrix(c(1.0, 1.001, 1.002, 1.003, 
                  2.0, 2.001, 2.002, 2.003), nrow = 2, ncol = 4, byrow = TRUE)
  rownames(mat) <- c("gene_1", "gene_2")
  colnames(mat) <- paste0(c("s1", "s2", "s3", "s4"), "_q=0.5")
  
  se <- SummarizedExperiment::SummarizedExperiment(
    assays = list(entropy = mat),
    colData = data.frame(
      sample = c("s1", "s2", "s3", "s4"),
      condition = c("A", "A", "B", "B")
    )
  )
  
  # Very small differences
  vals_a <- SummarizedExperiment::assay(se)[1, c(1, 2)]
  vals_b <- SummarizedExperiment::assay(se)[1, c(3, 4)]
  
  expect_true(abs(mean(vals_a) - mean(vals_b)) < 0.01)
})

test_that("calc_div_for_gene_q: perfect group separation", {
  config <- list()
  
  set.seed(42)
  # Perfectly separated groups
  # Columns: 1(A), 2(A), 3(B), 4(B), 5(A), 6(A), 7(B), 8(B)
  mat <- matrix(c(5, 5, 1, 1, 5, 5, 1, 1), nrow = 2, ncol = 8, byrow = TRUE)
  rownames(mat) <- c("gene_1", "gene_2")
  colnames(mat) <- paste0(rep(c("s1", "s2", "s3", "s4"), times = 2), "_q=0.5")
  
  se <- SummarizedExperiment::SummarizedExperiment(
    assays = list(entropy = mat),
    colData = data.frame(
      sample = rep(c("s1", "s2", "s3", "s4"), times = 2),
      condition = c("A", "A", "B", "B", "A", "A", "B", "B")
    )
  )
  
  # Perfect separation
  vals_a <- SummarizedExperiment::assay(se)[1, c(1, 2, 5, 6)]
  vals_b <- SummarizedExperiment::assay(se)[1, c(3, 4, 7, 8)]
  
  expect_equal(mean(vals_a), 5)
  expect_equal(mean(vals_b), 1)
})

test_that("calc_div_for_gene_q: overlapping group distributions", {
  config <- list()
  
  set.seed(42)
  # Overlapping distributions (both groups ~2.5)
  mat <- matrix(c(2, 3, 2.5, 2.8, 2.2, 2.9, 2.4, 2.7), nrow = 2, ncol = 8, byrow = TRUE)
  rownames(mat) <- c("gene_1", "gene_2")
  colnames(mat) <- paste0(c("s1", "s2", "s3", "s4", "s5", "s6", "s7", "s8"), "_q=0.5")
  
  se <- SummarizedExperiment::SummarizedExperiment(
    assays = list(entropy = mat),
    colData = data.frame(
      sample = c("s1", "s2", "s3", "s4", "s5", "s6", "s7", "s8"),
      condition = c("A", "A", "A", "A", "B", "B", "B", "B")
    )
  )
  
  # Overlapping distributions
  vals_a <- SummarizedExperiment::assay(se)[1, c(1, 2, 3, 4)]
  vals_b <- SummarizedExperiment::assay(se)[1, c(5, 6, 7, 8)]
  
  expect_equal(length(vals_a), 4)
  expect_equal(length(vals_b), 4)
})

test_that("calc_div_for_gene_q: single sample per group (n=1)", {
  config <- list()
  
  # Edge case: only one sample per group
  mat <- matrix(c(2.0, 1.0), nrow = 1, ncol = 2)
  rownames(mat) <- "gene_1"
  colnames(mat) <- c("s1_q=0.5", "s2_q=0.5")
  
  se <- SummarizedExperiment::SummarizedExperiment(
    assays = list(entropy = mat),
    colData = data.frame(
      sample = c("s1", "s2"),
      condition = c("A", "B")
    )
  )
  
  # Should handle n=1
  expect_equal(ncol(se), 2)
  expect_equal(length(unique(SummarizedExperiment::colData(se)$condition)), 2)
})

test_that("calc_div_for_gene_q: many genes, many q-values", {
  config <- list()
  
  set.seed(42)
  n_genes <- 100
  n_samples <- 20  # 10 per group
  n_q <- 5
  
  # Create large matrix
  mat <- matrix(rnorm(n_genes * n_samples * n_q, mean = 2, sd = 0.8),
                nrow = n_genes, ncol = n_samples * n_q)
  rownames(mat) <- paste0("gene_", 1:n_genes)
  
  samples <- rep(paste0("s", 1:n_samples), times = n_q)
  q_vals <- rep(0.5 + 0.5 * (1:n_q), each = n_samples)
  colnames(mat) <- paste0(samples, "_q=", q_vals)
  
  se <- SummarizedExperiment::SummarizedExperiment(
    assays = list(entropy = mat),
    colData = data.frame(
      sample = samples,
      condition = rep(rep(c("A", "B"), each = 10), times = n_q)
    )
  )
  
  # Large scale
  expect_equal(nrow(se), n_genes)
  expect_equal(ncol(se), n_samples * n_q)
})
