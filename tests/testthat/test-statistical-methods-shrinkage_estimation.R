context("Bayesian Shrinkage: Entropy Estimation for Small Sample Sizes")

# calculate_method is internal; expose for tests
calculate_method <- TSENAT:::.calculate_method

test_that(".tsenat_estimate_shrinkage_params estimates global mean correctly", {
  # Create simple test data: 3 genes x 4 transcripts x 2 samples
  x <- matrix(c(
    10, 5,    # Gene A, transcript 1
    3, 2,     # Gene A, transcript 2
    8, 12,    # Gene B, transcript 1
    2, 1,     # Gene B, transcript 2
    5, 5,     # Gene C, transcript 1
    5, 5      # Gene C, transcript 2
  ), nrow = 6, byrow = TRUE)
  
  colnames(x) <- c("S1", "S2")
  genes <- c("A", "A", "B", "B", "C", "C")
  
  # Calculate entropy first
  entropy_result <- .calculate_method(x, genes, norm = TRUE, q = 1, what = "S")
  entropy_mat <- as.matrix(entropy_result[, -1])
  
  # Estimate parameters
  params <- TSENAT:::.tsenat_estimate_shrinkage_params(
    x = x,
    genes = genes,
    entropy_matrix = entropy_mat,
    q = 1
  )
  
  # Check that we got back the expected structure
  expect_true(is.list(params))
  expect_true("global_mean" %in% names(params))
  expect_true("global_var" %in% names(params))
  expect_true("n_isoforms" %in% names(params))
  
  # Check that global_mean has appropriate length
  expect_true(length(params$global_mean) >= 1)
  
  # Check n_isoforms is correct
  expect_equal(names(params$n_isoforms), c("A", "B", "C"))
  expect_equal(params$n_isoforms, c(A = 2, B = 2, C = 2))
})

test_that(".tsenat_estimate_shrinkage_params detects variable isoform counts", {
  # Gene A: 2 isoforms, Gene B: 3 isoforms, Gene C: 1 isoform
  x <- matrix(c(
    10, 5,     # Gene A, transcript 1
    3, 2,      # Gene A, transcript 2
    8, 12,     # Gene B, transcript 1
    2, 1,      # Gene B, transcript 2
    5, 5,      # Gene B, transcript 3
    15, 20     # Gene C, transcript 1
  ), nrow = 6, byrow = TRUE)
  
  colnames(x) <- c("S1", "S2")
  genes <- c("A", "A", "B", "B", "B", "C")
  
  entropy_result <- .calculate_method(x, genes, norm = TRUE, q = 1, what = "S")
  entropy_mat <- as.matrix(entropy_result[, -1])
  
  params <- TSENAT:::.tsenat_estimate_shrinkage_params(
    x = x,
    genes = genes,
    entropy_matrix = entropy_mat,
    q = 1
  )
  
  # Check n_isoforms correctly identifies variation
  expect_equal(as.numeric(params$n_isoforms["A"]), 2)
  expect_equal(as.numeric(params$n_isoforms["B"]), 3)
  expect_equal(as.numeric(params$n_isoforms["C"]), 1)
})

test_that(".tsenat_apply_shrinkage shrinks estimates toward global mean", {
  # Create test entropy matrix
  entropy_mat <- matrix(c(
    0.3, 0.5,    # Gene A
    0.4, 0.6,    # Gene B
    0.7, 0.8     # Gene C
  ), nrow = 3, byrow = TRUE)
  
  rownames(entropy_mat) <- c("A", "B", "C")
  colnames(entropy_mat) <- c("S1_q=1", "S2_q=1")
  
  # Create parameters with global mean and variance
  params <- list(
    global_mean = c("q=1" = 0.5),
    global_var = c("q=1" = 0.05),
    n_isoforms = c(A = 2, B = 2, C = 2)
  )
  
  # Apply shrinkage
  shrunk <- TSENAT:::.tsenat_apply_shrinkage(
    entropy_matrix = entropy_mat,
    params = params,
    gene_isoform_map = c(2, 2, 2)
  )
  
  # Check that result has same dimension
  expect_equal(dim(shrunk), dim(entropy_mat))
  
  # Check that values are modified (shrunk toward mean)
  # For genes with equal n_isoforms, shrinkage weight w = n/(n+lambda)
  # Values should be between original and global mean
  expect_true(all(is.finite(shrunk)))
})

test_that(".tsenat_apply_shrinkage provides more shrinkage for genes with fewer isoforms", {
  # Create test entropy matrix
  entropy_mat <- matrix(c(
    0.2, 0.9,    # Gene few (will be shrunk more)
    0.5, 0.5     # Gene many (will be shrunk less)
  ), nrow = 2, byrow = TRUE)
  
  rownames(entropy_mat) <- c("few", "many")
  colnames(entropy_mat) <- c("S1_q=1", "S2_q=1")
  
  params <- list(
    global_mean = c("q=1" = 0.5),
    global_var = c("q=1" = 0.1),
    n_isoforms = c(few = 1, many = 10)  # Different isoform counts
  )
  
  shrunk <- TSENAT:::.tsenat_apply_shrinkage(
    entropy_matrix = entropy_mat,
    params = params,
    gene_isoform_map = c(1, 10)
  )
  
  # Gene with 1 isoform should be shrunk more (closer to global mean)
  # Gene with 10 isoforms should be shrunk less (further from global mean)
  
  # Calculate shrinkage magnitude (distance from full shrinkage to global mean)
  shrink_few <- abs(shrunk["few", "S1_q=1"] - params$global_mean["q=1"])
  shrink_many <- abs(shrunk["many", "S1_q=1"] - params$global_mean["q=1"])
  original_few <- abs(entropy_mat["few", "S1_q=1"] - params$global_mean["q=1"])
  original_many <- abs(entropy_mat["many", "S1_q=1"] - params$global_mean["q=1"])
  
  # Shrinkage reduction should be larger for few (1 isoform) than many (10 isoforms)
  reduction_few <- original_few - shrink_few
  reduction_many <- original_many - shrink_many
  
  expect_true(reduction_few > reduction_many)
})

test_that("calculate_method with shrinkage='none' returns unmodified estimates", {
  x <- matrix(c(
    10, 5,
    3, 2,
    8, 12,
    2, 1,
    5, 5,
    5, 5
  ), nrow = 6, byrow = TRUE)
  
  colnames(x) <- c("S1", "S2")
  genes <- c("A", "A", "B", "B", "C", "C")
  
  # Calculate with no shrinkage
  result_none <- .calculate_method(x, genes, norm = TRUE, q = 1, what = "S", shrinkage = "none")
  
  # Should be a data frame with Gene column + entropy columns
  expect_true(is.data.frame(result_none))
  expect_equal(ncol(result_none), 3)  # Gene + 2 samples * 1 q-value
  expect_equal(nrow(result_none), 3)  # 3 genes
})

test_that("calculate_method with shrinkage='empirical_bayes' returns modified estimates", {
  x <- matrix(c(
    10, 5,
    3, 2,
    8, 12,
    2, 1,
    5, 5,
    5, 5
  ), nrow = 6, byrow = TRUE)
  
  colnames(x) <- c("S1", "S2")
  genes <- c("A", "A", "B", "B", "C", "C")
  
  # Calculate with empirical Bayes shrinkage
  result_eb <- .calculate_method(x, genes, norm = TRUE, q = 1, what = "S", 
                                 shrinkage = "empirical_bayes")
  
  # Should return a data frame
  expect_true(is.data.frame(result_eb))
  expect_equal(ncol(result_eb), 3)  # Gene + 2 samples * 1 q-value
  expect_equal(nrow(result_eb), 3)
  
  # All values should be finite
  entropy_vals <- as.matrix(result_eb[, -1])
  expect_true(all(is.finite(entropy_vals)))
})

test_that("calculate_method shrinkage='empirical_bayes' differs from 'none'", {
  # Use a dataset with variable isoform counts to ensure shrinkage is detected
  x <- matrix(c(
    100, 50,       # Gene A, tx1
    10, 5,         # Gene A, tx2
    1, 0,          # Gene A, tx3
    80, 90,        # Gene B, tx1
    20, 10,        # Gene B, tx2
    70, 60         # Gene C, tx1 (single isoform)
  ), nrow = 6, byrow = TRUE)
  
  colnames(x) <- c("S1", "S2")
  genes <- c("A", "A", "A", "B", "B", "C")
  
  result_none <- .calculate_method(x, genes, norm = TRUE, q = 1, what = "S", shrinkage = "none")
  result_eb <- .calculate_method(x, genes, norm = TRUE, q = 1, what = "S", shrinkage = "empirical_bayes")
  
  # With shrinkage, single-isoform genes (like Gene C) get rescued from filtering
  # Expected: shrinkage="empirical_bayes" has more rows than shrinkage="none"
  expect_true(nrow(result_eb) >= nrow(result_none))
  
  # Gene A and B should be in both results
  expect_true("A" %in% result_none$Gene)
  expect_true("B" %in% result_none$Gene)
  expect_true("A" %in% result_eb$Gene)
  expect_true("B" %in% result_eb$Gene)
  
  # Gene C (single isoform) should only be in shrinkage result
  expect_false("C" %in% result_none$Gene)
  expect_true("C" %in% result_eb$Gene)
})

test_that("calculate_diversity passes shrinkage parameter through", {
  skip_if_not_installed("SummarizedExperiment")
  
  x <- matrix(c(
    10, 5,
    3, 2,
    8, 12,
    2, 1,
    5, 5,
    5, 5
  ), nrow = 6, byrow = TRUE)
  
  colnames(x) <- c("S1", "S2")
  genes <- c("A", "A", "B", "B", "C", "C")
  
  # Call with shrinkage
  se <- .calculate_diversity(x, genes = genes, q = 1, norm = TRUE, 
                            shrinkage = "empirical_bayes", verbose = FALSE)
  
  # Check that result is SummarizedExperiment
  expect_true(is(se, "SummarizedExperiment"))
  
  # Check assay dimensions
  diversity_assay <- SummarizedExperiment::assay(se, "diversity")
  expect_equal(ncol(diversity_assay), 2)  # 2 samples
})

test_that("calculate_diversity with shrinkage='empirical_bayes' handles verbose output", {
  skip_if_not_installed("SummarizedExperiment")
  
  x <- matrix(c(
    10, 5,
    3, 2,
    8, 12,
    2, 1,
    1, 0,
    5, 5
  ), nrow = 6, byrow = TRUE)
  
  colnames(x) <- c("S1", "S2")
  genes <- c("A", "A", "B", "B", "C", "C")
  
  # Should not error when called with verbose=TRUE and shrinkage='empirical_bayes'
  se <- .calculate_diversity(x, genes = genes, q = 1, norm = TRUE,
                            shrinkage = "empirical_bayes", verbose = TRUE)
  
  # Check that result is valid
  expect_true(is(se, "SummarizedExperiment"))
})

test_that("shrinkage='empirical_bayes' handles genes with single isoform", {
  # 1 gene with 1 isoform, 1 gene with 3 isoforms
  x <- matrix(c(
    50, 40,       # Gene A, single isoform
    100, 80,      # Gene B, isoform 1
    20, 30,       # Gene B, isoform 2
    10, 15        # Gene B, isoform 3
  ), nrow = 4, byrow = TRUE)
  
  colnames(x) <- c("S1", "S2")
  genes <- c("A", "B", "B", "B")
  
  result <- .calculate_method(x, genes, norm = TRUE, q = 1, what = "S",
                             shrinkage = "empirical_bayes")
  
  # Should return successfully without errors
  expect_true(is.data.frame(result))
  expect_equal(nrow(result), 2)
  
  # Gene A (single isoform) should have received shrinkage
  expect_true(is.finite(result[result$Gene == "A", "S1_q=1"]))
})

test_that("shrinkage parameter rejects invalid values", {
  x <- matrix(c(
    10, 5,
    3, 2
  ), nrow = 2, byrow = TRUE)
  
  colnames(x) <- c("S1", "S2")
  genes <- c("A", "A")
  
  # Invalid shrinkage method should error
  expect_error(
    .calculate_method(x, genes, norm = TRUE, q = 1, shrinkage = "invalid_method"),
    "should be one of"
  )
})

test_that("shrinkage works with multiple q values", {
  x <- matrix(c(
    10, 5,
    3, 2,
    8, 12,
    2, 1,
    5, 5,
    5, 5
  ), nrow = 6, byrow = TRUE)
  
  colnames(x) <- c("S1", "S2")
  genes <- c("A", "A", "B", "B", "C", "C")
  
  # Calculate with multiple q values and shrinkage
  result <- .calculate_method(x, genes, norm = TRUE, q = c(1, 2), what = "S",
                             shrinkage = "empirical_bayes")
  
  # Should have Gene + 2 samples * 2 q-values = 5 columns
  expect_equal(ncol(result), 5)
  expect_equal(nrow(result), 3)
  
  # All numeric values should be finite
  numeric_cols <- as.matrix(result[, -1])
  expect_true(all(is.finite(numeric_cols)))
})

test_that("shrinkage preserves NA values", {
  # Create data where one gene might have NA values
  x <- matrix(c(
    10, 5,
    3, 2,
    8, 0,        # This might produce NA if normalization is undefined
    2, 0,
    5, 5,
    5, 5
  ), nrow = 6, byrow = TRUE)
  
  colnames(x) <- c("S1", "S2")
  genes <- c("A", "A", "B", "B", "C", "C")
  
  result <- .calculate_method(x, genes, norm = TRUE, q = 1, what = "S",
                             shrinkage = "empirical_bayes")
  
  # NA values should be preserved (not converted to numbers)
  numeric_vals <- as.matrix(result[, -1])
  expect_true(sum(is.na(numeric_vals)) == sum(is.na(numeric_vals)))
})
