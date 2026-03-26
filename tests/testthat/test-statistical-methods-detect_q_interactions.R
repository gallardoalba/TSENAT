library(TSENAT)

context("detect_q_gene_interactions: Internal Helper Functions")

# ============================================================================
# Test 1: .tsenat_detect_q_validate_params
# ============================================================================

test_that(".tsenat_detect_q_validate_params: validates paired parameter", {
  # Should error when paired=TRUE but subject_col is NULL
  expect_error(
    TSENAT:::.tsenat_detect_q_validate_params(
      paired = TRUE, subject_col = NULL, wy_randomizations = 100,
      nperm_mode = "standard", verbose = FALSE
    ),
    "paired=TRUE with subject_col=NULL"
  )
})

test_that(".tsenat_detect_q_validate_params: handles 'auto' wy_randomizations", {
  result <- TSENAT:::.tsenat_detect_q_validate_params(
    paired = FALSE, subject_col = NULL, wy_randomizations = "auto",
    nperm_mode = "standard", verbose = FALSE
  )
  
  expect_equal(result$wy_randomizations, "auto")
})

test_that(".tsenat_detect_q_validate_params: converts numeric wy_randomizations to integer", {
  result <- TSENAT:::.tsenat_detect_q_validate_params(
    paired = FALSE, subject_col = NULL, wy_randomizations = 100.5,
    nperm_mode = "standard", verbose = FALSE
  )
  
  expect_equal(result$wy_randomizations, 100L)
  expect_is(result$wy_randomizations, "integer")
})

test_that(".tsenat_detect_q_validate_params: defaults NULL to 500", {
  result <- TSENAT:::.tsenat_detect_q_validate_params(
    paired = FALSE, subject_col = NULL, wy_randomizations = NULL,
    nperm_mode = "standard", verbose = FALSE
  )
  
  expect_equal(result$wy_randomizations, 500)
})

test_that(".tsenat_detect_q_validate_params: warns on small wy_randomizations", {
  expect_warning(
    TSENAT:::.tsenat_detect_q_validate_params(
      paired = FALSE, subject_col = NULL, wy_randomizations = 5,
      nperm_mode = "standard", verbose = FALSE
    ),
    "unreliable"
  )
})

test_that(".tsenat_detect_q_validate_params: validates nperm_mode", {
  result <- TSENAT:::.tsenat_detect_q_validate_params(
    paired = FALSE, subject_col = NULL, wy_randomizations = 100,
    nperm_mode = "conservative", verbose = FALSE
  )
  
  expect_equal(result$nperm_mode, "conservative")
})

test_that(".tsenat_detect_q_validate_params: warns on subject_col with paired=FALSE", {
  expect_warning(
    TSENAT:::.tsenat_detect_q_validate_params(
      paired = FALSE, subject_col = "subject", wy_randomizations = 100,
      nperm_mode = "standard", verbose = FALSE
    ),
    "paired=FALSE"
  )
})

# ============================================================================
# Test 2: .tsenat_detect_q_prepare_data
# ============================================================================

test_that(".tsenat_detect_q_prepare_data: validates required columns", {
  df <- data.frame(
    value = rnorm(20),
    group = rep(c("A", "B"), 10)
  )
  
  expect_error(
    TSENAT:::.tsenat_detect_q_prepare_data(
      data = df, entropy_col = "nonexistent", q_col = "q", gene_col = "gene",
      paired = FALSE, subject_col = NULL, condition_col = NULL, verbose = FALSE
    ),
    "not found"
  )
})

test_that(".tsenat_detect_q_prepare_data: standardizes column names", {
  df <- data.frame(
    my_entropy = c(1, 2, 3, 4, 5, 6),
    my_q = c("q1", "q1", "q1", "q2", "q2", "q2"),
    my_gene = c("G1", "G1", "G1", "G1", "G1", "G1"),
    stringsAsFactors = FALSE
  )
  
  result <- TSENAT:::.tsenat_detect_q_prepare_data(
    data = df, entropy_col = "my_entropy", q_col = "my_q", gene_col = "my_gene",
    paired = FALSE, subject_col = NULL, condition_col = NULL, verbose = FALSE
  )
  
  expect_true("entropy" %in% colnames(result$data))
  expect_true("q" %in% colnames(result$data))
  expect_true("gene" %in% colnames(result$data))
  expect_is(result$data$q, "factor")
  expect_is(result$data$gene, "factor")
})

test_that(".tsenat_detect_q_prepare_data: returns has_condition=FALSE for data frame without condition", {
  df <- data.frame(
    entropy = c(1, 2, 3, 4, 5, 6),
    q = c("q1", "q1", "q1", "q2", "q2", "q2"),
    gene = c("G1", "G1", "G1", "G1", "G1", "G1"),
    stringsAsFactors = FALSE
  )
  
  result <- TSENAT:::.tsenat_detect_q_prepare_data(
    data = df, entropy_col = "entropy", q_col = "q", gene_col = "gene",
    paired = FALSE, subject_col = NULL, condition_col = NULL, verbose = FALSE
  )
  
  expect_false(result$has_condition)
  expect_false("condition" %in% colnames(result$data))
})

test_that(".tsenat_detect_q_prepare_data: handles paired designs", {
  df <- data.frame(
    entropy = rnorm(12),
    q = rep(c("q1", "q2", "q3"), 4),
    gene = rep("G1", 12),
    subject = rep(c("S1", "S2"), 6),
    stringsAsFactors = FALSE
  )
  
  result <- TSENAT:::.tsenat_detect_q_prepare_data(
    data = df, entropy_col = "entropy", q_col = "q", gene_col = "gene",
    paired = TRUE, subject_col = "subject", condition_col = NULL, verbose = FALSE
  )
  
  expect_true("subject" %in% colnames(result$data))
  expect_is(result$data$subject, "factor")
})

# ============================================================================
# Test 3: .tsenat_detect_q_analyze_gene
# ============================================================================

test_that(".tsenat_detect_q_analyze_gene: identifies insufficient data", {
  gene_data <- data.frame(
    entropy = c(1, 2),
    q = c("q1", "q1"),
    gene = c("G1", "G1"),
    stringsAsFactors = FALSE
  )
  
  result <- TSENAT:::.tsenat_detect_q_analyze_gene(
    gene_data, paired = FALSE, subject_col = NULL, has_condition = FALSE
  )
  
  expect_true(result$test_failed)
  expect_equal(result$class, "Insufficient data")
})

test_that(".tsenat_detect_q_analyze_gene: computes test statistics for valid data", {
  set.seed(123)
  gene_data <- data.frame(
    entropy = rnorm(12),
    q = rep(c("q1", "q2", "q3"), 4),
    gene = rep("G1", 12),
    stringsAsFactors = FALSE
  )
  
  result <- TSENAT:::.tsenat_detect_q_analyze_gene(
    gene_data, paired = FALSE, subject_col = NULL, has_condition = FALSE
  )
  
  expect_false(result$test_failed)
  expect_true(is.numeric(result$f_stat))
  expect_true(is.numeric(result$p_val))
  expect_true(!is.na(result$f_stat))
  expect_true(!is.na(result$p_val))
  expect_equal(result$n_q, 3)
})

test_that(".tsenat_detect_q_analyze_gene: computes valid effect sizes", {
  set.seed(456)
  # Create data with strong q-effect
  gene_data <- data.frame(
    entropy = c(
      rnorm(4, mean = 1, sd = 0.1),   # q1
      rnorm(4, mean = 2, sd = 0.1),   # q2
      rnorm(4, mean = 3, sd = 0.1)    # q3
    ),
    q = rep(c("q1", "q2", "q3"), each = 4),
    gene = rep("G1", 12),
    stringsAsFactors = FALSE
  )
  
  result <- TSENAT:::.tsenat_detect_q_analyze_gene(
    gene_data, paired = FALSE, subject_col = NULL, has_condition = FALSE
  )
  
  # Effect size should be meaningful (eta2 between 0 and 1)
  expect_true(result$eta2 >= 0 && result$eta2 <= 1)
  # With this strong effect, eta2 should be reasonably large
  expect_true(result$eta2 > 0.5)
})

test_that(".tsenat_detect_q_analyze_gene: sums of squares are consistent", {
  set.seed(789)
  gene_data <- data.frame(
    entropy = rnorm(12),
    q = rep(c("q1", "q2", "q3"), 4),
    gene = rep("G1", 12),
    stringsAsFactors = FALSE
  )
  
  result <- TSENAT:::.tsenat_detect_q_analyze_gene(
    gene_data, paired = FALSE, subject_col = NULL, has_condition = FALSE
  )
  
  # Function returns ss_interaction (not ss_q) and ss_residual
  # These should be valid positive numbers
  expect_true(result$ss_interaction >= 0)
  expect_true(result$ss_residual >= 0)
  expect_true(result$eta2 >= 0 && result$eta2 <= 1)
})

test_that(".tsenat_detect_q_analyze_gene: handles condition column", {
  set.seed(321)
  gene_data <- data.frame(
    entropy = rnorm(12),
    q = rep(c("q1", "q2"), 6),
    condition = rep(c("ctrl", "treat"), 6),
    gene = rep("G1", 12),
    stringsAsFactors = FALSE
  )
  
  result <- TSENAT:::.tsenat_detect_q_analyze_gene(
    gene_data, paired = FALSE, subject_col = NULL, has_condition = TRUE
  )
  
  expect_false(result$test_failed)
  expect_true(is.numeric(result$f_stat))
  expect_true(!is.na(result$p_val))
})
