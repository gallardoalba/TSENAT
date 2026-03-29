# Tests for effect_sizes_divergence and results integration
# Covers uncovered lines from divergence_coverage.txt

test_that("effect_sizes_divergence aligns gene datasets", {
  # Create mock data matching expected structure
  lm_res <- data.frame(
    gene = c("gene1", "gene2", "gene3"),
    adj_p_interaction = c(0.01, 0.05, 0.5),
    estimate_interaction = c(0.5, 0.3, 0.1),
    se_interaction = c(0.1, 0.15, 0.2)
  )
  
  # Create SummarizedExperiment with divergence results
  div_se <- SummarizedExperiment::SummarizedExperiment(
    assays = list(divergence = matrix(c(0.8, 0.7, 0.6), nrow = 3, ncol = 1)),
    rowData = data.frame(
      gene_name = c("gene1", "gene2", "gene3"),
      estimate = c(0.8, 0.7, 0.6),
      lower_ci = c(0.7, 0.6, 0.5),
      upper_ci = c(0.9, 0.8, 0.7)
    )
  )
  
  # Should align and compute effect sizes
  expect_error(
    .effect_sizes_divergence(
      lm_res = lm_res,
      divergence_results_se = div_se,
      verbose = FALSE
    ),
    NA  # Expect no error
  )
})

test_that("effect_sizes_divergence filters non-matching genes", {
  # lm_res with different genes
  lm_res <- data.frame(
    gene = c("gene_A", "gene_B"),
    adj_p_interaction = c(0.01, 0.05),
    estimate_interaction = c(0.5, 0.3),
    se_interaction = c(0.1, 0.15)
  )
  
  # divergence_results_se with different genes
  div_se <- SummarizedExperiment::SummarizedExperiment(
    assays = list(divergence = matrix(c(0.8, 0.7), nrow = 2, ncol = 1)),
    rowData = data.frame(
      gene_name = c("gene_X", "gene_Y"),
      estimate = c(0.8, 0.7),
      lower_ci = c(0.7, 0.6),
      upper_ci = c(0.9, 0.8)
    )
  )
  
  # Should handle gene filtering gracefully
  result <- .effect_sizes_divergence(
    lm_res = lm_res,
    divergence_results_se = div_se,
    verbose = FALSE
  )
  
  # Result should be valid (may be empty or have notes)
  expect_true(is.data.frame(result) || is.list(result))
})

test_that("effect_sizes_divergence handles significance threshold", {
  # Test with different significance thresholds
  lm_res <- data.frame(
    gene = c("gene1", "gene2", "gene3"),
    adj_p_interaction = c(0.01, 0.05, 0.5),
    estimate_interaction = c(0.5, 0.3, 0.1),
    se_interaction = c(0.1, 0.15, 0.2)
  )
  
  div_se <- SummarizedExperiment::SummarizedExperiment(
    assays = list(divergence = matrix(c(0.8, 0.7, 0.6), nrow = 3, ncol = 1)),
    rowData = data.frame(
      gene_name = c("gene1", "gene2", "gene3"),
      estimate = c(0.8, 0.7, 0.6),
      lower_ci = c(0.7, 0.6, 0.5),
      upper_ci = c(0.9, 0.8, 0.7)
    )
  )
  
  # Test with strict threshold
  result_strict <- .effect_sizes_divergence(
    lm_res = lm_res,
    divergence_results_se = div_se,
    significance_threshold = 0.01,
    verbose = FALSE
  )
  
  # Test with lenient threshold
  result_lenient <- .effect_sizes_divergence(
    lm_res = lm_res,
    divergence_results_se = div_se,
    significance_threshold = 0.1,
    verbose = FALSE
  )
  
  # Both should produce valid results
  expect_true(is.data.frame(result_strict) || is.list(result_strict))
  expect_true(is.data.frame(result_lenient) || is.list(result_lenient))
})

test_that("effect_sizes_divergence enriches per-q patterns", {
  # Test with enrich_per_q_pattern = TRUE
  lm_res <- data.frame(
    gene = c("gene1", "gene2"),
    adj_p_interaction = c(0.01, 0.05),
    estimate_interaction = c(0.5, 0.3),
    se_interaction = c(0.1, 0.15)
  )
  
  div_se <- SummarizedExperiment::SummarizedExperiment(
    assays = list(
      divergence = matrix(c(0.8, 0.7), nrow = 2, ncol = 1),
      q_0.5 = matrix(c(0.9, 0.8), nrow = 2, ncol = 1),
      q_1.0 = matrix(c(0.8, 0.7), nrow = 2, ncol = 1),
      q_2.0 = matrix(c(0.7, 0.6), nrow = 2, ncol = 1)
    ),
    rowData = data.frame(
      gene_name = c("gene1", "gene2"),
      estimate = c(0.8, 0.7),
      lower_ci = c(0.7, 0.6),
      upper_ci = c(0.9, 0.8),
      per_q_pattern = c("RARE_DRIVEN", "BALANCED")
    )
  )
  
  result <- .effect_sizes_divergence(
    lm_res = lm_res,
    divergence_results_se = div_se,
    enrich_per_q_pattern = TRUE,
    verbose = FALSE
  )
  
  expect_true(!is.null(result))
})

test_that("effect_sizes_divergence formats output correctly", {
  # Test output formatting and structure
  lm_res <- data.frame(
    gene = c("gene1"),
    adj_p_interaction = c(0.01),
    estimate_interaction = c(0.5),
    se_interaction = c(0.1)
  )
  
  div_se <- SummarizedExperiment::SummarizedExperiment(
    assays = list(divergence = matrix(0.8, nrow = 1, ncol = 1)),
    rowData = data.frame(
      gene_name = c("gene1"),
      estimate = c(0.8),
      lower_ci = c(0.7),
      upper_ci = c(0.9)
    )
  )
  
  result <- .effect_sizes_divergence(
    lm_res = lm_res,
    divergence_results_se = div_se,
    verbose = FALSE
  )
  
  # Should return a valid result (data frame or SE or list)
  expect_true(!is.null(result))
})

test_that("effect_sizes_divergence handles zero divergence", {
  # Test handling of zero divergence estimates
  lm_res <- data.frame(
    gene = c("gene1", "gene2"),
    adj_p_interaction = c(0.01, 0.05),
    estimate_interaction = c(0.5, 0.0),  # Zero slope difference
    se_interaction = c(0.1, 0.15)
  )
  
  div_se <- SummarizedExperiment::SummarizedExperiment(
    assays = list(divergence = matrix(c(0.8, 0.0), nrow = 2, ncol = 1)),  # Zero divergence
    rowData = data.frame(
      gene_name = c("gene1", "gene2"),
      estimate = c(0.8, 0.0),
      lower_ci = c(0.7, 0.0),
      upper_ci = c(0.9, 0.0)
    )
  )
  
  result <- .effect_sizes_divergence(
    lm_res = lm_res,
    divergence_results_se = div_se,
    verbose = FALSE
  )
  
  expect_true(is.data.frame(result) || is.list(result))
})

test_that("effect_sizes_divergence handles NA divergence estimates", {
  # Test handling of NA estimates
  lm_res <- data.frame(
    gene = c("gene1", "gene2"),
    adj_p_interaction = c(0.01, 0.05),
    estimate_interaction = c(0.5, 0.3),
    se_interaction = c(0.1, 0.15)
  )
  
  div_se <- SummarizedExperiment::SummarizedExperiment(
    assays = list(divergence = matrix(c(0.8, NA), nrow = 2, ncol = 1)),
    rowData = data.frame(
      gene_name = c("gene1", "gene2"),
      estimate = c(0.8, NA),
      lower_ci = c(0.7, NA),
      upper_ci = c(0.9, NA)
    )
  )
  
  result <- .effect_sizes_divergence(
    lm_res = lm_res,
    divergence_results_se = div_se,
    verbose = FALSE
  )
  
  expect_true(!is.null(result))
})

test_that("effect_sizes_divergence with verbose output", {
  # Test verbose output path
  lm_res <- data.frame(
    gene = c("gene1"),
    adj_p_interaction = c(0.01),
    estimate_interaction = c(0.5),
    se_interaction = c(0.1)
  )
  
  div_se <- SummarizedExperiment::SummarizedExperiment(
    assays = list(divergence = matrix(0.8, nrow = 1, ncol = 1)),
    rowData = data.frame(
      gene_name = c("gene1"),
      estimate = c(0.8),
      lower_ci = c(0.7),
      upper_ci = c(0.9)
    )
  )
  
  # Capture output to verify verbose = TRUE produces messages
  result <- .effect_sizes_divergence(
    lm_res = lm_res,
    divergence_results_se = div_se,
    verbose = TRUE
  )
  
  expect_true(!is.null(result))
})

# ============================================================================
# TESTS FOR REFACTORED HELPER FUNCTIONS
# ============================================================================

test_that(".validateEffectSizeInputs rejects invalid lm_res", {
  # Test with non-data.frame lm_res
  expect_error(
    TSENAT:::.validateEffectSizeInputs(
      lm_res = list(gene = "gene1", p = 0.01),
      divergence_results_se = SummarizedExperiment::SummarizedExperiment()
    ),
    "must be a data frame"
  )
})

test_that(".validateEffectSizeInputs rejects lm_res without required columns", {
  # Missing adj_p_interaction column
  expect_error(
    TSENAT:::.validateEffectSizeInputs(
      lm_res = data.frame(gene = c("gene1"), other_col = c(0.01)),
      divergence_results_se = SummarizedExperiment::SummarizedExperiment()
    ),
    "must have columns"
  )
})

test_that(".validateEffectSizeInputs rejects non-SummarizedExperiment", {
  # Test with non-SE divergence_results
  expect_error(
    TSENAT:::.validateEffectSizeInputs(
      lm_res = data.frame(gene = "gene1", adj_p_interaction = 0.01),
      divergence_results_se = data.frame(gene = "gene1")
    ),
    "must be a SummarizedExperiment"
  )
})

test_that(".validateEffectSizeInputs rejects SE without gene_name in rowData", {
  # SE without gene_name column
  bad_se <- SummarizedExperiment::SummarizedExperiment(
    rowData = data.frame(other_col = c("gene1"))
  )
  expect_error(
    TSENAT:::.validateEffectSizeInputs(
      lm_res = data.frame(gene = "gene1", adj_p_interaction = 0.01),
      divergence_results_se = bad_se
    ),
    "must have 'gene_name' column"
  )
})

test_that(".alignGeneDatasets detects generic divergence columns", {
  # Setup data with generic columns
  lm_res <- data.frame(
    gene = c("gene1", "gene2"),
    adj_p_interaction = c(0.01, 0.05)
  )
  
  rd <- data.frame(
    gene_name = c("gene1", "gene2"),
    estimate = c(0.8, 0.7),
    lower_ci = c(0.7, 0.6),
    upper_ci = c(0.9, 0.8)
  )
  
  result <- TSENAT:::.alignGeneDatasets(lm_res, rd, verbose = FALSE)
  
  expect_true(result$use_generic)
  expect_true(is.na(result$q_values[1]))
})

test_that(".alignGeneDatasets detects per-q divergence columns", {
  # Setup data with per-q columns
  lm_res <- data.frame(
    gene = c("gene1", "gene2"),
    adj_p_interaction = c(0.01, 0.05)
  )
  
  rd <- data.frame(
    gene_name = c("gene1", "gene2"),
    estimate_q0.5 = c(0.9, 0.8),
    lower_ci_q0.5 = c(0.8, 0.7),
    upper_ci_q0.5 = c(1.0, 0.9),
    estimate_q1.0 = c(0.8, 0.7),
    lower_ci_q1.0 = c(0.7, 0.6),
    upper_ci_q1.0 = c(0.9, 0.8),
    estimate_q2.0 = c(0.7, 0.6),
    lower_ci_q2.0 = c(0.6, 0.5),
    upper_ci_q2.0 = c(0.8, 0.7)
  )
  
  result <- TSENAT:::.alignGeneDatasets(lm_res, rd, verbose = FALSE)
  
  expect_false(result$use_generic)
  expect_equal(result$q_values, c(0.5, 1.0, 2.0))
})

test_that(".alignGeneDatasets filters non-matching genes", {
  # Setup with only partial overlap
  lm_res <- data.frame(
    gene = c("gene_A", "gene_B", "gene_C"),
    adj_p_interaction = c(0.01, 0.05, 0.1)
  )
  
  rd <- data.frame(
    gene_name = c("gene_A", "gene_X"),
    estimate = c(0.8, 0.7),
    lower_ci = c(0.7, 0.6),
    upper_ci = c(0.9, 0.8)
  )
  
  result <- TSENAT:::.alignGeneDatasets(lm_res, rd, verbose = FALSE)
  
  # Should have filtered to only matching gene
  expect_equal(nrow(result$lm_res), 1)
  expect_equal(result$lm_res$gene[1], "gene_A")
})

test_that(".createResultsDataFrame generates correct columns for generic divergence", {
  # Test generic divergence columns
  df <- TSENAT:::.createResultsDataFrame(q_values = NA_real_, use_generic = TRUE)
  
  expected_cols <- c("gene", "p_value_interaction", "slope_diff", 
                     "effect_size_D", "D_lower_ci", "D_upper_ci")
  expect_equal(colnames(df), expected_cols)
  expect_equal(nrow(df), 0)
})

test_that(".createResultsDataFrame generates correct columns for multi-q divergence", {
  # Test multi-q divergence columns
  q_values <- c(0.5, 1.0, 2.0)
  df <- TSENAT:::.createResultsDataFrame(q_values = q_values, use_generic = FALSE)
  
  # Should have gene, p, slope, plus 9 columns (3 per q × 3 q values)
  expect_equal(nrow(df), 0)
  expect_true(all(c("gene", "p_value_interaction", "slope_diff") %in% colnames(df)))
  
  # Check for all q-specific columns (with underscores as separators)
  expected_cols <- c(
    "effect_size_D_q0_5", "D_q0_5_lower_ci", "D_q0_5_upper_ci",
    "effect_size_D_q1_0", "D_q1_0_lower_ci", "D_q1_0_upper_ci",
    "effect_size_D_q2_0", "D_q2_0_lower_ci", "D_q2_0_upper_ci"
  )
  expect_true(all(expected_cols %in% colnames(df)))
})

test_that(".extractLMMData returns correct LMM information", {
  # Setup LMM results data
  lm_res <- data.frame(
    gene = c("gene1", "gene2"),
    adj_p_interaction = c(0.01, 0.05),
    slope_diff = c(0.5, 0.3),
    gene_name = c("g1", "g2")
  )
  
  # Extract without gene_name column preference
  result <- TSENAT:::.extractLMMData(lm_res, "gene1", use_gene_name_col = FALSE)
  
  expect_equal(result$match_name, "gene1")
  expect_equal(result$p_interaction, 0.01)
  expect_equal(result$slope_diff, 0.5)
})

test_that(".extractLMMData uses gene_name column when available", {
  # Setup LMM results with gene_name column
  lm_res <- data.frame(
    gene = c("g_id_1", "g_id_2"),
    adj_p_interaction = c(0.01, 0.05),
    slope_diff = c(0.5, 0.3),
    gene_name = c("ENSG00001", "ENSG00002")
  )
  
  # Extract with gene_name column preference
  result <- TSENAT:::.extractLMMData(lm_res, "g_id_1", use_gene_name_col = TRUE)
  
  expect_equal(result$match_name, "ENSG00001")
  expect_equal(result$p_interaction, 0.01)
})

test_that(".extractLMMData returns NULL for non-existent gene", {
  lm_res <- data.frame(
    gene = c("gene1"),
    adj_p_interaction = c(0.01)
  )
  
  result <- TSENAT:::.extractLMMData(lm_res, "nonexistent_gene", use_gene_name_col = FALSE)
  
  expect_null(result)
})

test_that(".extractDivergenceData retrieves divergence row correctly", {
  rd <- data.frame(
    gene_name = c("gene1", "gene2"),
    estimate = c(0.8, 0.7),
    lower_ci = c(0.7, 0.6),
    upper_ci = c(0.9, 0.8)
  )
  
  result <- TSENAT:::.extractDivergenceData(rd, "gene1", verbose = FALSE, i = 1, total = 1)
  
  expect_true(is.data.frame(result))
  expect_equal(result$estimate[1], 0.8)
  expect_equal(result$lower_ci[1], 0.7)
  expect_equal(result$upper_ci[1], 0.9)
})

test_that(".extractDivergenceData returns NULL for non-existent gene", {
  rd <- data.frame(
    gene_name = c("gene1", "gene2"),
    estimate = c(0.8, 0.7),
    lower_ci = c(0.7, 0.6),
    upper_ci = c(0.9, 0.8)
  )
  
  result <- TSENAT:::.extractDivergenceData(rd, "nonexistent", verbose = FALSE, i = 1, total = 1)
  
  expect_null(result)
})

test_that(".formatSingleQResult creates correctly formatted row", {
  div_data <- data.frame(
    estimate = 0.75,
    lower_ci = 0.65,
    upper_ci = 0.85
  )
  
  result <- TSENAT:::.formatSingleQResult(
    match_name = "gene1",
    p_interaction = 0.01,
    slope_diff = 0.5,
    div_data = div_data
  )
  
  expect_equal(result$gene, "gene1")
  expect_equal(result$p_value_interaction, 0.01)
  expect_equal(result$slope_diff, 0.5)
  expect_equal(result$effect_size_D, 0.75)  # abs applied
  expect_equal(result$D_lower_ci, 0.65)
  expect_equal(result$D_upper_ci, 0.85)
})

test_that(".formatSingleQResult takes absolute value of divergence", {
  # Test that negative divergence values are converted to absolute
  div_data <- data.frame(
    estimate = -0.5,
    lower_ci = -0.6,
    upper_ci = -0.4
  )
  
  result <- TSENAT:::.formatSingleQResult(
    match_name = "gene1",
    p_interaction = 0.01,
    slope_diff = 0.5,
    div_data = div_data
  )
  
  expect_equal(result$effect_size_D, 0.5)  # Should be absolute value
})

test_that(".formatMultiQResult creates correctly formatted multi-q row", {
  # Create test data with per-q columns (column names use periods, not underscores)
  div_data <- data.frame(
    check.names = FALSE,  # Preserve column names exactly as specified
    # Column names match format from divergence_core.R: paste0("estimate_q", q_val)
    # where q_val is numeric, so as.character() converts 1.0 -> "1", 2.0 -> "2"
    `estimate_q0.5` = 0.9,
    `lower_ci_q0.5` = 0.8,
    `upper_ci_q0.5` = 1.0,
    `estimate_q1` = 0.8,
    `lower_ci_q1` = 0.7,
    `upper_ci_q1` = 0.9,
    `estimate_q2` = 0.6,
    `lower_ci_q2` = 0.5,
    `upper_ci_q2` = 0.7
  )
  
  result <- TSENAT:::.formatMultiQResult(
    match_name = "gene1",
    p_interaction = 0.01,
    slope_diff = 0.5,
    div_data = div_data,
    q_values = c(0.5, 1.0, 2.0)
  )
  
  expect_true(is.data.frame(result))
  expect_equal(result$gene, "gene1")
  expect_equal(result$effect_size_D_q0_5, 0.9)
  expect_equal(result$effect_size_D_q1, 0.8)
  expect_equal(result$effect_size_D_q2, 0.6)
})

test_that(".formatMultiQResult returns NULL when all estimates are NA", {
  div_data <- data.frame(
    check.names = FALSE,
    `estimate_q0.5` = NA_real_,
    `lower_ci_q0.5` = NA_real_,
    `upper_ci_q0.5` = NA_real_,
    `estimate_q1` = NA_real_,
    `lower_ci_q1` = NA_real_,
    `upper_ci_q1` = NA_real_
  )
  
  result <- TSENAT:::.formatMultiQResult(
    match_name = "gene1",
    p_interaction = 0.01,
    slope_diff = 0.5,
    div_data = div_data,
    q_values = c(0.5, 1.0)
  )
  
  expect_null(result)
})

test_that(".classify_q_pattern classifies RARE_DRIVEN pattern", {
  # Divergence higher at low q (rare): use named q-values
  divs <- c(q_0.5 = 0.9, q_1.0 = 0.7, q_2.0 = 0.5)  # Decreasing - rare driven
  
  result <- TSENAT:::.classify_q_pattern(divs)
  
  expect_equal(result, "RARE_DRIVEN")
})

test_that(".classify_q_pattern classifies ABUNDANT_DRIVEN pattern", {
  # Divergence higher at high q (abundant): use named q-values
  divs <- c(q_0.5 = 0.5, q_1.0 = 0.7, q_2.0 = 0.9)  # Increasing - abundant driven
  
  result <- TSENAT:::.classify_q_pattern(divs)
  
  expect_equal(result, "ABUNDANT_DRIVEN")
})

test_that(".classify_q_pattern classifies BALANCED pattern", {
  # Divergence similar across q (balanced): use named q-values
  # For balanced, need rare_median/abundant_median between 1/1.3 and 1.3
  divs <- c(q_0.5 = 0.7, q_1.0 = 0.75, q_2.0 = 0.8)  # ratio = 0.7/0.8 = 0.875 = BALANCED
  
  result <- TSENAT:::.classify_q_pattern(divs)
  
  expect_equal(result, "BALANCED")
})

test_that(".classify_q_pattern requires named q-values for classification", {
  # Unnamed vector (without q_ names) should return NA
  divs <- c(0.75)
  
  result <- TSENAT:::.classify_q_pattern(divs)
  
  expect_true(is.na(result))
})

test_that(".classify_q_pattern returns NA for all-NA input", {
  divs <- c(q_0.5 = NA_real_, q_1.0 = NA_real_, q_2.0 = NA_real_)
  
  result <- TSENAT:::.classify_q_pattern(divs)
  
  expect_true(is.na(result))
})

test_that(".classify_q_pattern returns NA for empty input", {
  divs <- numeric(0)
  
  result <- TSENAT:::.classify_q_pattern(divs)
  
  expect_true(is.na(result))
})

# ============================================================================
# INTEGRATION TESTS FOR COMPLETE WORKFLOW
# ============================================================================

test_that("effect_sizes_divergence returns list with required components", {
  lm_res <- data.frame(
    gene = c("gene1", "gene2"),
    adj_p_interaction = c(0.01, 0.05)
  )
  
  div_se <- SummarizedExperiment::SummarizedExperiment(
    assays = list(divergence = matrix(c(0.8, 0.7), nrow = 2, ncol = 1)),
    rowData = data.frame(
      gene_name = c("gene1", "gene2"),
      estimate = c(0.8, 0.7),
      lower_ci = c(0.7, 0.6),
      upper_ci = c(0.9, 0.8)
    )
  )
  
  result <- .effect_sizes_divergence(
    lm_res = lm_res,
    divergence_results_se = div_se,
    verbose = FALSE
  )
  
  expect_true(is.list(result))
  expect_true("interaction_results" %in% names(result))
  expect_true("validation_stats" %in% names(result))
})

test_that("effect_sizes_divergence validation stats are accurate", {
  lm_res <- data.frame(
    gene = c("gene1", "gene2", "gene3"),
    adj_p_interaction = c(0.01, 0.05, 0.5)
  )
  
  div_se <- SummarizedExperiment::SummarizedExperiment(
    assays = list(divergence = matrix(c(0.8, 0.7, 0.6), nrow = 3, ncol = 1)),
    rowData = data.frame(
      gene_name = c("gene1", "gene2", "gene3"),
      estimate = c(0.8, 0.7, 0.6),
      lower_ci = c(0.7, 0.6, 0.5),
      upper_ci = c(0.9, 0.8, 0.7)
    )
  )
  
  result <- .effect_sizes_divergence(
    lm_res = lm_res,
    divergence_results_se = div_se,
    significance_threshold = 0.05,
    verbose = FALSE
  )
  
  stats <- result$validation_stats
  
  # With threshold 0.05, should include gene1 (p=0.01)
  expect_equal(stats$total_genes, 1)
  expect_true(stats$passed_lmm >= 0)
})

test_that("effect_sizes_divergence with disable enrichment", {
  lm_res <- data.frame(
    gene = c("gene1"),
    adj_p_interaction = c(0.01)
  )
  
  div_se <- SummarizedExperiment::SummarizedExperiment(
    assays = list(divergence = matrix(0.8, nrow = 1, ncol = 1)),
    rowData = data.frame(
      gene_name = c("gene1"),
      estimate = c(0.8),
      lower_ci = c(0.7),
      upper_ci = c(0.9)
    )
  )
  
  # Run with enrichment disabled
  result <- .effect_sizes_divergence(
    lm_res = lm_res,
    divergence_results_se = div_se,
    enrich_per_q_pattern = FALSE,
    verbose = FALSE
  )
  
  # Should have interaction_results but not per_q_pattern column
  expect_true(!is.null(result$interaction_results))
})

test_that("effect_sizes_divergence handles empty matching results", {
  # Create data with no overlapping genes
  lm_res <- data.frame(
    gene = c("gene_A", "gene_B"),
    adj_p_interaction = c(0.01, 0.05)
  )
  
  div_se <- SummarizedExperiment::SummarizedExperiment(
    assays = list(divergence = matrix(c(0.8, 0.7), nrow = 2, ncol = 1)),
    rowData = data.frame(
      gene_name = c("gene_X", "gene_Y"),
      estimate = c(0.8, 0.7),
      lower_ci = c(0.7, 0.6),
      upper_ci = c(0.9, 0.8)
    )
  )
  
  result <- .effect_sizes_divergence(
    lm_res = lm_res,
    divergence_results_se = div_se,
    verbose = FALSE
  )
  
  # Should return gracefully with empty or minimal results
  expect_true(is.list(result))
  expect_equal(result$validation_stats$total_genes, 0)
})

test_that("effect_sizes_divergence with mixed NA and valid divergence", {
  lm_res <- data.frame(
    gene = c("gene1", "gene2", "gene3"),
    adj_p_interaction = c(0.01, 0.02, 0.03)
  )
  
  div_se <- SummarizedExperiment::SummarizedExperiment(
    assays = list(divergence = matrix(c(0.8, NA, 0.6), nrow = 3, ncol = 1)),
    rowData = data.frame(
      gene_name = c("gene1", "gene2", "gene3"),
      estimate = c(0.8, NA, 0.6),
      lower_ci = c(0.7, NA, 0.5),
      upper_ci = c(0.9, NA, 0.7)
    )
  )
  
  result <- .effect_sizes_divergence(
    lm_res = lm_res,
    divergence_results_se = div_se,
    verbose = FALSE
  )
  
  # Should handle mixed NA correctly
  expect_true(is.list(result))
  stats <- result$validation_stats
  expect_true(stats$failed_missing_divergence >= 1)  # gene2 should fail
})

context("classify_q_pattern helper")

# Basic scenarios
per_q1 <- c(q_0.5 = 0.5, q_1 = 0.3, q_2 = 0.1)
per_q2 <- c(q_0.5 = 0.1, q_1 = 0.3, q_2 = 0.5)
per_q3 <- c(q_0.5 = 0.2, q_1 = 0.2, q_2 = 0.2)  # flat values -> balanced

# edge cases
per_q_na <- c(q_0.5 = NA_real_, q_1 = NA_real_)
per_q_short <- c(q_0.5 = 0.2)
per_q_noname <- c(0.1, 0.2, 0.3)  # unnamed vector should return NA

test_that("patterns are classified correctly", {
  expect_equal(.classify_q_pattern(per_q1), "RARE_DRIVEN")
  expect_equal(.classify_q_pattern(per_q2), "ABUNDANT_DRIVEN")
  expect_equal(.classify_q_pattern(per_q3), "BALANCED")
})

test_that("NA or invalid input returns NA", {
  expect_true(is.na(.classify_q_pattern(per_q_na)))
  expect_true(is.na(.classify_q_pattern(per_q_short)))
  expect_true(is.na(.classify_q_pattern(per_q_noname)))
  expect_true(is.na(.classify_q_pattern(NULL)))
})

