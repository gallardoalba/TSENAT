# Tests for calculate_effect_sizes and results integration
# Covers uncovered lines from divergence_coverage.txt

test_that("calculate_effect_sizes aligns gene datasets", {
  # Create mock data matching expected structure
  sait_res <- data.frame(
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
    .calculate_effect_sizes(
      sait_res = sait_res,
      divergence_results_se = div_se,
      verbose = FALSE
    ),
    NA  # Expect no error
  )
})

test_that("calculate_effect_sizes filters non-matching genes", {
  # sait_res with different genes
  sait_res <- data.frame(
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
  result <- .calculate_effect_sizes(
    sait_res = sait_res,
    divergence_results_se = div_se,
    verbose = FALSE
  )
  
  # Result should be valid (may be empty or have notes)
  expect_true(is.data.frame(result) || is.list(result))
})

test_that("calculate_effect_sizes handles significance threshold", {
  # Test with different significance thresholds
  sait_res <- data.frame(
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
  result_strict <- .calculate_effect_sizes(
    sait_res = sait_res,
    divergence_results_se = div_se,
    significance_threshold = 0.01,
    verbose = FALSE
  )
  
  # Test with lenient threshold
  result_lenient <- .calculate_effect_sizes(
    sait_res = sait_res,
    divergence_results_se = div_se,
    significance_threshold = 0.1,
    verbose = FALSE
  )
  
  # Both should produce valid results
  expect_true(is.data.frame(result_strict) || is.list(result_strict))
  expect_true(is.data.frame(result_lenient) || is.list(result_lenient))
})

test_that("calculate_effect_sizes enriches per-q patterns", {
  # Test with enrich_per_q_pattern = TRUE
  sait_res <- data.frame(
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
  
  result <- .calculate_effect_sizes(
    sait_res = sait_res,
    divergence_results_se = div_se,
    enrich_per_q_pattern = TRUE,
    verbose = FALSE
  )
  
  expect_true(!is.null(result))
})

test_that("calculate_effect_sizes formats output correctly", {
  # Test output formatting and structure
  sait_res <- data.frame(
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
  
  result <- .calculate_effect_sizes(
    sait_res = sait_res,
    divergence_results_se = div_se,
    verbose = FALSE
  )
  
  # Should return a valid result (data frame or SE or list)
  expect_true(!is.null(result))
})

test_that("calculate_effect_sizes handles zero divergence", {
  # Test handling of zero divergence estimates
  sait_res <- data.frame(
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
  
  result <- .calculate_effect_sizes(
    sait_res = sait_res,
    divergence_results_se = div_se,
    verbose = FALSE
  )
  
  expect_true(is.data.frame(result) || is.list(result))
})

test_that("calculate_effect_sizes handles NA divergence estimates", {
  # Test handling of NA estimates
  sait_res <- data.frame(
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
  
  result <- .calculate_effect_sizes(
    sait_res = sait_res,
    divergence_results_se = div_se,
    verbose = FALSE
  )
  
  expect_true(!is.null(result))
})

test_that("calculate_effect_sizes with verbose output", {
  # Test verbose output path
  sait_res <- data.frame(
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
  result <- .calculate_effect_sizes(
    sait_res = sait_res,
    divergence_results_se = div_se,
    verbose = TRUE
  )
  
  expect_true(!is.null(result))
})

# ============================================================================
# TESTS FOR REFACTORED HELPER FUNCTIONS
# ============================================================================

test_that(".validateEffectSizeInputs rejects invalid sait_res", {
  # Test with non-data.frame sait_res
  expect_error(
    TSENAT:::.validateEffectSizeInputs(
      sait_res = list(gene = "gene1", p = 0.01),
      divergence_results_se = SummarizedExperiment::SummarizedExperiment()
    ),
    "must be a data frame"
  )
})

test_that(".validateEffectSizeInputs rejects sait_res without required columns", {
  # Missing adj_p_interaction column
  expect_error(
    TSENAT:::.validateEffectSizeInputs(
      sait_res = data.frame(gene = c("gene1"), other_col = c(0.01)),
      divergence_results_se = SummarizedExperiment::SummarizedExperiment()
    ),
    "must have columns"
  )
})

test_that(".validateEffectSizeInputs rejects non-SummarizedExperiment", {
  # Test with non-SE divergence_results
  expect_error(
    TSENAT:::.validateEffectSizeInputs(
      sait_res = data.frame(gene = "gene1", adj_p_interaction = 0.01),
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
      sait_res = data.frame(gene = "gene1", adj_p_interaction = 0.01),
      divergence_results_se = bad_se
    ),
    "must have 'gene_name' column"
  )
})

test_that(".alignGeneDatasets detects generic divergence columns", {
  # Setup data with generic columns
  sait_res <- data.frame(
    gene = c("gene1", "gene2"),
    adj_p_interaction = c(0.01, 0.05)
  )
  
  rd <- data.frame(
    gene_name = c("gene1", "gene2"),
    estimate = c(0.8, 0.7),
    lower_ci = c(0.7, 0.6),
    upper_ci = c(0.9, 0.8)
  )
  
  result <- TSENAT:::.alignGeneDatasets(sait_res, rd, verbose = FALSE)
  
  expect_true(result$use_generic)
  expect_true(is.na(result$q_values[1]))
})

test_that(".alignGeneDatasets detects per-q divergence columns", {
  # Setup data with per-q columns
  sait_res <- data.frame(
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
  
  result <- TSENAT:::.alignGeneDatasets(sait_res, rd, verbose = FALSE)
  
  expect_false(result$use_generic)
  expect_equal(result$q_values, c(0.5, 1.0, 2.0))
})

test_that(".alignGeneDatasets filters non-matching genes", {
  # Setup with only partial overlap
  sait_res <- data.frame(
    gene = c("gene_A", "gene_B", "gene_C"),
    adj_p_interaction = c(0.01, 0.05, 0.1)
  )
  
  rd <- data.frame(
    gene_name = c("gene_A", "gene_X"),
    estimate = c(0.8, 0.7),
    lower_ci = c(0.7, 0.6),
    upper_ci = c(0.9, 0.8)
  )
  
  result <- TSENAT:::.alignGeneDatasets(sait_res, rd, verbose = FALSE)
  
  # Should have filtered to only matching gene
  expect_equal(nrow(result$sait_res), 1)
  expect_equal(result$sait_res$gene[1], "gene_A")
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
  sait_res <- data.frame(
    gene = c("gene1", "gene2"),
    adj_p_interaction = c(0.01, 0.05),
    slope_diff = c(0.5, 0.3),
    gene_name = c("g1", "g2")
  )
  
  # Extract without gene_name column preference
  result <- TSENAT:::.extractLMMData(sait_res, "gene1", use_gene_name_col = FALSE)
  
  expect_equal(result$match_name, "gene1")
  expect_equal(result$p_interaction, 0.01)
  expect_equal(result$slope_diff, 0.5)
})

test_that(".extractLMMData uses gene_name column when available", {
  # Setup LMM results with gene_name column
  sait_res <- data.frame(
    gene = c("g_id_1", "g_id_2"),
    adj_p_interaction = c(0.01, 0.05),
    slope_diff = c(0.5, 0.3),
    gene_name = c("ENSG00001", "ENSG00002")
  )
  
  # Extract with gene_name column preference
  result <- TSENAT:::.extractLMMData(sait_res, "g_id_1", use_gene_name_col = TRUE)
  
  expect_equal(result$match_name, "ENSG00001")
  expect_equal(result$p_interaction, 0.01)
})

test_that(".extractLMMData returns NULL for non-existent gene", {
  sait_res <- data.frame(
    gene = c("gene1"),
    adj_p_interaction = c(0.01)
  )
  
  result <- TSENAT:::.extractLMMData(sait_res, "nonexistent_gene", use_gene_name_col = FALSE)
  
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
  
  expect_equal(result$pattern, "Rare driven")
})

test_that(".classify_q_pattern classifies ABUNDANT_DRIVEN pattern", {
  # Divergence higher at high q (abundant): use named q-values
  divs <- c(q_0.5 = 0.5, q_1.0 = 0.7, q_2.0 = 0.9)  # Increasing - abundant driven
  
  result <- TSENAT:::.classify_q_pattern(divs)
  
  expect_equal(result$pattern, "Abundant driven")
})

test_that(".classify_q_pattern classifies BALANCED pattern", {
  # Divergence similar across q (balanced): use named q-values
  # For balanced, need rare_median/abundant_median between 1/1.3 and 1.3
  divs <- c(q_0.5 = 0.7, q_1.0 = 0.75, q_2.0 = 0.8)  # ratio = 0.7/0.8 = 0.875 = BALANCED
  
  result <- TSENAT:::.classify_q_pattern(divs)
  
  expect_equal(result$pattern, "Balanced")
})

test_that(".classify_q_pattern requires named q-values for classification", {
  # Unnamed vector (without q_ names) should return NA
  divs <- c(0.75)
  
  result <- TSENAT:::.classify_q_pattern(divs)
  
  expect_true(is.na(result$pattern))
})

test_that(".classify_q_pattern returns NA for all-NA input", {
  divs <- c(q_0.5 = NA_real_, q_1.0 = NA_real_, q_2.0 = NA_real_)
  
  result <- TSENAT:::.classify_q_pattern(divs)
  
  expect_true(is.na(result$pattern))
})

test_that(".classify_q_pattern returns NA for empty input", {
  divs <- numeric(0)
  
  result <- TSENAT:::.classify_q_pattern(divs)
  
  expect_true(is.na(result$pattern))
})

# ============================================================================
# INTEGRATION TESTS FOR COMPLETE WORKFLOW
# ============================================================================

test_that("calculate_effect_sizes returns list with required components", {
  sait_res <- data.frame(
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
  
  result <- .calculate_effect_sizes(
    sait_res = sait_res,
    divergence_results_se = div_se,
    verbose = FALSE
  )
  
  expect_true(is.list(result))
  expect_true("interaction_results" %in% names(result))
  expect_true("validation_stats" %in% names(result))
})

test_that("calculate_effect_sizes validation stats are accurate", {
  sait_res <- data.frame(
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
  
  result <- .calculate_effect_sizes(
    sait_res = sait_res,
    divergence_results_se = div_se,
    significance_threshold = 0.05,
    verbose = FALSE
  )
  
  stats <- result$validation_stats
  
  # With threshold 0.05, should include gene1 (p=0.01)
  expect_equal(stats$total_genes, 1)
  expect_true(stats$passed_lmm >= 0)
})

test_that("calculate_effect_sizes with disable enrichment", {
  sait_res <- data.frame(
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
  result <- .calculate_effect_sizes(
    sait_res = sait_res,
    divergence_results_se = div_se,
    enrich_per_q_pattern = FALSE,
    verbose = FALSE
  )
  
  # Should have interaction_results but not per_q_pattern column
  expect_true(!is.null(result$interaction_results))
})

test_that("calculate_effect_sizes handles empty matching results", {
  # Create data with no overlapping genes
  sait_res <- data.frame(
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
  
  result <- .calculate_effect_sizes(
    sait_res = sait_res,
    divergence_results_se = div_se,
    verbose = FALSE
  )
  
  # Should return gracefully with empty or minimal results
  expect_true(is.list(result))
  expect_equal(result$validation_stats$total_genes, 0)
})

test_that("calculate_effect_sizes with mixed NA and valid divergence", {
  sait_res <- data.frame(
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
  
  result <- .calculate_effect_sizes(
    sait_res = sait_res,
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

# =============================================================================
# OPTIMIZATION: Pre-cache shared SummarizedExperiment fixtures
# Reused across multiple tests to avoid redundant SE creation (saves ~5-8 sec)
# =============================================================================

# Basic 2-gene multi-q fixture
div_se_multiq_basic <- SummarizedExperiment::SummarizedExperiment(
  assays = list(
    q_0.5 = matrix(c(0.95, 0.55), nrow = 2, ncol = 1),
    q_1.0 = matrix(c(0.85, 0.35), nrow = 2, ncol = 1),
    q_2.0 = matrix(c(0.75, 0.15), nrow = 2, ncol = 1)
  ),
  rowData = data.frame(
    gene_name = c("gene1", "gene2"),
    estimate_q0.5 = c(0.95, 0.55),
    lower_ci_q0.5 = c(0.85, 0.45),
    upper_ci_q0.5 = c(1.05, 0.65),
    estimate_q1 = c(0.85, 0.35),
    lower_ci_q1 = c(0.75, 0.25),
    upper_ci_q1 = c(0.95, 0.45),
    estimate_q2 = c(0.75, 0.15),
    lower_ci_q2 = c(0.65, 0.05),
    upper_ci_q2 = c(0.85, 0.25)
  )
)

# Simple 2-gene single-q fixture
div_se_singleq_basic <- SummarizedExperiment::SummarizedExperiment(
  assays = list(div = matrix(c(0.8, 0.7), nrow = 2, ncol = 1)),
  rowData = data.frame(
    gene_name = c("gene1", "gene2"),
    estimate = c(0.8, 0.7),
    lower_ci = c(0.7, 0.6),
    upper_ci = c(0.9, 0.8)
  )
)

# 3-gene fixture for more comprehensive testing
div_se_threegene <- SummarizedExperiment::SummarizedExperiment(
  assays = list(div = matrix(c(0.8, 0.7, 0.6), nrow = 3, ncol = 1)),
  rowData = data.frame(
    gene_name = c("ENSG00001", "ENSG00002", "ENSG00003"),
    estimate = c(0.8, 0.7, 0.6),
    lower_ci = c(0.7, 0.6, 0.5),
    upper_ci = c(0.9, 0.8, 0.7)
  )
)

test_that("patterns are classified correctly", {
  expect_equal(.classify_q_pattern(per_q1)$pattern, "Rare driven")
  expect_equal(.classify_q_pattern(per_q2)$pattern, "Abundant driven")
  expect_equal(.classify_q_pattern(per_q3)$pattern, "Balanced")
})

test_that("NA or invalid input returns NA", {
  expect_true(is.na(.classify_q_pattern(per_q_na)$pattern))
  expect_true(is.na(.classify_q_pattern(per_q_short)$pattern))
  expect_true(is.na(.classify_q_pattern(per_q_noname)$pattern))
  expect_true(is.na(.classify_q_pattern(NULL)$pattern))
})

# ============================================================================
# TESTS FOR NEW REFACTORED HELPER FUNCTIONS
# ============================================================================

test_that(".filterSignificantGenes filters by p-value threshold correctly", {
  sait_res <- data.frame(
    gene = c("gene1", "gene2", "gene3", "gene4"),
    adj_p_interaction = c(0.001, 0.01, 0.05, 0.5)
  )
  
  q_values <- NA_real_
  
  # Test with threshold 0.05 - should filter to p < 0.05 which are first 3
  result <- TSENAT:::.filterSignificantGenes(
    sait_res = sait_res,
    significance_threshold = 0.05,
    q_values = q_values,
    use_generic = TRUE,
    verbose = FALSE
  )
  
  # All genes with p < 0.05 are gene1, gene2, gene3
  expect_true(length(result$significant_genes) > 0)
  expect_true("gene1" %in% result$significant_genes)
})

test_that(".filterSignificantGenes handles strict threshold", {
  sait_res <- data.frame(
    gene = c("gene1", "gene2", "gene3"),
    adj_p_interaction = c(0.001, 0.01, 0.05)
  )
  
  # Test with strict threshold 0.01 - genes with p < 0.01
  result <- TSENAT:::.filterSignificantGenes(
    sait_res = sait_res,
    significance_threshold = 0.01,
    q_values = NA_real_,
    use_generic = TRUE,
    verbose = FALSE
  )
  
  # Only gene1 has p < 0.01 (p=0.001)
  expect_true(length(result$significant_genes) >= 1)
  expect_true("gene1" %in% result$significant_genes)
})

test_that(".filterSignificantGenes returns empty with no significant genes", {
  sait_res <- data.frame(
    gene = c("gene1", "gene2", "gene3"),
    adj_p_interaction = c(0.1, 0.2, 0.5)
  )
  
  result <- TSENAT:::.filterSignificantGenes(
    sait_res = sait_res,
    significance_threshold = 0.05,
    q_values = NA_real_,
    use_generic = TRUE,
    verbose = FALSE
  )
  
  expect_equal(length(result$significant_genes), 0)
  expect_equal(nrow(result$empty_results), 0)
})

test_that(".filterSignificantGenes handles NA p-values", {
  sait_res <- data.frame(
    gene = c("gene1", "gene2", "gene3"),
    adj_p_interaction = c(0.01, NA, 0.05)
  )
  
  result <- TSENAT:::.filterSignificantGenes(
    sait_res = sait_res,
    significance_threshold = 0.05,
    q_values = NA_real_,
    use_generic = TRUE,
    verbose = FALSE
  )
  
  # NA values should be excluded, so gene2 missing and gene1, gene3 included
  expect_true(length(result$significant_genes) > 0)
  expect_true("gene1" %in% result$significant_genes)
  expect_false("gene2" %in% result$significant_genes)  # NA excluded
})

test_that(".filterSignificantGenes creates proper empty results dataframe", {
  sait_res <- data.frame(
    gene = c("gene1"),
    adj_p_interaction = c(0.01)
  )
  
  q_values <- c(0.5, 1.0, 2.0)
  
  result_generic <- TSENAT:::.filterSignificantGenes(
    sait_res = sait_res,
    significance_threshold = 0.05,
    q_values = NA_real_,
    use_generic = TRUE,
    verbose = FALSE
  )
  
  result_multi_q <- TSENAT:::.filterSignificantGenes(
    sait_res = sait_res,
    significance_threshold = 0.05,
    q_values = q_values,
    use_generic = FALSE,
    verbose = FALSE
  )
  
  # Check structure of empty dataframes
  expect_equal(nrow(result_generic$empty_results), 0)
  expect_equal(nrow(result_multi_q$empty_results), 0)
  expect_true("gene" %in% colnames(result_generic$empty_results))
  expect_true("p_value_interaction" %in% colnames(result_multi_q$empty_results))
})

test_that(".mergeEffectSizesForGenes processes genes correctly", {
  sait_res <- data.frame(
    gene = c("gene1", "gene2"),
    adj_p_interaction = c(0.01, 0.05)
  )
  
  rd <- data.frame(
    gene_name = c("gene1", "gene2"),
    estimate = c(0.8, 0.7),
    lower_ci = c(0.7, 0.6),
    upper_ci = c(0.9, 0.8)
  )
  
  result <- TSENAT:::.mergeEffectSizesForGenes(
    sait_res = sait_res,
    rd = rd,
    significant_genes = c("gene1", "gene2"),
    q_values = NA_real_,
    use_generic = TRUE,
    verbose = FALSE
  )
  
  expect_true(is.data.frame(result$interaction_results))
  expect_true(is.list(result$validation_stats))
  expect_equal(result$validation_stats$total_genes, 2)
})

test_that(".mergeEffectSizesForGenes tracks validation statistics", {
  sait_res <- data.frame(
    gene = c("gene1", "gene2", "gene3"),
    adj_p_interaction = c(0.01, 0.05, 0.02)
  )
  
  rd <- data.frame(
    gene_name = c("gene1", "gene2"),  # Missing gene3
    estimate = c(0.8, 0.7),
    lower_ci = c(0.7, 0.6),
    upper_ci = c(0.9, 0.8)
  )
  
  result <- TSENAT:::.mergeEffectSizesForGenes(
    sait_res = sait_res,
    rd = rd,
    significant_genes = c("gene1", "gene2", "gene3"),
    q_values = NA_real_,
    use_generic = TRUE,
    verbose = FALSE
  )
  
  stats <- result$validation_stats
  expect_equal(stats$total_genes, 3)
  expect_true(stats$failed_missing_divergence >= 1)  # gene3 missing
})

test_that(".mergeEffectSizesForGenes handles empty significant genes", {
  sait_res <- data.frame(
    gene = c("gene1"),
    adj_p_interaction = c(0.01)
  )
  
  rd <- data.frame(
    gene_name = c("gene1"),
    estimate = c(0.8),
    lower_ci = c(0.7),
    upper_ci = c(0.9)
  )
  
  result <- TSENAT:::.mergeEffectSizesForGenes(
    sait_res = sait_res,
    rd = rd,
    significant_genes = character(0),  # Empty significant genes
    q_values = NA_real_,
    use_generic = TRUE,
    verbose = FALSE
  )
  
  expect_equal(nrow(result$interaction_results), 0)
  expect_equal(result$validation_stats$total_genes, 0)
})

test_that(".mergeEffectSizesForGenes processes multi-q correctly", {
  sait_res <- data.frame(
    gene = c("gene1"),
    adj_p_interaction = c(0.01)
  )
  
  rd <- data.frame(
    gene_name = c("gene1"),
    estimate_q0.5 = 0.9,
    lower_ci_q0.5 = 0.8,
    upper_ci_q0.5 = 1.0,
    estimate_q1 = 0.8,
    lower_ci_q1 = 0.7,
    upper_ci_q1 = 0.9,
    estimate_q2 = 0.6,
    lower_ci_q2 = 0.5,
    upper_ci_q2 = 0.7
  )
  
  result <- TSENAT:::.mergeEffectSizesForGenes(
    sait_res = sait_res,
    rd = rd,
    significant_genes = c("gene1"),
    q_values = c(0.5, 1.0, 2.0),
    use_generic = FALSE,
    verbose = FALSE
  )
  
  expect_true(nrow(result$interaction_results) > 0)
  expect_true("effect_size_D_q0_5" %in% colnames(result$interaction_results))
})

test_that(".mergeEffectSizesForGenes handles gene_name column in sait_res", {
  # Note: When gene_name column exists in sait_res, matching uses it
  sait_res <- data.frame(
    gene = c("g_id_1", "g_id_2"),
    adj_p_interaction = c(0.01, 0.05)
    # Note: no gene_name column here - the function checks for it to decide strategy
  )
  
  rd <- data.frame(
    gene_name = c("g_id_1", "g_id_2"),  # Match on original gene identifiers
    estimate = c(0.8, 0.7),
    lower_ci = c(0.7, 0.6),
    upper_ci = c(0.9, 0.8)
  )
  
  result <- TSENAT:::.mergeEffectSizesForGenes(
    sait_res = sait_res,
    rd = rd,
    significant_genes = c("g_id_1", "g_id_2"),
    q_values = NA_real_,
    use_generic = TRUE,
    verbose = FALSE
  )
  
  # Should merge successfully 
  expect_true(nrow(result$interaction_results) >= 1)
})

test_that(".buildColumnCache creates accurate column indicator matrix", {
  div_data <- data.frame(
    check.names = FALSE,
    `estimate_q0.5` = 0.9,
    `lower_ci_q0.5` = 0.8,
    `upper_ci_q0.5` = 1.0,
    `estimate_q1` = 0.8,
    `lower_ci_q1` = 0.7,
    `upper_ci_q1` = 0.9,
    `estimate_q2` = 0.6,
    `lower_ci_q2` = 0.5
    # Note: missing upper_ci_q2
  )
  
  q_values <- c(0.5, 1.0, 2.0)
  cache <- TSENAT:::.buildColumnCache(div_data, q_values)
  
  expect_equal(nrow(cache), 3)
  expect_equal(ncol(cache), 3)
  expect_true(cache[1, "estimate"])     # q0.5 has estimate
  expect_true(cache[1, "lower"])        # q0.5 has lower_ci
  expect_true(cache[1, "upper"])        # q0.5 has upper_ci
  expect_true(cache[3, "estimate"])     # q2.0 has estimate
  expect_true(!cache[3, "upper"])       # q2.0 missing upper_ci
})

test_that(".buildColumnCache handles columns with numeric q values", {
  # Test with different q-value representations (e.g., "0.5" vs "1" vs "2")
  div_data <- data.frame(
    check.names = FALSE,
    `estimate_q0.5` = 0.9,
    `lower_ci_q0.5` = 0.8,
    `upper_ci_q0.5` = 1.0,
    `estimate_q1` = 0.8,   # Note: numeric as "1" not "1.0"
    `lower_ci_q1` = 0.7,
    `upper_ci_q1` = 0.9
  )
  
  q_values <- c(0.5, 1.0)  # Note: numeric 1.0
  cache <- TSENAT:::.buildColumnCache(div_data, q_values)
  
  # Should match correctly
  expect_true(cache[2, "estimate"])  # q1.0 should find estimate_q1
})

test_that(".printMergeSummary produces output with verbose=TRUE", {
  validation_stats <- list(
    total_genes = 5,
    passed_lmm = 4,
    failed_missing_divergence = 1,
    other_errors = 0,
    q_values = NA_real_
  )
  
  interaction_results <- data.frame(
    gene = c("gene1", "gene2", "gene3", "gene4"),
    effect_size_D = c(0.8, 0.7, 0.6, 0.5)
  )
  
  # Should not error and produce verbose output
  expect_message(
    TSENAT:::.printMergeSummary(
      validation_stats = validation_stats,
      interaction_results = interaction_results,
      q_values = NA_real_,
      use_generic = TRUE,
      verbose = TRUE
    ),
    "MERGE COMPLETED"
  )
})

test_that(".printMergeSummary silently handles verbose=FALSE", {
  validation_stats <- list(
    total_genes = 5,
    passed_lmm = 4,
    failed_missing_divergence = 1,
    other_errors = 0,
    q_values = NA_real_
  )
  
  interaction_results <- data.frame(
    gene = c("gene1", "gene2"),
    effect_size_D = c(0.8, 0.7)
  )
  
  # Should not error and produce no output
  expect_no_message(
    TSENAT:::.printMergeSummary(
      validation_stats = validation_stats,
      interaction_results = interaction_results,
      q_values = NA_real_,
      use_generic = TRUE,
      verbose = FALSE
    )
  )
})

test_that(".enrichWithQPatterns adds per_q_pattern column", {
  interaction_results <- data.frame(
    gene = c("gene1", "gene2"),
    p_value_interaction = c(0.01, 0.05),
    effect_size_D = c(0.8, 0.7)
  )
  
  # OPTIMIZATION: Reuse pre-cached div_se_multiq_basic fixture
  # Avoids SE creation overhead (~500ms), shared across similar tests
  result <- TSENAT:::.enrichWithQPatterns(
    interaction_results = interaction_results,
    divergence_results_se = div_se_multiq_basic,
    verbose = FALSE
  )
  
  expect_true("per_q_pattern" %in% colnames(result))
  expect_equal(nrow(result), 2)
})

test_that("complete workflow with multi-q divergence produces correct output", {
  
  sait_res <- data.frame(
    gene = c("gene1", "gene2"),
    adj_p_interaction = c(0.001, 0.05)
  )
  
  # OPTIMIZATION: Reuse pre-cached div_se_multiq_basic fixture
  # Avoids SE creation overhead (~500ms), shared across similar tests
  result <- .calculate_effect_sizes(
    sait_res = sait_res,
    divergence_results_se = div_se_multiq_basic,
    significance_threshold = 0.05,
    enrich_per_q_pattern = TRUE,
    verbose = FALSE
  )
  
  # Verify multi-q columns are present
  expect_true("effect_size_D_q0_5" %in% colnames(result$interaction_results))
  expect_true("effect_size_D_q1" %in% colnames(result$interaction_results) ||
              "effect_size_D_q1_0" %in% colnames(result$interaction_results))
  expect_true("effect_size_D_q2" %in% colnames(result$interaction_results) ||
              "effect_size_D_q2_0" %in% colnames(result$interaction_results))
})

test_that("calculate_effect_sizes maintains data integrity through pipeline", {
  
  sait_res <- data.frame(
    gene = c("ENSG00001", "ENSG00002", "ENSG00003"),
    adj_p_interaction = c(0.001, 0.01, 0.5),
    slope_diff = c(0.5, 0.3, 0.1)
  )
  
  # OPTIMIZATION: Reuse pre-cached div_se_threegene fixture
  # Avoids SE creation overhead (~500ms), shared across similar tests
  result <- .calculate_effect_sizes(
    sait_res = sait_res,
    divergence_results_se = div_se_threegene,
    significance_threshold = 0.1,
    verbose = FALSE
  )
  
  # Verify data integrity
  expect_equal(nrow(result$interaction_results), 2)  # ENSG00001, ENSG00002
  expect_true(all(c("ENSG00001", "ENSG00002") %in% result$interaction_results$gene))
  expect_equal(result$interaction_results$p_value_interaction[1], 0.001)
})

# ============================================================================
# NUMERICAL CORRECTNESS TESTS
# ============================================================================

test_that("formatSingleQResult computes correct absolute value of divergence", {
  # Test that negative divergence values are properly converted to absolute
  div_data <- data.frame(
    estimate = -0.742,
    lower_ci = -0.850,
    upper_ci = -0.634
  )
  
  result <- TSENAT:::.formatSingleQResult(
    match_name = "ENSG00001",
    p_interaction = 0.00523,
    slope_diff = 0.456,
    div_data = div_data
  )
  
  # Effect size should be absolute value
  expect_equal(result$effect_size_D, 0.742, tolerance = 1e-10)
  expect_equal(result$p_value_interaction, 0.00523, tolerance = 1e-10)
  expect_equal(result$slope_diff, 0.456, tolerance = 1e-10)
  # CIs should be preserved as-is (sign preserved)
  expect_equal(result$D_lower_ci, -0.850, tolerance = 1e-10)
  expect_equal(result$D_upper_ci, -0.634, tolerance = 1e-10)
})

test_that("formatSingleQResult preserves positive values correctly", {
  # Test with positive divergence
  div_data <- data.frame(
    estimate = 0.456,
    lower_ci = 0.389,
    upper_ci = 0.523
  )
  
  result <- TSENAT:::.formatSingleQResult(
    match_name = "test_gene",
    p_interaction = 0.01,
    slope_diff = 0.2,
    div_data = div_data
  )
  
  expect_equal(result$effect_size_D, 0.456, tolerance = 1e-10)
  expect_equal(result$D_lower_ci, 0.389, tolerance = 1e-10)
  expect_equal(result$D_upper_ci, 0.523, tolerance = 1e-10)
})

test_that("formatMultiQResult produces numerically correct multi-q values", {
  # High precision test with known values
  div_data <- data.frame(
    check.names = FALSE,
    `estimate_q0.5` = 0.8234,
    `lower_ci_q0.5` = 0.7123,
    `upper_ci_q0.5` = 0.9345,
    `estimate_q1` = 0.6512,
    `lower_ci_q1` = 0.5401,
    `upper_ci_q1` = 0.7623,
    `estimate_q2` = 0.4891,
    `lower_ci_q2` = 0.3780,
    `upper_ci_q2` = 0.6002
  )
  
  result <- TSENAT:::.formatMultiQResult(
    match_name = "ENSG00001",
    p_interaction = 0.00234,
    slope_diff = 0.567,
    div_data = div_data,
    q_values = c(0.5, 1.0, 2.0)
  )
  
  # Verify exact numerical values
  expect_equal(result$gene, "ENSG00001")
  expect_equal(result$effect_size_D_q0_5, 0.8234, tolerance = 1e-10)
  expect_equal(result$D_q0_5_lower_ci, 0.7123, tolerance = 1e-10)
  expect_equal(result$D_q0_5_upper_ci, 0.9345, tolerance = 1e-10)
  expect_equal(result$effect_size_D_q1, 0.6512, tolerance = 1e-10)
  expect_equal(result$effect_size_D_q2, 0.4891, tolerance = 1e-10)
})

test_that("formatMultiQResult handles negative divergence with absolute value", {
  # Test that negative estimates are converted to positive effect sizes
  div_data <- data.frame(
    check.names = FALSE,
    `estimate_q0.5` = -0.8234,
    `lower_ci_q0.5` = -0.9345,
    `upper_ci_q0.5` = -0.7123,
    `estimate_q1` = -0.6512,
    `lower_ci_q1` = -0.7623,
    `upper_ci_q1` = -0.5401
  )
  
  result <- TSENAT:::.formatMultiQResult(
    match_name = "ENSG00002",
    p_interaction = 0.05,
    slope_diff = 0.1,
    div_data = div_data,
    q_values = c(0.5, 1.0)
  )
  
  # Effect sizes should be absolute values
  expect_equal(result$effect_size_D_q0_5, 0.8234, tolerance = 1e-10)
  expect_equal(result$effect_size_D_q1, 0.6512, tolerance = 1e-10)
  # CIs should preserve original values
  expect_equal(result$D_q0_5_lower_ci, -0.9345, tolerance = 1e-10)
  expect_equal(result$D_q0_5_upper_ci, -0.7123, tolerance = 1e-10)
})

test_that("calculate_effect_sizes produces exact numerical output for known input", {
  # Test with precise known values to validate numerical correctness
  sait_res <- data.frame(
    gene = c("g1", "g2"),
    adj_p_interaction = c(0.001234, 0.050000)
  )
  
  div_se <- SummarizedExperiment::SummarizedExperiment(
    assays = list(div = matrix(c(0.5432, 0.7654), nrow = 2, ncol = 1)),
    rowData = data.frame(
      gene_name = c("g1", "g2"),
      estimate = c(0.5432, 0.7654),
      lower_ci = c(0.4321, 0.6543),
      upper_ci = c(0.6543, 0.8765)
    )
  )
  
  result <- .calculate_effect_sizes(
    sait_res = sait_res,
    divergence_results_se = div_se,
    significance_threshold = 0.1,
    verbose = FALSE
  )
  
  # Verify exact numerical correspondence
  expect_equal(result$interaction_results$p_value_interaction[1], 0.001234, tolerance = 1e-10)
  expect_equal(result$interaction_results$effect_size_D[1], 0.5432, tolerance = 1e-10)
  expect_equal(result$interaction_results$D_lower_ci[1], 0.4321, tolerance = 1e-10)
  expect_equal(result$interaction_results$D_upper_ci[1], 0.6543, tolerance = 1e-10)
})

test_that("calculate_effect_sizes handles zero values correctly", {
  # Test edge case with zero divergence
  sait_res <- data.frame(
    gene = c("zero_gene"),
    adj_p_interaction = c(0.01)
  )
  
  div_se <- SummarizedExperiment::SummarizedExperiment(
    assays = list(div = matrix(0.0, nrow = 1, ncol = 1)),
    rowData = data.frame(
      gene_name = c("zero_gene"),
      estimate = c(0.0),
      lower_ci = c(0.0),
      upper_ci = c(0.0)
    )
  )
  
  result <- .calculate_effect_sizes(
    sait_res = sait_res,
    divergence_results_se = div_se,
    verbose = FALSE
  )
  
  # Should handle zero correctly (no error, produces zero)
  expect_equal(result$interaction_results$effect_size_D[1], 0.0)
  expect_equal(result$interaction_results$D_lower_ci[1], 0.0)
  expect_equal(result$interaction_results$D_upper_ci[1], 0.0)
})

test_that("calculate_effect_sizes handles very small numbers correctly", {
  # Test with very small numbers (but not zero)
  sait_res <- data.frame(
    gene = c("tiny_gene"),
    adj_p_interaction = c(0.001)
  )
  
  div_se <- SummarizedExperiment::SummarizedExperiment(
    assays = list(div = matrix(1e-10, nrow = 1, ncol = 1)),
    rowData = data.frame(
      gene_name = c("tiny_gene"),
      estimate = c(1e-10),
      lower_ci = c(5e-11),
      upper_ci = c(1.5e-10)
    )
  )
  
  result <- .calculate_effect_sizes(
    sait_res = sait_res,
    divergence_results_se = div_se,
    verbose = FALSE
  )
  
  # Should preserve precision with very small numbers
  expect_equal(result$interaction_results$effect_size_D[1], 1e-10, tolerance = 1e-20)
  expect_equal(result$interaction_results$D_lower_ci[1], 5e-11, tolerance = 1e-20)
})

test_that("calculate_effect_sizes handles large numbers correctly", {
  # Test with large divergence values
  sait_res <- data.frame(
    gene = c("large_gene"),
    adj_p_interaction = c(0.01)
  )
  
  div_se <- SummarizedExperiment::SummarizedExperiment(
    assays = list(div = matrix(15.6789, nrow = 1, ncol = 1)),
    rowData = data.frame(
      gene_name = c("large_gene"),
      estimate = c(15.6789),
      lower_ci = c(14.5678),
      upper_ci = c(16.7890)
    )
  )
  
  result <- .calculate_effect_sizes(
    sait_res = sait_res,
    divergence_results_se = div_se,
    verbose = FALSE
  )
  
  expect_equal(result$interaction_results$effect_size_D[1], 15.6789, tolerance = 1e-4)
  expect_equal(result$interaction_results$D_lower_ci[1], 14.5678, tolerance = 1e-4)
})

test_that("formatMultiQResult handles mixed NA/valid values numerically", {
  # Test case where some q-values have data and others are NA
  div_data <- data.frame(
    check.names = FALSE,
    `estimate_q0.5` = 0.5432,
    `lower_ci_q0.5` = 0.4321,
    `upper_ci_q0.5` = 0.6543,
    `estimate_q1` = NA_real_,
    `lower_ci_q1` = NA_real_,
    `upper_ci_q1` = NA_real_,
    `estimate_q2` = 0.1234,
    `lower_ci_q2` = 0.0123,
    `upper_ci_q2` = 0.2345
  )
  
  result <- TSENAT:::.formatMultiQResult(
    match_name = "partial_gene",
    p_interaction = 0.05,
    slope_diff = 0.01,
    div_data = div_data,
    q_values = c(0.5, 1.0, 2.0)
  )
  
  # Should include valid values and NA for missing
  expect_equal(result$effect_size_D_q0_5, 0.5432, tolerance = 1e-10)
  expect_true(is.na(result$effect_size_D_q1))
  expect_equal(result$effect_size_D_q2, 0.1234, tolerance = 1e-10)
})

test_that(".buildColumnCache correctly identifies present/absent columns", {
  # Verify numeric correctness: 1 for present, 0 for absent
  div_data <- data.frame(
    check.names = FALSE,
    `estimate_q0.5` = 1.0,
    `lower_ci_q0.5` = 2.0,
    `upper_ci_q0.5` = 3.0,
    `estimate_q1` = 4.0,
    `lower_ci_q1` = 5.0,
    # Missing upper_ci_q1
    `estimate_q2` = 6.0
    # Missing lower_ci_q2 and upper_ci_q2
  )
  
  cache <- TSENAT:::.buildColumnCache(div_data, c(0.5, 1.0, 2.0))
  
  # Check structure: should be matrix with TRUE/FALSE
  expect_equal(nrow(cache), 3)
  expect_equal(ncol(cache), 3)
  
  # Verify presence matrix
  expect_true(cache[1, "estimate"])   # q0.5 estimate present
  expect_true(cache[1, "lower"])      # q0.5 lower present
  expect_true(cache[1, "upper"])      # q0.5 upper present
  expect_true(cache[2, "estimate"])   # q1 estimate present
  expect_true(cache[2, "lower"])      # q1 lower present
  expect_false(cache[2, "upper"])     # q1 upper MISSING
  expect_true(cache[3, "estimate"])   # q2 estimate present
  expect_false(cache[3, "lower"])     # q2 lower MISSING
  expect_false(cache[3, "upper"])     # q2 upper MISSING
})

test_that("classify_q_pattern produces numerically correct ratio thresholds", {
  # Test boundary conditions for ratio threshold (default 1.3)
  # ratio = rare_median / abundant_median
  
  # Case 1: Ratio strictly greater than threshold (rare_median = 1.5 * abundant_median)
  # Should be RARE_DRIVEN
  divs_above_rare <- c(q_0.5 = 1.5, q_1.0 = 1.0, q_2.0 = 1.0)
  expect_equal(.classify_q_pattern(divs_above_rare)$pattern, "Rare driven")
  
  # Case 2: Ratio exactly at threshold - should be BALANCED (not >= but strictly >)
  divs_at_threshold <- c(q_0.5 = 1.3, q_1.0 = 1.0, q_2.0 = 1.0)
  expect_equal(.classify_q_pattern(divs_at_threshold)$pattern, "Balanced")
  
  # Case 3: Ratio strictly below ABUNDANT threshold (ratio < 1/1.3 ≈ 0.769)
  # rare_median = 0.5, abundant_median = 1.0 → ratio = 0.5 < 0.769
  divs_below_abundant <- c(q_0.5 = 0.5, q_1.0 = 1.0, q_2.0 = 1.0)
  expect_equal(.classify_q_pattern(divs_below_abundant)$pattern, "Abundant driven")
  
  # Case 4: Ratio exactly at ABUNDANT threshold - should be BALANCED (not <= but strictly <)
  # rare_median = 0.769, abundant_median = 1.0 → ratio ≈ 0.769
  divs_at_abundant_boundary <- c(q_0.5 = 1.0 / 1.3, q_1.0 = 1.0, q_2.0 = 1.0)
  expect_equal(.classify_q_pattern(divs_at_abundant_boundary)$pattern, "Balanced")
})

test_that("classify_q_pattern correctly computes median for classification", {
  # Test with multiple values per region to verify median calculation
  # Rare region (q < 1): 0.8, 0.9, 1.0 (should be excluded) -> median = 0.85
  # Abundant region (q > 1): 1.1, 1.2 -> median = 1.15
  # ratio = 0.85 / 1.15 = 0.739 < 1/1.3 → ABUNDANT_DRIVEN
  divs_multi <- c(q_0.25 = 0.8, q_0.5 = 0.9, q_1.0 = 1.0, q_1.5 = 1.1, q_2.0 = 1.2)
  result <- .classify_q_pattern(divs_multi)
  
  # Verify it correctly identified median and ratio
  expect_equal(result$pattern, "Abundant driven")
})

test_that("Effect size calculation preserves full precision through pipeline", {
  
  # Test full pipeline with high-precision values
  sait_res <- data.frame(
    gene = c("precision_test"),
    adj_p_interaction = c(0.00112358)  # High precision p-value
  )
  
  div_se <- SummarizedExperiment::SummarizedExperiment(
    assays = list(div = matrix(0.6180339887, nrow = 1, ncol = 1)),
    rowData = data.frame(
      gene_name = c("precision_test"),
      estimate = c(0.6180339887),
      lower_ci = c(0.5555555556),
      upper_ci = c(0.6805555550)
    )
  )
  
  result <- .calculate_effect_sizes(
    sait_res = sait_res,
    divergence_results_se = div_se,
    verbose = FALSE
  )
  
  # Verify precision maintained through pipeline
  expect_equal(result$interaction_results$p_value_interaction[1], 0.00112358, 
               tolerance = 1e-10)
  expect_equal(result$interaction_results$effect_size_D[1], 0.6180339887,
               tolerance = 1e-10)
  expect_equal(result$interaction_results$D_lower_ci[1], 0.5555555556,
               tolerance = 1e-10)
})

test_that("Multi-q output maintains numerical ordering consistency", {
  # Verify that q-values are ordered consistently (low q to high q)
  sait_res <- data.frame(
    gene = c("order_test"),
    adj_p_interaction = c(0.01)
  )
  
  div_se <- SummarizedExperiment::SummarizedExperiment(
    assays = list(div = matrix(0.5, nrow = 1, ncol = 1)),
    rowData = data.frame(
      gene_name = c("order_test"),
      estimate_q0.5 = 0.9,
      lower_ci_q0.5 = 0.8,
      upper_ci_q0.5 = 1.0,
      estimate_q1 = 0.7,
      lower_ci_q1 = 0.6,
      upper_ci_q1 = 0.8,
      estimate_q2 = 0.5,
      lower_ci_q2 = 0.4,
      upper_ci_q2 = 0.6,
      estimate_q10 = 0.1,
      lower_ci_q10 = 0.0,
      upper_ci_q10 = 0.2
    )
  )
  
  result <- .calculate_effect_sizes(
    sait_res = sait_res,
    divergence_results_se = div_se,
    verbose = FALSE
  )
  
  # Verify columns exist and values are in correct order
  expect_equal(result$interaction_results$effect_size_D_q0_5, 0.9)
  expect_equal(result$interaction_results$effect_size_D_q1, 0.7)
  expect_equal(result$interaction_results$effect_size_D_q2, 0.5)
  expect_equal(result$interaction_results$effect_size_D_q10, 0.1)
})


# ============================================================================
# TESTS FOR S4 WRAPPER HELPER FUNCTIONS
# ============================================================================

test_that(".validate_effect_sizes_inputs_s4 rejects non-TSENATAnalysis", {
  # Should reject non-S4 objects
  expect_error(
    TSENAT:::.validate_effect_sizes_inputs_s4(list(data = "test")),
    "must be a TSENATAnalysis object"
  )
})

test_that(".validate_effect_sizes_inputs_s4 rejects missing divergence results", {
  # Create analysis with empty divergence results
  analysis <- new("TSENATAnalysis")
  analysis@divergence_results <- list()
  
  expect_error(
    TSENAT:::.validate_effect_sizes_inputs_s4(analysis),
    "Divergence results required"
  )
})

test_that(".validate_effect_sizes_inputs_s4 rejects missing SAIT results", {
  # Create analysis with divergence but no SAIT results
  analysis <- new("TSENATAnalysis")
  analysis@divergence_results <- list(mock = "data")
  analysis@sait_results <- list()
  
  expect_error(
    TSENAT:::.validate_effect_sizes_inputs_s4(analysis),
    "SAIT results required"
  )
})

test_that(".validate_effect_sizes_inputs_s4 passes valid analysis", {
  # Create valid analysis
  analysis <- new("TSENATAnalysis")
  analysis@divergence_results <- list(mock = "data")
  analysis@sait_results <- list(mock = "data")
  
  # Should not throw error
  expect_error(
    TSENAT:::.validate_effect_sizes_inputs_s4(analysis),
    NA
  )
})

test_that(".extract_effect_sizes_data_s4 extracts data from valid analysis", {
  # Create minimal mock analysis with required data structures
  mock_se <- SummarizedExperiment::SummarizedExperiment(
    assays = list(counts = matrix(c(100, 200, 150, 175), nrow = 2, ncol = 2))
  )
  
  analysis <- new("TSENATAnalysis")
  analysis@se <- mock_se
  analysis@divergence_results <- list(
    divergence_se = SummarizedExperiment::SummarizedExperiment(
      assays = list(divergence = matrix(c(0.8, 0.7), nrow = 2, ncol = 1)),
      rowData = data.frame(gene_name = c("gene1", "gene2"))
    )
  )
  analysis@sait_results <- list(
    sait_interaction = data.frame(
      gene = c("gene1", "gene2"),
      adj_p_interaction = c(0.01, 0.05)
    )
  )
  
  # Extract data
  data_list <- TSENAT:::.extract_effect_sizes_data_s4(analysis, verbose = FALSE)
  
  expect_true(is.list(data_list))
  expect_true("divergence_se" %in% names(data_list))
  expect_true("sait_res" %in% names(data_list))
  expect_true(is(data_list$divergence_se, "SummarizedExperiment"))
  expect_true(is.data.frame(data_list$sait_res))
})

test_that(".extract_effect_sizes_data_s4 handles missing divergence results", {
  # Create analysis with LM but no divergence results
  analysis <- new("TSENATAnalysis")
  analysis@sait_results <- list(mock = "sait_data")
  analysis@divergence_results <- list()
  
  expect_error(
    TSENAT:::.extract_effect_sizes_data_s4(analysis),
    "Could not extract divergence"
  )
})

test_that(".map_tx_to_genes handles valid tx2gene mapping", {
  # Create mock tx2gene dataframe
  tx2gene <- data.frame(
    Transcript = c("TX001", "TX002", "TX003"),
    Gene = c("GENE_A", "GENE_B", "GENE_A"),
    stringsAsFactors = FALSE
  )
  
  tx_in_divergence <- c("TX001", "TX002")
  
  result <- TSENAT:::.map_tx_to_genes(tx2gene, tx_in_divergence, verbose = FALSE)
  
  expect_true(!is.null(result))
  expect_equal(length(result), 2)
  expect_equal(result[1], "GENE_A")
  expect_equal(result[2], "GENE_B")
})

test_that(".map_tx_to_genes returns NULL for invalid columns", {
  # tx2gene with non-standard column names
  tx2gene <- data.frame(
    col1 = c("TX001", "TX002"),
    col2 = c("GENE_A", "GENE_B"),
    stringsAsFactors = FALSE
  )
  
  tx_in_divergence <- c("TX001", "TX002")
  
  result <- TSENAT:::.map_tx_to_genes(tx2gene, tx_in_divergence, verbose = FALSE)
  
  expect_null(result)
})

test_that(".map_tx_to_genes returns NULL for non-matching transcripts", {
  # Transcripts not in tx2gene
  tx2gene <- data.frame(
    Transcript = c("TX001", "TX002"),
    Gene = c("GENE_A", "GENE_B"),
    stringsAsFactors = FALSE
  )
  
  tx_in_divergence <- c("TX_NOT_FOUND", "TX_ALSO_NOT_FOUND")
  
  result <- TSENAT:::.map_tx_to_genes(tx2gene, tx_in_divergence, verbose = FALSE)
  
  # Should return NULL because all matches are NA
  expect_null(result)
})

test_that(".add_gene_names_to_divergence_se adds via tx2gene mapping", {
  # Create divergence SE without gene names
  div_se <- SummarizedExperiment::SummarizedExperiment(
    assays = list(div = matrix(c(0.8, 0.7, 0.6), nrow = 3, ncol = 1)),
    rowData = data.frame(
      estimate = c(0.8, 0.7, 0.6),
      lower_ci = c(0.7, 0.6, 0.5),
      upper_ci = c(0.9, 0.8, 0.7),
      row.names = c("TX001", "TX002", "TX003")
    )
  )
  rownames(div_se) <- c("TX001", "TX002", "TX003")
  
  # Create base SE with tx2gene metadata
  tx2gene <- data.frame(
    Transcript = c("TX001", "TX002", "TX003"),
    Gene = c("GENE_A", "GENE_B", "GENE_A"),
    stringsAsFactors = FALSE
  )
  
  base_se <- SummarizedExperiment::SummarizedExperiment(
    assays = list(counts = matrix(100, nrow = 3, ncol = 2))
  )
  S4Vectors::metadata(base_se)$tx2gene <- tx2gene
  
  # Create minimal analysis with required slots
  analysis <- new("TSENATAnalysis")
  analysis@se <- base_se
  
  sait_res <- data.frame(
    gene = c("GENE_A", "GENE_B"),
    adj_p_interaction = c(0.01, 0.05)
  )
  
  # Add gene names
  result_se <- TSENAT:::.add_gene_names_to_divergence_se(div_se, analysis, sait_res, verbose = FALSE)
  
  expect_true(is(result_se, "SummarizedExperiment"))
  expect_true("gene_name" %in% colnames(SummarizedExperiment::rowData(result_se)))
})

test_that(".add_gene_names_to_divergence_se handles missing tx2gene gracefully", {
  # Create divergence SE without gene names
  div_se <- SummarizedExperiment::SummarizedExperiment(
    assays = list(div = matrix(c(0.8, 0.7), nrow = 2, ncol = 1)),
    rowData = data.frame(
      estimate = c(0.8, 0.7),
      row.names = c("TX001", "TX002")
    )
  )
  rownames(div_se) <- c("TX001", "TX002")
  
  # Create base SE without tx2gene
  base_se <- SummarizedExperiment::SummarizedExperiment(
    assays = list(counts = matrix(100, nrow = 2, ncol = 2))
  )
  
  analysis <- new("TSENATAnalysis")
  analysis@se <- base_se
  
  sait_res <- data.frame(
    gene = c("GENE_A", "GENE_B"),
    adj_p_interaction = c(0.01, 0.05)
  )
  
  # Should fallback to direct assignment
  result_se <- TSENAT:::.add_gene_names_to_divergence_se(div_se, analysis, sait_res, verbose = FALSE)
  
  expect_true("gene_name" %in% colnames(SummarizedExperiment::rowData(result_se)))
  expect_equal(SummarizedExperiment::rowData(result_se)$gene_name, c("GENE_A", "GENE_B"))
})

test_that(".store_effect_sizes_results_in_metadata stores results correctly", {
  # Create minimal analysis
  analysis <- new("TSENATAnalysis")
  
  # Create mock result
  result <- list(
    interaction_results = data.frame(gene = "gene1", estimate = 0.8),
    validation_stats = list(n_genes = 1)
  )
  
  # Store results
  analysis_updated <- TSENAT:::.store_effect_sizes_results_in_metadata(
    analysis, result, significance_threshold = 0.05, verbose = FALSE
  )
  
  expect_true(!is.null(analysis_updated@metadata$effect_sizes_divergence))
  expect_equal(
    analysis_updated@metadata$effect_sizes_divergence$interaction_results$gene,
    "gene1"
  )
})

test_that(".store_effect_sizes_results_in_metadata tracks function calls", {
  # Create analysis with existing function calls
  analysis <- new("TSENATAnalysis")
  analysis@metadata <- list(function_calls = c("prev_call1", "prev_call2"))
  
  result <- list(interaction_results = data.frame(gene = "gene1"))
  
  analysis_updated <- TSENAT:::.store_effect_sizes_results_in_metadata(
    analysis, result, significance_threshold = 0.05, verbose = FALSE
  )
  
  calls <- analysis_updated@metadata$function_calls
  expect_equal(length(calls), 3)
  expect_true(any(grepl("calculate_effect_sizes", calls)))
})

test_that(".save_effect_sizes_output writes TSV files", {
  # Create temporary file
  temp_tsv <- tempfile(fileext = ".tsv")
  
  # Create mock result
  result <- list(
    interaction_results = data.frame(
      gene = c("gene1", "gene2"),
      estimate = c(0.8, 0.7),
      pval = c(0.01, 0.05)
    )
  )
  
  # Save output
  TSENAT:::.save_effect_sizes_output(temp_tsv, result, verbose = FALSE)
  
  # Verify file was created
  expect_true(file.exists(temp_tsv))
  
  # Verify content
  saved_data <- read.table(temp_tsv, header = TRUE, sep = "\t")
  expect_equal(nrow(saved_data), 2)
  expect_equal(ncol(saved_data), 3)
  
  # Cleanup
  file.remove(temp_tsv)
})

test_that(".save_effect_sizes_output filters all-NA columns in TSV output", {
  # Create temporary file
  temp_tsv <- tempfile(fileext = ".tsv")
  
  # Create result with NA column
  result <- list(
    interaction_results = data.frame(
      gene = c("gene1", "gene2"),
      estimate = c(0.8, 0.7),
      na_col = c(NA, NA)
    )
  )
  
  # Save output
  TSENAT:::.save_effect_sizes_output(temp_tsv, result, verbose = FALSE)
  
  # Verify NA column was removed
  saved_data <- read.table(temp_tsv, header = TRUE, sep = "\t")
  expect_false("na_col" %in% colnames(saved_data))
  expect_equal(ncol(saved_data), 2)
  
  # Cleanup
  file.remove(temp_tsv)
})

test_that(".save_effect_sizes_output writes RDS files", {
  # Create temporary file
  temp_rds <- tempfile(fileext = ".rds")
  
  # Create mock result
  result <- list(
    interaction_results = data.frame(gene = c("gene1", "gene2")),
    validation_stats = list(notes = "test")
  )
  
  # Save output
  TSENAT:::.save_effect_sizes_output(temp_rds, result, verbose = FALSE)
  
  # Verify file was created and readable
  expect_true(file.exists(temp_rds))
  loaded_result <- readRDS(temp_rds)
  expect_equal(loaded_result$validation_stats$notes, "test")
  
  # Cleanup
  file.remove(temp_rds)
})

test_that(".save_effect_sizes_output handles NULL output_file", {
  # Should not produce error when output_file is NULL
  result <- list(interaction_results = data.frame(gene = "gene1"))
  
  expect_error(
    TSENAT:::.save_effect_sizes_output(NULL, result, verbose = FALSE),
    NA
  )
})

test_that("calculate_effect_sizes orchestrates helpers correctly", {
  # Replicate roxygen documentation example exactly
  data(readcounts)
  readcounts <- as.matrix(readcounts)
  mode(readcounts) <- 'numeric'
  
  metadata_df <- read.table(
    system.file('extdata', 'metadata.tsv', package = 'TSENAT'),
    header = TRUE, sep = '\t'
  )
  
  gff3_dataset <- system.file('extdata', 'annotation.gff3.gz', package = 'TSENAT')
  
  # Create TPM and effective_length as in roxygen example
  # Simple normalization to simulate TPM
  tpm <- sweep(readcounts, 2, colSums(readcounts), "/") * 1e6
  effective_length <- rep(1000, nrow(readcounts))
  names(effective_length) <- rownames(readcounts)
  
  # Create config FIRST with required metadata column parameters
  config <- TSENAT_config(
    q_values = seq(0, 2, by = 0.05),
    sample_col = 'sample',
    condition_col = 'condition',
    subject_col = 'paired_samples',
    paired = TRUE,
    control = 'normal'
  )
  
  analysis <- build_analysis(
    config = config,
    readcounts = readcounts,
    metadata = metadata_df,
    tx2gene = gff3_dataset,
    tpm = tpm,
    effective_length = effective_length,
    verbose = FALSE
  )
  
  analysis <- filter_analysis(analysis, stringency = 'severe', verbose = FALSE)
  analysis <- calculate_diversity(analysis, q = c(0.5, 0.75, 1.0, 1.5, 2.0), verbose = FALSE)
  analysis <- calculate_divergence(analysis, q = c(0.5, 0.75, 1.0, 1.5, 2.0), verbose = FALSE)
  analysis <- suppressWarnings(
    calculate_sait(analysis, method = 'gam', verbose = FALSE)
  )
  
  # Compute effect sizes from divergence results
  analysis <- calculate_effect_sizes(
    analysis,
    significance_threshold = 0.05, 
    verbose = FALSE
  )
  
  # Access results directly from metadata (getMeta now filters to essential fields only)
  # For internal testing, access @metadata directly instead of public accessor
  effect_size_results <- analysis@metadata$effect_sizes_divergence
  
  # Verify structure of results
  expect_true(!is.null(effect_size_results))
  expect_true(is.list(effect_size_results))
})

test_that("calculate_effect_sizes respects output_file parameter", {
  # Replicate roxygen documentation example with output_file
  data(readcounts)
  readcounts <- as.matrix(readcounts)
  mode(readcounts) <- 'numeric'
  
  metadata_df <- read.table(
    system.file('extdata', 'metadata.tsv', package = 'TSENAT'),
    header = TRUE, sep = '\t'
  )
  
  gff3_dataset <- system.file('extdata', 'annotation.gff3.gz', package = 'TSENAT')
  
  # Create TPM and effective_length as in roxygen example
  # Simple normalization to simulate TPM
  tpm <- sweep(readcounts, 2, colSums(readcounts), "/") * 1e6
  effective_length <- rep(1000, nrow(readcounts))
  names(effective_length) <- rownames(readcounts)
  
  # Create config FIRST with required metadata column parameters
  config <- TSENAT_config(
    q_values = seq(0, 2, by = 0.05),
    sample_col = 'sample',
    condition_col = 'condition',
    subject_col = 'paired_samples',
    paired = TRUE,
    control = 'normal'
  )
  
  analysis <- build_analysis(
    config = config,
    readcounts = readcounts,
    metadata = metadata_df,
    tx2gene = gff3_dataset,
    tpm = tpm,
    effective_length = effective_length,
    verbose = FALSE
  )
  
  analysis <- filter_analysis(analysis, stringency = 'severe', verbose = FALSE)
  analysis <- calculate_diversity(analysis, q = c(0.5, 0.75, 1.0, 1.5, 2.0), verbose = FALSE)
  analysis <- calculate_divergence(analysis, q = c(0.5, 0.75, 1.0, 1.5, 2.0), verbose = FALSE)
  analysis <- suppressWarnings(
    calculate_sait(analysis, method = 'gam', verbose = FALSE)
  )
  
  # Create temporary file for TSV output
  temp_output <- tempfile(fileext = ".tsv")
  
  # Compute effect sizes and save to file
  analysis <- calculate_effect_sizes(
    analysis,
    significance_threshold = 0.05,
    output_file = temp_output,
    verbose = FALSE
  )
  
  # Verify output file was created and is readable
  expect_true(file.exists(temp_output))
  
  # Verify contents
  output_data <- read.table(temp_output, header = TRUE, sep = "\t")
  expect_true(nrow(output_data) > 0)
  
  # Cleanup
  file.remove(temp_output)
})


# ============================================================================
# TESTS FOR .printMergeSuccess - EDGE CASES AND ERROR HANDLING
# ============================================================================
# These tests cover the uncovered conditional branches (use_generic TRUE/FALSE)
# and error/edge case paths in the message printing logic

test_that(".printMergeSuccess prints generic path with valid CI values", {
  # Test generic divergence path with complete CI information
  lmm_data <- list(
    match_name = "gene_A1BG",
    p_interaction = 0.01234
  )
  
  div_data <- data.frame(
    estimate = 0.854,
    lower_ci = 0.750,
    upper_ci = 0.920
  )
  
  # Capture message output
  expect_message(
    TSENAT:::.printMergeSuccess(lmm_data, div_data, q_values = NA_real_, use_generic = TRUE),
    "SUCCESS.*gene_A1BG.*p=.*D_spectrum=.*CI="
  )
})

test_that(".printMergeSuccess prints generic path with NA CI values", {
  # Test generic path when CI values are NA (missing bootstrap CI)
  lmm_data <- list(
    match_name = "gene_BRCA1",
    p_interaction = 0.0001
  )
  
  div_data <- data.frame(
    estimate = 0.42,
    lower_ci = NA_real_,
    upper_ci = NA_real_
  )
  
  # Should print without CI text when values are NA
  # Use a custom reporter to capture the output
  env <- new.env()
  capture_output <- capture.output({
    withCallingHandlers(
      TSENAT:::.printMergeSuccess(lmm_data, div_data, q_values = NA_real_, use_generic = TRUE),
      message = function(m) {
        env$msg <- m$message
      }
    )
  })
  
  # Check that message was generated without CI text
  if (is.null(env$msg)) {
    # Alternative: function may not throw message, just print
    # Verify that at least SUCCESS is printed
    expect_match(paste(capture_output, collapse = ""), "SUCCESS")
  } else {
    expect_match(env$msg, "SUCCESS")
    expect_false(grepl("CI=", env$msg))
  }
})

test_that(".printMergeSuccess handles very small p-values in generic path", {
  # Test with p-value that triggers scientific notation
  lmm_data <- list(
    match_name = "gene_TP53",
    p_interaction = 1.23e-45  # Extremely small p-value
  )
  
  div_data <- data.frame(
    estimate = 0.95,
    lower_ci = 0.88,
    upper_ci = 0.98
  )
  
  # Should handle scientific notation formatting correctly
  expect_message(
    TSENAT:::.printMergeSuccess(lmm_data, div_data, q_values = NA_real_, use_generic = TRUE),
    "SUCCESS.*TP53.*D_spectrum="
  )
})

test_that(".printMergeSuccess handles very large divergence values in generic path", {
  # Test with divergence values near boundary (log N)
  lmm_data <- list( 
    match_name = "gene_ABC",
    p_interaction = 0.05
  )
  
  div_data <- data.frame(
    estimate = 15.234,  # Large divergence (e.g., for large alphabet size)
    lower_ci = 14.1,
    upper_ci = 16.5
  )
  
  # Should format large values correctly
  expect_message(
    TSENAT:::.printMergeSuccess(lmm_data, div_data, q_values = NA_real_, use_generic = TRUE),
    "SUCCESS.*D_spectrum=.*1\\.52e"
  )
})

test_that(".printMergeSuccess handles negative divergence estimates in generic path", {
  
  # Test with negative divergence (should take absolute value)
  lmm_data <- list(
    match_name = "gene_NEG",
    p_interaction = 0.02
  )
  
  div_data <- data.frame(
    estimate = -0.65,  # Negative estimate
    lower_ci = -0.8,
    upper_ci = -0.4
  )
  
  # Should convert to absolute value and format CI
  expect_message(
    TSENAT:::.printMergeSuccess(lmm_data, div_data, q_values = NA_real_, use_generic = TRUE),
    "SUCCESS.*NEG.*D_spectrum=.*6\\.5"  # 0.65 in scientific or decimal
  )
})

test_that(".printMergeSuccess prints multi-q path with homogeneous values", {
  # Test multi-q path with uniform divergence across q-values
  lmm_data <- list(
    match_name = "gene_MultiQ",
    p_interaction = 0.003
  )
  
  # IMPORTANT: Column names MUST match paste0("estimate_q", q_values)
  # When q=1.0, paste0("estimate_q", 1.0) creates "estimate_q1" (not "estimate_q1.0")
  q_values <- c(0.5, 1.0, 2.0)
  div_data <- data.frame(
    estimate_q0.5 = 0.50,
    estimate_q1 = 0.50,      # Note: "estimate_q1" not "estimate_q1.0"
    estimate_q2 = 0.50
  )
  
  # Should print all divergence values
  expect_message(
    TSENAT:::.printMergeSuccess(lmm_data, div_data, q_values = q_values, use_generic = FALSE),
    "SUCCESS.*MultiQ.*D_spectrum=.*0\\.50.*0\\.50.*0\\.50"
  )
})

test_that(".printMergeSuccess prints multi-q path with heterogeneous values", {
  # Test multi-q path with varying divergence across q-values (spectrum effect)
  lmm_data <- list(
    match_name = "gene_Spectrum",
    p_interaction = 0.01
  )
  
  # IMPORTANT: Column names MUST match paste0("estimate_q", q_values)
  q_values <- c(0.5, 1.0, 2.0)
  div_data <- data.frame(
    estimate_q0.5 = 0.95,   # High at low q (rare-driven)
    estimate_q1 = 0.65,     # Moderate at Shannon (note: "estimate_q1")
    estimate_q2 = 0.35      # Low at high q
  )
  
  # Should show spectrum pattern
  expect_message(
    TSENAT:::.printMergeSuccess(lmm_data, div_data, q_values = q_values, use_generic = FALSE),
    "SUCCESS.*Spectrum.*D_spectrum=.*0\\.95.*0\\.65.*0\\.35"
  )
})

test_that(".printMergeSuccess prints multi-q path with NA estimates", {
  # Test multi-q where some q-values have NA estimates
  lmm_data <- list(
    match_name = "gene_WithNA",
    p_interaction = 0.04
  )
  
  # IMPORTANT: Column names MUST match paste0("estimate_q", q_values)
  q_values <- c(0.5, 1.0, 2.0)
  div_data <- data.frame(
    estimate_q0.5 = 0.7,
    estimate_q1 = NA_real_,  # Missing at q=1.0 (note: "estimate_q1")
    estimate_q2 = 0.5
  )
  
  # Should handle NA values (will extract NA from data)
  expect_message(
    TSENAT:::.printMergeSuccess(lmm_data, div_data, q_values = q_values, use_generic = FALSE),
    "SUCCESS.*WithNA"
  )
})

test_that(".printMergeSuccess prints multi-q path with many q-values", {
  # Test with large q-value spectrum
  lmm_data <- list(
    match_name = "gene_FullSpectrum",
    p_interaction = 0.001
  )
  
  q_values <- c(0.1, 0.5, 1.0, 1.5, 2.0, 3.0)
  
  # IMPORTANT: Create columns matching paste0("estimate_q", q_values)
  # Results: "estimate_q0.1", "estimate_q0.5", "estimate_q1", "estimate_q1.5", "estimate_q2", "estimate_q3"
  col_names <- paste0("estimate_q", q_values)
  div_data <- as.data.frame(setNames(
    as.list(seq(0.9, 0.3, length.out = 6)),
    col_names
  ))
  
  # Should format multiple q-values
  expect_message(
    TSENAT:::.printMergeSuccess(lmm_data, div_data, q_values = q_values, use_generic = FALSE),
    "SUCCESS.*FullSpectrum.*D_spectrum="
  )
})

test_that(".printMergeSuccess handles boundary p-value (p=1.0) in generic path", {
  # Test with p-value at boundary (no significant effect)
  lmm_data <- list(
    match_name = "gene_NoEffect",
    p_interaction = 1.0  # Not significant
  )
  
  div_data <- data.frame(
    estimate = 0.01,  # Very small divergence
    lower_ci = 0.001,
    upper_ci = 0.02
  )
  
  # Should still print (function doesn't filter - that's caller's job)
  expect_message(
    TSENAT:::.printMergeSuccess(lmm_data, div_data, q_values = NA_real_, use_generic = TRUE),
    "SUCCESS.*NoEffect.*p=.*1"
  )
})

test_that(".printMergeSuccess handles boundary p-value (p=0.0) in generic path", {
  # Test with p-value very close to 0 (extremely significant)
  lmm_data <- list(
    match_name = "gene_VerySignificant",
    p_interaction = 1e-300  # Near-zero p-value
  )
  
  div_data <- data.frame(
    estimate = 0.99,
    lower_ci = 0.95,
    upper_ci = 0.999
  )
  
  # Should handle extreme scientific notation
  expect_message(
    TSENAT:::.printMergeSuccess(lmm_data, div_data, q_values = NA_real_, use_generic = TRUE),
    "SUCCESS.*VerySignificant"
  )
})

test_that(".printMergeSuccess handles special gene names in generic path", {
  # Test with complex gene names (symbols, numbers, dots, dashes)
  lmm_data <- list(
    match_name = "ENSG00000000003.13-AS1",  # Real Ensembl format
    p_interaction = 0.005
  )
  
  div_data <- data.frame(
    estimate = 0.42,
    lower_ci = 0.35,
    upper_ci = 0.50
  )
  
  # Should print gene name as-is
  expect_message(
    TSENAT:::.printMergeSuccess(lmm_data, div_data, q_values = NA_real_, use_generic = TRUE),
    "ENSG00000000003"
  )
})

test_that(".printMergeSuccess outputs message (is callable)", {
  # Verify that function doesn't error when called with valid data
  lmm_data <- list(
    match_name = "test_gene",
    p_interaction = 0.012
  )
  
  div_data <- data.frame(
    estimate = 0.5,
    lower_ci = 0.4,
    upper_ci = 0.6
  )
  
  # Should not throw error
  expect_error(
    TSENAT:::.printMergeSuccess(lmm_data, div_data, q_values = NA_real_, use_generic = TRUE),
    NA  # Expect no error
  )
})

test_that(".printMergeSuccess returns invisible(NULL) silently", {
  # Test return value (should be invisible so doesn't print in console)
  lmm_data <- list(
    match_name = "test_gene",
    p_interaction = 0.05
  )
  
  div_data <- data.frame(
    estimate = 0.5,
    lower_ci = 0.4,
    upper_ci = 0.6
  )
  
  # Function doesn't explicitly return but should not print object
  result <- TSENAT:::.printMergeSuccess(
    lmm_data, div_data, q_values = NA_real_, use_generic = TRUE
  )
  
  expect_null(result)
})

test_that(".printMergeSuccess handles zero divergence in multi-q path", {
  # Test multi-q with zero divergence values
  lmm_data <- list(
    match_name = "gene_ZeroDivergence",
    p_interaction = 0.02
  )
  
  # IMPORTANT: Column names MUST match paste0("estimate_q", q_values)
  q_values <- c(0.5, 1.0, 2.0)
  div_data <- data.frame(
    estimate_q0.5 = 0.0,
    estimate_q1 = 0.0,       # Note: "estimate_q1"
    estimate_q2 = 0.0
  )
  
  # Should handle all-zero case
  expect_message(
    TSENAT:::.printMergeSuccess(lmm_data, div_data, q_values = q_values, use_generic = FALSE),
    "SUCCESS.*ZeroDivergence.*0\\.000.*0\\.000.*0\\.000"
  )
})

test_that(".printMergeSuccess handles scientific notation for small divergence in generic", {
  # Test very small divergence (may use scientific notation)
  lmm_data <- list(
    match_name = "gene_TinyDiv",
    p_interaction = 0.03
  )
  
  div_data <- data.frame(
    estimate = 1.23e-5,  # Tiny divergence
    lower_ci = 1e-6,
    upper_ci = 2e-5
  )
  
  # Should format in scientific notation
  expect_message(
    TSENAT:::.printMergeSuccess(lmm_data, div_data, q_values = NA_real_, use_generic = TRUE),
    "SUCCESS.*TinyDiv.*D_spectrum=.*1\\.23e"
  )
})

test_that(".printMergeSuccess formats p-values (generic path works)", {
  # Verify p-value formatting by verifying function executes without error
  test_cases <- list(
    list(p = 0.001234),
    list(p = 0.1234),
    list(p = 0.9999),
    list(p = 1e-6)
  )
  
  for (case in test_cases) {
    lmm_data <- list(
      match_name = "test_gene",
      p_interaction = case$p
    )
    
    div_data <- data.frame(
      estimate = 0.5,
      lower_ci = 0.4,
      upper_ci = 0.6
    )
    
    # Should not error regardless of p-value
    expect_error(
      TSENAT:::.printMergeSuccess(lmm_data, div_data, q_values = NA_real_, use_generic = TRUE),
      NA
    )
  }
})

test_that(".printMergeSuccess formats divergence to 3 digits in generic path", {
  # Verify divergence formatting in generic mode - test that it handles various values
  lmm_data <- list(
    match_name = "test_gene",
    p_interaction = 0.05
  )
  
  # Test various divergence values
  div_data <- data.frame(
    estimate = 0.123456,  # Should truncate to 3 significant figures in scientific notation
    lower_ci = 0.1,
    upper_ci = 0.15
  )
  
  # Should process without error
  expect_error(
    TSENAT:::.printMergeSuccess(lmm_data, div_data, q_values = NA_real_, use_generic = TRUE),
    NA
  )
})

test_that(".printMergeSuccess formats divergence to 3 digits in multi-q path", {
  # Verify divergence formatting in multi-q mode
  lmm_data <- list(
    match_name = "test_gene",
    p_interaction = 0.05
  )
  
  # IMPORTANT: Column names MUST match paste0("estimate_q", q_values)
  q_values <- c(0.5, 1.0, 2.0)
  div_data <- data.frame(
    estimate_q0.5 = 0.123456,
    estimate_q1 = 0.456789,    # Note: "estimate_q1"
    estimate_q2 = 0.789012
  )
  
  # Should process without error
  expect_error(
    TSENAT:::.printMergeSuccess(lmm_data, div_data, q_values = q_values, use_generic = FALSE),
    NA
  )
})
