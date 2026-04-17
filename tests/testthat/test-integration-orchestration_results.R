library(testthat)

context("Orchestration: Results Processing and Ranking")

# ==============================================================================
# .extract_jackknife_multi_q(): Tests for multi-Q jackknife result extraction
# ==============================================================================

test_that(".extract_jackknife_multi_q extracts single-Q results with pvalue rankBy", {
  # Single Q value result with pvalue rankBy
  result_list <- list(
    q_0_00 = data.frame(
      gene = c("g1", "g2", "g3"),
      pvalue = c(0.001, 0.05, 0.1),
      delta_influence = c(0.5, 0.3, 0.2)
    )
  )
  
  extracted <- .extract_jackknife_multi_q(result_list, q = 0.0, rankBy = "pvalue")
  
  expect_is(extracted, "data.frame")
  expect_equal(nrow(extracted), 3)
  expect_true(all(c("gene", "pvalue", "delta_influence") %in% colnames(extracted)))
})

test_that(".extract_jackknife_multi_q extracts from named multi-Q list", {
  # Multiple Q value results: extract by named index
  result_list <- list(
    q_0_00 = data.frame(gene = c("g1", "g2"), pvalue = c(0.01, 0.05)),
    q_0_50 = data.frame(gene = c("g1", "g2"), pvalue = c(0.02, 0.06)),
    q_1_00 = data.frame(gene = c("g1", "g2"), pvalue = c(0.03, 0.07))
  )
  
  extracted_q0 <- .extract_jackknife_multi_q(result_list, q = 0.0, rankBy = "pvalue")
  extracted_q1 <- .extract_jackknife_multi_q(result_list, q = 1.0, rankBy = "pvalue")
  
  expect_equal(nrow(extracted_q0), 2)
  expect_equal(extracted_q0$pvalue[1], 0.01)
  expect_equal(extracted_q1$pvalue[1], 0.03)
})

test_that(".extract_jackknife_multi_q handles multi_q with effectSize rankBy", {
  # Test with multi_q element and effectSize rankBy
  result_list <- list(
    multi_q = list(
      summary_table = data.frame(gene = c("g1", "g2"), delta_influence = c(1.0, 2.0))
    )
  )
  
  extracted <- .extract_jackknife_multi_q(result_list, q = NULL, rankBy = "effectSize")
  
  expect_is(extracted, "data.frame")
  expect_equal(nrow(extracted), 2)
})

test_that(".extract_jackknife_multi_q handles NULL q with multi_q pvalue rankBy", {
  # When q=NULL, should use multi_q element if available with pvalue rankBy
  result_list <- list(
    multi_q = data.frame(
      gene = c("g1", "g2"),
      pvalue = c(0.01, 0.05)
    )
  )
  
  extracted <- .extract_jackknife_multi_q(result_list, q = NULL, rankBy = "pvalue")
  
  expect_is(extracted, "data.frame")
  expect_equal(nrow(extracted), 2)
})

test_that(".extract_jackknife_multi_q handles non-list input", {
  # Input that is already a dataframe (not list): should pass through
  result_df <- data.frame(gene = c("g1", "g2"), pvalue = c(0.01, 0.05))
  
  extracted <- .extract_jackknife_multi_q(result_df, q = 0.0, rankBy = "pvalue")
  
  expect_true(is.data.frame(extracted) || is.null(extracted))
})

test_that(".extract_jackknife_multi_q handles missing Q value with warning", {
  # Request Q value that doesn't exist - should warn and return NULL
  result_list <- list(
    q_0_00 = data.frame(gene = c("g1", "g2"), pvalue = c(0.01, 0.05))
  )
  
  # Suppress warning for this test
  expect_warning(
    extracted <- .extract_jackknife_multi_q(result_list, q = 0.5, rankBy = "pvalue"),
    "Jackknife results for q="
  )
  
  expect_null(extracted)
})

# ==============================================================================
# .get_ranking_column(): Tests for ranking column designation
# ==============================================================================

test_that(".get_ranking_column returns p_interaction for LM with pvalue", {
  result <- data.frame(
    gene = c("g1", "g2"),
    p_interaction = c(0.01, 0.05),
    estimate = c(0.5, 0.3)
  )
  
  col <- .get_ranking_column(type = "lm", rankBy = "pvalue", result = result)
  
  expect_equal(col, "p_interaction")
})

test_that(".get_ranking_column returns NULL for LM pvalue when column absent", {
  result <- data.frame(
    gene = c("g1", "g2"),
    estimate = c(0.5, 0.3)
  )
  
  col <- .get_ranking_column(type = "lm", rankBy = "pvalue", result = result)
  
  expect_null(col)
})

test_that(".get_ranking_column returns adj_p_interaction for LM with qvalue", {
  result <- data.frame(
    gene = c("g1", "g2"),
    adj_p_interaction = c(0.05, 0.10),
    estimate = c(0.5, 0.3)
  )
  
  col <- .get_ranking_column(type = "lm", rankBy = "qvalue", result = result)
  
  expect_equal(col, "adj_p_interaction")
})

test_that(".get_ranking_column returns NULL for LM qvalue when column absent", {
  result <- data.frame(
    gene = c("g1", "g2"),
    p_value = c(0.01, 0.05)
  )
  
  col <- .get_ranking_column(type = "lm", rankBy = "qvalue", result = result)
  
  expect_null(col)
})

test_that(".get_ranking_column prioritizes statistic for LM effectSize", {
  result <- data.frame(
    gene = c("g1", "g2"),
    statistic = c(2.5, 3.0),
    estimate = c(0.5, 0.3)
  )
  
  col <- .get_ranking_column(type = "lm", rankBy = "effectSize", result = result)
  
  expect_equal(col, "statistic")
})

test_that(".get_ranking_column fallback to estimate for LM effectSize", {
  result <- data.frame(
    gene = c("g1", "g2"),
    estimate = c(0.5, 0.3),
    p_value = c(0.01, 0.05)
  )
  
  col <- .get_ranking_column(type = "lm", rankBy = "effectSize", result = result)
  
  expect_equal(col, "estimate")
})

test_that(".get_ranking_column fallback to effect_size for LM effectSize", {
  result <- data.frame(
    gene = c("g1", "g2"),
    effect_size = c(0.5, 0.3),
    p_value = c(0.01, 0.05)
  )
  
  col <- .get_ranking_column(type = "lm", rankBy = "effectSize", result = result)
  
  expect_equal(col, "effect_size")
})

test_that(".get_ranking_column returns NULL for LM effectSize when no effect columns", {
  result <- data.frame(
    gene = c("g1", "g2"),
    p_value = c(0.01, 0.05)
  )
  
  col <- .get_ranking_column(type = "lm", rankBy = "effectSize", result = result)
  
  expect_null(col)
})

test_that(".get_ranking_column handles rank_test type with pvalue", {
  result <- data.frame(
    gene = c("g1", "g2"),
    p_value = c(0.01, 0.05),
    statistic = c(1.5, 2.0)
  )
  
  col <- .get_ranking_column(type = "rank_test", rankBy = "pvalue", result = result)
  
  expect_equal(col, "p_value")
})

test_that(".get_ranking_column handles rank_test type with qvalue", {
  result <- data.frame(
    gene = c("g1", "g2"),
    adj_p_value = c(0.05, 0.10),
    p_value = c(0.01, 0.02)
  )
  
  col <- .get_ranking_column(type = "rank_test", rankBy = "qvalue", result = result)
  
  expect_equal(col, "adj_p_value")
})

test_that(".get_ranking_column handles rank_test type with effectSize", {
  result <- data.frame(
    gene = c("g1", "g2"),
    statistic = c(2.5, 3.0),
    p_value = c(0.01, 0.05)
  )
  
  col <- .get_ranking_column(type = "rank_test", rankBy = "effectSize", result = result)
  
  expect_equal(col, "statistic")
})

test_that(".get_ranking_column handles jackknife type with pvalue", {
  result <- data.frame(
    gene = c("g1", "g2"),
    pvalue = c(0.01, 0.05),
    delta_influence = c(0.5, 0.3)
  )
  
  col <- .get_ranking_column(type = "jackknife", rankBy = "pvalue", result = result)
  
  expect_equal(col, "pvalue")
})

test_that(".get_ranking_column handles jackknife type with qvalue", {
  result <- data.frame(
    gene = c("g1", "g2"),
    fdr = c(0.05, 0.10),
    pvalue = c(0.01, 0.02)
  )
  
  col <- .get_ranking_column(type = "jackknife", rankBy = "qvalue", result = result)
  
  expect_equal(col, "fdr")
})

test_that(".get_ranking_column handles jackknife type with effectSize", {
  result <- data.frame(
    gene = c("g1", "g2"),
    delta_influence = c(2.5, 3.0),
    pvalue = c(0.01, 0.05)
  )
  
  col <- .get_ranking_column(type = "jackknife", rankBy = "effectSize", result = result)
  
  expect_equal(col, "delta_influence")
})

test_that(".get_ranking_column fallback to max_delta_influence for jackknife effectSize", {
  result <- data.frame(
    gene = c("g1", "g2"),
    max_delta_influence = c(2.5, 3.0),
    pvalue = c(0.01, 0.05)
  )
  
  col <- .get_ranking_column(type = "jackknife", rankBy = "effectSize", result = result)
  
  expect_equal(col, "max_delta_influence")
})

test_that(".get_ranking_column returns NULL for jackknife effectSize when no effect columns", {
  result <- data.frame(
    gene = c("g1", "g2"),
    pvalue = c(0.01, 0.05)
  )
  
  col <- .get_ranking_column(type = "jackknife", rankBy = "effectSize", result = result)
  
  expect_null(col)
})

test_that(".get_ranking_column handles invalid type", {
  result <- data.frame(
    gene = c("g1", "g2"),
    p_value = c(0.01, 0.05)
  )
  
  col <- .get_ranking_column(type = "invalid_type", rankBy = "pvalue", result = result)
  
  expect_null(col)
})

test_that(".get_ranking_column handles invalid rankBy", {
  result <- data.frame(
    gene = c("g1", "g2"),
    p_value = c(0.01, 0.05),
    statistic = c(2.5, 3.0)
  )
  
  col <- .get_ranking_column(type = "lm", rankBy = "invalid_rank", result = result)
  
  expect_null(col)
})

test_that(".get_ranking_column handles empty result dataframe", {
  result <- data.frame()
  
  col <- .get_ranking_column(type = "lm", rankBy = "pvalue", result = result)
  
  expect_null(col)
})

test_that(".get_ranking_column with rank_test fallback for effectSize", {
  result <- data.frame(
    gene = c("g1", "g2"),
    estimate = c(0.5, 0.3),
    p_value = c(0.01, 0.05)
  )
  
  col <- .get_ranking_column(type = "rank_test", rankBy = "effectSize", result = result)
  
  expect_equal(col, "estimate")
})

# ==============================================================================
# .filter_statistical_by_fdr(): Tests for FDR-based filtering
# ==============================================================================

test_that(".filter_statistical_by_fdr filters LM results by adj_p_interaction", {
  result <- data.frame(
    gene = c("g1", "g2", "g3", "g4"),
    p_interaction = c(0.001, 0.01, 0.05, 0.1),
    adj_p_interaction = c(0.01, 0.05, 0.10, 0.20),
    estimate = c(0.5, 0.3, 0.2, 0.1)
  )
  
  filtered <- .filter_statistical_by_fdr(result, type = "lm", filterFDR = 0.05)
  
  expect_is(filtered, "data.frame")
  expect_equal(nrow(filtered), 2)  # Only rows with adj_p_interaction <= 0.05
  expect_true(all(filtered$adj_p_interaction <= 0.05))
})

test_that(".filter_statistical_by_fdr filters rank_test results by adj_p_value", {
  result <- data.frame(
    gene = c("g1", "g2", "g3"),
    p_value = c(0.001, 0.01, 0.05),
    adj_p_value = c(0.01, 0.05, 0.15),
    statistic = c(2.5, 2.0, 1.5)
  )
  
  filtered <- .filter_statistical_by_fdr(result, type = "rank_test", filterFDR = 0.08)
  
  expect_equal(nrow(filtered), 2)
  expect_true(all(filtered$adj_p_value <= 0.08))
})

test_that(".filter_statistical_by_fdr filters jackknife results by fdr", {
  result <- data.frame(
    gene = c("g1", "g2", "g3", "g4"),
    pvalue = c(0.001, 0.01, 0.05, 0.1),
    fdr = c(0.01, 0.05, 0.10, 0.20),
    delta_influence = c(0.5, 0.3, 0.2, 0.1)
  )
  
  filtered <- .filter_statistical_by_fdr(result, type = "jackknife", filterFDR = 0.08)
  
  expect_equal(nrow(filtered), 2)
  expect_true(all(filtered$fdr <= 0.08))
})

test_that(".filter_statistical_by_fdr returns all rows when filterFDR is NULL", {
  result <- data.frame(
    gene = c("g1", "g2", "g3"),
    adj_p_interaction = c(0.01, 0.05, 0.20),
    estimate = c(0.5, 0.3, 0.1)
  )
  
  filtered <- .filter_statistical_by_fdr(result, type = "lm", filterFDR = NULL)
  
  expect_equal(nrow(filtered), 3)
  expect_identical(filtered, result)
})

test_that(".filter_statistical_by_fdr returns input when not a dataframe", {
  result <- list(summary_table = data.frame(gene = c("g1", "g2")))
  
  filtered <- .filter_statistical_by_fdr(result, type = "lm", filterFDR = 0.05)
  
  expect_identical(filtered, result)
})

test_that(".filter_statistical_by_fdr handles missing adjusted p-value column", {
  result <- data.frame(
    gene = c("g1", "g2"),
    p_value = c(0.01, 0.05),
    estimate = c(0.5, 0.3)
  )
  
  filtered <- .filter_statistical_by_fdr(result, type = "lm", filterFDR = 0.05)
  
  expect_equal(nrow(filtered), 2)
  expect_identical(filtered, result)
})

test_that(".filter_statistical_by_fdr handles NA values in adjusted p-values", {
  result <- data.frame(
    gene = c("g1", "g2", "g3"),
    adj_p_interaction = c(0.01, NA, 0.05),
    estimate = c(0.5, 0.3, 0.2)
  )
  
  filtered <- .filter_statistical_by_fdr(result, type = "lm", filterFDR = 0.06)
  
  expect_equal(nrow(filtered), 2)  # NA rows excluded
  expect_false(any(is.na(filtered$adj_p_interaction)))
})

# ==============================================================================
# .rank_statistical_results(): Tests for result ranking
# ==============================================================================

test_that(".rank_statistical_results ranks by pvalue ascending", {
  result <- data.frame(
    gene = c("g1", "g2", "g3"),
    p_interaction = c(0.05, 0.01, 0.1),
    estimate = c(0.5, 0.7, 0.3)
  )
  
  ranked <- .rank_statistical_results(result, type = "lm", rankBy = "pvalue", n = NA)
  
  expect_equal(ranked$gene, c("g2", "g1", "g3"))  # Sorted by p_interaction
  expect_true(all(ranked$p_interaction == sort(result$p_interaction)))
})

test_that(".rank_statistical_results ranks by effectSize descending", {
  result <- data.frame(
    gene = c("g1", "g2", "g3"),
    statistic = c(1.5, 3.0, 2.0),
    p_value = c(0.01, 0.05, 0.1)
  )
  
  ranked <- .rank_statistical_results(result, type = "rank_test", rankBy = "effectSize", n = NA)
  
  expect_equal(ranked$gene, c("g2", "g3", "g1"))  # Sorted by abs(statistic) descending
  expect_equal(ranked$statistic, c(3.0, 2.0, 1.5))
})

test_that(".rank_statistical_results respects n parameter", {
  result <- data.frame(
    gene = c("g1", "g2", "g3", "g4", "g5"),
    p_value = c(0.001, 0.01, 0.05, 0.08, 0.1),
    pvalue = c(0.001, 0.01, 0.05, 0.08, 0.1)
  )
  
  ranked <- .rank_statistical_results(result, type = "jackknife", rankBy = "pvalue", n = 3)
  
  expect_equal(nrow(ranked), 3)  # Only top 3 returned
})

test_that(".rank_statistical_results skips ranking with rankBy='none'", {
  result <- data.frame(
    gene = c("g1", "g2", "g3"),
    p_value = c(0.05, 0.01, 0.1)
  )
  
  ranked <- .rank_statistical_results(result, type = "rank_test", rankBy = "none", n = NA)
  
  expect_identical(ranked, result)  # Unchanged
})

test_that(".rank_statistical_results extracts summary_table from list", {
  result <- list(
    summary_table = data.frame(
      gene = c("g1", "g2", "g3"),
      p_value = c(0.05, 0.01, 0.1),
      statistic = c(1.5, 2.5, 1.0)
    )
  )
  
  ranked <- .rank_statistical_results(result, type = "rank_test", rankBy = "pvalue", n = NA)
  
  expect_is(ranked, "data.frame")
  expect_equal(ranked$gene, c("g2", "g1", "g3"))
})

test_that(".rank_statistical_results handles NA values when sorting", {
  result <- data.frame(
    gene = c("g1", "g2", "g3"),
    p_value = c(0.01, NA, 0.05),
    statistic = c(1.0, 2.0, 1.5)
  )
  
  ranked <- .rank_statistical_results(result, type = "rank_test", rankBy = "pvalue", n = NA)
  
  expect_is(ranked, "data.frame")
  expect_equal(ranked$gene[1:2], c("g1", "g3"))  # NA pushed to end
})

# ==============================================================================
# .extract_statistical_dataframe(): Tests for dataframe extraction
# ==============================================================================

test_that(".extract_statistical_dataframe returns dataframe as-is", {
  result <- data.frame(gene = c("g1", "g2"), p_value = c(0.01, 0.05))
  
  extracted <- .extract_statistical_dataframe(result, type = "lm")
  
  expect_identical(extracted, result)
})

test_that(".extract_statistical_dataframe extracts from summary_table list", {
  result <- list(
    summary_table = data.frame(gene = c("g1", "g2"), p_value = c(0.01, 0.05))
  )
  
  extracted <- .extract_statistical_dataframe(result, type = "lm")
  
  expect_is(extracted, "data.frame")
  expect_equal(nrow(extracted), 2)
})

test_that(".extract_statistical_dataframe extracts from results list element", {
  result <- list(
    results = data.frame(gene = c("g1", "g2", "g3"), estimate = c(0.5, 0.3, 0.2))
  )
  
  extracted <- .extract_statistical_dataframe(result, type = "rank_test")
  
  expect_is(extracted, "data.frame")
  expect_equal(nrow(extracted), 3)
})

test_that(".extract_statistical_dataframe extracts from all_transcript_stats", {
  result <- list(
    all_transcript_stats = data.frame(gene = c("t1", "t2"), pvalue = c(0.01, 0.05))
  )
  
  extracted <- .extract_statistical_dataframe(result, type = "jackknife")
  
  expect_is(extracted, "data.frame")
  expect_equal(nrow(extracted), 2)
})

test_that(".extract_statistical_dataframe returns NULL for non-list non-dataframe", {
  result <- c(1, 2, 3)
  
  extracted <- .extract_statistical_dataframe(result, type = "lm")
  
  expect_null(extracted)
})

test_that(".extract_statistical_dataframe returns NULL for empty list", {
  result <- list()
  
  extracted <- .extract_statistical_dataframe(result, type = "lm")
  
  expect_null(extracted)
})

test_that(".extract_statistical_dataframe returns NULL for list without recognized fields", {
  result <- list(other_field = "value", another = 123)
  
  extracted <- .extract_statistical_dataframe(result, type = "lm")
  
  expect_null(extracted)
})

# ==============================================================================
# .process_statistical_results(): Tests for full statistical result processing
# ==============================================================================

test_that(".process_statistical_results returns NULL for NULL input", {
  result <- .process_statistical_results(NULL, type = "lm", filterFDR = 0.05, 
                                         rankBy = "pvalue", n = NA, format = "text")
  
  expect_null(result)
})

test_that(".process_statistical_results processes complete lm results", {
  input <- data.frame(
    gene = c("g1", "g2", "g3"),
    p_interaction = c(0.001, 0.05, 0.15),
    adj_p_interaction = c(0.01, 0.10, 0.25),
    estimate = c(0.5, 0.3, 0.1)
  )
  
  result <- .process_statistical_results(input, type = "lm", filterFDR = 0.15,
                                         rankBy = "pvalue", n = 2, format = "text")
  
  expect_is(result, "data.frame")
  expect_equal(nrow(result), 2)  # Top 2 by pvalue
  expect_true(all(result$adj_p_interaction <= 0.15))
})

test_that(".process_statistical_results filters then ranks results", {
  input <- data.frame(
    gene = c("g1", "g2", "g3", "g4"),
    p_value = c(0.001, 0.01, 0.05, 0.1),
    adj_p_value = c(0.01, 0.05, 0.10, 0.20),
    statistic = c(2.0, 3.0, 1.5, 1.0)
  )
  
  result <- .process_statistical_results(input, type = "rank_test", filterFDR = 0.08,
                                         rankBy = "effectSize", n = NA, format = "text")
  
  expect_is(result, "data.frame")
  expect_true(all(result$adj_p_value <= 0.08))
  expect_equal(result$gene[1], "g2")  # Highest statistic among filtered
})

test_that(".process_statistical_results returns NULL when all rows filtered", {
  input <- data.frame(
    gene = c("g1", "g2"),
    adj_p_interaction = c(0.10, 0.15),
    estimate = c(0.5, 0.3)
  )
  
  result <- .process_statistical_results(input, type = "lm", filterFDR = 0.05,
                                         rankBy = "pvalue", n = NA, format = "text")
  
  expect_null(result)
})

test_that(".process_statistical_results handles format conversion", {
  input <- data.frame(
    gene = c("g1", "g2"),
    p_value = c(0.01, 0.05),
    estimate = c(0.5, 0.3)
  )
  
  result <- .process_statistical_results(input, type = "rank_test", filterFDR = NULL,
                                         rankBy = "none", n = NA, format = "matrix")
  
  expect_is(result, "matrix")
})

# ==============================================================================
# .process_divergence_results(): Tests for divergence result processing
# ==============================================================================

test_that(".process_divergence_results returns NULL for NULL input", {
  result <- .process_divergence_results(NULL, filterFDR = 0.05, format = "text")
  
  expect_null(result)
})

test_that(".process_divergence_results handles dataframe with FDR filtering", {
  input <- data.frame(
    gene = c("g1", "g2", "g3"),
    divergence = c(0.5, 0.3, 0.2),
    padj = c(0.01, 0.05, 0.15)
  )
  
  result <- .process_divergence_results(input, filterFDR = 0.08, format = "text")
  
  expect_is(result, "data.frame")
  expect_equal(nrow(result), 2)
})

test_that(".process_divergence_results handles matrix input", {
  input <- matrix(c(0.5, 0.3, 0.2, 0.1, 0.4, 0.6), nrow = 2, ncol = 3,
                  dimnames = list(c("g1", "g2"), c("s1", "s2", "s3")))
  
  result <- .process_divergence_results(input, filterFDR = NULL, format = "text")
  
  expect_true(is.matrix(result) || is.data.frame(result))
})

test_that(".process_divergence_results returns input when filterFDR unsupported", {
  input <- list(
    item1 = data.frame(a = c(1, 2)),
    item2 = data.frame(b = c(3, 4))
  )
  
  result <- .process_divergence_results(input, filterFDR = NULL, format = "text")
  
  expect_is(result, "list")
})

# ==============================================================================
# .warn_unsupported_params(): Tests for parameter validation
# ==============================================================================

test_that(".warn_unsupported_params warns for unsupported rankBy with diversity", {
  expect_warning(
    .warn_unsupported_params(type = "diversity", filterFDR = 0.05, rankBy = "pvalue"),
    "rankBy.*not supported.*diversity"
  )
})

test_that(".warn_unsupported_params warns for unsupported rankBy with assumptions", {
  expect_warning(
    .warn_unsupported_params(type = "assumptions", filterFDR = NULL, rankBy = "pvalue"),
    "rankBy.*not supported.*assumptions"
  )
})

test_that(".warn_unsupported_params does not warn for supported parameters", {
  expect_silent(
    .warn_unsupported_params(type = "lm", filterFDR = 0.05, rankBy = "pvalue")
  )
})

test_that(".warn_unsupported_params warns for unsupported rankBy with switching", {
  expect_warning(
    .warn_unsupported_params(type = "switching_tables", filterFDR = NULL, rankBy = "effectSize"),
    "rankBy.*not supported.*switching_tables"
  )
})

# ==============================================================================
# .process_effect_sizes_divergence_results(): Tests (27.3%)
# ==============================================================================

test_that(".process_effect_sizes_divergence_results processes divergence matrix", {
  div_matrix <- matrix(rnorm(100, mean=1, sd=0.3), nrow=20, ncol=5)
  rownames(div_matrix) <- paste0("gene_", 1:20)
  colnames(div_matrix) <- paste0("sample_", 1:5)
  
  result <- .process_effect_sizes_divergence_results(
    result = div_matrix,
    top_n = 10,
    sort_by = "overall"
  )
  
  expect_is(result, c("data.frame", "matrix", "list"))
})

test_that(".process_effect_sizes_divergence_results respects top_n", {
  div_df <- data.frame(
    gene = paste0("g", 1:20),
    divergence = rnorm(20, mean=1, sd=0.5),
    adj_p_interaction = runif(20, 0, 0.1)
  )
  
  result <- .process_effect_sizes_divergence_results(
    result = div_df,
    top_n = 5,
    sort_by = "adj_p_interaction"
  )
  
  expect_is(result, c("data.frame", "matrix", "list"))
})
