# ============================================================================
# Integration tests derived from the S4-layer audit (s4_scripts.md, section 17)
#
# These tests exercise the public S4 API end-to-end so that the critical
# config-resolution and dimensional-semantics bugs stay fixed:
#   1. config$paired / config$bootstrap must reach the divergence core
#   2. calculate_assumptions(q = 1) must analyze EXACTLY q = 1
#   3. calculate_divergence(q = 0) must keep q = 0 (no silent replacement)
#   4. `[` subsetting must keep divergence q-value columns intact
#   5. Subsetting must expose an explicit staleness mechanism
# ============================================================================

# Self-contained factory: multi-isoform toy data with paired structure
.make_integration_analysis <- function(n_genes = 8, n_samples = 8, seed = 42) {
  set.seed(seed)
  n_tx <- n_genes * 10L
  counts <- pmax(matrix(rpois(n_tx * n_samples, lambda = 40), nrow = n_tx,
                        ncol = n_samples), 50)
  tx_names <- paste0("tx", seq_len(n_tx))
  sample_names <- paste0("s", seq_len(n_samples))
  rownames(counts) <- tx_names
  colnames(counts) <- sample_names

  gene_ids <- paste0("g", rep(seq_len(n_genes), each = n_tx / n_genes))
  row_data <- S4Vectors::DataFrame(
    transcript_id = tx_names,
    gene_id = gene_ids,
    row.names = tx_names
  )
  col_data <- S4Vectors::DataFrame(
    condition = rep(c("control", "treatment"), each = n_samples / 2),
    subject = rep(seq_len(n_samples / 2), 2),
    paired_samples = paste0("pair", rep(seq_len(n_samples / 2), 2)),
    row.names = sample_names
  )

  se <- SummarizedExperiment::SummarizedExperiment(
    assays = list(counts = counts),
    rowData = row_data,
    colData = col_data
  )
  S4Vectors::metadata(se)$tx2gene <- data.frame(
    Transcript = tx_names, Gene = gene_ids, stringsAsFactors = FALSE
  )

  TSENAT::TSENATAnalysis(se = se, config = list(
    condition_col = "condition",
    control_group = "control",
    nthreads = 1,
    verbose = FALSE
  ))
}

test_that("Integration: config$paired and config$bootstrap reach the divergence core", {
  analysis <- .make_integration_analysis()
  analysis <- calculate_diversity(analysis, q = c(0.5, 1, 1.5), min_valid_frac = 0)

  # Set via config (NOT via explicit arguments)
  analysis@config$paired <- TRUE
  analysis@config$subject_col <- "paired_samples"
  analysis@config$bootstrap <- TRUE

  analysis <- calculate_divergence(analysis, q = c(0.5, 1), nboot = 25)

  expect_true(isTRUE(analysis@metadata$last_divergence$paired))
  expect_true(isTRUE(analysis@metadata$last_divergence$bootstrap))

  # The bootstrap config must have been propagated into the result metadata
  div_se <- analysis@divergence_results$divergence_se
  expect_equal(div_se@metadata$bootstrap_config$nboot, 25)
})

test_that("Integration: calculate_assumptions uses exactly the requested q", {
  analysis <- .make_integration_analysis()
  analysis <- calculate_diversity(analysis, q = c(0.5, 1, 1.5), min_valid_frac = 0)

  analysis <- calculate_assumptions(analysis, q = 1)

  expect_equal(analysis@metadata$rankbased_assumptions$q_value_tested, 1)

  # Requesting an unavailable q must be an error, not a silent fallback
  expect_error(
    calculate_assumptions(analysis, q = 0.7),
    "not found in diversity results"
  )
})

test_that("Integration: calculate_divergence preserves q = 0", {
  analysis <- .make_integration_analysis()
  analysis <- calculate_diversity(analysis, q = c(0.5, 1), min_valid_frac = 0)

  analysis <- calculate_divergence(analysis, q = 0)

  div_se <- analysis@divergence_results$divergence_se
  q_cols <- SummarizedExperiment::colData(div_se)$q_value
  expect_true(any(q_cols == 0))
  expect_true(isTRUE(analysis@metadata$last_divergence$q == 0))
})

test_that("Integration: subsetting keeps divergence q-value columns intact", {
  analysis <- .make_integration_analysis()
  analysis <- calculate_diversity(analysis, q = c(0.5, 1, 1.5), min_valid_frac = 0)
  analysis <- calculate_divergence(analysis, q = c(0.5, 1, 1.5))

  n_q <- ncol(analysis@divergence_results$divergence_se)

  sub_obj <- analysis[1:6, 1:2]

  expect_equal(ncol(sub_obj@se), 2)
  # Divergence columns are q-VALUES, not samples: all must be preserved
  expect_equal(ncol(sub_obj@divergence_results$divergence_se), n_q)
  expect_equal(nrow(sub_obj@divergence_results$divergence_se), 6)
})

test_that("Integration: subsetting exposes an explicit staleness mechanism", {
  analysis <- .make_integration_analysis()
  analysis <- calculate_diversity(analysis, q = c(0.5, 1), min_valid_frac = 0)

  sub_obj <- analysis[1:6, ]

  expect_true(isTRUE(sub_obj@metadata$subset_applied))
  expect_equal(sub_obj@metadata$subset_info$genes, 6)
  expect_true(all(c("sait_results", "jackknife_results", "rank_test_results",
                    "divergence_results") %in% sub_obj@metadata$stale_results))
})
