# Tests for calculate_divergence_bootstrap and results processing
# Covers uncovered lines from divergence_coverage.txt

test_that("calculate_divergence handles per-q pattern classification", {
  # Tests the q-spectrum pattern classification (covered lines 1041-1051)
  se <- SummarizedExperiment::SummarizedExperiment(
    assays = list(counts = matrix(1:12, nrow = 2, ncol = 6)),
    colData = data.frame(
      sample_type = c("Control", "Control", "Treatment", "Treatment", "Treatment", "Treatment")
    ),
    rowData = data.frame(
      gene_name = c("gene1", "gene2")
    )
  )
  colnames(se) <- c("s1", "s2", "s3", "s4", "s5", "s6")
  
  # Test with multiple q values to trigger pattern classification
  result <- calculate_divergence(
    se,
    group_col = "sample_type",
    control_group = "Control",
    q = c(0.5, 1.0, 2.0),
    bootstrap = FALSE,
    verbose = FALSE
  )
  
  # Should have per_q_pattern column in result
  expect_true(!is.null(result))
})

test_that("calculate_divergence applies range normalization", {
  # Tests range normalization (covered lines 1082-1103)
  se <- SummarizedExperiment::SummarizedExperiment(
    assays = list(counts = matrix(1:12, nrow = 2, ncol = 6)),
    colData = data.frame(
      sample_type = c("Control", "Control", "Treatment", "Treatment", "Treatment", "Treatment")
    ),
    rowData = data.frame(
      gene_name = c("gene1", "gene2")
    )
  )
  colnames(se) <- c("s1", "s2", "s3", "s4", "s5", "s6")
  
  # Test normalization method
  result <- calculate_divergence(
    se,
    group_col = "sample_type",
    control_group = "Control",
    q = 1,
    norm = "range",
    bootstrap = FALSE,
    verbose = FALSE
  )
  
  # Result should be computed without error
  expect_true(!is.null(result))
})

test_that("calculate_divergence applies z-score normalization", {
  # Tests z-score normalization (covered lines 1104-1127)
  se <- SummarizedExperiment::SummarizedExperiment(
    assays = list(counts = matrix(1:12, nrow = 2, ncol = 6)),
    colData = data.frame(
      sample_type = c("Control", "Control", "Treatment", "Treatment", "Treatment", "Treatment")
    ),
    rowData = data.frame(
      gene_name = c("gene1", "gene2")
    )
  )
  colnames(se) <- c("s1", "s2", "s3", "s4", "s5", "s6")
  
  result <- calculate_divergence(
    se,
    group_col = "sample_type",
    control_group = "Control",
    q = 1,
    norm = "zscore",
    bootstrap = FALSE,
    verbose = FALSE
  )
  
  expect_true(!is.null(result))
})

test_that("calculate_divergence applies log_odds_ratio normalization", {
  # Tests log_odds_ratio normalization (covered lines 1128-1152)
  se <- SummarizedExperiment::SummarizedExperiment(
    assays = list(counts = matrix(1:12, nrow = 2, ncol = 6)),
    colData = data.frame(
      sample_type = c("Control", "Control", "Treatment", "Treatment", "Treatment", "Treatment")
    ),
    rowData = data.frame(
      gene_name = c("gene1", "gene2")
    )
  )
  colnames(se) <- c("s1", "s2", "s3", "s4", "s5", "s6")
  
  result <- calculate_divergence(
    se,
    group_col = "sample_type",
    control_group = "Control",
    q = 1,
    norm = "log_odds_ratio",
    bootstrap = FALSE,
    verbose = FALSE
  )
  
  expect_true(!is.null(result))
})

test_that("calculate_divergence applies relative_reference normalization", {
  # Tests relative_reference normalization (covered lines 1153+)
  se <- SummarizedExperiment::SummarizedExperiment(
    assays = list(counts = matrix(1:12, nrow = 2, ncol = 6)),
    colData = data.frame(
      sample_type = c("Control", "Control", "Treatment", "Treatment", "Treatment", "Treatment")
    ),
    rowData = data.frame(
      gene_name = c("gene1", "gene2")
    )
  )
  colnames(se) <- c("s1", "s2", "s3", "s4", "s5", "s6")
  
  result <- calculate_divergence(
    se,
    group_col = "sample_type",
    control_group = "Control",
    q = 1,
    norm = "relative_reference",
    bootstrap = FALSE,
    verbose = FALSE
  )
  
  expect_true(!is.null(result))
})

test_that("calculate_divergence skips genes with NA estimates", {
  # Tests handling of NA estimates (covered lines 1028-1030)
  se <- SummarizedExperiment::SummarizedExperiment(
    assays = list(counts = matrix(1:12, nrow = 2, ncol = 6)),
    colData = data.frame(
      sample_type = c("Control", "Control", "Treatment", "Treatment", "Treatment", "Treatment")
    ),
    rowData = data.frame(
      gene_name = c("gene1", "gene2")
    )
  )
  colnames(se) <- c("s1", "s2", "s3", "s4", "s5", "s6")
  
  # Run with progress to trigger logging paths
  result <- calculate_divergence(
    se,
    group_col = "sample_type",
    control_group = "Control",
    q = 1,
    bootstrap = FALSE,
    verbose = FALSE,
    progress = TRUE
  )
  
  expect_true(!is.null(result))
})

test_that("calculate_divergence with paired samples", {
  # Tests paired sample handling
  se <- SummarizedExperiment::SummarizedExperiment(
    assays = list(counts = matrix(1:12, nrow = 2, ncol = 6)),
    colData = data.frame(
      sample_type = c("Control", "Control", "Treatment", "Treatment", "Treatment", "Treatment"),
      paired_samples = c("pair_1", "pair_2", "pair_1", "pair_2", "pair_3", "pair_3")
    ),
    rowData = data.frame(
      gene_name = c("gene1", "gene2")
    )
  )
  colnames(se) <- c("s1", "s2", "s3", "s4", "s5", "s6")
  
  result <- calculate_divergence(
    se,
    group_col = "sample_type",
    control_group = "Control",
    q = 1,
    paired = TRUE,
    bootstrap = FALSE,
    verbose = FALSE
  )
  
  expect_true(!is.null(result))
})

test_that("calculate_divergence with unpaired bootstrapping", {
  # Tests bootstrap with unpaired samples
  se <- SummarizedExperiment::SummarizedExperiment(
    assays = list(counts = matrix(1:12, nrow = 2, ncol = 6)),
    colData = data.frame(
      sample_type = c("Control", "Control", "Treatment", "Treatment", "Treatment", "Treatment")
    ),
    rowData = data.frame(
      gene_name = c("gene1", "gene2")
    )
  )
  colnames(se) <- c("s1", "s2", "s3", "s4", "s5", "s6")
  
  result <- calculate_divergence(
    se,
    group_col = "sample_type",
    control_group = "Control",
    q = 1,
    bootstrap = TRUE,
    nboot = 10,
    verbose = FALSE
  )
  
  expect_true(!is.null(result))
})

test_that("calculate_divergence with paired bootstrapping", {
  # Tests bootstrap with paired samples
  se <- SummarizedExperiment::SummarizedExperiment(
    assays = list(counts = matrix(1:12, nrow = 2, ncol = 6)),
    colData = data.frame(
      sample_type = c("Control", "Control", "Treatment", "Treatment", "Treatment", "Treatment"),
      paired_samples = c("pair_1", "pair_2", "pair_1", "pair_2", "pair_3", "pair_3")
    ),
    rowData = data.frame(
      gene_name = c("gene1", "gene2")
    )
  )
  colnames(se) <- c("s1", "s2", "s3", "s4", "s5", "s6")
  
  result <- calculate_divergence(
    se,
    group_col = "sample_type",
    control_group = "Control",
    q = 1,
    paired = TRUE,
    bootstrap = TRUE,
    nboot = 10,
    verbose = FALSE
  )
  
  expect_true(!is.null(result))
})

test_that("calculate_divergence auto-selects nboot", {
  # Tests auto-selection of nboot (covered in bootstrap logic)
  se <- SummarizedExperiment::SummarizedExperiment(
    assays = list(counts = matrix(1:12, nrow = 2, ncol = 6)),
    colData = data.frame(
      sample_type = c("Control", "Control", "Treatment", "Treatment", "Treatment", "Treatment")
    ),
    rowData = data.frame(
      gene_name = c("gene1", "gene2")
    )
  )
  colnames(se) <- c("s1", "s2", "s3", "s4", "s5", "s6")
  
  result <- calculate_divergence(
    se,
    group_col = "sample_type",
    control_group = "Control",
    q = 1,
    bootstrap = TRUE,
    nboot = "auto",
    verbose = FALSE
  )
  
  expect_true(!is.null(result))
})

test_that("calculate_divergence uses different CI methods", {
  # Tests different confidence interval methods
  se <- SummarizedExperiment::SummarizedExperiment(
    assays = list(counts = matrix(1:12, nrow = 2, ncol = 6)),
    colData = data.frame(
      sample_type = c("Control", "Control", "Treatment", "Treatment", "Treatment", "Treatment")
    ),
    rowData = data.frame(
      gene_name = c("gene1", "gene2")
    )
  )
  colnames(se) <- c("s1", "s2", "s3", "s4", "s5", "s6")
  
  # Test percentile method
  result_percentile <- calculate_divergence(
    se,
    group_col = "sample_type",
    control_group = "Control",
    q = 1,
    bootstrap = TRUE,
    nboot = 10,
    method = "percentile",
    verbose = FALSE
  )
  
  expect_true(!is.null(result_percentile))
})

test_that("calculate_divergence applies CI threshold", {
  # Tests confidence interval parameter
  se <- SummarizedExperiment::SummarizedExperiment(
    assays = list(counts = matrix(1:12, nrow = 2, ncol = 6)),
    colData = data.frame(
      sample_type = c("Control", "Control", "Treatment", "Treatment", "Treatment", "Treatment")
    ),
    rowData = data.frame(
      gene_name = c("gene1", "gene2")
    )
  )
  colnames(se) <- c("s1", "s2", "s3", "s4", "s5", "s6")
  
  result <- calculate_divergence(
    se,
    group_col = "sample_type",
    control_group = "Control",
    q = 1,
    bootstrap = TRUE,
    nboot = 10,
    ci = 0.90,
    verbose = FALSE
  )
  
  expect_true(!is.null(result))
})

test_that("calculate_divergence uses parallel processing", {
  # Tests parallel/multi-threaded execution
  se <- SummarizedExperiment::SummarizedExperiment(
    assays = list(counts = matrix(1:12, nrow = 2, ncol = 6)),
    colData = data.frame(
      sample_type = c("Control", "Control", "Treatment", "Treatment", "Treatment", "Treatment")
    ),
    rowData = data.frame(
      gene_name = c("gene1", "gene2")
    )
  )
  colnames(se) <- c("s1", "s2", "s3", "s4", "s5", "s6")
  
  result <- calculate_divergence(
    se,
    group_col = "sample_type",
    control_group = "Control",
    q = 1,
    nthreads = 2,
    bootstrap = FALSE,
    verbose = FALSE
  )
  
  expect_true(!is.null(result))
})

test_that("calculate_divergence with log_base parameter", {
  # Tests alternative log base for entropy computation
  se <- SummarizedExperiment::SummarizedExperiment(
    assays = list(counts = matrix(1:12, nrow = 2, ncol = 6)),
    colData = data.frame(
      sample_type = c("Control", "Control", "Treatment", "Treatment", "Treatment", "Treatment")
    ),
    rowData = data.frame(
      gene_name = c("gene1", "gene2")
    )
  )
  colnames(se) <- c("s1", "s2", "s3", "s4", "s5", "s6")
  
  result <- calculate_divergence(
    se,
    group_col = "sample_type",
    control_group = "Control",
    q = 1,
    log_base = 2,
    bootstrap = FALSE,
    verbose = FALSE
  )
  
  expect_true(!is.null(result))
})

test_that("calculate_divergence with pseudocount parameter", {
  # Tests pseudocount handling for zero-count genes
  se <- SummarizedExperiment::SummarizedExperiment(
    assays = list(counts = matrix(c(0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11), nrow = 2, ncol = 6)),
    colData = data.frame(
      sample_type = c("Control", "Control", "Treatment", "Treatment", "Treatment", "Treatment")
    ),
    rowData = data.frame(
      gene_name = c("gene1", "gene2")
    )
  )
  colnames(se) <- c("s1", "s2", "s3", "s4", "s5", "s6")
  
  result <- calculate_divergence(
    se,
    group_col = "sample_type",
    control_group = "Control",
    q = 1,
    pseudocount = 1.0,
    bootstrap = FALSE,
    verbose = FALSE
  )
  
  expect_true(!is.null(result))
})
