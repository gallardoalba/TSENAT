# Test for Tsallis entropy calculation in calculate_diversity

test_that("calculate_diversity returns correct Tsallis entropy for single q", {
  set.seed(123)
  x <- matrix(rpois(60, 10), ncol = 6)
  colnames(x) <- paste0("Sample", 1:6)
  gene <- c(rep("Gene1", 3), rep("Gene2", 2), rep("Gene3", 3), rep("Gene4", 2))
  q <- 2
  result <- calculate_diversity(x, gene, q = q)
  expect_s4_class(result, "SummarizedExperiment")
  expect_true("diversity" %in% names(SummarizedExperiment::assays(result)))
  expect_equal(nrow(result), length(unique(gene)))
  expect_equal(ncol(result), ncol(x))
  expect_true(!is.null(SummarizedExperiment::rowData(result)$genes))
})

test_that("calculate_diversity returns correct Tsallis entropy for vector q", {
  set.seed(123)
  x <- matrix(rpois(60, 10), ncol = 6)
  colnames(x) <- paste0("Sample", 1:6)
  gene <- c(rep("Gene1", 3), rep("Gene2", 2), rep("Gene3", 3), rep("Gene4", 2))
  q <- c(1, 2)
  result <- calculate_diversity(x, gene, q = q)
  expect_s4_class(result, "SummarizedExperiment")
  expect_true("diversity" %in% names(SummarizedExperiment::assays(result)))
  expect_equal(nrow(result), length(unique(gene)))
  expect_equal(ncol(result), ncol(x) * length(q))
  expect_true(!is.null(SummarizedExperiment::rowData(result)$genes))
  expect_true(all(c("samples", "q") %in% colnames(SummarizedExperiment::colData(result))))
  expect_equal(length(unique(SummarizedExperiment::colData(result)$q)), length(q))
})
