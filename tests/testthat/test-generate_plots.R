context("Generate plots helpers")

test_that("plot_diversity_density returns ggplot object", {
  skip_if_not_installed("SummarizedExperiment")
  skip_if_not_installed("ggplot2")
  skip_if_not_installed("tidyr")
  skip_if_not_installed("dplyr")

  library(SummarizedExperiment)
  library(ggplot2)
  library(tidyr)
  library(dplyr)

  # construct minimal SummarizedExperiment
  mat <- matrix(runif(20), nrow = 5, ncol = 4)
  colnames(mat) <- c("S1_N", "S2_T", "S3_N", "S4_T")
  rownames(mat) <- paste0("G", 1:5)
  se <- SummarizedExperiment::SummarizedExperiment(assays = list(diversity = mat), rowData = DataFrame(genes = rownames(mat)), colData = DataFrame(samples = colnames(mat)))

  p <- plot_diversity_density(se)
  expect_s3_class(p, "gg")
})


test_that("plot_mean_violin returns ggplot object", {
  skip_if_not_installed("SummarizedExperiment")
  skip_if_not_installed("ggplot2")
  skip_if_not_installed("tidyr")
  skip_if_not_installed("dplyr")

  library(SummarizedExperiment)
  library(ggplot2)
  library(tidyr)
  library(dplyr)

  mat <- matrix(runif(20), nrow = 5, ncol = 4)
  colnames(mat) <- c("S1_N", "S2_T", "S3_N", "S4_T")
  rownames(mat) <- paste0("G", 1:5)
  se <- SummarizedExperiment::SummarizedExperiment(assays = list(diversity = mat), rowData = DataFrame(genes = rownames(mat)), colData = DataFrame(samples = colnames(mat)))

  p <- plot_mean_violin(se)
  expect_s3_class(p, "gg")
})


test_that("plot_ma returns ggplot object", {
  skip_if_not_installed("ggplot2")
  library(ggplot2)

  df <- data.frame(Gene = paste0("G", 1:10), A_mean = runif(10), B_mean = runif(10), log2_fold_change = rnorm(10), adjusted_p_values = runif(10))
  p <- plot_ma(df)
  expect_s3_class(p, "gg")
})

testthat::test_that("plot_top_transcripts returns ggplot for synthetic data", {
  skip_on_cran()
  set.seed(42)
  counts <- matrix(rpois(3 * 8, lambda = 20), nrow = 3)
  rownames(counts) <- paste0("tx", 1:3)
  colnames(counts) <- paste0("S", 1:8)
  samples <- c(rep("Normal", 4), rep("Tumor", 4))
  # create simple tx2gene mapping
  tx2 <- data.frame(Transcript = rownames(counts), Gen = rep("GENE1", 3), stringsAsFactors = FALSE)

  p <- plot_top_transcripts(counts, gene = "GENE1", samples = samples, tx2gene = tx2, top_n = 2)
  testthat::expect_s3_class(p, "ggplot")
})
