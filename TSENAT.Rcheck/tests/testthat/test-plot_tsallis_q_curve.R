library(testthat)

context("plot_tsallis_q_curve")

test_that("plot_tsallis_q_curve returns a ggplot object", {
  if (!requireNamespace("SummarizedExperiment", quietly = TRUE)) skip("SummarizedExperiment required")
  if (!requireNamespace("ggplot2", quietly = TRUE)) skip("ggplot2 required")
  if (!requireNamespace("tidyr", quietly = TRUE)) skip("tidyr required")
  if (!requireNamespace("dplyr", quietly = TRUE)) skip("dplyr required")

  set.seed(1)
  # small synthetic dataset: 30 transcripts, 3 samples
  readcounts <- matrix(rpois(30 * 3, lambda = 10), nrow = 30, ncol = 3)
  colnames(readcounts) <- c("S1_N", "S2_T", "S3_N")
  genes <- rep(paste0("G", 1:10), length.out = nrow(readcounts))

  p <- plot_tsallis_q_curve(readcounts, genes, q_values = seq(0.01, 0.05, by = 0.01),
                            group_pattern = "_N$", group_names = c("Normal", "Tumor"))

  expect_true(inherits(p, "ggplot"))
})
