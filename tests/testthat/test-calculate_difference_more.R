context("calculate_difference additional tests")

library(SummarizedExperiment)

test_that("calculate_difference accepts SummarizedExperiment and uses sample_type from colData", {
  # build simple SE with 3 genes and 8 samples (4 vs 4)
  mat <- matrix(runif(3 * 8), nrow = 3)
  rownames(mat) <- c("g1","g2","g3")
  colnames(mat) <- paste0("S", 1:8)
  colData_df <- S4Vectors::DataFrame(sample_type = c(rep("Healthy",4), rep("Pathogenic",4)), row.names = colnames(mat))
  se <- SummarizedExperiment(assays = S4Vectors::SimpleList(counts = mat), colData = colData_df)

  res <- calculate_difference(se, samples = NULL, control = "Healthy", method = "mean", test = "wilcoxon")
  expect_true(is.data.frame(res))
  expect_true("raw_p_values" %in% colnames(res) || "adjusted_p_values" %in% colnames(res))
})

test_that("calculate_difference errors on invalid assayno for SummarizedExperiment", {
  mat <- matrix(runif(2 * 4), nrow = 2)
  colnames(mat) <- paste0("S", 1:4)
  colData_df <- S4Vectors::DataFrame(sample_type = c("A","A","B","B"), row.names = colnames(mat))
  se <- SummarizedExperiment(assays = S4Vectors::SimpleList(a = mat), colData = colData_df)
  expect_error(calculate_difference(se, samples = NULL, control = "A", assayno = 2), "Invalid 'assayno'|Column count doesn't match length")
})

test_that("Genes with insufficient observations are reported with NA p-values (small group)", {
  # create a larger data.frame (20 samples) where second gene has many NAs
  set.seed(123)
  samples <- c(rep("A", 10), rep("B", 10))
  # build matrix for three genes; g2 will be mostly NA in group A
  vals_g1 <- rnorm(20, 5, 1)
  vals_g2 <- c(rep(NA, 12), rnorm(8, 2, 0.5))
  vals_g3 <- rnorm(20, 7, 1)
  df <- data.frame(Genes = c("g1","g2","g3"), stringsAsFactors = FALSE)
  df <- cbind(df, as.data.frame(rbind(vals_g1, vals_g2, vals_g3)))
  colnames(df)[-1] <- paste0("S", seq_len(20))
  # call calculate_difference; suppress low-sample wilcoxon warning for this case
  res <- suppressWarnings(calculate_difference(df, samples = samples, control = "A", method = "mean", test = "wilcoxon"))
  expect_true(is.data.frame(res))
  # find g2 row and check NA p-values
  row_g2 <- res[res$genes == "g2", , drop = FALSE]
  expect_true(nrow(row_g2) == 1)
  expect_true(is.na(row_g2$raw_p_values) || is.na(row_g2$adjusted_p_values) )
})
