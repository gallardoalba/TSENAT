context("Additional calculate_diversity tests: input types and edge cases")

library(SummarizedExperiment)

test_that("calculate_diversity supports data.frame input", {
  x_df <- as.data.frame(matrix(c(1,2,3,4,5,6), ncol=2))
  colnames(x_df) <- c("A","B")
  genes <- c("g1","g1","g2")
  res <- calculate_diversity(x_df, genes, q = 1)
  expect_s4_class(res, "SummarizedExperiment")
  expect_equal(ncol(res), ncol(as.matrix(x_df)))
})

test_that("calculate_diversity handles tximport-style list with tpm flag", {
  counts <- matrix(c(10,0,5,2,3,1), ncol = 2)
  abundance <- counts / rowSums(counts) * 1e6
  txlist <- list(counts = counts, abundance = abundance, length = NULL, countsFromAbundance = NULL)
  genes <- c("g1","g1","g2")
  # default uses counts
  res_counts <- calculate_diversity(txlist, genes, q = 1, tpm = FALSE)
  expect_s4_class(res_counts, "SummarizedExperiment")
  # using tpm should switch to abundance
  res_ab <- calculate_diversity(txlist, genes, q = 1, tpm = TRUE)
  expect_s4_class(res_ab, "SummarizedExperiment")
})

test_that("calculate_diversity handles object with class DGEList (list with class)", {
  counts <- matrix(c(1,2,3,4,5,6), ncol = 2)
  dgel <- list(counts = counts)
  class(dgel) <- "DGEList"
  genes <- c("g1","g1","g2")
  res <- calculate_diversity(dgel, genes, q = 2)
  expect_s4_class(res, "SummarizedExperiment")
})

test_that("calculate_diversity uses metadata$readcounts and tx2gene from a SummarizedExperiment", {
  rc <- matrix(c(5,5,0,1,2,3), ncol = 2)
  rownames(rc) <- c("tx1","tx2","tx3")
  tx2gene <- data.frame(Transcript = rownames(rc), Gen = c("g1","g1","g2"), stringsAsFactors = FALSE)
  se <- SummarizedExperiment(assays = SimpleList(dummy = matrix(0, nrow = 3, ncol = 2)))
  S4Vectors::metadata(se)$readcounts <- rc
  S4Vectors::metadata(se)$tx2gene <- tx2gene
  res <- calculate_diversity(se, genes = NULL, q = 1)
  expect_s4_class(res, "SummarizedExperiment")
  expect_equal(nrow(res), length(unique(tx2gene$Gen)))
})

test_that("calculate_diversity errors for invalid assay number in SummarizedExperiment", {
  se <- SummarizedExperiment(assays = SimpleList(a = matrix(1, nrow=2, ncol=2)))
  genes <- c("g1","g2")
  expect_error(calculate_diversity(se, genes, assayno = 2), "Please provide a valid assay number")
})

# Defensive / error cases

test_that("calculate_diversity errors on non-numeric input", {
  x <- matrix(letters[1:6], ncol = 2)
  genes <- c("g1","g1","g2")
  expect_error(calculate_diversity(x, genes), "Input data  must be numeric")
})

test_that("calculate_diversity errors on NA values", {
  x <- matrix(c(1, NA, 3, 4, 5, 6), ncol = 2)
  genes <- c("g1","g1","g2")
  expect_error(calculate_diversity(x, genes), "The data contains NA")
})

test_that("calculate_diversity errors when genes length mismatches rows", {
  x <- matrix(1:6, ncol = 2)
  genes <- c("g1","g2") # wrong length
  expect_error(calculate_diversity(x, genes), "The number of rows is not equal to the given gene set")
})

test_that("calculate_diversity errors for invalid q values (<=0)", {
  x <- matrix(1:6, ncol = 2)
  genes <- c("g1","g1","g2")
  expect_error(calculate_diversity(x, genes, q = 0), "Argument 'q' must be numeric and greater than 0")
})

test_that("calculate_diversity returns hill numbers when what='D'", {
  x <- matrix(c(1,2,3,4,5,6), ncol = 2)
  genes <- c("g1","g1","g2")
  res <- calculate_diversity(x, genes, q = 2, what = "D")
  expect_true("hill" %in% names(SummarizedExperiment::assays(res)))
})


