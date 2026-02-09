context("filter_se more tests")

library(SummarizedExperiment)

test_that("filter_se errors on non-SE input", {
  expect_error(filter_se(1:10), "must be a SummarizedExperiment")
})

test_that("filter_se uses specified assay by name and falls back", {
  mat1 <- matrix(c(0,6,7,2,8,9), nrow = 3)
  mat2 <- matrix(1:6, nrow = 3)
  se <- SummarizedExperiment(assays = list(counts = mat1, other = mat2))
  res <- filter_se(se, min_count = 5, min_samples = 1, assay_name = "counts", verbose = FALSE)
  expect_s4_class(res, "SummarizedExperiment")
  # using missing assay name warns and uses first
  expect_warning(filter_se(se, min_count = 5, min_samples = 1, assay_name = "nope", verbose = FALSE))
})

test_that("filter_se respects min_samples cap and returns empty SE when none kept", {
  mat <- matrix(0, nrow = 5, ncol = 3)
  se <- SummarizedExperiment(assays = list(counts = mat))
  expect_warning(res <- filter_se(se, min_count = 0, min_samples = 5, verbose = FALSE))
  # all zeros and min_samples > available leads to zero kept and a warning
  expect_equal(nrow(SummarizedExperiment::assay(res,1)), 0)
})

test_that("filter_se subsets metadata readcounts and tx2gene", {
  mat <- matrix(c(0,6,7,2,8,9), nrow = 3)
  rownames(mat) <- paste0("tx", 1:3)
  se <- SummarizedExperiment(assays = list(counts = mat))
  # attach metadata
  S4Vectors::metadata(se)$readcounts <- mat
  S4Vectors::metadata(se)$tx2gene <- data.frame(Transcript = rownames(mat), Gene = c("g1","g1","g2"), stringsAsFactors = FALSE)
  res <- filter_se(se, min_count = 5, min_samples = 1, verbose = FALSE)
  md <- S4Vectors::metadata(res)
  expect_true(is.null(md$readcounts) == FALSE)
  expect_true(is.null(md$tx2gene) == FALSE)
  # tx2gene rows should subset to remaining transcripts
  expect_true(all(md$tx2gene$Transcript %in% rownames(SummarizedExperiment::assay(res))))
})
