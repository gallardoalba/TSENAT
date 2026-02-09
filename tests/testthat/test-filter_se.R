testthat::test_that("filter_se keeps expected transcripts and updates metadata", {
  # create small toy dataset to avoid dependency on external data
  set.seed(1)
  n_tx <- 12L
  n_samps <- 4L
  genes <- paste0("G", sprintf("%03d", seq_len(n_tx)))
  readcounts <- matrix(sample(0:100, n_tx * n_samps, replace = TRUE), nrow = n_tx, ncol = n_samps)
  colnames(readcounts) <- paste0("S", seq_len(n_samps))
  tx2gene_df <- data.frame(
    Transcript = paste0("tx", seq_len(n_tx)),
    Gene = genes,
    stringsAsFactors = FALSE
  )
  readcounts <- map_tx_to_readcounts(readcounts, tx2gene_df)

  se <- build_se(tx2gene_df, readcounts, genes)
  se_f <- filter_se(se, min_count = 5L, min_samples = 5L, verbose = FALSE)

  testthat::expect_s4_class(se_f, "SummarizedExperiment")
  testthat::expect_true(nrow(se_f) < nrow(se))
  md <- S4Vectors::metadata(se_f)
  testthat::expect_true(is.null(md$readcounts) || nrow(md$readcounts) == nrow(se_f))
  if (!is.null(md$tx2gene)) testthat::expect_true(nrow(md$tx2gene) <= nrow(se))
})

testthat::test_that("filter_se warns when assay removes all rows", {
  # toy dataset matching the first test
  set.seed(1)
  n_tx <- 12L
  n_samps <- 4L
  genes <- paste0("G", sprintf("%03d", seq_len(n_tx)))
  readcounts <- matrix(sample(0:100, n_tx * n_samps, replace = TRUE), nrow = n_tx, ncol = n_samps)
  colnames(readcounts) <- paste0("S", seq_len(n_samps))
  tx2gene_df <- data.frame(
    Transcript = paste0("tx", seq_len(n_tx)),
    Gene = genes,
    stringsAsFactors = FALSE
  )
  readcounts <- map_tx_to_readcounts(readcounts, tx2gene_df)

  se <- build_se(tx2gene_df, readcounts, genes)
  # choose thresholds that remove all rows
  testthat::expect_warning(sf <- filter_se(se, min_count = 1e6, min_samples = 1, verbose = FALSE))
  testthat::expect_s4_class(sf, "SummarizedExperiment")
  testthat::expect_true(nrow(sf) == 0)
})

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
