testthat::test_that("build_se constructs SummarizedExperiment from tx2gene data.frame", {
  # small toy dataset
  set.seed(2)
  n_tx <- 10L
  n_samps <- 3L
  genes <- paste0("G", seq_len(n_tx))
  readcounts <- matrix(sample(0:50, n_tx * n_samps, replace = TRUE), nrow = n_tx, ncol = n_samps)
  colnames(readcounts) <- paste0("S", seq_len(n_samps))
  tx2gene_df <- data.frame(Transcript = paste0("tx", seq_len(n_tx)), Gene = genes, stringsAsFactors = FALSE)
  readcounts <- map_tx_to_readcounts(readcounts, tx2gene_df)

  se <- build_se(tx2gene_df, readcounts, genes)

  testthat::expect_s4_class(se, "SummarizedExperiment")
  testthat::expect_true("tx2gene" %in% names(S4Vectors::metadata(se)))
  testthat::expect_true("readcounts" %in% names(S4Vectors::metadata(se)))
  testthat::expect_equal(nrow(se), nrow(readcounts))
  testthat::expect_equal(dim(SummarizedExperiment::assay(se)), dim(readcounts))
  testthat::expect_equal(nrow(SummarizedExperiment::rowData(se)), nrow(readcounts))
})

testthat::test_that("build_se accepts tx2gene as data.frame and custom assay_name", {
  # toy dataset
  set.seed(3)
  n_tx <- 8L
  n_samps <- 2L
  genes <- paste0("G", seq_len(n_tx))
  readcounts <- matrix(sample(0:30, n_tx * n_samps, replace = TRUE), nrow = n_tx, ncol = n_samps)
  colnames(readcounts) <- paste0("S", seq_len(n_samps))
  tx2gene_df <- data.frame(Transcript = paste0("tx", seq_len(n_tx)), Gene = genes, stringsAsFactors = FALSE)
  readcounts <- map_tx_to_readcounts(readcounts, tx2gene_df)

  se2 <- build_se(tx2gene_df, readcounts, genes, assay_name = "mycounts")
  testthat::expect_s4_class(se2, "SummarizedExperiment")
  testthat::expect_true("tx2gene" %in% names(S4Vectors::metadata(se2)))
  testthat::expect_true("readcounts" %in% names(S4Vectors::metadata(se2)))
  testthat::expect_true("mycounts" %in% names(SummarizedExperiment::assays(se2)))
})

testthat::test_that("build_se errors on missing tx2gene path and mismatched genes length", {
  # toy dataset for mismatch test
  set.seed(4)
  n_tx <- 6L
  n_samps <- 2L
  genes <- paste0("G", seq_len(n_tx))
  readcounts <- matrix(sample(0:20, n_tx * n_samps, replace = TRUE), nrow = n_tx, ncol = n_samps)
  tx2gene_df <- data.frame(Transcript = paste0("tx", seq_len(n_tx)), Gene = genes, stringsAsFactors = FALSE)
  readcounts <- map_tx_to_readcounts(readcounts, tx2gene_df)

  testthat::expect_error(build_se("this_file_does_not_exist.tsv", readcounts, genes))
  testthat::expect_error(build_se(tx2gene_df, readcounts, genes[-1]))
})
