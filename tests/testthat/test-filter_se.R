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
