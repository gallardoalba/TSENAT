testthat::test_that("filter_se keeps expected transcripts and updates metadata", {
  pkgload::load_all()
  tx2gene_tsv <- system.file("extdata", "tx2gene.tsv", package = "TSENAT")
  data("tcga_brca_luma_dataset", package = "TSENAT", envir = globalenv())
  genes <- tcga_brca_luma_dataset[, 1]
  readcounts <- as.matrix(tcga_brca_luma_dataset[, -1])
  readcounts <- map_tx_to_readcounts(readcounts, tx2gene_tsv)

  se <- build_se(tx2gene_tsv, readcounts, genes)
  se_f <- filter_se(se, min_count = 5L, min_samples = 5L, verbose = FALSE)

  testthat::expect_s4_class(se_f, "SummarizedExperiment")
  testthat::expect_true(nrow(se_f) < nrow(se))
  md <- S4Vectors::metadata(se_f)
  testthat::expect_true(is.null(md$readcounts) || nrow(md$readcounts) == nrow(se_f))
  if (!is.null(md$tx2gene)) testthat::expect_true(nrow(md$tx2gene) <= nrow(se))
})

testthat::test_that("filter_se warns when assay removes all rows", {
  pkgload::load_all()
  tx2gene_tsv <- system.file("extdata", "tx2gene.tsv", package = "TSENAT")
  data("tcga_brca_luma_dataset", package = "TSENAT", envir = globalenv())
  genes <- tcga_brca_luma_dataset[, 1]
  readcounts <- as.matrix(tcga_brca_luma_dataset[, -1])
  readcounts <- map_tx_to_readcounts(readcounts, tx2gene_tsv)

  se <- build_se(tx2gene_tsv, readcounts, genes)
  # choose thresholds that remove all rows
  testthat::expect_warning(sf <- filter_se(se, min_count = 1e6, min_samples = 1, verbose = FALSE))
  testthat::expect_s4_class(sf, "SummarizedExperiment")
  testthat::expect_true(nrow(sf) == 0)
})
