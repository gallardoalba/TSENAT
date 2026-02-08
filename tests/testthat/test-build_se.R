testthat::test_that("build_se constructs SummarizedExperiment from path", {
  pkgload::load_all()
  tx2gene_tsv <- system.file("extdata", "tx2gene.tsv", package = "TSENAT")
  data("tcga_brca_luma_dataset", package = "TSENAT", envir = globalenv())
  genes <- tcga_brca_luma_dataset[, 1]
  readcounts <- as.matrix(tcga_brca_luma_dataset[, -1])
  readcounts <- map_tx_to_readcounts(readcounts, tx2gene_tsv)

  se <- build_se(tx2gene_tsv, readcounts, genes)

  testthat::expect_s4_class(se, "SummarizedExperiment")
  testthat::expect_true("tx2gene" %in% names(S4Vectors::metadata(se)))
  testthat::expect_true("readcounts" %in% names(S4Vectors::metadata(se)))
  testthat::expect_equal(nrow(se), nrow(readcounts))
  testthat::expect_equal(dim(SummarizedExperiment::assay(se)), dim(readcounts))
  testthat::expect_equal(nrow(SummarizedExperiment::rowData(se)), nrow(readcounts))
})

testthat::test_that("build_se accepts tx2gene as data.frame and custom assay_name", {
  pkgload::load_all()
  tx2gene_tsv <- system.file("extdata", "tx2gene.tsv", package = "TSENAT")
  tx2gene_df <- utils::read.table(tx2gene_tsv, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
  data("tcga_brca_luma_dataset", package = "TSENAT", envir = globalenv())
  genes <- tcga_brca_luma_dataset[, 1]
  readcounts <- as.matrix(tcga_brca_luma_dataset[, -1])
  readcounts <- map_tx_to_readcounts(readcounts, tx2gene_tsv)

  se2 <- build_se(tx2gene_df, readcounts, genes, assay_name = "mycounts")
  testthat::expect_s4_class(se2, "SummarizedExperiment")
  testthat::expect_true("tx2gene" %in% names(S4Vectors::metadata(se2)))
  testthat::expect_true("readcounts" %in% names(S4Vectors::metadata(se2)))
  testthat::expect_true("mycounts" %in% names(SummarizedExperiment::assays(se2)))
})

testthat::test_that("build_se errors on missing tx2gene path and mismatched genes length", {
  pkgload::load_all()
  tx2gene_tsv <- system.file("extdata", "tx2gene.tsv", package = "TSENAT")
  data("tcga_brca_luma_dataset", package = "TSENAT", envir = globalenv())
  genes <- tcga_brca_luma_dataset[, 1]
  readcounts <- as.matrix(tcga_brca_luma_dataset[, -1])
  readcounts <- map_tx_to_readcounts(readcounts, tx2gene_tsv)

  testthat::expect_error(build_se("this_file_does_not_exist.tsv", readcounts, genes))
  testthat::expect_error(build_se(tx2gene_tsv, readcounts, genes[-1]))
})
