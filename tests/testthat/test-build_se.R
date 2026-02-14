testthat::test_that("build_se constructs SummarizedExperiment from tx2gene data.frame", {
    # small toy dataset
    set.seed(2)
    n_tx <- 10L
    n_samps <- 3L
    genes <- paste0("G", seq_len(n_tx))
    readcounts <- matrix(sample(0:50, n_tx * n_samps, replace = TRUE), nrow = n_tx, ncol = n_samps)
    colnames(readcounts) <- paste0("S", seq_len(n_samps))
    tx2gene_df <- data.frame(Transcript = paste0("tx", seq_len(n_tx)), Gene = genes, stringsAsFactors = FALSE)
    rownames(readcounts) <- tx2gene_df$Transcript

    se <- build_se(readcounts, tx2gene_df)

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
    rownames(readcounts) <- tx2gene_df$Transcript

    se2 <- build_se(readcounts, tx2gene_df, assay_name = "mycounts")
    testthat::expect_s4_class(se2, "SummarizedExperiment")
    testthat::expect_true("tx2gene" %in% names(S4Vectors::metadata(se2)))
    testthat::expect_true("readcounts" %in% names(S4Vectors::metadata(se2)))
    testthat::expect_true("mycounts" %in% names(SummarizedExperiment::assays(se2)))
})

testthat::test_that("build_se errors on missing tx2gene path and unmatched transcript IDs", {
    # toy dataset
    set.seed(4)
    n_tx <- 6L
    n_samps <- 2L
    genes <- paste0("G", seq_len(n_tx))
    readcounts <- matrix(sample(0:20, n_tx * n_samps, replace = TRUE), nrow = n_tx, ncol = n_samps)
    colnames(readcounts) <- paste0("S", seq_len(n_samps))
    tx2gene_df <- data.frame(Transcript = paste0("tx", seq_len(n_tx)), Gene = genes, stringsAsFactors = FALSE)
    rownames(readcounts) <- tx2gene_df$Transcript

    testthat::expect_error(build_se("this_file_does_not_exist.tsv", readcounts))
    # Test with mismatched transcript IDs
    rownames(readcounts) <- paste0("wrong_tx", seq_len(n_tx))
    testthat::expect_error(build_se(readcounts, tx2gene_df))
})

context("build_se extra edge-case tests")

library(SummarizedExperiment)

# tx2gene provided as a file path should be read and stored in metadata
test_that("build_se accepts tx2gene as a file path and preserves metadata", {
    tx2 <- data.frame(Transcript = paste0("tx", 1:3), Gene = c("g1", "g1", "g2"), stringsAsFactors = FALSE)
    tf <- tempfile(fileext = ".tsv")
    utils::write.table(tx2, file = tf, sep = "\t", row.names = FALSE, quote = FALSE)

    rc <- matrix(c(5, 2, 9, 1, 0, 4), nrow = 3, byrow = FALSE)
    rownames(rc) <- tx2$Transcript

    se <- build_se(rc, tf, assay_name = "counts")
    md <- S4Vectors::metadata(se)
    expect_true(is.data.frame(md$tx2gene))
    expect_equal(md$tx2gene$Transcript, tx2$Transcript)
    expect_equal(S4Vectors::metadata(se)$readcounts, rc)
})

# readcounts can be a numeric data.frame and will be converted to matrix
test_that("build_se accepts numeric data.frame readcounts", {
    tx2 <- data.frame(Transcript = paste0("t", 1:2), Gene = c("g1", "g2"), stringsAsFactors = FALSE)
    rc_df <- data.frame(S1 = c(1, 3), S2 = c(2, 4))
    rownames(rc_df) <- tx2$Transcript

    se <- build_se(rc_df, tx2)
    expect_s4_class(se, "SummarizedExperiment")
    expect_true(is.matrix(SummarizedExperiment::assay(se, "counts")))
    expect_equal(SummarizedExperiment::assay(se, "counts"), as.matrix(rc_df))
})

# when readcounts has no rownames, errors should occur
test_that("build_se works when readcounts has no rownames", {
    tx2 <- data.frame(Transcript = paste0("tx", 1:2), Gene = c("g1", "g2"), stringsAsFactors = FALSE)
    rc <- matrix(c(1, 2, 3, 4), nrow = 2)
    rownames(rc) <- tx2$Transcript

    se <- build_se(rc, tx2)
    rd <- SummarizedExperiment::rowData(se)
    expect_equal(as.character(rd$genes), c("g1", "g2"))
    expect_equal(rownames(rd), rownames(rc))
})

# invalid tx2gene type should error
test_that("build_se errors on invalid tx2gene argument type", {
    rc <- matrix(1:4, nrow = 2)
    rownames(rc) <- c("tx1", "tx2")
    expect_error(build_se(rc, 12345), "'tx2gene' must be a path or a data.frame")
})

context("build_se additional tests")

library(SummarizedExperiment)

test_that("build_se accepts tx2gene data.frame and numeric matrix and sets metadata", {
    tx2 <- data.frame(Transcript = paste0("tx", 1:3), Gene = c("g1", "g1", "g2"), stringsAsFactors = FALSE)
    rc <- matrix(c(10, 0, 5, 2, 3, 1), nrow = 3)
    rownames(rc) <- tx2$Transcript
    se <- build_se(rc, tx2, assay_name = "counts")
    expect_s4_class(se, "SummarizedExperiment")
    md <- S4Vectors::metadata(se)
    expect_true(!is.null(md$tx2gene))
    expect_true(!is.null(md$readcounts))
    expect_equal(SummarizedExperiment::assayNames(se), "counts")
    expect_equal(SummarizedExperiment::rowData(se)$genes, c("g1", "g1", "g2"))
})

test_that("build_se errors on missing tx2gene file path", {
    rc <- matrix(1:4, nrow = 2)
    rownames(rc) <- c("tx1", "tx2")
    expect_error(build_se(rc, "/nonexistent/path.tsv"), "tx2gene file not found")
})

test_that("build_se errors on non-numeric readcounts or mismatched transcript IDs", {
    tx2 <- data.frame(Transcript = paste0("t", 1:2), Gene = c("g1", "g2"), stringsAsFactors = FALSE)
    rc_bad <- matrix(letters[1:4], nrow = 2)
    rownames(rc_bad) <- tx2$Transcript
    expect_error(build_se(rc_bad, tx2), "readcounts' must be a numeric")

    rc <- matrix(1:6, nrow = 3)
    rownames(rc) <- c("t1", "t2", "t3")
    expect_error(build_se(rc, tx2), "Some transcript IDs in readcounts were not found")
})
