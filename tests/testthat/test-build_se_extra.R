context("build_se extra edge-case tests")

library(SummarizedExperiment)

# tx2gene provided as a file path should be read and stored in metadata
test_that("build_se accepts tx2gene as a file path and preserves metadata", {
    tx2 <- data.frame(Transcript = paste0("tx", 1:3), Gene = c("g1", "g1", "g2"), stringsAsFactors = FALSE)
    tf <- tempfile(fileext = ".tsv")
    utils::write.table(tx2, file = tf, sep = "\t", row.names = FALSE, quote = FALSE)

    rc <- matrix(c(5, 2, 9, 1, 0, 4), nrow = 3, byrow = FALSE)
    rownames(rc) <- tx2$Transcript
    genes <- c("g1", "g1", "g2")

    se <- build_se(tf, rc, genes, assay_name = "counts")
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
    genes <- c("g1", "g2")

    se <- build_se(tx2, rc_df, genes)
    expect_s4_class(se, "SummarizedExperiment")
    expect_true(is.matrix(SummarizedExperiment::assay(se, "counts")))
    expect_equal(SummarizedExperiment::assay(se, "counts"), as.matrix(rc_df))
})

# when readcounts has no rownames, rowData should have no rownames
test_that("build_se works when readcounts has no rownames", {
    tx2 <- data.frame(Transcript = paste0("tx", 1:2), Gene = c("g1", "g2"), stringsAsFactors = FALSE)
    rc <- matrix(c(1, 2, 3, 4), nrow = 2)
    genes <- c("g1", "g2")

    se <- build_se(tx2, rc, genes)
    rd <- SummarizedExperiment::rowData(se)
    expect_equal(as.character(rd$genes), genes)
    # rownames should be NULL when input rc has no rownames
    expect_true(is.null(rownames(rd)))
})

# invalid tx2gene type should error
test_that("build_se errors on invalid tx2gene argument type", {
    rc <- matrix(1:4, nrow = 2)
    genes <- c("g1", "g2")
    expect_error(build_se(12345, rc, genes), "'tx2gene_tsv' must be a path or a data.frame")
})
