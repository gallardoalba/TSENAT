context("build_se additional tests")

library(SummarizedExperiment)

test_that("build_se accepts tx2gene data.frame and numeric matrix and sets metadata", {
    tx2 <- data.frame(Transcript = paste0("tx", 1:3), Gene = c("g1", "g1", "g2"), stringsAsFactors = FALSE)
    rc <- matrix(c(10, 0, 5, 2, 3, 1), nrow = 3)
    rownames(rc) <- tx2$Transcript
    genes <- c("g1", "g1", "g2")
    se <- build_se(tx2, rc, genes, assay_name = "counts")
    expect_s4_class(se, "SummarizedExperiment")
    md <- S4Vectors::metadata(se)
    expect_true(!is.null(md$tx2gene))
    expect_true(!is.null(md$readcounts))
    expect_equal(SummarizedExperiment::assayNames(se), "counts")
    expect_equal(SummarizedExperiment::rowData(se)$genes, genes)
})

test_that("build_se errors on missing tx2gene file path", {
    rc <- matrix(1:4, nrow = 2)
    genes <- c("g1", "g2")
    expect_error(build_se("/nonexistent/path.tsv", rc, genes), "tx2gene file not found")
})

test_that("build_se errors on non-numeric readcounts or mismatched genes length", {
    tx2 <- data.frame(Transcript = paste0("t", 1:2), Gene = c("g1", "g2"), stringsAsFactors = FALSE)
    rc_bad <- matrix(letters[1:4], nrow = 2)
    genes <- c("g1", "g2")
    expect_error(build_se(tx2, rc_bad, genes), "readcounts' must be a numeric")

    rc <- matrix(1:6, nrow = 3)
    genes_short <- c("g1", "g2")
    expect_error(build_se(tx2, rc, genes_short), "Length of 'genes' must equal nrow")
})
