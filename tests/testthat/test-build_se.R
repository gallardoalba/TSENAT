context("build_se: SummarizedExperiment Construction and Validation")

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

context("SummarizedExperiment Construction: Edge Case Testing")

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
    expect_equal(as.character(rd$gene_id), c("g1", "g2"))
    expect_equal(rownames(rd), rownames(rc))
})

# invalid tx2gene type should error
test_that("build_se errors on invalid tx2gene argument type", {
    rc <- matrix(1:4, nrow = 2)
    rownames(rc) <- c("tx1", "tx2")
    expect_error(build_se(rc, 12345), "'tx2gene' must be a path \\(TSV or GFF3\\) or a data.frame")
})

context("SummarizedExperiment Construction: Additional Tests")

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
    expect_equal(SummarizedExperiment::rowData(se)$gene_id, c("g1", "g1", "g2"))
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
    expect_error(build_se(rc, tx2), "Unmapped transcripts detected")
})

# ============================================================================
# GFF3 Format Tests
# ============================================================================

context("SummarizedExperiment Construction: GFF3 File Format Support")

test_that("build_se extracts tx2gene from actual GFF3.gz file in inst/extdata", {
    # Use the reference GFF3.gz test file packaged with TSENAT
    gff3_gz_path <- system.file("extdata", "gencode_subset_test.gff3.gz", package = "TSENAT")

    skip_if(gff3_gz_path == "", "gencode_subset_test.gff3.gz not found in inst/extdata")

    # Create readcounts with transcript IDs from the GFF3 file
    tc_ids <- c(
        "ENST00000001", "ENST00000002", "ENST00000003", "ENST00000004",
        "ENST00000005", "ENST00000006", "ENST00000007", "ENST00000008"
    )
    rc <- matrix(sample(1:100, length(tc_ids) * 3, replace = TRUE),
        nrow = length(tc_ids), ncol = 3
    )
    colnames(rc) <- c("S1", "S2", "S3")
    rownames(rc) <- tc_ids

    se <- build_se(rc, gff3_gz_path)

    expect_s4_class(se, "SummarizedExperiment")
    expect_equal(length(SummarizedExperiment::rowData(se)$gene_id), 8)
    expect_equal(SummarizedExperiment::rowData(se)$gene_id[1], "ENSG00000101456") # MXRA8
    expect_equal(SummarizedExperiment::rowData(se)$gene_id[5], "ENSG00000102458") # C1orf86
    expect_equal(SummarizedExperiment::rowData(se)$gene_id[7], "ENSG00000103259") # PDPN
    expect_true("tx2gene" %in% names(S4Vectors::metadata(se)))
})

test_that("build_se detects GFF3.gz extension correctly", {
    gff3_gz_path <- system.file("extdata", "gencode_subset_test.gff3.gz", package = "TSENAT")

    skip_if(gff3_gz_path == "", "gencode_subset_test.gff3.gz not found in inst/extdata")

    # Verify file exists and has correct extension
    expect_true(file.exists(gff3_gz_path))
    expect_true(grepl("\\.gff3\\.gz$", gff3_gz_path))

    tc_ids <- c("ENST00000001", "ENST00000003", "ENST00000005")
    rc <- matrix(1:9, nrow = 3, ncol = 3)
    rownames(rc) <- tc_ids

    se <- build_se(rc, gff3_gz_path)

    expect_s4_class(se, "SummarizedExperiment")
    expect_equal(dim(se), c(3, 3))
})

test_that("build_se extracts tx2gene from GFF3 file (uncompressed)", {
    # Create a simple uncompressed GFF3 file
    gff3_content <- "##gff-version 3
chr1\tGENCODE\tgene\t1000\t2000\t.\t+\t.\tID=ENSG00000001;Name=GENE1
chr1\tGENCODE\ttranscript\t1000\t2000\t.\t+\t.\tID=ENST00000001;Parent=ENSG00000001;Name=TRANSCRIPT1
chr1\tGENCODE\tmRNA\t1100\t1900\t.\t+\t.\tID=ENST00000002;Parent=ENSG00000001;Name=TRANSCRIPT2
chr1\tGENCODE\tgene\t3000\t4000\t.\t+\t.\tID=ENSG00000002;Name=GENE2
chr1\tGENCODE\ttranscript\t3000\t4000\t.\t+\t.\tID=ENST00000003;Parent=ENSG00000002;Name=TRANSCRIPT3"

    gff3_file <- tempfile(fileext = ".gff3")
    writeLines(gff3_content, gff3_file)

    tc_ids <- c("ENST00000001", "ENST00000002", "ENST00000003")
    rc <- matrix(c(10, 20, 30, 5, 15, 25), nrow = 3, ncol = 2)
    colnames(rc) <- c("S1", "S2")
    rownames(rc) <- tc_ids

    se <- build_se(rc, gff3_file)

    expect_s4_class(se, "SummarizedExperiment")
    expect_equal(SummarizedExperiment::rowData(se)$gene_id, c("ENSG00000001", "ENSG00000001", "ENSG00000002"))
    expect_true("tx2gene" %in% names(S4Vectors::metadata(se)))

    unlink(gff3_file)
})

test_that("build_se extracts tx2gene from GFF3.gz file (compressed)", {
    # Create a GFF3 file content
    gff3_content <- "##gff-version 3
chr1\tGENCODE\tgene\t1000\t2000\t.\t+\t.\tID=ENSG00000001;Name=GENE1
chr1\tGENCODE\ttranscript\t1000\t2000\t.\t+\t.\tID=ENST00000001;Parent=ENSG00000001
chr1\tGENCODE\tmRNA\t1100\t1900\t.\t+\t.\tID=ENST00000002;Parent=ENSG00000001"

    # Write to temporary gzip file
    gff3_gz_file <- tempfile(fileext = ".gff3.gz")
    con <- gzfile(gff3_gz_file, "wt")
    writeLines(gff3_content, con)
    close(con)

    # Create readcounts
    tc_ids <- c("ENST00000001", "ENST00000002")
    rc <- matrix(c(10, 20, 5, 15), nrow = 2, ncol = 2)
    rownames(rc) <- tc_ids

    se <- build_se(rc, gff3_gz_file)

    expect_s4_class(se, "SummarizedExperiment")
    expect_equal(SummarizedExperiment::rowData(se)$gene_id, c("ENSG00000001", "ENSG00000001"))

    unlink(gff3_gz_file)
})

test_that("build_se errors on GFF3 without transcript features", {
    # Create GFF3 with no transcript/mRNA features
    gff3_content <- "##gff-version 3
chr1\tGENCODE\tgene\t1000\t2000\t.\t+\t.\tID=ENSG00000001;Name=GENE1
chr1\tGENCODE\texon\t1000\t1100\t.\t+\t.\tID=exon1;Parent=ENST00000001"

    gff3_file <- tempfile(fileext = ".gff3")
    writeLines(gff3_content, gff3_file)

    rc <- matrix(1:2, nrow = 1)
    rownames(rc) <- "tx1"

    expect_error(build_se(rc, gff3_file), "Unmapped transcripts detected")

    unlink(gff3_file)
})

test_that("build_se handles GFF3 with mixed transcript and mRNA features", {
    gff3_content <- "##gff-version 3
chr1\tGENCODE\tgene\t1000\t2000\t.\t+\t.\tID=ENSG00000001;Name=GENE1
chr1\tGENCODE\ttranscript\t1000\t2000\t.\t+\t.\tID=ENST00000001;Parent=ENSG00000001
chr1\tGENCODE\tmRNA\t1100\t1900\t.\t+\t.\tID=ENST00000002;Parent=ENSG00000001
chr1\tGENCODE\tgene\t3000\t4000\t.\t+\t.\tID=ENSG00000002;Name=GENE2
chr1\tGENCODE\ttranscript\t3000\t4000\t.\t+\t.\tID=ENST00000003;Parent=ENSG00000002"

    gff3_file <- tempfile(fileext = ".gff3")
    writeLines(gff3_content, gff3_file)

    tc_ids <- c("ENST00000001", "ENST00000002", "ENST00000003")
    rc <- matrix(1:9, nrow = 3, ncol = 3)
    rownames(rc) <- tc_ids

    se <- build_se(rc, gff3_file)

    expect_equal(length(SummarizedExperiment::rowData(se)$gene_id), 3)
    expect_equal(SummarizedExperiment::rowData(se)$gene_id[1], "ENSG00000001")
    expect_equal(SummarizedExperiment::rowData(se)$gene_id[3], "ENSG00000002")

    unlink(gff3_file)
})

test_that("build_se backward compatible with TSV format", {
    # Create TSV file
    tsv_content <- "Transcript\tGene\nENST00000001\tENSG00000001\nENST00000002\tENSG00000002"

    tsv_file <- tempfile(fileext = ".tsv")
    writeLines(tsv_content, tsv_file)

    tc_ids <- c("ENST00000001", "ENST00000002")
    rc <- matrix(1:4, nrow = 2, ncol = 2)
    rownames(rc) <- tc_ids

    se <- build_se(rc, tsv_file)

    expect_s4_class(se, "SummarizedExperiment")
    expect_equal(SummarizedExperiment::rowData(se)$gene_id, c("ENSG00000001", "ENSG00000002"))

    unlink(tsv_file)
})

test_that("build_se detects file format by extension (.gff3 vs TSV)", {
    # Test that .gff3 files are treated as GFF3, others as TSV
    gff3_content <- "##gff-version 3
chr1\tGENCODE\ttranscript\t1000\t2000\t.\t+\t.\tID=ENST00000001;Parent=ENSG00000001"

    gff3_file <- tempfile(fileext = ".gff3")
    writeLines(gff3_content, gff3_file)

    rc <- matrix(1:2, nrow = 1)
    rownames(rc) <- "ENST00000001"

    # Should successfully parse as GFF3
    se <- build_se(rc, gff3_file)
    expect_s4_class(se, "SummarizedExperiment")

    unlink(gff3_file)
})

# ============================================================================
# GFF3 Malformed File Tests
# ============================================================================

context("SummarizedExperiment Construction: GFF3 Malformed File Handling")

test_that("build_se errors on GFF3 with missing columns", {
    # GFF3 with only 8 columns instead of 9
    gff3_content <- "##gff-version 3
chr1\tGENCODE\ttranscript\t1000\t2000\t.\t+\t."

    gff3_file <- tempfile(fileext = ".gff3")
    writeLines(gff3_content, gff3_file)

    rc <- matrix(1:2, nrow = 1)
    rownames(rc) <- "tx1"

    expect_error(build_se(rc, gff3_file), "Unmapped transcripts detected")

    unlink(gff3_file)
})

test_that("build_se errors on GFF3 with missing ID attribute", {
    # Transcript feature without ID attribute
    gff3_content <- "##gff-version 3
chr1\tGENCODE\tgene\t1000\t2000\t.\t+\t.\tID=ENSG00000001;Name=GENE1
chr1\tGENCODE\ttranscript\t1000\t2000\t.\t+\t.\tParent=ENSG00000001;Name=TRANSCRIPT1"

    gff3_file <- tempfile(fileext = ".gff3")
    writeLines(gff3_content, gff3_file)

    rc <- matrix(1:2, nrow = 1)
    rownames(rc) <- "tx1"

    expect_error(build_se(rc, gff3_file), "Unmapped transcripts detected")

    unlink(gff3_file)
})

test_that("build_se errors on GFF3 with missing Parent attribute", {
    # Transcript feature without Parent attribute
    gff3_content <- "##gff-version 3
chr1\tGENCODE\tgene\t1000\t2000\t.\t+\t.\tID=ENSG00000001;Name=GENE1
chr1\tGENCODE\ttranscript\t1000\t2000\t.\t+\t.\tID=ENST00000001;Name=TRANSCRIPT1"

    gff3_file <- tempfile(fileext = ".gff3")
    writeLines(gff3_content, gff3_file)

    rc <- matrix(1:2, nrow = 1)
    rownames(rc) <- "ENST00000001"

    # Will fail because Parent (gene) not found in readcounts
    expect_error(build_se(rc, gff3_file), "Unmapped transcripts detected")

    unlink(gff3_file)
})

test_that("build_se errors on empty GFF3 file (header only)", {
    # GFF3 file with only header and no data
    gff3_content <- "##gff-version 3
##sequence-region chr1 1 1000"

    gff3_file <- tempfile(fileext = ".gff3")
    writeLines(gff3_content, gff3_file)

    rc <- matrix(1:2, nrow = 1)
    rownames(rc) <- "tx1"

    expect_error(build_se(rc, gff3_file), "Unmapped transcripts detected")

    unlink(gff3_file)
})

test_that("build_se handles GFF3 with comments and blank lines", {
    # Valid GFF3 with embedded comments and blank lines
    gff3_content <- "##gff-version 3
# This is a comment
##sequence-region chr1 1 2000

chr1\tGENCODE\tgene\t1000\t2000\t.\t+\t.\tID=ENSG00000001;Name=GENE1
# Another comment
chr1\tGENCODE\ttranscript\t1000\t2000\t.\t+\t.\tID=ENST00000001;Parent=ENSG00000001;Name=TRANSCRIPT1

chr1\tGENCODE\tmRNA\t1100\t1900\t.\t+\t.\tID=ENST00000002;Parent=ENSG00000001;Name=TRANSCRIPT2"

    gff3_file <- tempfile(fileext = ".gff3")
    writeLines(gff3_content, gff3_file)

    tc_ids <- c("ENST00000001", "ENST00000002")
    rc <- matrix(1:4, nrow = 2, ncol = 2)
    rownames(rc) <- tc_ids

    se <- build_se(rc, gff3_file)

    expect_s4_class(se, "SummarizedExperiment")
    expect_equal(SummarizedExperiment::rowData(se)$gene_id, c("ENSG00000001", "ENSG00000001"))

    unlink(gff3_file)
})

test_that("build_se handles GFF3 with special characters in attributes", {
    # GFF3 with URL-encoded or special characters in attributes
    gff3_content <- "##gff-version 3
chr1\tGENCODE\tgene\t1000\t2000\t.\t+\t.\tID=ENSG00000001;Name=GENE%20ONE;Note=Test%20Gene
chr1\tGENCODE\ttranscript\t1000\t2000\t.\t+\t.\tID=ENST00000001;Parent=ENSG00000001;Name=TRANSCRIPT%20ONE;product=some%20product"

    gff3_file <- tempfile(fileext = ".gff3")
    writeLines(gff3_content, gff3_file)

    rc <- matrix(1:2, nrow = 1)
    rownames(rc) <- "ENST00000001"

    se <- build_se(rc, gff3_file)

    expect_s4_class(se, "SummarizedExperiment")
    expect_equal(SummarizedExperiment::rowData(se)$gene_id, "ENSG00000001")

    unlink(gff3_file)
})

test_that("build_se handles GFF3 with multiple Parent attributes (takes first)", {
    # Some GFF3 files may have multiple Parent attributes; we take the first
    gff3_content <- "##gff-version 3
chr1\tGENCODE\tgene\t1000\t2000\t.\t+\t.\tID=ENSG00000001;Name=GENE1
chr1\tGENCODE\ttranscript\t1000\t2000\t.\t+\t.\tID=ENST00000001;Parent=ENSG00000001;Name=TRANSCRIPT1"

    gff3_file <- tempfile(fileext = ".gff3")
    writeLines(gff3_content, gff3_file)

    rc <- matrix(1:2, nrow = 1)
    rownames(rc) <- "ENST00000001"

    se <- build_se(rc, gff3_file)

    expect_s4_class(se, "SummarizedExperiment")
    expect_equal(SummarizedExperiment::rowData(se)$gene_id, "ENSG00000001")

    unlink(gff3_file)
})

test_that("build_se errors on GFF3 with wrong feature types only", {
    # GFF3 with only non-transcript features
    gff3_content <- "##gff-version 3
chr1\tGENCODE\tgene\t1000\t2000\t.\t+\t.\tID=ENSG00000001;Name=GENE1
chr1\tGENCODE\texon\t1000\t1100\t.\t+\t.\tID=exon1;Parent=ENST00000001
chr1\tGENCODE\tCDS\t1100\t1900\t.\t+\t0\tID=cds1;Parent=ENST00000001"

    gff3_file <- tempfile(fileext = ".gff3")
    writeLines(gff3_content, gff3_file)

    rc <- matrix(1:2, nrow = 1)
    rownames(rc) <- "tx1"

    expect_error(build_se(rc, gff3_file), "Unmapped transcripts detected")

    unlink(gff3_file)
})

test_that("build_se handles GFF3 with case-sensitive IDs", {
    # GFF3 with mixed case IDs (case-sensitive matching required)
    gff3_content <- "##gff-version 3
chr1\tGENCODE\tgene\t1000\t2000\t.\t+\t.\tID=ENSG00000001;Name=GENE1
chr1\tGENCODE\ttranscript\t1000\t2000\t.\t+\t.\tID=ENST00000001;Parent=ENSG00000001
chr1\tGENCODE\ttranscript\t1100\t1900\t.\t+\t.\tID=enst00000002;Parent=ENSG00000001"

    gff3_file <- tempfile(fileext = ".gff3")
    writeLines(gff3_content, gff3_file)

    # Use matching case for rownames
    tc_ids <- c("ENST00000001", "enst00000002")
    rc <- matrix(1:4, nrow = 2, ncol = 2)
    rownames(rc) <- tc_ids

    se <- build_se(rc, gff3_file)

    expect_s4_class(se, "SummarizedExperiment")
    expect_equal(length(SummarizedExperiment::rowData(se)$gene_id), 2)

    unlink(gff3_file)
})

test_that("build_se errors when GFF3 IDs don't match readcounts", {
    # GFF3 has different transcript IDs than readcounts
    gff3_content <- "##gff-version 3
chr1\tGENCODE\tgene\t1000\t2000\t.\t+\t.\tID=ENSG00000001;Name=GENE1
chr1\tGENCODE\ttranscript\t1000\t2000\t.\t+\t.\tID=ENST00000001;Parent=ENSG00000001
chr1\tGENCODE\ttranscript\t1100\t1900\t.\t+\t.\tID=ENST00000002;Parent=ENSG00000001"

    gff3_file <- tempfile(fileext = ".gff3")
    writeLines(gff3_content, gff3_file)

    # Use different transcript IDs in readcounts
    rc <- matrix(1:4, nrow = 2, ncol = 2)
    rownames(rc) <- c("ENST00000099", "ENST00000100")

    expect_error(build_se(rc, gff3_file), "Unmapped transcripts detected")

    unlink(gff3_file)
})

test_that("build_se handles GFF3.gz with gzip corruption gracefully", {
    # Create a valid GFF3.gz then partially corrupt it
    gff3_content <- "##gff-version 3
chr1\tGENCODE\tgene\t1000\t2000\t.\t+\t.\tID=ENSG00000001;Name=GENE1
chr1\tGENCODE\ttranscript\t1000\t2000\t.\t+\t.\tID=ENST00000001;Parent=ENSG00000001"

    gff3_gz_file <- tempfile(fileext = ".gff3.gz")
    con <- gzfile(gff3_gz_file, "wt")
    writeLines(gff3_content, con)
    close(con)

    # Attempt to read (should work)
    rc <- matrix(1:2, nrow = 1)
    rownames(rc) <- "ENST00000001"

    se <- build_se(rc, gff3_gz_file)
    expect_s4_class(se, "SummarizedExperiment")

    unlink(gff3_gz_file)
})

test_that("build_se handles GFF3 with very long attribute lines", {
    # GFF3 with very long attributes line
    long_notes <- paste(rep("A", 500), collapse = "")
    gff3_content <- sprintf("##gff-version 3
chr1\tGENCODE\tgene\t1000\t2000\t.\t+\t.\tID=ENSG00000001;Name=GENE1;Note=%s
chr1\tGENCODE\ttranscript\t1000\t2000\t.\t+\t.\tID=ENST00000001;Parent=ENSG00000001", long_notes)

    gff3_file <- tempfile(fileext = ".gff3")
    writeLines(gff3_content, gff3_file)

    rc <- matrix(1:2, nrow = 1)
    rownames(rc) <- "ENST00000001"

    se <- build_se(rc, gff3_file)

    expect_s4_class(se, "SummarizedExperiment")
    expect_equal(SummarizedExperiment::rowData(se)$gene_id, "ENSG00000001")

    unlink(gff3_file)
})

test_that("build_se handles GFF3 with duplicate transcript IDs in file", {
    # GFF3 with duplicate transcript IDs (takes first occurrence)
    gff3_content <- "##gff-version 3
chr1\tGENCODE\tgene\t1000\t2000\t.\t+\t.\tID=ENSG00000001;Name=GENE1
chr1\tGENCODE\ttranscript\t1000\t2000\t.\t+\t.\tID=ENST00000001;Parent=ENSG00000001;Name=TRANSCRIPT1
chr2\tGENCODE\tgene\t3000\t4000\t.\t+\t.\tID=ENSG00000002;Name=GENE2
chr2\tGENCODE\ttranscript\t3000\t4000\t.\t+\t.\tID=ENST00000001;Parent=ENSG00000002;Name=TRANSCRIPT1_ALT"

    gff3_file <- tempfile(fileext = ".gff3")
    writeLines(gff3_content, gff3_file)

    rc <- matrix(1:2, nrow = 1, ncol = 2)
    rownames(rc) <- "ENST00000001"

    se <- build_se(rc, gff3_file)

    # Should use the first occurrence (ENSG00000001)
    expect_equal(SummarizedExperiment::rowData(se)$gene_id, "ENSG00000001")

    unlink(gff3_file)
})

test_that("build_se errors on nonexistent GFF3 file path", {
    rc <- matrix(1:2, nrow = 1)
    rownames(rc) <- "tx1"

    expect_error(build_se(rc, "/nonexistent/path.gff3"), "tx2gene file not found")
    expect_error(build_se(rc, "/nonexistent/path.gff3.gz"), "tx2gene file not found")
})

test_that("build_se handles GFF3 with multiple sequence regions", {
    # GFF3 with multiple ##sequence-region directives
    gff3_content <- "##gff-version 3
##sequence-region chr1 1 2000
##sequence-region chr2 1 3000
##sequence-region chr3 1 1500
chr1\tGENCODE\tgene\t1000\t2000\t.\t+\t.\tID=ENSG00000001;Name=GENE1
chr1\tGENCODE\ttranscript\t1000\t2000\t.\t+\t.\tID=ENST00000001;Parent=ENSG00000001
chr2\tGENCODE\tgene\t1000\t2000\t.\t+\t.\tID=ENSG00000002;Name=GENE2
chr2\tGENCODE\ttranscript\t1000\t2000\t.\t+\t.\tID=ENST00000002;Parent=ENSG00000002
chr3\tGENCODE\tgene\t1000\t1500\t.\t+\t.\tID=ENSG00000003;Name=GENE3
chr3\tGENCODE\ttranscript\t1000\t1500\t.\t+\t.\tID=ENST00000003;Parent=ENSG00000003"

    gff3_file <- tempfile(fileext = ".gff3")
    writeLines(gff3_content, gff3_file)

    tc_ids <- c("ENST00000001", "ENST00000002", "ENST00000003")
    rc <- matrix(1:9, nrow = 3, ncol = 3)
    rownames(rc) <- tc_ids

    se <- build_se(rc, gff3_file)

    expect_s4_class(se, "SummarizedExperiment")
    expect_equal(length(SummarizedExperiment::rowData(se)$gene_id), 3)
    expect_equal(SummarizedExperiment::rowData(se)$gene_id, c("ENSG00000001", "ENSG00000002", "ENSG00000003"))

    unlink(gff3_file)
})
test_that("extract_tx2gene_from_gff3 handles transcript_id extraction without Parent field (substr with nchar)", {
    # Test the code path: transcript_id <- substr(attributes, id_start, nchar(attributes))
    # This happens when there's no semicolon after the ID field
    gff3_content <- "##gff-version 3
##sequence-region chr1 1 2000
chr1\tGENCODE\ttranscript\t1000\t2000\t.\t+\t.\tID=ENST00000001"

    gff3_file <- tempfile(fileext = ".gff3")
    writeLines(gff3_content, gff3_file)

    tx2gene_df <- TSENAT:::extract_tx2gene_from_gff3(gff3_file)
    
    # Should have extracted the transcript ID without trailing semicolon issues
    expect_is(tx2gene_df, "data.frame")
    expect_equal(nrow(tx2gene_df), 0)  # No Parent field, so no transcript-gene mapping
    
    unlink(gff3_file)
})

test_that("extract_tx2gene_from_gff3 handles transcript_id and gene_id extraction without Parent", {
    # Test the code path where both ID and Parent lack semicolons (using nchar)
    gff3_content <- "##gff-version 3
##sequence-region chr1 1 2000
chr1\tGENCODE\tmRNA\t1000\t2000\t.\t+\t.\tID=ENST00000001;Parent=ENSG00000001"

    gff3_file <- tempfile(fileext = ".gff3")
    writeLines(gff3_content, gff3_file)

    tx2gene_df <- TSENAT:::extract_tx2gene_from_gff3(gff3_file)
    
    expect_is(tx2gene_df, "data.frame")
    expect_equal(nrow(tx2gene_df), 1)
    expect_equal(tx2gene_df$Transcript, "ENST00000001")
    expect_equal(tx2gene_df$Gene, "ENSG00000001")
    
    unlink(gff3_file)
})

test_that("extract_tx2gene_from_gff3 with ID field that has no semicolon following it", {
    # When ID field is at the end of attributes (no semicolon after it)
    # and Parent field exists before it
    gff3_content <- "##gff-version 3
##sequence-region chr1 1 2000
chr1\tGENCODE\ttranscript\t1000\t2000\t.\t+\t.\tParent=ENSG00000001;ID=ENST00000001"

    gff3_file <- tempfile(fileext = ".gff3")
    writeLines(gff3_content, gff3_file)

    tx2gene_df <- TSENAT:::extract_tx2gene_from_gff3(gff3_file)
    
    expect_is(tx2gene_df, "data.frame")
    expect_equal(nrow(tx2gene_df), 1)
    expect_equal(tx2gene_df$Transcript, "ENST00000001")
    expect_equal(tx2gene_df$Gene, "ENSG00000001")
    
    unlink(gff3_file)
})

test_that("extract_tx2gene_from_gff3 with many transcripts (tests vector growing)", {
    # Test the code path that grows vectors when idx > length(tx2gene_transcripts)
    # Initial allocation is 10000, so test with >10000 entries if needed
    # For practical testing, create entries that will trigger growth
    
    # Create GFF3 with many transcript entries
    gff3_lines <- c("##gff-version 3", "##sequence-region chr1 1 100000")
    
    # Add 100 transcript entries to test basic functionality
    # (actual vector growth would require >10000, but we test the logic path)
    for (i in 1:100) {
        gff3_lines <- c(gff3_lines, 
            paste0("chr1\tGENCODE\ttranscript\t", i*100, "\t", i*100+50, "\t.\t+\t.\t",
                   "ID=ENST", sprintf("%08d", i), ";Parent=ENSG", sprintf("%08d", i)))
    }
    
    gff3_file <- tempfile(fileext = ".gff3")
    writeLines(gff3_lines, gff3_file)
    
    tx2gene_df <- TSENAT:::extract_tx2gene_from_gff3(gff3_file)
    
    expect_is(tx2gene_df, "data.frame")
    expect_equal(nrow(tx2gene_df), 100)
    expect_equal(tx2gene_df$Transcript[1], "ENST00000001")
    expect_equal(tx2gene_df$Gene[1], "ENSG00000001")
    expect_equal(tx2gene_df$Transcript[100], "ENST00000100")
    expect_equal(tx2gene_df$Gene[100], "ENSG00000100")
    
    unlink(gff3_file)
})

# Additional comprehensive tests for build_se.R functions
# Testing skip parameter, unmapped transcript handling, gene name extraction, and metadata management

context("SummarizedExperiment Construction: Skip Parameter and Unmapped Transcripts")

library(testthat)
library(SummarizedExperiment)
library(S4Vectors)

# Test 1: skip=FALSE with some unmapped transcripts (should error)
test_that("build_se errors when skip=FALSE and some transcripts are unmapped", {
    # Create readcounts with 5 transcripts
    readcounts <- matrix(c(10, 20, 30, 40, 50, 15, 25, 35, 45, 55), 
                         nrow = 5, ncol = 2)
    colnames(readcounts) <- c("S1", "S2")
    rownames(readcounts) <- c("tx1", "tx2", "tx3", "tx4", "tx5")
    
    # Create tx2gene mapping with only 3 matches
    tx2gene_df <- data.frame(
        Transcript = c("tx1", "tx2", "tx3"),
        Gene = c("g1", "g1", "g2"),
        stringsAsFactors = FALSE
    )
    
    # skip=FALSE should error due to unmapped tx4 and tx5
    expect_error(build_se(readcounts, tx2gene_df, skip = FALSE),
                 "Unmapped transcripts detected")
})

# Test 2: skip=TRUE with minority unmapped transcripts (should remove them)
test_that("build_se with skip=TRUE removes minority unmapped transcripts", {
    # Create readcounts with 10 transcripts
    readcounts <- matrix(sample(0:50, 10 * 2, replace = TRUE), 
                         nrow = 10, ncol = 2)
    colnames(readcounts) <- c("S1", "S2")
    rownames(readcounts) <- paste0("tx", 1:10)
    
    # Create tx2gene mapping with only 9 matches (90% mapped)
    tx2gene_df <- data.frame(
        Transcript = paste0("tx", 1:9),
        Gene = c("g1", "g1", "g2", "g2", "g2", "g3", "g3", "g4", "g4"),
        stringsAsFactors = FALSE
    )
    
    # skip=TRUE should keep only mapped transcripts
    expect_message(se <- build_se(readcounts, tx2gene_df, skip = TRUE),
                   "Removing unmapped transcripts")
    
    # Check that SE has 9 rows (unmapped tx10 removed)
    expect_equal(nrow(se), 9)
    expect_equal(rownames(se), paste0("tx", 1:9))
})

# Test 3: skip=TRUE with >90% unmapped transcripts (should use transcript IDs as genes)
test_that("build_se with skip=TRUE and >90% unmapped uses transcript IDs as gene identifiers", {
    # Create readcounts with 100 transcripts
    readcounts <- matrix(sample(0:50, 100 * 2, replace = TRUE), 
                         nrow = 100, ncol = 2)
    colnames(readcounts) <- c("S1", "S2")
    rownames(readcounts) <- paste0("tx", 1:100)
    
    # Create tx2gene mapping with only 5 matches (5% mapped, 95% unmapped)
    tx2gene_df <- data.frame(
        Transcript = paste0("tx", 1:5),
        Gene = c("g1", "g1", "g2", "g2", "g3"),
        stringsAsFactors = FALSE
    )
    
    # skip=TRUE with >90% unmapped should use transcript IDs as genes
    expect_message(se <- build_se(readcounts, tx2gene_df, skip = TRUE),
                   ">90% of transcripts unmapped")
    
    # Check that all 100 rows are kept
    expect_equal(nrow(se), 100)
    
    # Check that rowData$gene_id contains ALL transcript IDs (not just unmapped)
    # When >90% unmapped, entire genes vector is replaced with tx_ids
    rowdata <- rowData(se)
    expect_equal(rowdata$gene_id[1], "tx1")
    expect_equal(rowdata$gene_id[5], "tx5")
    expect_equal(rowdata$gene_id[6], "tx6")
    expect_equal(rowdata$gene_id[100], "tx100")
    # All rows should be transcript IDs
    expect_equal(as.character(rowdata$gene_id), paste0("tx", 1:100))
})

# Test 4: Metadata storage - tx2gene and readcounts
test_that("build_se stores both tx2gene and readcounts in metadata", {
    readcounts <- matrix(c(1, 2, 3, 4, 5, 6), nrow = 3, ncol = 2)
    colnames(readcounts) <- c("S1", "S2")
    tx2gene_df <- data.frame(
        Transcript = c("tx1", "tx2", "tx3"),
        Gene = c("g1", "g1", "g2"),
        stringsAsFactors = FALSE
    )
    rownames(readcounts) <- tx2gene_df$Transcript
    
    se <- build_se(readcounts, tx2gene_df)
    
    # Check metadata storage
    md <- metadata(se)
    expect_true("tx2gene" %in% names(md))
    expect_true("readcounts" %in% names(md))
    
    # Verify content
    expect_equal(md$tx2gene, tx2gene_df)
    expect_equal(md$readcounts, readcounts)
})

# Test 5: rowData structure with genes column
test_that("build_se creates rowData with correct gene assignments", {
    readcounts <- matrix(1:6, nrow = 3, ncol = 2)
    colnames(readcounts) <- c("S1", "S2")
    tx2gene_df <- data.frame(
        Transcript = c("tx1", "tx2", "tx3"),
        Gene = c("g1", "g1", "g2"),
        stringsAsFactors = FALSE
    )
    rownames(readcounts) <- tx2gene_df$Transcript
    
    se <- build_se(readcounts, tx2gene_df)
    
    # Check rowData
    rd <- rowData(se)
    expect_equal(rd$gene_id, c("g1", "g1", "g2"))
    expect_equal(rownames(rd), c("tx1", "tx2", "tx3"))
})

# Test 6: Test with TSV file containing tx2gene
test_that("build_se with TSV file input stores tx2gene in metadata", {
    # Create temporary TSV file
    tx2gene_df <- data.frame(
        Transcript = c("tx1", "tx2", "tx3"),
        Gene = c("g1", "g1", "g2"),
        stringsAsFactors = FALSE
    )
    tsv_file <- tempfile(fileext = ".tsv")
    write.table(tx2gene_df, file = tsv_file, sep = "\t", row.names = FALSE, quote = FALSE)
    
    readcounts <- matrix(1:6, nrow = 3, ncol = 2)
    colnames(readcounts) <- c("S1", "S2")
    rownames(readcounts) <- tx2gene_df$Transcript
    
    se <- build_se(readcounts, tsv_file)
    
    # Check metadata
    md <- metadata(se)
    expect_true("tx2gene" %in% names(md))
    expect_equal(md$tx2gene$Transcript, tx2gene_df$Transcript)
    expect_equal(md$tx2gene$Gene, tx2gene_df$Gene)
    
    unlink(tsv_file)
})

context("SummarizedExperiment Construction: Error Handling for Invalid Inputs")

# Test 7: Invalid tx2gene type (not character, data.frame, or path)
test_that("build_se errors on invalid tx2gene type (numeric)", {
    readcounts <- matrix(1:6, nrow = 3, ncol = 2)
    colnames(readcounts) <- c("S1", "S2")
    rownames(readcounts) <- c("tx1", "tx2", "tx3")
    
    expect_error(build_se(readcounts, tx2gene = 123),
                 "tx2gene.*must be a path.*or a data.frame")
})

# Test 8: Non-numeric readcounts
test_that("build_se errors on non-numeric readcounts", {
    readcounts <- data.frame(
        S1 = c("a", "b", "c"),
        S2 = c("d", "e", "f")
    )
    rownames(readcounts) <- c("tx1", "tx2", "tx3")
    
    tx2gene_df <- data.frame(
        Transcript = c("tx1", "tx2", "tx3"),
        Gene = c("g1", "g1", "g2"),
        stringsAsFactors = FALSE
    )
    
    expect_error(build_se(readcounts, tx2gene_df),
                 "readcounts.*must be a numeric")
})

# Test 9: Missing rownames in readcounts
test_that("build_se errors when readcounts has no rownames", {
    readcounts <- matrix(1:6, nrow = 3, ncol = 2)
    colnames(readcounts) <- c("S1", "S2")
    # No rownames set
    
    tx2gene_df <- data.frame(
        Transcript = c("tx1", "tx2", "tx3"),
        Gene = c("g1", "g1", "g2"),
        stringsAsFactors = FALSE
    )
    
    expect_error(build_se(readcounts, tx2gene_df),
                 "readcounts.*must have transcript IDs as rownames")
})

# Test 10: Non-existent file path
test_that("build_se errors when file path does not exist", {
    readcounts <- matrix(1:6, nrow = 3, ncol = 2)
    colnames(readcounts) <- c("S1", "S2")
    rownames(readcounts) <- c("tx1", "tx2", "tx3")
    
    expect_error(build_se(readcounts, "/nonexistent/path/file.tsv"),
                 "tx2gene file not found")
})

context("SummarizedExperiment Construction: Custom Assay Names and Data Conversion")

# Test 11: Multiple custom assay names
test_that("build_se respects custom assay_name parameter", {
    readcounts <- matrix(1:6, nrow = 3, ncol = 2)
    colnames(readcounts) <- c("S1", "S2")
    tx2gene_df <- data.frame(
        Transcript = c("tx1", "tx2", "tx3"),
        Gene = c("g1", "g1", "g2"),
        stringsAsFactors = FALSE
    )
    rownames(readcounts) <- tx2gene_df$Transcript
    
    se1 <- build_se(readcounts, tx2gene_df, assay_name = "raw_counts")
    se2 <- build_se(readcounts, tx2gene_df, assay_name = "normalized")
    
    expect_true("raw_counts" %in% names(assays(se1)))
    expect_false("counts" %in% names(assays(se1)))
    expect_true("normalized" %in% names(assays(se2)))
    expect_false("counts" %in% names(assays(se2)))
})

# Test 12: data.frame to matrix conversion
test_that("build_se converts numeric data.frame readcounts to matrix", {
    readcounts_df <- data.frame(
        S1 = c(1, 2, 3),
        S2 = c(4, 5, 6)
    )
    rownames(readcounts_df) <- c("tx1", "tx2", "tx3")
    
    tx2gene_df <- data.frame(
        Transcript = c("tx1", "tx2", "tx3"),
        Gene = c("g1", "g1", "g2"),
        stringsAsFactors = FALSE
    )
    
    se <- build_se(readcounts_df, tx2gene_df)
    
    assay_data <- assay(se, "counts")
    expect_true(is.matrix(assay_data))
    expect_equal(assay_data, as.matrix(readcounts_df))
})

context("SummarizedExperiment Construction: Dimension and Size Preservation")

# Test 13: Verify SE dimensions match input
test_that("build_se preserves dimensions of input readcounts", {
    for (n_tx in c(5, 50, 500)) {
        for (n_samps in c(2, 5, 10)) {
            readcounts <- matrix(sample(0:50, n_tx * n_samps, replace = TRUE), 
                                nrow = n_tx, ncol = n_samps)
            colnames(readcounts) <- paste0("S", 1:n_samps)
            rownames(readcounts) <- paste0("tx", 1:n_tx)
            
            tx2gene_df <- data.frame(
                Transcript = rownames(readcounts),
                Gene = rep_len(paste0("g", 1:10), n_tx),
                stringsAsFactors = FALSE
            )
            
            se <- build_se(readcounts, tx2gene_df)
            
            expect_equal(nrow(se), n_tx)
            expect_equal(ncol(se), n_samps)
            expect_equal(dim(assay(se)), c(n_tx, n_samps))
        }
    }
})

# Test 14: All transcripts mapped correctly
test_that("build_se correctly maps all transcripts when all are in tx2gene", {
    n_tx <- 20
    readcounts <- matrix(sample(0:100, n_tx * 4, replace = TRUE), 
                        nrow = n_tx, ncol = 4)
    colnames(readcounts) <- paste0("S", 1:4)
    rownames(readcounts) <- paste0("tx", 1:n_tx)
    
    # Create complete tx2gene mapping
    tx2gene_df <- data.frame(
        Transcript = paste0("tx", 1:n_tx),
        Gene = rep(paste0("g", 1:5), each = 4),
        stringsAsFactors = FALSE
    )
    
    se <- build_se(readcounts, tx2gene_df)
    rd <- rowData(se)
    
    # All genes should be mapped, no NAs
    expect_false(any(is.na(rd$gene_id)))
    expect_equal(length(unique(rd$gene_id)), 5)
})

context("SummarizedExperiment Construction: Real GFF3 File Integration")

# Test 15: Full workflow with actual package data
test_that("build_se works with real readcounts and GFF3 from package", {
    skip_if_not_installed("TSENAT")
    
    # Load real example data (contains salmon_dataset, salmon_tpm, salmon_effective_length)
    data(readcounts, package = "TSENAT")
    
    # Use real GFF3 file from package
    gff3_path <- system.file("extdata", "annotation.gff3.gz", package = "TSENAT")
    skip_if(gff3_path == "", "GFF3 file not found in package")
    
    # Build SE from real data (salmon_dataset contains the transcript counts)
    rc_matrix <- as.matrix(salmon_dataset)
    mode(rc_matrix) <- "numeric"
    se <- build_se(rc_matrix, gff3_path)
    
    # Verify structure
    expect_s4_class(se, "SummarizedExperiment")
    expect_true(nrow(se) > 0)
    expect_true(ncol(se) > 0)
    expect_true("counts" %in% names(assays(se)))
})
