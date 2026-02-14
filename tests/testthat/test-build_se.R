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
    expect_error(build_se(rc, 12345), "'tx2gene' must be a path \\(TSV or GFF3\\) or a data.frame")
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

# ============================================================================
# GFF3 Format Tests
# ============================================================================

context("build_se GFF3 file format support")

test_that("build_se extracts tx2gene from actual GFF3.gz file in inst/extdata", {
    # Use the reference GFF3.gz test file packaged with TSENAT
    gff3_gz_path <- system.file("extdata", "gencode_subset_test.gff3.gz", package = "TSENAT")
    
    skip_if(gff3_gz_path == "", "gencode_subset_test.gff3.gz not found in inst/extdata")
    
    # Create readcounts with transcript IDs from the GFF3 file
    tc_ids <- c("ENST00000001", "ENST00000002", "ENST00000003", "ENST00000004", 
                "ENST00000005", "ENST00000006", "ENST00000007", "ENST00000008")
    rc <- matrix(sample(1:100, length(tc_ids) * 3, replace = TRUE), 
                 nrow = length(tc_ids), ncol = 3)
    colnames(rc) <- c("S1", "S2", "S3")
    rownames(rc) <- tc_ids
    
    se <- build_se(rc, gff3_gz_path)
    
    expect_s4_class(se, "SummarizedExperiment")
    expect_equal(length(SummarizedExperiment::rowData(se)$genes), 8)
    expect_equal(SummarizedExperiment::rowData(se)$genes[1], "ENSG00000101456")  # MXRA8
    expect_equal(SummarizedExperiment::rowData(se)$genes[5], "ENSG00000102458")  # C1orf86
    expect_equal(SummarizedExperiment::rowData(se)$genes[7], "ENSG00000103259")  # PDPN
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
    expect_equal(SummarizedExperiment::rowData(se)$genes, c("ENSG00000001", "ENSG00000001", "ENSG00000002"))
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
    expect_equal(SummarizedExperiment::rowData(se)$genes, c("ENSG00000001", "ENSG00000001"))

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

    expect_error(build_se(rc, gff3_file), "No transcript-gene mappings found")

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

    expect_equal(length(SummarizedExperiment::rowData(se)$genes), 3)
    expect_equal(SummarizedExperiment::rowData(se)$genes[1], "ENSG00000001")
    expect_equal(SummarizedExperiment::rowData(se)$genes[3], "ENSG00000002")

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
    expect_equal(SummarizedExperiment::rowData(se)$genes, c("ENSG00000001", "ENSG00000002"))

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

context("build_se GFF3 malformed file handling")

test_that("build_se errors on GFF3 with missing columns", {
    # GFF3 with only 8 columns instead of 9
    gff3_content <- "##gff-version 3
chr1\tGENCODE\ttranscript\t1000\t2000\t.\t+\t."

    gff3_file <- tempfile(fileext = ".gff3")
    writeLines(gff3_content, gff3_file)

    rc <- matrix(1:2, nrow = 1)
    rownames(rc) <- "tx1"

    expect_error(build_se(rc, gff3_file), "No transcript-gene mappings found")

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

    expect_error(build_se(rc, gff3_file), "No transcript-gene mappings found")

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
    expect_error(build_se(rc, gff3_file), "No transcript-gene mappings found")

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

    expect_error(build_se(rc, gff3_file), "No transcript-gene mappings found")

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
    expect_equal(SummarizedExperiment::rowData(se)$genes, c("ENSG00000001", "ENSG00000001"))

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
    expect_equal(SummarizedExperiment::rowData(se)$genes, "ENSG00000001")

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
    expect_equal(SummarizedExperiment::rowData(se)$genes, "ENSG00000001")

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

    expect_error(build_se(rc, gff3_file), "No transcript-gene mappings found")

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
    expect_equal(length(SummarizedExperiment::rowData(se)$genes), 2)

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

    expect_error(build_se(rc, gff3_file), "Some transcript IDs in readcounts were not found")

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
    expect_equal(SummarizedExperiment::rowData(se)$genes, "ENSG00000001")

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
    expect_equal(SummarizedExperiment::rowData(se)$genes, "ENSG00000001")

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
    expect_equal(length(SummarizedExperiment::rowData(se)$genes), 3)
    expect_equal(SummarizedExperiment::rowData(se)$genes, c("ENSG00000001", "ENSG00000002", "ENSG00000003"))

    unlink(gff3_file)
})
