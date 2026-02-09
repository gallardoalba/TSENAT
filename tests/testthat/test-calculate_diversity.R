context("Main diversity calculation")

test_that(
    "calculate_diversity supports q as a vector and returns correct metadata",
    {
        x <- matrix(c(0, 0, 5, 4, 1, 2, 2, 2, 2, 2), ncol = 2)
        colnames(x) <- c("Sample1", "Sample2")
        gene <- c("Gene1", "Gene1", "Gene1", "Gene1", "Gene1")
        qvec <- c(1.1, 1.5, 2)
        result <- calculate_diversity(x, gene, norm = TRUE, q = qvec)
        # Assay should have columns for each q and sample
        assay_names <- colnames(SummarizedExperiment::assay(result))
        for (qi in qvec) {
            expect_true(any(grepl(paste0("q=", qi), assay_names)))
        }
        # Metadata should contain the q vector
        expect_equal(S4Vectors::metadata(result)$q, qvec)
    }
)
test_that("calculate_diversity passes q parameter for Tsallis entropy", {
    x <- matrix(c(0, 0, 5, 4, 1, 2, 2, 2, 2, 2), ncol = 2)
    colnames(x) <- c("Sample1", "Sample2")
    gene <- c("Gene1", "Gene1", "Gene1", "Gene1", "Gene1")
    # Calculate with q = 2
    result_q2 <- calculate_diversity(x, gene, norm = TRUE, q = 2)
    # Calculate with q = 1.5
    result_q15 <- calculate_diversity(x, gene, norm = TRUE, q = 1.5)
    # Extract values
    val_q2 <- SummarizedExperiment::assay(result_q2)[1, 1]
    val_q15 <- SummarizedExperiment::assay(result_q15)[1, 1]
    expect_true(is.numeric(val_q2))
    expect_true(is.numeric(val_q15))
    expect_false(is.na(val_q2))
    expect_false(is.na(val_q15))
    expect_false(abs(val_q2 - val_q15) < 1e-8)
})


test_that("calculate_diversity supports data.frame input", {
  x_df <- as.data.frame(matrix(c(1,2,3,4,5,6), ncol=2))
  colnames(x_df) <- c("A","B")
  genes <- c("g1","g1","g2")
  res <- calculate_diversity(x_df, genes, q = 1)
  expect_s4_class(res, "SummarizedExperiment")
  expect_equal(ncol(res), ncol(as.matrix(x_df)))
})

test_that("calculate_diversity handles tximport-style list with tpm flag", {
  counts <- matrix(c(10,0,5,2,3,1), ncol = 2)
    colnames(counts) <- c("S1", "S2")
  abundance <- counts / rowSums(counts) * 1e6
    colnames(abundance) <- colnames(counts)
  txlist <- list(counts = counts, abundance = abundance, length = NULL, countsFromAbundance = NULL)
  genes <- c("g1","g1","g2")
  # default uses counts
  res_counts <- calculate_diversity(txlist, genes, q = 1, tpm = FALSE)
  expect_s4_class(res_counts, "SummarizedExperiment")
  # using tpm should switch to abundance
  res_ab <- calculate_diversity(txlist, genes, q = 1, tpm = TRUE)
  expect_s4_class(res_ab, "SummarizedExperiment")
})

test_that("calculate_diversity handles object with class DGEList (list with class)", {
  counts <- matrix(c(1,2,3,4,5,6), ncol = 2)
    colnames(counts) <- c("S1", "S2")
  dgel <- list(counts = counts)
  class(dgel) <- "DGEList"
  genes <- c("g1","g1","g2")
  res <- calculate_diversity(dgel, genes, q = 2)
  expect_s4_class(res, "SummarizedExperiment")
})

test_that("calculate_diversity uses metadata$readcounts and tx2gene from a SummarizedExperiment", {
  rc <- matrix(c(5,5,0,1,2,3), ncol = 2)
  rownames(rc) <- c("tx1","tx2","tx3")
    colnames(rc) <- c("S1", "S2")
  tx2gene <- data.frame(Transcript = rownames(rc), Gen = c("g1","g1","g2"), stringsAsFactors = FALSE)
  se <- SummarizedExperiment(assays = S4Vectors::SimpleList(dummy = matrix(0, nrow = 3, ncol = 2)))
  S4Vectors::metadata(se)$readcounts <- rc
  S4Vectors::metadata(se)$tx2gene <- tx2gene
  res <- calculate_diversity(se, genes = NULL, q = 1)
  expect_s4_class(res, "SummarizedExperiment")
  expect_true(all(SummarizedExperiment::rowData(res)$genes %in% unique(tx2gene$Gen)))
})

test_that("calculate_diversity errors for invalid assay number in SummarizedExperiment", {
  se <- SummarizedExperiment(assays = S4Vectors::SimpleList(a = matrix(1, nrow=2, ncol=2)))
  genes <- c("g1","g2")
  expect_error(calculate_diversity(se, genes, assayno = 2), "Please provide a valid assay number")
})

# Defensive / error cases

test_that("calculate_diversity errors on non-numeric input", {
  x <- matrix(letters[1:6], ncol = 2)
  genes <- c("g1","g1","g2")
  expect_error(calculate_diversity(x, genes), "Input data  must be numeric")
    colnames(x) <- c("S1", "S2")
})

test_that("calculate_diversity errors on NA values", {
  x <- matrix(c(1, NA, 3, 4, 5, 6), ncol = 2)
  genes <- c("g1","g1","g2")
  expect_error(calculate_diversity(x, genes), "The data contains NA")
    colnames(x) <- c("S1", "S2")
})

test_that("calculate_diversity errors when genes length mismatches rows", {
  x <- matrix(1:6, ncol = 2)
  genes <- c("g1","g2") # wrong length
  expect_error(calculate_diversity(x, genes), "The number of rows is not equal to the given gene set")
    colnames(x) <- c("S1", "S2")
})

test_that("calculate_diversity errors for invalid q values (<=0)", {
  x <- matrix(1:6, ncol = 2)
  genes <- c("g1","g1","g2")
  expect_error(calculate_diversity(x, genes, q = 0), "Argument 'q' must be numeric and greater than 0")
})

test_that("calculate_diversity returns hill numbers when what='D'", {
  x <- matrix(c(1,2,3,4,5,6), ncol = 2)
  colnames(x) <- c("S1", "S2")
  genes <- c("g1","g1","g2")
  res <- calculate_diversity(x, genes, q = 2, what = "D")
  expect_true("hill" %in% names(SummarizedExperiment::assays(res)))
})

library(SummarizedExperiment)

test_that("tpm logical on non-list gives informative message", {
  x <- matrix(1:6, ncol = 2)
  colnames(x) <- c("S1", "S2")
  genes <- c("g1","g1","g2")
  expect_message(calculate_diversity(x, genes, tpm = TRUE, verbose = TRUE), "tpm as a logical argument is only interpreted")
})

test_that("tximport-style list missing counts errors", {
  txbad <- list(a = 1, b = 2, c = 3)
  genes <- c("g1","g1","g2")
  expect_error(calculate_diversity(txbad, genes), "cannot find any expression data")
})

test_that("SummarizedExperiment metadata tx2gene with non-standard columns is accepted", {
  rc <- matrix(c(5,1,0,2,3,4), ncol = 2)
  colnames(rc) <- c("S1","S2")
  rownames(rc) <- c("tx1","tx2","tx3")
  tx2gene <- data.frame(txid = rownames(rc), geneid = c("gA","gA","gB"), stringsAsFactors = FALSE)
  se <- SummarizedExperiment(assays = S4Vectors::SimpleList(dummy = matrix(0, nrow = 3, ncol = 2)))
  S4Vectors::metadata(se)$readcounts <- rc
  S4Vectors::metadata(se)$tx2gene <- tx2gene
  res <- calculate_diversity(se, genes = NULL, q = 1)
  expect_s4_class(res, "SummarizedExperiment")
  expect_true(all(SummarizedExperiment::rowData(res)$genes %in% unique(tx2gene$geneid)))
})

test_that("All-zero counts lead to empty result (no valid genes)", {
  x <- matrix(0, nrow = 4, ncol = 3)
  colnames(x) <- paste0("S", 1:3)
  genes <- c("g1","g1","g2","g3")
  res <- calculate_diversity(x, genes, q = 1)
  expect_s4_class(res, "SummarizedExperiment")
  expect_equal(nrow(res), 0)
})
