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
    x_df <- as.data.frame(matrix(c(1, 2, 3, 4, 5, 6), ncol = 2))
    colnames(x_df) <- c("A", "B")
    genes <- c("g1", "g1", "g2")
    res <- calculate_diversity(x_df, genes, q = 1)
    expect_s4_class(res, "SummarizedExperiment")
    expect_equal(ncol(res), ncol(as.matrix(x_df)))
})

test_that("calculate_diversity handles tximport-style list with tpm flag", {
    counts <- matrix(c(10, 0, 5, 2, 3, 1), ncol = 2)
    colnames(counts) <- c("S1", "S2")
    abundance <- counts / rowSums(counts) * 1e6
    colnames(abundance) <- colnames(counts)
    txlist <- list(counts = counts, abundance = abundance, length = NULL, countsFromAbundance = NULL)
    genes <- c("g1", "g1", "g2")
    # default uses counts
    res_counts <- calculate_diversity(txlist, genes, q = 1, tpm = FALSE)
    expect_s4_class(res_counts, "SummarizedExperiment")
    # using tpm should switch to abundance
    res_ab <- calculate_diversity(txlist, genes, q = 1, tpm = TRUE)
    expect_s4_class(res_ab, "SummarizedExperiment")
})

test_that("calculate_diversity handles object with class DGEList (list with class)", {
    counts <- matrix(c(1, 2, 3, 4, 5, 6), ncol = 2)
    colnames(counts) <- c("S1", "S2")
    dgel <- list(counts = counts)
    class(dgel) <- "DGEList"
    genes <- c("g1", "g1", "g2")
    res <- calculate_diversity(dgel, genes, q = 2)
    expect_s4_class(res, "SummarizedExperiment")
})

test_that("calculate_diversity uses metadata$readcounts and tx2gene from a SummarizedExperiment", {
    rc <- matrix(c(5, 5, 0, 1, 2, 3), ncol = 2)
    rownames(rc) <- c("tx1", "tx2", "tx3")
    colnames(rc) <- c("S1", "S2")
    tx2gene <- data.frame(Transcript = rownames(rc), Gen = c("g1", "g1", "g2"), stringsAsFactors = FALSE)
    se <- SummarizedExperiment(assays = S4Vectors::SimpleList(dummy = matrix(0, nrow = 3, ncol = 2)))
    S4Vectors::metadata(se)$readcounts <- rc
    S4Vectors::metadata(se)$tx2gene <- tx2gene
    res <- calculate_diversity(se, genes = NULL, q = 1)
    expect_s4_class(res, "SummarizedExperiment")
    expect_true(all(SummarizedExperiment::rowData(res)$genes %in% unique(tx2gene$Gen)))
})

test_that("calculate_diversity errors for invalid assay number in SummarizedExperiment", {
    se <- SummarizedExperiment(assays = S4Vectors::SimpleList(a = matrix(1, nrow = 2, ncol = 2)))
    genes <- c("g1", "g2")
    expect_error(calculate_diversity(se, genes, assayno = 2), "Please provide a valid assay number")
})

# Defensive / error cases

test_that("calculate_diversity errors on non-numeric input", {
    x <- matrix(letters[1:6], ncol = 2)
    genes <- c("g1", "g1", "g2")
    expect_error(calculate_diversity(x, genes), "Input data  must be numeric")
    colnames(x) <- c("S1", "S2")
})

test_that("calculate_diversity errors on NA values", {
    x <- matrix(c(1, NA, 3, 4, 5, 6), ncol = 2)
    genes <- c("g1", "g1", "g2")
    expect_error(calculate_diversity(x, genes), "The data contains NA")
    colnames(x) <- c("S1", "S2")
})

test_that("calculate_diversity errors when genes length mismatches rows", {
    x <- matrix(1:6, ncol = 2)
    genes <- c("g1", "g2") # wrong length
    expect_error(calculate_diversity(x, genes), "The number of rows is not equal to the given gene set")
    colnames(x) <- c("S1", "S2")
})

test_that("calculate_diversity errors for invalid q values (<=0)", {
    x <- matrix(1:6, ncol = 2)
    genes <- c("g1", "g1", "g2")
    expect_error(calculate_diversity(x, genes, q = 0), "Argument 'q' must be numeric and greater than 0")
})

test_that("calculate_diversity returns hill numbers when what='D'", {
    x <- matrix(c(1, 2, 3, 4, 5, 6), ncol = 2)
    colnames(x) <- c("S1", "S2")
    genes <- c("g1", "g1", "g2")
    res <- calculate_diversity(x, genes, q = 2, what = "D")
    expect_true("hill" %in% names(SummarizedExperiment::assays(res)))
})

library(SummarizedExperiment)

test_that("tpm logical on non-list gives informative message", {
    x <- matrix(1:6, ncol = 2)
    colnames(x) <- c("S1", "S2")
    genes <- c("g1", "g1", "g2")
    expect_message(calculate_diversity(x, genes, tpm = TRUE, verbose = TRUE), "tpm as a logical argument is only interpreted")
})

test_that("tximport-style list missing counts errors", {
    txbad <- list(a = 1, b = 2, c = 3)
    genes <- c("g1", "g1", "g2")
    expect_error(calculate_diversity(txbad, genes), "cannot find any expression data")
})

test_that("SummarizedExperiment metadata tx2gene with non-standard columns is accepted", {
    rc <- matrix(c(5, 1, 0, 2, 3, 4), ncol = 2)
    colnames(rc) <- c("S1", "S2")
    rownames(rc) <- c("tx1", "tx2", "tx3")
    tx2gene <- data.frame(txid = rownames(rc), geneid = c("gA", "gA", "gB"), stringsAsFactors = FALSE)
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
    genes <- c("g1", "g1", "g2", "g3")
    res <- calculate_diversity(x, genes, q = 1)
    expect_s4_class(res, "SummarizedExperiment")
    expect_equal(nrow(res), 0)
})


### Additional Tsallis/entropy tests (merged)

test_that("calculate_diversity returns correct Tsallis entropy for single q", {
    set.seed(123)
    x <- matrix(rpois(60, 10), ncol = 6)
    colnames(x) <- paste0("Sample", 1:6)
    gene <- c(rep("Gene1", 3), rep("Gene2", 2), rep("Gene3", 3), rep("Gene4", 2))
    q <- 2
    result <- calculate_diversity(x, gene, q = q)
    expect_s4_class(result, "SummarizedExperiment")
    expect_true("diversity" %in% names(SummarizedExperiment::assays(result)))
    expect_equal(nrow(result), length(unique(gene)))
    expect_equal(ncol(result), ncol(x))
    expect_true(!is.null(SummarizedExperiment::rowData(result)$genes))
})


## calculate_tsallis_entropy low-level tests

test_that("calculate_tsallis_entropy computes correct q=1 Shannon entropy", {
    counts <- c(100, 50, 25, 25)
    total <- sum(counts)
    p <- counts / total
    expected <- -sum(p[p > 0] * log(p[p > 0]))
    result <- calculate_tsallis_entropy(counts, q = 1)
    expect_true(abs(result - expected) < 0.5)
})

test_that("calculate_tsallis_entropy requires q > 0", {
    counts <- c(10, 20, 0, 15)
    expect_error(calculate_tsallis_entropy(counts, q = 0))
})

test_that("calculate_tsallis_entropy handles uniform distribution", {
    counts <- c(25, 25, 25, 25)
    result_q05 <- calculate_tsallis_entropy(counts, q = 0.5)
    expect_true(is.numeric(result_q05))
    expect_true(!is.na(result_q05))
    expect_true(result_q05 > 0)
})

test_that("calculate_tsallis_entropy returns 0 for single taxon", {
    counts <- c(100, 0, 0)
    result <- calculate_tsallis_entropy(counts, q = 1.5)
    expect_equal(result, 0, tolerance = 1e-6)
})

test_that("calculate_tsallis_entropy increases with diversity", {
    uniform <- c(50, 50, 50, 50)
    uneven <- c(100, 40, 5, 5)
    entropy_uniform <- calculate_tsallis_entropy(uniform, q = 1)
    entropy_uneven <- calculate_tsallis_entropy(uneven, q = 1)
    expect_true(entropy_uniform > entropy_uneven)
})

test_that("calculate_tsallis_entropy is invariant to scale", {
    counts1 <- c(10, 20, 30)
    counts2 <- c(100, 200, 300)
    result1 <- calculate_tsallis_entropy(counts1, q = 1)
    result2 <- calculate_tsallis_entropy(counts2, q = 1)
    expect_equal(result1, result2, tolerance = 1e-10)
})

test_that("calculate_tsallis_entropy handles different q values (q > 0)", {
    counts <- c(100, 50, 30, 20)
    q_values <- c(0.1, 0.5, 1, 2, 3)
    results <- sapply(q_values, function(q) {
        calculate_tsallis_entropy(counts, q = q)
    })
    expect_length(results, 5)
    expect_true(all(!is.na(results)))
    expect_true(all(results >= 0))
})

test_that("calculate_tsallis_entropy returns numeric scalar", {
    counts <- c(10, 20, 15, 5)
    result <- calculate_tsallis_entropy(counts, q = 1.2)
    expect_is(result, "numeric")
    expect_length(result, 1)
})

test_that("calculate_tsallis_entropy handles zero-sum and q=1 correctly", {
    x_uniform <- c(1, 1, 1)
    s_unif <- calculate_tsallis_entropy(x_uniform, q = c(0.5, 1, 2), norm = TRUE, what = "S")
    expect_equal(as.numeric(s_unif), rep(1, 3))
    x_zero <- c(0, 0, 0)
    s_zero <- calculate_tsallis_entropy(x_zero, q = c(0.5, 1, 2), norm = TRUE, what = "S")
    expect_true(all(is.na(s_zero)))
    x <- c(10, 5, 0)
    p <- x / sum(x)
    sh <- -sum(ifelse(p > 0, p * log(p), 0))
    expected_D1 <- exp(sh)
    D1 <- calculate_tsallis_entropy(x, q = 1, what = "D")
    expect_equal(as.numeric(D1), expected_D1)
})

test_that("calculate_diversity rejects non-positive q", {
    mat <- matrix(1, nrow = 3, ncol = 2)
    genes <- letters[1:3]
    expect_error(calculate_diversity(mat, genes = genes, q = 0), "q")
})

test_that("calculate_diversity returns correct Tsallis entropy for vector q", {
    set.seed(123)
    x <- matrix(rpois(60, 10), ncol = 6)
    colnames(x) <- paste0("Sample", 1:6)
    gene <- c(rep("Gene1", 3), rep("Gene2", 2), rep("Gene3", 3), rep("Gene4", 2))
    q <- c(1, 2)
    result <- calculate_diversity(x, gene, q = q)
    expect_s4_class(result, "SummarizedExperiment")
    expect_true("diversity" %in% names(SummarizedExperiment::assays(result)))
    expect_equal(nrow(result), length(unique(gene)))
    expect_equal(ncol(result), ncol(x) * length(q))
    expect_true(!is.null(SummarizedExperiment::rowData(result)$genes))
    expect_true(all(c(
        "samples",
        "q"
    ) %in% colnames(SummarizedExperiment::colData(result))))
    expect_equal(
        length(unique(SummarizedExperiment::colData(result)$q)),
        length(q)
    )
})
test_that("calculate_diversity sets se_assay_mat when called directly with matrix input", {
    # This test covers: if (!exists("se_assay_mat")) { se_assay_mat <- x }
    # When calculate_diversity is called directly (not from within another function),
    # se_assay_mat should not exist initially, so it gets assigned from input x
    # Use multiple transcripts per gene to avoid NaN from single-isoform normalization
    x <- matrix(c(1, 2, 3, 4, 5, 6, 7, 8), nrow = 4, ncol = 2)
    colnames(x) <- c("S1", "S2")
    genes <- c("g1", "g1", "g2", "g2")  # 2 transcripts per gene
    
    # Make sure se_assay_mat doesn't exist in parent environment
    # The function should create it internally from the input matrix
    result <- calculate_diversity(x, genes, q = 1.5)
    
    # Verify the result is correct
    expect_s4_class(result, "SummarizedExperiment")
    expect_equal(nrow(result), 2)  # 2 genes
    expect_equal(ncol(result), 2)  # 2 samples
    expect_true("diversity" %in% names(SummarizedExperiment::assays(result)))
})

test_that("calculate_diversity with numeric matrix directly uses input as se_assay_mat", {
    # Additional test to verify that direct matrix input is properly handled
    # when se_assay_mat doesn't pre-exist
    # Use at least 2 transcripts per gene to avoid NaN from single-isoform normalization
    set.seed(42)
    x <- matrix(rpois(24, lambda = 5), nrow = 6, ncol = 4)
    colnames(x) <- c("SA", "SB", "SC", "SD")
    genes <- c("g1", "g1", "g2", "g2", "g3", "g3")  # 2 transcripts per gene
    
    result <- calculate_diversity(x, genes, norm = TRUE, q = 2)
    
    expect_s4_class(result, "SummarizedExperiment")
    expect_equal(ncol(result), 4)
    expect_equal(rownames(SummarizedExperiment::assay(result)), unique(genes))
})

test_that("calculate_diversity handles matrix input with multiple q values correctly", {
    # Tests the se_assay_mat assignment path with multi-value q vector
    # Use multiple transcripts per gene to avoid NaN from single-isoform normalization
    x <- matrix(c(10, 20, 15, 5, 8, 12, 6, 9, 11, 7, 4, 13), nrow = 6, ncol = 2)
    colnames(x) <- c("Sample1", "Sample2")
    genes <- c("Gene1", "Gene1", "Gene2", "Gene2", "Gene3", "Gene3")
    
    q_vec <- c(0.5, 1, 1.5, 2)
    result <- calculate_diversity(x, genes, q = q_vec)
    
    # Should have columns for each combination of sample and q value
    expect_equal(ncol(result), length(colnames(x)) * length(q_vec))
    expect_equal(nrow(result), length(unique(genes)))
    
    # Verify metadata has correct q values
    expect_equal(S4Vectors::metadata(result)$q, q_vec)
})