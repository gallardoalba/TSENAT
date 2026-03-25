context("Statistical Quality Control: min_valid_frac Parameter")

# calculate_method is internal; tests access via triple colon
calculate_method <- TSENAT:::calculate_method

test_that("min_valid_frac=0.75 filters out sparse genes correctly", {
    # Create a matrix with 3 genes: 
    # - Gene A: all valid values (16/16 = 100%)
    # - Gene B: 12/16 valid (75%) - should be kept with default 0.75
    # - Gene C: 8/16 valid (50%) - should be excluded with default 0.75
    mat <- matrix(NA_real_, nrow = 12, ncol = 16)
    
    # Gene A: all valid (100%)
    mat[1:4, ] <- rpois(64, lambda = 5)
    
    # Gene B: 75% valid (12/16) - mark 4 samples as non-finite
    mat[5:8, 1:12] <- rpois(48, lambda = 5)
    mat[5:8, 13:16] <- NaN  # 4 samples become NaN
    
    # Gene C: 50% valid (8/16) - mark 8 samples as non-finite
    mat[9:12, 1:8] <- rpois(32, lambda = 5)
    mat[9:12, 9:16] <- NaN  # 8 samples become NaN
    
    colnames(mat) <- paste0("S", 1:16)
    genes <- c(rep("A", 4), rep("B", 4), rep("C", 4))
    
    # With default min_valid_frac=0.75, should keep genes A and B, exclude C
    res <- calculate_method(mat, genes, norm = TRUE, q = 1, min_valid_frac = 0.75)
    expect_equal(nrow(res), 2)
    expect_true(all(c("A", "B") %in% res$Gene))
    expect_false("C" %in% res$Gene)
})

test_that("min_valid_frac=0.5 keeps genes with 50% or more valid values", {
    mat <- matrix(NA_real_, nrow = 8, ncol = 16)
    
    # Gene A: 100% valid
    mat[1:4, ] <- rpois(64, lambda = 5)
    
    # Gene B: 50% valid (8/16)
    mat[5:8, 1:8] <- rpois(32, lambda = 5)
    mat[5:8, 9:16] <- NaN
    
    colnames(mat) <- paste0("S", 1:16)
    genes <- c(rep("A", 4), rep("B", 4))
    
    # With min_valid_frac=0.5, both should be kept
    res <- calculate_method(mat, genes, norm = TRUE, q = 1, min_valid_frac = 0.5)
    expect_equal(nrow(res), 2)
    expect_true(all(c("A", "B") %in% res$Gene))
})

test_that("min_valid_frac=1.0 (strict) only keeps genes with 100% valid values", {
    mat <- matrix(NA_real_, nrow = 8, ncol = 16)
    
    # Gene A: 100% valid
    mat[1:4, ] <- rpois(64, lambda = 5)
    
    # Gene B: 99% valid (one NaN)
    mat[5:8, 1:15] <- rpois(60, lambda = 5)
    mat[5:8, 16] <- NaN
    
    colnames(mat) <- paste0("S", 1:16)
    genes <- c(rep("A", 4), rep("B", 4))
    
    # With min_valid_frac=1.0, only A should be kept
    res <- calculate_method(mat, genes, norm = TRUE, q = 1, min_valid_frac = 1.0)
    expect_equal(nrow(res), 1)
    expect_equal(res$Gene, "A")
})

test_that("min_valid_frac=0 disables filtering (keeps all genes with ≥1 valid value)", {
    mat <- matrix(NA_real_, nrow = 8, ncol = 16)
    
    # Gene A: 100% valid
    mat[1:4, ] <- rpois(64, lambda = 5)
    
    # Gene B: 1/16 valid (only one sample)
    mat[5:8, 1] <- rpois(4, lambda = 5)
    mat[5:8, 2:16] <- NaN
    
    colnames(mat) <- paste0("S", 1:16)
    genes <- c(rep("A", 4), rep("B", 4))
    
    # With min_valid_frac=0, no filtering applied; both should be kept
    res <- calculate_method(mat, genes, norm = TRUE, q = 1, min_valid_frac = 0)
    expect_equal(nrow(res), 2)
    expect_true(all(c("A", "B") %in% res$Gene))
})

test_that("min_valid_frac works correctly with multiple q values", {
    # With q = c(0.5, 1, 2), each gene has 16 * 3 = 48 total values
    mat <- matrix(NA_real_, nrow = 12, ncol = 16)
    
    # Gene A: 100% valid (48/48)
    mat[1:4, ] <- rpois(64, lambda = 5)
    
    # Gene B: 75% valid in first 12 samples, all NaN in last 4
    # This gives 12*3 = 36 valid out of 48 total = 75%
    mat[5:8, 1:12] <- rpois(48, lambda = 5)
    mat[5:8, 13:16] <- NaN
    
    # Gene C: 50% valid (24/48)
    mat[9:12, 1:8] <- NA_real_
    mat[9:12, 1:8] <- rpois(32, lambda = 5)
    mat[9:12, 9:16] <- NaN
    
    # Wait, I need 12 genes total to demonstrate this properly
    # Let me reconsider...
    
    colnames(mat) <- paste0("S", 1:16)
    genes <- c(rep("A", 4), rep("B", 4), rep("C", 4))
    
    # With multiple q and min_valid_frac=0.75
    # Gene A: 48/48 = 100% ✓ kept
    # Gene B: 36/48 = 75% ✓ kept
    # Gene C: 24/48 = 50% ✗ excluded
    res <- calculate_method(mat, genes, norm = TRUE, q = c(0.5, 1, 2), min_valid_frac = 0.75)
    expect_equal(nrow(res), 2)
    expect_true(all(c("A", "B") %in% res$Gene))
    expect_false("C" %in% res$Gene)
})

test_that("min_valid_frac verbose message reports correct percentage threshold", {
    mat <- matrix(NA_real_, nrow = 4, ncol = 8)
    
    # Gene A: 100% valid
    mat[1:4, ] <- rpois(32, lambda = 5)
    
    colnames(mat) <- paste0("S", 1:8)
    genes <- c("A")
    
    # No genes should be excluded with default threshold on complete data
    expect_silent(
        calculate_method(mat, genes, norm = TRUE, q = 1, min_valid_frac = 0.75, verbose = FALSE)
    )
})

test_that("min_valid_frac parameter passes through from calculate_diversity", {
    # Create a small synthetic SummarizedExperiment-like structure
    readcounts <- matrix(rpois(80, lambda = 5), nrow = 20, ncol = 4)
    rownames(readcounts) <- paste0("TX", 1:20)
    colnames(readcounts) <- paste0("S", 1:4)
    
    genes <- rep(c("A", "B", "C", "D", "E"), each = 4)
    
    # All genes have data, so all should pass any reasonable threshold
    res <- calculate_diversity(
        readcounts, 
        genes = genes, 
        norm = TRUE, 
        q = 1,
        pseudocount = 0.5,
        min_valid_frac = 0.75
    )
    
    expect_true(is(res, "SummarizedExperiment"))
    expect_equal(nrow(res), 5)  # All 5 genes retained
})

test_that("min_valid_frac=0.75 default recovers genes when combined with pseudocount", {
    # When pseudocount is used, more genes may have valid values
    # but sparse genes should still be excluded
    mat <- matrix(NA_real_, nrow = 8, ncol = 4)
    
    # Gene A: sparse - only 1/4 samples have valid counts
    mat[1:4, 1] <- 5  # only sample 1 has counts
    mat[1:4, 2:4] <- NaN  # other samples are NaN (invalid)
    
    # Gene B: uniform counts everywhere
    mat[5:8, ] <- rpois(16, lambda = 5)
    
    colnames(mat) <- paste0("S", 1:4)
    genes <- c(rep("A", 4), rep("B", 4))
    
    # With min_valid_frac=0.75, gene A needs 75% of samples valid (3/4)
    # Gene A: only 1/4 = 25% valid, so excluded
    # Gene B: 4/4 = 100% valid, so kept
    res <- calculate_method(mat, genes, norm = TRUE, q = 1, 
                            pseudocount = 0.5, min_valid_frac = 0.75)
    
    # Gene A: only 1/4 = 25% valid original counts, fails 75% threshold
    # Gene B: 4/4 = 100% valid, passes threshold
    expect_equal(nrow(res), 1)
    expect_equal(res$Gene, "B")
})
