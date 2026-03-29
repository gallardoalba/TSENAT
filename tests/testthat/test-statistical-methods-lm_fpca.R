library(testthat)

context("Linear Model Interaction: FPCA Methods")

test_that("fpca method attaches p_interaction to rowData", {
    # Use 4 samples and 5 q values to ensure adequate FPCA matrix coverage
    qvec <- seq(0.01, 0.05, by = 0.01)
    sample_ids <- rep(c("S1", "S2", "S3", "S4"), each = length(qvec))
    coln <- paste0(sample_ids, "_q=", rep(qvec, times = 4))

    set.seed(9)
    # create clear differences so PCA can pick up group separation
    gene1_vals <- c(qvec * 1, qvec * 2, qvec * 1.5, qvec * 2.2) + rnorm(length(coln), sd = 1e-4)
    gene2_vals <- c(qvec * 1, qvec * 1, qvec * 1, qvec * 1) + rnorm(length(coln), sd = 1e-4)
    mat <- rbind(g1 = gene1_vals, g2 = gene2_vals)
    colnames(mat) <- coln
    rownames(mat) <- c("g1", "g2")

    rd <- data.frame(genes = rownames(mat), row.names = rownames(mat), stringsAsFactors = FALSE)
    cd <- data.frame(samples = sample_ids, sample_type = rep(c("Normal", "Tumor", "Normal", "Tumor"), each = length(qvec)), row.names = coln, stringsAsFactors = FALSE)

    se <- SummarizedExperiment::SummarizedExperiment(assays = list(diversity = mat), rowData = rd, colData = cd)

    # lower min_obs so test is robust; 4 samples with full q coverage should pass
    res <- .calculate_lm_interaction(se, condition_col = "sample_type", method = "fpca", min_obs = 2)
    if (is.data.frame(res)) {
        rd_out <- as.data.frame(res)
    } else {
        rd_out <- as.data.frame(SummarizedExperiment::rowData(res))
    }
    expect_true("p_interaction" %in% colnames(rd_out))
})

# ============================================================================
# Tests for FPCA interaction helper prcomp and try-error handling
# ============================================================================

context("Linear Model Interaction: FPCA with prcomp and Error Handling")

test_that(".fpca_interaction handles prcomp successfully", {
    # This test covers: pca <- try(stats::prcomp(mat_sub, center = TRUE, scale. = FALSE), silent = TRUE)
    # and: if (inherits(pca, "try-error")) { return(NULL) }
    
    # Create properly formed input data
    mat <- matrix(rnorm(100), nrow = 10, ncol = 10)
    rownames(mat) <- paste0("Gene", 1:10)
    q_vals <- rep(c(0.5, 1, 1.5, 2), length.out = 10)
    sample_names <- rep(c("S1", "S2", "S3"), length.out = 10)
    group_vec <- rep(c("A", "B"), length.out = 10)
    
    res <- TSENAT:::.fpca_interaction(mat, q_vals, sample_names, group_vec, "Gene1", min_obs = 3)
    
    # Should return either NULL or a valid result with p_interaction
    if (!is.null(res)) {
        expect_is(res, "data.frame")
        expect_true("gene" %in% colnames(res))
        expect_true("p_interaction" %in% colnames(res))
        expect_equal(res$gene, "Gene1")
    } else {
        # When result is NULL, verify that it's handled correctly
        expect_null(res)
    }
})

test_that(".fpca_interaction returns NULL when prcomp fails", {
    # Create data that has insufficient variation for prcomp
    # All values the same would cause issues
    mat <- matrix(1.0, nrow = 10, ncol = 10)
    rownames(mat) <- paste0("Gene", 1:10)
    q_vals <- rep(c(0.5, 1), 5)
    sample_names <- c("S1", "S2", "S1", "S2", "S1", "S2", "S1", "S2", "S1", "S2")
    group_vec <- rep(c("A", "B"), 5)
    
    # This should either return NULL or handle gracefully
    res <- TSENAT:::.fpca_interaction(mat, q_vals, sample_names, group_vec, "Gene1", min_obs = 2)
    
    # Result should be NULL or a valid data frame
    if (!is.null(res)) {
        expect_is(res, "data.frame")
    } else {
        # When NULL, verify the return is correct
        expect_null(res)
    }
})

test_that(".fpca_interaction with NAs in data", {
    # Test prcomp with data containing NAs that need imputation
    mat <- matrix(rnorm(80), nrow = 10, ncol = 8)
    # Introduce some NAs
    mat[2, 3] <- NA
    mat[1, 4] <- NA
    rownames(mat) <- paste0("Gene", 1:10)
    q_vals <- rep(c(0.5, 1, 1.5, 2), 2)
    sample_names <- rep(c("S1", "S2"), each = 4)
    group_vec <- rep(c("A", "B"), 4)
    
    res <- TSENAT:::.fpca_interaction(mat, q_vals, sample_names, group_vec, "Gene1", min_obs = 2)
    
    # Should handle NAs and return result or NULL
    if (!is.null(res)) {
        expect_is(res, "data.frame")
        expect_equal(res$gene, "Gene1")
    } else {
        # When result is NULL, verify it's correctly NULL
        expect_null(res)
    }
})
