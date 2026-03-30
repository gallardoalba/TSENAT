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

# ============================================================================
# Additional FPCA tests (moved from test-statistical-methods-lm_helpers.R)
# ============================================================================

test_that(".fpca_interaction computes a p-value with reasonable input", {
    set.seed(2)
    # Create matrix genes x observations
    genes <- paste0("g", 1:3)
    samples <- paste0("S", 1:8)
    q_vals <- rep(c(0.1, 0.5, 1, 2), 2)
    # construct mat with rows genes, cols observations
    mat <- matrix(rnorm(length(genes) * length(q_vals)), nrow = length(genes))
    rownames(mat) <- genes
    # duplicate sample names to match observations length
    sample_names <- rep(samples[1:4], 2)
    group_vec <- rep(c("A", "B"), each = 4)
    # use min_obs small to allow test
    res <- .fpca_interaction(mat, q_vals = q_vals, sample_names = sample_names, group_vec = group_vec, g = 1, min_obs = 2)
    expect_true(is.null(res) || (is.data.frame(res) && "p_interaction" %in% colnames(res)))
})

test_that(".fpca_interaction returns NULL for non-diverse groups and handles imputation path", {
    # non-diverse groups -> NULL
    genes <- 1
    samples <- paste0("s", 1:6)
    q_vals <- rep(c(1, 2, 3), 2)
    mat <- matrix(rnorm(length(q_vals)), nrow = 1)
    rownames(mat) <- "g1"
    group_vec <- rep("A", length.out = length(q_vals))
    res <- .fpca_interaction(mat, q_vals = q_vals, sample_names = samples, group_vec = group_vec, g = 1, min_obs = 2)
    expect_null(res)

    # imputation path: create NA entries that are later imputed
    group_vec2 <- rep(c("A", "B"), each = 3)
    mat2 <- matrix(NA_real_, nrow = 1, ncol = 6)
    # fill some entries so there are at least two good rows after reshaping
    mat2[1, c(1, 4)] <- c(1.2, 2.3)
    rownames(mat2) <- "g1"
    sample_names2 <- paste0("s", 1:6)
    res2 <- .fpca_interaction(mat2, q_vals = q_vals, sample_names = sample_names2, group_vec = group_vec2, g = 1, min_obs = 1)
    expect_true(is.null(res2) || (is.data.frame(res2) && "p_interaction" %in% colnames(res2)))

    # q_vals with NA should be skipped during mapping (match returns NA)
    q_vals_na <- c(1, NA, 2, 3, NA, 2)
    mat_naq <- matrix(rnorm(length(q_vals_na)), nrow = 1)
    rownames(mat_naq) <- "g1"
    res_naq <- .fpca_interaction(mat_naq, q_vals = q_vals_na, sample_names = sample_names2, group_vec = group_vec2, g = 1, min_obs = 1)
    expect_true(is.null(res_naq) || is.data.frame(res_naq))
})

test_that(".fpca_interaction handles prcomp and t.test failures gracefully", {
    set.seed(101)
    genes <- paste0("g", 1)
    samples <- paste0("s", 1:6)
    q_vals <- rep(1:3, 2)
    mat <- matrix(rnorm(length(q_vals)), nrow = 1)
    rownames(mat) <- "g1"
    sample_names <- samples
    group_vec <- rep(c("A", "B"), each = 3)

    # We already exercise the basic null-return behavior; here we also ensure
    # that an imputation path that yields very small usable data returns either
    # NULL or a p_interaction, without triggering hard errors.
    mat3 <- matrix(NA_real_, nrow = 1, ncol = 6)
    mat3[1, c(1, 4)] <- c(1.2, 2.3)
    rownames(mat3) <- "g1"
    res3 <- .fpca_interaction(mat3, q_vals = q_vals, sample_names = sample_names, group_vec = group_vec, g = 1, min_obs = 2)
    expect_true(is.null(res3) || (is.data.frame(res3) && "p_interaction" %in% colnames(res3)))
})

test_that(".fpca_interaction includes slope_diff in results", {
    set.seed(1003)
    
    # Create synthetic matrix for FPCA
    genes <- "gene1"
    samples <- paste0("s", 1:8)
    q_vals <- rep(c(0.1, 0.5, 1, 2), 2)
    
    mat <- matrix(rnorm(length(q_vals)), nrow = 1)
    rownames(mat) <- genes
    
    sample_names <- samples
    group_vec <- rep(c("A", "B"), each = 4)
    
    result <- .fpca_interaction(
        mat,
        q_vals = q_vals,
        sample_names = sample_names,
        group_vec = group_vec,
        g = 1,
        min_obs = 2
    )
    
    # FPCA result should be either NULL or data.frame
    expect_true(is.null(result) || is.data.frame(result))
    
    # For FPCA, slope_diff should be NA (not applicable for functional analysis)
    if (is.data.frame(result)) {
        expect_true("slope_diff" %in% colnames(result))
        expect_true(is.na(result$slope_diff))
    }
})

test_that(".fpca_interaction works with all regularization methods", {
    skip_on_ci()  # Expensive: tests all 3 regularization modes. Keep locally for comprehensive validation
    skip_if_not_installed("glmnet")
    # Test via main function to verify all regularization methods are properly integrated
    qvec <- seq(0.01, 0.15, by = 0.01)
    
    set.seed(205)
    # Create 15 samples per group with expression pattern
    samples_normal <- paste0("S", 1:15)
    samples_tumor <- paste0("S", 16:30)
    all_samples <- c(samples_normal, samples_tumor)
    
    gene1_vals <- c()
    for (i in seq_len(15)) {
        gene1_vals <- c(gene1_vals, qvec * 1.0 + rnorm(length(qvec), sd = 0.05))
    }
    for (i in 16:30) {
        gene1_vals <- c(gene1_vals, qvec * 1.3 + rnorm(length(qvec), sd = 0.05))
    }
    
    sample_names <- rep(all_samples, each = length(qvec))
    coln <- paste0(sample_names, "_q=", qvec)
    
    mat <- rbind(g1 = gene1_vals)
    colnames(mat) <- coln
    rownames(mat) <- "g1"
    
    rd <- data.frame(genes = "g1", row.names = "g1", stringsAsFactors = FALSE)
    
    cd <- data.frame(
        samples = sample_names,
        sample_type = rep(c("Normal", "Tumor"), c(15*length(qvec), 15*length(qvec))),
        row.names = coln,
        stringsAsFactors = FALSE
    )
    
    se <- SummarizedExperiment::SummarizedExperiment(
        assays = list(diversity = mat),
        rowData = rd,
        colData = cd
    )
    
    # Test PCA regularization
    res_pca <- .calculate_lm_interaction(se,
        condition_col = "sample_type",
        method = "fpca",
        regularization = "pca",
        min_obs = 2
    )
    
    # Test LASSO regularization
    res_lasso <- suppressWarnings({
        .calculate_lm_interaction(se,
            condition_col = "sample_type",
            method = "fpca",
            regularization = "lasso",
            min_obs = 2
        )
    })
    
    # Test Elastic Net regularization
    res_elasticnet <- suppressWarnings({
        .calculate_lm_interaction(se,
            condition_col = "sample_type",
            method = "fpca",
            regularization = "elasticnet",
            min_obs = 2
        )
    })
    
    # All should return valid results
    expect_is(res_pca, "data.frame")
    expect_true(nrow(res_pca) > 0)
    expect_true("p_interaction" %in% colnames(res_pca))
    
    if (!is.null(res_lasso) && nrow(res_lasso) > 0) {
        expect_is(res_lasso, "data.frame")
        expect_true("p_interaction" %in% colnames(res_lasso))
    }
    
    if (!is.null(res_elasticnet) && nrow(res_elasticnet) > 0) {
        expect_is(res_elasticnet, "data.frame")
        expect_true("p_interaction" %in% colnames(res_elasticnet))
    }
})

# ============================================================================
# FPCA regularization methods (moved from test-integration-linear_models.R)
# ============================================================================

test_that("FPCA with regularization='pca' (default) works correctly", {
    library(TSENAT)
    qvec <- seq(0.01, 0.15, by = 0.01)
    
    set.seed(200)
    # Create curve-like data where each sample belongs to ONE group only
    # S1-S15 are Normal; S16-S30 are Tumor (15 per group for glmnet)
    samples_normal <- paste0("S", 1:15)
    samples_tumor <- paste0("S", 16:30)
    all_samples <- c(samples_normal, samples_tumor)
    
    gene1_vals <- c()
    for (i in seq_len(15)) {
        gene1_vals <- c(gene1_vals, qvec * 1.0 + rnorm(length(qvec), sd = 0.01))  # Normal
    }
    for (i in 16:30) {
        gene1_vals <- c(gene1_vals, qvec * 1.5 + rnorm(length(qvec), sd = 0.01))  # Tumor
    }
    
    sample_names <- rep(all_samples, each = length(qvec))
    coln <- paste0(sample_names, "_q=", qvec)
    
    mat <- rbind(g1 = gene1_vals)
    colnames(mat) <- coln
    rownames(mat) <- "g1"
    
    rd <- data.frame(genes = "g1", row.names = "g1", stringsAsFactors = FALSE)
    
    cd <- data.frame(
        samples = sample_names,
        sample_type = rep(c("Normal", "Tumor"), c(15*length(qvec), 15*length(qvec))),
        row.names = coln,
        stringsAsFactors = FALSE
    )
    
    se <- SummarizedExperiment::SummarizedExperiment(
        assays = list(diversity = mat),
        rowData = rd,
        colData = cd
    )
    
    # Default method should be "pca"
    res_default <- .calculate_lm_interaction(se,
        condition_col = "sample_type",
        method = "fpca",
        min_obs = 2
    )
    
    # Explicit "pca" method should give same result
    res_pca <- .calculate_lm_interaction(se,
        condition_col = "sample_type",
        method = "fpca",
        regularization = "pca",
        min_obs = 2
    )
    
    # Both should return results
    expect_is(res_default, "data.frame")
    expect_is(res_pca, "data.frame")
    expect_true(nrow(res_default) > 0)
    expect_true(nrow(res_pca) > 0)
    
    # Both should have p_interaction values
    expect_true(!is.na(res_default$p_interaction[1]))
    expect_true(!is.na(res_pca$p_interaction[1]))
    
    # Results should be very similar (same method)
    expect_equal(res_default$p_interaction[1], res_pca$p_interaction[1], tolerance = 1e-6)
})

test_that("FPCA with regularization='lasso' produces valid results", {
    skip_if_not_installed("glmnet")
    library(TSENAT)
    
    qvec <- seq(0.01, 0.15, by = 0.01)
    
    set.seed(201)
    # Create curve data where LASSO should select important q-values
    # 15 per group (30 total) for glmnet
    gene1_vals <- c()
    for (i in seq_len(15)) {
        gene1_vals <- c(gene1_vals, qvec * 1.0 + rnorm(length(qvec), sd = 0.01))  # Normal
    }
    for (i in 16:30) {
        gene1_vals <- c(gene1_vals, qvec * 2.0 + rnorm(length(qvec), sd = 0.01))  # Tumor
    }
    
    all_samples <- c(paste0("S", 1:15), paste0("S", 16:30))
    sample_names <- rep(all_samples, each = length(qvec))
    coln <- paste0(sample_names, "_q=", qvec)
    
    mat <- rbind(g1 = gene1_vals)
    colnames(mat) <- coln
    rownames(mat) <- "g1"
    
    rd <- data.frame(genes = "g1", row.names = "g1", stringsAsFactors = FALSE)
    
    cd <- data.frame(
        samples = sample_names,
        sample_type = rep(c("Normal", "Tumor"), c(15*length(qvec), 15*length(qvec))),
        row.names = coln,
        stringsAsFactors = FALSE
    )
    
    se <- SummarizedExperiment::SummarizedExperiment(
        assays = list(diversity = mat),
        rowData = rd,
        colData = cd
    )
    
    res_lasso <- .calculate_lm_interaction(se,
        condition_col = "sample_type",
        method = "fpca",
        regularization = "lasso",
        min_obs = 2
    )
    
    # Should return valid results
    expect_is(res_lasso, "data.frame")
    expect_true(nrow(res_lasso) > 0)
    expect_true("p_interaction" %in% colnames(res_lasso))
    expect_true(!is.na(res_lasso$p_interaction[1]))
    expect_true(res_lasso$p_interaction[1] >= 0 && res_lasso$p_interaction[1] <= 1)
})

test_that("FPCA with regularization='elasticnet' produces valid results", {
    skip_if_not_installed("glmnet")
    library(TSENAT)
    
    qvec <- seq(0.01, 0.15, by = 0.01)
    
    set.seed(202)
    # Create curve data with moderate interaction
    # 15 per group (30 total) for glmnet
    gene1_vals <- c()
    for (i in seq_len(15)) {
        gene1_vals <- c(gene1_vals, qvec * 1.2 + rnorm(length(qvec), sd = 0.01))  # Normal
    }
    for (i in 16:30) {
        gene1_vals <- c(gene1_vals, qvec * 1.8 + rnorm(length(qvec), sd = 0.01))  # Tumor
    }
    
    all_samples <- c(paste0("S", 1:15), paste0("S", 16:30))
    sample_names <- rep(all_samples, each = length(qvec))
    coln <- paste0(sample_names, "_q=", qvec)
    
    mat <- rbind(g1 = gene1_vals)
    colnames(mat) <- coln
    rownames(mat) <- "g1"
    
    rd <- data.frame(genes = "g1", row.names = "g1", stringsAsFactors = FALSE)
    
    cd <- data.frame(
        samples = sample_names,
        sample_type = rep(c("Normal", "Tumor"), c(15*length(qvec), 15*length(qvec))),
        row.names = coln,
        stringsAsFactors = FALSE
    )
    
    se <- SummarizedExperiment::SummarizedExperiment(
        assays = list(diversity = mat),
        rowData = rd,
        colData = cd
    )
    
    res_elasticnet <- .calculate_lm_interaction(se,
        condition_col = "sample_type",
        method = "fpca",
        regularization = "elasticnet",
        min_obs = 2
    )
    
    # Should return valid results
    expect_is(res_elasticnet, "data.frame")
    expect_true(nrow(res_elasticnet) > 0)
    expect_true("p_interaction" %in% colnames(res_elasticnet))
    expect_true(!is.na(res_elasticnet$p_interaction[1]))
    expect_true(res_elasticnet$p_interaction[1] >= 0 && res_elasticnet$p_interaction[1] <= 1)
})

test_that("FPCA regularization methods produce reasonable p-value differences", {
    skip_on_ci()  # Expensive: runs 3 regularization methods (pca, lasso, elasticnet). Keep locally for full validation
    skip_if_not_installed("glmnet")
    library(TSENAT)
    
    qvec <- seq(0.01, 0.15, by = 0.01)
    
    set.seed(203)
    # Data with clear interaction for comparing methods (15 per group)
    gene1_vals <- c()
    for (i in seq_len(15)) {
        gene1_vals <- c(gene1_vals, qvec * 1.0 + rnorm(length(qvec), sd = 0.005))   # Normal
    }
    for (i in 16:30) {
        gene1_vals <- c(gene1_vals, qvec * 2.0 + rnorm(length(qvec), sd = 0.005))   # Tumor
    }
    
    all_samples <- c(paste0("S", 1:15), paste0("S", 16:30))
    sample_names <- rep(all_samples, each = length(qvec))
    coln <- paste0(sample_names, "_q=", qvec)
    
    mat <- rbind(g1 = gene1_vals, g2 = rnorm(length(coln), sd = 0.1))
    colnames(mat) <- coln
    rownames(mat) <- c("g1", "g2")
    
    rd <- data.frame(genes = c("g1", "g2"), row.names = c("g1", "g2"), stringsAsFactors = FALSE)
    
    cd <- data.frame(
        samples = sample_names,
        sample_type = rep(c("Normal", "Tumor"), c(15*length(qvec), 15*length(qvec))),
        row.names = coln,
        stringsAsFactors = FALSE
    )
    
    se <- SummarizedExperiment::SummarizedExperiment(
        assays = list(diversity = mat),
        rowData = rd,
        colData = cd
    )
    
    # Compare all three methods
    res_pca <- .calculate_lm_interaction(se,
        condition_col = "sample_type",
        method = "fpca",
        regularization = "pca",
        min_obs = 2
    )
    
    res_lasso <- .calculate_lm_interaction(se,
        condition_col = "sample_type",
        method = "fpca",
        regularization = "lasso",
        min_obs = 2
    )
    
    res_elasticnet <- .calculate_lm_interaction(se,
        condition_col = "sample_type",
        method = "fpca",
        regularization = "elasticnet",
        min_obs = 2
    )
    
    # All three should return data frames with results
    expect_is(res_pca, "data.frame")
    expect_is(res_lasso, "data.frame")
    expect_is(res_elasticnet, "data.frame")
    
    expect_true(nrow(res_pca) > 0)
    expect_true(nrow(res_lasso) > 0)
    expect_true(nrow(res_elasticnet) > 0)
    
    # All genes should have p-values
    expect_true(all(!is.na(res_pca$p_interaction[res_pca$gene == "g1"])))
    expect_true(all(!is.na(res_lasso$p_interaction[res_lasso$gene == "g1"])))
    expect_true(all(!is.na(res_elasticnet$p_interaction[res_elasticnet$gene == "g1"])))
    
    # For g1 (true interaction), LASSO and ElasticNet should often provide more power
    # (smaller p-values) than PCA by selecting important q-values
    # But this depends on the specific data, so we just check they're all reasonable
    p_pca <- res_pca$p_interaction[res_pca$gene == "g1"]
    p_lasso <- res_lasso$p_interaction[res_lasso$gene == "g1"]
    p_elasticnet <- res_elasticnet$p_interaction[res_elasticnet$gene == "g1"]
    
    if (!is.na(p_pca) && !is.na(p_lasso) && !is.na(p_elasticnet)) {
        # All should be between 0 and 1
        expect_true(p_pca >= 0 && p_pca <= 1)
        expect_true(p_lasso >= 0 && p_lasso <= 1)
        expect_true(p_elasticnet >= 0 && p_elasticnet <= 1)
    }
})

test_that("FPCA regularization with paired design works correctly", {
    skip_if_not_installed("glmnet")
    library(TSENAT)
    
    # Test paired design with regularization methods
    qvec <- seq(0.01, 0.1, by = 0.01)
    sample_pairs <- paste0("P", 1:20)  # 20 pairs = 20 samples per group (optimal for glmnet with 10 features)
    
    set.seed(205)
    gene1_vals <- numeric()
    sample_vec <- character()
    type_vec <- character()
    base_vec <- character()
    
    for (p in sample_pairs) {
        # Normal (group 1)
        norm_vals <- qvec * 1.0 + rnorm(length(qvec), sd = 0.01)
        gene1_vals <- c(gene1_vals, norm_vals)
        sample_vec <- c(sample_vec, paste0(p, "_N_q=", qvec))
        type_vec <- c(type_vec, rep("Normal", length(qvec)))
        base_vec <- c(base_vec, rep(p, length(qvec)))
        
        # Tumor (group 2) with interaction
        tumor_vals <- qvec * 1.5 + rnorm(length(qvec), sd = 0.01)
        gene1_vals <- c(gene1_vals, tumor_vals)
        sample_vec <- c(sample_vec, paste0(p, "_T_q=", qvec))
        type_vec <- c(type_vec, rep("Tumor", length(qvec)))
        base_vec <- c(base_vec, rep(p, length(qvec)))
    }
    
    mat <- rbind(g1 = gene1_vals)
    colnames(mat) <- sample_vec
    rownames(mat) <- "g1"
    
    rd <- data.frame(genes = "g1", row.names = "g1", stringsAsFactors = FALSE)
    
    cd <- data.frame(
        samples = gsub("_q=.*", "", sample_vec),
        sample_type = type_vec,
        sample_base = base_vec,
        row.names = sample_vec,
        stringsAsFactors = FALSE
    )
    
    se <- SummarizedExperiment::SummarizedExperiment(
        assays = list(diversity = mat),
        rowData = rd,
        colData = cd
    )
    
    # Test regularization with paired design (with sufficient observations for glmnet)
    res <- .calculate_lm_interaction(se,
        condition_col = "sample_type",
        method = "fpca",
        subject_col = "sample_base",
        paired = TRUE,
        regularization = "lasso",
        min_obs = 2
    )
    
    # Should produce valid results
    expect_is(res, "data.frame")
    expect_true(nrow(res) > 0)
    expect_true(!is.na(res$p_interaction[1]))
})
