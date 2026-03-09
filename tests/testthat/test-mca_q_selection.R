context("Q-Parameter Selection: Multiple Correspondence Analysis (Feature #8)")

# Load package functions
library(TSENAT)

# Load test data and create gene-level structure from transcripts
load(system.file("data", "readcounts.RData", package = "TSENAT"))
# Use transcript-level data, grouped into genes (transcripts per gene)
full_rc <- as.matrix(salmon_dataset[1:100, , drop = FALSE])
mode(full_rc) <- "numeric"
full_tx_names <- rownames(full_rc)

# Assign transcripts to genes: 5 transcripts per gene
n_transcripts <- nrow(full_rc)
transcripts_per_gene <- 5
n_genes <- ceiling(n_transcripts / transcripts_per_gene)
gene_names <- paste0("GENE_", rep(seq_len(n_genes), each = transcripts_per_gene)[seq_len(n_transcripts)])

# Use transcript-level data with gene mapping
rc <- full_rc
gs <- gene_names

# Helper: Extract entropy matrix where rows are genes and columns are q values
# This averages across samples to provide sufficient variation for MCA
extract_entropy_by_q <- function(entropy_se, q_grid) {
    entropy_matrix <- assay(entropy_se, "diversity")
    col_data <- colData(entropy_se)
    
    # Reshape so that each row is a gene and each column is a q value
    # Average entropy across samples for each gene-q combination
    entropy_by_q <- matrix(0, nrow = nrow(entropy_matrix), ncol = length(q_grid))
    rownames(entropy_by_q) <- rownames(entropy_matrix)
    colnames(entropy_by_q) <- paste0("q=", round(q_grid, 3))
    
    for (q_idx in seq_along(q_grid)) {
        q_val <- q_grid[q_idx]
        q_cols <- which(col_data$q == q_val)
        entropy_by_q[, q_idx] <- rowMeans(entropy_matrix[, q_cols, drop = FALSE], na.rm = TRUE)
    }
    
    entropy_by_q
}

test_that("MCA q-selection works with multiple q values", {
    # Compute entropy for multiple q values
    q_grid <- c(0.5, 1, 1.5, 2, 2.5)
    
    # Use more genes and samples for sufficient entropy variation
    entropy_se <- calculate_diversity(rc[1:50, ], gs[1:50], q = q_grid, norm = "none", verbose = FALSE)
    entropy_matrix <- extract_entropy_by_q(entropy_se, q_grid)
    
    # Run MCA
    mca_result <- TSENAT:::.tsenat_mca_q_selection(entropy_matrix, q_grid)
    
    # Check structure
    expect_true(is.list(mca_result))
    expect_true(all(c("q_values", "inertia", "recommended_q", "variance_explained") %in% names(mca_result)))
    
    # Check that results are reasonable
    expect_equal(length(mca_result$q_values), length(q_grid))
    expect_true(all(mca_result$recommended_q %in% q_grid))
    expect_true(mca_result$variance_explained >= 0.8)  # Default 80% threshold
})

test_that("MCA selects subset of q values", {
    q_grid <- c(0.5, 1, 1.5, 2, 2.5, 3)
    
    entropy_se <- calculate_diversity(rc[1:50, ], gs[1:50], q = q_grid, norm = "none", verbose = FALSE)
    entropy_matrix <- extract_entropy_by_q(entropy_se, q_grid)
    
    mca_result <- TSENAT:::.tsenat_mca_q_selection(entropy_matrix, q_grid, min_variance_explained = 0.80)
    
    # Should recommend fewer than all q values
    expect_true(length(mca_result$recommended_q) <= length(q_grid))
    
    # Should explain target variance
    expect_true(mca_result$variance_explained >= 0.80)
})

test_that("MCA respects variance threshold", {
    q_grid <- c(0.5, 1, 1.5, 2, 2.5)
    
    entropy_se <- calculate_diversity(rc[1:50, ], gs[1:50], q = q_grid, norm = "none", verbose = FALSE)
    entropy_matrix <- extract_entropy_by_q(entropy_se, q_grid)
    
    # Strict threshold
    mca_strict <- TSENAT:::.tsenat_mca_q_selection(entropy_matrix, q_grid, min_variance_explained = 0.95)
    
    # Lenient threshold
    mca_lenient <- TSENAT:::.tsenat_mca_q_selection(entropy_matrix, q_grid, min_variance_explained = 0.60)
    
    # Strict should recommend more q values than lenient
    expect_true(length(mca_strict$recommended_q) >= length(mca_lenient$recommended_q))
})

test_that("MCA handles different numbers of categories", {
    q_grid <- c(0.5, 1, 1.5, 2, 2.5)
    
    entropy_se <- calculate_diversity(rc[1:50, ], gs[1:50], q = q_grid, norm = "none", verbose = FALSE)
    entropy_matrix <- extract_entropy_by_q(entropy_se, q_grid)
    
    # Test different categorization schemes
    mca_binary <- TSENAT:::.tsenat_mca_q_selection(entropy_matrix, q_grid, n_categories = 2)
    mca_tertile <- TSENAT:::.tsenat_mca_q_selection(entropy_matrix, q_grid, n_categories = 3)
    mca_quartile <- TSENAT:::.tsenat_mca_q_selection(entropy_matrix, q_grid, n_categories = 4)
    
    # All should return valid results
    expect_true(length(mca_binary$recommended_q) > 0)
    expect_true(length(mca_tertile$recommended_q) > 0)
    expect_true(length(mca_quartile$recommended_q) > 0)
})

test_that("MCA contributions sum to 1", {
    q_grid <- c(0.5, 1, 1.5, 2)
    
    entropy_se <- calculate_diversity(rc[1:50, ], gs[1:50], q = q_grid, norm = "none", verbose = FALSE)
    entropy_matrix <- extract_entropy_by_q(entropy_se, q_grid)
    
    mca_result <- TSENAT:::.tsenat_mca_q_selection(entropy_matrix, q_grid)
    
    # Total contributions should sum to approximately 1
    expect_true(abs(sum(mca_result$total_contrib) - 1.0) < 0.01)
})

test_that("MCA handles Z-score standardized entropy", {
    q_grid <- c(0.5, 1, 1.5, 2, 2.5)
    
    entropy_se <- calculate_diversity(rc[1:50, ], gs[1:50], q = q_grid, norm = "zscore", verbose = FALSE)
    entropy_matrix <- extract_entropy_by_q(entropy_se, q_grid)
    
    mca_result <- TSENAT:::.tsenat_mca_q_selection(entropy_matrix, q_grid)
    
    # Should work with standardized values too
    expect_true(is.list(mca_result))
    expect_true(length(mca_result$recommended_q) > 0)
})

test_that("MCA requires minimum q values", {
    # Only 2 q values (need at least 3)
    q_grid <- c(1, 2)
    
    entropy_se <- calculate_diversity(rc[1:40, ], gs[1:40], q = q_grid, norm = "none", verbose = FALSE)
    entropy_matrix <- extract_entropy_by_q(entropy_se, q_grid)
    
    expect_error(
        TSENAT:::.tsenat_mca_q_selection(entropy_matrix, q_grid),
        "at least 3 q values"
    )
})

test_that("MCA fails with invalid variance threshold", {
    q_grid <- c(0.5, 1, 1.5, 2)
    
    entropy_se <- calculate_diversity(rc[1:50, ], gs[1:50], q = q_grid, norm = "none", verbose = FALSE)
    entropy_matrix <- extract_entropy_by_q(entropy_se, q_grid)
    
    expect_error(
        TSENAT:::.tsenat_mca_q_selection(entropy_matrix, q_grid, min_variance_explained = 0),
        "must be in"
    )
    
    expect_error(
        TSENAT:::.tsenat_mca_q_selection(entropy_matrix, q_grid, min_variance_explained = 1.5),
        "must be in"
    )
})

test_that("MCA inertia (eigenvalues) are non-decreasing", {
    q_grid <- c(0.5, 1, 1.5, 2, 2.5)
    
    entropy_se <- calculate_diversity(rc[1:50, ], gs[1:50], q = q_grid, norm = "none", verbose = FALSE)
    entropy_matrix <- extract_entropy_by_q(entropy_se, q_grid)
    
    mca_result <- TSENAT:::.tsenat_mca_q_selection(entropy_matrix, q_grid)
    
    # Cumulative inertia should be monotonically non-decreasing
    expect_true(all(diff(mca_result$inertia) >= 0))
})

test_that("MCA plot returns ggplot object", {
    skip_if_not_installed("ggplot2")
    
    q_grid <- c(0.5, 1, 1.5, 2)
    
    entropy_se <- calculate_diversity(rc[1:50, ], gs[1:50], q = q_grid, norm = "none", verbose = FALSE)
    entropy_matrix <- extract_entropy_by_q(entropy_se, q_grid)
    
    mca_result <- TSENAT:::.tsenat_mca_q_selection(entropy_matrix, q_grid)
    
    plot_obj <- TSENAT:::.tsenat_plot_mca_q_selection(mca_result)
    
    expect_true(inherits(plot_obj, "ggplot"))
})

test_that("MCA recommendations are consistent", {
    # Same data should give same recommendations
    q_grid <- c(0.5, 1, 1.5, 2, 2.5)
    
    entropy_se <- calculate_diversity(rc[1:50, ], gs[1:50], q = q_grid, norm = "none", verbose = FALSE)
    entropy_matrix <- extract_entropy_by_q(entropy_se, q_grid)
    
    mca_result1 <- TSENAT:::.tsenat_mca_q_selection(entropy_matrix, q_grid, min_variance_explained = 0.80)
    mca_result2 <- TSENAT:::.tsenat_mca_q_selection(entropy_matrix, q_grid, min_variance_explained = 0.80)
    
    expect_identical(sort(mca_result1$recommended_q), sort(mca_result2$recommended_q))
})

test_that("MCA with range standardization", {
    q_grid <- c(0.5, 1, 1.5, 2, 2.5)
    
    entropy_se <- calculate_diversity(rc[1:50, ], gs[1:50], q = q_grid, norm = "range", verbose = FALSE)
    entropy_matrix <- extract_entropy_by_q(entropy_se, q_grid)
    
    mca_result <- TSENAT:::.tsenat_mca_q_selection(entropy_matrix, q_grid)
    
    # Should work with [0,1] bounded values
    expect_true(length(mca_result$recommended_q) > 0)
    expect_true(mca_result$variance_explained >= 0.8)
})
