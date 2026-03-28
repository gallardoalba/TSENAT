# ============================================================================
# TESTS FOR DIVERGENCE HELPER FUNCTIONS
# Tests for R/calculate_divergence_helpers.R helper functions
# ============================================================================

context("Divergence Helper Functions")

# INPUT VALIDATION HELPERS
# ============================================================================

test_that(".tsenat_validate_norm_parameter coerces logical to character", {
    expect_equal(.tsenat_validate_norm_parameter(TRUE), "range")
    expect_equal(.tsenat_validate_norm_parameter(FALSE), "none")
    expect_equal(.tsenat_validate_norm_parameter("zscore"), "zscore")
})

test_that(".tsenat_validate_norm_parameter rejects invalid values", {
    expect_error(.tsenat_validate_norm_parameter("invalid"),
                 "should be one of")
})

test_that(".tsenat_validate_and_sort_q_values sorts and validates", {
    result <- .tsenat_validate_and_sort_q_values(c(2, 0.5, 1))
    expect_equal(result, c(0.5, 1, 2))
})

test_that(".tsenat_validate_and_sort_q_values rejects negative q", {
    expect_error(.tsenat_validate_and_sort_q_values(c(1, -0.5)),
                 "q parameter must be >= 0")
})

test_that(".tsenat_validate_se_input accepts SummarizedExperiment", {
    se <- SummarizedExperiment::SummarizedExperiment(
        assays = list(counts = matrix(1:10, 2, 5))
    )
    expect_true(.tsenat_validate_se_input(se))
})

test_that(".tsenat_validate_se_input rejects non-SE objects", {
    expect_error(.tsenat_validate_se_input(data.frame(x = 1:5)),
                 "SummarizedExperiment")
})

# GENE COLUMN IDENTIFICATION
# ============================================================================

test_that(".tsenat_identify_gene_column finds gene_name column", {
    rd <- S4Vectors::DataFrame(gene_name = c("GENE1", "GENE2"))
    se <- SummarizedExperiment::SummarizedExperiment(
        assays = list(counts = matrix(1:10, 2, 5)),
        rowData = rd
    )
    expect_equal(.tsenat_identify_gene_column(se), "gene_name")
})

test_that(".tsenat_identify_gene_column falls back to gene_id", {
    rd <- S4Vectors::DataFrame(gene_id = c("ENSG001", "ENSG002"))
    se <- SummarizedExperiment::SummarizedExperiment(
        assays = list(counts = matrix(1:10, 2, 5)),
        rowData = rd
    )
    expect_equal(.tsenat_identify_gene_column(se), "gene_id")
})

test_that(".tsenat_identify_gene_column returns NA when no gene columns", {
    rd <- S4Vectors::DataFrame(other_col = c("A", "B"))
    se <- SummarizedExperiment::SummarizedExperiment(
        assays = list(counts = matrix(1:10, 2, 5)),
        rowData = rd
    )
    expect_true(is.na(.tsenat_identify_gene_column(se)))
})

# GENE LIST EXTRACTION
# ============================================================================

test_that(".tsenat_extract_gene_list extracts unique genes from gene_name", {
    rd <- S4Vectors::DataFrame(gene_name = c("GENE1", "GENE1", "GENE2"))
    se <- SummarizedExperiment::SummarizedExperiment(
        assays = list(counts = matrix(1:15, 3, 5)),
        rowData = rd
    )
    result <- .tsenat_extract_gene_list(se, "gene_name")
    expect_equal(sort(result), c("GENE1", "GENE2"))
})

test_that(".tsenat_extract_gene_list uses rownames when no gene column", {
    se <- SummarizedExperiment::SummarizedExperiment(
        assays = list(counts = matrix(1:10, 2, 5)),
        rowData = S4Vectors::DataFrame(x = 1:2)
    )
    rownames(se) <- c("GENE1", "GENE2")
    result <- .tsenat_extract_gene_list(se, NA_character_)
    expect_equal(result, c("GENE1", "GENE2"))
})

test_that(".tsenat_extract_gene_list rejects empty gene lists", {
    se <- SummarizedExperiment::SummarizedExperiment(
        assays = list(counts = matrix(nrow = 0, ncol = 0))
    )
    expect_error(.tsenat_extract_gene_list(se, NA_character_),
                 "gene identifiers")
})

# PARALLEL CONFIGURATION
# ============================================================================

test_that(".tsenat_configure_parallel auto-detects cores", {
    result <- .tsenat_configure_parallel(NULL, 10)
    expect_true(result$nthreads >= 1)
    expect_true(is.logical(result$use_parallel))
})

test_that(".tsenat_configure_parallel uses specified threads", {
    result <- .tsenat_configure_parallel(2, 10)
    expect_equal(result$nthreads, 2L)
})

test_that(".tsenat_configure_parallel decides parallel correctly", {
    result_seq <- .tsenat_configure_parallel(1, 3)
    expect_false(result_seq$use_parallel)
    
    result_par <- .tsenat_configure_parallel(2, 10)
    expect_true(result_par$use_parallel)
})

test_that(".tsenat_configure_parallel rejects invalid threads", {
    expect_error(.tsenat_configure_parallel(-1, 10),
                 "positive integer")
    expect_error(.tsenat_configure_parallel("invalid", 10),
                 "positive integer")
})

# GENE PROCESSING HELPERS
# ============================================================================

test_that(".tsenat_compute_aggregate_counts sums transcript counts", {
    counts_matrix <- matrix(c(1, 2, 3, 4, 5, 6), 2, 3)
    rd <- S4Vectors::DataFrame(gene_name = c("GENE1", "GENE1"))
    
    se <- SummarizedExperiment::SummarizedExperiment(
        assays = list(counts = counts_matrix),
        rowData = rd
    )
    
    result <- .tsenat_compute_aggregate_counts(se, "GENE1", "gene_name", rd)
    expect_equal(result, c(3, 7, 11))  # colSums of the two rows
})

test_that(".tsenat_compute_aggregate_counts returns NULL when gene not found", {
    counts_matrix <- matrix(1:6, 2, 3)
    rd <- S4Vectors::DataFrame(gene_name = c("GENE1", "GENE2"))
    
    se <- SummarizedExperiment::SummarizedExperiment(
        assays = list(counts = counts_matrix),
        rowData = rd
    )
    
    result <- .tsenat_compute_aggregate_counts(se, "MISSING", "gene_name", rd)
    expect_null(result)
})

test_that(".tsenat_extract_group_counts_gene splits by group", {
    counts <- c(10, 20, 30, 40, 50)
    groups <- c("A", "A", "B", "B", "B")
    
    result <- .tsenat_extract_group_counts_gene(counts, groups, "A")
    expect_equal(result$control, c(10, 20))
    expect_equal(result$treatment, c(30, 40, 50))
})

test_that(".tsenat_extract_group_counts_gene errors on length mismatch", {
    counts <- c(10, 20, 30)
    groups <- c("A", "B")
    
    expect_error(.tsenat_extract_group_counts_gene(counts, groups, "A"),
                 "Length mismatch")
})

test_that(".tsenat_make_error_result creates proper error structure", {
    result <- .tsenat_make_error_result("GENE1", c(0.5, 1, 2), "Test error", 1.5)
    
    expect_equal(result$gene_name, "GENE1")
    expect_equal(result$error, "Test error")
    expect_equal(result$computation_time_sec, 1.5)
    expect_equal(length(result$results_per_q), 3)
    expect_true(is.na(result$results_per_q[[1]]$estimate))
})

test_that(".tsenat_bootstrap_build_args builds arguments correctly", {
    x <- rnorm(10)
    y <- rnorm(10)
    
    result <- .tsenat_bootstrap_build_args(x, y, 1.0, 100, 0.95, "percentile",
                                     exp(1), 0.5, "GENE1", 42, NULL)
    
    expect_equal(result$x, x)
    expect_equal(result$y, y)
    expect_equal(result$q, 1.0)
    expect_equal(result$paired, FALSE)
    expect_null(result$pair_ids)
})

test_that(".tsenat_bootstrap_build_args includes pair_ids when provided", {
    x <- rnorm(10)
    y <- rnorm(10)
    pair_ids <- c(1, 1, 2, 2, 3)
    
    result <- .tsenat_bootstrap_build_args(x, y, 1.0, 100, 0.95, "percentile",
                                     exp(1), 0.5, "GENE1", 42, pair_ids)
    
    expect_equal(result$pair_ids, pair_ids)
    expect_equal(result$paired, TRUE)
})

# RESULTS COMPILATION HELPERS
# ============================================================================

test_that(".tsenat_initialize_matrices creates proper structure", {
    result <- .tsenat_initialize_matrices(5, c(0.5, 1, 2))
    
    expect_equal(nrow(result$assay), 5)
    expect_equal(ncol(result$assay), 3)
    expect_equal(nrow(result$rowData), 5)
    expect_true("gene_name" %in% colnames(result$rowData))
    expect_true("error" %in% colnames(result$rowData))
    expect_true("estimate_q0.5" %in% colnames(result$rowData))
})

test_that(".tsenat_initialize_matrices creates all q columns", {
    q_vals <- c(0.5, 1, 1.5, 2)
    result <- .tsenat_initialize_matrices(3, q_vals)
    
    for (q in q_vals) {
        expect_true(paste0("estimate_q", q) %in% colnames(result$rowData))
        expect_true(paste0("lower_ci_q", q) %in% colnames(result$rowData))
        expect_true(paste0("upper_ci_q", q) %in% colnames(result$rowData))
        expect_true(paste0("ci_width_q", q) %in% colnames(result$rowData))
        expect_true(paste0("method_q", q) %in% colnames(result$rowData))
        expect_true(paste0("nboot_q", q) %in% colnames(result$rowData))
    }
})

# NORMALIZATION HELPERS
# ============================================================================

test_that(".tsenat_normalize_range_matrix scales to [0,1]", {
    assay <- matrix(c(0, 5, 10, 1, 6, 11), 2, 3)
    row_data <- data.frame(
        estimate_q1 = c(0, 5),
        lower_ci_q1 = c(-1, 4),
        upper_ci_q1 = c(1, 6),
        stringsAsFactors = FALSE
    )
    
    result <- .tsenat_normalize_range_matrix(assay, row_data, c(1))
    
    expect_true(all(result$assay >= 0, na.rm = TRUE))
    expect_true(all(result$assay <= 1, na.rm = TRUE))
    expect_equal(result$assay[1, 1], 0)
    expect_equal(result$assay[2, ncol(result$assay)], 1)
})

test_that(".tsenat_divergence_normalize_zscore normalizes each column", {
    assay <- matrix(c(1:6), 2, 3)
    row_data <- data.frame(
        estimate_q1 = c(1, 4),
        lower_ci_q1 = c(2, 5),
        upper_ci_q1 = c(3, 6),
        stringsAsFactors = FALSE
    )
    
    result <- .tsenat_divergence_normalize_zscore(assay, row_data, c(1))
    
    # Z-score normalization is applied per column, so check individual columns
    col1 <- na.omit(result$assay[, 1])
    col2 <- na.omit(result$assay[, 2])
    col3 <- na.omit(result$assay[, 3])
    
    # Each column should have mean ~0 and sd ~1
    expect_true(abs(mean(col1)) < 0.01)
    expect_true(abs(sd(col1) - 1) < 0.01)
    expect_true(abs(mean(col2)) < 0.01)
    expect_true(abs(sd(col2) - 1) < 0.01)
    expect_true(abs(mean(col3)) < 0.01)
    expect_true(abs(sd(col3) - 1) < 0.01)
})

# APPLY NORMALIZATION DISPATCHER
# ============================================================================

test_that(".tsenat_normalize_divergence_matrix handles all methods", {
    assay <- matrix(1:6, 2, 3)
    row_data <- data.frame(
        estimate_q1 = 1:2,
        lower_ci_q1 = 1:2,
        upper_ci_q1 = 1:2
    )
    
    # Test each method
    result_none <- .tsenat_normalize_divergence_matrix(assay, row_data, c(1), "none")
    expect_equal(result_none$assay, assay)
    
    result_range <- .tsenat_normalize_divergence_matrix(assay, row_data, c(1), "range")
    expect_true(all(result_range$assay >= 0, na.rm = TRUE))
    
    result_zscore <- .tsenat_normalize_divergence_matrix(assay, row_data, c(1), "zscore")
    expect_true(!is.null(result_zscore))
})

# COMPUTE SUMMARY STATISTICS
# ============================================================================

test_that(".tsenat_generate_summary creates summary", {
    row_data_df <- data.frame(
        gene_name = c("G1", "G2", "G3"),
        error = c(NA_character_, NA_character_, "error"),
        stringsAsFactors = FALSE
    )
    
    result <- .tsenat_generate_summary(elapsed = 10, num_genes = 3, 
                                            num_errors = 1, row_data_df)
    
    expect_equal(result$successful, 2)
    expect_equal(result$failed, 1)
    expect_equal(result$total_elapsed, 10)
    expect_true(is.null(result$failed_details) || nrow(result$failed_details) <= 1)
})

# SE CONSTRUCTION
# ============================================================================

test_that(".tsenat_construct_result_se creates SummarizedExperiment", {
    assay <- matrix(0.1, 3, 2, dimnames = list(NULL, c("q_0.5", "q_1")))
    row_data <- data.frame(
        gene_name = c("G1", "G2", "G3"),
        error = c(NA, NA, NA),
        computation_time_sec = c(1, 1.5, 1.2),
        stringsAsFactors = FALSE
    )
    
    result <- .tsenat_construct_result_se(
        assay, row_data, c(0.5, 1),
        elapsed = 5, nboot = 1000, ci = 0.95, method = "percentile",
        norm = "none", use_parallel = FALSE, num_genes = 3, num_errors = 0
    )
    
    expect_true(methods::is(result, "SummarizedExperiment"))
    expect_equal(nrow(result), 3)
    expect_equal(ncol(result), 2)
    expect_true("divergence" %in% names(SummarizedExperiment::assays(result)))
})

test_that(".tsenat_construct_result_se preserves metadata", {
    assay <- matrix(0.1, 2, 1)
    row_data <- data.frame(
        gene_name = c("G1", "G2"),
        error = c(NA, NA),
        computation_time_sec = c(1, 1)
    )
    
    result <- .tsenat_construct_result_se(
        assay, row_data, c(1),
        elapsed = 2, nboot = 500, ci = 0.95, method = "bca",
        norm = "range", use_parallel = TRUE, num_genes = 2, num_errors = 0
    )
    
    meta <- S4Vectors::metadata(result)
    expect_equal(meta$bootstrap_config$nboot, 500)
    expect_equal(meta$bootstrap_config$method, "bca")
    expect_equal(meta$normalization, "range")
    expect_equal(meta$computation_mode, "parallel")
})
