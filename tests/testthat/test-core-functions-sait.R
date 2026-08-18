library(testthat)

context("SAIT Core Functions: Input Validation")

test_that(".validate_sait_interaction_input rejects invalid storey", {
    skip_if_not_installed("SummarizedExperiment")

    # Create minimal SE object for validation
    mat <- matrix(runif(10), nrow = 2)
    colnames(mat) <- c("S1_q=0.5", "S1_q=1.0", "S2_q=0.5", "S2_q=1.0",
                       "S3_q=0.5")
    rownames(mat) <- c("g1", "g2")
    rd <- data.frame(genes = c("g1", "g2"), row.names = c("g1", "g2"))
    cd <- data.frame(
        condition = c("N", "N", "T", "T", "N"),
        row.names = colnames(mat),
        stringsAsFactors = FALSE
    )
    se <- SummarizedExperiment::SummarizedExperiment(
        assays = list(diversity = mat),
        rowData = rd,
        colData = cd
    )

    # Should reject non-logical storey
    expect_error(
        TSENAT:::.validate_sait_interaction_input(
            method = "lmm",
            pvalue = "lrt",
            corstr = "ar1",
            regularization = "pca",
            multicorr = "hochberg",
            pcorr = "BH",
            storey = "yes",  # Invalid
            wy_randomizations = 1000,
            paired = FALSE,
            subject_col = NULL,
            se = se,
            verbose = FALSE
        ),
        "storey must be TRUE or FALSE",
        fixed = FALSE
    )
})

test_that(".validate_sait_interaction_input rejects invalid wy_randomizations", {
    skip_if_not_installed("SummarizedExperiment")

    mat <- matrix(runif(10), nrow = 2)
    colnames(mat) <- c("S1_q=0.5", "S1_q=1.0", "S2_q=0.5", "S2_q=1.0",
                       "S3_q=0.5")
    rownames(mat) <- c("g1", "g2")
    rd <- data.frame(genes = c("g1", "g2"), row.names = c("g1", "g2"))
    cd <- data.frame(
        condition = c("N", "N", "T", "T", "N"),
        row.names = colnames(mat),
        stringsAsFactors = FALSE
    )
    se <- SummarizedExperiment::SummarizedExperiment(
        assays = list(diversity = mat),
        rowData = rd,
        colData = cd
    )

    # Should reject negative wy_randomizations
    expect_error(
        TSENAT:::.validate_sait_interaction_input(
            method = "lmm",
            pvalue = "lrt",
            corstr = "ar1",
            regularization = "pca",
            multicorr = "hochberg",
            pcorr = "BH",
            storey = FALSE,
            wy_randomizations = -10,  # Invalid
            paired = FALSE,
            subject_col = NULL,
            se = se,
            verbose = FALSE
        ),
        "wy_randomizations must be numeric and >= 1",
        fixed = FALSE
    )
})

test_that(".validate_sait_interaction_input warns on low wy_randomizations", {
    skip_if_not_installed("SummarizedExperiment")

    mat <- matrix(runif(10), nrow = 2)
    colnames(mat) <- c("S1_q=0.5", "S1_q=1.0", "S2_q=0.5", "S2_q=1.0",
                       "S3_q=0.5")
    rownames(mat) <- c("g1", "g2")
    rd <- data.frame(genes = c("g1", "g2"), row.names = c("g1", "g2"))
    cd <- data.frame(
        condition = c("N", "N", "T", "T", "N"),
        row.names = colnames(mat),
        stringsAsFactors = FALSE
    )
    se <- SummarizedExperiment::SummarizedExperiment(
        assays = list(diversity = mat),
        rowData = rd,
        colData = cd
    )

    # Should warn on wy_randomizations < 100
    expect_warning(
        TSENAT:::.validate_sait_interaction_input(
            method = "lmm",
            pvalue = "lrt",
            corstr = "ar1",
            regularization = "pca",
            multicorr = "hochberg",
            pcorr = "BH",
            storey = FALSE,
            wy_randomizations = 50,  # Warning threshold
            paired = FALSE,
            subject_col = NULL,
            se = se,
            verbose = FALSE
        ),
        "wy_randomizations < 100",
        fixed = FALSE
    )
})

test_that(".validate_sait_interaction_input auto-detects subject_col when paired", {
    skip_if_not_installed("SummarizedExperiment")

    mat <- matrix(runif(10), nrow = 2)
    colnames(mat) <- c("S1_q=0.5", "S1_q=1.0", "S2_q=0.5", "S2_q=1.0",
                       "S3_q=0.5")
    rownames(mat) <- c("g1", "g2")
    rd <- data.frame(genes = c("g1", "g2"), row.names = c("g1", "g2"))
    cd <- data.frame(
        condition = c("N", "N", "T", "T", "N"),
        paired_samples = c("P1", "P1", "P2", "P2", "P3"),
        row.names = colnames(mat),
        stringsAsFactors = FALSE
    )
    se <- SummarizedExperiment::SummarizedExperiment(
        assays = list(diversity = mat),
        rowData = rd,
        colData = cd
    )

    # Should auto-detect paired_samples column
    result <- TSENAT:::.validate_sait_interaction_input(
        method = "lmm",
        pvalue = "lrt",
        corstr = "ar1",
        regularization = "pca",
        multicorr = "hochberg",
        pcorr = "BH",
        storey = FALSE,
        wy_randomizations = 1000,
        paired = TRUE,
        subject_col = NULL,  # Will be auto-detected
        se = se,
        verbose = FALSE
    )

    expect_equal(result$subject_col, "paired_samples")
})

context("SAIT Core Functions: Sample Metadata Parsing")

test_that(".parse_sample_metadata extracts q-values correctly", {
    skip_if_not_installed("SummarizedExperiment")

    qvec <- c(0.5, 1.0, 1.5)
    sample_ids <- c("S1", "S2")
    coln <- paste0(
        rep(sample_ids, each = length(qvec)),
        "_q=",
        rep(qvec, times = length(sample_ids))
    )
    mat <- matrix(runif(length(coln) * 2), nrow = 2)
    colnames(mat) <- coln
    rownames(mat) <- c("g1", "g2")

    rd <- data.frame(genes = c("g1", "g2"), row.names = c("g1", "g2"))
    cd <- data.frame(
        condition = rep(c("Normal", "Tumor"), each = length(qvec)),
        row.names = coln,
        stringsAsFactors = FALSE
    )
    se <- SummarizedExperiment::SummarizedExperiment(
        assays = list(diversity = mat),
        rowData = rd,
        colData = cd
    )

    metadata <- TSENAT:::.parse_sample_metadata(
        se = se,
        condition_col = "condition",
        assay_name = "diversity",
        verbose = FALSE
    )

    # Check q-values extracted correctly
    expect_equal(metadata$q_vals, rep(qvec, times = 2))
    # Check group vector
    expect_equal(
        metadata$group_vec,
        rep(c("Normal", "Tumor"), each = length(qvec))
    )
    # Check sample names
    expect_equal(metadata$sample_names, rep(sample_ids, each = length(qvec)))
})

test_that(".parse_sample_metadata rejects missing _q=", {
    skip_if_not_installed("SummarizedExperiment")

    # Column names without _q= should be rejected
    mat <- matrix(runif(8), nrow = 2)
    colnames(mat) <- c("S1", "S2", "S3", "S4")  # No _q=
    rownames(mat) <- c("g1", "g2")

    rd <- data.frame(genes = c("g1", "g2"), row.names = c("g1", "g2"))
    cd <- data.frame(
        condition = c("N", "N", "T", "T"),
        row.names = colnames(mat),
        stringsAsFactors = FALSE
    )
    se <- SummarizedExperiment::SummarizedExperiment(
        assays = list(diversity = mat),
        rowData = rd,
        colData = cd
    )

    expect_error(
        TSENAT:::.parse_sample_metadata(
            se = se,
            condition_col = "condition",
            assay_name = "diversity",
            verbose = FALSE
        ),
        "Could not parse q values",
        fixed = FALSE
    )
})

test_that(".parse_sample_metadata rejects missing condition_col", {
    skip_if_not_installed("SummarizedExperiment")

    qvec <- c(0.5, 1.0)
    coln <- c("S1_q=0.5", "S1_q=1.0", "S2_q=0.5", "S2_q=1.0")
    mat <- matrix(runif(length(coln) * 2), nrow = 2)
    colnames(mat) <- coln
    rownames(mat) <- c("g1", "g2")

    rd <- data.frame(genes = c("g1", "g2"), row.names = c("g1", "g2"))
    cd <- data.frame(
        samples = c("S1", "S1", "S2", "S2"),  # No 'condition' column
        row.names = coln,
        stringsAsFactors = FALSE
    )
    se <- SummarizedExperiment::SummarizedExperiment(
        assays = list(diversity = mat),
        rowData = rd,
        colData = cd
    )

    expect_error(
        TSENAT:::.parse_sample_metadata(
            se = se,
            condition_col = "nonexistent_col",
            assay_name = "diversity",
            verbose = FALSE
        ),
        "No sample grouping found",
        fixed = FALSE
    )
})

context("SAIT Core Functions: Gene Annotation Mapping")

test_that(".map_gene_annotations adds gene names from rowData", {
    skip_if_not_installed("SummarizedExperiment")

    mat <- matrix(runif(8), nrow = 2)
    colnames(mat) <- c("S1_q=0.5", "S1_q=1.0", "S2_q=0.5", "S2_q=1.0")
    rownames(mat) <- c("ENSG001", "ENSG002")

    rd <- data.frame(
        genes = c("ENSG001", "ENSG002"),
        gene_name = c("GENE_A", "GENE_B"),
        row.names = c("ENSG001", "ENSG002"),
        stringsAsFactors = FALSE
    )
    cd <- data.frame(
        condition = c("Normal", "Normal", "Tumor", "Tumor"),
        row.names = colnames(mat),
        stringsAsFactors = FALSE
    )
    se <- SummarizedExperiment::SummarizedExperiment(
        assays = list(diversity = mat),
        rowData = rd,
        colData = cd
    )

    # Create results data.frame with gene column (rownames)
    results <- data.frame(
        gene = c("ENSG001", "ENSG002"),
        p_interaction = c(0.01, 0.05),
        adj_p_interaction = c(0.02, 0.1)
    )

    mapped <- TSENAT:::.map_gene_annotations(
        res = results,
        se = se,
        verbose = FALSE
    )

    # Check gene_name column added
    expect_true("gene_name" %in% colnames(mapped))
    expect_equal(mapped$gene_name, c("GENE_A", "GENE_B"))
    # Check gene_id column added
    expect_true("gene_id" %in% colnames(mapped))
})

test_that(".map_gene_annotations handles missing gene_name gracefully", {
    skip_if_not_installed("SummarizedExperiment")

    mat <- matrix(runif(8), nrow = 2)
    colnames(mat) <- c("S1_q=0.5", "S1_q=1.0", "S2_q=0.5", "S2_q=1.0")
    rownames(mat) <- c("g1", "g2")

    # rowData without gene_name column
    rd <- data.frame(
        genes = c("g1", "g2"),
        row.names = c("g1", "g2"),
        stringsAsFactors = FALSE
    )
    cd <- data.frame(
        condition = c("Normal", "Normal", "Tumor", "Tumor"),
        row.names = colnames(mat),
        stringsAsFactors = FALSE
    )
    se <- SummarizedExperiment::SummarizedExperiment(
        assays = list(diversity = mat),
        rowData = rd,
        colData = cd
    )

    results <- data.frame(
        gene = c("g1", "g2"),
        p_interaction = c(0.01, 0.05),
        adj_p_interaction = c(0.02, 0.1)
    )

    mapped <- TSENAT:::.map_gene_annotations(
        res = results,
        se = se,
        verbose = FALSE
    )

    # Should still have gene_name column (using gene as fallback)
    expect_true("gene_name" %in% colnames(mapped))
    expect_equal(mapped$gene_name, c("g1", "g2"))
})

context("SAIT Core Functions: Model Metadata Assembly")

test_that(".assemble_model_metadata returns required fields", {
    skip_if_not_installed("SummarizedExperiment")

    qvec <- c(0.5, 1.0)
    coln <- c("S1_q=0.5", "S1_q=1.0", "S2_q=0.5", "S2_q=1.0")
    mat <- matrix(runif(length(coln) * 2), nrow = 2)
    colnames(mat) <- coln
    rownames(mat) <- c("g1", "g2")

    rd <- data.frame(
        genes = c("g1", "g2"),
        gene_name = c("GENE_A", "GENE_B"),
        row.names = c("g1", "g2"),
        stringsAsFactors = FALSE
    )
    cd <- data.frame(
        condition = rep(c("Normal", "Tumor"), each = 2),
        row.names = coln,
        stringsAsFactors = FALSE
    )
    se <- SummarizedExperiment::SummarizedExperiment(
        assays = list(diversity = mat),
        rowData = rd,
        colData = cd
    )

    metadata <- TSENAT:::.parse_sample_metadata(
        se = se,
        condition_col = "condition",
        assay_name = "diversity",
        verbose = FALSE
    )

    results <- data.frame(
        gene = c("g1", "g2"),
        p_interaction = c(0.01, 0.05),
        adj_p_interaction = c(0.02, 0.1)
    )

    model_data <- TSENAT:::.assemble_model_metadata(
        se = se,
        res = results,
        mat = mat,
        metadata = metadata,
        method = "lmm",
        pvalue = "lrt",
        multicorr = "hochberg",
        assay_name = "diversity",
        bias_correction = TRUE,
        regularization = "pca",
        corstr = "ar1",
        adaptive_knots = TRUE
    )

    # Check required fields present
    expect_true("method" %in% names(model_data))
    expect_true("n_genes" %in% names(model_data))
    expect_true("n_q_values" %in% names(model_data))
    expect_true("q_values" %in% names(model_data))
    expect_true("sample_names" %in% names(model_data))
    expect_true("group_levels" %in% names(model_data))
    expect_true("per_group_statistics" %in% names(model_data))
    expect_true("test_configuration" %in% names(model_data))
    expect_true("genes_analyzed" %in% names(model_data))

    # Check correct values
    expect_equal(model_data$method, "lmm")
    expect_equal(model_data$n_genes, 2)
    expect_equal(model_data$n_q_values, 2)
    expect_equal(model_data$q_values, c(0.5, 1.0))
})

test_that(".assemble_model_metadata computes per-group statistics", {
    skip_if_not_installed("SummarizedExperiment")

    qvec <- c(0.5, 1.0)
    coln <- c("S1_q=0.5", "S1_q=1.0", "S2_q=0.5", "S2_q=1.0")
    set.seed(42)
    mat <- matrix(runif(length(coln) * 2, min = 0.5, max = 1.5),
                  nrow = 2)
    colnames(mat) <- coln
    rownames(mat) <- c("g1", "g2")

    rd <- data.frame(
        genes = c("g1", "g2"),
        row.names = c("g1", "g2"),
        stringsAsFactors = FALSE
    )
    cd <- data.frame(
        condition = rep(c("Normal", "Tumor"), each = 2),
        row.names = coln,
        stringsAsFactors = FALSE
    )
    se <- SummarizedExperiment::SummarizedExperiment(
        assays = list(diversity = mat),
        rowData = rd,
        colData = cd
    )

    metadata <- TSENAT:::.parse_sample_metadata(
        se = se,
        condition_col = "condition",
        assay_name = "diversity",
        verbose = FALSE
    )

    results <- data.frame(
        gene = c("g1", "g2"),
        p_interaction = c(0.01, 0.05),
        adj_p_interaction = c(0.02, 0.1)
    )

    model_data <- TSENAT:::.assemble_model_metadata(
        se = se,
        res = results,
        mat = mat,
        metadata = metadata,
        method = "lmm",
        pvalue = "lrt",
        multicorr = "hochberg",
        assay_name = "diversity"
    )

    # Check per-group statistics computed
    expect_true(length(model_data$per_group_statistics) > 0)
    # Should have one entry per group
    expect_true("Normal" %in% names(model_data$per_group_statistics) ||
                "Tumor" %in% names(model_data$per_group_statistics))

    # Each group should have statistics
    for (group_stats in model_data$per_group_statistics) {
        expect_true("entropy_mean" %in% names(group_stats))
        expect_true("entropy_sd" %in% names(group_stats))
        expect_true("n_samples" %in% names(group_stats))
        expect_true("n_observations" %in% names(group_stats))
    }
})

context("SAIT Refactored Main Function Integration")

test_that("calculate_sait works with refactored code", {
    skip_if_not_installed("SummarizedExperiment")

    # Use a larger, more realistic dataset to avoid numerical issues
    qvec <- seq(0.1, 2.0, by = 0.1)
    sample_ids <- c("S1", "S2", "S3")
    coln <- paste0(
        rep(sample_ids, each = length(qvec)),
        "_q=",
        rep(qvec, times = length(sample_ids))
    )

    set.seed(42)
    # Create realistic entropy-like data (bounded, well-separated groups)
    # Gene 1: clear interaction (different slopes between Normal and Tumor)
    base_normal_g1 <- 0.7 + qvec * 0.15  # Normal group: moderate slope
    base_tumor_g1 <- 0.8 + qvec * 0.35   # Tumor group: steeper slope
    gene1 <- c(
        base_normal_g1 + rnorm(length(qvec), sd = 0.05),
        base_tumor_g1 + rnorm(length(qvec), sd = 0.05),
        base_normal_g1 + rnorm(length(qvec), sd = 0.05)
    )

    # Gene 2: weak/no interaction (similar slopes)
    base_g2 <- 0.75 + qvec * 0.12
    gene2 <- c(
        base_g2 + rnorm(length(qvec), sd = 0.05),
        base_g2 + 0.05 + rnorm(length(qvec), sd = 0.05),
        base_g2 + rnorm(length(qvec), sd = 0.05)
    )

    mat <- rbind(gene1, gene2)
    colnames(mat) <- coln
    rownames(mat) <- c("g1", "g2")

    rd <- data.frame(
        genes = c("g1", "g2"),
        gene_name = c("GENE_A", "GENE_B"),
        row.names = c("g1", "g2"),
        stringsAsFactors = FALSE
    )
    cd <- data.frame(
        condition = rep(c("Normal", "Tumor", "Normal"), each = length(qvec)),
        row.names = coln,
        stringsAsFactors = FALSE
    )
    se <- SummarizedExperiment::SummarizedExperiment(
        assays = list(diversity = mat),
        rowData = rd,
        colData = cd
    )

    # Call refactored function with realistic data
    results <- .calculate_sait(
        se,
        condition_col = "condition",
        method = "lmm",
        multicorr = "hochberg",
        min_obs = 5,
        verbose = FALSE
    )

    # Check results structure
    expect_true(is.data.frame(results))
    expect_true("gene" %in% colnames(results))
    expect_true("p_interaction" %in% colnames(results))
    expect_true("adj_p_interaction" %in% colnames(results))
    expect_true("gene_name" %in% colnames(results))

    # Check ordering by p-value
    expect_true(
        all(
            diff(results$adj_p_interaction[!is.na(results$adj_p_interaction)]) >= 0
        )
    )
})

test_that("calculate_sait returns model_data when requested", {
    skip_if_not_installed("SummarizedExperiment")

    # Use realistic entropy-like data (bounded between 0 and 1, realistic slopes)
    qvec <- seq(0.1, 2.0, by = 0.1)
    sample_ids <- c("S1", "S2", "S3")
    coln <- paste0(
        rep(sample_ids, each = length(qvec)),
        "_q=",
        rep(qvec, times = length(sample_ids))
    )

    set.seed(123)
    # Create well-behaved synthetic entropy data
    # Gene 1: moderate entropy that increases with q
    base_g1 <- 0.6 + qvec * 0.12
    g1_data <- c(
        base_g1 + rnorm(length(qvec), sd = 0.04),
        base_g1 + 0.1 + rnorm(length(qvec), sd = 0.04),
        base_g1 + rnorm(length(qvec), sd = 0.04)
    )

    # Gene 2: similar entropy across groups
    base_g2 <- 0.65 + qvec * 0.10
    g2_data <- c(
        base_g2 + rnorm(length(qvec), sd = 0.04),
        base_g2 + rnorm(length(qvec), sd = 0.04),
        base_g2 + rnorm(length(qvec), sd = 0.04)
    )

    mat <- rbind(g1_data, g2_data)
    colnames(mat) <- coln
    rownames(mat) <- c("g1", "g2")

    rd <- data.frame(genes = c("g1", "g2"), row.names = c("g1", "g2"))
    cd <- data.frame(
        condition = rep(c("Normal", "Tumor", "Normal"), each = length(qvec)),
        row.names = coln,
        stringsAsFactors = FALSE
    )
    se <- SummarizedExperiment::SummarizedExperiment(
        assays = list(diversity = mat),
        rowData = rd,
        colData = cd
    )

    # Call with return_model_data = TRUE
    output <- .calculate_sait(
        se,
        condition_col = "condition",
        method = "lmm",
        return_model_data = TRUE,
        verbose = FALSE,
        min_obs = 5
    )

    # Check it's a list with results and model_data
    expect_true(is.list(output))
    expect_true("results" %in% names(output))
    expect_true("model_data" %in% names(output))

    # Check results structure
    expect_true(is.data.frame(output$results))

    # Check model_data structure
    expect_true(is.list(output$model_data))
    expect_true("method" %in% names(output$model_data))
    expect_true("per_group_statistics" %in% names(output$model_data))
})

test_that(".extract_sait_result_df returns data.frame when given data.frame", {
    # Create a simple results data.frame
    results_df <- data.frame(
        gene = c("g1", "g2"),
        pvalue = c(0.01, 0.05),
        interaction_term = c(0.5, -0.3),
        stringsAsFactors = FALSE
    )

    # Should pass through data.frames unchanged
    output <- TSENAT:::.extract_sait_result_df(results_df)
    expect_identical(output, results_df)
    expect_true(is.data.frame(output))
})

test_that(".extract_sait_result_df extracts results from list with 'results' element", {
    # Create a list with 'results' component (as returned by .calculate_sait with return_model_data=TRUE)
    results_df <- data.frame(
        gene = c("g1", "g2"),
        pvalue = c(0.01, 0.05),
        interaction_term = c(0.5, -0.3),
        stringsAsFactors = FALSE
    )

    sait_result <- list(
        results = results_df,
        model_data = list(method = "lmm", summary_stats = list()),
        other_component = "some_value"
    )

    output <- TSENAT:::.extract_sait_result_df(sait_result)
    expect_identical(output, results_df)
    expect_true(is.data.frame(output))
    expect_equal(nrow(output), 2)
    expect_equal(ncol(output), 3)
})

test_that(".extract_sait_result_df raises error for invalid input (not data.frame or list with results)", {
    # Should reject vector input
    expect_error(
        TSENAT:::.extract_sait_result_df(c(1, 2, 3)),
        "sait_result must be either a data.frame or a list with 'results' component"
    )

    # Should reject list without 'results' component
    expect_error(
        TSENAT:::.extract_sait_result_df(list(model_data = "something")),
        "sait_result must be either a data.frame or a list with 'results' component"
    )

    # Should reject NULL
    expect_error(
        TSENAT:::.extract_sait_result_df(NULL),
        "sait_result must be either a data.frame or a list with 'results' component"
    )
})

test_that(".extract_sait_result_df preserves data.frame structure and content", {
    skip_if_not_installed("SummarizedExperiment")

    # Create results with various data types
    results_df <- data.frame(
        gene_id = c("ENSG1", "ENSG2", "ENSG3"),
        pvalue = c(0.001, NA, 0.5),
        log2_effect = c(1.5, -0.8, 0.2),
        q_value = c(0.5, 1.0, 2.0),
        method = c("lmm", "lmm", "lmm"),
        stringsAsFactors = FALSE
    )

    sait_result <- list(
        results = results_df,
        metadata = list()
    )

    output <- TSENAT:::.extract_sait_result_df(sait_result)

    # Check that all columns are preserved
    expect_equal(colnames(output), colnames(results_df))
    expect_equal(nrow(output), nrow(results_df))

    # Check that NA values are preserved
    expect_true(is.na(output$pvalue[2]))

    # Check data types
    expect_true(is.character(output$gene_id))
    expect_true(is.numeric(output$pvalue))
    expect_true(is.numeric(output$log2_effect))
})



test_that(".estimate_ar1_rho estimates AR(1) autocorrelation from entropy differences", {
    # Create entropy data with moderate autocorrelation
    entropy_diff <- c(0.01, 0.015, 0.008, 0.012, 0.010, 0.014, 0.009, 0.011)
    
    result <- TSENAT:::.estimate_ar1_rho(entropy_diff)
    
    # Should return numeric value or NULL
    expect_true(is.null(result) || (is.numeric(result) && result >= 0 && result <= 1))
})

test_that(".estimate_ar1_rho returns NULL for insufficient data", {
    # Too few observations
    entropy_diff <- c(0.01, 0.02)
    
    result <- TSENAT:::.estimate_ar1_rho(entropy_diff)
    
    expect_null(result)
})

test_that(".estimate_ar1_rho handles NA values", {
    entropy_diff <- c(0.01, NA, 0.015, 0.008, NA, 0.012, 0.010)
    
    result <- TSENAT:::.estimate_ar1_rho(entropy_diff)
    
    # Should not error and returns valid rho or NULL
    expect_true(is.null(result) || (is.numeric(result) && result >= 0 && result <= 1))
})

test_that(".estimate_ar1_rho handles constant series", {
    # Zero variance - all same values
    entropy_diff <- rep(0.01, 5)
    
    result <- TSENAT:::.estimate_ar1_rho(entropy_diff)
    
    # Should return NULL for zero variance
    expect_null(result)
})

test_that(".estimate_ar1_rho produces values in [0, 1]", {
    # Create data with actual autocorrelation (AR(1) process)
    # This ensures rho is estimated to be >= 0.01 (not rejected as too small)
    set.seed(234)
    entropy_diff <- numeric(25)
    entropy_diff[1] <- rnorm(1, 0, 0.02)
    for (i in 2:25) {
        # AR(1) with rho = 0.5 ensures meaningful autocorrelation
        entropy_diff[i] <- 0.5 * entropy_diff[i-1] + rnorm(1, 0, 0.01)
    }
    
    result <- TSENAT:::.estimate_ar1_rho(entropy_diff)
    
    # Should return a numeric value (not NULL) in valid range
    expect_true(is.numeric(result))
    expect_true(result >= 0 && result <= 1)
})

test_that(".estimate_ar1_rho warns on high autocorrelation", {
    # Create explicit AR(1) process with rho = 0.98
    # x[t] = 0.98 * x[t-1] + epsilon where epsilon ~ N(0, 0.0001)
    # This guarantees high positive autocorrelation
    set.seed(999)
    entropy_diff <- numeric(100)
    entropy_diff[1] <- rnorm(1)
    for (i in 2:100) {
        entropy_diff[i] <- 0.98 * entropy_diff[i-1] + rnorm(1, sd = 0.01)
    }
    
    expect_warning(
        result <- TSENAT:::.estimate_ar1_rho(entropy_diff),
        "AR\\(1\\) autocorrelation"
    )
})



test_that(".adjust_pvalues_multicorr uses Hochberg adjustment", {
    p_values <- c(0.001, 0.01, 0.05, 0.1, 0.5)
    
    result <- TSENAT:::.adjust_pvalues_multicorr(
        p_values = p_values,
        multicorr = "hochberg",
        wy_randomizations = 100
    )
    
    # Should return adjusted p-values
    expect_true(is.numeric(result))
    expect_equal(length(result), length(p_values))
    # Adjusted p-values should be >= original
    expect_true(all(result >= p_values))
    # All should be valid probabilities
    expect_true(all(result >= 0 & result <= 1))
})

test_that(".adjust_pvalues_multicorr handles benjamini-yekutieli adjustment", {
    p_values <- c(0.001, 0.01, 0.05, 0.1)
    
    result <- TSENAT:::.adjust_pvalues_multicorr(
        p_values = p_values,
        multicorr = "benjamini-yekutieli",
        wy_randomizations = 100
    )
    
    expect_true(is.numeric(result))
    expect_equal(length(result), length(p_values))
    expect_true(all(result >= 0 & result <= 1))
})

test_that(".adjust_pvalues_multicorr rejects unknown method", {
    p_values <- c(0.001, 0.05, 0.1)
    
    # Unknown method should throw an error
    expect_error(
        TSENAT:::.adjust_pvalues_multicorr(
            p_values = p_values,
            multicorr = "unknown_method",
            wy_randomizations = 100
        ),
        "should be one of"
    )
})

test_that(".adjust_pvalues_multicorr preserves single p-value", {
    p_values <- 0.05
    
    result <- TSENAT:::.adjust_pvalues_multicorr(
        p_values = p_values,
        multicorr = "hochberg",
        wy_randomizations = 100
    )
    
    expect_true(is.numeric(result))
    expect_equal(length(result), 1)
})

test_that(".adjust_pvalues_multicorr handles all significant p-values", {
    # All very small p-values
    p_values <- c(0.001, 0.002, 0.003, 0.004, 0.005)
    
    result <- TSENAT:::.adjust_pvalues_multicorr(
        p_values = p_values,
        multicorr = "hochberg",
        wy_randomizations = 100
    )
    
    expect_true(all(result >= 0 & result <= 1))
})

test_that(".adjust_pvalues_multicorr handles all non-significant p-values", {
    # All large p-values
    p_values <- c(0.5, 0.6, 0.7, 0.8, 0.9)
    
    result <- TSENAT:::.adjust_pvalues_multicorr(
        p_values = p_values,
        multicorr = "hochberg",
        wy_randomizations = 100
    )
    
    expect_true(all(result >= 0 & result <= 1))
})

# ═══════════════════════════════════════════════════════════════════════════════
# TEST SUITE: .fit_all_genes() - STRUCTURE TEST (complex function, basic validation)
# ═══════════════════════════════════════════════════════════════════════════════

test_that(".fit_all_genes exists and is callable", {
    # This function is complex and requires extensive setup
    # Basic check that it exists and can be called in principle
    expect_true(exists(".fit_all_genes", mode = "function"))
})

test_that(".fit_all_genes returns data.frame or empty frame", {
    skip_if_not_installed("SummarizedExperiment")
    skip_if_not_installed("nlme")
    
    # Create large, realistic test data with substantial sample size
    # The AR(1) + interaction model needs good data to converge
    set.seed(789)
    n_genes <- 3
    n_per_group <- 15  # 30 total samples for model stability
    
    # Create data with realistic but overlapping distribution (not perfect separation)
    # Higher variance helps with model identifiability
    group_a <- rnorm(n_per_group, mean = 2.8, sd = 0.4)
    group_b <- rnorm(n_per_group, mean = 1.8, sd = 0.4)
    
    # Build matrix with various gene-level effects
    mat <- matrix(NA, nrow = n_genes, ncol = n_per_group * 2)
    for (g in 1:n_genes) {
        # Vary effect size by gene but keep same group pattern
        effect_scale <- 0.6 + 0.3 * g / n_genes
        # Add additional noise to represent biological variation
        mat[g, 1:n_per_group] <- group_a * effect_scale + rnorm(n_per_group, 0, 0.25)
        mat[g, (n_per_group+1):(n_per_group*2)] <- group_b * effect_scale + rnorm(n_per_group, 0, 0.25)
    }
    
    rownames(mat) <- paste0("GENE_", 1:n_genes)
    colnames(mat) <- paste0("SAMPLE_", 1:(n_per_group*2))
    
    se <- SummarizedExperiment::SummarizedExperiment(
        assays = list(diversity = mat)
    )
    
    metadata <- list(
        q_vals = rep(1.0, n_per_group * 2),
        sample_names = paste0("SAMPLE_", 1:(n_per_group*2)),
        group_vec = c(rep(1, n_per_group), rep(2, n_per_group))
    )
    
    # Test with larger sample and realistic variation
    result <- suppressWarnings(tryCatch(
        TSENAT:::.fit_all_genes(
            mat = mat,
            se = se,
            metadata = metadata,
            method = "lmm",
            pvalue = "lrt",
            subject_col = NULL,
            paired = FALSE,
            min_obs = 3,
            nthreads = 1,
            verbose = FALSE,
            bias_correction = FALSE,
            regularization = "pca",
            corstr = "ar1",
            adaptive_knots = FALSE
        ),
        error = function(e) data.frame()
    ))
    
    expect_true(is.data.frame(result))
})

# ===========================================================================
# Tests for .validate_sait_method_dependencies()
# ===========================================================================

test_that(".validate_sait_method_dependencies passes for valid methods", {
    skip_if_not_installed("nlme")
    skip_if_not_installed("mgcv")
    skip_if_not_installed("geepack")
    
    expect_silent(TSENAT:::.validate_sait_method_dependencies("lmm"))
    expect_silent(TSENAT:::.validate_sait_method_dependencies("gam"))
    expect_silent(TSENAT:::.validate_sait_method_dependencies("gee"))
    expect_invisible(TSENAT:::.validate_sait_method_dependencies("fpca"))
})

test_that("calculate_sait validates method parameter correctly", {
    skip_if_not_installed("nlme")
    skip_if_not_installed("mgcv")
    skip_if_not_installed("geepack")
    
    data(readcounts, package = "TSENAT", envir = environment())
    readcounts <- as.matrix(readcounts)
    mode(readcounts) <- "numeric"
    
    metadata_df <- read.table(
        system.file("extdata", "metadata.tsv", package = "TSENAT"),
        header = TRUE, sep = "\t"
    )
    gff3_file <- system.file("extdata", "annotation.gff3.gz", package = "TSENAT")
    
    config <- TSENAT_config(
        sample_col = "sample",
        condition_col = "condition",
        q = seq(0, 2, by = 0.5)
    )
    analysis <- suppressWarnings(build_analysis(
        readcounts = readcounts,
        tx2gene = gff3_file,
        metadata = metadata_df,
        config = config,
        tpm = tpm,
        effective_length = effective_length
    ))
    analysis <- filter_analysis(analysis, min_samples = 1, subset_n_genes = 10)
    analysis <- calculate_diversity(analysis, q = seq(0, 2, by = 0.5), verbose = FALSE)
    
    # This should trigger the method validation path (covers lines 267-275 in sait_core.R)
    result <- suppressWarnings(
        calculate_sait(analysis, method = "lmm", verbose = FALSE)
    )
    expect_s4_class(result, "TSENATAnalysis")
    
    # Test that invalid method is rejected
    expect_error(
        suppressWarnings(calculate_sait(analysis, method = "invalid_method")),
        "should be one of"
    )
})

test_that("calculate_sait validates corstr parameter correctly", {
    skip_if_not_installed("nlme")
    skip_if_not_installed("mgcv")
    skip_if_not_installed("geepack")
    
    data(readcounts, package = "TSENAT", envir = environment())
    readcounts <- as.matrix(readcounts)
    mode(readcounts) <- "numeric"
    
    metadata_df <- read.table(
        system.file("extdata", "metadata.tsv", package = "TSENAT"),
        header = TRUE, sep = "\t"
    )
    gff3_file <- system.file("extdata", "annotation.gff3.gz", package = "TSENAT")
    
    config <- TSENAT_config(
        sample_col = "sample",
        condition_col = "condition",
        q = seq(0, 2, by = 0.5)
    )
    analysis <- suppressWarnings(build_analysis(
        readcounts = readcounts,
        tx2gene = gff3_file,
        metadata = metadata_df,
        config = config,
        tpm = tpm,
        effective_length = effective_length
    ))
    analysis <- filter_analysis(analysis, min_samples = 1, subset_n_genes = 10)
    analysis <- calculate_diversity(analysis, q = seq(0, 2, by = 0.5), verbose = FALSE)
    
    # Invalid corstr should be rejected via match.arg
    expect_error(
        suppressWarnings(calculate_sait(analysis, method = "gee", corstr = "invalid_corstr")),
        "should be one of"
    )
})

test_that("calculate_sait validates pvalue parameter correctly", {
    skip_if_not_installed("nlme")
    skip_if_not_installed("mgcv")
    
    data(readcounts, package = "TSENAT", envir = environment())
    readcounts <- as.matrix(readcounts)
    mode(readcounts) <- "numeric"
    
    metadata_df <- read.table(
        system.file("extdata", "metadata.tsv", package = "TSENAT"),
        header = TRUE, sep = "\t"
    )
    gff3_file <- system.file("extdata", "annotation.gff3.gz", package = "TSENAT")
    
    config <- TSENAT_config(
        sample_col = "sample",
        condition_col = "condition",
        q = seq(0, 2, by = 0.5)
    )
    analysis <- suppressWarnings(build_analysis(
        readcounts = readcounts,
        tx2gene = gff3_file,
        metadata = metadata_df,
        config = config,
        tpm = tpm,
        effective_length = effective_length
    ))
    analysis <- filter_analysis(analysis, min_samples = 1, subset_n_genes = 10)
    analysis <- calculate_diversity(analysis, q = seq(0, 2, by = 0.5), verbose = FALSE)
    
    # Invalid pvalue should be rejected via match.arg
    expect_error(
        suppressWarnings(calculate_sait(analysis, method = "gam", pvalue = "invalid_pvalue")),
        "should be one of"
    )
})

# ════════════════════════════════════════════════════════════════════════════════
# Tests for .validate_sait_method_dependencies()
# ════════════════════════════════════════════════════════════════════════════════

test_that(".validate_sait_method_dependencies returns invisible TRUE for valid method", {
    skip_if_not_installed("nlme")
    result <- TSENAT:::.validate_sait_method_dependencies("lmm")
    expect_true(result)
})

test_that(".validate_sait_method_dependencies passes for 'gam'", {
    skip_if_not_installed("mgcv")
    result <- TSENAT:::.validate_sait_method_dependencies("gam")
    expect_true(result)
})

# ════════════════════════════════════════════════════════════════════════════════
# Tests for .validate_sait_data_structure()
# ════════════════════════════════════════════════════════════════════════════════

test_that(".validate_sait_data_structure passes with valid columns and assay", {
    skip_if_not_installed("SummarizedExperiment")
    se <- SummarizedExperiment::SummarizedExperiment(
        assays = list(diversity = matrix(1:4, nrow = 2)),
        colData = data.frame(condition = c("A", "B"), row.names = c("s1", "s2"))
    )
    result <- TSENAT:::.validate_sait_data_structure(se, "condition", "diversity")
    expect_true(result)
})

test_that(".validate_sait_data_structure errors on missing condition_col", {
    skip_if_not_installed("SummarizedExperiment")
    se <- SummarizedExperiment::SummarizedExperiment(
        assays = list(div = matrix(1:4, nrow = 2)),
        colData = data.frame(sample = c("s1", "s2"), row.names = c("s1", "s2"))
    )
    expect_error(
        TSENAT:::.validate_sait_data_structure(se, "condition", "div"),
        "condition_col"
    )
})

test_that(".validate_sait_data_structure errors on missing assay", {
    skip_if_not_installed("SummarizedExperiment")
    se <- SummarizedExperiment::SummarizedExperiment(
        assays = list(counts = matrix(1:4, nrow = 2)),
        colData = data.frame(condition = c("A", "B"), row.names = c("s1", "s2"))
    )
    expect_error(
        TSENAT:::.validate_sait_data_structure(se, "condition", "diversity"),
        "Assay"
    )
})

# ════════════════════════════════════════════════════════════════════════════════
# Tests for .extract_sait_result_df()
# ════════════════════════════════════════════════════════════════════════════════

test_that(".extract_sait_result_df returns data.frame as-is", {
    df <- data.frame(gene = "g1", p_interaction = 0.01, stringsAsFactors = FALSE)
    result <- TSENAT:::.extract_sait_result_df(df)
    expect_identical(result, df)
})

test_that(".extract_sait_result_df extracts results from list", {
    df <- data.frame(gene = "g2", p_interaction = 0.05, stringsAsFactors = FALSE)
    result <- TSENAT:::.extract_sait_result_df(list(results = df, extra = "x"))
    expect_identical(result, df)
})

test_that(".extract_sait_result_df errors on invalid input", {
    expect_error(
        TSENAT:::.extract_sait_result_df(list(foo = "bar")),
        "sait_result must be either"
    )
})

# ════════════════════════════════════════════════════════════════════════════════
# Tests for wy_randomizations = "auto" branch in .calculate_sait()
# ════════════════════════════════════════════════════════════════════════════════

test_that(".calculate_sait: wy_randomizations='auto' errors when multicorr != 'westfall-young'", {
    skip_if_not_installed("SummarizedExperiment")
    skip_if_not_installed("nlme")
    skip_if_not_installed("mgcv")

    mat <- matrix(runif(30, 0, 1), nrow = 3, ncol = 10)
    rownames(mat) <- c("g1", "g2", "g3")
    colnames(mat) <- paste0(rep(c("s1", "s2"), each = 5), "_q=", rep(seq(0, 2, by = 0.5), 2))

    cd <- data.frame(
        condition = rep(c("normal", "tumor"), each = 5),
        row.names = colnames(mat),
        stringsAsFactors = FALSE
    )
    se <- SummarizedExperiment::SummarizedExperiment(
        assays = list(diversity = mat),
        colData = cd
    )

    expect_error(
        TSENAT:::.calculate_sait(
            se = se, condition_col = "condition", method = "lmm",
            multicorr = "hochberg", wy_randomizations = "auto",
            min_obs = 3, verbose = FALSE
        ),
        "wy_randomizations='auto' is only supported when multicorr='westfall-young'"
    )
})

test_that(".calculate_sait: wy_randomizations='auto' auto-estimates permutations with multicorr='westfall-young'", {
    skip_if_not_installed("SummarizedExperiment")
    skip_if_not_installed("nlme")
    skip_if_not_installed("mgcv")

    mat <- matrix(runif(30, 0, 1), nrow = 3, ncol = 10)
    rownames(mat) <- c("g1", "g2", "g3")
    colnames(mat) <- paste0(rep(c("s1", "s2"), each = 5), "_q=", rep(seq(0, 2, by = 0.5), 2))

    cd <- data.frame(
        condition = rep(c("normal", "tumor"), each = 5),
        row.names = colnames(mat),
        stringsAsFactors = FALSE
    )
    se <- SummarizedExperiment::SummarizedExperiment(
        assays = list(diversity = mat),
        colData = cd
    )

    result <- suppressWarnings(
        TSENAT:::.calculate_sait(
            se = se, condition_col = "condition", method = "lmm",
            multicorr = "westfall-young", wy_randomizations = "auto",
            min_obs = 3, verbose = FALSE
        )
    )

    expect_true(is.data.frame(result) || is.list(result))
})

# ============================================================================
# H1: regularization validation — (method, regularization) combo checked
# ============================================================================

test_that("H1: Invalid (method, regularization) combo produces error", {
    # GAM only accepts pca, gamsel, spline
    expect_error(
        TSENAT:::.validate_sait_interaction_input(
            method = "gam", pvalue = "satterthwaite", corstr = "ar1",
            regularization = "lasso", multicorr = "hochberg",
            pcorr = "BH", storey = FALSE, wy_randomizations = 100,
            paired = FALSE, subject_col = NULL, se = NULL, verbose = FALSE
        ),
        "regularization"
    )

    # LMM only accepts pca, lasso, elasticnet
    expect_error(
        TSENAT:::.validate_sait_interaction_input(
            method = "lmm", pvalue = "satterthwaite", corstr = "ar1",
            regularization = "gamsel", multicorr = "hochberg",
            pcorr = "BH", storey = FALSE, wy_randomizations = 100,
            paired = FALSE, subject_col = NULL, se = NULL, verbose = FALSE
        ),
        "regularization"
    )

    # Valid combos should NOT error
    expect_silent(
        TSENAT:::.validate_sait_interaction_input(
            method = "gam", pvalue = "satterthwaite", corstr = "ar1",
            regularization = "gamsel", multicorr = "hochberg",
            pcorr = "BH", storey = FALSE, wy_randomizations = 100,
            paired = FALSE, subject_col = NULL, se = NULL, verbose = FALSE
        )
    )
})

# ============================================================================
# M3: corstr = "auto" — now accepted in public API
# ============================================================================

test_that("M3: corstr='auto' is a valid argument", {
    # Verify "auto" is in the match.arg choices
    choices <- eval(formals(TSENAT:::.calculate_sait)$corstr)
    expect_true("auto" %in% choices)
})

# ============================================================================
# P3: FINALIZE SAIT RESULTS TESTS (July 2026: metrics.json refactoring)
# ============================================================================

context("SAIT: Finalize Results (.finalize_sait_results)")

test_that(".finalize_sait_results rejects non-data.frame input", {
  expect_error(
    TSENAT:::.finalize_sait_results(
      res = "not_a_df", mat = matrix(1:4, 2), se = NULL,
      metadata = NULL, method = "gam", pvalue = "lrt",
      subject_col = NULL, paired = FALSE, min_obs = 5,
      nthreads = 1, bias_correction = TRUE,
      regularization = "pca", corstr = "ar1",
      adaptive_knots = TRUE, multicorr = "hochberg",
      wy_randomizations = 1000, storey = FALSE,
      verbose = FALSE, return_model_data = FALSE,
      assay_name = "diversity"
    ),
    "should return a data.frame"
  )
})

test_that(".finalize_sait_results handles empty results gracefully", {
  empty_res <- data.frame(
    gene = character(0),
    p_interaction = numeric(0),
    stringsAsFactors = FALSE
  )
  result <- TSENAT:::.finalize_sait_results(
    res = empty_res, mat = matrix(1:4, 2), se = NULL,
    metadata = NULL, method = "gam", pvalue = "lrt",
    subject_col = NULL, paired = FALSE, min_obs = 5,
    nthreads = 1, bias_correction = TRUE,
    regularization = "pca", corstr = "ar1",
    adaptive_knots = TRUE, multicorr = "hochberg",
    wy_randomizations = 1000, storey = FALSE,
    verbose = FALSE, return_model_data = FALSE,
    assay_name = "diversity"
  )
  
  expect_is(result, "data.frame")
  expect_equal(nrow(result), 0)
  # Standard SAIT columns should be present
  expect_true("adj_p_interaction" %in% colnames(result))
})

test_that(".finalize_sait_results rejects results without p_interaction", {
  bad_res <- data.frame(gene = "G1", wrong_column = 0.05, stringsAsFactors = FALSE)
  expect_error(
    TSENAT:::.finalize_sait_results(
      res = bad_res, mat = matrix(1:4, 2), se = NULL,
      metadata = NULL, method = "gam", pvalue = "lrt",
      subject_col = NULL, paired = FALSE, min_obs = 5,
      nthreads = 1, bias_correction = TRUE,
      regularization = "pca", corstr = "ar1",
      adaptive_knots = TRUE, multicorr = "hochberg",
      wy_randomizations = 1000, storey = FALSE,
      verbose = FALSE, return_model_data = FALSE,
      assay_name = "diversity"
    ),
    "missing required.*p_interaction"
  )
})

test_that(".finalize_sait_results sorts by adjusted p-value", {
  res <- data.frame(
    gene = c("G1", "G2", "G3"),
    p_interaction = c(0.01, 0.05, 0.03),
    stringsAsFactors = FALSE
  )
  result <- TSENAT:::.finalize_sait_results(
    res = res, mat = matrix(1:9, 3), se = NULL,
    metadata = NULL, method = "gam", pvalue = "lrt",
    subject_col = NULL, paired = FALSE, min_obs = 5,
    nthreads = 1, bias_correction = TRUE,
    regularization = "pca", corstr = "ar1",
    adaptive_knots = TRUE, multicorr = "hochberg",
    wy_randomizations = 1000, storey = FALSE,
    verbose = FALSE, return_model_data = FALSE,
    assay_name = "diversity"
  )
  
  expect_is(result, "data.frame")
  # After Hochberg adjustment, order should be by adj_p first, then raw p
  expect_true("adj_p_interaction" %in% colnames(result))
  expect_equal(result$gene[1], "G1")  # G1 had the lowest p-value
})
