library(testthat)

context("LM Core Functions: Input Validation")

test_that(".tsenat_validate_lm_interaction_input rejects invalid storey", {
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
        TSENAT:::.tsenat_validate_lm_interaction_input(
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

test_that(".tsenat_validate_lm_interaction_input rejects invalid wy_randomizations", {
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
        TSENAT:::.tsenat_validate_lm_interaction_input(
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

test_that(".tsenat_validate_lm_interaction_input warns on low wy_randomizations", {
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
        TSENAT:::.tsenat_validate_lm_interaction_input(
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

test_that(".tsenat_validate_lm_interaction_input auto-detects subject_col when paired", {
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
    result <- TSENAT:::.tsenat_validate_lm_interaction_input(
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

context("LM Core Functions: Sample Metadata Parsing")

test_that(".tsenat_parse_sample_metadata extracts q-values correctly", {
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

    metadata <- TSENAT:::.tsenat_parse_sample_metadata(
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

test_that(".tsenat_parse_sample_metadata rejects missing _q=", {
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
        TSENAT:::.tsenat_parse_sample_metadata(
            se = se,
            condition_col = "condition",
            assay_name = "diversity",
            verbose = FALSE
        ),
        "Could not parse q values",
        fixed = FALSE
    )
})

test_that(".tsenat_parse_sample_metadata rejects missing condition_col", {
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
        TSENAT:::.tsenat_parse_sample_metadata(
            se = se,
            condition_col = "nonexistent_col",
            assay_name = "diversity",
            verbose = FALSE
        ),
        "No sample grouping found",
        fixed = FALSE
    )
})

context("LM Core Functions: Gene Annotation Mapping")

test_that(".tsenat_map_gene_annotations adds gene names from rowData", {
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

    mapped <- TSENAT:::.tsenat_map_gene_annotations(
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

test_that(".tsenat_map_gene_annotations handles missing gene_name gracefully", {
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

    mapped <- TSENAT:::.tsenat_map_gene_annotations(
        res = results,
        se = se,
        verbose = FALSE
    )

    # Should still have gene_name column (using gene as fallback)
    expect_true("gene_name" %in% colnames(mapped))
    expect_equal(mapped$gene_name, c("g1", "g2"))
})

context("LM Core Functions: Model Metadata Assembly")

test_that(".tsenat_assemble_model_metadata returns required fields", {
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

    metadata <- TSENAT:::.tsenat_parse_sample_metadata(
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

    model_data <- TSENAT:::.tsenat_assemble_model_metadata(
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

test_that(".tsenat_assemble_model_metadata computes per-group statistics", {
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

    metadata <- TSENAT:::.tsenat_parse_sample_metadata(
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

    model_data <- TSENAT:::.tsenat_assemble_model_metadata(
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

context("LM Refactored Main Function Integration")

test_that("calculate_lm_interaction works with refactored code", {
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
    results <- .calculate_lm_interaction(
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

test_that("calculate_lm_interaction returns model_data when requested", {
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
    output <- .calculate_lm_interaction(
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
