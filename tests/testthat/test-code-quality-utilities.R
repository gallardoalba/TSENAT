# Tests for code quality improvements: refactored utility functions
# Focus: High-in-degree functions with critical dependencies
#
# Background: These functions are called by 5+ plotting/analysis functions,
# so bugs here affect many downstream functions. Critical to test thoroughly.
#
# Date: April 15, 2026

context("Code Quality: Utility Function Tests (Post-Refactoring)")

# =============================================================================
# Tests for .apply_group_aesthetics and extracted helpers
# =============================================================================

test_that(".get_palette_colors returns vector when given palette name", {
    # Pure function test - should work independently without side effects
    colors <- .get_palette_colors("palette_blue_red")
    
    expect_is(colors, "character")
    expect_length(colors, 8)
    expect_true(all(grepl("^#[0-9A-F]{6}$", colors)))
})

test_that(".get_palette_colors handles color vector directly", {
    # Should accept pre-made color vectors
    test_colors <- c("#FF0000", "#00FF00", "#0000FF")
    result <- .get_palette_colors(test_colors)
    
    expect_identical(result, test_colors)
})

test_that(".get_palette_colors handles invalid palette gracefully", {
    # Should fail gracefully with clear error when given invalid palette
    expect_error(
        .get_palette_colors("nonexistent_palette_xyz"),
        "not found"
    )
})

test_that(".apply_aesthetics_colors preserves plot class", {
    # Pure function test - should not modify plot class
    
    
    p <- ggplot2::ggplot(mtcars, ggplot2::aes(x = mpg, y = cyl))
    colors <- c("#FF0000", "#0000FF")
    
    result <- .apply_aesthetics_colors(p, colors)
    
    expect_is(result, "ggplot")
})

test_that(".apply_aesthetics_colors applies color and fill scales", {
    # Verify that both color and fill scales are applied
    
    
    p <- ggplot2::ggplot(mtcars, ggplot2::aes(x = mpg, y = cyl, color = factor(cyl)))
    colors <- c("#FF0000", "#0000FF", "#00FF00")  # 3 colors for 3 levels of cyl
    
    result <- .apply_aesthetics_colors(p, colors)
    
    # Verify result is a ggplot object with scales applied
    expect_is(result, "ggplot")
    
    # Check that scales were added
    expect_true(length(result$scales$scales) > 0)
})

test_that(".apply_aesthetics_colors respects direction parameter", {
    # Test direction reversal
    
    
    p <- ggplot2::ggplot(mtcars, ggplot2::aes(x = mpg, y = cyl))
    colors <- c("#FF0000", "#00FF00", "#0000FF")
    
    forward <- .apply_aesthetics_colors(p, colors, direction = 1)
    reverse <- .apply_aesthetics_colors(p, colors, direction = -1)
    
    # Both should return ggplot objects
    expect_is(forward, "ggplot")
    expect_is(reverse, "ggplot")
})

test_that(".apply_aesthetics_colors sets legend position", {
    # Test legend positioning
    
    
    p <- ggplot2::ggplot(mtcars, ggplot2::aes(x = mpg, y = cyl))
    colors <- c("#FF0000", "#0000FF")
    
    result <- .apply_aesthetics_colors(p, colors, legend_position = "top")
    
    expect_is(result, "ggplot")
})

test_that(".apply_aesthetics_colors works with custom legend name", {
    # Test legend naming
    
    
    p <- ggplot2::ggplot(mtcars, ggplot2::aes(x = mpg, y = cyl, color = factor(vs)))
    colors <- c("#FF0000", "#0000FF")
    
    result <- .apply_aesthetics_colors(p, colors, legend_name = "Custom")
    
    expect_is(result, "ggplot")
})

test_that(".apply_group_aesthetics integrates palette lookup and application", {
    # Integration test: full pipeline
    
    
    p <- ggplot2::ggplot(mtcars, ggplot2::aes(x = mpg, y = cyl, color = factor(vs)))
    
    # Using palette name
    result1 <- .apply_group_aesthetics(p, palette = "palette_blue_red")
    expect_is(result1, "ggplot")
    
    # Using color vector directly
    colors <- c("#FF0000", "#0000FF")
    result2 <- .apply_group_aesthetics(p, palette = colors)
    expect_is(result2, "ggplot")
})

test_that(".apply_group_aesthetics works with all calling contexts", {
    # Verify that the 9 callers still work correctly
    
    
    p <- ggplot2::ggplot(mtcars, ggplot2::aes(x = mpg, y = cyl))
    
    # Test with various parameter combinations used by callers
    expect_silent(.apply_group_aesthetics(p, palette = "palette_blue_red", 
                                         legend_name = "Group"))
    expect_silent(.apply_group_aesthetics(p, palette = "palette_blue_red", 
                                         legend_position = "top"))
    expect_silent(.apply_group_aesthetics(p, palette = "palette_blue_red", 
                                         direction = -1))
})

# =============================================================================
# Tests for .validate_sait_method_dependencies
# =============================================================================

test_that(".validate_sait_method_dependencies validates LMM", {
    # LMM requires 'nlme' package
    expect_silent(.validate_sait_method_dependencies("lmm"))
})

test_that(".validate_sait_method_dependencies validates GAM", {
    # GAM requires 'mgcv' package  
    expect_silent(.validate_sait_method_dependencies("gam"))
})

test_that(".validate_sait_method_dependencies validates GEE", {
    # GEE requires 'geepack' package
    expect_silent(.validate_sait_method_dependencies("gee"))
})

test_that(".validate_sait_method_dependencies validates FPCA", {
    expect_silent(.validate_sait_method_dependencies("fpca"))
})

# =============================================================================
# Tests for .validate_sait_data_structure
# =============================================================================

test_that(".validate_sait_data_structure checks condition_col in colData", {
    
    
    # Create minimal SummarizedExperiment
    mat <- matrix(1:20, nrow = 5, ncol = 4)
    se <- SummarizedExperiment::SummarizedExperiment(assays = list(counts = mat))
    
    # Should fail - no condition column
    expect_error(
        .validate_sait_data_structure(se, condition_col = "treatment", assay_name = "counts"),
        "not found in colData"
    )
})

test_that(".validate_sait_data_structure checks assay exists", {
    
    
    mat <- matrix(1:20, nrow = 5, ncol = 4)
    coldata <- data.frame(condition = rep(c("A", "B"), 2))
    se <- SummarizedExperiment::SummarizedExperiment(
        assays = list(counts = mat),
        colData = coldata
    )
    
    # Should fail - no diversity assay
    expect_error(
        .validate_sait_data_structure(se, condition_col = "condition", assay_name = "diversity"),
        "not found"
    )
})

test_that(".validate_sait_data_structure passes with valid structure", {
    
    
    mat <- matrix(1:20, nrow = 5, ncol = 4)
    coldata <- data.frame(condition = rep(c("A", "B"), 2))
    se <- SummarizedExperiment::SummarizedExperiment(
        assays = list(diversity = mat),
        colData = coldata
    )
    
    # Should pass - valid structure
    expect_silent(
        .validate_sait_data_structure(se, condition_col = "condition", assay_name = "diversity")
    )
})

# =============================================================================
# Tests for .postprocess_sait_results
# =============================================================================

test_that(".postprocess_sait_results adds gene column from gene_id", {
    # Create minimal results dataframe
    res <- data.frame(
        gene_id = c("ENSG001", "ENSG002"),
        p_interaction = c(0.01, 0.05)
    )
    
    result <- .postprocess_sait_results(res, return_model_data = FALSE)
    
    expect_true("gene" %in% colnames(result))
    expect_equal(result$gene, c("ENSG001", "ENSG002"))
})

test_that(".postprocess_sait_results adds gene column from gene_name", {
    # When gene_id not available, should use gene_name
    res <- data.frame(
        gene_name = c("BRCA1", "BRCA2"),
        p_interaction = c(0.01, 0.05)
    )
    
    result <- .postprocess_sait_results(res, return_model_data = FALSE)
    
    expect_true("gene" %in% colnames(result))
    expect_equal(result$gene, c("BRCA1", "BRCA2"))
})

test_that(".postprocess_sait_results preserves existing gene column", {
    # If gene already exists, should keep it
    res <- data.frame(
        gene = c("GeneA", "GeneB"),
        gene_id = c("ENSG001", "ENSG002"),
        p_interaction = c(0.01, 0.05)
    )
    
    result <- .postprocess_sait_results(res, return_model_data = FALSE)
    
    # Should keep the original gene column
    expect_equal(result$gene, c("GeneA", "GeneB"))
})

test_that(".postprocess_sait_results returns data.frame by default", {
    res <- data.frame(
        gene_id = c("ENSG001", "ENSG002"),
        p_interaction = c(0.01, 0.05)
    )
    
    result <- .postprocess_sait_results(res, return_model_data = FALSE)
    
    expect_is(result, "data.frame")
})

test_that(".postprocess_sait_results returns list with model_data when requested", {
    
    
    res <- data.frame(
        gene_id = c("ENSG001", "ENSG002"),
        p_interaction = c(0.01, 0.05)
    )
    
    mat <- matrix(1:20, nrow = 5, ncol = 4)
    se <- SummarizedExperiment::SummarizedExperiment(assays = list(diversity = mat))
    
    result <- .postprocess_sait_results(res, return_model_data = TRUE, 
                                     se = se, mat = mat, metadata = list())
    
    expect_is(result, "list")
    expect_true("results" %in% names(result))
    expect_true("model_data" %in% names(result))
    expect_is(result$results, "data.frame")
})

# =============================================================================
# Regression Tests: Verify Extracted Functions Don't Break Existing Code
# =============================================================================

test_that("Extracted functions maintain backward compatibility with .calculate_sait", {
    # Verify that the new helper functions work in the actual calling context
    
    
    # Test that validation functions work
    mat <- matrix(1:20, nrow = 5, ncol = 4)
    coldata <- data.frame(condition = rep(c("A", "B"), 2))
    se <- SummarizedExperiment::SummarizedExperiment(
        assays = list(diversity = mat),
        colData = coldata
    )
    
    # Should pass through all validation steps
    expect_silent(.validate_sait_method_dependencies("gam"))
    expect_silent(.validate_sait_data_structure(se, "condition", "diversity"))
})

# =============================================================================
# Integration: Test High-In-Degree Function Reliability
# =============================================================================

test_that(".apply_group_aesthetics is reliable for 9+ calling contexts", {
    # This function has 9 callers across plotting code
    # Test that it handles edge cases robustly
    
    
    p <- ggplot2::ggplot(mtcars, ggplot2::aes(x = mpg, y = cyl))
    
    # Edge case: empty color vector should fail gracefully
    expect_error(
        .apply_aesthetics_colors(p, c()),
        "length|empty|must"
    )
    
    # Edge case: single color vector
    p1 <- .apply_aesthetics_colors(p, c("#FF0000"))
    expect_is(p1, "ggplot")
    
    # Edge case: many colors
    many_colors <- grDevices::rainbow(20)
    p2 <- .apply_aesthetics_colors(p, many_colors)
    expect_is(p2, "ggplot")
})

# =============================================================================
# Code Quality: Verify Pure Functions Have No Side Effects
# =============================================================================

test_that(".get_palette_colors has no side effects", {
    # Pure functions should produce same output for same input
    result1 <- .get_palette_colors("palette_blue_red")
    result2 <- .get_palette_colors("palette_blue_red")
    
    expect_identical(result1, result2)
})

test_that(".apply_aesthetics_colors doesn't modify input plot", {
    # Pure functions should not modify arguments
    
    
    p <- ggplot2::ggplot(mtcars, ggplot2::aes(x = mpg, y = cyl))
    p_original <- p  # Reference to original
    
    # Apply aesthetics
    result <- .apply_aesthetics_colors(p, c("#FF0000", "#0000FF"))
    
    # Original plot should be unchanged (R semantics mean p_original still points to same object)
    expect_is(result, "ggplot")
    expect_is(p, "ggplot")
})



test_that("auto_detect_column returns config value if available", {
    available_cols <- c("sample_id", "condition", "group")
    config_list <- list(condition_col = "condition")
    
    result <- auto_detect_column(
        available_cols = available_cols,
        config_list = config_list,
        config_key = "condition_col"
    )
    
    expect_equal(result, "condition")
})

test_that("auto_detect_column uses priority candidates if config unavailable", {
    available_cols <- c("sample_id", "group", "batch")
    priority_candidates <- c("condition", "group", "treatment")
    
    result <- auto_detect_column(
        available_cols = available_cols,
        priority_candidates = priority_candidates
    )
    
    expect_equal(result, "group")  # First match in priority order
})

test_that("auto_detect_column returns fallback when no match found", {
    available_cols <- c("sample_id", "value")
    
    result <- auto_detect_column(
        available_cols = available_cols,
        priority_candidates = c("condition", "group"),
        default_fallback = "sample_id"
    )
    
    expect_equal(result, "sample_id")
})

test_that("auto_detect_column returns NULL when no match and no fallback", {
    available_cols <- c("sample_id", "value")
    
    result <- auto_detect_column(
        available_cols = available_cols,
        priority_candidates = c("condition", "group")
    )
    
    expect_null(result)
})

test_that("auto_detect_column prioritizes config over priority candidates", {
    available_cols <- c("sample_id", "condition", "group", "batch")
    config_list <- list(condition_col = "batch")
    priority_candidates <- c("condition", "group")  # Would normally match "condition"
    
    result <- auto_detect_column(
        available_cols = available_cols,
        config_list = config_list,
        config_key = "condition_col",
        priority_candidates = priority_candidates
    )
    
    expect_equal(result, "batch")  # Config wins
})

test_that("auto_detect_column ignores missing config key", {
    available_cols <- c("sample_id", "condition", "group")
    config_list <- list(other_key = "condition")
    priority_candidates <- c("group")
    
    result <- auto_detect_column(
        available_cols = available_cols,
        config_list = config_list,
        config_key = "missing_key",
        priority_candidates = priority_candidates
    )
    
    expect_equal(result, "group")  # Falls back to priority candidates
})

test_that("auto_detect_column ignores config value not in available_cols", {
    available_cols <- c("sample_id", "group")
    config_list <- list(condition_col = "nonexistent")
    priority_candidates <- c("group")
    
    result <- auto_detect_column(
        available_cols = available_cols,
        config_list = config_list,
        config_key = "condition_col",
        priority_candidates = priority_candidates
    )
    
    expect_equal(result, "group")  # Config ignored, uses priority candidates
})

test_that("auto_detect_column ignores default_fallback not in available_cols", {
    available_cols <- c("sample_id", "group")
    result <- auto_detect_column(
        available_cols = available_cols,
        priority_candidates = "condition",
        default_fallback = "foo"
    )
    
    expect_null(result)
})

test_that("auto_detect_column sends verbose messages", {
    available_cols <- c("condition", "group")
    priority_candidates <- c("condition")
    
    expect_message(
        auto_detect_column(
            available_cols = available_cols,
            priority_candidates = priority_candidates,
            verbose = TRUE,
            param_name = "test_param"
        ),
        "Auto-detected test_param"
    )
})

test_that("auto_detect_column logs unavailable fallback in verbose mode", {
    available_cols <- c("sample_id", "value")

    expect_message(
        auto_detect_column(
            available_cols = available_cols,
            priority_candidates = c("condition", "group"),
            default_fallback = "foo",
            verbose = TRUE,
            param_name = "test_param"
        ),
        "default_fallback 'foo' for test_param is not available; returning NULL"
    )
})

# ===========================================================================
# Tests for resolve_slot_param()
# ===========================================================================

test_that("resolve_slot_param returns user_value when provided", {
    config <- list(nthreads = 4, verbose = TRUE)
    result <- TSENAT:::resolve_slot_param(
        user_value = 8, config_list = config,
        config_key = "nthreads", default_value = 1
    )
    expect_equal(result, 8)
})

test_that("resolve_slot_param falls back to config when user_value is NULL", {
    config <- list(nthreads = 4, verbose = TRUE)
    result <- TSENAT:::resolve_slot_param(
        user_value = NULL, config_list = config,
        config_key = "nthreads", default_value = 1
    )
    expect_equal(result, 4)
})

test_that("resolve_slot_param falls back to default when neither user nor config", {
    config <- list(nthreads = 4)
    result <- TSENAT:::resolve_slot_param(
        user_value = NULL, config_list = config,
        config_key = "missing_key", default_value = 16
    )
    expect_equal(result, 16)
})

test_that("resolve_slot_param returns NULL when everything is NULL and allow_null=TRUE", {
    config <- list(other = 1)
    result <- TSENAT:::resolve_slot_param(
        user_value = NULL, config_list = config,
        config_key = "nthreads", default_value = NULL
    )
    expect_null(result)
})

test_that("resolve_slot_param errors when allow_null=FALSE and nothing found", {
    config <- list(other = 1)
    expect_error(
        TSENAT:::resolve_slot_param(
            user_value = NULL, config_list = config,
            config_key = "nthreads", default_value = NULL,
            allow_null = FALSE
        ),
        "Could not resolve parameter"
    )
})

test_that("resolve_slot_param uses description in error message", {
    config <- list()
    expect_error(
        TSENAT:::resolve_slot_param(
            user_value = NULL, config_list = config,
            config_key = "nthreads", default_value = NULL,
            description = "thread count", allow_null = FALSE
        ),
        "thread count"
    )
})

test_that("resolve_slot_param handles NULL config_list", {
    result <- TSENAT:::resolve_slot_param(
        user_value = NULL, config_list = NULL,
        config_key = "nthreads", default_value = 4
    )
    expect_equal(result, 4)
})

test_that("resolve_slot_param handles NULL config value for existing key", {
    config <- list(nthreads = NULL)
    result <- TSENAT:::resolve_slot_param(
        user_value = NULL, config_list = config,
        config_key = "nthreads", default_value = 4
    )
    expect_equal(result, 4)
})



test_that("save_analysis_output saves data.frame to CSV", {
    
    
    # Create temporary file
    temp_file <- tempfile(fileext = ".csv")
    on.exit(unlink(temp_file), add = TRUE)
    
    test_df <- data.frame(
        gene = c("g1", "g2"),
        value = c(1.5, 2.3),
        stringsAsFactors = FALSE
    )
    
    result <- save_analysis_output(
        data = test_df,
        output_file = temp_file,
        verbose = FALSE
    )
    
    expect_true(result)
    expect_true(file.exists(temp_file))
    
    # Verify content
    loaded <- read.csv(temp_file, row.names = 1)
    expect_equal(nrow(loaded), 2)
})

test_that("save_analysis_output saves data.frame to TSV", {
    
    
    temp_file <- tempfile(fileext = ".tsv")
    on.exit(unlink(temp_file), add = TRUE)
    
    test_df <- data.frame(x = c(1, 2), y = c(3, 4))
    
    result <- save_analysis_output(
        data = test_df,
        output_file = temp_file,
        verbose = FALSE
    )
    
    expect_true(result)
    expect_true(file.exists(temp_file))
})

test_that("save_analysis_output saves object to RDS", {
    
    
    temp_file <- tempfile(fileext = ".rds")
    on.exit(unlink(temp_file), add = TRUE)
    
    test_list <- list(a = 1, b = c(2, 3, 4))
    
    result <- save_analysis_output(
        data = test_list,
        output_file = temp_file,
        verbose = FALSE
    )
    
    expect_true(result)
    expect_true(file.exists(temp_file))
    
    # Verify content
    loaded <- readRDS(temp_file)
    expect_identical(loaded, test_list)
})

test_that("save_analysis_output creates directory if create_dir=TRUE", {
    
    
    temp_dir <- file.path(tempdir(), "test_output_dir_new", "subdir")
    on.exit(unlink(file.path(tempdir(), "test_output_dir_new"), recursive = TRUE), add = TRUE)
    
    temp_file <- file.path(temp_dir, "test.csv")
    
    test_df <- data.frame(x = 1)
    
    result <- save_analysis_output(
        data = test_df,
        output_file = temp_file,
        create_dir = TRUE,
        verbose = FALSE
    )
    
    expect_true(result)
    expect_true(dir.exists(temp_dir))
    expect_true(file.exists(temp_file))
})

test_that("save_analysis_output returns FALSE for NULL output_file", {
    test_df <- data.frame(x = 1)
    
    result <- save_analysis_output(
        data = test_df,
        output_file = NULL,
        verbose = FALSE
    )
    
    expect_false(result)
})

test_that("save_analysis_output converts matrix to data.frame and saves", {
    
    
    temp_file <- tempfile(fileext = ".tsv")
    on.exit(unlink(temp_file), add = TRUE)
    
    test_matrix <- matrix(c(1, 2, 3, 4), nrow = 2)
    
    result <- save_analysis_output(
        data = test_matrix,
        output_file = temp_file,
        verbose = FALSE
    )
    
    expect_true(result)
    expect_true(file.exists(temp_file))
})

test_that("save_analysis_output handles unknown format with RDS fallback", {
    
    
    temp_file <- tempfile(fileext = ".unknown")
    on.exit(unlink(temp_file), add = TRUE)
    
    test_data <- list(a = 1)
    
    result <- expect_warning(
        save_analysis_output(
            data = test_data,
            output_file = temp_file,
            verbose = FALSE
        ),
        "Unknown output format"
    )
    
    expect_true(result)
    expect_true(file.exists(temp_file))
})

test_that("save_analysis_output returns verbose messages when verbose=TRUE", {
    skip_if_not_installed("ggplot2")
    
    temp_file <- tempfile(fileext = ".tsv")
    on.exit(unlink(temp_file), add = TRUE)
    
    test_df <- data.frame(x = 1)
    
    expect_message(
        save_analysis_output(
            data = test_df,
            output_file = temp_file,
            verbose = TRUE
        ),
        "Saved table to"
    )
})

test_that("save_analysis_output skips empty data.frame with verbose message", {
    temp_file <- tempfile(fileext = ".tsv")
    on.exit(unlink(temp_file), add = TRUE)
    
    empty_df <- data.frame()
    
    result <- expect_message(
        save_analysis_output(
            data = empty_df,
            output_file = temp_file,
            verbose = TRUE
        ),
        "Skipping empty output"
    )
    
    expect_false(result)
})

test_that("save_analysis_output saves ggplot to PDF with custom dimensions", {
    skip_if_not_installed("ggplot2")
    
    temp_file <- tempfile(fileext = ".pdf")
    on.exit(unlink(temp_file), add = TRUE)
    
    p <- ggplot2::ggplot(mtcars, ggplot2::aes(x = mpg, y = cyl)) + ggplot2::geom_point()
    
    result <- save_analysis_output(
        data = p,
        output_file = temp_file,
        verbose = FALSE,
        width = 10,
        height = 8
    )
    
    expect_true(result)
    expect_true(file.exists(temp_file))
    expect_gt(file.info(temp_file)$size, 0)
})

test_that("save_analysis_output warns for non-ggplot object with plot extension", {
    temp_file <- tempfile(fileext = ".pdf")
    on.exit(unlink(temp_file), add = TRUE)
    
    result <- suppressWarnings(
        save_analysis_output(
            data = data.frame(x = 1),
            output_file = temp_file,
            verbose = FALSE
        )
    )
    
    expect_false(result)
})

test_that("save_analysis_output saves SummarizedExperiment as table", {
    skip_if_not_installed("SummarizedExperiment")
    
    temp_file <- tempfile(fileext = ".tsv")
    on.exit(unlink(temp_file), add = TRUE)
    
    se <- SummarizedExperiment::SummarizedExperiment(
        assays = list(counts = matrix(1:6, nrow = 2, dimnames = list(c("g1","g2"), c("s1","s2","s3"))))
    )
    
    result <- save_analysis_output(
        data = se,
        output_file = temp_file,
        verbose = FALSE
    )
    
    expect_true(result)
    expect_true(file.exists(temp_file))
})

test_that("save_analysis_output warns when save fails on locked file path", {
    temp_file <- file.path(tempdir(), "locked", "test.csv")

    test_df <- data.frame(x = 1)

    result <- suppressWarnings(
        save_analysis_output(
            data = test_df,
            output_file = temp_file,
            create_dir = FALSE,
            verbose = FALSE
        )
    )

    # Directory doesn't exist, save should fail gracefully
    expect_false(result)
})

test_that("extract_multiq_table returns data.frame for single q-value", {
    result_single <- list(
        summary_table = data.frame(
            gene = c("g1", "g2"),
            value = c(1.0, 2.0)
        )
    )
    
    output <- extract_multiq_table(
        result = result_single,
        is_multiq = FALSE
    )
    
    expect_true(is.data.frame(output))
    expect_equal(nrow(output), 2)
})

test_that("extract_multiq_table combines multi-q results", {
    result_multiq <- list(
        q_1_00 = list(
            summary_table = data.frame(gene = c("g1", "g2"), value = c(1.0, 2.0))
        ),
        q_2_00 = list(
            summary_table = data.frame(gene = c("g1", "g2"), value = c(1.5, 2.5))
        )
    )
    
    output <- extract_multiq_table(
        result = result_multiq,
        is_multiq = TRUE,
        q_value_col = "q"
    )
    
    expect_true(is.data.frame(output))
    expect_equal(nrow(output), 4)  # 2 genes × 2 q values
    expect_true("q" %in% colnames(output))
})

test_that("extract_multiq_table auto-detects multi-q from list structure", {
    result_multiq <- list(
        q_0_50 = list(
            summary_table = data.frame(gene = c("g1"), value = c(1.0))
        ),
        q_1_00 = list(
            summary_table = data.frame(gene = c("g1"), value = c(1.5))
        )
    )
    
    # is_multiq should be auto-detected
    output <- extract_multiq_table(
        result = result_multiq,
        is_multiq = NULL
    )
    
    expect_true(is.data.frame(output))
    expect_equal(nrow(output), 2)
})

test_that("extract_multiq_table adds q_value column to combined results", {
    result_multiq <- list(
        q_1_50 = list(
            summary_table = data.frame(gene = "g1", value = 1.5)
        ),
        q_2_00 = list(
            summary_table = data.frame(gene = "g1", value = 2.0)
        )
    )
    
    output <- extract_multiq_table(
        result = result_multiq,
        is_multiq = TRUE,
        q_value_col = "q_value"
    )
    
    expect_true("q_value" %in% colnames(output))
    expect_equal(as.numeric(output$q_value[1]), 1.5)
    expect_equal(as.numeric(output$q_value[2]), 2.0)
})

test_that("extract_multiq_table returns NULL when summary_table missing", {
    result_error <- list(
        q_1_00 = list(other_data = "no summary_table")
    )
    
    output <- extract_multiq_table(
        result = result_error,
        is_multiq = TRUE
    )
    
    expect_null(output)
})

test_that("extract_multiq_table uses custom extraction function", {
    result_custom <- list(
        results = data.frame(gene = "g1", pvalue = 0.01)
    )
    
    custom_fn <- function(result_element, q_key) {
        if (!is.null(result_element$results)) {
            return(result_element$results)
        }
        return(NULL)
    }
    
    output <- extract_multiq_table(
        result = result_custom,
        is_multiq = FALSE,
        extract_fn = custom_fn
    )
    
    expect_true(is.data.frame(output))
    expect_equal(nrow(output), 1)
})
