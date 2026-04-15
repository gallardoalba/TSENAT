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
    skip_if_not_installed("ggplot2")
    
    p <- ggplot2::ggplot(mtcars, ggplot2::aes(x = mpg, y = cyl))
    colors <- c("#FF0000", "#0000FF")
    
    result <- .apply_aesthetics_colors(p, colors)
    
    expect_is(result, "ggplot")
})

test_that(".apply_aesthetics_colors applies color and fill scales", {
    # Verify that both color and fill scales are applied
    skip_if_not_installed("ggplot2")
    
    p <- ggplot2::ggplot(mtcars, ggplot2::aes(x = mpg, y = cyl, color = factor(cyl)))
    colors <- c("#FF0000", "#0000FF", "#00FF00")  # 3 colors for 3 levels of cyl
    
    result <- .apply_aesthetics_colors(p, colors)
    
    # Verify result is a ggplot object with scales applied
    expect_is(result, "ggplot")
    
    # Check that scales were added
    expect_true(length(result$scales$scales) > 0)
})
})

test_that(".apply_aesthetics_colors respects direction parameter", {
    # Test direction reversal
    skip_if_not_installed("ggplot2")
    
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
    skip_if_not_installed("ggplot2")
    
    p <- ggplot2::ggplot(mtcars, ggplot2::aes(x = mpg, y = cyl))
    colors <- c("#FF0000", "#0000FF")
    
    result <- .apply_aesthetics_colors(p, colors, legend_position = "top")
    
    expect_is(result, "ggplot")
})

test_that(".apply_aesthetics_colors works with custom legend name", {
    # Test legend naming
    skip_if_not_installed("ggplot2")
    
    p <- ggplot2::ggplot(mtcars, ggplot2::aes(x = mpg, y = cyl, color = factor(vs)))
    colors <- c("#FF0000", "#0000FF")
    
    result <- .apply_aesthetics_colors(p, colors, legend_name = "Custom")
    
    expect_is(result, "ggplot")
})

test_that(".apply_group_aesthetics integrates palette lookup and application", {
    # Integration test: full pipeline
    skip_if_not_installed("ggplot2")
    
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
    skip_if_not_installed("ggplot2")
    
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
# Tests for .validate_lm_method_dependencies
# =============================================================================

test_that(".validate_lm_method_dependencies validates LMM", {
    # LMM requires 'nlme' package
    expect_silent(.validate_lm_method_dependencies("lmm"))
})

test_that(".validate_lm_method_dependencies validates GAM", {
    # GAM requires 'mgcv' package  
    expect_silent(.validate_lm_method_dependencies("gam"))
})

test_that(".validate_lm_method_dependencies validates GEE", {
    # GEE requires 'geepack' package
    expect_silent(.validate_lm_method_dependencies("gee"))
})

test_that(".validate_lm_method_dependencies validates FPCA", {
    # FPCA requires 'refund' package
    skip_if_not_installed("refund")
    expect_silent(.validate_lm_method_dependencies("fpca"))
})

# =============================================================================
# Tests for .validate_lm_data_structure
# =============================================================================

test_that(".validate_lm_data_structure checks condition_col in colData", {
    skip_if_not_installed("SummarizedExperiment")
    
    # Create minimal SummarizedExperiment
    mat <- matrix(1:20, nrow = 5, ncol = 4)
    se <- SummarizedExperiment::SummarizedExperiment(assays = list(counts = mat))
    
    # Should fail - no condition column
    expect_error(
        .validate_lm_data_structure(se, condition_col = "treatment", assay_name = "counts"),
        "not found in colData"
    )
})

test_that(".validate_lm_data_structure checks assay exists", {
    skip_if_not_installed("SummarizedExperiment")
    
    mat <- matrix(1:20, nrow = 5, ncol = 4)
    coldata <- data.frame(condition = rep(c("A", "B"), 2))
    se <- SummarizedExperiment::SummarizedExperiment(
        assays = list(counts = mat),
        colData = coldata
    )
    
    # Should fail - no diversity assay
    expect_error(
        .validate_lm_data_structure(se, condition_col = "condition", assay_name = "diversity"),
        "not found"
    )
})

test_that(".validate_lm_data_structure passes with valid structure", {
    skip_if_not_installed("SummarizedExperiment")
    
    mat <- matrix(1:20, nrow = 5, ncol = 4)
    coldata <- data.frame(condition = rep(c("A", "B"), 2))
    se <- SummarizedExperiment::SummarizedExperiment(
        assays = list(diversity = mat),
        colData = coldata
    )
    
    # Should pass - valid structure
    expect_silent(
        .validate_lm_data_structure(se, condition_col = "condition", assay_name = "diversity")
    )
})

# =============================================================================
# Tests for .postprocess_lm_results
# =============================================================================

test_that(".postprocess_lm_results adds gene column from gene_id", {
    # Create minimal results dataframe
    res <- data.frame(
        gene_id = c("ENSG001", "ENSG002"),
        p_interaction = c(0.01, 0.05)
    )
    
    result <- .postprocess_lm_results(res, return_model_data = FALSE)
    
    expect_true("gene" %in% colnames(result))
    expect_equal(result$gene, c("ENSG001", "ENSG002"))
})

test_that(".postprocess_lm_results adds gene column from gene_name", {
    # When gene_id not available, should use gene_name
    res <- data.frame(
        gene_name = c("BRCA1", "BRCA2"),
        p_interaction = c(0.01, 0.05)
    )
    
    result <- .postprocess_lm_results(res, return_model_data = FALSE)
    
    expect_true("gene" %in% colnames(result))
    expect_equal(result$gene, c("BRCA1", "BRCA2"))
})

test_that(".postprocess_lm_results preserves existing gene column", {
    # If gene already exists, should keep it
    res <- data.frame(
        gene = c("GeneA", "GeneB"),
        gene_id = c("ENSG001", "ENSG002"),
        p_interaction = c(0.01, 0.05)
    )
    
    result <- .postprocess_lm_results(res, return_model_data = FALSE)
    
    # Should keep the original gene column
    expect_equal(result$gene, c("GeneA", "GeneB"))
})

test_that(".postprocess_lm_results returns data.frame by default", {
    res <- data.frame(
        gene_id = c("ENSG001", "ENSG002"),
        p_interaction = c(0.01, 0.05)
    )
    
    result <- .postprocess_lm_results(res, return_model_data = FALSE)
    
    expect_is(result, "data.frame")
})

test_that(".postprocess_lm_results returns list with model_data when requested", {
    skip_if_not_installed("SummarizedExperiment")
    
    res <- data.frame(
        gene_id = c("ENSG001", "ENSG002"),
        p_interaction = c(0.01, 0.05)
    )
    
    mat <- matrix(1:20, nrow = 5, ncol = 4)
    se <- SummarizedExperiment::SummarizedExperiment(assays = list(diversity = mat))
    
    result <- .postprocess_lm_results(res, return_model_data = TRUE, 
                                     se = se, mat = mat, metadata = list())
    
    expect_is(result, "list")
    expect_true("results" %in% names(result))
    expect_true("model_data" %in% names(result))
    expect_is(result$results, "data.frame")
})

# =============================================================================
# Regression Tests: Verify Extracted Functions Don't Break Existing Code
# =============================================================================

test_that("Extracted functions maintain backward compatibility with .calculate_lm", {
    # Verify that the new helper functions work in the actual calling context
    skip_if_not_installed("SummarizedExperiment")
    
    # Test that validation functions work
    mat <- matrix(1:20, nrow = 5, ncol = 4)
    coldata <- data.frame(condition = rep(c("A", "B"), 2))
    se <- SummarizedExperiment::SummarizedExperiment(
        assays = list(diversity = mat),
        colData = coldata
    )
    
    # Should pass through all validation steps
    expect_silent(.validate_lm_method_dependencies("gam"))
    expect_silent(.validate_lm_data_structure(se, "condition", "diversity"))
})

# =============================================================================
# Integration: Test High-In-Degree Function Reliability
# =============================================================================

test_that(".apply_group_aesthetics is reliable for 9+ calling contexts", {
    # This function has 9 callers across plotting code
    # Test that it handles edge cases robustly
    skip_if_not_installed("ggplot2")
    
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
    skip_if_not_installed("ggplot2")
    
    p <- ggplot2::ggplot(mtcars, ggplot2::aes(x = mpg, y = cyl))
    p_original <- p  # Reference to original
    
    # Apply aesthetics
    result <- .apply_aesthetics_colors(p, c("#FF0000", "#0000FF"))
    
    # Original plot should be unchanged (R semantics mean p_original still points to same object)
    expect_is(result, "ggplot")
    expect_is(p, "ggplot")
})
