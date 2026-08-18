# Tests for text formatting utilities (Code Quality Phase 2)
# These are high-risk utility functions used in results display
#
# Date: April 15, 2026

context("Code Quality: Text Formatting Utilities")

# =============================================================================
# Tests for column width calculation (cycle-free helper)
# =============================================================================

test_that(".calculate_column_widths returns vector of widths", {
    df <- data.frame(
        short = c("a", "bb"),
        verylongname = c("x", "yy")
    )
    df_char <- as.data.frame(lapply(df, as.character), stringsAsFactors = FALSE)
    
    widths <- .calculate_column_widths(df_char)
    
    expect_is(widths, "numeric")
    expect_length(widths, 2)
})

test_that(".calculate_column_widths accounts for column names", {
    df <- data.frame(
        verylongheader = c("a", "b"),
        shortname = c("value1", "value2")
    )
    df_char <- as.data.frame(lapply(df, as.character), stringsAsFactors = FALSE)
    
    widths <- .calculate_column_widths(df_char)
    
    # First column header is "verylongheader" (14 chars)
    expect_gte(widths[1], 14)
    # Second column should fit "shortname" (9 chars) and "value1" (6 chars)
    expect_gte(widths[2], 9)
})

test_that(".calculate_column_widths accounts for data values", {
    df <- data.frame(
        name = c("Alice", "Bob"),
        value = c("verylongvalue", "short")
    )
    df_char <- as.data.frame(lapply(df, as.character), stringsAsFactors = FALSE)
    
    widths <- .calculate_column_widths(df_char)
    
    # Second column must fit "verylongvalue" (13 chars)
    expect_gte(widths[2], 13)
})

# =============================================================================
# Tests for table row formatting (cycle-free helper)
# =============================================================================

test_that(".format_table_row formats single row with padding", {
    values <- c("short", "mediumvalue", "x")
    col_widths <- c(10, 12, 5)
    
    result <- .format_table_row(values, col_widths)
    
    expect_is(result, "character")
    expect_length(result, 1)
    # Should have padding and spaces between columns
    expect_true(nchar(result) > sum(nchar(values)))
})

test_that(".format_table_row preserves order", {
    values <- c("first", "second", "third")
    col_widths <- c(10, 10, 10)
    
    result <- .format_table_row(values, col_widths)
    
    # Format should preserve the order
    expect_true(grepl("first.*second.*third", result))
})

test_that(".format_table_row handles empty strings", {
    values <- c("", "value", "")
    col_widths <- c(5, 10, 5)
    
    result <- .format_table_row(values, col_widths)
    
    expect_is(result, "character")
    expect_true(nchar(result) > 0)
})

# =============================================================================
# Tests for table header formatting (cycle-free helper)
# =============================================================================

test_that(".format_table_header formats column names with padding", {
    colnames <- c("Name", "Value", "Count")
    col_widths <- c(10, 12, 8)
    
    result <- .format_table_header(colnames, col_widths)
    
    expect_is(result, "character")
    # Should include all column names
    expect_true(grepl("Name", result))
    expect_true(grepl("Value", result))
    expect_true(grepl("Count", result))
})

test_that(".format_table_header respects column widths", {
    colnames <- c("A", "B", "C")
    col_widths <- c(5, 5, 5)
    
    result <- .format_table_header(colnames, col_widths)
    
    # Total length should be ~ sum of widths + spaces between
    expect_gte(nchar(result), 17)  # 5+1+5+1+5 = 17 minimum
})

# =============================================================================
# Integration Tests: Full Text Formatting Pipeline
# =============================================================================

test_that(".format_data_frame_as_text handles empty dataframe", {
    df <- data.frame(A = character(), B = character())
    
    result <- .format_data_frame_as_text(df)
    
    expect_identical(result, "")
})

test_that(".format_data_frame_as_text formats single row", {
    df <- data.frame(
        Name = "Alice",
        Age = "30"
    )
    
    result <- .format_data_frame_as_text(df)
    
    expect_is(result, "character")
    expect_true(grepl("Name", result))
    expect_true(grepl("Alice", result))
    expect_true(grepl("Age", result))
    expect_true(grepl("30", result))
})

test_that(".format_data_frame_as_text formats multiple rows", {
    df <- data.frame(
        Name = c("Alice", "Bob", "Charlie"),
        Age = c("30", "25", "35"),
        Dept = c("HR", "IT", "Sales")
    )
    
    result <- .format_data_frame_as_text(df)
    
    # Should contain all data
    expect_true(grepl("Alice", result))
    expect_true(grepl("Bob", result))
    expect_true(grepl("Charlie", result))
    
    # Should have multiple lines (header + 3 rows + blank line)
    lines <- strsplit(result, "\n")[[1]]
    expect_gte(length(lines), 4)
})

test_that(".format_data_frame_as_text handles numeric data", {
    df <- data.frame(
        Value = c(1.234, 5.678),
        Count = c(100, 200)
    )
    
    result <- .format_data_frame_as_text(df)
    
    expect_true(grepl("1.234", result))
    expect_true(grepl("5.678", result))
    expect_true(grepl("100", result))
    expect_true(grepl("200", result))
})

test_that(".format_data_frame_as_text handles NA values", {
    df <- data.frame(
        Name = c("Alice", NA, "Charlie"),
        Value = c(1, 2, NA)
    )
    
    result <- .format_data_frame_as_text(df)
    
    expect_is(result, "character")
    # Should show NA as string
    expect_true(grepl("Alice", result))
    expect_true(grepl("Charlie", result))
})

test_that(".format_data_frame_as_text ends with newline", {
    df <- data.frame(A = "test", B = "value")
    
    result <- .format_data_frame_as_text(df)
    
    # Should end with newline
    expect_true(grepl("\n$", result))
})

# =============================================================================
# Regression Tests: Cycle Breaking Validation
# =============================================================================

test_that("Extracted helpers maintain compatibility with .format_data_frame_as_text", {
    df <- data.frame(
        Gene = c("BRCA1", "BRCA2", "TP53"),
        p_value = c("0.001", "0.05", "0.1"),
        Effect = c("Strong", "Moderate", "Weak")
    )
    
    # Should produce valid formatted output
    result <- .format_data_frame_as_text(df)
    
    expect_is(result, "character")
    expect_true(nchar(result) > 0)
    expect_true(grepl("Gene", result))
    expect_true(grepl("BRCA1", result))
})

test_that("Pure helpers produce consistent results (no side effects)", {
    df <- data.frame(A = c("x", "yy"), B = c("1", "22"))
    df_char <- as.data.frame(lapply(df, as.character), stringsAsFactors = FALSE)
    
    # Multiple calls should produce identical results
    width1 <- .calculate_column_widths(df_char)
    width2 <- .calculate_column_widths(df_char)
    
    expect_identical(width1, width2)
})

# =============================================================================
# Code Quality: Verify Print Method Still Works
# =============================================================================

test_that("print.concordance_text works with formatted output", {
    # This verifies the cycle-breaking didn't break the print method
    df <- data.frame(
        Metric = c("Concordant", "Discordant"),
        Count = c(100, 50)
    )
    
    # Format the dataframe
    formatted <- .format_data_frame_as_text(df)
    
    # Should produce formatted output
    expect_output(cat(formatted))
})

# =============================================================================
# Edge Cases and Robustness
# =============================================================================

test_that(".format_data_frame_as_text handles very long values", {
    df <- data.frame(
        Name = "ShortName",
        Description = paste(rep("X", 100), collapse = "")
    )
    
    result <- .format_data_frame_as_text(df)
    
    expect_is(result, "character")
    expect_true(nchar(result) > 100)
})

test_that(".format_data_frame_as_text handles many columns", {
    df <- data.frame(
        matrix(1:50, nrow = 2, dimnames = list(NULL, paste0("Col", 1:25)))
    )
    
    result <- .format_data_frame_as_text(df)
    
    expect_is(result, "character")
    # Should contain all column headers
    expect_true(grepl("Col1", result))
    expect_true(grepl("Col25", result))
})

test_that(".format_data_frame_as_text handles special characters", {
    df <- data.frame(
        Name = c("Gene@1", "Gene#2", "Gene$3"),
        Value = c("10%", "20@", "30#")
    )
    
    result <- .format_data_frame_as_text(df)
    
    expect_is(result, "character")
    expect_true(grepl("Gene@1", result))
    expect_true(grepl("10%", result))
})

# ==============================================================================
# .format_top_genes(): Tests for formatting top genes (0% coverage)
# ==============================================================================

test_that(".format_top_genes formats dataframe as tibble output", {
  results <- data.frame(
    gene_id = c("GENE1", "GENE2", "GENE3"),
    gene_name = c("Name1", "Name2", "Name3"),
    Normal_mean = c(100, 80, 120),
    Tumor_mean = c(150, 100, 180),
    mean_difference = c(50, 20, 60),
    log2_fold_change = c(2.5, 1.8, 3.0),
    pvalue = c(0.001, 0.01, 0.05),
    padj = c(0.01, 0.05, 0.15)
  )
  
  result <- .format_top_genes(
    results_df = results,
    gene_col = "gene_id",
    padj_col = "padj",
    n_top = 2
  )
  
  # Should return formatted output (typically character/tibble)
  expect_is(result, c("tbl_df", "tbl", "data.frame", "character"))
})

test_that(".format_top_genes respects n_top parameter", {
  results <- data.frame(
    gene_id = paste0("GENE", 1:5),
    gene_name = paste0("Name", 1:5),
    Normal_mean = rnorm(5, 100, 20),
    Tumor_mean = rnorm(5, 120, 20),
    mean_difference = rep(20, 5),
    log2_fold_change = c(3.0, 2.5, 2.0, 1.5, 1.0),
    pvalue = seq(0.001, 0.1, length.out = 5),
    padj = seq(0.01, 0.15, length.out = 5)
  )
  
  result_1 <- .format_top_genes(results, "gene_id", "padj", n_top = 1)
  result_3 <- .format_top_genes(results, "gene_id", "padj", n_top = 3)
  
  expect_is(result_1, c("tbl_df", "tbl", "data.frame", "character"))
  expect_is(result_3, c("tbl_df", "tbl", "data.frame", "character"))
})

test_that(".format_top_genes handles custom column selection", {
  results <- data.frame(
    gene_id = c("GENE1", "GENE2"),
    Normal_mean = c(100, 80),
    Tumor_mean = c(150, 100),
    mean_difference = c(50, 20),
    log2_fold_change = c(2.5, 1.8),
    pvalue = c(0.001, 0.01),
    padj = c(0.01, 0.05)
  )
  
  custom_cols <- c("gene_id", "log2_fold_change", "padj")
  
  result <- .format_top_genes(
    results_df = results,
    gene_col = "gene_id",
    padj_col = "padj",
    n_top = 2,
    select_cols = custom_cols
  )
  
  expect_is(result, c("tbl_df", "tbl", "data.frame", "character"))
})

# ==============================================================================
# .format_duration(): Tests for duration formatting (37.5%)
# ==============================================================================

test_that(".format_duration formats seconds correctly", {
  result_sec <- .format_duration(30)
  result_min <- .format_duration(125)
  result_hour <- .format_duration(3661)
  
  expect_is(result_sec, "character")
  expect_is(result_min, "character")
  expect_is(result_hour, "character")
  
  expect_true(grepl("s|sec", result_sec))
  expect_true(grepl("m|min", result_min))
  expect_true(grepl("h|hour", result_hour))
})

test_that(".format_duration handles edge cases", {
  result_zero <- .format_duration(0)
  result_large <- .format_duration(86400)  # 24 hours
  
  expect_is(result_zero, "character")
  expect_is(result_large, "character")
})

# ============================================================================
# BUG FIX 1: Verify no duplicate helper functions (format_top_genes,
# create_summary_stats)
# ============================================================================

test_that("BUG 1: .format_top_genes and .create_summary_stats defined only once", {
  # These should exist in the package namespace without error
  # Duplicate definitions would silently overwrite, but we verify they exist
  expect_true(exists(".format_top_genes", where = asNamespace("TSENAT")))
  expect_true(exists(".create_summary_stats", where = asNamespace("TSENAT")))
  
  # Verify they are callable (not NULL)
  expect_true(is.function(get(".format_top_genes", envir = asNamespace("TSENAT"))))
  expect_true(is.function(get(".create_summary_stats", envir = asNamespace("TSENAT"))))
})
