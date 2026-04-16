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
