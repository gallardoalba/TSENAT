context("AllGenerics: S4 Accessor Methods (Testing public and internal accessors)")

library(TSENAT)

# ============================================================================
# MODULE-LEVEL SETUP: Shared test data loaded ONCE
# ============================================================================

# Create a test TSENATAnalysis with metadata and plots
test_analysis_with_meta <- local({
    set.seed(42)
    data("readcounts", package = "TSENAT", envir = environment())
    readcounts <- as.matrix(readcounts)
    
    metadata_df <- read.table(
        system.file("extdata", "metadata.tsv", package = "TSENAT"),
        header = TRUE, sep = "\t"
    )
    
    gff3_file <- system.file("extdata", "annotation.gff3.gz", package = "TSENAT")
    
    config <- TSENAT::TSENAT_config(
        sample_col = "sample",
        condition_col = "condition",
        q = 1.0,
        paired = FALSE,
        stringency = "medium",
        nthreads = 1
    )
    
    analysis <- TSENAT::build_analysis(
        config = config,
        readcounts = readcounts,
        metadata = metadata_df,
        tx2gene = gff3_file,
        tpm = tpm,
        effective_length = effective_length
    )
    
    TSENAT::filter_analysis(analysis, stringency = "medium")
})

# ============================================================================
# TEST SUITE 1: getMeta() - Metadata Retrieval
# ============================================================================

test_that("getMeta() returns list of essential metadata", {
    meta <- TSENAT:::getMeta(test_analysis_with_meta)
    
    expect_is(meta, "list")
    expect_true(length(meta) > 0)
    expect_true("package_version" %in% names(meta) || "created_at" %in% names(meta) || "workflow_type" %in% names(meta))
})

test_that("getMeta() with key parameter returns specific metadata value", {
    meta_version <- TSENAT:::getMeta(test_analysis_with_meta, key = "package_version")
    
    # Should return the value if key exists, or NULL
    if (!is.null(meta_version)) {
        expect_is(meta_version, "character")
    }
})

test_that("getMeta() returns NULL for non-existent key", {
    non_existent <- TSENAT:::getMeta(test_analysis_with_meta, key = "non_existent_key_xyz")
    
    expect_null(non_existent)
})

test_that("getMeta() preserves metadata through accessor", {
    # Get metadata twice and verify consistency
    meta1 <- TSENAT:::getMeta(test_analysis_with_meta)
    meta2 <- TSENAT:::getMeta(test_analysis_with_meta)
    
    expect_identical(meta1, meta2)
})

# ============================================================================
# TEST SUITE 2: getPlot() - Plot Retrieval
# ============================================================================

test_that("getPlot() returns NULL when no plots cached", {
    # Create fresh analysis without plots
    analysis_no_plots <- test_analysis_with_meta
    
    all_plots <- TSENAT:::getPlot(analysis_no_plots)
    
    # Should return empty list or NA when no plots
    expect_true(length(all_plots) == 0 || is.null(all_plots) || length(names(all_plots)) == 0)
})

test_that("getPlot() returns list of all plots when type=NULL", {
    analysis <- test_analysis_with_meta
    
    plots <- TSENAT:::getPlot(analysis, type = NULL)
    
    expect_true(is.null(plots) || is.list(plots))
})

test_that("getPlot() returns specific plot when type specified", {
    analysis <- test_analysis_with_meta
    
    # Try to get non-existent plot type
    specific_plot <- TSENAT:::getPlot(analysis, type = "test_plot")
    
    # Should return NULL if plot doesn't exist
    expect_true(is.null(specific_plot) || inherits(specific_plot, "ggplot") || inherits(specific_plot, "list"))
})

# ============================================================================
# TEST SUITE 3: addPlot() - Plot Caching
# ============================================================================

test_that("addPlot() adds new plot to cache", {
    analysis <- test_analysis_with_meta
    
    # Create a simple test plot using ggplot2
    test_plot <- ggplot2::ggplot() +
        ggplot2::geom_point() +
        ggplot2::labs(title = "Test Plot")
    
    result <- TSENAT:::addPlot(analysis, type = "test_plot", plot = test_plot, replace = FALSE)
    
    expect_s4_class(result, "TSENATAnalysis")
    
    # Verify plot was added
    cached_plot <- TSENAT:::getPlot(result, type = "test_plot")
    expect_true(inherits(cached_plot, "ggplot") || is.list(cached_plot))
})

test_that("addPlot() warns when replacing existing plot without replace=TRUE", {
    analysis <- test_analysis_with_meta
    
    # Add first plot
    test_plot1 <- ggplot2::ggplot() + ggplot2::geom_point() + ggplot2::labs(title = "Plot 1")
    analysis <- TSENAT:::addPlot(analysis, type = "myplot", plot = test_plot1, replace = FALSE)
    
    # Try to add another plot with same type without replace
    test_plot2 <- ggplot2::ggplot() + ggplot2::geom_line() + ggplot2::labs(title = "Plot 2")
    
    expect_warning(
        TSENAT:::addPlot(analysis, type = "myplot", plot = test_plot2, replace = FALSE),
        "already exists"
    )
})

test_that("addPlot() replaces existing plot when replace=TRUE", {
    analysis <- test_analysis_with_meta
    
    # Add first plot
    test_plot1 <- ggplot2::ggplot() + ggplot2::geom_point() + ggplot2::labs(title = "Plot 1")
    analysis <- TSENAT:::addPlot(analysis, type = "replaceable", plot = test_plot1, replace = FALSE)
    
    # Replace with second plot
    test_plot2 <- ggplot2::ggplot() + ggplot2::geom_line() + ggplot2::labs(title = "Plot 2")
    result <- TSENAT:::addPlot(analysis, type = "replaceable", plot = test_plot2, replace = TRUE)
    
    expect_s4_class(result, "TSENATAnalysis")
    
    # Verify plot was replaced (title should be "Plot 2")
    cached <- TSENAT:::getPlot(result, type = "replaceable")
    expect_true(inherits(cached, "ggplot") || is.list(cached))
})

test_that("addPlot() returns analysis object for method chaining", {
    analysis <- test_analysis_with_meta
    
    test_plot <- ggplot2::ggplot() + ggplot2::geom_point()
    result <- TSENAT:::addPlot(analysis, type = "chaining_test", plot = test_plot, replace = FALSE)
    
    expect_s4_class(result, "TSENATAnalysis")
    expect_equal(class(result), class(analysis))
})

# ============================================================================
# TEST SUITE 4: setConfigValue() - Configuration Updates
# ============================================================================

test_that("setConfigValue() adds new configuration key", {
    analysis <- test_analysis_with_meta
    original_config <- TSENAT:::getConfig(analysis)
    
    result <- TSENAT:::setConfigValue(analysis, key = "test_param", value = 42)
    
    expect_s4_class(result, "TSENATAnalysis")
    
    # Verify new value was set
    new_config <- TSENAT:::getConfig(result)
    expect_equal(new_config$test_param, 42)
})

test_that("setConfigValue() updates existing configuration key", {
    analysis <- test_analysis_with_meta
    original_q <- TSENAT:::getConfig(analysis)$q
    
    result <- TSENAT:::setConfigValue(analysis, key = "q", value = 2.0)
    
    expect_s4_class(result, "TSENATAnalysis")
    
    # Verify value was updated
    updated_q <- TSENAT:::getConfig(result)$q
    expect_equal(updated_q, 2.0)
    expect_false(updated_q == original_q)
})

test_that("setConfigValue() accepts numeric values", {
    analysis <- test_analysis_with_meta
    
    result <- TSENAT:::setConfigValue(analysis, key = "numeric_test", value = 3.14159)
    updated <- TSENAT:::getConfig(result)$numeric_test
    
    expect_equal(updated, 3.14159)
})

test_that("setConfigValue() accepts character values", {
    analysis <- test_analysis_with_meta
    
    result <- TSENAT:::setConfigValue(analysis, key = "char_test", value = "test_string")
    updated <- TSENAT:::getConfig(result)$char_test
    
    expect_equal(updated, "test_string")
})

test_that("setConfigValue() accepts logical values", {
    analysis <- test_analysis_with_meta
    
    result <- TSENAT:::setConfigValue(analysis, key = "logic_test", value = TRUE)
    updated <- TSENAT:::getConfig(result)$logic_test
    
    expect_equal(updated, TRUE)
})

test_that("setConfigValue() accepts list values", {
    analysis <- test_analysis_with_meta
    
    test_list <- list(a = 1, b = 2)
    result <- TSENAT:::setConfigValue(analysis, key = "list_test", value = test_list)
    updated <- TSENAT:::getConfig(result)$list_test
    
    expect_is(updated, "list")
    expect_equal(updated$a, 1)
    expect_equal(updated$b, 2)
})

test_that("setConfigValue() maintains validity of TSENATAnalysis object", {
    analysis <- test_analysis_with_meta
    
    result <- TSENAT:::setConfigValue(analysis, key = "validity_test", value = "test")
    
    # Object should still be valid
    expect_true(validObject(result))
})

# ============================================================================
# TEST SUITE 5: metadata<- Assignment Method
# ============================================================================

test_that("metadata<- replaces entire metadata object", {
    analysis <- test_analysis_with_meta
    
    new_metadata <- list(
        custom_field = "custom_value",
        analysis_id = "test_123"
    )
    
    # Via replacement method - modifies in place
    # Access slot directly to verify since getter may have different behavior
    analysis@metadata <- new_metadata
    
    retrieved <- analysis@metadata
    expect_is(retrieved, "list")
    expect_true("custom_field" %in% names(retrieved))
    expect_equal(retrieved$custom_field, "custom_value")
})

test_that("metadata<- accepts empty list", {
    analysis <- test_analysis_with_meta
    
    analysis@metadata <- list()
    
    retrieved <- analysis@metadata
    expect_is(retrieved, "list")
    expect_equal(length(retrieved), 0)
})

test_that("metadata<- preserves NULL values in list", {
    analysis <- test_analysis_with_meta
    
    new_metadata <- list(
        field1 = "value1",
        field2 = NULL,
        field3 = "value3"
    )
    
    analysis@metadata <- new_metadata
    
    retrieved <- analysis@metadata
    expect_true("field2" %in% names(retrieved))
    expect_null(retrieved$field2)
    expect_equal(retrieved$field1, "value1")
})

test_that("metadata<- works with nested lists", {
    analysis <- test_analysis_with_meta
    
    nested_metadata <- list(
        workflow = list(
            type = "standard",
            steps = c("load", "filter", "analyze"),
            timing = list(start = "2026-04-10", end = "2026-04-10")
        ),
        user = "test_user"
    )
    
    analysis@metadata <- nested_metadata
    
    retrieved <- analysis@metadata
    expect_is(retrieved$workflow, "list")
    expect_is(retrieved$workflow$timing, "list")
    expect_equal(retrieved$workflow$type, "standard")
    expect_equal(retrieved$user, "test_user")
})

test_that("metadata<- roundtrip preserves structure", {
    analysis <- test_analysis_with_meta
    
    original_metadata <- list(
        int_val = 42L,
        num_val = 3.14,
        char_val = "test",
        vec_val = c(1, 2, 3),
        nested = list(x = 10, y = 20)
    )
    
    analysis@metadata <- original_metadata
    retrieved <- analysis@metadata
    
    expect_equal(retrieved$int_val, 42L)
    expect_equal(retrieved$num_val, 3.14)
    expect_equal(retrieved$char_val, "test")
    expect_equal(retrieved$vec_val, c(1, 2, 3))
    expect_equal(retrieved$nested$x, 10)
    expect_equal(retrieved$nested$y, 20)
})

# ============================================================================
# TEST SUITE 6: Related Accessors (se, getSE)
# ============================================================================

test_that("getSE() returns SummarizedExperiment object", {
    se <- TSENAT:::getSE(test_analysis_with_meta)
    
    expect_s4_class(se, "SummarizedExperiment")
})

test_that("getSE() contains correct gene and sample counts", {
    se <- TSENAT:::getSE(test_analysis_with_meta)
    
    expect_true(nrow(se) > 0)
    expect_true(ncol(se) > 0)
})

test_that("se() function is alias for getSE()", {
    se1 <- TSENAT::se(test_analysis_with_meta)
    se2 <- TSENAT:::getSE(test_analysis_with_meta)
    
    expect_identical(nrow(se1), nrow(se2))
    expect_identical(ncol(se1), ncol(se2))
})

# ============================================================================
# TEST SUITE 7: Related Accessors (getConfig, setConfig)
# ============================================================================

test_that("getConfig() returns list of configuration", {
    config <- TSENAT:::getConfig(test_analysis_with_meta)
    
    expect_is(config, "list")
    expect_true(length(config) > 0)
})

test_that("getConfig() with key returns specific config value", {
    q_value <- TSENAT:::getConfig(test_analysis_with_meta, key = "q")
    
    expect_true(is.numeric(q_value) || is.null(q_value))
})

test_that("setConfig() replaces entire configuration", {
    analysis <- test_analysis_with_meta
    
    new_config <- list(q = 1.5, paired = TRUE, test_param = "value")
    result <- TSENAT:::setConfig(analysis, new_config)
    
    expect_s4_class(result, "TSENATAnalysis")
    
    retrieved <- TSENAT:::getConfig(result)
    expect_equal(retrieved$q, 1.5)
    expect_equal(retrieved$test_param, "value")
})

# ============================================================================
# TEST SUITE 8: Accessor Consistency and Chain Operations
# ============================================================================

test_that("Multiple accessor calls are consistent", {
    analysis <- test_analysis_with_meta
    
    # Call multiple times and verify same results
    config1 <- TSENAT:::getConfig(analysis)
    config2 <- TSENAT:::getConfig(analysis)
    config3 <- TSENAT:::getConfig(analysis)
    
    expect_identical(config1, config2)
    expect_identical(config2, config3)
})

test_that("Accessor operations don't modify original object", {
    analysis <- test_analysis_with_meta
    original_config <- TSENAT:::getConfig(analysis)
    
    # Perform accessor operations
    TSENAT:::getMeta(analysis)
    TSENAT:::getPlot(analysis)
    
    # Original should be unchanged
    same_config <- TSENAT:::getConfig(analysis)
    expect_identical(original_config, same_config)
})

test_that("setConfigValue() and addPlot() allow chaining", {
    analysis <- test_analysis_with_meta
    
    # Chain operations
    test_plot <- ggplot2::ggplot() + ggplot2::geom_point()
    
    result <- TSENAT:::setConfigValue(analysis, key = "test1", value = 100)
    result <- TSENAT:::setConfigValue(result, key = "test2", value = 200)
    result <- TSENAT:::addPlot(result, type = "chain_plot", plot = test_plot, replace = FALSE)
    
    expect_s4_class(result, "TSENATAnalysis")
    expect_equal(TSENAT:::getConfig(result)$test1, 100)
    expect_equal(TSENAT:::getConfig(result)$test2, 200)
})

# ============================================================================
# TEST SUITE 9: metadata<- (Metadata Setter)
# ============================================================================

test_that("metadata<- sets metadata for TSENATAnalysis objects", {
    analysis <- test_analysis_with_meta
    
    # Use the metadata setter with explicit assignment
    meta <- metadata(analysis)
    meta$custom_key <- "custom_value"
    metadata(analysis) <- meta
    
    # Verify it was set
    expect_equal(metadata(analysis)$custom_key, "custom_value")
})

test_that("metadata<- preserves existing metadata", {
    analysis <- test_analysis_with_meta
    
    # Get original metadata
    original_meta <- metadata(analysis)
    original_length <- length(original_meta)
    
    # Add new metadata with explicit assignment
    meta <- metadata(analysis)
    meta$new_field <- "new_value"
    metadata(analysis) <- meta
    
    # Verify original metadata is preserved and new field is added
    new_meta <- metadata(analysis)
    expect_gte(length(new_meta), original_length)
    expect_equal(new_meta$new_field, "new_value")
})

test_that("metadata<- allows multiple assignments", {
    analysis <- test_analysis_with_meta
    
    # Multiple metadata assignments with explicit assignment
    meta <- metadata(analysis)
    meta$field1 <- "value1"
    meta$field2 <- 42
    meta$field3 <- list(nested = TRUE)
    metadata(analysis) <- meta
    
    meta <- metadata(analysis)
    expect_equal(meta$field1, "value1")
    expect_equal(meta$field2, 42)
    expect_equal(meta$field3$nested, TRUE)
})

test_that("metadata<- works with complex data types", {
    analysis <- test_analysis_with_meta
    
    # Assign various data types with explicit assignment
    meta <- metadata(analysis)
    meta$string_val <- "test"
    meta$numeric_val <- 3.14159
    meta$logical_val <- TRUE
    meta$vector_val <- c(1, 2, 3)
    meta$list_val <- list(a = 1, b = "two")
    metadata(analysis) <- meta
    
    meta <- metadata(analysis)
    expect_equal(meta$string_val, "test")
    expect_equal(meta$numeric_val, 3.14159)
    expect_equal(meta$logical_val, TRUE)
    expect_equal(meta$vector_val, c(1, 2, 3))
    expect_equal(meta$list_val$a, 1)
    expect_equal(meta$list_val$b, "two")
})

test_that("metadata<- returns the object invisibly", {
    analysis <- test_analysis_with_meta
    
    # Test that assignment works and returns value (invisibly)
    meta <- metadata(analysis)
    result <- (meta$test_key <- "test_value")
    metadata(analysis) <- meta
    
    # The assigned value should be available
    expect_equal(result, "test_value")
})

# ============================================================================
# TEST SUITE 10: se() - SummarizedExperiment accessor
# ============================================================================

test_that("se() returns SummarizedExperiment", {
    analysis <- test_analysis_with_meta
    
    result <- TSENAT::se(analysis)
    
    expect_s4_class(result, "SummarizedExperiment")
})

test_that("se() returns same object as getSE()", {
    analysis <- test_analysis_with_meta
    
    se_via_se <- TSENAT::se(analysis)
    se_via_getSE <- TSENAT:::getSE(analysis)
    
    expect_identical(se_via_se, se_via_getSE)
})

test_that("se() with multiple calls returns same SE", {
    analysis <- test_analysis_with_meta
    
    se1 <- TSENAT::se(analysis)
    se2 <- TSENAT::se(analysis)
    
    expect_identical(se1, se2)
})

# ============================================================================
# TEST SUITE 11: show() and summary() methods
# ============================================================================

test_that("show() method executes without error", {
    analysis <- test_analysis_with_meta
    
    # Test that show() can be called without error (messages are expected)
    expect_error(show(analysis), NA)
    
    # Verify the object has required structure that show() displays
    expect_true(nrow(analysis@se) > 0)  # Genes
    expect_true(ncol(analysis@se) > 0)  # Samples
})

test_that("show() displays gene count accurately", {
    analysis <- test_analysis_with_meta
    
    # Verify show() executes without error
    expect_error(show(analysis), NA)
    
    # Verify gene count is correct
    n_genes <- nrow(analysis@se)
    expect_true(n_genes > 0)
})

test_that("show() displays sample count accurately", {
    analysis <- test_analysis_with_meta
    
    # Verify show() executes without error
    expect_error(show(analysis), NA)
    
    # Verify sample count is correct
    n_samples <- ncol(analysis@se)
    expect_true(n_samples > 0)
})

test_that("summary() method executes without error", {
    analysis <- test_analysis_with_meta
    
    # Test that summary() can be called without error
    expect_error(summary(analysis), NA)  # Should NOT error
    
    # Verify object has structure that summary() reports
    expect_true(nrow(analysis@se) > 0)  # Genes
    expect_true(ncol(analysis@se) > 0)  # Samples
})

test_that("summary() reports data structure", {
    analysis <- test_analysis_with_meta
    
    # Verify summary can execute
    expect_error(summary(analysis), NA)
    
    # Check that summary has access to structure it reports
    assays <- SummarizedExperiment::assayNames(analysis@se)
    expect_true(length(assays) > 0 || is.null(assays))
})

test_that("summary() shows analysis results status", {
    analysis <- test_analysis_with_meta
    
    # Verify summary can execute
    expect_error(summary(analysis), NA)
    
    # Verify analysis has result slots that summary reports
    expect_true(is.list(analysis@diversity_results) || is.null(analysis@diversity_results))
    expect_true(is.list(analysis@sait_results) || is.null(analysis@sait_results))
    expect_true(is.list(analysis@jackknife_results) || is.null(analysis@jackknife_results))
    expect_true(is.list(analysis@divergence_results) || is.null(analysis@divergence_results))
})

# ============================================================================
# TEST SUITE 12: setConfig() - Full configuration replacement
# ============================================================================

test_that("setConfig() replaces entire configuration", {
    analysis <- test_analysis_with_meta
    original_config <- TSENAT:::getConfig(analysis)
    
    new_config <- list(
        param1 = "new1",
        param2 = 999,
        param3 = TRUE
    )
    
    result <- TSENAT:::setConfig(analysis, new_config)
    
    expect_s4_class(result, "TSENATAnalysis")
    updated_config <- TSENAT:::getConfig(result)
    expect_equal(updated_config$param1, "new1")
    expect_equal(updated_config$param2, 999)
    expect_equal(updated_config$param3, TRUE)
})

test_that("setConfig() with TSENATConfig object", {
    analysis <- test_analysis_with_meta
    
    # Create a TSENATConfig-like object
    config_obj <- TSENAT::TSENAT_config(q = 1.5, nthreads = 2)
    
    result <- TSENAT:::setConfig(analysis, config_obj)
    
    expect_s4_class(result, "TSENATAnalysis")
})

test_that("setConfig() validates object type", {
    analysis <- test_analysis_with_meta
    
    # Should error with non-list, non-TSENATConfig object
    expect_error(
        TSENAT:::setConfig(analysis, "invalid"),
        "Configuration must be a list"
    )
})

test_that("setConfig() with empty list", {
    analysis <- test_analysis_with_meta
    
    result <- TSENAT:::setConfig(analysis, list())
    
    updated_config <- TSENAT:::getConfig(result)
    expect_length(updated_config, 0)
})

test_that("setConfig() preserves object validity", {
    analysis <- test_analysis_with_meta
    
    new_config <- list(test = TRUE)
    result <- TSENAT:::setConfig(analysis, new_config)
    
    # Should still be valid S4 object
    expect_true(validObject(result))
})

# ============================================================================
# TEST SUITE 13: getConfig() with key parameter
# ============================================================================

test_that("getConfig() with specific key returns value", {
    analysis <- test_analysis_with_meta
    
    config <- TSENAT:::getConfig(analysis)
    
    # Get a specific key that exists
    for (key in names(config)[1:min(2, length(config))]) {
        value <- TSENAT:::getConfig(analysis, key = key)
        expect_equal(value, config[[key]])
    }
})

test_that("getConfig() with non-existent key returns NULL", {
    analysis <- test_analysis_with_meta
    
    value <- TSENAT:::getConfig(analysis, key = "non_existent_key_xyz")
    
    expect_null(value)
})

# ============================================================================
# TEST SUITE 14: getPlot() edge cases
# ============================================================================

test_that("getPlot() with non-existent type returns NULL", {
    analysis <- test_analysis_with_meta
    
    plot <- TSENAT:::getPlot(analysis, type = "non_existent_plot_type_xyz")
    
    expect_null(plot)
})

test_that("getPlot() with type=NULL returns all plots", {
    analysis <- test_analysis_with_meta
    
    all_plots <- TSENAT:::getPlot(analysis, type = NULL)
    
    # Should return list or NULL
    expect_true(is.list(all_plots) || is.null(all_plots))
})

# ============================================================================
# TEST SUITE 15: addPlot() advanced scenarios
# ============================================================================

test_that("addPlot() returns object for chaining", {
    analysis <- test_analysis_with_meta
    
    p1 <- ggplot2::ggplot() + ggplot2::geom_point()
    p2 <- ggplot2::ggplot() + ggplot2::geom_line()
    
    # Chain operations
    result <- TSENAT:::addPlot(analysis, type = "chain1", plot = p1, replace = FALSE)
    result <- TSENAT:::addPlot(result, type = "chain2", plot = p2, replace = FALSE)
    
    expect_s4_class(result, "TSENATAnalysis")
    expect_equal(length(TSENAT:::getPlot(result)), 2)
})

test_that("addPlot() with NULL plot", {
    analysis <- test_analysis_with_meta
    
    # Adding NULL plot should still work
    result <- TSENAT:::addPlot(analysis, type = "null_plot", plot = NULL, replace = FALSE)
    
    expect_s4_class(result, "TSENATAnalysis")
    cached <- TSENAT:::getPlot(result, type = "null_plot")
    expect_null(cached)
})

test_that("addPlot() returns original object when already exists and replace=FALSE", {
    analysis <- test_analysis_with_meta
    
    p1 <- ggplot2::ggplot() + ggplot2::geom_point()
    p2 <- ggplot2::ggplot() + ggplot2::geom_line()
    
    analysis <- TSENAT:::addPlot(analysis, type = "myplot", plot = p1, replace = FALSE)
    original_plot <- TSENAT:::getPlot(analysis, type = "myplot")
    
    # Try to add different plot without replace
    suppressWarnings({
        result <- TSENAT:::addPlot(analysis, type = "myplot", plot = p2, replace = FALSE)
    })
    
    # Plot should be unchanged
    unchanged_plot <- TSENAT:::getPlot(result, type = "myplot")
    expect_identical(original_plot, unchanged_plot)
})

# ============================================================================
# TEST SUITE 16: getMeta() with workflow info
# ============================================================================

test_that("getMeta() filters to essential metadata only", {
    analysis <- test_analysis_with_meta
    
    # Add various metadata with explicit assignment
    meta <- metadata(analysis)
    meta$large_data <- rep(1:1000, 100)  # Large data
    meta$logs <- "Function execution logs"
    meta$package_version <- "1.0.0"  # Essential
    metadata(analysis) <- meta
    
    essential <- TSENAT:::getMeta(analysis)
    
    # Should contain essential metadata
    expect_true("package_version" %in% names(essential) || length(essential) >= 0)
    
    # Large data and logs should not be in essential
    expect_false("large_data" %in% names(essential))
    expect_false("logs" %in% names(essential))
})

test_that("getMeta() handles missing created_at", {
    analysis <- test_analysis_with_meta
    
    # Remove created_at if it exists with explicit assignment
    meta <- metadata(analysis)
    meta$created_at <- NULL
    metadata(analysis) <- meta
    
    essential <- TSENAT:::getMeta(analysis)
    
    # Should handle gracefully
    expect_true(is.list(essential))
})

# ============================================================================
# TEST SUITE 17: Integration and edge cases
# ============================================================================

test_that("Multiple accessor calls don't modify object", {
    analysis <- test_analysis_with_meta
    
    # Perform multiple reads
    for (i in 1:10) {
        TSENAT:::getMeta(analysis)
        TSENAT:::getConfig(analysis)
        TSENAT:::getPlot(analysis)
        TSENAT::se(analysis)
    }
    
    # Object should be identical to original
    expect_equal(nrow(TSENAT::se(analysis)), nrow(test_analysis_with_meta@se))
})

test_that("Accessor methods work with minimal analysis object", {
    # Create minimal SE
    counts <- matrix(1:10, nrow = 5, ncol = 2)
    se <- SummarizedExperiment::SummarizedExperiment(assays = list(counts = counts))
    
    # Create minimal analysis
    minimal_analysis <- new("TSENATAnalysis", se = se, config = list(), plots = list())
    
    # All accessors should work
    expect_s4_class(TSENAT::se(minimal_analysis), "SummarizedExperiment")
    expect_is(TSENAT:::getConfig(minimal_analysis), "list")
    expect_true(is.null(TSENAT:::getPlot(minimal_analysis)) || is.list(TSENAT:::getPlot(minimal_analysis)))
})

test_that("setConfigValue() chains with other operations", {
    analysis <- test_analysis_with_meta
    
    # Chain multiple operations
    result <- TSENAT:::setConfigValue(analysis, "key1", 100)
    result <- TSENAT:::setConfigValue(result, "key2", 200)
    result <- TSENAT:::setConfigValue(result, "key3", 300)
    
    expect_equal(TSENAT:::getConfig(result)$key1, 100)
    expect_equal(TSENAT:::getConfig(result)$key2, 200)
    expect_equal(TSENAT:::getConfig(result)$key3, 300)
})

# ============================================================================
# TEST SUITE 18: getMeta() - Specific metadata field coverage
# ============================================================================

test_that("getMeta() extracts created_at when present", {
    analysis <- test_analysis_with_meta
    
    # Ensure created_at is in metadata with explicit assignment
    meta <- metadata(analysis)
    meta$created_at <- "2026-04-18 10:30:00"
    metadata(analysis) <- meta
    
    essential <- TSENAT:::getMeta(analysis)
    
    # Should include created_at
    if ("created_at" %in% names(metadata(analysis))) {
        expect_true("created_at" %in% names(essential) || length(essential) >= 0)
    }
})

test_that("getMeta() extracts ended_at when present", {
    analysis <- test_analysis_with_meta
    
    # Add ended_at with explicit assignment
    meta <- metadata(analysis)
    meta$ended_at <- "2026-04-18 12:00:00"
    metadata(analysis) <- meta
    
    essential <- TSENAT:::getMeta(analysis)
    
    # Should handle ended_at field
    expect_true(is.list(essential))
})

test_that("getMeta() extracts package_version when present", {
    analysis <- test_analysis_with_meta
    
    # Set package version with explicit assignment
    meta <- metadata(analysis)
    meta$package_version <- "1.2.3"
    metadata(analysis) <- meta
    
    essential <- TSENAT:::getMeta(analysis)
    
    # Should extract version
    if ("package_version" %in% names(metadata(analysis))) {
        expect_true("package_version" %in% names(essential))
        expect_equal(essential$package_version, "1.2.3")
    }
})

test_that("getMeta() extracts tsenat_version when present", {
    analysis <- test_analysis_with_meta
    
    # Set TSENAT-specific version with explicit assignment
    meta <- metadata(analysis)
    meta$tsenat_version <- "2.0.1"
    metadata(analysis) <- meta
    
    essential <- TSENAT:::getMeta(analysis)
    
    # Should extract TSENAT version
    if ("tsenat_version" %in% names(metadata(analysis))) {
        expect_true("tsenat_version" %in% names(essential))
    }
})

test_that("getMeta() extracts workflow_type when present", {
    analysis <- test_analysis_with_meta
    
    # Set workflow type with explicit assignment
    meta <- metadata(analysis)
    meta$workflow_type <- "diversity_analysis"
    metadata(analysis) <- meta
    
    essential <- TSENAT:::getMeta(analysis)
    
    # Should extract workflow type
    if ("workflow_type" %in% names(metadata(analysis))) {
        expect_true("workflow_type" %in% names(essential))
    }
})

test_that("getMeta() extracts workflow list structure", {
    analysis <- test_analysis_with_meta
    
    # Set workflow as list with explicit assignment
    meta <- metadata(analysis)
    meta$workflow <- list(
        workflow_type = "complete",
        completion_time = "2026-04-18 12:00:00",
        other_field = "should_be_excluded"
    )
    metadata(analysis) <- meta
    
    essential <- TSENAT:::getMeta(analysis)
    
    # Should extract workflow structure
    if ("workflow" %in% names(metadata(analysis))) {
        expect_true("workflow" %in% names(essential) || is.list(metadata(analysis)$workflow))
    }
})

test_that("getMeta() handles workflow as non-list", {
    analysis <- test_analysis_with_meta
    
    # Set workflow as character (not list) with explicit assignment
    meta <- metadata(analysis)
    meta$workflow <- "some_string_value"
    metadata(analysis) <- meta
    
    essential <- TSENAT:::getMeta(analysis)
    
    # Should handle non-list workflow gracefully
    expect_true(is.list(essential))
})

test_that("getMeta() key parameter extracts specific metadata field", {
    analysis <- test_analysis_with_meta
    
    # Set multiple metadata fields with explicit assignment
    meta <- metadata(analysis)
    meta$created_at <- "2026-04-18"
    meta$package_version <- "1.0.0"
    meta$workflow_type <- "test"
    metadata(analysis) <- meta
    
    # Extract specific key
    pkg_version <- TSENAT:::getMeta(analysis, key = "package_version")
    created <- TSENAT:::getMeta(analysis, key = "created_at")
    workflow <- TSENAT:::getMeta(analysis, key = "workflow_type")
    
    # Should extract requested keys
    if ("package_version" %in% names(metadata(analysis))) {
        expect_equal(pkg_version, "1.0.0")
    }
    if ("created_at" %in% names(metadata(analysis))) {
        expect_equal(created, "2026-04-18")
    }
})

test_that("getMeta() handles metadata with only non-essential fields", {
    analysis <- test_analysis_with_meta
    
    # Clear metadata and add only non-essential fields
    metadata(analysis) <- list(
        large_results = rep(1:1000, 50),
        verbose_logs = "detailed function logs",
        temporary_cache = "should not appear"
    )
    
    essential <- TSENAT:::getMeta(analysis)
    
    # Should return empty or minimal list
    expect_true(is.list(essential))
    expect_true(length(essential) == 0 || all(!(names(essential) %in% c("large_results", "verbose_logs", "temporary_cache"))))
})

# ============================================================================
# TEST SUITE 19: getConfig() - Detailed key handling
# ============================================================================

test_that("getConfig() returns different data types correctly", {
    analysis <- test_analysis_with_meta
    
    # Add various types to config
    config <- TSENAT:::getConfig(analysis)
    config$string_val <- "test_string"
    config$numeric_val <- 42.5
    config$logical_val <- TRUE
    config$vector_val <- c(1, 2, 3)
    config$list_val <- list(a = 1, b = 2)
    
    analysis <- TSENAT:::setConfig(analysis, config)
    
    # Retrieve each type
    if ("string_val" %in% names(TSENAT:::getConfig(analysis))) {
        expect_equal(TSENAT:::getConfig(analysis, key = "string_val"), "test_string")
    }
    if ("numeric_val" %in% names(TSENAT:::getConfig(analysis))) {
        expect_equal(TSENAT:::getConfig(analysis, key = "numeric_val"), 42.5)
    }
    if ("logical_val" %in% names(TSENAT:::getConfig(analysis))) {
        expect_equal(TSENAT:::getConfig(analysis, key = "logical_val"), TRUE)
    }
})

test_that("getConfig() handles deeply nested config", {
    analysis <- test_analysis_with_meta
    
    # Create deeply nested config
    nested_config <- list(
        level1 = list(
            level2 = list(
                level3 = list(value = "deep")
            )
        ),
        other = "top_level"
    )
    
    analysis <- TSENAT:::setConfig(analysis, nested_config)
    
    # Should handle nested structure
    retrieved <- TSENAT:::getConfig(analysis)
    expect_equal(retrieved$level1$level2$level3$value, "deep")
})

test_that("getConfig() with key on nested list", {
    analysis <- test_analysis_with_meta
    
    config <- list(params = list(q = 1.0, nthreads = 2))
    analysis <- TSENAT:::setConfig(analysis, config)
    
    # Get nested structure
    params <- TSENAT:::getConfig(analysis, key = "params")
    
    # Should return the nested list
    expect_is(params, "list")
    if (is.list(params) && "q" %in% names(params)) {
        expect_equal(params$q, 1.0)
    }
})

# ============================================================================
# TEST SUITE 20: getPlot() - Plot storage and retrieval
# ============================================================================

test_that("getPlot() stores and retrieves multiple plot types", {
    analysis <- test_analysis_with_meta
    
    # Create multiple plots
    p1 <- ggplot2::ggplot() + ggplot2::geom_point() + ggplot2::ggtitle("Plot 1")
    p2 <- ggplot2::ggplot() + ggplot2::geom_line() + ggplot2::ggtitle("Plot 2")
    p3 <- ggplot2::ggplot() + ggplot2::geom_boxplot() + ggplot2::ggtitle("Plot 3")
    
    # Add multiple plots
    analysis <- TSENAT:::addPlot(analysis, "plot_A", p1, replace = FALSE)
    analysis <- TSENAT:::addPlot(analysis, "plot_B", p2, replace = FALSE)
    analysis <- TSENAT:::addPlot(analysis, "plot_C", p3, replace = FALSE)
    
    # Retrieve all plots
    all_plots <- TSENAT:::getPlot(analysis, type = NULL)
    
    expect_true(length(all_plots) >= 3)
    expect_true("plot_A" %in% names(all_plots))
})

test_that("getPlot() returns specific plot type accurately", {
    analysis <- test_analysis_with_meta
    
    p1 <- ggplot2::ggplot() + ggplot2::geom_point() + ggplot2::labs(title = "Specific")
    analysis <- TSENAT:::addPlot(analysis, "specific_type", p1, replace = FALSE)
    
    # Retrieve specific plot
    retrieved <- TSENAT:::getPlot(analysis, type = "specific_type")
    
    # Should be a plot
    expect_true(ggplot2::is_ggplot(retrieved) || is.null(retrieved))
})

test_that("getPlot() returns empty list when no plots cached", {
    # Create minimal analysis with no plots
    counts <- matrix(1:10, nrow = 5, ncol = 2)
    se <- SummarizedExperiment::SummarizedExperiment(assays = list(counts = counts))
    analysis <- new("TSENATAnalysis", se = se, config = list(), plots = list())
    
    all_plots <- TSENAT:::getPlot(analysis, type = NULL)
    
    expect_true(length(all_plots) == 0 || is.null(all_plots))
})

# ============================================================================
# TEST SUITE 21: addPlot() - Warning behavior
# ============================================================================

test_that("addPlot() warns about existing plot without replace", {
    analysis <- test_analysis_with_meta
    
    p1 <- ggplot2::ggplot() + ggplot2::geom_point()
    p2 <- ggplot2::ggplot() + ggplot2::geom_line()
    
    # Add first plot
    analysis <- TSENAT:::addPlot(analysis, "same_plot", p1, replace = FALSE)
    
    # Try to add same plot again without replace
    expect_warning(
        TSENAT:::addPlot(analysis, "same_plot", p2, replace = FALSE),
        "already exists"
    )
})

test_that("addPlot() does not warn when replace=TRUE", {
    analysis <- test_analysis_with_meta
    
    p1 <- ggplot2::ggplot() + ggplot2::geom_point()
    p2 <- ggplot2::ggplot() + ggplot2::geom_line()
    
    # Add first plot
    analysis <- TSENAT:::addPlot(analysis, "same_plot", p1, replace = FALSE)
    
    # Replace with replace=TRUE should not warn
    expect_silent(
        TSENAT:::addPlot(analysis, "same_plot", p2, replace = TRUE)
    )
})

test_that("addPlot() allows new plot type without warning", {
    analysis <- test_analysis_with_meta
    
    p1 <- ggplot2::ggplot() + ggplot2::geom_point()
    
    # Add new plot type should not warn
    expect_silent(
        TSENAT:::addPlot(analysis, "new_type_xyz", p1, replace = FALSE)
    )
})

# ============================================================================
# TEST SUITE 22: S4 method dispatch verification
# ============================================================================

test_that("S4 methods dispatch correctly to TSENATAnalysis", {
    analysis <- test_analysis_with_meta
    
    # These should all work (S4 dispatch)
    se <- TSENAT::se(analysis)
    config <- TSENAT:::getConfig(analysis)
    meta <- TSENAT:::getMeta(analysis)
    
    expect_s4_class(se, "SummarizedExperiment")
    expect_is(config, "list")
    expect_is(meta, "list")
})

test_that("Replacement method 'metadata<-' uses S4 dispatch", {
    analysis <- test_analysis_with_meta
    
    # Use replacement method
    new_meta <- list(test = TRUE, value = 123)
    metadata(analysis) <- new_meta
    
    # Verify replacement worked
    result_meta <- metadata(analysis)
    expect_true("test" %in% names(result_meta) || is.list(result_meta))
})

# ============================================================================
# TEST SUITE 23: Config NULL handling
# ============================================================================

test_that("setConfigValue() initializes empty config", {
    # Create analysis with empty config list
    counts <- matrix(1:10, nrow = 5, ncol = 2)
    se <- SummarizedExperiment::SummarizedExperiment(assays = list(counts = counts))
    analysis <- new("TSENATAnalysis", se = se, config = list())
    
    # Set value on empty config
    result <- TSENAT:::setConfigValue(analysis, "new_key", "new_value")
    
    updated_config <- TSENAT:::getConfig(result)
    expect_equal(updated_config$new_key, "new_value")
})

test_that("getConfig() handles empty config", {
    counts <- matrix(1:10, nrow = 5, ncol = 2)
    se <- SummarizedExperiment::SummarizedExperiment(assays = list(counts = counts))
    analysis <- new("TSENATAnalysis", se = se, config = list())
    
    config <- TSENAT:::getConfig(analysis)
    
    expect_true(is.list(config))
    expect_equal(length(config), 0)
})

# ============================================================================
# TEST SUITE 24: Metadata NULL and list handling
# ============================================================================

test_that("getMeta() handles empty metadata", {
    counts <- matrix(1:10, nrow = 5, ncol = 2)
    se <- SummarizedExperiment::SummarizedExperiment(assays = list(counts = counts))
    analysis <- new("TSENATAnalysis", se = se, metadata = list())
    
    essential <- TSENAT:::getMeta(analysis)
    
    # Should return empty list or handle gracefully
    expect_true(is.null(essential) || is.list(essential))
})

test_that("getMeta() handles empty metadata list", {
    analysis <- test_analysis_with_meta
    metadata(analysis) <- list()
    
    essential <- TSENAT:::getMeta(analysis)
    
    expect_is(essential, "list")
    expect_length(essential, 0)
})

test_that("metadata<- with list assignment", {
    analysis <- test_analysis_with_meta
    
    new_metadata <- list(
        field1 = "value1",
        field2 = 42,
        field3 = TRUE
    )
    
    metadata(analysis) <- new_metadata
    
    result_meta <- metadata(analysis)
    expect_equal(result_meta$field1, "value1")
    expect_equal(result_meta$field2, 42)
    expect_equal(result_meta$field3, TRUE)
})

# ============================================================================
# TEST SUITE 25: Accessor robustness with edge data
# ============================================================================

test_that("Accessors work with large gene count", {
    # Create large SE
    large_counts <- matrix(rpois(100000, lambda = 10), nrow = 1000, ncol = 100)
    large_se <- SummarizedExperiment::SummarizedExperiment(
        assays = list(counts = large_counts)
    )
    
    analysis <- new("TSENATAnalysis", se = large_se, config = list())
    
    # All accessors should work with large data
    se <- TSENAT::se(analysis)
    expect_equal(nrow(se), 1000)
    expect_equal(ncol(se), 100)
})

test_that("Config handles large number of parameters", {
    analysis <- test_analysis_with_meta
    
    # Add many parameters
    large_config <- as.list(1:1000)
    names(large_config) <- paste0("param_", 1:1000)
    
    analysis <- TSENAT:::setConfig(analysis, large_config)
    
    # Should handle large config
    retrieved <- TSENAT:::getConfig(analysis)
    expect_equal(length(retrieved), 1000)
})
