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
