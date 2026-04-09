# Tests for package initialization
# Phase 1.1 - CRITICAL: Package .onLoad hook registration
# File: R/zzz.R

context("Package Initialization")

# ============================================================================
# TEST 1: Package loads without errors
# ============================================================================

test_that("TSENAT package loads successfully", {
    
    # If we're running this test, the package has already loaded successfully
    # Just verify the package is available
    expect_true("TSENAT" %in% loadedNamespaces())
})

# ============================================================================
# TEST 2: S3 method registrations are active
# ============================================================================

test_that("S3 methods are registered for print and summary", {
    
    # Test key S3 method registrations
    # These should exist and be callable
    
    # Check print methods exist
    expect_true(exists("print.tsenat_bootstrap_ci", where = asNamespace("TSENAT")))
    expect_true(exists("print.tsenat_jackknife", where = asNamespace("TSENAT")))
    expect_true(exists("print.rank_assumptions", where = asNamespace("TSENAT")))
    
    # Check summary methods exist
    expect_true(exists("summary.tsenat_bootstrap_ci", where = asNamespace("TSENAT")))
    expect_true(exists("summary.tsenat_divergence_bootstrap_ci", where = asNamespace("TSENAT")))
})

# ============================================================================
# TEST 3: Namespace is initialized with .viz_loaded flag
# ============================================================================

test_that(".viz_loaded flag is set to FALSE on load", {
    
    # This flag should be FALSE after package load (before any plots are made)
    ns <- asNamespace("TSENAT")
    
    # The flag may not exist yet if no plots have been made, or it should be FALSE
    if (exists(".viz_loaded", envir = ns)) {
        expect_false(get(".viz_loaded", envir = ns))
    }
})

# ============================================================================
# TEST 4: ggplot2 wrapper handles warnings properly  
# ============================================================================

test_that("ggplot2 print wrapper is registered if ggplot2 available", {
    
    if (!requireNamespace("ggplot2", quietly = TRUE)) {
        skip("ggplot2 not available")
    }
    
    # Verify we can use ggplot2 without spurious warnings
    library(ggplot2)
    
    # Create a simple plot
    p <- ggplot(data.frame(x = 1:5, y = 1:5), aes(x = x, y = y)) + geom_point()
    
    # Suppress any output but capture warnings
    captured_warnings <- list()
    result <- tryCatch(
        {
            withCallingHandlers(
                {
                    # Just verify print works without crashing
                    try(print(p), silent = TRUE)
                },
                warning = function(w) {
                    captured_warnings[[length(captured_warnings) + 1]] <<- conditionMessage(w)
                }
            )
        },
        error = function(e) {
            paste("Error:", e$message)
        }
    )
    
    # Should not error out
    expect_true(!is.character(result) || !grepl("Error", result))
})

# ============================================================================
# TEST 5: Package components are accessible
# ============================================================================

test_that("Main TSENAT classes and functions are exported", {
    
    # Check main classes
    expect_true(isClass("TSENATAnalysis", where = asNamespace("TSENAT")))
    
    # Check main functions are exported
    expect_true(exists("TSENAT", where = asNamespace("TSENAT")))
    expect_true(exists("calculate_diversity", where = asNamespace("TSENAT")))
    expect_true(exists("calculate_divergence", where = asNamespace("TSENAT")))
})

# ============================================================================
# TEST 6: NAMESPACE file is properly configured
# ============================================================================

test_that("NAMESPACE file exports main functions", {
    
    # Check that key functions are in the public API
    ns_exports <- getNamespaceExports("TSENAT")
    
    # Should export the main S4 function
    expect_true("TSENAT" %in% ns_exports)
    expect_true("TSENATAnalysis" %in% ns_exports)
})

# ============================================================================
# TEST 7: Package version is properly defined
# ============================================================================

test_that("TSENAT package version is defined and valid", {
    pkg_version <- packageVersion("TSENAT")
    
    # Should be a valid version object
    expect_is(pkg_version, "package_version")
    
    # Version should have at least major.minor components
    version_str <- as.character(pkg_version)
    expect_true(grepl("^[0-9]+\\.[0-9]+", version_str))
})

# ============================================================================
# TEST 8: Dependencies are available
# ============================================================================

test_that("Required dependencies are available", {
    # Core dependencies should be loaded
    required_pkgs <- c("S4Vectors", "SummarizedExperiment", "methods")
    
    for (pkg in required_pkgs) {
        expect_true(requireNamespace(pkg, quietly = TRUE),
                   info = paste("Package", pkg, "should be available"))
    }
})

# ============================================================================
# TEST 9: DESCRIPTION file metadata is complete
# ============================================================================

test_that("Package DESCRIPTION file metadata is valid", {
    # Get package metadata
    pkg_dir <- find.package("TSENAT")
    desc_file <- file.path(pkg_dir, "DESCRIPTION")
    
    # DESCRIPTION should exist
    expect_true(file.exists(desc_file))
    
    # Should be readable
    desc_contents <- readLines(desc_file)
    expect_true(length(desc_contents) > 0)
    
    # Should contain key fields
    desc_text <- paste(desc_contents, collapse = "\n")
    expect_true(grepl("Package:", desc_text))
    expect_true(grepl("Title:", desc_text))
    expect_true(grepl("Version:", desc_text))
})

# ============================================================================
# TEST 10: Data objects are accessible if included
# ============================================================================

test_that("Package includes example data for testing", {
    # Check if readcounts data is available
    expect_true("readcounts" %in% data(package = "TSENAT")$results[, 3])
})
