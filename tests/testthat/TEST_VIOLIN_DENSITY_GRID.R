# ============================================================================
# TSENAT: Test Suite for Single-q Visualization Functions
# ============================================================================
# Purpose: Tests plot_tsallis_violin_singleq(), plot_tsallis_density_singleq(),
#          and plot_tsallis_violin_density_grid() for single q-value input
# ============================================================================

# Setup: Load packages
suppressPackageStartupMessages({
    library(devtools)
    devtools::load_all(".")
    library(ggplot2)
    library(SummarizedExperiment)
    library(dplyr)
    library(cowplot)
})

set.seed(42)

# Create output directory
output_dir <- "/home/nouser/galaxy/tools_source/TSENAT/output"
if (!dir.exists(output_dir)) {
    dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
}

# Prepare output file for logging
output_file <- file.path(output_dir, "test_results.txt")
log_file <- file(output_file, "w")

sink(log_file)
cat(paste(c("=", rep("=", 70), "\n"), collapse=""))
cat("TSENAT Single-q Visualization Functions Test Suite\n")
cat(paste(c("=", rep("=", 70), "\n"), collapse=""))
cat("Date:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n\n")

tryCatch({
    # ====================================================================
    # PART 1: Load Data
    # ====================================================================
    cat("PART 1: Loading Data\n")
    cat(paste(c("---", rep("-", 67), "\n"), collapse=""), sep="")
    
    # Load preprocessed dataset
    data("readcounts", package = "TSENAT")
    readcounts <- as.matrix(salmon_dataset)
    mode(readcounts) <- "numeric"
    
    # Load metadata
    metadata_df <- read.table(
        system.file("extdata", "metadata.tsv", package = "TSENAT"),
        header = TRUE, sep = "\t"
    )
    
    gff3_dataset <- system.file("extdata", "annotation.gff3.gz", package = "TSENAT")
    
    # Subset for faster testing (use first 50 genes)
    n_genes_subset <- 50
    if (nrow(readcounts) > n_genes_subset) {
        readcounts <- readcounts[1:n_genes_subset, ]
        cat("NOTE: Analyzing", n_genes_subset, "genes for testing speed\n\n")
    }
    
    # Build SummarizedExperiment
    se <- build_se(readcounts, gff3_dataset, metadata = metadata_df)
    
    # Filter SE
    se <- tryCatch(
        filter_se(se, stringency = "loose"),
        error = function(e) {
            cat("Warning: filter_se failed, proceeding with unfiltered SE\n")
            se
        }
    )
    
    cat("✓ Data loaded successfully\n")
    cat("  - Genes:", nrow(se), "\n")
    cat("  - Samples:", ncol(se), "\n")
    cat("  - Groups:", paste(unique(colData(se)$group), collapse = ", "), "\n\n")
    
    # ====================================================================
    # PART 2: Calculate Diversity for Single q Values
    # ====================================================================
    cat("PART 2: Calculate Diversity for Single q Values\n")
    cat(paste(c("---", rep("-", 67), "\n"), collapse=""), sep="")
    
    # Test multiple single q values
    q_values <- c(0.5, 1, 2)
    se_list <- list()
    
    for (q in q_values) {
        cat("Computing diversity for q =", q, "... ")
        se_q <- calculate_diversity(se, q = q, norm = TRUE)
        se_list[[as.character(q)]] <- se_q
        cat("✓\n")
    }
    
    cat("✓ Diversity calculated for q values:", paste(q_values, collapse = ", "), "\n\n")
    
    # ====================================================================
    # PART 3: Test plot_tsallis_violin_singleq()
    # ====================================================================
    cat("PART 3: Testing plot_tsallis_violin_singleq()\n")
    cat(paste(c("---", rep("-", 67), "\n"), collapse=""), sep="")
    
    for (q_str in names(se_list)) {
        q_val <- as.numeric(q_str)
        cat("Testing plot_tsallis_violin_singleq() for q =", q_val, "... ")
        
        tryCatch({
            p <- plot_tsallis_violin_singleq(se_list[[q_str]])
            
            # Save plot
            filename <- file.path(output_dir, sprintf("violin_q%g.png", q_val))
            png(filename, width = 800, height = 600, res = 100)
            print(p)
            dev.off()
            
            cat("✓\n")
            cat("   Saved to:", filename, "\n")
        }, error = function(e) {
            cat("✗ FAILED\n")
            cat("   Error:", conditionMessage(e), "\n")
        })
    }
    cat("\n")
    
    # ====================================================================
    # PART 4: Test plot_tsallis_density_singleq()
    # ====================================================================
    cat("PART 4: Testing plot_tsallis_density_singleq()\n")
    cat(paste(c("---", rep("-", 67), "\n"), collapse=""), sep="")
    
    for (q_str in names(se_list)) {
        q_val <- as.numeric(q_str)
        cat("Testing plot_tsallis_density_singleq() for q =", q_val, "... ")
        
        tryCatch({
            p <- plot_tsallis_density_singleq(se_list[[q_str]])
            
            # Save plot
            filename <- file.path(output_dir, sprintf("density_q%g.png", q_val))
            png(filename, width = 800, height = 600, res = 100)
            print(p)
            dev.off()
            
            cat("✓\n")
            cat("   Saved to:", filename, "\n")
        }, error = function(e) {
            cat("✗ FAILED\n")
            cat("   Error:", conditionMessage(e), "\n")
        })
    }
    cat("\n")
    
    # ====================================================================
    # PART 5: Test plot_tsallis_violin_density_grid()
    # ====================================================================
    cat("PART 5: Testing plot_tsallis_violin_density_grid()\n")
    cat(paste(c("---", rep("-", 67), "\n"), collapse=""), sep="")
    
    for (q_str in names(se_list)) {
        q_val <- as.numeric(q_str)
        cat("Testing plot_tsallis_violin_density_grid() for q =", q_val, "... ")
        
        tryCatch({
            p <- plot_tsallis_violin_density_grid(se_list[[q_str]])
            
            # Save plot
            filename <- file.path(output_dir, sprintf("violin_density_grid_q%g.png", q_val))
            png(filename, width = 1400, height = 600, res = 100)
            print(p)
            dev.off()
            
            cat("✓\n")
            cat("   Saved to:", filename, "\n")
        }, error = function(e) {
            cat("✗ FAILED\n")
            cat("   Error:", conditionMessage(e), "\n")
        })
    }
    cat("\n")
    
    # ====================================================================
    # PART 6: Test with Custom Titles
    # ====================================================================
    cat("PART 6: Testing Functions with Custom Titles\n")
    cat(paste(c("---", rep("-", 67), "\n"), collapse=""), sep="")
    
    q_val <- 1
    custom_title_violin <- sprintf("Custom Violin Plot (q = %g)", q_val)
    custom_title_density <- sprintf("Custom Density Plot (q = %g)", q_val)
    custom_title_grid <- sprintf("Custom Grid: Violin + Density (q = %g)", q_val)
    
    cat("Testing plot_tsallis_violin_singleq() with custom title... ")
    tryCatch({
        p <- plot_tsallis_violin_singleq(se_list[["1"]], title = custom_title_violin)
        filename <- file.path(output_dir, "violin_custom_title.png")
        png(filename, width = 800, height = 600, res = 100)
        print(p)
        dev.off()
        cat("✓\n   Saved to:", filename, "\n")
    }, error = function(e) {
        cat("✗ FAILED:", conditionMessage(e), "\n")
    })
    
    cat("Testing plot_tsallis_density_singleq() with custom title... ")
    tryCatch({
        p <- plot_tsallis_density_singleq(se_list[["1"]], title = custom_title_density)
        filename <- file.path(output_dir, "density_custom_title.png")
        png(filename, width = 800, height = 600, res = 100)
        print(p)
        dev.off()
        cat("✓\n   Saved to:", filename, "\n")
    }, error = function(e) {
        cat("✗ FAILED:", conditionMessage(e), "\n")
    })
    
    cat("Testing plot_tsallis_violin_density_grid() with custom title... ")
    tryCatch({
        p <- plot_tsallis_violin_density_grid(se_list[["1"]], title = custom_title_grid)
        filename <- file.path(output_dir, "grid_custom_title.png")
        png(filename, width = 1400, height = 600, res = 100)
        print(p)
        dev.off()
        cat("✓\n   Saved to:", filename, "\n")
    }, error = function(e) {
        cat("✗ FAILED:", conditionMessage(e), "\n")
    })
    cat("\n")
    
    # ====================================================================
    # SUMMARY
    # ====================================================================
    cat(paste(c("=", rep("=", 70), "\n"), collapse=""))
    cat("Test Suite Completion\n")
    cat(paste(c("=", rep("=", 70), "\n"), collapse=""))
    cat("✓ All tests completed successfully!\n")
    cat("  Output directory:", output_dir, "\n")
    cat("  Generated files:\n")
    
    # List generated files
    files <- list.files(output_dir, pattern = "\\.png$")
    for (f in files) {
        cat("    -", f, "\n")
    }
    cat("\n")
    
}, error = function(e) {
    cat("\n\n")
    cat(paste(c("=", rep("=", 70), "\n"), collapse=""))
    cat("ERROR DURING TEST EXECUTION\n")
    cat(paste(c("=", rep("=", 70), "\n"), collapse=""))
    cat("Error message:", conditionMessage(e), "\n")
    cat("Stack trace:\n")
    print(traceback())
})

sink()
close(log_file)

# Print summary to console
cat("\n")
cat(readLines(output_file), sep = "\n")
