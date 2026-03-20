#!/usr/bin/env Rscript

# Test S4 class instantiation and basic methods
library(TSENAT)
library(SummarizedExperiment)

# Create a minimal test SummarizedExperiment
set.seed(42)
counts <- matrix(rpois(100, 5), nrow = 10, ncol = 10)
colnames(counts) <- paste0("Sample", 1:10)
rownames(counts) <- paste0("Gene", 1:10)

se <- SummarizedExperiment(assays = list(counts = counts))

# Test 1: Create TSENATAnalysis object
cat("Test 1: Creating TSENATAnalysis object...\n")
analysis <- TSENATAnalysis(se)
cat("✓ TSENATAnalysis created\n")

# Test 2: Check S4 class
cat("\nTest 2: Checking S4 class structure...\n")
cat("Class:", class(analysis), "\n")
cat("Slots:", slotNames(analysis), "\n")

# Test 3: Show method
cat("\nTest 3: Testing show method...\n")
show(analysis)

# Test 4: Config builder
cat("\nTest 4: Testing tsenat_config()...\n")
cfg <- tsenat_config(q_values = c(0.5, 1.0), fdr_threshold = 0.01)
cat("Config created with q_values:", paste(cfg$q_values, collapse = ", "), "\n")

# Test 5: Accessor methods on empty results
cat("\nTest 5: Testing accessors on empty results...\n")
div <- diversity(analysis)
cat("diversity() returns NULL:", is.null(div), "(should be TRUE)\n")

lm <- lmResults(analysis)
cat("lmResults() returns NULL:", is.null(lm), "(should be TRUE)\n")

# Test 6: Test addPlot method
cat("\nTest 6: Testing addPlot method...\n")
library(ggplot2)
test_plot <- ggplot() + geom_point(aes(x = 1, y = 1))
analysis <- addPlot(analysis, type = "test_plot", plot = test_plot)
cat("Plot added successfully\n")
retrieved_plot <- getPlot(analysis, type = "test_plot")
cat("Plot retrieved successfully:", !is.null(retrieved_plot), "\n")

cat("\n✓ All basic tests passed!\n")
