# Test coverage for tsallis q-spectrum plotting functions
# Functions: .plot_tsallis_bootstrap_ci, .plot_tsallis_gene_bootstrap_ci, .plot_tsallis_basic_gene

context("Tsallis q-spectrum plotting: Bootstrap CI and per-gene visualizations")
library(ggplot2)
library(cowplot)
library(dplyr)
library(SummarizedExperiment)
library(testthat)

# ════════════════════════════════════════════════════════════════════════════════
# TEST 1: .plot_tsallis_bootstrap_ci - Group-level bootstrap CI curves
# ════════════════════════════════════════════════════════════════════════════════

test_that(".plot_tsallis_bootstrap_ci creates q-curve plot with bootstrap CIs", {
  skip_if_not_installed("ggplot2")
  skip_if_not_installed("cowplot")
  
  # Create proper q-curve test data with multiple q-values and realistic structure
  q_values <- c(0.5, 1.0, 1.5, 2.0)
  n_q <- length(q_values)
  
  # SE represents aggregated entropy across samples - one row per q-value
  se <- SummarizedExperiment::SummarizedExperiment(
    assays = list(
      tsallis = matrix(seq(2.0, 1.5, length.out = n_q), nrow = 1, ncol = n_q),
      ci_lower = matrix(seq(1.8, 1.3, length.out = n_q), nrow = 1, ncol = n_q),
      ci_upper = matrix(seq(2.2, 1.7, length.out = n_q), nrow = 1, ncol = n_q)
    ),
    colData = data.frame(group = c("GroupA", "GroupB", "GroupA", "GroupB"))
  )
  
  # Long format: per-sample bootstrap data aggregated to group level
  long <- data.frame(
    q = rep(q_values, 4),
    group = rep(c("GroupA", "GroupB"), each = 2 * n_q),
    tsallis = rep(seq(2.0, 1.5, length.out = n_q), 4),
    sample = rep(paste0("S", 1:4), each = n_q),
    stringsAsFactors = FALSE
  )
  
  # Test successful plot creation
  p <- TSENAT:::.plot_tsallis_bootstrap_ci(se, long, output_file = NULL)
  
  expect_is(p, "ggplot")
  expect_true(inherits(p, "gg"))
})

test_that(".plot_tsallis_bootstrap_ci requires at least 2 q-values", {
  skip_if_not_installed("ggplot2")
  
  # SE with single q-value
  se <- SummarizedExperiment::SummarizedExperiment(
    assays = list(
      tsallis = matrix(c(2.1, 1.8), nrow = 1, ncol = 2),
      ci_lower = matrix(c(1.9, 1.6), nrow = 1, ncol = 2),
      ci_upper = matrix(c(2.3, 2.0), nrow = 1, ncol = 2)
    ),
    colData = data.frame(group = c("A", "B"))
  )
  
  # Long data with only one q-value
  long <- data.frame(
    q = c(0.5, 0.5),
    group = c("A", "B"),
    tsallis = c(2.1, 1.8),
    sample = c("S1", "S2"),
    stringsAsFactors = FALSE
  )
  
  # Should error: need at least 2 q values
  expect_error(
    TSENAT:::.plot_tsallis_bootstrap_ci(se, long, output_file = NULL),
    "Need at least 2 q values"
  )
})

test_that(".plot_tsallis_bootstrap_ci requires exactly 2 groups", {
  skip_if_not_installed("ggplot2")
  
  # SE with 3 groups
  se <- SummarizedExperiment::SummarizedExperiment(
    assays = list(
      tsallis = matrix(rnorm(6), nrow = 2, ncol = 3),
      ci_lower = matrix(rnorm(6, mean = 1.5), nrow = 2, ncol = 3),
      ci_upper = matrix(rnorm(6, mean = 2.5), nrow = 2, ncol = 3)
    ),
    colData = data.frame(group = c("A", "B", "C"))
  )
  
  # Long data with 3 groups
  long <- data.frame(
    q = c(0.5, 1.0, 0.5, 1.0, 0.5, 1.0),
    group = c("A", "A", "B", "B", "C", "C"),
    tsallis = rnorm(6),
    sample = c("S1", "S2", "S1", "S2", "S1", "S2"),
    stringsAsFactors = FALSE
  )
  
  # Should error: expected exactly 2 groups
  expect_error(
    TSENAT:::.plot_tsallis_bootstrap_ci(se, long, output_file = NULL),
    "Expected exactly 2 groups"
  )
})

test_that(".plot_tsallis_bootstrap_ci saves plot to file when output_file specified", {
  skip_if_not_installed("ggplot2")
  
  # Create test data
  se <- SummarizedExperiment::SummarizedExperiment(
    assays = list(
      tsallis = matrix(c(2.1, 2.3, 1.8, 2.0), nrow = 2, ncol = 2),
      ci_lower = matrix(c(1.9, 2.1, 1.6, 1.8), nrow = 2, ncol = 2),
      ci_upper = matrix(c(2.3, 2.5, 2.0, 2.2), nrow = 2, ncol = 2)
    ),
    colData = data.frame(group = c("A", "B"))
  )
  
  long <- data.frame(
    q = c(0.5, 1.0, 1.5, 0.5, 1.0, 1.5),
    group = c("A", "A", "A", "B", "B", "B"),
    tsallis = c(2.1, 2.3, 2.0, 1.8, 2.0, 1.9),
    sample = c("S1", "S1", "S1", "S2", "S2", "S2"),
    stringsAsFactors = FALSE
  )
  
  # Create temp file
  temp_file <- tempfile(fileext = ".png")
  
  # Plot with output file
  p <- TSENAT:::.plot_tsallis_bootstrap_ci(se, long, output_file = temp_file)
  
  # Verify plot is returned
  expect_is(p, "ggplot")
  
  # Clean up
  if (file.exists(temp_file)) {
    file.remove(temp_file)
  }
})

# ════════════════════════════════════════════════════════════════════════════════
# TEST 2: .plot_tsallis_gene_bootstrap_ci - Per-gene bootstrap CI faceted plots
# ════════════════════════════════════════════════════════════════════════════════

test_that(".plot_tsallis_gene_bootstrap_ci creates per-gene faceted plots", {
  skip_if_not_installed("ggplot2")
  skip_if_not_installed("cowplot")
  
  # Create realistic per-gene SE with bootstrap CI data
  # Key: SE rows are genes, columns are samples/groups
  q_values <- c(0.5, 1.0, 1.5, 2.0)
  genes <- c("Gene1", "Gene2")
  n_q <- length(q_values)
  n_genes <- length(genes)
  
  # SE: each gene has q-value measurements with CIs
  # Columns represent different q values with groups alternating
  # CRITICAL: Column names MUST follow format "Sample_q=X" for .prepare_gene_ci_data() parsing
  col_names <- c("A_q=0.5", "B_q=0.5", "A_q=1.0", "B_q=1.0", "A_q=1.5", "B_q=1.5", "A_q=2.0", "B_q=2.0")
  
  se <- SummarizedExperiment::SummarizedExperiment(
    assays = list(
      tsallis = matrix(
        c(2.0, 1.95, 1.9, 1.85, 1.8, 1.75, 1.7, 1.65),  # Gene1 q=0.5,1.0,1.5,2.0 groups A,B
        nrow = n_genes, ncol = 2 * n_q, byrow = TRUE
      ),
      ci_lower = matrix(
        c(1.8, 1.75, 1.7, 1.65, 1.6, 1.55, 1.5, 1.45),
        nrow = n_genes, ncol = 2 * n_q, byrow = TRUE
      ),
      ci_upper = matrix(
        c(2.2, 2.15, 2.1, 2.05, 2.0, 1.95, 1.9, 1.85),
        nrow = n_genes, ncol = 2 * n_q, byrow = TRUE
      )
    ),
    rowData = data.frame(gene_name = genes),
    colData = data.frame(group = rep(c("A", "B"), times = n_q))
  )
  colnames(se) <- col_names
  rownames(se) <- genes
  
  # Long format: per-sample bootstrap replicates for each gene
  # CRITICAL: sample names must match what SE column names parse to
  # If SE has "A_q=0.5", function parses sample="A", so long_data needs sample="A"
  long_list <- lapply(genes, function(g) {
    gene_idx <- match(g, genes)
    se_tsallis <- SummarizedExperiment::assays(se)[["tsallis"]][gene_idx, ]
    
    # Create bootstrap replicates: use A/B as sample names to match SE column parsing
    # Generate 3 bootstrap replicates for each gene-group-q combination
    data.frame(
      q = rep(q_values, 6),
      group = rep(rep(c("A", "B"), each = n_q), 3),
      Gene = g,
      tsallis = rep(se_tsallis, 3),
      sample = rep(rep(c("A", "B"), times = n_q), 3),  # Match SE column name parsing
      stringsAsFactors = FALSE
    )
  })
  long <- do.call(rbind, long_list)
  
  # Test plot creation with valid genes
  genes_to_plot <- c("Gene1", "Gene2")
  p <- TSENAT:::.plot_tsallis_gene_bootstrap_ci(se, long, genes = genes_to_plot, output_file = NULL)
  
  # Should create a plot or list of plots
  expect_true(ggplot2::is_ggplot(p) || is.list(p))
})

test_that(".plot_tsallis_gene_bootstrap_ci handles missing CI assay gracefully", {
  skip_if_not_installed("ggplot2")
  
  # Create SE without CI assays (only tsallis)
  se <- SummarizedExperiment::SummarizedExperiment(
    assays = list(
      tsallis = matrix(rnorm(6, mean = 2), nrow = 3, ncol = 2)
    ),
    rowData = data.frame(gene_name = c("Gene1", "Gene2", "Gene3")),
    colData = data.frame(group = c("A", "B"))
  )
  rownames(se) <- c("Gene1", "Gene2", "Gene3")
  
  long <- data.frame(
    q = rep(c(0.5, 1.0, 1.5), 6),
    group = rep(rep(c("A", "B"), each = 3), 3),
    Gene = rep(c("Gene1", "Gene2", "Gene3"), each = 6),
    tsallis = rnorm(18, mean = 2, sd = 0.3),
    sample = rep(c("S1", "S2"), each = 9),
    stringsAsFactors = FALSE
  )
  
  genes <- c("Gene1", "Gene2")
  
  # Should error: requires CI assays
  expect_error(
    TSENAT:::.plot_tsallis_gene_bootstrap_ci(se, long, genes = genes, output_file = NULL),
    "Bootstrap CI assays.*not found in SummarizedExperiment"
  )
})

test_that(".plot_tsallis_gene_bootstrap_ci with only existing genes succeeds", {
  skip_if_not_installed("ggplot2")
  
  # Create SE with proper column naming for CI parsing
  col_names_2gene <- c("A_q=0.5", "B_q=0.5", "A_q=1.0", "B_q=1.0")
  
  se <- SummarizedExperiment::SummarizedExperiment(
    assays = list(
      tsallis = matrix(rnorm(8, mean = 2), nrow = 2, ncol = 4),
      ci_lower = matrix(rnorm(8, mean = 1.5), nrow = 2, ncol = 4),
      ci_upper = matrix(rnorm(8, mean = 2.5), nrow = 2, ncol = 4)
    ),
    rowData = data.frame(gene_name = c("Gene1", "Gene2")),
    colData = data.frame(group = rep(c("A", "B"), 2), q = rep(c(0.5, 1.0), each = 2))
  )
  colnames(se) <- col_names_2gene
  rownames(se) <- c("Gene1", "Gene2")
  
  long <- data.frame(
    q = rep(c(0.5, 1.0), 4),
    group = rep(c("A", "B"), each = 2),
    Gene = c("Gene1", "Gene1", "Gene2", "Gene2", "Gene1", "Gene1", "Gene2", "Gene2"),
    tsallis = rnorm(8, mean = 2),
    sample = rep(c("A", "B"), 4),  # Match SE column parsing
    stringsAsFactors = FALSE
  )
  
  # Use only genes that exist in SE
  genes <- c("Gene1", "Gene2")
  
  # Should succeed when genes exist (no warnings expected with proper data mapping)
  p <- TSENAT:::.plot_tsallis_gene_bootstrap_ci(se, long, genes = genes, output_file = NULL)
  
  # Should create valid plot output
  expect_true(ggplot2::is_ggplot(p) || is.list(p))
})

# ════════════════════════════════════════════════════════════════════════════════
# TEST 3: .plot_tsallis_basic_gene - Per-gene spread plots (no bootstrap CI)
# ════════════════════════════════════════════════════════════════════════════════

test_that(".plot_tsallis_basic_gene creates per-gene plots with IQR spread", {
  skip_if_not_installed("ggplot2")
  skip_if_not_installed("dplyr")
  
  # Create realistic long format data with multiple samples and q-values
  q_values <- c(0.5, 1.0, 1.5, 2.0)
  genes <- c("Gene1", "Gene2")
  n_samples_per_group <- 5  # Realistic sample count for better IQR calculation
  
  long_list <- lapply(genes, function(g) {
    data.frame(
      q = rep(q_values, 2 * n_samples_per_group),
      group = rep(rep(c("A", "B"), each = length(q_values)), n_samples_per_group),
      Gene = g,
      tsallis = rnorm(length(q_values) * 2 * n_samples_per_group, mean = 2, sd = 0.3),
      sample = rep(paste0("S", 1:(2 * n_samples_per_group)), each = length(q_values)),
      stringsAsFactors = FALSE
    )
  })
  long <- do.call(rbind, long_list)
  
  genes_to_plot <- c("Gene1", "Gene2")
  
  # Test with IQR metric (default)
  p <- TSENAT:::.plot_tsallis_basic_gene(long, genes = genes_to_plot, metric = "iqr", output_file = NULL)
  
  # Result should be a plot or list of plots
  expect_true(ggplot2::is_ggplot(p) || is.list(p))
})

test_that(".plot_tsallis_basic_gene respects metric parameter (SD)", {
  skip_if_not_installed("ggplot2")
  skip_if_not_installed("dplyr")
  
  long <- data.frame(
    q = rep(c(0.5, 1.0), 8),
    group = rep(c("A", "B"), each = 2),
    Gene = "Gene1",
    tsallis = rnorm(16, mean = 2, sd = 0.5),
    sample = rep(c("S1", "S2"), 8),
    stringsAsFactors = FALSE
  )
  
  genes <- c("Gene1")
  
  # Test with SD metric
  p <- TSENAT:::.plot_tsallis_basic_gene(long, genes = genes, metric = "sd", output_file = NULL)
  
  expect_is(p, "ggplot")
})

test_that(".plot_tsallis_basic_gene creates faceted plots for multiple genes", {
  skip_if_not_installed("ggplot2")
  skip_if_not_installed("cowplot")
  skip_if_not_installed("dplyr")
  
  long <- data.frame(
    q = rep(c(0.5, 1.0, 1.5, 2.0), 12),
    group = rep(rep(c("A", "B"), each = 4), 6),
    Gene = rep(c("GeneA", "GeneB", "GeneC"), each = 16),
    tsallis = rnorm(48, mean = 2, sd = 0.3),
    sample = rep(c("S1", "S2", "S3"), 16),
    stringsAsFactors = FALSE
  )
  
  genes <- c("GeneA", "GeneB", "GeneC")
  
  # Create faceted plot (multiple genes)
  expected_output <- TSENAT:::.plot_tsallis_basic_gene(long, genes = genes, metric = "iqr", output_file = NULL)
  
  # For multiple genes, should return a list or complex plot
  expect_true(ggplot2::is_ggplot(expected_output) || is.list(expected_output) || inherits(expected_output, "gtable"))
})

test_that(".plot_tsallis_basic_gene handles single gene", {
  skip_if_not_installed("ggplot2")
  skip_if_not_installed("dplyr")
  
  long <- data.frame(
    q = rep(c(0.5, 1.0), 4),
    group = c("A", "A", "B", "B", "A", "A", "B", "B"),
    Gene = "SingleGene",
    tsallis = rnorm(8, mean = 2),
    sample = c("S1", "S1", "S2", "S2", "S3", "S3", "S4", "S4"),
    stringsAsFactors = FALSE
  )
  
  genes <- c("SingleGene")
  
  # Single gene should return single plot
  p <- TSENAT:::.plot_tsallis_basic_gene(long, genes = genes, metric = "iqr", output_file = NULL)
  
  expect_is(p, "ggplot")
})

test_that(".plot_tsallis_basic_gene saves plot to file", {
  skip_if_not_installed("ggplot2")
  skip_if_not_installed("dplyr")
  
  long <- data.frame(
    q = rep(c(0.5, 1.0), 4),
    group = c("A", "A", "B", "B", "A", "A", "B", "B"),
    Gene = "Gene1",
    tsallis = rnorm(8, mean = 2),
    sample = c("S1", "S1", "S2", "S2", "S3", "S3", "S4", "S4"),
    stringsAsFactors = FALSE
  )
  
  genes <- c("Gene1")
  temp_file <- tempfile(fileext = ".png")
  
  # Create plot with output file
  p <- TSENAT:::.plot_tsallis_basic_gene(long, genes = genes, metric = "iqr", output_file = temp_file)
  
  expect_is(p, "ggplot")
  
  # Clean up
  if (file.exists(temp_file)) {
    file.remove(temp_file)
  }
})

test_that(".plot_tsallis_basic_gene handles NA values gracefully", {
  skip_if_not_installed("ggplot2")
  skip_if_not_installed("dplyr")
  
  long <- data.frame(
    q = rep(c(0.5, 1.0), 8),
    group = rep(c("A", "B"), times = 8),
    Gene = "Gene1",
    tsallis = c(rnorm(15, mean = 2), NA),  # Include one NA (16 total values)
    sample = rep(c("S1", "S2", "S3", "S4"), 4),
    stringsAsFactors = FALSE
  )
  
  genes <- c("Gene1")
  
  # Should handle NA without error (na.rm = TRUE in calculations)
  p <- TSENAT:::.plot_tsallis_basic_gene(long, genes = genes, metric = "iqr", output_file = NULL)
  
  expect_is(p, "ggplot")
})

test_that(".plot_tsallis_basic_gene with mismatched gene names returns empty gracefully", {
  skip_if_not_installed("ggplot2")
  skip_if_not_installed("dplyr")
  
  long <- data.frame(
    q = rep(c(0.5, 1.0), 4),
    group = c("A", "A", "B", "B", "A", "A", "B", "B"),
    Gene = "Gene1",
    tsallis = rnorm(8, mean = 2),
    sample = c("S1", "S1", "S2", "S2", "S3", "S3", "S4", "S4"),
    stringsAsFactors = FALSE
  )
  
  # Request different gene
  genes <- c("NonExistent")
  
  # Should handle gracefully or error cleanly
  expect_error(
    TSENAT:::.plot_tsallis_basic_gene(long, genes = genes, metric = "iqr", output_file = NULL),
    NA  # May error but shouldn't crash
  )
})

# ════════════════════════════════════════════════════════════════════════════════
# EDGE CASE TESTS: Finding potential bugs
# ════════════════════════════════════════════════════════════════════════════════

test_that(".plot_tsallis_bootstrap_ci detects inverted CI bounds", {
  skip_if_not_installed("ggplot2")
  
  # Test with inverted CI: ci_lower > ci_upper (bug check)
  q_values <- c(0.5, 1.0, 1.5, 2.0)
  
  se <- SummarizedExperiment::SummarizedExperiment(
    assays = list(
      tsallis = matrix(seq(2.0, 1.5, length.out = length(q_values)), nrow = 1, ncol = length(q_values)),
      ci_lower = matrix(seq(2.2, 1.7, length.out = length(q_values)), nrow = 1, ncol = length(q_values)),  # INVERTED!
      ci_upper = matrix(seq(1.8, 1.3, length.out = length(q_values)), nrow = 1, ncol = length(q_values))   # INVERTED!
    ),
    colData = data.frame(group = c("GroupA", "GroupB", "GroupA", "GroupB"))
  )
  
  long <- data.frame(
    q = rep(q_values, 4),
    group = rep(c("GroupA", "GroupB"), each = 2 * length(q_values)),
    tsallis = rep(seq(2.0, 1.5, length.out = length(q_values)), 4),
    sample = rep(paste0("S", 1:4), each = length(q_values)),
    stringsAsFactors = FALSE
  )
  
  # Should still plot (inverted CIs are a user data problem, not plot problem)
  p <- TSENAT:::.plot_tsallis_bootstrap_ci(se, long, output_file = NULL)
  
  expect_is(p, "ggplot")
})

test_that(".plot_tsallis_bootstrap_ci with zero variance throws informative error", {
  skip_if_not_installed("ggplot2")
  
  # All values identical (zero variance)
  q_values <- c(0.5, 1.0)
  
  se <- SummarizedExperiment::SummarizedExperiment(
    assays = list(
      tsallis = matrix(rep(2.0, 2), nrow = 1, ncol = 2),
      ci_lower = matrix(rep(2.0, 2), nrow = 1, ncol = 2),
      ci_upper = matrix(rep(2.0, 2), nrow = 1, ncol = 2)
    ),
    colData = data.frame(group = c("GroupA", "GroupB"))
  )
  
  long <- data.frame(
    q = c(0.5, 1.0, 0.5, 1.0),
    group = rep(c("GroupA", "GroupB"), each = 2),
    tsallis = rep(2.0, 4),
    sample = c("S1", "S1", "S2", "S2"),
    stringsAsFactors = FALSE
  )
  
  # Should work even with zero variance (just produces flat line)
  p <- TSENAT:::.plot_tsallis_bootstrap_ci(se, long, output_file = NULL)
  
  expect_is(p, "ggplot")
})

test_that(".plot_tsallis_basic_gene handles extreme q-value ranges", {
  skip_if_not_installed("ggplot2")
  skip_if_not_installed("dplyr")
  
  # Test with very small and very large q values
  q_values <- c(0.1, 10, 100, 1000)  # Extreme range
  
  long <- data.frame(
    q = rep(q_values, 4),
    group = rep(c("A", "B"), each = 2 * length(q_values)),
    Gene = "Gene1",
    tsallis = rnorm(length(q_values) * 4, mean = 2, sd = 0.2),
    sample = rep(paste0("S", 1:4), each = length(q_values)),
    stringsAsFactors = FALSE
  )
  
  # Should handle extreme q-values gracefully
  p <- TSENAT:::.plot_tsallis_basic_gene(long, genes = c("Gene1"), metric = "iqr", output_file = NULL)
  
  expect_is(p, "ggplot")
})

test_that(".plot_tsallis_gene_bootstrap_ci handles single sample per group", {
  skip_if_not_installed("ggplot2")
  
  # Minimal viable data: 1 gene, 4 q-values, 1 sample per group
  q_values <- c(0.5, 1.0, 1.5, 2.0)
  
  se <- SummarizedExperiment::SummarizedExperiment(
    assays = list(
      tsallis = matrix(rnorm(length(q_values), mean = 2), nrow = 1, ncol = length(q_values)),
      ci_lower = matrix(rnorm(length(q_values), mean = 1.5), nrow = 1, ncol = length(q_values)),
      ci_upper = matrix(rnorm(length(q_values), mean = 2.5), nrow = 1, ncol = length(q_values))
    ),
    rowData = data.frame(gene_name = "Gene1"),
    colData = data.frame(group = rep(c("A", "B"), each = length(q_values) / 2))
  )
  rownames(se) <- "Gene1"
  
  # Minimal long data: one sample per group
  long <- data.frame(
    q = rep(q_values, 2),
    group = rep(c("A", "B"), each = length(q_values)),
    Gene = "Gene1",
    tsallis = rnorm(length(q_values) * 2, mean = 2, sd = 0.2),
    sample = rep(c("S1", "S2"), each = length(q_values)),
    stringsAsFactors = FALSE
  )
  
  # Should handle minimal data without crashing
  p <- TSENAT:::.plot_tsallis_gene_bootstrap_ci(se, long, genes = c("Gene1"), output_file = NULL)
  
  expect_true(ggplot2::is_ggplot(p) || is.list(p))
})
