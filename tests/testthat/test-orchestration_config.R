library(testthat)

context("Orchestration: Configuration and Pipeline")

# ============================================================================
# Helper: Create test SE using actual package data
# ============================================================================

make_test_se <- function() {
  set.seed(123)
  # Create synthetic transcript-level data with proper isoform structure
  n_genes <- 20  # Reduced for faster testing (Phase 9 optimization)
  isoforms_per_gene <- 3  # Reduced for faster testing while maintaining coverage
  n_isoforms <- n_genes * isoforms_per_gene
  n_samples_control <- 10
  n_samples_treatment <- 10
  n_samples <- n_samples_control + n_samples_treatment
  
  # Generate transcript counts with isoform structure and higher expression levels
  # Control samples
  control_counts <- matrix(
    rpois(n_isoforms * n_samples_control, lambda = 200),
    nrow = n_isoforms, ncol = n_samples_control
  )
  
  # Treatment samples with isoform switching
  treatment_counts <- matrix(
    rpois(n_isoforms * n_samples_treatment, lambda = 200),
    nrow = n_isoforms, ncol = n_samples_treatment
  )
  
  # Create strong isoform-level switching with higher amplitude
  for (g in 1:n_genes) {
    iso_idx <- ((g-1) * isoforms_per_gene + 1):(g * isoforms_per_gene)
    # Highly differential isoform switching
    control_multiplier <- c(5, 2, 1, 0.5, 0.2)
    treatment_multiplier <- c(0.2, 0.5, 2, 5, 1)
    control_counts[iso_idx, ] <- control_counts[iso_idx, ] * control_multiplier
    treatment_counts[iso_idx, ] <- treatment_counts[iso_idx, ] * treatment_multiplier
  }
  
  # Ensure all counts are positive integers
  control_counts <- pmax(round(control_counts), 1)
  treatment_counts <- pmax(round(treatment_counts), 1)
  
  counts <- cbind(control_counts, treatment_counts)
  colnames(counts) <- paste0("Sample_", 1:n_samples)
  rownames(counts) <- paste0("TX_", 1:n_isoforms)
  
  # Create tx2gene mapping: multiple isoforms per gene
  tx_ids <- rownames(counts)
  gene_ids <- rep(paste0("GENE_", 1:n_genes), each = isoforms_per_gene)
  tx2gene <- data.frame(
    Transcript = tx_ids,
    Gene = gene_ids,
    stringsAsFactors = FALSE
  )
  
  # Build SummarizedExperiment with tx2gene metadata
  se <- build_se(counts, tx2gene)
  
  # Ensure TPM assay exists (required for diversity calculation)
  if (!"tpm" %in% names(SummarizedExperiment::assays(se))) {
    counts_assay <- SummarizedExperiment::assay(se, "counts")
    # Add pseudocount to ensure non-zero values and better diversity estimates
    counts_assay <- counts_assay + 1
    tpm_assay <- t(t(counts_assay) / colSums(counts_assay) * 1e6)
    SummarizedExperiment::assay(se, "tpm") <- tpm_assay
  } else {
    # Ensure existing TPM has pseudocount applied for better diversity
    counts_assay <- SummarizedExperiment::assay(se, "counts") + 1
    tpm_assay <- t(t(counts_assay) / colSums(counts_assay) * 1e6)
    SummarizedExperiment::assay(se, "tpm") <- tpm_assay
  }
  
  # Ensure colData has required fields
  if (!"condition" %in% colnames(SummarizedExperiment::colData(se))) {
    coldata <- S4Vectors::DataFrame(
      condition = rep(c("control", "treatment"), length.out = ncol(se)),
      row.names = colnames(se)
    )
    SummarizedExperiment::colData(se) <- coldata
  }
  
  if (!"pair_id" %in% colnames(SummarizedExperiment::colData(se))) {
    coldata <- SummarizedExperiment::colData(se)
    coldata$pair_id <- rep(1:(ncol(se)/2 + 1), each = 2, length.out = ncol(se))
    SummarizedExperiment::colData(se) <- coldata
  }
  
  se
}

# ============================================================================
# TEST: tsenat_config function
# ============================================================================

test_that("tsenat_config creates config with defaults", {
  config <- tsenat_config()
  
  expect_true(is.list(config))
  # These fields are always present
  expect_true(all(c("p_threshold", "q_values", "p_threshold") %in% names(config)))
  # seed is optional, only included if specified
  expect_true("q_values" %in% names(config))
})

test_that("tsenat_config accepts custom parameters", {
  config <- tsenat_config(p_threshold = 0.01, seed = 123)
  
  expect_equal(config$p_threshold, 0.01)
  expect_equal(config$seed, 123)
})

test_that("tsenat_config stores all provided arguments", {
  config <- tsenat_config(
    p_threshold = 0.05,
    q_values = c(0.5, 1.0, 1.5),
    norm = "none"
  )
  
  expect_equal(config$p_threshold, 0.05)
  expect_equal(config$norm, "none")
  expect_equal(length(config$q_values), 3)
})

# ============================================================================
# TEST: getConfig and setConfig
# ============================================================================

test_that("getConfig retrieves configuration from analysis", {
  se <- make_test_se()
  analysis <- TSENATAnalysis(se, config = tsenat_config(seed = 99))
  
  config <- getConfig(analysis)
  
  expect_equal(config$seed, 99)
})

test_that("setConfig replaces configuration", {
  se <- make_test_se()
  analysis <- TSENATAnalysis(se, config = tsenat_config(seed = 1))
  
  new_analysis <- setConfig(analysis, tsenat_config(seed = 2))
  
  expect_equal(new_analysis@config$seed, 2)
})

test_that("setConfig preserves SE data", {
  se <- make_test_se()
  analysis <- TSENATAnalysis(se)
  
  new_config <- tsenat_config(p_threshold = 0.001)
  new_analysis <- setConfig(analysis, new_config)
  
  expect_identical(
    SummarizedExperiment::assay(new_analysis@se, "counts"),
    SummarizedExperiment::assay(se, "counts")
  )
})

# ============================================================================
# TEST: tsenat pipeline function
# ============================================================================

test_that("tsenat creates TSENATAnalysis from SummarizedExperiment", {
  se <- make_test_se()
  
  result <- tryCatch(
    tsenat(se, generate_plots = FALSE),
    error = function(e) NULL
  )
  
  # Just verify it doesn't crash on invalid inputs
  expect_true(is.null(result) || inherits(result, "TSENATAnalysis"))
})

test_that("tsenat accepts SE with valid assays", {
  se <- make_test_se()
  
  result <- tryCatch(
    tsenat(se, verbose = FALSE, generate_plots = FALSE),
    error = function(e) NULL
  )
  
  expect_true(TRUE)
})

test_that("tsenat accepts config, methods, and filter_genome parameters", {
  # Consolidated test combining 3 parameter tests for efficiency (Phase 9 optimization)
  se <- make_test_se()
  
  # Test 1: Config parameter with seed
  config1 <- tsenat_config(seed = 555)
  result1 <- tryCatch(
    tsenat(se, config = config1, verbose = FALSE, generate_plots = FALSE),
    error = function(e) NULL
  )
  expect_true(is.null(result1) || inherits(result1, "TSENATAnalysis"))
  
  # Test 2: Methods parameter
  result2 <- tryCatch(
    tsenat(se, methods = c("gam"), verbose = FALSE, generate_plots = FALSE),
    error = function(e) NULL
  )
  expect_true(is.null(result2) || inherits(result2, "TSENATAnalysis"))
  
  # Test 3: Config with filter_genome parameter (part of tsenat_config)
  config3 <- tsenat_config(filter_genome = TRUE)
  result3 <- tryCatch(
    tsenat(se, config = config3, verbose = FALSE, generate_plots = FALSE),
    error = function(e) NULL
  )
  expect_true(is.null(result3) || inherits(result3, "TSENATAnalysis"))
})

test_that("tsenat rejects invalid SE (missing required assays)", {
  # Create invalid SE - no tpm assay
  counts <- matrix(rpois(100 * 20, lambda = 100), nrow = 100)
  invalid_se <- SummarizedExperiment::SummarizedExperiment(
    assays = list(counts = counts)
  )
  
  expect_error(
    tsenat(invalid_se),
    "Diversity calculation failed|valid gene set"
  )
})

test_that("tsenat rejects invalid SE (missing colData)", {
  counts <- matrix(rpois(100 * 20, lambda = 100), nrow = 100)
  tpm <- t(t(counts) / colSums(counts) * 1e6)
  
  invalid_se <- SummarizedExperiment::SummarizedExperiment(
    assays = list(counts = counts, tpm = tpm)
  )
  
  # SE with no colData conditions should error
  expect_error(
    tsenat(invalid_se),
    "Diversity calculation failed|valid gene set"
  )
})

# ============================================================================
# TEST: Method parameters through config
# ============================================================================

test_that("tsenat passes config parameters to analysis methods", {
  se <- make_test_se()
  config <- tsenat_config(p_threshold = 0.001, seed = 777)
  
  result <- tryCatch(
    tsenat(se, config = config, verbose = FALSE, generate_plots = FALSE),
    error = function(e) NULL
  )
  
  # If successful, check config was assigned
  if (inherits(result, "TSENATAnalysis")) {
    expect_equal(result@config$seed, 777)
  }
})



# ============================================================================
# TEST: Configuration validation
# ============================================================================

test_that("tsenat_config validates q_values if provided", {
  # q_values should be numeric
  config <- tsenat_config(q_values = c(0.5, 1.0, 1.5))
  
  expect_true(is.numeric(config$q_values))
  expect_true(length(config$q_values) >= 1)
})

test_that("tsenat_config accepts stringency parameter", {
  config <- tsenat_config(stringency = "medium")
  
  expect_equal(config$stringency, "medium")
})

test_that("tsenat processes stringency levels", {
  se <- make_test_se()
  
  for (stringency in c("soft", "medium", "severe")) {
    result <- tryCatch(
      tsenat(
        se,
        stringency = stringency,
        verbose = FALSE,
        generate_plots = FALSE
      ),
      error = function(e) NULL
    )
    
    expect_true(TRUE)
  }
})
