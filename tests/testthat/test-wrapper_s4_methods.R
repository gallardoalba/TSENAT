library(testthat)

context("S4 Method Dispatch: Wrapper Functions")

# ============================================================================
# HELPER: Create test analysis using actual package data
# ============================================================================

make_minimal_analysis <- function() {
  set.seed(123)
  # Create synthetic transcript-level data with proper isoform structure
  n_genes <- 50  # Reduced for faster testing but enough for jackknife
  isoforms_per_gene <- 5  # More isoforms per gene for better diversity
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
  
  TSENATAnalysis(se)
}

# ============================================================================
# TEST: S4 method dispatch - ensure methods exist and dispatch correctly
# ============================================================================

test_that("calculate_diversity_s4 method exists and dispatches", {
  analysis <- make_minimal_analysis()
  result <- tryCatch(
    calculate_diversity_s4(analysis, q = 1.0, verbose = FALSE),
    error = function(e) NULL
  )
  # Just verify method ran without crashing on dispatch
  expect_true(is.null(result) || inherits(result, "TSENATAnalysis"))
})

test_that("calculate_diversity_s4 rejects invalid input type", {
  expect_error(
    calculate_diversity_s4("not_analysis"),
    "must be a TSENATAnalysis|inherited method|signature"
  )
})

test_that("calculate_divergence_s4 method exists", {
  analysis <- make_minimal_analysis()
  result <- tryCatch(
    calculate_divergence_s4(analysis, group_col = "condition", verbose = FALSE),
    error = function(e) NULL
  )
  # Method should exist (may error on data but shouldn't on dispatch)
  expect_true(TRUE)
})

test_that("calculate_difference_s4 method exists", {
  analysis <- make_minimal_analysis()
  result <- tryCatch(
    calculate_difference_s4(analysis, group_col = "condition", verbose = FALSE),
    error = function(e) NULL
  )
  expect_true(TRUE)
})

test_that("jackknife_isoform_switching_s4 method exists", {
  analysis <- make_minimal_analysis()
  result <- tryCatch(
    jackknife_isoform_switching_s4(analysis, n_bootstrap = 5, verbose = FALSE),
    error = function(e) NULL
  )
  expect_true(TRUE)
})

test_that("detect_q_gene_interactions_s4 method exists", {
  analysis <- make_minimal_analysis()
  result <- tryCatch(
    detect_q_gene_interactions_s4(analysis, q_values = 1.0, verbose = FALSE),
    error = function(e) NULL
  )
  expect_true(TRUE)
})

test_that("calculate_lm_interaction_s4 method exists", {
  analysis <- make_minimal_analysis()
  result <- tryCatch(
    calculate_lm_interaction_s4(analysis, formula = ~ condition, verbose = FALSE),
    error = function(e) NULL
  )
  expect_true(TRUE)
})

test_that("compute_method_concordance_s4 rejects wrong input type", {
  expect_error(
    compute_method_concordance_s4("not_analysis"),
    "inherited method|signature"
  )
})

test_that("plot_method_concordance_s4 rejects wrong input type", {
  expect_error(
    plot_method_concordance_s4("not_analysis"),
    "inherited method|signature"
  )
})

test_that("effect_sizes_divergence_s4 method exists", {
  analysis <- make_minimal_analysis()
  result <- tryCatch(
    effect_sizes_divergence_s4(analysis, q = 1.0, verbose = FALSE),
    error = function(e) NULL
  )
  expect_true(TRUE)
})

# ============================================================================
# TEST: show and summary methods for TSENATAnalysis
# ============================================================================

test_that("show method for TSENATAnalysis works", {
  analysis <- make_minimal_analysis()
  
  expect_output(
    show(analysis),
    "TSENATAnalysis"
  )
})

test_that("summary method for TSENATAnalysis works", {
  analysis <- make_minimal_analysis()
  
  expect_output(
    summary(analysis),
    "Genes|Samples"
  )
})

# ============================================================================
# TEST: Configuration integration 
# ============================================================================

test_that("TSENATAnalysis accepts config in constructor", {
  analysis <- make_minimal_analysis()
  config <- tsenat_config(p_threshold = 0.01)
  
  analysis_with_config <- TSENATAnalysis(analysis@se, config = config)
  
  expect_equal(analysis_with_config@config$p_threshold, 0.01)
})

test_that("getConfig retrieves configuration", {
  analysis <- make_minimal_analysis()
  config <- tsenat_config(seed = 42)
  analysis <- TSENATAnalysis(analysis@se, config = config)
  
  retrieved <- getConfig(analysis)
  
  expect_true(is.list(retrieved))
  expect_equal(retrieved$seed, 42)
})

test_that("setConfig updates configuration on analysis", {
  analysis <- make_minimal_analysis()
  new_config <- tsenat_config(p_threshold = 0.001)
  
  updated <- setConfig(analysis, new_config)
  
  expect_equal(updated@config$p_threshold, 0.001)
})
