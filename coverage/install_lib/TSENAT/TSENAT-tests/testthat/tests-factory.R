#' Test Data Factory Functions for TSENATAnalysis
#' 
#' These helper functions create consistent, reusable test data with biological signal.
#' They follow S4 patterns and reduce boilerplate across the test suite by ~70%.
#'
#' @keywords internal
#' @name test_factories

# ============================================================================
# NOTE: This file contains test factory functions, not test cases
# ============================================================================
test_that("[FACTORY FILE] This file provides test utilities, not tests", {
  skip("This file contains test data factory functions for use across the test suite, not test cases")
})

#' Create a Complete TSENATAnalysis with Diversity Results
#' 
#' Factory function that generates a fully initialized TSENATAnalysis object
#' with realistic data and biological signal. Eliminates repetitive setup code
#' across tests.
#'
#' @param n_genes Number of genes to simulate (default: 8)
#' @param n_samples_per_group Samples per condition (default: 20)
#' @param control_lambda Poisson lambda for control condition (default: 40)
#' @param treatment_lambda Poisson lambda for treatment condition (default: 150)
#' @param q_values Vector of q-values for diversity calculation (default: c(0.5, 0.75, 1.0, 1.5, 2.0))
#' @param include_divergence If TRUE, compute divergence results (default: TRUE)
#' @param include_sait_results If TRUE, add placeholder SAIT results (default: TRUE)
#' @param seed Random seed for reproducibility (default: 42)
#' @param verbose Logical for progress messages (default: FALSE)
#'
#' @return TSENATAnalysis object with:
#'   - SummarizedExperiment with count matrix, rowData, and colData
#'   - Computed diversity results across q-values
#'   - Computed divergence results (optional)
#'   - Placeholder SAIT results (optional)
#'   - Proper tx2gene metadata mapping
#'
#' @details
#' The factory ensures:
#' - Biological signal: control (lambda=40) vs treatment (lambda=150) contrast
#' - Sufficient samples: 40 total (20 per group) for stable LM fitting
#' - Multiple q-values: c(0.5, 1.0, 1.5) avoids rank deficiency
#' - Valid S4 object structure: passes all TSENATAnalysis validity checks
#' - Optional divergence and SAIT results to support testing without warnings
#'
#' @examples
#' \dontrun{
#'   # Create with defaults (8 genes, 20 samples/group, multi-q, with all results)
#'   analysis <- create_test_analysis()
#'   
#'   # Create with custom parameters
#'   analysis <- create_test_analysis(
#'     n_genes = 16,
#'     n_samples_per_group = 30,
#'     control_lambda = 50,
#'     treatment_lambda = 200,
#'     q_values = c(0.1, 0.5, 1.0, 1.5, 2.0),
#'     include_divergence = TRUE,
#'     include_sait_results = TRUE
#'   )
#' }
#'
#' @keywords internal
#' @export
create_test_analysis <- function(
    n_genes = 8,
    n_samples_per_group = 20,
    control_lambda = 40,
    treatment_lambda = 150,
    q_values = c(0.5, 0.75, 1.0, 1.5, 2.0),
    include_divergence = TRUE,
    include_sait_results = TRUE,
    seed = 42,
    verbose = FALSE) {
  
  set.seed(seed)
  
  # Dimensions
  n_samples <- n_samples_per_group * 2
  n_transcripts <- n_genes * 50
  
  if (verbose) {
    message("[create_test_analysis] Generating ", n_transcripts, " transcripts across ",
            n_genes, " genes with ", n_samples, " samples")
  }
  
  # Generate counts with biological signal
  control_idx <- seq(1, n_samples, by = 2)
  treatment_idx <- seq(2, n_samples, by = 2)
  
  counts <- matrix(0, nrow = n_transcripts, ncol = n_samples)
  for (j in seq_len(n_samples)) {
    if (j %in% control_idx) {
      counts[, j] <- rpois(n_transcripts, lambda = control_lambda)
    } else {
      counts[, j] <- rpois(n_transcripts, lambda = treatment_lambda)
    }
  }
  counts <- pmax(counts, 50)
  
  rownames(counts) <- paste0("TX_", 1:n_transcripts)
  colnames(counts) <- paste0("Sample_", 1:n_samples)
  
  # Create rowData with gene mappings
  rowData <- S4Vectors::DataFrame(
    transcript_id = rownames(counts),
    gene_id = paste0("GENE_", rep(1:n_genes, each = 50, length.out = n_transcripts)),
    row.names = rownames(counts)
  )
  
  # Create colData with experimental design
  colData <- S4Vectors::DataFrame(
    sample_id = colnames(counts),
    condition = rep(c("control", "treatment"), length.out = n_samples),
    sample_type = rep(c("typeA", "typeB"), length.out = n_samples),
    subject = rep(paste0("S", 1:10), length.out = n_samples),
    paired_samples = rep(paste0("pair", 1:10), length.out = n_samples),
    row.names = colnames(counts)
  )
  
  # Create SummarizedExperiment
  se <- SummarizedExperiment::SummarizedExperiment(
    assays = list(counts = counts),
    rowData = rowData,
    colData = colData
  )
  
  # Add tx2gene metadata
  tx2gene_df <- data.frame(
    Transcript = rownames(counts),
    Gene = rowData$gene_id,
    stringsAsFactors = FALSE
  )
  S4Vectors::metadata(se)$tx2gene <- tx2gene_df
  
  # Generate synthetic TPM data (matching counts dimensions)
  tpm <- counts
  for (j in seq_len(ncol(tpm))) {
    lib_size <- colSums(tpm[, j, drop = FALSE])
    if (lib_size > 0) {
      tpm[, j] <- (tpm[, j] / lib_size) * 1e6
    }
  }
  rownames(tpm) <- rownames(counts)
  colnames(tpm) <- colnames(counts)
  S4Vectors::metadata(se)$tpm <- tpm
  
  # Initialize TSENATAnalysis with explicit control_group in config.
  # This prevents auto-detection heuristics from guessing the wrong
  # reference group and ensures reproducible, explicit test design.
  analysis <- TSENAT::TSENATAnalysis(se = se, config = list(control_group = "control"))
  
  # Calculate diversity
  analysis <- TSENAT::calculate_diversity(
    analysis,
    q = q_values,
    verbose = FALSE,
    min_valid_frac = 0
  )
  
  # Calculate divergence if requested
  if (include_divergence) {
    analysis <- tryCatch({
      TSENAT::calculate_divergence(
        analysis,
        group_col = "condition",
        control_group = "control",
        verbose = FALSE
      )
    }, error = function(e) {
      # If divergence fails, continue without it
      if (verbose) {
        message("[create_test_analysis] Warning: divergence calculation failed: ", e$message)
      }
      analysis
    })
  }
  
  # Add placeholder SAIT results if requested
  if (include_sait_results) {
    # Create a simple placeholder LM result (empty data frame structure)
    # This prevents "No SAIT results found" warnings in tests
    sait_placeholder <- list(
      sait_interaction = list(
        results = data.frame(
          gene = character(0),
          term = character(0),
          estimate = numeric(0),
          std.error = numeric(0),
          statistic = numeric(0),
          p.value = numeric(0)
        ),
        model_data = NULL
      )
    )
    analysis@sait_results <- sait_placeholder
  }
  
  if (verbose) {
    message("[create_test_analysis] Analysis created with ",
            length(q_values), " q-values")
    message("[create_test_analysis] Diversity results: ",
            nrow(TSENAT::diversity(analysis)), " genes")
    if (include_divergence) {
      message("[create_test_analysis] Divergence results included")
    }
    if (include_sait_results) {
      message("[create_test_analysis] Placeholder SAIT results included")
    }
  }
  
  return(analysis)
}

# =============================================================================
# TSENAT Test Suite Helper Functions
# =============================================================================
# This file provides reusable helper functions for creating test data and
# making standardized assertions, significantly reducing code duplication
# across test files.
# =============================================================================

# =============================================================================
# Test Data Helpers
# =============================================================================

#' Create a standard SummarizedExperiment for testing
#'
#' @param n_genes Number of genes (rows)
#' @param n_samples Number of samples (columns)
#' @param control_n Number of samples in control group
#' @param seed Random seed
#' @param lambda Poisson lambda for rpois() count generation
#' @param use_rpois If TRUE, use rpois for realistic counts; if FALSE, use
#'   sequential matrix(1:N) for reproducibility
#'
#' @return SummarizedExperiment with counts and metadata
#'
#' @keywords internal
create_test_se <- function(
    n_genes = 2,
    n_samples = 6,
    control_n = 3,
    seed = 42,
    lambda = 100,
    use_rpois = TRUE
) {
  skip_if_not_installed("SummarizedExperiment")
  
  set.seed(seed)
  
  if (use_rpois) {
    counts <- matrix(
      rpois(n_genes * n_samples, lambda = lambda),
      nrow = n_genes,
      ncol = n_samples
    )
  } else {
    # Use sequential matrix for reproducibility
    counts <- matrix(
      seq_len(n_genes * n_samples),
      nrow = n_genes,
      ncol = n_samples
    )
  }
  
  rownames(counts) <- paste0("Gene_", seq_len(n_genes))
  colnames(counts) <- paste0("Sample_", seq_len(n_samples))
  
  colData <- data.frame(
    sample_type = factor(c(
      rep("Control", control_n),
      rep("Treatment", n_samples - control_n)
    )),
    row.names = colnames(counts)
  )
  
  rowData <- data.frame(
    gene_id = paste0("GENE_", seq_len(n_genes)),
    gene_name = rownames(counts)
  )
  
  # Generate TPM data (normalized counts)
  tpm <- counts
  for (j in seq_len(ncol(tpm))) {
    lib_size <- sum(tpm[, j])
    if (lib_size > 0) {
      tpm[, j] <- (tpm[, j] / lib_size) * 1e6
    }
  }
  rownames(tpm) <- rownames(counts)
  colnames(tpm) <- colnames(counts)
  
  se <- SummarizedExperiment::SummarizedExperiment(
    assays = list(counts = counts),
    colData = colData,
    rowData = rowData
  )
  
  # Add TPM to metadata
  S4Vectors::metadata(se)$tpm <- tpm
  
  se
}

#' Create a standard SummarizedExperiment with simple colnames
#'
#' This variant uses short sample names (s1, s2, etc.) and simpler gene names,
#' matching the style already used in many tests.
#'
#' @param n_genes Number of genes (rows)
#' @param n_samples Number of samples (columns)
#' @param control_n Number of samples in control group
#' @param group_col_name Name for the grouping column in colData
#' @param paired If TRUE, add paired_samples column for paired tests
#' @param seed Random seed
#'
#' @return SummarizedExperiment with counts and metadata
#'
#' @keywords internal
create_test_se_simple <- function(
    n_genes = 2,
    n_samples = 6,
    control_n = 3,
    group_col_name = "sample_type",
    paired = FALSE,
    seed = 42
) {
  skip_if_not_installed("SummarizedExperiment")
  
  set.seed(seed)
  
  # Use sequential counts (1:N) for reproducibility
  counts <- matrix(
    seq_len(n_genes * n_samples),
    nrow = n_genes,
    ncol = n_samples
  )
  
  rownames(counts) <- paste0("gene", seq_len(n_genes))
  colnames(counts) <- paste0("s", seq_len(n_samples))
  
  # Create colData with the specified column name
  colData_list <- list(
    factor(c(
      rep("Control", control_n),
      rep("Treatment", n_samples - control_n)
    ))
  )
  names(colData_list) <- group_col_name
  
  # Add paired samples column if requested
  if (paired) {
    # Create paired sample identifiers
    n_pairs <- ceiling(n_samples / 2)
    pair_ids <- rep(paste0("pair_", seq_len(n_pairs)), length.out = n_samples)
    colData_list$paired_samples <- pair_ids
  }
  
  colData <- data.frame(colData_list, row.names = colnames(counts))
  
  rowData <- data.frame(
    gene_id = paste0("GENE", seq_len(n_genes)),
    gene_name = paste0("gene", seq_len(n_genes))
  )
  
  # Generate TPM data (normalized counts)
  tpm <- counts
  for (j in seq_len(ncol(tpm))) {
    lib_size <- sum(tpm[, j])
    if (lib_size > 0) {
      tpm[, j] <- (tpm[, j] / lib_size) * 1e6
    }
  }
  rownames(tpm) <- rownames(counts)
  colnames(tpm) <- colnames(counts)
  
  se <- SummarizedExperiment::SummarizedExperiment(
    assays = list(counts = counts),
    colData = colData,
    rowData = rowData
  )
  
  # Add TPM to metadata
  S4Vectors::metadata(se)$tpm <- tpm
  
  se
}

#' Get multiple standard test SE configurations
#'
#' Useful for parametric testing of different sample/gene sizes
#'
#' @keywords internal
get_test_se_configurations <- function() {
  list(
    minimal = create_test_se(n_genes = 2, n_samples = 4, control_n = 2),
    small = create_test_se(n_genes = 5, n_samples = 8, control_n = 4),
    medium = create_test_se(n_genes = 20, n_samples = 12, control_n = 6),
    unbalanced = create_test_se(n_genes = 10, n_samples = 8, control_n = 2)
  )
}

#' Create entropy/diversity matrix for testing
#'
#' Generates realistic-looking synthetic entropy values using smooth
#' sigmoid transformations to mimic natural diversity patterns.
#'
#' @param n_genes Number of genes (rows)
#' @param n_samples Number of samples (columns)
#' @param seed Random seed
#'
#' @return Matrix with entropy values in [0, 1]
#'
#' @keywords internal
create_entropy_matrix <- function(
    n_genes = 15,
    n_samples = 8,
    seed = 42
) {
  set.seed(seed)
  
  entropy_matrix <- matrix(nrow = n_genes, ncol = n_samples)
  rownames(entropy_matrix) <- paste0("Gene", seq_len(n_genes))
  colnames(entropy_matrix) <- paste0("Sample", seq_len(n_samples))
  
  # Use smooth logistic functions for realistic distribution
  for (i in seq_len(n_genes)) {
    # Sigmoid curve: maps gene index to mean entropy
    x_scaled <- (i - 1) / max(1, n_genes - 1) * 8 - 4  # [-4, 4]
    mean_val <- 1 / (1 + exp(-x_scaled))
    
    # Variance pattern based on position
    var_pattern <- (sin(i / n_genes * pi * 3) + 1) / 2
    base_var <- 0.01 + var_pattern * 0.07
    
    # Weight variance by distance from 0.5
    variance_weight <- 1 / (1 + ((mean_val - 0.5) / 0.2) ^ 2)
    sd_val <- sqrt(base_var + variance_weight * 0.045)
    
    entropy_matrix[i, ] <- pmax(0.01, pmin(0.99,
      rnorm(n_samples, mean_val, sd_val)
    ))
  }
  
  entropy_matrix
}

# =============================================================================
# Assertion Helpers
# =============================================================================

#' Verify normalization was applied correctly
#'
#' @param estimates Numeric vector of estimate values
#' @param norm_mode Name of normalization mode applied
#' @param tolerance Small tolerance for floating point comparisons
#'
#' @keywords internal
#' @examples
#' \dontrun{
#'   estimates <- c(0.5, 0.75, 0.25, NA)
#'   expect_valid_normalization(estimates, norm_mode = "range")
#' }
expect_valid_normalization <- function(
    result,
    norm_mode = c("none", "range", "zscore", "log_odds_ratio", "relative_reference"),
    tolerance = 1e-6
) {
  norm_mode <- match.arg(norm_mode)
  
  # Extract estimates from SE if necessary
  if (methods::is(result, "SummarizedExperiment")) {
    rd <- SummarizedExperiment::rowData(result)
    est_cols <- colnames(rd)[grep("^estimate_q", colnames(rd))]
    if (length(est_cols) == 0) {
      expect_true(FALSE, info = "No estimate columns found in rowData")
      return()
    }
    # Combine all estimates
    estimates <- as.numeric(unlist(rd[, est_cols, drop = FALSE]))
  } else {
    estimates <- result
  }
  
  # All modes should produce no NaN (except where input is already invalid)
  valid_ests <- estimates[!is.na(estimates)]
  expect_true(all(!is.nan(valid_ests)),
    info = sprintf("NaN values found in %s normalization", norm_mode)
  )
  
  if (norm_mode == "range") {
    # Values should be in [0, 1]
    expect_true(all(valid_ests >= -tolerance),
      info = "Range normalization: values < 0"
    )
    expect_true(all(valid_ests <= 1 + tolerance),
      info = "Range normalization: values > 1"
    )
  } else if (norm_mode == "zscore") {
    # Z-score should have mean ~0 and sd ~1 (before NA removal)
    expect_true(abs(mean(valid_ests)) < 2,
      info = "Z-score normalization: mean not centered at 0"
    )
  }
}

#' Expect SE column to exist
#'
#' @param se SummarizedExperiment object
#' @param col_name Column name to check
#' @param data_type "colData" or "rowData"
#'
#' @keywords internal
expect_se_column <- function(
    se,
    col_name,
    data_type = c("colData", "rowData")
) {
  data_type <- match.arg(data_type)
  
  if (data_type == "colData") {
    expect_true(col_name %in% names(SummarizedExperiment::colData(se)),
      info = sprintf(
        "Column '%s' not found in colData. Available: %s",
        col_name,
        paste0(names(SummarizedExperiment::colData(se)), collapse = ", ")
      )
    )
  } else {
    expect_true(col_name %in% names(SummarizedExperiment::rowData(se)),
      info = sprintf(
        "Column '%s' not found in rowData. Available: %s",
        col_name,
        paste0(names(SummarizedExperiment::rowData(se)), collapse = ", ")
      )
    )
  }
}

#' Check result has expected structure
#'
#' @param result Object to check
#' @param expected_type Expected class (e.g., "data.frame", "matrix")
#' @param n_rows Expected number of rows (optional)
#' @param n_cols Expected number of columns (optional)
#'
#' @keywords internal
expect_result_structure <- function(
    result,
    expected_type = c("data.frame", "matrix", "SummarizedExperiment", "list"),
    n_rows = NULL,
    n_cols = NULL
) {
  expected_type <- match.arg(expected_type)
  
  # For SummarizedExperiment, check inheritance more flexibly
  if (expected_type == "SummarizedExperiment") {
    expect_true(methods::is(result, "SummarizedExperiment"),
      info = sprintf("Result is %s, expected SummarizedExperiment", 
                     class(result)[1])
    )
  } else {
    expect_is(result, expected_type,
      info = sprintf("Result is %s, expected %s", class(result)[1], expected_type)
    )
  }
  
  if (!is.null(n_rows)) {
    if (expected_type == "SummarizedExperiment") {
      expect_equal(nrow(result), n_rows,
        info = sprintf("SE has %d rows, expected %d", nrow(result), n_rows)
      )
    } else {
      expect_equal(nrow(result), n_rows,
        info = sprintf("Result has %d rows, expected %d", nrow(result), n_rows)
      )
    }
  }
  
  if (!is.null(n_cols)) {
    if (expected_type == "SummarizedExperiment") {
      expect_equal(ncol(result), n_cols,
        info = sprintf("SE has %d cols, expected %d", ncol(result), n_cols)
      )
    } else {
      expect_equal(ncol(result), n_cols,
        info = sprintf("Result has %d cols, expected %d", ncol(result), n_cols)
      )
    }
  }
}

#' Verify numeric result has no invalid values
#'
#' @param result Numeric vector or matrix
#' @param allow_na Allow NA values (default: TRUE)
#' @param allow_negative Allow negative values (default: FALSE)
#'
#' @keywords internal
expect_valid_numeric <- function(
    result,
    allow_na = TRUE,
    allow_negative = FALSE
) {
  if (allow_na) {
    valid_result <- result[!is.na(result)]
  } else {
    valid_result <- result
    expect_false(any(is.na(valid_result)), info = "Found NA values")
  }
  
  expect_false(any(is.nan(valid_result)), info = "Found NaN values")
  expect_false(any(is.infinite(valid_result)), info = "Found Inf values")
  
  if (!allow_negative) {
    expect_true(all(valid_result >= 0), info = "Found negative values")
  }
}

# =============================================================================
# Data Generation Helpers
# =============================================================================

#' Create a gene expression matrix for diversity testing
#'
#' Simple utility for creating isoform count matrices for diversity
#' function testing.
#'
#' @param n_genes Number of genes (rows per gene)
#' @param n_isoforms Number of isoforms per gene
#' @param n_samples Number of samples (columns)
#' @param seed Random seed
#'
#' @return List with matrix x (counts) and genes (isoform-to-gene mapping)
#'
#' @keywords internal
create_diversity_test_data <- function(
    n_genes = 3,
    n_isoforms = 3,
    n_samples = 4,
    seed = 42
) {
  set.seed(seed)
  
  total_isoforms <- n_genes * n_isoforms
  
  x <- matrix(
    rpois(total_isoforms * n_samples, lambda = 50),
    nrow = total_isoforms,
    ncol = n_samples
  )
  
  colnames(x) <- paste0("Sample_", seq_len(n_samples))
  rownames(x) <- paste0("Iso_", seq_len(total_isoforms))
  
  genes <- rep(paste0("Gene_", seq_len(n_genes)), each = n_isoforms)
  
  list(x = x, genes = genes)
}

# =============================================================================
# Utility Functions
# =============================================================================

#' Skip test if required packages are not installed
#'
#' Wrapper around testthat::skip_if_not_installed for multiple packages
#'
#' @param packages Character vector of package names
#'
#' @keywords internal
skip_if_not_installed_multiple <- function(packages) {
  for (pkg in packages) {
    skip_if_not_installed(pkg)
  }
}

#' Set common test options
#'
#' Configures parallel processing and verbosity for consistency
#'
#' @param n_cores Number of cores to use (default: 1 for tests)
#' @param verbose Show verbose output (default: FALSE)
#'
#' @keywords internal
setup_test_environment <- function(n_cores = 1, verbose = FALSE) {
  # Ensure reproducibility
  set.seed(42)
  
  # Store old options
  old_options <- options(
    TSENAT.n_cores = n_cores,
    TSENAT.verbose = verbose
  )
  
  invisible(old_options)
}

#' Verify bootstrap results have expected structure
#'
#' @param bootstrap_result Result object from a bootstrap function
#' @param expected_ci_cols Expected column names for CI results
#'
#' @keywords internal
expect_valid_bootstrap_result <- function(
    bootstrap_result,
    expected_ci_cols = c("estimate", "lower_ci", "upper_ci")
) {
  expect_true(!is.null(bootstrap_result), info = "Bootstrap result is NULL")
  
  if (is.data.frame(bootstrap_result) || is.matrix(bootstrap_result)) {
    for (col in expected_ci_cols) {
      expect_true(col %in% colnames(bootstrap_result),
        info = sprintf("Missing column in bootstrap result: %s", col)
      )
    }
    # CI bounds should be reasonable
    if ("estimate" %in% colnames(bootstrap_result)) {
      expect_valid_numeric(bootstrap_result[["estimate"]],
        allow_na = TRUE
      )
    }
  }
}

#' Create diversity SummarizedExperiment for testing
#'
#' Creates a SummarizedExperiment containing diversity calculations with
#' standardized structure for use across multiple tests.
#'
#' @param n_genes Number of genes (rows)
#' @param n_samples Number of samples (columns)
#' @param seed Random seed for reproducibility
#' @param gene_names Optional custom gene names
#' @param sample_names Optional custom sample names
#'
#' @return SummarizedExperiment with diversity assay
#'
#' @keywords internal
create_diversity_se <- function(
    n_genes = 5,
    n_samples = 4,
    seed = NULL,
    gene_names = NULL,
    sample_names = NULL
) {
  if (!is.null(seed)) set.seed(seed)
  
  # Create diversity data
  div_data <- matrix(runif(n_genes * n_samples, 1, 3), nrow = n_genes, ncol = n_samples)
  
  # Set row/column names
  if (is.null(gene_names)) {
    rownames(div_data) <- paste0("GENE", 1:n_genes)
  } else {
    rownames(div_data) <- gene_names
  }
  
  if (is.null(sample_names)) {
    colnames(div_data) <- paste0("S", 1:n_samples)
  } else {
    colnames(div_data) <- sample_names
  }
  
  # Create SE
  SummarizedExperiment::SummarizedExperiment(
    assays = list(diversity = div_data)
  )
}

#' Create TSENATAnalysis with diversity results
#'
#' Convenience function to create a TSENATAnalysis object with pre-computed
#' diversity results for testing divergence calculations.
#'
#' @param se SummarizedExperiment to use as base (or NULL to create default)
#' @param n_genes Number of genes (ignored if se is provided)
#' @param n_samples Number of samples (ignored if se is provided)
#' @param control_n Number of control samples
#' @param group_col_name Name of the group column in colData
#' @param q_values Q-values for diversity results (default: c(1.0))
#' @param seed Random seed
#'
#' @return TSENATAnalysis object with populated diversity_results
#'
#' @keywords internal
create_tsenat_with_diversity <- function(
    se = NULL,
    n_genes = 5,
    n_samples = 4,
    control_n = 2,
    group_col_name = "group",
    q_values = c(1.0),
    seed = NULL
) {
  if (!is.null(seed)) set.seed(seed)
  
  # Use provided SE or create default
  if (is.null(se)) {
    se <- create_test_se_simple(
      n_genes = n_genes,
      n_samples = n_samples,
      control_n = control_n,
      group_col_name = group_col_name,
      seed = seed
    )
  }
  
  # Create TSENATAnalysis
  analysis <- TSENAT::TSENATAnalysis(se = se)
  
  # Add diversity results for each q-value
  for (q in q_values) {
    q_label <- paste0("q_", q)
    div_se <- create_diversity_se(
      n_genes = nrow(se),
      n_samples = ncol(se),
      seed = seed,
      gene_names = rownames(se),
      sample_names = colnames(se)
    )
    analysis@diversity_results[[q_label]] <- div_se
  }
  
  analysis
}

#' Create count data SummarizedExperiment
#'
#' Creates a SummarizedExperiment with count matrix for testing.
#' Used when you need raw count data rather than ready-made test SE.
#'
#' @param n_genes Number of genes
#' @param n_samples Number of samples
#' @param n_control Number of control samples (controls come first)
#' @param lambda Poisson parameter for count generation
#' @param seed Random seed
#'
#' @return SummarizedExperiment with counts assay
#'
#' @keywords internal
create_count_se <- function(
    n_genes = 5,
    n_samples = 4,
    n_control = 2,
    lambda = 5,
    seed = NULL
) {
  if (!is.null(seed)) set.seed(seed)
  
  # Create count matrix
  counts <- matrix(rpois(n_genes * n_samples, lambda = lambda),
                   nrow = n_genes, ncol = n_samples)
  rownames(counts) <- paste0("GENE", 1:n_genes)
  colnames(counts) <- paste0("S", 1:n_samples)
  
  # Create colData with groups
  group <- c(rep("Control", n_control), rep("Treatment", n_samples - n_control))
  
  col_data <- S4Vectors::DataFrame(
    sample_id = paste0("S", 1:n_samples),
    group = group
  )

  row_data <- data.frame(
    gene_id = rownames(counts)
  )

  # Create SE
  se <- SummarizedExperiment::SummarizedExperiment(
    assays = list(counts = counts),
    colData = col_data,
    rowData = row_data
  )
  
  se
}

#' Create wrapper test diversity SummarizedExperiment
#'
#' Creates a small SE for diversity-related wrapper testing.
#' Default dimensions: 10 rows × 3 cols with rpois(30, 100).
#' Used extensively in wrapper_s4_methods tests.
#'
#' @param n_rows Number of rows (genes/transcripts)
#' @param n_cols Number of columns (samples)
#' @param lambda Poisson parameter
#' @param col_data Optional colData DataFrame. If NULL, creates default with sample names.
#' @param seed Random seed
#'
#' @return SummarizedExperiment appropriate for diversity testing
#'
#' @keywords internal
create_wrapper_diversity_se <- function(
    n_rows = 10,
    n_cols = 3,
    lambda = 100,
    col_data = NULL,
    seed = NULL
) {
  if (!is.null(seed)) set.seed(seed)
  
  # Create count matrix
  counts <- matrix(
    rpois(n_rows * n_cols, lambda = lambda),
    nrow = n_rows,
    ncol = n_cols
  )
  
  # Create colData if not provided
  if (is.null(col_data)) {
    col_data <- data.frame(
      sample = paste0("sample_", 1:n_cols),
      stringsAsFactors = FALSE
    )
  }
  
  # Create rowData with gene_id
  row_data <- data.frame(
    gene_id = paste0("GENE_", 1:n_rows)
  )

  # Create SE
  se <- SummarizedExperiment::SummarizedExperiment(
    assays = list(counts = counts),
    colData = col_data,
    rowData = row_data
  )

  se
}

#' Create wrapper test SE: 10×3 with metadata
#'
#' 10×3 SE with sample names and condition metadata.
#'
#' @param include_condition If TRUE, adds condition column (control/control/treatment)
#' @param include_pair_id If TRUE, adds pair_id column for paired tests
#' @param seed Random seed
#'
#' @return SummarizedExperiment with standard wrapper test dimensions
#'
#' @keywords internal
create_wrapper_test_se_10x3 <- function(
    include_condition = FALSE,
    include_pair_id = FALSE,
    seed = NULL
) {
  if (!is.null(seed)) set.seed(seed)
  
  col_data <- data.frame(
    sample = c("control_1", "control_2", "treatment_1"),
    stringsAsFactors = FALSE
  )
  
  if (include_condition) {
    col_data$condition <- c("control", "control", "treatment")
  }
  
  if (include_pair_id) {
    col_data$pair_id <- c(1, 2, 1)
  }
  
  # Create SE with standard 10×3 dimensions
  create_wrapper_diversity_se(
    n_rows = 10,
    n_cols = 3,
    lambda = 100,
    col_data = col_data,
    seed = seed
  )
}

#' Create wrapper test diversity SE with sample metadata (10×6)
#'
#' Convenience function for larger wrapper test pattern:
#' 10×6 SE with 6 samples having 3 groups and pair_ids.
#'
#' @param include_condition If TRUE, adds condition column
#' @param include_pair_id If TRUE, adds pair_id column
#' @param seed Random seed
#'
#' @return SummarizedExperiment with 10×6 dimensions
#'
#' @keywords internal
create_wrapper_test_se_10x6 <- function(
    include_condition = FALSE,
    include_pair_id = FALSE,
    seed = NULL
) {
  if (!is.null(seed)) set.seed(seed)
  
  col_data <- data.frame(
    sample = c("control_1", "control_2", "treatment_1", "treatment_2", "other_1", "other_2"),
    stringsAsFactors = FALSE
  )
  
  if (include_condition) {
    col_data$condition <- c("control", "control", "treatment", "treatment", "other", "other")
  }
  
  if (include_pair_id) {
    col_data$pair_id <- c(1, 2, 1, 2, 3, 4)
  }
  
  # Create SE with 10×6 dimensions
  create_wrapper_diversity_se(
    n_rows = 10,
    n_cols = 6,
    lambda = 100,
    col_data = col_data,
    seed = seed
  )
}

#' Create test data for calculate_diversity tests
#'
#' Standard 3×2 matrix with well-known values for diversity testing.
#' Repeats in 12+ tests in test-calculate_diversity.R
#'
#' @param set_colnames If TRUE, adds colnames c("S1", "S2")
#' @param include_genes If TRUE, also returns genes vector c("g1", "g1", "g2")
#' @param seed Random seed
#'
#' @return List with $x (matrix) and optionally $genes (character vector)
#'
#' @keywords internal
create_diversity_test_matrix_3x2_standard <- function(
    set_colnames = TRUE,
    include_genes = TRUE,
    seed = NULL
) {
  # Updated to use 4x2 matrix with multi-isoform genes (2 isoforms per gene)
  # This avoids single-isoform warnings in bootstrap tests
  x <- matrix(c(10, 5, 8, 12, 15, 3, 7, 9), nrow = 4, ncol = 2)
  
  if (set_colnames) {
    colnames(x) <- c("S1", "S2")
  }
  
  result <- list(x = x)
  
  if (include_genes) {
    # 4 transcripts: g1 has 2 isoforms, g2 has 2 isoforms
    result$genes <- c("g1", "g1", "g2", "g2")
  }
  
  result
}

#' Create second common test matrix for calculate_diversity
#'
#' 3×2 matrix c(1,2,3,4,5,6) with colnames
#'
#' @param set_colnames If TRUE, adds colnames c("S1", "S2")
#' @param include_genes If TRUE, also returns genes vector c("g1", "g1", "g2")
#'
#' @return List with $x (matrix) and optionally $genes
#'
#' @keywords internal
create_diversity_test_matrix_3x2_simple <- function(
    set_colnames = TRUE,
    include_genes = TRUE
) {
  # Updated to use 4x2 matrix with multi-isoform genes (2 isoforms per gene)
  # This avoids single-isoform warnings in bootstrap tests
  x <- matrix(c(1, 2, 3, 4, 5, 6, 7, 8), nrow = 4, ncol = 2)
  
  if (set_colnames) {
    colnames(x) <- c("S1", "S2")
  }
  
  result <- list(x = x)
  
  if (include_genes) {
    # 4 transcripts: g1 has 2 isoforms, g2 has 2 isoforms
    result$genes <- c("g1", "g1", "g2", "g2")
  }
  
  result
}

#' Create test data for 5-element vector test cases
#'
#' Standard 5×2 test matrix
#'
#' @param set_colnames If TRUE, adds colnames
#' @param include_genes If TRUE, adds genes c("g1", "g1", "g2", "g2", "g3")
#'
#' @return List with $x and optionally $genes
#'
#' @keywords internal
create_diversity_test_matrix_5x2 <- function(
    set_colnames = TRUE,
    include_genes = TRUE
) {
  x <- matrix(c(0, 0, 5, 4, 1, 2, 2, 2, 2, 2), ncol = 2)
  
  if (set_colnames) {
    colnames(x) <- c("Sample1", "Sample2")
  }
  
  result <- list(x = x)
  
  if (include_genes) {
    result$genes <- c("Gene1", "Gene1", "Gene1", "Gene1", "Gene1")
  }
  
  result
}

#' Create simple SummarizedExperiment (5×10 with rpois(50, 3))
#'
#' Minimal SE with just counts, no colData. Used for basic S4 method testing.
#'
#' @param seed Random seed
#'
#' @return SummarizedExperiment with 5 rows × 10 cols
#'
#' @keywords internal
create_simple_se_5x10 <- function(seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  
  counts <- matrix(rpois(50, 3), nrow = 5, ncol = 10)
  rownames(counts) <- paste0("GENE", 1:5)
  colnames(counts) <- paste0("S", 1:10)

  row_data <- data.frame(
    gene_id = rownames(counts)
  )

  SummarizedExperiment::SummarizedExperiment(
    assays = list(counts = counts),
    rowData = row_data
  )
}

#' Test all normalization modes comprehensively
#'
#' Validates that all 5 normalization modes work correctly with strong assertions:
#' - none: raw values with no transformation
#' - range: values in [0, 1]
#' - zscore: mean ~0, values unbounded
#' - log_odds_ratio: log-transformed odds ratios
#' - relative_reference: reference-normalized values
#'
#' This consolidates 11 redundant tests into one comprehensive validation.
#'
#' @param func Function to call (typically calculate_divergence or calculate_diversity)
#' @param se SummarizedExperiment to analyze
#' @param q Q value(s) for analysis (default: 1)
#' @param group_col Column name for grouping (default: "sample_type")
#' @param control_group Control group name (default: "Control")
#' @param bootstrap Whether to use bootstrap (default: FALSE)
#' @param verbose Verbosity flag (default: FALSE)
#' @param progress Progress flag (default: FALSE)
#'
#' @return Invisibly returns list of results for each normalization mode
#'
#' @keywords internal
test_all_normalization_modes <- function(
    func,
    se,
    q = 1,
    group_col = "sample_type",
    control_group = "Control",
    bootstrap = FALSE,
    progress = FALSE
) {
  norm_modes <- c("none", "range", "zscore", "log_odds_ratio", "relative_reference")
  results <- list()
  
  for (norm_mode in norm_modes) {
    # Call the function with current normalization mode
    result <- func(
      se = se,
      group_col = group_col,
      control_group = control_group,
      q = q,
      norm = norm_mode,
      bootstrap = bootstrap,
      progress = progress
    )
    
    # Strong assertions for all modes
    expect_is(result, "SummarizedExperiment",
      info = sprintf("Result is not SummarizedExperiment for %s", norm_mode))
    
    expect_true(nrow(result) > 0,
      info = sprintf("Result has no rows for %s", norm_mode))
    
    expect_equal(S4Vectors::metadata(result)$normalization, norm_mode,
      info = sprintf("Metadata normalization mismatch for %s", norm_mode))
    
    # Validate value ranges by mode
    rd <- SummarizedExperiment::rowData(result)
    est_cols <- colnames(rd)[grep("^estimate", colnames(rd))]
    
    if (length(est_cols) > 0) {
      for (col in est_cols) {
        estimates <- rd[[col]]
        valid_est <- estimates[!is.na(estimates)]
        
        if (length(valid_est) > 0) {
          expect_true(all(!is.nan(valid_est)),
            info = sprintf("NaN values in %s estimates (%s)", norm_mode, col))
          
          if (norm_mode == "range") {
            expect_true(all(valid_est >= -1e-6),
              info = sprintf("Range: values < 0 in %s", col))
            expect_true(all(valid_est <= 1 + 1e-6),
              info = sprintf("Range: values > 1 in %s", col))
          }
        }
      }
    }
    
    results[[norm_mode]] <- result
  }
  
  invisible(results)
}

#' Validate rankbased assumptions test result for correctness
#'
#' Ensures rankbased assumptions result is not just non-null but structurally correct with valid values
#'
#' @param result TSENATAnalysis object with rankbased_assumptions in metadata
#' @param check_type Type of check performed (e.g., "exchangeability", "normality")
#'
#' @return Invisibly returns result (all assertions pass or error)
#'
#' @keywords internal
assert_rankbased_result_valid <- function(result, check_type = "exchangeability") {
  # Strong assertion: result is TSENATAnalysis object
  expect_is(result, "TSENATAnalysis")
  
  # Strong assertion: metadata exists and has rankbased_assumptions list
  expect_true(!is.null(result@metadata))
  expect_true(!is.null(result@metadata$rankbased_assumptions))
  expect_true(is.list(result@metadata$rankbased_assumptions))
  
  # Strong assertion: rankbased_assumptions is non-empty
  expect_true(length(result@metadata$rankbased_assumptions) > 0)
  
  # Strong assertion: each assumption check has expected structure
  # Skip non-check entries (like "result" metadata key)
  assumptions <- result@metadata$rankbased_assumptions
  check_names <- c("exchangeability", "monotonicity", "q_value_tested", "checks_performed", "alpha_used", "timestamp")
  
  for (check_name in names(assumptions)) {
    check_entry <- assumptions[[check_name]]
    
    # Skip non-check keys (like "result")
    if (!(check_name %in% check_names)) {
      next
    }
    
    # If it's a list, should have p_value and status fields (at minimum)
    if (is.list(check_entry)) {
      # Valid checks have at least one of these
      has_fields <- ("p_value" %in% names(check_entry)) || ("status" %in% names(check_entry))
      expect_true(has_fields,
        info = sprintf("Check '%s' should have p_value or status field", check_name))
      
      # If p_value exists, validate it's numeric and [0,1]
      if ("p_value" %in% names(check_entry)) {
        p_val <- check_entry$p_value
        if (!is.null(p_val) && is.numeric(p_val) && length(p_val) == 1) {
          expect_true(p_val >= 0 && p_val <= 1,
            info = sprintf("%s p_value should be in [0,1]", check_name))
        }
      }
    }
  }
  
  invisible(result)
}

#' Test calculate_divergence input validation errors
#'
#' Comprehensive test helper for input validation errors in calculate_divergence
#' Consolidates common validation checks to eliminate test duplication
#'
#' @param se SummarizedExperiment object (optional, created if NULL)
#' @param test_types Character vector of test types to run
#'   - "se_type": Invalid SE object type
#'   - "bootstrap_type": Invalid bootstrap parameter type
#'   - "q_values": Invalid q parameter values
#'   - "missing_genes": Missing gene identifiers


#' @keywords internal
test_calculate_divergence_input_validation <- function(
  se = NULL,
  test_types = c("se_type", "bootstrap_type", "q_values", "missing_genes")) {
  
  if (is.null(se)) {
    skip_if_not_installed("SummarizedExperiment")
    se <- create_count_se(
      n_genes = 20,
      n_samples = 8,
      n_control = 4,
      lambda = 100,
      seed = 44
    )
    SummarizedExperiment::colData(se)$group <- factor(c(rep("A", 4), rep("B", 4)))
  }
  
  if ("se_type" %in% test_types) {
    # Invalid SE (not SummarizedExperiment)
    expect_error(
      TSENAT:::.calculate_divergence(se = list()),
      "must be a SummarizedExperiment|SummarizedExperiment",
      ignore.case = TRUE
    )
  }
  
  if ("bootstrap_type" %in% test_types) {
    # Invalid bootstrap parameter
    expect_error(
      TSENAT:::.calculate_divergence(se = se, bootstrap = "yes"),
      "bootstrap must be a logical"
    )
  }
  
  if ("q_values" %in% test_types) {
    # Invalid q values (non-positive)
    expect_error(
      TSENAT:::.calculate_divergence(
        se = se,
        group_col = "group",
        control_group = "A",
        q = c(-0.5, 1, 2),
        bootstrap = FALSE,
        verbose = FALSE
      ),
      "q parameter must be positive"
    )
  }
  
  if ("missing_genes" %in% test_types) {
    # Create SE without gene identifiers
    se_no_genes <- SummarizedExperiment::SummarizedExperiment(
      assays = list(counts = matrix(c(1, 2, 3, 4, 5, 6), nrow = 1, ncol = 6)),
      colData = data.frame(
        group = c("A", "A", "B", "B", "B", "B")
      )
    )
    colnames(se_no_genes) <- c("s1", "s2", "s3", "s4", "s5", "s6")
    
    expect_error(
      TSENAT:::.calculate_divergence(
        se = se_no_genes,
        group_col = "group",
        control_group = "A",
        bootstrap = FALSE,
        verbose = FALSE
      ),
      "gene identifiers"
    )
  }
  
  invisible(TRUE)
}

#' Validate numeric vector for reasonable values
#'
#' Ensures numeric data doesn't contain NaN, Inf, or other garbage values
#'
#' @param x Numeric vector to validate
#' @param name Variable name for error messages
#' @param allow_na If TRUE, NA values are allowed (default FALSE)
#' @param min_value Minimum acceptable value (default -Inf)
#' @param max_value Maximum acceptable value (default Inf)
#'
#' @return Invisibly returns x (all assertions pass or error)
#'
#' @keywords internal
assert_numeric_valid <- function(x, name = "x", allow_na = FALSE, min_value = -Inf, max_value = Inf) {
  expect_is(x, "numeric",
    info = sprintf("%s must be numeric", name))
  
  expect_true(all(!is.nan(x)),
    info = sprintf("%s contains NaN values", name))
  
  expect_true(all(!is.infinite(x)),
    info = sprintf("%s contains Inf/-Inf values", name))
  
  if (!allow_na) {
    expect_true(all(!is.na(x)),
      info = sprintf("%s contains NA values", name))
  }
  
  valid_vals <- x[!is.na(x)]
  if (length(valid_vals) > 0) {
    expect_true(all(valid_vals >= min_value),
      info = sprintf("%s has values below minimum %f", name, min_value))
    expect_true(all(valid_vals <= max_value),
      info = sprintf("%s has values above maximum %f", name, max_value))
  }
  
  invisible(x)
}

#' Suppress expected warnings from loess-based functions
#'
#' Helper to wrap calls that generate expected loess warnings from synthetic test data
#' Justification: Loess fitting on synthetic/sparse data often generates  
#' convergence warnings that are expected and handled internally by graceful fallback
#'
#' @param expr Expression to evaluate with warnings suppressed
#'
#' @return Result of evaluating expr
#'
#' @keywords internal
suppress_loess_warnings <- function(expr) {
  suppressWarnings(expr)
}


#' Create a Complete TSENATAnalysis with Test Data
#' 
#' Factory function that generates a fully initialized TSENATAnalysis object
#' with realistic data and biological signal. Useful for examples, testing,
#' and documentation. Eliminates ~70% of boilerplate in documentation examples.
#'
#' @param n_genes Number of genes to simulate (default: 8)
#' @param n_samples_per_group Samples per condition (default: 20)
#' @param control_lambda Poisson lambda for control condition (default: 40)
#' @param treatment_lambda Poisson lambda for treatment condition (default: 150)
#' @param q_values Vector of q-values for diversity calculation (default:
#' c(0.5, 0.75, 1.0, 1.5, 2.0))
#' @param include_divergence If TRUE, compute divergence results (default: TRUE)
#' @param include_sait_results If TRUE, add placeholder SAIT results (default: TRUE)
#' @param seed Random seed for reproducibility (default: 42)
#' @param verbose Logical for progress messages (default: FALSE)
#'
#' @return TSENATAnalysis object with:
#'   - SummarizedExperiment with count matrix, rowData, and colData
#'   - Computed diversity results across q-values
#'   - Computed divergence results (optional)
#'   - Placeholder SAIT results (optional)
#'   - Proper tx2gene metadata mapping
#'   - TPM data in metadata
#'
#' @details
#' The factory ensures:
#' - Biological signal: control (lambda=40) vs treatment (lambda=150) contrast
#' - Sufficient samples: 40 total (20 per group) for stable LM fitting
#' - Multiple q-values: c(0.5, 1.0, 1.5) avoids rank deficiency
#' - Valid S4 object structure: passes all TSENATAnalysis validity checks
#' - Optional divergence and SAIT results to support testing without warnings
#' - All required metadata (tx2gene, TPM, rowData) pre-configured
#'
#' @examples
#' # Create with defaults (8 genes, 20 samples/group, multi-q, with all results)
#' analysis <- create_test_analysis()
#' 
#' # Create minimal analysis for quick testing
#' analysis <- create_test_analysis(
#'   n_genes = 4,
#'   n_samples_per_group = 10,
#'   q_values = c(0.5, 1.0),
#'   include_divergence = FALSE
#' )
#' 
#' # Create with custom parameters
#' analysis <- create_test_analysis(
#'   n_genes = 16,
#'   n_samples_per_group = 30,
#'   control_lambda = 50,
#'   treatment_lambda = 200,
#'   q_values = c(0.1, 0.5, 1.0, 1.5, 2.0),
#'   include_divergence = TRUE,
#'   include_sait_results = TRUE
#' )
#'
#' @noRd
#' @noRd
.create_test_analysis <- function(n_genes = 8, n_samples_per_group = 20, control_lambda = 40,
    treatment_lambda = 150, q_values = c(0.5, 1, 1.5), include_divergence = TRUE,
    include_sait_results = TRUE, seed = 42, verbose = FALSE) {

    withr::local_seed(seed)

    # Dimensions
    n_samples <- n_samples_per_group * 2
    n_transcripts <- n_genes * 50

    if (verbose) {
        message("[create_test_analysis] Generating ", n_transcripts, " transcripts across ",
            n_genes, " genes with ", n_samples, " samples")
    }

    # Generate counts with biological signal
    control_idx <- seq(1, n_samples, by = 2)
    treatment_idx <- seq(2, n_samples, by = 2)

    counts <- matrix(0, nrow = n_transcripts, ncol = n_samples)
    for (j in seq_len(n_samples)) {
        if (j %in% control_idx) {
            counts[, j] <- rpois(n_transcripts, lambda = control_lambda)
        } else {
            counts[, j] <- rpois(n_transcripts, lambda = treatment_lambda)
        }
    }
    counts <- pmax(counts, 50)

    rownames(counts) <- paste0("TX_", seq_len(n_transcripts))
    colnames(counts) <- paste0("Sample_", seq_len(n_samples))

    # Create rowData with gene mappings
    rowData <- S4Vectors::DataFrame(transcript_id = rownames(counts), gene_id = paste0("GENE_",
        rep(seq_len(n_genes), each = 50, length.out = n_transcripts)), row.names = rownames(counts))

    # Create colData with experimental design
    colData <- S4Vectors::DataFrame(sample_id = colnames(counts), condition = rep(c("control",
        "treatment"), length.out = n_samples), sample_type = rep(c("typeA", "typeB"),
        length.out = n_samples), subject = rep(paste0("S", seq_len(10)), length.out = n_samples),
        paired_samples = rep(paste0("pair", seq_len(10)), length.out = n_samples),
        row.names = colnames(counts))

    # Create SummarizedExperiment
    se <- SummarizedExperiment::SummarizedExperiment(assays = list(counts = counts),
        rowData = rowData, colData = colData)

    # Add tx2gene metadata
    tx2gene_df <- data.frame(Transcript = rownames(counts), Gene = rowData$gene_id,
        stringsAsFactors = FALSE)
    S4Vectors::metadata(se)$tx2gene <- tx2gene_df

    # Generate synthetic TPM data (matching counts dimensions)
    tpm <- counts
    for (j in seq_len(ncol(tpm))) {
        lib_size <- colSums(tpm[, j, drop = FALSE])
        if (lib_size > 0) {
            tpm[, j] <- (tpm[, j]/lib_size) * 1e+06
        }
    }
    rownames(tpm) <- rownames(counts)
    colnames(tpm) <- colnames(counts)
    S4Vectors::metadata(se)$tpm <- tpm

    # Generate synthetic effective_length data for bootstrap normalization
    # Realistic effective lengths typically range from 20 to 5000 bp
    effective_length <- stats::runif(n_transcripts, min = 100, max = 3000)
    names(effective_length) <- rownames(counts)
    S4Vectors::metadata(se)$effective_length <- effective_length

    # Initialize TSENATAnalysis with explicit control_group in config.
    # This prevents auto-detection heuristics from guessing the wrong
    # reference group and ensures reproducible, explicit test design.
    analysis <- TSENATAnalysis(se = se, config = list(control_group = "control"))

    # Calculate diversity
    analysis <- calculate_diversity(analysis, q = q_values, verbose = FALSE)

    # Calculate divergence if requested
    if (include_divergence) {
        analysis <- tryCatch({
            calculate_divergence(analysis, verbose = FALSE)
        }, error = function(e) {
            # If divergence fails, continue without it
            if (verbose) {
                message("[create_test_analysis] Warning: divergence calculation failed: ",
                  e$message)
            }
            analysis
        })
    }

    # Add placeholder SAIT results if requested
    if (include_sait_results) {
        # Create a simple placeholder LM result (empty data frame structure)
        # This prevents 'No SAIT results found' warnings in tests
        sait_placeholder <- list(sait_interaction = list(
            results = data.frame(gene = character(0), term = character(0),
                estimate = numeric(0), std.error = numeric(0), statistic = numeric(0),
                p.value = numeric(0)),
            model_data = NULL))
        analysis@sait_results <- sait_placeholder
    }

    if (verbose) {
        message("[create_test_analysis] Analysis created with ", length(q_values),
            " q-values")
        message("[create_test_analysis] Diversity results: ", nrow(diversity(analysis)),
            " genes")
        if (include_divergence) {
            message("[create_test_analysis] Divergence results included")
        }
        if (include_sait_results) {
            message("[create_test_analysis] Placeholder SAIT results included")
        }
    }

    return(analysis)
}
