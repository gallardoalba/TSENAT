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
    gene_name = rownames(counts)
  )
  
  SummarizedExperiment::SummarizedExperiment(
    assays = list(counts = counts),
    colData = colData,
    rowData = rowData
  )
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
    gene_name = paste0("gene", seq_len(n_genes))
  )
  
  SummarizedExperiment::SummarizedExperiment(
    assays = list(counts = counts),
    colData = colData,
    rowData = rowData
  )
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
  
  # Create SE
  se <- SummarizedExperiment::SummarizedExperiment(
    assays = list(counts = counts),
    colData = col_data
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
  
  # Create SE
  se <- SummarizedExperiment::SummarizedExperiment(
    assays = list(counts = counts),
    colData = col_data
  )
  
  se
}

#' Create wrapper test diversity SE with sample metadata
#'
#' Convenience function for the most common wrapper test pattern:
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
  
  SummarizedExperiment::SummarizedExperiment(
    assays = list(counts = counts)
  )
}
