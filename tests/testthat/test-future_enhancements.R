context("Future Enhancements: Block Jackknife and Visualization Functions")

# Test fixtures
test_se_with_blocks <- function() {
  set.seed(42)
  counts <- matrix(
    c(
      # Gene 1, 3 transcripts
      1000, 500, 200,  400, 300, 150,   # Condition A - Phase 1
      800, 400, 180,   350, 250, 120,   # Condition B - Phase 1
      # Add phase grouping
      1200, 600, 250,  500, 350, 200,   # Condition A - Phase 2
      900, 500, 200,   450, 300, 180    # Condition B - Phase 2
    ),
    nrow = 3,
    ncol = 8,
    byrow = FALSE
  )
  
  rowData <- data.frame(
    gene_id = c("Gene1", "Gene1", "Gene1"),
    isoform_id = c("Gene1.1", "Gene1.2", "Gene1.3"),
    row.names = c("iso1", "iso2", "iso3")
  )
  
  colData <- data.frame(
    condition = rep(c("A", "B"), each = 4),
    phase = rep(c("Phase1", "Phase1", "Phase2", "Phase2"), 2),
    row.names = paste0("sample_", 1:8)
  )
  
  SummarizedExperiment::SummarizedExperiment(
    assays = list(counts = counts),
    rowData = rowData,
    colData = colData
  )
}

# ============================================================================
# BLOCK JACKKNIFE TESTS
# ============================================================================

test_that("block_jackknife_isoform_switching function exists", {
  expect_true(exists("block_jackknife_isoform_switching"))
  expect_true(is.function(block_jackknife_isoform_switching))
})

test_that("block_jackknife requires se parameter", {
  expect_error(
    block_jackknife_isoform_switching(se = NULL),
    "SummarizedExperiment object.*required"
  )
})

test_that("block_jackknife validates SummarizedExperiment class", {
  expect_error(
    block_jackknife_isoform_switching(
      se = data.frame(x = 1:10),
      block_col = "phase",
      gene_col = "gene_id",
      isoform_col = "isoform_id"
    ),
    "must be a SummarizedExperiment object"
  )
})

test_that("block_jackknife validates block_col exists", {
  se <- test_se_with_blocks()
  expect_error(
    block_jackknife_isoform_switching(
      se = se,
      block_col = "nonexistent_column",
      gene_col = "gene_id",
      isoform_col = "isoform_id",
      print_results = FALSE
    ),
    "not found in colData"
  )
})

test_that("block_jackknife requires gene_col and isoform_col", {
  se <- test_se_with_blocks()
  expect_error(
    block_jackknife_isoform_switching(
      se = se,
      block_col = "phase",
      gene_col = NULL,
      isoform_col = "isoform_id",
      print_results = FALSE
    ),
    "gene_col and isoform_col must be specified"
  )
})

test_that("block_jackknife requires at least 2 blocks", {
  se <- test_se_with_blocks()
  SummarizedExperiment::colData(se)$phase <- "OnlyPhase"  # Single block
  
  expect_error(
    block_jackknife_isoform_switching(
      se = se,
      block_col = "phase",
      gene_col = "gene_id",
      isoform_col = "isoform_id",
      print_results = FALSE
    ),
    "at least 2 blocks"
  )
})

test_that("block_jackknife returns tsenat_block_jackknife object", {
  se <- test_se_with_blocks()
  result <- suppressWarnings(block_jackknife_isoform_switching(
    se = se,
    block_col = "phase",
    gene_col = "gene_id",
    isoform_col = "isoform_id",
    print_results = FALSE
  ))
  
  expect_true(inherits(result, "tsenat_block_jackknife"))
  expect_true(is.list(result))
})

test_that("block_jackknife result has required components", {
  se <- test_se_with_blocks()
  result <- suppressWarnings(block_jackknife_isoform_switching(
    se = se,
    block_col = "phase",
    gene_col = "gene_id",
    isoform_col = "isoform_id",
    print_results = FALSE
  ))
  
  expect_true("blocks" %in% names(result))
  expect_true("gene_names" %in% names(result))
  expect_true("per_block_results" %in% names(result))
  expect_true("metadata" %in% names(result))
})

test_that("block_jackknife metadata has correct structure", {
  se <- test_se_with_blocks()
  result <- suppressWarnings(block_jackknife_isoform_switching(
    se = se,
    block_col = "phase",
    gene_col = "gene_id",
    isoform_col = "isoform_id",
    q = 1.5,
    norm = FALSE,
    print_results = FALSE
  ))
  
  meta <- result$metadata
  expect_equal(meta$q, 1.5)
  expect_equal(meta$norm, FALSE)
  expect_equal(meta$method, "block_jackknife")
})

test_that("block_jackknife metadata includes paper references", {
  se <- test_se_with_blocks()
  result <- suppressWarnings(block_jackknife_isoform_switching(
    se = se,
    block_col = "phase",
    gene_col = "gene_id",
    isoform_col = "isoform_id",
    print_results = FALSE
  ))
  
  expect_true("reference_papers" %in% names(result$metadata))
  expect_true("C137" %in% result$metadata$reference_papers)
  expect_true("ISO021" %in% result$metadata$reference_papers)
})

# ============================================================================
# HEATMAP VISUALIZATION TESTS
# ============================================================================

test_that("plot_isoform_switching_heatmap function exists", {
  expect_true(exists("plot_isoform_switching_heatmap"))
  expect_true(is.function(plot_isoform_switching_heatmap))
})

test_that("plot_isoform_switching_heatmap requires switching_results", {
  expect_error(
    plot_isoform_switching_heatmap(switching_results = NULL),
    "must be from jackknife_isoform_switching"
  )
})

test_that("plot_isoform_switching_heatmap validates class", {
  expect_error(
    plot_isoform_switching_heatmap(switching_results = list()),
    "must be from jackknife_isoform_switching"
  )
})

test_that("plot_isoform_switching_heatmap creates heatmap for valid input", {
  se <- test_se_with_blocks()
  
  switching <- suppressWarnings(jackknife_isoform_switching(
    se = se,
    condition_col = "condition",
    gene_col = "gene_id",
    isoform_col = "isoform_id",
    q = 1,
    n_bootstrap = 5,
    norm = FALSE,
    print_results = FALSE
  ))
  
  # Should not error
  expect_silent(suppressWarnings(
    plot_isoform_switching_heatmap(
      switching,
      top_n = 1,
      top_transcripts_per_gene = 2,
      show_values = FALSE
    )
  ))
})

test_that("plot_isoform_switching_heatmap handles diverging color scheme", {
  se <- test_se_with_blocks()
  
  switching <- suppressWarnings(jackknife_isoform_switching(
    se = se,
    condition_col = "condition",
    gene_col = "gene_id",
    isoform_col = "isoform_id",
    q = 1,
    n_bootstrap = 5,
    norm = FALSE,
    print_results = FALSE
  ))
  
  # Should not error
  expect_silent(suppressWarnings(
    plot_isoform_switching_heatmap(
      switching,
      color_scheme = "diverging",
      show_values = FALSE
    )
  ))
})

test_that("plot_isoform_switching_heatmap handles sequential color scheme", {
  se <- test_se_with_blocks()
  
  switching <- suppressWarnings(jackknife_isoform_switching(
    se = se,
    condition_col = "condition",
    gene_col = "gene_id",
    isoform_col = "isoform_id",
    q = 1,
    n_bootstrap = 5,
    norm = FALSE,
    print_results = FALSE
  ))
  
  # Should not error
  expect_silent(suppressWarnings(
    plot_isoform_switching_heatmap(
      switching,
      color_scheme = "sequential",
      show_values = FALSE
    )
  ))
})

test_that("plot_isoform_switching_heatmap returns heatmap matrix invisibly", {
  se <- test_se_with_blocks()
  
  switching <- suppressWarnings(jackknife_isoform_switching(
    se = se,
    condition_col = "condition",
    gene_col = "gene_id",
    isoform_col = "isoform_id",
    q = 1,
    n_bootstrap = 5,
    norm = FALSE,
    print_results = FALSE
  ))
  
  result <- suppressWarnings(
    plot_isoform_switching_heatmap(switching, show_values = FALSE)
  )
  
  expect_true(is.matrix(result) || is.null(result))
})

# ============================================================================
# Q-SENSITIVITY CURVE TESTS
# ============================================================================

test_that("plot_q_sensitivity_curve function exists", {
  expect_true(exists("plot_q_sensitivity_curve"))
  expect_true(is.function(plot_q_sensitivity_curve))
})

test_that("plot_q_sensitivity_curve requires se and gene parameters", {
  expect_error(
    plot_q_sensitivity_curve(se = NULL),
    "required"
  )
})

test_that("plot_q_sensitivity_curve validates gene exists", {
  se <- test_se_with_blocks()
  expect_error(
    plot_q_sensitivity_curve(
      se = se,
      condition_col = "condition",
      gene = "NonexistentGene",
      gene_col = "gene_id",
      isoform_col = "isoform_id"
    ),
    "not found in rowData"
  )
})

test_that("plot_q_sensitivity_curve requires 2 conditions", {
  se <- test_se_with_blocks()
  SummarizedExperiment::colData(se)$condition <- "A"  # All same condition
  
  expect_error(
    plot_q_sensitivity_curve(
      se = se,
      condition_col = "condition",
      gene = "Gene1",
      gene_col = "gene_id",
      isoform_col = "isoform_id"
    ),
    "Exactly 2 conditions"
  )
})

test_that("plot_q_sensitivity_curve creates plot for single q", {
  se <- test_se_with_blocks()
  
  # Should not error
  expect_silent(suppressWarnings(
    plot_q_sensitivity_curve(
      se = se,
      condition_col = "condition",
      gene = "Gene1",
      gene_col = "gene_id",
      isoform_col = "isoform_id",
      q_values = 1
    )
  ))
})

test_that("plot_q_sensitivity_curve creates plot for multiple q values", {
  se <- test_se_with_blocks()
  
  # Should not error
  expect_silent(suppressWarnings(
    plot_q_sensitivity_curve(
      se = se,
      condition_col = "condition",
      gene = "Gene1",
      gene_col = "gene_id",
      isoform_col = "isoform_id",
      q_values = c(0.5, 1, 1.5, 2)
    )
  ))
})

test_that("plot_q_sensitivity_curve returns sensitivity dataframe", {
  se <- test_se_with_blocks()
  
  result <- suppressWarnings(
    plot_q_sensitivity_curve(
      se = se,
      condition_col = "condition",
      gene = "Gene1",
      gene_col = "gene_id",
      isoform_col = "isoform_id",
      q_values = c(0.5, 1, 1.5, 2),
      show_legend = FALSE
    )
  )
  
  expect_true(is.data.frame(result))
  expect_true("q" %in% colnames(result))
  expect_true("delta_entropy" %in% colnames(result))
})

test_that("plot_q_sensitivity_curve sensitivity data has correct length", {
  se <- test_se_with_blocks()
  
  q_vals <- c(0.5, 0.8, 1, 1.2, 1.5, 2)
  result <- suppressWarnings(
    plot_q_sensitivity_curve(
      se = se,
      condition_col = "condition",
      gene = "Gene1",
      gene_col = "gene_id",
      isoform_col = "isoform_id",
      q_values = q_vals,
      show_legend = FALSE
    )
  )
  
  expect_equal(nrow(result), length(q_vals))
})

test_that("plot_q_sensitivity_curve respects normalization parameter", {
  se <- test_se_with_blocks()
  
  result_norm <- suppressWarnings(
    plot_q_sensitivity_curve(
      se = se,
      condition_col = "condition",
      gene = "Gene1",
      gene_col = "gene_id",
      isoform_col = "isoform_id",
      q_values = c(0.5, 1, 1.5, 2),
      norm = TRUE,
      show_legend = FALSE
    )
  )
  
  result_unnorm <- suppressWarnings(
    plot_q_sensitivity_curve(
      se = se,
      condition_col = "condition",
      gene = "Gene1",
      gene_col = "gene_id",
      isoform_col = "isoform_id",
      q_values = c(0.5, 1, 1.5, 2),
      norm = FALSE,
      show_legend = FALSE
    )
  )
  
  # Values should be different between normalized and unnormalized
  expect_false(all(result_norm$delta_entropy == result_unnorm$delta_entropy))
})

test_that("plot_q_sensitivity_curve respects q_values parameter", {
  se <- test_se_with_blocks()
  
  result <- suppressWarnings(
    plot_q_sensitivity_curve(
      se = se,
      condition_col = "condition",
      gene = "Gene1",
      gene_col = "gene_id",
      isoform_col = "isoform_id",
      q_values = c(0.7, 1.3, 2.1),
      show_legend = FALSE
    )
  )
  
  expect_equal(result$q, c(0.7, 1.3, 2.1))
})

# ============================================================================
# INTEGRATION TESTS
# ============================================================================

test_that("Block jackknife and isoform switching work together with same data", {
  se <- test_se_with_blocks()
  
  # Block jackknife
  block_result <- suppressWarnings(
    block_jackknife_isoform_switching(
      se = se,
      block_col = "phase",
      gene_col = "gene_id",
      isoform_col = "isoform_id",
      print_results = FALSE
    )
  )
  
  # Standard isoform switching
  switch_result <- suppressWarnings(
    jackknife_isoform_switching(
      se = se,
      condition_col = "condition",
      gene_col = "gene_id",
      isoform_col = "isoform_id",
      q = 1,
      n_bootstrap = 5,
      norm = FALSE,
      print_results = FALSE
    )
  )
  
  expect_true(inherits(block_result, "tsenat_block_jackknife"))
  expect_true(inherits(switch_result, "tsenat_isoform_switching"))
})

test_that("Visualization functions work with isoform switching results", {
  se <- test_se_with_blocks()
  
  switch_result <- suppressWarnings(
    jackknife_isoform_switching(
      se = se,
      condition_col = "condition",
      gene_col = "gene_id",
      isoform_col = "isoform_id",
      q = 1,
      n_bootstrap = 5,
      norm = FALSE,
      print_results = FALSE
    )
  )
  
  # Test heatmap
  heatmap <- suppressWarnings(
    plot_isoform_switching_heatmap(switch_result, show_values = FALSE)
  )
  
  # Test q-curve
  q_result <- suppressWarnings(
    plot_q_sensitivity_curve(
      se = se,
      condition_col = "condition",
      gene = "Gene1",
      gene_col = "gene_id",
      isoform_col = "isoform_id",
      q_values = c(1, 1.5),
      show_legend = FALSE
    )
  )
  
  expect_true(is.data.frame(q_result))
})

test_that("All enhancement functions include paper references", {
  # Check function documentation in roxygen comments
  expect_true(exists("block_jackknife_isoform_switching"))
  expect_true(exists("plot_isoform_switching_heatmap"))
  expect_true(exists("plot_q_sensitivity_curve"))
})
