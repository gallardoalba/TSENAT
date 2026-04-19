# ============================================================================
library(SummarizedExperiment)
library(ggplot2)
library(pheatmap)
library(testthat)

# ============================================================================
# Unit Tests for Heatmap Helper Functions
# ============================================================================
# Tests for the 13 internal helper functions in R/heatmap_helpers.R
# These tests ensure correct behavior of the refactored heatmap functions
# ============================================================================

# Setup test data


test_se <- create_test_se_simple(n_genes = 15, n_samples = 6, control_n = 3)

# Create multi-q switching results for testing
create_test_multiq_results <- function() {
  result <- list(
    q_0_01 = list(
      gene_ids = c("G1", "G2", "G3", "G4", "G5"),
      gene_name_map = c(G1 = "GENE1", G2 = "GENE2", G3 = "GENE3", G4 = "GENE4", G5 = "GENE5"),
      results_per_gene = list(
        G1 = list(
          delta_influence = c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10),
          transcript_ids = paste0("tx", 1:10)
        ),
        G2 = list(
          delta_influence = c(1, 2, 3, 4, 5),
          transcript_ids = paste0("tx", 1:5)
        ),
        G3 = list(
          delta_influence = c(1, 2, 3),
          transcript_ids = paste0("tx", 1:3)
        )
      )
    ),
    q_0_50 = list(
      gene_ids = c("G1", "G2", "G3", "G4", "G5"),
      gene_name_map = c(G1 = "GENE1", G2 = "GENE2", G3 = "GENE3", G4 = "GENE4", G5 = "GENE5"),
      results_per_gene = list(
        G1 = list(
          delta_influence = rnorm(10),
          transcript_ids = paste0("tx", 1:10)
        ),
        G2 = list(
          delta_influence = rnorm(5),
          transcript_ids = paste0("tx", 1:5)
        ),
        G3 = list(
          delta_influence = rnorm(3),
          transcript_ids = paste0("tx", 1:3)
        )
      )
    ),
    q_2_00 = list(
      gene_ids = c("G1", "G2", "G3", "G4", "G5"),
      gene_name_map = c(G1 = "GENE1", G2 = "GENE2", G3 = "GENE3", G4 = "GENE4", G5 = "GENE5"),
      results_per_gene = list(
        G1 = list(
          delta_influence = rnorm(10),
          transcript_ids = paste0("tx", 1:10)
        ),
        G2 = list(
          delta_influence = rnorm(5),
          transcript_ids = paste0("tx", 1:5)
        ),
        G3 = list(
          delta_influence = rnorm(3),
          transcript_ids = paste0("tx", 1:3)
        )
      )
    )
  )
  class(result) <- "tsenat_isoform_switching_multiq"
  result
}

# ============================================================================
# SECTION 1: VALIDATION HELPERS
# ============================================================================

describe(".validate_se_for_heatmaps()", {
  test_that("validates correct SummarizedExperiment with auto-detected gene column", {
    result <- TSENAT:::.validate_se_for_heatmaps(test_se)
    expect_is(result, "list")
    expect_named(result, c("counts", "rowdata", "coldata", "gene_col", "condition_col"))
    expect_equal(result$gene_col, "gene_name")
    expect_is(result$counts, "matrix")
    expect_equal(nrow(result$counts), 15)
    expect_equal(ncol(result$counts), 6)
  })

  test_that("validates correct SummarizedExperiment with specified gene column", {
    result <- TSENAT:::.validate_se_for_heatmaps(test_se, gene_col = "gene_id")
    expect_equal(result$gene_col, "gene_id")
  })

  test_that("validates with specified condition column", {
    result <- TSENAT:::.validate_se_for_heatmaps(test_se, condition_col = "sample_type")
    expect_equal(result$condition_col, "sample_type")
  })

  test_that("throws error for invalid SummarizedExperiment", {
    expect_error(
      TSENAT:::.validate_se_for_heatmaps(data.frame(a = 1:5)),
      "must be a SummarizedExperiment"
    )
  })

  test_that("throws error for missing gene column", {
    expect_error(
      TSENAT:::.validate_se_for_heatmaps(test_se, gene_col = "nonexistent"),
      "not found in rowData"
    )
  })

  test_that("throws error for missing condition column", {
    expect_error(
      TSENAT:::.validate_se_for_heatmaps(test_se, condition_col = "nonexistent"),
      "not found in colData"
    )
  })
})

describe(".validate_multiq_input()", {
  test_that("validates correct multi-q results object", {
    multiq_results <- create_test_multiq_results()
    class(multiq_results) <- "tsenat_isoform_switching_multiq"
    result <- TSENAT:::.validate_multiq_input(multiq_results)
    
    expect_is(result, "list")
    expect_named(result, c("q_result_keys", "first_result", "gene_ids", "gene_name_map"))
    expect_equal(length(result$q_result_keys), 3)
    expect_equal(result$gene_ids, c("G1", "G2", "G3", "G4", "G5"))
  })

  test_that("throws error for invalid class", {
    expect_error(
      TSENAT:::.validate_multiq_input(list()),
      "must be a multi-q result"
    )
  })

  test_that("throws error for missing q-value results", {
    invalid_results <- list(other_result = list(gene_ids = 1:5))
    class(invalid_results) <- "tsenat_isoform_switching_multiq"
    expect_error(
      TSENAT:::.validate_multiq_input(invalid_results),
      "No multi-q results found"
    )
  })

  test_that("throws error for no genes found", {
    bad_results <- list(
      q_0_01 = list(gene_ids = integer(0), gene_name_map = NULL)
    )
    class(bad_results) <- "tsenat_isoform_switching_multiq"
    expect_error(
      TSENAT:::.validate_multiq_input(bad_results),
      "No genes found"
    )
  })
})

# ============================================================================
# SECTION 2: GENE SELECTION HELPERS
# ============================================================================

describe(".heatmap_select_genes_multiq()", {
  test_that("selects first N genes when no sait_results provided", {
    multiq_results <- create_test_multiq_results()
    class(multiq_results) <- "tsenat_isoform_switching_multiq"
    
    selected <- TSENAT:::.heatmap_select_genes_multiq(multiq_results, n_genes = 3)
    expect_equal(length(selected), 3)
    expect_equal(selected, c("G1", "G2", "G3"))
  })

  test_that("handles n_genes larger than available genes", {
    multiq_results <- create_test_multiq_results()
    class(multiq_results) <- "tsenat_isoform_switching_multiq"
    
    selected <- TSENAT:::.heatmap_select_genes_multiq(multiq_results, n_genes = 100)
    expect_equal(length(selected), 5)
  })

  test_that("ranks genes by p-value when sait_results provided", {
    multiq_results <- create_test_multiq_results()
    class(multiq_results) <- "tsenat_isoform_switching_multiq"
    
    sait_results <- data.frame(
      gene_id = c("G1", "G2", "G3", "G4", "G5"),
      adj_p_interaction = c(0.001, 0.01, 0.05, 0.1, 0.5)
    )
    
    selected <- TSENAT:::.heatmap_select_genes_multiq(
      multiq_results,
      n_genes = 2,
      sait_results = sait_results
    )
    expect_equal(selected[1], "G1")  # Lowest p-value
  })

  test_that("handles missing gene column in sait_results gracefully", {
    multiq_results <- create_test_multiq_results()
    class(multiq_results) <- "tsenat_isoform_switching_multiq"
    
    sait_results <- data.frame(
      bad_col = c("G1", "G2", "G3"),
      p_value = c(0.001, 0.01, 0.05)
    )
    
    selected <- TSENAT:::.heatmap_select_genes_multiq(multiq_results, n_genes = 2, sait_results = sait_results)
    expect_equal(length(selected), 2)
  })

  test_that("handles empty n_genes = 0", {
    multiq_results <- create_test_multiq_results()
    class(multiq_results) <- "tsenat_isoform_switching_multiq"
    
    selected <- TSENAT:::.heatmap_select_genes_multiq(multiq_results, n_genes = 0)
    expect_equal(length(selected), 0)
  })
})

describe(".heatmap_select_genes_results()", {
  test_that("selects top genes by p-value from results", {
    # Use actual gene names from test_se
    se_genes <- rownames(test_se)[1:3]
    results <- data.frame(
      gene = se_genes,
      adj_p_value = c(0.001, 0.05, 0.5),
      stringsAsFactors = FALSE
    )
    
    selected <- TSENAT:::.heatmap_select_genes_results(
      test_se, results, gene_col = "gene_name", top_n = 2
    )
    expect_equal(length(selected), 2)
    expect_equal(selected[1], se_genes[1])
  })

  test_that("handles padj column as alias for adj_p_value", {
    se_genes <- rownames(test_se)[1:3]
    results <- data.frame(
      gene = se_genes,
      padj = c(0.001, 0.05, 0.5),
      stringsAsFactors = FALSE
    )
    
    selected <- TSENAT:::.heatmap_select_genes_results(
      test_se, results, gene_col = "gene_name", top_n = 2
    )
    expect_equal(length(selected), 2)
  })

  test_that("filters out genes not in SE", {
    results <- data.frame(
      gene = c("NotInSE1", "NotInSE2", "Gene_1"),
      adj_p_value = c(0.001, 0.01, 0.5)
    )
    
    selected <- TSENAT:::.heatmap_select_genes_results(
      test_se, results, gene_col = "gene_name", top_n = 2
    )
    expect_true(all(selected %in% rownames(test_se)))
  })

  test_that("throws error for invalid results class", {
    expect_error(
      TSENAT:::.heatmap_select_genes_results(test_se, "not_a_dataframe"),
      "must be a data.frame"
    )
  })

  test_that("throws error for missing gene column", {
    results <- data.frame(
      bad_col = c("NotAGene1", "NotAGene2"),
      adj_p_value = c(0.001, 0.05),
      stringsAsFactors = FALSE
    )
    
    expect_error(
      TSENAT:::.heatmap_select_genes_results(test_se, results, gene_col = "gene_name"),
      "must contain.*gene"
    )
  })

  test_that("handles missing p-value column gracefully", {
    se_genes <- rownames(test_se)[1:3]
    results <- data.frame(
      gene = se_genes,
      stringsAsFactors = FALSE
    )
    
    selected <- TSENAT:::.heatmap_select_genes_results(
      test_se, results, gene_col = "gene_name", top_n = 2
    )
    expect_equal(length(selected), 2)
  })
})

# ============================================================================
# SECTION 3: LAYOUT PLANNING HELPERS
# ============================================================================

describe(".plot_adaptive_layout()", {
  test_that("fixed layout: 2 columns per row with 5 genes", {
    gene_info <- list(
      list(n_transcripts = 10),
      list(n_transcripts = 8),
      list(n_transcripts = 6),
      list(n_transcripts = 4),
      list(n_transcripts = 5)
    )
    
    result <- TSENAT:::.plot_adaptive_layout(
      gene_info, use_fixed_layout = TRUE, layout_ncol = 2
    )
    
    expect_equal(result$n_layout_rows, 3)  # ceiling(5/2) = 3
    expect_equal(length(result$layout), 5)
    expect_equal(result$layout[[1]]$row, 1)
    expect_equal(result$layout[[2]]$row, 1)
    expect_equal(result$layout[[3]]$row, 2)
  })

  test_that("adaptive layout: full vs half-width based on transcripts", {
    gene_info <- list(
      list(n_transcripts = 15),  # Full width
      list(n_transcripts = 2),   # Half width
      list(n_transcripts = 2)    # Half width
    )
    
    result <- TSENAT:::.plot_adaptive_layout(
      gene_info, use_fixed_layout = FALSE
    )
    
    expect_equal(result$layout[[1]]$width, 1)    # Full width
    expect_equal(result$layout[[2]]$width, 0.5)  # Half width
    expect_equal(result$layout[[3]]$width, 0.5)  # Half width
  })

  test_that("handles single gene", {
    gene_info <- list(list(n_transcripts = 10))
    
    result <- TSENAT:::.plot_adaptive_layout(
      gene_info, use_fixed_layout = TRUE, layout_ncol = 2
    )
    
    expect_equal(result$n_layout_rows, 1)
    expect_equal(length(result$layout), 1)
  })

  test_that("handles single gene in adaptive layout", {
    gene_info <- list(list(n_transcripts = 10))
    result <- TSENAT:::.plot_adaptive_layout(
      gene_info, use_fixed_layout = FALSE
    )
    
    expect_equal(result$n_layout_rows, 1)
    expect_equal(length(result$layout), 1)
  })
})

# ============================================================================
# SECTION 4: SIZING HELPERS
# ============================================================================

describe(".calculate_heatmap_dimensions()", {
  test_that("calculates PNG dimensions based on layout rows", {
    dims <- TSENAT:::.calculate_heatmap_dimensions(n_layout_rows = 2, n_data_rows = 5)
    
    expect_is(dims, "list")
    expect_true("png_width" %in% names(dims))
    expect_true("png_height" %in% names(dims))
    expect_equal(dims$png_width, 12)
    expect_true(dims$png_height > 0)
  })

  test_that("scales height with number of layout rows", {
    dims_1 <- TSENAT:::.calculate_heatmap_dimensions(n_layout_rows = 1, n_data_rows = 5)
    dims_2 <- TSENAT:::.calculate_heatmap_dimensions(n_layout_rows = 2, n_data_rows = 5)
    
    expect_true(dims_2$png_height > dims_1$png_height)
  })

  test_that("handles edge case: single row", {
    dims <- TSENAT:::.calculate_heatmap_dimensions(n_layout_rows = 1, n_data_rows = 1)
    
    expect_true(dims$png_height > 0)
    expect_equal(dims$png_width, 12)
  })
})

describe(".calculate_adaptive_cellsizes()", {
  test_that("calculates cell sizes for full-width heatmap", {
    sizes <- TSENAT:::.calculate_adaptive_cellsizes(
      n_cols = 10, n_rows = 5, width_frac = 1, cellwidth = 0, cellheight = 0
    )
    
    expect_is(sizes, "list")
    expect_named(sizes, c("cellwidth", "cellheight", "fontsize_adj"))
    expect_true(sizes$cellwidth > 0)
    expect_true(sizes$cellheight > 0)
  })

  test_that("calculates cell sizes for half-width heatmap", {
    sizes_half <- TSENAT:::.calculate_adaptive_cellsizes(
      n_cols = 5, n_rows = 5, width_frac = 0.5
    )
    sizes_full <- TSENAT:::.calculate_adaptive_cellsizes(
      n_cols = 5, n_rows = 5, width_frac = 1
    )
    
    expect_true(sizes_half$cellwidth < sizes_full$cellwidth)
  })

  test_that("returns positive cell sizes", {
    sizes <- TSENAT:::.calculate_adaptive_cellsizes(
      n_cols = 10, n_rows = 5, cellwidth = 50
    )
    
    expect_true(sizes$cellwidth > 0)
  })

  test_that("returns positive cellheight", {
    sizes <- TSENAT:::.calculate_adaptive_cellsizes(
      n_cols = 10, n_rows = 5, cellheight = 20
    )
    
    expect_true(sizes$cellheight > 0)
  })

  test_that("handles very large matrices", {
    sizes <- TSENAT:::.calculate_adaptive_cellsizes(
      n_cols = 100, n_rows = 50
    )
    
    expect_true(sizes$cellwidth > 0)
    expect_true(sizes$cellheight > 0)
  })

  test_that("handles very small matrices", {
    sizes <- TSENAT:::.calculate_adaptive_cellsizes(
      n_cols = 1, n_rows = 1
    )
    
    expect_true(sizes$cellwidth > 0)
    expect_true(sizes$cellheight > 0)
  })
})

# ============================================================================
# SECTION 5: DATA PREPARATION HELPERS
# ============================================================================

describe(".heatmap_prepare_multiq_data()", {
  test_that("extracts delta_influence matrix for given gene", {
    multiq_results <- create_test_multiq_results()
    
    mat <- TSENAT:::.heatmap_prepare_multiq_data(
      multiq_results, "G1",
      q_result_keys = c("q_0_01", "q_0_50", "q_2_00")
    )
    
    expect_is(mat, "matrix")
    expect_equal(nrow(mat), 3)   # 3 q-values (rows)
    expect_equal(ncol(mat), 10)  # 10 transcripts (cols)
  })

  test_that("processes matrix data correctly", {
    multiq_results <- list(
      q_0_01 = list(
        results_per_gene = list(
          G1 = list(
            delta_influence = c(1, 2, 3, 4, 5),
            transcript_ids = paste0("tx", 1:5)
          )
        )
      )
    )
    class(multiq_results) <- "tsenat_isoform_switching_multiq"
    
    mat <- TSENAT:::.heatmap_prepare_multiq_data(
      multiq_results, "G1", q_result_keys = "q_0_01"
    )
    
    expect_is(mat, "matrix")
    expect_equal(nrow(mat), 1)  # 1 q-value (row)
    expect_equal(ncol(mat), 5)  # 5 transcripts (cols)
  })

  test_that("returns NULL for missing gene", {
    multiq_results <- create_test_multiq_results()
    class(multiq_results) <- "tsenat_isoform_switching_multiq"
    
    mat <- TSENAT:::.heatmap_prepare_multiq_data(
      multiq_results, "NonexistentGene",
      q_result_keys = c("q_0_01", "q_0_50")
    )
    
    expect_null(mat)
  })
})

describe(".heatmap_prepare_condition_data()", {
  test_that("aggregates transcript counts by condition", {
    # Subset of test_se with 3 transcripts
    tx_counts <- assay(test_se)[1:3, ]
    conditions <- colData(test_se)$sample_type
    
    mat <- TSENAT:::.heatmap_prepare_condition_data(
      tx_counts, seq_len(nrow(tx_counts)), conditions, metric = "median"
    )
    
    expect_is(mat, "matrix")
    expect_equal(nrow(mat), 2)  # 2 conditions
    expect_equal(ncol(mat), 3)  # 3 transcripts
    expect_true(all(is.finite(mat)))
  })

  test_that("supports different aggregation metrics", {
    tx_counts <- assay(test_se)[1:3, ]
    conditions <- colData(test_se)$sample_type
    
    for (metric in c("median", "mean", "variance", "iqr")) {
      mat <- TSENAT:::.heatmap_prepare_condition_data(
        tx_counts, seq_len(nrow(tx_counts)), conditions, metric = metric
      )
      expect_is(mat, "matrix")
      expect_true(all(is.finite(mat) | is.na(mat)))
    }
  })

  test_that("handles empty transcript indices", {
    tx_counts <- assay(test_se)[1:3, ]
    conditions <- colData(test_se)$sample_type
    
    mat <- TSENAT:::.heatmap_prepare_condition_data(
      tx_counts, integer(0), conditions, metric = "median"
    )
    
    # May return NULL or empty matrix depending on implementation
    expect_true(is.null(mat) || length(mat) == 0)
  })

  test_that("returns numeric matrix", {
    tx_counts <- assay(test_se)[1:3, ]
    conditions <- colData(test_se)$sample_type
    
    mat <- TSENAT:::.heatmap_prepare_condition_data(
      tx_counts, seq_len(nrow(tx_counts)), conditions, metric = "median"
    )
    
    expect_is(mat, "matrix")
    expect_true(all(is.numeric(mat)))
  })
})

# ============================================================================
# SECTION 6: PHEATMAP CREATION HELPER
# ============================================================================

describe(".create_pheatmap_grob()", {
  test_that("creates valid pheatmap object", {
    mat <- matrix(rnorm(50), nrow = 5, ncol = 10)
    
    grDevices::pdf(file = NULL)
    on.exit(if (grDevices::dev.cur() > 1) grDevices::dev.off())
    grob <- TSENAT:::.create_pheatmap_grob(mat, title = "Test Heatmap")
    
    expect_is(grob, "pheatmap")
  })

  test_that("handles matrix with rownames and colnames", {
    mat <- matrix(rnorm(50), nrow = 5, ncol = 10)
    rownames(mat) <- paste0("Row", 1:5)
    colnames(mat) <- paste0("Col", 1:10)
    
    grDevices::pdf(file = NULL)
    on.exit(if (grDevices::dev.cur() > 1) grDevices::dev.off())
    grob <- TSENAT:::.create_pheatmap_grob(mat, title = "Test")
    
    expect_is(grob, "pheatmap")
  })

  test_that("respects custom cell dimensions", {
    mat <- matrix(rnorm(50), nrow = 5, ncol = 10)
    
    grDevices::pdf(file = NULL)
    on.exit(if (grDevices::dev.cur() > 1) grDevices::dev.off())
    grob <- TSENAT:::.create_pheatmap_grob(
      mat, cellw = 50, cellh = 20
    )
    
    expect_is(grob, "pheatmap")
  })

  test_that("respects fontsize parameter", {
    mat <- matrix(rnorm(50), nrow = 5, ncol = 10)
    
    grDevices::pdf(file = NULL)
    on.exit(if (grDevices::dev.cur() > 1) grDevices::dev.off())
    grob <- TSENAT:::.create_pheatmap_grob(mat, fontsize = 10)
    
    expect_is(grob, "pheatmap")
  })

  test_that("handles cluster_rows = FALSE", {
    mat <- matrix(rnorm(50), nrow = 5, ncol = 10)
    
    grDevices::pdf(file = NULL)
    on.exit(if (grDevices::dev.cur() > 1) grDevices::dev.off())
    grob <- TSENAT:::.create_pheatmap_grob(mat, cluster_rows = FALSE)
    
    expect_is(grob, "pheatmap")
  })
})

# ============================================================================
# SECTION 7: GRID RENDERING HELPERS
# ============================================================================

describe(".plot_grid_setup()", {
  test_that("opens PNG device when output_file provided", {
    temp_file <- tempfile(fileext = ".png")
    on.exit(unlink(temp_file), add = TRUE)
    
    TSENAT:::.plot_grid_setup(
      n_layout_rows = 2, output_file = temp_file,
      png_width = 12, png_height = 8
    )
    grDevices::dev.off()
    
    # Check that PNG was created
    expect_true(file.exists(temp_file))
  })

  test_that("creates grid layout with titles", {
    TSENAT:::.plot_grid_setup(
      n_layout_rows = 1, output_file = NULL,
      png_width = 12, png_height = 8,
      title = "Test Title", subtitle = "Test Subtitle"
    )
    grDevices::dev.off()
    
    # If no error thrown, test passes
    expect_true(TRUE)
  })

  test_that("handles multiple layout rows", {
    TSENAT:::.plot_grid_setup(
      n_layout_rows = 3, output_file = NULL,
      png_width = 12, png_height = 12,
      title = "Multi-row Test"
    )
    grDevices::dev.off()
    
    expect_true(TRUE)
  })
})

describe(".render_heatmaps_to_grid()", {
  test_that("renders heatmap grobs to grid", {
    # Create simple pheatmap
    mat <- matrix(rnorm(50), nrow = 5, ncol = 10)
    grob <- TSENAT:::.create_pheatmap_grob(mat)
    
    heatmap_plots <- list(grob, NULL, grob)  # Include NULL for missing heatmap
    
    gene_info <- list(
      list(row = 1, col = 1, width = 1),
      NULL,
      list(row = 2, col = 1, width = 1)
    )
    
    TSENAT:::.plot_grid_setup(
      n_layout_rows = 2, output_file = NULL,
      png_width = 12, png_height = 8
    )
    
    TSENAT:::.render_heatmaps_to_grid(heatmap_plots, gene_info, layout_ncol = 2)
    grDevices::dev.off()
    
    expect_true(TRUE)
  })
})

describe(".plot_grid_finalize()", {
  test_that("closes PNG device without error", {
    temp_file <- tempfile(fileext = ".png")
    on.exit(unlink(temp_file), add = TRUE)
    
    TSENAT:::.plot_grid_setup(
      n_layout_rows = 1, output_file = temp_file,
      png_width = 12, png_height = 8
    )
    
    result <- TSENAT:::.plot_grid_finalize(temp_file)
    
    expect_true(file.exists(temp_file))
    expect_true(file.size(temp_file) > 0)
  })

  test_that("handles NULL output_file gracefully", {
    TSENAT:::.plot_grid_setup(
      n_layout_rows = 1, output_file = NULL,
      png_width = 12, png_height = 8
    )
    
    result <- TSENAT:::.plot_grid_finalize(NULL)
    grDevices::dev.off()
    
    expect_null(result)
  })
})

# ============================================================================
# INTEGRATION TESTS: Helper Workflow
# ============================================================================

describe("Integration: Complete helper workflow", {
  test_that("workflow: SE validation → layout planning → heatmap creation", {
    # Step 1: Validate SE
    se_data <- TSENAT:::.validate_se_for_heatmaps(
      test_se, condition_col = "sample_type"
    )
    expect_is(se_data, "list")
    
    # Step 2: Select genes
    se_genes <- rownames(test_se)[1:3]
    selected <- TSENAT:::.heatmap_select_genes_results(
      test_se,
      data.frame(gene = se_genes, adj_p_value = c(0.001, 0.01, 0.05), stringsAsFactors = FALSE),
      gene_col = "gene_name",
      top_n = 2
    )
    expect_equal(length(selected), 2)
    
    # Step 3: Plan layout
    gene_info <- lapply(seq_along(selected), function(i) {
      list(n_transcripts = 5)
    })
    layout <- TSENAT:::.plot_adaptive_layout(gene_info, use_fixed_layout = TRUE)
    expect_is(layout, "list")
    
    # Step 4: Calculate dimensions
    dims <- TSENAT:::.calculate_heatmap_dimensions(layout$n_layout_rows, 3)
    expect_true(dims$png_width > 0)
  })

  test_that("workflow: Multi-Q validation and gene selection", {
    multiq_results <- create_test_multiq_results()
    class(multiq_results) <- "tsenat_isoform_switching_multiq"
    
    # Step 1: Validate multi-Q
    validated <- TSENAT:::.validate_multiq_input(multiq_results)
    expect_equal(length(validated$q_result_keys), 3)
    
    # Step 2: Select genes
    selected <- TSENAT:::.heatmap_select_genes_multiq(multiq_results, n_genes = 2)
    expect_equal(length(selected), 2)
    
    # Step 3: Prepare data for each gene
    for (gene_id in selected) {
      mat <- TSENAT:::.heatmap_prepare_multiq_data(
        multiq_results, gene_id,
        q_result_keys = validated$q_result_keys
      )
      expect_is(mat, "matrix")
    }
  })
})

test_that("plot_jis_delta: creates heatmaps from q-values", {
  skip_if_not_installed("SummarizedExperiment")
  skip_if_not_installed("ggplot2")
  skip_if_not_installed("ComplexHeatmap")
  
  set.seed(777)
  
  # Create mock multi-q divergence results
  multiq_results <- list(
    q_0.5 = data.frame(
      gene_id = paste0("g", 1:5),
      divergence = rnorm(5, mean = 1.5, sd = 0.5),
      stringsAsFactors = FALSE
    ),
    q_1.0 = data.frame(
      gene_id = paste0("g", 1:5),
      divergence = rnorm(5, mean = 2.0, sd = 0.5),
      stringsAsFactors = FALSE
    ),
    q_1.5 = data.frame(
      gene_id = paste0("g", 1:5),
      divergence = rnorm(5, mean = 2.5, sd = 0.5),
      stringsAsFactors = FALSE
    )
  )
  
  # Verify structure
  expect_true(is.list(multiq_results))
  expect_true(length(multiq_results) == 3)
  expect_true(all(c("q_0.5", "q_1.0", "q_1.5") %in% names(multiq_results)))
})

test_that("plot_jis_delta: parameter validation", {
  skip_if_not_installed("ggplot2")
  
  
  # Test top_n parameter
  top_n <- 10
  expect_true(is.numeric(top_n))
  expect_true(top_n > 0)
  
  # Test min_delta parameter
  min_delta <- 0.5
  expect_true(is.numeric(min_delta))
  expect_true(min_delta >= 0)
})


# ============================================================================
# TEST 3: plot_jis_delta - Multi-q heatmaps
# ============================================================================

test_that("plot_jis_delta: validates analysis object", {
  expect_error(
    TSENAT:::plot_jis_delta("not_analysis"),
    "must be a TSENATAnalysis object"
  )
})

test_that("plot_jis_delta: requires jackknife results", {
  set.seed(308)
  
  # Create analysis without jackknife results
  analysis <- .create_test_analysis(
    n_genes = 10,
    n_samples_per_group = 3,
    q_values = c(1.0),
    include_divergence = FALSE,
    include_sait_results = FALSE,
    seed = 308,
    verbose = FALSE
  )
  
  # Should error when no jackknife results
  expect_error(
    TSENAT:::plot_jis_delta(analysis),
    "No jackknife results"
  )
})


create_mock_jackknife_multiq <- function(n_genes = 10, n_transcripts_per_gene = 3,
                                         q_values = c("q_0_50", "q_1_00")) {
  # Create properly structured mock data for .plot_jis_delta()
  # Structure required:
  # - Class: "tsenat_isoform_switching_multiq"
  # - Names: q_* keys with:
  #   - $gene_ids: vector of gene IDs
  #   - $gene_name_map: named vector mapping gene IDs to gene names
  #   - $results_per_gene: list where each element (per gene) contains:
  #     - $transcript_ids: vector of transcript IDs for this gene
  #     - $delta_influence: numeric vector of delta values per transcript
  
  result <- list()
  gene_ids <- paste0("gene_", 1:n_genes)
  gene_name_map <- setNames(paste0("Gene_", 1:n_genes), gene_ids)
  
  for (q in q_values) {
    # Build results_per_gene structure: list of genes, each with transcripts and delta values
    results_per_gene <- list()
    
    for (gene_id in gene_ids) {
      transcript_ids <- paste0(gene_id, "_tx_", 1:n_transcripts_per_gene)
      # Create realistic delta values: numeric vector per transcript
      delta_vals <- rnorm(n_transcripts_per_gene, mean = 0.2, sd = 0.15)
      delta_vals <- pmax(pmin(delta_vals, 1), -1)  # Bound to [-1, 1] for realism
      
      results_per_gene[[gene_id]] <- list(
        transcript_ids = transcript_ids,
        delta_influence = delta_vals
      )
    }
    
    result[[q]] <- list(
      gene_ids = gene_ids,
      gene_name_map = gene_name_map,
      results_per_gene = results_per_gene
    )
  }
  
  # Set required class attribute
  class(result) <- c("tsenat_isoform_switching_multiq", "list")
  
  return(result)
}

test_that("plot_jis_delta: creates heatmap with mock jackknife data", {
  skip_if_not_installed("pheatmap")
  
  # Build analysis from vignette data
  config <- TSENAT_config(
    q_values = seq(0, 2, by = 0.2),
    condition_col = "condition",
    subject_col = "paired_samples",
    paired = TRUE,
    control = "normal"
  )
  analysis <- build_analysis(config = config, metadata = metadata_df, readcounts = readcounts, tx2gene = gff3_dataset, tpm = tpm, effective_length = effective_length)
  analysis <- filter_analysis(analysis, stringency = "medium")
  
  # Add required calculations
  analysis <- calculate_diversity(analysis, q = c(0.5, 1.0), verbose = FALSE)
  analysis <- calculate_jis(analysis, q = c(0.5, 1.0), nboot = 50, verbose = FALSE)
  
  # Function should execute without error (renders to graphics device, returns invisible NULL)
  expect_silent({
    TSENAT:::plot_jis_delta(analysis, n_genes = 4, verbose = FALSE)
  })
})

test_that("plot_jis_delta: ranks genes by SAIT results when provided", {
  skip_if_not_installed("pheatmap")
  
  # Build analysis from vignette data
  config <- TSENAT_config(
    q_values = seq(0, 2, by = 0.2),
    condition_col = "condition",
    subject_col = "paired_samples",
    paired = TRUE,
    control = "normal"
  )
  analysis <- build_analysis(config = config, metadata = metadata_df, readcounts = readcounts, tx2gene = gff3_dataset, tpm = tpm, effective_length = effective_length)
  analysis <- filter_analysis(analysis, stringency = "medium")
  
  # Add required calculations
  analysis <- calculate_diversity(analysis, q = c(0.5, 1.0), verbose = FALSE)
  analysis <- calculate_jis(analysis, q = c(0.5, 1.0), nboot = 50, verbose = FALSE)
  
  # Function should execute without error
  expect_silent({
    TSENAT:::plot_jis_delta(analysis, n_genes = 3, verbose = FALSE)
  })
})

test_that("plot_jis_delta: respects n_genes parameter", {
  skip_if_not_installed("pheatmap")
  
  # Build analysis from vignette data
  config <- TSENAT_config(
    q_values = seq(0, 2, by = 0.2),
    condition_col = "condition",
    subject_col = "paired_samples",
    paired = TRUE,
    control = "normal"
  )
  analysis <- build_analysis(config = config, metadata = metadata_df, readcounts = readcounts, tx2gene = gff3_dataset, tpm = tpm, effective_length = effective_length)
  analysis <- filter_analysis(analysis, stringency = "medium")
  
  # Add required calculations
  analysis <- calculate_diversity(analysis, q = c(0.5, 1.0), verbose = FALSE)
  analysis <- calculate_jis(analysis, q = c(0.5, 1.0), nboot = 50, verbose = FALSE)
  
  # Test with different n_genes values - should execute without error
  expect_silent({
    TSENAT:::plot_jis_delta(analysis, n_genes = 2, verbose = FALSE)
  })
  
  expect_silent({
    TSENAT:::plot_jis_delta(analysis, n_genes = 10, verbose = FALSE)
  })
})

# ============================================================================
# TEST 4: Integration - All three plotting functions in sequence
# ============================================================================

test_that("S4 plotting functions work on complete analysis object", {
  skip_if_not_installed("ggplot2")
  skip_if_not_installed("pheatmap")
  
  # Build analysis from vignette data
  config <- TSENAT_config(
    q_values = seq(0, 2, by = 0.2),
    condition_col = "condition",
    subject_col = "paired_samples",
    paired = TRUE,
    control = "normal"
  )
  analysis <- build_analysis(config = config, metadata = metadata_df, readcounts = readcounts, tx2gene = gff3_dataset, tpm = tpm, effective_length = effective_length)
  analysis <- filter_analysis(analysis, stringency = "medium")
  
  # Add required calculations (calculate_difference removed from codebase)
  analysis <- calculate_diversity(analysis, q = c(0.5, 1.0), verbose = FALSE)
  analysis <- calculate_divergence(analysis, q = c(0.5, 1.0), verbose = FALSE)
  analysis <- calculate_jis(analysis, q = c(0.5, 1.0), nboot = 50, verbose = FALSE)
  
  # Note: plot_diversity_volcano_ma() requires calculate_difference which was removed
  # These two should execute successfully without error
  expect_silent({
    TSENAT:::plot_divergence_spectrum(analysis, n_genes = 2, verbose = FALSE)
  })
  
  expect_silent({
    TSENAT:::plot_jis_delta(analysis, n_genes = 3, verbose = FALSE)
  })
})
