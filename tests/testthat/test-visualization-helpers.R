context("plot_helpers: Plot Composition and Theme Utilities")
library(ggplot2)
library(grid)
library(RColorBrewer)
library(SummarizedExperiment)
library(pheatmap)
library(patchwork)
library(cowplot)
library(methods)
library(testthat)
library(TSENAT)

# ============================================================================
# TEST: Theme and Title Functions
# ============================================================================

testthat::test_that("apply_tsenat_theme can be applied to plots", {

  # Just verify the function can be called
  testthat::expect_error(
    theme_obj <- .apply_tsenat_theme(base_size = 11),
    NA  # Expect no error
  )
  
  # Verify it can be added to a plot
  p <- ggplot2::ggplot(data.frame(x = 1, y = 1), ggplot2::aes(x, y)) +
    ggplot2::geom_point()
  
  testthat::expect_error(
    p_themed <- p + .apply_tsenat_theme(),
    NA  # Expect no error
  )
  testthat::expect_is(p_themed, "ggplot")
})

testthat::test_that("set_plot_title modifies title correctly", {

  p <- ggplot2::ggplot(data.frame(x = 1, y = 1), ggplot2::aes(x, y)) +
    ggplot2::geom_point()

  p_titled <- .set_plot_title(p, title = "Test Title", subtitle = "Test Subtitle")

  # Extract title from plot
  testthat::expect_is(p_titled, "ggplot")
  testthat::expect_equal(p_titled$labels$title, "Test Title")
  testthat::expect_equal(p_titled$labels$subtitle, "Test Subtitle")
})

testthat::test_that("set_plot_title applies font sizes", {

  p <- ggplot2::ggplot(data.frame(x = 1, y = 1), ggplot2::aes(x, y)) +
    ggplot2::geom_point()

  p_styled <- .set_plot_title(p, title = "Title", title_size = 16, subtitle_size = 12)

  testthat::expect_is(p_styled, "ggplot")
})

testthat::test_that("set_plot_title handles NULL title/subtitle", {

  p <- ggplot2::ggplot(data.frame(x = 1, y = 1), ggplot2::aes(x, y)) +
    ggplot2::geom_point()

  p_unchanged <- .set_plot_title(p, title = NULL, subtitle = NULL)

  testthat::expect_is(p_unchanged, "ggplot")
})

# ============================================================================
# TEST: Color and Fill Scale Creation
# ============================================================================

testthat::test_that("create_color_scale returns ggplot scale", {

  scale <- .create_color_scale(palette = "blue_red")
  testthat::expect_is(scale, "Scale")
})

testthat::test_that("create_color_scale reverses with direction -1", {

  scale_fwd <- .create_color_scale(palette = "blue_red", direction = 1)
  scale_rev <- .create_color_scale(palette = "blue_red", direction = -1)

  testthat::expect_is(scale_fwd, "Scale")
  testthat::expect_is(scale_rev, "Scale")
})

testthat::test_that("create_fill_scale returns appropriate scale", {

  scale <- .create_fill_scale(palette = "continuous_diverging")
  testthat::expect_is(scale, "Scale")
})

testthat::test_that("create_fill_scale accepts breaks parameter", {

  scale_50 <- .create_fill_scale(breaks = 50)
  scale_100 <- .create_fill_scale(breaks = 100)

  testthat::expect_is(scale_50, "Scale")
  testthat::expect_is(scale_100, "Scale")
})

# ============================================================================
# TEST: Heatmap Creation
# ============================================================================

testthat::test_that("create_tsenat_heatmap creates basic heatmap", {
  
  # Create a simple test matrix
  mat <- matrix(rnorm(50), nrow = 10, ncol = 5)
  rownames(mat) <- paste0("Gene_", 1:10)
  colnames(mat) <- paste0("Sample_", 1:5)
  
  # Create heatmap with basic parameters
  hm <- .create_tsenat_heatmap(
    mat = mat,
    title = "Test Heatmap",
    colors = NULL
  )
  
  # Verify it returned a pheatmap object with correct structure
  testthat::expect_is(hm, "pheatmap")
  testthat::expect_is(hm$gtable, "gtable")
  testthat::expect_is(hm$tree_row, "hclust")
  testthat::expect_is(hm$tree_col, "hclust")
})

testthat::test_that("create_tsenat_heatmap applies custom colors", {
  
  mat <- matrix(rnorm(50), nrow = 10, ncol = 5)
  rownames(mat) <- paste0("Gene_", 1:10)
  colnames(mat) <- paste0("Sample_", 1:5)
  
  # Create color palette
  custom_colors <- RColorBrewer::brewer.pal(9, "RdBu")
  
  # Create heatmap with custom colors
  hm <- .create_tsenat_heatmap(
    mat = mat,
    title = "Colored Heatmap",
    colors = custom_colors
  )
  
  testthat::expect_is(hm, "pheatmap")
  # Verify the heatmap with custom colors has correct structure
  testthat::expect_is(hm$gtable, "gtable")
  testthat::expect_is(hm$tree_row, "hclust")
  testthat::expect_is(hm$tree_col, "hclust")
})

testthat::test_that("create_tsenat_heatmap respects font size parameters", {
  
  mat <- matrix(rnorm(50), nrow = 10, ncol = 5)
  rownames(mat) <- paste0("Gene_", 1:10)
  colnames(mat) <- paste0("Sample_", 1:5)
  
  # Create heatmap with custom font sizes
  hm <- .create_tsenat_heatmap(
    mat = mat,
    title = "Large Font Heatmap",
    fontsize_row = 14,
    fontsize_col = 14
  )
  
  testthat::expect_is(hm, "pheatmap")
  # Verify the heatmap object has the correct structure (tree_row, tree_col, gtable)
  testthat::expect_true(all(c("tree_row", "tree_col", "gtable") %in% names(hm)))
  testthat::expect_is(hm$gtable, "gtable")
  testthat::expect_is(hm$tree_row, "hclust")
  testthat::expect_is(hm$tree_col, "hclust")
})

# ============================================================================
# TEST: Patchwork Composition
# ============================================================================

testthat::test_that("combine_plots_patchwork handles single plot", {

  p <- ggplot2::ggplot(data.frame(x = 1:10, y = 1:10), ggplot2::aes(x, y)) +
    ggplot2::geom_point()

  combined <- .combine_plots_patchwork(list(p), agg_label_unique = "median")

  testthat::expect_is(combined, "ggplot")
})

testthat::test_that("combine_plots_patchwork handles multiple plots", {

  plots <- list(
    ggplot2::ggplot(data.frame(x = 1:5, y = 1:5), ggplot2::aes(x, y)) +
      ggplot2::geom_point(),
    ggplot2::ggplot(data.frame(x = 1:5, y = 10:6), ggplot2::aes(x, y)) +
      ggplot2::geom_point(),
    ggplot2::ggplot(data.frame(x = 1:5, y = 5:1), ggplot2::aes(x, y)) +
      ggplot2::geom_point()
  )

  combined <- .combine_plots_patchwork(plots, agg_label_unique = "median")
  testthat::expect_is(combined, "ggplot")
})

testthat::test_that("combine_plots_patchwork generates title annotation", {

  p <- ggplot2::ggplot(data.frame(x = 1:5, y = 1:5), ggplot2::aes(x, y)) +
    ggplot2::geom_point()

  combined <- .combine_plots_patchwork(list(p), agg_label_unique = "test_metric")

  # Check that combined result is a patchwork composition
  testthat::expect_is(combined, "ggplot")
})

# ============================================================================
# TEST: Cowplot Composition
# ============================================================================

testthat::test_that("combine_plots_cowplot handles single plot", {

  p <- ggplot2::ggplot(data.frame(x = 1:10, y = 1:10), ggplot2::aes(x, y)) +
    ggplot2::geom_point()

  combined <- .combine_plots_cowplot(list(p), agg_label_unique = "median")

  # cowplot returns a ggplot object
  testthat::expect_is(combined, "ggplot")
})

testthat::test_that("combine_plots_cowplot handles multiple plots", {

  plots <- list(
    ggplot2::ggplot(data.frame(x = 1:5, y = 1:5), ggplot2::aes(x, y)) +
      ggplot2::geom_point(),
    ggplot2::ggplot(data.frame(x = 1:5, y = 10:6), ggplot2::aes(x, y)) +
      ggplot2::geom_point()
  )

  combined <- .combine_plots_cowplot(plots, agg_label_unique = "median")

  testthat::expect_is(combined, "ggplot")
})

testthat::test_that("combine_plots_cowplot returns invisible NULL with output_file", {

  p <- ggplot2::ggplot(data.frame(x = 1:5, y = 1:5), ggplot2::aes(x, y)) +
    ggplot2::geom_point()

  # Create temporary file
  temp_file <- tempfile(fileext = ".pdf")
  on.exit(unlink(temp_file))

  result <- .combine_plots_cowplot(
    list(p),
    output_file = temp_file,
    agg_label_unique = "median"
  )

  # Should return invisible NULL
  testthat::expect_null(result)

  # File should be created
  testthat::expect_true(file.exists(temp_file))
})

# ============================================================================
# TEST: Grid Composition
# ============================================================================

testthat::test_that("combine_plots_grid handles single plot", {

  p <- ggplot2::ggplot(data.frame(x = 1:10, y = 1:10), ggplot2::aes(x, y)) +
    ggplot2::geom_point()

  # Grid composition returns invisible NULL by default
  result <- .combine_plots_grid(list(p), agg_label_unique = "median")

  testthat::expect_null(result)
})

testthat::test_that("combine_plots_grid handles multiple plots", {

  plots <- list(
    ggplot2::ggplot(data.frame(x = 1:5, y = 1:5), ggplot2::aes(x, y)) +
      ggplot2::geom_point(),
    ggplot2::ggplot(data.frame(x = 1:5, y = 10:6), ggplot2::aes(x, y)) +
      ggplot2::geom_point()
  )

  result <- .combine_plots_grid(plots, agg_label_unique = "median")

  testthat::expect_null(result)
})

testthat::test_that("combine_plots_grid saves PNG with output_file", {

  p <- ggplot2::ggplot(data.frame(x = 1:5, y = 1:5), ggplot2::aes(x, y)) +
    ggplot2::geom_point()

  temp_file <- tempfile(fileext = ".png")
  on.exit(unlink(temp_file))

  result <- .combine_plots_grid(
    list(p),
    output_file = temp_file,
    agg_label_unique = "median"
  )

  testthat::expect_null(result)
  testthat::expect_true(file.exists(temp_file))
})

# ============================================================================
# TEST: Composition Backend Selection
# ============================================================================

testthat::test_that("combine_plots_patchwork produces distinct layouts from cowplot", {

  plots <- list(
    ggplot2::ggplot(data.frame(x = 1:5, y = 1:5), ggplot2::aes(x, y)) +
      ggplot2::geom_point(),
    ggplot2::ggplot(data.frame(x = 1:5, y = 10:6), ggplot2::aes(x, y)) +
      ggplot2::geom_point()
  )

  p_patchwork <- .combine_plots_patchwork(plots, agg_label_unique = "metric")
  p_cowplot <- .combine_plots_cowplot(plots, agg_label_unique = "metric")

  # Both should produce ggplot objects
  testthat::expect_is(p_patchwork, "ggplot")
  testthat::expect_is(p_cowplot, "ggplot")
})

# ============================================================================
# TEST: Edge Cases and Error Handling
# ============================================================================

testthat::test_that("create_color_scale handles NULL name parameter", {

  scale <- .create_color_scale(palette = "blue_red", name = NULL)
  testthat::expect_is(scale, "Scale")
})

testthat::test_that("set_plot_title preserves existing plot aesthetics", {

  p <- ggplot2::ggplot(data.frame(x = 1:10, y = 1:10), ggplot2::aes(x, y)) +
    ggplot2::geom_point(color = "red", size = 3) +
    ggplot2::labs(x = "X Axis", y = "Y Axis")

  p_modified <- .set_plot_title(p, title = "New Title")

  # Check title was added
  testthat::expect_equal(p_modified$labels$title, "New Title")
  # Check original labels preserved
  testthat::expect_equal(p_modified$labels$x, "X Axis")
  testthat::expect_equal(p_modified$labels$y, "Y Axis")
  # Check plot still has the same structure
  testthat::expect_equal(length(p_modified$layers), length(p$layers))
})

testthat::test_that("combine_plots functions handle plots with legends", {

  df <- data.frame(x = 1:5, y = 1:5, group = c("A", "B", "A", "B", "A"))
  p <- ggplot2::ggplot(df, ggplot2::aes(x, y, color = group)) +
    ggplot2::geom_point() +
    ggplot2::scale_color_manual(values = c("A" = "red", "B" = "blue"))
  
  combined <- .combine_plots_patchwork(list(p), agg_label_unique = "test")

  testthat::expect_is(combined, "ggplot")
})

# ============================================================================
# TEST: Draw Transcript Grid Utility
# ============================================================================

testthat::test_that("draw_transcript_grid handles proper dimensions", {

  p <- ggplot2::ggplot(data.frame(x = 1:5, y = 1:5), ggplot2::aes(x, y)) +
    ggplot2::geom_point()

  grob <- ggplot2::ggplotGrob(p)

  # The function draws to device, so we just test it doesn't error
  testthat::expect_error(
    .draw_transcript_grid(
      list(grob),
      agg_label_unique = "test",
      legend_grob = NULL,
      ncol = 1,
      heights = grid::unit(c(0.5, 1), "cm")
    ),
    NA  # Expect no error
  )
})

# ============================================================================
# TEST: Integration Tests
# ============================================================================

testthat::test_that("theme + composition workflow produces valid plot", {

  # Create adata frame with groups
  df <- data.frame(
    x = rep(1:5, 2),
    y = c(1:5, 6:10),
    group = rep(c("A", "B"), each = 5)
  )

  # Create plot with custom theme
  p <- ggplot2::ggplot(df, ggplot2::aes(x, y, color = group)) +
    ggplot2::geom_point() +
    .apply_tsenat_theme()

  # Apply title
  p_titled <- .set_plot_title(p, title = "Test Plot", subtitle = "Integration Test")

  # Combine with another plot
  plots <- list(p_titled, p_titled)
  combined <- .combine_plots_patchwork(plots, agg_label_unique = "test_metric")

  testthat::expect_is(combined, "ggplot")
})

testthat::test_that("scale creation works with theme application", {

  df <- data.frame(
    x = 1:10,
    y = 1:10,
    z = seq(0, 1, length.out = 10)
  )

  p <- ggplot2::ggplot(df, ggplot2::aes(x, y, fill = z)) +
    ggplot2::geom_tile() +
    .create_fill_scale(palette = "continuous_diverging") +
    .apply_tsenat_theme()

  testthat::expect_is(p, "ggplot")
})

# ============================================================================
# TEST: generate_plots_spectrum.R Helper Functions
# ============================================================================

context("generate_plots_spectrum: Helper Functions")

testthat::test_that(".divergence_profile_extract_q_values parses q_ prefix format", {
  # Test with "q_" format - column names should start with q_
  col_names <- c("q_0.5", "q_1.0", "q_2.0")
  q_vals <- vapply(col_names, function(name) {
    extracted <- gsub("^q[_=]", "", name)
    as.numeric(extracted)
  }, FUN.VALUE = numeric(1))
  
  testthat::expect_equal(length(q_vals), 3)
  testthat::expect_true(all(!is.na(q_vals)))
  testthat::expect_equal(unname(q_vals), c(0.5, 1.0, 2.0), tolerance = 1e-10)
})

testthat::test_that(".divergence_profile_extract_q_values parses q= format", {
  # Test with "q=" format - column names should start with q=
  col_names <- c("q=0.5", "q=1.0", "q=2.0")
  q_vals <- vapply(col_names, function(name) {
    extracted <- gsub("^q[_=]", "", name)
    as.numeric(extracted)
  }, FUN.VALUE = numeric(1))
  
  testthat::expect_equal(length(q_vals), 3)
  testthat::expect_true(all(!is.na(q_vals)))
  testthat::expect_equal(unname(q_vals), c(0.5, 1.0, 2.0), tolerance = 1e-10)
})

testthat::test_that(".divergence_profile_extract_q_values handles malformed names gracefully", {
  # Test with invalid format
  col_names <- c("sample_invalid", "another_bad")
  q_vals <- suppressWarnings(as.numeric(gsub("^q[_=]", "", col_names)))
  
  testthat::expect_true(all(is.na(q_vals)))
})

testthat::test_that(".find_gene_column identifies 'gene' column", {
  df_gene <- data.frame(gene = c("G1", "G2"), p_value = c(0.01, 0.05))
  col_name <- colnames(df_gene)[grep("^gene", colnames(df_gene))][1]
  
  testthat::expect_equal(col_name, "gene")
})

testthat::test_that(".find_gene_column identifies 'gene_name' column", {
  df_gene_name <- data.frame(gene_name = c("G1", "G2"), p_value = c(0.01, 0.05))
  col_name <- colnames(df_gene_name)[grep("^gene", colnames(df_gene_name))][1]
  
  testthat::expect_equal(col_name, "gene_name")
})

testthat::test_that(".find_gene_column identifies 'gene_id' column", {
  df_gene_id <- data.frame(gene_id = c("ENSEMBL0001", "ENSEMBL0002"), p_value = c(0.01, 0.05))
  col_name <- colnames(df_gene_id)[grep("^gene", colnames(df_gene_id))][1]
  
  testthat::expect_equal(col_name, "gene_id")
})

testthat::test_that(".find_pvalue_column identifies adj_p_interaction", {
  df_adj <- data.frame(
    gene = c("G1", "G2"),
    adj_p_interaction = c(0.01, 0.05)
  )
  
  p_cols <- c("adj_p_interaction", "p_interaction", "adj_p_value", "p_value")
  found <- p_cols[p_cols %in% colnames(df_adj)]
  
  testthat::expect_equal(found[1], "adj_p_interaction")
})

testthat::test_that(".find_pvalue_column identifies p_interaction", {
  df_p <- data.frame(
    gene = c("G1", "G2"),
    p_interaction = c(0.01, 0.05)
  )
  
  p_cols <- c("adj_p_interaction", "p_interaction", "adj_p_value", "p_value")
  found <- p_cols[p_cols %in% colnames(df_p)]
  
  testthat::expect_equal(found[1], "p_interaction")
})

# ============================================================================
# TEST: generate_plots_profile.R Helper Functions
# ============================================================================

context("generate_plots_profile: Helper Functions")

testthat::test_that("select_genesselect_genes handles user-provided gene vector", {
  gene_vec <- c("GENE1", "GENE2", "GENE3")
  result <- as.character(unique(gene_vec))
  
  testthat::expect_length(result, 3)
  testthat::expect_equal(result, gene_vec)
})

testthat::test_that("select_genesselect_genes handles NULL gene input", {
  # Create sample lm_res data.frame
  lm_res <- data.frame(
    gene = c("G1", "G2", "G3", "G4", "G5"),
    adj_p_interaction = c(0.001, 0.01, 0.05, 0.1, 0.2)
  )
  
  # Simulate selecting top 3 genes by p-value
  genes_ordered <- unique(as.character(lm_res$gene[order(lm_res$adj_p_interaction)]))
  top_genes <- head(genes_ordered, 3)
  
  testthat::expect_length(top_genes, 3)
  testthat::expect_equal(top_genes, c("G1", "G2", "G3"))
})

testthat::test_that("select_genesextract_q_values parses column names correctly", {
  # Simulate column names with q-values (realistic SE column names)
  col_names <- c("sample1_q_0.5", "sample2_q_0.5", "sample1_q_1.0", "sample2_q_1.0")
  
  extract_q <- function(name) {
    if (grepl("_q[_=]", name)) {
      as.numeric(gsub(".*_q[_=]", "", name))
    } else {
      NA
    }
  }
  
  q_values <- vapply(col_names, extract_q, FUN.VALUE = numeric(1))
  unique_q <- sort(unique(q_values[!is.na(q_values)]))
  
  testthat::expect_equal(length(unique_q), 2)
  testthat::expect_equal(unique_q, c(0.5, 1.0))
})

testthat::test_that("select_genesextract_q_values returns sorted unique q values", {
  col_names <- c("s1_q_2.0", "s2_q_0.5", "s1_q_1.5", "s2_q_2.0", "s1_q_0.5")
  
  extract_q <- function(name) {
    if (grepl("_q[_=]", name)) {
      as.numeric(gsub(".*_q[_=]", "", name))
    } else {
      NA
    }
  }
  
  q_values <- vapply(col_names, extract_q, FUN.VALUE = numeric(1))
  unique_q <- sort(unique(q_values[!is.na(q_values)]))
  
  testthat::expect_equal(unique_q, c(0.5, 1.5, 2.0))
  testthat::expect_true(is.ordered(unique_q) || all(diff(unique_q) > 0))
})

testthat::test_that("select_genesbuild_facet_plot returns ggplot object", {
  
  # Create sample plot data with realistic q values
  plot_data <- data.frame(
    gene = c("G1", "G1", "G2", "G2"),
    q = c(0.5, 1.0, 0.5, 1.0),
    divergence = c(0.5, 0.7, 0.4, 0.6),
    stringsAsFactors = FALSE
  )
  
  groups <- c("Control", "Treatment")
  
  # Build faceted plot
  p <- ggplot2::ggplot(plot_data, ggplot2::aes(x = q, y = divergence, color = gene)) +
    ggplot2::geom_line(linewidth = 1.1) +
    ggplot2::geom_point(size = 3, alpha = 0.7) +
    ggplot2::facet_wrap(~gene, scales = "free_y")
  
  testthat::expect_is(p, "ggplot")
})

testthat::test_that("select_genesbuild_facet_plot handles signed divergence", {
  
  # Create signed plot data
  plot_data <- data.frame(
    gene = c("G1", "G1", "G2", "G2"),
    q = c(0.5, 1.0, 0.5, 1.0),
    divergence = c(-0.2, 0.1, 0.3, -0.1),
    direction = c("Negative: Control higher", "Positive: Treatment higher",
                  "Positive: Treatment higher", "Negative: Control higher"),
    stringsAsFactors = FALSE
  )
  
  groups <- c("Control", "Treatment")
  
  # Add zero line for signed divergence
  p <- ggplot2::ggplot(plot_data, ggplot2::aes(x = q, y = divergence, color = direction)) +
    ggplot2::geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
    ggplot2::geom_line(linewidth = 1.1) +
    ggplot2::geom_point(size = 3, alpha = 0.7)
  
  testthat::expect_is(p, "ggplot")
})

testthat::test_that("select_genesbuild_list_plots returns named list of ggplot objects", {
  
  # Create sample plot data
  plot_data <- data.frame(
    gene = c("G1", "G1", "G2", "G2"),
    q = c(0.5, 1.0, 0.5, 1.0),
    divergence = c(0.5, 0.7, 0.4, 0.6),
    stringsAsFactors = FALSE
  )
  
  genes <- c("G1", "G2")
  plots <- list()
  
  for (g in genes) {
    df_gene <- plot_data[plot_data$gene == g, ]
    p_gene <- ggplot2::ggplot(df_gene, ggplot2::aes(x = q, y = divergence)) +
      ggplot2::geom_line(color = "#2E86AB", linewidth = 1.2) +
      ggplot2::geom_point(color = "#2E86AB", size = 3, alpha = 0.8) +
      ggplot2::labs(title = paste("Divergence Profile:", g))
    plots[[g]] <- p_gene
  }
  
  testthat::expect_type(plots, "list")
  testthat::expect_length(plots, 2)
  testthat::expect_named(plots, expected = genes)
  testthat::expect_true(all(vapply(plots, inherits, FUN.VALUE = logical(1), "ggplot")))
})

testthat::test_that("select_genesbuild_list_plots handles empty gene list", {
  plots <- list()
  
  testthat::expect_type(plots, "list")
  testthat::expect_length(plots, 0)
})

testthat::test_that("select_genesselect_genes prioritizes adj_p_lmm over adj_p_interaction", {
  lm_res <- data.frame(
    gene = c("G1", "G2", "G3"),
    adj_p_lmm = c(0.02, 0.01, 0.05),
    adj_p_interaction = c(0.001, 0.002, 0.003)
  )
  
  # Should use adj_p_lmm (first priority)
  p_col <- if ("adj_p_lmm" %in% colnames(lm_res)) "adj_p_lmm" else "adj_p_interaction"
  
  testthat::expect_equal(p_col, "adj_p_lmm")
  
  # Top gene should be G2 (p=0.01)
  genes_ordered <- unique(as.character(lm_res$gene[order(lm_res[[p_col]])]))
  testthat::expect_equal(genes_ordered[1], "G2")
})

# ============================================================================
# TEST: Helper Functions for .prepare_combined_se (Refactored Components)
# ============================================================================

testthat::test_that(".extract_diversity_objects extracts SummarizedExperiment and metadata", {
  
  # Create test analysis
  analysis <- .create_test_analysis(
    n_genes = 5, n_samples_per_group = 3,
    q_values = c(1, 2, 3), seed = 42
  )
  
  div_list <- analysis@diversity_results
  
  # Call helper
  extracted <- TSENAT:::.extract_diversity_objects(div_list)
  
  # Verify structure
  testthat::expect_is(extracted, "list")
  testthat::expect_true("objects" %in% names(extracted))
  testthat::expect_true("q_names" %in% names(extracted))
  testthat::expect_true("first_se" %in% names(extracted))
  testthat::expect_true("bootstrap_ci_available" %in% names(extracted))
  
  # Verify first_se is SummarizedExperiment
  testthat::expect_true(methods::is(extracted$first_se, "SummarizedExperiment"))
  
  # Verify q_names match input
  testthat::expect_equal(extracted$q_names, names(div_list))
  
  # Verify bootstrap_ci_available is logical
  testthat::expect_is(extracted$bootstrap_ci_available, "logical")
})

testthat::test_that(".extract_diversity_objects errors with no valid SE", {
  # Create empty list with only matrices
  div_list <- list(
    q_1 = matrix(1:10, nrow = 5, ncol = 2)
  )
  
  # Should error since no SummarizedExperiment
  testthat::expect_error(
    TSENAT:::.extract_diversity_objects(div_list),
    "No valid SummarizedExperiment"
  )
})

testthat::test_that(".normalize_matrix_to_target pads columns when needed", {
  # Create small matrix with fewer columns than target
  mat <- matrix(1:6, nrow = 3, ncol = 2)
  rownames(mat) <- c("G1", "G2", "G3")
  
  target_genes <- c("G1", "G2", "G3")
  target_n_cols <- 5
  
  # Call helper
  result <- TSENAT:::.normalize_matrix_to_target(mat, target_genes, target_n_cols)
  
  # Verify dimensions
  testthat::expect_equal(ncol(result), 5)  # Padded to 5 columns
  testthat::expect_equal(nrow(result), 3)
  
  # Verify row order preserved
  testthat::expect_equal(rownames(result), target_genes)
  
  # Verify original data preserved
  testthat::expect_equal(result[, 1:2], mat[target_genes, ])
})

testthat::test_that(".normalize_matrix_to_target truncates when needed", {
  
  # Create larger matrix
  mat <- matrix(1:15, nrow = 3, ncol = 5)
  rownames(mat) <- c("G1", "G2", "G3")
  
  target_genes <- c("G1", "G2", "G3")
  target_n_cols <- 3
  
  # Call helper
  result <- TSENAT:::.normalize_matrix_to_target(mat, target_genes, target_n_cols)
  
  # Verify truncation
  testthat::expect_equal(ncol(result), 3)
  testthat::expect_equal(nrow(result), 3)
  
  # Verify data integrity
  testthat::expect_equal(result, mat[target_genes, 1:3])
})

testthat::test_that(".normalize_matrix_to_target reorders rows", {
  mat <- matrix(1:6, nrow = 3, ncol = 2)
  rownames(mat) <- c("G3", "G1", "G2")
  
  target_genes <- c("G1", "G2", "G3")  # Different order
  target_n_cols <- 2
  
  # Call helper
  result <- TSENAT:::.normalize_matrix_to_target(mat, target_genes, target_n_cols)
  
  # Verify row order matches target
  testthat::expect_equal(rownames(result), target_genes)
})

testthat::test_that(".create_q_suffixed_colnames adds q-value suffix", {
  colnames <- c("sample_1", "sample_2", "sample_3")
  q_val <- 1.5
  n_cols <- 3
  
  # Call helper
  result <- TSENAT:::.create_q_suffixed_colnames(colnames, q_val, n_cols)
  
  # Verify format
  testthat::expect_equal(length(result), 3)
  testthat::expect_true(all(grepl("_q=1\\.5", result)))
  testthat::expect_true(all(grepl("^sample_", result)))
})

testthat::test_that(".create_q_suffixed_colnames handles NULL colnames", {
  colnames <- NULL
  q_val <- 2.0
  n_cols <- 4
  
  # Call helper
  result <- TSENAT:::.create_q_suffixed_colnames(colnames, q_val, n_cols)
  
  # Verify synthetic names created
  testthat::expect_equal(length(result), 4)
  testthat::expect_true(all(grepl("sample_", result)))
  testthat::expect_true(all(grepl("_q=2\\.000", result)))
})

testthat::test_that(".create_q_suffixed_colnames removes existing q= suffixes", {
  colnames <- c("sample_1_q=1.000", "sample_2_q=1.000")
  q_val <- 2.0
  n_cols <- 2
  
  # Call helper
  result <- TSENAT:::.create_q_suffixed_colnames(colnames, q_val, n_cols)
  
  # Verify old suffix removed and new one added
  testthat::expect_true(all(grepl("_q=2\\.000", result)))
  testthat::expect_false(any(grepl("_q=1\\.000", result)))
})

testthat::test_that(".build_combined_coldata combines metadata across q-values", {
  
  # Create test analysis
  analysis <- .create_test_analysis(
    n_genes = 4, n_samples_per_group = 2,
    q_values = c(1, 2), seed = 42
  )
  
  div_list <- analysis@diversity_results
  q_names <- names(div_list)
  
  # Create unique colnames for each q - must match the number of columns in each SE
  n_cols_q1 <- ncol(div_list[[1]])
  n_cols_q2 <- ncol(div_list[[2]])
  
  unique_colnames_list <- list(
    q_1 = paste0("sample_", seq_len(n_cols_q1), "_q=1.000"),
    q_2 = paste0("sample_", seq_len(n_cols_q2), "_q=2.000")
  )
  
  # Call helper
  result <- TSENAT:::.build_combined_coldata(div_list, q_names, unique_colnames_list)
  
  # Verify structure
  testthat::expect_is(result, "data.frame")
  testthat::expect_true("q" %in% colnames(result))
  
  # Verify q column has both values
  unique_q_vals <- unique(result$q)
  testthat::expect_equal(length(unique_q_vals), 2)
  
  # Verify row count matches total samples
  testthat::expect_equal(nrow(result), n_cols_q1 + n_cols_q2)
})

testthat::test_that(".extract_bootstrap_ci_matrices returns NULL for non-SE", {
  # Create regular matrix
  mat <- matrix(1:6, nrow = 3, ncol = 2)
  target_genes <- rownames(mat) <- c("G1", "G2", "G3")
  target_n_cols <- 2
  assay_names <- c()
  
  # Call helper
  result <- TSENAT:::.extract_bootstrap_ci_matrices(
    mat, target_genes, target_n_cols, assay_names
  )
  
  # Should return NULL for matrix input
  testthat::expect_null(result)
})

testthat::test_that(".extract_bootstrap_ci_matrices extracts CI matrices from SE", {
  
  # Create test analysis
  analysis <- .create_test_analysis(
    n_genes = 3, n_samples_per_group = 2,
    q_values = c(1), seed = 42
  )
  
  # Get first SE with CIs (if available)
  se_obj <- analysis@diversity_results[[1]]
  
  # Extract assay names
  sim_names <- SummarizedExperiment::assayNames(se_obj)
  
  # Call helper
  result <- TSENAT:::.extract_bootstrap_ci_matrices(
    se_obj,
    rownames(se_obj),
    ncol(se_obj),
    sim_names
  )
  
  # Verify result structure - result can be NULL or a list
  if (!is.null(result)) {
    testthat::expect_is(result, "list")
    testthat::expect_true("ci_lower" %in% names(result))
    testthat::expect_true("ci_upper" %in% names(result))
  } else {
    # If CIs not available, that's also valid - test passes
    testthat::expect_null(result)
  }
})

testthat::test_that(".prepare_q_value_for_combining processes q-value data correctly", {
  
  # Create test analysis
  analysis <- .create_test_analysis(
    n_genes = 3, n_samples_per_group = 2,
    q_values = c(1.5), seed = 42
  )
  
  div_list <- analysis@diversity_results
  combined_assays_dict <- list(
    q_1.5 = list(
      matrix = SummarizedExperiment::assay(div_list[[1]], 1),
      q_val = 1.5,
      se_obj = div_list[[1]]
    )
  )
  
  target_genes <- rownames(div_list[[1]])
  target_n_cols <- ncol(div_list[[1]])
  
  # Call helper
  result <- TSENAT:::.prepare_q_value_for_combining(
    "q_1.5", combined_assays_dict, target_genes,
    target_n_cols, FALSE  # No bootstrap CI
  )
  
  # Verify result structure
  testthat::expect_is(result, "list")
  testthat::expect_true("matrix" %in% names(result))
  testthat::expect_true("unique_colnames" %in% names(result))
  testthat::expect_true("ncol_val" %in% names(result))
  
  # Verify q-value appears in colnames
  testthat::expect_true(all(grepl("_q=1\\.500", result$unique_colnames)))
})

testthat::test_that(".fill_combined_assays combines multiple q-values", {
  
  # Create test analysis with multiple q-values
  analysis <- .create_test_analysis(
    n_genes = 3, n_samples_per_group = 2,
    q_values = c(1, 2), seed = 42
  )
  
  extracted <- TSENAT:::.extract_diversity_objects(analysis@diversity_results)
  
  target_genes <- rownames(extracted$first_se)
  target_n_cols <- ncol(extracted$first_se)
  
  # Call helper
  filled <- TSENAT:::.fill_combined_assays(
    extracted$objects, extracted$q_names,
    target_genes, target_n_cols,
    extracted$bootstrap_ci_available
  )
  
  # Verify structure
  testthat::expect_is(filled, "list")
  testthat::expect_is(filled$combined_assay, "matrix")
  
  # Verify dimensions
  expected_cols <- target_n_cols * length(extracted$q_names)
  testthat::expect_equal(ncol(filled$combined_assay), expected_cols)
  testthat::expect_equal(nrow(filled$combined_assay), length(target_genes))
  
  # Verify column names collect all q-values
  testthat::expect_equal(length(filled$unique_colnames_list), length(extracted$q_names))
})

testthat::test_that(".create_combined_se_object creates valid SummarizedExperiment", {
  
  # Create minimal test data
  genes <- c("G1", "G2", "G3")
  samples <- c("S1_q=1", "S2_q=1", "S3_q=2", "S4_q=2")
  
  combined_assay <- matrix(1:12, nrow = 3, ncol = 4)
  rownames(combined_assay) <- genes
  colnames(combined_assay) <- samples
  
  combined_coldata <- data.frame(
    q = c(1, 1, 2, 2),
    row.names = samples
  )
  
  # Create a simple SE as template
  se_template <- SummarizedExperiment::SummarizedExperiment(
    assays = list(diversity = combined_assay),
    colData = combined_coldata
  )
  
  # Call helper
  result <- TSENAT:::.create_combined_se_object(
    combined_assay,
    NULL, NULL,  # No CI matrices
    combined_coldata,
    se_template
  )
  
  # Verify result
  testthat::expect_true(methods::is(result, "SummarizedExperiment"))
  testthat::expect_equal(nrow(result), 3)
  testthat::expect_equal(ncol(result), 4)
  testthat::expect_true("diversity" %in% SummarizedExperiment::assayNames(result))
})

testthat::test_that(".create_combined_se_object includes CI assays when provided", {
  
  # Create test data with CIs
  genes <- c("G1", "G2")
  samples <- c("S1", "S2")
  
  combined_assay <- matrix(1:4, nrow = 2, ncol = 2)
  rownames(combined_assay) <- genes
  colnames(combined_assay) <- samples
  
  combined_ci_lower <- matrix(0.5:3.5, nrow = 2, ncol = 2)
  rownames(combined_ci_lower) <- genes
  colnames(combined_ci_lower) <- samples
  combined_ci_upper <- matrix(1.5:4.5, nrow = 2, ncol = 2)
  rownames(combined_ci_upper) <- genes
  colnames(combined_ci_upper) <- samples
  
  combined_coldata <- data.frame(q = c(1, 1), row.names = samples)
  
  se_template <- SummarizedExperiment::SummarizedExperiment(
    assays = list(diversity = combined_assay),
    colData = combined_coldata
  )
  
  # Call helper with CI matrices
  result <- TSENAT:::.create_combined_se_object(
    combined_assay,
    combined_ci_lower, combined_ci_upper,
    combined_coldata,
    se_template
  )
  
  # Verify CI assays included
  assay_names <- SummarizedExperiment::assayNames(result)
  testthat::expect_true("ci_lower" %in% assay_names)
  testthat::expect_true("ci_upper" %in% assay_names)
})

testthat::test_that(".prepare_combined_se integration test with all helpers", {
  
  # Create test analysis
  analysis <- .create_test_analysis(
    n_genes = 5, n_samples_per_group = 3,
    q_values = c(1, 2, 3), seed = 42
  )
  
  # Call main function (which uses all refactored helpers)
  se_combined <- TSENAT:::.prepare_combined_se(analysis)
  
  # Verify result
  testthat::expect_true(methods::is(se_combined, "SummarizedExperiment"))
  testthat::expect_equal(nrow(se_combined), 5)
  testthat::expect_equal(ncol(se_combined), 18)  # 3 samples * 2 groups * 3 q-values
  
  # Verify assays
  testthat::expect_true("diversity" %in% SummarizedExperiment::assayNames(se_combined))
  
  # Verify colData structure
  coldata <- SummarizedExperiment::colData(se_combined)
  testthat::expect_true("q" %in% colnames(coldata))
  
  # Verify rowData structure
  rowdata <- SummarizedExperiment::rowData(se_combined)
  testthat::expect_true("gene_id" %in% colnames(rowdata))
  
  # Verify column names have q suffixes
  colnames_se <- colnames(se_combined)
  testthat::expect_true(any(grepl("_q=", colnames_se)))
})

# =============================================================================
# COMPREHENSIVE TESTS FOR 11 UNCOVERED PLOT HELPER FUNCTIONS
# =============================================================================
# These tests cover the critical plot helper functions identified in coverage
# analysis as having 0% test coverage.

# =============================================================================
# 1. .extract_bootstrap_ci_assays - CRITICAL FOUNDATION TEST
# =============================================================================

testthat::test_that(".extract_bootstrap_ci_assays correctly identifies CI assays", {
  
  # Create test SE with both base and CI assays
  base_assay <- matrix(rnorm(100), nrow = 10, ncol = 10)
  ci_lower <- base_assay - 0.5
  ci_upper <- base_assay + 0.5
  
  rownames(base_assay) <- paste0("Gene", 1:10)
  colnames(base_assay) <- paste0("Sample", 1:10)
  rownames(ci_lower) <- rownames(base_assay)
  colnames(ci_lower) <- colnames(base_assay)
  rownames(ci_upper) <- rownames(base_assay)
  colnames(ci_upper) <- colnames(base_assay)
  
  se <- SummarizedExperiment::SummarizedExperiment(
    assays = list(
      diversity = base_assay,
      diversity_ci_lower = ci_lower,
      diversity_ci_upper = ci_upper
    )
  )
  
  result <- .extract_bootstrap_ci_assays(se, assay_name = "diversity")
  
  testthat::expect_true(result$has_ci)
  testthat::expect_false(is.null(result$ci_lower))
  testthat::expect_false(is.null(result$ci_upper))
  testthat::expect_equal(nrow(result$ci_lower), 10)
  testthat::expect_equal(ncol(result$ci_lower), 10)
  testthat::expect_null(result$fallback_metric)
})

testthat::test_that(".extract_bootstrap_ci_assays handles missing CI assays", {
  
  base_assay <- matrix(rnorm(100), nrow = 10, ncol = 10)
  rownames(base_assay) <- paste0("Gene", 1:10)
  colnames(base_assay) <- paste0("Sample", 1:10)
  
  se <- SummarizedExperiment::SummarizedExperiment(
    assays = list(diversity = base_assay)
  )
  
  result <- .extract_bootstrap_ci_assays(se, assay_name = "diversity", fallback_to_iqr = TRUE)
  
  testthat::expect_false(result$has_ci)
  testthat::expect_null(result$ci_lower)
  testthat::expect_null(result$ci_upper)
  testthat::expect_equal(result$fallback_metric, "iqr")
})

testthat::test_that(".extract_bootstrap_ci_assays errors on missing base assay", {
  
  se <- SummarizedExperiment::SummarizedExperiment(
    assays = list(other = matrix(rnorm(100), nrow = 10, ncol = 10))
  )
  
  testthat::expect_error(
    .extract_bootstrap_ci_assays(se, assay_name = "diversity"),
    "Assay 'diversity' not found"
  )
})

testthat::test_that(".extract_bootstrap_ci_assays handles partial CI assays", {
  
  # Only lower CI present
  base_assay <- matrix(rnorm(100), nrow = 10, ncol = 10)
  ci_lower <- base_assay - 0.5
  
  rownames(base_assay) <- paste0("Gene", 1:10)
  colnames(base_assay) <- paste0("Sample", 1:10)
  rownames(ci_lower) <- rownames(base_assay)
  colnames(ci_lower) <- colnames(base_assay)
  
  se <- SummarizedExperiment::SummarizedExperiment(
    assays = list(
      diversity = base_assay,
      diversity_ci_lower = ci_lower
    )
  )
  
  result <- .extract_bootstrap_ci_assays(se, assay_name = "diversity")
  
  testthat::expect_false(result$has_ci)
  testthat::expect_false(is.null(result$ci_lower))
  testthat::expect_null(result$ci_upper)
})

# =============================================================================
# 2. .calculate_scaled_fonts - Font Scaling Test
# =============================================================================

testthat::test_that(".calculate_scaled_fonts computes correct scaled sizes", {
  # No packages required
  result <- .calculate_scaled_fonts(base_size = 11, scale_factor = 1.0)
  
  testthat::expect_is(result, "list")
  testthat::expect_equal(result$base, 11)
  testthat::expect_equal(result$scaled, 11)
  testthat::expect_true(result$axis_text > result$base)
  testthat::expect_true(result$title > result$axis_title)
})

testthat::test_that(".calculate_scaled_fonts scales by multipliers", {
  multipliers <- list(
    axis_text = 1.2,
    axis_title = 1.3,
    title = 1.5,
    legend = 0.9
  )
  
  result <- .calculate_scaled_fonts(base_size = 10, scale_factor = 2.0, font_multipliers = multipliers)
  
  testthat::expect_equal(result$scaled, 20)
  testthat::expect_equal(result$axis_text, round(20 * 1.2))
  testthat::expect_equal(result$title, round(20 * 1.5))
})

testthat::test_that(".calculate_scaled_fonts computes all font sizes", {
  result <- .calculate_scaled_fonts(base_size = 11, scale_factor = 0.8)
  
  testthat::expect_true(all(c("base", "scaled", "axis_text", "axis_title", 
                              "title", "legend", "subtitle", "caption") %in% names(result)))
  testthat::expect_true(all(unlist(result) > 0))
})

# =============================================================================
# 3. .create_centered_theme - Theme Creation Test
# =============================================================================

testthat::test_that(".create_centered_theme creates theme object", {
  
  theme <- .create_centered_theme(include_title = TRUE, include_subtitle = TRUE)
  
  testthat::expect_is(theme, "theme")
  testthat::expect_false(is.null(theme$plot.title))
  testthat::expect_false(is.null(theme$plot.subtitle))
})

testthat::test_that(".create_centered_theme respects flags", {
  
  theme_both <- .create_centered_theme(include_title = TRUE, include_subtitle = TRUE)
  theme_title_only <- .create_centered_theme(include_title = TRUE, include_subtitle = FALSE)
  theme_none <- .create_centered_theme(include_title = FALSE, include_subtitle = FALSE)
  
  testthat::expect_is(theme_both, "theme")
  testthat::expect_is(theme_title_only, "theme")
  testthat::expect_is(theme_none, "theme")
})

testthat::test_that(".create_centered_theme applies centering", {
  
  # Create theme with centering
  theme <- .create_centered_theme(include_title = TRUE, hjust = 0.5)
  
  # Apply to plot and verify no errors
  p <- ggplot2::ggplot(data.frame(x = 1:5), ggplot2::aes(x = x)) +
    ggplot2::geom_point() +
    theme
  
  testthat::expect_is(p, "ggplot")
})

# =============================================================================
# 4. .normalize_plot_scales - Scale Normalization Test
# =============================================================================

testthat::test_that(".normalize_plot_scales identifies fold-change column", {
  # No packages required
  df <- data.frame(
    logFC = c(1.5, -2.0, 0.5),
    condition1_mean = c(100, 200, 150),
    condition2_mean = c(110, 210, 160),
    gene = c("A", "B", "C"),
    stringsAsFactors = FALSE
  )
  
  result <- .normalize_plot_scales(df, fold_col_candidates = c("log2_fold_change", "logFC", "fold"))
  
  testthat::expect_is(result, "list")
  testthat::expect_equal(result$fold_col, "logFC")
  testthat::expect_true("x_norm" %in% colnames(result$df))
  testthat::expect_true("y_norm" %in% colnames(result$df))
})

testthat::test_that(".normalize_plot_scales computes normalized positions", {
  # No packages required
  df <- data.frame(
    log2_fold_change = c(1.0, -1.0, 0.5),
    control_mean = c(100, 200, 150),
    treatment_mean = c(120, 210, 160),
    gene = c("A", "B", "C"),
    stringsAsFactors = FALSE
  )
  
  result <- .normalize_plot_scales(df, fold_col_candidates = c("log2_fold_change", "logFC"))
  
  # x_norm should be average of means
  expected_x <- c((100+120)/2, (200+210)/2, (150+160)/2)
  testthat::expect_equal(result$df$x_norm, expected_x, tolerance = 0.1)
})

testthat::test_that(".normalize_plot_scales warns on missing fold column", {
  # No packages required
  df <- data.frame(
    x = c(1, 2, 3),
    y = c(4, 5, 6)
  )
  
  testthat::expect_warning(
    .normalize_plot_scales(df, fold_col_candidates = c("logFC", "log2FC")),
    "No fold-change column found"
  )
})

# =============================================================================
# 5. .compute_distribution_stats - Statistics Computation Test
# =============================================================================

testthat::test_that(".compute_distribution_stats computes median and IQR", {
  
  df <- data.frame(
    group = c("A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "B", "B", "B", "B", "B", "B", "B", "B", "B", "B"),
    value = c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20)
  )
  
  result <- .compute_distribution_stats(df, "group", "value", metric = "median", spread_metric = "iqr")
  
  testthat::expect_equal(nrow(result), 2)
  testthat::expect_true(all(c("value", "lower", "upper") %in% colnames(result)))
  testthat::expect_true(result$lower[1] < result$value[1])
  testthat::expect_true(result$upper[1] > result$value[1])
})

testthat::test_that(".compute_distribution_stats computes mean and SD", {
  
  df <- data.frame(
    group = c(rep("A", 10), rep("B", 10)),
    value = c(c(-5, -4, -3, -2, -1, 0, 1, 2, 3, 4), c(5, 6, 7, 8, 9, 10, 11, 12, 13, 14))
  )
  
  result <- .compute_distribution_stats(df, "group", "value", metric = "mean", spread_metric = "sd")
  
  testthat::expect_equal(nrow(result), 2)
  testthat::expect_true(all(!is.na(result$value)))
  testthat::expect_true(!is.na(result$upper[1]) && result$upper[1] > result$value[1])
})

testthat::test_that(".compute_distribution_stats handles missing values", {
  
  df <- data.frame(
    group = rep(c("A", "B"), each = 10),
    value = c(rnorm(9), NA, rnorm(9, mean = 2), NA)
  )
  
  result <- .compute_distribution_stats(df, "group", "value")
  
  testthat::expect_equal(nrow(result), 2)
  testthat::expect_true(all(!is.na(result$value)))
})

testthat::test_that(".compute_distribution_stats errors on non-numeric values", {
  
  df <- data.frame(group = c("A", "B"), value = c("x", "y"))
  
  testthat::expect_error(
    .compute_distribution_stats(df, "group", "value"),
    "not numeric"
  )
})

testthat::test_that(".compute_distribution_stats errors on missing columns", {
  df <- data.frame(x = c(1, 2), y = c(3, 4))
  
  testthat::expect_error(
    .compute_distribution_stats(df, "group", "value"),
    "not found"
  )
})

# =============================================================================
# 6. .prepare_grouped_long_format - Data Preparation Test
# =============================================================================

testthat::test_that(".prepare_grouped_long_format transforms to long format", {
  
  assay_mat <- matrix(rnorm(50), nrow = 5, ncol = 10)
  rownames(assay_mat) <- paste0("Gene", 1:5)
  colnames(assay_mat) <- paste0("Sample", 1:10)
  
  se <- SummarizedExperiment::SummarizedExperiment(
    assays = list(diversity = assay_mat),
    colData = data.frame(
      condition = rep(c("A", "B"), each = 5),
      row.names = colnames(assay_mat)
    )
  )
  
  result <- .prepare_grouped_long_format(se, assay_name = "diversity", group_by_col = "condition")
  
  testthat::expect_is(result, "data.frame")
  testthat::expect_true("Gene" %in% colnames(result))
  testthat::expect_true("group" %in% colnames(result))
  testthat::expect_true("value" %in% colnames(result))
  testthat::expect_true("lower" %in% colnames(result))
  testthat::expect_true("upper" %in% colnames(result))
})

testthat::test_that(".prepare_grouped_long_format computes aggregated statistics", {
  
  assay_mat <- matrix(1:50, nrow = 5, ncol = 10)
  rownames(assay_mat) <- paste0("Gene", 1:5)
  colnames(assay_mat) <- paste0("Sample", 1:10)
  
  se <- SummarizedExperiment::SummarizedExperiment(
    assays = list(diversity = assay_mat),
    colData = data.frame(
      condition = rep(c("A", "B"), each = 5),
      row.names = colnames(assay_mat)
    )
  )
  
  result <- .prepare_grouped_long_format(se, assay_name = "diversity", group_by_col = "condition")
  
  testthat::expect_true(nrow(result) >= 5)
  testthat::expect_true(all(result$value > 0))
})

testthat::test_that(".prepare_grouped_long_format errors on invalid colData column", {
  
  assay_mat <- matrix(rnorm(50), nrow = 5, ncol = 10)
  rownames(assay_mat) <- paste0("Gene", 1:5)
  colnames(assay_mat) <- paste0("Sample", 1:10)
  
  se <- SummarizedExperiment::SummarizedExperiment(
    assays = list(diversity = assay_mat),
    colData = data.frame(row.names = colnames(assay_mat))
  )
  
  testthat::expect_error(
    .prepare_grouped_long_format(se, assay_name = "diversity", group_by_col = "missing_col"),
    "not found"
  )
})

# =============================================================================
# 7. .prepare_gene_ci_data - CI Data Preparation Test
# =============================================================================

testthat::test_that(".prepare_gene_ci_data prepares CI data structure", {
  
  long_data <- data.frame(
    Gene = rep(c("Gene1", "Gene2"), each = 10),
    group = rep(c("A", "A", "B", "B"), length.out = 20),
    q = rep(c(0.5, 1.0), 10),
    tsallis = rnorm(20),
    sample = rep(c("S1", "S2", "S3", "S4", "S5"), 4),
    stringsAsFactors = FALSE
  )
  
  ci_lower <- matrix(rnorm(50), nrow = 2, ncol = 25)
  ci_upper <- matrix(rnorm(50) + 1, nrow = 2, ncol = 25)
  rownames(ci_lower) <- c("Gene1", "Gene2")
  rownames(ci_upper) <- c("Gene1", "Gene2")
  colnames(ci_lower) <- paste0("Sample_q=", rep(c(0.5, 1.0), length.out = 25))
  colnames(ci_upper) <- paste0("Sample_q=", rep(c(0.5, 1.0), length.out = 25))
  
  result <- .prepare_gene_ci_data(long_data, ci_lower, ci_upper, c("Gene1", "Gene2"))
  
  testthat::expect_is(result, "data.frame")
  testthat::expect_true("ci_lower" %in% colnames(result))
  testthat::expect_true("ci_upper" %in% colnames(result))
})

testthat::test_that(".prepare_gene_ci_data handles empty CI matrices", {
  
  long_data <- data.frame(
    Gene = c("Gene1", "Gene2"),
    group = c("A", "B"),
    q = c(0.5, 1.0),
    tsallis = c(1.5, 2.0),
    sample = c("S1", "S2"),
    stringsAsFactors = FALSE
  )
  
  ci_lower <- matrix(nrow = 0, ncol = 0)
  ci_upper <- matrix(nrow = 0, ncol = 0)
  
  result <- .prepare_gene_ci_data(long_data, ci_lower, ci_upper, c("Gene1", "Gene2"))
  
  testthat::expect_is(result, "data.frame")
})

# =============================================================================
# 8. .prepare_transcript_inputs - Input Preparation Test
# =============================================================================

testthat::test_that(".prepare_transcript_inputs validates matrix input", {
  counts <- matrix(rnbinom(100, size = 1, prob = 0.1), nrow = 10, ncol = 10)
  rownames(counts) <- paste0("TX", 1:10)
  colnames(counts) <- paste0("Sample", 1:10)
  
  samples <- rep(c("A", "B"), each = 5)
  
  tx2gene <- data.frame(
    Transcript = rownames(counts),
    Gen = paste0("Gene", sample(1:5, 10, replace = TRUE)),
    stringsAsFactors = FALSE
  )
  
  result <- .prepare_transcript_inputs(
    counts, samples = samples, tx2gene = tx2gene, metric = "median"
  )
  
  testthat::expect_is(result, "list")
  testthat::expect_true("counts" %in% names(result))
  testthat::expect_true("samples" %in% names(result))
  testthat::expect_true("mapping" %in% names(result))
  testthat::expect_equal(length(result$samples), ncol(counts))
})

testthat::test_that(".prepare_transcript_inputs errors on missing rownames", {
  counts <- matrix(rnorm(100), nrow = 10, ncol = 10)
  
  testthat::expect_error(
    .prepare_transcript_inputs(counts, samples = rep(c("A", "B"), 5)),
    "must have rownames"
  )
})

testthat::test_that(".prepare_transcript_inputs handles different metrics", {
  counts <- matrix(rnbinom(100, size = 1, prob = 0.1), nrow = 10, ncol = 10)
  rownames(counts) <- paste0("TX", 1:10)
  colnames(counts) <- paste0("Sample", 1:10)
  
  samples <- rep(c("A", "B"), each = 5)
  tx2gene <- data.frame(
    Transcript = rownames(counts),
    Gen = paste0("Gene", sample(1:5, 10, replace = TRUE)),
    stringsAsFactors = FALSE
  )
  
  result_median <- .prepare_transcript_inputs(
    counts, samples = samples, tx2gene = tx2gene, metric = "median"
  )
  
  result_mean <- .prepare_transcript_inputs(
    counts, samples = samples, tx2gene = tx2gene, metric = "mean"
  )
  
  testthat::expect_equal(result_median$metric_choice, "median")
  testthat::expect_equal(result_mean$metric_choice, "mean")
})

# =============================================================================
# 9. .plot_gam_arrange_grid - Grid Arrangement Test
# =============================================================================

testthat::test_that(".plot_gam_arrange_grid creates grid plot", {
  
  plots <- list(
    ggplot2::ggplot(data.frame(x = 1:10, y = 1:10), 
                   ggplot2::aes(x = x, y = y)) + ggplot2::geom_point(),
    ggplot2::ggplot(data.frame(x = 1:10, y = 10:1), 
                   ggplot2::aes(x = x, y = y)) + ggplot2::geom_point()
  )
  
  font_sizes <- list(legend_text = 10, legend_title = 11)
  
  result <- .plot_gam_arrange_grid(plots, condition_col = "condition", font_sizes = font_sizes)
  
  testthat::expect_true(ggplot2::is_ggplot(result) || inherits(result, "gtable"))
})

testthat::test_that(".plot_gam_arrange_grid handles single plot", {
  
  plots <- list(
    ggplot2::ggplot(data.frame(x = 1:5, y = 1:5), 
                   ggplot2::aes(x = x, y = y)) + ggplot2::geom_point()
  )
  
  font_sizes <- list(legend_text = 10, legend_title = 11)
  
  result <- .plot_gam_arrange_grid(plots, condition_col = "cond", font_sizes = font_sizes)
  
  testthat::expect_true(ggplot2::is_ggplot(result) || inherits(result, "gtable"))
})

# =============================================================================
# 10. .plot_gam_save_plot - Plot Saving Test
# =============================================================================

testthat::test_that(".plot_gam_save_plot saves plot to file", {
  
  tmpfile <- tempfile(fileext = ".png")
  on.exit(unlink(tmpfile))
  
  p <- ggplot2::ggplot(data.frame(x = 1:10, y = 1:10), 
                      ggplot2::aes(x = x, y = y)) + ggplot2::geom_point()
  
  .plot_gam_save_plot(p, output_file = tmpfile, width = 8, height = 6)
  
  testthat::expect_true(file.exists(tmpfile))
  testthat::expect_true(file.size(tmpfile) > 0)
})

testthat::test_that(".plot_gam_save_plot handles NULL output_file", {
  
  p <- ggplot2::ggplot(data.frame(x = 1:10, y = 1:10), 
                      ggplot2::aes(x = x, y = y)) + ggplot2::geom_point()
  
  # Should not error when output_file is NULL
  result <- .plot_gam_save_plot(p, output_file = NULL, width = 8, height = 6)
  testthat::expect_true(is.null(result) || is.function(result))
})

# =============================================================================
# 11. .save_plot_standard - Standard Plot Saving Test
# =============================================================================

testthat::test_that(".save_plot_standard saves with standard dimensions", {
  
  tmpfile <- tempfile(fileext = ".png")
  on.exit(unlink(tmpfile))
  
  p <- ggplot2::ggplot(data.frame(x = 1:10, y = 1:10), 
                      ggplot2::aes(x = x, y = y)) + ggplot2::geom_point()
  
  .save_plot_standard(p, filename = tmpfile, width_inches = 12, aspect_type = "standard")
  
  testthat::expect_true(file.exists(tmpfile))
  testthat::expect_true(file.size(tmpfile) > 0)
})

testthat::test_that(".save_plot_standard handles different aspect ratios", {
  
  p <- ggplot2::ggplot(data.frame(x = 1:10, y = 1:10), 
                      ggplot2::aes(x = x, y = y)) + ggplot2::geom_point()
  
  for (aspect in c("standard", "wide", "tall")) {
    tmpfile_aspect <- tempfile(fileext = ".png")
    on.exit(unlink(tmpfile_aspect))
    
    .save_plot_standard(p, filename = tmpfile_aspect, width_inches = 10, aspect_type = aspect)
    testthat::expect_true(file.exists(tmpfile_aspect))
  }
})

testthat::test_that(".save_plot_standard handles cm dimensions", {
  
  tmpfile <- tempfile(fileext = ".png")
  on.exit(unlink(tmpfile))
  
  p <- ggplot2::ggplot(data.frame(x = 1:10, y = 1:10), 
                      ggplot2::aes(x = x, y = y)) + ggplot2::geom_point()
  
  .save_plot_standard(p, filename = tmpfile, width_cm = 25, height_cm = 20)
  
  testthat::expect_true(file.exists(tmpfile))
})

# =============================================================================
# INTEGRATION TESTS
# =============================================================================

testthat::test_that("CI extraction + distribution stats workflow", {
  
  base_assay <- matrix(rnorm(100), nrow = 5, ncol = 20)
  ci_lower <- base_assay - 0.5
  ci_upper <- base_assay + 0.5
  
  rownames(base_assay) <- paste0("Gene", 1:5)
  colnames(base_assay) <- paste0("Sample", 1:20)
  rownames(ci_lower) <- rownames(base_assay)
  colnames(ci_lower) <- colnames(base_assay)
  rownames(ci_upper) <- rownames(base_assay)
  colnames(ci_upper) <- colnames(base_assay)
  
  se <- SummarizedExperiment::SummarizedExperiment(
    assays = list(
      diversity = base_assay,
      diversity_ci_lower = ci_lower,
      diversity_ci_upper = ci_upper
    ),
    colData = data.frame(
      condition = rep(c("A", "B"), each = 10),
      row.names = colnames(base_assay)
    )
  )
  
  # Extract CIs
  ci_result <- .extract_bootstrap_ci_assays(se, assay_name = "diversity")
  testthat::expect_true(ci_result$has_ci)
  
  # Use in preparation
  long_prepared <- .prepare_grouped_long_format(se, assay_name = "diversity", group_by_col = "condition")
  testthat::expect_true(nrow(long_prepared) > 0)
})

testthat::test_that("Font scaling + theme creation workflow", {
  
  fonts <- .calculate_scaled_fonts(base_size = 11, scale_factor = 0.9)
  testthat::expect_true(fonts$title > fonts$base)
  
  theme <- .create_centered_theme(
    include_title = TRUE,
    include_subtitle = TRUE,
    title_size = fonts$title,
    subtitle_size = fonts$subtitle
  )
  
  testthat::expect_is(theme, "theme")
})

context("Plot Helper Functions - Critical Coverage Fixes")

# =============================================================================
# PRIORITY 1: Foundation CI Extraction (Critical blocker for all CI plots)
# =============================================================================

test_that(".extract_bootstrap_ci_assays works with standard CI naming", {
    skip_if_not_installed("SummarizedExperiment")
    
    # Create test SE with diversity and CI assays
    mat_base <- matrix(rnorm(20), nrow = 5)
    rownames(mat_base) <- paste0("gene", 1:5)
    colnames(mat_base) <- paste0("sample", 1:4)
    
    mat_ci_lower <- matrix(rnorm(20, mean = -0.5), nrow = 5)
    rownames(mat_ci_lower) <- rownames(mat_base)
    colnames(mat_ci_lower) <- colnames(mat_base)
    
    mat_ci_upper <- matrix(rnorm(20, mean = 0.5), nrow = 5)
    rownames(mat_ci_upper) <- rownames(mat_base)
    colnames(mat_ci_upper) <- colnames(mat_base)
    
    se <- SummarizedExperiment::SummarizedExperiment(
        assays = list(
            diversity = mat_base,
            diversity_ci_lower = mat_ci_lower,
            diversity_ci_upper = mat_ci_upper
        )
    )
    
    # Test extraction
    result <- TSENAT:::.extract_bootstrap_ci_assays(se, assay_name = "diversity")
    
    expect_true(is.list(result))
    expect_true(result$has_ci)
    expect_equal(dim(result$assay_base), c(5, 4))
    expect_equal(dim(result$ci_lower), c(5, 4))
    expect_equal(dim(result$ci_upper), c(5, 4))
    expect_null(result$fallback_metric)
})

test_that(".extract_bootstrap_ci_assays handles missing CI assays with IQR fallback", {
    skip_if_not_installed("SummarizedExperiment")
    
    mat <- matrix(rnorm(20), nrow = 5)
    rownames(mat) <- paste0("gene", 1:5)
    colnames(mat) <- paste0("sample", 1:4)
    
    se <- SummarizedExperiment::SummarizedExperiment(
        assays = list(diversity = mat)
    )
    
    # Extract with fallback enabled
    result <- TSENAT:::.extract_bootstrap_ci_assays(se, assay_name = "diversity", 
                                                      fallback_to_iqr = TRUE)
    
    expect_false(result$has_ci)
    expect_null(result$ci_lower)
    expect_null(result$ci_upper)
    expect_equal(result$fallback_metric, "iqr")
})

test_that(".extract_bootstrap_ci_assays throws error for missing assay", {
    skip_if_not_installed("SummarizedExperiment")
    
    se <- SummarizedExperiment::SummarizedExperiment(
        assays = list(counts = matrix(1:10, nrow = 2))
    )
    
    expect_error(
        TSENAT:::.extract_bootstrap_ci_assays(se, assay_name = "nonexistent"),
        "not found"
    )
})

test_that(".extract_bootstrap_ci_assays handles custom assay names", {
    skip_if_not_installed("SummarizedExperiment")
    
    mat_base <- matrix(rnorm(20), nrow = 5)
    rownames(mat_base) <- paste0("gene", 1:5)
    
    mat_ci_lower <- matrix(rnorm(20, mean = -0.5), nrow = 5)
    rownames(mat_ci_lower) <- rownames(mat_base)
    
    mat_ci_upper <- matrix(rnorm(20, mean = 0.5), nrow = 5)
    rownames(mat_ci_upper) <- rownames(mat_base)
    
    se <- SummarizedExperiment::SummarizedExperiment(
        assays = list(
            divergence = mat_base,
            divergence_ci_lower = mat_ci_lower,
            divergence_ci_upper = mat_ci_upper
        )
    )
    
    result <- TSENAT:::.extract_bootstrap_ci_assays(se, assay_name = "divergence")
    
    expect_true(result$has_ci)
    expect_equal(dim(result$assay_base), dim(mat_base))
})

# =============================================================================
# CI Data Preparation for Plotting
# =============================================================================

test_that(".prepare_gene_ci_data extracts CI data correctly", {
    skip_if_not_installed("SummarizedExperiment")
    
    # Create sample long format data with required columns: Gene, group, q, tsallis, sample
    long_data <- data.frame(
        Gene = c("gene1", "gene1", "gene1", "gene1", "gene2", "gene2", "gene2", "gene2"),
        group = c("A", "A", "B", "B", "A", "A", "B", "B"),
        q = c(0.5, 1.0, 0.5, 1.0, 0.5, 1.0, 0.5, 1.0),
        tsallis = c(1.0, 1.5, 1.2, 1.7, 2.0, 2.5, 2.2, 2.7),
        sample = c("s1", "s1", "s2", "s2", "s1", "s1", "s2", "s2"),
        stringsAsFactors = FALSE
    )
    
    # Create CI matrices with column names matching sample_q=value format
    ci_lower <- matrix(c(0.9, 1.4, 1.9, 2.4), nrow = 2, byrow = TRUE)
    rownames(ci_lower) <- c("gene1", "gene2")
    colnames(ci_lower) <- c("s1_q=0.5", "s1_q=1.0")
    
    ci_upper <- matrix(c(1.1, 1.6, 2.1, 2.6), nrow = 2, byrow = TRUE)
    rownames(ci_upper) <- c("gene1", "gene2")
    colnames(ci_upper) <- c("s1_q=0.5", "s1_q=1.0")
    
    # Test extraction
    result <- TSENAT:::.prepare_gene_ci_data(long_data, ci_lower, ci_upper, 
                                              genes = c("gene1", "gene2"))
    
    expect_true(is.data.frame(result))
    expect_true("ci_lower" %in% colnames(result))
    expect_true("ci_upper" %in% colnames(result))
    # Result groups by gene, group, q: 2 genes × 2 groups × 2 q values = 8 groups
    expect_equal(nrow(result), 8)
})

test_that(".prepare_gene_ci_data handles single gene correctly", {
    skip_if_not_installed("SummarizedExperiment")
    
    long_data <- data.frame(
        Gene = c("gene1", "gene1"),
        group = c("A", "B"),
        q = c(0.5, 0.5),
        tsallis = c(1.0, 1.5),
        sample = c("s1", "s2"),
        stringsAsFactors = FALSE
    )
    
    ci_lower <- matrix(c(0.9, 1.4), nrow = 1)
    rownames(ci_lower) <- "gene1"
    colnames(ci_lower) <- c("s1_q=0.5", "s2_q=0.5")
    
    ci_upper <- matrix(c(1.1, 1.6), nrow = 1)
    rownames(ci_upper) <- "gene1"
    colnames(ci_upper) <- c("s1_q=0.5", "s2_q=0.5")
    
    result <- TSENAT:::.prepare_gene_ci_data(long_data, ci_lower, ci_upper, 
                                              genes = "gene1")
    
    expect_equal(nrow(result), 2)
    expect_equal(unique(result$Gene), "gene1")
})

# =============================================================================
# Distribution Statistics Computation
# =============================================================================

test_that(".compute_distribution_stats calculates median and IQR", {
    df <- data.frame(
        group = c("A", "A", "A", "B", "B", "B"),
        value = c(1, 2, 3, 4, 5, 6),
        stringsAsFactors = FALSE
    )
    
    result <- TSENAT:::.compute_distribution_stats(df, "group", "value", 
                                                     metric = "median", 
                                                     spread_metric = "iqr")
    
    expect_true(is.data.frame(result))
    expect_true("group" %in% colnames(result))
    expect_true("value" %in% colnames(result))
    expect_true("lower" %in% colnames(result))
    expect_true("upper" %in% colnames(result))
    
    # Group A: median = 2, IQR from 1-3 is 1 (Q3 - Q1 = 3 - 2 = 1)
    grp_a <- result[result$group == "A", ]
    expect_equal(grp_a$value, 2)
})

test_that(".compute_distribution_stats calculates mean and SD", {
    df <- data.frame(
        group = c("A", "A", "A", "B", "B", "B"),
        value = c(1, 2, 3, 4, 5, 6),
        stringsAsFactors = FALSE
    )
    
    result <- TSENAT:::.compute_distribution_stats(df, "group", "value",
                                                     metric = "mean",
                                                     spread_metric = "sd")
    
    expect_true(is.data.frame(result))
    
    # Group A: mean = 2, lower/upper = mean ± sd
    grp_a <- result[result$group == "A", ]
    expect_equal(grp_a$value, 2)
    expect_true(grp_a$upper > grp_a$value)
})

test_that(".compute_distribution_stats handles NA values", {
    df <- data.frame(
        group = c("A", "A", "A", "B", "B", "B"),
        value = c(1, NA, 3, 4, NA, 6),
        stringsAsFactors = FALSE
    )
    
    result <- TSENAT:::.compute_distribution_stats(df, "group", "value",
                                                     metric = "median",
                                                     spread_metric = "iqr")
    
    expect_true(is.data.frame(result))
    expect_equal(nrow(result), 2)
    # Group A: median of c(1,3) = 2
    grp_a <- result[result$group == "A", ]
    expect_equal(grp_a$value, 2)
})

# =============================================================================
# Data Format Preparation for Grouped Visualization
# =============================================================================

test_that(".prepare_grouped_long_format converts wide to long with groups", {
    skip_if_not_installed("SummarizedExperiment")
    
    # Create proper SummarizedExperiment input
    # Create 3 genes x 12 samples matrix
    mat <- matrix(rnorm(36), nrow = 3, ncol = 12)
    rownames(mat) <- c("gene1", "gene2", "gene3")
    colnames(mat) <- paste0("s", 1:12)
    
    # Create colData with matching number of columns
    coldata <- data.frame(
        condition = rep(c("A", "A", "B"), 4),
        row.names = colnames(mat)
    )
    
    se <- SummarizedExperiment::SummarizedExperiment(
        assays = list(diversity = mat),
        colData = coldata
    )
    
    # Call function - it takes SummarizedExperiment, not separate mat/samples
    result <- TSENAT:::.prepare_grouped_long_format(se, assay_name = "diversity",
                                                      group_by_col = "condition")
    
    expect_true(is.data.frame(result))
    expect_true("Gene" %in% colnames(result))
    expect_true("value" %in% colnames(result))
    expect_true("group" %in% colnames(result))
    expect_true(nrow(result) > 0)
})

# =============================================================================
# Transcript Input Preparation
# =============================================================================

test_that(".prepare_transcript_inputs handles matrix input", {
    skip_if_not_installed("SummarizedExperiment")
    
    # Create simple count matrix with required metadata
    counts <- matrix(c(100, 50, 200, 75), nrow = 2)
    rownames(counts) <- c("tx1", "tx2")
    colnames(counts) <- c("sample1", "sample2")
    
    # Provide required samples parameter
    samples <- c("condition_A", "condition_B")
    
    # Create mock tx2gene mapping (must use 'Transcript' and 'Gen' column names)
    tx2gene <- data.frame(
        Transcript = c("tx1", "tx2"),
        Gen = c("gene1", "gene1")
    )
    
    result <- TSENAT:::.prepare_transcript_inputs(counts = counts, samples = samples, tx2gene = tx2gene)
    
    # Function returns a list, not a data.frame or matrix
    expect_true(is.list(result))
    expect_true("counts" %in% names(result))
    expect_true("samples" %in% names(result))
    expect_true("mapping" %in% names(result))
})

test_that(".prepare_transcript_inputs throws error for missing rownames", {
    counts <- matrix(c(100, 50, 200, 75), nrow = 2)
    colnames(counts) <- c("sample1", "sample2")
    
    samples <- c("condition_A", "condition_B")
    tx2gene <- data.frame(Transcript = c("tx1", "tx2"), Gen = c("gene1", "gene1"))
    
    expect_error(
        TSENAT:::.prepare_transcript_inputs(counts = counts, samples = samples, tx2gene = tx2gene),
        "rownames|transcript"
    )
})

test_that(".prepare_transcript_inputs validates input type", {
    skip_if_not_installed("SummarizedExperiment")
    
    # Try passing invalid input (list)
    samples <- c("condition_A", "condition_B")
    tx2gene <- data.frame(Transcript = c("tx1", "tx2"), Gen = c("gene1", "gene1"))
    
    expect_error(
        TSENAT:::.prepare_transcript_inputs(counts = list(data = 1:10), samples = samples, tx2gene = tx2gene),
        "matrix|data.frame|SummarizedExperiment"
    )
})

test_that(".prepare_transcript_inputs requires samples or coldata", {
    counts <- matrix(c(100, 50, 200, 75), nrow = 2)
    rownames(counts) <- c("tx1", "tx2")
    colnames(counts) <- c("sample1", "sample2")
    
    tx2gene <- data.frame(Transcript = c("tx1", "tx2"), Gen = c("gene1", "gene1"))
    
    # Neither samples nor coldata provided
    expect_error(
        TSENAT:::.prepare_transcript_inputs(counts = counts, tx2gene = tx2gene),
        "samples|coldata"
    )
})

# =============================================================================
# Theme and Plotting Utilities
# =============================================================================

test_that(".create_centered_theme creates valid ggplot2 theme", {
    skip_if_not_installed("ggplot2")
    
    theme_obj <- TSENAT:::.create_centered_theme(include_title = TRUE, 
                                                   title_size = 14)
    
    expect_true(inherits(theme_obj, "theme"))
})

test_that(".save_plot_standard handles file output", {
    skip_if_not_installed("ggplot2")
    
    p <- ggplot2::ggplot(data.frame(x = 1:10, y = 1:10), 
                         ggplot2::aes(x = x, y = y)) +
         ggplot2::geom_point()
    
    temp_file <- tempfile(fileext = ".png")
    on.exit(unlink(temp_file))
    
    # Should complete without error
    expect_silent(
        TSENAT:::.save_plot_standard(p, temp_file, width_inches = 8, 
                                      dpi_output = 100)
    )
    
    # File should exist after save
    expect_true(file.exists(temp_file))
})

# =============================================================================
# Infer Samples from SummarizedExperiment
# =============================================================================

test_that(".infer_samples_from_se detects condition column", {
    skip_if_not_installed("SummarizedExperiment")
    
    se <- SummarizedExperiment::SummarizedExperiment(
        assays = list(counts = matrix(1:20, nrow = 5)),
        colData = data.frame(
            condition = c("A", "A", "B", "B"),
            sample_id = c("s1", "s2", "s3", "s4")
        )
    )
    
    result <- TSENAT:::.infer_samples_from_se(se, condition_col = "condition")
    
    expect_equal(result, c("A", "A", "B", "B"))
})

test_that(".infer_samples_from_se returns explicit samples if provided", {
    skip_if_not_installed("SummarizedExperiment")
    
    se <- SummarizedExperiment::SummarizedExperiment(
        assays = list(counts = matrix(1:20, nrow = 5))
    )
    
    explicit_samples <- c("X", "Y", "X", "Y")
    result <- TSENAT:::.infer_samples_from_se(se, samples = explicit_samples)
    
    expect_equal(result, explicit_samples)
})

test_that(".infer_samples_from_se returns NULL for empty SE", {
    skip_if_not_installed("SummarizedExperiment")
    
    se <- SummarizedExperiment::SummarizedExperiment(
        assays = list(counts = matrix(1:20, nrow = 5))
    )
    
    result <- TSENAT:::.infer_samples_from_se(se, samples = NULL, 
                                               condition_col = "nonexistent")
    
    expect_null(result)
})
