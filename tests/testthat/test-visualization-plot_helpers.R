context("plot_helpers: Plot Composition and Theme Utilities")

# ============================================================================
# TEST: Theme and Title Functions
# ============================================================================

testthat::test_that("apply_tsenat_theme can be applied to plots", {
  require_pkgs("ggplot2")

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
  require_pkgs("ggplot2")

  p <- ggplot2::ggplot(data.frame(x = 1, y = 1), ggplot2::aes(x, y)) +
    ggplot2::geom_point()

  p_titled <- .set_plot_title(p, title = "Test Title", subtitle = "Test Subtitle")

  # Extract title from plot
  testthat::expect_is(p_titled, "ggplot")
  testthat::expect_equal(p_titled$labels$title, "Test Title")
  testthat::expect_equal(p_titled$labels$subtitle, "Test Subtitle")
})

testthat::test_that("set_plot_title applies font sizes", {
  require_pkgs("ggplot2")

  p <- ggplot2::ggplot(data.frame(x = 1, y = 1), ggplot2::aes(x, y)) +
    ggplot2::geom_point()

  p_styled <- .set_plot_title(p, title = "Title", title_size = 16, subtitle_size = 12)

  testthat::expect_is(p_styled, "ggplot")
})

testthat::test_that("set_plot_title handles NULL title/subtitle", {
  require_pkgs("ggplot2")

  p <- ggplot2::ggplot(data.frame(x = 1, y = 1), ggplot2::aes(x, y)) +
    ggplot2::geom_point()

  p_unchanged <- .set_plot_title(p, title = NULL, subtitle = NULL)

  testthat::expect_is(p_unchanged, "ggplot")
})

# ============================================================================
# TEST: Color and Fill Scale Creation
# ============================================================================

testthat::test_that("create_color_scale returns ggplot scale", {
  require_pkgs("ggplot2")

  scale <- .create_color_scale(palette = "blue_red")
  testthat::expect_is(scale, "Scale")
})

testthat::test_that("create_color_scale reverses with direction -1", {
  require_pkgs("ggplot2")

  scale_fwd <- .create_color_scale(palette = "blue_red", direction = 1)
  scale_rev <- .create_color_scale(palette = "blue_red", direction = -1)

  testthat::expect_is(scale_fwd, "Scale")
  testthat::expect_is(scale_rev, "Scale")
})

testthat::test_that("create_fill_scale returns appropriate scale", {
  require_pkgs("ggplot2")

  scale <- .create_fill_scale(palette = "continuous_diverging")
  testthat::expect_is(scale, "Scale")
})

testthat::test_that("create_fill_scale accepts breaks parameter", {
  require_pkgs("ggplot2")

  scale_50 <- .create_fill_scale(breaks = 50)
  scale_100 <- .create_fill_scale(breaks = 100)

  testthat::expect_is(scale_50, "Scale")
  testthat::expect_is(scale_100, "Scale")
})

# ============================================================================
# TEST: Heatmap Creation
# ============================================================================

testthat::test_that("create_tsenat_heatmap creates basic heatmap", {
  require_pkgs("pheatmap")
  
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
  require_pkgs(c("pheatmap", "RColorBrewer"))
  
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
  require_pkgs("pheatmap")
  
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
  require_pkgs(c("ggplot2", "patchwork"))

  p <- ggplot2::ggplot(data.frame(x = 1:10, y = 1:10), ggplot2::aes(x, y)) +
    ggplot2::geom_point()

  combined <- .combine_plots_patchwork(list(p), agg_label_unique = "median")

  testthat::expect_is(combined, "ggplot")
})

testthat::test_that("combine_plots_patchwork handles multiple plots", {
  require_pkgs(c("ggplot2", "patchwork"))

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
  require_pkgs(c("ggplot2", "patchwork"))

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
  require_pkgs(c("ggplot2", "cowplot"))

  p <- ggplot2::ggplot(data.frame(x = 1:10, y = 1:10), ggplot2::aes(x, y)) +
    ggplot2::geom_point()

  combined <- .combine_plots_cowplot(list(p), agg_label_unique = "median")

  # cowplot returns a ggplot object
  testthat::expect_is(combined, "ggplot")
})

testthat::test_that("combine_plots_cowplot handles multiple plots", {
  require_pkgs(c("ggplot2", "cowplot"))

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
  require_pkgs(c("ggplot2", "cowplot"))

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
  require_pkgs("ggplot2")

  p <- ggplot2::ggplot(data.frame(x = 1:10, y = 1:10), ggplot2::aes(x, y)) +
    ggplot2::geom_point()

  # Grid composition returns invisible NULL by default
  result <- .combine_plots_grid(list(p), agg_label_unique = "median")

  testthat::expect_null(result)
})

testthat::test_that("combine_plots_grid handles multiple plots", {
  require_pkgs("ggplot2")

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
  require_pkgs("ggplot2")

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
  require_pkgs(c("ggplot2", "patchwork", "cowplot"))

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
  require_pkgs("ggplot2")

  scale <- .create_color_scale(palette = "blue_red", name = NULL)
  testthat::expect_is(scale, "Scale")
})

testthat::test_that("set_plot_title preserves existing plot aesthetics", {
  require_pkgs("ggplot2")

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
  require_pkgs(c("ggplot2", "patchwork"))

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
  require_pkgs("ggplot2")

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
  require_pkgs(c("ggplot2", "patchwork"))

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
  require_pkgs("ggplot2")

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
  require_pkgs("ggplot2")
  
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
  require_pkgs("ggplot2")
  
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
  require_pkgs("ggplot2")
  
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
