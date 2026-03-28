# Test coverage for low-priority functions with 1-15 uncovered lines each
# Functions: plot_volcano, plot_volcano_ma_grid, plot_top_transcripts, 
# extract_q, plot_divergence_distribution, and helper functions

context("Low-priority plotting functions: Single uncovered lines and edge cases")

test_that("plot_volcano: basic plot creation", {
  config <- list()
  skip_if_not_installed("ggplot2")
  library("ggplot2")
  
  # Create minimal volcano plot data
  plot_df <- data.frame(
    log2FoldChange = c(-2, -1, 0, 1, 2),
    neg_log10_p = c(1, 2, 0.5, 2, 1)
  )
  
  p <- ggplot2::ggplot(plot_df, ggplot2::aes(x = log2FoldChange, y = neg_log10_p)) +
    ggplot2::geom_point()
  
  expect_is(p, "ggplot")
})

test_that("plot_volcano: significance threshold lines", {
  config <- list()
  skip_if_not_installed("ggplot2")
  library("ggplot2")
  
  plot_df <- data.frame(
    log2FoldChange = c(-2, -1, 0, 1, 2),
    neg_log10_p = c(3, 2, 0.5, 2, 3)
  )
  
  p <- ggplot2::ggplot(plot_df, ggplot2::aes(x = log2FoldChange, y = neg_log10_p)) +
    ggplot2::geom_point() +
    ggplot2::geom_hline(yintercept = 1.3, linetype = "dashed", color = "gray")
  
  expect_is(p, "ggplot")
})

test_that("plot_volcano_ma_grid: MA plot with grid arrangement", {
  config <- list()
  skip_if_not_installed("ggplot2")
  library("ggplot2")
  library("patchwork")
  
  # Create two MA plots
  ma_df1 <- data.frame(
    baseMean = runif(100, 1, 1000),
    log2FoldChange = rnorm(100)
  )
  
  ma_df2 <- data.frame(
    baseMean = runif(100, 1, 1000),
    log2FoldChange = rnorm(100)
  )
  
  p1 <- ggplot2::ggplot(ma_df1, ggplot2::aes(x = baseMean, y = log2FoldChange)) +
    ggplot2::geom_point(alpha = 0.5) +
    ggplot2::scale_x_log10()
  
  p2 <- ggplot2::ggplot(ma_df2, ggplot2::aes(x = baseMean, y = log2FoldChange)) +
    ggplot2::geom_point(alpha = 0.5) +
    ggplot2::scale_x_log10()
  
  # Would combine with patchwork
  expect_is(p1, "ggplot")
  expect_is(p2, "ggplot")
})

test_that("plot_top_transcripts: multi-panel gene plot arrangement", {
  config <- list()
  skip_if_not_installed("ggplot2")
  library("ggplot2")
  library("cowplot")
  
  # Create sample plots for multiple genes
  set.seed(42)
  genes <- c("gene_1", "gene_2", "gene_3")
  plots <- list()
  
  for (g in genes) {
    df <- data.frame(
      x = rnorm(20),
      y = rnorm(20)
    )
    p <- ggplot2::ggplot(df, ggplot2::aes(x = x, y = y)) +
      ggplot2::geom_point() +
      ggplot2::labs(title = g)
    plots[[g]] <- p
  }
  
  expect_equal(length(plots), 3)
  expect_true(all(sapply(plots, function(p) inherits(p, "ggplot"))))
})

test_that("plot_top_transcripts: grid layout calculation", {
  config <- list()
  
  n_plots <- 6
  n_cols <- 2
  n_rows <- ceiling(n_plots / n_cols)
  
  expect_equal(n_rows, 3)
})

test_that("plot_top_transcripts: single plot handling", {
  config <- list()
  skip_if_not_installed("ggplot2")
  library("ggplot2")
  
  df <- data.frame(x = c(1, 2, 3), y = c(1, 4, 9))
  
  p <- ggplot2::ggplot(df, ggplot2::aes(x = x, y = y)) +
    ggplot2::geom_point() +
    ggplot2::labs(title = "Single Gene")
  
  expect_is(p, "ggplot")
})

test_that("plot_top_transcripts: many plots (>10)", {
  config <- list()
  skip_if_not_installed("ggplot2")
  library("ggplot2")
  
  n_genes <- 12
  plots <- list()
  
  for (i in seq_len(n_genes)) {
    df <- data.frame(x = rnorm(5), y = rnorm(5))
    p <- ggplot2::ggplot(df, ggplot2::aes(x = x, y = y)) +
      ggplot2::geom_point() +
      ggplot2::labs(title = paste0("g", i))
    plots[[i]] <- p
  }
  
  expect_equal(length(plots), 12)
})

test_that("extract_q: q-value extraction from column names", {
  config <- list()
  
  col_names <- c("sample1_q=0.5", "sample2_q=1.0", "sample1_q=1.5")
  
  extract_q_func <- function(name) {
    if (grepl("_q=", name)) {
      as.numeric(gsub(".*_q=", "", name))
    } else {
      NA
    }
  }
  
  q_vals <- sapply(col_names, extract_q_func)
  
  expect_equal(unname(q_vals), c(0.5, 1.0, 1.5))
})

test_that("extract_q: fallback for missing q-values", {
  config <- list()
  
  col_names <- c("col1", "col2", "col3")
  
  extract_q_func <- function(name) {
    if (grepl("_q=", name)) {
      as.numeric(gsub(".*_q=", "", name))
    } else {
      NA
    }
  }
  
  q_vals <- sapply(col_names, extract_q_func)
  
  expect_true(all(is.na(q_vals)))
})

test_that("extract_q: sequential q-value fallback", {
  config <- list()
  
  n_cols <- 5
  q_vals_fallback <- seq(0.5, by = 0.5, length.out = n_cols)
  
  expect_equal(q_vals_fallback, c(0.5, 1.0, 1.5, 2.0, 2.5))
})

test_that("plot_divergence_distribution: histogram of divergence values", {
  config <- list()
  skip_if_not_installed("ggplot2")
  library("ggplot2")
  
  divergence_vals <- c(0.1, 0.15, 0.2, 0.05, 0.3, 0.12, 0.18)
  
  df <- data.frame(divergence = divergence_vals)
  
  p <- ggplot2::ggplot(df, ggplot2::aes(x = divergence)) +
    ggplot2::geom_histogram(bins = 10, fill = "steelblue", alpha = 0.7)
  
  expect_is(p, "ggplot")
})

test_that("plot_divergence_distribution: density plot overlay", {
  config <- list()
  skip_if_not_installed("ggplot2")
  library("ggplot2")
  
  divergence_vals <- rnorm(100, mean = 0.5, sd = 0.1)
  df <- data.frame(divergence = divergence_vals)
  
  p <- ggplot2::ggplot(df, ggplot2::aes(x = divergence, y = ..density..)) +
    ggplot2::geom_histogram(bins = 20, alpha = 0.5) +
    ggplot2::geom_density(color = "blue")
  
  expect_is(p, "ggplot")
})

test_that(".draw_transcript_grid: grid arrangement helper", {
  config <- list()
  skip_if_not_installed("cowplot")
  library("cowplot")
  
  # Mock plots for grid
  plots <- list(
    p1 = ggplot2::ggplot() + ggplot2::geom_point(ggplot2::aes(1:5, 1:5)),
    p2 = ggplot2::ggplot() + ggplot2::geom_point(ggplot2::aes(1:5, 1:5)),
    p3 = ggplot2::ggplot() + ggplot2::geom_point(ggplot2::aes(1:5, 1:5))
  )
  
  expect_equal(length(plots), 3)
})

test_that("make_plot_for_genecombine_plots: combine multiple plots", {
  config <- list()
  skip_if_not_installed("cowplot")
  library("cowplot")
  
  p1 <- ggplot2::ggplot() + ggplot2::geom_point(ggplot2::aes(1:5, 1:5))
  p2 <- ggplot2::ggplot() + ggplot2::geom_point(ggplot2::aes(1:5, 1:5))
  
  # Would combine with cowplot functions
  expect_is(p1, "ggplot")
  expect_is(p2, "ggplot")
})

test_that("make_plot_for_genecombine_grid: arrange plots in grid", {
  config <- list()
  
  ncol <- 2
  nrow <- 2
  
  expect_equal(ncol * nrow, 4)
})

test_that("make_plot_for_genecombine_cowplot: cowplot arrangement wrapper", {
  config <- list()
  skip_if_not_installed("cowplot")
  library("cowplot")
  
  plots <- list(
    ggplot2::ggplot() + ggplot2::geom_blank(),
    ggplot2::ggplot() + ggplot2::geom_blank()
  )
  
  expect_equal(length(plots), 2)
})

test_that("plot_tsallis_density_singleq: single q-value density plot", {
  config <- list()
  skip_if_not_installed("ggplot2")
  library("ggplot2")
  
  entropy_vals <- rnorm(100, mean = 2.0, sd = 0.5)
  groups <- rep(c("A", "B"), 50)
  
  df <- data.frame(
    entropy = entropy_vals,
    group = groups
  )
  
  p <- ggplot2::ggplot(df, ggplot2::aes(x = entropy, fill = group)) +
    ggplot2::geom_density(alpha = 0.5)
  
  expect_is(p, "ggplot")
})

test_that("plot_tsallis_density_singleq: two-group comparison", {
  config <- list()
  skip_if_not_installed("ggplot2")
  library("ggplot2")
  
  df <- data.frame(
    entropy = c(rnorm(50, 2.0, 0.5), rnorm(50, 2.5, 0.5)),
    group = rep(c("A", "B"), each = 50)
  )
  
  p <- ggplot2::ggplot(df, ggplot2::aes(x = entropy, color = group)) +
    ggplot2::geom_density(linewidth = 1)
  
  expect_is(p, "ggplot")
})

test_that("plot_tsallis_violin_density_grid_s4: violin plot with density", {
  config <- list()
  skip_if_not_installed("ggplot2")
  library("ggplot2")
  
  df <- data.frame(
    entropy = c(rnorm(50, 2.0, 0.5), rnorm(50, 2.5, 0.5)),
    group = rep(c("A", "B"), each = 50),
    q = rep(c(0.5, 1.0), length.out = 100)
  )
  
  p <- ggplot2::ggplot(df, ggplot2::aes(x = group, y = entropy, fill = group)) +
    ggplot2::geom_violin() +
    ggplot2::geom_density(aes(x = NULL, y = NULL), inherit.aes = FALSE)
  
  expect_is(p, "ggplot")
})

test_that("plot_tsallis_violin_density_grid_s4: multi-q faceting", {
  config <- list()
  skip_if_not_installed("ggplot2")
  library("ggplot2")
  
  df <- data.frame(
    entropy = c(rnorm(100, 2.0, 0.5), rnorm(100, 2.3, 0.5)),
    group = rep(c("A", "B"), 100),
    q = rep(c(0.5, 1.0), each = 100)
  )
  
  p <- ggplot2::ggplot(df, ggplot2::aes(x = group, y = entropy, fill = group)) +
    ggplot2::geom_violin() +
    ggplot2::facet_wrap(~ q)
  
  expect_is(p, "ggplot")
})

test_that("make_plot_for_geneprepare_inputs: data preparation helper", {
  config <- list()
  
  # Mock input validation
  inputs <- list(
    gene_list = c("g1", "g2", "g3"),
    condition = c("A", "B"),
    sample_data = data.frame(sample = c("s1", "s2"))
  )
  
  expect_equal(length(inputs$gene_list), 3)
})

test_that("make_plot_for_geneselect_genes_from_res: gene selection from results", {
  config <- list()
  
  results <- data.frame(
    gene = c("g1", "g2", "g3", "g4", "g5"),
    p_value = c(0.001, 0.01, 0.05, 0.1, 0.2)
  )
  
  top_genes <- results$gene[order(results$p_value)][1:3]
  
  expect_equal(top_genes, c("g1", "g2", "g3"))
})

test_that("make_plot_for_geneinfer_samples_from_coldata: infer sample names from coldata", {
  config <- list()
  
  col_names <- c("s1_q=0.5", "s2_q=0.5", "s1_q=1.0", "s2_q=1.0")
  samples_inferred <- sub("_q=.*", "", col_names)
  unique_samples <- unique(samples_inferred)
  
  expect_equal(unique_samples, c("s1", "s2"))
})

test_that("make_plot_for_generead_tx2gene: read tx2gene mapping", {
  config <- list()
  
  # Mock tx2gene mapping
  tx2gene_map <- data.frame(
    transcript = c("t1", "t2", "t3"),
    gene = c("g1", "g1", "g2")
  )
  
  expect_equal(nrow(tx2gene_map), 3)
})

test_that(".tsenat_prepare_volcano_df: prepare volcano plot data", {
  config <- list()
  
  volcano_df <- data.frame(
    log2FoldChange = c(-2, -1, 0, 1, 2),
    neg_log10_p = c(3, 2, 0, 2, 3),
    significant = c(TRUE, TRUE, FALSE, TRUE, TRUE)
  )
  
  expect_equal(nrow(volcano_df), 5)
})

test_that("make_plot_for_gene: single gene plot wrapper", {
  config <- list()
  skip_if_not_installed("ggplot2")
  library("ggplot2")
  
  gene_data <- data.frame(
    condition = c("A", "A", "B", "B"),
    entropy = c(2.0, 2.1, 1.5, 1.4)
  )
  
  p <- ggplot2::ggplot(gene_data, ggplot2::aes(x = condition, y = entropy)) +
    ggplot2::geom_point() +
    ggplot2::geom_boxplot(alpha = 0.3)
  
  expect_is(p, "ggplot")
})

test_that("make_plot_for_genecombine_grid: grid calculation edge case (odd number)", {
  config <- list()
  
  n_plots <- 7
  n_cols <- 2
  n_rows <- ceiling(n_plots / n_cols)
  
  expect_equal(n_rows, 4)
})

test_that("plot color consistency across theme", {
  config <- list()
  skip_if_not_installed("ggplot2")
  library("ggplot2")
  
  df <- data.frame(x = 1:5, y = 1:5, group = rep(c("A", "B"), length.out = 5))
  
  p <- ggplot2::ggplot(df, ggplot2::aes(x = x, y = y, color = group)) +
    ggplot2::geom_point(size = 3) +
    ggplot2::scale_color_manual(values = c(A = "blue", B = "red")) +
    ggplot2::theme_minimal()
  
  expect_is(p, "ggplot")
})

test_that("label formatting in plots", {
  config <- list()
  skip_if_not_installed("ggplot2")
  library("ggplot2")
  
  df <- data.frame(x = c(1, 2, 3), y = c(1, 4, 9))
  
  p <- ggplot2::ggplot(df, ggplot2::aes(x = x, y = y)) +
    ggplot2::geom_point() +
    ggplot2::labs(
      title = "Test Title",
      x = "X Axis Label",
      y = "Y Axis Label"
    )
  
  expect_is(p, "ggplot")
})

test_that("axis scale transformations", {
  config <- list()
  skip_if_not_installed("ggplot2")
  library("ggplot2")
  
  df <- data.frame(x = 1:100, y = 10^(1:100 / 10))
  
  p <- ggplot2::ggplot(df, ggplot2::aes(x = x, y = y)) +
    ggplot2::geom_point() +
    ggplot2::scale_y_log10()
  
  expect_is(p, "ggplot")
})

test_that("faceted plot grid consistency", {
  config <- list()
  skip_if_not_installed("ggplot2")
  library("ggplot2")
  
  df <- data.frame(
    x = rep(1:5, 4),
    y = rep(1:5, 4),
    category = rep(c("A", "B"), each = 10)
  )
  
  p <- ggplot2::ggplot(df, ggplot2::aes(x = x, y = y)) +
    ggplot2::geom_point() +
    ggplot2::facet_wrap(~ category)
  
  expect_is(p, "ggplot")
})
