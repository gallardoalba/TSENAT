#!/usr/bin/env Rscript
# Debug script for plot_top_transcripts spacing issues
# Based on TEMPLATE_TESTS.R data preparation

suppressPackageStartupMessages({
    library(devtools)
    devtools::load_all(".")
    library(ggplot2)
    library(SummarizedExperiment)
    library(dplyr)
})

set.seed(42)

message("=== Data Preparation (from TEMPLATE_TESTS.R) ===\n")

# Load data
message("Loading readcounts data...")
data(readcounts)
counts_data <- as.matrix(salmon_dataset)
mode(counts_data) <- "numeric"

# Load TPM data
message("Loading TPM data...")
tpm_data <- as.matrix(salmon_tpm)
mode(tpm_data) <- "numeric"

# Subset for faster analysis
n_genes_subset <- 50
if (!is.null(n_genes_subset) && nrow(counts_data) > n_genes_subset) {
  counts_data <- counts_data[1:n_genes_subset, ]
  tpm_data <- tpm_data[1:n_genes_subset, ]
  message(sprintf("Analyzing %d genes for speed.\n", n_genes_subset))
}

# Load metadata
message("Loading metadata...")
metadata_df <- read.table(
    system.file("extdata", "metadata.tsv", package = "TSENAT"),
    header = TRUE, sep = "\t"
)

gff3_dataset <- system.file("extdata", "annotation.gff3.gz", package = "TSENAT")

# Build SummarizedExperiment
message("Building SummarizedExperiment...")
se <- build_se(counts_data, gff3_dataset, metadata = metadata_df)

# Add TPM as an assay
tpm_subset <- tpm_data[rownames(se), colnames(se)]
assay(se, "tpm", withDimnames = FALSE) <- tpm_subset

message(sprintf("SE dimensions: %d x %d\n", nrow(se), ncol(se)))

# Filter
message("Filtering SE...")
se <- tryCatch(
  filter_se(se, stringency = "loose"),
  error = function(e) se
)

message("=== Generating plot with plot_top_transcripts ===\n")

# Get available genes
available_genes <- unique(as.character(rowData(se)$gene_name))
message(sprintf("Available genes: %s\n", paste(head(available_genes, 5), collapse=", ")))

# Select first 4 genes to plot
genes_to_plot <- head(available_genes, 4)
message(sprintf("Plotting genes: %s\n", paste(genes_to_plot, collapse=", ")))

# Generate the plot with specific genes
message("Creating plot with top_n=4...")
p <- plot_top_transcripts(
    se, 
    gene = genes_to_plot,
    top_n = 4,
    metric = "median"
)

message("Plot generated successfully!")
message(sprintf("Plot class: %s\n", paste(class(p), collapse=", ")))

# Save outputs
dir.create("/home/nouser/galaxy/tools_source/TSENAT/output", showWarnings = FALSE)

output_png <- "/home/nouser/galaxy/tools_source/TSENAT/output/debug_plot_spacing.png"
output_pdf <- "/home/nouser/galaxy/tools_source/TSENAT/output/debug_plot_spacing.pdf"

message(sprintf("Saving to PNG: %s", output_png))
ggplot2::ggsave(output_png, p, width = 14, height = 10, dpi = 100)

message(sprintf("Saving to PDF: %s\n", output_pdf))
ggplot2::ggsave(output_pdf, p, width = 14, height = 10)

message("\n=== Creating sample-based heatmap (plot_multiq style) ===\n")

# Extract count data and create per-gene heatmaps with sample grouping
genes_to_heatmap <- genes_to_plot
counts_mat <- as.matrix(SummarizedExperiment::assay(se))
samples <- as.character(colData(se)$sample_type)

# Create heatmaps for each gene showing transcript x sample expression
require_pkgs(c("pheatmap", "grid"))

output_heatmap_path <- "/home/nouser/galaxy/tools_source/TSENAT/output/debug_plot_sample_heatmaps.png"

# Build heatmap data for each gene
heatmap_list <- list()
gene_names_with_data <- character(0)

for (gene_name in genes_to_heatmap) {
    # Get transcripts for this gene
    gene_txs <- rownames(se)[rowData(se)$gene_name == gene_name]
    
    if (length(gene_txs) == 0) next
    
    # Extract expression matrix for this gene's transcripts x samples
    expr_mat <- counts_mat[gene_txs, , drop = FALSE]
    expr_log <- log2(expr_mat + 1)  # Log transform
    
    # Add gene name to list
    gene_names_with_data <- c(gene_names_with_data, gene_name)
    
    # Create pheatmap for this gene
    p <- pheatmap::pheatmap(
        expr_log,
        main = gene_name,
        cluster_rows = FALSE,
        cluster_cols = FALSE,
        display_numbers = FALSE,
        na_col = "lightgray",
        color = grDevices::colorRampPalette(c("#4575B4", "#FFFFFF", "#D73027"))(100),
        cellwidth = 60,
        cellheight = 60,
        fontsize = 26,
        fontsize_row = 26,
        fontsize_col = 26,
        fontsize_number = 20,
        margins = c(11, 180),
        show_rownames = TRUE,
        show_colnames = TRUE,
        silent = TRUE
    )
    
    heatmap_list[[length(heatmap_list) + 1]] <- p
}

if (length(heatmap_list) > 0) {
    message(sprintf("Creating combined heatmap for %d genes...", length(heatmap_list)))
    
    # Combine into multi-panel figure using grid
    n_genes <- length(heatmap_list)
    n_cols <- 2
    n_rows <- ceiling(n_genes / n_cols)
    
    # Calculate PNG dimensions
    heatmap_height <- 9 * n_rows + 2 * (n_rows - 1) + 5
    
    png(filename = output_heatmap_path, width = 28, height = heatmap_height, units = "in", res = 72)
    
    # Create title and layout structure (similar to plot_multiq_delta_influence_heatmaps)
    grid::grid.newpage()
    grid::pushViewport(grid::viewport(layout = grid::grid.layout(n_rows + 1, n_cols)))
    
    # Title row
    title_vp <- grid::viewport(layout.pos.row = 1, layout.pos.col = 1:n_cols)
    grid::pushViewport(title_vp)
    grid::grid.text("Sample-based Transcript Expression Heatmaps", x = 0.5, y = 0.5,
                    gp = grid::gpar(fontsize = 40, fontface = "bold"))
    grid::upViewport()
    
    # Plot rows
    for (i in seq_along(heatmap_list)) {
        plot_row_idx <- ((i - 1) %/% n_cols) + 2
        plot_col_idx <- ((i - 1) %% n_cols) + 1
        vp <- grid::viewport(layout.pos.row = plot_row_idx, layout.pos.col = plot_col_idx)
        grid::pushViewport(vp)
        grid::grid.draw(heatmap_list[[i]]$gtable)
        grid::upViewport()
    }
    
    grid::upViewport()
    grDevices::dev.off()
    
    message(sprintf("Sample heatmaps saved to: %s", output_heatmap_path))
} else {
    message("No genes with data found for heatmap creation")
}

message("\n=== SUCCESS: All plots saved ===")
message(sprintf("Transcript plot PNG: %s", output_png))
message(sprintf("Transcript plot PDF: %s", output_pdf))
message(sprintf("Sample heatmaps PNG: %s", output_heatmap_path))


