#!/usr/bin/env Rscript

# Test script for unused plotting functions:
# - plot_multi_q_spectrum
# - plot_q_sensitivity_curve

library(devtools)
load_all('.')

library(SummarizedExperiment)
library(S4Vectors)

# Create output directory if it doesn't exist
output_dir <- "output"
if (!dir.exists(output_dir)) {
  dir.create(output_dir, showWarnings = FALSE)
}

cat("========================================\n")
cat("Testing: plot_multi_q_spectrum\n")
cat("========================================\n\n")

# Generate mock LMM interaction results for plot_multi_q_spectrum
set.seed(42)
n_genes <- 8

# Create interaction results with effect sizes across q-values
interaction_results <- data.frame(
  gene = paste0("GENE", 1:n_genes),
  gene_id = paste0("ENS", 1:n_genes),
  effect_size_D_q0.5 = rnorm(n_genes, mean = 0.5, sd = 0.3),
  effect_size_D_q1.0 = rnorm(n_genes, mean = 0.6, sd = 0.3),
  effect_size_D_q1.5 = rnorm(n_genes, mean = 0.55, sd = 0.3),
  effect_size_D_q2.0 = rnorm(n_genes, mean = 0.4, sd = 0.3),
  per_q_pattern = rep(NA, n_genes),
  stringsAsFactors = FALSE
)

# Add per_q_pattern as comma-separated values for each gene
interaction_results$per_q_pattern <- apply(
  interaction_results[, c("effect_size_D_q0.5", "effect_size_D_q1.0", 
                          "effect_size_D_q1.5", "effect_size_D_q2.0")], 
  1, 
  function(x) paste(x, collapse = ",")
)

# Create the lmm_results structure
lmm_results <- list(
  interaction_results = interaction_results
)

cat("Created mock LMM results with", nrow(interaction_results), "genes\n")
cat("Columns:", paste(colnames(interaction_results), collapse = ", "), "\n\n")

# Test plot_multi_q_spectrum
png_file_1 <- file.path(output_dir, "plot_multi_q_spectrum.png")
png(png_file_1, width = 1200, height = 400)

cat("Generating plot_multi_q_spectrum...\n")
plot_multi_q_spectrum(lmm_results, n_genes = 5)
mtext("Top 5 Genes - Multi-Q Spectrum", outer = TRUE, side = 3, line = -0.5, cex = 1.2)

dev.off()
cat("✓ Saved to:", png_file_1, "\n\n")

# =====================================
cat("========================================\n")
cat("Testing: plot_q_sensitivity_curve\n")
cat("========================================\n\n")

# Create sample SummarizedExperiment for plot_q_sensitivity_curve
set.seed(123)
n_transcripts <- 6
n_samples <- 8

# Create counts matrix: transcripts × samples
counts <- matrix(
  rnbinom(n_transcripts * n_samples, size = 10, prob = 0.5),
  nrow = n_transcripts,
  ncol = n_samples
)

# Add reasonable count levels
counts <- counts + 50  # Baseline counts

rownames(counts) <- paste0("TX", 1:n_transcripts)
colnames(counts) <- paste0("S", 1:n_samples)

# Create row data with gene and isoform info
rowData_df <- DataFrame(
  transcript_id = paste0("TX", 1:n_transcripts),
  gene = rep(c("GENE_A", "GENE_B", "GENE_C"), each = 2),
  gene_id = rep(c("ENS0001", "ENS0002", "ENS0003"), each = 2),
  stringsAsFactors = FALSE
)

# Create column data with conditions
colData_df <- DataFrame(
  sample = colnames(counts),
  condition = c(rep("Control", 4), rep("Treatment", 4)),
  sample_type = c(rep("WT", 4), rep("KO", 4)),
  stringsAsFactors = FALSE
)

# Create SE object
se <- SummarizedExperiment(
  assays = list(counts = counts),
  rowData = rowData_df,
  colData = colData_df
)

cat("Created SummarizedExperiment with:\n")
cat("  -", nrow(se), "transcripts (3 genes × 2 isoforms each)\n")
cat("  -", ncol(se), "samples (4 Control + 4 Treatment)\n\n")

# Generate plots for multiple genes
png_file_2 <- file.path(output_dir, "plot_q_sensitivity_curve_grid.png")

# Test for 3 different genes with multi-panel layout
png(png_file_2, width = 1200, height = 400)
par(mfrow = c(1, 3), mar = c(4, 4, 3, 2), oma = c(0, 0, 2, 0))

gene_list <- c("GENE_A", "GENE_B", "GENE_C")

for (gene in gene_list) {
  cat("Plotting sensitivity curve for:", gene, "\n")
  
  plot_q_sensitivity_curve(
    se = se,
    condition_col = "condition",
    gene = gene,
    gene_col = "gene",
    isoform_col = "transcript_id",
    q_values = c(0.5, 1.0, 1.5, 2.0),
    norm = TRUE,
    main = gene,
    show_legend = FALSE
  )
}

mtext("Q-Sensitivity Curves: Entropy Changes Across Diversity Scales", 
      outer = TRUE, side = 3, line = 0.5, cex = 1.2, font = 2)

dev.off()

cat("\n✓ Saved to:", png_file_2, "\n")

# Also create individual plots for each gene
for (gene in gene_list) {
  png_file_ind <- file.path(output_dir, paste0("plot_q_sensitivity_", gene, ".png"))
  png(png_file_ind, width = 600, height = 500)
  
  plot_q_sensitivity_curve(
    se = se,
    condition_col = "condition",
    gene = gene,
    gene_col = "gene",
    isoform_col = "transcript_id",
    q_values = c(0.5, 1.0, 1.5, 2.0, 2.5),
    norm = TRUE,
    main = paste0("Q-Sensitivity: ", gene),
    show_legend = TRUE
  )
  
  dev.off()
  cat("✓ Saved to:", png_file_ind, "\n")
}

cat("\n========================================\n")
cat("All tests completed successfully!\n")
cat("Output files saved to:", output_dir, "/\n")
cat("========================================\n")
