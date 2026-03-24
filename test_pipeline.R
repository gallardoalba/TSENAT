#!/usr/bin/env Rscript
library(TSENAT)
library(SummarizedExperiment)
set.seed(42)
n_genes <- 10
n_isoforms_per_gene <- 3
n_isoforms <- n_genes * n_isoforms_per_gene
n_samples_per_group <- 5
n_samples <- n_samples_per_group * 2

control_counts <- matrix(rpois(n_isoforms * n_samples_per_group, lambda = 35),
                         nrow = n_isoforms, ncol = n_samples_per_group)
treatment_counts <- matrix(rpois(n_isoforms * n_samples_per_group, lambda = 35),
                           nrow = n_isoforms, ncol = n_samples_per_group)
counts <- cbind(control_counts, treatment_counts)
rownames(counts) <- paste0("TX_", 1:n_isoforms)
colnames(counts) <- paste0("Sample_", 1:n_samples)
se <- SummarizedExperiment(assays = list(counts = counts))
S4Vectors::metadata(se)$tx2gene <- data.frame(
  Transcript = rownames(counts),
  Gene = rep(paste0("GENE_", 1:n_genes), each = n_isoforms_per_gene))

cat("colnames(se):\n")
print(colnames(se))

SummarizedExperiment::colData(se) <- S4Vectors::DataFrame(
  condition = rep(c("control", "treatment"), each = n_samples_per_group),
  row.names = colnames(se))

cat("\ncolData(se):\n")
print(colData(se))

analysis <- TSENATAnalysis(se)
analysis <- calculate_diversity_s4(analysis, q = 1.0, verbose = FALSE)
analysis <- calculate_divergence_s4(analysis, q = 1.0, verbose = FALSE)

cat("Before calculate_lm_interaction_s4:\n")
cat("lm_results exists:", !is.null(analysis@lm_results$lm_interaction), "\n")

suppressWarnings(
  analysis <- calculate_lm_interaction_s4(analysis,
    condition_col = "condition", verbose = FALSE)
)

cat("After calculate_lm_interaction_s4:\n")
cat("lm_results exists:", !is.null(analysis@lm_results$lm_interaction), "\n")

cat("\nAll tests completed successfully!\n")
