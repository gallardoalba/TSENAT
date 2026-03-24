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

SummarizedExperiment::colData(se) <- S4Vectors::DataFrame(
  condition = rep(c("control", "treatment"), each = n_samples_per_group),
  row.names = colnames(se))

analysis <- TSENATAnalysis(se)
analysis <- calculate_diversity_s4(analysis, q = 1.0, verbose = FALSE)
analysis <- calculate_divergence_s4(analysis, q = 1.0, verbose = FALSE)

suppressWarnings(
  analysis <- calculate_lm_interaction_s4(analysis,
    condition_col = "condition", verbose = FALSE)
)

# Debug effect_sizes call
cat("Checking effect_sizes_divergence inputs:\n")
cat("lm_results exists:", !is.null(analysis@lm_results$lm_interaction), "\n")
cat("lm_results rows:", nrow(analysis@lm_results$lm_interaction), "\n")
cat("lm_results cols:", ncol(analysis@lm_results$lm_interaction), "\n")
cat("\nFirst few LM results:\n")
print(head(analysis@lm_results$lm_interaction[, 1:5]))

cat("\nDivergence SE dims:", dim(analysis@divergence_results$divergence_se), "\n")

cat("\nTrying effect_sizes_divergence_s4:\n")
tryCatch({
  analysis <- effect_sizes_divergence_s4(analysis, verbose = TRUE)
  cat("SUCCESS!\n")
}, error = function(e) {
  cat("ERROR:", e$message, "\n")
})
