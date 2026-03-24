#!/usr/bin/env Rscript
# Use devtools::load_all to load the development package
devtools::load_all("/home/nouser/galaxy/tools_source/TSENAT")

library(SummarizedExperiment)

# Use the EXACT setup from test-wrapper_s4_effect_size_divergence.R
set.seed(456)

n_transcripts <- 100  # Larger than my 30
n_genes <- 20         # Size matches
n_samples <- 12       # Larger than my 10

counts <- matrix(
  rpois(n_transcripts * n_samples, lambda = 50),  # lambda=50, not 35
  nrow = n_transcripts, ncol = n_samples
)
counts <- pmax(counts, 5)  # Ensure minimum count of 5

rownames(counts) <- paste0("TX_", 1:n_transcripts)
colnames(counts) <- paste0("Sample_", 1:n_samples)

rowData <- S4Vectors::DataFrame(
  transcript_id = rownames(counts),
  gene_id = paste0("GENE_", rep(1:n_genes, length.out = n_transcripts)),
  row.names = rownames(counts)
)

colData <- S4Vectors::DataFrame(
  sample_id = colnames(counts),
  condition = rep(c("A", "B"), length.out = n_samples),
  row.names = colnames(counts)
)

se <- SummarizedExperiment::SummarizedExperiment(
  assays = list(counts = counts),
  rowData = rowData,
  colData = colData
)

tx2gene_df <- data.frame(
  Transcript = rownames(counts),
  Gene = rowData$gene_id,
  stringsAsFactors = FALSE
)
S4Vectors::metadata(se)$tx2gene <- tx2gene_df

analysis <- TSENATAnalysis(se = se)

# Calculate diversity
cat("Calculating diversity...\n")
analysis <- calculate_diversity_s4(
  analysis,
  q = 1.0,
  verbose = FALSE,
  min_valid_frac = 0
)

# Calculate divergence
cat("Calculating divergence...\n")
analysis <- calculate_divergence_s4(
  analysis,
  verbose = FALSE
)

# Calculate LM interaction (needed for effect sizes)
cat("Calculating LM interaction...\n")
analysis <- suppressWarnings(
  calculate_lm_interaction_s4(
    analysis,
    verbose = FALSE
  )
)

cat("\nCalling effect_sizes_divergence_s4...\n")
result <- tryCatch({
  analysis <- effect_sizes_divergence_s4(analysis, verbose = TRUE)
  cat("\nSUCCESS! effect_sizes_divergence_s4 completed\n")
  TRUE
}, error = function(e) {
  cat("\nERROR:", e$message, "\n")
  FALSE
})

if (!result) {
  cat("Test FAILED\n")
  quit(save = "no", status = 1)
} else {
  cat("Test PASSED\n")
  quit(save = "no", status = 0)
}
