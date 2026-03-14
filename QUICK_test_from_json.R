#!/usr/bin/env Rscript
# Load precomputed results from CSV and JSON

suppressPackageStartupMessages({
  library(devtools)
  devtools::load_all(".")
  library(SummarizedExperiment)
  library(jsonlite)
  library(ggplot2)
})

cat("\n════════════════════════════════════════════════════════════\n")
cat("QUICK TEST: Load Results from CSV & JSON, Then Test Plotting\n")
cat("════════════════════════════════════════════════════════════\n\n")

output_dir <- "/home/nouser/galaxy/tools_source/TSENAT/output"

# STEP 1: Load precomputed results from CSV and JSON
cat("[1] Loading precomputed results...\n")

model_json_file <- file.path(output_dir, "model_data.json")
results_csv_file <- file.path(output_dir, "lm_interaction_results.csv")

if (!file.exists(model_json_file) || !file.exists(results_csv_file)) {
  cat("✗ Files not found. Run TEST_return_model_data.R first.\n")
  quit(status = 1)
}

# Load metadata from JSON
model_data_loaded <- read_json(model_json_file)

# Load results from CSV (cleaner for sorting and inspection)
results_df <- read.csv(results_csv_file, stringsAsFactors = FALSE)

# Sort by adj_p_interaction (most significant first)
results_df <- results_df[order(results_df$adj_p_interaction), ]

# Simulate calculate_lm_interaction output with return_model_data = TRUE
lm_res <- list(
    results = results_df,
    model_data = model_data_loaded$model_metadata
)

cat("  ✓ Loaded model_data.json (metadata)\n")
cat("  ✓ Loaded lm_interaction_results.csv (results)\n")
cat("  - Results rows:", nrow(results_df), "\n")
cat("  - Model groups:", paste(unlist(lm_res$model_data$group_levels), collapse=", "), "\n")
cat("  - Sorted by adj_p_interaction (ascending)\n\n")

# Display top genes
cat("  Top genes by significance:\n")
for (i in 1:min(10, nrow(results_df))) {
  cat(sprintf("    %d. %s (adj_p = %.3e)\n", i, results_df$gene[i], results_df$adj_p_interaction[i]))
}
cat("\n")

# STEP 2: Quickly rebuild diversity SE (needed for plot_lm_interaction_gam)
cat("[2] Rebuilding diversity SE (minimal steps)...\n")

data(readcounts)
readcounts <- as.matrix(salmon_dataset)
mode(readcounts) <- "numeric"

tpm_data <- as.matrix(salmon_tpm)
mode(tpm_data) <- "numeric"

metadata_df <- read.table(
    system.file("extdata", "metadata.tsv", package = "TSENAT"),
    header = TRUE, sep = "\t"
)

gff3_dataset <- system.file("extdata", "annotation.gff3.gz", package = "TSENAT")

se <- build_se(readcounts, gff3_dataset, metadata = metadata_df)
tpm_subset <- tpm_data[rownames(se), colnames(se)]
assay(se, "tpm", withDimnames = FALSE) <- tpm_subset
se <- filter_se(se, stringency = "medium")

ts_se <- calculate_diversity(
  se, q = seq(0.1, 2, by = 0.05), norm = TRUE, 
  metadata = metadata_df
)

cat("  ✓ Diversity SE rebuilt\n\n")

# STEP 2b: Subset diversity SE to only include genes in results
cat("[2b] Subsetting diversity SE to match results CSV...\n")

# After calculate_diversity, rownames(ts_se) are gene names
# results_df$gene column also contains gene names
# So we can match directly!
gene_names_in_results <- results_df$gene
gene_names_in_se <- rownames(ts_se)

cat("  Genes in results CSV:", length(unique(gene_names_in_results)), "\n")
cat("  Genes in diversity SE:", length(gene_names_in_se), "\n")

# Find genes that exist in both
available_genes <- gene_names_in_se[gene_names_in_se %in% gene_names_in_results]
cat("  Genes available in both:", length(available_genes), "\n")

if (length(available_genes) > 0) {
    # Subset SE to only genes in results
    ts_se <- ts_se[available_genes, ]
    cat("  ✓ Subsetted diversity SE to", nrow(ts_se), "genes\n")
} else {
    cat("  WARNING: No genes from results found in diversity SE!\n")
    cat("  First genes in results:", paste(head(gene_names_in_results, 5), collapse=", "), "\n")
    cat("  First genes in SE rownames:", paste(head(gene_names_in_se, 5), collapse=", "), "\n")
}

cat("\n")

# STEP 3: Recreate model_data in proper R format
cat("[3] Recreate model data...\n")

# Extract model_data from lm_res (simulated return_model_data output)
model_data <- list(
  method = lm_res$model_data$method[[1]],
  n_genes = lm_res$model_data$n_genes[[1]],
  n_q_values = lm_res$model_data$n_q_values[[1]],
  q_values = unlist(lm_res$model_data$q_values),
  sample_names = unlist(lm_res$model_data$sample_names),
  group_levels = unlist(lm_res$model_data$group_levels),
  genes_analyzed = unlist(lm_res$model_data$genes_analyzed),
  test_configuration = lm_res$model_data$test_configuration
)

cat("  ✓ Model data prepared\n\n")

# STEP 4: Test plotting with loaded results
cat("[4] Testing plot_lm_interaction_gam with loaded results...\n")

# Filter results to only include genes that are in the SE
results_df_available <- lm_res$results[lm_res$results$gene %in% rownames(ts_se), ]
results_df_available <- results_df_available[order(results_df_available$adj_p_interaction), ]

cat("  Results available in SE:", nrow(results_df_available), "of", nrow(results_df), "\n")

# Show top genes
if (nrow(results_df_available) > 0) {
    cat("  Top genes by adj_p_interaction:\n")
    for (i in 1:min(10, nrow(results_df_available))) {
        cat(sprintf("    %d. %s (adj_p = %.3e)\n", i, results_df_available$gene[i], results_df_available$adj_p_interaction[i]))
    }
    cat("\n")
    
    # Get top 10 genes for plotting
    top_genes <- results_df_available$gene[1:10]
} else {
    cat("  No genes from results found in subsetted SE!\n")
    top_genes <- character(0)
}

if (length(top_genes) > 0) {
    tryCatch({
      plots <- plot_lm_interaction_gam(
          se = ts_se,
          lm_res = results_df_available,
          sample_type_col = "sample_type",
          genes = top_genes,
          model_data = model_data,
          palette = "Set1"
      )
      
      if (is.list(plots)) {
        cat("  ✓ Generated", length(plots), "plots\n")
        
        # Save plots
        for (i in seq_along(plots)) {
          gene_id <- names(plots)[i]
          output_file <- file.path(output_dir, sprintf("quick_test_plot_%d.png", i))
          ggplot2::ggsave(output_file, plots[[i]], width = 10, height = 6, dpi = 300)
          cat(sprintf("    ✓ Saved: %s\n", basename(output_file)))
        }
      } else {
        cat("  ✓ Generated single plot\n")
        output_file <- file.path(output_dir, "quick_test_plot_single.png")
        ggplot2::ggsave(output_file, plots, width = 10, height = 6, dpi = 300)
        cat(sprintf("    ✓ Saved: %s\n", basename(output_file)))
      }
    }, error = function(e) {
      cat("  ✗ Error:", e$message, "\n")
    })
} else {
    cat("  ✗ No genes available for plotting\n")
}

cat("\n════════════════════════════════════════════════════════════\n")
cat("COMPLETE: Loaded results from JSON and generated plots\n")
cat("════════════════════════════════════════════════════════════\n\n")
