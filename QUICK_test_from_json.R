#!/usr/bin/env Rscript
# Load precomputed results from CSV and JSON

suppressPackageStartupMessages({
  library(devtools)
  devtools::load_all(".")
  library(SummarizedExperiment)
  library(jsonlite)
  library(ggplot2)
  library(cowplot)
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
cat("  - Results rows:", nrow(lm_res$results), "\n")
cat("  - Model groups:", paste(unlist(lm_res$model_data$group_levels), collapse=", "), "\n")
cat("  - Sorted by adj_p_interaction (ascending)\n\n")

# Display top genes from loaded results
cat("  Top genes by significance:\n")
for (i in 1:min(10, nrow(lm_res$results))) {
  cat(sprintf("    %d. %s (adj_p = %.3e)\n", i, lm_res$results$gene[i], lm_res$results$adj_p_interaction[i]))
}
cat("\n")

# STEP 2: Quickly rebuild diversity SE (needed for plot_lm_interaction_gam)
cat("[2] Rebuilding diversity SE (minimal steps)...\n")

# Load the TSENAT example datasets
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

# STEP 3: Generate plots
cat("[3] Generating plots with plot_lm_interaction_gam...\n")
cat("  Note: Function automatically extracts model_data and handles gene matching\n")
cat("  Note: grid=TRUE arranges plots automatically\n\n")

# Generate plots with error handling
# Function will automatically:
# - Extract model_data from lm_res if needed
# - Match genes between SE and results
# - Subset both to only available genes
# - If grid=TRUE, arrange in grid and return single combined plot
tryCatch({
    combined_plot <- plot_lm_interaction_gam(
        se = ts_se,
        lm_res = lm_res,
        n_top = 4  # Pass full lm_res list with $results and $model_data
    )

    if (!is.null(combined_plot)) {
        cat("  ✓ Generated combined plot grid (2 cols × 3 rows)\n\n")
        
        # Save combined grid plot
        output_file <- file.path(output_dir, "quick_test_plot_grid.png")
        ggplot2::ggsave(output_file, combined_plot, width = 14, height = 12, dpi = 300)
        cat(sprintf("  ✓ Saved combined grid: %s\n", basename(output_file)))
        
    } else {
        cat("  ⚠ plot_lm_interaction_gam returned NULL, skipping save\n")
    }
}, error = function(e) {
    cat("  ✗ Error generating plots:", e$message, "\n")
    # Uncomment for debugging:
    # traceback()
})

cat("\n════════════════════════════════════════════════════════════\n")
cat("COMPLETE: Loaded results from JSON and generated plots\n")
cat("════════════════════════════════════════════════════════════\n\n")
