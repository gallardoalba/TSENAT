# TSENAT: Tsallis Entropy Analysis Toolbox
# Jackknife Diagnostics for Isoform Switching
# 
# This script contains the complete analysis workflow from the test_jackknife.Rmd vignette

# =============================================================================
# Setup and Installation
# =============================================================================

# Install TSENAT if not already installed
# if (!require("BiocManager", quietly = TRUE))
#     install.packages("BiocManager")
# BiocManager::install("TSENAT")

setwd("~/galaxy/tools_source/TSENAT")

# Install devtools if not already installed
if (!require("devtools", quietly = TRUE))
  install.packages("devtools")

# Load the package from the local folder
devtools::load_all(".")

# Now you can use TSENAT
library(ggplot2)
library(SummarizedExperiment)

# Load packages
suppressPackageStartupMessages({
    library(ggplot2)
    library(SummarizedExperiment)
})

# devtools::install("~/galaxy/tools_source/TSENAT", dependencies = TRUE)
library(TSENAT)

# =============================================================================
# Data Loading and Preprocessing
# =============================================================================

# Load example dataset
data(readcounts)
readcounts <- as.matrix(salmon_dataset)
mode(readcounts) <- "numeric"

# Load sample metadata
metadata_df <- read.table(
    system.file("extdata", "metadata.tsv", package = "TSENAT"),
    header = TRUE, sep = "\t"
)

# Get the GFF3.gz annotation file path
gff3_dataset <- system.file("extdata", "annotation.gff3.gz", package = "TSENAT")

# =============================================================================
# Build SummarizedExperiment
# =============================================================================

# Build a SummarizedExperiment from readcounts + GFF3.gz annotation
se <- build_se(readcounts, gff3_dataset, metadata = metadata_df)

# Examine the structure of the SummarizedExperiment
cat("SummarizedExperiment dimensions:\n")
cat("  Transcripts (rows):", nrow(se), "\n")
cat("  Samples (columns):", ncol(se), "\n")
cat("  Assay (counts):", nrow(assay(se, "counts")), "×", ncol(assay(se, "counts")), "\n")
cat("\nRowData columns:", paste(colnames(rowData(se)), collapse = ", "), "\n")
cat("Sample metadata available:", exists("metadata_df"), "\n\n")

# Show rowData with gene information
rd <- as.data.frame(rowData(se)[1:5, ])

# =============================================================================
# Filter lowly-expressed transcripts
# =============================================================================

# Filter lowly-expressed transcripts and report counts
# Use stringency = "medium" for balanced filtering
se <- filter_se(se, stringency = "medium")

# =============================================================================
# Compute Tsallis entropy
# =============================================================================

# Compute Tsallis entropy for a sequence of values (normalized)
qvec <- seq(0.1, 2, by = 0.05)
ts_se <- calculate_diversity(se,
    q = qvec, norm = TRUE, metadata = metadata_df
)

# =============================================================================
# Plot q-curve
# =============================================================================

# Plot overall q-curve (mean entropy across all genes)
p_qcurve <- plot_tsallis_q_curve(ts_se)
print(p_qcurve)

# =============================================================================
# Quality Control: Sample Influence Assessment
# =============================================================================

# Perform multi-q sample influence analysis
sample_qc <- m_estimate(ts_se, samples = "sample_type", 
                        loss_type = "huber", 
                        q_combine_method = "mean",
                        influence_threshold = 0.75,
                        paired = TRUE)

# Display QC results
base_display_cols <- c("Sample", "Condition", "Proportion_Affected", 
                       "Genes_Affected", "Distance_from_Centroid", "Status")
qc_table <- sample_qc[, base_display_cols]

# =============================================================================
# Linear-model interaction test
# =============================================================================

# Test for interactions between q and sample groups
if (requireNamespace("mgcv", quietly = TRUE)) {
    lm_res <- calculate_lm_interaction(ts_se,
        method = "gam",
        paired = TRUE,
        multicorr = "hochberg",
        nthreads = 3,
        verbose = TRUE
    )
}

# =============================================================================
# Plot q-curve profiles for top genes
# =============================================================================

# Plot q-curve profile for LINC03040
plot_target <- plot_tsallis_gene_profile(ts_se, gene = "LINC03040")
print(plot_target)

# Plot q-curve profile for PNPT1
plot_target <- plot_tsallis_gene_profile(ts_se, gene = "PNPT1")
print(plot_target)

# Plot q-curve profile for CXCL12
plot_target <- plot_tsallis_gene_profile(ts_se, gene = "CXCL12")
print(plot_target)

# =============================================================================
# Effect size analysis across q-spectrum
# =============================================================================

if (requireNamespace("mgcv", quietly = TRUE) && exists("lm_res") && nrow(lm_res) > 0) {
  # Set rownames of lm_res to gene names (required by calculate_divergence)
  rownames(lm_res) <- lm_res$gene
  
  # Calculate Tsallis divergence with bootstrap CIs
  divergence_results_se <- calculate_divergence(
      se = se,
      group_col = "sample_type",
      control_group = "normal",
      q = qvec,
      bootstrap = TRUE,
      nboot = 50,
      ci = 0.95,
      method = "percentile"
  )
}

# Merge LMM results with divergence effect sizes
if (exists("divergence_results_se") && nrow(divergence_results_se) > 0) {
    eff_res <- effect_sizes_divergence(
        lm_res = lm_res,
        divergence_results_se = divergence_results_se,
        significance_threshold = 0.05,
        enrich_per_q_pattern = TRUE,
        verbose = FALSE
    )
}

# =============================================================================
# Visualization: Divergence distribution
# =============================================================================

if (exists("eff_res") && 
    !is.null(eff_res) && 
    !is.null(eff_res$interaction_results) &&
    nrow(eff_res$interaction_results) > 0) {
    
    p <- plot_divergence_distribution(eff_res$interaction_results, 
        threshold = 0.1)
    print(p)
}

# =============================================================================
# Classify gene pattern types
# =============================================================================

# Extract classification table
classification_tbl <- extract_classification_table(eff_res, top_n = 10, sort_by = 'padj')

# =============================================================================
# Visualize Q-Spectrum Curve
# =============================================================================

if (exists("divergence_results_se")) {
    plot_gene_q_spectrum(divergence_results_se, target_gene = "LINC03040", verbose = FALSE)
}

# =============================================================================
# Compare q-spectra across multiple genes
# =============================================================================

plot_generated <- FALSE

if (exists("eff_res") && 
    !is.null(eff_res) &&
    !is.null(eff_res$interaction_results) &&
    nrow(eff_res$interaction_results) > 0) {
  
  # Check if patchwork is available (required for plot_multi_gene_q_spectrum)
  if (!requireNamespace("patchwork", quietly = TRUE)) {
    cat("Note: The patchwork package is required for plot_multi_gene_q_spectrum().\n")
    cat("Install it with: install.packages('patchwork')\n")
    cat("Skipping multi-gene q-spectrum plot for now.\n\n")
  } else {
    tryCatch({
      p_multi <- plot_multi_gene_q_spectrum(
        eff_res = eff_res, 
        n_genes = 9, 
        ncol = 3,
        verbose = TRUE
      )
      
      if (!is.null(p_multi)) {
        print(p_multi)
        plot_generated <- TRUE
      }
    }, error = function(e) {
      # Print error message for debugging
      cat("Note: plot_multi_gene_q_spectrum() encountered an error:\n")
      cat(conditionMessage(e), "\n")
      cat("This function may not be available in the current TSENAT version.\n")
    })
  }
}

# =============================================================================
# Analysis 1: Basic Jackknife Stability Analysis
# =============================================================================

# Select genes with significant q×condition interaction
# significant_lm_genes <- lm_res$gene[lm_res$adj_p_interaction < 0.05]

#if (length(significant_lm_genes) == 0) {
#  cat("No genes with FDR < 0.05 in LM test. Using top LM genes for demonstration.\n")
#  significant_lm_genes <- lm_res$gene[1:min(3, nrow(lm_res))]
#}

#cat("Genes with significant q×condition interaction:", length(significant_lm_genes), "\n\n")

# Find the first LM-significant gene with multiple transcripts
# selected_gene <- NULL
# sample_col <- 1
 
# for (gene_candidate in lm_res$gene[lm_res$adj_p_interaction < 0.05]) {
#   gene_counts_test <- assays(se)$counts[rowData(se)$gene_id == gene_candidate, sample_col]
#   gene_counts_test <- gene_counts_test[gene_counts_test > 0]
# 
#   if (length(gene_counts_test) >= 2) {
#     selected_gene <- gene_candidate
#     gene_counts <- gene_counts_test
#     break
#   }
# }
# 
# Fallback to any gene with 2+ transcripts if needed
# if (is.null(selected_gene)) {
#   cat("Finding a gene with multiple isoforms for demonstration...\n\n")
#   gene_list <- unique(rowData(se)$gene_id)
#   for (gene_candidate in gene_list[1:min(20, length(gene_list))]) {
#     gene_counts_test <- assays(se)$counts[rowData(se)$gene_id == gene_candidate, sample_col]
#     gene_counts_test <- gene_counts_test[gene_counts_test > 0]
# 
#     if (length(gene_counts_test) >= 2) {
#       selected_gene <- gene_candidate
#       gene_counts <- gene_counts_test
#       break
#     }
#   }
# }


# =============================================================================
# Analysis 2: Multi-Q Sensitivity Analysis
# =============================================================================

# Run jackknife for multiple q values on the same gene
# q_values <- c(0.5, 1, 1.5, 2)
# jack_multiq <- jackknife_tsallis_entropy(
#   x = gene_counts,
#   q = q_values,
#   norm = TRUE,
#   print_results = FALSE
# )
# 
# # Create summary table
# q_summary <- data.frame(
#   q = q_values,
#   entropy = sapply(jack_multiq, function(x) x$estimate),
#   jackknife_se = sapply(jack_multiq, function(x) x$jackknife_se),
#   n_outliers = sapply(jack_multiq, function(x) length(x$outlier_indices))
# )
# 
# cat("Multi-Q Summary for", selected_gene, ":\n")
# print(knitr::kable(q_summary, digits = 4))
# 
# cat("\nInterpretation:\n")
# cat("- q=0.5: Emphasizes rare transcripts\n")
# cat("- q=1: Shannon entropy (balanced)\n")
# cat("- q=1.5: Slight emphasis on abundant transcripts\n")
# cat("- q=2: Strong emphasis on abundant transcripts\n")
# cat("- Notice how entropy changes and outlier count varies with q\n")

# =============================================================================
# Analysis 3: Condition-Stratified Isoform Switching
# =============================================================================
# 
# # Identify the condition column name
# condition_col <- "sample_type"
# 
# # Identify the paired samples column if available
# pair_col <- NULL
# coldata_cols <- colnames(colData(se))
# if ("paired_samples" %in% coldata_cols) {
#   pair_col <- "paired_samples"
# }
# 
# # Map gene names to Ensembl IDs for jackknife compatibility
# rd <- rowData(se)
# gene_name_to_id <- setNames(as.character(rd$gene_id), as.character(rd$gene_name))
# gene_name_to_id <- gene_name_to_id[!is.na(names(gene_name_to_id))]
# gene_name_to_id <- gene_name_to_id[!duplicated(names(gene_name_to_id))]
# 
# # Create a copy of lm_res with gene IDs for jackknife
# lm_res_for_jackknife <- lm_res
# lm_res_for_jackknife$gene <- unname(gene_name_to_id[as.character(lm_res$gene)])
# lm_res_for_jackknife <- lm_res_for_jackknife[!is.na(lm_res_for_jackknife$gene), ]
# 
# n_lm_sig <- sum(lm_res_for_jackknife$adj_p_interaction < 0.05)
# 
# # Run isoform switching analysis
# switching_results <- jackknife_isoform_switching(
#   se = se,
#   condition_col = condition_col,
#   pair_col = pair_col,
#   gene_col = "gene_id",
#   isoform_col = "transcript_id",
#   q = 1,
#   norm = TRUE,
#   n_bootstrap = 100,
#   print_results = FALSE,
#   threshold = 90,
#   top_n = n_lm_sig,
#   lm_results = lm_res_for_jackknife,
#   lm_p_threshold = 0.05,
#   use_lm_fdr = TRUE
# )
# 
# # Prepare results for display
# summary_df <- switching_results$summary_table
# 
# # Sort by LM test adjusted p-value (most significant first)
# lm_pval_map <- setNames(lm_res$adj_p_interaction, lm_res$gene)
# summary_df$lm_adj_pval <- lm_pval_map[summary_df$gene_name]
# summary_df <- summary_df[order(summary_df$lm_adj_pval, na.last = TRUE), ]

# =============================================================================
# Analysis 3b: Multi-Q Switching Validation
# =============================================================================

# Identify the paired samples column if available
pair_col <- NULL
coldata_cols <- colnames(colData(se))
if ("paired_samples" %in% coldata_cols) {
  pair_col <- "paired_samples"
}

# Identify the condition column name
condition_col <- "sample_type"

# Map gene names to Ensembl IDs for jackknife compatibility
rd <- rowData(se)
gene_name_to_id <- setNames(as.character(rd$gene_id), as.character(rd$gene_name))
gene_name_to_id <- gene_name_to_id[!is.na(names(gene_name_to_id))]
gene_name_to_id <- gene_name_to_id[!duplicated(names(gene_name_to_id))]

# Create a copy of lm_res with gene IDs for jackknife
lm_res_for_jackknife <- lm_res
lm_res_for_jackknife$gene <- unname(gene_name_to_id[as.character(lm_res$gene)])
lm_res_for_jackknife <- lm_res_for_jackknife[!is.na(lm_res_for_jackknife$gene), ]

n_lm_sig <- sum(lm_res_for_jackknife$adj_p_interaction < 0.05)

# Multi-Q analysis: Test switching at different q values
q_test_values <- c(0.01, 0.5, 1.0, 1.5, 2.0)

# Use top 4 LM-significant genes (by adj p-value)
lm_res_sorted <- lm_res_for_jackknife[order(lm_res_for_jackknife$adj_p_interaction), ]
top_genes_for_comparison <- head(lm_res_sorted$gene, 4)

# Single call for all q values at once (vectorized approach)
multi_q_results <- jackknife_isoform_switching(
  se = se,
  condition_col = condition_col,
  pair_col = pair_col,
  gene_col = "gene_id",
  isoform_col = "transcript_id",
  q = q_test_values,
  norm = TRUE,
  n_bootstrap = 100,
  print_results = TRUE,
  threshold = 90,
  top_n = n_lm_sig,
  lm_results = lm_res_for_jackknife,
  lm_p_threshold = 0.05,
  use_lm_fdr = TRUE
)

# =============================================================================
# Analysis 4: Visualization - Multi-Q Heatmaps
# =============================================================================

# Prepare data for combined multi-Q heatmap
all_gene_matrices <- list()
all_gene_info <- list()

for (gene_idx in 1:min(6, length(top_genes_for_comparison))) {
  gene_id <- top_genes_for_comparison[gene_idx]
  gene_name <- lm_res_sorted[lm_res_sorted$gene == gene_id, "gene_name"]
  if (length(gene_name) == 0) gene_name <- gene_id
  
  # Collect delta_influence for all transcripts across all q-values
  heatmap_data <- NULL
  
  # Iterate over the multi_q_results keys (q_0_01, q_0_50, etc.)
  for (q_key in names(multi_q_results)) {
    if (!is.null(multi_q_results[[q_key]]) && 
        !is.null(multi_q_results[[q_key]]$results_per_gene) &&
        gene_id %in% names(multi_q_results[[q_key]]$results_per_gene)) {
      gene_res <- multi_q_results[[q_key]]$results_per_gene[[gene_id]]
      if (!is.null(gene_res$delta_influence)) {
        # Extract q value from key name (e.g., "q_1_00" -> "q_1.00")
        col_display_name <- gsub("_", ".", q_key)
        delta_vals <- gene_res$delta_influence
        
        # Replace Inf and NaN with NA for clean handling
        delta_vals[is.infinite(delta_vals) | is.nan(delta_vals)] <- NA
        
        if (is.null(heatmap_data)) {
          # Initialize with transcript IDs
          heatmap_data <- data.frame(
            transcript = gene_res$transcript_ids,
            stringsAsFactors = FALSE
          )
        }
        
        # Add column for this q-value (with Inf/NaN as NA)
        # Match length to existing data rows
        n_to_add <- min(length(delta_vals), nrow(heatmap_data))
        heatmap_data[[col_display_name]] <- c(delta_vals[1:n_to_add], rep(NA_real_, nrow(heatmap_data) - n_to_add))
      }
    }
  }
  
  if (!is.null(heatmap_data) && nrow(heatmap_data) > 0) {
    # Convert to matrix for heatmap
    heatmap_matrix <- as.matrix(heatmap_data[, -1])
    rownames(heatmap_matrix) <- heatmap_data$transcript
    
    # Remove transcripts with exclusively no valid data (all NA)
    valid_transcript_rows <- rowSums(!is.na(heatmap_matrix)) > 0
    heatmap_matrix <- heatmap_matrix[valid_transcript_rows, ]
    heatmap_data <- heatmap_data[valid_transcript_rows, ]
    
    if (nrow(heatmap_matrix) == 0) {
      next
    }
    
    # Cap outliers for remaining finite values
    finite_vals <- heatmap_matrix[is.finite(heatmap_matrix)]
    if (length(finite_vals) > 0) {
      cap_val <- quantile(abs(finite_vals), 0.95)
      # Cap finite values that exceed the 95th percentile
      heatmap_matrix[is.finite(heatmap_matrix) & abs(heatmap_matrix) > cap_val] <- 
        sign(heatmap_matrix[is.finite(heatmap_matrix) & abs(heatmap_matrix) > cap_val]) * cap_val
    }
    
    # Transpose: q-values as rows, transcripts as columns
    heatmap_matrix <- t(heatmap_matrix)
    
    # Store matrix and gene info for combined plot
    all_gene_matrices[[gene_idx]] <- heatmap_matrix
    all_gene_info[[gene_idx]] <- list(
      gene_id = gene_id,
      gene_name = gene_name,
      n_transcripts = ncol(heatmap_matrix)
    )
  }
}

# Create individual heatmaps for each gene
if (length(all_gene_matrices) > 0 && requireNamespace("pheatmap", quietly = TRUE)) {
  for (gene_idx in seq_along(all_gene_matrices)) {
    if (is.null(all_gene_matrices[[gene_idx]]) || is.null(all_gene_info[[gene_idx]])) {
      next
    }
    
    mat <- all_gene_matrices[[gene_idx]]
    gene_info <- all_gene_info[[gene_idx]]
    gene_name <- gene_info$gene_name
    
    # Apply outlier capping
    finite_vals <- abs(as.vector(mat[is.finite(mat)]))
    if (length(finite_vals) == 0) {
      percentile_95 <- 1  # Default if no finite values
    } else {
      percentile_95 <- quantile(finite_vals, 0.95)
    }
    
    mat_viz <- mat
    if (length(finite_vals) > 0) {
      mat_viz[is.finite(mat_viz) & abs(mat_viz) > percentile_95] <- 
        sign(mat_viz[is.finite(mat_viz) & abs(mat_viz) > percentile_95]) * percentile_95
    }
    
    # Create pheatmap
    p <- pheatmap::pheatmap(
      mat_viz,
      main = gene_name,
      cluster_rows = FALSE,
      cluster_cols = (ncol(mat_viz) > 1),
      display_numbers = FALSE,
      na_col = "lightgray",
      color = colorRampPalette(c("#4575B4", "#FFFFFF", "#D73027"))(100),
      cellwidth = 65,
      cellheight = 65,
      fontsize = 16,
      fontsize_row = 16,
      fontsize_col = 16,
      fontsize_number = 13,
      margins = c(11, 10),
      show_rownames = TRUE,
      show_colnames = TRUE,
      silent = TRUE
    )
    
    print(p)
  }
}

# Debug commands
print(dim(mat_viz))           # Check dimensions
print(range(mat_viz, na.rm=TRUE))  # Check value range
print(sum(is.na(mat_viz)))    # Count NAs
print(sum(!is.finite(mat_viz)))  # Count non-finite values

# =============================================================================
# Basic Isoform Switching Heatmap
# =============================================================================
# 
# # Create heatmap showing switching patterns in top LM-significant genes
# top_n_genes <- min(4, nrow(lm_res))
# top_lm_genes_names <- lm_res$gene_name[1:top_n_genes]
# 
# # Create a filtered version of switching results
# switching_filtered <- switching_results
# switching_filtered$summary_table <- switching_results$summary_table[
#   switching_results$summary_table$gene_name %in% top_lm_genes_names,
# ]
# 
# # Filter the per-gene results as well
# top_lm_genes_ids <- switching_results$summary_table[
#   switching_results$summary_table$gene_name %in% top_lm_genes_names,
# ]$gene
# 
# switching_filtered$results_per_gene <- switching_results$results_per_gene[
#   names(switching_results$results_per_gene) %in% top_lm_genes_ids
# ]
# 
# # Check if there's any switching to visualize
# max_switching <- max(abs(switching_filtered$summary_table$max_delta_influence), na.rm = TRUE)
# 
# if (max_switching > 0.001) {
#   tryCatch({
#     suppressWarnings({
#       max_tx_per_gene <- max(sapply(switching_filtered$results_per_gene, 
#                                     function(x) length(x$transcript_ids)))
#       
#       plot_isoform_switching_heatmap(
#         switching_filtered,
#         top_n = nrow(switching_filtered$summary_table),
#         top_transcripts_per_gene = max_tx_per_gene,
#         color_scheme = "diverging",
#         show_values = TRUE,
#         main = "Isoform Switching in Top LM-Significant Genes\n(Red=up in first condition, Blue=up in second)"
#       )
#     })
#   }, error = function(e) {
#     cat("Note: Heatmap visualization encountered an issue.\n")
#     cat("This may occur if selected genes have minimal switching patterns.\n")
#     cat("Summary statistics are still available.\n")
#   })
# } else {
#   cat("Note: No substantial isoform switching detected in selected genes.\n")
# }

# =============================================================================
# Analysis 5: Q-Parameter Sensitivity for Specific Gene
# =============================================================================
# 
# # Select a gene with significant q×condition interaction
# significant_lm_genes_names <- lm_res$gene_name[lm_res$adj_p_interaction < 0.05]
# 
# if (length(significant_lm_genes_names) > 0) {
#   # Find the first LM-significant gene that exists in the filtered SE dataset
#   selected_gene_for_q <- NULL
#   genes_in_se <- unique(rowData(se)$gene_id)
#   
#   # Get gene IDs that match the selected gene names
#   for (lm_gene_name in significant_lm_genes_names) {
#     # Find the gene ID for this gene name
#     matching_gene_id <- rowData(se)$gene_id[
#       which(rowData(se)$gene_name == lm_gene_name)[1]
#     ]
#     if (!is.na(matching_gene_id)) {
#       selected_gene_for_q <- matching_gene_id
#       break
#     }
#   }
#   
#   if (!is.null(selected_gene_for_q)) {
#     cat("Analyzing Q-sensitivity for gene:", selected_gene_for_q, "\n")
#     cat("This gene shows a significant q×condition interaction pattern\n\n")
#     
#     tryCatch({
#       # Plot q-sensitivity for this gene
#       suppressWarnings({
#         plot_q_sensitivity_curve(
#           se = se,
#           condition_col = condition_col,
#           gene = selected_gene_for_q,
#           gene_col = "gene_id",
#           isoform_col = "transcript_id",
#           q_values = c(0.5, 0.8, 1, 1.2, 1.5, 1.8, 2),
#           norm = TRUE,
#           show_legend = TRUE
#         )
#       })
#       
#       cat("\nQ-Sensitivity Curve Interpretation for", selected_gene_for_q, ":\n")
#       cat("- Red/colored lines: Entropy for first and second conditions\n")
#       cat("- The separation between these lines shows condition-specific diversity\n")
#       cat("- Variation with q reveals scale-dependent condition effects\n")
#       cat("- Gray shaded area: Recommended q ∈ [0.5, 2] range\n")
#       cat("- If lines diverge: Effect consistent across abundance scales\n")
#       cat("- If lines converge: Effect is scale-dependent\n")
#     }, error = function(e) {
#       cat("Note: Could not generate q-sensitivity plot for this gene\n")
#       cat("This can occur if the gene has insufficient transcript variety\n")
#       cat("Error:", conditionMessage(e), "\n")
#     })
#   } else {
#     cat("Note: No LM-significant genes found in the filtered dataset\n")
#     cat("The filtered SE may contain different genes than the full dataset used for LM test\n")
#   }
# } else {
#   cat("Note: No genes with significant q×condition interaction detected\n")
#   cat("This can occur with real data where scale-dependent effects are subtle\n")
# }

# =============================================================================
# Analysis 6: Block Jackknife for Grouped Data
# =============================================================================
# 
# # Filter to LM-significant genes
# significant_lm_genes_names <- lm_res$gene_name[lm_res$adj_p_interaction < 0.05]
# 
# if (length(significant_lm_genes_names) == 0) {
#   cat("Note: No genes reached FDR < 0.05 in LM test.\n")
#   cat("Using top LM-significant genes for block jackknife analysis.\n\n")
#   significant_lm_genes_names <- lm_res$gene_name[1:min(3, nrow(lm_res))]
# }
# 
# cat("Block jackknife analysis for LM-significant genes:\n")
# cat("Genes analyzed:", paste(head(significant_lm_genes_names, min(3, length(significant_lm_genes_names))), collapse = ", "), "\n")
# cat("Total LM-significant genes:", length(significant_lm_genes_names), "\n\n")
# 
# # Identify available grouping columns
# group_cols <- colnames(colData(se))
# 
# # Use pairing information if available
# block_column <- if ("paired_samples" %in% group_cols) {
#   "paired_samples"
# } else if ("pair" %in% tolower(group_cols)) {
#   grep("pair", colnames(colData(se)), ignore.case = TRUE, value = TRUE)[1]
# } else if ("batch" %in% tolower(group_cols)) {
#   grep("batch", colnames(colData(se)), ignore.case = TRUE, value = TRUE)[1]
# } else {
#   condition_col  # Fallback to condition column
# }
# 
# cat("Using column '", block_column, "' for block jackknife grouping\n", sep = "")
# cat("Groups:", paste(unique(colData(se)[[block_column]]), collapse = ", "), "\n\n")
# 
# # Block jackknife analysis interpretation
# if (length(significant_lm_genes_names) > 0) {
#   first_lm_gene_name <- significant_lm_genes_names[1]
#   
#   cat("\nBlock Jackknife Results for", first_lm_gene_name, ":\n")
#   cat("- Analyzes how entropy changes when entire", block_column, "groups are removed\n")
#   cat("- Useful for assessing group/pairing/batch effects\n")
#   cat("- Shows how each group contributes to overall diversity in this gene\n")
# }
# 
# cat("\nBlock Jackknife Interpretation:\n")
# cat("- Groups with high influence warrant further investigation\n")
# cat("- Helps distinguish block/batch/pairing effects from biological signal\n")
# cat("- Essential for validating paired design results\n")
# 
# # =============================================================================
# # Analysis 7: Summary Statistics and Key Findings
# # =============================================================================
# 
# # Filter switching results to LM-significant genes
# significant_lm_genes_names <- lm_res$gene_name[lm_res$adj_p_interaction < 0.05]
# 
# if (length(significant_lm_genes_names) == 0) {
#   cat("Note: No genes reached FDR < 0.05 in LM test\n")
#   cat("Showing results for all analyzed genes\n\n")
#   summary_df <- switching_results$summary_table
# } else {
#   summary_df <- switching_results$summary_table[
#     switching_results$summary_table$gene_name %in% significant_lm_genes_names,
#   ]
#   summary_df <- summary_df[order(-summary_df$max_delta_influence), ]
# }
# 
# cat("Summary of Jackknife Results for LM-Significant Genes:\n\n")
# 
# if (nrow(summary_df) > 0) {
#   cat("Summary Table:\n")
#   summary_display <- summary_df[1:min(10, nrow(summary_df)),
#                                 c("gene", "gene_name", "n_transcripts", "n_switching_transcripts",
#                                   "max_delta_influence", "n_fdr_significant")]
#   print(summary_display)
#   
#   cat("\n\nKey Findings:\n")
#   cat("Total LM-significant genes analyzed:", nrow(summary_df), "\n")
#   cat("Genes with switching transcripts:", sum(summary_df$n_switching_transcripts > 0), "\n")
#   cat("Genes with FDR-significant switching (FDR < 0.05):",
#       sum(summary_df$n_fdr_significant > 0), "\n\n")
#   
#   # Show top genes
#   top_genes <- summary_df[summary_df$n_switching_transcripts > 0, ]
#   if (nrow(top_genes) > 0) {
#     cat("Top genes by switching magnitude:\n")
#     for (i in 1:min(5, nrow(top_genes))) {
#       gene <- top_genes$gene_name[i]
#       n_switch <- top_genes$n_switching_transcripts[i]
#       n_sig <- top_genes$n_fdr_significant[i]
#       max_delta <- top_genes$max_delta_influence[i]
#       
#       cat("  ", gene, ": ", n_switch, " switching transcripts,",
#           n_sig, " FDR-significant, max delta = ", round(max_delta, 3), "\n", sep = "")
#     }
#   }
#   
#   # Access transcripts with significant FDR
#   cat("\n\nTranscripts with FDR < 0.05 in LM-Significant Genes:\n")
#   sig_transcript_genes <- switching_results$summary_table[
#     switching_results$summary_table$gene_name %in% significant_lm_genes_names,
#   ]$gene
#   
#   sig_transcripts <- switching_results$all_transcript_stats[
#     switching_results$all_transcript_stats$fdr < 0.05 &
#     switching_results$all_transcript_stats$gene %in% sig_transcript_genes,
#   ]
#   
#   if (nrow(sig_transcripts) > 0) {
#     sig_display <- sig_transcripts[1:min(20, nrow(sig_transcripts)),
#                                    c("gene", "transcript_id", "delta_influence",
#                                      "pvalue", "fdr")]
#     print(sig_display)
#     cat("\nTotal:", nrow(sig_transcripts), "transcripts with FDR < 0.05\n")
#   } else {
#     cat("None found at FDR < 0.05\n")
#     cat("(This is common with real data showing subtle switching effects)\n")
#   }
# } else {
#   cat("No LM-significant genes with switching detected\n")
# }

# =============================================================================
# Session Information
# =============================================================================

sessionInfo()
