#' Merge LMM Interaction Results with Tsallis Divergence Effect Sizes
#'
#' Combines LMM interaction test p-values with pre-computed Tsallis divergence
#' effect sizes and bootstrap confidence intervals across ALL q values. This is a 
#' **data merger**, not a model fitter--all statistical computation happens upstream in:
#' - `.calculate_lm_interaction()` -> LMM p-values
#' - `.calculate_divergence()` -> Divergence estimates and CIs for multiple q
#'
#' This function merges the results into a single data frame for downstream
#' interpretation. When multiple q values are present, effect sizes are computed
#' for each q to capture the full biological spectrum (rare->abundant isoforms).
#'
#' **Architecture:**
#' ```
#' Input 1: LMM results (from calculate_lm_interaction)
#'   - gene names
#'   - adj_p_interaction values
#'
#' Input 2: Divergence SE (from calculate_divergence)
#'   - gene names
#'   - divergence estimates and bootstrap CIs for each q value
#'
#' Output: Merged data frame
#'   - gene name
#'   - statistical significance (p-value)
#'   - effect magnitude for EACH q: D_q, lower_ci_q, upper_ci_q
#' ```
#'
#' @param lm_res A data frame of LMM interaction test results from `.calculate_lm_interaction()`,
#'   with columns: `gene` (character, gene name), `adj_p_interaction` (numeric, multiple-test adjusted p-value).
#'   Genes with adj_p_interaction below `significance_threshold` are included.
#'
#' @param divergence_results_se A SummarizedExperiment from `.calculate_divergence()`,
#'   containing rowData with columns: `gene_name` and either:
#'   - Generic: `estimate`, `lower_ci`, `upper_ci` (single q-value results), OR
#'   - Per-q: `estimate_q*`, `lower_ci_q*`, `upper_ci_q*` (multiple q-values)
#'   All divergence-related data is self-contained in this object.
#'
#' @param significance_threshold Numeric; p-value threshold for filtering significant
#'   genes (default: 0.05). Only genes with adj_p_interaction < threshold are included.
#'
#' @param enrich_per_q_pattern Logical; if TRUE (default), adds a 'per_q_pattern' column
#'   to the output data frame containing comma-separated divergence values across the
#'   q spectrum for each gene. This column enables visualization and classification
#'   of whether treatment effects are driven by rare (low-q) or abundant (high-q)
#'   isoforms. Set to FALSE to reduce output size if this annotation is not needed.
#'
#' @param verbose Logical; if TRUE, print detailed validation and merge statistics
#'   to console (default: TRUE). Shows counts of passed, skipped, and failed genes.
#'
#' @return A list with two elements:
#'   \describe{
#'     \item{\code{interaction_results}}{Data frame (genes * columns) with merged results:
#'       - `gene`: Gene name (character)
#'       - `p_value_interaction`: LMM adjusted p-value for q:group interaction
#'       - `slope_diff`: q:group interaction slope coefficient (if present in lm_res)
#'       - For EACH q value found: 
#'         - `effect_size_D_q*`: Absolute Tsallis divergence at q
#'         - `D_q*_lower_ci`: Bootstrap lower confidence bound
#'         - `D_q*_upper_ci`: Bootstrap upper confidence bound
#'     }
#'     \item{\code{validation_stats}}{List with merge quality metrics:
#'       - `total_genes`: Total significant genes from LMM
#'       - `passed_lmm`: Successfully merged with divergence data
#'       - `failed_missing_divergence`: Missing or NA divergence estimate
#'       - `other_errors`: Other processing failures
#'       - `q_values`: Numeric vector of q-values processed
#'     }
#'   }
#'
#' @details
#' **Interpretation of Effect Sizes:**
#'
#' Tsallis divergence D_q quantifies the information-theoretic distance between
#' control and treatment isoform distributions at each q-value:
#' - D_q > 0.05: Small effect size
#' - D_q > 0.10: Medium effect size (meaningful biological significance)
#' - D_q > 0.20: Large effect size
#'
#' **Different q-values capture different biological scales:**
#' - q=0.5: Rare (low-abundance) isoforms dominate
#' - q=1.0: Shannon entropy (balanced across abundances)
#' - q=2.0: Common (high-abundance) isoforms dominate
#'
#' Examining the divergence spectrum across q reveals whether treatment effects
#' are driven by rare transcripts (high D at low q) or abundant transcripts (high D at high q).
#'
#' For paired designs, divergence is computed separately within each pair,
#' then averaged to account for pairing structure.
#'
#' **Database References:**
#' - Papers I002-I004: Tsallis divergence mathematical foundation
#' - Papers C016: Bootstrap CI computation respecting data structure
#' - Papers S197: Quality filtering and effect size thresholds
#'

#' @noRd

.effect_sizes_divergence <- function(
    lm_res,
    divergence_results_se,
    significance_threshold = 0.05,
    enrich_per_q_pattern = TRUE,
    verbose = FALSE) {

  # =========================================================================
  # INPUT VALIDATION
  # =========================================================================

  .validateEffectSizeInputs(lm_res, divergence_results_se)
  rd <- SummarizedExperiment::rowData(divergence_results_se)

  # =========================================================================
  # ALIGN DATASETS
  # =========================================================================
  
  alignment_result <- .alignGeneDatasets(lm_res, rd, verbose)
  lm_res <- alignment_result$lm_res
  divergence_genes <- alignment_result$divergence_genes
  q_values <- alignment_result$q_values
  use_generic <- alignment_result$use_generic

  # =========================================================================
  # FILTERING & INITIALIZATION
  # =========================================================================

  # Filter to genes with valid p-values and adj_p_interaction < threshold
  valid_p_idx <- !is.na(lm_res$adj_p_interaction)
  significant_idx <- valid_p_idx & (lm_res$adj_p_interaction < significance_threshold)
  significant_genes <- lm_res$gene[significant_idx]

  if (verbose) {
    message("\n**Filtering effect size analysis to significant genes:**")
    message("- Genes with valid p-values: ", sum(valid_p_idx))
    message("- Genes with adj_p_interaction <", significance_threshold, ": ", 
        length(significant_genes))
  }

  # Create empty results data frame
  interaction_results <- .createResultsDataFrame(q_values, use_generic)

  if (length(significant_genes) == 0) {
    if (verbose) {
      message("No genes with significant q*group interaction detected.")
    }
    return(list(
      interaction_results = interaction_results,
      validation_stats = list(
        total_genes = 0,
        passed_lmm = 0,
        failed_missing_divergence = 0,
        other_errors = 0,
        q_values = q_values
      )
    ))
  }

  # =========================================================================
  # MERGE RESULTS
  # =========================================================================

  # For each significant gene, combine LMM p-value with divergence effect sizes
  validation_stats <- list(
    total_genes = length(significant_genes),
    passed_lmm = 0,
    failed_missing_divergence = 0,
    other_errors = 0,
    q_values = q_values
  )

  if (verbose) {
    message("\nMerging LMM results with divergence effect sizes...")
  }

  # Get matching strategy info
  use_gene_name_col <- "gene_name" %in% colnames(lm_res)
  
  if (verbose) {
    message("[effect_sizes_divergence] Gene name matching strategy:")
    message("  - gene_name column in lm_res:", use_gene_name_col)
    if (use_gene_name_col) {
      message("  - lm_res$gene (first 5):", paste(head(lm_res$gene, 5), collapse=", "))
      message("  - lm_res$gene_name (first 5):", paste(head(lm_res$gene_name, 5), collapse=", "))
    } else {
      message("  - lm_res$gene (first 5):", paste(head(lm_res$gene, 5), collapse=", "))
    }
    message("  - divergence gene_name (first 5):", paste(head(rd$gene_name, 5), collapse=", "))
  }

  if (verbose) {
    message("\n[effect_sizes_divergence] MERGE STARTING")
    message("  - significant_genes count:", length(significant_genes))
    message("  - lm_res rows:", nrow(lm_res))
    message("  - divergence rowData rows:", nrow(rd))
    message("  - use_gene_name_col:", use_gene_name_col)
  }

  # OPTIMIZATION: Pre-allocate list to collect results instead of using rbind in loop
  # This avoids O(n²) behavior of repeated rbind operations
  result_list <- vector("list", length(significant_genes))
  debug_messages <- character(0)

  for (i in seq_along(significant_genes)) {
    gene_id <- significant_genes[i]
    
    # Extract LMM data for this gene
    lmm_data <- .extractLMMData(lm_res, gene_id, use_gene_name_col)
    if (is.null(lmm_data)) {
      validation_stats$other_errors <- validation_stats$other_errors + 1
      next
    }

    # Collect debug info for first few genes (batch print after loop)
    if (verbose && i <= min(3, length(significant_genes))) {
      debug_messages <- c(debug_messages, 
        sprintf("  [Gene %d] gene_id='%s' match_name='%s'", i, gene_id, lmm_data$match_name))
    }

    # Extract divergence data for this gene
    div_idx <- which(rd$gene_name == lmm_data$match_name)
    if (length(div_idx) == 0) {
      div_idx <- which(rownames(rd) == lmm_data$match_name)
    }
    
    if (length(div_idx) == 0) {
      validation_stats$failed_missing_divergence <- validation_stats$failed_missing_divergence + 1
      next
    }
    
    div_data <- as.data.frame(rd[div_idx[1], , drop = FALSE])

    if (verbose && i <= min(3, length(significant_genes))) {
      debug_messages[length(debug_messages)] <- paste0(debug_messages[length(debug_messages)], 
        " -> found 1 row(s)")
    }

    # Format result row
    if (use_generic) {
      # Single q result
      if (is.na(div_data$estimate[1])) {
        validation_stats$failed_missing_divergence <- validation_stats$failed_missing_divergence + 1
        if (verbose) {
          debug_messages <- c(debug_messages, 
            sprintf("  [Skipped] %s - divergence estimate is NA", lmm_data$match_name))
        }
        next
      }
      
      new_row <- .formatSingleQResult(
        match_name = lmm_data$match_name,
        p_interaction = lmm_data$p_interaction,
        slope_diff = lmm_data$slope_diff,
        div_data = div_data
      )
    } else {
      # Per-q results
      new_row <- .formatMultiQResult(
        match_name = lmm_data$match_name,
        p_interaction = lmm_data$p_interaction,
        slope_diff = lmm_data$slope_diff,
        div_data = div_data,
        q_values = q_values
      )
      
      if (is.null(new_row)) {
        validation_stats$failed_missing_divergence <- validation_stats$failed_missing_divergence + 1
        if (verbose) {
          debug_messages <- c(debug_messages, 
            sprintf("  [Skipped] %s - all divergence estimates are NA", lmm_data$match_name))
        }
        next
      }
    }

    # Store in list (OPTIMIZATION: avoid rbind in loop)
    result_list[[i]] <- new_row
    validation_stats$passed_lmm <- validation_stats$passed_lmm + 1

    if (verbose && i <= min(3, length(significant_genes))) {
      debug_messages <- c(debug_messages, 
        sprintf("  [SUCCESS] %s - p=%.3e, D_spectrum=[...]", 
          lmm_data$match_name, lmm_data$p_interaction))
    }
  }

  # OPTIMIZATION: Replace rbind loop with do.call(rbind) - converts O(n²) to O(n)
  result_list <- result_list[!vapply(result_list, is.null, logical(1))]
  
  if (length(result_list) > 0) {
    interaction_results <- do.call(rbind, result_list)
    rownames(interaction_results) <- NULL  # Reset rownames
  } else {
    interaction_results <- .createResultsDataFrame(q_values, use_generic)
  }

  # Print batched debug messages after loop completes (OPTIMIZATION: move verbose logging outside loop)
  if (verbose && length(debug_messages) > 0) {
    for (msg in debug_messages) {
      message(msg)
    }
  }

  # =========================================================================
  # SUMMARY
  # =========================================================================

  # Print merge results and statistics
  .printMergeSummary(validation_stats, interaction_results, q_values, use_generic, verbose)

  # =========================================================================
  # ENRICH RESULTS: Add per_q_pattern column
  # =========================================================================

  if (enrich_per_q_pattern && nrow(interaction_results) > 0) {
    interaction_results <- .enrichWithQPatterns(
      interaction_results,
      divergence_results_se,
      verbose
    )
  }

  return(list(
    interaction_results = interaction_results,
    validation_stats = validation_stats
  ))
}


# ============================================================================
# INTERNAL HELPER FUNCTIONS
# ============================================================================


#' @noRd
.validateEffectSizeInputs <- function(lm_res, divergence_results_se) {
  if (!is.data.frame(lm_res)) {
    stop("lm_res must be a data frame")
  }

  if (!("gene" %in% colnames(lm_res)) || !("adj_p_interaction" %in% colnames(lm_res))) {
    stop("lm_res must have columns 'gene' and 'adj_p_interaction'")
  }

  if (!methods::is(divergence_results_se, "SummarizedExperiment")) {
    stop("divergence_results_se must be a SummarizedExperiment from .calculate_divergence()")
  }

  rd <- SummarizedExperiment::rowData(divergence_results_se)
  if (!("gene_name" %in% colnames(rd))) {
    stop("divergence_results_se rowData must have 'gene_name' column")
  }
}



#' @noRd
.alignGeneDatasets <- function(lm_res, rd, verbose) {
  # Filter lm_res to include only genes present in divergence_results_se
  
  # Extract gene identifiers from divergence_results_se
  if ("gene_name" %in% colnames(rd)) {
    divergence_genes <- as.character(rd$gene_name)
  } else {
    divergence_genes <- as.character(rownames(rd))
  }
  
  # Get lm_res gene identifiers (prefer gene column, fall back to rownames)
  if ("gene" %in% colnames(lm_res)) {
    lm_res_genes <- as.character(lm_res$gene)
  } else {
    lm_res_genes <- as.character(rownames(lm_res))
    if (length(lm_res_genes) == 0 || all(is.na(lm_res_genes))) {
      stop("lm_res must have either a 'gene' column or valid gene names in rownames")
    }
  }
  
  # Find matching genes and filter lm_res
  matching_idx <- lm_res_genes %in% divergence_genes
  n_before_filter <- nrow(lm_res)
  n_after_filter <- sum(matching_idx)
  
  # Always filter lm_res, even if no matches (results in empty dataframe)
  lm_res <- lm_res[matching_idx, , drop = FALSE]
  
  if (n_after_filter > 0) {
    if (verbose) {
      message("[effect_sizes_divergence] Gene alignment:")
      message("  - lm_res before filtering: ", n_before_filter, " genes")
      message("  - lm_res after filtering: ", n_after_filter, " genes")
      message("  - Genes filtered out: ", n_before_filter - n_after_filter)
    }
  } else {
    if (verbose) {
      message("[effect_sizes_divergence] No matching genes found between lm_res and divergence_results_se")
    }
  }
  
  # Auto-detect available q values from per-q columns
  estimate_cols <- grep("^estimate_q", colnames(rd), value = TRUE)
  
  if (length(estimate_cols) == 0) {
    # Check for generic columns (single q result)
    if (!all(c("estimate", "lower_ci", "upper_ci") %in% colnames(rd))) {
      stop("divergence_results_se rowData must have either:\n",
           "  - Generic columns: 'estimate', 'lower_ci', 'upper_ci', OR\n",
           "  - Per-q columns: 'estimate_q*', 'lower_ci_q*', 'upper_ci_q*'")
    }
    q_values <- NA_real_
    use_generic <- TRUE
    if (verbose) {
      message("[effect_sizes_divergence] Using generic divergence columns (single q-value results)")
    }
  } else {
    # Extract q values from column names
    q_values <- as.numeric(sub("estimate_q", "", estimate_cols))
    q_values <- sort(q_values)  # Sort for consistent output
    use_generic <- FALSE
    
    if (verbose) {
      message("[effect_sizes_divergence] Detected per-q columns for q values: ",
          paste(q_values, collapse = ", "))
    }
  }
  
  return(list(
    lm_res = lm_res,
    divergence_genes = divergence_genes,
    q_values = q_values,
    use_generic = use_generic
  ))
}



#' @noRd
.createResultsDataFrame <- function(q_values, use_generic) {
  interaction_results <- data.frame(
    gene = character(0),
    p_value_interaction = numeric(0),
    slope_diff = numeric(0),
    stringsAsFactors = FALSE
  )
  
  # Add per-q effect size columns
  if (use_generic) {
    interaction_results$effect_size_D <- numeric(0)
    interaction_results$D_lower_ci <- numeric(0)
    interaction_results$D_upper_ci <- numeric(0)
  } else {
    for (q_val in q_values) {
      # Format q value consistently (e.g., "0.5" -> "0_5", "1.0" -> "1_0")
      q_label <- gsub("\\.", "_", sprintf("%.1f", q_val))
      interaction_results[[paste0("effect_size_D_q", q_label)]] <- numeric(0)
      interaction_results[[paste0("D_q", q_label, "_lower_ci")]] <- numeric(0)
      interaction_results[[paste0("D_q", q_label, "_upper_ci")]] <- numeric(0)
    }
  }
  
  return(interaction_results)
}



#' @noRd
.extractLMMData <- function(lm_res, gene_id, use_gene_name_col) {
  lmm_row <- lm_res[lm_res$gene == gene_id, ]
  if (nrow(lmm_row) == 0) {
    return(NULL)
  }

  # Determine the gene name to use for matching against divergence_results_se
  match_name <- if (use_gene_name_col && !is.na(lmm_row$gene_name[1])) {
    lmm_row$gene_name[1]
  } else {
    gene_id
  }

  p_interaction <- lmm_row$adj_p_interaction[1]
  slope_diff <- if ("slope_diff" %in% colnames(lmm_row)) {
    lmm_row$slope_diff[1]
  } else {
    NA_real_
  }
  
  return(list(
    match_name = match_name,
    p_interaction = p_interaction,
    slope_diff = slope_diff
  ))
}



#' @noRd
.extractDivergenceData <- function(rd, match_name, verbose, i, total) {
  # Get divergence info from rowData
  # Try to match by gene_name first (most reliable), then by rownames
  div_row <- rd[rd$gene_name == match_name, ]
  
  if (nrow(div_row) == 0) {
    # If match_name is an ID and we have rownames, try matching rownames
    div_row <- rd[rownames(rd) == match_name, ]
  }

  if (nrow(div_row) == 0) {
    if (verbose && i > min(3, total)) {
      # Only show skipped messages for genes after the debug ones
      message("  [Skipped] ", match_name, " - divergence data not found")
    }
    return(NULL)
  }

  return(as.data.frame(div_row[1, , drop = FALSE]))
}



#' @noRd
.formatSingleQResult <- function(match_name, p_interaction, slope_diff, div_data) {
  data.frame(
    gene = match_name,
    p_value_interaction = p_interaction,
    slope_diff = slope_diff,
    effect_size_D = abs(div_data$estimate[1]),
    D_lower_ci = div_data$lower_ci[1],
    D_upper_ci = div_data$upper_ci[1],
    stringsAsFactors = FALSE
  )
}



#' @noRd
.formatMultiQResult <- function(match_name, p_interaction, slope_diff, div_data, q_values) {
  any_valid <- FALSE
  new_row <- data.frame(
    gene = match_name,
    p_value_interaction = p_interaction,
    slope_diff = slope_diff,
    stringsAsFactors = FALSE
  )
  
  for (q_val in q_values) {
    # Column names in rowData use direct numeric representation
    # (e.g., q=0.5 -> "estimate_q0.5", q=1.0 -> "estimate_q1", q=1.5 -> "estimate_q1.5")
    q_str <- as.character(q_val)
    estimate_col <- paste0("estimate_q", q_str)
    lower_col <- paste0("lower_ci_q", q_str)
    upper_col <- paste0("upper_ci_q", q_str)
    
    # For output columns, convert decimal to underscore (e.g., "1.5" -> "1_5")
    q_label <- gsub("\\.", "_", q_str)
    
    # Safely check if column exists and extract value
    if (estimate_col %in% colnames(div_data)) {
      estimate_val <- div_data[[estimate_col]][1]
    } else {
      estimate_val <- NA_real_
    }
    
    if (lower_col %in% colnames(div_data)) {
      lower_val <- div_data[[lower_col]][1]
    } else {
      lower_val <- NA_real_
    }
    
    if (upper_col %in% colnames(div_data)) {
      upper_val <- div_data[[upper_col]][1]
    } else {
      upper_val <- NA_real_
    }
    
    if (!is.na(estimate_val) && is.finite(estimate_val)) {
      any_valid <- TRUE
      new_row[[paste0("effect_size_D_q", q_label)]] <- abs(estimate_val)
      new_row[[paste0("D_q", q_label, "_lower_ci")]] <- lower_val
      new_row[[paste0("D_q", q_label, "_upper_ci")]] <- upper_val
    } else {
      # Set to NA for this q (handles NA, NaN, Inf cases)
      new_row[[paste0("effect_size_D_q", q_label)]] <- NA_real_
      new_row[[paste0("D_q", q_label, "_lower_ci")]] <- NA_real_
      new_row[[paste0("D_q", q_label, "_upper_ci")]] <- NA_real_
    }
  }
  
  if (!any_valid) {
    return(NULL)
  }
  
  return(new_row)
}



#' @noRd
.printMergeSuccess <- function(lmm_data, div_data, q_values, use_generic) {
  if (use_generic) {
    div_val <- abs(div_data$estimate[1])
    ci_text <- if (!is.na(div_data$lower_ci[1])) {
      sprintf(" CI=[%.4f, %.4f]", div_data$lower_ci[1], div_data$upper_ci[1])
    } else {
      ""
    }
    div_summary <- format(div_val, scientific = TRUE, digits = 3)
  } else {
    # Show summary across q values
    div_vals <- vapply(q_values, function(q) {
      estimate_col <- paste0("estimate_q", q)
      div_data[[estimate_col]][1]
    }, FUN.VALUE = numeric(1))
    div_summary <- paste(sprintf("%.3f", abs(div_vals)), collapse = ", ")
    ci_text <- ""
  }
  
  message("  [SUCCESS] ", lmm_data$match_name, " - p=", 
      format(lmm_data$p_interaction, digits = 3), ", D_spectrum=[", 
      div_summary, "]", ci_text)
}



#' @noRd
.printMergeSummary <- function(validation_stats, interaction_results, q_values, use_generic, verbose) {
  if (!verbose) {
    return(invisible(NULL))
  }
  
  message("\n[effect_sizes_divergence] MERGE COMPLETED\n",
          "  - Total significant genes: ", validation_stats$total_genes, "\n",
          "  - Passed merge: ", validation_stats$passed_lmm, "\n",
          "  - Failed (missing divergence): ", validation_stats$failed_missing_divergence, "\n",
          "  - Other errors: ", validation_stats$other_errors, "\n",
          "  - interaction_results rows: ", nrow(interaction_results))

  if (nrow(interaction_results) > 0) {
    message("\n**Effect Size Distribution Across q Values:**")
    
    if (use_generic) {
      div_col <- "effect_size_D"
      message("- Mean D:", round(mean(interaction_results[[div_col]], na.rm = TRUE), 4))
      message("- Median D:", round(median(interaction_results[[div_col]], na.rm = TRUE), 4))
      message("- Range: [", 
          round(min(interaction_results[[div_col]], na.rm = TRUE), 4), ", ",
          round(max(interaction_results[[div_col]], na.rm = TRUE), 4), "]")
    } else {
      for (q_val in q_values) {
        q_label <- gsub("\\.", "_", as.character(q_val))
        div_col <- paste0("effect_size_D_q", q_label)
        if (div_col %in% colnames(interaction_results)) {
          valid_vals <- interaction_results[[div_col]][!is.na(interaction_results[[div_col]])]
          if (length(valid_vals) > 0) {
            message("- q=", q_val, ": mean=", round(mean(valid_vals, na.rm = TRUE), 4),
                ", median=", round(median(valid_vals, na.rm = TRUE), 4))
          }
        }
      }
    }
    
    message("- Interpretation: D > 0.05 = small, D > 0.1 = medium, D > 0.2 = large")
  }
  
  invisible(NULL)
}



#' @noRd
.enrichWithQPatterns <- function(interaction_results, divergence_results_se, verbose) {
  # Extract gene names from divergence_results_se rowData and assay matrix
  div_assay <- SummarizedExperiment::assay(divergence_results_se)
  div_rd <- as.data.frame(SummarizedExperiment::rowData(divergence_results_se))

  if (nrow(div_assay) > 0 && ncol(div_assay) > 0) {
    # Map genes from divergence_results_se
    div_gene_names <- if ("gene_name" %in% colnames(div_rd)) {
      div_rd$gene_name
    } else {
      rownames(div_assay)
    }

    # Create per_q_pattern column: classify divergence patterns
    # RARE_DRIVEN = divergence higher at low q (rare isoforms drive changes)
    # ABUNDANT_DRIVEN = divergence higher at high q (abundant isoforms drive changes)  
    # BALANCED = similar divergence across diversity scales
    per_q_patterns <- character(nrow(interaction_results))
    for (i in seq_len(nrow(interaction_results))) {
      gene_name <- interaction_results$gene[i]
      gene_idx <- which(div_gene_names == gene_name)

      if (length(gene_idx) > 0) {
        # Get divergence values for this gene across q values
        divs <- div_assay[gene_idx[1], ]
        
        # Create named vector for classify_q_pattern
        # Column names in divs should be like "q_0.01", "q_0.5", "q_1.0", etc.
        per_q_patterns[i] <- .classify_q_pattern(divs)
        
        # If classification failed, return "UNCLASSIFIED"
        if (is.na(per_q_patterns[i])) {
          per_q_patterns[i] <- "UNCLASSIFIED"
        }
      }
    }
    interaction_results$per_q_pattern <- per_q_patterns
  }
  
  return(interaction_results)
}



#' @noRd
.classify_q_pattern <- function(per_q_divs, ratio_threshold = 1.3) {
  # Classify q-value divergence pattern based on median divergence in rare vs abundant regions
  # per_q_divs: named vector where names are like "q_0.01", "q_0.5", "q_1.0", "q_2.0"
  #             or unnamed numeric vector (classification by position)
  # ratio_threshold: factor for classifying patterns
  #   - ratio > threshold → RARE_DRIVEN (rare isoforms drive effect)
  #   - ratio < 1/threshold → ABUNDANT_DRIVEN (abundant isoforms drive effect)
  #   - 1/threshold <= ratio <= threshold → BALANCED (all regions equally important)
  
  # Input validation
  if (length(per_q_divs) == 0 || all(is.na(per_q_divs))) {
    return(NA_character_)
  }
  
  # Check if input can be converted to numeric
  # as.numeric() naturally produces NAs for non-numeric values—no wrapping needed
  divs_numeric <- as.numeric(per_q_divs)
  
  if (all(is.na(divs_numeric))) {
    return(NA_character_)
  }
  
  # Check if input has names
  has_names <- !is.null(names(per_q_divs)) && length(names(per_q_divs)) > 0
  
  # If it has names, they must all be valid q_ names or all be NA (which is invalid)
  if (has_names) {
    # Can't have any NA names
    if (any(is.na(names(per_q_divs)))) {
      return(NA_character_)
    }
    
    # Check if names start with "q_"
    all_q_names <- all(grepl("^q_", names(per_q_divs)))
    
    if (!all_q_names) {
      # Names exist but don't match q_ pattern - invalid
      return(NA_character_)
    }
  }
  
  # Detect if input has named q-values: check if names start with "q_"
  has_q_names <- has_names && all(grepl("^q_", names(per_q_divs)))
  
  # Case 1: Named vector with q-value names (e.g., "q_0.01", "q_0.5", "q_1.0")
  if (has_q_names) {
    # For named q-values, require at least 2 values for classification
    if (length(per_q_divs) < 2) {
      return(NA_character_)
    }
    
    # Extract q values from names and classify by ratio
    # Try to extract numeric part after "q_" prefix  
    # Handle both "q_0.01" and "q_0_01" formats
    name_parts <- sub("q_", "", names(per_q_divs))
    # Replace underscore with dot for parsing
    name_parts_normalized <- gsub("_", ".", name_parts)
    # as.numeric() naturally produces NAs for non-numeric values
    q_values <- as.numeric(name_parts_normalized)
    
    # Check if all q values parsed successfully (all non-NA)
    if (all(is.na(q_values))) {
      # If q name extraction failed completely, return NA
      return(NA_character_)
    }
    
    # Separate into rare (q < 1) and abundant (q > 1) regions
    # q < 1: emphasizes rare (low-probability) isoforms
    # q > 1: emphasizes abundant (high-probability) isoforms
    # q = 1: neutral (Shannon entropy) - excluded from both regions
    rare_divs <- divs_numeric[!is.na(q_values) & q_values < 1.0]
    abundant_divs <- divs_numeric[!is.na(q_values) & q_values > 1.0]
    
    # Remove NA values
    valid_rare <- rare_divs[!is.na(rare_divs)]
    valid_abundant <- abundant_divs[!is.na(abundant_divs)]
    
    # Need at least one valid value in EACH region for classification
    if (length(valid_rare) == 0 || length(valid_abundant) == 0) {
      # Insufficient data in regions for classification
      return(NA_character_)
    }
    
    # Calculate median divergence in each region
    rare_median <- median(valid_rare, na.rm = TRUE)
    abundant_median <- median(valid_abundant, na.rm = TRUE)
    
    # Avoid division by zero or NA
    if (!is.na(rare_median) && !is.na(abundant_median) && abundant_median != 0) {
      # Calculate ratio
      ratio <- rare_median / abundant_median
      
      # Classify based on ratio (with NA check)
      if (!is.na(ratio) && ratio > ratio_threshold) {
        return("RARE_DRIVEN")
      } else if (!is.na(ratio) && ratio < 1 / ratio_threshold) {
        return("ABUNDANT_DRIVEN")
      } else {
        return("BALANCED")
      }
    }
    
    # If median calculation failed, return NA
    return(NA_character_)
  }
  
  # Case 2: Unnamed vector
  # Unnamed vectors lack semantic meaning (no q-value labels) → return NA
  # This requires proper q-value names for meaningful classification
  return(NA_character_)
}


