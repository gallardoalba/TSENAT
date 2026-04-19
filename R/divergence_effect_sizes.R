#' Merge LMM Interaction Results with Tsallis Divergence Effect Sizes
#'
#' Combines LMM interaction test p-values with pre-computed Tsallis divergence
#' effect sizes and bootstrap confidence intervals across ALL q values. This
#' is a
#' **data merger**, not a model fitter--all statistical computation happens
#' upstream in:
#' - `.calculate_sait()` -> SAIT/LMM p-values
#' - `.calculate_divergence()` -> Divergence estimates and CIs for multiple q
#'
#' This function merges the results into a single data frame for downstream
#' interpretation. When multiple q values are present, effect sizes are computed
#' for each q to capture the full biological spectrum (rare->abundant isoforms).
#'
#' **Architecture:**
#' ```
#' Input 1: LMM results (from calculate_sait_interaction)
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
#' @param sait_res A data frame of LMM interaction test results from
#' `.calculate_sait()`,
#' with columns: `gene` (character, gene name), `adj_p_interaction`
#' (numeric, multiple-test adjusted p-value).
#'   Genes with adj_p_interaction below `significance_threshold` are included.
#'
#' @param divergence_results_se A SummarizedExperiment from
#' `.calculate_divergence()`,
#'   containing rowData with columns: `gene_name` and either:
#'   - Generic: `estimate`, `lower_ci`, `upper_ci` (single q-value results), OR
#'   - Per-q: `estimate_q*`, `lower_ci_q*`, `upper_ci_q*` (multiple q-values)
#'   All divergence-related data is self-contained in this object.
#'
#' @param significance_threshold Numeric; p-value threshold for filtering
#' significant
#' genes (default: 0.05). Only genes with adj_p_interaction < threshold are
#' included.
#'
#' @param enrich_per_q_pattern Logical; if TRUE (default), adds a
#' 'per_q_pattern' column
#' to the output data frame containing comma-separated divergence values
#' across the
#' q spectrum for each gene. This column enables visualization and
#' classification
#' of whether treatment effects are driven by rare (low-q) or abundant
#' (high-q)
#' isoforms. Set to FALSE to reduce output size if this annotation is not
#' needed.
#'
#' @param verbose Logical; if TRUE, print detailed validation and merge
#' statistics
#' to console (default: TRUE). Shows counts of passed, skipped, and failed
#' genes.
#'
#' @return A list with two elements:
#'   \describe{
#'     \item{\code{interaction_results}}{Data frame (genes * columns) with 
#' merged results:
#'       - `gene`: Gene name (character)
#'       - `p_value_interaction`: LMM adjusted p-value for q:group interaction
#' - `slope_diff`: q:group interaction slope coefficient (if present in
#' sait_res)
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
#' are driven by rare transcripts (high D at low q) or abundant transcripts
#' (high D at high q).
#'
#' For paired designs, divergence is computed separately within each pair,
#' then averaged to account for pairing structure.
#'
#' **Database References:**
#' - Papers I002-I004: Tsallis divergence mathematical foundation
#' - Papers Ramsay (2005), Springer Series in Statistics: Bootstrap CI computation respecting data structure
#' - Papers S197: Quality filtering and effect size thresholds
#'
#' @noRd
.calculate_effect_sizes <- function(sait_res, divergence_results_se, significance_threshold = 0.05,
    enrich_per_q_pattern = TRUE, verbose = FALSE) {
    .validateEffectSizeInputs(sait_res, divergence_results_se)
    rd <- SummarizedExperiment::rowData(divergence_results_se)
    align_list <- .alignGeneDatasets(sait_res, rd, verbose)
    filter_list <- .filterSignificantGenes(align_list$sait_res, significance_threshold,
        align_list$q_values, align_list$use_generic, verbose)
    if (length(filter_list$significant_genes) == 0) {
        if (verbose)
            message("No genes with significant q*group interaction detected.")
        return(list(interaction_results = filter_list$empty_results, validation_stats = list(total_genes = 0,
            passed_lmm = 0, failed_missing_divergence = 0, other_errors = 0, q_values = align_list$q_values)))
    }
    merge_list <- .mergeEffectSizesForGenes(align_list$sait_res, rd, filter_list$significant_genes,
        align_list$q_values, align_list$use_generic, verbose)
    .printMergeSummary(merge_list$validation_stats, merge_list$interaction_results,
        align_list$q_values, align_list$use_generic, verbose)
    if (enrich_per_q_pattern && nrow(merge_list$interaction_results) > 0) {
        merge_list$interaction_results <- .enrichWithQPatterns(merge_list$interaction_results,
            divergence_results_se, verbose)
    }
    return(list(interaction_results = merge_list$interaction_results, validation_stats = merge_list$validation_stats))
}


# ============================================================================
# REFACTORED HELPER FUNCTIONS
# ============================================================================

#' @noRd
.filterSignificantGenes <- function(sait_res, significance_threshold, q_values, use_generic,
    verbose) {
    # Filter to genes with valid p-values and adj_p_interaction < threshold
    valid_p_idx <- !is.na(sait_res$adj_p_interaction)
    significant_idx <- valid_p_idx & (sait_res$adj_p_interaction < significance_threshold)
    significant_genes <- sait_res$gene[significant_idx]

    if (verbose) {
        message("\n**Filtering effect size analysis to significant genes:**")
        message("- Genes with valid p-values: ", sum(valid_p_idx))
        message("- Genes with adj_p_interaction <", significance_threshold, ": ",
            length(significant_genes))
    }

    # Create empty results data frame
    empty_results <- .createResultsDataFrame(q_values, use_generic)

    return(list(significant_genes = significant_genes, empty_results = empty_results))
}


#' @noRd
.mergeEffectSizesForGenes <- function(sait_res, rd, significant_genes, q_values, use_generic,
    verbose) {
    validation_stats <- list(total_genes = length(significant_genes), passed_lmm = 0,
        failed_missing_divergence = 0, other_errors = 0, q_values = q_values)

    if (verbose) {
        message("\nMerging LMM results with divergence effect sizes...")
    }

    # Get matching strategy info
    use_gene_name_col <- "gene_name" %in% colnames(sait_res)

    # OPTIMIZATION: Consolidate verbose logging into a single condition block
    if (verbose) {
        message("[calculate_effect_sizes] Gene name matching strategy:")
        message("  - gene_name column in sait_res:", use_gene_name_col)
        if (use_gene_name_col) {
            message("  - sait_res$gene (first 5):", paste(head(sait_res$gene, 5), collapse = ", "))
            message("  - sait_res$gene_name (first 5):", paste(head(sait_res$gene_name,
                5), collapse = ", "))
        } else {
            message("  - sait_res$gene (first 5):", paste(head(sait_res$gene, 5), collapse = ", "))
        }
        message("  - divergence gene_name (first 5):", paste(head(rd$gene_name, 5),
            collapse = ", "))
        message("\n[calculate_effect_sizes] MERGE STARTING")
        message("  - significant_genes count:", length(significant_genes))
        message("  - sait_res rows:", nrow(sait_res))
        message("  - divergence rowData rows:", nrow(rd))
        message("  - use_gene_name_col:", use_gene_name_col)
    }

    # OPTIMIZATION: Pre-compute gene name index map ONCE instead of searching
    # per gene
    gene_idx_map <- match(significant_genes, rd$gene_name)
    missing_idx <- which(is.na(gene_idx_map))
    if (length(missing_idx) > 0) {
        rowname_idx <- match(significant_genes[missing_idx], rownames(rd))
        gene_idx_map[missing_idx] <- rowname_idx
    }

    if (verbose) {
        n_found <- sum(!is.na(gene_idx_map))
        message("[calculate_effect_sizes] OPTIMIZATION: Pre-computed gene index map")
        message("  - Genes found: ", n_found, "/", length(significant_genes))
    }

    # OPTIMIZATION: Pre-allocate and pre-compute all q-metadata ONCE for all
    # genes
    result_list <- vector("list", length(significant_genes))
    debug_messages <- character(0)
    verbose_threshold <- if (verbose)
        min(3, length(significant_genes)) else 0

    # Pre-compute q-labels that will be used for all genes (not per-gene)
    q_strs_precomp <- NULL
    q_labels_precomp <- NULL
    estimate_cols_precomp <- NULL
    lower_cols_precomp <- NULL
    upper_cols_precomp <- NULL

    if (!use_generic) {
        q_strs_precomp <- as.character(q_values)
        q_labels_precomp <- gsub("\\.", "_", q_strs_precomp)
        estimate_cols_precomp <- paste0("estimate_q", q_strs_precomp)
        lower_cols_precomp <- paste0("lower_ci_q", q_strs_precomp)
        upper_cols_precomp <- paste0("upper_ci_q", q_strs_precomp)
    }

    for (i in seq_along(significant_genes)) {
        gene_id <- significant_genes[i]

        # Extract LMM data for this gene
        lmm_data <- .extractLMMData(sait_res, gene_id, use_gene_name_col)
        if (is.null(lmm_data)) {
            validation_stats$other_errors <- validation_stats$other_errors + 1
            next
        }

        # Collect debug info for first few genes (batch print after loop)
        if (verbose && i <= verbose_threshold) {
            debug_messages <- c(debug_messages, sprintf("  [Gene %d] gene_id='%s' match_name='%s'",
                i, gene_id, lmm_data$match_name))
        }

        # OPTIMIZATION: Use pre-computed index instead of searching
        div_idx <- gene_idx_map[i]

        if (is.na(div_idx)) {
            validation_stats$failed_missing_divergence <- validation_stats$failed_missing_divergence +
                1
            next
        }

        div_data <- as.data.frame(rd[div_idx, , drop = FALSE])

        if (verbose && i <= verbose_threshold) {
            debug_messages[length(debug_messages)] <- paste0(debug_messages[length(debug_messages)],
                " -> found 1 row(s)")
        }

        # Format result row
        if (use_generic) {
            # Single q result
            if (is.na(div_data$estimate[1])) {
                validation_stats$failed_missing_divergence <- validation_stats$failed_missing_divergence +
                  1
                if (verbose) {
                  debug_messages <- c(debug_messages, sprintf("  [Skipped] %s - divergence estimate is NA",
                    lmm_data$match_name))
                }
                next
            }

            new_row <- .formatSingleQResult(match_name = lmm_data$match_name, p_interaction = lmm_data$p_interaction,
                slope_diff = lmm_data$slope_diff, div_data = div_data)
        } else {
            # OPTIMIZATION: Reuse pre-computed columns if multi-q
            if (!use_generic) {
                col_cache <- .buildColumnCache(div_data, q_values)
                new_row <- .formatMultiQResult(match_name = lmm_data$match_name,
                  p_interaction = lmm_data$p_interaction, slope_diff = lmm_data$slope_diff,
                  div_data = div_data, q_values = q_values, col_cache = col_cache,
                  q_strs = q_strs_precomp, q_labels = q_labels_precomp, estimate_cols = estimate_cols_precomp,
                  lower_cols = lower_cols_precomp, upper_cols = upper_cols_precomp)
            } else {
                new_row <- .formatSingleQResult(match_name = lmm_data$match_name,
                  p_interaction = lmm_data$p_interaction, slope_diff = lmm_data$slope_diff,
                  div_data = div_data)
            }

            if (is.null(new_row)) {
                validation_stats$failed_missing_divergence <- validation_stats$failed_missing_divergence +
                  1
                if (verbose) {
                  debug_messages <- c(debug_messages, sprintf("  [Skipped] %s - all divergence estimates are NA",
                    lmm_data$match_name))
                }
                next
            }
        }

        # Store in list (OPTIMIZATION: avoid rbind in loop)
        result_list[[i]] <- new_row
        validation_stats$passed_lmm <- validation_stats$passed_lmm + 1

        if (verbose && i <= verbose_threshold) {
            debug_messages <- c(debug_messages, sprintf("  [SUCCESS] %s - p=%.3e, D_spectrum=[...]",
                lmm_data$match_name, lmm_data$p_interaction))
        }
    }

    # OPTIMIZATION: Replace rbind loop with do.call(rbind) - converts O(n²) to
    # O(n)
    result_list <- result_list[!vapply(result_list, is.null, logical(1))]

    if (length(result_list) > 0) {
        interaction_results <- do.call(rbind, result_list)
        rownames(interaction_results) <- NULL  # Reset rownames
    } else {
        interaction_results <- .createResultsDataFrame(q_values, use_generic)
    }

    # Print batched debug messages after loop completes
    if (verbose && length(debug_messages) > 0) {
        for (msg in debug_messages) {
            message(msg)
        }
    }

    return(list(interaction_results = interaction_results, validation_stats = validation_stats))
}


# ============================================================================
# INTERNAL HELPER FUNCTIONS
# ============================================================================

# OPTIMIZATION: Pre-compute column existence for all q-values This runs ONCE
# per gene instead of per q-value, reducing column checks from gene_count ×
# q_count to gene_count × 1
#' @noRd
.buildColumnCache <- function(div_data, q_values) {
    # Build a matrix showing which columns exist for which q-values
    # OPTIMIZATION: Vectorized operations - pre-compute all column names once
    col_cache <- matrix(FALSE, nrow = length(q_values), ncol = 3)
    colnames(col_cache) <- c("estimate", "lower", "upper")

    q_strs <- as.character(q_values)
    col_names <- colnames(div_data)

    # Single pass with vectorized operations
    estimate_names <- paste0("estimate_q", q_strs)
    lower_names <- paste0("lower_ci_q", q_strs)
    upper_names <- paste0("upper_ci_q", q_strs)

    col_cache[, "estimate"] <- !is.na(match(estimate_names, col_names))
    col_cache[, "lower"] <- !is.na(match(lower_names, col_names))
    col_cache[, "upper"] <- !is.na(match(upper_names, col_names))

    return(col_cache)
}

#' @noRd
.validateEffectSizeInputs <- function(sait_res, divergence_results_se) {
    if (!is.data.frame(sait_res)) {
        stop("sait_res must be a data frame")
    }

    if (!("gene" %in% colnames(sait_res)) || !("adj_p_interaction" %in% colnames(sait_res))) {
        stop("sait_res must have columns 'gene' and 'adj_p_interaction'")
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
.alignGeneDatasets <- function(sait_res, rd, verbose) {
    # Filter sait_res to include only genes present in divergence_results_se

    # OPTIMIZATION: Consolidate string conversions in single pass
    divergence_genes <- as.character(if ("gene_name" %in% colnames(rd)) rd$gene_name else rownames(rd))
    sait_res_genes <- as.character(if ("gene" %in% colnames(sait_res)) sait_res$gene else rownames(sait_res))

    if (length(sait_res_genes) == 0 || all(is.na(sait_res_genes))) {
        stop("sait_res must have either a 'gene' column or valid gene names in rownames")
    }

    # Find matching genes and filter sait_res
    matching_idx <- sait_res_genes %in% divergence_genes
    n_before_filter <- nrow(sait_res)
    n_after_filter <- sum(matching_idx)

    # Always filter sait_res, even if no matches (results in empty dataframe)
    sait_res <- sait_res[matching_idx, , drop = FALSE]

    if (n_after_filter > 0) {
        if (verbose) {
            message("[calculate_effect_sizes] Gene alignment:")
            message("  - sait_res before filtering: ", n_before_filter, " genes")
            message("  - sait_res after filtering: ", n_after_filter, " genes")
            message("  - Genes filtered out: ", n_before_filter - n_after_filter)
        }
    } else {
        if (verbose) {
            message("[calculate_effect_sizes] No matching genes found between sait_res and divergence_results_se")
        }
    }

    # Auto-detect available q values from per-q columns
    estimate_cols <- grep("^estimate_q", colnames(rd), value = TRUE)

    if (length(estimate_cols) == 0) {
        # Check for generic columns (single q result)
        if (!all(c("estimate", "lower_ci", "upper_ci") %in% colnames(rd))) {
            stop("divergence_results_se rowData must have either:\n", "  - Generic columns: 'estimate', 'lower_ci', 'upper_ci', OR\n",
                "  - Per-q columns: 'estimate_q*', 'lower_ci_q*', 'upper_ci_q*'")
        }
        q_values <- NA_real_
        use_generic <- TRUE
        if (verbose) {
            message("[calculate_effect_sizes] Using generic divergence columns (single q-value results)")
        }
    } else {
        # Extract q values from column names
        q_values <- as.numeric(sub("estimate_q", "", estimate_cols))
        q_values <- sort(q_values)  # Sort for consistent output
        use_generic <- FALSE

        if (verbose) {
            message("[calculate_effect_sizes] Detected per-q columns for q values: ",
                paste(q_values, collapse = ", "))
        }
    }

    return(list(sait_res = sait_res, divergence_genes = divergence_genes, q_values = q_values,
        use_generic = use_generic))
}



#' @noRd
.createResultsDataFrame <- function(q_values, use_generic) {
    # OPTIMIZATION: Build entire df at once with vectorized column creation
    if (use_generic) {
        interaction_results <- data.frame(gene = character(0), p_value_interaction = numeric(0),
            slope_diff = numeric(0), effect_size_D = numeric(0), D_lower_ci = numeric(0),
            D_upper_ci = numeric(0), stringsAsFactors = FALSE)
    } else {
        q_labels <- gsub("\\.", "_", sprintf("%.1f", q_values))
        estimate_cols <- paste0("effect_size_D_q", q_labels)
        lower_cols <- paste0("D_q", q_labels, "_lower_ci")
        upper_cols <- paste0("D_q", q_labels, "_upper_ci")

        # All columns at once using setNames
        col_list <- setNames(replicate(length(c(estimate_cols, lower_cols, upper_cols)),
            numeric(0)), c(estimate_cols, lower_cols, upper_cols))
        interaction_results <- data.frame(gene = character(0), p_value_interaction = numeric(0),
            slope_diff = numeric(0), stringsAsFactors = FALSE)
        interaction_results <- cbind(interaction_results, as.data.frame(col_list))
    }

    return(interaction_results)
}



#' @noRd
.extractLMMData <- function(sait_res, gene_id, use_gene_name_col) {
    lmm_row <- sait_res[sait_res$gene == gene_id, ]
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

    return(list(match_name = match_name, p_interaction = p_interaction, slope_diff = slope_diff))
}



#' @noRd
.extractDivergenceData <- function(rd, match_name, verbose, i, total) {
    # Get divergence info from rowData Try to match by gene_name first (most
    # reliable), then by rownames
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
    data.frame(gene = match_name, p_value_interaction = p_interaction, slope_diff = slope_diff,
        effect_size_D = abs(div_data$estimate[1]), D_lower_ci = div_data$lower_ci[1],
        D_upper_ci = div_data$upper_ci[1], stringsAsFactors = FALSE)
}



#' @noRd
.formatMultiQResult <- function(match_name, p_interaction, slope_diff, div_data,
    q_values, col_cache = NULL, q_strs = NULL, q_labels = NULL, estimate_cols = NULL,
    lower_cols = NULL, upper_cols = NULL) {
    # OPTIMIZATION: Accept pre-computed q-labels parameter for reuse across
    # genes
    if (is.null(col_cache)) {
        col_cache <- .buildColumnCache(div_data, q_values)
    }

    if (is.null(q_strs)) {
        q_strs <- as.character(q_values)
        q_labels <- gsub("\\.", "_", q_strs)
        estimate_cols <- paste0("estimate_q", q_strs)
        lower_cols <- paste0("lower_ci_q", q_strs)
        upper_cols <- paste0("upper_ci_q", q_strs)
    }

    any_valid <- FALSE
    new_row <- data.frame(gene = match_name, p_value_interaction = p_interaction,
        slope_diff = slope_diff, stringsAsFactors = FALSE)

    output_cols_est <- paste0("effect_size_D_q", q_labels)
    output_cols_lower <- paste0("D_q", q_labels, "_lower_ci")
    output_cols_upper <- paste0("D_q", q_labels, "_upper_ci")

    for (q_idx in seq_along(q_values)) {
        # Use pre-computed column names to avoid paste0() in loop
        estimate_val <- if (col_cache[q_idx, "estimate"])
            div_data[[estimate_cols[q_idx]]][1] else NA_real_
        lower_val <- if (col_cache[q_idx, "lower"])
            div_data[[lower_cols[q_idx]]][1] else NA_real_
        upper_val <- if (col_cache[q_idx, "upper"])
            div_data[[upper_cols[q_idx]]][1] else NA_real_

        if (!is.na(estimate_val) && is.finite(estimate_val)) {
            any_valid <- TRUE
            new_row[[output_cols_est[q_idx]]] <- abs(estimate_val)
            new_row[[output_cols_lower[q_idx]]] <- lower_val
            new_row[[output_cols_upper[q_idx]]] <- upper_val
        } else {
            new_row[[output_cols_est[q_idx]]] <- NA_real_
            new_row[[output_cols_lower[q_idx]]] <- NA_real_
            new_row[[output_cols_upper[q_idx]]] <- NA_real_
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

    message("  [SUCCESS] ", lmm_data$match_name, " - p=", format(lmm_data$p_interaction,
        digits = 3), ", D_spectrum=[", div_summary, "]", ci_text)
}



#' @noRd
.printMergeSummary <- function(validation_stats, interaction_results, q_values, use_generic,
    verbose) {
    if (!verbose) {
        return(invisible(NULL))
    }

    message("\n[calculate_effect_sizes] MERGE COMPLETED\n", "  - Total significant genes: ",
        validation_stats$total_genes, "\n", "  - Passed merge: ", validation_stats$passed_lmm,
        "\n", "  - Failed (missing divergence): ", validation_stats$failed_missing_divergence,
        "\n", "  - Other errors: ", validation_stats$other_errors, "\n", "  - interaction_results rows: ",
        nrow(interaction_results))

    if (nrow(interaction_results) > 0) {
        message("\n**Effect Size Distribution Across q Values:**")

        if (use_generic) {
            div_col <- "effect_size_D"
            message("- Mean D:", round(mean(interaction_results[[div_col]], na.rm = TRUE),
                4))
            message("- Median D:", round(median(interaction_results[[div_col]], na.rm = TRUE),
                4))
            message("- Range: [", round(min(interaction_results[[div_col]], na.rm = TRUE),
                4), ", ", round(max(interaction_results[[div_col]], na.rm = TRUE),
                4), "]")
        } else {
            # OPTIMIZATION: Pre-compute q-value labels once instead of in loop
            q_labels <- gsub("\\.", "_", as.character(q_values))
            div_cols <- paste0("effect_size_D_q", q_labels)

            for (i in seq_along(q_values)) {
                if (div_cols[i] %in% colnames(interaction_results)) {
                  valid_vals <- interaction_results[[div_cols[i]]][!is.na(interaction_results[[div_cols[i]]])]
                  if (length(valid_vals) > 0) {
                    message("- q=", q_values[i], ": mean=", round(mean(valid_vals,
                      na.rm = TRUE), 4), ", median=", round(median(valid_vals, na.rm = TRUE),
                      4))
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

        # OPTIMIZATION: Pre-compute gene index map ONCE with vectorized match()
        gene_idx_map <- match(interaction_results$gene, div_gene_names)

        # OPTIMIZATION: Vectorized classification - classify all genes at once
        per_q_patterns <- character(nrow(interaction_results))
        rare_median_vals <- numeric(nrow(interaction_results))
        abundant_median_vals <- numeric(nrow(interaction_results))
        q_ratio_vals <- numeric(nrow(interaction_results))

        for (i in seq_len(nrow(interaction_results))) {
            gene_idx <- gene_idx_map[i]

            if (!is.na(gene_idx)) {
                # Get divergence values for this gene across q values
                divs <- div_assay[gene_idx, ]
                # .classify_q_pattern now returns list with pattern + metrics
                class_result <- .classify_q_pattern(divs)
                per_q_patterns[i] <- class_result$pattern
                rare_median_vals[i] <- class_result$rare_median
                abundant_median_vals[i] <- class_result$abundant_median
                q_ratio_vals[i] <- class_result$ratio

                # If classification failed, mark as UNCLASSIFIED
                if (is.na(per_q_patterns[i])) {
                  per_q_patterns[i] <- "UNCLASSIFIED"
                }
            } else {
                per_q_patterns[i] <- "UNCLASSIFIED"
                rare_median_vals[i] <- NA_real_
                abundant_median_vals[i] <- NA_real_
                q_ratio_vals[i] <- NA_real_
            }
        }
        interaction_results$per_q_pattern <- per_q_patterns
        interaction_results$div_rare_median <- rare_median_vals
        interaction_results$div_abundant_median <- abundant_median_vals
        interaction_results$q_ratio <- q_ratio_vals
    }

    return(interaction_results)
}


#' @noRd
.classify_q_pattern <- function(per_q_divs, ratio_threshold = 1.3) {
    # Classify q-value divergence pattern based on median divergence Returns
    # list with: pattern (classification), rare_median, abundant_median, ratio
    # OPTIMIZATION: Early exit for invalid inputs
    if (length(per_q_divs) == 0 || all(is.na(per_q_divs))) {
        return(list(pattern = NA_character_, rare_median = NA_real_, abundant_median = NA_real_,
            ratio = NA_real_))
    }

    divs_numeric <- as.numeric(per_q_divs)
    if (all(is.na(divs_numeric))) {
        return(list(pattern = NA_character_, rare_median = NA_real_, abundant_median = NA_real_,
            ratio = NA_real_))
    }

    # OPTIMIZATION: Consolidate name validation in single pass
    has_names <- !is.null(names(per_q_divs)) && length(names(per_q_divs)) > 0
    if (has_names && any(is.na(names(per_q_divs)))) {
        return(list(pattern = NA_character_, rare_median = NA_real_, abundant_median = NA_real_,
            ratio = NA_real_))
    }

    # OPTIMIZATION: Quick first-element check before all(grepl())
    pattern_names <- if (has_names)
        names(per_q_divs) else character(0)
    has_q_names <- has_names && (length(pattern_names) > 0 && substr(pattern_names[1],
        1, 2) == "q_")

    # Case 1: Named vector with q-value names
    if (has_q_names) {
        if (length(per_q_divs) < 2) {
            return(list(pattern = NA_character_, rare_median = NA_real_, abundant_median = NA_real_,
                ratio = NA_real_))
        }

        # OPTIMIZATION: Single pass to extract q-values and classify
        name_parts <- sub("q_", "", pattern_names)
        # Replace underscores with dots for parsing
        name_parts_normalized <- gsub("_", ".", name_parts)
        q_values <- as.numeric(name_parts_normalized)

        if (all(is.na(q_values))) {
            return(list(pattern = NA_character_, rare_median = NA_real_, abundant_median = NA_real_,
                ratio = NA_real_))
        }

        # OPTIMIZATION: Vectorized filtering in single pass
        rare_idx <- !is.na(q_values) & q_values < 1
        abundant_idx <- !is.na(q_values) & q_values > 1

        valid_rare <- divs_numeric[rare_idx & !is.na(divs_numeric)]
        valid_abundant <- divs_numeric[abundant_idx & !is.na(divs_numeric)]

        # Need at least one valid value in EACH region
        if (length(valid_rare) == 0 || length(valid_abundant) == 0) {
            return(list(pattern = NA_character_, rare_median = NA_real_, abundant_median = NA_real_,
                ratio = NA_real_))
        }

        # Calculate medians
        rare_median <- median(valid_rare, na.rm = TRUE)
        abundant_median <- median(valid_abundant, na.rm = TRUE)

        # OPTIMIZATION: Single ratio calculation
        if (!is.na(rare_median) && !is.na(abundant_median) && abundant_median !=
            0) {
            ratio <- rare_median/abundant_median
            if (!is.na(ratio)) {
                pattern <- if (ratio > ratio_threshold) {
                  "Rare driven"
                } else if (ratio < 1/ratio_threshold) {
                  "Abundant driven"
                } else {
                  "Balanced"
                }
                return(list(pattern = pattern, rare_median = rare_median, abundant_median = abundant_median,
                  ratio = ratio))
            }
        }
        return(list(pattern = NA_character_, rare_median = NA_real_, abundant_median = NA_real_,
            ratio = NA_real_))
    }

    # Case 2: Unnamed vector
    return(list(pattern = NA_character_, rare_median = NA_real_, abundant_median = NA_real_,
        ratio = NA_real_))
}

