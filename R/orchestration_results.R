# ============================================================================
# RESULT ACCESSOR FUNCTIONS (EXPORTED)
# ============================================================================

#' Extract analysis results from TSENATAnalysis object
#'
#' Provides flexible access to diversity, divergence, and statistical test results
#' with options for ranking, filtering, and format conversion.
#'
#' @param analysis \code{TSENATAnalysis} object containing computed results.
#' @param type \code{character}. \strong{Required.} Type of results to extract:
#'   'diversity', 'divergence', 'lm', 'jackknife', 'rank_test', 'effect_sizes_divergence',
#'   'assumptions', or 'switching_tables'.
#' @param q \code{numeric}. For diversity results, optionally return results for 
#'   a specific q-value only. When specified, returns a single SummarizedExperiment 
#'   for that q-value instead of the full list. Default: NULL (return all q-values 
#'   as list). Example: q = 1.0 returns only the q=1.0 results.
#' @param rankBy \code{character}. For LM/Jackknife results, ranking method:
#'   'none' (default), 'pvalue', 'effectSize', or 'qvalue'.
#'   Applies to statistical test results. Default: 'none'.
#' @param n \code{integer}. Return top N features/genes ranked by rankBy.
#'   Use NA (default) to return all results. Requires rankBy != 'none'.
#' @param filterFDR \code{numeric}. FDR threshold for significance filtering
#'   (0.0-1.0). Only results with adjusted p-value <= filterFDR retained.
#'   Default: NULL (no filtering).
#' @param format \code{character}. Output format: 'auto' (sensible default for type),
#'   'list', 'dataframe', or 'matrix'. Default: 'auto'.
#' @param display_table \code{logical}. For diversity results with display_table=TRUE,
#'   returns a formatted table showing diversity values across selected q-values
#'   for each gene. Default: FALSE (returns SummarizedExperiment or list).
#' @param n_genes \code{integer}. Number of genes to display in diversity tables
#'   when display_table=TRUE. Default: 4.
#' @param q_values_table \code{numeric}. Vector of q-values to include in diversity
#'   table display when display_table=TRUE. Default: c(0, 0.5, 1.0, 1.5, 2.0).
#' @param top_n \code{integer}. For effect_sizes_divergence, return top N genes ranked by sort_by.
#'   When specified, results are sorted by sort_by column and limited to top N rows.
#'   Default: NULL (return all results). Use NA to return all.
#' @param sort_by \code{character}. For effect_sizes_divergence, column name to sort by.
#'   Common choices: 'adj_p_interaction' (p-value, ascending), 'Mean_Divergence' (descending).
#'   Default: 'adj_p_interaction' (most significant first).
#' @return 
#'   - For diversity with q=NULL: A named list of SummarizedExperiment objects, one per q-value
#'   - For diversity with q specified: A single SummarizedExperiment for that q-value
#'   - For divergence: A SummarizedExperiment (rows=genes, columns=q-values), data.frame, or other format depending on divergence computation method
#'   - For lm/jackknife: A data.frame or list based on type and format
#'   - For pairwise: A data.frame with pairwise comparison difference metrics
#'   - For effect_sizes_divergence: A list containing effect size divergence results with components like interaction_results
#'   - For assumptions: A list containing rank-based assumption checks (exchangeability, monotonicity, consistency) and optional method-specific diagnostics (gam_metrics, gee_metrics, lmm_metrics, fpca_metrics)
#'   - For switching_tables: A list containing gene switching comparison tables
#'   Returns NULL if requested result type not computed or no results pass filtering.
#'
#' @details
#' This function provides flexible access to all computed results with ranking,
#' filtering, and format conversion. Compatible with DESeq2/edgeR design patterns
#' for familiar result extraction workflows.
#'
#' **Lazy Computation for switching_tables:**
#' When requesting \code{type = "switching_tables"}, the function automatically
#' computes and caches the tables if they don't exist yet but the prerequisites
#' do (LM and jackknife results). This eliminates the need for a separate
#' \code{prepare_gene_switching_tables_s4()} call - simply request the results
#' and they will be computed on-demand.
#'
#' @examples
#' # Load example data
#' data(readcounts, package = 'TSENAT')
#'
#' # Create TSENATAnalysis from count matrix
#' config <- TSENAT_config(
#'   q = 1.0,
#'   condition_col = 'group'
#' )
#' se <- SummarizedExperiment::SummarizedExperiment(
#'   assays = list(counts = readcounts),
#'   colData = data.frame(
#'     group = rep(c('A', 'B'), length.out = ncol(readcounts))
#'   )
#' )
#' analysis <- TSENATAnalysis(se = se, config = config)
#' analysis <- calculate_diversity(analysis)
#'
#' # Get all diversity results (list of SummarizedExperiment objects, one per q)
#' div_all <- results(analysis, type = 'diversity')
#'
#' # Get diversity for specific q-value (single SummarizedExperiment)
#' div_q1 <- results(analysis, type = 'diversity', q = 1.0)
#'
#' # Get results ranked by p-value, top 20 genes
#' # Using accessor function instead of @ slot access
#' top_lm <- results(analysis, type = 'lm', rankBy = 'pvalue', n = 20)
#'
#' # Get pairwise results (e.g., differential diversity metrics between conditions)
#' pairwise_diff <- results(analysis, type = 'pairwise')
#'
#' # Get effect size results, top 6 genes by p-value (most significant first)
#' top_effect_sizes <- results(analysis, type = 'effect_sizes_divergence', 
#'                              top_n = 6, sort_by = 'adj_p_interaction')
#'
#' # Get effect sizes sorted by mean divergence (largest effect sizes first)
#' large_effects <- results(analysis, type = 'effect_sizes_divergence',
#'                          top_n = 10, sort_by = 'Mean_Divergence')
#'
#' # Get switching tables - automatically computed if prerequisites exist
#' # (no need to call prepare_gene_switching_tables_s4 separately)
#' switching <- results(analysis, type = 'switching_tables')
#'
#' @rdname results
#' @export
results <- function(analysis, type, q = NULL, rankBy = "none", 
                       n = NA, filterFDR = NULL, format = "auto", display_table = FALSE,
                       n_genes = 4, q_values_table = c(0, 0.5, 1.0, 1.5, 2.0),
                       top_n = NULL, sort_by = "adj_p_interaction", sample = NULL) {
    # Validate parameters
    .validate_results_params(analysis, type, rankBy, format, filterFDR)
    
    # Extract result based on type
    result <- .extract_result_by_type(analysis, type)
    
    if (is.null(result)) {
        return(NULL)
    }
    
    # Warn about unsupported parameter combinations
    .warn_unsupported_params(type, filterFDR, rankBy)
    
    # Route to type-specific processor
    switch(type,
        diversity = .process_diversity_results(result, q, display_table, analysis, 
                                                n_genes, q_values_table, sample),
        divergence = .process_divergence_results(result, filterFDR, format),
        pairwise = .process_pairwise_results(result, filterFDR, format),
        lm = ,
        jackknife = ,
        rank_test = .process_statistical_results(result, type, filterFDR, rankBy, n, format),
        effect_sizes_divergence = .process_effect_sizes_divergence_results(result, top_n, sort_by, display_table, analysis),
        assumptions = result,
        switching_tables = .process_switching_tables_results(result, display_table),
        metadata = result,
        result
    )
}



# ============================================================================
# HELPER FUNCTIONS FOR RESULTS ACCESSOR
# ============================================================================

#' Internal: Validate results accessor parameters
#'
#' @noRd
.validate_results_params <- function(analysis, type, rankBy, format, filterFDR) {
    if (!is(analysis, "TSENATAnalysis")) {
        stop("'analysis' must be a TSENATAnalysis object", call. = FALSE)
    }
    
    valid_rank_methods <- c("none", "pvalue", "qvalue", "effectSize")
    if (!rankBy %in% valid_rank_methods) {
        stop("'rankBy' must be one of: ", paste(valid_rank_methods, collapse = ", "), 
             call. = FALSE)
    }
    
    valid_formats <- c("auto", "list", "dataframe", "matrix")
    if (!format %in% valid_formats) {
        stop("'format' must be one of: ", paste(valid_formats, collapse = ", "), 
             call. = FALSE)
    }
    
    if (!is.null(filterFDR) && (filterFDR < 0 || filterFDR > 1)) {
        stop("'filterFDR' must be between 0 and 1 or NULL", call. = FALSE)
    }
}

# ============================================================================
# HELPER: Get diversity result for specific q-value
# ============================================================================
.get_diversity_q_value <- function(result, q) {
    q_key <- NULL
    supported_decimals <- c(3, 2, 1, 0)
    
    for (decimals in supported_decimals) {
        candidate_key <- paste0("q_", formatC(q, format = "f", digits = decimals))
        if (candidate_key %in% names(result)) {
            q_key <- candidate_key
            break
        }
    }
    
    if (is.null(q_key)) {
        q_formatted <- formatC(q, format = "f", digits = 1)
        q_char_underscore <- paste0("q_", gsub("\\.", "_", q_formatted))
        if (q_char_underscore %in% names(result)) {
            q_key <- q_char_underscore
        }
    }
    
    if (is.null(q_key) && q == as.integer(q)) {
        candidate_key <- paste0("q_", as.integer(q))
        if (candidate_key %in% names(result)) {
            q_key <- candidate_key
        }
    }
    
    if (is.null(q_key)) {
        available_q <- vapply(names(result), function(x) {
            numeric_part <- sub("^q_", "", x)
            numeric_part <- gsub("_", ".", numeric_part)
            as.numeric(numeric_part)
        }, numeric(1))
        available_q <- available_q[!is.na(available_q)]
        available_q <- sort(unique(available_q))
        stop("Q-value ", q, " not found in results. Available q-values: ", 
             paste(available_q, collapse = ", "), call. = FALSE)
    }
    
    result[[q_key]]
}

# ============================================================================
# HELPER: Extract diversity table as a data.frame
# ============================================================================
.extract_diversity_table <- function(analysis, result, q, n_genes, q_values_table, sample = NULL) {
    all_div_results <- if (is.null(q)) analysis@diversity_results else list(result)
    
    if (length(all_div_results) == 0) {
        return(NULL)
    }
    
    first_se <- all_div_results[[1]]
    first_sample <- colnames(SummarizedExperiment::assay(first_se))[1]
    
    # Use specified sample if provided, validate it exists
    if (!is.null(sample)) {
        sample_names <- colnames(SummarizedExperiment::assay(first_se))
        if (!sample %in% sample_names) {
            stop("Sample '", sample, "' not found. Available samples: ",
                 paste(sample_names, collapse = ", "), call. = FALSE)
        }
        first_sample <- sample
    }
    
    # Build table as data.frame
    n_show <- min(n_genes, nrow(SummarizedExperiment::assay(first_se)))
    gene_names <- rownames(SummarizedExperiment::assay(first_se))[seq_len(n_show)]
    
    # Initialize with gene names
    table_df <- data.frame(Gene = gene_names, stringsAsFactors = FALSE)
    
    # Add columns for each q-value
    for (q_val in q_values_table) {
        q_name <- paste0("q_", sprintf("%.3f", q_val))
        col_name <- paste0("q_", sprintf("%.1f", q_val))
        
        if (q_name %in% names(all_div_results)) {
            mat <- SummarizedExperiment::assay(all_div_results[[q_name]])
            values <- mat[seq_len(n_show), first_sample]
            table_df[[col_name]] <- values
        }
    }
    
    # Add sample and n_genes as attributes
    attr(table_df, "sample") <- first_sample
    attr(table_df, "n_genes_total") <- nrow(SummarizedExperiment::assay(first_se))
    
    table_df
}

# ============================================================================
# HELPER: Display diversity table
# ============================================================================
.display_diversity_table <- function(analysis, result, q, n_genes, q_values_table, sample = NULL) {
    all_div_results <- if (is.null(q)) analysis@diversity_results else list(result)
    
    if (length(all_div_results) == 0) {
        return()
    }
    
    first_se <- all_div_results[[1]]
    first_sample <- colnames(SummarizedExperiment::assay(first_se))[1]
    
    # Use specified sample if provided, validate it exists
    if (!is.null(sample)) {
        sample_names <- colnames(SummarizedExperiment::assay(first_se))
        if (!sample %in% sample_names) {
            stop("Sample '", sample, "' not found. Available samples: ",
                 paste(sample_names, collapse = ", "), call. = FALSE)
        }
        first_sample <- sample
    }
    
    table_lines <- character()
    table_lines <- c(table_lines, "\n[results] Tsallis entropy across q-spectrum")
    table_lines <- c(table_lines, sprintf("[results] Sample: %s", first_sample))
    table_lines <- c(table_lines, sprintf("[results] Gene count: %d (showing %d)\n", 
                nrow(SummarizedExperiment::assay(first_se)), 
                min(n_genes, nrow(SummarizedExperiment::assay(first_se)))))
    
    header <- sprintf("%-15s", "Gene")
    for (q_val in q_values_table) {
        header <- paste0(header, sprintf("%12s", paste0("q=", sprintf("%.1f", q_val))))
    }
    table_lines <- c(table_lines, header)
    
    for (gene_idx in seq_len(min(n_genes, nrow(SummarizedExperiment::assay(first_se))))) {
        gene_name <- rownames(SummarizedExperiment::assay(first_se))[gene_idx]
        row_str <- sprintf("%-15s", gene_name)
        
        for (q_val in q_values_table) {
            q_name <- paste0("q_", sprintf("%.3f", q_val))
            if (q_name %in% names(all_div_results)) {
                mat <- SummarizedExperiment::assay(all_div_results[[q_name]])
                if (gene_idx <= nrow(mat)) {
                    val <- mat[gene_idx, first_sample]
                    row_str <- paste0(row_str, sprintf("%12.5f", val))
                }
            }
        }
        table_lines <- c(table_lines, row_str)
    }
    
    message(paste(table_lines, collapse = "\n"))
}

# ============================================================================
# HELPER: Extract result by type from analysis object
# ============================================================================
.extract_result_by_type <- function(analysis, type) {
    switch(type, 
        diversity = if (length(analysis@diversity_results) > 0) analysis@diversity_results else NULL,
        divergence = if (length(analysis@divergence_results) > 0) analysis@divergence_results else NULL,
        lm = if (length(analysis@lm_results) > 0) {
            if ("lm_interaction" %in% names(analysis@lm_results)) {
                analysis@lm_results$lm_interaction
            } else {
                analysis@lm_results
            }
        } else NULL,
        jackknife = .extract_jackknife_result(analysis),
        rank_test = if (!is.null(analysis@rank_test_results) && "rank_test" %in% names(analysis@rank_test_results)) {
            analysis@rank_test_results$rank_test
        } else NULL,
        pairwise = if (length(analysis@pairwise_results) > 0) analysis@pairwise_results else NULL,
        effect_sizes_divergence = .get_metadata_field(analysis, "effect_sizes_divergence"),
        assumptions = {
            meta <- .get_metadata_field(analysis, "rankbased_assumptions")
            if (!is.null(meta) && is.list(meta) && !is.null(meta$result)) {
                attr(meta$result, "checks")
            } else {
                NULL
            }
        },
        switching_tables = .extract_or_compute_switching_tables(analysis),
        metadata = analysis@metadata,
        stop("Unknown result type: '", type, "'. Must be one of: ", 
             "diversity, divergence, lm, jackknife, rank_test, pairwise, effect_sizes_divergence, assumptions, switching_tables, metadata", 
             call. = FALSE)
    )
}

# ============================================================================
# HELPER: Metadata accessor
# ============================================================================
.get_metadata_field <- function(analysis, field) {
    if (is.null(analysis@metadata)) {
        return(NULL)
    }
    analysis@metadata[[field]]
}

.set_metadata_field <- function(analysis, field, value) {
    if (is.null(analysis@metadata)) {
        analysis@metadata <- list()
    }
    analysis@metadata[[field]] <- value
    analysis
}

# ============================================================================
# HELPER: Extract jackknife result
# ============================================================================
.extract_jackknife_result <- function(analysis) {
    if (length(analysis@jackknife_results) == 0) {
        return(NULL)
    }
    
    jk_res <- NULL
    if (is.data.frame(analysis@jackknife_results)) {
        jk_res <- analysis@jackknife_results
    } else if ("results" %in% names(analysis@jackknife_results)) {
        jk_res <- analysis@jackknife_results$results
    } else if ("ci" %in% names(analysis@jackknife_results)) {
        jk_res <- analysis@jackknife_results$ci
    } else {
        jk_res <- analysis@jackknife_results
    }
    
    jk_res
}

# ============================================================================
# HELPER: Extract or compute switching tables
# ============================================================================
.extract_or_compute_switching_tables <- function(analysis) {
    existing_tables <- .get_metadata_field(analysis, "switching_tables")
    if (!is.null(existing_tables)) {
        return(existing_tables)
    }
    
    if (length(analysis@lm_results) == 0 || length(analysis@jackknife_results) == 0) {
        return(NULL)
    }
    
    tryCatch({
        lm_results_list <- analysis@lm_results
        lm_res <- if (!is.null(lm_results_list$lm_interaction)) {
            if (is.data.frame(lm_results_list$lm_interaction$results)) {
                lm_results_list$lm_interaction$results
            } else if (is.data.frame(lm_results_list$lm_interaction)) {
                lm_results_list$lm_interaction
            } else NULL
        } else if (is.data.frame(lm_results_list)) {
            lm_results_list
        } else NULL
        
        if (is.null(lm_res)) {
            return(NULL)
        }
        
        jk_list <- analysis@jackknife_results
        q_key_pattern <- "^q_[0-9]+_[0-9]{2}$"
        q_keyed <- jk_list[grep(q_key_pattern, names(jk_list))]
        
        if (length(q_keyed) == 0) {
            return(NULL)
        }
        
        computed_tables <- .prepare_gene_switching_tables(
            lm_res = lm_res, 
            multi_q_results = q_keyed,
            verbose = FALSE
        )
        
        analysis <- .set_metadata_field(analysis, "switching_tables", computed_tables)
        computed_tables
    }, error = function(e) {
        NULL
    })
}

# ============================================================================
# HELPER: Process diversity results
# ============================================================================
.process_diversity_results <- function(result, q, display_table, analysis, n_genes, 
                                       q_values_table, sample = NULL) {
    if (!is.null(q)) {
        result <- .get_diversity_q_value(result, q)
    }
    
    if (display_table) {
        .display_diversity_table(analysis, result, q, n_genes, q_values_table, sample)
        # Return the table as a data.frame instead of the full list
        table_df <- .extract_diversity_table(analysis, result, q, n_genes, q_values_table, sample)
        return(invisible(table_df))
    }
    
    result
}

# ============================================================================
# HELPER: Warn about unsupported parameters
# ============================================================================
.warn_unsupported_params <- function(type, filterFDR, rankBy) {
    if (type == "switching_tables" && (!is.null(filterFDR) || rankBy != "none")) {
        warning("rankBy and filterFDR are not supported for type='switching_tables'. ",
                "Ignoring these parameters. Switching tables are automatically pre-computed with optimal ",
                "ranking and filtering.", call. = FALSE)
    } else if (rankBy != "none" && !type %in% c("lm", "jackknife", "rank_test")) {
        warning("rankBy='", rankBy, "' is not supported for type='", type, "'. ",
                "Ignoring rankBy parameter. rankBy is only supported for types: ",
                "'lm', 'jackknife', 'rank_test'.", call. = FALSE)
    }
}

# ============================================================================
# HELPER: Handle jackknife multi-q extraction
# ============================================================================
.extract_jackknife_multi_q <- function(jk_res, q, rankBy) {
    if (is.list(jk_res) && !is.data.frame(jk_res)) {
        if (!is.null(q)) {
            q_char_underscore <- sprintf("q_%s", gsub("\\.", "_", sprintf("%.2f", q)))
            q_char_dot <- sprintf("q_%.2f", q)
            
            q_char <- if (q_char_underscore %in% names(jk_res)) q_char_underscore 
                      else if (q_char_dot %in% names(jk_res)) q_char_dot 
                      else NULL
            
            if (!is.null(q_char)) {
                jk_q_result <- jk_res[[q_char]]
                if (rankBy %in% c("pvalue", "qvalue")) {
                    jk_res <- jk_q_result
                } else if (is.list(jk_q_result) && "summary_table" %in% names(jk_q_result)) {
                    jk_res <- jk_q_result$summary_table
                } else if (is.data.frame(jk_q_result)) {
                    jk_res <- jk_q_result
                } else {
                    jk_res <- jk_q_result
                }
            } else {
                warning("Jackknife results for q=", q, " not found. Available q-values: ",
                        paste(grep("^q_", names(jk_res), value = TRUE), collapse = ", "),
                        call. = FALSE)
                jk_res <- NULL
            }
        } else if ("multi_q" %in% names(jk_res)) {
            jk_multi <- jk_res[["multi_q"]]
            if (rankBy %in% c("pvalue", "qvalue")) {
                jk_res <- jk_multi
            } else if (is.list(jk_multi) && "summary_table" %in% names(jk_multi)) {
                jk_res <- jk_multi$summary_table
            } else if (is.data.frame(jk_multi)) {
                jk_res <- jk_multi
            } else {
                jk_res <- jk_multi
            }
        } else if ("summary_table" %in% names(jk_res) && !(rankBy %in% c("pvalue", "qvalue"))) {
            jk_res <- jk_res$summary_table
        }
    }
    jk_res
}

# ============================================================================
# HELPER: Filter statistical results by FDR
# ============================================================================
.filter_statistical_by_fdr <- function(result, type, filterFDR) {
    if (is.null(filterFDR) || !is.data.frame(result)) {
        return(result)
    }
    
    padj_col <- switch(type,
        lm = if ("adj_p_interaction" %in% colnames(result)) "adj_p_interaction" else NULL,
        rank_test = if ("adj_p_value" %in% colnames(result)) "adj_p_value" else NULL,
        jackknife = if ("fdr" %in% colnames(result)) "fdr" 
                    else if ("delta_fdr" %in% colnames(result)) "delta_fdr" 
                    else NULL,
        NULL
    )
    
    if (is.null(padj_col)) {
        return(result)
    }
    
    result[!is.na(result[[padj_col]]) & result[[padj_col]] <= filterFDR, , drop = FALSE]
}

# ============================================================================
# HELPER: Rank and subset statistical results
# ============================================================================
.rank_statistical_results <- function(result, type, rankBy, n) {
    if (rankBy == "none") {
        return(result)
    }
    
    if (is.list(result) && !is.data.frame(result)) {
        if (type == "jackknife") {
            if (rankBy %in% c("pvalue", "qvalue") && "all_transcript_stats" %in% names(result)) {
                result <- result$all_transcript_stats
            } else if ("summary_table" %in% names(result) && is.data.frame(result$summary_table)) {
                result <- result$summary_table
            }
        } else {
            if ("summary_table" %in% names(result) && is.data.frame(result$summary_table)) {
                result <- result$summary_table
            } else if ("results" %in% names(result) && is.data.frame(result$results)) {
                result <- result$results
            }
        }
    }
    
    if (!is.data.frame(result)) {
        stop("Ranking requires data.frame results. Type '", type, "' returned different format.",
             call. = FALSE)
    }
    
    rank_col <- .get_ranking_column(type, rankBy, result)
    
    if (is.null(rank_col)) {
        warning("Column for rankBy='", rankBy, "' not found in results. Skipping ranking.",
               call. = FALSE)
        return(result)
    }
    
    idx <- if (rankBy == "effectSize") {
        order(abs(result[[rank_col]]), decreasing = TRUE, na.last = TRUE)
    } else {
        order(result[[rank_col]], na.last = TRUE)
    }
    result <- result[idx, , drop = FALSE]
    
    if (!is.na(n) && n > 0) {
        n <- min(n, nrow(result))
        result <- result[seq_len(n), , drop = FALSE]
    }
    
    result
}

# ============================================================================
# HELPER: Get ranking column name for statistical results
# ============================================================================
.get_ranking_column <- function(type, rankBy, result) {
    switch(type,
        lm = switch(rankBy,
            pvalue = if ("p_interaction" %in% colnames(result)) "p_interaction" else NULL,
            qvalue = if ("adj_p_interaction" %in% colnames(result)) "adj_p_interaction" else NULL,
            effectSize = if ("statistic" %in% colnames(result)) "statistic"
                        else if ("estimate" %in% colnames(result)) "estimate"
                        else if ("effect_size" %in% colnames(result)) "effect_size"
                        else NULL,
            NULL
        ),
        rank_test = switch(rankBy,
            pvalue = if ("p_value" %in% colnames(result)) "p_value" else NULL,
            qvalue = if ("adj_p_value" %in% colnames(result)) "adj_p_value" else NULL,
            effectSize = if ("statistic" %in% colnames(result)) "statistic"
                        else if ("estimate" %in% colnames(result)) "estimate"
                        else NULL,
            NULL
        ),
        jackknife = switch(rankBy,
            pvalue = if ("pvalue" %in% colnames(result)) "pvalue" else NULL,
            qvalue = if ("fdr" %in% colnames(result)) "fdr" else NULL,
            effectSize = if ("delta_influence" %in% colnames(result)) "delta_influence"
                        else if ("max_delta_influence" %in% colnames(result)) "max_delta_influence"
                        else NULL,
            NULL
        ),
        NULL
    )
}

# ============================================================================
# HELPER: Process statistical results (LM/Jackknife/RankTest)
# ============================================================================
.process_statistical_results <- function(result, type, filterFDR, rankBy, n, format) {
    if (is.null(result)) {
        return(NULL)
    }
    
    result <- .extract_statistical_dataframe(result, type)
    if (is.null(result)) {
        return(NULL)
    }
    
    result <- .filter_statistical_by_fdr(result, type, filterFDR)
    if (nrow(result) == 0) {
        return(NULL)
    }
    
    result <- .rank_statistical_results(result, type, rankBy, n)
    
    if (format != "auto") {
        result <- .convert_result_format(result, format, type)
    }
    
    result
}

# ============================================================================
# HELPER: Extract data.frame from statistical results
# ============================================================================
.extract_statistical_dataframe <- function(result, type) {
    if (is.data.frame(result)) {
        return(result)
    }
    
    if (!is.list(result)) {
        return(NULL)
    }
    
    if ("summary_table" %in% names(result) && is.data.frame(result$summary_table)) {
        result$summary_table
    } else if ("results" %in% names(result) && is.data.frame(result$results)) {
        result$results
    } else if ("all_transcript_stats" %in% names(result) && is.data.frame(result$all_transcript_stats)) {
        result$all_transcript_stats
    } else {
        NULL
    }
}

# ============================================================================
# HELPER: Process divergence results
# ============================================================================
.process_divergence_results <- function(result, filterFDR, format) {
    if (is.null(result)) {
        return(NULL)
    }
    
    if (methods::is(result, "SummarizedExperiment")) {
        result <- SummarizedExperiment::assay(result, "divergence")
    } else if (is.list(result) && length(result) > 0) {
        for (i in seq_along(result)) {
            if (methods::is(result[[i]], "SummarizedExperiment")) {
                result <- SummarizedExperiment::assay(result[[i]], "divergence")
                break
            }
        }
    }
    
    if (is.data.frame(result) || is.matrix(result)) {
        if (!is.null(filterFDR) && is.data.frame(result)) {
            padj_col <- if ("padj" %in% colnames(result)) "padj" else NULL
            if (!is.null(padj_col)) {
                result <- result[!is.na(result[[padj_col]]) & result[[padj_col]] <= filterFDR, , drop = FALSE]
                if (nrow(result) == 0) return(NULL)
            }
        }
        
        if (format != "auto") {
            result <- .convert_result_format(result, format, "divergence")
        }
    }
    
    result
}

# ============================================================================
# HELPER: Process pairwise results
# ============================================================================
.process_pairwise_results <- function(result, filterFDR, format) {
    if (is.null(result)) {
        return(NULL)
    }
    
    if (is.list(result)) {
        if ("difference" %in% names(result)) {
            result <- result$difference
        } else if (length(result) > 0) {
            result <- result[[1]]
        } else {
            return(NULL)
        }
    }
    
    if (is.data.frame(result)) {
        if (!is.null(filterFDR)) {
            padj_col <- if ("padj" %in% colnames(result)) "padj"
                       else if ("adj_p_value" %in% colnames(result)) "adj_p_value"
                       else NULL
            if (!is.null(padj_col)) {
                result <- result[!is.na(result[[padj_col]]) & result[[padj_col]] <= filterFDR, , drop = FALSE]
                if (nrow(result) == 0) return(NULL)
            }
        }
        
        if (format != "auto") {
            result <- .convert_result_format(result, format, "pairwise")
        }
    }
    
    result
}

# ============================================================================
# ============================================================================
# HELPER: Display effect sizes divergence table
# ============================================================================
.display_effect_sizes_table <- function(results_df, top_n) {
    if (is.null(results_df) || nrow(results_df) == 0) {
        return()
    }
    
    table_lines <- character()
    table_lines <- c(table_lines, "\n[results] Top genes by linear model significance with effect sizes and q-spectrum patterns")
    table_lines <- c(table_lines, "[results] Ranked by statistical significance (ascending P-values, Benjamini-Hochberg q-value < 0.05)\n")
    
    # Build header - select key columns for display
    header <- sprintf("%-20s", "Gene")
    header <- paste0(header, sprintf("%16s", "Slope Diff"))
    header <- paste0(header, sprintf("%16s", "Effect Size"))
    header <- paste0(header, sprintf("%18s", "LM p-value"))
    table_lines <- c(table_lines, header)
    
    # Build rows
    for (i in seq_len(min(top_n, nrow(results_df)))) {
        row_data <- results_df[i, ]
        
        # Gene name or ID
        gene_name <- if ("gene" %in% colnames(results_df)) {
            row_data$gene[1]
        } else if ("gene_name" %in% colnames(results_df)) {
            row_data$gene_name[1]
        } else if ("name" %in% colnames(results_df)) {
            row_data$name[1]
        } else {
            rownames(results_df)[i]
        }
        
        row_str <- sprintf("%-20s", gene_name)
        
        # Slope difference (main effect size metric)
        if ("slope_diff" %in% colnames(results_df)) {
            row_str <- paste0(row_str, sprintf("%16.4f", row_data$slope_diff[1]))
        } else {
            row_str <- paste0(row_str, sprintf("%16s", "N/A"))
        }
        
        # Effect size (use q=0.1 if available, otherwise q=0.05)
        effect_size_val <- NA
        if ("effect_size_D_q0_1" %in% colnames(results_df)) {
            effect_size_val <- row_data$effect_size_D_q0_1[1]
        } else if ("effect_size_D_q0_05" %in% colnames(results_df)) {
            effect_size_val <- row_data$effect_size_D_q0_05[1]
        } else if ("effect_size_D_q0_01" %in% colnames(results_df)) {
            effect_size_val <- row_data$effect_size_D_q0_01[1]
        }
        
        if (!is.na(effect_size_val)) {
            row_str <- paste0(row_str, sprintf("%16.4f", effect_size_val))
        } else {
            row_str <- paste0(row_str, sprintf("%16s", "N/A"))
        }
        
        # P-value (check for both column names)
        pval_col <- if ("adj_p_interaction" %in% colnames(results_df)) {
            "adj_p_interaction"
        } else if ("p_value_interaction" %in% colnames(results_df)) {
            "p_value_interaction"
        } else {
            NULL
        }
        
        if (!is.null(pval_col)) {
            pval <- row_data[[pval_col]][1]
            if (pval < 0.001) {
                pval_str <- sprintf("%.1e", pval)
            } else {
                pval_str <- sprintf("%.4f", pval)
            }
            row_str <- paste0(row_str, sprintf("%18s", pval_str))
        } else {
            row_str <- paste0(row_str, sprintf("%18s", "N/A"))
        }
        
        table_lines <- c(table_lines, row_str)
    }
    
    message(paste(table_lines, collapse = "\n"))
}

# ============================================================================
# HELPER: Process effect_sizes_divergence results with sorting and filtering
# ============================================================================
.process_effect_sizes_divergence_results <- function(result, top_n = NULL, sort_by = "adj_p_interaction", 
                                                      display_table = FALSE, analysis = NULL) {
    if (is.null(result)) {
        return(NULL)
    }
    
    # If result is already a data.frame, process directly (check this BEFORE is.list!)
    if (is.data.frame(result)) {
        results_df <- result
        
        if (!is.null(top_n) && !is.na(top_n) && top_n > 0) {
            # Verify sort_by column exists
            if (!(sort_by %in% colnames(results_df))) {
                stop("Column '", sort_by, "' not found in results. ",
                     "Available columns: ", paste(colnames(results_df), collapse = ", "),
                     call. = FALSE)
            }
            
            # Determine sort direction
            decreasing <- !grepl("p_value|pvalue|padj|adj_p", sort_by, ignore.case = TRUE)
            
            # Sort and limit
            order_idx <- order(results_df[[sort_by]], na.last = TRUE, decreasing = decreasing)
            results_df <- results_df[order_idx, , drop = FALSE]
            results_df <- head(results_df, top_n)
        }
        
        # Display table if requested
        if (display_table) {
            display_top_n <- if (!is.null(top_n) && !is.na(top_n)) top_n else nrow(results_df)
            .display_effect_sizes_table(results_df, display_top_n)
        }
        
        return(if (display_table) invisible(results_df) else results_df)
    }
    
    # If result is a list, try to extract the main results data.frame
    if (is.list(result)) {
        # Look for common data.frame names in the list
        if ("results" %in% names(result) && is.data.frame(result$results)) {
            results_df <- result$results
        } else if ("interaction_results" %in% names(result) && is.data.frame(result$interaction_results)) {
            results_df <- result$interaction_results
        } else if (length(result) == 1 && is.data.frame(result[[1]])) {
            results_df <- result[[1]]
        } else {
            # No data.frame found, return as-is
            return(result)
        }
        
        original_results_df <- results_df
        
        # Apply sorting and limiting if top_n specified
        if (!is.null(top_n) && !is.na(top_n) && top_n > 0) {
            # Verify sort_by column exists
            if (!(sort_by %in% colnames(results_df))) {
                stop("Column '", sort_by, "' not found in results. ",
                     "Available columns: ", paste(colnames(results_df), collapse = ", "),
                     call. = FALSE)
            }
            
            # Determine sort direction: ascending for p-values, descending for divergence/effect sizes
            decreasing <- !grepl("p_value|pvalue|padj|adj_p", sort_by, ignore.case = TRUE)
            
            # Sort results
            order_idx <- order(results_df[[sort_by]], na.last = TRUE, decreasing = decreasing)
            results_df <- results_df[order_idx, , drop = FALSE]
            
            # Limit to top_n
            results_df <- head(results_df, top_n)
        }
        
        # Display table if requested
        if (display_table) {
            display_top_n <- if (!is.null(top_n) && !is.na(top_n)) top_n else nrow(results_df)
            .display_effect_sizes_table(results_df, display_top_n)
        }
        
        # Return the processed data.frame (not the list)
        return(if (display_table) invisible(results_df) else results_df)
    }
    
    # Return as-is if it's some other type
    result
}

# ============================================================================
# Switching Tables Result Processing
# ============================================================================

#' Process switching_tables results with optional display
#'
#' @keywords internal
#' @noRd
.process_switching_tables_results <- function(result, display_table = FALSE) {
    if (is.null(result)) {
        return(NULL)
    }
    
    # Display formatted tables if requested
    if (display_table && is.list(result) && !is.null(result$comparison_tables)) {
        .display_switching_tables(result)
    }
    
    # Return the original result (invisibly if displayed)
    if (display_table) invisible(result) else result
}

#' Display switching tables in formatted output
#'
#' Displays transcript switching tables for each gene in a human-readable format,
#' with transcript delta-influence values across q-spectrum and direction consistency.
#'
#' @keywords internal
#' @noRd
.display_switching_tables <- function(switching_tables_result) {
    if (is.null(switching_tables_result)) {
        return()
    }
    
    comparison_tables <- switching_tables_result$comparison_tables
    gene_headers <- switching_tables_result$gene_headers
    q_metadata <- switching_tables_result$q_metadata
    
    if (is.null(comparison_tables) || length(comparison_tables) == 0) {
        return()
    }
    
    # Print each gene's switching table
    for (i in seq_along(comparison_tables)) {
        if (is.null(comparison_tables[[i]])) {
            next
        }
        
        # Print gene header
        message("")
        message(paste0("Gene: ", gene_headers[i]))
        
        table_data <- comparison_tables[[i]]
        
        # Get q-value metadata for this gene
        q_info <- q_metadata[[i]]
        q_values_available <- q_info$q_values_available
        q_key_to_value <- q_info$q_key_to_value
        
        # Build header with column names
        header_parts <- "Transcript"
        
        for (q_key in q_values_available) {
            q_val <- q_key_to_value[[q_key]]
            q_formatted <- sprintf("%.2f", q_val)
            header_parts <- c(header_parts, paste0("q=", q_formatted))
        }
        
        header_parts <- c(header_parts, "Direction Consistency")
        
        # Print header
        # Calculate column widths: transcript gets 20 chars, each q-value gets 10, consistency gets 20
        header_line <- sprintf("%-20s", header_parts[1])
        for (j in 2:(length(header_parts)-1)) {
            header_line <- paste0(header_line, sprintf("%10s", header_parts[j]))
        }
        header_line <- paste0(header_line, sprintf("  %-20s", header_parts[length(header_parts)]))
        message(header_line)
        
        # Print data rows
        for (row_idx in seq_len(nrow(table_data))) {
            tx_name <- table_data$transcript[row_idx]
            consistency <- table_data$Consistency[row_idx]
            
            if (is.na(consistency)) {
                consistency <- "Unknown"
            }
            
            row_line <- sprintf("%-20s", tx_name)
            
            # Add delta_influence values for each q
            for (q_key in q_values_available) {
                delta_val <- table_data[[q_key]][row_idx]
                if (is.na(delta_val)) {
                    row_line <- paste0(row_line, sprintf("%10s", "N/A"))
                } else {
                    row_line <- paste0(row_line, sprintf("%10.3f", delta_val))
                }
            }
            
            row_line <- paste0(row_line, sprintf("  %-20s", consistency))
            
            message(row_line)
        }
    }
}

