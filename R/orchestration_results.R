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
#'   'assumptions', 'concordance', or 'switching_tables'.
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
#' @param format \code{character}. Output format: 'text' (pre-formatted character vector,
#'   default) or 'list' (structured components). For diversity results with format='se', 
#'   returns the SummarizedExperiment object directly (useful for downstream processing 
#'   with other packages like SplicingFactory). For concordance results, 'text' returns 
#'   pre-formatted character vector for display, while 'list' returns structured components 
#'   (summary_table, agreement_dist, high_conf, etc.) suitable for custom display.
#'   Default: 'text'.
#' @param n_genes \code{integer}. Number of genes to display in diversity results.
#'   Default: 4.
#' @param q_values_table \code{numeric}. Vector of q-values to include in diversity
#'   results. Default: c(0, 0.5, 1.0, 1.5, 2.0).
#' @param top_n \code{integer}. For effect_sizes_divergence, return top N genes ranked by sort_by.
#'   When specified, results are sorted by sort_by column and limited to top N rows.
#'   Default: NULL (return all results). Use NA to return all.
#' @param sort_by \code{character}. For effect_sizes_divergence, column name to sort by.
#'   Common choices: 'adj_p_interaction' (p-value, ascending), 'Mean_Divergence' (descending).
#'   Default: 'adj_p_interaction' (most significant first).
#' @return 
#'   - For diversity with q=NULL: A named list of SummarizedExperiment objects, one per q-value
#'   - For diversity with q specified and format='auto' or 'table': A data.frame table for display
#'   - For diversity with q specified and format='se': The SummarizedExperiment object directly (useful for SplicingFactory)
#'   - For divergence: A SummarizedExperiment (rows=genes, columns=q-values), data.frame, or other format depending on divergence computation method
#'   - For lm/jackknife: A data.frame or list based on type and format
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
#' # Get diversity for specific q-value (single SummarizedExperiment table for display)
#' div_q1 <- results(analysis, type = 'diversity', q = 1.0)
#'
#' # Get diversity for specific q-value as SummarizedExperiment for downstream processing (e.g., SplicingFactory)
#' div_q1_se <- results(analysis, type = 'diversity', q = 1.0, format = 'se')
#'
#' # Get results ranked by p-value, top 20 genes
#' # Using accessor function instead of @ slot access
#' top_lm <- results(analysis, type = 'lm', rankBy = 'pvalue', n = 20)
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
                       n = NA, filterFDR = NULL, format = "auto",
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
        diversity = .process_diversity_results(result, q, analysis, 
                                                n_genes, q_values_table, sample, format),
        divergence = .process_divergence_results(result, filterFDR, format),
        lm = ,
        jackknife = ,
        rank_test = .process_statistical_results(result, type, filterFDR, rankBy, n, format),
        effect_sizes_divergence = .process_effect_sizes_divergence_results(result, top_n, sort_by, analysis),
        assumptions = .process_assumptions_results(result, format = format),
        switching_tables = .process_switching_tables_results(result),
        concordance = .process_concordance_results(result, format = format),
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
    
    valid_formats <- c("auto", "list", "dataframe", "matrix", "se", "table")
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
# HELPER: Extract result by type from analysis object
# ============================================================================
.extract_result_by_type <- function(analysis, type) {
    switch(type, 
        diversity = if (length(analysis@diversity_results) > 0) analysis@diversity_results else NULL,
        divergence = if (length(analysis@divergence_results) > 0) analysis@divergence_results else NULL,
        concordance = .get_metadata_field(analysis, "method_concordance"),
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
             "diversity, divergence, lm, jackknife, rank_test, effect_sizes_divergence, assumptions, switching_tables, concordance, metadata", 
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
.process_diversity_results <- function(result, q, analysis, n_genes, 
                                       q_values_table, sample = NULL, format = "auto") {
    # Handle format = "se": Return SummarizedExperiment for downstream processing
    # Useful for passing to other packages like SplicingFactory
    if (format == "se" && !is.null(q)) {
        result_se <- .get_diversity_q_value(result, q)
        return(result_se)
    }
    
    # Default (format = "auto" or "table"): Return formatted table for display
    if (!is.null(q)) {
        result <- .get_diversity_q_value(result, q)
    }
    
    # Return formatted table without automatic display
    table_df <- .extract_diversity_table(analysis, result, q, n_genes, q_values_table, sample)
    table_df
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
# HELPER: Display effect sizes divergence table
# ============================================================================
# ============================================================================
# HELPER: Process effect_sizes_divergence results with sorting and filtering
# ============================================================================
.process_effect_sizes_divergence_results <- function(result, top_n = NULL, sort_by = "adj_p_interaction", 
                                                      analysis = NULL) {
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
        
        # Return the processed data frame (silently)
        return(results_df)
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
        
        # Return the processed data.frame (silently)
        return(results_df)
    }
    
    # Return as-is if it's some other type
    result
}

# ============================================================================
# Assumptions Result Processing
# ============================================================================

#' Process assumptions results into formatted table
#'
#' Converts nested assumption check list into formatted data frame
#' with characteristic, test, result, and interpretation columns.
#'
#' @keywords internal
#' @noRd
.process_assumptions_results <- function(result, format = "text") {
    if (is.null(result)) {
        return(NULL)
    }
    
    # Check if result is a list of assumption checks
    if (!is.list(result)) {
        return(result)
    }
    
    # Build table as data frame
    rows <- list()
    
    # Helper to format p-values and test statistics
    format_value <- function(x) {
        if (is.null(x)) return("N/A")
        # Handle vectors: extract first element
        if (length(x) > 1) x <- x[1]
        if (is.na(x)) return("N/A")
        if (is.numeric(x)) {
            if (x < 0.001) return(sprintf("%.0e", x))
            if (x < 0.01) return(sprintf("%.4f", x))
            return(sprintf("%.3f", x))
        }
        return(as.character(x))
    }
    
    # Helper to capitalize first letter and remove underscores
    capitalize_first <- function(x) {
        x_clean <- gsub("_", " ", x)
        paste0(toupper(substring(x_clean, 1, 1)), substring(x_clean, 2))
    }
    
    # ========== Rank-based tests (always present) ==========
    
    # Exchangeability (core assumption for rank-based tests)
    if (!is.null(result$exchangeability)) {
        check <- result$exchangeability
        status_val <- if (!is.null(check$status)) check$status else "unknown"
        rows[[length(rows) + 1]] <- list(
            Characteristic = "Exchangeability",
            Test = "Permutation test",
            Result = paste0("p=", format_value(check$p_value)),
            Interpretation = paste0(toupper(substring(status_val, 1, 1)), substring(status_val, 2))
        )
    }
    
    # Monotonicity (rank ordering consistency)
    if (!is.null(result$monotonicity)) {
        check <- result$monotonicity
        r_val <- format_value(check$mean_correlation)
        interp <- if (!is.null(check$mean_correlation)) {
            if (abs(check$mean_correlation) < 0.3) "Heterogeneous" else "Homogeneous"
        } else "Unknown"
        rows[[length(rows) + 1]] <- list(
            Characteristic = "Monotonicity",
            Test = "Spearman ρ",
            Result = paste0("r=", r_val),
            Interpretation = interp
        )
    }
    
    # Consistency (replicate agreement across samples)
    if (!is.null(result$consistency)) {
        check <- result$consistency
        # Consistency stores: kendall_w and icc_simplified (not w_statistic and icc)
        w_val <- format_value(check$kendall_w)
        icc_val <- format_value(check$icc_simplified)
        interp <- if (!is.na(check$icc_simplified)) {
            if (check$icc_simplified < 0.5) "Low" else "Moderate"
        } else "Unknown"
        rows[[length(rows) + 1]] <- list(
            Characteristic = "Consistency",
            Test = "Kendall's W / ICC",
            Result = paste0("W=", w_val, ", ICC=", icc_val),
            Interpretation = interp
        )
    }
    
    # ========== GAM-specific tests (corrected key names) ==========
    
    if (!is.null(result$gam_metrics)) {
        gam_metrics <- result$gam_metrics
        
        # Concurvity: actual key is gam_metrics$concurvity, not concurvity_index
        if (!is.null(gam_metrics$concurvity)) {
            metric <- gam_metrics$concurvity
            has_error <- isTRUE(metric$error)
            # Key is overall_concurvity, not concurvity_index
            val_valid <- !is.null(metric$overall_concurvity) && !is.na(metric$overall_concurvity)
            if (has_error || !val_valid) {
                result_str <- "NA"
            } else {
                result_str <- sprintf("%.3f", metric$overall_concurvity)
            }
            interp <- if (!val_valid || has_error) "Unknown" else if (metric$overall_concurvity < 0.5) "Low" else "High"
            rows[[length(rows) + 1]] <- list(
                Characteristic = "Concurvity",
                Test = "Smooth collinearity",
                Result = result_str,
                Interpretation = interp
            )
        }
        
        # Effective DoF: actual key is edf with edf_ratio field
        if (!is.null(gam_metrics$edf)) {
            metric <- gam_metrics$edf
            has_error <- isTRUE(metric$error)
            val_valid <- !is.null(metric$edf_ratio) && !is.na(metric$edf_ratio)
            if (has_error || !val_valid) {
                result_str <- "NA"
            } else {
                result_str <- sprintf("%.3f", metric$edf_ratio)
            }
            interp <- if (!val_valid || has_error) "Unknown" else if (metric$edf_ratio < 0.1) "Over-smoothed" else if (metric$edf_ratio > 0.9) "Under-smoothed" else "Adequate"
            rows[[length(rows) + 1]] <- list(
                Characteristic = "EDF Ratio",
                Test = "Smoothing",
                Result = result_str,
                Interpretation = interp
            )
        }
        
        # Non-linearity: actual key is nonlinearity with r2_improvement_percent field
        if (!is.null(gam_metrics$nonlinearity)) {
            metric <- gam_metrics$nonlinearity
            has_error <- isTRUE(metric$error)
            val_valid <- !is.null(metric$r2_improvement_percent) && !is.na(metric$r2_improvement_percent)
            if (has_error || !val_valid) {
                result_str <- "NA"
            } else {
                result_str <- sprintf("%.1f%%", metric$r2_improvement_percent)
            }
            interp <- if (!val_valid || has_error) "Unknown" else if (metric$r2_improvement_percent < 1) "Use linear" else "Use GAM"
            rows[[length(rows) + 1]] <- list(
                Characteristic = "Non-linearity",
                Test = "ΔR² vs LM",
                Result = result_str,
                Interpretation = interp
            )
        }
        
        # Basis Adequacy: actual key is basis_adequacy with optimal_basis_dimension field
        if (!is.null(gam_metrics$basis_adequacy)) {
            metric <- gam_metrics$basis_adequacy
            has_error <- isTRUE(metric$error)
            val_valid <- !is.null(metric$optimal_basis_dimension) && !is.na(metric$optimal_basis_dimension)
            if (has_error || !val_valid) {
                result_str <- "NA"
            } else {
                result_str <- sprintf("k=%d", metric$optimal_basis_dimension)
            }
            interp <- if (!val_valid || has_error) "Unknown" else "Adequate"
            rows[[length(rows) + 1]] <- list(
                Characteristic = "Basis Dimension",
                Test = "Spline basis",
                Result = result_str,
                Interpretation = interp
            )
        }
    }
    
    # ========== GEE-specific tests ==========
    
    if (!is.null(result$gee_metrics)) {
        gee_metrics <- result$gee_metrics
        
        # Extract details from each GEE metric (correlation_fit, cluster_variation, etc.)
        for (metric_name in names(gee_metrics)) {
            if (metric_name == "consolidated") next  # Skip consolidated summary
            
            metric <- gee_metrics[[metric_name]]
            if (!is.list(metric)) next
            
            # Format result: just extract the essential value from details
            result_val <- if (!is.null(metric$details)) {
                # Extract just the numeric value if possible
                gsub("<.*?>|\\s+\\(.*\\)", "", metric$details)  # Remove HTML/parenthetical info
            } else "N/A"
            
            test_val <- if (!is.null(metric$method)) sub("^[^:]*:\\s*", "", metric$method) else capitalize_first(metric_name)
            interp_val <- if (!is.null(metric$status)) gsub("✓|⚠|✗", "", metric$status) else "unknown"
            
            row_item <- list(
                Characteristic = capitalize_first(metric_name),
                Test = capitalize_first(test_val),
                Result = substr(result_val, 1, 50),  # Truncate to 50 chars max
                Interpretation = paste0(toupper(substring(interp_val, 1, 1)), substring(interp_val, 2))
            )
            rows[[length(rows) + 1]] <- row_item
        }
    }
    
    # ========== LMM-specific tests ==========
    
    if (!is.null(result$lmm_metrics)) {
        lmm_metrics <- result$lmm_metrics
        
        # Extract details from each LMM metric (variance_components, normality, etc.)
        for (metric_name in names(lmm_metrics)) {
            if (metric_name == "consolidated") next  # Skip consolidated summary
            
            metric <- lmm_metrics[[metric_name]]
            if (!is.list(metric)) next
            
            # Format result: just extract the essential value from details
            result_val <- if (!is.null(metric$details)) {
                gsub("<.*?>|\\s+\\(.*\\)", "", metric$details)  # Remove HTML/parenthetical info
            } else "N/A"
            
            test_val <- if (!is.null(metric$method)) sub("^[^:]*:\\s*", "", metric$method) else capitalize_first(metric_name)
            interp_val <- if (!is.null(metric$status)) gsub("✓|⚠|✗", "", metric$status) else "unknown"
            
            row_item <- list(
                Characteristic = capitalize_first(metric_name),
                Test = capitalize_first(test_val),
                Result = substr(result_val, 1, 50),  # Truncate to 50 chars max
                Interpretation = paste0(toupper(substring(interp_val, 1, 1)), substring(interp_val, 2))
            )
            rows[[length(rows) + 1]] <- row_item
        }
    }
    
    # ========== FPCA-specific tests ==========
    
    if (!is.null(result$fpca_metrics)) {
        fpca_metrics <- result$fpca_metrics
        
        # Extract details from each FPCA metric (variance_adequacy, bootstrap_stability, etc.)
        for (metric_name in names(fpca_metrics)) {
            if (metric_name == "consolidated") next  # Skip consolidated summary
            
            metric <- fpca_metrics[[metric_name]]
            if (!is.list(metric)) next
            
            # Format result: just extract the essential value from details
            result_val <- if (!is.null(metric$details)) {
                gsub("<.*?>|\\s+\\(.*\\)", "", metric$details)  # Remove HTML/parenthetical info
            } else "N/A"
            
            test_val <- if (!is.null(metric$method)) sub("^[^:]*:\\s*", "", metric$method) else capitalize_first(metric_name)
            interp_val <- if (!is.null(metric$status)) gsub("✓|⚠|✗", "", metric$status) else "unknown"
            
            row_item <- list(
                Characteristic = capitalize_first(metric_name),
                Test = capitalize_first(test_val),
                Result = substr(result_val, 1, 50),  # Truncate to 50 chars max
                Interpretation = paste0(toupper(substring(interp_val, 1, 1)), substring(interp_val, 2))
            )
            rows[[length(rows) + 1]] <- row_item
        }
    }
    
    # Convert list of rows to data frame
    if (length(rows) == 0) {
        return(NULL)
    }
    
    df <- do.call(rbind, lapply(rows, as.data.frame, stringsAsFactors = FALSE))
    rownames(df) <- NULL
    
    # Return format based on format parameter
    if (format == "list") {
        # Return structured list
        return(list(
            assumptions_table = df,
            raw_result = result
        ))
    }
    
    # Default (format = "text"): Return formatted text
    output_lines <- c("\nRank-Based Test Assumptions\n")
    output_lines <- c(output_lines, .format_data_frame_as_text(df))
    
    # Combine all lines into single text string
    formatted_text <- paste(output_lines, collapse = "")
    
    # Add custom class for printing
    class(formatted_text) <- c("assumptions_text", "character")
    formatted_text
}

# ============================================================================
# Custom Print Method for assumptions_text
# ============================================================================

#' Print method for assumptions_text class
#'
#' Properly displays formatted assumptions text by interpreting newline characters
#'
#' @param x An object of class \code{assumptions_text}
#' @param ... Additional arguments (unused)
#'
#' @return Invisibly returns the input object
#'
#' @exportS3Method base::print
print.assumptions_text <- function(x, ...) {
    cat(x)
    invisible(x)
}

# ============================================================================
# Switching Tables Result Processing
# ============================================================================

#' Process switching_tables results
#'
#' @keywords internal
#' @noRd
.process_switching_tables_results <- function(result) {
    if (is.null(result)) {
        return(NULL)
    }
    
    # Result is already formatted as a list of data frames per gene
    # Simply return it
    result
}

#' Display switching tables in formatted output
#'
#' Displays transcript switching tables for each gene in a human-readable format,
#' with transcript delta-influence values across q-spectrum and direction consistency.
#'
#' @keywords internal
#' @noRd
# ============================================================================
# HELPER: Process concordance results
# ============================================================================

#' @noRd
.process_concordance_results <- function(result, format = "text") {
    if (is.null(result)) {
        return(NULL)
    }
    
    # If result is a plain data frame, return as-is
    if (is.data.frame(result)) {
        return(result)
    }
    
    # If result is a list with concordance data, build comprehensive summary
    if (!is.list(result) || is.null(result$comparison_df)) {
        return(result)
    }
    
    comparison_df <- result$comparison_df
    spearman_rho <- result$spearman_rho %||% NA
    high_conf <- result$high_conf %||% data.frame()
    agreement_table <- result$agreement_table %||% table()
    
    n_total <- nrow(comparison_df)
    
    # Calculate agreement statistics
    both_sig <- sum(comparison_df$agreement == "Both significant", na.rm = TRUE)
    lm_only <- sum(comparison_df$agreement == "LM only", na.rm = TRUE)
    rank_only <- sum(comparison_df$agreement == "Rank test only", na.rm = TRUE)
    neither <- sum(comparison_df$agreement == "Neither significant", na.rm = TRUE)
    
    # Calculate rates
    concordance_rate <- if (n_total > 0) (both_sig / n_total) * 100 else 0
    discordance_rate <- if (n_total > 0) ((lm_only + rank_only) / n_total) * 100 else 0
    
    # ====== 1. SUMMARY METRICS TABLE ======
    summary_table <- data.frame(
        Metric = c(
            "Total genes compared",
            "Spearman correlation (p-values)",
            "Both methods significant (p < 0.05)",
            "LM only significant",
            "Rank test only significant",
            "Neither significant",
            "Concordance rate",
            "Discordance rate"
        ),
        Value = c(
            paste0(n_total),
            paste0("rho = ", sprintf("%.4f", spearman_rho)),
            paste0(both_sig, " (", sprintf("%.1f%%", (both_sig/n_total)*100), ")"),
            paste0(lm_only, " (", sprintf("%.1f%%", (lm_only/n_total)*100), ")"),
            paste0(rank_only, " (", sprintf("%.1f%%", (rank_only/n_total)*100), ")"),
            paste0(neither, " (", sprintf("%.1f%%", (neither/n_total)*100), ")"),
            paste0(sprintf("%.1f%%", concordance_rate)),
            paste0(sprintf("%.1f%%", discordance_rate))
        ),
        stringsAsFactors = FALSE
    )
    
    # ====== 2. AGREEMENT DISTRIBUTION TABLE ======
    agreement_dist <- data.frame(
        "Agreement Category" = c(
            "Both significant",
            "LM only",
            "Rank test only",
            "Neither significant"
        ),
        "Number of Genes" = c(both_sig, lm_only, rank_only, neither),
        "Percentage" = c(
            sprintf("%.1f%%", (both_sig/n_total)*100),
            sprintf("%.1f%%", (lm_only/n_total)*100),
            sprintf("%.1f%%", (rank_only/n_total)*100),
            sprintf("%.1f%%", (neither/n_total)*100)
        ),
        stringsAsFactors = FALSE,
        check.names = FALSE
    )
    
    # ====== 3. HIGH-CONFIDENCE GENES TABLE ======
    high_conf_table <- NULL
    if (!is.null(high_conf) && nrow(high_conf) > 0) {
        high_conf_table <- data.frame(
            Gene = high_conf$gene,
            "LM adj p" = sapply(high_conf$padj_lm, function(x) {
                if (x < 1e-50) sprintf("%.2e", x) else sprintf("%.3e", x)
            }),
            "Rank test adj p" = sapply(high_conf$padj_rank, function(x) {
                if (x < 1e-50) sprintf("%.2e", x) else sprintf("%.3e", x)
            }),
            "LM Effect" = sprintf("%.1f%%", high_conf$effect_lm * 100),
            "Rank test η²" = sprintf("%.3f", high_conf$effect_rank),
            stringsAsFactors = FALSE,
            check.names = FALSE
        )
    }
    
    # ====== 4. ALL GENES TABLE ======
    all_genes_table <- data.frame(
        Gene = comparison_df$gene,
        "LM p" = comparison_df$p_lm,
        "LM adj p" = comparison_df$padj_lm,
        "Rank test p" = comparison_df$p_rank,
        "Rank test adj p" = comparison_df$padj_rank,
        "LM Effect" = comparison_df$effect_lm,
        "Rank test η²" = comparison_df$effect_rank,
        Agreement = comparison_df$agreement,
        stringsAsFactors = FALSE,
        check.names = FALSE
    )
    
    # ====== BUILD FORMATTED TEXT OUTPUT ======
    output_lines <- c()
    
    output_lines <- c(output_lines, "\nGlobal Concordance Metrics: LM vs Scheirer-Ray-Hare Methods\n")
    output_lines <- c(output_lines, .format_data_frame_as_text(summary_table))
    
    output_lines <- c(output_lines, "\nMethod Agreement Distribution\n")
    output_lines <- c(output_lines, .format_data_frame_as_text(agreement_dist))
    
    if (!is.null(high_conf_table) && nrow(high_conf_table) > 0) {
        n_hc <- nrow(high_conf_table)
        output_lines <- c(output_lines, 
            paste0("\nRobust Entropic Order Index Interactions: High-Confidence Genes Detected by Both Methods (n=", n_hc, ", ranked by statistical significance)\n"))
        output_lines <- c(output_lines, .format_data_frame_as_text(high_conf_table))
    }
    
    output_lines <- c(output_lines, "\nAll Genes with Agreement Classification\n")
    output_lines <- c(output_lines, .format_data_frame_as_text(all_genes_table))
    
    # Combine all lines into single text string
    formatted_text <- paste(output_lines, collapse = "")
    
    # Return based on format parameter
    if (format == "list") {
        # Return structured list with all components
        return(list(
            summary_table = summary_table,
            agreement_dist = agreement_dist,
            high_conf = high_conf,
            high_conf_table = high_conf_table,
            all_genes_table = all_genes_table,
            spearman_rho = spearman_rho,
            concordance_rate = concordance_rate,
            discordance_rate = discordance_rate
        ))
    }
    
    # Default (format = "text"): Return formatted text
    # Add custom class for printing
    class(formatted_text) <- c("concordance_text", "character")
    formatted_text
}

# ============================================================================
# HELPER: Format data frame as aligned text
# ============================================================================

#' @noRd
.format_data_frame_as_text <- function(df) {
    if (nrow(df) == 0) {
        return("")
    }
    
    # Convert to character for formatting
    df_char <- as.data.frame(lapply(df, as.character), stringsAsFactors = FALSE)
    
    # Calculate column widths
    col_widths <- sapply(seq_len(ncol(df_char)), function(j) {
        max(nchar(colnames(df_char)[j]), max(nchar(df_char[[j]])))
    })
    
    # Format header
    header <- paste(
        mapply(function(name, width) {
            sprintf("%-*s", width, name)
        }, colnames(df_char), col_widths),
        collapse = " "
    )
    
    # Format rows
    rows <- apply(df_char, 1, function(row) {
        paste(
            mapply(function(val, width) {
                sprintf("%-*s", width, val)
            }, row, col_widths),
            collapse = " "
        )
    })
    
    # Combine and add newlines
    text_lines <- c(header, rows, "")
    paste(text_lines, collapse = "\n")
}

# ============================================================================
# DISPLAY METHOD: Print method for concordance_text class
# ============================================================================

#' @exportS3Method base::print
print.concordance_text <- function(x, ...) {
    cat(x)
    invisible(x)
}

