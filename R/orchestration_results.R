# ============================================================================
# RESULT ACCESSOR FUNCTIONS (EXPORTED)
# ============================================================================

#' Extract analysis results from TSENATAnalysis object
#'
#' Provides flexible access to diversity, divergence, and statistical test results
#' with options for ranking, filtering, and format conversion.
#'
#' @param analysis \code{TSENATAnalysis} object containing computed results.
#' @param type \code{character}. Type of results to extract:
#'   'diversity', 'divergence', 'lm', 'jackknife', 'rank_test', 'effect_sizes_divergence',
#'   or 'switching_tables'. Default: 'diversity'.
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
#'
#' @return 
#'   - For diversity with q=NULL: A named list of SummarizedExperiment objects, one per q-value
#'   - For diversity with q specified: A single SummarizedExperiment for that q-value
#'   - For divergence: A SummarizedExperiment (rows=genes, columns=q-values), data.frame, or other format depending on divergence computation method
#'   - For lm/jackknife: A data.frame or list based on type and format
#'   - For pairwise: A data.frame with pairwise comparison difference metrics
#'   - For effect_sizes_divergence: A list containing effect size divergence results with components like interaction_results
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
#' # Get switching tables - automatically computed if prerequisites exist
#' # (no need to call prepare_gene_switching_tables_s4 separately)
#' switching <- results(analysis, type = 'switching_tables')
#'
#' @rdname TSENATAnalysis-methods
#' @export
results <- function(analysis, type = "diversity", q = NULL, rankBy = "none", 
                       n = NA, filterFDR = NULL, format = "auto", display_table = FALSE,
                       n_genes = 4, q_values_table = c(0, 0.5, 1.0, 1.5, 2.0)) {
    if (!is(analysis, "TSENATAnalysis")) {
        stop("'analysis' must be a TSENATAnalysis object", call. = FALSE)
    }

    # Validate parameters
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

    # Extract result based on type
    result <- switch(type, 
        diversity = if (length(analysis@diversity_results) > 0) analysis@diversity_results else NULL,
        divergence = if (length(analysis@divergence_results) > 0) analysis@divergence_results else NULL,
        lm = if (length(analysis@lm_results) > 0) {
            # LM results are stored as list(lm_interaction = data.frame(...))
            if ("lm_interaction" %in% names(analysis@lm_results)) {
                analysis@lm_results$lm_interaction
            } else {
                analysis@lm_results
            }
        } else NULL,
        jackknife = if (length(analysis@jackknife_results) > 0) {
            # Jackknife results are stored as list, extract data.frame or summary_table if present
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
            
            # Handle multi-q jackknife: if we have nested q-value lists, extract specific q if specified
            if (is.list(jk_res) && !is.data.frame(jk_res)) {
                if (!is.null(q)) {
                    # User specified a q-value, extract that specific q-value result
                    # Use consistent %.2f format matching jackknife storage (see s4_functions_jis.R line 466)
                    q_char_underscore <- sprintf("q_%s", gsub("\\.", "_", sprintf("%.2f", q)))
                    q_char_dot <- sprintf("q_%.2f", q)  # Fallback for dot format
                    
                    q_char <- NULL
                    if (q_char_underscore %in% names(jk_res)) {
                        q_char <- q_char_underscore
                    } else if (q_char_dot %in% names(jk_res)) {
                        q_char <- q_char_dot
                    }
                    
                    if (!is.null(q_char)) {
                        jk_q_result <- jk_res[[q_char]]
                        # If rankBy is pvalue/qvalue, keep full structure for ranking logic to select table
                        # Otherwise, extract summary_table as default
                        if (rankBy %in% c("pvalue", "qvalue")) {
                            # Keep as list for ranking logic to select all_transcript_stats
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
                    # For multi-q jackknife without specific q parameter, try to use the multi_q result
                    jk_multi <- jk_res[["multi_q"]]
                    # If rankBy is pvalue/qvalue, keep full structure; otherwise extract summary_table
                    if (rankBy %in% c("pvalue", "qvalue")) {
                        jk_res <- jk_multi
                    } else if (is.list(jk_multi) && "summary_table" %in% names(jk_multi)) {
                        jk_res <- jk_multi$summary_table
                    } else if (is.data.frame(jk_multi)) {
                        jk_res <- jk_multi
                    } else {
                        # Keep original multi_q structure for now
                        jk_res <- jk_multi
                    }
                } else if ("summary_table" %in% names(jk_res)) {
                    # Single-q jackknife with summary_table
                    # Only extract if not ranking by pvalue/qvalue
                    if (rankBy %in% c("pvalue", "qvalue")) {
                        # Keep as list to allow ranking logic to select all_transcript_stats
                    } else {
                        jk_res <- jk_res$summary_table
                    }
                }
            }
            jk_res
        } else NULL,
        rank_test = if (!is.null(analysis@rank_test_results) && "rank_test" %in% names(analysis@rank_test_results)) {
            analysis@rank_test_results$rank_test
        } else NULL,
        pairwise = if (length(analysis@pairwise_results) > 0) analysis@pairwise_results else NULL,
        effect_sizes_divergence = S4Vectors::metadata(analysis)$effect_sizes_divergence,
        switching_tables = {
            # Lazy computation: compute tables if not already cached but prerequisites exist
            existing_tables <- metadata(analysis)$switching_tables
            if (!is.null(existing_tables)) {
                existing_tables
            } else if (length(analysis@lm_results) > 0 && length(analysis@jackknife_results) > 0) {
                # Prerequisites exist but tables not computed yet - compute and cache
                tryCatch({
                    # Compute tables using base function
                    lm_results_list <- analysis@lm_results
                    if (!is.null(lm_results_list$lm_interaction)) {
                        if (is.data.frame(lm_results_list$lm_interaction$results)) {
                            lm_res <- lm_results_list$lm_interaction$results
                        } else if (is.data.frame(lm_results_list$lm_interaction)) {
                            lm_res <- lm_results_list$lm_interaction
                        } else {
                            NULL
                        }
                    } else if (is.data.frame(lm_results_list)) {
                        lm_res <- lm_results_list
                    } else {
                        NULL
                    }
                    
                    if (!is.null(lm_res)) {
                        # Extract multi-q results
                        jk_list <- analysis@jackknife_results
                        q_key_pattern <- "^q_[0-9]+_[0-9]{2}$"
                        q_keyed <- jk_list[grep(q_key_pattern, names(jk_list))]
                        
                        if (length(q_keyed) > 0) {
                            multi_q_results <- q_keyed
                            
                            # Call base function
                            computed_tables <- .prepare_gene_switching_tables(
                                lm_res = lm_res, 
                                multi_q_results = multi_q_results,
                                verbose = FALSE
                            )
                            
                            # Cache in analysis metadata using the proper replacement method
                            meta <- metadata(analysis)
                            meta$switching_tables <- computed_tables
                            metadata(analysis) <- meta
                            computed_tables
                        } else {
                            NULL
                        }
                    } else {
                        NULL
                    }
                }, error = function(e) {
                    # If computation fails, return NULL silently
                    NULL
                })
            } else {
                NULL
            }
        },
        stop("Unknown result type: '", type, "'. Must be one of: ", 
             "diversity, divergence, lm, jackknife, rank_test, pairwise, effect_sizes_divergence, switching_tables", call. = FALSE)
    )

    if (is.null(result)) {
        return(NULL)
    }

    # Check for incompatible rankBy usage (type-specific validation)
    # Effect sizes divergence: more specific message
    if (type == "effect_sizes_divergence" && (!is.null(filterFDR) || rankBy != "none")) {
        warning("rankBy and filterFDR are not supported for type='effect_sizes_divergence'. ",
                "Ignoring these parameters. Access metadata directly for advanced filtering: ",
                "analysis@metadata$effect_sizes_divergence",
                call. = FALSE)
    } else if (type == "switching_tables" && (!is.null(filterFDR) || rankBy != "none")) {
        warning("rankBy and filterFDR are not supported for type='switching_tables'. ",
                "Ignoring these parameters. Switching tables are automatically pre-computed with optimal ",
                "ranking and filtering.",
                call. = FALSE)
    } else if (rankBy != "none" && !type %in% c("lm", "jackknife", "rank_test")) {
        # Generic warning for other unsupported types
        warning("rankBy='", rankBy, "' is not supported for type='", type, "'. ",
                "Ignoring rankBy parameter. rankBy is only supported for types: ",
                "'lm', 'jackknife', 'rank_test'.",
                call. = FALSE)
    }

    # === DIVERSITY RESULTS HANDLING ===
    if (type == "diversity") {
        # Filter by q-value if specified
        if (!is.null(q)) {
            # result is a list of SummarizedExperiment objects, one per q-value
            # Try multiple precision levels for robustness (matches diversity() accessor pattern)
            q_key <- NULL
            supported_decimals <- c(3, 2, 1, 0)  # Try 3 decimals first, then fewer
            
            for (decimals in supported_decimals) {
                candidate_key <- paste0("q_", formatC(q, format = "f", digits = decimals))
                if (candidate_key %in% names(result)) {
                    q_key <- candidate_key
                    break
                }
            }
            
            # If still not found, try underscore format (q_X_XXX instead of q_X.XXX)
            if (is.null(q_key)) {
                q_formatted <- formatC(q, format = "f", digits = 1)  # Default to 1 decimal
                q_char_underscore <- paste0("q_", gsub("\\.", "_", q_formatted))
                if (q_char_underscore %in% names(result)) {
                    q_key <- q_char_underscore
                }
            }
            
            # Last resort: try direct character conversion with various formats
            if (is.null(q_key)) {
                # Try as integer if q is whole number
                if (q == as.integer(q)) {
                    candidate_key <- paste0("q_", as.integer(q))
                    if (candidate_key %in% names(result)) {
                        q_key <- candidate_key
                    }
                }
            }
            
            if (is.null(q_key)) {
                # Format error message extracting just the numeric part
                available_q <- vapply(names(result), function(x) {
                    numeric_part <- sub("^q_", "", x)
                    numeric_part <- gsub("_", ".", numeric_part)
                    as.numeric(numeric_part)
                }, numeric(1))
                available_q <- available_q[!is.na(available_q)]
                available_q <- sort(unique(available_q))
                stop("Q-value ", q, " not found in results. Available q-values: ", 
                     paste(available_q, collapse = ", "), 
                     call. = FALSE)
            }
            
            # Return single SE for specified q-value
            result <- result[[q_key]]
        }
        
        # Display table if requested (diversity type)
        if (display_table) {
            # Get all diversity results as a list to access all q-values
            all_div_results <- if (is.null(q)) analysis@diversity_results else list(result)
            
            if (length(all_div_results) > 0) {
                first_se <- all_div_results[[1]]
                first_sample <- colnames(SummarizedExperiment::assay(first_se))[1]
                
                # Build table as character vector then message once
                table_lines <- character()
                table_lines <- c(table_lines, "\n[results] Tsallis entropy across q-spectrum")
                table_lines <- c(table_lines, sprintf("[results] Sample: %s", first_sample))
                table_lines <- c(table_lines, sprintf("[results] Gene count: %d (showing %d)\n", 
                            nrow(SummarizedExperiment::assay(first_se)), 
                            min(n_genes, nrow(SummarizedExperiment::assay(first_se)))))
                
                # Build table header
                header <- sprintf("%-15s", "Gene")
                for (q_val in q_values_table) {
                    header <- paste0(header, sprintf("%12s", paste0("q=", sprintf("%.1f", q_val))))
                }
                table_lines <- c(table_lines, header)
                
                # Display genes
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
        }
        
        return(result)
    }

    # === LM/STATISTICAL RESULTS HANDLING ===
    if (type %in% c("lm", "jackknife", "rank_test")) {
        
        # Filter by FDR if specified - type-specific column detection
        if (!is.null(filterFDR) && is.data.frame(result)) {
            padj_col <- if (type == "lm") {
                # LM results use adj_p_interaction
                if ("adj_p_interaction" %in% colnames(result)) "adj_p_interaction"
                else NULL
            } else if (type == "rank_test") {
                # rank_test results use adj_p_value
                if ("adj_p_value" %in% colnames(result)) "adj_p_value"
                else NULL
            } else if (type == "jackknife") {
                # Jackknife results use fdr (summary_table) or delta_fdr (per-transcript)
                if ("fdr" %in% colnames(result)) "fdr"
                else if ("delta_fdr" %in% colnames(result)) "delta_fdr"
                else NULL
            } else {
                NULL
            }
            
            if (!is.null(padj_col)) {
                result <- result[!is.na(result[[padj_col]]) & result[[padj_col]] <= filterFDR, , drop = FALSE]
                if (nrow(result) == 0) {
                    return(NULL)  # No results pass filter
                }
            }
        }

        # Rank and subset if requested
        if (rankBy != "none") {
            # For jackknife, check if we need to switch to transcript-level data for pvalue/qvalue ranking
            # This is needed because initial extraction returns summary_table by default
            if (type == "jackknife") {
                # If we have list structure with all_transcript_stats (not yet extracted)
                if (is.list(result) && !is.data.frame(result)) {
                    if (rankBy %in% c("pvalue", "qvalue") && "all_transcript_stats" %in% names(result)) {
                        result <- result$all_transcript_stats
                    } else if ("summary_table" %in% names(result) && is.data.frame(result$summary_table)) {
                        result <- result$summary_table
                    }
                }
                # Note: If result is already a data.frame, keep it as-is and proceed to ranking
                # (it should have the appropriate columns based on what was extracted)
            }
            
            # For non-jackknife list results, try to extract data.frame
            if (!is.data.frame(result) && is.list(result)) {
                # Try to extract data.frame from nested list
                if ("summary_table" %in% names(result) && is.data.frame(result$summary_table)) {
                    result <- result$summary_table
                } else if ("results" %in% names(result) && is.data.frame(result$results)) {
                    result <- result$results
                }
            }
            
            if (!is.data.frame(result)) {
                stop("Ranking requires data.frame results. Type '", type, "' returned different format.",
                     call. = FALSE)
            }
            
            # Determine ranking column based on result type with type-specific column detection
            # LM: p_interaction / adj_p_interaction / statistic (effect) + estimate + effect_size
            # rank_test: p_value / adj_p_value / statistic or estimate
            # jackknife: pvalue / fdr / delta_influence + max_delta_influence
            rank_col <- if (type == "lm") {
                # LM results
                switch(rankBy,
                    pvalue = if ("p_interaction" %in% colnames(result)) "p_interaction"
                            else NULL,
                    qvalue = if ("adj_p_interaction" %in% colnames(result)) "adj_p_interaction"
                            else NULL,
                    effectSize = if ("statistic" %in% colnames(result)) "statistic"
                                else if ("estimate" %in% colnames(result)) "estimate"
                                else if ("effect_size" %in% colnames(result)) "effect_size"
                                else NULL,
                    NULL
                )
            } else if (type == "rank_test") {
                # rank_test results
                switch(rankBy,
                    pvalue = if ("p_value" %in% colnames(result)) "p_value"
                            else NULL,
                    qvalue = if ("adj_p_value" %in% colnames(result)) "adj_p_value"
                            else NULL,
                    effectSize = if ("statistic" %in% colnames(result)) "statistic"
                                else if ("estimate" %in% colnames(result)) "estimate"
                                else NULL,
                    NULL
                )
            } else if (type == "jackknife") {
                # Jackknife (JIS) results
                switch(rankBy,
                    pvalue = if ("pvalue" %in% colnames(result)) "pvalue"
                            else NULL,
                    qvalue = if ("fdr" %in% colnames(result)) "fdr"
                            else NULL,
                    effectSize = if ("delta_influence" %in% colnames(result)) "delta_influence"
                                else if ("max_delta_influence" %in% colnames(result)) "max_delta_influence"
                                else NULL,
                    NULL
                )
            } else {
                NULL
            }
            
            if (is.null(rank_col)) {
                warning("Column for rankBy='", rankBy, "' not found in results. Skipping ranking.",
                       call. = FALSE)
            } else {
                # Sort by absolute value of effect sizes, or by p-value directly
                if (rankBy == "effectSize") {
                    idx <- order(abs(result[[rank_col]]), decreasing = TRUE, na.last = TRUE)
                } else {
                    idx <- order(result[[rank_col]], na.last = TRUE)
                }
                result <- result[idx, , drop = FALSE]
                
                # Subset to top N if specified
                if (!is.na(n) && n > 0) {
                    n <- min(n, nrow(result))
                    result <- result[seq_len(n), , drop = FALSE]
                }
            }
        }

        # Convert to requested format
        if (format != "auto") {
            result <- .convert_result_format(result, format, type)
        }
        
        return(result)
    }

    # === DIVERGENCE RESULTS HANDLING ===
    if (type == "divergence") {
        # Extract actual divergence values if result is a SummarizedExperiment
        if (methods::is(result, "SummarizedExperiment")) {
            result <- SummarizedExperiment::assay(result, "divergence")
        }
        # If result is a list with SummarizedExperiment (e.g., divergence_se element)
        else if (is.list(result) && length(result) > 0) {
            # Try to extract SummarizedExperiment from list
            for (i in seq_along(result)) {
                if (methods::is(result[[i]], "SummarizedExperiment")) {
                    result <- SummarizedExperiment::assay(result[[i]], "divergence")
                    break
                }
            }
        }
        
        # Filtering and format conversion for data.frame results
        if (is.data.frame(result) || is.matrix(result)) {
            if (!is.null(filterFDR)) {
                if (is.data.frame(result)) {
                    padj_col <- if ("padj" %in% colnames(result)) "padj" else NULL
                    if (!is.null(padj_col)) {
                        result <- result[!is.na(result[[padj_col]]) & result[[padj_col]] <= filterFDR, , drop = FALSE]
                        if (nrow(result) == 0) return(NULL)
                    }
                }
            }
            
            if (format != "auto") {
                result <- .convert_result_format(result, format, type)
            }
        }
        
        return(result)
    }

    # === PAIRWISE RESULTS HANDLING ===
    if (type == "pairwise") {
        # Pairwise results are stored as list(difference = data.frame(...))
        # Extract the difference component directly
        if (is.list(result)) {
            if ("difference" %in% names(result)) {
                result <- result$difference
            } else if (length(result) > 0) {
                # Fallback: return first element if difference not found
                result <- result[[1]]
            } else {
                return(NULL)
            }
        }
        
        # Apply filtering and formatting
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
                result <- .convert_result_format(result, format, type)
            }
        }
        
        return(result)
    }

    # === DIVERGENCE RESULTS HANDLING ===
    if (type == "divergence") {
        # Extract actual divergence values if result is a SummarizedExperiment
        if (methods::is(result, "SummarizedExperiment")) {
            result <- SummarizedExperiment::assay(result, "divergence")
        }
        # If result is a list with SummarizedExperiment (e.g., divergence_se element)
        else if (is.list(result) && length(result) > 0) {
            # Try to extract SummarizedExperiment from list
            for (i in seq_along(result)) {
                if (methods::is(result[[i]], "SummarizedExperiment")) {
                    result <- SummarizedExperiment::assay(result[[i]], "divergence")
                    break
                }
            }
        }
        
        # Filtering and format conversion for data.frame results
        if (is.data.frame(result) || is.matrix(result)) {
            if (!is.null(filterFDR)) {
                if (is.data.frame(result)) {
                    padj_col <- if ("padj" %in% colnames(result)) "padj" else NULL
                    if (!is.null(padj_col)) {
                        result <- result[!is.na(result[[padj_col]]) & result[[padj_col]] <= filterFDR, , drop = FALSE]
                        if (nrow(result) == 0) return(NULL)
                    }
                }
            }
            
            if (format != "auto") {
                result <- .convert_result_format(result, format, type)
            }
        }
        
        return(result)
    }

    return(result)
}
