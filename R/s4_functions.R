# S4 Wrapper Functions for TSENAT Pipeline These functions provide
# S4-integrated alternatives to existing analysis functions. They extract input
# from TSENATAnalysis slots, run analysis, and store results back to
# appropriate slots.



#' Calculate LM interactions and store in TSENATAnalysis
#'
#' @param analysis \code{TSENATAnalysis} object.
#' @param fdr_threshold \code{numeric}. FDR cutoff for significance.
#'   Default: 0.05.
#' @param formula \code{formula} or NULL. Reserved for future use.
#' @param condition_col \code{character} or \code{NULL}. Column name in colData identifying sample conditions.
#'   If NULL, reads from \code{@config$condition_col} or auto-detects common column names.
#' @param method \code{character}. Statistical method (e.g., 'lmm', 'gam', 'gee').
#'   If NULL, uses method from @config$method or defaults to 'lmm'.
#' @param paired \code{logical}. Whether to use paired design. Default: FALSE.
#'   If not specified, reads from \code{@config$paired} if available.
#' @param subject_col \code{character} or \code{NULL}. Column name identifying subject IDs for paired designs.
#'   If NULL, reads from \code{@config$subject_col} if available.
#' @param nthreads \code{numeric} or \code{NULL}. Number of CPU threads for parallel processing.
#'   If NULL, reads from \code{@config$nthreads} (or defaults to NULL, letting base function decide).
#' @param multicorr \code{character} or \code{NULL}. Multiple comparison correction method.
#'   Options: 'hochberg', 'westfall-young', 'benjamini-yekutieli'.
#'   If NULL, uses method from @config or base function defaults.
#' @param corstr \code{character} or \code{NULL}. Correlation structure for GEE models.
#'   Options: 'ar1', 'exchangeable', 'independence'.
#'   If NULL, uses method from @config or base function defaults.
#' @param pcorr \code{character} or \code{NULL}. P-value correction method.
#'   Default: 'BH' (Benjamini-Hochberg).
#'   If NULL, reads from \code{@config$pcorr} if available.
#' @param verbose \code{logical}. Print progress messages. Default: FALSE.
#' @param return_model_data \code{logical}. Return model data for visualization. Default: TRUE.
#' @param output_file \code{character} or \code{NULL}. Optional file path to save results.
#'   Supported formats: .rds (for S4 objects), .tsv, .csv, .txt (for tables). Default: NULL (no file output).
#' @param ... Additional arguments passed to the base LM function,
#'   including: pvalue, min_obs, assay_name, bias_correction, regularization, storey, wy_randomizations, adaptive_knots, etc.
#'
#' @return Modified TSENATAnalysis with results in @lm_results$lm_interaction.
#'
#' @details
#' Extracts diversity results from @diversity_results (prerequisite),
#' combines across q-values into single SummarizedExperiment,
#' then runs \code{.calculate_lm_interaction()}.
#' 
#' **Parameter Priority Resolution:**
#' \itemize{
#'   \item \code{nthreads}: Priority: explicit > @config > NULL
#' }
#' 
#' Parameters are resolved in priority order:
#' 1. Explicit arguments passed to function
#' 2. Values from analysis@config (if present)
#' 3. Function defaults
#'
#' @examples
#' # Create test analysis with appropriate sample structure for paired design
#' # (requires lme4 package for LMM fitting; uses synthetic data for reproducibility)
#' set.seed(42)
#' 
#' # Create transcript-level counts with biological signal
#' # Note: Use adequate complexity (transcripts/genes, samples, expression) 
#' # to avoid filtering away all genes during diversity computation
#' n_genes <- 50
#' n_transcripts_per_gene <- 30
#' n_transcripts <- n_genes * n_transcripts_per_gene
#' n_samples <- 16  # 8 subjects \u00D7 2 conditions (paired design)
#' 
#' # Generate counts with clear biological signal
#' control_idx <- seq(1, n_samples, by = 2)
#' treatment_idx <- seq(2, n_samples, by = 2)
#' 
#' counts <- matrix(0, nrow = n_transcripts, ncol = n_samples)
#' for (j in seq_len(n_samples)) {
#'   lambda <- if (j %in% control_idx) 100 else 180
#'   counts[, j] <- rpois(n_transcripts, lambda = lambda)
#' }
#' counts <- pmax(counts, 50)  # Ensure minimum expression
#' 
#' rownames(counts) <- paste0('TX_', seq_len(n_transcripts))
#' colnames(counts) <- paste0('Sample_', seq_len(n_samples))
#' 
#' # Create rowData with gene mapping (tx2gene structure)
#' rowdata <- data.frame(
#'   transcript_id = rownames(counts),
#'   gene_id = rep(paste0('GENE_', 1:n_genes), 
#'                 each = n_transcripts_per_gene),
#'   row.names = rownames(counts)
#' )
#' 
#' # Create colData with paired design metadata
#' coldata <- data.frame(
#'   sample_id = colnames(counts),
#'   condition = rep(c('control', 'treatment'), 
#'                   length.out = n_samples),
#'   subject = rep(paste0('Subject_', 1:8), 
#'                 length.out = n_samples),
#'   row.names = colnames(counts)
#' )
#' 
#' # Build SummarizedExperiment
#' se <- SummarizedExperiment::SummarizedExperiment(
#'   assays = list(counts = counts),
#'   rowData = S4Vectors::DataFrame(rowdata),
#'   colData = S4Vectors::DataFrame(coldata)
#' )
#' 
#' # Add tx2gene metadata for gene-level aggregation
#' S4Vectors::metadata(se)$tx2gene <- 
#'   data.frame(Transcript = rowdata$transcript_id,
#'              Gene = rowdata$gene_id)
#' 
#' # Initialize TSENATAnalysis
#' analysis <- TSENATAnalysis(se = se, config = list())
#' 
#' # Compute diversity (prerequisite for LM interaction analysis)
#' analysis <- calculate_diversity_s4(
#'   analysis, 
#'   q = c(0.5, 1.0, 1.5), 
#'   verbose = FALSE
#' )
#' 
#' # Calculate q \u00D7 condition interactions using GAM
#' # (avoids convergence issues with LMM)
#' analysis <- calculate_lm_interaction_s4(
#'   analysis,
#'   condition_col = 'condition',
#'   method = 'gam',
#'   verbose = FALSE
#' )
#' 
#' # View interaction test results
#' head(lmResults(analysis, 'lm_interaction'))
#'
#' @export
#' @importFrom utils write.table
calculate_lm_interaction_s4 <- function(analysis, fdr_threshold = NULL, formula = NULL,
    condition_col = NULL, method = NULL, paired = FALSE, subject_col = NULL, nthreads = NULL,
    multicorr = NULL, corstr = NULL, pcorr = NULL, verbose = FALSE, return_model_data = TRUE,
    output_file = NULL, ...) {
    if (!is(analysis, "TSENATAnalysis")) {
        stop("'analysis' must be a TSENATAnalysis object", call. = FALSE)
    }

    # OPTIMIZATION: Clear memoization cache for new analysis Ensures fresh
    # calculations for new dataset
    .clear_lm_helper_cache()

    # Check prerequisites
    if (length(analysis@diversity_results) == 0) {
        stop("Diversity results required. Run calculate_diversity_s4() first.", call. = FALSE)
    }

    # Ensure colData consistency: If diversity results have colData, sync to
    # analysis@se This maintains pairing structure through the analysis
    # pipeline
    if (length(analysis@diversity_results) > 0) {
        # Get colData from first diversity result (should be consistent across
        # q-values)
        first_diversity_se <- analysis@diversity_results[[1]]
        if (is(first_diversity_se, "SummarizedExperiment") && ncol(first_diversity_se) >
            0) {
            diversity_coldata <- SummarizedExperiment::colData(first_diversity_se)
            if (!is.null(diversity_coldata) && nrow(diversity_coldata) == ncol(analysis@se)) {
                # If diversity results have richer colData, use that
                if (!is.null(diversity_coldata) && ncol(diversity_coldata) > 0) {
                  SummarizedExperiment::colData(analysis@se) <- diversity_coldata
                }
            }
        }
    }

    # =========================================================================
    # PARAMETER EXTRACTION FROM @config using centralized helpers
    # =========================================================================

    # Auto-detect condition_col if not provided
    if (is.null(condition_col)) {
        cd_cols <- colnames(colData(analysis@se))
        condition_col <- auto_detect_column(cd_cols, config_list = analysis@config,
            config_key = "condition_col", priority_candidates = c("condition", "sample_type",
                "group", "treatment"), default_fallback = NULL, verbose = verbose,
            param_name = "condition_col")
    }

    # Resolve remaining parameters using centralized handler
    method <- resolve_slot_param(method, analysis@config, "method", "lmm")
    subject_col <- resolve_slot_param(subject_col, analysis@config, "subject_col",
        NULL)
    multicorr <- resolve_slot_param(multicorr, analysis@config, "multicorr", NULL)
    corstr <- resolve_slot_param(corstr, analysis@config, "corstr", NULL)
    pcorr <- resolve_slot_param(pcorr, analysis@config, "pcorr", "BH")
    nthreads <- resolve_slot_param(nthreads, analysis@config, "nthreads", NULL)

    # Paired design (logical flag)
    if (!paired && "paired" %in% names(analysis@config)) {
        paired <- analysis@config$paired
    }

    # Validate that required parameters are available
    if (is.null(condition_col)) {
        cd_cols <- colnames(colData(analysis@se))
        if (length(cd_cols) > 0) {
            message("condition_col not specified. Available columns: ", paste(cd_cols,
                collapse = ", "))
        } else {
            message("condition_col not specified and colData is empty. ", "Will be determined by .calculate_lm_interaction().")
        }
    }

    # Extract q-values from diversity_results keys (format: 'q_0.5', 'q_1.0',
    # etc.)
    q_keys <- names(analysis@diversity_results)

    # =========================================================================
    # APPROACH: Combine existing per-q diversity results NO recomputation -
    # just stack the already-computed SEs horizontally Add _q= suffix to column
    # names to distinguish q-values
    # =========================================================================

    # Combine per-q SEs into single multi-q SE for calculate_lm_interaction
    diversity_combined <- tryCatch({
        # Extract assay matrices and colData from each per-q SE
        assay_list <- list()
        coldata_list <- list()
        rowdata_first <- NULL

        for (key in sort(q_keys)) {
            se <- analysis@diversity_results[[key]]
            if (!is(se, "SummarizedExperiment")) {
                stop("[calculate_lm_interaction_s4] Diversity result ", key, " is not a SummarizedExperiment",
                  call. = FALSE)
            }

            # Extract q-value from key (e.g., 'q_0.100' -> 0.1)
            q_val_str <- sub("^q_", "", key)
            q_val <- as.numeric(q_val_str)

            # Get assay matrix (should be named 'diversity')
            if (length(SummarizedExperiment::assays(se)) == 0) {
                stop("[calculate_lm_interaction_s4] Diversity SE ", key, " has no assays",
                  call. = FALSE)
            }
            assay_matrix <- as.matrix(SummarizedExperiment::assay(se))

            # Add _q= suffix to ASSAY column names ONLY to distinguish q-values
            # DO NOT add suffix to colData column names (rbind requires
            # identical column names across all DFrames)
            assay_colnames_with_q <- paste0(colnames(assay_matrix), "_q=", q_val_str)
            colnames(assay_matrix) <- assay_colnames_with_q

            assay_list[[key]] <- assay_matrix

            # Get colData - IMPORTANT: colData column names must remain the
            # same across all q-values Only update rownames to match new assay
            # column names
            coldata <- SummarizedExperiment::colData(se)
            rownames(coldata) <- assay_colnames_with_q
            coldata_list[[key]] <- coldata

            # Save rowData from first SE (same genes for all q-values)
            if (is.null(rowdata_first)) {
                rowdata_first <- SummarizedExperiment::rowData(se)
            }
        }

        # Combine assays horizontally (columns from different q-values)
        combined_assay <- do.call(cbind, assay_list)

        # Combine colData vertically (each q-value's samples with q-suffix in
        # rownames)
        combined_coldata <- do.call(rbind, coldata_list)

        # Ensure colnames of assay match rownames of colData
        colnames(combined_assay) <- rownames(combined_coldata)

        # Create combined SE
        se_combined <- SummarizedExperiment::SummarizedExperiment(assays = list(diversity = combined_assay),
            colData = combined_coldata, rowData = rowdata_first)

        se_combined
    }, error = function(e) {
        stop("[calculate_lm_interaction_s4] Failed to combine diversity results: ",
            conditionMessage(e), call. = FALSE)
    })

    # Build arguments for calculate_lm_interaction
    args <- list(se = diversity_combined)

    if (!is.null(method)) {
        args$method <- method
    }

    # Add extracted parameters
    if (!is.null(condition_col)) {
        args$condition_col <- condition_col
    }

    if (paired) {
        args$paired <- paired
    }

    if (!is.null(subject_col)) {
        args$subject_col <- subject_col
    }

    if (!is.null(nthreads)) {
        args$nthreads <- nthreads
    }

    if (!is.null(multicorr)) {
        args$multicorr <- multicorr
    }

    if (!is.null(corstr)) {
        args$corstr <- corstr
    }

    if (!is.null(pcorr)) {
        args$pcorr <- pcorr
    }

    if (verbose) {
        args$verbose <- verbose
    }

    # Request model_data for plotting compatibility
    args$return_model_data <- return_model_data

    # Merge with additional args (which may override these values)
    args <- c(args, list(...))

    # Run LM analysis
    result <- tryCatch({
        do.call(.calculate_lm_interaction, args)
    }, error = function(e) {
        stop("lm_interaction calculation failed:\n", conditionMessage(e), call. = FALSE)
    })

    # =========================================================================
    # HANDLE RESULT: May be list($results, $model_data) or just data.frame
    # =========================================================================
    lm_results_df <- result
    model_data <- NULL

    # Check if result is a list with $results and $model_data
    # (return_model_data=TRUE)
    if (is.list(result) && "results" %in% names(result)) {
        lm_results_df <- result$results
        model_data <- result$model_data
    }

    # =========================================================================
    # CRITICAL VALIDATION: Ensure result has required structure
    # =========================================================================
    # The result MUST be a data.frame with columns 'gene' and
    # 'adj_p_interaction' to be compatible with downstream functions like
    # .effect_sizes_divergence()
    if (!is.data.frame(lm_results_df)) {
        stop("[calculate_lm_interaction_s4] Result from calculate_lm_interaction must be a data.frame, ",
            "but got: ", class(lm_results_df), ". This indicates an unexpected change in the underlying function.",
            call. = FALSE)
    }

    # Check if result is empty (0 rows OR 0 columns)
    if (nrow(lm_results_df) == 0 || ncol(lm_results_df) == 0) {
        warning("[calculate_lm_interaction_s4] Result is empty (", nrow(lm_results_df),
            " rows, ", ncol(lm_results_df), " columns). ", "This can occur with: low sample counts per condition, ",
            "insufficient signal, or model convergence issues. ", "Try: increasing samples, using higher lambda for data generation, ",
            "or checking colData grouping structure.", call. = FALSE)
        # Return empty results gracefully rather than error
        analysis@lm_results <- list(lm_interaction = data.frame())
        return(analysis)
    }

    required_cols <- c("gene", "adj_p_interaction")
    missing_cols <- setdiff(required_cols, colnames(lm_results_df))
    if (length(missing_cols) > 0) {
        stop("[calculate_lm_interaction_s4] Result missing required columns: ", paste(missing_cols,
            collapse = ", "), ". Available columns: ", paste(colnames(lm_results_df),
            collapse = ", "), call. = FALSE)
    }

    # Store results with model_data
    if (is.list(analysis@lm_results) && "lm_interaction" %in% names(analysis@lm_results)) {
        # Preserve other lm_results components
        analysis@lm_results$lm_interaction <- lm_results_df
        if (!is.null(model_data)) {
            analysis@lm_results$lm_interaction_model_data <- model_data
        }
    } else {
        analysis@lm_results <- list(lm_interaction = lm_results_df)
        if (!is.null(model_data)) {
            analysis@lm_results$lm_interaction_model_data <- model_data
        }
    }

    # Track metadata
    analysis@metadata$function_calls <- c(analysis@metadata$function_calls, "calculate_lm_interaction")

    # Save if output_file provided
    if (!is.null(output_file)) {
        lm_data <- NULL

        # Try multiple ways to extract data
        if (!is.null(analysis@lm_results)) {
            if (is.list(analysis@lm_results) && "lm_interaction" %in% names(analysis@lm_results)) {
                lm_int <- analysis@lm_results$lm_interaction
                if (is.data.frame(lm_int)) {
                  lm_data <- lm_int
                } else if (is.list(lm_int) && "results" %in% names(lm_int)) {
                  lm_data <- lm_int$results
                }
            }
        }

        if (!is.null(lm_data) && is.data.frame(lm_data)) {
            save_analysis_output(lm_data, output_file, object = analysis, verbose = verbose,
                func_name = "calculate_lm_interaction_s4")
        }
    }

    analysis
}

# ============================================================================
# JACKKNIFE WRAPPER
# ============================================================================

#' Jackknife resampling with confidence intervals
#'
#' @param analysis \code{TSENATAnalysis} object.
#' @param q \code{numeric}. Q-value(s) for jackknife. Default: 1.0.
#' @param norm \code{logical}. Normalization flag. Default: NULL (uses @config$norm or TRUE).
#' @param log_base \code{numeric}. Logarithm base for entropy normalization. Default: NULL (uses e).
#' @param seed \code{numeric} or \code{NULL}. Random seed for reproducibility. Default: NULL.
#' @param top_n \code{numeric}. Number of top outlier samples to report. Default: 5.
#' @param pseudocount \code{numeric}. Pseudocount value for count regularization. Default: NULL (uses @config$pseudocount or 0).
#' @param verbose \code{logical}. Print jackknife results summary. Default: FALSE.
#' @param nthreads \code{numeric} or \code{NULL}. Number of CPU threads for parallel processing.
#'   If NULL, reads from \code{@config$nthreads} (or defaults to 1).
#'   If > 1 and multiple q-values provided, uses parallel PSOCK cluster.
#' @param output_file \code{character} or \code{NULL}. Optional file path to save results.
#'   Supported formats: .tsv, .csv, .txt (for jackknife results table with estimates, influence, outliers),
#'   .rds (for entire S4 object). Default: NULL (no file output).
#' @param ... Additional arguments passed to the base function.
#'
#' @return Modified TSENATAnalysis with jackknife results in @jackknife_results.
#'
#' @details
#' Requires diversity results to exist first. Will error if
#' \code{calculate_diversity_s4()} has not been run.
#'
#' **Parameter Priority Resolution:**
#' \describe{
#'   \item{nthreads}{Priority: explicit > \code{@config$nthreads} > 1}
#' }
#'
#' @examples
#' # Load example data (matching TSENAT.Rmd workflow)
#' data(readcounts)
#' readcounts <- as.matrix(salmon_dataset)
#' mode(readcounts) <- 'numeric'
#' metadata_df <- read.table(
#'   system.file('extdata', 'metadata.tsv', package = 'TSENAT'),
#'   header = TRUE, sep = '\t'
#' )
#' gff3_dataset <- system.file('extdata', 'annotation.gff3.gz', package = 'TSENAT')
#' 
#' # Build analysis from vignette data
#' analysis <- build_analysis_s4(readcounts = readcounts, tx2gene = gff3_dataset, metadata = metadata_df,
#'   tpm = salmon_tpm, effective_length = salmon_effective_length)
#'
#' # Filter low-abundance genes (required for reliable jackknife estimates)
#' analysis <- filter_analysis_s4(analysis, stringency = 'severe')
#' 
#' # Compute diversity first (required for jackknife)
#' analysis <- calculate_diversity_s4(analysis, q = c(0.5, 1.0, 1.5),
#'   verbose = FALSE)
#' 
#' # Run jackknife estimation
#' analysis <- jackknife_entropy_outliers_s4(analysis, q = c(0.5, 1.0, 1.5),
#'   verbose = FALSE)
#' # Check jackknife results
#' names(jackKnife(analysis))
#'
#' @export
#' @importFrom utils write.table
jackknife_entropy_outliers_s4 <- function(analysis, q = NULL, norm = NULL, log_base = NULL,
    seed = NULL, top_n = NULL, verbose = FALSE, nthreads = NULL, pseudocount = NULL,
    output_file = NULL, ...) {
    if (!is(analysis, "TSENATAnalysis")) {
        stop("'analysis' must be a TSENATAnalysis object", call. = FALSE)
    }

    # Check prerequisite: diversity must be calculated
    if (length(analysis@diversity_results) == 0) {
        stop("Diversity results required. Run calculate_diversity_s4() first.", call. = FALSE)
    }

    # PARAMETER EXTRACTION using utility function
    q <- resolve_slot_param(q, analysis@config, "q_values", 1)
    norm <- resolve_slot_param(norm, analysis@config, "norm", TRUE)
    log_base <- resolve_slot_param(log_base, analysis@config, "log_base", exp(1))
    seed <- resolve_slot_param(seed, analysis@config, "seed", NULL)
    top_n <- resolve_slot_param(top_n, analysis@config, "top_n", 5)
    nthreads <- resolve_slot_param(nthreads, analysis@config, "nthreads", 1)
    pseudocount <- resolve_slot_param(pseudocount, analysis@config, "pseudocount",
        0)

    # Ensure q is numeric
    if (!is.numeric(q)) {
        stop("'q' must be numeric", call. = FALSE)
    }

    for (q_val in q) {
        # Verify diversity was computed for this q-value (prerequisite for
        # jackknife)
        div_key <- paste0("q_", formatC(q_val, format = "f", digits = 3))
        if (!(div_key %in% names(analysis@diversity_results))) {
            stop("Diversity not calculated for q=", q_val, ". Run calculate_diversity_s4(analysis, q=",
                q_val, ") first.", call. = FALSE)
        }

        # Run jackknife - pass the COUNTS (not diversity values!) to jackknife
        # function Jackknife stability analysis requires raw count data, not
        # pre-computed diversity
        tryCatch({
            # Extract counts matrix from SummarizedExperiment
            counts_matrix <- SummarizedExperiment::assay(analysis@se, "counts")

            result <- .jackknife_entropy_outliers(x = counts_matrix, q = q_val, norm = norm,
                log_base = log_base, seed = seed, top_n = top_n, pseudocount = pseudocount,
                verbose = verbose, nthreads = nthreads, ...)

            # Store with key 'q_X.XXX' (consistent 3 decimal formatting)
            jk_key <- paste0("q_", formatC(q_val, format = "f", digits = 3))
            analysis@jackknife_results[[jk_key]] <- result

            # Track metadata
            analysis@metadata$function_calls <- c(analysis@metadata$function_calls,
                paste0("jackknife_tsallis_entropy[q=", q_val, "]"))
        }, error = function(e) {
            stop("Jackknife computation failed for q=", q_val, ":\n", e$message,
                call. = FALSE)
        })
    }

    # Save if output_file provided
    if (!is.null(output_file)) {
        # Create directory if it doesn't exist
        output_dir <- dirname(output_file)
        if (output_dir != "." && !dir.exists(output_dir)) {
            dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
        }

        if (grepl("\\.tsv$|\\.csv$|\\.txt$", tolower(output_file))) {
            # Convert jackknife results to data.frame for text output
            tryCatch({
                # Extract all jackknife results and convert to data.frame
                all_results <- list()
                for (jk_key in names(analysis@jackknife_results)) {
                  jk_result <- analysis@jackknife_results[[jk_key]]

                  # Handle both single result and list of results
                  if (inherits(jk_result, "tsenat_jackknife")) {
                    # Single result - wrap in list for uniform processing
                    jk_result <- list(jk_result)
                    names(jk_result) <- "gene1"
                  }

                  if (inherits(jk_result, "tsenat_jackknife_list")) {
                    # Convert list of results to data.frame
                    df_list <- lapply(names(jk_result), function(gene_name) {
                      res <- jk_result[[gene_name]]
                      data.frame(gene = gene_name, estimate = res$estimate, jackknife_se = res$jackknife_se,
                        n_transcripts = res$n_transcripts, n_outliers = length(res$outlier_indices),
                        outlier_indices = paste(res$outlier_indices, collapse = ","),
                        outlier_threshold = res$outlier_threshold, outlier_cutoff_value = res$outlier_cutoff_value,
                        q_value = res$q, normalized = res$norm, stringsAsFactors = FALSE)
                    })

                    df <- do.call(rbind, df_list)
                    df$q_key <- jk_key  # Add q-value key for multi-q results
                    rownames(df) <- NULL
                    all_results[[jk_key]] <- df
                  }
                }

                # Combine all results into single data.frame
                if (length(all_results) > 0) {
                  output_df <- do.call(rbind, all_results)
                  rownames(output_df) <- NULL

                  # Determine separator based on file extension
                  sep <- if (grepl("\\.csv$", tolower(output_file)))
                    "," else "\t"

                  write.table(output_df, file = output_file, sep = sep, quote = FALSE,
                    row.names = FALSE)

                  if (verbose) {
                    message("[jackknife_entropy_outliers_s4] Results saved to ",
                      output_file)
                  }
                }
            }, error = function(e) {
                warning("[jackknife_entropy_outliers_s4] Could not write jackknife results to file: ",
                  conditionMessage(e), call. = FALSE)
            })
        } else {
            # Default to RDS for S4 object
            saveRDS(analysis, file = output_file)
            if (verbose) {
                message("[jackknife_entropy_outliers_s4] Analysis object saved to ",
                  output_file)
            }
        }
    }
    analysis
}

# ============================================================================

#' Calculate divergence metrics and store in TSENATAnalysis
#'
#' @param analysis \code{TSENATAnalysis} object.
#' @param q \code{numeric}. Q-value for divergence.
#'   If NULL, uses first q_value from @config$q_values if available, else defaults to 1.0.
#' @param verbose \code{logical}. Print progress messages. Default: TRUE.
#' @param nthreads \code{numeric} or \code{NULL}. Number of CPU threads for parallel processing.
#'   If NULL, reads from \code{@config$nthreads} (or defaults to 1).
#' @param output_file \code{character} or \code{NULL}. Optional file path to save results.
#'   Supported formats: .rds (for S4 objects), .tsv, .csv, .txt (for tables). Default: NULL (no file output).
#' @param control_group \code{character} or \code{NULL}. Control group identifier for divergence comparison.
#'   If NULL, reads from \code{@config$control_group} if available.
#' @param paired \code{logical}. Whether to use paired design. Default: FALSE.
#'   If not specified, reads from \code{@config$paired} if available.
#' @param method \code{character} or \code{NULL}. Statistical method for divergence calculation.
#'   If NULL, reads from \code{@config$method} if available.
#' @param bootstrap \code{logical}. Whether to compute bootstrap confidence intervals. Default: FALSE.
#'   If not specified, reads from \code{@config$bootstrap} if available.
#' @param nboot \code{numeric} or \code{NULL}. Number of bootstrap replicates. Default: NULL.
#'   If NULL, reads from \code{@config$nboot} if available.
#' @param seed \code{numeric} or \code{NULL}. Random seed for reproducibility. Default: NULL.
#'   If NULL, reads from \code{@config$seed} if available.
#' @param progress \code{logical}. Show progress bar during computation. Default: FALSE.
#' @param ... Additional arguments passed to the base divergence function.
#'
#' @return Modified TSENATAnalysis with divergence metrics in @divergence_results
#'   (stored as list of data.frames or matrices).
#'
#' @details
#' Requires diversity results from calculate_diversity_s4() as prerequisite.
#' 
#' Parameters are resolved in priority order:
#' 1. Explicit arguments passed to function
#' 2. Values from analysis@config (if present)
#' 3. Function defaults
#'
#' @examples
#' # Load example data (matching TSENAT.Rmd workflow)
#' data(readcounts)
#' readcounts <- as.matrix(salmon_dataset)
#' mode(readcounts) <- 'numeric'
#' metadata_df <- read.table(
#'   system.file('extdata', 'metadata.tsv', package = 'TSENAT'),
#'   header = TRUE, sep = '\t'
#' )
#' gff3_dataset <- system.file('extdata', 'annotation.gff3.gz', package = 'TSENAT')
#' 
#' # Build analysis from vignette data and create manageable subset
#' analysis <- build_analysis_s4(readcounts = readcounts, tx2gene = gff3_dataset, metadata = metadata_df,
#'   tpm = salmon_tpm, effective_length = salmon_effective_length)
#' # Use 200+ genes to ensure diversity filtering doesn't remove all genes
#' analysis <- filter_analysis_s4(analysis, min_samples = 1, subset_n_genes = 200)
#' 
#' # Compute diversity first (required for divergence)
#' analysis <- calculate_diversity_s4(analysis, q = c(0.5, 1.0, 1.5), verbose = FALSE)
#' 
#' # Calculate divergence across q-values
#' analysis <- calculate_divergence_s4(analysis, verbose = FALSE)
#' 
#' # Check divergence results  
#' head(divergence(analysis))
#'
#' @export
#' @importFrom utils write.table
calculate_divergence_s4 <- function(analysis, q = NULL, verbose = TRUE, nthreads = NULL,
    output_file = NULL, control_group = NULL, paired = FALSE, method = NULL, bootstrap = FALSE,
    nboot = NULL, seed = NULL, progress = FALSE, ...) {
    if (!is(analysis, "TSENATAnalysis")) {
        stop("'analysis' must be a TSENATAnalysis object", call. = FALSE)
    }

    if (nrow(analysis@se) == 0) {
        stop("SummarizedExperiment in @se is empty", call. = FALSE)
    }

    # Check if diversity available
    if (length(analysis@diversity_results) == 0) {
        stop("Diversity results required. Run calculate_diversity_s4() first.", call. = FALSE)
    }

    # PARAMETER EXTRACTION using utility functions
    q <- resolve_slot_param(q, analysis@config, "q_values", 1)
    control_group <- resolve_slot_param(control_group, analysis@config, "control_group",
        NULL)
    method <- resolve_slot_param(method, analysis@config, "method", "percentile")
    nthreads <- resolve_slot_param(nthreads, analysis@config, "nthreads", 1)
    nboot <- resolve_slot_param(nboot, analysis@config, "nboot", NULL)
    seed <- resolve_slot_param(seed, analysis@config, "seed", NULL)

    # Replace q=0 with q=0.01 for practical approximation (q=0 divergence
    # always returns 0, which is mathematically correct but uninformative)
    if (is.vector(q)) {
        q[q == 0] <- 0.01
    } else if (is.numeric(q) && q == 0) {
        q <- 0.01
    }

    # Paired design flag - ensure it's always a logical value
    if (!isTRUE(paired) && !isFALSE(paired)) {
        if ("paired" %in% names(analysis@config)) {
            paired <- analysis@config$paired
        } else {
            paired <- FALSE
        }
        if (!is.logical(paired) || is.na(paired))
            paired <- FALSE
    }

    # Bootstrap parameters - ensure it's always a logical value
    if (!isTRUE(bootstrap) && !isFALSE(bootstrap)) {
        if ("bootstrap" %in% names(analysis@config)) {
            bootstrap <- analysis@config$bootstrap
        } else {
            bootstrap <- FALSE
        }
        if (!is.logical(bootstrap) || is.na(bootstrap))
            bootstrap <- FALSE
    }

    # Sanitize method if needed
    if (!is.null(method)) {
        method <- as.character(method[1])
    }
    if (is.null(method) || is.na(method)) {
        method <- "percentile"
    }

    # Run divergence calculation with extracted parameters Include progress
    # explicitly to prevent NA boolean operations
    args <- list(se = analysis@se, q = q, verbose = verbose, nthreads = nthreads,
        progress = progress  # Use parameter value, not hardcoded
)

    # Add parameters if they are not NULL/FALSE
    if (!is.null(control_group)) {
        args$control_group <- control_group
    }

    # Use isTRUE to safely handle NA values
    if (isTRUE(paired)) {
        args$paired <- paired
    }

    if (!is.null(method)) {
        args$method <- method
    }

    # Use isTRUE to safely handle NA values
    if (isTRUE(bootstrap)) {
        args$bootstrap <- bootstrap
        # Add bootstrap parameters if bootstrap is TRUE
        if (!is.null(nboot)) {
            args$nboot <- nboot
        }
        if (!is.null(seed)) {
            args$seed <- seed
        }
    }

    # Merge with additional args (which may override values)
    args <- c(args, list(...))

    # Run divergence calculation
    result <- tryCatch({
        do.call(.calculate_divergence, args)
    }, error = function(e) {
        # More detailed error handling
        stop("Divergence calculation failed [NA boolean likely in: ", paste0(names(args),
            collapse = ", "), "]:\n", e$message, call. = FALSE)
    })

    # =========================================================================
    # VALIDATE AND STORE RESULTS
    # =========================================================================
    # Store results in list format (required by @divergence_results slot)
    if (is.null(result)) {
        warning("[calculate_divergence_s4] Result is NULL. Check .calculate_divergence() output.",
            call. = FALSE)
        analysis@divergence_results <- list()
    } else if (is(result, "SummarizedExperiment")) {
        # Wrap SummarizedExperiment in list with 'divergence_se' key
        analysis@divergence_results <- list(divergence_se = result)
    } else if (is.list(result)) {
        analysis@divergence_results <- result
    } else if (is.data.frame(result) || is.matrix(result)) {
        # Wrap single result in list
        analysis@divergence_results <- list(main = result)
    } else {
        warning("[calculate_divergence_s4] Unexpected result type: ", class(result),
            ". Wrapping in list.", call. = FALSE)
        analysis@divergence_results <- list(result = result)
    }

    # Track metadata
    analysis@metadata$function_calls <- c(analysis@metadata$function_calls, paste0("calculate_divergence[q=",
        q, "]"))

    # Save if output_file provided
    if (!is.null(output_file)) {
        # Create directory if it doesn't exist
        output_dir <- dirname(output_file)
        if (output_dir != "." && !dir.exists(output_dir)) {
            dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
        }

        if (grepl("\\.tsv$|\\.csv$|\\.txt$", tolower(output_file))) {
            # Write divergence results as text table
            tryCatch({
                if (length(analysis@divergence_results) > 0) {
                  # Extract and combine all divergence results
                  div_list <- analysis@divergence_results

                  # Handle wrapped structure (divergence_se key) or direct
                  # SummarizedExperiments
                  if ("divergence_se" %in% names(div_list)) {
                    write_data <- as.data.frame(SummarizedExperiment::assay(div_list$divergence_se))
                  } else if (is(div_list[[1]], "SummarizedExperiment")) {
                    # Combine multiple SE objects if present
                    data_list <- lapply(div_list, function(se) {
                      as.data.frame(SummarizedExperiment::assay(se, 1))
                    })
                    write_data <- do.call(cbind, data_list)
                  } else {
                    write_data <- as.data.frame(div_list)
                  }

                  # Determine separator
                  sep <- if (grepl("\\.csv$", tolower(output_file)))
                    "," else "\t"
                  write.table(write_data, file = output_file, sep = sep, quote = FALSE,
                    row.names = TRUE)
                } else {
                  warning("[calculate_divergence_s4] No divergence results to save",
                    call. = FALSE)
                }
            }, error = function(e) {
                warning("[calculate_divergence_s4] Could not write divergence results to file: ",
                  conditionMessage(e), call. = FALSE)
            })

            # If bootstrap was used, also save bootstrap results to separate
            # file
            if (isTRUE(bootstrap)) {
                tryCatch({
                  bootstrap_file <- gsub("\\.(tsv|csv|txt)$", "_bootstrap.\\1", tolower(output_file),
                    ignore.case = TRUE)
                  if (basename(bootstrap_file) == basename(output_file)) {
                    # Fallback: if replacement didn't work, append before
                    # extension
                    bootstrap_file <- sub("(\\.[^.]+)$", "_bootstrap\\1", output_file)
                  }

                  if (isTRUE(verbose)) {
                    message("[calculate_divergence_s4] Attempting to save bootstrap results to: ",
                      bootstrap_file)
                  }

                  # Extract bootstrap results from rowData
                  if (length(analysis@divergence_results) > 0) {
                    div_list <- analysis@divergence_results

                    # Handle wrapped structure (divergence_se key)
                    se <- if ("divergence_se" %in% names(div_list)) {
                      div_list$divergence_se
                    } else if (is(div_list[[1]], "SummarizedExperiment")) {
                      div_list[[1]]
                    } else {
                      NULL
                    }

                    if (!is.null(se) && is(se, "SummarizedExperiment")) {
                      rd <- SummarizedExperiment::rowData(se)

                      # Extract bootstrap-related columns: gene_name + all
                      # CI/estimate/nboot columns
                      bootstrap_cols <- c("gene_name", grep("^estimate_q|^lower_ci_q|^upper_ci_q|^ci_width_q|^nboot_q|^method_q",
                        colnames(rd), value = TRUE))

                      # Include computation_time_sec and error if present
                      if ("computation_time_sec" %in% colnames(rd)) {
                        bootstrap_cols <- c(bootstrap_cols, "computation_time_sec")
                      }
                      if ("error" %in% colnames(rd)) {
                        bootstrap_cols <- c(bootstrap_cols, "error")
                      }

                      # Filter to only existing columns
                      bootstrap_cols <- bootstrap_cols[bootstrap_cols %in% colnames(rd)]

                      if (length(bootstrap_cols) > 1) {
                        # More than just gene_name
                        bootstrap_data <- as.data.frame(rd[, bootstrap_cols])
                        sep <- if (grepl("\\.csv$", tolower(bootstrap_file)))
                          "," else "\t"
                        write.table(bootstrap_data, file = bootstrap_file, sep = sep,
                          quote = FALSE, row.names = FALSE)

                        if (isTRUE(verbose)) {
                          message("[calculate_divergence_s4] Saved bootstrap results to: ",
                            bootstrap_file)
                        }
                      }
                    }
                  }
                }, error = function(e) {
                  warning("[calculate_divergence_s4] Could not write bootstrap results to file: ",
                    conditionMessage(e), call. = FALSE)
                })
            }
        } else {
            # Default to RDS for S4 object
            saveRDS(analysis, file = output_file)
        }

    }

    analysis
}



# ============================================================================
# CALCULATE DIFFERENCE WRAPPER
# ============================================================================

#' Calculate Difference Between Control and Treatment Groups (S4 Wrapper)
#'
#' S4 wrapper that operates on TSENATAnalysis objects to calculate differences
#' between control and treatment groups. Uses diversity results from \code{@diversity_results} 
#' slot (from \code{calculate_diversity_s4()}) and stores results in the \code{lm_results} slot.
#'
#' @param analysis A \code{TSENATAnalysis} object with diversity results in \code{@diversity_results}.
#' @param q \code{numeric}. Q-value to use. If NULL, uses first diversity result or q=1.0.
#' @param control Character string specifying the control group identifier. If \code{NULL},
#'   attempts to retrieve from \code{analysis@config$control}.
#' @param condition_col \code{character} or \code{NULL}. Column name in colData identifying sample conditions.
#'   If NULL, reads from \code{@config$condition_col} or auto-detects.
#' @param method \code{character}. Difference calculation method. Default: 'mean'.
#'   If NULL, reads from \code{@config$method} if available.
#' @param test \code{character}. Statistical test type. Default: 'wilcoxon'.
#'   If NULL, reads from \code{@config$test} if available.
#' @param randomizations \code{numeric}. Number of randomizations. Default: 100.
#'   If NULL, reads from \code{@config$randomizations} if available.
#' @param pcorr \code{character}. P-value correction method. Default: 'BH'.
#'   If NULL, reads from \code{@config$pcorr} if available.
#' @param assayno \code{numeric}. Assay number to use. Default: 1.
#'   If NULL, reads from \code{@config$assayno} if available.
#' @param verbose \code{logical}. Print progress messages. Default: TRUE.
#'   If not specified, reads from \code{@config$verbose} if available.
#' @param paired \code{logical}. Whether data is paired. Default: FALSE.
#'   If not specified, reads from \code{@config$paired} if available.
#' @param exact \code{logical}. Use exact test. Default: FALSE.
#'   If not specified, reads from \code{@config$exact} if available.
#' @param pseudocount \code{numeric}. Pseudocount for normalization. Default: 0.
#'   If NULL, reads from \code{@config$pseudocount} if available.
#' @param nthreads \code{numeric} or \code{NULL}. Number of CPU threads for parallel processing.
#'   If NULL, reads from \code{@config$nthreads} (or defaults to 1).
#' @param seed \code{numeric} or \code{NULL}. Random seed. Default: NULL.
#'   If NULL, reads from \code{@config$seed} if available.
#' @param robust_loss_type \code{character}. Robust regression loss type. Default: 'huber'.
#'   If NULL, reads from \code{@config$robust_loss_type} if available.
#' @param robust_scale_method \code{character}. Robust scaling method. Default: 'mad'.
#'   If NULL, reads from \code{@config$robust_scale_method} if available.
#' @param output_file \code{character} or \code{NULL}. Optional file path to save results.
#'   Supported formats: .rds (for S4 objects), .tsv, .csv, .txt (for tables). Default: NULL (no file output).
#' @param ... Additional arguments passed to the base function.
#'
#' @return Returns the modified \code{analysis} object invisibly with results stored in
#'   \code{analysis@lm_results$difference}.
#'
#' @details
#' **IMPORTANT:** Requires diversity results to exist first via \code{calculate_diversity_s4()}.
#' This wrapper extracts the diversity SummarizedExperiment from \code{@diversity_results},
#' not the raw input data in \code{@se}. This ensures you're comparing diversity values
#' between control and treatment groups, not raw abundance data.
#'
#' **Parameter resolution priority** (explicit > @config > auto-detect > error):
#' \itemize{
#'   \item \code{control}: Uses explicit arg, else \code{@config$control}, else error
#'   \item \code{nthreads}: Uses explicit arg, else \code{@config$nthreads}, else 1
#'   \item \code{condition_col} (sample grouping): Uses \code{@config$condition_col},
#'     else auto-detects from colData columns: 'group', 'sample_type', 'condition'
#' }
#'
#' @importFrom utils write.table
#' @seealso
#' \code{\link{calculate_diversity_s4}} for computing diversity.
#'
#' @examples
#' # Load example data (matching TSENAT.Rmd workflow)
#' data(readcounts)
#' readcounts <- as.matrix(salmon_dataset)
#' mode(readcounts) <- 'numeric'
#' metadata_df <- read.table(
#'   system.file('extdata', 'metadata.tsv', package = 'TSENAT'),
#'   header = TRUE, sep = '\t'
#' )
#' gff3_dataset <- system.file('extdata', 'annotation.gff3.gz', package = 'TSENAT')
#' 
#' # Build analysis from vignette data and create small subset
#' analysis <- build_analysis_s4(readcounts = readcounts, tx2gene = gff3_dataset, metadata = metadata_df,
#'   tpm = salmon_tpm, effective_length = salmon_effective_length)
#' analysis <- filter_analysis_s4(analysis, min_samples = 1, subset_n_genes = 200)
#' analysis <- calculate_diversity_s4(analysis, q = c(0.5, 1.0, 1.5), verbose = FALSE)
#' result <- calculate_difference_s4(analysis, control = 'normal', verbose = FALSE)
#'
#' @export
calculate_difference_s4 <- function(analysis, control = NULL, q = NULL, condition_col = NULL,
    method = NULL, test = NULL, randomizations = NULL, pcorr = NULL, assayno = NULL,
    verbose = NULL, paired = FALSE, exact = FALSE, pseudocount = NULL, nthreads = NULL,
    seed = NULL, robust_loss_type = NULL, robust_scale_method = NULL, output_file = NULL,
    ...) {
    if (!is(analysis, "TSENATAnalysis")) {
        stop("'analysis' must be a TSENATAnalysis object", call. = FALSE)
    }

    # Check prerequisites: diversity results must exist
    if (length(analysis@diversity_results) == 0) {
        stop("Diversity results required. Run calculate_diversity_s4() first.", call. = FALSE)
    }

    # Priority 1: Use explicit parameter Priority 2: Use @config$control
    if (is.null(control)) {
        if ("control" %in% names(analysis@config)) {
            control <- analysis@config$control
        } else {
            stop("'control' must be specified (control group identifier) or set in config",
                call. = FALSE)
        }
    }

    # Determine which diversity result to use
    if (is.null(q)) {
        div_keys <- names(analysis@diversity_results)
        if (length(div_keys) == 0) {
            stop("No diversity results found in @diversity_results", call. = FALSE)
        }
        diversity_se <- analysis@diversity_results[[div_keys[1]]]
        q_used <- sub("^q_", "", div_keys[1])
    } else {
        q_key <- paste0("q_", formatC(q, format = "f", digits = 3))
        if (!(q_key %in% names(analysis@diversity_results))) {
            stop("Diversity not calculated for q=", q, ". Available: ", paste(names(analysis@diversity_results),
                collapse = ", "), call. = FALSE)
        }
        diversity_se <- analysis@diversity_results[[q_key]]
        q_used <- q
    }

    # Determine condition column to use
    condition_col <- resolve_slot_param(condition_col, analysis@config, "condition_col",
        NULL)

    # Resolve remaining parameters using centralized handler
    method <- resolve_slot_param(method, analysis@config, "method", "mean")
    test <- resolve_slot_param(test, analysis@config, "test", "wilcoxon")
    randomizations <- resolve_slot_param(randomizations, analysis@config, "randomizations",
        100)
    pcorr <- resolve_slot_param(pcorr, analysis@config, "pcorr", "BH")
    assayno <- resolve_slot_param(assayno, analysis@config, "assayno", 1)
    verbose <- resolve_slot_param(verbose, analysis@config, "verbose", TRUE)
    pseudocount <- resolve_slot_param(pseudocount, analysis@config, "pseudocount",
        0)
    nthreads <- resolve_slot_param(nthreads, analysis@config, "nthreads", 1)
    robust_loss_type <- resolve_slot_param(robust_loss_type, analysis@config, "robust_loss_type",
        "huber")
    robust_scale_method <- resolve_slot_param(robust_scale_method, analysis@config,
        "robust_scale_method", "mad")

    # Parameters with logical defaults (check config if FALSE)
    if (!paired && "paired" %in% names(analysis@config)) {
        paired <- analysis@config$paired
    }
    if (!exact && "exact" %in% names(analysis@config)) {
        exact <- analysis@config$exact
    }

    # Optional parameters (may be NULL)
    seed <- resolve_slot_param(seed, analysis@config, "seed", NULL)

    # Run difference calculation on diversity results Note: diversity_se and
    # its colData are already prepared by calculate_diversity_s4
    result <- tryCatch({
        .calculate_difference(x = diversity_se, condition_col = condition_col, control = control,
            method = method, test = test, randomizations = randomizations, pcorr = pcorr,
            assayno = assayno, verbose = verbose, paired = paired, exact = exact,
            pseudocount = pseudocount, nthreads = nthreads, seed = seed, robust_loss_type = robust_loss_type,
            robust_scale_method = robust_scale_method, ...)
    }, error = function(e) {
        stop("Difference calculation failed:\n", e$message, call. = FALSE)
    })

    # Store in lm_results under 'difference' key
    if (is.list(analysis@lm_results)) {
        analysis@lm_results$difference <- result
    } else {
        analysis@lm_results <- list(difference = result)
    }

    # Track metadata
    analysis@metadata$function_calls <- c(analysis@metadata$function_calls, paste0("calculate_difference_s4[q=",
        q_used, ", control=", control, "]"))

    # Save if output_file provided (using centralized output handler)
    if (!is.null(output_file)) {
        diff_data <- if (!is.null(analysis@lm_results$difference$results)) {
            analysis@lm_results$difference$results
        } else {
            as.data.frame(analysis@lm_results$difference)
        }
        save_analysis_output(diff_data, output_file, object = analysis, verbose = verbose,
            func_name = "calculate_difference_s4")
    }

    analysis
}

# ============================================================================
# TEST RANKBASED ASSUMPTIONS WRAPPER
# ============================================================================

#' Test rank-based method assumptions in TSENATAnalysis
#'
#' @param analysis \code{TSENATAnalysis} object with diversity results stored
#'   in \code{@diversity_results}.
#' @param q \code{numeric}. Q-value(s) to extract from diversity results.
#'   If NULL, uses the first available diversity result or q=1.0.
#' @param checks \code{character}. Which assumptions to test. Default includes:
#'   'exchangeability', 'monotonicity', 'consistency'.
#' @param alpha \code{numeric}. Significance level for tests (default: 0.05).
#' @param ... Additional arguments (for future extensibility).
#'
#' @return Modified TSENATAnalysis object with assumption test results stored
#'   in \code{@metadata$rankbased_assumptions}.
#'
#' @details
#' This wrapper calls \code{.test_rankbased_assumptions()} on diversity data
#' extracted from the analysis object. Results include:
#'
#' \describe{
#'   \item{exchangeability}{Permutation test for temporal/spatial ordering effects}
#'   \item{monotonicity}{Spearman correlation stability across rows}
#'   \item{consistency}{Kendall's W concordance and ICC across samples}
#' }
#'
#' **Data Extraction Priority:**
#' 1. If q specified: uses diversity result for that q-value
#' 2. If q NULL: uses first available diversity result
#' 3. If no diversity results: extracts from cached combined result (\code{@metadata$diversity_combined})
#'
#' @examples
#' # Load example data (matching TSENAT.Rmd workflow)
#' data(readcounts)
#' readcounts <- as.matrix(salmon_dataset)
#' mode(readcounts) <- 'numeric'
#' metadata_df <- read.table(
#'   system.file('extdata', 'metadata.tsv', package = 'TSENAT'),
#'   header = TRUE, sep = '\t'
#' )
#' gff3_dataset <- system.file('extdata', 'annotation.gff3.gz', package = 'TSENAT')
#' 
#' # Build analysis from vignette data and create small subset
#' analysis <- build_analysis_s4(readcounts = readcounts, tx2gene = gff3_dataset, metadata = metadata_df,
#'   tpm = salmon_tpm, effective_length = salmon_effective_length)
#' analysis <- filter_analysis_s4(analysis, min_samples = 1, subset_n_genes = 200)
#' analysis <- calculate_diversity_s4(analysis, q = c(0.5, 1.0, 1.5), verbose = FALSE)
#' analysis <- test_rankbased_assumptions_s4(analysis, q = 1.0)
#' # Access results using getMeta S4 accessor
#' names(getMeta(analysis, 'rankbased_assumptions'))
#'
#' @export
#' @rdname test_rankbased_assumptions_s4
setGeneric("test_rankbased_assumptions_s4", function(analysis, q = NULL, checks = c("exchangeability",
    "monotonicity", "consistency"), alpha = 0.05, ...) {
    standardGeneric("test_rankbased_assumptions_s4")
})

#' @rdname test_rankbased_assumptions_s4
setMethod("test_rankbased_assumptions_s4", signature(analysis = "TSENATAnalysis"),
    function(analysis, q = NULL, checks = c("exchangeability", "monotonicity", "consistency"),
        alpha = 0.05, ...) {

        # Validate inputs
        if (!methods::is(analysis, "TSENATAnalysis")) {
            stop("'analysis' must be a TSENATAnalysis object", call. = FALSE)
        }

        # Extract diversity data
        diversity_data <- NULL
        q_used <- q

        # If q is specified, try to get that specific q-value
        if (!is.null(q)) {
            q_key <- if (nchar(as.character(q)) > 3) {
                paste0("q_", round(q, 1))
            } else {
                paste0("q_", q)
            }

            if (q_key %in% names(analysis@diversity_results)) {
                div_se <- analysis@diversity_results[[q_key]]
                diversity_data <- assay(div_se, "diversity")
            }
        }

        # If q is NULL and multiple diversity results exist, combine all
        # q-values
        if (is.null(diversity_data) && is.null(q) && length(analysis@diversity_results) >
            1) {
            entropy_list <- lapply(analysis@diversity_results, function(se) {
                mat <- assay(se, "diversity")
                if (!is.matrix(mat)) {
                  mat <- as.matrix(mat)
                }
                return(mat)
            })

            # Use complete case analysis: keep only genes present in ALL
            # q-value matrices This is mathematically sound for rank-based
            # tests (Friedman) and follows best practices per scholarly
            # literature (Springer Handbook, Permutation Tests)
            all_genes <- lapply(entropy_list, rownames)
            common_genes <- Reduce(intersect, all_genes)

            # Subset all matrices to common genes in same order
            entropy_list <- lapply(entropy_list, function(mat) {
                mat[common_genes, , drop = FALSE]
            })

            # Combine all matrices column-wise (genes x all samples across
            # q-values)
            diversity_data <- do.call(cbind, entropy_list)
            # Keep natural column names from cbind to preserve structure
            q_used <- "all"
        }

        # If q is NULL and only one result, use it
        if (is.null(diversity_data) && is.null(q) && length(analysis@diversity_results) ==
            1) {
            div_se <- analysis@diversity_results[[1]]
            diversity_data <- assay(div_se, "diversity")
            q_used <- .extract_q_from_key(names(analysis@diversity_results)[1])
        }

        # Fallback: use first diversity result
        if (is.null(diversity_data) && length(analysis@diversity_results) > 0) {
            div_se <- analysis@diversity_results[[1]]
            diversity_data <- assay(div_se, "diversity")
            if (is.null(q_used) || is.na(q_used)) {
                q_used <- .extract_q_from_key(names(analysis@diversity_results)[1])
            }
        }

        # Check that we have data
        if (is.null(diversity_data)) {
            stop("No diversity results found in analysis object. ", "Run calculate_diversity_s4() first.",
                call. = FALSE)
        }

        # Ensure we have a matrix
        if (!is.matrix(diversity_data)) {
            diversity_data <- as.matrix(diversity_data)
        }

        # Run assumptions test
        result <- tryCatch({
            .test_rankbased_assumptions(data = diversity_data, checks = checks, alpha = alpha)
        }, error = function(e) {
            stop("Rankbased assumptions test failed:\n", e$message, call. = FALSE)
        })

        # Store results
        analysis@metadata$rankbased_assumptions <- list(result = result, q_value_tested = q_used,
            checks_performed = checks, alpha_used = alpha, timestamp = Sys.time())

        # Track function call
        analysis@metadata$function_calls <- c(analysis@metadata$function_calls, paste0("test_rankbased_assumptions_s4[q=",
            q_used, "]"))

        analysis
    })

# Helper function to extract q-value from key

.extract_q_from_key <- function(key) {
    # Extract numeric part from 'q_X.X' format
    as.numeric(sub("^q_", "", key))
}

#' Plot Volcano and MA Grid from Differential Analysis Results (S4 Wrapper)
#'
#' S4 wrapper that extracts differential analysis results from a TSENATAnalysis
#' object and creates side-by-side volcano and MA plots for comparing control
#' and treatment groups.
#'
#' @param analysis \code{TSENATAnalysis} object with calculated differences
#'   (typically via \code{\link{calculate_difference_s4}}).
#' @param x_col \code{character}. Column name for x-axis in MA plot.
#'   Default: NULL (uses mean_difference if available, else mean fold-change).
#' @param padj_col \code{character}. Column name for adjusted p-values.
#'   Default: 'padj' (the standard column name from calculate_difference).
#' @param label_thresh \code{numeric}. P-value threshold for labeling top genes.
#'   Genes with adjusted p-value below this threshold are labeled.
#'   Default: 0.1.
#' @param sig_alpha \code{numeric}. Significance threshold for coloring significant
#'   differences. Points with adjusted p-value below sig_alpha are highlighted.
#'   Default: 0.05.
#' @param top_n \code{integer}. Number of top genes (by significance) to label
#'   in volcano plot. Default: 5.
#' @param title_volcano \code{character}. Title for volcano plot.
#'   Default: NULL (no title).
#' @param title_ma \code{character}. Title for MA plot.
#'   Default: 'Tsallis-based MA plot'.
#' @param verbose \code{logical}. Print status messages. Default: FALSE.
#' @param output_file \code{character} or \code{NULL}. Optional file path to save the plot.
#'   Default: NULL (no file output).
#' @param width \code{numeric}. Width of the output plot in inches (default: 12).
#'   Only used if output_file is provided.
#' @param height \code{numeric}. Height of the output plot in inches (default: 7.2).
#'   Only used if output_file is provided.
#' @param ... Additional arguments passed to the base plotting function.
#'
#' @return
#' Invisibly returns a cowplot grid object containing both volcano and MA plots
#' combined side-by-side. If the plot cannot be created, returns NULL invisibly.
#'
#' @details
#' This wrapper extracts the difference results data frame from
#' \code{analysis@lm_results$difference} and passes it to the base
#' \code{.plot_volcano_ma_grid()} function.
#'
#' **Required Data:**
#' \itemize{
#'   \item Differential analysis must be computed via \code{calculate_difference_s4()}
#'   \item Results are stored in \code{analysis@lm_results$difference}
#' }
#'
#' **Expected Columns in Difference Results:**
#' \itemize{
#'   \item \code{genes} or \code{gene_id}: Gene identifiers
#'   \item \code{Normal_mean}, \code{Tumor_mean}: Group means (or equivalent controls/treatments)
#'   \item \code{mean_difference}: Calculated difference between groups
#'   \item \code{log2_fold_change}: Log2 fold-change values
#'   \item \code{raw_p_values} or \code{pvalue}: Un-adjusted p-values
#'   \item \code{adjusted_p_values} or \code{padj}: Adjusted p-values (default column used)
#' }
#'
#' **Volcano Plot Features:**
#' \itemize{
#'   \item X-axis: log2 fold-change or mean difference
#'   \item Y-axis: -log10(adjusted p-value)
#'   \item Top significant genes labeled
#'   \item Points colored by significance threshold
#' }
#'
#' **MA Plot Features:**
#' \itemize{
#'   \item X-axis: Average expression level (A)
#'   \item Y-axis: Log2 fold-change (M)
#'   \item Loess curve showing trend
#'   \item Significant changes highlighted
#' }
#'
#' @examples
#' # Load example data (matching TSENAT.Rmd workflow)
#' data(readcounts)
#' readcounts <- as.matrix(salmon_dataset)
#' mode(readcounts) <- 'numeric'
#' metadata_df <- read.table(
#'   system.file('extdata', 'metadata.tsv', package = 'TSENAT'),
#'   header = TRUE, sep = '\t'
#' )
#' gff3_dataset <- system.file('extdata', 'annotation.gff3.gz', package = 'TSENAT')
#' 
#' # Build analysis from vignette data and create small subset
#' analysis <- build_analysis_s4(readcounts = readcounts, tx2gene = gff3_dataset, metadata = metadata_df,
#'   tpm = salmon_tpm, effective_length = salmon_effective_length)
#' analysis <- filter_analysis_s4(analysis, min_samples = 1, subset_n_genes = 200)
#' analysis <- calculate_diversity_s4(analysis, q = c(0.5, 1.0, 1.5), verbose = FALSE)
#' analysis <- calculate_difference_s4(analysis, control = 'normal', verbose = FALSE)
#'   
#' # Plot volcano and MA plots
#' p <- plot_volcano_ma_grid_s4(analysis, sig_alpha = 0.05, top_n = 3)
#' print(p)
#'
#' @seealso
#' \code{\link{calculate_difference_s4}} for computing differential analysis.
#'
#' @export
plot_volcano_ma_grid_s4 <- function(analysis, x_col = NULL, padj_col = "padj", label_thresh = 0.1,
    sig_alpha = 0.05, top_n = 5, title_volcano = NULL, title_ma = "Tsallis-based MA plot",
    verbose = FALSE, output_file = NULL, width = 12, height = 7.2, ...) {

    # Load visualization dependencies (ggplot2, cowplot, etc.)
    .load_visualization_deps()

    # Extract verbose parameter if not provided
    verbose <- resolve_slot_param(verbose, analysis@config, "verbose", FALSE)

    # Validate input
    if (!is(analysis, "TSENATAnalysis")) {
        stop("'analysis' must be a TSENATAnalysis object", call. = FALSE)
    }

    # Extract difference results from S4 object
    if (is.null(analysis@lm_results) || !is.list(analysis@lm_results)) {
        stop("No LM results found in analysis@lm_results. ", "Run calculate_difference_s4() first.",
            call. = FALSE)
    }

    if (!("difference" %in% names(analysis@lm_results))) {
        stop("Difference results not found in analysis@lm_results$difference. ",
            "Run calculate_difference_s4() first.", call. = FALSE)
    }

    diff_df <- analysis@lm_results$difference

    if (!is.data.frame(diff_df) || nrow(diff_df) == 0) {
        stop("Difference results are empty or not a data frame", call. = FALSE)
    }

    # Auto-detect column names - padj_col with fallbacks, x_col for effect size
    actual_padj_col <- auto_detect_column(colnames(diff_df), analysis@config, "padj_col",
        c("padj", "adjusted_p_values", "pvalue"), verbose = verbose, param_name = "padj_col")

    if (is.null(x_col)) {
        x_col <- auto_detect_column(colnames(diff_df), analysis@config, "x_col",
            c("mean_difference", "log2_fold_change", "effect_size"), verbose = verbose,
            param_name = "x_col")
    }

    plot_obj <- tryCatch({
        .plot_volcano_ma_grid(diff_df = diff_df, x_col = x_col, padj_col = actual_padj_col,
            label_thresh = label_thresh, sig_alpha = sig_alpha, top_n = top_n, title_volcano = title_volcano,
            title_ma = title_ma, ...)
    }, error = function(e) {
        stop("[plot_volcano_ma_grid_s4]", conditionMessage(e), call. = FALSE)
    })

    if (verbose) {
        message("[plot_volcano_ma_grid_s4] Plot created successfully")
    }

    # Save plot to file if requested
    if (!is.null(output_file)) {
        save_analysis_output(plot_obj, output_file, object = analysis, verbose = verbose,
            func_name = "plot_volcano_ma_grid_s4", width = width, height = height)
    }

    return(invisible(plot_obj))
}

# ============================================================================
# CONCORDANCE WRAPPER - Compute Method Concordance (GAM vs Friedman/KW)
# ============================================================================

#' Compute concordance between two analysis methods in TSENATAnalysis
#'
#' @param analysis \code{TSENATAnalysis} object with LM results (e.g., GAM).
#' @param gam_method \code{character}. Key for GAM/interaction results in \code{@lm_results}.
#'   Default: 'q_interactions' (results from \code{rank_test_q_condition_s4})
#' @param friedman_method \code{character}. Key for Friedman/rank-based results in \code{@lm_results}.
#'   Default: 'rankbased' (results from \code{test_rankbased_assumptions_s4})
#' @param gam_results \code{data.frame} or \code{NULL}. Optional GAM results data frame to store
#'   in the analysis object. If provided, automatically stored in \code{@lm_results} under the
#'   key specified by \code{gam_method}. Useful for importing external results or results 
#'   computed outside the S4 wrapper. Default: NULL (use existing results in analysis).
#' @param verbose \code{logical}. Print progress messages (default: FALSE).
#' @param output_file \code{character} or \code{NULL}. Optional file path to save results.
#'   Supported formats: .rds (for S4 objects). Default: NULL (no file output).
#' @param ... Additional arguments for future extensibility.
#'
#' @return Modified TSENATAnalysis object with concordance results stored in:
#'   \code{@metadata$method_concordance}:
#'   \describe{
#'     \item{comparison_df}{Data frame comparing results from both methods}
#'     \item{spearman_rho}{Spearman correlation between adjusted p-values}
#'     \item{high_confidence}{Genes with strong agreement}
#'     \item{agreement_table}{Contingency table of significant/non-significant calls}
#'     \item{gam_method}{Method name used for GAM analysis}
#'     \item{friedman_method}{Method name used for Friedman analysis}
#'     \item{timestamp}{When concordance was computed}
#'   }
#'
#' @details
#' Compares results from two different statistical methods (typically GAM for continuous 
#' and Friedman/Kruskal-Wallis for rank-based analysis) on the same data. Identifies:
#' - Genes significant in both methods (high confidence)
#' - Genes detected by one method only (potential false positives or method-specific signal)
#' - Spearman correlation of p-values (overall agreement trends)
#'
#' @examples
#' # Load example data (matching TSENAT.Rmd workflow)
#' data(readcounts)
#' readcounts <- as.matrix(salmon_dataset)
#' mode(readcounts) <- 'numeric'
#' metadata_df <- read.table(
#'   system.file('extdata', 'metadata.tsv', package = 'TSENAT'),
#'   header = TRUE, sep = '\t'
#' )
#' gff3_dataset <- system.file('extdata', 'annotation.gff3.gz', package = 'TSENAT')
#' 
#' # Build analysis from vignette data (matching TSENAT.Rmd workflow)
#' analysis <- build_analysis_s4(readcounts = readcounts, tx2gene = gff3_dataset, metadata = metadata_df,
#'   tpm = salmon_tpm, effective_length = salmon_effective_length)
#' 
#' # Configure analysis parameters (best practice for reproducibility)
#' config <- tsenat_config(
#'   condition_col = 'condition',
#'   subject_col = 'paired_samples',
#'   paired = TRUE,
#'   control = 'normal',
#'   metadata = metadata_df
#' )
#' analysis <- setConfig(analysis, config)
#' 
#' analysis <- filter_analysis_s4(analysis, stringency = 'severe')
#' analysis <- calculate_diversity_s4(analysis, q = c(0.5, 1.0, 1.5), verbose = FALSE)
#' analysis <- calculate_divergence_s4(analysis, q = c(0.5, 1.0, 1.5), verbose = FALSE)
#' analysis <- calculate_lm_interaction_s4(analysis, method = 'gam', verbose = FALSE)
#' # Note: compute_method_concordance_s4 requires results from both
#' # rank_test_q_condition_s4 and test_rankbased_assumptions_s4
#'
#' @aliases compute_method_concordance_s4
#' @export
setGeneric("compute_method_concordance_s4", function(analysis, ...) {
    standardGeneric("compute_method_concordance_s4")
})

#' @rdname compute_method_concordance_s4
setMethod("compute_method_concordance_s4", "TSENATAnalysis", function(analysis, gam_method = "q_interactions",
    friedman_method = "rankbased", gam_results = NULL, verbose = FALSE, output_file = NULL) {

    # ===================================================================
    # VALIDATION
    # ===================================================================

    if (!is(analysis, "TSENATAnalysis")) {
        stop("'analysis' must be a TSENATAnalysis object", call. = FALSE)
    }

    # If gam_results provided, store them in lmResults automatically
    if (!is.null(gam_results)) {
        if (!is.data.frame(gam_results)) {
            stop("gam_results must be a data.frame", call. = FALSE)
        }
        # Store GAM results in analysis@lm_results
        if (is.null(analysis@lm_results)) {
            analysis@lm_results <- list()
        }
        analysis@lm_results[[gam_method]] <- gam_results
        if (verbose) {
            message("[compute_method_concordance_s4] Stored GAM results as '", gam_method,
                "'")
        }
    }

    if (is.null(analysis@lm_results)) {
        stop("No LM results found in analysis@lm_results. Run rank_test_q_condition_s4() first.",
            call. = FALSE)
    }

    # Check for required methods
    if (!(gam_method %in% names(analysis@lm_results))) {
        available_methods <- paste(names(analysis@lm_results), collapse = ", ")
        stop("GAM method '", gam_method, "' not found in LM results. ", "Available: ",
            available_methods, call. = FALSE)
    }

    if (!(friedman_method %in% names(analysis@lm_results))) {
        available_methods <- paste(names(analysis@lm_results), collapse = ", ")
        stop("Friedman method '", friedman_method, "' not found in LM results. ",
            "Available: ", available_methods, call. = FALSE)
    }

    # Extract results
    gam_results_final <- analysis@lm_results[[gam_method]]
    friedman_results <- analysis@lm_results[[friedman_method]]

    # Validate they're data frames
    if (!is.data.frame(gam_results_final)) {
        stop("GAM results ('", gam_method, "') must be a data.frame", call. = FALSE)
    }

    if (!is.data.frame(friedman_results)) {
        stop("Friedman results ('", friedman_method, "') must be a data.frame", call. = FALSE)
    }

    # ===================================================================
    # COMPUTE CONCORDANCE
    # ===================================================================

    if (verbose) {
        message("[compute_method_concordance_s4] Computing concordance between ",
            gam_method, " and ", friedman_method)
    }

    # Call the standard function
    concordance_result <- tryCatch({
        .compute_method_concordance(gam_results_final, friedman_results)
    }, error = function(e) {
        stop("[compute_method_concordance_s4]", conditionMessage(e), call. = FALSE)
    })

    # =================================================================== STORE
    # RESULTS
    # ===================================================================

    analysis@metadata$method_concordance <- list(comparison_df = concordance_result$comparison_df,
        spearman_rho = concordance_result$spearman_rho, high_confidence = concordance_result$high_conf,
        agreement_table = concordance_result$agreement_table, gam_method = gam_method,
        friedman_method = friedman_method, timestamp = Sys.time())

    # Track function call
    analysis@metadata$function_calls <- c(analysis@metadata$function_calls, paste0("compute_method_concordance_s4[",
        gam_method, " vs ", friedman_method, "]"))

    if (verbose) {
        message("[compute_method_concordance_s4] Concordance computed successfully")
        if (!is.na(concordance_result$spearman_rho)) {
            message("[compute_method_concordance_s4] Spearman corr = ", round(concordance_result$spearman_rho,
                3))
        }
    }

    # =================================================================== SAVE
    # TO FILE (if output_file provided)
    # ===================================================================

    if (!is.null(output_file)) {
        if (verbose) {
            message("[compute_method_concordance_s4] Writing results to: ", output_file)
        }
        saveRDS(analysis, file = output_file)
    }

    analysis
})

#' Plot Global Divergence q-Curve Across All Genes (S4 Wrapper)
#'
#' S4 wrapper that extracts divergence results from a TSENATAnalysis object
#' and visualizes the average Tsallis divergence across all genes (or specified
#' genes) as a function of q-value.
#'
#' @param analysis \code{TSENATAnalysis} object with divergence results
#'   (typically via \code{\link{calculate_divergence_s4}}).
#' @param gene \code{character}. Optional specific gene name to plot.
#'   If NULL, plots global divergence curve (aggregated across all genes).
#' @param n_genes \code{integer}. Number of top genes to plot when showing
#'   multi-gene spectra. Default is 4. Genes are sorted by p-value significance.
#' @param ncol \code{integer}. Number of columns in grid layout for multi-gene
#'   plots. Default is 2. Number of rows is automatically calculated.
#' @param metric \code{character}. Summary statistic for global curve:
#'   'median' (default) or 'mean'. Only used when gene = NULL.
#' @param variability_metric \code{character}. Error bar type for global curve:
#'   'iqr' (default) or 'sd'. Only used when gene = NULL.
#' @param use_pvalue_ranking \code{logical}. If TRUE, uses LM results to rank and
#'   display top n_genes by p-value significance. If FALSE (default), plots global
#'   divergence curve when gene = NULL. Default is FALSE.
#' @param output_file \code{character}. Optional file path to save the plot.
#'   If NULL, plot is returned but not saved.
#' @param width \code{numeric}. Plot width in inches. Default is 10.
#' @param height \code{numeric}. Plot height in inches. Default is 6.
#' @param verbose \code{logical}. Print status messages. Default is TRUE.
#' @param ... Additional arguments passed to the underlying plotting function.
#'
#' @return
#' Invisibly returns the file path if saved, otherwise the ggplot object.
#' If the plot cannot be created (missing data, ggplot2 not available),
#' returns NULL invisibly with an informative message.
#'
#' @details
#' This wrapper extracts the divergence SummarizedExperiment from
#' \code{analysis@divergence_results} and optionally the LM results from
#' \code{analysis@lm_results$lm_interaction} to pass to the base function.
#'
#' **Data Requirements:**
#' \itemize{
#'   \item Divergence must be computed via \code{calculate_divergence_s4()}
#'   \item \code{@divergence_results$divergence_se} or direct divergence SE
#' }
#'
#' **Modes:**
#' \itemize{
#'   \item \strong{Global mode} (gene = NULL): Shows median/mean divergence
#'         across all genes with variability bands
#'   \item \strong{Gene-specific mode} (gene specified): Shows divergence
#'         spectrum for a single named gene
#'   \item \strong{Top genes mode} (gene = NULL, lm_res provided): Shows
#'         top n_genes by significance
#' }
#'
#' @examples
#' # Load example data (matching TSENAT.Rmd workflow)
#' data(readcounts)
#' readcounts <- as.matrix(salmon_dataset)
#' mode(readcounts) <- 'numeric'
#' metadata_df <- read.table(
#'   system.file('extdata', 'metadata.tsv', package = 'TSENAT'),
#'   header = TRUE, sep = '\t'
#' )
#' gff3_dataset <- system.file('extdata', 'annotation.gff3.gz', package = 'TSENAT')
#' 
#' # Build analysis from vignette data and create small subset
#' analysis <- build_analysis_s4(readcounts = readcounts, tx2gene = gff3_dataset, metadata = metadata_df,
#'   tpm = salmon_tpm, effective_length = salmon_effective_length)
#' analysis <- filter_analysis_s4(analysis, min_samples = 1, subset_n_genes = 200)
#' analysis <- calculate_diversity_s4(analysis, q = c(0.5, 1, 1.5), verbose = FALSE)
#' analysis <- calculate_divergence_s4(analysis, q = c(0.5, 1, 1.5), verbose = FALSE)
#' p_global <- plot_divergence_spectrum_s4(analysis)
#' print(p_global)
#'
#' @seealso
#' \code{\link{calculate_divergence_s4}} for computing divergence.
#'
#' @export
plot_divergence_spectrum_s4 <- function(analysis, gene = NULL, n_genes = 4, ncol = 2,
    metric = c("median", "mean"), variability_metric = c("iqr", "sd"), use_pvalue_ranking = FALSE,
    output_file = NULL, width = 12, height = NULL, verbose = TRUE, ...) {

    # Load visualization dependencies (ggplot2, cowplot, etc.)
    .load_visualization_deps()

    # Extract verbose parameter if not provided
    verbose <- resolve_slot_param(verbose, analysis@config, "verbose", TRUE)

    # Validate input
    if (!is(analysis, "TSENATAnalysis")) {
        stop("'analysis' must be a TSENATAnalysis object", call. = FALSE)
    }

    # Match metric and variability_metric arguments
    metric <- match.arg(metric)
    variability_metric <- match.arg(variability_metric)

    # Extract divergence SE from analysis object
    if (is.null(analysis@divergence_results)) {
        stop("Divergence results not found in analysis@divergence_results. ", "Run calculate_divergence_s4() first.",
            call. = FALSE)
    }

    # Handle both direct SE and wrapped 'divergence_se' key
    divergence_results_se <- if (is.list(analysis@divergence_results) && "divergence_se" %in%
        names(analysis@divergence_results)) {
        analysis@divergence_results$divergence_se
    } else if (is(analysis@divergence_results, "SummarizedExperiment")) {
        analysis@divergence_results
    } else {
        stop("Invalid divergence_results structure. Expected SummarizedExperiment or list with 'divergence_se' key",
            call. = FALSE)
    }

    if (nrow(divergence_results_se) == 0 || ncol(divergence_results_se) == 0) {
        stop("Divergence SummarizedExperiment is empty", call. = FALSE)
    }

    # Extract LM results for lm_res parameter (optional) Only use for ranking
    # if use_pvalue_ranking = TRUE
    lm_res <- NULL
    if (use_pvalue_ranking && !is.null(analysis@lm_results) && is.list(analysis@lm_results)) {
        if ("lm_interaction" %in% names(analysis@lm_results)) {
            lm_res <- analysis@lm_results$lm_interaction
        } else if (length(analysis@lm_results) > 0) {
            lm_res <- analysis@lm_results[[1]]
        }
    }

    # Validate LM results if using multi-gene mode
    if (is.null(gene) && use_pvalue_ranking && !is.null(lm_res)) {
        if (!is.data.frame(lm_res) || nrow(lm_res) == 0) {
            if (verbose) {
                message("Note: Invalid LM results. Plotting global curve without gene ranking.")
            }
            lm_res <- NULL
        }
    }

    # Calculate height if not provided (based on grid layout)
    if (is.null(height)) {
        n_rows <- ceiling(n_genes/ncol)
        height <- 3 + (3.5 * n_rows)  # 3' base + 3.5' per row
    }

    # Create the plot using base function
    p <- tryCatch({
        .plot_divergence_spectrum(divergence_results_se = divergence_results_se,
            gene = gene, lm_res = lm_res, n_genes = n_genes, ncol = ncol, metric = metric,
            variability_metric = variability_metric, ...)
    }, error = function(e) {
        if (verbose) {
            message("plot_divergence_spectrum failed: ", e$message)
        }
        return(NULL)
    })

    # If plot creation failed, return NULL invisibly
    if (is.null(p)) {
        return(invisible(NULL))
    }

    # Save plot to file if requested
    if (!is.null(output_file)) {
        save_analysis_output(p, output_file, object = analysis, verbose = verbose,
            func_name = "plot_divergence_spectrum_s4", width = width, height = height)
    }

    # Return file path if saved, otherwise return plot
    invisible(p)
}

# ============================================================================
# PLOT WRAPPER - Plot Method Concordance Comparison
# ============================================================================

#' Plot method concordance results from TSENATAnalysis
#'
#' @param analysis \code{TSENATAnalysis} object with computed method concordance
#'   (from \code{compute_method_concordance_s4()}).
#' @param verbose \code{logical}. Print progress messages. Default: FALSE
#'
#' @return A ggplot/cowplot object showing:
#'   \describe{
#'     \item{Panel 1}{Scatter plot of -log10(p-values) comparing methods}
#'     \item{Panel 2}{Histogram of p-value distributions by method}
#'   }
#'
#' @details
#' Creates visualization of method concordance including:
#' - Comparison of significance across two methods (with color-coded agreement)
#' - P-value distribution histograms for both methods
#' - Significance threshold lines at p < 0.05
#'
#' Requires that \code{compute_method_concordance_s4()} has already been run
#' to populate \code{@metadata$method_concordance}.
#'
#' @examples
#' # Load example data (matching TSENAT.Rmd workflow)
#' data(readcounts)
#' readcounts <- as.matrix(salmon_dataset)
#' mode(readcounts) <- 'numeric'
#' metadata_df <- read.table(
#'   system.file('extdata', 'metadata.tsv', package = 'TSENAT'),
#'   header = TRUE, sep = '\t'
#' )
#' gff3_dataset <- system.file('extdata', 'annotation.gff3.gz', package = 'TSENAT')
#' 
#' # Build analysis from vignette data and create small subset
#' analysis <- build_analysis_s4(readcounts = readcounts, tx2gene = gff3_dataset, metadata = metadata_df,
#'   tpm = salmon_tpm, effective_length = salmon_effective_length)
#' analysis <- filter_analysis_s4(analysis, min_samples = 1, subset_n_genes = 200)
#' # Note: compute_method_concordance_s4 requires additional LM and Friedman results
#' # For demo, we show that plot_method_concordance_s4 needs pre-computed concordance
#'
#' @aliases plot_method_concordance_s4
#' @export
setGeneric("plot_method_concordance_s4", function(analysis, verbose = FALSE) {
    standardGeneric("plot_method_concordance_s4")
})

#' @rdname plot_method_concordance_s4
setMethod("plot_method_concordance_s4", "TSENATAnalysis", function(analysis, verbose = FALSE) {

    # Load visualization dependencies (ggplot2, cowplot, etc.)
    .load_visualization_deps()

    # Validate that concordance results exist
    if (is.null(analysis@metadata$method_concordance)) {
        stop("[plot_method_concordance_s4] No concordance results found in @metadata.\n",
            "  Please run compute_method_concordance_s4() first.")
    }

    concordance_results <- analysis@metadata$method_concordance

    # Extract comparison dataframe
    comparison_df <- concordance_results$comparison_df

    if (is.null(comparison_df) || nrow(comparison_df) == 0) {
        stop("[plot_method_concordance_s4] Concordance comparison_df is empty or missing.")
    }

    if (verbose) {
        message("[plot_method_concordance_s4] Plotting concordance for ", nrow(comparison_df),
            " genes")
        message("[plot_method_concordance_s4] Methods compared: ", concordance_results$gam_method,
            " vs ", concordance_results$friedman_method)
    }

    # Call standard plotting function
    plot_obj <- .plot_method_concordance(comparison_df)

    if (verbose) {
        message("[plot_method_concordance_s4] Plot generated successfully")
    }

    return(plot_obj)
})

# ============================================================================
# OPTIMIZATION: Consolidated object extraction helper
# ============================================================================
# This helper consolidates redundant fallback extraction patterns into a single
# source of truth, reducing code duplication and improving maintainability.
#' @noRd
.extract_object_with_fallbacks <- function(obj, expected_class, key_name = NULL,
    verbose = FALSE) {
    # Single source of extraction logic for common pattern: Try direct class
    # match, then named list access, then list[1]

    if (is(obj, expected_class)) {
        if (verbose) {
            message("[extract_object] Found object via direct class match: ", expected_class)
        }
        return(obj)
    }

    if (is.list(obj)) {
        # Try named access first
        if (!is.null(key_name) && key_name %in% names(obj)) {
            if (verbose) {
                message("[extract_object] Found object via key: ", key_name)
            }
            return(obj[[key_name]])
        }

        # Fall back to first element
        if (length(obj) > 0) {
            if (verbose) {
                message("[extract_object] Using first element of list")
            }
            return(obj[[1]])
        }
    }

    if (verbose) {
        message("[extract_object] Could not extract object of class ", expected_class)
    }
    return(NULL)
}

#' Compute Effect Sizes from Divergence Results (S4 Wrapper)
#'
#' S4 wrapper for \code{.effect_sizes_divergence()} that extracts divergence and LM
#' results directly from a TSENATAnalysis object.
#'
#' @param analysis \code{TSENATAnalysis}. An S4 object containing divergence results
#'   (from \code{calculate_divergence_s4()}) and LM interaction results
#'   (from \code{calculate_lm_interaction_s4()}).
#'
#' @param significance_threshold \code{numeric}. Adjusted p-value threshold for
#'   filtering significant genes (default: 0.05).
#'
#' @param enrich_per_q_pattern \code{logical}. If TRUE, enriches results with per-q
#'   divergence patterns (default: TRUE).
#'
#' @param verbose \code{logical}. If TRUE, print diagnostic messages (default: TRUE).
#'
#' @param output_file \code{character} or \code{NULL}. Optional file path to save results.
#'   Supported formats: .rds (for S4 objects). Default: NULL (no file output).
#'
#' @param ... Additional arguments passed to the base function.
#'
#' @return Modified TSENATAnalysis with effect size results stored via
#'   \code{metadata(analysis)$effect_sizes_divergence}. Returns the analysis object visibly
#'   to support piping and method chaining.
#'
#' @details
#' **Workflow steps:**
#' \describe{
#'   \item{Extracting}{Divergence results (via accessors) and LM interaction results}
#'   \item{Computing}{Effect sizes using standard \code{.effect_sizes_divergence()} function}
#'   \item{Storing}{Results as list with \code{interaction_results} (data.frame) and
#'     \code{validation_stats}}
#'   \item{Tracking}{Function call in analysis metadata}
#' }
#'
#' **Parameter resolution priority** (explicit > metadata > default):
#' \itemize{
#'   \item \code{significance_threshold}: Uses explicit arg, else \code{metadata(analysis)$significance_threshold},
#'     else 0.05
#'   \item \code{enrich_per_q_pattern}: Uses explicit arg, else \code{metadata(analysis)$enrich_per_q_pattern},
#'     else TRUE
#'   \item \code{verbose}: Uses explicit arg, else \code{metadata(analysis)$verbose}, else TRUE
#' }
#'
#' Results are accessed via: \code{metadata(analysis)$effect_sizes_divergence}
#'
#' @examples
#' # Setup: Create test analysis with divergence and LM interaction results
#' data(readcounts)
#' readcounts <- as.matrix(salmon_dataset)
#' mode(readcounts) <- 'numeric'
#' metadata_df <- read.table(
#'   system.file('extdata', 'metadata.tsv', package = 'TSENAT'),
#'   header = TRUE, sep = '\t'
#' )
#' gff3_dataset <- system.file('extdata', 'annotation.gff3.gz', package = 'TSENAT')
#' analysis <- build_analysis_s4(readcounts = readcounts, tx2gene = gff3_dataset, metadata = metadata_df,
#'   tpm = salmon_tpm, effective_length = salmon_effective_length)
#' 
#' # Configure analysis parameters (best practice for reproducibility)
#' config <- tsenat_config(
#'   condition_col = 'condition',
#'   subject_col = 'paired_samples',
#'   paired = TRUE,
#'   control = 'normal',
#'   metadata = metadata_df
#' )
#' analysis <- setConfig(analysis, config)
#' 
#' analysis <- filter_analysis_s4(analysis, stringency = 'severe')
#' analysis <- calculate_diversity_s4(analysis, q = c(0.5, 1.0, 1.5), verbose = FALSE)
#' analysis <- calculate_divergence_s4(analysis, q = c(0.5, 1.0, 1.5), verbose = FALSE)
#' analysis <- calculate_lm_interaction_s4(analysis, method = 'gam', verbose = FALSE)
#'
#' # Compute effect sizes from divergence results
#' analysis <- effect_sizes_divergence_s4(analysis,
#'   significance_threshold = 0.05, verbose = FALSE)
#'
#' # Access results using metadata accessor
#' effect_size_results <- getMeta(analysis, 'effect_sizes_divergence')
#'
#' # View structure of results
#' str(effect_size_results, max.level = 1)
#'
#' @seealso
#' \code{\link{calculate_divergence_s4}} for divergence wrapper,
#' \code{\link{calculate_lm_interaction_s4}} for LM interaction wrapper
#'
#' @export
#' @importFrom methods is
#' @importFrom utils write.table
effect_sizes_divergence_s4 <- function(analysis, significance_threshold = 0.05, enrich_per_q_pattern = TRUE,
    verbose = FALSE, output_file = NULL, ...) {

    # Extract verbose parameter if not provided
    verbose <- resolve_slot_param(verbose, getConfig(analysis), "verbose", FALSE)

    # =========================================================================
    # INPUT VALIDATION
    # =========================================================================
    if (!is(analysis, "TSENATAnalysis")) {
        stop("'analysis' must be a TSENATAnalysis object", call. = FALSE)
    }

    # Check for required results
    if (length(divRes(analysis)) == 0) {
        stop("Divergence results required. Run calculate_divergence_s4() first.",
            call. = FALSE)
    }

    if (is.null(lmRes(analysis)) || length(lmRes(analysis)) == 0) {
        stop("LM results required. Run calculate_lm_interaction_s4() first.", call. = FALSE)
    }

    # =========================================================================
    # PARAMETER EXTRACTION using utility functions
    significance_threshold <- resolve_slot_param(significance_threshold, getConfig(analysis),
        "significance_threshold", 0.05)
    enrich_per_q_pattern <- resolve_slot_param(enrich_per_q_pattern, getConfig(analysis),
        "enrich_per_q_pattern", TRUE)
    verbose <- resolve_slot_param(verbose, getConfig(analysis), "verbose", FALSE)

    # =========================================================================
    # EXTRACT RESULTS FROM ANALYSIS OBJECT
    # =========================================================================

    # OPTIMIZATION: Use consolidated extraction helper instead of repeated
    # fallback chains Extract divergence SE
    analysis_divres <- divRes(analysis)
    divergence_se <- .extract_object_with_fallbacks(analysis_divres, "SummarizedExperiment",
        key_name = "divergence_se", verbose = verbose)

    if (is.null(divergence_se) || !is(divergence_se, "SummarizedExperiment")) {
        stop("Could not extract divergence SummarizedExperiment from divergence results",
            call. = FALSE)
    }

    # Extract LM results
    analysis_lmres <- lmRes(analysis)
    lm_res <- .extract_object_with_fallbacks(analysis_lmres, "data.frame", key_name = "lm_interaction",
        verbose = verbose)

    if (is.null(lm_res) || !is.data.frame(lm_res)) {
        stop("Could not extract LM results data.frame from LM results", call. = FALSE)
    }

    if (verbose) {
        message("[effect_sizes_divergence_s4] Extracted:")
        message("  - Divergence SE: ", paste(dim(divergence_se), collapse = " x "))
        message("  - LM results: ", nrow(lm_res), " genes")
    }

    # =========================================================================
    # COMPUTE EFFECT SIZES
    # =========================================================================
    if (verbose) {
        message("[effect_sizes_divergence_s4] Computing effect sizes...")
    }

    # =========================================================================
    # ENSURE DIVERGENCE SE HAS GENE NAMES FOR EFFECT SIZES
    # =========================================================================
    # The effect_sizes_divergence base function requires gene_name in rowData
    # Create a proper mapping that expands genes to match transcript-level rows

    rd <- rowData(divergence_se)
    if (is.null(rd) || !("gene_name" %in% colnames(rd))) {
        if (is.null(rd)) {
            rd <- DataFrame(row.names = rownames(divergence_se))
        }

        gene_names_added <- FALSE

        # Try tx2gene mapping first
        base_se <- getSE(analysis)
        tx2gene <- metadata(base_se)$tx2gene
        if (!is.null(tx2gene) && nrow(tx2gene) > 0 && ncol(tx2gene) >= 2) {
            # Find transcript and gene columns by name
            col_names <- tolower(colnames(tx2gene))
            tx_col <- NULL
            gene_col <- NULL

            # Find transcript column (first priority: Transcript, second: Tx)
            for (cn in colnames(tx2gene)) {
                if (tolower(cn) %in% c("transcript", "tx")) {
                  tx_col <- cn
                  break
                }
            }

            # Find gene column (first priority: Gene, second: Gen)
            for (cn in colnames(tx2gene)) {
                if (tolower(cn) %in% c("gene", "gen")) {
                  gene_col <- cn
                  break
                }
            }

            if (!is.null(tx_col) && !is.null(gene_col)) {
                # Use tx2gene to map rownames to genes
                tx_in_divergence <- rownames(divergence_se)
                tx_vector <- as.character(tx2gene[[tx_col]])
                gene_vector <- as.character(tx2gene[[gene_col]])

                # Find matching indices
                match_idx <- match(tx_in_divergence, tx_vector)

                # Get gene names using matched indices
                gene_names_expanded <- gene_vector[match_idx]

                # Check if all transcripts have valid mappings (no NAs)
                if (all(!is.na(gene_names_expanded))) {
                  rd$gene_name <- gene_names_expanded
                  gene_names_added <- TRUE
                  if (verbose) {
                    message("[effect_sizes_divergence_s4] Added gene_name to rowData from tx2gene mapping")
                  }
                }
            }
        }

        # If tx2gene didn't work, try direct assignment (for same nrow case)
        if (!gene_names_added && "gene" %in% colnames(lm_res)) {
            if (nrow(lm_res) == nrow(divergence_se)) {
                rd$gene_name <- as.character(lm_res$gene)
                gene_names_added <- TRUE
                if (verbose) {
                  message("[effect_sizes_divergence_s4] Added gene_name to rowData via direct assignment")
                }
            }
        }

        # Store updated rowData
        rowData(divergence_se) <- rd

        # Validate that gene_name was successfully added without NAs
        rd_final <- rowData(divergence_se)
        if (!("gene_name" %in% colnames(rd_final)) || any(is.na(rd_final$gene_name))) {
            stop("[effect_sizes_divergence_s4] Failed to add valid gene_name column to divergence_se rowData. ",
                "Ensure tx2gene metadata is properly set with matching transcript IDs.",
                call. = FALSE)
        }
    }

    result <- tryCatch({
        .effect_sizes_divergence(lm_res = lm_res, divergence_results_se = divergence_se,
            significance_threshold = significance_threshold, enrich_per_q_pattern = enrich_per_q_pattern,
            verbose = verbose, ...)
    }, error = function(e) {
        stop("effect_sizes_divergence", e$message, call. = FALSE)
    })

    # =========================================================================
    # STORE RESULTS IN ANALYSIS OBJECT
    # =========================================================================
    # Store result list in metadata (not in dedicated slot since none exists)
    current_meta <- getMeta(analysis)
    if (is.null(current_meta)) {
        analysis@metadata <- list()
    }

    analysis@metadata$effect_sizes_divergence <- result

    # Track function call
    current_calls <- getMeta(analysis, "function_calls")
    if (is.null(current_calls)) {
        current_calls <- character(0)
    }
    analysis@metadata$function_calls <- c(current_calls, paste0("effect_sizes_divergence[threshold=",
        significance_threshold, "]"))

    if (verbose) {
        message("[effect_sizes_divergence_s4] Results stored in metadata")
        if (!is.null(result$interaction_results)) {
            message("  - Effect size results: ", nrow(result$interaction_results),
                " genes")
        }
    }

    # Save if output_file provided
    if (!is.null(output_file)) {
        if (grepl("\\.tsv$|\\.csv$|\\.txt$", tolower(output_file))) {
            # Write effect size results as text table, filtering out all-NA
            # columns
            results_df <- as.data.frame(result$interaction_results)

            # Filter out columns that are all NA
            all_na_cols <- colnames(results_df)[vapply(results_df, function(x) all(is.na(x)),
                FUN.VALUE = logical(1))]
            if (length(all_na_cols) > 0) {
                results_df <- results_df[, !colnames(results_df) %in% all_na_cols]
                if (verbose) {
                  message("[effect_sizes_divergence_s4] Removed all-NA columns: ",
                    paste(all_na_cols, collapse = ", "))
                }
            }

            write.table(results_df, file = output_file, sep = "\t", quote = FALSE,
                row.names = FALSE)
        } else {
            # Default to RDS for S4 object
            saveRDS(analysis, file = output_file)
        }

    }

    analysis
}

#' Plot Top Transcripts from TSENATAnalysis Object
#'
#' S4 wrapper for \code{.plot_top_transcripts()} that extracts data directly from
#' a TSENATAnalysis object. Automatically retrieves the SummarizedExperiment and
#' LM results for visualizing transcript abundance across conditions.
#'
#' @param analysis \code{TSENATAnalysis}. An S4 object containing a processed
#'   SummarizedExperiment and optional LM interaction results.
#'
#' @param gene \code{character} or \code{NULL}. Gene identifier(s) to plot. If a vector 
#'   of multiple genes is provided, plots all of them. If NULL, automatically selects 
#'   the top genes from LM results based on \code{top_n} parameter (genes with lowest p-values).
#'   Default: NULL (auto-extract from lm_results).
#'
#' @param condition_col \code{character}. Column name in colData(se) specifying
#'   group assignments (default: 'sample_type').
#'
#' @param top_n \code{numeric}. Number of top transcripts to display for each
#'   condition (default: 3).
#'
#' @param output_file \code{character} or \code{NULL}. Path to save the plot as
#'   a PNG file. If NULL, saves to a temporary location (default: NULL).
#'
#' @param metric \code{character}. Method for ranking transcripts within genes.
#'   One of 'median', 'mean', 'variance', or 'iqr' (default: 'median').
#'
#' @param width \code{numeric} or \code{NULL}. Output image width in inches. 
#'   If NULL, automatically calculated based on number of genes (default: ~13 inches per column).
#'
#' @param height \code{numeric} or \code{NULL}. Output image height in inches.
#'   If NULL, automatically calculated based on number of genes (default: ~10 inches per row + headers).
#'
#' @param fontsize \code{numeric}. Base font size for heatmap titles and labels 
#'   (default: 16pt). Automatically scaled for readability.
#'
#' @param cellwidth \code{numeric}. Width of individual heatmap cells in pixels.
#'   If 0 (default), uses adaptive sizing based on data dimensions and layout.
#'   Set > 0 to override dynamic sizing.
#'
#' @param cellheight \code{numeric}. Height of individual heatmap cells in pixels.
#'   If 0 (default), uses adaptive sizing based on data dimensions and layout.
#'   Set > 0 to override dynamic sizing.
#'
#' @param layout_ncol \code{numeric}. Number of heatmaps per row in fixed layout
#'   (default: 2). If NULL, uses adaptive layout based on transcript counts.
#'
#' @param use_tpm \code{logical}. If \code{TRUE}, uses TPM (Transcripts Per Million) 
#'   from metadata instead of raw counts (default: FALSE). TPM is normalized for sequencing 
#'   depth and is recommended for comparing expression across samples. Requires TPM data 
#'   in metadata from `build_analysis_s4()` or `.build_se()` with `tpm` parameter. 
#'   Raises error if TPM not available and `use_tpm = TRUE`.
#'
#' @param verbose \code{logical}. If \code{TRUE}, print diagnostic messages
#'   during plotting (default: FALSE).
#'
#' @param output_file \code{character} or \code{NULL}. Optional file path to save the plot.
#'   Supported formats: .pdf, .png, .jpg. Default: NULL (no file output).
#' @param ... Additional arguments passed to the base plotting function.
#'
#' @return Invisibly returns the output file path (if `output_file` provided), or invisible(NULL) 
#'   if rendering to active graphics device. Graphics are rendered to the active grid device 
#'   for capture during vignette compilation.
#'
#' @details
#' This wrapper extracts the following from \code{analysis}:
#' \describe{
#'   \item{SummarizedExperiment}{From \code{analysis@se} containing transcript counts}
#'   \item{LM results}{From \code{analysis@lm_results$lm_interaction} for gene selection}
#' }
#'
#' If no gene is specified, the function automatically selects the top gene from
#' the LM results (lowest p-value). This simplifies visualization of genes with
#' significant q x condition interaction effects.
#'
#' @examples
#' # Plot 6: Top transcripts across groups
#' data(readcounts)
#' readcounts <- as.matrix(salmon_dataset)
#' mode(readcounts) <- 'numeric'
#' metadata_df <- read.table(
#'   system.file('extdata', 'metadata.tsv', package = 'TSENAT'),
#'   header = TRUE, sep = '\t'
#' )
#' gff3_dataset <- system.file('extdata', 'annotation.gff3.gz', package = 'TSENAT')
#' analysis <- build_analysis_s4(readcounts = readcounts, tx2gene = gff3_dataset, metadata = metadata_df,
#'   tpm = salmon_tpm, effective_length = salmon_effective_length)
#' 
#' # Configure analysis parameters
#' config <- tsenat_config(
#'   condition_col = 'condition',
#'   subject_col = 'paired_samples',
#'   paired = TRUE,
#'   control = 'normal',
#'   metadata = metadata_df
#' )
#' analysis <- setConfig(analysis, config)
#' 
#' analysis <- filter_analysis_s4(analysis, stringency = 'severe')
#' analysis <- calculate_diversity_s4(analysis, q = c(0.5, 1, 1.5), verbose = FALSE)
#' analysis <- calculate_lm_interaction_s4(analysis, method = 'gam', verbose = FALSE)
#' plot_file <- plot_top_transcripts_s4(analysis, top_n = 3)
#'
#' @seealso
#' \code{\link{TSENATAnalysis}} for object structure
#'
#' @export
#' @importFrom methods is
plot_top_transcripts_s4 <- function(analysis, gene = NULL, condition_col = NULL,
    top_n = 4, output_file = NULL, metric = c("median", "mean", "variance", "iqr"),
    use_tpm = TRUE, width = NULL, height = NULL, fontsize = 16, cellwidth = 0, cellheight = 0,
    layout_ncol = 2, verbose = FALSE, ...) {

    # Load visualization dependencies (ggplot2, cowplot, pheatmap, etc.)
    .load_visualization_deps()

    # Extract verbose parameter if not provided
    verbose <- resolve_slot_param(verbose, analysis@config, "verbose", FALSE)

    # =========================================================================
    # INPUT VALIDATION
    # =========================================================================
    if (!is(analysis, "TSENATAnalysis")) {
        stop("analysis must be a TSENATAnalysis object", call. = FALSE)
    }

    se <- analysis@se
    if (!inherits(se, "SummarizedExperiment")) {
        stop("analysis@se must be a SummarizedExperiment object", call. = FALSE)
    }

    # =========================================================================
    # AUTO-DETECT condition_col
    # =========================================================================
    if (is.null(condition_col)) {
        cd_cols <- colnames(colData(se))
        condition_col <- auto_detect_column(cd_cols, analysis@config, "condition_col",
            c("condition", "sample_type", "group", "treatment"), verbose = verbose,
            param_name = "condition_col")
    }

    # =========================================================================
    # EXTRACT LM RESULTS (for gene ranking if not specified)
    # =========================================================================
    lm_results_df <- NULL
    if (is.null(gene) && !is.null(analysis@lm_results)) {
        if ("lm_interaction" %in% names(analysis@lm_results)) {
            # Extract results data.frame from list structure
            if (is.data.frame(analysis@lm_results$lm_interaction)) {
                lm_results_df <- analysis@lm_results$lm_interaction
            } else if (is.list(analysis@lm_results$lm_interaction) && "results" %in%
                names(analysis@lm_results$lm_interaction)) {
                lm_results_df <- analysis@lm_results$lm_interaction$results
            }

            # Auto-select top gene from LM results
            if (!is.null(lm_results_df) && nrow(lm_results_df) > 0) {
                # Find p-value and gene columns
                p_col <- auto_detect_column(colnames(lm_results_df), analysis@config,
                  "p_col", c("p_interaction", "padj", "pvalue", "p.value", "p_value"),
                  verbose = FALSE, param_name = "p_col")

                gene_col <- auto_detect_column(colnames(lm_results_df), analysis@config,
                  "gene_col", c("gene", "gene_name", "gene_id"), verbose = FALSE,
                  param_name = "gene_col")

                if (!is.null(p_col) && p_col %in% colnames(lm_results_df)) {
                  # Get top genes (sorted by p-value, select top_n)
                  top_indices <- order(lm_results_df[[p_col]])[seq_len(min(top_n,
                    nrow(lm_results_df)))]
                  gene <- as.character(lm_results_df[top_indices, gene_col])

                  if (verbose) {
                    message("[plot_top_transcripts_s4] Auto-selected top ", length(gene),
                      " genes from LM results")
                    message("  ", paste(gene, collapse = ", "))
                  }
                }
            }
        }
    }

    if (is.null(gene)) {
        stop("[plot_top_transcripts_s4] No gene specified and cannot auto-detect from LM results. ",
            "Provide gene explicitly.", call. = FALSE)
    }

    # =========================================================================
    # CALL BASE FUNCTION
    # =========================================================================
    if (verbose) {
        if (length(gene) > 1) {
            message("[plot_top_transcripts_s4] Calling .plot_top_transcripts() for genes: ",
                paste(gene, collapse = ", "))
        } else {
            message("[plot_top_transcripts_s4] Calling .plot_top_transcripts() for gene: ",
                gene)
        }
    }

    plot_file <- tryCatch({
        .plot_top_transcripts(se = se, gene = gene, condition_col = condition_col,
            res = lm_results_df, top_n = top_n, output_file = output_file, metric = metric[1],
            use_tpm = use_tpm, width = width, height = height, fontsize = fontsize,
            cellwidth = cellwidth, cellheight = cellheight, layout_ncol = layout_ncol,
            ...)
    }, error = function(e) {
        stop("[plot_top_transcripts_s4]", conditionMessage(e), call. = FALSE)
    })

    # Optionally save to file if output_file provided
    if (!is.null(output_file)) {
        if (verbose) {
            message("[plot_top_transcripts_s4] Plot saved to: ", output_file)
        }
    }

    # Always return the plot object (ggplot) Knitr will auto-manage figure
    # rendering
    invisible(plot_file)
}

#' Plot Tsallis Divergence Effect Size Distribution (S4 Wrapper)
#'
#' S4 wrapper that extracts effect size results from a TSENATAnalysis object
#' and generates a histogram visualization of Tsallis divergence effect sizes
#' across genes.
#'
#' @param analysis \code{TSENATAnalysis} object with effect sizes computed
#'   (typically via \code{\link{effect_sizes_divergence_s4}}).
#' @param threshold \code{numeric}. Effect size threshold for visual marking
#'   in the plot. Default is 0.1 (information-theoretic significance level).
#' @param output_file \code{character}. Optional file path to save the plot.
#'   If NULL, plot is returned but not saved.
#' @param width \code{numeric}. Plot width in inches. Default is 10.
#' @param height \code{numeric}. Plot height in inches. Default is 6.
#' @param verbose \code{logical}. Print status messages. Default is TRUE.
#' @param ... Additional arguments passed to the underlying plotting function.
#'
#' @return
#' Invisibly returns the file path if saved, otherwise the ggplot object.
#' If the plot cannot be created (missing data, ggplot2 not available),
#' returns NULL invisibly with an informative message.
#'
#' @details
#' This wrapper extracts the interaction results (with effect size columns)
#' from \code{analysis@metadata$effect_sizes_divergence$interaction_results}
#' and passes them to the base \code{.plot_divergence_distribution()} function.
#'
#' The function visualizes the distribution of effect sizes using the median
#' q-value's divergence (typically around q=1.0, close to Shannon entropy).
#' A red dashed line marks the information-theoretic significance threshold.
#'
#' **Data Requirements:**
#' \itemize{
#'   \item Effect sizes must be computed via \code{effect_sizes_divergence_s4()}
#'   \item \code{@metadata$effect_sizes_divergence$interaction_results} must
#'         contain columns matching pattern \code{effect_size_D_q*}
#' }
#'
#' @examples
#' # Plot 2: Distribution of effect sizes across genes
#' data(readcounts)
#' readcounts <- as.matrix(salmon_dataset)
#' mode(readcounts) <- 'numeric'
#' metadata_df <- read.table(
#'   system.file('extdata', 'metadata.tsv', package = 'TSENAT'),
#'   header = TRUE, sep = '\t'
#' )
#' gff3_dataset <- system.file('extdata', 'annotation.gff3.gz', package = 'TSENAT')
#' analysis <- build_analysis_s4(readcounts = readcounts, tx2gene = gff3_dataset, metadata = metadata_df,
#'   tpm = salmon_tpm, effective_length = salmon_effective_length)
#' 
#' # Configure analysis parameters
#' config <- tsenat_config(
#'   condition_col = 'condition',
#'   subject_col = 'paired_samples',
#'   paired = TRUE,
#'   control = 'normal',
#'   metadata = metadata_df
#' )
#' analysis <- setConfig(analysis, config)
#' 
#' analysis <- filter_analysis_s4(analysis, stringency = 'severe')
#' analysis <- calculate_diversity_s4(analysis, q = c(0.5, 1, 1.5), verbose = FALSE)
#' analysis <- calculate_divergence_s4(analysis, q = c(0.5, 1, 1.5), verbose = FALSE)
#' analysis <- calculate_lm_interaction_s4(analysis, method = 'gam', verbose = FALSE)
#' analysis <- effect_sizes_divergence_s4(analysis, verbose = FALSE)
#' p_dist <- plot_divergence_distribution_s4(analysis, verbose = FALSE)
#' print(p_dist)
#'
#' @seealso
#' \code{\link{effect_sizes_divergence_s4}} for computing effect sizes.
#'
#' @export
plot_divergence_distribution_s4 <- function(analysis, threshold = 0.1, output_file = NULL,
    width = 12, height = 6, verbose = TRUE, ...) {

    # Load visualization dependencies (ggplot2, cowplot, etc.)
    .load_visualization_deps()

    # Extract verbose parameter if not provided
    verbose <- resolve_slot_param(verbose, analysis@config, "verbose", TRUE)

    # Validate input
    if (!is(analysis, "TSENATAnalysis")) {
        stop("'analysis' must be a TSENATAnalysis object", call. = FALSE)
    }

    # Extract effect sizes from metadata
    if (is.null(analysis@metadata$effect_sizes_divergence)) {
        stop("Effect sizes not found in analysis@metadata$effect_sizes_divergence. ",
            "Run effect_sizes_divergence_s4() first.", call. = FALSE)
    }

    effect_sizes <- analysis@metadata$effect_sizes_divergence

    # Extract interaction results (with effect size columns)
    if (!is.list(effect_sizes) || is.null(effect_sizes$interaction_results)) {
        stop("Invalid effect size structure. Expected @metadata$effect_sizes_divergence$interaction_results",
            call. = FALSE)
    }

    interaction_results <- effect_sizes$interaction_results

    if (!is.data.frame(interaction_results) || nrow(interaction_results) == 0) {
        stop("interaction_results must be a non-empty data frame", call. = FALSE)
    }

    # Create the plot using base function
    p <- tryCatch({
        .plot_divergence_distribution(interaction_results = interaction_results,
            threshold = threshold, ...)
    }, error = function(e) {
        if (verbose) {
            message("plot_divergence_distribution failed: ", e$message)
        }
        return(NULL)
    })

    # If plot creation failed, return NULL invisibly
    if (is.null(p)) {
        return(invisible(NULL))
    }

    # Save to file if requested
    if (!is.null(output_file)) {
        save_analysis_output(p, output_file, object = analysis, verbose = verbose,
            func_name = "plot_divergence_distribution_s4", width = width, height = height)
    }

    # Always return the plot object (not file path) Knitr will auto-manage
    # figure rendering File is saved separately if output_file provided
    invisible(p)
}


#' Prepare Gene Switching Tables from TSENATAnalysis Object
#'
#' S4 wrapper for \code{.prepare_gene_switching_tables()} that extracts results
#' directly from a TSENATAnalysis object. Automatically retrieves LM results and
#' jackknife switching results from the analysis object slots.
#'
#' @param analysis \code{TSENATAnalysis}. An S4 object containing completed
#'   LM interaction and jackknife isoform switching analyses.
#'
#' @param n_top_genes \code{numeric} or \code{NULL}. Number of top genes
#'   (by adjusted p-value) to include in output tables. If \code{NULL},
#'   all genes with significant LM results are included.
#'
#' @param n_transcripts_per_gene \code{numeric}. Maximum number of transcripts
#'   to display per gene in the output tables (default: 10).
#'
#' @param verbose \code{logical}. If \code{TRUE}, print diagnostic messages
#'   during table preparation.
#'
#' @param ... Additional arguments passed to the base function.
#'
#' @return A list containing:
#'   \describe{
#'     \item{\code{$summary_table}}{Gene-level summary with LM p-values and
#'           significant q-values}
#'     \item{\code{$transcript_tables}}{Named list of data.frames, one per gene,
#'           showing transcript-level switching metrics}
#'     \item{\code{$q_vector}}{Vector of q-values analyzed}
#'   }
#'
#' @details
#' This function extracts the following from \code{analysis}:
#' \describe{
#'   \item{LM results}{From \code{analysis@lm_results$lm_interaction$results}}
#'   \item{Jackknife results}{From \code{analysis@jackknife_results} or
#'         extracted from the switching analysis metadata}
#' }
#'
#' The wrapper automatically handles column detection and parameter extraction,
#' providing a simplified interface compared to the base function.
#'
#' @examples
#' data(readcounts)
#' readcounts <- as.matrix(salmon_dataset)
#' mode(readcounts) <- 'numeric'
#' metadata_df <- read.table(
#'   system.file('extdata', 'metadata.tsv', package = 'TSENAT'),
#'   header = TRUE, sep = '\t'
#' )
#' gff3_dataset <- system.file('extdata', 'annotation.gff3.gz', package = 'TSENAT')
#' analysis <- build_analysis_s4(readcounts = readcounts, tx2gene = gff3_dataset, metadata = metadata_df,
#'   tpm = salmon_tpm, effective_length = salmon_effective_length)
#' 
#' # Configure analysis parameters
#' config <- tsenat_config(
#'   condition_col = 'condition',
#'   subject_col = 'paired_samples',
#'   paired = TRUE,
#'   control = 'normal',
#'   metadata = metadata_df
#' )
#' analysis <- setConfig(analysis, config)
#' 
#' analysis <- filter_analysis_s4(analysis, stringency = 'severe')
#' analysis <- calculate_diversity_s4(analysis, q = c(0.5, 1.0, 1.5), verbose = FALSE)
#' analysis <- calculate_lm_interaction_s4(analysis, method = 'gam', verbose = FALSE)
#' analysis <- jackknife_isoform_switching_s4(analysis, n_bootstrap = 50,
#'   verbose = FALSE)
#' tables <- prepare_gene_switching_tables_s4(analysis)
#' head(tables$summary_df)
#'
#' @param output_file \code{character} or \code{NULL}. Optional file path to save results.
#'   Supported formats: .rds (for S4 objects), .tsv, .csv, .txt (for tables). Default: NULL (no file output).
#'
#' @export
#' @importFrom methods is
#' @importFrom utils write.table
prepare_gene_switching_tables_s4 <- function(analysis, n_top_genes = NULL, n_transcripts_per_gene = 10,
    verbose = FALSE, output_file = NULL, ...) {

    # Extract verbose parameter if not provided
    verbose <- resolve_slot_param(verbose, analysis@config, "verbose", FALSE)

    # Validation
    if (!is(analysis, "TSENATAnalysis")) {
        stop("analysis must be a TSENATAnalysis object")
    }

    # Extract LM results
    if (verbose)
        message("Extracting LM results from analysis object...")

    lm_results_list <- analysis@lm_results
    if (is.null(lm_results_list) || length(lm_results_list) == 0) {
        stop("No LM results found in analysis@lm_results. Run calculate_lm_interaction_s4() first.")
    }

    # Try to get lm_interaction results first, fallback to first available
    if (!is.null(lm_results_list$lm_interaction)) {
        if (is.data.frame(lm_results_list$lm_interaction$results)) {
            lm_res <- lm_results_list$lm_interaction$results
        } else if (is.data.frame(lm_results_list$lm_interaction)) {
            lm_res <- lm_results_list$lm_interaction
        } else {
            stop("Cannot extract LM results from analysis@lm_results$lm_interaction")
        }
    } else if (is.data.frame(lm_results_list)) {
        lm_res <- lm_results_list
    } else {
        stop("Cannot find LM results data.frame in analysis@lm_results")
    }

    if (verbose)
        message("  [OK] Extracted LM results with ", nrow(lm_res), " genes")

    # Extract jackknife/switching results
    if (verbose)
        message("Extracting jackknife switching results from analysis object...")

    jackknife_results_list <- analysis@jackknife_results
    if (is.null(jackknife_results_list) || length(jackknife_results_list) == 0) {
        stop("No jackknife results found in analysis@jackknife_results. Run jackknife_isoform_switching_s4() first.")
    }

    # Check if results are stored under 'multi_q' key
    # (jackknife_isoform_switching_multiq class)
    if ("multi_q" %in% names(jackknife_results_list)) {
        multi_q_object <- jackknife_results_list[["multi_q"]]

        # If it's a tsenat_isoform_switching_multiq object, use it directly
        if (inherits(multi_q_object, "tsenat_isoform_switching_multiq")) {
            multi_q_results <- multi_q_object
            if (verbose) {
                message("  [OK] Found multi-q results under 'multi_q' key with ",
                  length(multi_q_results), " q-values")
            }
        } else {
            stop("Element at jackknife_results$multi_q is not a tsenat_isoform_switching_multiq object")
        }
    } else {
        # Fallback: look for q-keyed results directly Extract only the q-keyed
        # results (filter by pattern 'q_X_XX')
        q_key_pattern <- "^q_[0-9]+_[0-9]{2}$"
        q_keyed_results <- jackknife_results_list[grep(q_key_pattern, names(jackknife_results_list))]

        if (length(q_keyed_results) == 0) {
            stop("No q-keyed jackknife results found in analysis@jackknife_results. ",
                "Expected keys in format 'q_X_XX' (e.g., 'q_0_01', 'q_1_00') or 'multi_q'.")
        }

        # Wrap q-keyed results as a multi_q object for consistency
        multi_q_results <- q_keyed_results
        if (verbose) {
            message("  [OK] Found ", length(multi_q_results), " q-keyed results")
        }
    }

    if (verbose)
        message("  [OK] Extracted jackknife results with ", length(multi_q_results),
            " q-values")

    # Call base function with extracted parameters
    if (verbose)
        message("Calling .prepare_gene_switching_tables()...")

    result <- .prepare_gene_switching_tables(lm_res = lm_res, multi_q_results = multi_q_results,
        n_top_genes = n_top_genes, n_transcripts_per_gene = n_transcripts_per_gene,
        verbose = verbose, ...)

    if (verbose)
        message("[OK] Gene switching tables prepared successfully")

    # Track function call in metadata
    analysis@metadata$function_calls <- c(analysis@metadata$function_calls, "prepare_gene_switching_tables_s4")
    analysis@metadata$function_timestamps <- c(analysis@metadata$function_timestamps,
        as.character(Sys.time()))

    # Save if output_file provided
    if (!is.null(output_file)) {
        if (grepl("\\.tsv$|\\.csv$|\\.txt$", tolower(output_file))) {
            # Write as TSV/CSV
            write.table(result, file = output_file, sep = "\t", quote = FALSE, row.names = TRUE)
        } else {
            # Default to RDS for arbitrary objects
            saveRDS(result, file = output_file)
        }

    }

    return(result)
}

#' Plot Multi-Q Delta Influence Heatmaps from TSENATAnalysis Object
#'
#' S4 wrapper for \code{.plot_multiq_delta_influence_heatmaps()} that extracts results
#' directly from a TSENATAnalysis object. Automatically retrieves jackknife switching
#' results from the analysis object slots.
#'
#' @param analysis \code{TSENATAnalysis}. An S4 object containing completed
#'   jackknife isoform switching analysis across multiple q-values.
#'
#' @param n_genes \code{numeric}. Number of top genes to display in heatmaps
#'   (default: 4). Genes are ranked by LM p-values if available, otherwise
#'   by order of appearance in results.
#'
#' @param lm_results \code{data.frame} or \code{NULL}. Optional LM interaction
#'   results for ranking genes (default: NULL). If NULL, attempts to extract from
#'   \code{analysis@lm_results$lm_interaction}.
#'
#' @param verbose \code{logical}. If \code{TRUE}, print diagnostic messages
#'   during plot generation (default: FALSE).
#'
#' @param output_file \code{character} or \code{NULL}. Optional file path to save the plot.
#'   Default: NULL (no file output).
#'
#' @param ... Additional arguments passed to the base function.
#'
#' @return A file path (character) to the saved heatmap PNG file, invisibly.
#'
#' @details
#' This function extracts the following from \code{analysis}:
#' \describe{
#'   \item{Jackknife results}{From \code{analysis@jackknife_results}, which should
#'         contain multi-q switching results keyed by q-value (e.g., 'q_1.00')}
#'   \item{LM results}{From \code{analysis@lm_results$lm_interaction} if not
#'         explicitly provided, for ranking genes by significance}
#' }
#'
#' The wrapper automatically handles parameter extraction and provides a simplified
#' interface compared to the base function.
#'
#' @examples
#' # Plot 5: Multi-q delta influence (isoform switching) heatmaps
#' data(readcounts)
#' readcounts <- as.matrix(salmon_dataset)
#' mode(readcounts) <- 'numeric'
#' metadata_df <- read.table(
#'   system.file('extdata', 'metadata.tsv', package = 'TSENAT'),
#'   header = TRUE, sep = '\t'
#' )
#' gff3_dataset <- system.file('extdata', 'annotation.gff3.gz', package = 'TSENAT')
#' analysis <- build_analysis_s4(readcounts = readcounts, tx2gene = gff3_dataset, metadata = metadata_df,
#'   tpm = salmon_tpm, effective_length = salmon_effective_length)
#' 
#' # Configure analysis parameters
#' config <- tsenat_config(
#'   condition_col = 'condition',
#'   subject_col = 'paired_samples',
#'   paired = TRUE,
#'   control = 'normal',
#'   metadata = metadata_df
#' )
#' analysis <- setConfig(analysis, config)
#' 
#' analysis <- filter_analysis_s4(analysis, stringency = 'severe')
#' analysis <- calculate_diversity_s4(analysis, q = c(0.5, 1, 1.5), verbose = FALSE)
#' analysis <- calculate_divergence_s4(analysis, q = c(0.5, 1, 1.5), verbose = FALSE)
#' analysis <- calculate_lm_interaction_s4(analysis, method = 'gam', verbose = FALSE)
#' analysis <- jackknife_isoform_switching_s4(analysis, q = c(0.5, 1, 1.5),
#'   n_bootstrap = 50, verbose = FALSE)
#' heatmap_file <- plot_multiq_delta_influence_heatmaps_s4(analysis, n_genes = 2)
#'
#' @seealso
#' \code{\link{jackknife_isoform_switching_s4}} for computing switching results
#'
#' @export
#' @importFrom methods is
plot_multiq_delta_influence_heatmaps_s4 <- function(analysis, n_genes = 4, lm_results = NULL,
    verbose = FALSE, output_file = NULL, ...) {

    # Load visualization dependencies (ggplot2, cowplot, pheatmap, etc.)
    .load_visualization_deps()

    # Extract verbose parameter if not provided
    verbose <- resolve_slot_param(verbose, analysis@config, "verbose", FALSE)

    # Validation
    if (!is(analysis, "TSENATAnalysis")) {
        stop("analysis must be a TSENATAnalysis object", call. = FALSE)
    }

    if (verbose)
        message("Extracting jackknife switching results from analysis object...")

    # Extract jackknife/switching results
    jackknife_results_list <- analysis@jackknife_results
    if (is.null(jackknife_results_list) || length(jackknife_results_list) == 0) {
        stop("No jackknife results found in analysis@jackknife_results. ", "Run jackknife_isoform_switching_s4() first.",
            call. = FALSE)
    }

    # Check for multi-q result (stored under 'multi_q' key when multiple
    # q-values provided)
    if ("multi_q" %in% names(jackknife_results_list)) {
        switching_results <- jackknife_results_list$multi_q
        if (verbose) {
            message("  [OK] Found multi-q result with class: ", class(switching_results)[1])
        }
    } else {
        # Fallback: use all results as list (for single or multiple q-values)
        switching_results <- jackknife_results_list
        if (verbose) {
            message("  [OK] Using individual q-value results (", length(switching_results),
                " q-values)")
        }
    }

    if (verbose) {
        message("  Q-values: ", paste(names(switching_results), collapse = ", "))
    }

    # Extract LM results if not provided
    if (is.null(lm_results)) {
        if (verbose)
            message("Extracting LM results from analysis@lm_results...")

        lm_results_list <- analysis@lm_results
        if (!is.null(lm_results_list)) {
            if (!is.null(lm_results_list$lm_interaction)) {
                if (is.data.frame(lm_results_list$lm_interaction$results)) {
                  lm_results <- lm_results_list$lm_interaction$results
                } else if (is.data.frame(lm_results_list$lm_interaction)) {
                  lm_results <- lm_results_list$lm_interaction
                }
            }

            if (!is.null(lm_results)) {
                if (verbose)
                  message("  [OK] Extracted LM results with ", nrow(lm_results),
                    " genes")
            } else if (verbose) {
                message("  LM results not found; genes will be ranked by appearance")
            }
        }
    }

    if (verbose)
        message("Calling .plot_multiq_delta_influence_heatmaps()...")

    # Call base function with extracted parameters Note: output_file parameter
    # can be used to save heatmap as PNG file
    result <- .plot_multiq_delta_influence_heatmaps(switching_results = switching_results,
        n_genes = n_genes, lm_results = lm_results, verbose = verbose, output_file = output_file,
        ...)

    if (verbose) {
        message("[OK] Heatmap plot generated successfully")
    }

    invisible(result)
}

#' Plot GAM q-curves from TSENATAnalysis object
#'
#' S4 wrapper that accepts a TSENATAnalysis object and generates GAM q-curve plots
#' @param analysis \code{TSENATAnalysis} object with diversity and LM interaction results.
#' @param n_top \code{integer}. Number of top genes (by adjusted p-value) to plot 
#'   (default: 6). Only used if genes = NULL.
#'
#' @param genes \code{character} vector. Optional specific gene names to plot. 
#'   If provided, these genes are plotted directly regardless of significance.
#'   If NULL (default), top n_top significant genes are selected.
#'
#' @param condition_col \code{character}. Column name in colData(se) specifying 
#'   group assignments for samples. If NULL, attempts to auto-detect from @config
#'   (looks for 'condition_col' or 'condition'). If still NULL, defaults to 'sample_type'.
#'
#' @param sig_alpha \code{numeric}. Significance threshold for adjusted p-values 
#'   (default: 0.05). Only used if genes = NULL; filters lm_res to significant 
#'   genes before selecting top n.
#'
#' @param assay_name \code{character}. Name of the assay in se to extract 
#'   (default: 'diversity').
#'
#' @param output_file \code{character} or \code{NULL}. Optional file path to save the plot.
#'   Default: NULL (no file output).
#'
#' @param width \code{numeric}. Width of the plot in inches (default: 12).
#'   Only used if output_file is not NULL.
#'
#' @param height \code{numeric} or \code{NULL}. Height of the plot in inches.
#'   Default: NULL (automatically calculated based on width and aspect ratio).
#'   Only used if output_file is not NULL.
#'
#' @param ... Additional arguments passed to the base function.
#'
#' @return A single \code{ggplot} object with all selected genes arranged in a 
#'   grid layout. Can be saved with \code{ggplot2::ggsave()}.
#'
#' @details
#' This wrapper automatically:
#' 1. Extracts SummarizedExperiment from \code{@se} slot
#' 2. Extracts LM results from \code{@lm_results$lm_interaction} slot
#' 3. Detects condition_col from \code{@config} or uses default
#' 4. Calls \code{.plot_lm_interaction_gam()} with extracted parameters
#'
#' **Parameter Resolution (condition_col):**
#' \enumerate{
#'   \item Explicit \code{condition_col} parameter (highest priority)
#'   \item \code{@config$condition_col} if available
#'   \item \code{@config$condition} if available  
#'   \item Default: 'sample_type'
#' }
#'
#' @seealso
#' \code{\link{calculate_lm_interaction_s4}} for running LM analysis on TSENATAnalysis.
#'
#' @examples
#' # Plot 3: GAM q-curves for genes with q-by-condition interactions
#' data(readcounts)
#' readcounts <- as.matrix(salmon_dataset)
#' mode(readcounts) <- 'numeric'
#' metadata_df <- read.table(
#'   system.file('extdata', 'metadata.tsv', package = 'TSENAT'),
#'   header = TRUE, sep = '\t'
#' )
#' gff3_dataset <- system.file('extdata', 'annotation.gff3.gz', package = 'TSENAT')
#' analysis <- build_analysis_s4(readcounts = readcounts, tx2gene = gff3_dataset, metadata = metadata_df,
#'   tpm = salmon_tpm, effective_length = salmon_effective_length)
#' 
#' # Configure analysis parameters
#' config <- tsenat_config(
#'   condition_col = 'condition',
#'   subject_col = 'paired_samples',
#'   paired = TRUE,
#'   control = 'normal',
#'   metadata = metadata_df
#' )
#' analysis <- setConfig(analysis, config)
#' 
#' analysis <- filter_analysis_s4(analysis, stringency = 'severe')
#' analysis <- calculate_diversity_s4(analysis, q = seq(0.2, 2.5, by = 0.15), verbose = FALSE)
#' analysis <- calculate_lm_interaction_s4(analysis, method = 'gam', verbose = FALSE)
#' 
#' p_gam <- plot_lm_interaction_gam_s4(analysis, n_top = 2, sig_alpha = 0.15)
#' print(p_gam)
#'
#' @export
plot_lm_interaction_gam_s4 <- function(analysis, n_top = 6, genes = NULL, condition_col = NULL,
    sig_alpha = 0.05, assay_name = "diversity", output_file = NULL, width = 12, height = NULL,
    ...) {
    # Load visualization dependencies (ggplot2, cowplot, mgcv, etc.)
    .load_visualization_deps()

    # =========================================================================
    # INPUT VALIDATION
    # =========================================================================
    if (!is(analysis, "TSENATAnalysis")) {
        stop("'analysis' must be a TSENATAnalysis object", call. = FALSE)
    }

    # Check that LM results exist
    if (is.null(analysis@lm_results) || is.null(analysis@lm_results$lm_interaction)) {
        stop("[plot_lm_interaction_gam_s4] No LM interaction results found in @lm_results$lm_interaction. ",
            "Run calculate_lm_interaction_s4() first.", call. = FALSE)
    }

    lm_res <- analysis@lm_results$lm_interaction

    if (!is.data.frame(lm_res)) {
        stop("[plot_lm_interaction_gam_s4] @lm_results$lm_interaction must be a data.frame",
            call. = FALSE)
    }

    # Check that diversity results exist (needed for SE reconstruction)
    if (length(analysis@diversity_results) == 0) {
        stop("[plot_lm_interaction_gam_s4] No diversity results found in @diversity_results. ",
            "Run calculate_diversity_s4() first.", call. = FALSE)
    }

    # =========================================================================
    # AUTO-DETECT condition_col FROM @config IF NOT PROVIDED
    # =========================================================================
    if (is.null(condition_col)) {
        cd_cols <- colnames(colData(analysis@se))
        condition_col <- auto_detect_column(cd_cols, analysis@config, "condition_col",
            c("condition", "sample_type", "group", "treatment"), verbose = verbose,
            param_name = "condition_col")
    }

    # Validate that condition_col exists in colData
    if (!(condition_col %in% colnames(colData(analysis@se)))) {
        stop("[plot_lm_interaction_gam_s4] Specified condition_col='", condition_col,
            "' not found in colData. Available columns: ", paste(colnames(colData(analysis@se)),
                collapse = ", "), call. = FALSE)
    }

    # =========================================================================
    # RECONSTRUCT COMBINED DIVERSITY SE FOR PLOTTING (Same approach as in
    # calculate_lm_interaction_s4)
    # =========================================================================
    # Extract q-values from diversity_results keys
    q_keys <- names(analysis@diversity_results)
    q_computed <- as.numeric(sub("^q_", "", q_keys))

    # Reconstruct combined diversity SE with all q-values
    diversity_combined <- tryCatch({
        verbose <- if ("verbose" %in% names(analysis@config)) {
            analysis@config$verbose
        } else {
            FALSE
        }

        .calculate_diversity(x = analysis@se, q = sort(q_computed), norm = TRUE,
            verbose = verbose, bootstrap = FALSE)
    }, error = function(e) {
        stop("[plot_lm_interaction_gam_s4] Failed to reconstruct diversity SE:\n",
            conditionMessage(e), call. = FALSE)
    })

    # =========================================================================
    # EXTRACT model_data FROM STORED RESULTS
    # =========================================================================
    model_data <- NULL
    if ("lm_interaction_model_data" %in% names(analysis@lm_results)) {
        model_data <- analysis@lm_results$lm_interaction_model_data
    }

    # =========================================================================
    # CALCULATE HEIGHT IF NOT PROVIDED
    # =========================================================================
    if (is.null(height)) {
        # Estimate number of genes to be plotted
        if (!is.null(genes)) {
            n_genes_plot <- length(genes)
        } else {
            # Count significant genes
            if ("adj_p_interaction" %in% colnames(lm_res)) {
                sig_genes <- lm_res$adj_p_interaction <= sig_alpha
            } else if ("p_interaction" %in% colnames(lm_res)) {
                sig_genes <- lm_res$p_interaction <= sig_alpha
            } else {
                sig_genes <- rep(TRUE, nrow(lm_res))
            }
            n_genes_plot <- min(sum(sig_genes), n_top)
        }
        # Calculate height: 2 rows per 3-gene group, ~3.5 inches per row
        n_rows <- ceiling(n_genes_plot/2)
        height <- 2 + (3.5 * n_rows)
    }

    # =========================================================================
    # CALL plot_lm_interaction_gam WITH RECONSTRUCTED DIVERSITY SE
    # =========================================================================
    result <- tryCatch({
        .plot_lm_interaction_gam(se = diversity_combined, lm_res = lm_res, condition_col = condition_col,
            n_top = n_top, genes = genes, sig_alpha = sig_alpha, assay_name = assay_name,
            model_data = model_data, output_file = output_file, width = width, height = height,
            ...)
    }, error = function(e) {
        stop("[plot_lm_interaction_gam_s4]", conditionMessage(e), call. = FALSE)
    })

    # =========================================================================
    # RETURN PLOT
    # =========================================================================
    # Track that plotting occurred
    if (is.list(analysis@metadata)) {
        analysis@metadata$function_calls <- c(analysis@metadata$function_calls, paste0("plot_lm_interaction_gam_s4[n_top=",
            n_top, ", condition_col=", condition_col, "]"))
    }

    # Save plot to file if requested
    if (!is.null(output_file)) {
        save_analysis_output(result, output_file, object = analysis, verbose = verbose,
            func_name = "plot_lm_interaction_gam_s4", width = width, height = height)
    }

    # Return the plot object directly (not the analysis object)
    result
}

#' Jackknife isoform switching analysis on TSENATAnalysis object
#'
#' S4 wrapper that accepts a TSENATAnalysis object and performs jackknife-based
#' isoform switching detection. Automatically extracts the SummarizedExperiment
#' and metadata from object slots.
#'
#' @param analysis \code{TSENATAnalysis} object containing:
#'   \itemize{
#'     \item \code{@se}: SummarizedExperiment with count data
#'     \item \code{@config}: Configuration metadata
#'   }
#'
#' @param condition_col \code{character}. Column name in colData(se) specifying 
#'   group assignments (default: 'sample_type'). If NULL, attempts auto-detection.
#'
#' @param subject_col \code{character}. Optional column for paired/repeated measures design.
#'   If provided, enables paired analysis. Default: NULL (unpaired).
#'
#' @param gene_col \code{character}. Column name in rowData(se) or metadata identifying genes.
#'   Default: 'gene'.
#'
#' @param isoform_col \code{character}. Column name in rowData(se) or metadata identifying 
#'   isoforms/transcripts. Default: 'transcript' or 'isoform'.
#'
#' @param q \code{numeric}. Tsallis entropy parameter(s) to analyze. Can be single value 
#'   or vector for multi-q analysis (default: 1).
#'
#' @param norm \code{logical}. Whether to use normalized diversity values 
#'   (default: TRUE).
#'
#' @param log_base \code{numeric}. Logarithm base for entropy calculations. Default: NULL (uses e).
#'
#' @param pseudocount \code{numeric}. Pseudocount value for count regularization. Default: NULL (uses @config$pseudocount or 0).
#'
#' @param threshold \code{numeric}. Percentile threshold for detecting transcript switching 
#'   (default: 90). Transcripts with delta_influence >= threshold percentile are classified 
#'   as 'switching'.
#'
#' @param n_bootstrap \code{integer}. Number of bootstrap resamples for confidence 
#'   intervals (default: 1000).
#'
#' @param lm_results \code{data.frame}. Optional LM interaction results to filter genes.
#'   If provided, only genes in lm_results are analyzed.
#'
#' @param lm_p_threshold \code{numeric}. P-value threshold for filtering genes from 
#'   lm_results (default: 0.05).
#'
#' @param use_lm_fdr \code{logical}. If TRUE, uses adjusted p-values from lm_results 
#'   (default: TRUE).
#'
#' @param verbose \code{logical}. Print progress messages (default: FALSE).
#' @param output_file \code{character} or \code{NULL}. Optional file path to save results.
#'   Supported formats: .tsv, .csv, .txt (for tables), or .rds (for S4 objects).
#'   When text format is specified, generates TWO files:
#'   \itemize{
#'     \item \code{output_file}: Gene-level summary (one row per gene with switching statistics)
#'     \item \code{output_file_transcripts.ext}: Transcript-level details (one row per transcript with p-values and FDR)
#'   }
#'   For .rds format, saves only the full analysis object.
#'   Default: NULL (no file output).
#' @param ... Additional arguments for future extensibility.
#'
#' @return \code{TSENATAnalysis} object with jackknife results stored in \code{@jackknife_results}
#'   slot. Results are keyed by q-value (e.g., 'q_1.00'). For multi-q analysis, multiple
#'   calls will accumulate results in the slot.
#'
#'   The analysis object is returned visibly to support method chaining:
#'   \preformatted{
#'     analysis <- jackknife_isoform_switching_s4(analysis, q = 0.5)
#'     analysis <- jackknife_isoform_switching_s4(analysis, q = 1.0)
#'   }
#'
#' @details
#' This wrapper automatically:
#' 1. Extracts SummarizedExperiment from \code{@se} slot
#' 2. Detects condition_col, gene_col, isoform_col from colData/rowData or @config
#' 3. Calls \code{.jackknife_isoform_switching()} with extracted parameters
#'
#' **Parameter Auto-Detection:**
#' \enumerate{
#'   \item \code{condition_col}: Uses explicit parameter, then @config, then 'sample_type'
#'   \item \code{gene_col}: Uses explicit parameter, then looks for 'gene' or 'Gene'
#'   \item \code{isoform_col}: Uses explicit parameter, then looks for 'transcript', 'isoform', or 'Isoform'
#' }
#'
#' @seealso \code{jackknife_isoform_switching} for the underlying implementation,
#' \code{\link{TSENATAnalysis}} for object structure.
#'
#' @examples
#' data(readcounts)
#' metadata_df <- read.table(system.file('extdata', 'metadata.tsv', package = 'TSENAT'),
#'                           header = TRUE, sep = '\t')
#' gff3_file <- system.file('extdata', 'annotation.gff3.gz', package = 'TSENAT')
#' analysis <- build_analysis_s4(readcounts = salmon_dataset, tx2gene = gff3_file, metadata = metadata_df,
#'                              tpm = salmon_tpm, effective_length = salmon_effective_length)
#' analysis <- filter_analysis_s4(analysis, min_samples = 1, subset_n_genes = 20, subset_n_samples = 8)
#' analysis <- calculate_diversity_s4(analysis, q = 1, verbose = FALSE)
#' 
#' @export
jackknife_isoform_switching_s4 <- function(analysis, condition_col = NULL, subject_col = NULL,
    gene_col = NULL, isoform_col = NULL, q = 1, norm = NULL, log_base = NULL, threshold = 90,
    n_bootstrap = 1000, pseudocount = NULL, lm_results = NULL, lm_p_threshold = 0.05,
    use_lm_fdr = TRUE, output_file = NULL, verbose = FALSE, ...) {
    # =========================================================================
    # INPUT VALIDATION
    # =========================================================================
    if (!is(analysis, "TSENATAnalysis")) {
        stop("'analysis' must be a TSENATAnalysis object", call. = FALSE)
    }

    se <- analysis@se

    if (!inherits(se, "SummarizedExperiment")) {
        stop("[jackknife_isoform_switching_s4] @se must be a SummarizedExperiment object",
            call. = FALSE)
    }

    # =========================================================================
    # AUTO-DETECT condition_col (using centralized helper)
    # =========================================================================
    cd_cols <- colnames(colData(se))
    condition_col <- auto_detect_column(cd_cols, config_list = analysis@config, config_key = "condition_col",
        priority_candidates = c("sample_type", "condition", "group", "sample_group"),
        default_fallback = if (length(cd_cols) > 0)
            cd_cols[1] else NULL, verbose = verbose, param_name = "condition_col")

    if (is.null(condition_col)) {
        stop("[jackknife_isoform_switching_s4] Cannot auto-detect condition_col.\n",
            "  Available colData columns: ", paste(cd_cols, collapse = ", "), "\n\n",
            "SOLUTION: Set @config$condition_col or pass explicit parameter\n", "  Example: analysis@config$condition_col <- 'sample_type'\n",
            "  Or:      jackknife_isoform_switching_s4(analysis, condition_col = 'sample_type')\n",
            call. = FALSE)
    }

    # =========================================================================
    # AUTO-DETECT gene_col AND isoform_col FROM rowData (using centralized
    # helper)
    # =========================================================================
    rd <- if (!is.null(rowData(se)) && nrow(rowData(se)) > 0) {
        rowData(se)
    } else {
        NULL
    }

    rd_cols <- if (!is.null(rd))
        colnames(rd) else character(0)

    # Detect gene_col
    if (is.null(gene_col)) {
        gene_col <- auto_detect_column(rd_cols, config_list = analysis@config, config_key = "gene_col",
            priority_candidates = c("gene_id", "gene", "Gene", "gene_name"), default_fallback = "gene",
            verbose = verbose, param_name = "gene_col")
    }

    # Detect isoform_col
    if (is.null(isoform_col)) {
        isoform_col <- auto_detect_column(rd_cols, config_list = analysis@config,
            config_key = "isoform_col", priority_candidates = c("transcript_id",
                "transcript", "isoform", "Isoform", "tx_id"), default_fallback = "transcript",
            verbose = verbose, param_name = "isoform_col")
    }

    # =========================================================================
    # EXTRACT LM_RESULTS IF PROVIDED VIA ANALYSIS OBJECT
    # =========================================================================
    if (is.null(lm_results) && !is.null(analysis@lm_results)) {
        # Try to extract LM results from analysis object
        if ("lm_interaction" %in% names(analysis@lm_results)) {
            lm_results <- analysis@lm_results$lm_interaction
            if (verbose) {
                message("[jackknife_isoform_switching_s4] Using LM interaction results from @lm_results")
            }
        }
    }

    # =========================================================================
    # CHECK Q-VALUES AVAILABILITY IN DIVERSITY RESULTS
    # =========================================================================
    if (!is.null(analysis@diversity_results) && length(analysis@diversity_results) >
        0) {
        available_q_keys <- names(analysis@diversity_results)
        available_q <- as.numeric(sub("^q_", "", available_q_keys))
        available_q <- sort(unique(available_q))

        # Check if requested q-values are available
        q_vals <- if (is.numeric(q))
            q else c(q)
        missing_q <- setdiff(q_vals, available_q)

        if (length(missing_q) > 0) {
            warning("[jackknife_isoform_switching_s4] Requested q-values not in @diversity_results:\n",
                "  Requested: ", paste(q_vals, collapse = ", "), "\n", "  Available: ",
                paste(available_q, collapse = ", "), "\n", "  Missing:   ", paste(missing_q,
                  collapse = ", "), "\n\n", "SOLUTION: Recompute diversity with all desired q-values before calling jackknife_isoform_switching_s4:\n",
                "  analysis <- calculate_diversity_s4(analysis, q = c(", paste(q_vals,
                  collapse = ", "), "), norm = TRUE)\n", "  analysis <- jackknife_isoform_switching_s4(analysis, q = c(",
                paste(q_vals, collapse = ", "), "))\n", call. = FALSE)
        }

        if (verbose && length(missing_q) == 0) {
            message("[jackknife_isoform_switching_s4] All requested q-values available in diversity results")
        }
    }

    # =========================================================================
    # RESOLVE PARAMETERS
    # =========================================================================
    norm <- resolve_slot_param(norm, analysis@config, "norm", TRUE)
    log_base <- resolve_slot_param(log_base, analysis@config, "log_base", exp(1))
    pseudocount <- resolve_slot_param(pseudocount, analysis@config, "pseudocount",
        0)

    # =========================================================================
    # VALIDATE n_bootstrap PARAMETER
    # =========================================================================
    if (!is.numeric(n_bootstrap) || length(n_bootstrap) != 1 || n_bootstrap < 1) {
        stop("'n_bootstrap' must be a positive integer", call. = FALSE)
    }

    if (n_bootstrap < 50) {
        warning("n_bootstrap = ", n_bootstrap, " is less than the recommended minimum of 50. ",
            "Results may be unreliable. Consider increasing to at least 50-100 for stable estimates.",
            call. = FALSE)
    }

    # =========================================================================
    # CALL BASE FUNCTION
    # =========================================================================
    result <- tryCatch({
        .jackknife_isoform_switching(se = se, condition_col = condition_col, subject_col = subject_col,
            gene_col = gene_col, isoform_col = isoform_col, q = q, norm = norm, log_base = log_base,
            threshold = threshold, n_bootstrap = n_bootstrap, pseudocount = pseudocount,
            verbose = verbose, lm_results = lm_results, lm_p_threshold = lm_p_threshold,
            use_lm_fdr = use_lm_fdr)
    }, error = function(e) {
        stop("[jackknife_isoform_switching_s4] Jackknife analysis failed:\\n", conditionMessage(e),
            call. = FALSE)
    })

    # =========================================================================
    # STORE RESULTS IN ANALYSIS OBJECT (OPTIMIZED - vectorized q-value storage)
    # =========================================================================
    # Ensure q is a vector and validate
    if (!is.numeric(q) && !(is.vector(q) && all(vapply(q, is.numeric, FUN.VALUE = logical(1))))) {
        stop("[jackknife_isoform_switching_s4] 'q' must be numeric or numeric vector",
            call. = FALSE)
    }
    q_vals <- if (is.numeric(q))
        q else as.numeric(c(q))

    # Check if result has multi-q class
    if (inherits(result, "tsenat_isoform_switching_multiq")) {
        # Multi-q result: store as-is (already has correct underscore-format
        # keys from base function)
        for (q_key in names(result)) {
            analysis@jackknife_results[[q_key]] <- result[[q_key]]
        }
        # Also store the entire multi-q result object for plotting function
        # access
        analysis@jackknife_results[["multi_q"]] <- result
        if (verbose) {
            message(sprintf("[jackknife_isoform_switching_s4] Stored multi-q result with keys: %s",
                paste(names(result), collapse = ", ")))
        }
    } else {
        # Store results for each q-value (vectorized - no explicit loop) Use
        # underscore format to match base function: paste0('q_', gsub('\\.',
        # '_', sprintf('%.2f', q)))
        q_keys <- paste0("q_", gsub("\\.", "_", sprintf("%.2f", q_vals)))

        # Store each result
        for (i in seq_along(q_keys)) {
            q_key <- q_keys[i]

            # Check if result is list with named q-values (should have
            # underscore format now)
            if (is.list(result) && q_key %in% names(result)) {
                analysis@jackknife_results[[q_key]] <- result[[q_key]]
            } else if (length(q_vals) == 1) {
                # Single q-value: store result directly
                analysis@jackknife_results[[q_key]] <- result
            } else {
                # Multiple q-values: log warning if not found
                warning("[jackknife_isoform_switching_s4] Result for q=", q_vals[i],
                  " (key: ", q_key, ") not found in base function result.", call. = FALSE)
            }

            if (verbose) {
                message(sprintf("[jackknife_isoform_switching_s4] Stored results for %s",
                  q_key))
            }
        }
    }

    # =========================================================================
    # TRACK FUNCTION CALL IN METADATA (OPTIMIZED - single paste())
    # =========================================================================
    if (is.list(analysis@metadata)) {
        # Pre-format all metadata in one call (more efficient)
        call_str <- sprintf("jackknife_isoform_switching_s4[q=%s, condition_col=%s]",
            paste(q_vals, collapse = ","), condition_col)

        analysis@metadata$function_calls <- c(analysis@metadata$function_calls, call_str)
        analysis@metadata$function_timestamps <- c(analysis@metadata$function_timestamps,
            as.character(Sys.time()))
    }

    # =========================================================================
    # SAVE RESULTS IF output_file PROVIDED
    # =========================================================================
    if (!is.null(output_file)) {
        output_dir <- dirname(output_file)
        if (output_dir != "." && !dir.exists(output_dir)) {
            dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
        }

        if (grepl("\\.tsv$|\\.csv$|\\.txt$", tolower(output_file))) {
            # Write as text table (extract from jackknife results)
            tryCatch({
                output_ext <- sub("^.*\\.", ".", tolower(output_file))
                output_base <- sub(paste0(output_ext, "$"), "", output_file)
                transcript_file <- paste0(output_base, "_transcripts", output_ext)

                # Extract gene-level and transcript-level tables using helper
                gene_data <- extract_multiq_table(result, is_multiq = inherits(result,
                  "tsenat_isoform_switching_multiq"), extract_fn = function(res,
                  q_key) {
                  if (!is.null(res$summary_table)) {
                    return(res$summary_table)
                  }
                  return(NULL)
                }, q_value_col = "q_value")

                transcript_data <- extract_multiq_table(result, is_multiq = inherits(result,
                  "tsenat_isoform_switching_multiq"), extract_fn = function(res,
                  q_key) {
                  if (!is.null(res$all_transcript_stats)) {
                    return(res$all_transcript_stats)
                  }
                  return(NULL)
                }, q_value_col = "q_value")

                # Write gene-level to original file
                if (!is.null(gene_data)) {
                  write.table(gene_data, file = output_file, sep = "\t", quote = FALSE,
                    row.names = FALSE)
                }

                # Write transcript-level to separate file
                if (!is.null(transcript_data)) {
                  write.table(transcript_data, file = transcript_file, sep = "\t",
                    quote = FALSE, row.names = FALSE)
                }
            }, error = function(e) {
                warning("[jackknife_isoform_switching_s4] Could not write jackknife results to file: ",
                  conditionMessage(e), call. = FALSE)
            })
        } else {
            # Default to RDS for S4 object
            saveRDS(analysis, file = output_file)
        }
    }

    # Return modified analysis object (visibly for method chaining as
    # documented)
    analysis
}

#' M-Estimation for Sample Quality (S4 Wrapper)
#'
#' S4 wrapper for \code{m_estimate} that performs robust M-estimation
#' on diversity results stored in a TSENATAnalysis object and stores results
#' back into the object.
#'
#' @param analysis \code{TSENATAnalysis} object with diversity results
#'   (typically via \code{\link{calculate_diversity_s4}}).
#' @param condition_col \code{character}. Column name in sample metadata indicating
#'   condition/sample grouping. Auto-detected from \code{@config$condition_col}
#'   if available.
#' @param loss_type \code{character}. Type of loss function: 'huber' (default),
#'   'tukey', or 'lsq'. Determines robustness vs efficiency trade-off.
#' @param scale \code{numeric}. Manual scale parameter. If NULL, estimated from data.
#' @param max_iter \code{integer}. Maximum iterations for M-estimation. Default: 50.
#' @param tol \code{numeric}. Convergence tolerance. Default: 1e-6.
#' @param paired \code{logical}. If TRUE, adjusts degrees of freedom for paired
#'   designs. Auto-detected from \code{@config$paired} if available.
#'   Default: FALSE.
#' @param pcorr \code{character}. P-value correction method. Default: 'BH' (Benjamini-Hochberg).
#' @param q_combine_method \code{character}. How to collapse multi-q results:
#'   'mean' (default) or 'median'.
#' @param influence_threshold \code{numeric}. Threshold for classifying samples
#'   as high-influence. Default: 0.75.
#' @param scale_method \code{character}. Scale estimation method: 'mad' (default),
#'   'proposal2', or 's-estimator'.
#' @param output_file \code{character} or \code{NULL}. Optional file path to save results.
#'   Supported formats: .rds (for S4 objects), .tsv, .csv, .txt (for tables). Default: NULL (no file output).
#' @param verbose \code{logical}. Print status messages. Default: TRUE.
#'
#' @return
#' Modified TSENATAnalysis object with M-estimation results stored in
#' \code{analysis@metadata$m_estimate_results}. Contains data frame with
#' influence scores, robustness weights, entropy statistics, and QC classifications.
#' Returns visibly to support method chaining and piping.
#'
#' @details
#' This wrapper extracts diversity results from \code{analysis@diversity_results},
#' performs robust M-estimation on diversity values (entropy), and stores results
#' in the analysis object metadata.
#'
#' **M-Estimation:** Robust regression technique that down-weights outliers based
#' on their residuals. Useful for detecting low-quality samples that show
#' unusual diversity patterns.
#'
#' **M-Estimation Results include:**
#' \itemize{
#'   \item \code{sample_influence}: How much each sample affects the overall fit
#'   \item \code{robustness_weight}: Down-weighting factor (lower = more outlying)
#'   \item \code{entropy_mean}: Average entropy for the sample
#'   \item \code{entropy_sd}: Entropy variability within the sample
#'   \item \code{Status}: QC Classification ('OK' or 'Flag for QC' based on influence_threshold)
#' }
#'
#' **Parameter resolution priority** (explicit > @config > error):
#' \itemize{
#'   \item \code{samples}: Uses explicit arg, else \code{@config$condition_col},
#'     else error. Note: despite parameter name 'samples', maps to condition grouping column
#' }
#'
#' **Data Requirements:**
#' \itemize{
#'   \item Diversity results must be computed via \code{calculate_diversity_s4()}
#'   \item Sample grouping column required in colData (auto-detected from @config$condition_col
#'     or via 'samples' parameter)
#' }
#'
#' @examples
#' # Create test analysis and compute M-estimation
#' data(readcounts)
#' readcounts <- as.matrix(salmon_dataset)
#' mode(readcounts) <- 'numeric'
#' metadata_df <- read.table(
#'   system.file('extdata', 'metadata.tsv', package = 'TSENAT'),
#'   header = TRUE, sep = '\t'
#' )
#' gff3_dataset <- system.file('extdata', 'annotation.gff3.gz', package = 'TSENAT')
#' analysis <- build_analysis_s4(readcounts = readcounts, tx2gene = gff3_dataset, metadata = metadata_df,
#'   tpm = salmon_tpm, effective_length = salmon_effective_length)
#' analysis <- filter_analysis_s4(analysis, min_samples = 1, subset_n_genes = 200)
#' analysis <- calculate_diversity_s4(analysis, q = c(0.5, 1.0, 1.5), verbose = FALSE)
#' analysis <- m_estimate_s4(
#'   analysis,
#'   condition_col = 'condition',
#'   loss_type = 'huber',
#'   verbose = FALSE
#' )
#'
#' @seealso
#' \code{\link{calculate_diversity_s4}} for computing diversity
#'
#' @export
#' @importFrom utils write.table
m_estimate_s4 <- function(analysis, condition_col = NULL, loss_type = "huber", scale = NULL,
    max_iter = 50, tol = 1e-06, paired = NULL, pcorr = "BH", q_combine_method = "mean",
    influence_threshold = 0.75, scale_method = "mad", output_file = NULL, verbose = FALSE) {

    # Validate input
    if (!is(analysis, "TSENATAnalysis")) {
        stop("'analysis' must be a TSENATAnalysis object", call. = FALSE)
    }

    # Check for diversity results
    if (is.null(analysis@diversity_results) || length(analysis@diversity_results) ==
        0) {
        stop("Diversity results not found. Run calculate_diversity_s4() first.",
            call. = FALSE)
    }

    # Validate that diversity_results is a properly structured named list
    if (!is.list(analysis@diversity_results) || is.null(names(analysis@diversity_results))) {
        stop("Diversity results must be a named list of SummarizedExperiment objects",
            call. = FALSE)
    }

    # Auto-detect condition_col if not provided
    if (is.null(condition_col)) {
        if ("condition_col" %in% names(analysis@config)) {
            condition_col <- analysis@config$condition_col
            if (is.null(condition_col) || !is.character(condition_col) || condition_col ==
                "") {
                stop("@config$condition_col must be a non-empty character value",
                  call. = FALSE)
            }
            if (verbose) {
                message(sprintf("Auto-detected 'condition_col' from config: %s",
                  condition_col))
            }
        } else {
            cd_cols <- colnames(SummarizedExperiment::colData(analysis@diversity_results[[1]]))
            stop("Sample grouping column not specified:\n", "  Available colData columns: ",
                paste(cd_cols, collapse = ", "), "\n\n", "SOLUTION: Set @config$condition_col or pass 'condition_col' parameter\n",
                "  Example: analysis@config$condition_col <- 'sample_type'\n", "  Or:      m_estimate_s4(analysis, condition_col = 'sample_type')\n",
                call. = FALSE)
        }
    } else if (!is.character(condition_col) || length(condition_col) != 1) {
        stop("'condition_col' must be a single character value", call. = FALSE)
    }

    # Auto-detect paired parameter from @config if not explicitly provided
    paired <- resolve_slot_param(paired, analysis@config, "paired", FALSE)

    if (!is.logical(paired) || length(paired) != 1) {
        stop("'paired' must be a single logical value (TRUE or FALSE)", call. = FALSE)
    }

    # Extract diversity results - get first SE to access sample metadata
    diversity_se <- analysis@diversity_results[[1]]

    if (is.null(diversity_se) || nrow(diversity_se) == 0) {
        stop("Diversity SummarizedExperiment is empty", call. = FALSE)
    }

    # Verify condition_col exists
    sample_info <- SummarizedExperiment::colData(diversity_se)
    if (!(condition_col %in% colnames(sample_info))) {
        stop("Column '", condition_col, "' not found in sample metadata.\n", "Available columns: ",
            paste(colnames(sample_info), collapse = ", "), "\n\n", "SOLUTION: Use a valid column name\n",
            "  Example: m_estimate_s4(analysis, condition_col = 'sample_type')\n",
            call. = FALSE)
    }

    # Combine all q-value diversity results into a single matrix (m_estimate
    # needs all diversity data in one SE)
    if (verbose) {
        message(sprintf("Combining %d q-value diversity results...", length(analysis@diversity_results)))
    }

    first_se <- analysis@diversity_results[[1]]
    combined_assay <- SummarizedExperiment::assay(first_se)
    combined_colnames <- colnames(first_se)

    # Add other q-values
    for (q_name in names(analysis@diversity_results)[-1]) {
        se_q <- analysis@diversity_results[[q_name]]
        combined_assay <- cbind(combined_assay, SummarizedExperiment::assay(se_q))
        combined_colnames <- c(combined_colnames, colnames(se_q))
    }

    # Update column names to reflect combined data
    colnames(combined_assay) <- combined_colnames
    single_colData <- SummarizedExperiment::colData(first_se)

    # Replicate colData for each q-value
    n_q_values <- length(analysis@diversity_results)
    combined_colData <- do.call(rbind, replicate(n_q_values, single_colData, simplify = FALSE))
    rownames(combined_colData) <- combined_colnames

    # Create combined SummarizedExperiment
    combined_se <- SummarizedExperiment::SummarizedExperiment(assays = list(counts = combined_assay),
        colData = combined_colData)

    # Run m_estimate
    if (verbose) {
        message("Running M-estimation on combined diversity...")
    }

    m_est_results <- tryCatch({
        result <- .m_estimate(x = combined_se, samples = condition_col, loss_type = loss_type,
            scale = scale, max_iter = max_iter, tol = tol, paired = paired, pcorr = pcorr,
            q_combine_method = q_combine_method, influence_threshold = influence_threshold,
            scale_method = scale_method)
        result
    }, error = function(e) {
        stop("M-estimation failed:\n", e$message, call. = FALSE)
    })

    # Store results in metadata
    analysis@metadata$m_estimate_results <- m_est_results

    # Save to output file if provided
    if (!is.null(output_file)) {
        if (!is.character(output_file) || length(output_file) != 1) {
            stop("'output_file' must be a character string (file path)", call. = FALSE)
        }
        # Use save_analysis_output for consistent handling (data frame for TSV,
        # RDS for S4)
        save_analysis_output(m_est_results, output_file, object = analysis, verbose = verbose,
            func_name = "m_estimate_s4")
    }

    # Track function call
    analysis@metadata$function_calls <- c(analysis@metadata$function_calls, paste0("m_estimate_s4[condition_col=",
        condition_col, ",loss_type=", loss_type, "]"))

    if (verbose) {
        message("M-estimation complete. Results stored in @metadata$m_estimate_results")
    }

    invisible(analysis)
}

#' Filter Low-Abundance Transcripts in a TSENATAnalysis Object
#'
#' S4 wrapper for \code{.filter_se()} that filters low-abundance transcripts
#' directly within a \code{TSENATAnalysis} object. This maintains the consistent
#' S4 workflow pattern where functions accept and return analysis objects.
#'
#' @param analysis A \code{TSENATAnalysis} S4 object containing the
#'   \code{SummarizedExperiment} to be filtered.
#'
#' @param min_tpm Numeric TPM threshold (default 1.0).
#'   Keeps transcripts with TPM >= \code{min_tpm} in >= \code{min_samples} samples.
#'   Ignored if \code{stringency} is specified.
#'
#' @param tpm_assay_name Character; name of assay containing TPM data (default: NULL).
#'   If NULL, searches for TPM assay automatically.
#'
#' @param min_samples Numeric. Minimum number of samples in which a transcript
#'   must be present (default: 5). Ignored if \code{stringency} is specified.
#'
#' @param stringency Character. Filtering stringency level: 'soft' (permissive),
#'   'medium' (balanced), or 'severe' (stringent). When specified, auto-estimates:
#'   \code{min_samples}, \code{min_tpm}, \code{min_tx_per_gene}, and \code{min_isoform_abundance}
#'   from data. Requires \code{pair_col} in colData for paired designs.
#'   User-provided values for any parameter override stringency defaults.
#'   Default: NULL (use explicit parameters).
#'
#' @param pair_col Character; column name in colData containing pair IDs for paired designs.
#'   Default: NULL (auto-detect if needed).
#'
#' @param min_tx_per_gene Integer minimum number of transcripts per gene required
#'   (default 2L). Single-transcript genes are always kept. Ignored if \code{stringency}
#'   is specified; when specified, automatically adjusted based on stringency level.
#'
#' @param min_isoform_abundance Numeric in [0, 1]; minimum relative abundance threshold
#'   for isoforms within each gene. Implements Soneson et al. (2016) filtering.
#'   Default behavior:
#'   - If \code{stringency} is specified: uses stringency-based default (soft: 0.01, medium: 0.05, severe: 0.15)
#'   - If \code{stringency} is NULL: uses default 0.05 (5%)
#'   - If explicitly provided: overrides any stringency default
#'   Set to 0 or NULL (post-stringency processing) to skip isoform-level filtering.
#'
#' @param assay_name Character; name or index of the assay to use for filtering
#'   (default: 'counts'). Deprecated: use \code{tpm_assay_name} instead.
#'
#' @param subset_n_genes Integer; optional number of genes to retain after filtering.
#'   If provided, genes are selected based on \code{subset_select_by}. Default: NULL.
#'
#' @param subset_genes Character vector; optional specific genes to retain after filtering.
#'   Default: NULL.
#'
#' @param subset_n_samples Integer; optional number of samples to retain after filtering.
#'   If provided, samples are selected (balanced by condition if available).
#'   Default: NULL.
#'
#' @param subset_samples Character vector; optional specific samples to retain after filtering.
#'   Default: NULL.
#'
#' @param subset_select_by Character; gene selection method for \code{subset_n_genes}:
#'   'variance' (highest variance), 'mean' (highest mean expression), or 'random'.
#'   Default: 'variance'.
#'
#' @param subset_seed Integer; random seed for reproducibility when \code{subset_select_by = 'random'}.
#'   Default: 42.
#'
#' @param subset_min_count Numeric; optional minimum count threshold applied during subsetting.
#'   Default: NULL.
#'
#' @param verbose Logical. If TRUE, print filtering progress and summary statistics
#'   (default: FALSE).
#'
#' @return Invisibly returns the modified \code{analysis} object with filtered
#'   \code{SummarizedExperiment} in the \code{@se} slot. The filtering operation
#'   modifies the analysis object in-place while maintaining all other slots
#'   (results, metadata, etc.).
#'
#' @details
#' This wrapper applies \code{.filter_se()} to the SummarizedExperiment within
#' the TSENATAnalysis object, optionally followed by subsetting parameters.
#' The filtering and subsetting operations are applied in sequence:
#'
#' 1. Extracts the SE from \code{analysis@se}
#' 2. Filters using \code{.filter_se()} with specified filtering parameters
#' 3. If any subset parameters are provided, applies gene/sample selection
#'    to select specific genes and/or samples
#' 4. Stores the filtered/subsetted SE back in \code{analysis@se}
#' 5. Returns the modified analysis object invisibly
#'
#' **Important:** Filtering should be performed BEFORE computing diversity,
#' divergence, or LM interaction results. If called after analysis results
#' have been computed, those results will be based on unfiltered data and
#' may not align with the filtered SE dimensions.
#'
#' @seealso
#' \code{\link{build_analysis_s4}} for creating a new analysis object
#'
#' @examples
#' # Create test analysis and filter
#' data(readcounts)
#' readcounts <- as.matrix(salmon_dataset)
#' mode(readcounts) <- 'numeric'
#' metadata_df <- read.table(
#'   system.file('extdata', 'metadata.tsv', package = 'TSENAT'),
#'   header = TRUE, sep = '\t'
#' )
#' gff3_dataset <- system.file('extdata', 'annotation.gff3.gz', package = 'TSENAT')
#' analysis <- build_analysis_s4(readcounts = readcounts, tx2gene = gff3_dataset, metadata = metadata_df,
#'   tpm = salmon_tpm, effective_length = salmon_effective_length)
#' analysis <- filter_analysis_s4(analysis, stringency = 'medium')
#'
#' @export
filter_analysis_s4 <- function(analysis, min_tpm = 1, tpm_assay_name = NULL, min_samples = 5L,
    stringency = NULL, pair_col = NULL, min_tx_per_gene = 2L, min_isoform_abundance = NULL,
    assay_name = "counts", subset_n_genes = NULL, subset_genes = NULL, subset_n_samples = NULL,
    subset_samples = NULL, subset_select_by = c("variance", "mean", "random"), subset_seed = 42,
    subset_min_count = NULL, verbose = FALSE) {
    # Validate input
    if (!inherits(analysis, "TSENATAnalysis")) {
        stop("analysis must be a TSENATAnalysis object", call. = FALSE)
    }

    # Extract SE from analysis
    se <- analysis@se

    # Apply filtering via .filter_se() with all parameters
    se_filtered <- .filter_se(se = se, min_tpm = min_tpm, tpm_assay_name = tpm_assay_name,
        min_samples = min_samples, stringency = stringency, pair_col = pair_col,
        min_tx_per_gene = min_tx_per_gene, min_isoform_abundance = min_isoform_abundance,
        assay_name = assay_name, verbose = verbose)

    # Store filtered SE back in analysis object
    analysis@se <- se_filtered

    # Apply optional subsetting after filtering
    has_subset_params <- !is.null(subset_n_genes) || !is.null(subset_genes) || !is.null(subset_n_samples) ||
        !is.null(subset_samples) || !is.null(subset_min_count)

    if (has_subset_params) {
        if (verbose) {
            message("Applying subset_analysis to filtered data...")
        }

        # Use match.arg to validate subset_select_by
        subset_select_by <- match.arg(subset_select_by)

        # Apply subsetting via .subset_analysis
        analysis <- .subset_analysis(analysis = analysis, n_genes = subset_n_genes,
            n_samples = subset_n_samples, genes = subset_genes, samples = subset_samples,
            select_by = subset_select_by, seed = subset_seed, min_count = subset_min_count,
            verbose = verbose)
    }

    # Return modified analysis object
    analysis
}

#' Build a Complete TSENATAnalysis Object
#'
#' Convenience wrapper that combines \code{.build_se()} and \code{TSENATAnalysis()}
#' into a single function call. This creates a complete analysis object ready for
#' Tsallis entropy computation and downstream analysis.
#'
#' @param readcounts A matrix or data.frame of transcript-level read counts with
#'   transcript IDs as row names and sample names as column names. Typically output
#'   from quantification tools (SALMON, kallisto, etc.). Optional when \code{salmon_dir}
#'   is provided (in which case readcounts are auto-loaded from quant.sf files).
#'
#' @param salmon_dir Optional character path to directory containing Salmon quantification
#'   output. Expected structure: \code{salmon_dir/sample_name/quant.sf}. When provided,
#'   automatically discovers and reads all quant.sf files. If both \code{readcounts} and
#'   \code{salmon_dir} are provided, \code{salmon_dir} takes precedence. Default: \code{NULL}.
#'
#' @param tx2gene Either:
#'   - A path to a GFF3 or GFF3.gz file containing transcript-to-gene mapping
#'   - A path to a TSV file with columns 'Transcript' and 'Gene'
#'   - A data.frame with transcript-to-gene mapping
#'
#' @param assay_name Character. Name for the assay (default: 'counts').
#'
#' @param metadata Optional data.frame with sample metadata. Should have sample
#'   names as row names and metadata columns (e.g., sample_type, condition, etc.).
#'
#' @param tpm Optional matrix of transcript-level TPM values. If provided, will be
#'   stored in the SummarizedExperiment. Same dimensions as readcounts required.
#'
#' @param effective_length Optional numeric vector of transcript effective lengths
#'   (e.g., from SALMON). Length should match nrow(readcounts).
#'
#' @param config Optional list of configuration parameters to store in the
#'   TSENATAnalysis object. Useful for tracking analysis parameters.
#'
#' @param skip Logical. If TRUE, allow unmapped transcripts (transcripts not found
#'   in tx2gene mapping) and remove them from analysis. If FALSE (default), stop with
#'   an error when unmapped transcripts are detected. Useful for handling data with
#'   transcript IDs that don't match the annotation file provided.
#'
#' @param verbose Logical. If TRUE, print informative messages during execution
#'   (e.g., Salmon sample discovery, progress on data loading). Default: TRUE.
#'
#' @details
#' When using \code{salmon_dir}, the function automatically:
#' \enumerate{
#'   \item Discovers all Salmon sample folders and quant.sf files
#'   \item Reads transcript counts (NumReads), TPM, and effective_length
#'   \item Extracts sample names from directory structure
#'   \item Creates count matrix ready for analysis
#' }
#'
#' The \code{salmon_dir} parameter provides a convenient alternative to manually
#' constructing the \code{readcounts} matrix, especially useful in Galaxy workflows.
#'
#' @return A \code{TSENATAnalysis} S4 object with:
#'   \item{@se}{The SummarizedExperiment containing transcript counts and metadata}
#'   \item{@config}{Analysis configuration (empty list or user-provided)}
#'   \item{@diversity_results}{Empty list (populated by calculate_diversity_s4())}
#'   \item{@divergence_results}{Empty list (populated by calculate_divergence_s4())}
#'   \item{@lm_results}{Empty list (populated by calculate_lm_interaction_s4())}
#'   \item{@jackknife_results}{Empty list (populated by jackknife functions)}
#'   \item{@plots}{Empty list (populated by plotting functions)}
#'   \item{@metadata}{Metadata with package version and creation timestamp}
#'
#' @details
#' This wrapper combines two steps into one:
#' \enumerate{
#'   \item Call \code{.build_se()} to create a SummarizedExperiment from transcript counts
#'   \item Wrap the result in \code{TSENATAnalysis()} to create the analysis object
#' }
#'
#' The returned object is ready for diversity analysis via \code{calculate_diversity_s4()}.
#'
#' If you need to inspect or filter the SummarizedExperiment before creating the
#' TSENATAnalysis object, call \code{.build_se()} and \code{TSENATAnalysis()} separately.
#'
#' @seealso
#' \code{\link{TSENATAnalysis}} for the S4 class structure
#' \code{\link{calculate_diversity_s4}} for computing Tsallis entropy
#'
#' @examples
#' # Create example transcript count data
#' set.seed(42)
#' n_genes <- 10
#' n_isoforms_per_gene <- 3
#' n_isoforms <- n_genes * n_isoforms_per_gene
#' n_samples <- 10
#'
#' # Generate count matrix
#' counts <- matrix(rpois(n_isoforms * n_samples, lambda = 20),
#'                  nrow = n_isoforms, ncol = n_samples)
#' rownames(counts) <- paste0('TX_', 1:n_isoforms)
#' colnames(counts) <- paste0('Sample_', 1:n_samples)
#'
#' # Create tx2gene mapping
#' tx2gene <- data.frame(
#'   Transcript = rownames(counts),
#'   Gene = rep(paste0('GENE_', 1:n_genes), each = n_isoforms_per_gene))
#'
#' # Create sample metadata
#' metadata <- data.frame(
#'   condition = rep(c('control', 'treatment'), each = 5),
#'   row.names = colnames(counts))
#'
#' # Build analysis object - use NAMED parameters to avoid confusion
#' # Method 1: With explicit tx2gene data.frame (most common)
#' analysis <- build_analysis_s4(
#'   readcounts = counts,
#'   tx2gene = tx2gene,
#'   metadata = metadata)
#'
#' # Verify the analysis object was created
#' analysis
#' print(dim(analysis))
#'
#' \donttest{
#' # Method 2: From Salmon quantification folder
#' # Requires directory structure like:
#' #   salmon_output/
#' #     sample1/quant.sf
#' #     sample2/quant.sf
#' #     ...
#' # 
#' # salmon_dir <- '/path/to/salmon/directory'
#' # 
#' # First create sample metadata matching Salmon sample names
#' # salmon_metadata <- data.frame(
#' #   condition = c('control', 'control', 'treatment', 'treatment'),
#' #   row.names = c('sample1', 'sample2', 'sample3', 'sample4')
#' # )
#' # 
#' # analysis_salmon <- build_analysis_s4(
#' #   salmon_dir = salmon_dir,
#' #   tx2gene = 'annotation.gff3.gz',  # Auto-parsed from GFF3
#' #   metadata = salmon_metadata
#' # )
#' #
#' # Method 3: Hybrid - Salmon counts with manual tx2gene
#' # analysis_hybrid <- build_analysis_s4(
#' #   salmon_dir = salmon_dir,
#' #   tx2gene = tx2gene,  # data.frame instead of file
#' #   metadata = salmon_metadata
#' # )
#' }
#'
#' @export
build_analysis_s4 <- function(readcounts = NULL, salmon_dir = NULL, tx2gene, assay_name = "counts",
    metadata = NULL, tpm = NULL, effective_length = NULL, config = list(), skip = FALSE,
    verbose = TRUE) {
    # Handle salmon_dir parameter - auto-load Salmon quantification data
    if (!is.null(salmon_dir)) {
        # Detect Salmon samples
        salmon_info <- .detect_salmon_samples(salmon_dir)

        if (verbose)
            message("[build_analysis_s4] Found ", salmon_info$count, " Salmon samples")

        # Validate Salmon sample names match metadata (if metadata provided)
        if (!is.null(metadata)) {
            metadata_samples <- rownames(metadata)
            salmon_samples <- salmon_info$sample_names

            # Check if all Salmon samples have corresponding metadata
            missing_in_metadata <- setdiff(salmon_samples, metadata_samples)
            missing_in_salmon <- setdiff(metadata_samples, salmon_samples)

            if (length(missing_in_metadata) > 0 || length(missing_in_salmon) > 0) {
                error_msg <- "[build_analysis_s4] Sample name mismatch between Salmon folder and metadata:\n"

                if (length(missing_in_metadata) > 0) {
                  error_msg <- paste0(error_msg, "  Salmon samples NOT in metadata (",
                    length(missing_in_metadata), "): ", paste(missing_in_metadata,
                      collapse = ", "), "\n")
                }

                if (length(missing_in_salmon) > 0) {
                  error_msg <- paste0(error_msg, "  Metadata samples NOT in Salmon folder (",
                    length(missing_in_salmon), "): ", paste(missing_in_salmon, collapse = ", "),
                    "\n")
                }

                error_msg <- paste0(error_msg, "\n  Suggestion: Ensure sample folder names in Salmon directory exactly match",
                  "\n  the row names of the metadata data.frame (case-sensitive)")

                stop(error_msg)
            }

            if (verbose)
                message("[build_analysis_s4] ✓ Sample names match metadata")
        }

        # Read Salmon quantification files
        salmon_data <- .read_salmon_samples(file_paths = salmon_info$file_paths,
            sample_names = salmon_info$sample_names, include_tpm = is.null(tpm),
            include_eff_length = is.null(effective_length), verbose = verbose)

        # Assign extracted data to parameters
        readcounts <- salmon_data$counts
        if (is.null(tpm))
            tpm <- salmon_data$tpm
        if (is.null(effective_length))
            effective_length <- salmon_data$effective_length
    }

    # Validate readcounts is now available
    if (is.null(readcounts)) {
        stop("[build_analysis_s4] Either 'readcounts' or 'salmon_dir' must be provided\n",
            "  readcounts: matrix/data.frame of transcript counts\n", "  salmon_dir: path to Salmon quantification output directory")
    }
    # Build SummarizedExperiment
    se <- .build_se(readcounts = readcounts, tx2gene = tx2gene, assay_name = assay_name,
        metadata = metadata, tpm = tpm, effective_length = effective_length, skip = skip,
        verbose = verbose)

    # Ensure sample_id column exists in colData (required by TSENATAnalysis)
    # OPTIMIZATION: Only add if not already present
    if (!"sample_id" %in% colnames(SummarizedExperiment::colData(se))) {
        SummarizedExperiment::colData(se)$sample_id <- colnames(se)
    }

    # Store metadata in config for later use (e.g., in
    # calculate_lm_interaction_s4)
    if (!is.null(metadata)) {
        config$metadata <- metadata
    }

    # Wrap in TSENATAnalysis
    analysis <- TSENATAnalysis(se = se, config = config)

    return(analysis)
}
