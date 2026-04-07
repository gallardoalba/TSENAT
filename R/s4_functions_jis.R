#' Jackknife isoform switching analysis on TSENATAnalysis object
#'
#' Wrapper around .jackknife_isoform_switching() that manages TSENATAnalysis
#' object. Identifies transcripts with significant isoform switching patterns
#' using jackknife resampling across samples to detect influential isoforms.
#'
#' Key Features:
#' \itemize{
#'   \item Jackknife resampling: Robust outlier detection across all samples
#'   \item Delta-influence metric: Measures how much each isoform drives phenotype
#'   \item Confidence intervals: Bootstrap-based uncertainty quantification
#'   \item Multi-q analysis: Tests across full q-spectrum simultaneously
#'   \item LM filtering: Optional restriction to genes with significant interactions
#'   \item Paired designs: Supports repeated measures/longitudinal data
#' }
#'
#' @param analysis \code{TSENATAnalysis} object containing:
#'   \itemize{
#'     \item \code{@se}: SummarizedExperiment with count data
#'     \item \code{@config}: Configuration metadata
#'   }
#'
#' @param condition_col \code{character}. Column name in colData(se) specifying 
#' group assignments (default: 'sample_type'). If NULL, attempts
#' auto-detection.
#'
#' @param subject_col \code{character}.  Optional column for 
#' paired/repeated measures design.
#'   If provided, enables paired analysis. Default: NULL (unpaired).
#'
#' @param gene_col \code{character}.  Column name in rowData(se) or 
#' metadata identifying genes.
#'   Default: 'gene'.
#'
#' @param isoform_col \code{character}.  Column name in rowData(se) or 
#' metadata identifying 
#'   isoforms/transcripts. Default: 'transcript' or 'isoform'.
#'
#' @param q \code{numeric}.  Tsallis entropy parameter(s) to analyze.
#'  Can be single value 
#'   or vector for multi-q analysis (default: c(0, 0.5, 1, 1.5, 2)).
#'   If NULL, uses @config$q_values.
#'
#' @param norm \code{logical}. Whether to use normalized diversity values 
#'   (default: TRUE).
#'
#' @param log_base \code{numeric}.  Logarithm base for  entropy calculations.
#'  Default:  NULL (uses e).
#'
#' @param pseudocount \code{numeric}.  Pseudocount value for 
#' count regularization.  Default:  NULL (uses @config$pseudocount or  0).
#'
#' @param threshold \code{numeric}.  Percentile threshold for 
#' detecting transcript switching 
#' (default: 90). Transcripts with delta_influence >= threshold percentile
#' are classified
#'   as 'switching'.
#'
#' @param n_bootstrap \code{integer}.  Number of bootstrap resamples for 
#' confidence 
#'   intervals (default: 1000).
#'
#' @param lm_results \code{data. frame}.
#'  Optional LM interaction results to filter genes.
#'   If provided, only genes in lm_results are analyzed.
#'
#' @param lm_p_threshold \code{numeric}.  P-value threshold for 
#' filtering genes from 
#'   lm_results (default: 0.05).
#'
#' @param use_lm_fdr \code{logical}.  If TRUE,
#'  uses adjusted p-values from lm_results 
#'   (default: TRUE).
#'
#' @param verbose \code{logical}. Print progress messages (default: FALSE).
#' @param output_file \code{character} or  \code{NULL}.
#'  Optional file path to save results.
#'   Supported formats: .tsv, .csv, .txt (for tables), or .rds (for S4 objects).
#'   When text format is specified, generates TWO files:
#'   \itemize{
#'     \item \code{output_file}:  Gene-level summary (one row per gene with 
#' switching statistics)
#'     \item \code{output_file_transcripts. ext}:
#'  Transcript-level details (one row per transcript with  p-values and  FDR)
#'   }
#'   For .rds format, saves only the full analysis object.
#'   Default: NULL (no file output).
#' @param ... Additional arguments for future extensibility.
#'
#' @details
#' **Parameter Resolution from Config**
#'
#' The following parameters are resolved using a three-level priority system:
#' \enumerate{
#'   \item User-provided argument (if not NULL)
#'   \item Value from \code{analysis@config} (if key exists)
#'   \item Function default value
#' }
#'
#' Affected parameters:
#' \itemize{
#'   \item \code{q}: Uses config's \code{q_values} if available, else c(0, 0.5, 1, 1.5, 2)
#'   \item \code{n_bootstrap}: Uses config's \code{n_bootstrap} if available, else 1000
#'   \item \code{threshold}: Uses config's \code{threshold} if available, else 90
#'   \item \code{lm_p_threshold}: Uses config's \code{lm_p_threshold} if available, else 0.05
#' }
#'
#' This allows setting defaults once in the config and reusing across multiple analyses.
#'
#' @return \code{TSENATAnalysis} object with 
#' jackknife results stored in \code{@jackknife_results}
#' slot. Results are keyed by q-value (e.g., 'q_1.00'). For multi-q
#' analysis, multiple
#'   calls will accumulate results in the slot.
#'
#'   The analysis object is returned visibly to support method chaining:
#'   \preformatted{
#'     analysis <- jackknife_isoform_switching_s4(analysis, q = 0.5)
#'     analysis <- jackknife_isoform_switching_s4(analysis, q = 1.0)
#'   }
#'
#' @details
#' **Mathematical Background:**
#' Delta-influence measures how much removing each sample changes entropy:
#' \preformatted{
#'   Delta = H_q(leave-one-out) - H_q(original)
#' }
#' High |Delta| for specific isoforms indicates those isoforms drive differences.
#' Identifies 'outlier samples' where isoforms contribute unusually much.
#'
#' **Example Use Case:**
#' Sample shows high Delta for isoform X → X has outsized importance in that sample\cr
#' Classifying transcripts as 'switching' if top percentile (e.g., 90th) Delta\cr
#' Reveals condition-specific isoforms crucial for phenotype determination.\cr
#'
#' **Automatic Setup:**
#' This wrapper automatically:
#' 1. Extracts SummarizedExperiment from \code{@se} slot
#' 2. Detects condition_col, gene_col, isoform_col from colData/rowData or
#' \code{@config}
#' 3. Calls \code{.jackknife_isoform_switching()} with extracted parameters
#'
#' **Parameter Auto-Detection:**
#' \enumerate{
#'   \item \code{condition_col}:  Uses explicit parameter,  then @config,
#'  then 'sample_type'
#'   \item \code{gene_col}:  Uses explicit parameter,  then looks for 
#' 'gene' or  'Gene'
#'   \item \code{isoform_col}:  Uses explicit parameter,  then looks for 
#' 'transcript',  'isoform',  or  'Isoform'
#' }
#'
#' @seealso \code{jackknife_isoform_switching} for 
#' the underlying implementation,
#' \code{\link{TSENATAnalysis}} for object structure.
#'
#' @examples
#' data(readcounts)
#' metadata_df <- read.table(system.file('extdata', 'metadata.tsv', package
#' = 'TSENAT'),
#'                           header = TRUE, sep = '\t')
#' gff3_file <- system.file('extdata', 'annotation.gff3.gz', package = 'TSENAT')
#' config <- tsenat_config(sample_col = 'sample', condition_col = 'condition')
#' analysis <- build_analysis_s4(readcounts = readcounts, tx2gene =
#' gff3_file, metadata = metadata_df, config = config,
#'                              tpm = tpm,
#'  effective_length = effective_length)
#' analysis <- filter_analysis_s4(analysis, min_samples = 1, subset_n_genes
#' = 20, subset_n_samples = 8)
#' analysis <- calculate_diversity_s4(analysis, q = 1)
#' 
#' @export
jackknife_isoform_switching_s4 <- function(analysis, condition_col = NULL, subject_col = NULL,
    gene_col = NULL, isoform_col = NULL, q = c(0, 0.5, 1, 1.5, 2), norm = NULL, log_base = NULL, threshold = 90,
    n_bootstrap = 1000, pseudocount = NULL, lm_results = NULL, lm_p_threshold = 0.05,
    use_lm_fdr = TRUE, output_file = NULL, verbose = FALSE, ...) {
    # Validate input and extract SummarizedExperiment
    se <- .validate_jis_input(analysis)

    # Auto-detect and validate column names
    col_info <- .detect_jis_columns(analysis, se, condition_col, gene_col, isoform_col,
        verbose)

    condition_col <- col_info$condition_col
    gene_col <- col_info$gene_col
    isoform_col <- col_info$isoform_col

    # Extract or validate LM results
    lm_results <- .extract_lm_results(analysis, lm_results, verbose)

    # Resolve and validate all parameters (including threshold and lm_p_threshold)
    params <- .resolve_and_validate_jis_params(q, norm, log_base, pseudocount, n_bootstrap,
        threshold, lm_p_threshold, analysis, verbose)

    # Call base jackknife function
    result <- tryCatch({
        .jackknife_isoform_switching(se = se, condition_col = condition_col, subject_col = subject_col,
            gene_col = gene_col, isoform_col = isoform_col, q = params$q, norm = params$norm,
            log_base = params$log_base, threshold = params$threshold, n_bootstrap = params$n_bootstrap,
            pseudocount = params$pseudocount, verbose = verbose, lm_results = lm_results,
            lm_p_threshold = params$lm_p_threshold, use_lm_fdr = use_lm_fdr)
    }, error = function(e) {
        stop("[jackknife_isoform_switching_s4] Jackknife analysis failed:\n", conditionMessage(e),
            call. = FALSE)
    })

    # Store results and update metadata
    analysis <- .store_jis_results(analysis, result, params$q_vals, condition_col,
        verbose)

    # Save results to file if requested
    if (!is.null(output_file)) {
        .save_jis_output(output_file, result, analysis, verbose)
    }

    # Return modified analysis object (visibly for method chaining as
    # documented)
    analysis
}

# ============================================================================
# HELPER FUNCTIONS FOR JACKKNIFE ISOFORM SWITCHING
# ============================================================================

#' Validate input analysis and extract SummarizedExperiment
#'
#' @param analysis TSENATAnalysis object
#'
#' @return SummarizedExperiment object
#'
#' @noRd
.validate_jis_input <- function(analysis) {
    if (!is(analysis, "TSENATAnalysis")) {
        stop("'analysis' must be a TSENATAnalysis object", call. = FALSE)
    }

    se <- analysis@se
    if (!inherits(se, "SummarizedExperiment")) {
        stop("[jackknife_isoform_switching_s4] @se must be a SummarizedExperiment object",
            call. = FALSE)
    }

    se
}

#' Auto-detect column names from rowData and colData
#'
#' @param analysis TSENATAnalysis object
#' @param se SummarizedExperiment
#' @param condition_col User-provided condition column
#' @param gene_col User-provided gene column
#' @param isoform_col User-provided isoform column
#' @param verbose Logical, print messages
#'
#' @return List with detected columns: condition_col, gene_col, isoform_col
#'
#' @noRd
.detect_jis_columns <- function(analysis, se, condition_col, gene_col, isoform_col,
    verbose) {
    cd_cols <- colnames(colData(se))

    # Use explicit condition_col if provided, otherwise auto-detect Note:
    # Always pass verbose=TRUE for condition_col to ensure users are aware of
    # auto-detection
    if (is.null(condition_col)) {
        condition_col <- auto_detect_column(cd_cols, config_list = analysis@config,
            config_key = "condition_col", priority_candidates = c("sample_type",
                "condition", "group", "sample_group"), default_fallback = if (length(cd_cols) >
                0)
                cd_cols[1] else NULL, verbose = FALSE, param_name = "condition_col")
    } else {
        # Validate explicit condition_col exists
        if (!condition_col %in% cd_cols) {
            stop("[jackknife_isoform_switching_s4] condition_col '", condition_col,
                "' not found in colData.\n", "  Available columns: ", paste(cd_cols,
                  collapse = ", "), call. = FALSE)
        }
    }

    if (is.null(condition_col)) {
        stop("[jackknife_isoform_switching_s4] Cannot auto-detect condition_col.\n",
            "  Available colData columns: ", paste(cd_cols, collapse = ", "), "\n\n",
            "SOLUTION: Set @config$condition_col or pass explicit parameter\n", call. = FALSE)
    }

    # Extract rowData columns
    rd <- if (!is.null(rowData(se)) && nrow(rowData(se)) > 0)
        rowData(se) else NULL
    rd_cols <- if (!is.null(rd))
        colnames(rd) else character(0)

    # Use explicit gene_col if provided, otherwise auto-detect
    if (is.null(gene_col)) {
        gene_col <- auto_detect_column(rd_cols, config_list = analysis@config, config_key = "gene_col",
            priority_candidates = c("gene_id", "gene", "Gene", "gene_name"), default_fallback = "gene",
            verbose = verbose, param_name = "gene_col")
    } else {
        # Validate explicit gene_col exists if rowData is present
        if (length(rd_cols) > 0 && !gene_col %in% rd_cols) {
            stop("[jackknife_isoform_switching_s4] gene_col '", gene_col, "' not found in rowData.\n",
                "  Available columns: ", paste(rd_cols, collapse = ", "), call. = FALSE)
        }
    }

    # Use explicit isoform_col if provided, otherwise auto-detect
    if (is.null(isoform_col)) {
        isoform_col <- auto_detect_column(rd_cols, config_list = analysis@config,
            config_key = "isoform_col", priority_candidates = c("transcript_id",
                "transcript", "isoform", "Isoform", "tx_id"), default_fallback = "transcript",
            verbose = verbose, param_name = "isoform_col")
    } else {
        # Validate explicit isoform_col exists if rowData is present
        if (length(rd_cols) > 0 && !isoform_col %in% rd_cols) {
            stop("[jackknife_isoform_switching_s4] isoform_col '", isoform_col, "' not found in rowData.\n",
                "  Available columns: ", paste(rd_cols, collapse = ", "), call. = FALSE)
        }
    }

    list(condition_col = condition_col, gene_col = gene_col, isoform_col = isoform_col)
}

#' Extract LM results from analysis object or use provided results
#'
#' @param analysis TSENATAnalysis object
#' @param lm_results Optional user-provided LM results
#' @param verbose Logical, print messages
#'
#' @return LM results data frame or NULL
#'
#' @noRd
.extract_lm_results <- function(analysis, lm_results, verbose) {
    if (is.null(lm_results) && !is.null(analysis@lm_results)) {
        if ("lm_interaction" %in% names(analysis@lm_results)) {
            lm_results <- analysis@lm_results$lm_interaction
            if (verbose) {
                message("[jackknife_isoform_switching_s4] Using LM interaction results")
            }
        }
    }
    lm_results
}

#' Validate q-values are available in diversity results
#'
#' @param analysis TSENATAnalysis object
#' @param q Target q-values
#' @param verbose Logical, print messages
#'
#' @noRd
.validate_diversity_q_values <- function(analysis, q, verbose) {
    if (is.null(analysis@diversity_results) || length(analysis@diversity_results) ==
        0) {
        return(invisible(NULL))
    }

    available_q_keys <- names(analysis@diversity_results)
    # Extract numeric q values from keys like 'q_1_00' -> 1.00
    available_q <- as.numeric(gsub("_", ".", sub("^q_", "", available_q_keys)))
    available_q <- sort(unique(available_q[!is.na(available_q)]))

    q_vals <- if (is.numeric(q))
        q else c(q)
    
    # Use intersection of requested and available q-values
    # This allows jackknife to work with any set of q-values from diversity
    usable_q <- intersect(q_vals, available_q)
    
    if (length(usable_q) == 0) {
        warning("[jackknife_isoform_switching_s4] No matching q-values found.\n",
                "  Requested: ", paste(q_vals, collapse = ", "), "\n",
                "  Available: ", paste(available_q, collapse = ", "), call. = FALSE)
    } else if (verbose) {
        message("[jackknife_isoform_switching_s4] Using q-values: ", paste(usable_q, collapse = ", "))
    }

    invisible(NULL)
}

#' Resolve and validate all numerical parameters
#'
#' @param q Q-values to analyze
#' @param norm Normalization flag
#' @param log_base Logarithm base
#' @param pseudocount Pseudocount for regularization
#' @param n_bootstrap Number of bootstrap resamples
#' @param analysis TSENATAnalysis object for config access
#' @param verbose Logical, print messages
#'
#' @return List with resolved parameters and q_vals (numeric vector)
#'
#' @noRd
.resolve_and_validate_jis_params <- function(q, norm = NULL, log_base = NULL, pseudocount = NULL, n_bootstrap = 1000,
    threshold = NULL, lm_p_threshold = NULL, analysis, verbose = FALSE) {
    # Resolve q: use config if available, otherwise use the provided value (which has function default)
    if (is.null(q) && "q_values" %in% names(analysis@config)) {
        q <- analysis@config$q_values
    }
    
    # Convert q to numeric vector
    q_vals <- if (is.numeric(q))
        q else as.numeric(q)

    # Resolve all parameters from config using standard resolver
    # Defaults match TSENAT.Rmd vignette usage
    norm <- resolve_slot_param(norm, analysis@config, "norm", TRUE)
    log_base <- resolve_slot_param(log_base, analysis@config, "log_base", exp(1))
    pseudocount <- resolve_slot_param(pseudocount, analysis@config, "pseudocount", 0)
    n_bootstrap <- resolve_slot_param(n_bootstrap, analysis@config, "n_bootstrap", 1000)
    threshold <- resolve_slot_param(threshold, analysis@config, "threshold", 90)
    lm_p_threshold <- resolve_slot_param(lm_p_threshold, analysis@config, "lm_p_threshold", 0.05)

    # Validate n_bootstrap
    if (!is.numeric(n_bootstrap) || length(n_bootstrap) != 1 || n_bootstrap < 1) {
        stop("'n_bootstrap' must be a positive integer", call. = FALSE)
    }

    if (n_bootstrap < 50) {
        warning("n_bootstrap = ", n_bootstrap, " is less than recommended minimum 50",
            call. = FALSE)
    }

    # Validate threshold
    if (!is.numeric(threshold) || length(threshold) != 1) {
        stop("'threshold' must be a single numeric value", call. = FALSE)
    }

    # Validate lm_p_threshold
    if (!is.numeric(lm_p_threshold) || length(lm_p_threshold) != 1) {
        stop("'lm_p_threshold' must be a single numeric value", call. = FALSE)
    }

    list(q = q, q_vals = q_vals, norm = norm, log_base = log_base, pseudocount = pseudocount,
        n_bootstrap = n_bootstrap, threshold = threshold, lm_p_threshold = lm_p_threshold)
}

#' Store jackknife results in analysis object
#'
#' @param analysis TSENATAnalysis object
#' @param result Jackknife result object
#' @param q_vals Numeric vector of q-values
#' @param condition_col Column name used
#' @param verbose Logical, print messages
#'
#' @return Modified TSENATAnalysis object
#'
#' @noRd
.store_jis_results <- function(analysis, result, q_vals, condition_col, verbose) {
    if (inherits(result, "tsenat_isoform_switching_multiq")) {
        # Multi-q result
        for (q_key in names(result)) {
            analysis@jackknife_results[[q_key]] <- result[[q_key]]
        }
        analysis@jackknife_results[["multi_q"]] <- result
        if (verbose) {
            message("[jackknife_isoform_switching_s4] Stored multi-q results")
        }
    } else {
        # Single or vector q-values
        q_keys <- paste0("q_", gsub("\\.", "_", sprintf("%.2f", q_vals)))

        for (i in seq_along(q_keys)) {
            q_key <- q_keys[i]
            if (is.list(result) && q_key %in% names(result)) {
                analysis@jackknife_results[[q_key]] <- result[[q_key]]
            } else if (length(q_vals) == 1) {
                analysis@jackknife_results[[q_key]] <- result
            } else {
                warning("[jackknife_isoform_switching_s4] Result for q=", q_vals[i],
                  " (key: ", q_key, ") not found", call. = FALSE)
            }

            if (verbose) {
                message("[jackknife_isoform_switching_s4] Stored results for ", q_key)
            }
        }
    }

    # Update metadata
    if (is.list(analysis@metadata)) {
        call_str <- sprintf("jackknife_isoform_switching_s4[q=%s, condition_col=%s]",
            paste(q_vals, collapse = ","), condition_col)
        analysis@metadata$function_calls <- c(analysis@metadata$function_calls, call_str)
        analysis@metadata$function_timestamps <- c(analysis@metadata$function_timestamps,
            as.character(Sys.time()))
    }

    analysis
}

#' Save jackknife results to file
#'
#' @param output_file File path for output
#' @param result Jackknife result object
#' @param analysis TSENATAnalysis object (for RDS output)
#' @param verbose Logical, print messages
#'
#' @noRd
.save_jis_output <- function(output_file, result, analysis, verbose) {
    output_dir <- dirname(output_file)
    if (output_dir != "." && !dir.exists(output_dir)) {
        dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
    }

    if (grepl("\\.tsv$|\\.csv$|\\.txt$", tolower(output_file))) {
        .write_jis_tables(output_file, result, verbose)
    } else {
        tryCatch({
            saveRDS(analysis, file = output_file)
            if (verbose) {
                message("[jackknife_isoform_switching_s4] Saved to ", output_file)
            }
        }, error = function(e) {
            warning("[jackknife_isoform_switching_s4] Failed to save: ", conditionMessage(e),
                call. = FALSE)
        })
    }
}

#' Write jackknife tables to text files
#'
#' @param output_file File path for output
#' @param result Jackknife result object
#' @param verbose Logical, print messages
#'
#' @noRd
.write_jis_tables <- function(output_file, result, verbose) {
    tryCatch({
        output_ext <- sub("^.*\\.", ".", tolower(output_file))
        output_base <- sub(paste0(output_ext, "$"), "", output_file)
        transcript_file <- paste0(output_base, "_transcripts", output_ext)

        # Extract gene-level table
        gene_data <- extract_multiq_table(result, is_multiq = inherits(result, "tsenat_isoform_switching_multiq"),
            extract_fn = function(res, q_key) {
                if (!is.null(res$summary_table))
                  res$summary_table else NULL
            }, q_value_col = "q_value")

        # Extract transcript-level table
        transcript_data <- extract_multiq_table(result, is_multiq = inherits(result,
            "tsenat_isoform_switching_multiq"), extract_fn = function(res, q_key) {
            if (!is.null(res$all_transcript_stats))
                res$all_transcript_stats else NULL
        }, q_value_col = "q_value")

        # Write gene-level table
        if (!is.null(gene_data)) {
            write.table(gene_data, file = output_file, sep = "\t", quote = FALSE,
                row.names = FALSE)
        }

        # Write transcript-level table
        if (!is.null(transcript_data)) {
            write.table(transcript_data, file = transcript_file, sep = "\t", quote = FALSE,
                row.names = FALSE)
        }
    }, error = function(e) {
        warning("[jackknife_isoform_switching_s4] Could not write results: ", conditionMessage(e),
            call. = FALSE)
    })
}
