
#' Calculate LM interactions and store in TSENATAnalysis
#'
#' @param analysis \code{TSENATAnalysis} object.
#' @param fdr_threshold \code{numeric}. FDR cutoff for significance.
#'   Default: 0.05.
#' @param formula \code{formula} or NULL. Reserved for future use.
#' @param condition_col \code{character} or  \code{NULL}.
#'  Column name in colData identifying sample conditions.
#'   If NULL,  reads from \code{@config$condition_col} or 
#' auto-detects common column names.
#' @param method \code{character}.  Statistical method (e. g. ,  'lmm',
#'  'gam',  'gee').
#'   If NULL, uses method from @config$method or defaults to 'lmm'.
#' @param paired \code{logical}. Whether to use paired design. Default: FALSE.
#'   If not specified, reads from \code{@config$paired} if available.
#' @param subject_col \code{character} or  \code{NULL}.
#'  Column name identifying subject IDs for  paired designs.
#'   If NULL, reads from \code{@config$subject_col} if available.
#' @param nthreads \code{numeric} or  \code{NULL}.  Number of CPU threads for 
#' parallel processing.
#'   If NULL,  reads from \code{@config$nthreads} (or  defaults to NULL,
#'  letting base function decide).
#' @param multicorr \code{character} or  \code{NULL}.
#'  Multiple comparison correction method.
#'   Options: 'hochberg', 'westfall-young', 'benjamini-yekutieli'.
#'   If NULL, uses method from @config or base function defaults.
#' @param corstr \code{character} or  \code{NULL}.  Correlation structure for 
#' GEE models.
#'   Options: 'ar1', 'exchangeable', 'independence'.
#'   If NULL, uses method from @config or base function defaults.
#' @param pcorr \code{character} or \code{NULL}. P-value correction method.
#'   Default: 'BH' (Benjamini-Hochberg).
#'   If NULL, reads from \code{@config$pcorr} if available.
#' @param verbose \code{logical}. Print progress messages. Default: FALSE.
#' @param return_model_data \code{logical}.  Return model data for 
#' visualization.  Default:  TRUE.
#' @param output_file \code{character} or  \code{NULL}.
#'  Optional file path to save results.
#' Supported formats: .rds (for S4 objects), .tsv, .csv, .txt (for tables).
#' Default: NULL (no file output).
#' @param ... Additional arguments passed to the base LM function,
#' including: pvalue, min_obs, assay_name, bias_correction, regularization,
#' storey, wy_randomizations, adaptive_knots, etc.
#'
#' @return Modified TSENATAnalysis with results in @lm_results$lm_interaction.
#'
#' @details
#' Extracts diversity results from @diversity_results (prerequisite),
#' combines across q-values into single SummarizedExperiment,
#' then runs \code{.calculate_lm()}.
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
#' # Note: requires lme4 package for LMM fitting; uses synthetic data
#' set.seed(42)
#' 
#' # Create transcript-level counts with biological signal
#' # Note: Use adequate complexity (transcripts/genes, samples, expression) 
#' # to avoid filtering away all genes during diversity computation
#' n_genes <- 50
#' n_transcripts_per_gene <- 30
#' n_transcripts <- n_genes * n_transcripts_per_gene
#' n_samples <- 16  # 8 subjects x 2 conditions (paired design)
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
#' analysis <- calculate_diversity(
#'   analysis, 
#'   q = c(0.5, 1.0, 1.5, 2.0, 2.5)
#' )
#' 
#' # Calculate q x condition interactions using GAM
#' analysis <- suppressWarnings(calculate_lm(
#'   analysis,
#'   condition_col = 'condition',
#'   method = 'gam'
#' ))
#' 
#' # View top interaction results using unified accessor (first 3 genes)
#' res <- results(analysis, type = "lm")
#' if (!is.null(res)) head(res, 3)
#'
#' @export
#' @importFrom utils write.table
calculate_lm <- function(analysis, fdr_threshold = NULL, formula = NULL,
    condition_col = NULL, method = "gam", paired = NULL, subject_col = NULL, nthreads = NULL,
    multicorr = NULL, corstr = NULL, pcorr = NULL, verbose = NULL, return_model_data = NULL,
    output_file = NULL, ...) {
    # Validate input
    if (!is(analysis, "TSENATAnalysis")) {
        stop("'analysis' must be a TSENATAnalysis object", call. = FALSE)
    }

    # Clear cache for fresh calculations
    .clear_lm_helper_cache()

    # Check prerequisites
    if (length(analysis@diversity_results) == 0) {
        stop("Diversity results required. Run calculate_diversity() first.", call. = FALSE)
    }

    # Sync colData from diversity results
    analysis <- .sync_coldata_from_diversity(analysis, verbose = verbose)

    # Resolve parameters from config first
    fdr_threshold <- resolve_slot_param(fdr_threshold, analysis@config, "fdr_threshold",
        NULL)
    formula <- resolve_slot_param(formula, analysis@config, "formula", NULL)
    output_file <- resolve_slot_param(output_file, analysis@config, "output_file",
        NULL)
    verbose <- resolve_slot_param(verbose, analysis@config, "verbose", FALSE)
    paired <- resolve_slot_param(paired, analysis@config, "paired", FALSE)
    return_model_data <- resolve_slot_param(return_model_data, analysis@config, "return_model_data",
        TRUE)

    # Extract and resolve remaining parameters Note: .extract_lm_params()
    # receives resolved paired value
    params <- .extract_lm_params(analysis, condition_col = condition_col, method = method,
        subject_col = subject_col, nthreads = nthreads, multicorr = multicorr, corstr = corstr,
        pcorr = pcorr, paired = paired, verbose = verbose)

    # Combine diversity results across q-values
    diversity_combined <- .combine_diversity_results_for_lm(analysis@diversity_results)

    # REQUIREMENT: Check that we have at least 5 unique q-values
    # ARIMA(1,1,0) differencing removes 1 observation per subject, leaving (n_q - 1) unique values
    # With 5 q-values: after ARIMA -> 4 unique q-values
    # This provides sufficient degrees of freedom for GAM spline fitting (k=3 or k=4 works with 4 unique values)
    q_values <- sort(as.numeric(unique(sub(".*q=", "", colnames(diversity_combined)))))
    if (length(q_values) < 5) {
        stop(sprintf("[calculate_lm] At least 5 unique q-values are required for interaction analysis. Current data has only %d unique q-value(s). Ensure diversity_results contains >=5 distinct q values. After ARIMA(1,1,0) differencing, this leaves sufficient degrees of freedom for GAM fitting.", 
                     length(q_values)), call. = FALSE)
    }

    # Build arguments for LM calculation
    args <- .build_lm_args(diversity_combined, params, return_model_data = return_model_data,
        verbose = verbose, ...)

    # Run LM analysis Phase 15: Catch errors gracefully - return empty results
    # instead of crashing
    result <- tryCatch({
        do.call(.calculate_lm, args)
    }, error = function(e) {
        # Return empty data.frame on error instead of stopping workflow
        warning("lm_interaction calculation failed:\n", conditionMessage(e), call. = FALSE)
        data.frame()
    })

    # Validate and extract results
    extracted <- .validate_and_extract_lm_result(result)
    if (nrow(extracted$results) == 0) {
        # Return early with empty results
        analysis@lm_results <- list(lm_interaction = data.frame())
        return(analysis)
    }

    # Store results in analysis object
    analysis <- .store_lm_results_in_analysis(analysis, extracted$results, extracted$model_data)

    # Save output if requested
    if (!is.null(output_file) && is.data.frame(extracted$results)) {
        save_analysis_output(extracted$results, output_file, object = analysis, verbose = verbose,
            func_name = "calculate_lm")
    }

    analysis
}

# Helper: Sync colData from diversity results to analysis@se @param analysis
# TSENATAnalysis object @param verbose Logical: print messages @keywords
# internal
.sync_coldata_from_diversity <- function(analysis, verbose = FALSE) {
    if (length(analysis@diversity_results) == 0) {
        return(analysis)
    }

    first_diversity_se <- analysis@diversity_results[[1]]
    if (!is(first_diversity_se, "SummarizedExperiment") || ncol(first_diversity_se) ==
        0) {
        return(analysis)
    }

    diversity_coldata <- SummarizedExperiment::colData(first_diversity_se)
    if (is.null(diversity_coldata) || nrow(diversity_coldata) != ncol(analysis@se)) {
        return(analysis)
    }

    if (ncol(diversity_coldata) > 0) {
        SummarizedExperiment::colData(analysis@se) <- diversity_coldata
    }

    analysis
}


# Helper: Extract and resolve LM interaction parameters from config @param
# analysis TSENATAnalysis object
.extract_lm_params <- function(analysis, condition_col = NULL, method = NULL, subject_col = NULL,
    nthreads = NULL, multicorr = NULL, corstr = NULL, pcorr = NULL, paired = FALSE,
    verbose = FALSE) {
    # Auto-detect condition_col if not provided Note: Always pass verbose=TRUE
    # for condition_col to ensure users are aware of auto-detection
    if (is.null(condition_col)) {
        cd_cols <- colnames(colData(analysis@se))
        condition_col <- auto_detect_column(cd_cols, config_list = analysis@config,
            config_key = "condition_col", priority_candidates = c("condition", "sample_type",
                "group", "treatment"), default_fallback = NULL, verbose = FALSE,
            param_name = "condition_col")
    }

    # Resolve remaining parameters using centralized handler
    method <- resolve_slot_param(method, analysis@config, "method", "gam")
    subject_col <- resolve_slot_param(subject_col, analysis@config, "subject_col",
        NULL)
    multicorr <- resolve_slot_param(multicorr, analysis@config, "multicorr", NULL)
    corstr <- resolve_slot_param(corstr, analysis@config, "corstr", NULL)
    pcorr <- resolve_slot_param(pcorr, analysis@config, "pcorr", "BH")
    nthreads <- resolve_slot_param(nthreads, analysis@config, "nthreads", NULL)

    # Note: 'paired' is already resolved by calling function
    # (calculate_lm) to avoid duplicate resolution. Use as-is.

    # Log condition_col info
    if (is.null(condition_col)) {
        cd_cols <- colnames(colData(analysis@se))
        if (length(cd_cols) > 0) {
            message("condition_col not specified. Available columns: ", paste(cd_cols,
                collapse = ", "))
        } else {
            message("condition_col not specified and colData is empty. ", "Will be determined by .calculate_lm().")
        }
    }

    list(condition_col = condition_col, method = method, subject_col = subject_col,
        nthreads = nthreads, multicorr = multicorr, corstr = corstr, pcorr = pcorr,
        paired = paired)
}


# Helper: Combine per-q diversity results into single SE
.combine_diversity_results_for_lm <- function(diversity_results) {
    tryCatch({
        assay_list <- list()
        coldata_list <- list()
        rowdata_first <- NULL
        q_keys <- sort(names(diversity_results))

        for (key in q_keys) {
            se <- diversity_results[[key]]
            if (!is(se, "SummarizedExperiment")) {
                stop("[combine_diversity_results_for_lm] Result ", key, " is not a SummarizedExperiment",
                  call. = FALSE)
            }

            # Extract q-value from key (e.g., 'q_0.100' -> '0.100')
            q_val_str <- sub("^q_", "", key)

            # Validate assays exist
            if (length(SummarizedExperiment::assays(se)) == 0) {
                stop("[combine_diversity_results_for_lm] SE ", key, " has no assays",
                  call. = FALSE)
            }

            # Get assay and add q-value suffix to column names
            assay_matrix <- as.matrix(SummarizedExperiment::assay(se))
            assay_colnames_with_q <- paste0(colnames(assay_matrix), "_q=", q_val_str)
            colnames(assay_matrix) <- assay_colnames_with_q
            assay_list[[key]] <- assay_matrix

            # Get colData with matching rownames
            coldata <- SummarizedExperiment::colData(se)
            rownames(coldata) <- assay_colnames_with_q
            coldata_list[[key]] <- coldata

            # Save rowData from first SE
            if (is.null(rowdata_first)) {
                rowdata_first <- SummarizedExperiment::rowData(se)
            }
        }

        # Combine assays and colData
        combined_assay <- do.call(cbind, assay_list)
        combined_coldata <- do.call(rbind, coldata_list)
        colnames(combined_assay) <- rownames(combined_coldata)

        # Create combined SE
        SummarizedExperiment::SummarizedExperiment(assays = list(diversity = combined_assay),
            colData = combined_coldata, rowData = rowdata_first)
    }, error = function(e) {
        stop("[combine_diversity_results_for_lm] Failed: ", conditionMessage(e),
            call. = FALSE)
    })
}


# Helper: Build argument list for calculate_lm_interaction
.build_lm_args <- function(diversity_se, params, return_model_data = TRUE, verbose = FALSE,
    ...) {
    args <- list(se = diversity_se)

    # Add resolved parameters if non-NULL
    if (!is.null(params$method)) {
        args$method <- params$method
    }
    if (!is.null(params$condition_col)) {
        args$condition_col <- params$condition_col
    }
    if (params$paired) {
        args$paired <- params$paired
    }
    if (!is.null(params$subject_col)) {
        args$subject_col <- params$subject_col
    }
    if (!is.null(params$nthreads)) {
        args$nthreads <- params$nthreads
    }
    if (!is.null(params$multicorr)) {
        args$multicorr <- params$multicorr
    }
    if (!is.null(params$corstr)) {
        args$corstr <- params$corstr
    }
    if (!is.null(params$pcorr)) {
        args$pcorr <- params$pcorr
    }

    args$return_model_data <- return_model_data
    if (verbose) {
        args$verbose <- verbose
    }

    # Merge with additional args
    c(args, list(...))
}


# Helper: Validate and extract LM results from raw output
.validate_and_extract_lm_result <- function(result) {
    lm_results_df <- result
    model_data <- NULL

    # Extract components if result is list with $results and $model_data
    if (is.list(result) && "results" %in% names(result)) {
        lm_results_df <- result$results
        model_data <- result$model_data
    }

    # Validate structure
    if (!is.data.frame(lm_results_df)) {
        stop("[validate_and_extract_lm_result] Result must be data.frame, got: ",
            class(lm_results_df), call. = FALSE)
    }

    # Check for empty results
    if (nrow(lm_results_df) == 0 || ncol(lm_results_df) == 0) {
        warning("[validate_and_extract_lm_result] Result is empty (", nrow(lm_results_df),
            " rows, ", ncol(lm_results_df), " columns). ", "This can occur with: low sample counts per condition, ",
            "insufficient signal, or model convergence issues.", call. = FALSE)
        return(list(results = data.frame(), model_data = NULL))
    }

    # Validate required columns
    required_cols <- c("gene", "adj_p_interaction")
    missing_cols <- setdiff(required_cols, colnames(lm_results_df))
    if (length(missing_cols) > 0) {
        stop("[validate_and_extract_lm_result] Missing columns: ", paste(missing_cols,
            collapse = ", "), ". Available: ", paste(colnames(lm_results_df), collapse = ", "),
            call. = FALSE)
    }

    list(results = lm_results_df, model_data = model_data)
}


# Helper: Store LM results in analysis object
.store_lm_results_in_analysis <- function(analysis, lm_results_df, model_data = NULL) {
    if (is.list(analysis@lm_results) && "lm_interaction" %in% names(analysis@lm_results)) {
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

    # Track function call
    analysis@metadata$function_calls <- c(analysis@metadata$function_calls, "calculate_lm")

    analysis
}


