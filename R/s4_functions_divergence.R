#' Calculate divergence metrics and store in TSENATAnalysis
#'
#' Wrapper around .calculate_divergence() that manages TSENATAnalysis object.
#' Calculates Tsallis divergence between experimental conditions for transcripts
#' across multiple q-values to detect condition-specific transcript remodeling.
#'
#' Key Features:
#' \itemize{
#'   \item Multi-q analysis: Divergence computed across full q-spectrum simultaneously
#'   \item Bootstrap confidence intervals: Quantify uncertainty in divergence estimates
#'   \item Multiple testing correction: Hochberg, Benjamini-Yekutieli, or no correction
#'   \item Paired designs: Supports paired/repeated measures via subject_col parameter
#'   \item Effect size reporting: Log-fold-change and confidence intervals per gene
#'   \item Flexible control group: Compare any condition vs. any other condition
#' }
#'
#' @param analysis \code{TSENATAnalysis} object.
#' @param q \code{numeric}. Q-value for divergence.
#' If NULL, uses first q_value from @config$q_values if available, else
#' defaults to 1.0.
#' @param verbose \code{logical}. Print progress messages. Default: TRUE.
#' @param nthreads \code{numeric} or  \code{NULL}.  Number of CPU threads for 
#' parallel processing.
#'   If NULL, reads from \code{@config$nthreads} (or defaults to 1).
#' @param output_file \code{character} or  \code{NULL}.
#'  Optional file path to save results.
#' Supported formats: .rds (for S4 objects), .tsv, .csv, .txt (for tables).
#' Default: NULL (no file output).
#' @param control_group \code{character} or  \code{NULL}.
#'  Control group identifier for  divergence comparison.
#'   If NULL, reads from \code{@config$control_group} if available.
#' @param paired \code{logical}. Whether to use paired design. Default: FALSE.
#'   If not specified, reads from \code{@config$paired} if available.
#' @param method \code{character} or  \code{NULL}.  Statistical method for 
#' divergence calculation.
#'   If NULL, reads from \code{@config$method} if available.
#' @param bootstrap \code{logical}.
#'  Whether to compute bootstrap confidence intervals.  Default:  FALSE.
#'   If not specified, reads from \code{@config$bootstrap} if available.
#' @param nboot \code{numeric} or  \code{NULL}.
#'  Number of bootstrap replicates.  Default:  NULL.
#'   If NULL, reads from \code{@config$nboot} if available.
#' @param seed \code{numeric} or  \code{NULL}.  Random seed for 
#' reproducibility.  Default:  NULL.
#'   If NULL, reads from \code{@config$seed} if available.
#' @param progress \code{logical}.  Show progress bar during computation.
#'  Default:  FALSE.
#' @param ... Additional arguments passed to the base divergence function.
#'
#' @return Modified TSENATAnalysis with divergence metrics in
#'   \code{@divergence_results} (stored as list of data.frames or matrices).
#'
#' @details
#' **Mathematical Background:**
#' Tsallis divergence D_q between two probability distributions:
#' \preformatted{
#'   D_q(P||Q) = (log_2(N) - entropy_q(P) + entropy_q(Q)) / (q - 1)
#' }
#' Measures how much transcript composition changes from control to condition.
#' Values near 0: Similar isoform composition; Large positive values: Major change.
#'
#' **Example Use Case:**
#' Control sample: All reads from dominant isoform (low entropy)\cr
#' Tumor sample: Reads spread across multiple isoforms (high entropy)\cr
#' Result: Large divergence indicating condition-specific isoform switching.\cr
#'
#' **Parameter Resolution:**
#' Parameters are resolved in priority order:
#' 1. Explicit arguments passed to function
#' 2. Values from analysis@config (if present)
#' 3. Function defaults
#'
#' Requires diversity results from calculate_diversity_s4() as prerequisite.
#'
#' @examples
#' # Load example data (matching TSENAT.Rmd workflow)
#' data(readcounts)
#' readcounts <- as.matrix(readcounts)
#' mode(readcounts) <- 'numeric'
#' metadata_df <- read.table(
#'   system.file('extdata', 'metadata.tsv', package = 'TSENAT'),
#'   header = TRUE, sep = '\t'
#' )
#' gff3_dataset <- system.file('extdata', 'annotation.gff3.gz', package =
#' 'TSENAT')
#' 
#' # Build analysis from vignette data and create manageable subset
#' analysis <- build_analysis_s4(readcounts = readcounts, tx2gene =
#' gff3_dataset, metadata = metadata_df,
#'   tpm = tpm, effective_length = effective_length)
#' # Use 200+ genes to ensure diversity filtering doesn't remove all genes
#' analysis <- filter_analysis_s4(analysis, min_samples = 1, subset_n_genes
#' = 200)
#' 
#' # Compute diversity first (required for divergence)
#' analysis <- calculate_diversity_s4(analysis, q = c(0.5, 1.0, 1.5),
#' verbose = FALSE)
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
    
    # Step 1: Validate input
    .validate_divergence_input(analysis)

    # Step 2: Resolve and process parameters
    params <- .resolve_divergence_parameters(q, control_group, method, nthreads, nboot, seed,
        paired, bootstrap, analysis)

    # Step 3: Build arguments for computation
    args <- .build_divergence_args(analysis, params, verbose, progress, ...)

    # Step 4: Run divergence calculation
    result <- tryCatch({
        do.call(.calculate_divergence, args)
    }, error = function(e) {
        stop("Divergence calculation failed:\n", e$message, call. = FALSE)
    })

    # Step 5: Store results
    analysis <- .store_divergence_results(analysis, result)

    # Step 6: Track metadata
    analysis@metadata$function_calls <- c(analysis@metadata$function_calls, paste0("calculate_divergence[q=",
        params$q, "]"))

    # Step 7: Write output file if requested
    if (!is.null(output_file)) {
        .write_divergence_output(analysis, output_file, verbose)

        # Write bootstrap results if applicable
        if (isTRUE(params$bootstrap)) {
            .write_divergence_bootstrap_output(analysis, output_file, verbose)
        }
    }

    analysis
}

# INTERNAL HELPERS FOR DIVERGENCE CALCULATION
# ============================================================================

#' Validate divergence calculation inputs
#'
#' @noRd
#' @noRd
.validate_divergence_input <- function(analysis) {
    if (!is(analysis, "TSENATAnalysis")) {
        stop("'analysis' must be a TSENATAnalysis object", call. = FALSE)
    }

    if (nrow(analysis@se) == 0) {
        stop("SummarizedExperiment in @se is empty", call. = FALSE)
    }

    if (length(analysis@diversity_results) == 0) {
        stop("Diversity results required. Run calculate_diversity_s4() first.", call. = FALSE)
    }

    invisible(NULL)
}

#' Resolve and process divergence parameters
#'
#' @noRd
#' @noRd
.resolve_divergence_parameters <- function(q, control_group, method, nthreads, nboot, seed,
    paired, bootstrap, analysis) {
    # Extract parameters using utility functions
    q <- resolve_slot_param(q, analysis@config, "q_values", 1)
    control_group <- resolve_slot_param(control_group, analysis@config, "control_group", NULL)
    method <- resolve_slot_param(method, analysis@config, "method", "percentile")
    nthreads <- resolve_slot_param(nthreads, analysis@config, "nthreads", 1)
    nboot <- resolve_slot_param(nboot, analysis@config, "nboot", NULL)
    seed <- resolve_slot_param(seed, analysis@config, "seed", NULL)

    # Replace q=0 with q=0.01 (q=0 always returns 0, which is uninformative)
    if (is.vector(q)) {
        q[q == 0] <- 0.01
    } else if (is.numeric(q) && q == 0) {
        q <- 0.01
    }

    # Ensure paired is logical
    if (!isTRUE(paired) && !isFALSE(paired)) {
        paired <- if ("paired" %in% names(analysis@config)) analysis@config$paired else FALSE
        if (!is.logical(paired) || is.na(paired)) paired <- FALSE
    }

    # Ensure bootstrap is logical
    if (!isTRUE(bootstrap) && !isFALSE(bootstrap)) {
        bootstrap <- if ("bootstrap" %in% names(analysis@config)) analysis@config$bootstrap else FALSE
        if (!is.logical(bootstrap) || is.na(bootstrap)) bootstrap <- FALSE
    }

    # Sanitize method
    if (!is.null(method)) {
        method <- as.character(method[1])
    }
    if (is.null(method) || is.na(method)) {
        method <- "percentile"
    }

    list(q = q, control_group = control_group, method = method, nthreads = nthreads,
        nboot = nboot, seed = seed, paired = paired, bootstrap = bootstrap)
}

#' Build arguments list for divergence calculation
#'
#' @noRd
#' @noRd
.build_divergence_args <- function(analysis, params, verbose, progress, ...) {
    args <- list(se = analysis@se, q = params$q, verbose = verbose, progress = progress)

    if (!is.null(params$control_group)) {
        args$control_group <- params$control_group
    }

    if (isTRUE(params$paired)) {
        args$paired <- params$paired
    }

    if (!is.null(params$method)) {
        args$method <- params$method
    }

    if (!is.null(params$nthreads)) {
        args$nthreads <- params$nthreads
    }

    if (isTRUE(params$bootstrap)) {
        args$bootstrap <- params$bootstrap
        if (!is.null(params$nboot)) {
            args$nboot <- params$nboot
        }
        if (!is.null(params$seed)) {
            args$seed <- params$seed
        }
    }

    # Merge with additional args (may override defaults)
    c(args, list(...))
}

#' Store divergence results in analysis object
#'
#' @noRd
#' @noRd
.store_divergence_results <- function(analysis, result) {
    if (is.null(result)) {
        warning("[calculate_divergence_s4] Result is NULL. Check .calculate_divergence() output.",
            call. = FALSE)
        analysis@divergence_results <- list()
    } else if (is(result, "SummarizedExperiment")) {
        analysis@divergence_results <- list(divergence_se = result)
    } else if (is.list(result)) {
        analysis@divergence_results <- result
    } else if (is.data.frame(result) || is.matrix(result)) {
        analysis@divergence_results <- list(main = result)
    } else {
        warning("[calculate_divergence_s4] Unexpected result type: ", class(result),
            ". Wrapping in list.", call. = FALSE)
        analysis@divergence_results <- list(result = result)
    }

    analysis
}

#' Extract divergence data for file writing
#'
#' @noRd
#' @noRd
.extract_divergence_write_data <- function(div_list) {
    if (length(div_list) == 0) {
        return(NULL)
    }

    if ("divergence_se" %in% names(div_list)) {
        return(as.data.frame(SummarizedExperiment::assay(div_list$divergence_se)))
    }

    if (is(div_list[[1]], "SummarizedExperiment")) {
        data_list <- lapply(div_list, function(se) {
            as.data.frame(SummarizedExperiment::assay(se, 1))
        })
        return(do.call(cbind, data_list))
    }

    as.data.frame(div_list)
}

#' Write divergence results to file
#'
#' @noRd
#' @noRd
.write_divergence_output <- function(analysis, output_file, verbose) {
    # Create output directory if needed
    output_dir <- dirname(output_file)
    if (output_dir != "." && !dir.exists(output_dir)) {
        dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
    }

    if (grepl("\\.tsv$|\\.csv$|\\.txt$", tolower(output_file))) {
        tryCatch({
            write_data <- .extract_divergence_write_data(analysis@divergence_results)
            if (!is.null(write_data)) {
                sep <- if (grepl("\\.csv$", tolower(output_file))) "," else "\t"
                write.table(write_data, file = output_file, sep = sep, quote = FALSE, row.names = TRUE)
            } else {
                warning("[calculate_divergence_s4] No divergence results to save", call. = FALSE)
            }
        }, error = function(e) {
            warning("[calculate_divergence_s4] Could not write divergence results: ",
                conditionMessage(e), call. = FALSE)
        })
    } else {
        saveRDS(analysis, file = output_file)
    }

    invisible(NULL)
}

#' Extract SummarizedExperiment from divergence results
#'
#' @noRd
#' @noRd
.extract_divergence_se <- function(div_list) {
    if ("divergence_se" %in% names(div_list)) {
        return(div_list$divergence_se)
    }

    if (is(div_list[[1]], "SummarizedExperiment")) {
        return(div_list[[1]])
    }

    NULL
}

#' Build bootstrap column list from rowData
#'
#' @noRd
#' @noRd
.build_bootstrap_cols <- function(rd) {
    cols <- c("gene_name", grep("^estimate_q|^lower_ci_q|^upper_ci_q|^ci_width_q|^nboot_q|^method_q",
        colnames(rd), value = TRUE))

    if ("computation_time_sec" %in% colnames(rd)) {
        cols <- c(cols, "computation_time_sec")
    }
    if ("error" %in% colnames(rd)) {
        cols <- c(cols, "error")
    }

    cols[cols %in% colnames(rd)]
}

#' Write bootstrap divergence results to file
#'
#' @noRd
#' @noRd
.write_divergence_bootstrap_output <- function(analysis, output_file, verbose) {
    # Generate bootstrap filename (preserving original case)
    bootstrap_file <- gsub("\\.(tsv|csv|txt)$", "_bootstrap.\\1", output_file,
        ignore.case = TRUE)
    if (basename(bootstrap_file) == basename(output_file)) {
        bootstrap_file <- sub("(\\.[^.]+)$", "_bootstrap\\1", output_file)
    }

    tryCatch({
        if (isTRUE(verbose)) {
            message("[calculate_divergence_s4] Attempting to save bootstrap results to: ",
                bootstrap_file)
        }

        # Extract SE and build column list
        se <- .extract_divergence_se(analysis@divergence_results)
        if (is.null(se) || !is(se, "SummarizedExperiment")) {
            return(invisible(NULL))
        }

        rd <- SummarizedExperiment::rowData(se)
        bootstrap_cols <- .build_bootstrap_cols(rd)

        # Write if columns available
        if (length(bootstrap_cols) > 1) {
            bootstrap_data <- as.data.frame(rd[, bootstrap_cols])
            sep <- if (grepl("\\.csv$", tolower(bootstrap_file))) "," else "\t"
            
            # Ensure parent directory exists
            parent_dir <- dirname(bootstrap_file)
            if (!dir.exists(parent_dir)) {
                dir.create(parent_dir, recursive = TRUE, showWarnings = FALSE)
            }
            
            write.table(bootstrap_data, file = bootstrap_file, sep = sep, quote = FALSE, row.names = FALSE)

            if (isTRUE(verbose)) {
                message("[calculate_divergence_s4] Saved bootstrap results to: ", bootstrap_file)
            }
        }
    }, error = function(e) {
        warning("[calculate_divergence_s4] Could not write bootstrap results: ",
            conditionMessage(e), call. = FALSE)
    })

    invisible(NULL)
}


