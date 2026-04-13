# ============================================================================
# JACKKNIFE WRAPPER
# ============================================================================

#' Jackknife resampling with confidence intervals
#'
#' @param analysis \code{TSENATAnalysis} object.
#' @param q \code{numeric}. Q-value(s) for jackknife estimation.
#'   If NULL, reads from \code{analysis@config$q}. 
#'   Default: c(0, 0.5, 1, 1.5, 2) (matches \code{calculate_jis} spectrum).
#' @param norm \code{logical}.  Normalization flag.  Default:
#'  NULL (uses @config$norm or  TRUE).
#' @param log_base \code{numeric}.  Logarithm base for  entropy normalization.
#'  Default:  NULL (uses e).
#' @param top_n \code{numeric}.  Number of top outlier samples to report.
#'  Default:  5.
#' @param pseudocount \code{numeric}.  Pseudocount value for 
#' count regularization.  Default:  NULL (uses @config$pseudocount or  0).
#' @param verbose \code{logical}.  Print jackknife results summary.  Default:
#'  FALSE.
#' @param nthreads \code{numeric} or  \code{NULL}.  Number of CPU threads for 
#' parallel processing.
#'   If NULL, reads from \code{@config$nthreads} (or defaults to 1).
#'   If > 1 and multiple q-values provided, uses parallel PSOCK cluster.
#' @param output_file \code{character} or  \code{NULL}.
#'  Optional file path to save results.
#' Supported formats: .tsv, .csv, .txt (for jackknife results table with
#' estimates, influence, outliers),
#'   .rds (for entire S4 object). Default: NULL (no file output).
#' @param ... Additional arguments passed to the base function.
#'
#' @return Modified TSENATAnalysis with jackknife results in @jackknife_results.
#'
#' @details
#' Requires diversity results to exist first. Will error if
#' \code{calculate_diversity()} has not been run.
#'
#' **Parameter Priority Resolution:**
#' \describe{
#'   \item{nthreads}{Priority: explicit > \code{@config$nthreads} > 1}
#' }
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
#' # TPM and effective_length REQUIRED for filter_analysis()
#' tpm <- matrix(runif(nrow(readcounts) * ncol(readcounts), 0.1, 10),
#'               nrow = nrow(readcounts), ncol = ncol(readcounts),
#'               dimnames = dimnames(readcounts))
#' effective_length <- matrix(100, nrow = nrow(readcounts), ncol = ncol(readcounts))
#' 
#' # Create config (metadata passed as explicit parameter to build_analysis)
#' config <- TSENAT_config(
#'   sample_col = 'sample',
#'   condition_col = 'condition',
#'   q = seq(0, 2, by = 0.05),
#'   paired = FALSE
#' )
#' 
#' # Build analysis from vignette data - metadata as explicit parameter
#' analysis <- build_analysis(
#'   readcounts = readcounts,
#'   metadata = metadata_df,
#'   tx2gene = gff3_dataset,
#'   config = config,
#'   tpm = tpm,
#'   effective_length = effective_length
#' )
#'
#' # Filter low-abundance genes (required for reliable jackknife estimates)
#' analysis <- filter_analysis(analysis, stringency = 'severe')
#' 
#' # Compute diversity first (required for jackknife)
#' analysis <- calculate_diversity(analysis, q = c(0.5, 1.0, 1.5))
#' 
#' # Run jackknife estimation
#' analysis <- calculate_jeo(analysis, q = c(0.5, 1.0, 1.5))
#' # Check jackknife results using unified accessor
#' jackknife_res <- results(analysis, type = 'jackknife')
#' if (!is.null(jackknife_res)) names(jackknife_res)
#'
#' @export
#' @importFrom utils write.table
calculate_jeo <- function(analysis, q = NULL, norm = NULL, log_base = NULL, top_n = NULL,
    verbose = NULL, nthreads = NULL, pseudocount = NULL, output_file = NULL, ...) {
    if (!is(analysis, "TSENATAnalysis")) {
        stop("'analysis' must be a TSENATAnalysis object", call. = FALSE)
    }

    # Check prerequisite: diversity must be calculated
    if (length(analysis@diversity_results) == 0) {
        stop("Diversity results required. Run calculate_diversity() first.", call. = FALSE)
    }

    # PARAMETER EXTRACTION using utility function
    q_source <- "explicit"
    q <- resolve_slot_param(q, analysis@config, "q", NULL)
    if (is.null(q)) {
        q <- c(0, 0.5, 1, 1.5, 2)
        q_source <- "default"
    } else if (!is.null(analysis@config$q) && identical(q, analysis@config$q)) {
        q_source <- "config"
    }

    # Show message about q-value(s) being used
    verbose_resolved <- resolve_slot_param(verbose, analysis@config, "verbose", FALSE)
    if (verbose_resolved && q_source != "explicit") {
        if (length(q) == 1) {
            message("[calculate_jeo] Using q = ", formatC(q, format = "f", digits = 3),
                " (", q_source, ")")
        } else {
            message("[calculate_jeo] Using q spectrum: ", paste(formatC(q, format = "f",
                digits = 3), collapse = ", "), " (", q_source, ")")
        }
    }

    norm <- resolve_slot_param(norm, analysis@config, "norm", TRUE)
    log_base <- resolve_slot_param(log_base, analysis@config, "log_base", exp(1))
    top_n <- resolve_slot_param(top_n, analysis@config, "top_n", 5)
    nthreads <- resolve_slot_param(nthreads, analysis@config, "nthreads", 1)
    pseudocount <- resolve_slot_param(pseudocount, analysis@config, "pseudocount",
        0)
    verbose <- resolve_slot_param(verbose, analysis@config, "verbose", FALSE)
    output_file <- resolve_slot_param(output_file, analysis@config, "output_file",
        NULL)

    # Ensure q is numeric
    if (!is.numeric(q)) {
        stop("'q' must be numeric", call. = FALSE)
    }

    for (q_val in q) {
        # Verify diversity was computed for this q-value (prerequisite for
        # jackknife)
        div_key <- paste0("q_", formatC(q_val, format = "f", digits = 3))
        if (!(div_key %in% names(analysis@diversity_results))) {
            stop("Diversity not calculated for q=", q_val, ". Run calculate_diversity(analysis, q=",
                q_val, ") first.", call. = FALSE)
        }

        # Run jackknife - pass the COUNTS (not diversity values!) to jackknife
        # function Jackknife stability analysis requires raw count data, not
        # pre-computed diversity
        tryCatch({
            # Extract counts matrix from SummarizedExperiment
            counts_matrix <- SummarizedExperiment::assay(analysis@se, "counts")

            result <- .calculate_jeo(x = counts_matrix, q = q_val, norm = norm, log_base = log_base,
                top_n = top_n, pseudocount = pseudocount, verbose = verbose, nthreads = nthreads,
                ...)

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
                    message("[calculate_jeo] Results saved to ", output_file)
                  }
                }
            }, error = function(e) {
                warning("[calculate_jeo] Could not write jackknife results to file: ",
                  conditionMessage(e), call. = FALSE)
            })
        } else {
            # Default to RDS for S4 object
            saveRDS(analysis, file = output_file)
            if (verbose) {
                message("[calculate_jeo] Analysis object saved to ", output_file)
            }
        }
    }
    analysis
}
