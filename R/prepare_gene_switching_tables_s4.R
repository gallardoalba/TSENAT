#' Prepare Gene Switching Tables from TSENATAnalysis Object
#'
#' S4 wrapper for \code{prepare_gene_switching_tables()} that extracts results
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
#' \dontrun{
#'   # After running full analysis pipeline
#'   analysis <- calculate_lm_interaction_s4(analysis, ...)
#'   analysis <- jackknife_isoform_switching_s4(analysis, ...)
#'
#'   # Prepare tables using S4 wrapper
#'   tables <- prepare_gene_switching_tables_s4(analysis, n_top_genes = 20)
#'
#'   # Access individual components
#'   summary_table <- tables$summary_table
#'   gene_tables <- tables$transcript_tables
#' }
#'
#' @export
#' @importFrom methods is
prepare_gene_switching_tables_s4 <- function(
    analysis,
    n_top_genes = NULL,
    n_transcripts_per_gene = 10,
    verbose = FALSE) {
  
  # Validation
  if (!is(analysis, "TSENATAnalysis")) {
    stop("analysis must be a TSENATAnalysis object")
  }
  
  # Extract LM results
  if (verbose) cat("Extracting LM results from analysis object...\n")
  
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
  
  if (verbose) cat("  ✓ Extracted LM results with", nrow(lm_res), "genes\n")
  
  # Extract jackknife/switching results
  if (verbose) cat("Extracting jackknife switching results from analysis object...\n")
  
  jackknife_results_list <- analysis@jackknife_results
  if (is.null(jackknife_results_list) || length(jackknife_results_list) == 0) {
    stop("No jackknife results found in analysis@jackknife_results. Run jackknife_isoform_switching_s4() first.")
  }
  
  # Check if results are stored under "multi_q" key (jackknife_isoform_switching_multiq class)
  if ("multi_q" %in% names(jackknife_results_list)) {
    multi_q_object <- jackknife_results_list[["multi_q"]]
    
    # If it's a tsenat_isoform_switching_multiq object, use it directly
    if (inherits(multi_q_object, "tsenat_isoform_switching_multiq")) {
      multi_q_results <- multi_q_object
      if (verbose) {
        cat("  ✓ Found multi-q results under 'multi_q' key with", 
            length(multi_q_results), "q-values\n")
      }
    } else {
      stop("Element at jackknife_results$multi_q is not a tsenat_isoform_switching_multiq object")
    }
  } else {
    # Fallback: look for q-keyed results directly
    # Extract only the q-keyed results (filter by pattern "q_X_XX")
    q_key_pattern <- "^q_[0-9]+_[0-9]{2}$"
    q_keyed_results <- jackknife_results_list[grep(q_key_pattern, names(jackknife_results_list))]
    
    if (length(q_keyed_results) == 0) {
      stop("No q-keyed jackknife results found in analysis@jackknife_results. ",
           "Expected keys in format 'q_X_XX' (e.g., 'q_0_01', 'q_1_00') or 'multi_q'.")
    }
    
    # Wrap q-keyed results as a multi_q object for consistency
    multi_q_results <- q_keyed_results
    if (verbose) {
      cat("  ✓ Found", length(multi_q_results), "q-keyed results\n")
    }
  }
  
  if (verbose) cat("  ✓ Extracted jackknife results with", length(multi_q_results), "q-values\n")
  
  # Call base function with extracted parameters
  if (verbose) cat("Calling prepare_gene_switching_tables()...\n")
  
  result <- prepare_gene_switching_tables(
    lm_res = lm_res,
    multi_q_results = multi_q_results,
    n_top_genes = n_top_genes,
    n_transcripts_per_gene = n_transcripts_per_gene,
    verbose = verbose
  )
  
  if (verbose) cat("✓ Gene switching tables prepared successfully\n")
  
  # Track function call in metadata
  analysis@metadata$function_calls <- c(analysis@metadata$function_calls, 
                                        "prepare_gene_switching_tables_s4")
  analysis@metadata$function_timestamps <- c(analysis@metadata$function_timestamps,
                                             as.character(Sys.time()))
  
  return(result)
}
