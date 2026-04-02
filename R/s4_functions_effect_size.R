#' Validate Effect Sizes Computation Inputs
#'
#' @keywords internal
#' @noRd
.validate_effect_sizes_inputs_s4 <- function(analysis) {
    if (!is(analysis, "TSENATAnalysis")) {
        stop("'analysis' must be a TSENATAnalysis object", call. = FALSE)
    }
    if (length(divRes(analysis)) == 0) {
        stop("Divergence results required. Run calculate_divergence_s4() first.", call. = FALSE)
    }
    if (is.null(lmRes(analysis)) || length(lmRes(analysis)) == 0) {
        stop("LM results required. Run calculate_lm_interaction_s4() first.", call. = FALSE)
    }
}

#' Extract Data for Effect Sizes Computation
#'
#' @keywords internal
#' @noRd
.extract_effect_sizes_data_s4 <- function(analysis, verbose = FALSE) {
    analysis_divres <- divRes(analysis)
    divergence_se <- .extract_object_with_fallbacks(analysis_divres, "SummarizedExperiment",
        key_name = "divergence_se", verbose = verbose)

    if (is.null(divergence_se) || !is(divergence_se, "SummarizedExperiment")) {
        stop("Could not extract divergence SummarizedExperiment from divergence results", call. = FALSE)
    }

    analysis_lmres <- lmRes(analysis)
    lm_res <- .extract_object_with_fallbacks(analysis_lmres, "data.frame", key_name = "lm_interaction",
        verbose = verbose)

    if (is.null(lm_res) || !is.data.frame(lm_res)) {
        stop("Could not extract LM results data.frame from LM results", call. = FALSE)
    }

    if (verbose) {
        message("[effect_sizes_divergence_s4] Extracted: Divergence SE ", paste(dim(divergence_se),
            collapse = " x "), ", LM results: ", nrow(lm_res), " genes")
    }

    list(divergence_se = divergence_se, lm_res = lm_res)
}

#' Add Gene Names to Divergence SE
#'
#' @keywords internal
#' @noRd
.add_gene_names_to_divergence_se <- function(divergence_se, analysis, lm_res, verbose = FALSE) {
    rd <- rowData(divergence_se)
    if (!is.null(rd) && "gene_name" %in% colnames(rd)) {
        return(divergence_se)
    }

    if (is.null(rd)) {
        rd <- DataFrame(row.names = rownames(divergence_se))
    }

    # Try tx2gene mapping
    base_se <- getSE(analysis)
    tx2gene <- metadata(base_se)$tx2gene
    if (!is.null(tx2gene) && nrow(tx2gene) > 0) {
        gene_names_expanded <- .map_tx_to_genes(tx2gene, rownames(divergence_se), verbose)
        if (!is.null(gene_names_expanded)) {
            rd$gene_name <- gene_names_expanded
            rowData(divergence_se) <- rd
            if (verbose) {
                message("[effect_sizes_divergence_s4] Added gene_name via tx2gene mapping")
            }
            return(divergence_se)
        }
    }

    # Fallback: direct assignment from lm_res
    if ("gene" %in% colnames(lm_res) && nrow(lm_res) == nrow(divergence_se)) {
        rd$gene_name <- as.character(lm_res$gene)
        rowData(divergence_se) <- rd
        if (verbose) {
            message("[effect_sizes_divergence_s4] Added gene_name via direct assignment")
        }
        return(divergence_se)
    }

    rd_final <- rowData(divergence_se)
    if (!("gene_name" %in% colnames(rd_final)) || any(is.na(rd_final$gene_name))) {
        stop("[effect_sizes_divergence_s4] Failed to add valid gene_name. ",
            "Ensure tx2gene metadata is properly set.", call. = FALSE)
    }
    divergence_se
}

#' Map Transcripts to Genes via tx2gene
#'
#' @keywords internal
#' @noRd
.map_tx_to_genes <- function(tx2gene, tx_in_divergence, verbose = FALSE) {
    tx_col <- gene_col <- NULL
    for (cn in colnames(tx2gene)) {
        if (tolower(cn) %in% c("transcript", "tx")) tx_col <- cn
        if (tolower(cn) %in% c("gene", "gen")) gene_col <- cn
    }

    if (is.null(tx_col) || is.null(gene_col)) return(NULL)

    match_idx <- match(tx_in_divergence, as.character(tx2gene[[tx_col]]))
    gene_names <- as.character(tx2gene[[gene_col]])[match_idx]

    if (all(!is.na(gene_names))) gene_names else NULL
}

#' Store Effect Sizes Results in Analysis Metadata
#'
#' @keywords internal
#' @noRd
.store_effect_sizes_results_in_metadata <- function(analysis, result, significance_threshold,
    verbose = FALSE) {
    current_meta <- getMeta(analysis)
    if (is.null(current_meta)) analysis@metadata <- list()
    analysis@metadata$effect_sizes_divergence <- result

    current_calls <- getMeta(analysis, "function_calls") %||% character(0)
    analysis@metadata$function_calls <- c(current_calls, paste0("effect_sizes_divergence[threshold=",
        significance_threshold, "]"))

    if (verbose && !is.null(result$interaction_results)) {
        message("[effect_sizes_divergence_s4] Stored results: ", nrow(result$interaction_results),
            " genes")
    }
    analysis
}

#' Save Effect Sizes Output to File
#'
#' @keywords internal
#' @noRd
.save_effect_sizes_output <- function(output_file, result, verbose = FALSE) {
    if (is.null(output_file)) return()

    if (grepl("\\.tsv$|\\.csv$|\\.txt$", tolower(output_file))) {
        results_df <- as.data.frame(result$interaction_results)
        all_na_cols <- colnames(results_df)[vapply(results_df, function(x) all(is.na(x)),
            FUN.VALUE = logical(1))]
        if (length(all_na_cols) > 0) {
            results_df <- results_df[, !colnames(results_df) %in% all_na_cols]
            if (verbose) {
                message("[effect_sizes_divergence_s4] Removed NA columns: ", paste(all_na_cols,
                  collapse = ", "))
            }
        }
        write.table(results_df, file = output_file, sep = "\t", quote = FALSE, row.names = FALSE)
    } else {
        saveRDS(result, file = output_file)
    }
}

#' Compute Effect Sizes from Divergence Results (S4 Wrapper)
#'
#' S4 wrapper for  \code{. effect_sizes_divergence()} that 
#' extracts divergence and  LM
#' results directly from a TSENATAnalysis object.
#'
#' @param analysis \code{TSENATAnalysis}.
#'  An S4 object containing divergence results
#'   (from \code{calculate_divergence_s4()}) and LM interaction results
#'   (from \code{calculate_lm_interaction_s4()}).
#'
#' @param significance_threshold \code{numeric}. Adjusted p-value threshold for
#'   filtering significant genes (default: 0.05).
#'
#' @param enrich_per_q_pattern \code{logical}.  If TRUE,
#'  enriches results with  per-q
#'   divergence patterns (default: TRUE).
#'
#' @param verbose \code{logical}.  If TRUE,
#'  print diagnostic messages (default:  TRUE).
#'
#' @param output_file \code{character} or  \code{NULL}.
#'  Optional file path to save results.
#'   Supported formats: .rds (for S4 objects). Default: NULL (no file output).
#'
#' @param ... Additional arguments passed to the base function.
#'
#' @return Modified TSENATAnalysis with effect size results stored via
#'   \code{metadata(analysis)$effect_sizes_divergence}.
#'  Returns the analysis object visibly
#'   to support piping and method chaining.
#'
#' @details
#' **Workflow steps:**
#' \describe{
#'   \item{Validating}{Input analysis object and required results}
#'   \item{Extracting}{Divergence SE and LM results from analysis slots}
#'   \item{Enriching}{Divergence SE with gene names via tx2gene or direct mapping}
#'   \item{Computing}{Effect sizes using base \code{.effect_sizes_divergence()}}
#'   \item{Storing}{Results in metadata with function call tracking}
#' }
#'
#' **Parameter resolution priority** (explicit > metadata > default):
#' \itemize{
#'   \item \code{significance_threshold}:  Uses explicit arg,
#'  else \code{metadata(analysis)$significance_threshold},
#'     else 0.05
#'   \item \code{enrich_per_q_pattern}:  Uses explicit arg,
#'  else \code{metadata(analysis)$enrich_per_q_pattern},
#'     else TRUE
#'   \item \code{verbose}:  Uses explicit arg,
#'  else \code{metadata(analysis)$verbose},  else TRUE
#' }
#'
#' Results are accessed via: \code{metadata(analysis)$effect_sizes_divergence}
#'
#' @examples
#' # Setup: Create test analysis with divergence and LM interaction results
#' data(readcounts)
#' readcounts <- as.matrix(readcounts)
#' mode(readcounts) <- 'numeric'
#' metadata_df <- read.table(
#'   system.file('extdata', 'metadata.tsv', package = 'TSENAT'),
#'   header = TRUE, sep = '\t'
#' )
#' gff3_dataset <- system.file('extdata', 'annotation.gff3.gz', package =
#' 'TSENAT')
#' analysis <- build_analysis_s4(readcounts = readcounts, tx2gene =
#' gff3_dataset, metadata = metadata_df,
#'   tpm = tpm, effective_length = effective_length)
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
#' analysis <- calculate_diversity_s4(analysis, q = c(0.5, 1.0, 1.5),
#' verbose = FALSE)
#' analysis <- calculate_divergence_s4(analysis, q = c(0.5, 1.0, 1.5),
#' verbose = FALSE)
#' analysis <- calculate_lm_interaction_s4(analysis, method = 'gam', verbose
#' = FALSE)
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

    verbose <- resolve_slot_param(verbose, getConfig(analysis), "verbose", FALSE)
    .validate_effect_sizes_inputs_s4(analysis)

    significance_threshold <- resolve_slot_param(significance_threshold, getConfig(analysis),
        "significance_threshold", 0.05)
    enrich_per_q_pattern <- resolve_slot_param(enrich_per_q_pattern, getConfig(analysis),
        "enrich_per_q_pattern", TRUE)
    verbose <- resolve_slot_param(verbose, getConfig(analysis), "verbose", FALSE)

    data_list <- .extract_effect_sizes_data_s4(analysis, verbose)
    divergence_se <- data_list$divergence_se
    lm_res <- data_list$lm_res

    divergence_se <- .add_gene_names_to_divergence_se(divergence_se, analysis, lm_res, verbose)

    if (verbose) {
        message("[effect_sizes_divergence_s4] Computing effect sizes...")
    }

    result <- tryCatch({
        .effect_sizes_divergence(lm_res = lm_res, divergence_results_se = divergence_se,
            significance_threshold = significance_threshold, enrich_per_q_pattern = enrich_per_q_pattern,
            verbose = verbose, ...)
    }, error = function(e) {
        stop("effect_sizes_divergence: ", e$message, call. = FALSE)
    })

    analysis <- .store_effect_sizes_results_in_metadata(analysis, result, significance_threshold,
        verbose)
    .save_effect_sizes_output(output_file, result, verbose)

    analysis
}