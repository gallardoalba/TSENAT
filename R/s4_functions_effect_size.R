#' Compute Effect Sizes from Divergence Results (S4 Wrapper)
#'
#' S4 wrapper for  \code{. calculate_effect_sizes()} that 
#' extracts divergence and  LM
#' results directly from a TSENATAnalysis object.
#'
#' @param analysis \code{TSENATAnalysis}.
#'  An S4 object containing divergence results
#'   (from \code{calculate_divergence()}) and LM interaction results
#'   (from \code{calculate_lm()}).
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
#'   \item{Computing}{Effect sizes using base \code{.calculate_effect_sizes()}}
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
#' 
#' # Configure analysis parameters (best practice for reproducibility)
#' config <- TSENAT_config(
#'   sample_col = 'sample',
#'   condition_col = 'condition',
#'   subject_col = 'paired_samples',
#'   paired = TRUE,
#'   control = 'normal',
#'   q = seq(0, 2, by = 0.5)  # Multiple q-values for LM interaction analysis (5 unique: 0, 0.5, 1, 1.5, 2)
#' )
#' analysis <- build_analysis(readcounts = readcounts, tx2gene =
#' gff3_dataset, metadata = metadata_df, config = config,
#'   tpm = tpm, effective_length = effective_length)
#' 
#' analysis <- filter_analysis(analysis, stringency = 'severe')
#' analysis <- calculate_diversity(analysis)
#' analysis <- calculate_divergence(analysis)
#' analysis <- suppressWarnings(calculate_lm(analysis, method = 'gam'))
#'
#' # Compute effect sizes from divergence results
#' analysis <- calculate_effect_sizes(analysis,
#'   significance_threshold = 0.05)
#'
#' # Access results using unified results accessor
#' effect_size_results <- results(analysis, type = 'effect_sizes_divergence')
#'
#' # View structure of results
#' str(effect_size_results, max.level = 1)
#'
#' @seealso
#' \code{\link{calculate_divergence}} for divergence wrapper,
#' \code{\link{calculate_lm}} for LM interaction wrapper
#'
#' @export
#' @importFrom methods is
#' @importFrom utils write.table
calculate_effect_sizes <- function(analysis, significance_threshold = NULL, enrich_per_q_pattern = NULL,
    verbose = NULL, output_file = NULL, ...) {

    verbose <- resolve_slot_param(verbose, getConfig(analysis), "verbose", FALSE)
    .validate_effect_sizes_inputs_s4(analysis)

    significance_threshold <- resolve_slot_param(significance_threshold, getConfig(analysis),
        "significance_threshold", 0.05)
    enrich_per_q_pattern <- resolve_slot_param(enrich_per_q_pattern, getConfig(analysis),
        "enrich_per_q_pattern", TRUE)
    output_file <- resolve_slot_param(output_file, getConfig(analysis), "output_file",
        NULL)

    data_list <- .extract_effect_sizes_data_s4(analysis, verbose)
    divergence_se <- data_list$divergence_se
    lm_res <- data_list$lm_res

    divergence_se <- .add_gene_names_to_divergence_se(divergence_se, analysis, lm_res,
        verbose)

    if (verbose) {
        message("[calculate_effect_sizes] Computing effect sizes...")
    }

    result <- tryCatch({
        .calculate_effect_sizes(lm_res = lm_res, divergence_results_se = divergence_se,
            significance_threshold = significance_threshold, enrich_per_q_pattern = enrich_per_q_pattern,
            verbose = verbose, ...)
    }, error = function(e) {
        stop("calculate_effect_sizes: ", e$message, call. = FALSE)
    })

    analysis <- .store_effect_sizes_results_in_metadata(analysis, result, significance_threshold,
        verbose)
    .save_effect_sizes_output(output_file, result, verbose)

    analysis
}

#' Validate Effect Sizes Computation Inputs
#'
#' @keywords internal
#' @noRd
.validate_effect_sizes_inputs_s4 <- function(analysis) {
    if (!is(analysis, "TSENATAnalysis")) {
        stop("'analysis' must be a TSENATAnalysis object", call. = FALSE)
    }
    # Check slots directly to avoid warnings from accessors
    if (length(analysis@divergence_results) == 0) {
        stop("Divergence results required. Run calculate_divergence() first.", call. = FALSE)
    }
    # Check for LM interaction results (exclude rank_test which belongs to
    # rankResults)
    lm_only <- analysis@lm_results
    if (is.list(lm_only) && "rank_test" %in% names(lm_only)) {
        lm_only$rank_test <- NULL
    }
    if (length(lm_only) == 0) {
        stop("LM results required. Run calculate_lm() first.", call. = FALSE)
    }
}

#' Extract Data for Effect Sizes Computation
#'
#' @keywords internal
#' @noRd
.extract_effect_sizes_data_s4 <- function(analysis, verbose = FALSE) {
    # Access slots directly to avoid accessor method warnings (validation
    # already checked slots exist)
    analysis_divres <- analysis@divergence_results
    divergence_se <- .extract_object_with_fallbacks(analysis_divres, "SummarizedExperiment",
        key_name = "divergence_se", verbose = verbose)

    if (is.null(divergence_se) || !is(divergence_se, "SummarizedExperiment")) {
        stop("Could not extract divergence SummarizedExperiment from divergence results",
            call. = FALSE)
    }

    # Filter out rank_test (rank test results) and access LM results
    analysis_lmres <- analysis@lm_results
    if (is.list(analysis_lmres) && "rank_test" %in% names(analysis_lmres)) {
        analysis_lmres$rank_test <- NULL
    }

    lm_res <- .extract_object_with_fallbacks(analysis_lmres, "data.frame", key_name = "lm_interaction",
        verbose = verbose)

    if (is.null(lm_res) || !is.data.frame(lm_res)) {
        stop("Could not extract LM results data.frame from LM results", call. = FALSE)
    }

    if (verbose) {
        message("[calculate_effect_sizes] Extracted: Divergence SE ", paste(dim(divergence_se),
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
        gene_names_expanded <- .map_tx_to_genes(tx2gene, rownames(divergence_se),
            verbose)
        if (!is.null(gene_names_expanded)) {
            rd$gene_name <- gene_names_expanded
            rowData(divergence_se) <- rd
            if (verbose) {
                message("[calculate_effect_sizes] Added gene_name via tx2gene mapping")
            }
            return(divergence_se)
        }
    }

    # Fallback: direct assignment from lm_res
    if ("gene" %in% colnames(lm_res) && nrow(lm_res) == nrow(divergence_se)) {
        rd$gene_name <- as.character(lm_res$gene)
        rowData(divergence_se) <- rd
        if (verbose) {
            message("[calculate_effect_sizes] Added gene_name via direct assignment")
        }
        return(divergence_se)
    }

    rd_final <- rowData(divergence_se)
    if (!("gene_name" %in% colnames(rd_final)) || any(is.na(rd_final$gene_name))) {
        stop("[calculate_effect_sizes] Failed to add valid gene_name. ", "Ensure tx2gene metadata is properly set.",
            call. = FALSE)
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
        if (tolower(cn) %in% c("transcript", "tx"))
            tx_col <- cn
        if (tolower(cn) %in% c("gene", "gen"))
            gene_col <- cn
    }

    if (is.null(tx_col) || is.null(gene_col))
        return(NULL)

    match_idx <- match(tx_in_divergence, as.character(tx2gene[[tx_col]]))
    gene_names <- as.character(tx2gene[[gene_col]])[match_idx]

    if (all(!is.na(gene_names)))
        gene_names else NULL
}

#' Store Effect Sizes Results in Analysis Metadata
#'
#' @keywords internal
#' @noRd
.store_effect_sizes_results_in_metadata <- function(analysis, result, significance_threshold,
    verbose = FALSE) {
    # Internal function: access object@metadata directly for large result
    # storage getMeta() will filter these out when users call it (returns only
    # essential metadata)

    if (is.null(analysis@metadata))
        analysis@metadata <- list()
    analysis@metadata$effect_sizes_divergence <- result

    current_calls <- analysis@metadata$function_calls %||% character(0)
    analysis@metadata$function_calls <- c(current_calls, paste0("calculate_effect_sizes[threshold=",
        significance_threshold, "]"))

    if (verbose && !is.null(result$interaction_results)) {
        message("[calculate_effect_sizes] Stored results: ", nrow(result$interaction_results),
            " genes")
    }
    analysis
}

#' Save Effect Sizes Output to File
#'
#' @keywords internal
#' @noRd
.save_effect_sizes_output <- function(output_file, result, verbose = FALSE) {
    if (is.null(output_file))
        return()

    if (grepl("\\.tsv$|\\.csv$|\\.txt$", tolower(output_file))) {
        results_df <- as.data.frame(result$interaction_results)
        all_na_cols <- colnames(results_df)[vapply(results_df, function(x) all(is.na(x)),
            FUN.VALUE = logical(1))]
        if (length(all_na_cols) > 0) {
            results_df <- results_df[, !colnames(results_df) %in% all_na_cols]
            if (verbose) {
                message("[calculate_effect_sizes] Removed NA columns: ", paste(all_na_cols,
                  collapse = ", "))
            }
        }
        write.table(results_df, file = output_file, sep = "\t", quote = FALSE, row.names = FALSE)
    } else {
        saveRDS(result, file = output_file)
    }
}


