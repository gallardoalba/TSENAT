#' Jackknife Isoform Switching Detection
#'
#' Detects isoform/transcript switches between conditions using
#' condition-stratified
#' jackknife analysis on Tsallis entropy. Identifies which transcripts change
#' importance between conditions (e.g., condition A vs B).
#'
#' @param se A SummarizedExperiment object with transcript-level counts.
#' @param condition_col Character: column name in colData for condition labels.
#' @param subject_col Character: optional column name for paired design
#' (subject/individual IDs).
#' @param gene_col Character: column name in rowData for gene IDs.
#' @param isoform_col Character: column name in rowData for
#' transcript/isoform IDs.
#' @param q Numeric: Tsallis entropy order (default 1 = Shannon entropy).
#' Can be a vector
#' for multi-q analysis (e.g., q = c(0.5, 1.0, 1.5, 2.0)); results will be
#' nested by q value.
#' @param norm Logical: normalize entropy to [0,1]? (default TRUE).
#' @param log_base Numeric: log base for entropy (default e).
#' @param pseudocount Numeric: pseudocount to add (default 0).
#' @param threshold Numeric: percentile for outlier detection on influences
#' (default 90).
#' @param n_bootstrap Numeric: number of bootstrap resamples (default 1000).
#' @param verbose Logical: print results and verbose output? (default TRUE).
#' @param lm_results Data frame: results from .calculate_lm()
#' with 'gene' column.
#' Can contain either gene names or gene IDs; function automatically maps
#' names to IDs
#' using rowData(se). Include 'p_interaction' and/or 'adj_p_interaction'
#' columns for
#' filtering genes by significance. When provided, only genes passing
#' lm_p_threshold are
#'   analyzed; all matching genes are included (top_n parameter removed).
#' @param lm_p_threshold Numeric: p-value threshold for LM gene filtering
#' (default 0.05).
#' @param use_lm_fdr Logical: use adjusted p-values from LM results if
#' available (default TRUE).
#'
#' @return If q is a single value, returns a list of class
#' tsenat_isoform_switching with:
#'   \describe{
#'     \item{results_per_gene}{named list of per-gene results}
#'     \item{summary_table}{data.frame with per-gene summary}
#'     \item{all_transcript_stats}{data.frame with all transcript statistics}
#'     \item{gene_names}{character vector of analyzed genes}
#'     \item{conditions}{character vector of the two conditions compared}
#'     \item{metadata}{list with analysis metadata}
#'   }
#'   
#' If q is a vector, returns a list of class tsenat_isoform_switching_multiq
#' where each
#' element is a complete tsenat_isoform_switching result for that q value.
#' Keys are
#'   formatted as 'q_X_XX' for ease of iteration (e.g., q_0_01, q_1_00, q_2_00).
#'
#' @section Sample Metadata Parameters (Unified Naming Convention):
#' TSENAT functions use consistent parameter names for sample grouping and
#' subject identification:
#' \itemize{
#'   \item{\code{condition_col}: Character string specifying the colData column 
#' containing sample group/condition labels (e.g., 'Normal', 'Tumor',
#' 'control', 'treatment').
#' Default: 'condition'. Required for identifying the two conditions to
#' compare.}
#'   \item{\code{subject_col}:  For paired/blocked designs,
#'  character string specifying 
#' the colData column with subject/individual/patient identifiers. Default:
#' NULL.}
#' }
#' All functions use \code{SummarizedExperiment: :
#' colData()} as the single source of truth
#' for sample metadata. This eliminates parameter fragmentation and improves
#' API discoverability
#' across the TSENAT package.
#'
#' @noRd
# MAIN FUNCTION - Refactored to ~45 lines using helper functions
.jackknife_isoform_switching <- function(se = NULL, condition_col = "condition",
    subject_col = NULL, gene_col = NULL, isoform_col = NULL, q = 1, norm = TRUE,
    log_base = exp(1), pseudocount = 0, threshold = 90, n_bootstrap = 1000, verbose = TRUE,
    lm_results = NULL, lm_p_threshold = 0.05, use_lm_fdr = TRUE) {
    # 1. Validate input
    conditions <- .jis_validate_input(se, condition_col, gene_col, isoform_col)

    # 2. Handle multiple q values
    if (is.numeric(q) && length(q) > 1) {
        q_params <- list(condition_col = condition_col, subject_col = subject_col,
            gene_col = gene_col, isoform_col = isoform_col, norm = norm, log_base = log_base,
            pseudocount = pseudocount, threshold = threshold, n_bootstrap = n_bootstrap,
            lm_results = lm_results, lm_p_threshold = lm_p_threshold, use_lm_fdr = use_lm_fdr)
        return(.jis_handle_multi_q(se, q, q_params, verbose))
    }

    # 3. Setup
    conditions <- sort(conditions)
    paired_info <- .setup_paired_design_jis(se, subject_col, condition_col)
    gene_id_to_name <- .build_gene_id_mapping(se, gene_col)
    gene_ids <- unique(rowData(se)[[gene_col]])

    # 4. Setup LM filtering
    lm_setup <- .jis_setup_lm_filtering(se, lm_results, lm_p_threshold, use_lm_fdr,
        gene_ids, gene_col)
    gene_ids <- lm_setup$filtered_genes
    lm_gene_mapping <- lm_setup$lm_gene_mapping

    # 5. Process genes
    gene_results <- .jis_process_all_genes(se, gene_ids, gene_col, isoform_col, condition_col,
        conditions, paired_info, q, norm, log_base, pseudocount, n_bootstrap, lm_gene_mapping)

    # 6. Build results
    summary_results <- .jis_build_summary_results(se, gene_results$results_per_gene,
        gene_results$all_pvalues, gene_col, lm_gene_mapping)

    # 7. Compile final result
    n_fdr_sig_total <- if (nrow(summary_results$all_transcript_stats) > 0)
        sum(summary_results$all_transcript_stats$fdr < 0.05, na.rm = TRUE) else 0

    metadata <- list(q = q, is_paired = paired_info$is_paired, subject_col = if (paired_info$is_paired) paired_info$subject_col else NULL,
        pair_info = paired_info$pair_info, norm = norm, log_base = log_base, pseudocount = pseudocount,
        threshold = threshold, n_transcripts_tested = nrow(summary_results$all_transcript_stats),
        n_fdr_significant = n_fdr_sig_total, lm_results_provided = !is.null(lm_results),
        lm_p_threshold = lm_p_threshold, lm_genes_filtered = lm_setup$lm_genes_filtered,
        gene_processing_log = gene_results$gene_processing_log)

    result_gene_names <- unname(gene_id_to_name[names(summary_results$results_per_gene)])
    if (length(result_gene_names) == 0)
        result_gene_names <- rep(NA_character_, length(summary_results$results_per_gene))

    result <- list(gene_names = names(summary_results$results_per_gene), gene_ids = names(summary_results$results_per_gene),
        gene_name_map = result_gene_names, conditions = conditions, results_per_gene = summary_results$results_per_gene,
        summary_table = summary_results$summary_table, all_transcript_stats = summary_results$all_transcript_stats,
        metadata = metadata)

    class(result) <- c("tsenat_isoform_switching", "list")

    # 8. Verbose output
    if (verbose)
        .jis_print_results(result, paired_info, lm_results, lm_p_threshold, lm_setup$lm_genes_filtered)

    return(invisible(result))
}

# HELPER FUNCTIONS - Internal Implementation Details
# ===================================================
#' Validate input parameters

#' @noRd
.jis_validate_input <- function(se, condition_col, gene_col, isoform_col) {
    if (is.null(se))
        stop("SummarizedExperiment object (se) is required")
    if (!inherits(se, "SummarizedExperiment"))
        stop("se must be a SummarizedExperiment object")
    if (!(condition_col %in% colnames(colData(se))))
        stop("condition_col '", condition_col, "' not found in colData")
    conditions <- unique(colData(se)[[condition_col]])
    if (length(conditions) != 2)
        stop("Exactly 2 conditions required, found: ", length(conditions))
    if (is.null(gene_col))
        stop("gene_col must be specified")
    if (is.null(isoform_col))
        stop("isoform_col must be specified")
    if (!(gene_col %in% colnames(rowData(se))))
        stop("gene_col '", gene_col, "' not found in rowData")
    if (!(isoform_col %in% colnames(rowData(se))))
        stop("isoform_col '", isoform_col, "' not found in rowData")
    sort(conditions)
}

#' Setup paired design if applicable

#' @noRd
.setup_paired_design_jis <- function(se, subject_col, condition_col) {
    if (is.null(subject_col))
        return(list(is_paired = FALSE, pair_info = NULL, subject_col = NULL))
    if (!(subject_col %in% colnames(colData(se))))
        stop("subject_col '", subject_col, "' not found in colData")
    pairs <- colData(se)[[subject_col]]
    conds <- colData(se)[[condition_col]]
    pair_conds <- table(pairs, conds)

    # BUG FIX #3: Validate paired design - check for unmatched subjects
    conditions <- unique(conds)
    matched_info <- list()
    for (condition in conditions) {
        subjects_in_cond <- names(which(pair_conds[, condition] > 0))
        matched_info[[condition]] <- subjects_in_cond
    }

    unmatched_subjects <- setdiff(matched_info[[1]], matched_info[[2]])
    if (length(unmatched_subjects) > 0 && length(unmatched_subjects) < length(matched_info[[1]])) {
        warning(sprintf("Paired design: %d/%d subjects unmatched in condition '%s'. ",
            length(unmatched_subjects), length(matched_info[[1]]), conditions[1]),
            "Only matched subjects will be analyzed.")
    }

    paired_ratio <- sum(apply(pair_conds > 0, 1, sum) == 2)/nrow(pair_conds)
    if (paired_ratio >= 0.5) {
        paired_indices <- apply(pair_conds > 0, 1, sum) == 2
        list(is_paired = TRUE, pair_info = list(subject_col = subject_col, n_pairs = sum(paired_indices),
            matched_pairs = names(which(paired_indices))), subject_col = subject_col)
    } else {
        list(is_paired = FALSE, pair_info = NULL, subject_col = NULL)
    }
}

#' Build gene ID to name mapping

#' @noRd
.build_gene_id_mapping <- function(se, gene_col) {
    rd_mapping <- rowData(se)
    gene_id_to_name <- character(0)
    if ("gene_name" %in% colnames(rd_mapping)) {
        gene_id_to_name <- setNames(as.character(rd_mapping$gene_name), as.character(rd_mapping[[gene_col]]))
        gene_id_to_name <- gene_id_to_name[!duplicated(names(gene_id_to_name))]
    }
    gene_id_to_name
}

#' Normalize pseudocount parameter

#' @noRd
.jis_normalize_pseudocount <- function(pseudocount, min_value = 1e-08) {
    if (pseudocount <= 0)
        min_value else pseudocount
}

#' Tsallis Entropy - C++ Wrapper [PRODUCTION]
#' 
#' Optimized C++ implementation for Tsallis entropy computation.
#' Requires compiled Rcpp code.
#' 
#' @noRd
.jis_tsallis_entropy_fast <- function(counts, q, norm, log_base, pseudocount, n_tx_fixed = NULL) {
    if (is.null(n_tx_fixed))
        n_tx_fixed <- -1L
    jis_tsallis_entropy_cpp(counts, q = q, normalize = norm, log_base = log_base,
        pseudocount = pseudocount, n_tx_fixed = n_tx_fixed)
}

#' Jackknife Influences - C++ Wrapper [PRODUCTION]
#' 
#' Optimized C++ implementation for jackknife influence computation.
#' Requires compiled Rcpp code.
#' 
#' @noRd
.jis_jackknife_influences_fast <- function(counts, q, norm, log_base, pseudocount,
    n_tx_fixed = NULL) {
    if (is.null(n_tx_fixed))
        n_tx_fixed <- -1L
    jis_jackknife_influences_cpp(counts, q = q, normalize = norm, log_base = log_base,
        pseudocount = pseudocount, n_tx_fixed = n_tx_fixed)
}

#' C++ Wrapper: Fast bootstrap delta statistics computation
#' 
#' Calls C++ implementation if available (Rcpp compiled) and translates
#' field names,
#' otherwise falls back to pure R
#' 
#' @noRd
.jis_bootstrap_delta_fast <- function(counts_A, counts_B, delta_influence, q = 1,
    norm = TRUE, log_base = exp(1), pseudocount = 0, n_bootstrap = 1000, confidence = 0.95,
    method = "percentile", n_transcripts = NULL) {
    tryCatch({
        if (exists("jis_bootstrap_delta_cpp", mode = "function")) {
            result_cpp <- jis_bootstrap_delta_cpp(counts_A, counts_B, delta_influence,
                q = q, normalize = norm, log_base = log_base, pseudocount = pseudocount,
                n_bootstrap = n_bootstrap, confidence = confidence, method = method)

            # Translate field names from C++ (p_value) to R (pvalue) convention
            result_R <- list(delta_influence = if (!is.null(result_cpp$delta_influence)) {
                result_cpp$delta_influence
            } else {
                delta_influence
            }, variance = result_cpp$variance, se = sqrt(result_cpp$variance), ci_lower = result_cpp$ci_lower,
                ci_upper = result_cpp$ci_upper, pvalue = result_cpp$p_value, effect_size = result_cpp$effect_size,
                ci_width = result_cpp$ci_width, relative_ci_width = result_cpp$relative_ci_width,
                power_assessment = NA)
            return(result_R)
        }

        # Fallback to R implementation
        .compute_delta_statistics(counts_A, counts_B, delta_influence, q = q, norm = norm,
            log_base = log_base, pseudocount = pseudocount, n_bootstrap = n_bootstrap,
            n_transcripts = n_transcripts)
    }, error = function(e) {
        # Fallback to R if C++ fails
        .compute_delta_statistics(counts_A, counts_B, delta_influence, q = q, norm = norm,
            log_base = log_base, pseudocount = pseudocount, n_bootstrap = n_bootstrap,
            n_transcripts = n_transcripts)
    })
}



#' Setup LM filtering and gene mapping

#' @noRd
.jis_setup_lm_filtering <- function(se, lm_results, lm_p_threshold, use_lm_fdr, gene_ids,
    gene_col) {
    if (is.null(lm_results))
        return(list(lm_gene_mapping = NULL, filtered_genes = gene_ids, lm_genes_filtered = 0))

    if (!("gene" %in% colnames(lm_results)))
        stop("lm_results must have 'gene' column")

    # Detect if genes are IDs or names
    sample_lm_genes <- lm_results$gene[seq_len(min(5, nrow(lm_results)))]
    genes_are_ids <- all(sample_lm_genes %in% gene_ids)

    # Map gene names to IDs if needed
    if (!genes_are_ids) {
        rd <- rowData(se)
        gene_name_to_id <- setNames(as.character(rd[[gene_col]]), as.character(rd$gene_name))
        gene_name_to_id <- gene_name_to_id[!is.na(names(gene_name_to_id))]
        gene_name_to_id <- gene_name_to_id[!duplicated(names(gene_name_to_id))]
        lm_results$gene <- unname(gene_name_to_id[as.character(lm_results$gene)])
        lm_results <- lm_results[!is.na(lm_results$gene), ]
    }

    if (nrow(lm_results) == 0)
        return(list(lm_gene_mapping = NULL, filtered_genes = gene_ids, lm_genes_filtered = 0))

    # Filter by p-value threshold
    p_col <- if (use_lm_fdr && "adj_p_interaction" %in% colnames(lm_results))
        "adj_p_interaction" else if ("p_interaction" %in% colnames(lm_results))
        "p_interaction" else NULL

    if (!is.null(p_col)) {
        sig_genes <- lm_results[lm_results[[p_col]] < lm_p_threshold, "gene"]
        filtered_genes <- intersect(gene_ids, sig_genes)
        lm_genes_filtered <- length(sig_genes)
    } else {
        filtered_genes <- gene_ids
        lm_genes_filtered <- nrow(lm_results)
    }

    list(lm_gene_mapping = lm_results, filtered_genes = filtered_genes, lm_genes_filtered = lm_genes_filtered)
}

#' Process all genes for isoform switching analysis

#' @noRd
.jis_process_all_genes <- function(se, gene_ids, gene_col, isoform_col, condition_col,
    conditions, paired_info, q, norm, log_base, pseudocount, n_bootstrap, lm_gene_mapping) {
    results_per_gene <- list()
    all_pvalues <- list()
    gene_processing_log <- data.frame(gene = character(), n_transcripts = numeric(),
        has_2_transcripts = logical(), in_results = logical(), stringsAsFactors = FALSE)

    for (gene in gene_ids) {
        gene_mask <- rowData(se)[[gene_col]] == gene
        gene_isos <- rowData(se)[gene_mask, isoform_col]

        if (length(gene_isos) < 2) {
            gene_processing_log <- rbind(gene_processing_log, data.frame(gene = gene,
                n_transcripts = length(gene_isos), has_2_transcripts = FALSE, in_results = FALSE))
            next
        }

        # Extract condition-specific counts
        counts_matrix <- assays(se)$counts[gene_mask, , drop = FALSE]
        cond_mask_A <- colData(se)[[condition_col]] == conditions[1]
        cond_mask_B <- colData(se)[[condition_col]] == conditions[2]
        counts_A_all <- counts_matrix[, cond_mask_A, drop = FALSE]
        counts_B_all <- counts_matrix[, cond_mask_B, drop = FALSE]

        # Handle paired design if applicable
        if (paired_info$is_paired) {
            pairs_A <- colData(se)[cond_mask_A, paired_info$subject_col]
            pairs_B <- colData(se)[cond_mask_B, paired_info$subject_col]
            matched_pairs <- intersect(pairs_A, pairs_B)
            match_idx_A <- match(matched_pairs, pairs_A, nomatch = NA)
            match_idx_B <- match(matched_pairs, pairs_B, nomatch = NA)
            valid_pairs <- !is.na(match_idx_A) & !is.na(match_idx_B)
            if (any(valid_pairs)) {
                cond_A_idx <- which(cond_mask_A)[match_idx_A[valid_pairs]]
                cond_B_idx <- which(cond_mask_B)[match_idx_B[valid_pairs]]
                counts_A <- counts_matrix[, cond_A_idx, drop = FALSE]
                counts_B <- counts_matrix[, cond_B_idx, drop = FALSE]
            } else {
                counts_A <- counts_A_all
                counts_B <- counts_B_all
            }
        } else {
            counts_A <- counts_A_all
            counts_B <- counts_B_all
        }

        # Compute delta influence
        n_tx_original <- nrow(counts_A)
        delta_influence <- .jis_jackknife_influences_fast(counts_A, q, norm, log_base,
            pseudocount, n_tx_original) - .jis_jackknife_influences_fast(counts_B,
            q, norm, log_base, pseudocount, n_tx_original)

        # Get bootstrap statistics
        delta_stats <- .jis_bootstrap_delta_fast(counts_A, counts_B, delta_influence,
            q = q, norm = norm, log_base = log_base, pseudocount = pseudocount, n_bootstrap = n_bootstrap,
            confidence = 0.95, method = "percentile", n_transcripts = nrow(counts_A))

        # Determine switching status
        switching_status <- ifelse(delta_influence > 0, "up", ifelse(delta_influence <
            0, "down", "neutral"))

        # Compute effect size
        max_abs_influence <- max(abs(c(.jis_jackknife_influences_fast(counts_A, q,
            norm, log_base, pseudocount, nrow(counts_A)), .jis_jackknife_influences_fast(counts_B,
            q, norm, log_base, pseudocount, nrow(counts_A)))), na.rm = TRUE)
        effect_size <- ifelse(max_abs_influence > 0, abs(delta_influence)/max_abs_influence,
            0)
        ci_width <- delta_stats$ci_upper - delta_stats$ci_lower
        relative_ci_width <- ifelse(abs(delta_influence) > 1e-10, ci_width/(abs(delta_influence) +
            1e-10), NA)

        # Build gene result
        gene_result <- list(gene_id = gene, transcript_ids = as.character(gene_isos),
            delta_influence = delta_influence, delta_se = delta_stats$se, delta_ci_lower = delta_stats$ci_lower,
            delta_ci_upper = delta_stats$ci_upper, delta_pvalue = delta_stats$pvalue,
            switching_status = switching_status, effect_size = effect_size, ci_width = ci_width,
            relative_ci_width = relative_ci_width, power_assessment = delta_stats$power_assessment)

        # Add LM results if available
        if (!is.null(lm_gene_mapping)) {
            lm_row <- lm_gene_mapping[lm_gene_mapping$gene == gene, ]
            if (nrow(lm_row) > 0) {
                gene_result$lm_p_interaction <- if ("p_interaction" %in% colnames(lm_row))
                  lm_row$p_interaction[1] else NA
                gene_result$lm_adj_p_interaction <- if ("adj_p_interaction" %in%
                  colnames(lm_row))
                  lm_row$adj_p_interaction[1] else NA
            }
        }

        results_per_gene[[gene]] <- gene_result
        gene_processing_log <- rbind(gene_processing_log, data.frame(gene = gene,
            n_transcripts = length(gene_isos), has_2_transcripts = TRUE, in_results = TRUE))

        # Record p-values
        for (i in seq_along(gene_isos)) {
            all_pvalues[[length(all_pvalues) + 1]] <- list(gene = gene, transcript = gene_isos[i],
                pvalue = delta_stats$pvalue[i])
        }
    }

    list(results_per_gene = results_per_gene, all_pvalues = all_pvalues, gene_processing_log = gene_processing_log)
}

#' Build summary results from per-gene analysis

#' @noRd
.jis_build_summary_results <- function(se, results_per_gene, all_pvalues, gene_col,
    lm_gene_mapping) {
    # Apply FDR correction first
    results_per_gene <- .jis_apply_fdr(results_per_gene, all_pvalues)

    # Build transcript-level statistics
    all_transcript_stats <- do.call(rbind, lapply(names(results_per_gene), function(gene) {
        res <- results_per_gene[[gene]]
        fdr_vals <- if (!is.null(res$delta_fdr) && length(res$delta_fdr) > 0)
            res$delta_fdr else rep(NA_real_, length(res$transcript_ids))
        pval_vals <- if (!is.null(res$delta_pvalue) && length(res$delta_pvalue) >
            0)
            res$delta_pvalue else rep(NA_real_, length(res$transcript_ids))

        df <- data.frame(gene = gene, transcript_id = res$transcript_ids, pvalue = pval_vals,
            fdr = fdr_vals, stringsAsFactors = FALSE)
        if (!is.null(res$lm_p_interaction))
            df$lm_p_interaction <- res$lm_p_interaction
        if (!is.null(res$lm_adj_p_interaction))
            df$lm_adj_p_interaction <- res$lm_adj_p_interaction
        df
    }))

    # Build gene-level summary
    summary_rows <- lapply(names(results_per_gene), function(gene) {
        res <- results_per_gene[[gene]]
        delta_vals <- res$delta_influence[is.finite(res$delta_influence)]
        gene_name <- gene
        if ("gene_name" %in% colnames(rowData(se))) {
            gn <- rowData(se)[which(rowData(se)[[gene_col]] == gene)[1], "gene_name"]
            if (!is.na(gn))
                gene_name <- as.character(gn)
        }
        data.frame(gene = gene, gene_name = gene_name, n_transcripts = length(res$transcript_ids),
            max_delta_influence = if (length(delta_vals) > 0)
                max(abs(delta_vals)) else 0, n_switching_transcripts = sum(res$switching_status != "neutral",
                na.rm = TRUE), n_fdr_significant = sum(res$delta_fdr < 0.05, na.rm = TRUE),
            stringsAsFactors = FALSE)
    })

    summary_table <- do.call(rbind, summary_rows)
    rownames(summary_table) <- NULL

    list(results_per_gene = results_per_gene, all_transcript_stats = all_transcript_stats,
        summary_table = summary_table)
}

#' Apply FDR correction

#' @noRd
.jis_apply_fdr <- function(results_per_gene, all_pvalues) {
    if (length(all_pvalues) == 0)
        return(invisible(results_per_gene))
    pvals_vec <- vapply(all_pvalues, function(x) x$pvalue, FUN.VALUE = numeric(1))
    fdr_vec <- p.adjust(pvals_vec, method = "BH")

    # OPTIMIZED: Vectorized FDR assignment (replaces loop with direct
    # assignment)
    for (i in seq_along(all_pvalues)) {
        all_pvalues[[i]]$fdr <- fdr_vec[i]
    }

    # OPTIMIZED: Pre-compute matches to avoid repeated lookups in loop
    genes <- vapply(all_pvalues, function(x) x$gene, FUN.VALUE = character(1))
    transcripts <- vapply(all_pvalues, function(x) x$transcript, FUN.VALUE = character(1))
    fdr_values <- fdr_vec

    # Single vectorized pass through genes
    gene_names_results <- names(results_per_gene)

    for (px_idx in seq_along(all_pvalues)) {
        px <- all_pvalues[[px_idx]]
        gene_idx <- match(px$gene, gene_names_results)
        if (!is.na(gene_idx)) {
            trans_ids <- results_per_gene[[gene_idx]]$transcript_ids
            iso_idx <- match(px$transcript, trans_ids)
            if (!is.na(iso_idx)) {
                if (is.null(results_per_gene[[gene_idx]]$delta_fdr)) {
                  results_per_gene[[gene_idx]]$delta_fdr <- rep(NA, length(trans_ids))
                }
                results_per_gene[[gene_idx]]$delta_fdr[iso_idx] <- px$fdr
            }
        }
    }
    invisible(results_per_gene)
}

#' Handle multi-q analysis

#' @noRd
.jis_handle_multi_q <- function(se, q, q_params, verbose) {
    results_list <- lapply(q, function(q_val) {
        do.call(.jackknife_isoform_switching, c(list(se = se, q = q_val, verbose = FALSE),
            q_params))
    })
    names(results_list) <- paste0("q_", gsub("\\.", "_", sprintf("%.2f", q)))
    for (gene_idx in seq_along(results_list[[1]]$results_per_gene)) {
        gene_name <- names(results_list[[1]]$results_per_gene)[gene_idx]
        transcript_ids <- results_list[[1]]$results_per_gene[[gene_name]]$transcript_ids
        delta_matrix <- matrix(NA_real_, nrow = length(transcript_ids), ncol = length(q))
        colnames(delta_matrix) <- paste0("q_", gsub("\\.", "_", sprintf("%.2f", q)))
        rownames(delta_matrix) <- transcript_ids
        for (q_idx in seq_along(results_list)) {
            if (!is.null(results_list[[q_idx]]$results_per_gene[[gene_name]])) {
                delta_matrix[, q_idx] <- results_list[[q_idx]]$results_per_gene[[gene_name]]$delta_influence
            }
        }
        consistency_results <- apply(delta_matrix, 1, function(x) {
            x_valid <- x[!is.na(x) & !is.infinite(x)]
            if (length(x_valid) >= 2) {
                # Ignore zeros when checking for consistency (0 is a
                # neutral/boundary point)
                x_nonzero <- x_valid[x_valid != 0]
                if (length(x_nonzero) == 0) {
                  "All zero"
                } else {
                  pos_count <- sum(x_nonzero > 0)
                  neg_count <- sum(x_nonzero < 0)
                  # Check if ALL non-zero values have the same sign
                  if (pos_count == length(x_nonzero))
                    "Consistent positive" else if (neg_count == length(x_nonzero))
                    "Consistent negative" else "Mixed directions"
                }
            } else if (length(x_valid) == 1)
                "Single q-value" else "No data"
        })
        for (q_idx in seq_along(results_list)) {
            if (!is.null(results_list[[q_idx]]$results_per_gene[[gene_name]])) {
                results_list[[q_idx]]$results_per_gene[[gene_name]]$direction_consistency <- consistency_results
            }
        }
    }
    class(results_list) <- c("tsenat_isoform_switching_multiq", "list")
    if (verbose) {
        message("Isoform Switching Analysis - Multi-Q Comparison")
        message("================================================")
        for (i in seq_along(results_list)) {
            message(sprintf("q = %s", q[i]))
            res <- results_list[[i]]
            if (!is.null(res$metadata)) {
                message(sprintf("  Genes analyzed:      %d", length(res$gene_names)))
                message(sprintf("  Transcripts tested:  %d", res$metadata$n_transcripts_tested))
                message(sprintf("  FDR-significant:     %d", res$metadata$n_fdr_significant))
                message(sprintf("  Genes with switching: %d", sum(res$summary_table$n_switching_transcripts >
                  0)))
            }
        }
        message("[OK] Access results$q_<value>$results_per_gene$<gene> for per-q, per-gene details")
        message("[OK] Compare q values to assess scale-dependent isoform switching patterns")
    }
    invisible(results_list)
}



#' Print isoform switching analysis results

#' @noRd
.jis_print_results <- function(result, paired_info, lm_results, lm_p_threshold, lm_genes_filtered) {
    message("Isoform Switching Analysis Results")
    message("===================================")
    message(sprintf("Conditions: '%s' vs. '%s'", result$conditions[1], result$conditions[2]))
    if (paired_info$is_paired && !is.null(paired_info$pair_info))
        message(sprintf("Design: PAIRED (%s) - %d matched pairs", paired_info$subject_col,
            paired_info$pair_info$n_pairs)) else message("Design: UNPAIRED")
    if (!is.null(lm_results))
        message(sprintf("LM filtering: Genes with p < %.6f (N = %d)", lm_p_threshold,
            lm_genes_filtered))
    message(sprintf("\nGenes analyzed: %d", length(result$gene_names)))
    message(sprintf("Total transcripts tested: %d", result$metadata$n_transcripts_tested))
    message(sprintf("FDR-significant transcripts (FDR<0.05): %d", result$metadata$n_fdr_significant))
    message("Summary Table:")
    message(paste(capture.output(str(result$summary_table)), collapse = "\n"))
    message("\n[OK] Use results$results_per_gene$'GeneName' to access per-gene switching details")
    message("[OK] Use results$summary_table for overview across genes")
    message("[OK] Use results$all_transcript_stats for FDR-corrected p-values per transcript")
    message("[OK] Use results$metadata$is_paired to check if paired design was applied")
    if (!is.null(lm_results))
        message("[OK] Access lm_p_interaction in each gene$lm_p_interaction for LM test results")
    message("")
}

#' Bootstrap Helper Function for Delta Statistics
#'
#' Internal function to compute bootstrap confidence intervals and p-values
#' for delta_influence (difference in Tsallis entropy influence between
#' conditions).
#' Computes per-transcript statistics from bootstrap resamples.
#'

#' @noRd

.compute_delta_statistics <- function(counts_A, counts_B, delta_influence, q = 1,
    norm = TRUE, log_base = exp(1), pseudocount = 0, n_bootstrap = 1000, confidence = 0.95,
    n_transcripts = NULL) {
    .calculate_tsallis <- function(counts, q, norm, log_base, pseudocount, n_transcripts_fixed) {
        if (is.vector(counts))
            counts <- t(as.matrix(counts))

        # Check for zero-sum columns BEFORE adding pseudocount to detect
        # samples with no expression
        raw_col_sums <- colSums(counts)
        with_zero_counts <- raw_col_sums == 0

        # BUG FIX #2: Use centralized pseudocount normalization
        pseudocount <- .jis_normalize_pseudocount(pseudocount)
        counts <- counts + pseudocount

        col_sums <- colSums(counts)
        # Avoid division by zero
        if (any(col_sums <= 0)) {
            return(rep(NA_real_, ncol(counts)))
        }

        # Avoid division by near-zero for zero-count columns by using raw
        # column sums intelligently For zero-count columns, set col_sums to 1
        # to avoid division by near-zero
        col_sums_safe <- pmax(col_sums, 1)

        # Normalize by column (sample): divide each column by its total
        p <- t(t(counts)/col_sums_safe)

        if (q == 1) {
            # Shannon entropy: H = -sum(p_i * log(p_i)) Use epsilon offset to
            # avoid log(0) while preserving probability conservation
            eps <- 1e-15
            if (log_base == exp(1)) {
                h <- -colSums(p * log(pmax(p, eps)))
            } else {
                h <- -colSums(p * log(pmax(p, eps)))/log(log_base)
            }
        } else {
            # Tsallis entropy: (1 - sum(p^q)) / (q - 1)
            p_q_sum <- colSums(p^q)
            h <- (1 - p_q_sum)/(q - 1)
        }

        h[!is.finite(h)] <- NA_real_
        # Set entropy to NA for zero-count columns
        h[with_zero_counts] <- NA_real_

        if (norm) {
            if (q == 1) {
                # Use FIXED n_transcripts if provided, otherwise use current
                # matrix dimensions
                n_tx <- if (!is.null(n_transcripts_fixed))
                  n_transcripts_fixed else nrow(counts)
                max_h <- log(n_tx)/log(log_base)
            } else {
                n_tx <- if (!is.null(n_transcripts_fixed))
                  n_transcripts_fixed else nrow(counts)
                # Tsallis: (1 - n^(1-q)) / (q - 1) is always positive by
                # definition CRITICAL: Use (q - 1), not (1/(1-q)), to get
                # correct sign
                max_h <- (1 - n_tx^(1 - q))/(q - 1)
            }
            if (!is.na(max_h) && is.finite(max_h) && max_h > 0) {
                h <- h/max_h
            }
        }
        return(h)
    }

    .jackknife_entropy <- function(counts, q, norm, log_base, pseudocount, n_transcripts_fixed) {
        if (is.vector(counts))
            counts <- matrix(counts, nrow = 1)
        n_tx <- nrow(counts)
        h_full <- .calculate_tsallis(counts, q, norm, log_base, pseudocount, n_transcripts_fixed)
        influences <- numeric(n_tx)

        for (i in seq_len(n_tx)) {
            h_leave_i <- .calculate_tsallis(counts[-i, , drop = FALSE], q, norm,
                log_base, pseudocount, n_transcripts_fixed)
            influences[i] <- mean(abs(h_full - h_leave_i), na.rm = TRUE)
        }
        list(full_entropy = h_full, influences = influences)
    }

    n_tx <- length(delta_influence)
    # Seed handling left to caller for Bioconductor compliance
    bootstrap_deltas_matrix <- matrix(nrow = n_bootstrap, ncol = n_tx)

    for (b in seq_len(n_bootstrap)) {
        idx_A <- sample(seq_len(ncol(counts_A)), size = ncol(counts_A), replace = TRUE)
        idx_B <- sample(seq_len(ncol(counts_B)), size = ncol(counts_B), replace = TRUE)

        boot_A <- counts_A[, idx_A, drop = FALSE]
        boot_B <- counts_B[, idx_B, drop = FALSE]

        jack_A <- .jackknife_entropy(boot_A, q, norm, log_base, pseudocount, n_transcripts)
        jack_B <- .jackknife_entropy(boot_B, q, norm, log_base, pseudocount, n_transcripts)

        # Compute per-transcript delta for this bootstrap sample
        bootstrap_deltas_matrix[b, ] <- jack_A$influences - jack_B$influences
    }

    # Compute per-transcript statistics
    alpha <- 1 - confidence
    # Use type=1 (nearest-rank) quantile method for consistency with C++
    # implementation
    ci_lower <- apply(bootstrap_deltas_matrix, 2, function(x) quantile(x, alpha/2,
        na.rm = TRUE, type = 1))
    ci_upper <- apply(bootstrap_deltas_matrix, 2, function(x) quantile(x, 1 - alpha/2,
        na.rm = TRUE, type = 1))

    # Compute two-tailed bootstrap p-values Two-tailed test: proportion of
    # bootstrap samples with |bootstrap_delta| >= |observed_delta|
    pvalues <- numeric(n_tx)
    for (i in seq_len(n_tx)) {
        obs_abs <- abs(delta_influence[i])
        boot_abs <- abs(bootstrap_deltas_matrix[, i])
        # Two-tailed p-value: proportion of bootstrap samples as or more
        # extreme than observed Handle case where all bootstrap values are
        # NA/NaN
        valid_boots <- !is.na(boot_abs) & is.finite(boot_abs)
        if (!is.na(obs_abs) && is.finite(obs_abs) && any(valid_boots)) {
            pvalues[i] <- mean(boot_abs[valid_boots] >= obs_abs, na.rm = FALSE)
            pvalues[i] <- max(pvalues[i], 1/n_bootstrap)  # Minimum p-value = 1/n_bootstrap
        } else {
            pvalues[i] <- NA_real_  # Return NA if insufficient data
        }
    }

    # Standard error per transcript
    se <- apply(bootstrap_deltas_matrix, 2, sd, na.rm = TRUE)

    # BUG FIX #1: Handle degenerate cases (zero-width CIs)
    ci_width <- ci_upper - ci_lower
    degenerate_idx <- which(ci_width < 1e-10)
    if (length(degenerate_idx) > 0) {
        warning(sprintf("Degenerate bootstrap distributions detected for %d transcript(s). ",
            length(degenerate_idx)), "CI width = 0. Consider increasing nboot or checking data quality.")
    }

    # IMPROVEMENT #2: Add power/sample size assessment
    power_assessment <- data.frame(n_effective = n_bootstrap * (1 - sum(is.na(bootstrap_deltas_matrix))/length(bootstrap_deltas_matrix)),
        avg_ci_width = mean(ci_width, na.rm = TRUE), min_recommended_nboot = NA_integer_)

    if (power_assessment$avg_ci_width > 0.05) {
        power_assessment$min_recommended_nboot <- ceiling(n_bootstrap * (power_assessment$avg_ci_width/0.05)^2)
    }

    return(list(ci_lower = as.numeric(ci_lower), ci_upper = as.numeric(ci_upper),
        pvalue = pvalues, se = se, power_assessment = power_assessment))
}
