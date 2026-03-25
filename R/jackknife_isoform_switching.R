#' Jackknife Isoform Switching Detection
#'
#' Detects isoform/transcript switches between conditions using condition-stratified
#' jackknife analysis on Tsallis entropy. Identifies which transcripts change
#' importance between conditions (e.g., condition A vs B).
#'
#' @param se A SummarizedExperiment object with transcript-level counts.
#' @param condition_col Character: column name in colData for condition labels.
#' @param subject_col Character: optional column name for paired design (subject/individual IDs).
#' @param gene_col Character: column name in rowData for gene IDs.
#' @param isoform_col Character: column name in rowData for transcript/isoform IDs.
#' @param q Numeric: Tsallis entropy order (default 1 = Shannon entropy). Can be a vector
#'   for multi-q analysis (e.g., q = c(0.5, 1.0, 1.5, 2.0)); results will be nested by q value.
#' @param norm Logical: normalize entropy to [0,1]? (default TRUE).
#' @param log_base Numeric: log base for entropy (default e).
#' @param pseudocount Numeric: pseudocount to add (default 0).
#' @param threshold Numeric: percentile for outlier detection on influences (default 90).
#' @param n_bootstrap Numeric: number of bootstrap resamples (default 1000).
#' @param print_results Logical: print results? (default TRUE).
#' @param verbose Logical: verbose output? (default FALSE).
#' @param lm_results Data frame: results from calculate_lm_interaction() with 'gene' column.
#'   Can contain either gene names or gene IDs; function automatically maps names to IDs
#'   using rowData(se). Include 'p_interaction' and/or 'adj_p_interaction' columns for
#'   filtering genes by significance. When provided, only genes passing lm_p_threshold are
#'   analyzed; all matching genes are included (top_n parameter removed).
#' @param lm_p_threshold Numeric: p-value threshold for LM gene filtering (default 0.05).
#' @param use_lm_fdr Logical: use adjusted p-values from LM results if available (default TRUE).
#'
#' @return If q is a single value, returns a list of class tsenat_isoform_switching with:
#'   \describe{
#'     \item{results_per_gene}{named list of per-gene results}
#'     \item{summary_table}{data.frame with per-gene summary}
#'     \item{all_transcript_stats}{data.frame with all transcript statistics}
#'     \item{gene_names}{character vector of analyzed genes}
#'     \item{conditions}{character vector of the two conditions compared}
#'     \item{metadata}{list with analysis metadata}
#'   }
#'   
#'   If q is a vector, returns a list of class tsenat_isoform_switching_multiq where each
#'   element is a complete tsenat_isoform_switching result for that q value. Keys are
#'   formatted as "q_X_XX" for ease of iteration (e.g., q_0_01, q_1_00, q_2_00).
#'
#' @section Sample Metadata Parameters (Unified Naming Convention):
#' TSENAT functions use consistent parameter names for sample grouping and subject identification:
#' \itemize{
#'   \item{\code{condition_col}: Character string specifying the colData column 
#'         containing sample group/condition labels (e.g., "Normal", "Tumor", "control", "treatment"). 
#'         Default: "condition". Required for identifying the two conditions to compare.}
#'   \item{\code{subject_col}: For paired/blocked designs, character string specifying 
#'         the colData column with subject/individual/patient identifiers. Default: NULL.}
#' }
#' All functions use \code{SummarizedExperiment::colData()} as the single source of truth
#' for sample metadata. This eliminates parameter fragmentation and improves API discoverability
#' across the TSENAT package.
#'
#' @keywords internal
#' @noRd
jackknife_isoform_switching <- function(
  se = NULL,
  condition_col = "condition",
  subject_col = NULL,
  gene_col = NULL,
  isoform_col = NULL,
  q = 1,
  norm = TRUE,
  log_base = exp(1),
  pseudocount = 0,
  threshold = 90,
  n_bootstrap = 1000,
  print_results = TRUE,
  verbose = FALSE,
  lm_results = NULL,
  lm_p_threshold = 0.05,
  use_lm_fdr = TRUE
) {
  
  # Input validation
  if (is.null(se)) {
    stop("SummarizedExperiment object (se) is required")
  }
  
  if (!inherits(se, "SummarizedExperiment")) {
    stop("se must be a SummarizedExperiment object")
  }
  
  # Handle multiple q values
  if (is.numeric(q) && length(q) > 1) {
    # Recursive call for each q value
    results_list <- lapply(q, function(q_val) {
      jackknife_isoform_switching(
        se = se,
        condition_col = condition_col,
        subject_col = subject_col,
        gene_col = gene_col,
        isoform_col = isoform_col,
        q = q_val,
        norm = norm,
        log_base = log_base,
        pseudocount = pseudocount,
        threshold = threshold,
        n_bootstrap = n_bootstrap,
        print_results = FALSE,
        verbose = verbose,
        lm_results = lm_results,
        lm_p_threshold = lm_p_threshold,
        use_lm_fdr = use_lm_fdr
      )
    })
    
    # Create named list with q values formatted as keys
    names(results_list) <- paste0("q_", gsub("\\.", "_", sprintf("%.2f", q)))
    
    # Add direction consistency to multi-q results
    # This analyzes each transcript's delta_influence pattern across q-values
    for (gene_idx in seq_along(results_list[[1]]$results_per_gene)) {
      gene_name <- names(results_list[[1]]$results_per_gene)[gene_idx]
      
      # Collect delta_influence for this gene across all q-values
      transcript_ids <- results_list[[1]]$results_per_gene[[gene_name]]$transcript_ids
      
      # Build matrix: rows = transcripts, cols = q-values
      delta_matrix <- matrix(NA_real_, nrow = length(transcript_ids), ncol = length(q))
      colnames(delta_matrix) <- paste0("q_", gsub("\\.", "_", sprintf("%.2f", q)))
      rownames(delta_matrix) <- transcript_ids
      
      for (q_idx in seq_along(results_list)) {
        if (!is.null(results_list[[q_idx]]$results_per_gene[[gene_name]])) {
          delta_matrix[, q_idx] <- results_list[[q_idx]]$results_per_gene[[gene_name]]$delta_influence
        }
      }
      
      # Compute direction consistency for each transcript
      consistency_results <- apply(delta_matrix, 1, function(x) {
        x_valid <- x[!is.na(x) & !is.infinite(x)]
        if (length(x_valid) >= 2) {
          pos_count <- sum(x_valid > 0)
          neg_count <- sum(x_valid < 0)
          zero_count <- sum(x_valid == 0)
          
          if (pos_count == length(x_valid)) {
            return("Consistent positive")
          } else if (neg_count == length(x_valid)) {
            return("Consistent negative")
          } else if (zero_count == length(x_valid)) {
            return("All zero")
          } else {
            return("Mixed directions")
          }
        } else if (length(x_valid) == 1) {
          return("Single q-value")
        } else {
          return("No data")
        }
      })
      
      # Add direction_consistency to each q-value result for this gene
      for (q_idx in seq_along(results_list)) {
        if (!is.null(results_list[[q_idx]]$results_per_gene[[gene_name]])) {
          results_list[[q_idx]]$results_per_gene[[gene_name]]$direction_consistency <- consistency_results
        }
      }
    }
    
    class(results_list) <- c("tsenat_isoform_switching_multiq", "list")
    
    # Optional printing for multi-q results
    if (print_results) {
      message("Isoform Switching Analysis - Multi-Q Comparison")
      message("================================================")
      for (i in seq_along(results_list)) {
        message(sprintf("q = %s", q[i]))
        res <- results_list[[i]]
        if (!is.null(res$metadata)) {
          message(sprintf("  Genes analyzed:      %d", length(res$gene_names)))
          message(sprintf("  Transcripts tested:  %d", res$metadata$n_transcripts_tested))
          message(sprintf("  FDR-significant:     %d", res$metadata$n_fdr_significant))
          message(sprintf("  Genes with switching: %d", sum(res$summary_table$n_switching_transcripts > 0)))
        }
      }
      message("[OK] Access results$q_<value>$results_per_gene$<gene> for per-q, per-gene details")
      message("[OK] Compare q values to assess scale-dependent isoform switching patterns")
    }
    
    return(invisible(results_list))
  }
  
  if (!(condition_col %in% colnames(colData(se)))) {
    stop("condition_col '", condition_col, "' not found in colData")
  }
  
  # Get conditions
  conditions <- unique(colData(se)[[condition_col]])
  if (length(conditions) != 2) {
    stop("Exactly 2 conditions required, found: ", length(conditions))
  }
  conditions <- sort(conditions)
  
  # Validate gene and isoform columns
  if (is.null(gene_col)) {
    stop("gene_col must be specified")
  }
  if (is.null(isoform_col)) {
    stop("isoform_col must be specified")
  }
  
  if (!(gene_col %in% colnames(rowData(se)))) {
    stop("gene_col '", gene_col, "' not found in rowData")
  }
  if (!(isoform_col %in% colnames(rowData(se)))) {
    stop("isoform_col '", isoform_col, "' not found in rowData")
  }
  
  # Handle paired design
  is_paired <- FALSE
  pair_info <- NULL
  
  if (!is.null(subject_col)) {
    if (!(subject_col %in% colnames(colData(se)))) {
      stop("subject_col '", subject_col, "' not found in colData")
    }
    
    # Check if pairing is valid (50% threshold)
    pairs <- colData(se)[[subject_col]]
    conds <- colData(se)[[condition_col]]
    pair_conds <- table(pairs, conds)
    paired_ratio <- sum(apply(pair_conds > 0, 1, sum) == 2) / nrow(pair_conds)
    
    if (paired_ratio >= 0.5) {
      is_paired <- TRUE
      paired_indices <- apply(pair_conds > 0, 1, sum) == 2
      pair_info <- list(
        subject_col = subject_col,
        n_pairs = sum(paired_indices),
        matched_pairs = names(which(paired_indices))
      )
    }
  }
  
  # Get gene list
  gene_ids <- unique(rowData(se)[[gene_col]])
  
  # Create mapping from gene_id to gene_name from rowData
  rd_mapping <- rowData(se)
  gene_id_to_name <- character(0)
  if ("gene_name" %in% colnames(rd_mapping)) {
    gene_id_to_name <- setNames(
      as.character(rd_mapping$gene_name),
      as.character(rd_mapping[[gene_col]])
    )
    # Remove duplicates, keeping first occurrence
    gene_id_to_name <- gene_id_to_name[!duplicated(names(gene_id_to_name))]
  }
  
  # Handle LM filtering with automatic gene name to ID mapping
  lm_gene_mapping <- NULL
  lm_genes_filtered <- 0
  
  if (!is.null(lm_results)) {
    if (!("gene" %in% colnames(lm_results))) {
      stop("lm_results must have 'gene' column")
    }
    
    # Check if lm_results contains gene names or gene IDs
    # First, check if the "gene" column matches gene IDs in se
    sample_lm_genes <- lm_results$gene[seq_len(min(5, nrow(lm_results)))]
    genes_are_ids <- all(sample_lm_genes %in% gene_ids)
    
    # If genes are names, map them to IDs
    if (!genes_are_ids) {
      # Map gene names to Ensembl IDs from rowData
      rd <- rowData(se)
      gene_name_to_id <- setNames(as.character(rd[[gene_col]]), as.character(rd$gene_name))
      gene_name_to_id <- gene_name_to_id[!is.na(names(gene_name_to_id))]
      gene_name_to_id <- gene_name_to_id[!duplicated(names(gene_name_to_id))]
      
      # Map genes in lm_results
      lm_results <- lm_results
      lm_results$gene <- unname(gene_name_to_id[as.character(lm_results$gene)])
      lm_results <- lm_results[!is.na(lm_results$gene), ]
      
      if (nrow(lm_results) == 0) {
        warning("No genes from lm_results matched in SummarizedExperiment rowData.")
        lm_gene_mapping <- NULL
      } else {
        lm_gene_mapping <- lm_results
      }
    } else {
      lm_gene_mapping <- lm_results
    }
    
    # Determine which p-value column to use
    if (!is.null(lm_gene_mapping)) {
      p_col <- NULL
      if (use_lm_fdr && "adj_p_interaction" %in% colnames(lm_gene_mapping)) {
        p_col <- "adj_p_interaction"
      } else if ("p_interaction" %in% colnames(lm_gene_mapping)) {
        p_col <- "p_interaction"
      } else {
        warning("lm_results missing p_interaction or adj_p_interaction columns")
      }
      
      if (!is.null(p_col)) {
        sig_genes <- lm_gene_mapping[lm_gene_mapping[[p_col]] < lm_p_threshold, "gene"]
        lm_genes_filtered <- length(sig_genes)
        
        if (length(sig_genes) == 0) {
          warning("No genes pass LM threshold p < ", lm_p_threshold,
                  ". Analyzing all genes.")
        } else {
          gene_ids <- intersect(gene_ids, sig_genes)
        }
      }
    }
  }
  
  # Note: top_n parameter removed - all LM-significant genes are analyzed
  
  # Initialize results
  results_per_gene <- list()
  all_pvalues <- list()
  all_fdr <- list()
  summary_rows <- list()
  
  # Track gene processing for debugging (if verbose mode enabled)
  gene_processing_log <- data.frame(
    gene = character(),
    n_transcripts = numeric(),
    has_2_transcripts = logical(),
    in_results = logical(),
    stringsAsFactors = FALSE
  )
  
  # Process each gene
  for (gene in gene_ids) {
    # Get transcript indices for this gene
    gene_mask <- rowData(se)[[gene_col]] == gene
    gene_isos <- rowData(se)[gene_mask, isoform_col]
    
    if (length(gene_isos) < 2) {
      gene_processing_log <- rbind(gene_processing_log, data.frame(
        gene = gene, n_transcripts = length(gene_isos),
        has_2_transcripts = FALSE, in_results = FALSE
      ))
      next # Skip single-transcript genes
    }
    
    # Mark that it passed 2-transcript check
    processing_check_passed <- TRUE
    
    # Get counts for this gene
    counts_matrix <- assays(se)$counts[gene_mask, , drop = FALSE]
    
    # Extract condition A and B counts
    cond_mask_A <- colData(se)[[condition_col]] == conditions[1]
    cond_mask_B <- colData(se)[[condition_col]] == conditions[2]
    
    counts_A_all <- counts_matrix[, cond_mask_A, drop = FALSE]
    counts_B_all <- counts_matrix[, cond_mask_B, drop = FALSE]
    
    # Handle paired design: use matched pairs only
    if (is_paired && !is.null(subject_col)) {
      pairs_A <- colData(se)[cond_mask_A, subject_col]
      pairs_B <- colData(se)[cond_mask_B, subject_col]
      
      matched_pairs <- intersect(pairs_A, pairs_B)
      
      cond_A_paired_idx <- which(cond_mask_A)[match(matched_pairs, pairs_A, nomatch = 0)]
      cond_B_paired_idx <- which(cond_mask_B)[match(matched_pairs, pairs_B, nomatch = 0)]
      
      if (length(cond_A_paired_idx) > 0 && length(cond_B_paired_idx) > 0) {
        counts_A <- counts_matrix[, cond_A_paired_idx, drop = FALSE]
        counts_B <- counts_matrix[, cond_B_paired_idx, drop = FALSE]
      } else {
        counts_A <- counts_A_all
        counts_B <- counts_B_all
      }
    } else {
      counts_A <- counts_A_all
      counts_B <- counts_B_all
    }
    
    # Calculate Tsallis entropy and influences
    # CRITICAL: Need to maintain consistent normalization across jackknife leave-one-out iterations
    n_tx_original <- nrow(counts_A)  # Store original transcript count
    
    .tsallis_entropy <- function(counts, q, norm, log_base, pseudocount, n_tx_fixed = NULL) {
      # Check for zero-sum columns BEFORE adding pseudocount to detect samples with no expression
      raw_col_sums <- colSums(counts)
      with_zero_counts <- raw_col_sums == 0
      
      # Ensure pseudocount is used to avoid log(0) and division by zero
      if (pseudocount == 0) {
        pseudocount <- 1e-8
      }
      counts <- counts + pseudocount
      
      col_sums <- colSums(counts)
      if (any(col_sums <= 0)) {
        return(rep(NA_real_, ncol(counts)))
      }
      
      # Return NA for samples with zero total counts (no expression for this gene)
      h_result <- rep(NA_real_, ncol(counts))
      if (all(with_zero_counts)) {
        return(h_result)
      }
      
      # Avoid division by zero for zero-count columns by using raw column sums intelligently
      # For zero-count columns, set col_sums to 1 to avoid division by near-zero
      # These samples will return NA anyway (no valid probability distribution)
      col_sums_safe <- pmax(col_sums, 1)
      
      p <- counts / col_sums_safe
      
      if (q == 1) {
        h <- if (log_base == exp(1)) {
          -colSums(p * log(p + 1e-100))
        } else {
          -colSums(p * log(p + 1e-100, log_base))
        }
      } else {
        p_q_sum <- colSums(p^q)
        h <- (1 / (1 - q)) * (1 - p_q_sum)
      }
      
      h[!is.finite(h)] <- NA_real_
      # Set entropy to NA for zero-count columns
      h[with_zero_counts] <- NA_real_
      
      if (norm) {
        # Use FIXED n_transcripts if provided (for consistent jackknife normalization)
        n_tx_use <- if (!is.null(n_tx_fixed)) n_tx_fixed else nrow(counts)
        if (q == 1) {
          max_h <- log(n_tx_use) / log(log_base)
        } else {
          # Use absolute value because the formula produces negative values for q<1 and q>1
          max_h <- abs((1 / (1 - q)) * (1 - n_tx_use^(1 - q)))
        }
        
        if (!is.na(max_h) && is.finite(max_h) && max_h > 0) {
          h <- h / max_h
        }
      }
      return(h)
    }
    
    .jackknife_influences <- function(counts, q, norm, log_base, pseudocount, n_tx_fixed = NULL) {
      h_full <- .tsallis_entropy(counts, q, norm, log_base, pseudocount, n_tx_fixed)
      n_tx <- nrow(counts)
      influences <- numeric(n_tx)
      
      if (n_tx < 2) {
        return(influences)
      }
      
      for (i in seq_len(n_tx)) {
        counts_leave_i <- counts[-i, , drop = FALSE]
        h_leave_i <- .tsallis_entropy(counts_leave_i, q, norm, log_base, pseudocount, n_tx_fixed)
        # Compute mean absolute difference across samples
        diffs <- abs(h_full - h_leave_i)
        influences[i] <- mean(diffs, na.rm = TRUE)
      }
      return(influences)
    }
    
    # Calculate influences for each condition with fixed transcript count
    # This ensures consistent normalization across all jackknife iterations
    influences_A <- .jackknife_influences(counts_A, q, norm, log_base, pseudocount, n_tx_original)
    influences_B <- .jackknife_influences(counts_B, q, norm, log_base, pseudocount, n_tx_original)
    
    # Calculate delta influence
    delta_influence <- influences_A - influences_B
    
    # Get original transcript count for consistent normalization during jackknife bootstrap
    n_tx <- nrow(counts_A)
    
    # Calculate bootstrap statistics
    delta_stats <- compute_delta_statistics(
      counts_A, counts_B, delta_influence,
      q = q, norm = norm, log_base = log_base,
      pseudocount = pseudocount, n_bootstrap = n_bootstrap,
      n_transcripts = n_tx  # Pass original transcript count for consistent normalization
    )
    
    # Determine switching status
    switching_status <- ifelse(delta_influence > 0, "up", ifelse(delta_influence < 0, "down", "neutral"))
    
    # Store gene result (delta_stats now contains per-transcript vectors)
    gene_result <- list(
      gene_id = gene,
      transcript_ids = as.character(gene_isos),
      delta_influence = delta_influence,
      delta_se = delta_stats$se,
      delta_ci_lower = delta_stats$ci_lower,
      delta_ci_upper = delta_stats$ci_upper,
      delta_pvalue = delta_stats$pvalue,
      switching_status = switching_status
    )
    
    # Add LM annotation if available
    if (!is.null(lm_gene_mapping)) {
      lm_row <- lm_gene_mapping[lm_gene_mapping$gene == gene, ]
      if (nrow(lm_row) > 0) {
        gene_result$lm_p_interaction <- if ("p_interaction" %in% colnames(lm_row)) {
          lm_row$p_interaction[1]
        } else {
          NA
        }
        gene_result$lm_adj_p_interaction <- if ("adj_p_interaction" %in% colnames(lm_row)) {
          lm_row$adj_p_interaction[1]
        } else {
          NA
        }
      }
    }
    
    results_per_gene[[gene]] <- gene_result
    
    # Log successful processing
    gene_processing_log <- rbind(gene_processing_log, data.frame(
      gene = gene, n_transcripts = length(gene_isos),
      has_2_transcripts = TRUE, in_results = TRUE
    ))
    
    # Collect p-values for global FDR correction (per-transcript)
    for (i in seq_along(gene_isos)) {
      all_pvalues[[length(all_pvalues) + 1]] <- list(
        gene = gene,
        transcript = gene_isos[i],
        pvalue = delta_stats$pvalue[i]  # Per-transcript p-value
      )
    }
  }
  
  # Global FDR correction (Benjamini-Hochberg)
  if (length(all_pvalues) > 0) {
    pvals_vec <- vapply(all_pvalues, function(x) x$pvalue, FUN.VALUE = numeric(1))
    fdr_vec <- p.adjust(pvals_vec, method = "BH")
    
    for (i in seq_along(all_pvalues)) {
      all_pvalues[[i]]$fdr <- fdr_vec[i]
    }
    
    # Add FDR to gene results
    for (px in all_pvalues) {
      gene_idx <- match(px$gene, names(results_per_gene))
      if (!is.na(gene_idx)) {
        iso_idx <- match(px$transcript, results_per_gene[[gene_idx]]$transcript_ids)
        if (!is.na(iso_idx)) {
          if (is.null(results_per_gene[[gene_idx]]$delta_fdr)) {
            results_per_gene[[gene_idx]]$delta_fdr <- rep(NA, length(results_per_gene[[gene_idx]]$transcript_ids))
          }
          results_per_gene[[gene_idx]]$delta_fdr[iso_idx] <- px$fdr
        }
      }
    }
  }
  
  # Build all_transcript_stats table
  all_transcript_stats <- do.call(rbind, lapply(names(results_per_gene), function(gene) {
    res <- results_per_gene[[gene]]
    df <- data.frame(
      gene = gene,
      transcript_id = res$transcript_ids,
      pvalue = ifelse(is.null(res$delta_pvalue), NA_real_, res$delta_pvalue),
      fdr = ifelse(is.null(res$delta_fdr), NA_real_, res$delta_fdr),
      stringsAsFactors = FALSE
    )
    
    if (!is.null(res$lm_p_interaction)) {
      df$lm_p_interaction <- res$lm_p_interaction
    }
    if (!is.null(res$lm_adj_p_interaction)) {
      df$lm_adj_p_interaction <- res$lm_adj_p_interaction
    }
    
    return(df)
  }))
  
  # Build summary table with gene names
  summary_rows <- lapply(names(results_per_gene), function(gene) {
    res <- results_per_gene[[gene]]
    n_tx <- length(res$transcript_ids)
    n_switching <- sum(res$switching_status != "neutral", na.rm = TRUE)
    # Only compute max_delta if there are finite values
    delta_vals <- res$delta_influence[is.finite(res$delta_influence)]
    max_delta <- if (length(delta_vals) > 0) max(abs(delta_vals)) else 0
    n_fdr_sig <- sum(res$delta_fdr < 0.05, na.rm = TRUE)
    
    # Get gene name from rowData if available
    gene_name <- gene  # Default to gene ID
    gene_rows <- which(rowData(se)[[gene_col]] == gene)
    if (length(gene_rows) > 0 && "gene_name" %in% colnames(rowData(se))) {
      gene_name_from_data <- rowData(se)[gene_rows[1], "gene_name"]
      if (!is.na(gene_name_from_data)) {
        gene_name <- as.character(gene_name_from_data)
      }
    }
    
    data.frame(
      gene = gene,
      gene_name = gene_name,
      n_transcripts = n_tx,
      max_delta_influence = max_delta,
      n_switching_transcripts = n_switching,
      n_fdr_significant = n_fdr_sig,
      stringsAsFactors = FALSE
    )
  })
  
  summary_table <- do.call(rbind, summary_rows)
  rownames(summary_table) <- NULL
  
  # Build metadata
  n_fdr_sig_total <- if (nrow(all_transcript_stats) > 0) {
    sum(all_transcript_stats$fdr < 0.05, na.rm = TRUE)
  } else {
    0
  }
  
  metadata <- list(
    q = q,
    is_paired = is_paired,
    subject_col = if (is_paired) subject_col else NULL,
    pair_info = pair_info,
    norm = norm,
    log_base = log_base,
    pseudocount = pseudocount,
    threshold = threshold,
    n_transcripts_tested = nrow(all_transcript_stats),
    n_fdr_significant = n_fdr_sig_total,
    lm_results_provided = !is.null(lm_results),
    lm_p_threshold = lm_p_threshold,
    lm_genes_filtered = lm_genes_filtered,
    gene_processing_log = gene_processing_log
  )
  
  # Create result object
  # Build gene_id and gene_name vectors for analyzed genes
  analyzed_genes <- names(results_per_gene)
  
  # gene_names currently contains gene IDs (from rowData gene_col)
  # Add explicit gene_id vector for clarity
  result_gene_ids <- analyzed_genes
  
  # Add gene_names vector (mapped from gene_id_to_name)
  result_gene_names <- unname(gene_id_to_name[analyzed_genes])
  if (length(result_gene_names) == 0) {
    result_gene_names <- rep(NA_character_, length(analyzed_genes))
  }
  
  result <- list(
    gene_names = analyzed_genes,  # Keep for backwards compatibility (actually gene IDs)
    gene_ids = result_gene_ids,   # Explicit gene IDs
    gene_name_map = result_gene_names,  # Mapping to gene symbols/names
    conditions = conditions,
    results_per_gene = results_per_gene,
    summary_table = summary_table,
    all_transcript_stats = all_transcript_stats,
    metadata = metadata
  )
  
  class(result) <- c("tsenat_isoform_switching", "list")
  
  # Print results if requested
  if (print_results) {
    message("Isoform Switching Analysis Results")
    message("===================================")
    message(sprintf("Conditions: '%s' vs. '%s'", conditions[1], conditions[2]))
    
    if (is_paired && !is.null(pair_info)) {
      message(sprintf("Design: PAIRED (%s) - %d matched pairs", subject_col, pair_info$n_pairs))
    } else {
      message("Design: UNPAIRED")
    }
    
    if (!is.null(lm_results)) {
      message(sprintf("LM filtering: Genes with p < %.6f (N = %d)", lm_p_threshold, lm_genes_filtered))
    }
    
    message(sprintf("\nGenes analyzed: %d", length(result$gene_names)))
    message(sprintf("Total transcripts tested: %d", result$metadata$n_transcripts_tested))
    message(sprintf("FDR-significant transcripts (FDR<0.05): %d", result$metadata$n_fdr_significant))
    
    message("Summary Table:")
    message(paste(capture.output(print(summary_table)), collapse = ""))
    
    message("\n[OK] Use results$results_per_gene$'GeneName' to access per-gene switching details")
    message("[OK] Use results$summary_table for overview across genes")
    message("[OK] Use results$all_transcript_stats for FDR-corrected p-values per transcript")
    message("[OK] Use results$metadata$is_paired to check if paired design was applied")
    
    if (!is.null(lm_results)) {
      message("[OK] Access lm_p_interaction in each gene$lm_p_interaction for LM test results")
    }
    message("")
  }
  
  return(invisible(result))
}
