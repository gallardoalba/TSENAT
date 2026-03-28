#' Computation Statistics Helpers for TSENAT Visualization
#'
#' This module provides centralized functions for data transformation and
#' statistical summary calculations used across all visualization functions.
#' These helpers enable consistent, reusable computations and improve
#' maintainability by removing duplication.
#'
#' @name compute_stats

#' @noRd
NULL

# ============================================================================
# DIVERSITY SPECTRUM COMPUTATION
# ============================================================================

#' Compute Diversity Spectrum Statistics
#'
#' Aggregates diversity measurements across q-values and groups.
#' Calculates median/mean and variability (IQR/SD) for each q-value.
#'
#' @param se A \code{SummarizedExperiment} with diversity assays.
#' @param q_values Numeric vector of q-values to compute (optional, auto-detect if NULL).
#' @param metric Character: "median" (default) or "mean" for central tendency.
#' @param variability_metric Character: "iqr" (default) or "sd" for spread.
#' @param condition_col Character: column name for grouping conditions (optional).
#'
#' @return Data frame with columns:
#'   - q: q-value
#'   - group: condition group (if condition_col provided)
#'   - central: median or mean divergence
#'   - spread: IQR or SD of divergence
#'   - count: number of valid measurements
#'

#' @noRd

.compute_diversity_spectrum <- function(se,
                                       q_values = NULL,
                                       metric = c("median", "mean"),
                                       variability_metric = c("iqr", "sd"),
                                       condition_col = NULL) {

  require_pkgs(c("SummarizedExperiment", "dplyr"))

  # Validate input
  if (!inherits(se, "SummarizedExperiment")) {
    stop("se must be a SummarizedExperiment object", call. = FALSE)
  }

  if (nrow(se) == 0 || ncol(se) == 0) {
    stop("SummarizedExperiment is empty", call. = FALSE)
  }

  # Match arguments
  metric <- match.arg(metric)
  variability_metric <- match.arg(variability_metric)

  # Prepare long format data
  long_data <- .prepare_tsallis_long(se,
    assay_name = "diversity",
    condition_col = condition_col
  )

  if (nrow(long_data) == 0) {
    stop("No valid diversity data found in SummarizedExperiment", call. = FALSE)
  }

  # Ensure q is numeric
  long_data$q <- as.numeric(as.character(long_data$q))

  # Compute statistics by group and q-value
  if (!is.null(condition_col) && condition_col %in% colnames(long_data)) {
    # Group by condition
    stats <- long_data %>%
      dplyr::group_by(group, q) %>%
      dplyr::summarise(
        central = if (metric == "median") {
          median(.data$tsallis, na.rm = TRUE)
        } else {
          mean(.data$tsallis, na.rm = TRUE)
        },
        spread = if (variability_metric == "iqr") {
          IQR(.data$tsallis, na.rm = TRUE)
        } else {
          sqrt(stats::var(.data$tsallis, na.rm = TRUE))
        },
        count = sum(!is.na(.data$tsallis)),
        .groups = "drop"
      )
  } else {
    # No grouping
    stats <- long_data %>%
      dplyr::group_by(q) %>%
      dplyr::summarise(
        central = if (metric == "median") {
          median(.data$tsallis, na.rm = TRUE)
        } else {
          mean(.data$tsallis, na.rm = TRUE)
        },
        spread = if (variability_metric == "iqr") {
          IQR(.data$tsallis, na.rm = TRUE)
        } else {
          sqrt(stats::var(.data$tsallis, na.rm = TRUE))
        },
        count = sum(!is.na(.data$tsallis)),
        .groups = "drop"
      )
  }

  return(stats)
}


# ============================================================================
# GENE FILTERING & RANKING
# ============================================================================

#' Select Top Genes by P-Value
#'
#' Ranks genes by statistical significance and selects top N.
#'
#' @param results Data frame with at least one p-value column.
#' @param p_col Character: column name for p-values
#'   ("adj_p_interaction", "p_interaction", "padj", "pvalue").
#' @param gene_col Character: column name for gene identifiers
#'   ("gene_id", "gene", "gene_name").
#' @param n_genes Integer: number of top genes to select (default: 4).
#'
#' @return Character vector of top gene IDs, sorted by p-value (smallest first).
#'

#' @noRd

.select_top_genes <- function(results,
                             p_col = NULL,
                             gene_col = NULL,
                             n_genes = 4) {

  require_pkgs("dplyr")

  if (!is.data.frame(results) || nrow(results) == 0) {
    stop("results must be a non-empty data frame", call. = FALSE)
  }

  # Auto-detect p-value column
  if (is.null(p_col)) {
    candidate_cols <- c("adj_p_interaction", "p_interaction", "padj", "pvalue")
    matched <- candidate_cols[candidate_cols %in% colnames(results)]
    if (length(matched) > 0) {
      p_col <- matched[1]
    } else {
      stop("Could not find p-value column. ",
        "Provide p_col explicitly.",
        call. = FALSE
      )
    }
  }

  # Auto-detect gene column
  if (is.null(gene_col)) {
    candidate_cols <- c("gene_id", "gene", "gene_name")
    matched <- candidate_cols[candidate_cols %in% colnames(results)]
    if (length(matched) > 0) {
      gene_col <- matched[1]
    } else {
      stop("Could not find gene column. ",
        "Provide gene_col explicitly.",
        call. = FALSE
      )
    }
  }

  # Select top genes
  top_genes <- results %>%
    dplyr::arrange(.data[[p_col]]) %>%
    dplyr::slice(seq_len(min(n_genes, nrow(results)))) %>%
    dplyr::pull(.data[[gene_col]])

  return(as.character(top_genes))
}

#' Filter Genes by Significance Threshold
#'
#' Selects genes with p-value below threshold.
#'
#' @param results Data frame with p-values and gene identifiers.
#' @param p_threshold Numeric: p-value cutoff (default: 0.05).
#' @param p_col Character: p-value column name (auto-detected if NULL).
#' @param gene_col Character: gene identifier column (auto-detected if NULL).
#'
#' @return Character vector of significant gene IDs.
#'

#' @noRd

.filter_genes_by_pvalue <- function(results,
                                   p_threshold = 0.05,
                                   p_col = NULL,
                                   gene_col = NULL) {

  require_pkgs("dplyr")

  if (!is.data.frame(results) || nrow(results) == 0) {
    stop("results must be a non-empty data frame", call. = FALSE)
  }

  # Auto-detect columns (same logic as select_top_genes)
  if (is.null(p_col)) {
    candidate_cols <- c("adj_p_interaction", "p_interaction", "padj", "pvalue")
    matched <- candidate_cols[candidate_cols %in% colnames(results)]
    if (length(matched) > 0) {
      p_col <- matched[1]
    } else {
      stop("Could not find p-value column", call. = FALSE)
    }
  }

  if (is.null(gene_col)) {
    candidate_cols <- c("gene_id", "gene", "gene_name")
    matched <- candidate_cols[candidate_cols %in% colnames(results)]
    if (length(matched) > 0) {
      gene_col <- matched[1]
    } else {
      stop("Could not find gene column", call. = FALSE)
    }
  }

  # Filter and return
  sig_genes <- results %>%
    dplyr::filter(.data[[p_col]] < p_threshold) %>%
    dplyr::arrange(.data[[p_col]]) %>%
    dplyr::pull(.data[[gene_col]])

  return(as.character(sig_genes))
}

# ============================================================================
# DATA VALIDATION & QUALITY CHECKS
# ============================================================================

#' Validate Diversity SummarizedExperiment
#'
#' Checks that SE has required structure for diversity visualization.
#'
#' @param se A \code{SummarizedExperiment} object.
#' @param check_metadata Logical: also validate metadata? (default: TRUE)
#'
#' @return Logical TRUE if valid, else error with message.
#'

#' @noRd

.validate_diversity_se <- function(se, check_metadata = TRUE) {

  if (!inherits(se, "SummarizedExperiment")) {
    stop("se must be a SummarizedExperiment object", call. = FALSE)
  }

  if (nrow(se) == 0) {
    stop("SummarizedExperiment has no rows (samples)", call. = FALSE)
  }

  if (ncol(se) == 0) {
    stop("SummarizedExperiment has no columns (genes)", call. = FALSE)
  }

  # Check for diversity assay
  assay_names <- SummarizedExperiment::assayNames(se)
  if (!("diversity" %in% assay_names)) {
    stop("Required 'diversity' assay not found. ",
      "Available: ", paste(assay_names, collapse = ", "),
      call. = FALSE
    )
  }

  # Check for valid data
  div_mat <- SummarizedExperiment::assay(se, "diversity")
  if (all(is.na(div_mat))) {
    stop("All diversity values are NA", call. = FALSE)
  }

  if (check_metadata) {
    # Check for at least one q-value
    meta <- S4Vectors::metadata(se)
    if (!("q" %in% names(meta)) || length(meta$q) == 0) {
      warning("q-values not found in SE metadata", call. = FALSE)
    }
  }

  return(TRUE)
}

#' Validate Results Data Frame for Gene Selection
#'
#' Checks that results DataFrame has required columns.
#'
#' @param results Data frame (LM results, effect sizes, etc.).
#' @param require_pvalue Logical: check for p-value column? (default: TRUE)
#'
#' @return Logical TRUE if valid, else error.
#'

#' @noRd

.validate_results_df <- function(results, require_pvalue = TRUE) {

  if (!is.data.frame(results)) {
    stop("results must be a data frame", call. = FALSE)
  }

  if (nrow(results) == 0) {
    stop("results data frame is empty", call. = FALSE)
  }

  # Check for gene column
  gene_cols <- c("gene_id", "gene", "gene_name")
  has_gene <- any(gene_cols %in% colnames(results))
  if (!has_gene) {
    stop("No gene identifier column found. ",
      "Expected one of: ", paste(gene_cols, collapse = ", "),
      call. = FALSE
    )
  }

  # Check for p-value column
  if (require_pvalue) {
    p_cols <- c("adj_p_interaction", "p_interaction", "padj", "pvalue")
    has_pval <- any(p_cols %in% colnames(results))
    if (!has_pval) {
      stop("No p-value column found. ",
        "Expected one of: ", paste(p_cols, collapse = ", "),
        call. = FALSE
      )
    }
  }

  return(TRUE)
}

# ============================================================================
# FORMATTING & UTILITY FUNCTIONS
# ============================================================================

#' Format P-Value for Display
#'
#' Converts p-value to formatted string (scientific or threshold).
#'
#' @param pval Numeric p-value.
#' @param threshold Numeric: cutoff for "< threshold" format (default: 0.001).
#' @param digits Integer: decimal places for scientific notation (default: 2).
#'
#' @return Character string formatted p-value.
#'

#' @noRd

.format_pvalue <- function(pval, threshold = 0.001, digits = 2) {

  if (is.na(pval)) {
    return("NA")
  }

  if (pval < threshold) {
    return(paste0("< ", threshold))
  }

  return(format(
    pval,
    scientific = TRUE,
    digits = digits
  ))
}

#' Format Q-Value Label
#'
#' Converts numeric q-value to display label (e.g., "q = 1.0").
#'
#' @param q_val Numeric q-value.
#' @param prefix Character: prefix for label (default: "q").
#'
#' @return Character string label.
#'

#' @noRd

.format_q_label <- function(q_val, prefix = "q") {
  if (is.na(q_val)) {
    return("NA")
  }
  return(sprintf("%s = %.2f", prefix, as.numeric(q_val)))
}

#' Format Label for Display
#'
#' Converts underscored/raw column names to readable labels.
#' Replaces underscores with spaces and formats capitalization.
#'
#' @param lbl Character: label to format (may contain underscores).
#'
#' @return Character string, properly capitalized.
#'

#' @noRd

.format_label <- function(lbl) {
  if (is.null(lbl)) {
    return(NULL)
  }
  s <- gsub("_", " ", lbl)
  s <- gsub("\\s+", " ", s)
  s <- trimws(s)
  s <- tolower(s)
  if (nchar(s) == 0) {
    return(s)
  }
  if (nchar(s) == 1) {
    return(toupper(s))
  }
  paste0(toupper(substr(s, 1, 1)), substr(s, 2, nchar(s)))
}

# ============================================================================
# TRANSCRIPT PLOTTING HELPERS: Data Preparation
# ============================================================================

#' Prepare Inputs for Transcript-Level Plotting
#'
#' Normalizes and validates counts, samples, and tx2gene mapping.
#' Creates aggregation function based on chosen metric.
#'
#' @param counts Matrix or data.frame with transcripts as rows, samples as columns.
#'   Can also be a \code{SummarizedExperiment}.
#' @param readcounts Character: name of assay in SE (if counts is SE). Default: NULL.
#' @param samples Character vector: sample group assignments (optional).
#' @param coldata Character/data.frame: sample metadata (optional).
#' @param condition_col Character: column name for grouping
#'   in coldata (default: "sample_type").
#' @param tx2gene data.frame/character: Transcript-to-gene mapping with columns
#'   "Transcript" and "Gen". Can be file path or data.frame.
#' @param res Optional data.frame with results (gene names and p-values).
#' @param top_n Integer: number of transcripts to select.
#' @param pseudocount Numeric: pseudocount for log transformation (default: 0).
#' @param output_file Character: file path for saving plot (optional).
#' @param metric Character: aggregation metric
#'   ("median" [default], "mean", "variance", "iqr").
#'
#' @return List with elements:
#'   - counts: normalized count matrix
#'   - samples: sample group assignments
#'   - mapping: tx2gene data.frame
#'   - metric_choice: chosen metric
#'   - agg_fun: aggregation function
#'   - agg_label_unique: metric label for display
#'   - top_n: number of transcripts
#'   - pseudocount: pseudocount value
#'   - output_file: output file path (if provided)
#'

#' @noRd

.prepare_transcript_inputs <- function(counts,
                                      readcounts = NULL,
                                      samples = NULL,
                                      coldata = NULL,
                                      condition_col = "sample_type",
                                      tx2gene = NULL,
                                      res = NULL,
                                      top_n = NULL,
                                      pseudocount = 0,
                                      output_file = NULL,
                                      metric = c("median", "mean", "variance", "iqr")) {

  require_pkgs(c("SummarizedExperiment", "S4Vectors"))

  # Handle SummarizedExperiment input
  if (inherits(counts, "SummarizedExperiment")) {
    se <- counts
    counts_mat <- .get_readcounts_from_se(se, readcounts)
    counts <- as.matrix(counts_mat)
    samples <- .infer_samples_from_se(se, samples, condition_col = condition_col)

    if (is.null(tx2gene)) {
      txres <- .get_tx2gene_from_se(se, counts)
      if (!is.null(txres) && !is.null(txres$mapping)) {
        mapping <- data.frame(
          Transcript = rownames(counts),
          Gen = as.character(txres$mapping),
          stringsAsFactors = FALSE
        )
        tx2gene <- mapping
      }
    }
  }

  # Validate counts
  if (!is.matrix(counts) && !is.data.frame(counts)) {
    stop("`counts` must be a matrix, data.frame, or SummarizedExperiment",
      call. = FALSE
    )
  }

  counts <- as.matrix(counts)
  if (is.null(rownames(counts))) {
    stop("`counts` must have rownames (transcript identifiers)", call. = FALSE)
  }

  # Infer samples from coldata if needed
  if (is.null(samples)) {
    if (!is.null(coldata)) {
      samples <- .infer_samples_from_coldata(coldata, counts, condition_col)
    } else {
      stop("Either 'samples' or 'coldata' must be provided", call. = FALSE)
    }
  }

  # Validate and normalize tx2gene
  if (is.null(tx2gene)) {
    stop("`tx2gene` must be provided", call. = FALSE)
  }
  mapping <- .read_tx2gene(tx2gene)

  if (length(samples) != ncol(counts)) {
    stop("Length of `samples` must equal columns in `counts`", call. = FALSE)
  }

  # Create aggregation function
  metric_choice <- match.arg(metric)
  agg_fun <- switch(metric_choice,
    median = function(x) stats::median(x, na.rm = TRUE),
    mean = function(x) base::mean(x, na.rm = TRUE),
    variance = function(x) stats::var(x, na.rm = TRUE),
    iqr = function(x) stats::IQR(x, na.rm = TRUE)
  )

  agg_label_metric <- if (metric_choice == "iqr") "IQR" else metric_choice
  agg_label <- sprintf("Transcript-level expression with metric %s", agg_label_metric)
  agg_label_unique <- agg_label

  list(
    counts = counts,
    samples = samples,
    mapping = mapping,
    metric_choice = metric_choice,
    agg_fun = agg_fun,
    agg_label_unique = agg_label_unique,
    top_n = top_n,
    pseudocount = pseudocount,
    output_file = output_file
  )
}

#' Read and Validate tx2gene Mapping
#'
#' Reads transcript-to-gene mapping from file or data.frame.
#' Validates required columns: "Transcript" and "Gen".
#'
#' @param tx2gene Character (file path) or data.frame mapping.
#'
#' @return data.frame with columns "Transcript" and "Gen".
#'

#' @noRd

.read_tx2gene <- function(tx2gene) {
  if (is.null(tx2gene)) {
    stop("`tx2gene` must be provided as file path or data.frame",
      call. = FALSE
    )
  }

  if (is.character(tx2gene) && length(tx2gene) == 1) {
    if (!file.exists(tx2gene)) {
      stop("tx2gene file not found: ", tx2gene, call. = FALSE)
    }
    mapping <- utils::read.delim(tx2gene, stringsAsFactors = FALSE, header = TRUE)
  } else if (is.data.frame(tx2gene)) {
    mapping <- tx2gene
  } else {
    stop("`tx2gene` must be file path or data.frame", call. = FALSE)
  }

  if (!all(c("Transcript", "Gen") %in% colnames(mapping))) {
    stop("tx2gene must have columns 'Transcript' and 'Gen'",
      call. = FALSE
    )
  }

  return(mapping)
}

#' Infer Samples from Column Metadata
#'
#' Extracts sample group assignments from coldata.
#' Aligns sample IDs from coldata to counts columns.
#'
#' @param coldata Character (file path) or data.frame with sample metadata.
#' @param counts Count matrix (for column name alignment).
#' @param condition_col Character: column name for grouping variable.
#'
#' @return Character vector of sample group assignments.
#'

#' @noRd

.infer_samples_from_coldata <- function(coldata, counts, condition_col) {
  if (is.character(coldata) && length(coldata) == 1) {
    if (!file.exists(coldata)) {
      stop("coldata file not found: ", coldata, call. = FALSE)
    }
    cdf <- utils::read.delim(coldata, header = TRUE, stringsAsFactors = FALSE)
  } else if (is.data.frame(coldata)) {
    cdf <- coldata
  } else {
    stop("`coldata` must be file path or data.frame", call. = FALSE)
  }

  # Try row-indexed matching first
  if (!is.null(rownames(cdf)) && all(colnames(counts) %in% rownames(cdf))) {
    return(as.character(cdf[colnames(counts), condition_col]))
  }

  # Try sample ID column matching
  sample_id_cols <- c("sample", "Sample", "sample_id", "id")
  sid <- intersect(sample_id_cols, colnames(cdf))

  if (length(sid) > 0) {
    sid <- sid[1]
    if (!all(colnames(counts) %in% as.character(cdf[[sid]]))) {
      stop("coldata sample ID column doesn't match counts columns",
        call. = FALSE
      )
    }
    row_ix <- match(colnames(counts), as.character(cdf[[sid]]))
    return(as.character(cdf[[condition_col]][row_ix]))
  }

  stop("Could not match coldata to counts. Provide rownames or sample ID column.",
    call. = FALSE
  )
}

#' Create Aggregation Function
#'
#' Builds an aggregation function based on chosen metric.
#'
#' @param metric Character: "median" (default), "mean", "variance", or "iqr".
#'
#' @return List with:
#'   - metric_choice: the selected metric
#'   - agg_fun: function that computes the metric
#'   - agg_label_unique: display label
#'

#' @noRd

.create_aggregation_function <- function(metric = c("median", "mean", "variance", "iqr")) {
  metric_choice <- match.arg(metric)

  agg_fun <- switch(metric_choice,
    median = function(x) stats::median(x, na.rm = TRUE),
    mean = function(x) base::mean(x, na.rm = TRUE),
    variance = function(x) stats::var(x, na.rm = TRUE),
    iqr = function(x) stats::IQR(x, na.rm = TRUE)
  )

  agg_label_metric <- if (metric_choice == "iqr") "IQR" else metric_choice
  agg_label <- sprintf("Transcript-level expression with metric %s", agg_label_metric)
  agg_label_unique <- agg_label

  list(
    metric_choice = metric_choice,
    agg_fun = agg_fun,
    agg_label_unique = agg_label_unique
  )
}

#' Build Long-Format Transcript Data
#'
#' Transforms wide count matrix to long-format data.frame
#' for ggplot visualization.
#'
#' @param gene_single Character: single gene identifier.
#' @param mapping data.frame: tx2gene mapping with Transcript and Gen columns.
#' @param counts Matrix: transcript count matrix.
#' @param samples Character vector: sample group assignments.
#' @param top_n Integer: limit to top N transcripts (optional).
#'
#' @return List with:
#'   - df_long: long-format data.frame (columns: tx, sample, expr, group)
#'   - txs: selected transcript identifiers
#'

#' @noRd

.build_transcript_long <- function(gene_single, mapping, counts, samples, top_n = NULL) {
  txs <- mapping$Transcript[mapping$Gen == gene_single]
  txs <- intersect(txs, rownames(counts))

  if (length(txs) == 0) {
    stop("No transcripts found for gene: ", gene_single, call. = FALSE)
  }

  if (!is.null(top_n)) {
    txs <- head(txs, top_n)
  }

  # Create long-format data
  mat <- counts[txs, , drop = FALSE]
  df_all <- as.data.frame(mat)
  df_all$tx <- rownames(mat)

  require_pkgs("tidyr")
  df_long <- tidyr::pivot_longer(df_all, -tx, names_to = "sample", values_to = "expr")
  df_long$group <- rep(samples, times = length(txs))

  list(df_long = df_long, txs = txs)
}

#' Aggregate Long-Format Transcript Data
#'
#' Summarizes expression by transcript and group.
#'
#' @param df_long Long-format data.frame from \code{build_transcript_long}.
#' @param agg_fun Function: aggregation function (e.g., median, mean).
#' @param pseudocount Numeric: pseudocount for log transformation.
#'
#' @return data.frame with columns: tx, group, expr, log2expr.
#'

#' @noRd

.aggregate_transcript_data <- function(df_long, agg_fun, pseudocount = 0) {
  df_summary <- stats::aggregate(expr ~ tx + group, data = df_long, FUN = agg_fun)
  df_summary$log2expr <- log2(df_summary$expr + pseudocount)
  df_summary$tx <- factor(df_summary$tx, levels = unique(df_summary$tx))
  return(df_summary)
}

#' Select Top Genes from Results  
#'
#' Extracts top genes from results by p-value.
#'
#' @param res data.frame: results with gene and p-value columns.
#' @param top_n Integer: number of genes to select.
#'
#' @return Character vector of top gene IDs.
#'

#' @noRd

.select_genes_from_results <- function(res, top_n) {
  if (is.null(res)) {
    stop("Either 'gene' or 'res' must be provided", call. = FALSE)
  }

  if (!("genes" %in% colnames(res))) {
    stop("res must contain a 'genes' column", call. = FALSE)
  }

  # Find p-value column (ordered by preference)
  p_cols <- c("padj", "adjusted_p_values", "pvalue", "raw_p_values")
  p_col <- intersect(p_cols, colnames(res))[1]

  if (is.na(p_col)) {
    # No p-value column; just return gene order
    ord <- seq_len(nrow(res))
  } else {
    ord <- order(res[[p_col]], na.last = NA)
  }

  genes_sel <- as.character(res$genes[ord])
  genes_sel <- unique(genes_sel)
  return(head(genes_sel, top_n))
}

# ============================================================================
# PLOT TSALLIS Q-CURVE HELPERS
# ============================================================================

#' Convert TSENATAnalysis to combined SummarizedExperiment
#'
#' @param analysis TSENATAnalysis object with diversity_results
#' @return SummarizedExperiment with combined assay across all q-values

#' @noRd
.tsenat_prepare_combined_se <- function(analysis) {
  require_pkgs(c("SummarizedExperiment", "S4Vectors"))
  
  div_list <- analysis@diversity_results
  
  # Extract first SE to get dimensions
  first_se <- NULL
  combined_assays_dict <- list()
  
  for (q_name in names(div_list)) {
    obj <- div_list[[q_name]]
    if (methods::is(obj, "SummarizedExperiment")) {
      mat <- SummarizedExperiment::assay(obj, 1)
      if (is.null(first_se)) {
        first_se <- obj
      }
    } else {
      mat <- as.matrix(obj)
    }
    
    q_val <- as.numeric(sub("^q_", "", q_name))
    combined_assays_dict[[q_name]] <- list(matrix = mat, q_val = q_val)
  }
  
  if (is.null(first_se)) {
    stop("No valid SummarizedExperiment found in analysis@diversity_results")
  }
  
  # Get dimensions
  target_genes <- rownames(first_se)
  target_n_cols <- ncol(first_se)
  target_n_qs <- length(combined_assays_dict)
  total_cols <- target_n_cols * target_n_qs
  
  # Create combined assay matrix
  combined_assay <- matrix(0, nrow = length(target_genes), ncol = total_cols)
  rownames(combined_assay) <- target_genes
  
  combined_coldata_list <- list()
  col_idx <- 1
  
  for (q_name in names(combined_assays_dict)) {
    mat <- combined_assays_dict[[q_name]]$matrix
    q_val <- combined_assays_dict[[q_name]]$q_val
    
    # Handle dimension mismatches
    if (ncol(mat) != target_n_cols) {
      if (ncol(mat) > target_n_cols) {
        mat <- mat[, seq_len(target_n_cols), drop = FALSE]
      } else {
        mat <- cbind(mat, matrix(0, nrow = nrow(mat), ncol = target_n_cols - ncol(mat)))
      }
    }
    
    # Reorder rows to match first_se
    mat <- mat[target_genes, , drop = FALSE]
    
    # Add q-value suffix to column names
    orig_colnames <- colnames(mat)
    if (is.null(orig_colnames)) {
      orig_colnames <- paste0("sample_", seq_len(ncol(mat)))
    }
    
    clean_colnames <- sub("_q=.*$", "", orig_colnames)
    if (is.na(clean_colnames[1]) || identical(clean_colnames, orig_colnames)) {
      clean_colnames <- orig_colnames
    }
    
    unique_colnames <- paste0(clean_colnames, "_q=", formatC(q_val, format = "f", digits = 3))
    
    # Fill in combined assay
    if (col_idx + ncol(mat) - 1 > total_cols) {
      stop("Dimension mismatch: ", col_idx, " to ", col_idx + ncol(mat) - 1,
           " exceeds total_cols=", total_cols)
    }
    
    for (i in seq_len(ncol(mat))) {
      combined_assay[, col_idx] <- mat[, i]
      col_idx <- col_idx + 1
    }
    
    # Build colData for this q-value
    if (is(div_list[[q_name]], "SummarizedExperiment")) {
      cd <- as.data.frame(SummarizedExperiment::colData(div_list[[q_name]]))
    } else {
      cd <- data.frame(row.names = unique_colnames)
    }
    cd$q <- q_val
    rownames(cd) <- unique_colnames
    combined_coldata_list[[q_name]] <- cd
  }
  
  # Combine colData
  combined_coldata_df <- do.call(rbind, combined_coldata_list)
  colnames(combined_assay) <- rownames(combined_coldata_df)
  
  # Get/create rowData
  rd_combined <- tryCatch({
    rd_temp <- SummarizedExperiment::rowData(first_se)
    if (!is.null(rd_temp) && nrow(rd_temp) > 0) {
      rd_temp
    } else {
      NULL
    }
  }, error = function(e) NULL)
  
  if (is.null(rd_combined) || nrow(rd_combined) == 0) {
    rd_combined <- data.frame(
      gene_id = rownames(combined_assay),
      row.names = rownames(combined_assay),
      stringsAsFactors = FALSE
    )
  }
  
  # Return combined SE
  SummarizedExperiment::SummarizedExperiment(
    assays = list(diversity = combined_assay),
    colData = combined_coldata_df,
    rowData = rd_combined
  )
}

#' Compute gene-level statistics (median +/- SD) by group and q-value
#'
#' @param long_data Long-format data frame with Gene, q, group, tsallis columns
#' @return Data frame with central tendency and spread by gene, group, q

#' @noRd
.tsenat_compute_gene_group_stats <- function(long_data) {
  require_pkgs("dplyr")
  
  long_data$qnum <- as.numeric(as.character(long_data$q))
  
  dplyr::summarise(
    dplyr::group_by(long_data, group, qnum),
    central = median(tsallis, na.rm = TRUE),
    spread = sqrt(stats::var(tsallis, na.rm = TRUE)),
    .groups = "drop"
  )
}

#' Aggregate bootstrap CI bounds by group and q-value
#'
#' @param se SummarizedExperiment with ci_lower and ci_upper assays
#' @param long Long-format data with group, q, sample, tsallis
#' @return Data frame with q, median, ci_lower, ci_upper, group

#' @noRd
.tsenat_bootstrap_aggregate_ci <- function(se, long) {
  require_pkgs(c("SummarizedExperiment", "dplyr"))
  
  ci_lower_mat <- SummarizedExperiment::assay(se, "ci_lower")
  ci_upper_mat <- SummarizedExperiment::assay(se, "ci_upper")
  
  sample_names <- colnames(ci_lower_mat)
  if (is.null(sample_names)) {
    sample_names <- paste0("Sample", seq_len(ncol(ci_lower_mat)))
  }
  
  groups <- unique(sort(long$group))
  unique_q <- sort(unique(long$q))
  
  plot_df <- data.frame(
    q = numeric(),
    median = numeric(),
    ci_lower = numeric(),
    ci_upper = numeric(),
    group = character(),
    stringsAsFactors = FALSE
  )
  
  for (group_val in groups) {
    for (q_val in unique_q) {
      group_q_data <- long %>%
        dplyr::filter(group == group_val, q == q_val)
      
      if (nrow(group_q_data) > 0) {
        median_val <- median(group_q_data$tsallis, na.rm = TRUE)
        
        group_samples <- unique(group_q_data$sample)
        all_ci_lower <- c()
        all_ci_upper <- c()
        
        for (samp in group_samples) {
          samp_idx <- which(sample_names == samp)
          if (length(samp_idx) > 0) {
            all_ci_lower <- c(all_ci_lower, mean(ci_lower_mat[, samp_idx], na.rm = TRUE))
            all_ci_upper <- c(all_ci_upper, mean(ci_upper_mat[, samp_idx], na.rm = TRUE))
          }
        }
        
        if (length(all_ci_lower) > 0) {
          ci_lower_final <- median(all_ci_lower, na.rm = TRUE)
          ci_upper_final <- median(all_ci_upper, na.rm = TRUE)
        } else {
          ci_lower_final <- median(ci_lower_mat, na.rm = TRUE)
          ci_upper_final <- median(ci_upper_mat, na.rm = TRUE)
        }
        
        plot_df <- rbind(plot_df, data.frame(
          q = q_val,
          median = median_val,
          ci_lower = ci_lower_final,
          ci_upper = ci_upper_final,
          group = group_val,
          stringsAsFactors = FALSE
        ))
      }
    }
  }
  
  plot_df
}

# ============================================================================
# GAM INTERACTION HELPERS
# ============================================================================

#' Build sample-to-group mapping from colData
#'
#' @param cdata SummarizedExperiment colData with sample metadata
#' @param condition_col Column name for group assignments
#' @return Named character vector: sample name -> group value

#' @noRd
.tsenat_prepare_sample_group_mapping <- function(cdata, condition_col) {
  coldata_rownames <- rownames(cdata)
  coldata_sample_names <- sub("_q=.*", "", coldata_rownames)
  
  unique_samples <- unique(coldata_sample_names)
  sample_to_group <- character(length(unique_samples))
  names(sample_to_group) <- unique_samples
  
  for (samp in unique_samples) {
    idx <- which(coldata_sample_names == samp)[1]
    sample_to_group[samp] <- as.character(cdata[[condition_col]][idx])
  }
  
  sample_to_group
}

#' Prepare long-format plot data for a single gene
#'
#' @param gene Gene ID to extract
#' @param mat Assay matrix (genes x samples*q)
#' @param sample_to_group Named vector mapping sample names to groups
#' @return Data frame with columns: sample, group, q, entropy (or NULL if invalid)

#' @noRd
.tsenat_plot_gam_prepare_gene_data <- function(gene, mat, sample_to_group) {
  if (!(gene %in% rownames(mat))) {
    return(NULL)
  }
  
  gene_vals <- mat[gene, ]
  col_names_full <- colnames(mat)
  
  # Parse column names: "Sample_q=value"
  col_sample_names <- sub("_q=.*", "", col_names_full)
  col_q_values <- as.numeric(sub(".*_q=", "", col_names_full))
  
  # Look up group for each column
  col_groups <- unname(sample_to_group[col_sample_names])
  
  if (any(is.na(col_groups))) {
    return(NULL)
  }
  
  # Build long-format data frame
  plot_df <- data.frame(
    sample = col_sample_names,
    group = col_groups,
    q = col_q_values,
    entropy = as.numeric(gene_vals),
    stringsAsFactors = FALSE
  )
  
  # Remove NA entries
  plot_df <- plot_df[!is.na(plot_df$entropy), , drop = FALSE]
  
  if (nrow(plot_df) == 0) {
    return(NULL)
  }
  
  plot_df
}

#' Fit GAM models per group and generate predictions
#'
#' @param plot_df Long-format data frame with sample, group, q, entropy
#' @return List with $plot_data and $pred_data data frames (or NULL if fitting fails)

#' @noRd
.tsenat_plot_gam_fit_group <- function(plot_df) {
  require_pkgs(c("mgcv", "dplyr"))
  
  unique_groups <- unique(plot_df$group)
  
  if (length(unique_groups) < 2) {
    return(NULL)
  }
  
  # Generate prediction grid
  q_range <- range(plot_df$q, na.rm = TRUE)
  if (!is.finite(q_range[1]) || !is.finite(q_range[2])) {
    return(NULL)
  }
  
  pred_q <- seq(q_range[1], q_range[2], length.out = 100)
  
  # Fit GAM and predict for each group
  pred_list <- list()
  for (gr in unique_groups) {
    subset_data <- subset(plot_df, group == gr)
    
    if (nrow(subset_data) < 3) {
      next
    }
    
    tryCatch(
      {
        k <- min(10, max(2, round(nrow(subset_data) / 2)))
        gam_fit <- mgcv::gam(entropy ~ s(q, k = k), data = subset_data)
        
        pred_data <- data.frame(q = pred_q)
        pred_vals <- stats::predict(gam_fit, newdata = pred_data, se.fit = TRUE)
        
        pred_list[[as.character(gr)]] <- data.frame(
          group = gr,
          q = pred_q,
          entropy_fit = pred_vals$fit,
          se = pred_vals$se.fit,
          stringsAsFactors = FALSE
        )
      },
      error = function(e) {
        # Silently skip failed fits
        NULL
      }
    )
  }
  
  if (length(pred_list) == 0) {
    return(NULL)
  }
  
  pred_df <- do.call(rbind, pred_list)
  
  # Ensure group is factor with consistent levels
  group_levels <- sort(unique(c(as.character(plot_df$group), as.character(pred_df$group))))
  plot_df$group <- factor(plot_df$group, levels = group_levels)
  pred_df$group <- factor(pred_df$group, levels = group_levels)
  
  list(plot_data = plot_df, pred_data = pred_df, group_levels = group_levels)
}

#' Select genes to plot based on significance
#'
#' @param lm_res Data frame with gene and p-value columns
#' @param genes Optional character vector of specific genes
#' @param n_top Number of top genes to select
#' @param sig_alpha Significance threshold
#' @return Character vector of gene IDs to plot (or NULL if none selected)

#' @noRd
.tsenat_plot_select_genes <- function(lm_res, genes = NULL, n_top = 6, sig_alpha = 0.05) {
  if (!is.null(genes)) {
    if (!is.character(genes)) {
      stop("genes must be a character vector of gene names", call. = FALSE)
    }
    return(genes)
  }
  
  # Select top genes by p-value
  if ("adj_p_interaction" %in% colnames(lm_res)) {
    sig_mask <- lm_res$adj_p_interaction <= sig_alpha
  } else if ("p_interaction" %in% colnames(lm_res)) {
    sig_mask <- lm_res$p_interaction <= sig_alpha
  } else {
    stop("lm_res must contain 'adj_p_interaction' or 'p_interaction' column", call. = FALSE)
  }
  
  sig_genes <- lm_res[sig_mask, , drop = FALSE]
  
  if (nrow(sig_genes) == 0) {
    return(NULL)
  }
  
  # Select top n
  sig_genes$gene[seq_len(min(n_top, nrow(sig_genes)))]
}
