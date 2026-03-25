# # Helpers for generating diversity plots (density, violin, MA)

#' @importFrom utils head
#' @importFrom grDevices png dev.off
if (getRversion() >= "2.15.1") {
    utils::globalVariables(
        c(
            "Gene",
            "diversity",
            "sample",
            "sample_q",
            "sample_type",
            "fold",
            "significant",
            "value",
            ".",
            "group",
            "tsallis",
            "q",
            "median",
            "IQR",
            "padj_num",
            "padj_clean",
            "xval",
            "label_flag",
            "genes",
            "mean",
            # transcript plotting globals
            "tx",
            "expr",
            "tx_cond",
            "log2expr",
            # generic plotting helpers
            "x",
            "y",
            "padj"
        )
    )

    require_pkgs <- function(pkgs) {
        missing <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
        if (length(missing)) {
            stop(
                sprintf(
                    "%s required",
                    paste(missing, collapse = ", ")
                )
            )
        }
        invisible(TRUE)
    }
}

# Small utility helpers to reduce repeated code paths when extracting
# samples, readcounts and tx->gene mappings from a
# SummarizedExperiment. Keeping these as focused helpers improves
# readability of the longer plotting functions below.
infer_samples_from_se <- function(se, samples = NULL, condition_col = "sample_type") {
    if (!is.null(samples)) {
        return(as.character(samples))
    }
    cd <- NULL
    try(cd <- SummarizedExperiment::colData(se), silent = TRUE)
    if (is.null(cd)) {
        return(NULL)
    }

    # Common column names to try
    candidates <- c(
        condition_col,
        "condition",
        "group",
        "sample_group",
        "sampleType",
        "class",
        "status",
        "phenotype"
    )
    for (nm in candidates) {
        if (nm %in% colnames(cd)) {
            return(as.character(cd[[nm]]))
        }
    }

    # Fallback: choose column with smallest >1 unique values
    uniq_counts <- vapply(cd, function(col) length(unique(na.omit(col))), integer(1))
    valid_cols <- names(uniq_counts[uniq_counts > 1])
    if (length(valid_cols) > 0) {
        bin_cols <- valid_cols[uniq_counts[valid_cols] == 2]
        pick <- if (length(bin_cols) > 0) bin_cols[1] else valid_cols[which.min(uniq_counts[valid_cols])]
        return(as.character(cd[[pick]]))
    }

    NULL
}

get_readcounts_from_se <- function(se, readcounts_arg = NULL) {
    # If user provided a readcounts object/path, accept it first
    if (!is.null(readcounts_arg)) {
        if (is.character(readcounts_arg) && length(readcounts_arg) == 1) {
            if (!file.exists(readcounts_arg)) stop("readcounts file not found: ", readcounts_arg)
            rc_df <- utils::read.delim(readcounts_arg, header = TRUE, stringsAsFactors = FALSE)
            if (!is.null(colnames(rc_df)) && ncol(rc_df) > 1) {
                counts <- as.matrix(rc_df[, -1, drop = FALSE])
                rownames(counts) <- rc_df[[1]]
            } else {
                counts <- as.matrix(rc_df)
            }
            return(counts)
        } else if (is.matrix(readcounts_arg) || is.data.frame(readcounts_arg)) {
            return(as.matrix(readcounts_arg))
        } else {
            stop("`readcounts` must be a matrix/data.frame or path to a file")
        }
    }

    md <- NULL
    try(md <- S4Vectors::metadata(se), silent = TRUE)
    if (!is.null(md) && !is.null(md$readcounts)) {
        return(as.matrix(md$readcounts))
    }

    assay_names <- SummarizedExperiment::assayNames(se)
    preferred_assays <- c("readcounts", "counts", "tx_counts", "counts_tx")
    chosen <- intersect(preferred_assays, assay_names)
    if (length(chosen) > 0) {
        return(as.matrix(SummarizedExperiment::assay(se, chosen[1])))
    }

    # fallback to first assay with a warning
    warning(
        "Using first assay from SummarizedExperiment to compute",
        " expression-based fold changes; ensure it contains",
        " transcript-level readcounts or provide metadata$readcounts"
    )
    as.matrix(SummarizedExperiment::assay(se))
}

get_tx2gene_from_se <- function(se, readcounts_mat = NULL) {
    md <- NULL
    try(md <- S4Vectors::metadata(se), silent = TRUE)
    # prefer explicit tx2gene in metadata
    if (!is.null(md) && !is.null(md$tx2gene) && is.data.frame(md$tx2gene)) {
        txmap <- md$tx2gene
        # attempt to find sensible columns
        tx_col <- if ("Transcript" %in% colnames(txmap)) "Transcript" else colnames(txmap)[1]
        gene_col <- if ("Gen" %in% colnames(txmap)) "Gen" else colnames(txmap)[2]
        return(list(
            type = "vector",
            mapping = as.character(txmap[[gene_col]][match(rownames(readcounts_mat), txmap[[tx_col]])])
        ))
    }

    # fallback: try rowData mapping
    rdata <- SummarizedExperiment::rowData(se)
    if (!is.null(rdata) && "genes" %in% colnames(rdata)) {
        genes_vec <- as.character(rdata$genes)
        if (!is.null(readcounts_mat) && length(genes_vec) == nrow(readcounts_mat)) {
            return(list(type = "vector", mapping = genes_vec))
        }
    }

    # last resort: use rownames of readcounts as gene identifiers
    if (!is.null(readcounts_mat)) {
        return(list(type = "vector", mapping = rownames(readcounts_mat)))
    }

    NULL
}

validate_control_in_samples <- function(control, samples) {
    uniq <- unique(samples)
    if (!is.null(control) && control %in% uniq) {
        return(control)
    }
    if ("Normal" %in% uniq) {
        return("Normal")
    }
    # fallback to first level and message
    chosen <- uniq[1]
    message(sprintf("`control` not found; using '%s' instead", chosen))
    chosen
}


# Core MA plotting implementation documentation moved to internal block
#' Plot MA using Tsallis-based fold changes
#'
#' Wrapper around `plot_ma(..., type = "tsallis")` for convenience and
#' clearer API separation.
#'
#' @param x Data.frame from `calculate_difference()`.
#' @param sig_alpha Numeric significance threshold for adjusted p-values (default: 0.05).
#' @param x_label Optional x-axis label passed to `plot_ma`.
#' @param y_label Optional y-axis label passed to `plot_ma`.
#' @param title Optional plot title passed to `plot_ma`.
#' @param ... Additional arguments passed to `plot_ma()`.
#' @return A `ggplot2` object representing the MA plot.
#' @noRd
#' @keywords internal
plot_ma_tsallis <- function(x, sig_alpha = 0.05, x_label = NULL, y_label = NULL, title = NULL, ...) {
    title_use <- title %||% "Tsallis-based MA plot"
    x_label_use <- x_label %||% "mean_difference"
    y_label_use <- y_label %||% "Log10 fold-change of entropy"
    .plot_ma_core(x, fc_df = NULL, sig_alpha = sig_alpha, x_label = x_label_use, y_label = y_label_use, title = title_use)
}



# Core MA plotting implementation used by wrappers above. Accepts a
# differential results `x` (data.frame) and an optional `fc_df` with
# fold-changes (genes as rownames or a `genes` column). Returns a
# `ggplot` MA-plot.
#' Core MA plotting implementation (internal)
#'
#' This is an internal helper used by `plot_ma_tsallis()`.
#' It is documented here for developers but is not exported.
#' @noRd
.plot_ma_core <- function(x,
                          fc_df = NULL,
                          diff_res = NULL,
                          sig_alpha = 0.05,
                          x_label = NULL,
                          y_label = NULL,
                          title = NULL,
                          ...) {
    require_pkgs(c("ggplot2"))

    df <- as.data.frame(x, stringsAsFactors = FALSE)
    # ensure a gene identifier column exists
    # Support both new "gene_id" and legacy "genes" column names
    if (!("genes" %in% colnames(df))) {
        if ("gene_id" %in% colnames(df)) {
            # Rename gene_id to genes for backward compatibility with plotting code
            df$genes <- df$gene_id
        } else if (!is.null(rownames(df))) {
            df$genes <- rownames(df)
        }
    }

    # If an external fc_df is provided, merge fold values into df
    if (!is.null(fc_df)) {
        fdf <- as.data.frame(fc_df, stringsAsFactors = FALSE)
        # Support both "gene_id" and "genes" column names
        if (!("genes" %in% colnames(fdf))) {
            if ("gene_id" %in% colnames(fdf)) {
                fdf$genes <- fdf$gene_id
            } else if (!is.null(rownames(fdf))) {
                fdf$genes <- rownames(fdf)
            }
        }
        if (!("log2_fold_change" %in% colnames(fdf))) {
            stop("Provided `fc_df` must contain 'log2_fold_change' column")
        }
        df <- merge(df, fdf[, c("genes", "log2_fold_change")], by = "genes", all.x = TRUE, suffixes = c("", ".fc"))
        # prefer fc_df fold values when present
        if ("log2_fold_change.fc" %in% colnames(df)) df$log2_fold_change <- ifelse(!is.na(df$log2_fold_change.fc), df$log2_fold_change.fc, df$log2_fold_change)
    }

    # Detect fold-change column
    fold_candidates <- c("log2_fold_change", "logFC", "fold", "estimate_interaction", "fold_change")
    fold_col <- intersect(fold_candidates, colnames(df))
    if (length(fold_col) == 0) stop("Could not find a fold-change column in input")
    fold_col <- fold_col[1]

    # Detect mean/average columns for x-axis
    mean_cols <- grep("_mean$|_median$", colnames(df), value = TRUE)
    if (length(mean_cols) >= 2) {
        # ensure the two chosen columns are consistent (both _mean or both _median)
        c1 <- mean_cols[1]
        c2 <- mean_cols[2]
        is_mean1 <- grepl("_mean$", c1)
        is_mean2 <- grepl("_mean$", c2)
        is_med1 <- grepl("_median$", c1)
        is_med2 <- grepl("_median$", c2)
        if (!((is_mean1 && is_mean2) || (is_med1 && is_med2))) {
            stop("Could not find two mean or two median columns")
        }
        xvals <- rowMeans(df[, mean_cols[seq_len(2)],
            drop = FALSE
        ], na.rm = TRUE)
        x_label <- x_label %||% paste0(mean_cols[1], " vs ", mean_cols[2])
    } else if (length(mean_cols) == 1) {
        xvals <- as.numeric(df[[mean_cols[1]]])
        x_label <- x_label %||% mean_cols[1]
    } else if ("mean" %in% colnames(df)) {
        xvals <- as.numeric(df$mean)
        x_label <- x_label %||% "Mean"
    } else {
        # fallback: use rank or index
        xvals <- seq_len(nrow(df))
        x_label <- x_label %||% "Index"
    }

    yvals <- as.numeric(df[[fold_col]])

    # p-value / adjusted p-value detection
    padj_candidates <- c("padj", "adjusted_p_values", "adj_p_value", "adj_p", "p.adjust")
    padj_col <- intersect(padj_candidates, colnames(df))
    padj_col <- if (length(padj_col)) padj_col[1] else NULL

    padj <- if (!is.null(padj_col)) as.numeric(df[[padj_col]]) else rep(1, length(yvals))
    padj[is.na(padj)] <- 1

    sig_flag <- ifelse(abs(yvals) > 0 & padj < sig_alpha, "significant", "non-significant")

    plot_df <- data.frame(genes = df$genes, x = xvals, y = yvals, padj = padj, significant = sig_flag, stringsAsFactors = FALSE)

    prep <- .tsenat_prepare_ma_plot_df(df, fold_col = fold_col, mean_cols = mean_cols, x_label = x_label, y_label = y_label)
    plot_df <- prep$plot_df
    x_label_formatted <- .tsenat_format_label(prep$x_label)
    y_label_raw <- prep$y_label %||% fold_col
    y_label_formatted <- .tsenat_format_label(y_label_raw)
    if (!is.null(y_label_formatted)) {
        y_label_formatted <- sub("\\blog2\\b", "log10", y_label_formatted, ignore.case = TRUE)
    }

    p <- ggplot2::ggplot(plot_df, ggplot2::aes(x = x, y = y, color = significant)) +
        ggplot2::geom_point(alpha = 0.75, size = 3.2) +
        ggplot2::scale_color_manual(
            values = .tsenat_significance_colors(),
            guide = "none"
        ) +
        ggplot2::labs(
            title = title %||% "MA plot: mean vs log10 fold-change",
            x = x_label_formatted,
            y = y_label_formatted
        ) +
        .tsenat_theme_base(base_size = 11) +
        ggplot2::theme(
            axis.title = ggplot2::element_text(face = "bold")
        )

    p
}


#' Plot Tsallis Entropy q-Curve
#'
#' Visualize Tsallis entropy (S_q) as a function of the diversity parameter q across sample groups.
#' Supports three modes: aggregate q-curves (default), gene-specific q-curves (when `gene` provided),
#' or bootstrap confidence interval bands.
#'
#' @param se A `SummarizedExperiment` returned by `calculate_diversity()` with diversity assay.
#'   For CI mode (bootstrap=TRUE), must contain pre-computed bootstrap confidence intervals.
#' @param assay_name Character; name of the assay to plot (default: "diversity").
#' @param condition_col Character; column name in colData indicating group/sample type
#'   (default: "sample_type"). Only used in aggregate and CI modes.
#' @param bootstrap Logical; if TRUE, plots bootstrap confidence interval bands for aggregate mode.
#'   Requires SE to contain pre-computed CI assays (ci_lower/ci_upper).
#'   Requires 2+ q values and exactly 2 groups (default: FALSE).
#' @param gene Character vector (optional); if provided, plot q-curves for specified gene(s).
#'   Overrides default aggregate behavior. When provided, uses median +/- SD for each gene.
#' @param lm_res Data frame (optional); gene interaction test results with `gene` column and
#'   p-value column. Accepts either:
#'   - Results from `calculate_lm_interaction()` (has `adj_p_interaction` or `p_interaction` columns)
#'   - Results from `detect_q_gene_interactions()` (has `adj_p_value` or `p_value` columns from Friedman/Wilcoxon tests)
#'   If provided (and `gene` is NULL), plots top `n_top` genes ranked by p-value.
#'   Useful for plotting significant genes from any interaction analysis.
#' @param n_top Integer or NULL; number of top genes to select from `lm_res` when `gene` is NULL
#'   (default: NULL). When NULL, defaults to showing the single most significant gene (n_top=1),
#'   providing a conservative view of the strongest effect. Set to a numeric value to show that many top genes.
#' @param output_file \code{character} or \code{NULL}. Optional file path to save the plot.
#'   Default: NULL (no file output).
#'
#' @return
#' **Aggregate mode (gene=NULL, lm_res=NULL)**:
#' - With bootstrap=FALSE: A ggplot object showing median entropy with IQR ribbons.
#' - With bootstrap=TRUE: A ggplot object with bootstrap confidence interval bands.
#'
#' **Gene-specific mode (gene or lm_res provided)**:
#' - Single gene: A ggplot object showing median entropy +/- SD for that gene.
#' - Multiple genes: A grid plot object arranged in 2 rows x 2 columns with a shared legend at the bottom.
#'   The legend appears once beneath the grid, avoiding repetition across subplots.
#'
#' @details
#' **Aggregate mode (default, gene=NULL, lm_res=NULL)**:
#' - Plots median Tsallis entropy +/- IQR across all genes for each group
#' - Works with any SummarizedExperiment from calculate_diversity()
#' - Supports single or multiple q values and any number of groups
#' - No CI data required for basic plots; bootstrap CIs optional
#'
#' **Gene-specific mode (gene or lm_res provided)**:
#' - Plots q-curve separately for each selected gene
#' - Shows median entropy +/- SD (variance) for each gene across q-values and groups
#' - When `lm_res` provided: automatically ranks genes and selects top `n_top` by p-value
#' - Single gene: returns a ggplot object; multiple genes: returns a grid plot (2 rows x 2 columns) with shared legend
#' - For multiple genes: legend appears once at the bottom of the grid to avoid repetition and save space
#' - Useful for highlighting specific genes of interest or significant discoveries
#' - Bootstrap mode not supported in gene-specific mode
#'
#' **Bootstrap CI mode (bootstrap=TRUE in aggregate mode)**:
#' - Displays bootstrap confidence interval bands for each group across q-values
#' - Requires exactly 2 groups for comparison
#' - Requires 2+ q values for q-curve visualization
#' - Requires pre-computed bootstrap CIs from `calculate_diversity(..., bootstrap=TRUE)`
#' - Produces ci_lower, ci_upper assays that properly propagate through entropy transformation
#'
#' @importFrom ggplot2 ggplot aes geom_line geom_ribbon geom_point theme_minimal
#'   scale_color_manual scale_fill_manual labs theme element_text annotate
#' @importFrom dplyr filter group_by summarise pull
#' @importFrom SummarizedExperiment assayNames assay colData rowData
#' @importFrom tidyr pivot_longer
#'
#' @examples
#' # Plot 7: Tsallis entropy q-curve (combined across all sample diversity)
#' analysis <- create_test_analysis(n_genes = 8, n_samples_per_group = 25,
#'   q_values = seq(0.1, 3, by = 0.1), seed = 123)
#' analysis <- calculate_diversity_s4(analysis, q = seq(0.1, 3, by = 0.1), verbose = FALSE)
#' p <- plot_tsallis_q_curve_s4(analysis)
#' if (!is.null(p)) print(p)
#'
#' @export
plot_tsallis_q_curve_s4 <- function(
  se,
  assay_name = "diversity",
  condition_col = "sample_type",
  bootstrap = FALSE,
  gene = NULL,
  lm_res = NULL,
  n_top = NULL,
  output_file = NULL
) {
  require_pkgs(c("ggplot2", "dplyr", "tidyr", "SummarizedExperiment", "cowplot"))
  
  # Handle TSENATAnalysis objects - extract diversity_results
  if (methods::is(se, "TSENATAnalysis")) {
    # Convert diversity_results list to combined SummarizedExperiment
    div_list <- se@diversity_results
    
    # Combine all q-values into one SE with single assay containing all q-value columns
    assay_list <- list()
    first_se <- NULL
    combined_assays_dict <- list()
    
    for (q_name in names(div_list)) {
      obj <- div_list[[q_name]]
      if (methods::is(obj, "SummarizedExperiment")) {
        mat <- SummarizedExperiment::assay(obj, 1)
        if (is.null(first_se)) {
          first_se <- obj  # Save first SE for rowData template
        }
      } else {
        mat <- as.matrix(obj)
      }
      
      # Store with q-value appended to column names
      q_val <- as.numeric(sub("^q_", "", q_name))
      combined_assays_dict[[q_name]] <- list(
        matrix = mat,
        q_val = q_val
      )
    }
    
    # Ensure first_se is not NULL
    if (is.null(first_se)) {
      stop("No valid SummarizedExperiment found in analysis@diversity_results")
    }
    
    # Get target dimensions
    target_genes <- rownames(first_se)
    target_n_cols <- ncol(first_se)
    target_n_qs <- length(combined_assays_dict)
    
    # Create combined assay and colData
    total_cols <- target_n_cols * target_n_qs
    combined_assay <- matrix(0, nrow = length(target_genes), ncol = total_cols)
    rownames(combined_assay) <- target_genes
    
    # Build combined colData
    combined_coldata_list <- list()
    col_idx <- 1
    
    for (q_name in names(combined_assays_dict)) {
      mat <- combined_assays_dict[[q_name]]$matrix
      q_val <- combined_assays_dict[[q_name]]$q_val
      
      # Safety check: ensure matrix has the expected number of columns
      if (ncol(mat) != target_n_cols) {
        message("[plot_tsallis_q_curve] q=", q_val, ": Expected ", target_n_cols, 
            " columns but got ", ncol(mat), ". Adjusting...")
        # If mat has more columns, take only first target_n_cols
        if (ncol(mat) > target_n_cols) {
          mat <- mat[, seq_len(target_n_cols), drop=FALSE]
        } else {
          # If mat has fewer columns, pad with zeros (shouldn't happen)
          mat <- cbind(mat, matrix(0, nrow=nrow(mat), ncol=target_n_cols-ncol(mat)))
        }
      }
      
      # Reorder to match first_se if needed
      mat <- mat[target_genes, , drop = FALSE]
      
      # Get original colnames
      orig_colnames <- colnames(mat)
      if (is.null(orig_colnames)) {
        orig_colnames <- paste0("sample_", seq_len(ncol(mat)))
      }
      
      # Strip any existing _q= suffix before re-adding it (avoid double suffixes)
      # The columns from diversity results may already have _q=X format
      clean_colnames <- sub("_q=.*$", "", orig_colnames)
      if (is.na(clean_colnames[1]) || identical(clean_colnames, orig_colnames)) {
        # If no _q pattern found, use originals as-is
        clean_colnames <- orig_colnames
      }
      
      # Add q-value suffix for uniqueness (use _q= format with 3 decimal precision)
      unique_colnames <- paste0(clean_colnames, "_q=", formatC(q_val, format="f", digits=3))
      
      # Fill in the combined assay
      if (col_idx + ncol(mat) - 1 > total_cols) {
        stop("Dimension mismatch in plot_tsallis_q_curve: ",
                    "Trying to assign to columns ", col_idx, " to ", col_idx + ncol(mat) - 1,
                    ", but combined_assay only has ", total_cols, " columns.\n",
                    "Matrix dimensions: ", nrow(mat), " x ", ncol(mat), "\n",
                    "target_genes: ", length(target_genes), ", target_n_cols: ", target_n_cols,
                    ", target_n_qs: ", target_n_qs)
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
    
    # Set colnames on combined assay to match colData
    combined_coldata_df <- do.call(rbind, combined_coldata_list)
    colnames(combined_assay) <- rownames(combined_coldata_df)
    
    # Ensure rowData is present (fallback to rownames as gene names if needed)
    rd_combined <- tryCatch({
      rd_temp <- SummarizedExperiment::rowData(first_se)
      if (!is.null(rd_temp) && nrow(rd_temp) > 0) {
        rd_temp
      } else {
        NULL
      }
    }, error = function(e) NULL)
    
    # If no rowData, create one with gene identifiers
    if (is.null(rd_combined) || nrow(rd_combined) == 0) {
      rd_combined <- data.frame(
        gene_id = rownames(combined_assay),
        row.names = rownames(combined_assay),
        stringsAsFactors = FALSE
      )
    }
    
    # Create combined SE
    se <- SummarizedExperiment::SummarizedExperiment(
      assays = list(diversity = combined_assay),
      colData = combined_coldata_df,
      rowData = rd_combined
    )
    
    # Use the single assay name
    assay_name <- "diversity"
  }
  
  # Validate input
  if (!methods::is(se, "SummarizedExperiment")) {
    stop("plot_tsallis_q_curve requires a SummarizedExperiment or TSENATAnalysis object")
  }
  
  if (!(assay_name %in% SummarizedExperiment::assayNames(se))) {
    stop("Assay '", assay_name, "' not found in SummarizedExperiment")
  }
  
  # =========================================================================
  # GENE-SPECIFIC MODE (when gene or lm_res is provided)
  # =========================================================================
  if (!is.null(gene) || !is.null(lm_res)) {
    long <- prepare_tsallis_long(se, assay_name = assay_name, condition_col = condition_col)
    
    # Diagnostic: check what columns were created
    if (!("Gene" %in% colnames(long))) {
      # Try to create Gene column if missing - fallback for robustness
      if ("gene" %in% colnames(long)) {
        long$Gene <- long$gene
        long$gene <- NULL
      } else {
        # Last resort: reconstruct from SE rownames
        se_rownames <- rownames(se)
        if (!is.null(se_rownames) && length(se_rownames) > 0) {
          # Each gene should appear the same number of times in long format
          n_per_gene <- nrow(long) / length(se_rownames)
          if (is.integer(n_per_gene) && n_per_gene > 0) {
            long$Gene <- rep(se_rownames, each = n_per_gene)
          } else {
            stop("prepare_tsallis_long did not return Gene column and reconstruction failed.\n",
                 "  long nrow=", nrow(long), ", se nrow=", length(se_rownames), "\n",
                 "  Available columns: ", paste(colnames(long), collapse = ", "))
          }
        } else {
          stop("prepare_tsallis_long did not return Gene column and SE has no rownames.\n",
               "  Available columns: ", paste(colnames(long), collapse = ", "))
        }
      }
    }
    
    # Resolve genes to plot
    if (is.null(gene)) {
      if (is.null(lm_res)) stop("Either 'gene' or 'lm_res' must be provided")
      if (!is.data.frame(lm_res) || !("gene" %in% colnames(lm_res))) stop("'lm_res' must be a data.frame with a 'gene' column")
      
      # Determine p-value column: handle both calculate_lm_interaction and detect_q_gene_interactions formats
      pcol <- NULL
      if ("adj_p_interaction" %in% colnames(lm_res)) {
        # calculate_lm_interaction format
        pcol <- "adj_p_interaction"
      } else if ("p_interaction" %in% colnames(lm_res)) {
        pcol <- "p_interaction"
      } else if ("adj_p_value" %in% colnames(lm_res)) {
        # detect_q_gene_interactions format (from Friedman/Wilcoxon)
        pcol <- "adj_p_value"
      } else if ("p_value" %in% colnames(lm_res)) {
        pcol <- "p_value"
      }
      
      if (is.null(pcol)) stop("'lm_res' must contain one of: 'adj_p_interaction', 'p_interaction', 'adj_p_value', or 'p_value' columns")
      
      genes_ordered <- unique(as.character(lm_res$gene[order(lm_res[[pcol]])]))
      # If n_top is NULL, default to top 1 gene
      n_genes_to_plot <- if (is.null(n_top)) 1 else n_top
      genes <- head(genes_ordered, n_genes_to_plot)
    } else {
      genes <- as.character(unlist(gene))
    }
    
    if (length(genes) == 0) stop("No genes selected for plotting")
    
    # Helper to build single plot for a gene (uses median and variance by default)
    make_plot_for_gene <- function(sel) {
      long_g <- long[as.character(long$Gene) == sel, , drop = FALSE]
      if (nrow(long_g) == 0) stop("Gene not found in assay: ", sel)
      long_g$qnum <- as.numeric(as.character(long_g$q))
      
      # Compute median +/- SD (variance)
      stats_df <- dplyr::summarise(dplyr::group_by(long_g, group, qnum),
        central = median(tsallis, na.rm = TRUE),
        spread = sqrt(stats::var(tsallis, na.rm = TRUE)),
        .groups = "drop"
      )
      
      # Build plot
      p <- ggplot2::ggplot() +
        .tsenat_theme_spectrum(base_size = 11)
      
      # Median +/- SD ribbon
      p <- p +
        ggplot2::geom_ribbon(data = stats_df, ggplot2::aes(x = qnum, ymin = central - spread, ymax = central + spread, fill = group), alpha = 0.2, inherit.aes = FALSE) +
        ggplot2::geom_line(data = stats_df, ggplot2::aes(x = qnum, y = central, color = group), linewidth = 1.3) +
        ggplot2::labs(title = sel, x = "q value", y = "Tsallis entropy", color = "Group", fill = "Group") +
        ggplot2::scale_color_manual(values = .tsenat_palette_blue_red(), name = "Group") + 
        ggplot2::scale_fill_manual(values = .tsenat_palette_blue_red(), name = "Group") +
        ggplot2::theme(plot.title = ggplot2::element_text(
          hjust = 0.5, size = 16, face = "bold",
          margin = ggplot2::margin(b = 10)
        ))
      p
    }
    
    # Return single ggplot for single gene, or a gridded arrangement for multiple genes
    if (length(genes) == 1) {
      return(make_plot_for_gene(genes))
    }
    
    plots <- lapply(genes, make_plot_for_gene)
    names(plots) <- genes
    
    # For multiple genes: create grid with shared legend at bottom using cowplot
    require_pkgs(c("cowplot", "gridExtra"))
    
    # Use cowplot::plot_grid for cleaner handling of shared legends
    # Extract legend from first plot with horizontal layout
    legend_obj <- cowplot::get_legend(
      plots[[1]] + ggplot2::theme(
        legend.position = "bottom",
        legend.direction = "horizontal",
        legend.justification = "center"
      )
    )
    
    # Remove legends from all plots
    plots_no_legend <- lapply(plots, function(p) {
      p + ggplot2::theme(legend.position = "none")
    })
    
    # Arrange plots in 2x2 grid without legend
    grid_with_plots <- do.call(cowplot::plot_grid, c(
      plots_no_legend,
      list(nrow = 2, ncol = 2, align = "hv", axis = "lrtb")
    ))
    
    # Create title and subtitle
    title_plot <- cowplot::ggdraw() + 
      cowplot::draw_label("Tsallis Entropy q-Curve Profile", 
                         fontface = "bold", size = 19, x = 0.5, y = 0.7) +
      cowplot::draw_label("Top genes ranked by statistical significance (Median +/- SD)", 
                         fontface = "italic", size = 15, x = 0.5, y = 0.35, color = "gray40")
    
    # Add legend at bottom
    grid_with_legend <- cowplot::plot_grid(
      title_plot,
      grid_with_plots,
      legend_obj,
      nrow = 3,
      rel_heights = c(0.12, 1, 0.08)
    )
    
    return(grid_with_legend)
  }
  
  # =========================================================================
  # AGGREGATE MODE (when gene and lm_res are NULL)
  # =========================================================================
  
  # =========================================================================
  # CONFIDENCE INTERVAL MODE (Bootstrap only)
  # =========================================================================
  has_bootstrap_ci <- "ci_lower" %in% SummarizedExperiment::assayNames(se) &&
                       "ci_upper" %in% SummarizedExperiment::assayNames(se)
  
  use_ci_mode <- FALSE
  
  if (bootstrap && has_bootstrap_ci) {
    use_ci_mode <- TRUE
  } else if (bootstrap && !has_bootstrap_ci) {
    # User requested bootstrap but not available
    warning("bootstrap=TRUE but bootstrap CIs not found in SE. ",
            "Available assays: ", paste(SummarizedExperiment::assayNames(se), collapse = ", "),
            "\n  Run calculate_diversity(..., bootstrap=TRUE) to generate CI data.",
            "\n  Falling back to basic (non-bootstrap) plot")
  }
  
  bootstrap <- use_ci_mode  # Update bootstrap flag for downstream code
  
  # =========================================================================
  # BASIC MODE (no bootstrap or bootstrap CIs not available)
  # =========================================================================
  if (!bootstrap) {
    long <- prepare_tsallis_long(se, assay_name = assay_name, condition_col = condition_col)
    y_label <- expression("Tsallis entropy (" * S[q] * ")")
    if (nrow(long) == 0) stop("No tsallis values found in SummarizedExperiment")
    
    # Ensure q is numeric
    long$q <- as.numeric(as.character(long$q))
    
    # Compute median and IQR at each q-value for each group
    stats_df <- dplyr::summarise(
      dplyr::group_by(long, group, q),
      median = median(tsallis, na.rm = TRUE),
      IQR = stats::IQR(tsallis, na.rm = TRUE),
      .groups = "drop"
    )
    
    # Create plot with IQR ribbons
    p <- ggplot2::ggplot(
      stats_df,
      ggplot2::aes(x = q, y = median, color = group, fill = group)
    ) +
      ggplot2::geom_line(linewidth = 1.3) +
      ggplot2::geom_ribbon(
        ggplot2::aes(ymin = median - IQR / 2, ymax = median + IQR / 2),
        alpha = 0.2, color = NA
      ) +
      ggplot2::scale_color_manual(values = .tsenat_palette_blue_red(), name = "Group") +
      ggplot2::scale_fill_manual(values = .tsenat_palette_blue_red(), name = "Group") +
      .tsenat_theme_base(base_size = 11) +
      ggplot2::labs(
        title = "Group Comparison: Tsallis Entropy Across Diversity Scales (q-spectrum)",
        subtitle = "Median +/- IQR across samples",
        x = "q value (diversity scale parameter)",
        y = y_label,
        color = "Group",
        fill = "Group"
      ) +
      ggplot2::theme(plot.subtitle = ggplot2::element_text(hjust = 0.5, face = "italic"))
    
    # If only one group, hide legend
    if (length(unique(long$group)) == 1) {
      p <- p + ggplot2::theme(legend.position = "none")
    }
    
    # Save to file if output_file is provided
    if (!is.null(output_file)) {
      ggplot2::ggsave(output_file, plot = p, width = 12, height = 7.2, dpi = 100, create.dir = TRUE)
    }
    
    return(p)
  }
  
  # =========================================================================
  # BOOTSTRAP MODE (with CI data)
  # =========================================================================
  long <- prepare_tsallis_long(se, assay_name = assay_name, condition_col = condition_col)
  if (nrow(long) == 0) {
    stop("No tsallis values found in SummarizedExperiment")
  }
  
  long$q <- as.numeric(as.character(long$q))
  unique_q <- sort(unique(long$q))
  n_q <- length(unique_q)
  
  if (n_q < 2) {
    stop("Need at least 2 q values for q-curve analysis")
  }
  
  # Extract group information
  if (!(condition_col %in% colnames(SummarizedExperiment::colData(se)))) {
    stop("'", condition_col, "' not found in colData")
  }
  
  groups <- unique(sort(long$group))
  if (length(groups) != 2) {
    stop("Expected exactly 2 groups for bootstrap comparison, found ", length(groups))
  }
  
  # Extract bootstrap CIs
  if (!("ci_lower" %in% SummarizedExperiment::assayNames(se) &&
        "ci_upper" %in% SummarizedExperiment::assayNames(se))) {
    stop("No bootstrap CI data found (ci_lower/ci_upper assays required)")
  }
  
  ci_lower_mat <- SummarizedExperiment::assay(se, "ci_lower")
  ci_upper_mat <- SummarizedExperiment::assay(se, "ci_upper")
  
  # Prepare plot data with CI bands
  # Get sample names and their indices in the CI matrices
  sample_names <- colnames(ci_lower_mat)
  if (is.null(sample_names)) {
    sample_names <- paste0("Sample", seq_len(ncol(ci_lower_mat)))
  }
  
  plot_df <- data.frame(
    q = numeric(),
    median = numeric(),
    ci_lower = numeric(),
    ci_upper = numeric(),
    group = character(),
    stringsAsFactors = FALSE
  )
  
  # For each q-value and group, extract median diversity and aggregate CI bounds
  for (group_val in groups) {
    for (q_val in unique_q) {
      # Filter data for this group and q-value
      group_q_data <- long %>%
        dplyr::filter(group == group_val, q == q_val)
      
      if (nrow(group_q_data) > 0) {
        median_val <- median(group_q_data$tsallis, na.rm = TRUE)
        
        # Get unique samples in this group for this q
        group_samples <- unique(group_q_data$sample)
        
        # For each sample, find the CI bounds (average across genes)
        all_ci_lower <- c()
        all_ci_upper <- c()
        
        for (samp in group_samples) {
          # Find column index for this sample in CI matrices
          samp_idx <- which(sample_names == samp)
          if (length(samp_idx) > 0) {
            # Average CI bounds across genes for this sample
            ci_lower_val_samp <- mean(ci_lower_mat[, samp_idx], na.rm = TRUE)
            ci_upper_val_samp <- mean(ci_upper_mat[, samp_idx], na.rm = TRUE)
            all_ci_lower <- c(all_ci_lower, ci_lower_val_samp)
            all_ci_upper <- c(all_ci_upper, ci_upper_val_samp)
          }
        }
        
        # Use median CI bounds across samples in the group
        if (length(all_ci_lower) > 0) {
          ci_lower_final <- median(all_ci_lower, na.rm = TRUE)
          ci_upper_final <- median(all_ci_upper, na.rm = TRUE)
        } else {
          # Fallback: use overall CI from all genes and samples
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
  
  # Create ggplot with CI bands
  p <- ggplot2::ggplot(
    plot_df,
    ggplot2::aes(x = q, y = median, color = group, fill = group)
  ) +
    ggplot2::geom_line(linewidth = 1.2) +
    ggplot2::geom_ribbon(
      ggplot2::aes(ymin = ci_lower, ymax = ci_upper),
      alpha = 0.15,
      color = NA
    ) +
    ggplot2::scale_color_manual(values = .tsenat_palette_blue_red(), name = "Group") +
    ggplot2::scale_fill_manual(values = .tsenat_palette_blue_red(), name = "Group") +
    .tsenat_theme_base(base_size = 11) +
    ggplot2::labs(
      title = "Group Comparison: Tsallis Entropy Across Diversity Scales (q-spectrum)",
      subtitle = "Median with 95% confidence intervals",
      x = "q value (diversity scale parameter)",
      y = expression("Tsallis entropy (" * S[q] * ")"),
      color = "Group",
      fill = "Group"
    ) +
    ggplot2::theme(
      plot.subtitle = ggplot2::element_text(hjust = 0.5, face = "italic")
    )
  
  if (length(groups) == 1) {
    p <- p + ggplot2::theme(legend.position = "none")
  }
  
  return(p)
}


#' Violin plot of Tsallis entropy for a single q value
#'
#' Creates a violin plot showing the distribution of Tsallis entropy for a specific q value,
#' with groups (conditions) displayed side by side.
#'
#' @param se A `SummarizedExperiment` returned by `calculate_diversity` containing
#'   entropy values at one or more q values.
#' @param q_value The specific q value to plot (numeric, e.g., 1, 2, 0.5).
#' @param assay_name Name of the assay to use (default: "diversity").
#' @param title Optional plot title. If NULL, auto-generated based on q value.
#'
#' @return A `ggplot2` object showing a violin plot with groups on the x-axis.
#'  
#' @noRd
#' @keywords internal
plot_tsallis_violin_singleq <- function(se, assay_name = "diversity", title = NULL) {
    require_pkgs(c("ggplot2", "tidyr", "dplyr"))
    
    # Try to extract q from SE metadata first (best source for single-q SE)
    q_val <- NA
    if (!is.null(S4Vectors::metadata(se)$q) && length(S4Vectors::metadata(se)$q) > 0) {
        q_vals <- unique(as.numeric(S4Vectors::metadata(se)$q))
        if (length(q_vals) > 0 && !all(is.na(q_vals))) {
            q_val <- q_vals[1]
        }
    }
    
    # Fallback: use prepare_tsallis_long for data transformation
    long <- prepare_tsallis_long(se, assay_name = assay_name)
    
    if (nrow(long) == 0) stop("No data found in the long format dataframe")
    
    # If still no q, extract from data
    if (is.na(q_val)) {
        q_values <- unique(long$q)
        q_values <- q_values[!is.na(q_values)]
        if (length(q_values) > 0) {
            q_val <- q_values[1]
        }
    }
    
    # Set title
    title_use <- title %||% sprintf("Violin plot: Tsallis entropy at q = %g", q_val)
    
    # Create violin plot
    ggplot2::ggplot(
        long,
        ggplot2::aes(x = group, y = tsallis, fill = group)
    ) +
        ggplot2::geom_violin(
            alpha = 0.5, width = 0.7,
            position = ggplot2::position_dodge(width = 0.8)
        ) +
        ggplot2::geom_boxplot(
            width = 0.2,
            position = ggplot2::position_dodge(width = 0.8),
            outlier.shape = NA,
            alpha = 0.8
        ) +
        .tsenat_theme_base(base_size = 11) +
        ggplot2::scale_fill_manual(values = .tsenat_palette_blue_red(), name = "Group", guide = "none") +
        ggplot2::labs(
            title = title_use,
            x = "Group",
            y = "Tsallis entropy",
            fill = "Group"
        ) +
        ggplot2::theme(
            axis.title = ggplot2::element_text(size = .tsenat_font_sizes$axis_title)
        )
}


#' Density plot of Tsallis entropy for a single q value
#'
#' Creates a density plot showing the distribution of Tsallis entropy for a specific q value,
#' with different groups (conditions) represented by different colors.
#'
#' @param se A `SummarizedExperiment` returned by `calculate_diversity` containing
#'   entropy values at one or more q values.
#' @param q_value The specific q value to plot (numeric, e.g., 1, 2, 0.5).
#' @param assay_name Name of the assay to use (default: "diversity").
#' @param title Optional plot title. If NULL, auto-generated based on q value.
#'
#' @return A `ggplot2` object showing a density plot colored by group.
#'
#' @noRd
#' @keywords internal
plot_tsallis_density_singleq <- function(se, assay_name = "diversity", title = NULL) {
    require_pkgs(c("ggplot2", "tidyr", "dplyr"))
    
    # Try to extract q from SE metadata first (best source for single-q SE)
    q_val <- NA
    if (!is.null(S4Vectors::metadata(se)$q) && length(S4Vectors::metadata(se)$q) > 0) {
        q_vals <- unique(as.numeric(S4Vectors::metadata(se)$q))
        if (length(q_vals) > 0 && !all(is.na(q_vals))) {
            q_val <- q_vals[1]
        }
    }
    
    # Fallback: use prepare_tsallis_long for data transformation
    long <- prepare_tsallis_long(se, assay_name = assay_name)
    
    if (nrow(long) == 0) stop("No data found in the long format dataframe")
    
    # If still no q, extract from data
    if (is.na(q_val)) {
        q_values <- unique(long$q)
        q_values <- q_values[!is.na(q_values)]
        if (length(q_values) > 0) {
            q_val <- q_values[1]
        }
    }
    
    # Set title
    title_use <- title %||% sprintf("Density plot: Tsallis entropy at q = %g", q_val)
    
    # Create density plot
    ggplot2::ggplot(
        long,
        ggplot2::aes(x = tsallis, color = group, fill = group)
    ) +
        ggplot2::geom_density(alpha = 0.3, linewidth = 1) +
        .tsenat_theme_base(base_size = 11) +
        ggplot2::scale_color_manual(values = .tsenat_palette_blue_red(), name = "Group") +
        ggplot2::scale_fill_manual(values = .tsenat_palette_blue_red(), name = "Group") +
        ggplot2::labs(
            title = title_use,
            x = "Tsallis entropy",
            y = "Density",
            color = "Group",
            fill = "Group"
        ) +
        ggplot2::theme(
            axis.title = ggplot2::element_text(size = .tsenat_font_sizes$axis_title)
        )
}


#' Combined Violin and Density Plot Grid for Single q Value
#'
#' Creates a side-by-side grid layout with a violin plot on the left and a density plot
#' on the right, both showing Tsallis entropy distribution for the q value in the provided
#' SummarizedExperiment (which should contain a single q value).
#'
#' @param se A `SummarizedExperiment` returned by `calculate_diversity` containing
#'   entropy values at a single q value.
#' @param assay_name Name of the assay to use (default: "diversity").
#' @param title Optional base title. If NULL, auto-generated based on q value.
#' @param output_file Character or NULL. Optional file path to save the plot as an image.
#'   If provided, the plot will be saved with appropriate dimensions.
#'   Default: NULL (no file output, only return object).
#'
#' @return A `ggplot2` object showing a 1x2 grid with violin plot on the left and
#'   density plot on the right.
#'
#' @export
#' @examples
#' # Plot 8: Violin and density plots of Tsallis entropy distribution
#' analysis <- create_test_analysis(n_genes = 8, n_samples_per_group = 25,
#'   q_values = seq(0.1, 3, by = 0.1), seed = 123)
#' analysis <- calculate_diversity_s4(analysis, q = seq(0.1, 3, by = 0.1), verbose = FALSE)
#' p <- plot_tsallis_violin_density_grid_s4(analysis)
#' if (!is.null(p)) print(p)
#'
plot_tsallis_violin_density_grid_s4 <- function(se, assay_name = "diversity", title = NULL, output_file = NULL) {
    # Require cowplot for grid arrangement
    if (!requireNamespace("cowplot", quietly = TRUE)) {
        stop("cowplot package required for plot_tsallis_violin_density_grid_s4()")
    }
    
    # Handle TSENATAnalysis objects - extract first diversity result
    if (methods::is(se, "TSENATAnalysis")) {
        if (length(se@diversity_results) == 0) {
            stop("No diversity results found in TSENATAnalysis object. Run calculate_diversity_s4() first.")
        }
        # Extract first diversity result
        se <- se@diversity_results[[1]]
    }
    
    # Try to extract q from SE metadata first (best source for single-q SE)
    q_val <- NA
    if (!is.null(S4Vectors::metadata(se)$q) && length(S4Vectors::metadata(se)$q) > 0) {
        q_vals <- unique(as.numeric(S4Vectors::metadata(se)$q))
        if (length(q_vals) > 0 && !all(is.na(q_vals))) {
            q_val <- q_vals[1]
        }
    }
    
    # Fallback: use prepare_tsallis_long for data transformation
    long <- prepare_tsallis_long(se, assay_name = assay_name)
    
    # If still no q, extract from data
    if (is.na(q_val)) {
        q_values <- unique(long$q)
        q_values <- q_values[!is.na(q_values)]
        if (length(q_values) > 0) {
            q_val <- q_values[1]
        }
    }
    
    # Generate base title
    base_title <- title %||% sprintf("Tsallis entropy at q = %g", q_val)
    
    # Create individual plots
    p_violin <- plot_tsallis_violin_singleq(
        se = se,
        assay_name = assay_name,
        title = "Violin"
    )
    
    p_density <- plot_tsallis_density_singleq(
        se = se,
        assay_name = assay_name,
        title = "Density"
    )
    
    # Arrange plots side by side: violin on left, density on right
    grid <- cowplot::plot_grid(
        p_violin,
        p_density,
        nrow = 1,
        ncol = 2,
        align = "h",
        axis = "b"
    )
    
    # Add overall title and subtitle above the grid
    grid_with_title <- cowplot::plot_grid(
        cowplot::ggdraw() + 
            cowplot::draw_label("Tsallis Entropy Distribution by Group", fontface = "bold", size = 19, x = 0.5, y = 0.75) +
            cowplot::draw_label("Violin and density plots across samples", fontface = "italic", size = 15, x = 0.5, y = 0.25, color = "gray40"),
        grid,
        nrow = 2,
        rel_heights = c(0.08, 1)
    )
    
    # Save to file if output_file is provided
    if (!is.null(output_file)) {
      ggplot2::ggsave(output_file, plot = grid_with_title, width = 12, height = 7.2, dpi = 100, create.dir = TRUE)
    }
    
    return(grid_with_title)
}


#' Volcano plot for differential results
#'
#' Create a volcano plot showing fold-change (x-axis) versus adjusted
#' p-value significance (y-axis). The function auto-detects a suitable x-axis
#' column if one is not provided and expects an adjusted p-value column for
#' significance coloring.
#'
#' @param diff_df Data.frame with differential expression results.
#'   Should contain p-values and optionally fold-change columns.
#' @param x_col Optional column name for the x-axis. If `NULL`, the function
#'   will try to auto-detect a suitable numeric column (excluding p-values).
#' @param padj_col Adjusted p-value column name (default: "padj").
#' @param label_thresh Fold-change threshold used to annotate points (default: 0.1).
#' @param sig_alpha Adjusted p-value cutoff for significance (default: 0.05).
#' @param top_n Number of top significant genes to label (default: 5).
#' @param title Optional plot title; if `NULL` a default title is used.
#'
#' @return A `ggplot2` object.
#' @noRd
#' @keywords internal
plot_volcano <- function(
  diff_df,
  x_col = NULL,
  padj_col = "padj",
  label_thresh = 0.1,
  sig_alpha = 0.05,
  top_n = 5,
  title = NULL
) {
    if (!requireNamespace("ggplot2", quietly = TRUE)) {
        stop("ggplot2 required")
    }

    prep_volcano <- .tsenat_prepare_volcano_df(diff_df = diff_df, x_col = x_col, padj_col = padj_col, label_thresh = label_thresh, sig_alpha = sig_alpha, title = title)
    df <- prep_volcano$df
    x_col <- prep_volcano$x_col
    padj_col <- prep_volcano$padj_col
    x_label_formatted <- prep_volcano$x_label_formatted
    padj_label_formatted <- prep_volcano$padj_label_formatted
    title_use <- prep_volcano$title_use

    p <- ggplot2::ggplot(
        df,
        ggplot2::aes(x = xval, y = -log10(padj), color = significant)
    ) +
        ggplot2::geom_point(alpha = 0.75, size = 3.4) +
        ggplot2::scale_color_manual(
            values = .tsenat_significance_colors(),
            guide = "none"
        ) +
        ggplot2::geom_hline(
            yintercept = -log10(sig_alpha),
            linetype = "dashed",
            color = "gray50"
        ) +
        ggplot2::geom_vline(
            xintercept = c(-label_thresh, label_thresh),
            linetype = "dashed",
            color = "gray50"
        ) +
        ggplot2::labs(
            title = title_use,
            x = x_label_formatted,
            y = paste0("-Log10(", padj_label_formatted, ")")
        ) +
        .tsenat_theme_base(base_size = 11) +
        ggplot2::theme(
            axis.title = ggplot2::element_text(face = "bold")
        )

    p
}


#' Combine Volcano and MA-Tsallis Plots in a Grid Layout
#'
#' Creates a side-by-side grid layout with a volcano plot on the left and an MA-Tsallis plot on the right.
#' Both plots are generated from differential analysis results data.
#'
#' @param diff_df Data.frame from differential analysis containing required columns for both volcano and MA plots.
#' @param x_col Column name for x-axis in volcano plot (e.g., "mean_difference"). Auto-detected if NULL.
#' @param padj_col Column name for adjusted p-values (default: "padj").
#' @param label_thresh Threshold for volcano plot labels (default: 0.1).
#' @param sig_alpha Numeric significance threshold for adjusted p-values (default: 0.05).
#' @param top_n Number of top genes to annotate in volcano plot (default: 5).
#' @param title_volcano Title for volcano plot. If NULL, auto-generated.
#' @param title_ma Title for MA plot (default: "Tsallis-based MA plot").
#' @param ... Additional arguments passed to plotting functions.
#'
#' @return A `ggplot2` object showing a 1x2 grid with volcano plot on the left and MA plot on the right.
#'
#' @examples
#' # Simulate differential analysis results
#' x <- data.frame(
#'   genes = paste0("g", seq_len(20)),
#'   mean_difference = rnorm(20, sd = 1),
#'   padj = runif(20, 1e-5, 0.1),
#'   log2_fold_change = rnorm(20, sd = 0.8)
#' )
#' # Placeholder: actual usage would require valid differential results
#' # plot_volcano_ma_grid(x, sig_alpha = 0.05)
#'
#' @keywords internal
#' @noRd
plot_volcano_ma_grid <- function(
  diff_df,
  x_col = NULL,
  padj_col = "padj",
  label_thresh = 0.1,
  sig_alpha = 0.05,
  top_n = 5,
  title_volcano = NULL,
  title_ma = "Tsallis-based MA plot",
  ...
) {
    # Require cowplot for grid arrangement
    if (!requireNamespace("cowplot", quietly = TRUE)) {
        stop("cowplot package required for plot_volcano_ma_grid()")
    }

    # Create volcano plot
    p_volcano <- plot_volcano(
        diff_df = diff_df,
        x_col = x_col,
        padj_col = padj_col,
        label_thresh = label_thresh,
        sig_alpha = sig_alpha,
        top_n = top_n,
        title = title_volcano
    )

    # Create MA plot
    p_ma <- plot_ma_tsallis(
        x = diff_df,
        sig_alpha = sig_alpha,
        title = title_ma,
        ...
    )

    # Arrange plots side by side: volcano on left, MA on right
    grid <- cowplot::plot_grid(
        p_volcano,
        p_ma,
        nrow = 1,
        ncol = 2,
        align = "h",
        axis = "b"
    )

    return(grid)
}


#' Internal helper to compute fill limits across multiple genes (not exported)
#' @noRd
.compute_transcript_fill_limits <- function(genes, mapping, counts, samples, top_n, agg_fun, pseudocount) {
    mins <- maxs <- c()
    for (g in genes) {
        txs <- mapping$Transcript[mapping$Gen == g]
        txs <- intersect(txs, rownames(counts))
        if (length(txs) == 0) next
        if (!is.null(top_n)) txs <- head(txs, top_n)
        mat <- counts[txs, , drop = FALSE]
        df_all <- as.data.frame(mat)
        df_all$tx <- rownames(mat)
        df_long <- tidyr::pivot_longer(df_all, -tx, names_to = "sample", values_to = "expr")
        df_long$group <- rep(samples, times = length(txs))
        df_summary <- stats::aggregate(expr ~ tx + group, data = df_long, FUN = agg_fun)
        df_summary$log2expr <- log2(df_summary$expr + pseudocount)
        mins <- c(mins, min(df_summary$log2expr, na.rm = TRUE))
        maxs <- c(maxs, max(df_summary$log2expr, na.rm = TRUE))
    }
    if (length(mins) == 0) stop("No transcripts found for provided genes")
    c(min(mins, na.rm = TRUE), max(maxs, na.rm = TRUE))
}

#' Internal helper to draw grid layout with title, plots, and legend using base grid
#' @noRd
.draw_transcript_grid <- function(grobs, title, legend_grob, ncol, heights, to_file = NULL) {
    # If no output file is provided and no graphics device is open, render to a
    # temporary pdf device so that plotting in non-interactive sessions does not
    # create `Rplots.pdf` in the working directory.
    temp_dev <- FALSE
    # Only open a temporary PDF device when:
    #  - caller did not supply an output file (`to_file` is NULL),
    #  - the session is non-interactive, and
    #  - no graphics device is currently open (dev.cur() == 1)
    if (is.null(to_file) && !interactive() && grDevices::dev.cur() == 1L) {
        tmp <- tempfile("TSENAT_plot_", fileext = ".pdf")
        grDevices::pdf(tmp)
        temp_dev <- TRUE
        # Ensure device is closed and temporary file removed on exit
        on.exit(
            {
                try(grDevices::dev.off(), silent = TRUE)
                if (file.exists(tmp)) unlink(tmp)
            },
            add = TRUE
        )
    }

    # Calculate number of rows needed for plots
    nrow_plots <- ceiling(length(grobs) / ncol)
    nrow_total <- 2 + nrow_plots  # title + plot rows + legend

    grid::grid.newpage()
    grid::pushViewport(
        grid::viewport(layout = grid::grid.layout(nrow_total, ncol, heights = heights))
    )
    # Title row
    vp_title <- grid::viewport(layout.pos.row = 1, layout.pos.col = seq_len(ncol))
    grid::pushViewport(vp_title)
    grid::grid.text("Transcript level expression", x = 0.5, y = 0.6, gp = grid::gpar(fontsize = 14, fontface = "bold"))
    grid::grid.text(paste0("Top genes with metric ", title), x = 0.5, y = 0.2, gp = grid::gpar(fontsize = 11, fontface = "italic", col = "gray40"))
    grid::upViewport()
    # Plot rows
    for (i in seq_along(grobs)) {
        plot_row_idx <- ((i - 1) %/% ncol) + 2
        plot_col_idx <- ((i - 1) %% ncol) + 1
        vp <- grid::viewport(layout.pos.row = plot_row_idx, layout.pos.col = plot_col_idx)
        grid::pushViewport(vp)
        grid::grid.draw(grobs[[i]])
        grid::upViewport()
    }
    # Legend row
    if (!is.null(legend_grob)) {
        vp_leg <- grid::viewport(layout.pos.row = nrow_total, layout.pos.col = seq_len(ncol))
        grid::pushViewport(vp_leg)
        grid::grid.draw(legend_grob)
        grid::upViewport()
    }
    grid::upViewport()
    # If caller supplied a file (caller likely opened a device), close it here.
    if (!is.null(to_file)) grDevices::dev.off()
    invisible(NULL)
}

## Internal helpers for `plot_top_transcripts` refactor
## Create per-gene plot and combine multiple gene plots into final output
.ptt_make_plot_for_gene <- function(gene_single, mapping, counts, samples, top_n, agg_fun, pseudocount, agg_label_unique, fill_limits = NULL, font_scale = 1.0) {
    require_pkgs(c("ggplot2", "tidyr"))
    built <- .ptt_build_tx_long(gene_single, mapping, counts, samples, NULL)
    df_summary <- .ptt_aggregate_df_long(built$df_long, agg_fun, pseudocount)
    .ptt_build_plot_from_summary(df_summary, agg_label_unique, fill_limits, font_scale = font_scale)
}

.ptt_combine_plots <- function(plots, output_file = NULL, agg_label_unique = NULL) {
    require_pkgs(c("ggplot2"))
    # Allow callers to pass a single character second argument as the
    # `agg_label_unique` for convenience (legacy test call patterns).
    if (is.null(agg_label_unique) && !is.null(output_file) && is.character(output_file) && length(output_file) == 1) {
        agg_label_unique <- output_file
        output_file <- NULL
    }
    if (requireNamespace("patchwork", quietly = TRUE)) {
        .ptt_combine_patchwork(plots, agg_label_unique)
    } else if (requireNamespace("cowplot", quietly = TRUE)) {
        .ptt_combine_cowplot(plots, output_file = output_file, agg_label_unique = agg_label_unique)
    } else {
        .ptt_combine_grid(plots, output_file = output_file, agg_label_unique = agg_label_unique)
    }
}

## Prepare and validate inputs for `plot_top_transcripts`
.ptt_prepare_inputs <- function(counts, readcounts = NULL, samples = NULL, coldata = NULL, condition_col = "sample_type", tx2gene = NULL, res = NULL, top_n = NULL, pseudocount = 0, output_file = NULL, metric = c("median", "mean", "variance", "iqr")) {
    # handle selecting genes from `res` is left to caller; this function focuses
    # on normalizing counts, samples and tx2gene mapping and preparing agg functions
    if (inherits(counts, "SummarizedExperiment")) {
        require_pkgs(c("SummarizedExperiment", "S4Vectors"))
        se <- counts
        counts_mat <- get_readcounts_from_se(se, readcounts)
        counts <- as.matrix(counts_mat)
        samples <- infer_samples_from_se(se, samples, condition_col = condition_col)

        if (is.null(tx2gene)) {
            txres <- get_tx2gene_from_se(se, counts)
            if (!is.null(txres) && !is.null(txres$mapping)) {
                mapping <- data.frame(Transcript = rownames(counts), Gen = as.character(txres$mapping), stringsAsFactors = FALSE)
                tx2gene <- mapping
            }
        }
    }

    if (!is.matrix(counts) && !is.data.frame(counts)) stop("`counts` must be a matrix or data.frame with transcripts as rownames")
    counts <- as.matrix(counts)
    if (is.null(rownames(counts))) stop("`counts` must have rownames corresponding to transcript identifiers")

    # derive samples from coldata if needed
    if (is.null(samples)) {
        if (!is.null(coldata)) {
            if (is.character(coldata) && length(coldata) == 1) {
                if (!file.exists(coldata)) stop("coldata file not found: ", coldata)
                cdf <- utils::read.delim(coldata, header = TRUE, stringsAsFactors = FALSE)
            } else if (is.data.frame(coldata)) {
                cdf <- coldata
            } else {
                stop("`coldata` must be a data.frame or path to a tab-delimited file")
            }

            if (!is.null(rownames(cdf)) && all(colnames(counts) %in% rownames(cdf))) {
                samples <- as.character(cdf[colnames(counts), condition_col])
            } else {
                sample_id_cols <- c("sample", "Sample", "sample_id", "id")
                sid <- intersect(sample_id_cols, colnames(cdf))
                if (length(sid) > 0) {
                    sid <- sid[1]
                    if (!all(colnames(counts) %in% as.character(cdf[[sid]]))) stop("coldata sample id column does not match column names of counts")
                    row_ix <- match(colnames(counts), as.character(cdf[[sid]]))
                    samples <- as.character(cdf[[condition_col]][row_ix])
                } else {
                    stop("Could not match `coldata` rows to `counts` columns. Provide `samples` or a row-named `coldata`.")
                }
            }
        } else {
            stop("Either 'samples' or 'coldata' must be provided to determine sample groups")
        }
    }

    # normalize tx2gene mapping
    if (is.null(tx2gene)) stop("`tx2gene` must be provided as a file path or data.frame (or include mapping in metadata of provided SummarizedExperiment)")
    if (is.character(tx2gene) && length(tx2gene) == 1) {
        if (!file.exists(tx2gene)) stop("tx2gene file not found: ", tx2gene)
        mapping <- utils::read.delim(tx2gene, stringsAsFactors = FALSE, header = TRUE)
    } else if (is.data.frame(tx2gene)) {
        mapping <- tx2gene
    } else {
        stop("`tx2gene` must be provided as a file path or data.frame (or include mapping in metadata of provided SummarizedExperiment)")
    }

    if (!all(c("Transcript", "Gen") %in% colnames(mapping))) stop("tx2gene must have columns 'Transcript' and 'Gen'")

    if (!requireNamespace("ggplot2", quietly = TRUE)) stop("ggplot2 required for plotting")

    if (!is.null(samples) && length(samples) != ncol(counts)) stop("Length of `samples` must equal number of columns in `counts`")

    metric_choice <- match.arg(metric)
    agg_fun <- switch(metric_choice,
        median = function(x) stats::median(x, na.rm = TRUE),
        mean = function(x) base::mean(x, na.rm = TRUE),
        variance = function(x) stats::var(x, na.rm = TRUE),
        iqr = function(x) stats::IQR(x, na.rm = TRUE)
    )
    agg_label_metric <- if (metric_choice == "iqr") "IQR" else metric_choice
    agg_label <- agg_label_metric
    .cnt <- as.integer(getOption("TSENAT.plot_top_counter", 0)) + 1L
    options(TSENAT.plot_top_counter = .cnt)
    agg_label_unique <- agg_label

    list(counts = counts, samples = samples, mapping = mapping, metric_choice = metric_choice, agg_fun = agg_fun, agg_label_unique = agg_label_unique, top_n = top_n, pseudocount = pseudocount, output_file = output_file)
}

#' Plot top transcripts for a gene
#' @param se A `SummarizedExperiment` with transcript counts as assay and gene information in rowData.
#'   Must have a "genes" column in rowData specifying which gene each transcript belongs to.
#' @param gene Character vector; gene symbol(s) to inspect. If NULL and `res` is provided, 
#'   top genes are selected by p-value.
#' @param condition_col Character; column name in colData(se) to use for sample grouping 
#'   (default: "sample_type").
#' @param res Optional result data.frame from differential/interaction analysis with gene identifiers and p-values.
#'   Supported sources:
#'   - `calculate_lm_interaction(..., return_model_data = TRUE)` returns a list with $results and $model_data
#'   - `calculate_lm_interaction(..., return_model_data = FALSE)` returns a data.frame with adj_p_interaction column
#'   - `detect_q_gene_interactions()` returns a data.frame with adj_p_value column (for Friedman/Kruskal-Wallis tests)
#'   If provided and `gene` is NULL, top genes are selected by adjusted p-value.
#' @param top_n Integer number of transcripts to show (default = 3). Use NULL to plot all transcripts for the gene.
#' @param output_file Optional file path to save the plot. If `NULL`, the `ggplot` object is returned.
#' @param metric Aggregation metric: "median", "mean", "variance", or "iqr" (default: "median").
#' @param width Output image width in inches. If NULL, automatically calculated based on number of genes (default: ~13 inches per column).
#' @param height Output image height in inches. If NULL, automatically calculated based on number of genes (default: ~10 inches per row + headers).
#' @return If `output_file` is `NULL`, returns a `ggplot` object. Otherwise saves to file and returns NULL invisibly.
#' @examples
#' library(SummarizedExperiment)
#' library(S4Vectors)
#' # Create example SummarizedExperiment
#' counts <- matrix(sample(1:100, 36, replace = TRUE), nrow = 6, ncol = 6)
#' rownames(counts) <- paste0("tx", 1:6)
#' rowData_df <- DataFrame(genes = rep(paste0("G", 1:3), each = 2))
#' colData_df <- DataFrame(sample_type = rep(c("Normal", "Tumor"), 3))
#' se <- SummarizedExperiment(assays = list(counts = counts), 
#'                           rowData = rowData_df, colData = colData_df)
#' # Plot top transcripts
#' plot_top_transcripts(se, gene = "G1", top_n = 2)
#' @keywords internal
#' @noRd
plot_top_transcripts <- function(
  se,
  gene = NULL,
  condition_col = "sample_type",
  res = NULL,
  top_n = 3,
  output_file = NULL,
  metric = c("median", "mean", "variance", "iqr"),
  width = NULL,
  height = NULL
) {
    require_pkgs(c("SummarizedExperiment", "S4Vectors"))
    
    # Validate input
    if (!inherits(se, "SummarizedExperiment")) {
        stop("se must be a SummarizedExperiment object", call. = FALSE)
    }
    
    # Extract components from SE
    counts <- as.matrix(SummarizedExperiment::assay(se))
    rd <- SummarizedExperiment::rowData(se)
    cd <- SummarizedExperiment::colData(se)
    
    # Identify gene column (prefer "genes", then "gene_name", then "gene_id")
    gene_col <- if ("genes" %in% colnames(rd)) {
        "genes"
    } else if ("gene_name" %in% colnames(rd)) {
        "gene_name"
    } else if ("gene_id" %in% colnames(rd)) {
        "gene_id"
    } else {
        stop("rowData(se) must contain a 'genes', 'gene_name', or 'gene_id' column", call. = FALSE)
    }
    
    # Build tx2gene mapping
    tx2gene <- data.frame(
        Transcript = rownames(counts),
        Gen = as.character(rd[[gene_col]]),
        stringsAsFactors = FALSE
    )
    
    # Extract sample groups
    if (!condition_col %in% colnames(cd)) {
        stop("Column '", condition_col, "' not found in colData(se)", call. = FALSE)
    }
    samples <- as.character(cd[[condition_col]])
    
    # Handle gene selection from results if needed
    if (is.null(gene) && !is.null(res)) {
        # Extract data.frame if res is a list (from calculate_lm_interaction with return_model_data = TRUE)
        if (is.list(res) && !is.data.frame(res) && "results" %in% names(res)) {
            res <- res$results
        }
        
        if (!is.data.frame(res)) {
            stop("res must be a data.frame or a list with $results component from calculate_lm_interaction() or similar analysis function", call. = FALSE)
        }
        
        # Normalize p-value column name for compatibility across sources:
        # - calculate_lm_interaction uses: adj_p_interaction
        # - detect_q_gene_interactions uses: adj_p_value
        # Standardize to a common column for downstream use
        if ("adj_p_interaction" %in% colnames(res) && "adj_p_value" %in% colnames(res) == FALSE) {
            # From calculate_lm_interaction (LMM/GAM method)
            colnames(res)[colnames(res) == "adj_p_interaction"] <- "adj_p_value"
        }
        
        # Find gene column in results
        res_gene_col <- if ("gene" %in% colnames(res)) {
            "gene"
        } else if ("genes" %in% colnames(res)) {
            "genes"
        } else if ("gene_id" %in% colnames(res)) {
            "gene_id"
        } else {
            stop("res must contain 'gene', 'genes', or 'gene_id' column", call. = FALSE)
        }
        
        # Find adjusted p-value column (supports multiple naming conventions)
        # Priority: adj_p_value (from both calculate_lm_interaction and detect_q_gene_interactions)
        #         padj (legacy support)
        #         adjusted_p_values (legacy support)
        p_col <- if ("adj_p_value" %in% colnames(res)) {
            "adj_p_value"
        } else if ("padj" %in% colnames(res)) {
            "padj"
        } else if ("adjusted_p_values" %in% colnames(res)) {
            "adjusted_p_values"
        } else {
            stop("res must contain adjusted p-value column. Expected: 'adj_p_value' (from calculate_lm_interaction or detect_q_gene_interactions), 'padj', or 'adjusted_p_values'", call. = FALSE)
        }
        
        # Sort by p-value and select top genes
        res_sorted <- res[order(res[[p_col]], na.last = NA), ]
        gene <- head(as.character(res_sorted[[res_gene_col]]), top_n)
        
        # Filter to genes present in SE
        se_genes <- unique(tx2gene$Gen)
        gene <- gene[gene %in% se_genes]
        
        if (length(gene) == 0) {
            stop("No genes from res found in rowData of SE", call. = FALSE)
        }
    }
    
    if (is.null(gene)) {
        stop("gene must be provided or derivable from res", call. = FALSE)
    }

    ## Prepare inputs and normalization via helper
    prep <- .ptt_prepare_inputs(
        counts = counts, 
        readcounts = NULL, 
        samples = samples, 
        coldata = NULL, 
        condition_col = condition_col, 
        tx2gene = tx2gene, 
        res = NULL,
        top_n = top_n, 
        pseudocount = 1e-6, 
        output_file = output_file, 
        metric = metric
    )

    # Extract prepared data
    counts <- prep$counts
    samples <- prep$samples
    mapping <- prep$mapping
    metric_choice <- prep$metric_choice
    agg_fun <- prep$agg_fun
    agg_label_unique <- prep$agg_label_unique
    top_n <- prep$top_n
    pseudocount <- prep$pseudocount
    output_file <- prep$output_file

    make_plot_for_gene <- function(gene_single, fill_limits = NULL) {
        .ptt_make_plot_for_gene(gene_single, mapping, counts, samples, top_n, agg_fun, pseudocount, agg_label_unique, fill_limits, font_scale = font_scale)
    }

    # Calculate dimensions and font scaling upfront
    n_cols <- 2
    n_rows <- ceiling(length(gene) / n_cols)
    
    # Standardized dimensions: 12 inches width with proportional height
    # Base: 6 inches per column width, 3 inches per row height, plus 2 inch margin
    calc_width <- if (is.null(width)) 12 else width
    calc_height <- if (is.null(height)) 2 + (3 * n_rows) else height
    
    # Calculate font scaling based on actual dimensions
    # Reference: 12x8 inches (96 sq in) uses base formula for 11pt fonts
    # For other sizes: font_scale = sqrt(area / 96)
    plot_area <- calc_width * calc_height
    font_scale <- sqrt(plot_area / 96)

    # Produce plots (single or multiple)
    if (length(gene) > 1) {
        fill_limits <- .compute_transcript_fill_limits(gene, mapping, counts, samples, top_n, agg_fun, pseudocount)

        plots <- lapply(seq_along(gene), function(i) {
            gname <- gene[i]
            pp <- make_plot_for_gene(gname, fill_limits = fill_limits)
            per_gene_title <- if (!is.na(gname) && nzchar(as.character(gname))) as.character(gname) else ""
            # Scale title font proportionally
            scaled_title_size <- 16 * font_scale
            pp <- pp + ggplot2::labs(title = per_gene_title) + ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5, size = scaled_title_size, face = "bold", margin = ggplot2::margin(b = 5)))
            pp
        })

        result_plot <- .ptt_combine_plots(plots, output_file = output_file, agg_label_unique = agg_label_unique)
    } else {
        result_plot <- make_plot_for_gene(gene)
    }

    if (!is.null(output_file)) {
        # Create directory if it doesn't exist
        output_dir <- dirname(output_file)
        if (!dir.exists(output_dir) && nzchar(output_dir) && output_dir != ".") {
            dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
        }
        
        # Use pre-calculated dimensions
        plot_width <- calc_width
        plot_height <- calc_height
        
        ggplot2::ggsave(output_file, result_plot, width = plot_width, height = plot_height, dpi = 100)
        invisible(output_file)
    } else {
        result_plot
    }
}

#' Plot GAM q-curves for top genes identified by FPCA/GAM interaction tests
#'
#' Visualizes smooth q-curve profiles (GAM fits) for selected genes from
#' `calculate_lm_interaction()` results. Useful for understanding which q-ranges (rare
#' vs. dominant isoforms) drive significant PC differences between groups.
#'
#' @param se A `SummarizedExperiment` with q-sequence diversity values
#'   (multiple q values per sample). Typically output from `calculate_diversity()`
#'   with multiple q (e.g., q = seq(0.1, 2, by = 0.1)).
#' @param lm_res A `data.frame` from `calculate_lm_interaction()` with columns
#'   `gene`, `p_interaction`, and `adj_p_interaction`. Can be from method="fpca"
#'   or method="gam".
#' @param condition_col Column name in `colData(se)` specifying group
#'   assignments for samples (e.g., "sample_type").
#' @param genes Optional character vector of specific gene names to plot. If provided,
#'   these genes are plotted directly regardless of significance or n_top. If NULL
#'   (default), the top n_top significant genes are selected.
#' @param n_top Number of top genes (by adjusted p-value) to plot (default: 6).
#'   Only used if genes = NULL.
#' @param sig_alpha Significance threshold for adjusted p-values (default: 0.05).
#'   Only used if genes = NULL; filters lm_res to significant genes before selecting top n.
#' @param assay_name Name of the assay in `se` to extract (default: "diversity").
#' @param model_data Required list from `calculate_lm_interaction(..., return_model_data = TRUE)$model_data`
#'   containing metadata (q_values, sample configuration, etc.). This is the preferred way to use
#'   this function as it ensures all visualizations are based on the exact analysis configuration.
#'
#' @return A single `ggplot` object with all selected genes arranged in a grid layout 
#'   (2 columns per row). Can be saved with `ggplot2::ggsave()`.
#'
#' @details
#' For each selected gene, this function:
#' 1. Extracts per-sample entropy values across all q values
#' 2. Fits GAM models: entropy ~ s(q, k=...) independently for each group
#' 3. Generates smooth predictions for visualization
#' 4. Overlays predicted curves for each group with a distinct color
#'
#' By providing `model_data` from `calculate_lm_interaction()`, the function can directly
#' access the q-values used in the original analysis for more accurate visualization.
#'
#' This complements FPCA by providing interpretable visualization of empirical
#' q-curve shape differences that drive PC-level significance.
#'
#' @examples
#' # Create example data with multiple q values
#' set.seed(123)
#' counts <- matrix(
#'   sample(1:100, 60, replace = TRUE),
#'   nrow = 15, ncol = 4
#' )
#' rownames(counts) <- paste0("tx_", 1:15)
#' colnames(counts) <- paste0("sample_", 1:4)
#' genes <- rep(paste0("gene_", 1:5), each = 3)
#' 
#' # Calculate diversity across multiple q values
#' se <- calculate_diversity(counts, genes = genes, q = seq(0.5, 2, by = 0.5), norm = TRUE)
#' 
#' # Add sample metadata
#' SummarizedExperiment::colData(se) <- S4Vectors::DataFrame(
#'   condition = rep(c('Normal', 'Tumor'), length.out = ncol(se)),
#'   row.names = colnames(se)
#' )
#' 
#' # Run linear model analysis with model_data  
#' lm_result <- calculate_lm_interaction(se, condition_col = "condition", method = "gam",
#'                                       return_model_data = TRUE)
#' 
#' # Plot GAM curves for top genes
#' if (nrow(lm_result$results) > 0) {
#'   grid_plot <- plot_lm_interaction_gam(se, lm_result$results, condition_col = "condition",
#'                                         n_top = 2, model_data = lm_result$model_data)
#' }
#'
#' @keywords internal
#' @noRd
#' @importFrom ggplot2 ggplot aes geom_line geom_point facet_wrap labs theme_minimal scale_color_brewer
#' @importFrom cowplot plot_grid
plot_lm_interaction_gam <- function(se, lm_res, condition_col = "sample_type", genes = NULL, n_top = 6,
    sig_alpha = 0.05, assay_name = "diversity", model_data = NULL, output_file = NULL, width = NULL, height = NULL) {

    require_pkgs(c("ggplot2", "mgcv", "SummarizedExperiment", "dplyr", "tidyr", "cowplot"))

    # Validate inputs
    if (!inherits(se, "SummarizedExperiment")) {
        stop("se must be a SummarizedExperiment", call. = FALSE)
    }

    # Handle flexible input: lm_res can be either:
    # 1. A data.frame with results (traditional usage)
    # 2. A list with $results and $model_data (return_model_data = TRUE format)
    if (is.list(lm_res) && !is.data.frame(lm_res)) {
        # lm_res is a list with components
        if ("results" %in% names(lm_res) && is.data.frame(lm_res$results)) {
            # Extract results and model_data from the list
            extracted_results <- lm_res$results
            
            # If model_data not provided, extract from lm_res
            if (is.null(model_data) && "model_data" %in% names(lm_res)) {
                model_data <- lm_res$model_data
            }
            
            lm_res <- extracted_results
        } else {
            stop("lm_res is a list but does not contain 'results' data.frame component",
                call. = FALSE)
        }
    }
    
    if (!is.data.frame(lm_res) || !("gene" %in% colnames(lm_res))) {
        stop("lm_res must be either:\n  1. A data.frame with 'gene' column from calculate_lm_interaction()\n  2. A list with $results and $model_data from return_model_data = TRUE",
            call. = FALSE)
    }

    if (nrow(lm_res) == 0) {
        stop("lm_res has no rows; calculate_lm_interaction() returned no genes", call. = FALSE)
    }
    
    # =========================================================================
    # VALIDATE DATA COMPATIBILITY (NEW: Issue #5 validation)
    # =========================================================================
    # Check that SE and LM results have compatible gene sets before processing
    data_validation <- validate_plot_data(
        se = se,
        lm_results = lm_res,
        stop_on_error = TRUE,
        verbose = FALSE  # Suppress verbose; we'll only see output if validation fails
    )
    
    # Validate and extract metadata from model_data
    if (is.null(model_data)) {
        stop("model_data is required. Provide it as a parameter or pass full lm_res list with $model_data component",
            call. = FALSE)
    }
    
    if (!is.list(model_data)) {
        stop("model_data must be a list from calculate_lm_interaction(..., return_model_data = TRUE)",
            call. = FALSE)
    }
    
    # Extract required metadata - handle both wrapped (from JSON) and unwrapped formats
    # JSON returns arrays: method = [["gam"]], need [[1]]
    # Direct list returns: method = "gam", no [[1]] needed
    q_values <- model_data$q_values
    if (is.null(q_values)) {
        stop("model_data must contain 'q_values' from the original analysis",
            call. = FALSE)
    }
    
    # Normalize q_values in case it's wrapped in list
    if (is.list(q_values) && length(q_values) == 1) {
        q_values <- unlist(q_values)
    } else {
        q_values <- unlist(q_values)
    }

    # Match and filter genes between SE and lm_res
    # After calculate_diversity, rownames(SE) are gene names
    # lm_res$gene column also contains gene names
    gene_names_in_results <- lm_res$gene
    gene_names_in_se <- rownames(se)

    # Find genes that exist in both
    available_genes <- gene_names_in_se[gene_names_in_se %in% gene_names_in_results]

    if (length(available_genes) == 0) {
        stop(sprintf("No genes from lm_res found in rownames(se). \n  Examples from lm_res: %s\n  Examples from SE: %s",
            paste(head(gene_names_in_results, 3), collapse=", "),
            paste(head(gene_names_in_se, 3), collapse=", ")),
            call. = FALSE)
    }

    # Subset SE to only genes that are in results
    se <- se[available_genes, ]

    # Subset results to only genes that are in SE
    lm_res <- lm_res[lm_res$gene %in% rownames(se), ]

    # Extract assay matrix and colData
    mat <- SummarizedExperiment::assay(se, assay_name)
    cdata <- SummarizedExperiment::colData(se)

    if (!condition_col %in% colnames(cdata)) {
        stop(sprintf("Column '%s' not found in colData(se)", condition_col), call. = FALSE)
    }

    # Build sample-to-group mapping
    # colData is duplicated for each q-value (one row per sample x q combination)
    # We need a unique mapping of sample name to group value
    coldata_rownames <- rownames(cdata)
    coldata_sample_names <- sub("_q=.*", "", coldata_rownames)
    
    # Get unique samples and their corresponding group values
    unique_samples <- unique(coldata_sample_names)
    sample_to_group <- character(length(unique_samples))
    names(sample_to_group) <- unique_samples
    
    for (samp in unique_samples) {
        idx <- which(coldata_sample_names == samp)[1]  # Get first occurrence
        sample_to_group[samp] <- as.character(cdata[[condition_col]][idx])
    }

    # Determine which genes to plot
    if (!is.null(genes)) {
        # User provided specific genes
        if (!is.character(genes)) {
            stop("genes must be a character vector of gene names", call. = FALSE)
        }
        top_genes <- genes
    } else {
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
            warning(sprintf("No genes significant at alpha = %g", sig_alpha), call. = FALSE)
            return(NULL)
        }

        # Select top n
        top_genes <- sig_genes$gene[seq_len(min(n_top, nrow(sig_genes)))]
    }

    # Create gene ID to display name mapping from lm_res
    gene_name_map <- setNames(lm_res$gene, lm_res$gene)  # default: use gene ID
    
    # If lm_res has a gene_name column (e.g., from return_model_data), use it
    if ("gene_name" %in% colnames(lm_res)) {
        gene_name_map <- setNames(lm_res$gene_name, lm_res$gene)
    }

    # Helper to build data.frame for a single gene using model_data
    make_gam_plot <- function(g, gene_display_name = NULL) {
        if (!(g %in% rownames(mat))) {
            warning(sprintf("Gene '%s' not found in assay", g), call. = FALSE)
            return(NULL)
        }

        # Use provided gene name, or look it up from mapping, or default to gene ID
        if (is.null(gene_display_name)) {
            if (g %in% names(gene_name_map)) {
                gene_display_name <- gene_name_map[[g]]
            } else {
                gene_display_name <- g
            }
        }
        
        # Extract data for this gene across all columns (samples x q-values)
        # CRITICAL: Extract gene_vals fresh for each gene!
        gene_vals <- mat[g, ]
        col_names_full <- colnames(mat)
        
        # Parse column names to extract sample and q-value
        # Column names are expected to be format: "Sample_q=value"
        col_sample_names <- sub("_q=.*", "", col_names_full)
        col_q_values <- as.numeric(sub(".*_q=", "", col_names_full))
        
        # Look up group for each column using the unique sample-to-group mapping
        col_groups <- unname(sample_to_group[col_sample_names])
        
        # Check for unmapped columns
        if (any(is.na(col_groups))) {
            unmapped_idx <- which(is.na(col_groups))
            warning(sprintf("Cannot map %d columns to groups for gene '%s'; columns not found in colData",
                length(unmapped_idx), g), call. = FALSE)
            return(NULL)
        }
        
        # Build long-format data.frame directly from column annotations
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
            warning(sprintf("No valid data for gene '%s'", g), call. = FALSE)
            return(NULL)
        }

        # Fit GAM per group
        unique_groups <- unique(plot_df$group)

        if (length(unique_groups) < 2) {
            warning(sprintf("Less than 2 groups for gene '%s'", g), call. = FALSE)
            return(NULL)
        }

        # Generate prediction grid
        q_range <- range(plot_df$q, na.rm = TRUE)
        
        # Check if q_range is valid (not all NA)
        if (!is.finite(q_range[1]) || !is.finite(q_range[2])) {
            warning(sprintf("Invalid q values for gene '%s'", g), call. = FALSE)
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
                    # Fit GAM with adaptive k (bases)
                    k <- min(10, max(2, round(nrow(subset_data) / 2)))
                    gam_fit <- mgcv::gam(entropy ~ s(q, k = k), data = subset_data)
                    
                    # Predict
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
                    warning(sprintf("GAM fit failed for gene '%s' group '%s': %s", g, gr, e$message),
                        call. = FALSE)
                }
            )
        }

        if (length(pred_list) == 0) {
            warning(sprintf("No GAM fits succeeded for gene '%s'", g), call. = FALSE)
            return(NULL)
        }

        pred_df <- do.call(rbind, pred_list)

        # CRITICAL FIX: Ensure group is a factor with consistent levels across both dataframes
        group_levels <- sort(unique(c(as.character(plot_df$group), as.character(pred_df$group))))
        plot_df$group <- factor(plot_df$group, levels = group_levels)
        pred_df$group <- factor(pred_df$group, levels = group_levels)

        # Create explicit color mapping for all groups
        # Use .tsenat_palette_blue_red() for harmonized blue-red color scheme
        palette_colors <- .tsenat_palette_blue_red()
        
        # Map each group to a color from the palette
        color_mapping <- c()
        for (i in seq_along(group_levels)) {
            # Cycle through palette if more groups than palette colors
            color_idx <- ((i - 1) %% length(palette_colors)) + 1
            color_mapping[group_levels[i]] <- palette_colors[color_idx]
        }

        # Create plot with explicit color scale
        # Make sure both geoms explicitly get color aesthetic
        p <- ggplot2::ggplot(plot_df, ggplot2::aes(x = q, y = entropy, color = group)) +
            ggplot2::geom_point(
                data = plot_df, 
                ggplot2::aes(x = q, y = entropy, color = group),
                alpha = 0.5, size = 2
            ) +
            ggplot2::geom_line(
                data = pred_df, 
                ggplot2::aes(x = q, y = entropy_fit, color = group, linetype = "GAM fit"), 
                linewidth = 1, alpha = 0.9
            ) +
            ggplot2::scale_color_manual(
                values = color_mapping, 
                name = condition_col,
                breaks = group_levels
            ) +
            ggplot2::scale_linetype_manual(values = c("GAM fit" = 1), name = "") +
            ggplot2::labs(
                x = "q parameter",
                y = "Tsallis entropy",
                title = ifelse(gene_display_name != g, sprintf("%s (%s)", gene_display_name, g), gene_display_name)
            ) +
            .tsenat_theme_spectrum(base_size = 11) +
            ggplot2::theme(
                legend.position = "none"
            )

        return(p)
    }

    # Generate plots for top genes
    plots <- list()
    for (g in top_genes) {
        p <- make_gam_plot(g)
        if (!is.null(p)) {
            plots[[g]] <- p
        }
    }

    if (length(plots) == 0) {
        warning("No valid plots generated", call. = FALSE)
        return(NULL)
    }

    # Arrange plots in a grid and return single combined plot
    n_plots <- length(plots)
    n_cols <- 2
    n_rows <- ceiling(n_plots / n_cols)
    
    # Add margins to plots for spacing, particularly between rows
    plots_with_margins <- lapply(seq_along(plots), function(i) {
        p <- plots[[i]]
        # Add larger bottom margin for plots in the first row to create space before second row
        if (i <= n_cols) {
            p <- p + ggplot2::theme(plot.margin = ggplot2::margin(b = 15, unit = "pt"))
        }
        p
    })
    
    # Extract legend from first plot
    legend <- cowplot::get_legend(plots[[1]] + 
        ggplot2::theme(legend.position = "bottom",
                      legend.title = ggplot2::element_text(size = .tsenat_font_sizes$legend_title),
                      legend.text = ggplot2::element_text(size = .tsenat_font_sizes$legend_text)))
    
    # Create grid without legends
    combined_plot <- cowplot::plot_grid(
        plotlist = plots_with_margins,
        nrow = n_rows,
        ncol = n_cols,
        align = "hv",
        axis = "lr"
    )
    
    # Add main title and subtitle above the grid
    title_plot <- cowplot::ggdraw() + 
        cowplot::draw_label("GAM q-curve: Top genes with group interaction", 
                           fontface = "bold", size = 20, x = 0.5, y = 0.80) +
        cowplot::draw_label("Fitted smooth curves by group", 
                           fontface = "italic", size = 16, x = 0.5, y = 0.25, color = "gray40")
    
    # Combine title, plots, and single legend at bottom
    final_plot <- cowplot::plot_grid(
        title_plot,
        combined_plot,
        legend,
        nrow = 3,
        rel_heights = c(0.12, 1, 0.08)
    )
    
    # Save to file if output_file is provided
    if (!is.null(output_file)) {
        # Use provided dimensions or defaults
        save_width <- if (is.null(width)) 12 else width
        save_height <- if (is.null(height)) 10.3 else height
        ggplot2::ggsave(output_file, plot = final_plot, width = save_width, height = save_height, dpi = 100, create.dir = TRUE)
    }
    
    return(final_plot)
}

# Helpers for plot_top_transcripts internals

.ptt_select_genes_from_res <- function(res, top_n) {
    if (is.null(res)) {
        stop("Either 'gene' or 'res' must be provided")
    }
    if (!("genes" %in% colnames(res))) {
        stop("Provided 'res' must contain a 'genes' column")
    }
    # Look for adjusted p-value columns in order of preference
    if ("padj" %in% colnames(res)) {
        ord <- order(res$padj, na.last = NA)
    } else if ("adjusted_p_values" %in% colnames(res)) {
        ord <- order(res$adjusted_p_values, na.last = NA)
    } else if ("pvalue" %in% colnames(res)) {
        ord <- order(res$pvalue, na.last = NA)
    } else if ("raw_p_values" %in% colnames(res)) {
        ord <- order(res$raw_p_values, na.last = NA)
    } else {
        ord <- seq_len(nrow(res))
    }
    genes_sel <- as.character(res$genes[ord])
    genes_sel <- unique(genes_sel)
    head(genes_sel, top_n)
}

.ptt_infer_samples_from_coldata <- function(coldata, counts, condition_col) {
    if (is.character(coldata) && length(coldata) == 1) {
        if (!file.exists(coldata)) {
            stop("coldata file not found: ", coldata)
        }
        cdf <- utils::read.delim(coldata, header = TRUE, stringsAsFactors = FALSE)
    } else if (is.data.frame(coldata)) {
        cdf <- coldata
    } else {
        stop("`coldata` must be a data.frame or path to a tab-delimited file")
    }

    if (!is.null(rownames(cdf)) && all(colnames(counts) %in% rownames(cdf))) {
        as.character(cdf[colnames(counts), condition_col])
    } else {
        sample_id_cols <- c("sample", "Sample", "sample_id", "id")
        sid <- intersect(sample_id_cols, colnames(cdf))
        if (length(sid) > 0) {
            sid <- sid[1]
            if (!all(colnames(counts) %in% as.character(cdf[[sid]]))) {
                stop("coldata sample id column does not match column names of counts")
            }
            row_ix <- match(colnames(counts), as.character(cdf[[sid]]))
            as.character(cdf[[condition_col]][row_ix])
        } else {
            stop("Could not match `coldata` rows to `counts` columns. Provide `samples` or a row-named `coldata`.")
        }
    }
}

.ptt_read_tx2gene <- function(tx2gene) {
    if (is.null(tx2gene)) {
        stop("`tx2gene` must be provided as a file path or data.frame (or include mapping in metadata of provided SummarizedExperiment)")
    }
    if (is.character(tx2gene) && length(tx2gene) == 1) {
        if (!file.exists(tx2gene)) {
            stop("tx2gene file not found: ", tx2gene)
        }
        mapping <- utils::read.delim(tx2gene, stringsAsFactors = FALSE, header = TRUE)
    } else if (is.data.frame(tx2gene)) {
        mapping <- tx2gene
    } else {
        stop("`tx2gene` must be provided as a file path or data.frame (or include mapping in metadata of provided SummarizedExperiment)")
    }
    if (!all(c("Transcript", "Gen") %in% colnames(mapping))) {
        stop("tx2gene must have columns 'Transcript' and 'Gen'")
    }
    mapping
}

.ptt_make_agg <- function(metric = c("median", "mean", "variance", "iqr")) {
    metric_choice <- match.arg(metric)
    agg_fun <- switch(metric_choice, median = function(x) stats::median(x, na.rm = TRUE),
        mean = function(x) base::mean(x, na.rm = TRUE), variance = function(x) {
            stats::var(x, na.rm = TRUE)
        }, iqr = function(x) stats::IQR(x, na.rm = TRUE))
    agg_label_metric <- if (metric_choice == "iqr") {
        "IQR"
    } else {
        metric_choice
    }
    agg_label <- sprintf("Transcript-level expression with metric %s", agg_label_metric)
    .cnt <- as.integer(getOption("TSENAT.plot_top_counter", 0)) + 1L
    options(TSENAT.plot_top_counter = .cnt)
    agg_label_unique <- agg_label
    list(metric_choice = metric_choice, agg_fun = agg_fun, agg_label_unique = agg_label_unique)
}

.ptt_build_tx_long <- function(gene_single, mapping, counts, samples, top_n) {
    txs <- mapping$Transcript[mapping$Gen == gene_single]
    txs <- intersect(txs, rownames(counts))
    if (length(txs) == 0) {
        stop("No transcripts found for gene: ", gene_single)
    }
    if (!is.null(top_n)) {
        txs <- head(txs, top_n)
    }
    mat <- counts[txs, , drop = FALSE]
    df_all <- as.data.frame(mat)
    df_all$tx <- rownames(mat)
    df_long <- tidyr::pivot_longer(df_all, -tx, names_to = "sample", values_to = "expr")
    df_long$group <- rep(samples, times = length(txs))
    list(df_long = df_long, txs = txs)
}

.ptt_aggregate_df_long <- function(df_long, agg_fun, pseudocount) {
    df_summary <- stats::aggregate(expr ~ tx + group, data = df_long, FUN = agg_fun)
    df_summary$log2expr <- log2(df_summary$expr + pseudocount)
    df_summary$tx <- factor(df_summary$tx, levels = unique(df_summary$tx))
    df_summary
}

.ptt_build_plot_from_summary <- function(df_summary, agg_label_unique, fill_limits = NULL, 
                                        font_scale = 1.0) {
    # Calculate font sizes proportionally to output dimensions
    # Reference: 12x8 inches (96 sq in) uses font_base=11
    # For other sizes, scale base font as: base_font = 11 * sqrt(area/96)
    # This ensures readability is maintained across different output sizes
    
    base_font <- 11 * font_scale
    y_axis_font <- 12 * font_scale
    x_axis_font <- 14 * font_scale
    title_font <- 16 * font_scale
    legend_font <- 9 * font_scale
    
    p <- ggplot2::ggplot(df_summary, ggplot2::aes(x = group, y = tx, fill = log2expr)) +
        ggplot2::geom_tile(color = "black", linewidth = 0.3, width = 0.95, height = 0.92) + 
        ggplot2::geom_vline(xintercept = 1.5, color = "white", linewidth = 1.5) +
        ggplot2::scale_x_discrete(expand = c(0, 0)) +
        ggplot2::scale_y_discrete(expand = c(0, 0)) +
        ggplot2::scale_fill_distiller(
            palette = "Blues",
            na.value = "lightgray", 
            limits = fill_limits,
            name = "log2(expr)"
        ) + 
        .tsenat_theme_base(base_size = base_font) +
        ggplot2::labs(title = agg_label_unique, x = NULL, y = NULL, fill = "log2(expr)") +
        ggplot2::theme(
            axis.text.y = ggplot2::element_text(size = y_axis_font, face = "plain"), 
            axis.text.x = ggplot2::element_text(size = x_axis_font),
            plot.title = ggplot2::element_text(size = title_font, hjust = 0.5, face = "bold"),
            legend.position = "bottom",
            legend.justification = "center",
            legend.key.width = ggplot2::unit(2, "cm"), 
            legend.text = ggplot2::element_text(size = legend_font),
            plot.margin = ggplot2::margin(4, 4, 4, 4)
        ) + 
        ggplot2::guides(fill = ggplot2::guide_colorbar(
            title.position = "top",
            barwidth = 10, 
            barheight = 0.5,
            title.theme = ggplot2::element_text(size = title_font)
        ))
    p
}

.ptt_combine_patchwork <- function(plots, agg_label_unique) {
    # Use 2 columns (2 genes per row) with controlled spacing between rows
    n_cols <- 2
    n_rows <- ceiling(length(plots) / n_cols)
    
    # Build rows of 2 plots each with spacing between columns
    plot_rows <- list()
    for (row in seq_len(n_rows)) {
        start_idx <- (row - 1) * n_cols + 1
        end_idx <- min(row * n_cols, length(plots))
        row_plots <- plots[start_idx:end_idx]
        # Add right margin to first plot to create column spacing
        if (length(row_plots) >= 1) {
            row_plots[[1]] <- row_plots[[1]] + ggplot2::theme(plot.margin = ggplot2::margin(r = 1.0, unit = "cm"))
        }
        # Use patchwork composition (| for horizontal) to avoid scale conflicts
        if (length(row_plots) == 1) {
            row_combined <- row_plots[[1]]
        } else if (length(row_plots) == 2) {
            row_combined <- row_plots[[1]] | row_plots[[2]]  # Horizontal with patchwork
        } else {
            row_combined <- Reduce(function(x, y) x | y, row_plots)
        }
        plot_rows[[row]] <- row_combined
    }
    
    # Combine rows with spacers between them
    combined_elements <- list()
    heights_spec <- c()
    
    for (i in seq_along(plot_rows)) {
        combined_elements[[length(combined_elements) + 1]] <- plot_rows[[i]]
        heights_spec <- c(heights_spec, 1)
        
        if (i < length(plot_rows)) {
            # Add spacer between rows
            spacer <- ggplot2::ggplot() + ggplot2::theme_void()
            combined_elements[[length(combined_elements) + 1]] <- spacer
            heights_spec <- c(heights_spec, 0.17)  # Space between rows (reduced by half)
        }
    }
    
    # Combine all elements
    combined_plots_section <- Reduce(`/`, combined_elements) +
        patchwork::plot_layout(heights = heights_spec)
    
    # Add title spacer above plots (use / for vertical, not | for horizontal)
    spacer <- ggplot2::ggplot() + ggplot2::theme_void()
    title_row <- spacer | patchwork::plot_spacer()
    
    # Combine: title row on top, plot grid below
    combined <- title_row / combined_plots_section +
        patchwork::plot_annotation(
            title = "Transcript level expression",
            subtitle = paste0("Top genes with metric ", agg_label_unique),
            theme = ggplot2::theme(
                plot.title = ggplot2::element_text(hjust = 0.5, size = .tsenat_font_sizes$title, face = "bold", margin = ggplot2::margin(t = 10, b = 10)),
                plot.subtitle = ggplot2::element_text(hjust = 0.5, size = .tsenat_font_sizes$subtitle, face = "italic", margin = ggplot2::margin(t = 5, b = 0.4)),
                legend.position = "bottom"
            )
        ) +
        patchwork::plot_layout(heights = c(0.045, 1), guides = "collect")
    combined
}

.ptt_combine_cowplot <- function(plots, output_file = NULL, agg_label_unique) {
    p_for_legend <- plots[[1]] + ggplot2::theme(legend.position = "bottom")
    legend <- cowplot::get_legend(p_for_legend)
    plots_nolegend <- lapply(plots, function(pp) pp + ggplot2::theme(legend.position = "none"))
    
    # Use 2 columns (2 genes per row), auto-calculate rows
    ncol <- 2
    nrow_val <- ceiling(length(plots_nolegend) / ncol)
    
    grid <- cowplot::plot_grid(plotlist = plots_nolegend, ncol = ncol, nrow = nrow_val, align = "hv")
    title_grob <- cowplot::ggdraw() + cowplot::draw_label("Transcript level expression", fontface = "bold",
        x = 0.5, hjust = 0.5, size = 18)
    subtitle_grob <- cowplot::ggdraw() + cowplot::draw_label(paste0("Top genes with metric ", agg_label_unique), fontface = "italic",
        x = 0.5, hjust = 0.5, size = 14, color = "gray40")
    # Add spacer between title and plots
    spacer_grob <- cowplot::ggdraw() + ggplot2::theme_void()
    result_plot <- cowplot::plot_grid(title_grob, subtitle_grob, spacer_grob, grid, legend, ncol = 1, rel_heights = c(0.05, 0.04,
        0.0015, 1, 0.08), align = "h", axis = "l")
    if (!is.null(output_file)) {
        ggplot2::ggsave(output_file, result_plot)
        invisible(NULL)
    }
    result_plot
}

.ptt_combine_grid <- function(plots, output_file = NULL, agg_label_unique) {
    plots_nolegend <- lapply(plots, function(pp) pp + ggplot2::theme(legend.position = "none"))
    grobs <- lapply(plots_nolegend, ggplot2::ggplotGrob)
    g_full <- ggplot2::ggplotGrob(plots[[1]])
    legend_idx <- which(vapply(g_full$grobs, function(x) x$name, character(1)) ==
        "guide-box")
    if (length(legend_idx)) {
        legend_grob <- g_full$grobs[[legend_idx[1]]]
    } else {
        legend_grob <- NULL
    }
    
    # Default to 2 columns (2 genes per row), adjust for smaller numbers
    ncol <- min(2, length(grobs))
    nrow <- ceiling(length(grobs) / ncol)
    
    # Create heights: title (0.5cm) + plot rows with gaps + legend (0.7cm)
    plot_heights <- list()
    for (i in seq_len(nrow)) {
        plot_heights[[length(plot_heights) + 1]] <- grid::unit(1, "null")
        if (i < nrow) {  # Add gap after each row except the last (reduced by half)
            plot_heights[[length(plot_heights) + 1]] <- grid::unit(0.17, "cm")
        }
    }
    # Combine all heights properly using do.call
    all_heights <- c(list(grid::unit(0.55, "cm")), plot_heights, list(grid::unit(0.7, "cm")))
    heights <- do.call(grid::unit.c, all_heights)
    
    if (!is.null(output_file)) {
        # Adjust PNG dimensions based on layout
        png_width <- 800 * ncol
        png_height <- 480 * nrow
        png(filename = output_file, width = png_width, height = png_height, res = 150)
        .draw_transcript_grid(grobs, agg_label_unique, legend_grob, ncol, heights, to_file = output_file)
        invisible(NULL)
    } else {
        .draw_transcript_grid(grobs, agg_label_unique, legend_grob, ncol, heights)
        invisible(NULL)
    }
}

#' @importFrom ggplot2 ggplot aes geom_col geom_point geom_line scale_y_continuous
#' @importFrom ggplot2 labs theme_minimal theme element_text geom_hline geom_vline
#' @importFrom ggplot2 scale_color_manual scale_shape_manual geom_tile scale_fill_gradient2
NULL

if (getRversion() >= "2.15.1") {
  utils::globalVariables(
    c(
      "dimension",
      "variable",
      "contribution",
      "dim1",
      "dim2",
      "type",
      "coord_x",
      "coord_y"
    )
  )
}



# Internal plot helpers

.tsenat_format_label <- function(lbl) {
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

.tsenat_prepare_ma_plot_df <- function(df, fold_col, mean_cols, x_label, y_label) {
    # Detect x-axis values
    if (length(mean_cols) >= 2) {
        xvals <- rowMeans(df[, mean_cols[seq_len(2)], drop = FALSE], na.rm = TRUE)
        x_label <- x_label %||% paste0(mean_cols[1], " vs ", mean_cols[2])
    } else if (length(mean_cols) == 1) {
        xvals <- as.numeric(df[[mean_cols[1]]])
        x_label <- x_label %||% mean_cols[1]
    } else if ("mean" %in% colnames(df)) {
        xvals <- as.numeric(df$mean)
        x_label <- x_label %||% "Mean"
    } else {
        xvals <- seq_len(nrow(df))
        x_label <- x_label %||% "Index"
    }

    yvals <- as.numeric(df[[fold_col]])

    padj_candidates <- c("adjusted_p_values", "adj_p_value", "adj_p", "padj", "p.adjust")
    padj_col <- intersect(padj_candidates, colnames(df))
    padj_col <- if (length(padj_col)) {
        padj_col[1]
    } else {
        NULL
    }

    padj <- if (!is.null(padj_col)) {
        as.numeric(df[[padj_col]])
    } else {
        rep(1, length(yvals))
    }
    padj[is.na(padj)] <- 1

    sig_flag <- ifelse(abs(yvals) > 0 & padj < 0.05, "significant", "non-significant")

    plot_df <- data.frame(genes = df$genes, x = xvals, y = yvals, padj = padj, significant = sig_flag,
        stringsAsFactors = FALSE)

    list(plot_df = plot_df, x_label = x_label, y_label = y_label)
}

.tsenat_prepare_volcano_df <- function(diff_df, x_col = NULL, padj_col = "adjusted_p_values",
    label_thresh = 0.1, sig_alpha = 0.05, title = NULL) {
    df <- as.data.frame(diff_df)
    cn <- colnames(df)

    # Auto-detect x-axis column if not specified
    if (is.null(x_col)) {
        diff_cols <- grep("_difference$", cn, value = TRUE, ignore.case = TRUE)
        if (length(diff_cols) > 0) {
            x_col <- diff_cols[1]
        } else {
            numeric_cols <- vapply(df, is.numeric, logical(1))
            p_cols <- grep("p_value|p.value", cn, ignore.case = TRUE)
            numeric_cols[p_cols] <- FALSE
            if (any(numeric_cols)) {
                x_col <- cn[which(numeric_cols)[1]]
            } else {
                stop("Could not find suitable column for x-axis. Specify 'x_col' explicitly.")
            }
        }
    }

    # Verify columns
    if (!(x_col %in% cn)) {
        stop(sprintf("Column '%s' not found in diff_df", x_col))
    }
    if (!(padj_col %in% cn)) {
        stop(sprintf("Column '%s' not found in diff_df", padj_col))
    }

    df$xval <- as.numeric(df[[x_col]])
    df$padj <- as.numeric(df[[padj_col]])
    df$padj[is.na(df$padj)] <- 1
    df$padj[df$padj <= 0] <- .Machine$double.xmin

    df$significant <- ifelse(abs(df$xval) >= label_thresh & df$padj < sig_alpha,
        "significant", "non-significant")
    df <- df[is.finite(df$xval) & is.finite(df$padj), ]

    if (nrow(df) == 0) {
        stop("No valid points to plot")
    }

    metric_label <- if (grepl("median", x_col, ignore.case = TRUE)) {
        "Median"
    } else if (grepl("mean", x_col, ignore.case = TRUE)) {
        "Mean"
    } else {
        "Value"
    }

    title_use <- title %||% "Volcano plot: fold-change vs significance"

    x_label_formatted <- .tsenat_format_label(x_col)
    padj_label_formatted <- .tsenat_format_label(padj_col)

    list(df = df, x_col = x_col, padj_col = padj_col, x_label_formatted = x_label_formatted,
        padj_label_formatted = padj_label_formatted, title_use = title_use)
}


#' Plot Tsallis Divergence Effect Size Distribution
#'
#' Generate a histogram visualization of Tsallis divergence effect sizes across genes,
#' showing the distribution of information-theoretic measures of isoform switching.
#'
#' @param interaction_results A data frame containing LMM results merged with per-q divergence estimates.
#'   Must contain columns matching the pattern `effect_size_D_q*` (e.g., `effect_size_D_q0_5`, `effect_size_D_q1_0`).
#'   Typically the result from [effect_sizes_divergence()].
#'
#' @param threshold Numeric. Effect size threshold for visual marking. Default is 0.1 (information-theoretic significance level).
#'
#' @return If ggplot2 is available and `interaction_results` contains valid data, returns a ggplot object.
#'   Otherwise returns NULL invisibly and prints an informative message.
#'
#' @details
#' The function visualizes the distribution of effect sizes using the median q-value's divergence
#' (typically around q=1.0, close to Shannon entropy). The red dashed line marks the default
#' information-theoretic significance threshold of D=0.1.
#'
#' @references
#' - Chanda et al. (2020). Information Theory in Computational Biology. *Entropy*, 22(6), 627.
#' - Tsallis, C. (1988). Possible Generalization of Boltzmann-Gibbs Statistics. *Journal of Statistical Physics*, 52(1), 479-487.
#'
#' @examples
#' # Create example interaction results with divergence effect sizes
#' set.seed(123)
#' interaction_results <- data.frame(
#'   gene = paste0("gene_", 1:20),
#'   effect_size_D_q0_5 = runif(20, 0, 0.3),
#'   effect_size_D_q1_0 = runif(20, 0, 0.2),
#'   effect_size_D_q1_5 = runif(20, 0, 0.25)
#' )
#' 
#' # Plot divergence distribution
#' plot_divergence_distribution(interaction_results, threshold = 0.1)
#'
#' @keywords internal
#' @noRd
plot_divergence_distribution <- function(interaction_results, threshold = 0.1) {
  
  # Check for ggplot2
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    message("ggplot2 package required for plotting. Please install: install.packages('ggplot2')")
    return(invisible(NULL))
  }
  
  # Validate input
  if (is.null(interaction_results) || nrow(interaction_results) == 0) {
    message("Plot not generated: interaction_results is empty or NULL.")
    return(invisible(NULL))
  }
  
  # Find per-q effect size columns
  effect_cols <- grep("^effect_size_D_q", colnames(interaction_results), value = TRUE)
  
  if (length(effect_cols) == 0) {
    message("Plot not generated: no per-q effect size columns found in interaction_results.")
    message("Expected columns like 'effect_size_D_q0_5', 'effect_size_D_q1_0', etc.")
    return(invisible(NULL))
  }
  
  # Use median q effect size for visualization
  median_idx <- ceiling(length(effect_cols) / 2)
  median_col <- effect_cols[median_idx]
  
  # Create visualization of effect size distribution
  p_effect <- ggplot2::ggplot(interaction_results, ggplot2::aes(x = .data[[median_col]])) +
    ggplot2::geom_histogram(binwidth = 0.02, fill = .tsenat_palette_blue_red()[1], alpha = 0.7, color = "black") +
    ggplot2::geom_vline(xintercept = threshold, linetype = "dashed", color = "red", linewidth = 1) +
    ggplot2::labs(
      title = expression("Distribution of Tsallis Divergence (" ~ D[q] ~ ") effect sizes across genes"),
      subtitle = "Information-theoretic measure respecting Tsallis multi-q entropy properties",
      x = bquote("Effect size (Tsallis Divergence" ~ D[q] ~ "; D >" ~ .(threshold) ~ "= meaningful information separation)"),
      y = "Number of genes",
      caption = paste("Red dashed line: D =", threshold, "filtering threshold (information-theoretic significance for q-dependent entropy)")
    ) +
    .tsenat_theme_base(base_size = 11) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(size = .tsenat_font_sizes$title, face = "bold", hjust = 0.5),
      plot.subtitle = ggplot2::element_text(face = "italic", size = .tsenat_font_sizes$subtitle, hjust = 0.5),
      panel.grid.major = ggplot2::element_line(color = "gray90")
    ) +
    ggplot2::annotate("text", x = threshold, y = Inf, 
                      label = paste("Information\nthreshold\n(D=", threshold, ")", sep = ""),
                      vjust = 1.5, hjust = -0.1, color = "red", size = 3.5)
  
  return(p_effect)
}




#' Plot Q-Spectrum Curves for Multiple Top Genes
#'
#' Creates a multi-panel grid comparing per-q divergence profiles across the top
#' N genes identified by LMM interaction analysis. Each panel shows the full q-spectrum
#' divergence curve with the gene name and adjusted p-value in the title.
#'
#' @param eff_res Output from effect size computation OR \code{NULL}.
#'   If provided, must contain `$interaction_results` with columns: gene, adj_p_interaction, per_q_pattern.
#'   If \code{NULL}, uses fallback with lm_res + divergence_results_se.
#'
#' @param lm_res (Optional) Data frame from LMM analysis with columns: gene, adj_p_interaction.
#'   Only used if eff_res is NULL. Must be provided for fallback mode.
#'
#' @param divergence_results_se (Optional) SummarizedExperiment from divergence calculation.
#'   Only used if eff_res is NULL. Must be provided for fallback mode.
#'
#' @param n_genes Integer; number of top genes to plot (default: 9). Genes are sorted by
#'   increasing adjusted p-value (most significant first).
#'
#' @param ncol Integer; number of columns in grid layout (default: 3). Number of rows is
#'   automatically calculated as ceiling(n_genes / ncol).
#'
#' @param verbose Logical; if TRUE, print diagnostic messages (default: TRUE).
#'
#' @param output_file Character or NULL. Optional file path to save the plot as an image.
#'   If provided, the plot will be saved with appropriate dimensions.
#'   Default: NULL (no file output, only return object).
#'
#' @return A ggplot2 object created via \code{patchwork} combining all gene panels,
#'   or NULL if gene data is unavailable. The function automatically handles ggplot2 grid
#'   creation and returns a print-ready object.
#'
#' @details
#' **Input Modes:**
#' - **Mode 1 (Primary)**: Pass eff_res directly (from effect_sizes_divergence)
#' - **Mode 2 (Fallback)**: Pass lm_res + divergence_results_se instead
#'
#' **Gene Filtering:**
#' Genes are ranked by decreasing statistical significance (increasing adj_p_interaction).
#' Only genes with complete per-q divergence data are included. If fewer than n_genes
#' have valid data, the function returns all available genes.
#'
#' **Plot Features:**
#' - Title shows: gene name and adjusted p-value (q-value format)
#' - Per-q divergence curve with point estimates and 95% bootstrap CI bands
#' - Vertical reference line at q=1 (Kullback-Leibler divergence point)
#' - Region labels: "Rare Isoforms" (q<1), "Balanced" (q~=1), "Abundant Isoforms" (q>1)
#' - All plots use consistent ggplot2 styling matching plot_q_spectrum
#'
#' @examples
#' # Plot 4: Multi-gene q-spectrum profiles
#' set.seed(42)
#' n_genes <- 8
#' n_isoforms_per_gene <- 3
#' n_isoforms <- n_genes * n_isoforms_per_gene
#' n_samples_per_group <- 20
#' n_samples <- n_samples_per_group * 2
#' 
#' control_counts <- matrix(rpois(n_isoforms * n_samples_per_group, lambda = 40),
#'   nrow = n_isoforms, ncol = n_samples_per_group)
#' treatment_counts <- matrix(rpois(n_isoforms * n_samples_per_group, lambda = 150),
#'   nrow = n_isoforms, ncol = n_samples_per_group)
#' counts <- cbind(control_counts, treatment_counts)
#' rownames(counts) <- paste0("TX_", 1:n_isoforms)
#' colnames(counts) <- paste0("Sample_", 1:n_samples)
#' 
#' se <- SummarizedExperiment::SummarizedExperiment(assays = list(counts = counts))
#' tx2gene_df <- data.frame(Transcript = rownames(counts),
#'   Gene = rep(paste0("GENE_", 1:n_genes), each = n_isoforms_per_gene))
#' S4Vectors::metadata(se)$tx2gene <- tx2gene_df
#' SummarizedExperiment::colData(se) <- S4Vectors::DataFrame(
#'   sample_id = paste0("Sample_", 1:n_samples),
#'   condition = rep(c("control", "treatment"), each = n_samples_per_group),
#'   row.names = colnames(se))
#' SummarizedExperiment::rowData(se)$transcript_id <- rownames(se)
#' SummarizedExperiment::rowData(se)$gene_id <- tx2gene_df$Gene[match(rownames(se),
#'   tx2gene_df$Transcript)]
#' 
#' analysis <- TSENATAnalysis(se)
#' analysis <- calculate_diversity_s4(analysis, q = c(0.5, 1, 1.5), verbose = FALSE)
#' analysis <- calculate_divergence_s4(analysis, q = c(0.5, 1, 1.5), verbose = FALSE)
#' analysis <- calculate_lm_interaction_s4(analysis,
#'   condition_col = "condition", verbose = FALSE)
#' analysis <- effect_sizes_divergence_s4(analysis, verbose = FALSE)
#' 
#' p <- plot_multi_gene_q_spectrum_s4(analysis, n_genes = 4, verbose = FALSE)
#' if (!is.null(p)) print(p)
#'
#' @seealso \code{\link{calculate_divergence_s4}} for computing divergence values.
#'
#' @export
plot_multi_gene_q_spectrum_s4 <- function(eff_res = NULL, 
                                           lm_res = NULL, 
                                           divergence_results_se = NULL,
                                           n_genes = 9, 
                                           ncol = 3, 
                                           verbose = TRUE,
                                           output_file = NULL) {
  
  require_pkgs(c("ggplot2", "patchwork", "SummarizedExperiment"))
  
  # Handle TSENATAnalysis S4 object
  if (methods::is(eff_res, "TSENATAnalysis")) {
    if (verbose) message("[plot_multi_gene_q_spectrum_s4] Detected TSENATAnalysis object, extracting lm_results and divergence_results...")
    
    # Extract lm_results (contains interaction results)
    if (length(eff_res@lm_results) > 0) {
      lm_res <- eff_res@lm_results[[1]]
      if (verbose) message("[plot_multi_gene_q_spectrum_s4] Extracted lm_results with ", nrow(lm_res), " rows")
    } else {
      stop("TSENATAnalysis object has no lm_results. Run calculate_lm_interaction_s4() first.")
    }
    
    # Extract first divergence result as divergence_results_se
    if (length(eff_res@diversity_results) > 0) {
      divergence_results_se <- eff_res@diversity_results[[1]]
      if (verbose) message("[plot_multi_gene_q_spectrum_s4] Extracted divergence_results with ", nrow(divergence_results_se), " rows")
    } else {
      stop("TSENATAnalysis object has no diversity_results. Run calculate_divergence_s4() first.")
    }
    
    # For S4 mode, set eff_res to NULL so we use lm_res + divergence_results_se
    eff_res <- NULL
  }
  
  # ============================================================================
  # Step 1: Extract gene data from either eff_res or (lm_res + divergence_results_se)
  # ============================================================================
  
  genes_to_plot <- NULL
  per_q_patterns <- NULL
  adj_p_values <- NULL
  
  # Mode 1: Try eff_res first
  if (!is.null(eff_res) && is.list(eff_res)) {
    if (!is.null(eff_res$interaction_results) && nrow(eff_res$interaction_results) > 0) {
      int_res <- eff_res$interaction_results
      
      # Check for required columns (flexible on p-value column name)
      has_gene <- "gene" %in% colnames(int_res)
      has_per_q <- "per_q_pattern" %in% colnames(int_res)
      has_p_adj <- "adj_p_interaction" %in% colnames(int_res)
      has_p_raw <- "p_value_interaction" %in% colnames(int_res)
      has_p_any <- has_p_adj || has_p_raw
      
      # Determine which p-value column to use
      p_col <- NA_character_
      if (has_p_adj) {
        p_col <- "adj_p_interaction"
      } else if (has_p_raw) {
        p_col <- "p_value_interaction"
      }
      
      if (has_gene && has_per_q && has_p_any) {
        # Sort by significance (lowest p-value first = most significant)
        int_res_sorted <- int_res[order(int_res[[p_col]], na.last = TRUE), ]
        
        # Get top n_genes
        int_res_subset <- head(int_res_sorted, n_genes)
        
        # Check if per_q_pattern has valid data
        valid_patterns <- !is.na(int_res_subset$per_q_pattern) & 
                          int_res_subset$per_q_pattern != "" &
                          int_res_subset$per_q_pattern != "NA"
        
        if (any(valid_patterns)) {
          genes_to_plot <- int_res_subset$gene[valid_patterns]
          per_q_patterns <- int_res_subset$per_q_pattern[valid_patterns]
          adj_p_values <- int_res_subset[[p_col]][valid_patterns]
          
          if (verbose) message(sprintf("[plot_multi_gene_q_spectrum_s4] Mode 1: Using eff_res with %s column (%d valid genes)", p_col, length(genes_to_plot)))
        } else {
          if (verbose) message("[plot_multi_gene_q_spectrum_s4] Mode 1 failed: per_q_pattern values are empty or invalid")
        }
      } else {
        if (verbose) {
          message("[plot_multi_gene_q_spectrum_s4] Mode 1 failed: Missing required columns")
          message("  - has 'gene':", has_gene)
          message("  - has 'per_q_pattern':", has_per_q)
          message("  - has 'adj_p_interaction':", has_p_adj)
          message("  - has 'p_value_interaction':", has_p_raw)
        }
      }
    } else {
      if (verbose) message("[plot_multi_gene_q_spectrum_s4] Mode 1 failed: eff_res$interaction_results is NULL or empty")
    }
  } else {
    if (verbose) message("[plot_multi_gene_q_spectrum_s4] Mode 1 failed: eff_res is NULL or not a list")
  }
  
  # Mode 2: Fallback to lm_res + divergence_results_se
  if (is.null(genes_to_plot) && !is.null(lm_res) && !is.null(divergence_results_se)) {
    if (nrow(lm_res) > 0 && nrow(divergence_results_se) > 0) {
      # Check required columns
      if (all(c("gene", "adj_p_interaction") %in% colnames(lm_res))) {
        div_rd <- as.data.frame(rowData(divergence_results_se))
        div_assay <- assay(divergence_results_se)
        
        # Get gene names from divergence rowData
        div_gene_names <- if ("gene_name" %in% colnames(div_rd)) {
          div_rd$gene_name
        } else {
          rownames(div_assay)
        }
        
        if (length(div_gene_names) > 0 && nrow(div_assay) > 0) {
          # Sort lm_res by significance
          lm_sorted <- lm_res[order(lm_res$adj_p_interaction, na.last = TRUE), ]
          top_genes <- head(lm_sorted$gene, n_genes)
          
          # Get indices in divergence_results_se
          gene_indices <- match(top_genes, div_gene_names)
          valid_idx <- !is.na(gene_indices)
          valid_genes <- top_genes[valid_idx]
          
          if (length(valid_genes) > 0) {
            genes_to_plot <- valid_genes
            adj_p_values <- lm_sorted$adj_p_interaction[seq_along(valid_genes)]
            
            # Extract per_q patterns from assay
            per_q_patterns <- character(length(valid_genes))
            for (i in seq_along(valid_genes)) {
              gene_idx <- which(div_gene_names == valid_genes[i])[1]
              if (!is.na(gene_idx)) {
                divs <- div_assay[gene_idx, ]
                per_q_patterns[i] <- paste(divs[!is.na(divs)], collapse = ",")
              }
            }
            
            if (verbose) message("[plot_multi_gene_q_spectrum_s4] Mode 2 (fallback): Using lm_res + divergence_results_se")
          }
        }
      }
    }
  }
  
  # ============================================================================
  # Step 2: Validate and parse gene data
  # ============================================================================
  
  if (is.null(genes_to_plot) || length(genes_to_plot) == 0) {
    if (verbose) {
      message("No valid genes to plot. Check input data:")
      message("  - eff_res provided:", !is.null(eff_res))
      if (!is.null(eff_res)) {
        message("  - eff_res$interaction_results exists:", !is.null(eff_res$interaction_results))
        if (!is.null(eff_res$interaction_results)) {
          message("  - Number of rows:", nrow(eff_res$interaction_results))
          message("  - Has 'per_q_pattern' column:", "per_q_pattern" %in% colnames(eff_res$interaction_results))
        }
      }
      message("  - lm_res provided:", !is.null(lm_res))
      message("  - divergence_results_se provided:", !is.null(divergence_results_se))
    }
    # Return NULL visibly (no invisible) for consistency
    return(NULL)
  }
  
  n_genes_actual <- length(genes_to_plot)
  if (verbose) message(sprintf("Plotting %d genes in %d-column grid", n_genes_actual, ncol))
  
  # ============================================================================
  # Step 3: Create individual q-spectrum plots for each gene
  # ============================================================================
  
  plot_list <- list()
  
  for (i in seq_along(genes_to_plot)) {
    gene_name <- genes_to_plot[i]
    pattern_str <- per_q_patterns[i]
    adj_p <- adj_p_values[i]
    
    # Parse per-q pattern string
    tryCatch({
      per_q_vals <- as.numeric(strsplit(pattern_str, ",")[[1]])
      
      if (length(per_q_vals) == 0 || all(is.na(per_q_vals))) {
        if (verbose) message(sprintf("  Skipping %s: no valid per-q values", gene_name))
        next
      }
      
      # Create plot using plot_q_spectrum logic
      q_vals <- seq(0.1, by = 0.05, length.out = length(per_q_vals))
      
      plot_df <- data.frame(
        q = q_vals,
        divergence = per_q_vals,
        stringsAsFactors = FALSE
      )
      
      # Create base plot
      p <- ggplot2::ggplot(plot_df, ggplot2::aes(x = q, y = divergence)) +
        .tsenat_theme_base(base_size = 11) +
        ggplot2::geom_line(color = "#4575B4", linewidth = 1.2) +
        ggplot2::geom_point(color = "#4575B4", size = 2.8, alpha = 0.8) +
        ggplot2::geom_vline(xintercept = 1, linetype = 3, color = "gray60", linewidth = 0.8, alpha = 0.7) +
        ggplot2::labs(
          title = sprintf("%s", gene_name),
          subtitle = sprintf("adj p = %.2e", adj_p),
          x = "q (Tsallis parameter)",
          y = "Tsallis Divergence D[q]"
        ) +
        ggplot2::theme(
          plot.subtitle = ggplot2::element_text(
            hjust = 0.5, size = 10, color = "gray40",
            margin = ggplot2::margin(b = 8)
          ),
          plot.margin = ggplot2::margin(t = 8, b = 8, l = 6, r = 6),
          panel.grid.major = ggplot2::element_line(color = "gray92", linewidth = 0.25),
          axis.text = ggplot2::element_text(size = .tsenat_font_sizes$axis_text),
          axis.title = ggplot2::element_text(size = .tsenat_font_sizes$axis_title, face = "plain")
        )
      
      plot_list[[i]] <- p
      
    }, error = function(e) {
      if (verbose) message(sprintf("  Failed to plot %s: %s", gene_name, e$message))
    })
  }
  
  # ============================================================================
  # Step 4: Combine plots into grid using patchwork
  # ============================================================================
  
  if (length(plot_list) == 0) {
    if (verbose) {
      message("No valid plots were created.\n",
              "This may occur if:\n",
              "  - per_q_pattern values cannot be parsed as numeric comma-separated strings\n",
              "  - All genes had parsing errors in tryCatch blocks\n",
              "  - Sample size or q-value count was too small")
    }
    # Return NULL visibly (no invisible) for consistency
    return(NULL)
  }
  
  nrow <- ceiling(length(plot_list) / ncol)
  
  # Build layout with spacers between rows to prevent overlap
  layout_plots <- list()
  for (row_idx in seq_len(nrow)) {
    row_start <- (row_idx - 1) * ncol + 1
    row_end <- min(row_idx * ncol, length(plot_list))
    row_plots <- plot_list[row_start:row_end]
    
    # Combine plots in this row horizontally
    if (length(row_plots) == 1) {
      row_combined <- row_plots[[1]]
    } else {
      row_combined <- Reduce(function(x, y) x + y, row_plots)
    }
    
    layout_plots[[length(layout_plots) + 1]] <- row_combined
    
    # Add spacer between rows (except after last row)
    if (row_idx < nrow) {
      layout_plots[[length(layout_plots) + 1]] <- patchwork::plot_spacer()
    }
  }
  
  # Combine all rows with spacers vertically
  combined_plot <- Reduce(function(x, y) x / y, layout_plots) +
                   patchwork::plot_layout(heights = c(rep(c(1, 0.1), nrow - 1), 1), guides = "collect") +
                   patchwork::plot_annotation(
                     title = "Tsallis Divergence q-Spectrum Profiles",
                     subtitle = "Per-q divergence curves for top-ranked genes",
                     theme = ggplot2::theme(
                       plot.title = ggplot2::element_text(hjust = 0.5, face = "bold", size = .tsenat_font_sizes$title, margin = ggplot2::margin(b = 8)),
                       plot.subtitle = ggplot2::element_text(hjust = 0.5, face = "italic", size = .tsenat_font_sizes$subtitle, color = "gray40", margin = ggplot2::margin(b = 12))
                     )
                   )
  
  if (verbose) message(sprintf("[OK] Multi-gene q-spectrum plot created with %d genes", length(plot_list)))
  
  # Save to file if output_file is provided
  if (!is.null(output_file)) {
    ggplot2::ggsave(output_file, plot = combined_plot, width = 12, height = 7.2, dpi = 100, create.dir = TRUE)
  }
  
  return(combined_plot)
}

#' Plot Tsallis Divergence Profiles Across q-Spectrum for Top Genes
#'
#' Visualizes how Tsallis divergence D_q varies across the q-spectrum for selected genes.
#' This reveals which diversity scales (rare vs. abundant isoforms) drive the observed
#' biological differences between groups.
#'
#' @param se A `SummarizedExperiment` returned by `calculate_diversity` with multiple q-values.
#' @param gene Optional character vector of gene names to plot. Can be a single gene name,
#'   a vector of names, or NULL. If NULL, uses top genes from lm_res.
#' @param lm_res Optional data.frame with columns: gene, p_value_interaction, 
#'   adj_p_lmm (or adj_p_interaction). If provided and gene=NULL, top genes are selected by significance.
#' @param readcounts Optional matrix of raw read counts (genes * transcripts) for computing
#'   true Tsallis divergence from isoform distributions. If NULL, uses entropy-based approximation.
#' @param tx2gene_map Optional data.frame mapping transcripts to genes (columns: "transcript", "gene").
#' @param group_col Character name of the column in colData(se) indicating group assignment.
#'   Default: "group".
#' @param n_top Integer. When using lm_res, plot top n_top genes. Default: 3.
#' @param assay_name Character name of the assay containing diversity measures. Default: "diversity".
#' @param arrange_type How to arrange multiple plots. Options: "facet" (faceted grid),
#'   "list" (named list). Default: "facet".
#' @param signed Logical. If TRUE, compute signed divergence (mean2 - mean1); 
#'   positive = group2 higher, negative = group1 higher. If FALSE (default), 
#'   absolute divergence. Signed divergence reveals directional differences.
#'
#' @return If arrange_type="facet", returns a single ggplot with facets by gene.
#'   If arrange_type="list", returns a named list of ggplot objects (one per gene).
#'   When signed=TRUE, negative values indicate group1 preference, positive indicate group2 preference.
#'
#' @details
#' The plot shows:
#' - X-axis: q value (from 0.1 to ~3, depending on SE)
#' - Y-axis: Tsallis divergence D_q between control and treatment groups
#' - Shape: Profile reveals which q-ranges drive divergence:
#'   - Low q (<1): rare isoform divergence
#'   - q~=1: KL divergence region
#'   - High q (>1): dominant isoform divergence
#'
#' @importFrom ggplot2 ggplot aes geom_line geom_point facet_wrap labs theme_minimal
#'   theme element_text scale_color_manual scale_size_manual
#' @importFrom dplyr group_by summarise
#' @importFrom SummarizedExperiment colData assay
#'
#' @examples
#' library(SummarizedExperiment)
#' set.seed(42)
#' # Create sample Tsallis divergence data
#' ts_se <- SummarizedExperiment(
#'   assays = list(divergence = matrix(rnorm(100, mean=2, sd=0.5), nrow=10, ncol=10)),
#'   rowData = data.frame(gene = paste0("gene_", 1:10)),
#'   colData = data.frame(condition = rep(c("A", "B"), 5))
#' )
#' # p <- plot_tsallis_divergence_profile(
#' #   ts_se, gene = c("gene_1", "gene_2")
#' # )
#'
#' @keywords internal
#' @noRd
plot_tsallis_divergence_profile <- function(se,
                                            gene = NULL,
                                            lm_res = NULL,
                                            readcounts = NULL,
                                            tx2gene_map = NULL,
                                            group_col = "group",
                                            n_top = 3,
                                            assay_name = "diversity",
                                            arrange_type = c("facet", "list"),
                                            signed = TRUE) {
    require_pkgs(c("ggplot2", "dplyr", "SummarizedExperiment"))
    arrange_type <- match.arg(arrange_type)

    # Validate SE
    if (!inherits(se, "SummarizedExperiment")) stop("se must be a SummarizedExperiment")
    if (!assay_name %in% names(SummarizedExperiment::assays(se))) {
        stop("Assay '", assay_name, "' not found in se")
    }

    # Determine which genes to plot
    if (is.null(gene)) {
        if (is.null(lm_res)) stop("Either 'gene' or 'lm_res' must be provided")
        if (!is.data.frame(lm_res)) stop("'lm_res' must be a data.frame")
        if (!("gene" %in% colnames(lm_res))) stop("'lm_res' must have a 'gene' column")
        
        # Use FDR-adjusted p-value if available, else raw p-value
        p_col <- if ("adj_p_lmm" %in% colnames(lm_res)) "adj_p_lmm" 
                 else if ("adj_p_interaction" %in% colnames(lm_res)) "adj_p_interaction"
                 else if ("p_value_interaction" %in% colnames(lm_res)) "p_value_interaction"
                 else stop("'lm_res' must contain p-value column")
        
        genes_ordered <- unique(as.character(lm_res$gene[order(lm_res[[p_col]])]))
        genes <- head(genes_ordered, n_top)
    } else {
        genes <- as.character(unlist(gene))
    }

    if (length(genes) == 0) stop("No genes selected for plotting")

    # Extract metadata
    col_data <- as.data.frame(SummarizedExperiment::colData(se))
    if (!group_col %in% colnames(col_data)) {
        stop("Column '", group_col, "' not found in colData(se)")
    }

    # Parse column names to extract q values
    col_names <- colnames(se)
    extract_q <- function(name) {
        if (grepl("_q=", name)) {
            as.numeric(gsub(".*_q=", "", name))
        } else {
            NA
        }
    }
    q_values <- vapply(col_names, extract_q, FUN.VALUE = numeric(1))
    unique_q <- sort(unique(q_values[!is.na(q_values)]))

    if (length(unique_q) < 2) {
        stop("SE must contain multiple q-values in column names (format: *_q=0.5)")
    }

    # Get group information
    groups <- unique(col_data[[group_col]])
    if (length(groups) != 2) {
        stop("Exactly 2 groups required in '", group_col, "' column; found: ", 
             paste(groups, collapse = ", "))
    }

    # Helper: calculate divergence for a single gene at a specific q
    calc_div_for_gene_q <- function(gene_name, q_val) {
        # Get columns matching this q value
        cols_q <- which(q_values == q_val)
        if (length(cols_q) == 0) return(NA)

        # Extract entropy data for this gene at this q
        diversity_matrix <- SummarizedExperiment::assay(se, assay_name)
        if (!gene_name %in% rownames(diversity_matrix)) return(NA)

        entropy_vals <- diversity_matrix[gene_name, cols_q]
        group_vals <- col_data[[group_col]][cols_q]

        # Compute divergence using calculate_tsallis_divergence_paired_gene if readcounts available
        if (!is.null(readcounts) && !is.null(tx2gene_map)) {
            # TIER 1: True Tsallis divergence from isoform distributions
            tryCatch({
                # Prepare gene subset data
                gene_cols <- cols_q
                entropy_pred <- data.frame(
                    entropy_pred = entropy_vals,
                    q = q_val,
                    group = group_vals,
                    stringsAsFactors = FALSE
                )

                div <- calculate_tsallis_divergence_paired_gene(
                    gene_name = gene_name,
                    gene_data = data.frame(
                        entropy = entropy_vals,
                        q = q_val,
                        group = group_vals,
                        sample = colnames(se)[cols_q],
                        stringsAsFactors = FALSE
                    ),
                    readcounts = readcounts,
                    tx2gene_map = tx2gene_map,
                    entropy_pred = entropy_pred,
                    group_levels = groups
                )
                # Return signed or absolute divergence based on parameter
                return(if (signed) div else abs(div))
            }, error = function(e) {
                return(NA)
            })
        }

        # TIER 2: Entropy-based approximation (fallback)
        # Divergence ~= difference in mean entropy between groups
        group1_vals <- entropy_vals[group_vals == groups[1]]
        group2_vals <- entropy_vals[group_vals == groups[2]]

        if (length(group1_vals) == 0 || length(group2_vals) == 0) return(NA)

        mean1 <- mean(group1_vals, na.rm = TRUE)
        mean2 <- mean(group2_vals, na.rm = TRUE)
        
        # Return signed or absolute divergence based on parameter
        if (signed) {
            div_approx <- mean2 - mean1  # Signed: positive if group2 > group1
        } else {
            div_approx <- max(0, abs(mean1 - mean2))  # Unsigned: always non-negative
        }

        return(div_approx)
    }

    # Calculate divergence for all genes * q combinations
    plot_data_list <- list()
    for (gene_name in genes) {
        divergences <- vapply(unique_q, function(q) calc_div_for_gene_q(gene_name, q), FUN.VALUE = numeric(1))
        df_gene <- data.frame(
            gene = gene_name,
            q = unique_q,
            divergence = divergences,
            stringsAsFactors = FALSE
        )
        plot_data_list[[gene_name]] <- df_gene
    }

    all_plot_data <- do.call(rbind, plot_data_list)
    rownames(all_plot_data) <- NULL

    # Remove NA divergences
    all_plot_data <- all_plot_data[!is.na(all_plot_data$divergence), ]

    if (nrow(all_plot_data) == 0) {
        stop("No valid divergence values computed; check readcounts/tx2gene_map or data structure")
    }

    # Add direction indicator for signed divergence
    if (signed) {
        all_plot_data$direction <- ifelse(all_plot_data$divergence > 0, 
                                         paste0("Positive: ", groups[2], " higher"),
                                         paste0("Negative: ", groups[1], " higher"))
    }

    # Build plot
    if (signed) {
        # Signed divergence plot with color-coded direction
        p <- ggplot2::ggplot(all_plot_data, ggplot2::aes(x = q, y = divergence, color = direction)) +
            ggplot2::geom_hline(yintercept = 0, linetype = "dashed", color = "gray50", linewidth = 0.8) +
            ggplot2::geom_line(linewidth = 1.1) +
            ggplot2::geom_point(size = 3, alpha = 0.7) +
            ggplot2::scale_color_manual(
                name = "Divergence Direction:",
                values = setNames(c("#E63946", "#1D3557"), 
                                 c(paste0("Negative: ", groups[1], " higher"),
                                   paste0("Positive: ", groups[2], " higher")))
            ) +
            ggplot2::labs(
                title = "Tsallis Divergence Profile: Directional (Signed)",
                x = "q value (diversity scale parameter)",
                y = "Divergence D[q] (Positive = Right Group Higher, Negative = Left Group Higher)"
            ) +
            .tsenat_theme_spectrum(base_size = 11) +
            ggplot2::theme(
                plot.title = ggplot2::element_text(
                    size = .tsenat_font_sizes$title
                ),
                axis.title = ggplot2::element_text(size = .tsenat_font_sizes$axis_title),
                axis.text = ggplot2::element_text(size = .tsenat_font_sizes$axis_text)
            )
    } else {
        # Absolute divergence plot (original)
        p <- ggplot2::ggplot(all_plot_data, ggplot2::aes(x = q, y = divergence, color = gene)) +
            ggplot2::geom_line(linewidth = 1.1) +
            ggplot2::geom_point(size = 3, alpha = 0.7) +
            ggplot2::labs(
                title = "Tsallis Divergence Profile Across q-Spectrum",
                x = "q value (diversity scale parameter)",
                y = "Divergence D[q] (Absolute)",
                color = "Gene"
            ) +
            .tsenat_theme_spectrum(base_size = 11) +
            ggplot2::theme(
                plot.title = ggplot2::element_text(
                    size = .tsenat_font_sizes$title
                ),
                axis.title = ggplot2::element_text(size = .tsenat_font_sizes$axis_title),
                axis.text = ggplot2::element_text(size = .tsenat_font_sizes$axis_text)
            )
    }

    # Conditionally apply faceting
    if (arrange_type == "facet" && length(unique(all_plot_data$gene)) > 1) {
        if (signed) {
            p <- p + ggplot2::facet_wrap(~gene, scales = "free_y") +
                ggplot2::theme(legend.position = "top")
        } else {
            p <- p + ggplot2::facet_wrap(~gene, scales = "free_y") +
                ggplot2::theme(legend.position = "top")
        }
    }

    # Return as list if requested
    if (arrange_type == "list") {
        plots <- list()
        for (gene_name in genes) {
            df_gene <- all_plot_data[all_plot_data$gene == gene_name, ]
            if (nrow(df_gene) > 0) {
                if (signed) {
                    p_gene <- ggplot2::ggplot(df_gene, ggplot2::aes(x = q, y = divergence, fill = direction, color = direction)) +
                        ggplot2::geom_hline(yintercept = 0, linetype = "dashed", color = "gray50", linewidth = 0.8) +
                        ggplot2::geom_line(linewidth = 1.2) +
                        ggplot2::geom_point(size = 3.5, alpha = 0.8) +
                        ggplot2::scale_color_manual(
                            name = "Divergence Direction:",
                            values = setNames(c("#E63946", "#1D3557"), 
                                             c(paste0("Negative: ", groups[1], " higher"),
                                               paste0("Positive: ", groups[2], " higher")))
                        ) +
                        ggplot2::scale_fill_manual(
                            name = "Divergence Direction:",
                            values = setNames(c("#E63946", "#1D3557"), 
                                             c(paste0("Negative: ", groups[1], " higher"),
                                               paste0("Positive: ", groups[2], " higher")))
                        ) +
                        ggplot2::labs(
                            title = paste("Divergence Profile (Signed):", gene_name),
                            x = "q value",
                            y = expression("Divergence D[q]"),
                            subtitle = paste0("Red = ", groups[1], " higher | Blue = ", groups[2], " higher")
                        ) +
                        .tsenat_theme_spectrum(base_size = 11) +
                        ggplot2::theme(
                            plot.title = ggplot2::element_text(hjust = 0.5, size = .tsenat_font_sizes$title, face = "bold"),
                            panel.grid.minor = ggplot2::element_blank(),
                            legend.position = "right"
                        )
                } else {
                    p_gene <- ggplot2::ggplot(df_gene, ggplot2::aes(x = q, y = divergence)) +
                        ggplot2::geom_line(color = "#2E86AB", linewidth = 1.2) +
                        ggplot2::geom_point(color = "#2E86AB", size = 3.5, alpha = 0.8) +
                        ggplot2::labs(
                            title = paste("Divergence Profile:", gene_name),
                            x = "q value",
                            y = "Divergence D[q] (Absolute)"
                        ) +
                        .tsenat_theme_spectrum(base_size = 11) +
                        ggplot2::theme(
                            plot.title = ggplot2::element_text(hjust = 0.5, size = .tsenat_font_sizes$title, face = "bold"),
                            panel.grid.minor = ggplot2::element_blank()
                    )
                }
                plots[[gene_name]] <- p_gene
            }
        }
        return(plots)
    }

    return(p)
}

#' Plot Global Divergence q-Curve Across All Genes
#'
#' Visualizes the average (mean/median) Tsallis divergence D_q across all genes
#' as a function of q-value. This provides a **global view** of which diversity scales
#' (rare vs. abundant isoforms) drive the most divergence on average across the dataset,
#' complementing gene-specific divergence profiles.
#'
#' @param divergence_results_se A `SummarizedExperiment` containing pre-computed divergence values.
#'   Rows = genes, columns = q-values. Column names should indicate q-values (e.g., "q_0.5", "q_1.0").
#' @param gene Optional character. If provided, plot divergence spectrum for this specific gene.
#'   If NULL, plot global divergence curve (aggregated across all genes).
#' @param lm_res Optional data.frame with columns for gene identifiers and p-values. Used to select
#'   top genes when gene = NULL and lm_res is provided. Default: NULL.
#' @param n_genes Integer; number of top genes to plot when showing multi-gene spectra (default: 4).
#'   Genes are sorted by p-value significance (lowest p-values first).
#' @param ncol Integer; number of columns in grid layout for multi-gene plots (default: 2).
#'   Number of rows is automatically calculated as ceiling(n_genes / ncol).
#' @param metric Character. Summary statistic for global curve: "median" or "mean". Default: "median".
#'   Only used when gene = NULL.
#' @param variability_metric Character. Error bar type for global curve: "sd" or "iqr". Default: "iqr".
#'   Only used when gene = NULL.
#'
#' @return A `ggplot` object. Gene-specific calls return a line plot.
#'   Global calls return an aggregated curve with variability bands.
#'
#' @details
#' **Gene-specific mode (gene provided)**:
#' - Extracts divergence values for the specified gene across all q-values
#' - Plots as a line chart with points
#' - Reveals whether this gene shows q-dependent divergence patterns
#'
#' **Global mode (gene = NULL)**:
#' - Aggregates divergence across all genes at each q-value
#' - Shows which diversity scales (q-values) drive the most divergence on average
#' - Useful for identifying dominant biological mechanisms (rare vs. abundant isoform driven)
#'
#' **Interpretation**: Compare with `plot_tsallis_q_curve` (entropy) to understand
#' the relationship between entropy changes and divergence patterns.
#'
#' @importFrom ggplot2 ggplot aes geom_line geom_point geom_ribbon labs theme_minimal element_text
#' @importFrom SummarizedExperiment assay
#'
#' @examples
#' # Create synthetic divergence data
#' set.seed(123)
#' divergence_matrix <- matrix(
#'   rnorm(80, mean = 0.5, sd = 0.1),
#'   nrow = 20, ncol = 4
#' )
#' rownames(divergence_matrix) <- paste0("gene_", 1:20)
#' colnames(divergence_matrix) <- c("q_0.5", "q_1.0", "q_1.5", "q_2.0")
#' divergence_se <- SummarizedExperiment::SummarizedExperiment(
#'   assays = list(divergence = divergence_matrix)
#' )
#' 
#' # Global divergence curve (all genes aggregated)
#' p_global <- plot_divergence_spectrum(divergence_se)
#' 
#' # Gene-specific divergence spectrum
#' p_gene <- plot_divergence_spectrum(divergence_se, gene = "gene_1")
#'
#' @keywords internal
#' @noRd
plot_divergence_spectrum <- function(divergence_results_se,
                                     gene = NULL,
                                     lm_res = NULL,
                                     n_genes = 4,
                                     ncol = 2,
                                     metric = c("median", "mean"),
                                     variability_metric = c("iqr", "sd")) {
    require_pkgs(c("ggplot2", "patchwork", "SummarizedExperiment"))
    metric <- match.arg(metric)
    variability_metric <- match.arg(variability_metric)
    
    # Validate input
    if (!inherits(divergence_results_se, "SummarizedExperiment")) {
        stop("divergence_results_se must be a SummarizedExperiment")
    }
    
    # Extract divergence matrix and get gene names
    div_mat <- SummarizedExperiment::assay(divergence_results_se, 1)
    if (is.null(div_mat) || ncol(div_mat) == 0) {
        stop("divergence_results_se has no assays or is empty")
    }
    
    # Get gene names from rowData if available
    rd <- SummarizedExperiment::rowData(divergence_results_se)
    gene_names <- if (!is.null(rd) && "gene_name" %in% colnames(rd)) {
        rd$gene_name
    } else {
        rownames(div_mat)
    }
    
    # Extract q values from column names
    col_names <- colnames(div_mat)
    extracted_q <- gsub(".*q[_=]?", "", col_names)
    q_vals <- as.numeric(extracted_q)
    
    if (all(is.na(q_vals))) {
        # Fallback: assume sequential q-values
        q_vals <- seq(0.5, by = 0.5, length.out = ncol(div_mat))
    }
    
    # Sort by q
    sort_idx <- order(q_vals)
    q_vals_sorted <- q_vals[sort_idx]
    div_mat_sorted <- div_mat[, sort_idx]
    
    # =========================================================================
    # Case 1: Gene-specific spectrum (single gene)
    # =========================================================================
    if (!is.null(gene)) {
        if (!gene %in% c(gene_names, rownames(div_mat_sorted))) {
            stop("Gene '", gene, "' not found in divergence_results_se")
        }
        
        # Find index of gene
        if (gene %in% gene_names) {
            gene_idx <- which(gene_names == gene)[1]
        } else {
            gene_idx <- which(rownames(div_mat_sorted) == gene)[1]
        }
        
        gene_div <- as.numeric(div_mat_sorted[gene_idx, ])
        
        plot_df <- data.frame(
            q = q_vals_sorted,
            divergence = gene_div,
            stringsAsFactors = FALSE
        )
        
        p <- ggplot2::ggplot(plot_df, ggplot2::aes(x = q, y = divergence)) +
            ggplot2::geom_line(color = "#4575B4", linewidth = 1.2) +
            ggplot2::geom_point(color = "#4575B4", size = 3.5, alpha = 0.8) +
            ggplot2::labs(
                title = paste("Divergence Spectrum:", gene),
                x = "q value (diversity scale parameter)",
                y = "Tsallis Divergence D_q"
            ) +
            .tsenat_theme_base(base_size = 11) +
            ggplot2::theme(
plot.title = ggplot2::element_text(size = .tsenat_font_sizes$title, face = "bold", hjust = 0.5)
            )
        
        return(p)
    }
    
    # =========================================================================
    # Case 2: Top N genes spectra (multi-gene faceted plot)
    # =========================================================================
    if (!is.null(lm_res)) {
        # Get top genes ranked by interaction p-value
        if (!is.data.frame(lm_res)) {
            stop("lm_res must be a data.frame with gene and p-value columns")
        }
        
        # Find gene identifier column
        gene_col <- if ("gene" %in% colnames(lm_res)) {
            "gene"
        } else if ("gene_name" %in% colnames(lm_res)) {
            "gene_name"
        } else if ("gene_id" %in% colnames(lm_res)) {
            "gene_id"
        } else {
            stop("lm_res must have 'gene', 'gene_name', or 'gene_id' column")
        }
        
        # Find p-value column
        p_col <- if ("adj_p_interaction" %in% colnames(lm_res)) {
            "adj_p_interaction"
        } else if ("p_interaction" %in% colnames(lm_res)) {
            "p_interaction"
        } else if ("adj_p_value" %in% colnames(lm_res)) {
            "adj_p_value"
        } else if ("p_value" %in% colnames(lm_res)) {
            "p_value"
        } else {
            stop("lm_res must have a p-value column (adj_p_interaction, p_interaction, etc.)")
        }
        
        # Sort by p-value and get top genes
        lm_sorted <- lm_res[order(lm_res[[p_col]], na.last = TRUE), ]
        top_genes_vec <- head(lm_sorted[[gene_col]], n_genes)
        
        # Match to divergence matrix
        gene_indices <- match(top_genes_vec, c(gene_names, rownames(div_mat_sorted)))
        gene_indices <- gene_indices[!is.na(gene_indices)]
        
        if (length(gene_indices) == 0) {
            warning("No genes from lm_res found in divergence_results_se. Falling back to global spectrum.")
            top_genes_vec <- NULL
        } else {
            # Build multi-gene plotting data frame
            plot_list <- list()
            
            for (i in seq_along(gene_indices)) {
                gene_idx <- gene_indices[i]
                gene_name_i <- gene_names[gene_idx]
                gene_div <- as.numeric(div_mat_sorted[gene_idx, ])
                p_val <- lm_sorted[[p_col]][i]
                
                plot_list[[i]] <- data.frame(
                    q = q_vals_sorted,
                    divergence = gene_div,
                    gene = gene_name_i,
                    p_value = p_val,
                    stringsAsFactors = FALSE
                )
            }
            
            multi_gene_df <- do.call(rbind, plot_list)
            
            # Extract confidence intervals from rowData for each gene at each q-value
            rd <- SummarizedExperiment::rowData(divergence_results_se)
            ci_data_list <- list()
            
            for (idx in seq_along(gene_indices)) {
                gene_idx <- gene_indices[idx]
                gene_name_i <- gene_names[gene_idx]
                
                # Extract CIs for this gene across all q-values
                ci_lower <- numeric(length(q_vals_sorted))
                ci_upper <- numeric(length(q_vals_sorted))
                
                for (j in seq_along(q_vals_sorted)) {
                    q_val <- q_vals_sorted[j]
                    lower_col <- paste0("lower_ci_q", q_val)
                    upper_col <- paste0("upper_ci_q", q_val)
                    
                    if (!is.null(rd) && lower_col %in% colnames(rd) && upper_col %in% colnames(rd)) {
                        ci_lower[j] <- rd[[lower_col]][gene_idx]
                        ci_upper[j] <- rd[[upper_col]][gene_idx]
                    } else {
                        ci_lower[j] <- NA_real_
                        ci_upper[j] <- NA_real_
                    }
                }
                
                ci_data_list[[idx]] <- data.frame(
                    q = q_vals_sorted,
                    lower = ci_lower,
                    upper = ci_upper,
                    gene = gene_name_i,
                    stringsAsFactors = FALSE
                )
            }
            
            if (length(ci_data_list) > 0) {
                ci_df <- do.call(rbind, ci_data_list)
                rownames(ci_df) <- NULL
            } else {
                ci_df <- NULL
            }
            
            # Sort genes by p-value for proper facet order (most significant first)
            gene_p_values <- multi_gene_df[!duplicated(multi_gene_df$gene), c("gene", "p_value")]
            gene_p_values <- gene_p_values[order(gene_p_values$p_value), ]
            gene_order <- gene_p_values$gene
            multi_gene_df$gene <- factor(multi_gene_df$gene, levels = gene_order)
            
            # Create faceted plot with individual gene confidence intervals
            p <- ggplot2::ggplot(multi_gene_df, ggplot2::aes(x = q, y = divergence)) +
                ggplot2::facet_wrap(~ gene, ncol = ncol, scales = "free_y")
            
            # Add CI ribbons if available
            if (!is.null(ci_df)) {
                ci_df$gene <- factor(ci_df$gene, levels = gene_order)
                p <- p + ggplot2::geom_ribbon(
                    data = ci_df,
                    ggplot2::aes(x = q, ymin = lower, ymax = upper),
                    inherit.aes = FALSE,
                    alpha = 0.15,
                    fill = "#4575B4",
                    color = NA
                )
            }
            
            p <- p +
                ggplot2::geom_line(color = "#4575B4", linewidth = 1.2, alpha = 0.8) +
                ggplot2::geom_point(color = "#4575B4", size = 3, alpha = 0.8) +
                ggplot2::labs(
                    title = "Divergence Spectra: Per-gene Comparisons",
                    subtitle = paste0("Ranked by interaction significance (", metric, ")"),
                    x = "q value (diversity scale parameter)",
                    y = expression("Divergence D[q]")
                ) +
                .tsenat_theme_base(base_size = 11) +
                ggplot2::theme(
                    plot.title = ggplot2::element_text(size = .tsenat_font_sizes$title, face = "bold", hjust = 0.5),
                    plot.subtitle = ggplot2::element_text(face = "italic", size = .tsenat_font_sizes$subtitle, hjust = 0.5),
                    panel.spacing = ggplot2::unit(1.5, "lines"),
                    strip.text = ggplot2::element_text(face = "bold", size = .tsenat_font_sizes$subtitle)
                )
            
            return(p)
        }
    }
    
    # =========================================================================
    # Case 3: Global divergence curve (all genes aggregated)
    # =========================================================================
    
    # Aggregate across genes at each q
    if (variability_metric == "iqr") {
        summary_stats <- data.frame(
            q = q_vals_sorted,
            central = apply(div_mat_sorted, 2, function(x) {
                if (metric == "median") median(x, na.rm = TRUE) else mean(x, na.rm = TRUE)
            }),
            spread = apply(div_mat_sorted, 2, function(x) {
                stats::IQR(x, na.rm = TRUE)
            }),
            stringsAsFactors = FALSE
        )
        spread_factor <- 1/2  # IQR/2 for symmetric ribbon
        spread_label <- "IQR"
    } else {  # sd
        summary_stats <- data.frame(
            q = q_vals_sorted,
            central = apply(div_mat_sorted, 2, function(x) {
                if (metric == "median") median(x, na.rm = TRUE) else mean(x, na.rm = TRUE)
            }),
            spread = apply(div_mat_sorted, 2, function(x) {
                sqrt(stats::var(x, na.rm = TRUE))
            }),
            stringsAsFactors = FALSE
        )
        spread_factor <- 1
        spread_label <- "SD"
    }
    
    metric_label <- if (metric == "median") "Median" else "Mean"
    p <- ggplot2::ggplot(summary_stats, ggplot2::aes(x = q, y = central)) +
        ggplot2::geom_ribbon(
            ggplot2::aes(ymin = central - spread * spread_factor, 
                        ymax = central + spread * spread_factor),
            alpha = 0.1,
            fill = "#4575B4",
            color = NA
        ) +
        ggplot2::geom_line(color = "#4575B4", linewidth = 1.3) +
        ggplot2::geom_point(color = "#4575B4", size = 3.5, alpha = 0.8) +
        ggplot2::labs(
            title = expression("Global Divergence Spectrum: Average " * D[q] * " Across All Genes"),
            x = "q value (diversity scale parameter)",
            y = expression("Divergence D[q]"),
            subtitle = paste0(metric_label, " +/- ", spread_label, " (", nrow(div_mat_sorted), " genes)")
        ) +
        .tsenat_theme_base(base_size = 11) +
        ggplot2::theme(
            plot.title = ggplot2::element_text(size = .tsenat_font_sizes$title, face = "bold", hjust = 0.5),
            plot.subtitle = ggplot2::element_text(face = "italic", size = .tsenat_font_sizes$subtitle, hjust = 0.5)
        )
    
    return(p)
}

#' Plot Divergence Spectrum Heatmaps (Multi-q Transcript Switching)
#'
#' Creates combined heatmap panels showing delta influence (transcript switching magnitude)
#' across multiple q-values (diversity scales) for selected genes. Each heatmap shows
#' how transcript importance differs between conditions (delta_influence) across the
#' q-spectrum from 0.01 to 2.0.
#'
#' @param switching_results A multi-q switching analysis result from
#'   \code{\link{jackknife_isoform_switching}(q = c(...))}. Must include
#'   gene_ids, gene_name_map, and per-gene results for each q-value.
#' @param n_genes Numeric: number of top genes to visualize (default 4).
#'   When \code{lm_results} is provided, selects the n genes with lowest p-values.
#'   Otherwise, selects the first n genes from results.
#' @param lm_results Optional data.frame from \code{\link{calculate_lm_interaction}()}
#'   containing gene interaction statistics. Should have columns for gene identifiers
#'   ('gene_name' or 'gene_id') and p-values ('p_interaction' or 'adj_p_interaction').
#'   If provided, genes are ranked by p-value significance for selection of top genes.
#' @param verbose Logical; if TRUE, print detailed validation report of heatmap data
#'   including which genes were included and any skipped due to insufficient data.
#'   Default: FALSE (no validation output).
#' @param cellwidth Numeric; width of heatmap cells in pixels (default: 35).
#'   Following pheatmap best practices for publication-quality heatmaps.
#'   Larger values (50+) make cells more visible but reduce number of visible transcripts.
#' @param cellheight Numeric; height of heatmap cells in pixels (default: 10.25).
#'   Following pheatmap best practices. Smaller values allow more q-values to be visible.
#' @param fontsize Numeric; font size in points for heatmap labels (default: 11).
#'   Following pheatmap best practices for publication-quality figures. Applies to
#'   row labels (q-values) and column labels (transcript IDs).
#'
#' @return Character path to saved PNG file containing the combined heatmaps.
#'   The plot is automatically saved to a temporary file and can be displayed
#'   in R Markdown with \code{knitr::include_graphics()}.
#'
#' @details
#' **Heatmap interpretation:**
#'
#' - \bold{Rows}: Different q-values from 0.01 (rare isoforms) to 2.0 (dominant isoforms)
#' - \bold{Columns}: Individual transcripts of each gene
#' - \bold{Color scale}: Blue (negative delta_influence) = transcript more important in second condition;
#'   Red (positive delta_influence) = transcript more important in first condition;
#'   White = no switching effect
#' - \bold{Intensity}: Darker colors indicate stronger switching magnitude
#'
#' **Layout:**
#' - Multiple genes displayed in separate panels (up to 2 per row)
#' - Panels combined into single PNG for reproducible visualization
#' - Outliers (>95th percentile) are capped for better color contrast
#' - NaN and Inf values treated as missing (light gray)
#'
#' **Use case:**
#' Identify whether transcript switching is consistent across diversity scales
#' or scale-dependent. Genes with similar patterns across q-values show robust
#' isoform shifts; genes with varying patterns across q indicate q-dependent
#' switching driven by rare vs. abundant isoforms.
#'
#' @examples
#' # Example: Create synthetic multi-q switching results
#' # For real analysis, use jackknife_isoform_switching() output
#' set.seed(123)
#' gene_names <- paste0("gene_", 1:4)
#' names(gene_names) <- 1:4
#' 
#' # Create multi-q results structure
#' q_values <- c(0.5, 1.0, 1.5)
#' switching_results <- structure(
#'   list(
#'     q_0_50 = list(
#'       gene_ids = 1:4,
#'       gene_name_map = gene_names,
#'       switching_results = list(
#'         list(gene_id = 1, delta_influence = matrix(rnorm(30, 0, 0.2), 5, 6)),
#'         list(gene_id = 2, delta_influence = matrix(rnorm(30, 0, 0.2), 5, 6)),
#'         list(gene_id = 3, delta_influence = matrix(rnorm(30, 0, 0.2), 5, 6)),
#'         list(gene_id = 4, delta_influence = matrix(rnorm(30, 0, 0.2), 5, 6))
#'       )
#'     )
#'   ),
#'   class = "tsenat_isoform_switching_multiq"
#' )
#' 
#' # Create heatmap visualization
#' plot_multiq_delta_influence_heatmaps(switching_results, n_genes = 2)
#'
#' @import grid
#' @import pheatmap
#' @importFrom grDevices png dev.off colorRampPalette
#' @keywords internal
#' @noRd
plot_multiq_delta_influence_heatmaps <- function(
    switching_results,
  n_genes = 4,
  lm_results = NULL,
  verbose = FALSE,
  cellwidth = 0,
  cellheight = 0,
  fontsize = 18,
  layout_ncol = 2,
  output_file = NULL) {
  # cellwidth, cellheight, fontsize follow pheatmap best practices:
  # - fontsize=18pt default for readable, large-format heatmaps (GLOBAL constant from .tsenat_font_sizes$heatmap_main)
  # - cellwidth=0, cellheight=0 (default) trigger dynamic sizing based on layout and data dimensions
  # - Dynamic sizing is aggressive: prioritizes visibility over whitespace
  # - Set cellwidth > 0 and cellheight > 0 to use fixed cell sizes instead and override dynamic sizing
  # - These are applied per individual heatmap in the grid layout
  # - layout_ncol: Fixed number of heatmaps per row (default: 2); set to NULL for adaptive layout based on transcript counts
  # Input validation
  if (!inherits(switching_results, "tsenat_isoform_switching_multiq")) {
    stop("switching_results must be a multi-q result from jackknife_isoform_switching()")
  }
  
  # Extract q-values from result names (q_0_01, q_0_50, etc.)
  q_result_keys <- names(switching_results)[grepl("^q_", names(switching_results))]
  if (length(q_result_keys) == 0) {
    stop("No multi-q results found in switching_results")
  }
  
  # Get the first result to extract gene information (all q values analyze same genes)
  first_q_key <- q_result_keys[1]
  first_result <- switching_results[[first_q_key]]
  
  # Extract gene IDs and names from first q-value result
  gene_ids <- first_result$gene_ids
  gene_name_map <- first_result$gene_name_map
  
  if (length(gene_ids) == 0) {
    stop("No genes found in switching_results")
  }
  
  # Select top N genes
  if (!is.null(lm_results)) {
    # Find p-value or adjusted p-value column for ranking
    p_col <- if ("adj_p_interaction" %in% colnames(lm_results)) {
      "adj_p_interaction"
    } else if ("p_interaction" %in% colnames(lm_results)) {
      "p_interaction"
    } else {
      NULL
    }
    
    # If we have a p-value column, rank genes by statistical significance
    if (!is.null(p_col)) {
      # Find which column in lm_results contains gene identifiers
      lm_col <- if ("gene_id" %in% colnames(lm_results)) {
        "gene_id"
      } else if ("gene_name" %in% colnames(lm_results)) {
        "gene_name"
      } else {
        NULL
      }
      
      # Match genes and sort by p-value
      if (!is.null(lm_col)) {
        matches <- match(gene_ids, lm_results[[lm_col]])
        p_values <- rep(Inf, length(gene_ids))
        matched_idx <- !is.na(matches)
        p_values[matched_idx] <- lm_results[[p_col]][matches[matched_idx]]
        
        # Sort by p-value (lowest p-values = most significant)
        gene_order <- order(p_values)
        gene_ids_sorted <- gene_ids[gene_order]
        top_genes_for_comparison <- gene_ids_sorted[seq_len(min(n_genes, length(gene_ids_sorted)))]
      } else {
        # Could not find gene identifier column, use first N genes
        top_genes_for_comparison <- gene_ids[seq_len(min(n_genes, length(gene_ids)))]
      }
    } else {
      # No p-value column found, use first N genes
      top_genes_for_comparison <- gene_ids[seq_len(min(n_genes, length(gene_ids)))]
    }
  } else {
    # No lm_results provided, use first N genes
    top_genes_for_comparison <- gene_ids[seq_len(min(n_genes, length(gene_ids)))]
  }
  
  # Prepare data for combined multi-Q heatmap
  all_gene_matrices <- list()
  all_gene_info <- list()
  data_validity_report <- list()  # Track validation for debugging
  
  for (gene_idx in seq_along(top_genes_for_comparison)) {
    gene_id <- top_genes_for_comparison[gene_idx]
    
    # Look up gene name from gene_name_map, use gene_id as fallback if NA
    gene_name_idx <- which(gene_ids == gene_id)[1]
    gene_name <- gene_id  # Default to gene_id
    if (!is.na(gene_name_idx) && !is.na(gene_name_map[gene_name_idx])) {
      gene_name <- gene_name_map[gene_name_idx]
    }
    
    # Initialize validation report for this gene
    validity_report <- list(
      gene_id = gene_id,
      gene_name = gene_name,
      has_heatmap_data = FALSE,
      has_valid_rows = FALSE,
      has_valid_cols = FALSE,
      has_valid_transcripts = FALSE,
      reason_skipped = NA_character_
    )
    
    # Collect delta_influence for all transcripts across all q-values
    heatmap_data <- NULL
    
    for (q_key in q_result_keys) {
      
      if (!is.null(switching_results[[q_key]]) && 
          !is.null(switching_results[[q_key]]$results_per_gene) &&
          gene_id %in% names(switching_results[[q_key]]$results_per_gene)) {
        gene_res <- switching_results[[q_key]]$results_per_gene[[gene_id]]
        if (!is.null(gene_res$delta_influence)) {
          # Extract q value from key: q_0_01 -> remove "q_" -> "0_01" -> replace "_" with "." -> "0.01"
          q_str_cleaned <- gsub("_", ".", gsub("^q_", "", q_key))
          q_num <- as.numeric(q_str_cleaned)
          col_name <- paste0("q_", sprintf("%.2f", q_num))
          delta_vals <- as.numeric(gene_res$delta_influence)
          
          # Replace Inf and NaN with NA for clean handling
          delta_vals[!is.finite(delta_vals)] <- NA
          
          if (is.null(heatmap_data)) {
            # Initialize with transcript IDs
            heatmap_data <- data.frame(
              transcript = as.character(gene_res$transcript_ids),
              stringsAsFactors = FALSE
            )
          }
          
          # Add column for this q-value (with Inf/NaN as NA)
          # Match number of rows and ensure alignment
          n_rows <- nrow(heatmap_data)
          if (length(delta_vals) < n_rows) {
            delta_vals <- c(delta_vals, rep(NA_real_, n_rows - length(delta_vals)))
          } else if (length(delta_vals) > n_rows) {
            delta_vals <- delta_vals[seq_len(n_rows)]
          }
          heatmap_data[[col_name]] <- as.numeric(delta_vals)
        }
      }
    }
    
    if (!is.null(heatmap_data) && nrow(heatmap_data) > 0 && ncol(heatmap_data) > 1) {
      validity_report$has_heatmap_data <- TRUE
      validity_report$has_valid_rows <- (nrow(heatmap_data) > 0)
      validity_report$has_valid_cols <- (ncol(heatmap_data) > 1)
      
      # Convert to matrix for heatmap (transcripts as rows, q-values as columns)
      heatmap_matrix <- as.matrix(heatmap_data[, -1, drop = FALSE])
      rownames(heatmap_matrix) <- heatmap_data$transcript
      # Explicitly preserve column names (q-value labels)
      colnames(heatmap_matrix) <- colnames(heatmap_data)[-1]
      
      # Remove transcripts with exclusively no valid data (all NA)
      valid_transcript_rows <- rowSums(!is.na(heatmap_matrix)) > 0
      if (sum(valid_transcript_rows) > 0) {
        heatmap_matrix <- heatmap_matrix[valid_transcript_rows, , drop = FALSE]
        
        # Skip if no valid transcripts remain or fewer than 1 q-value column
        if (nrow(heatmap_matrix) > 0 && ncol(heatmap_matrix) > 0) {
          validity_report$has_valid_transcripts <- TRUE
          
          # Cap outliers for remaining finite values
          finite_vals <- heatmap_matrix[is.finite(heatmap_matrix)]
          if (length(finite_vals) > 0) {
            # Convert to numeric to avoid type coercion issues
            fin_abs <- abs(as.numeric(finite_vals))
            fin_abs <- fin_abs[is.finite(fin_abs)]
            if (length(fin_abs) > 0) {
              cap_val <- as.numeric(quantile(fin_abs, probs = 0.95))
              # Cap finite values that exceed the 95th percentile
              abs_hm <- abs(heatmap_matrix)
              mask <- which(is.finite(abs_hm) & abs_hm > cap_val)
              if (length(mask) > 0) {
                heatmap_matrix[mask] <- sign(heatmap_matrix[mask]) * cap_val
              }
            }
          }
          
          # Transpose: q-values as rows, transcripts as columns
          # This converts from (transcripts x q-values) to (q-values x transcripts)
          heatmap_matrix <- t(heatmap_matrix)
          
          # Verify we have the expected dimensions (q-values as rows)
          if (nrow(heatmap_matrix) > 0 && ncol(heatmap_matrix) > 0) {
            # Store matrix and gene info for combined plot
            all_gene_matrices[[gene_idx]] <- heatmap_matrix
            all_gene_info[[gene_idx]] <- list(
              gene_id = gene_id,
              gene_name = gene_name,
              n_transcripts = ncol(heatmap_matrix)
            )
          } else {
            validity_report$reason_skipped <- "Matrix dimensions invalid after transpose"
          }
        } else {
          validity_report$reason_skipped <- paste0("No valid transcripts (nrow=", nrow(heatmap_matrix), ", ncol=", ncol(heatmap_matrix), ")")
        }
      } else {
        validity_report$reason_skipped <- "All transcript rows are all-NA"
      }
    } else {
      if (is.null(heatmap_data)) {
        validity_report$reason_skipped <- "No heatmap_data collected (no delta_influence found)"
      } else {
        validity_report$reason_skipped <- paste0("Insufficient data (nrow=", nrow(heatmap_data), ", ncol=", ncol(heatmap_data), ")")
      }
    }
    
    data_validity_report[[gene_idx]] <- validity_report
  }
  
  # Note: Rendering directly to active graphics device (like DESeq2)
  # This allows knitr to capture the output during vignette compilation
  
  if (length(all_gene_matrices) == 0) {
    warning("No valid heatmap data generated for any genes")
    # Suppress validation report output (only shown if verbose=TRUE)
    if (verbose) {
      message("\n=== Data Validation Report ===")
      for (i in seq_along(data_validity_report)) {
        report <- data_validity_report[[i]]
        message("\nGene #", i, ": ", report$gene_name, " (", report$gene_id, ")")
        message("  - Has heatmap data: ", report$has_heatmap_data)
        message("  - Has valid rows: ", report$has_valid_rows)
        message("  - Has valid cols: ", report$has_valid_cols)
        message("  - Has valid transcripts: ", report$has_valid_transcripts)
        if (!is.na(report$reason_skipped)) {
          message("  - Reason skipped: ", report$reason_skipped)
        }
      }
      message("\n================================\n")
    }
    return(invisible(NULL))
  }
  
  # Show validation report if verbose=TRUE (shows all genes, even those skipped)
  if (verbose) {
    message("\n=== Data Validation Report ===")
    for (i in seq_along(data_validity_report)) {
      report <- data_validity_report[[i]]
      status <- if (!is.na(report$reason_skipped)) "? SKIPPED" else "[OK] VALID"
      message("\nGene #", i, ": ", report$gene_name, " (", report$gene_id, ") - ", status)
      message("  - Has heatmap data: ", report$has_heatmap_data)
      message("  - Has valid rows: ", report$has_valid_rows)
      message("  - Has valid cols: ", report$has_valid_cols)
      message("  - Has valid transcripts: ", report$has_valid_transcripts)
      if (!is.na(report$reason_skipped)) {
        message("  - Reason: ", report$reason_skipped)
      }
    }
    message("\n================================\n")
  }
  
  if (!requireNamespace("pheatmap", quietly = TRUE)) {
    stop("pheatmap package required for this function. Install with: install.packages('pheatmap')")
  }
  
  tryCatch({
    # First, determine layout based on heatmap transcript counts
    # This is done BEFORE creating heatmaps so cellsizes can be adjusted
    n_total_genes <- length(top_genes_for_comparison)
    
    # Determine which genes have >5 transcripts (will use full-width rows)
    gene_layout <- list()  # Will store layout info for each gene
    n_layout_rows <- 0     # Count of actual rows needed
    
    if (length(all_gene_matrices) > 0) {
      # Determine layout strategy (fixed columns or adaptive)
      use_fixed_layout <- !is.null(layout_ncol) && layout_ncol > 0
      
      if (use_fixed_layout) {
        # FIXED LAYOUT: Force layout_ncol heatmaps per row
        n_cols <- as.integer(layout_ncol)
        i <- 1
        while (i <= n_total_genes) {
          for (col_pos in seq_len(n_cols)) {
            if (i <= n_total_genes) {
              gene_layout[[i]] <- list(row = n_layout_rows + 1, col = col_pos, width = 1/n_cols)
              i <- i + 1
            }
          }
          n_layout_rows <- n_layout_rows + 1
        }
      } else {
        # ADAPTIVE LAYOUT: Original logic based on transcript counts
        i <- 1
        while (i <= n_total_genes) {
          gene_idx <- i
          has_data_i <- !is.null(all_gene_matrices[[gene_idx]]) && nrow(all_gene_matrices[[gene_idx]]) > 0
          n_transcripts_i <- if (has_data_i) ncol(all_gene_matrices[[gene_idx]]) else 0
          
          # Check if next gene exists and has data
          has_next <- i < n_total_genes
          has_data_next <- has_next && !is.null(all_gene_matrices[[i + 1]]) && nrow(all_gene_matrices[[i + 1]]) > 0
          n_transcripts_next <- if (has_data_next) ncol(all_gene_matrices[[i + 1]]) else 0
          
          if (has_data_i && n_transcripts_i > 5) {
            # Gene with >5 transcripts: full-width row
            gene_layout[[gene_idx]] <- list(row = n_layout_rows + 1, col = 1, width = 1)
            n_layout_rows <- n_layout_rows + 1
            i <- i + 1
          } else if (has_data_i && n_transcripts_i <= 5 && has_data_next && n_transcripts_next <= 5) {
            # Two consecutive genes BOTH with <=5 transcripts: pair them
            gene_layout[[gene_idx]] <- list(row = n_layout_rows + 1, col = 1, width = 0.5)
            gene_layout[[i + 1]] <- list(row = n_layout_rows + 1, col = 2, width = 0.5)
            n_layout_rows <- n_layout_rows + 1
            i <- i + 2
          } else {
            # Single gene or last gene: full row
            gene_layout[[gene_idx]] <- list(row = n_layout_rows + 1, col = 1, width = 1)
            n_layout_rows <- n_layout_rows + 1
            i <- i + 1
          }
        }
      }
    } else {
      # No data: assign all genes to layouts (as placeholders)
      i <- 1
      while (i <= n_total_genes) {
        if (i < n_total_genes) {
          # Pair genes if possible
          gene_layout[[i]] <- list(row = n_layout_rows + 1, col = 1, width = 0.5)
          gene_layout[[i + 1]] <- list(row = n_layout_rows + 1, col = 2, width = 0.5)
          n_layout_rows <- n_layout_rows + 1
          i <- i + 2
        } else {
          # Last unpaired gene
          gene_layout[[i]] <- list(row = n_layout_rows + 1, col = 1, width = 1)
          n_layout_rows <- n_layout_rows + 1
          i <- i + 1
        }
      }
    }
    
    # Calculate heatmap dimensions BEFORE processing individual heatmaps
    # These values are needed for adaptive cell sizing calculations
    # Height scales proportionally with number of q-values (rows per heatmap)
    n_q_values <- length(q_result_keys)  # Number of q-values (rows per heatmap)
    # Base: 3 inches per layout row for ~5 q-values
    # Scale linearly: more q-values = taller heatmaps = more space needed
    height_per_layout_row <- 3 * (n_q_values / 5)  # Scales from 3" for 5 q-values
    gap_between_rows <- 1.5
    heatmap_height <- height_per_layout_row * n_layout_rows + gap_between_rows * (n_layout_rows - 1)
    
    # Create individual heatmaps for each gene and store as grobs
    heatmap_plots <- list()
    plot_gene_names <- character(0)
    genes_with_data <- integer(0)
    
    # Create plots for genes with data, track which genes have valid data
    for (gene_idx in seq_along(top_genes_for_comparison)) {
      if (is.null(all_gene_matrices[[gene_idx]]) || is.null(all_gene_info[[gene_idx]])) {
        # Mark this position with empty placeholder
        heatmap_plots[[gene_idx]] <- NULL
        next
      }
      
      genes_with_data <- c(genes_with_data, gene_idx)
      mat <- all_gene_matrices[[gene_idx]]
      gene_info <- all_gene_info[[gene_idx]]
      gene_name <- gene_info$gene_name
      gene_id <- gene_info$gene_id
      plot_gene_names <- c(plot_gene_names, gene_name)
      
      # Construct header text for this gene
      if (is.na(gene_name) || gene_name == "") {
        header_text <- gene_id
      } else {
        header_text <- gene_name
      }
      
      # Apply outlier capping (safely handle if all values are NA)
      finite_vals <- as.numeric(mat[is.finite(mat)])
      if (length(finite_vals) == 0) {
        percentile_95 <- 1  # Default if no finite values
      } else {
        percentile_95 <- as.numeric(quantile(abs(finite_vals), 0.95))
      }
      
      mat_viz <- mat
      if (length(finite_vals) > 0) {
        # Cap outliers using indices to avoid coercion warnings
        abs_mat <- abs(mat_viz)
        mask_idx <- which(is.finite(abs_mat) & abs_mat > percentile_95)
        if (length(mask_idx) > 0) {
          mat_viz[mask_idx] <- sign(mat_viz[mask_idx]) * percentile_95
        }
      }
      
      # Create pheatmap (returns a grob object)
      # Using best practices from pheatmap documentation:
      # - fontsize=13pt (standard for publication heatmaps)
      # - cellwidth/cellheight dynamically scaled based on matrix dimensions
      # - Color scale: diverging palette (blue-white-red) centered at zero
      
      # Calculate dynamic cell sizes based on AVAILABLE GRID SPACE
      # Simplified adaptive sizing: scale based on number of columns/rows
      n_cols_mat <- ncol(mat_viz)
      n_rows_mat <- nrow(mat_viz)
      
      # Get layout info for this gene
      layout_info <- gene_layout[[gene_idx]]
      heatmap_width_fraction <- if (!is.null(layout_info)) layout_info$width else 1
      
      # BASE cell sizes: SCALED FOR 1200px WIDTH (12 inches @ 100 DPI)
      # Reduced by 30% for proportional sizing
      base_cellwidth <- 35   # 50 * 0.7 for 30% reduction
      base_cellheight <- 29  # 42 * 0.7 for 30% reduction
      
      # For half-width heatmaps (2-per-row), calculate cellwidth to ensure ~32% plot width per heatmap (reduced by 20%)
      # Plot width: 1200px × 0.96 (viewport) × 0.4 (reduced from 0.5) = ~461px per heatmap
      # Minus margins (~40px left labels) and borders (~30px) = ~405px for heatmap cells
      # Divide by number of transcripts to get cellwidth
      if (heatmap_width_fraction < 1) {
        # Half-width: allocate ~32% of plot width per heatmap (20% reduction)
        # 1200px × 0.96 × 0.4 = 461px available
        # Minus margins (~40px left labels) and borders (~30px) = ~405px for cells
        available_width_px <- 1200 * 0.65 * heatmap_width_fraction - 40 - 30
        cellwidth_calc <- available_width_px / n_cols_mat
        # Allow cellwidth to vary based on transcript count (no artificial caps)
        adaptive_cellwidth <- max(15, cellwidth_calc)  # Minimum 15px to stay readable
      } else {
        # Full-width: use base scaling with expansion
        if (n_cols_mat > 8) {
          scale_factor_width <- 2.3
        } else {
          scale_factor_width <- 2.6
        }
        adaptive_cellwidth <- base_cellwidth * scale_factor_width
      }
      
      # Scale height based on number of rows (q-values) using linear formula
      # Continuously reduces cell height as rows increase
      scale_factor_height <- max(0.7, 1.15 - n_rows_mat * 0.03)
      
      adaptive_cellheight <- base_cellheight * scale_factor_height
      
      # If explicit cellwidth/cellheight provided and >0, use those; else use adaptive
      final_cellwidth <- if (cellwidth > 0) cellwidth else adaptive_cellwidth
      final_cellheight <- if (cellheight > 0) cellheight else adaptive_cellheight
      
    p <- pheatmap::pheatmap(
        mat_viz,
        main = header_text,
        cluster_rows = FALSE,
        cluster_cols = (ncol(mat_viz) > 1),
        display_numbers = FALSE,
        na_col = "lightgray",
        border_color = "black",           # BLACK borders for better color separation
        color = grDevices::colorRampPalette(c("#4575B4", "#FFFFFF", "#D73027"))(70),
        cellwidth = final_cellwidth,
        cellheight = final_cellheight,
        fontsize = fontsize * 0.7,        # Reduced by 30%
        fontsize_row = fontsize * 0.7,
        fontsize_col = fontsize * 0.7,
        fontsize_number = fontsize * 0.56, # Reduced by 30%
        margins = c(8, 10),               # REDUCED margins to maximize heatmap space
        show_rownames = TRUE,
        show_colnames = TRUE,
        silent = TRUE
      )
      
      heatmap_plots[[gene_idx]] <- p
    }
    
    # Pad with placeholder grobs for genes without data to maintain grid structure
    n_total_genes <- length(top_genes_for_comparison)
    for (gene_idx in seq_len(n_total_genes)) {
      if (is.null(heatmap_plots[[gene_idx]])) {
        # Create an empty placeholder grob
        placeholder_grob <- grid::gTree(
          children = grid::gList(
            grid::rectGrob(gp = grid::gpar(fill = "white", col = "lightgray", lwd = 2)),
            grid::textGrob("No data available", x = 0.5, y = 0.5, 
                          gp = grid::gpar(col = "gray50", fontsize = 18))
          )
        )
        heatmap_plots[[gene_idx]] <- placeholder_grob
      }
    }
    
    # Combine all panels into one figure using manual grid layout
    # Use mixed layout: full-width rows for large heatmaps, 2-per-row for small ones
    
    n_genes <- n_total_genes
    
    # Conditional rendering: PNG file if output_file provided, otherwise active device
    if (!is.null(output_file)) {
      # Render to PNG file when output_file is specified
      # PNG dimensions: standardized to 1200px width (12 inches @ 100 DPI) for consistency with other plots
      png_width <- 12   # 12 inches @ 100 DPI = 1200 pixels
      png_dpi <- 100
      
      # Scale heatmap_height proportionally: was 6 inches per row @ 150 DPI, now at 100 DPI
      # Adjust from (heatmap_height @ 150 DPI context) to (heatmap_height @ 100 DPI context)
      heatmap_height_scaled <- heatmap_height * (png_dpi / 150)
      
      grDevices::png(output_file, width = png_width, height = heatmap_height_scaled + 4, 
                     units = "in", res = png_dpi)
    }
    
    # Render directly to active graphics device (managed by knitr during vignette compilation)
    # or to PNG file if output_file was specified (opened above)
    grid::grid.newpage()
    
    # Calculate dynamic title/subtitle sizes based on number of rows
    # Scales proportionally: 15% increase per additional row layout
    # For typical 2-row layout: ~30% bigger than base (16pt * 1.30 = 20.8pt ≈ 21pt)
    title_fontsize <- 16 * (1 + 0.15 * n_layout_rows)
    subtitle_fontsize <- 12 * (1 + 0.15 * n_layout_rows)
    
    # Add main title 
    grid::grid.text("Delta Influence Across Diversity Scales", 
                    x = 0.5, y = 0.97, 
                    just = "top",
                    gp = grid::gpar(fontsize = title_fontsize, fontface = "bold"))
    
    # Add subtitle
    grid::grid.text("Jackknife weights across q-spectrum for selected genes", 
                    x = 0.5, y = 0.94, 
                    just = "top",
                    gp = grid::gpar(fontsize = subtitle_fontsize, fontface = "italic", col = "gray40"))
    
    # Create viewport layout with variable columns per row
    # Each actual content row is followed by a gap row
    n_grid_rows <- n_layout_rows * 2 - 1
    row_heights <- rep(c(1, 0.15), n_layout_rows)[seq_len(n_grid_rows)]
    
    # Use 3 columns as grid basis: column 1 (left), column 2 (gap), column 3 (right)
    # This allows proper spacing between heatmaps in 2-column layouts
    # MAXIMUM width (0.96) and height (0.85) to fill plot space with reduced title/subtitle spacing
    grid::pushViewport(grid::viewport(
      x = 0.5, 
      y = 0.48, 
      width = 0.96,
      height = 0.85,
      layout = grid::grid.layout(
        n_grid_rows, 
        3,
        heights = grid::unit(row_heights, "null"),
        widths = c(1, 0.12, 1),  # Col 1: left (1), gap (0.12), col 3: right (1)
        respect = FALSE
      )
    ))
    
    # Draw heatmaps using gene_layout positions
    for (plot_idx in seq_along(heatmap_plots)) {
      layout_info <- gene_layout[[plot_idx]]
      if (!is.null(layout_info)) {
        grid_row <- layout_info$row * 2 - 1  # Convert to grid row (accounting for gaps)
        grid_col <- layout_info$col
        width_frac <- layout_info$width
        
        if (width_frac == 1) {
          # Full width: span all columns (1, 2, 3)
          grid_col_start <- 1
          grid_col_end <- 3
        } else if (layout_info$col == 1) {
          # Left column: only column 1
          grid_col_start <- 1
          grid_col_end <- 1
        } else {
          # Right column: only column 3 (skip gap)
          grid_col_start <- 3
          grid_col_end <- 3
        }
        
        grid::pushViewport(grid::viewport(
          layout.pos.row = grid_row,
          layout.pos.col = grid_col_start:grid_col_end
        ))
        grid::grid.draw(heatmap_plots[[plot_idx]])
        grid::popViewport()
      }
    }
    
    grid::popViewport()
    
    # Close PNG device if it was opened
    if (!is.null(output_file)) {
      grDevices::dev.off()
      if (verbose) {
        message("Heatmap saved to: ", output_file)
      }
      return(output_file)
    } else {
      # Return invisible(NULL) - graphics are captured by knitr during vignette compilation
      return(invisible(NULL))
    }
    
  }, error = function(e) {
    if (!is.null(output_file)) {
      tryCatch(grDevices::dev.off(), silent = TRUE)
    }
    stop("Heatmap creation failed: ", e$message)
  })
}


