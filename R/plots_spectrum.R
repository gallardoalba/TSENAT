#' Plot Global Divergence q-Curve Across All Genes
#'
#' Visualizes the average (mean/median) Tsallis divergence D_q across all genes
#' as a function of q-value. This provides a **global view** of which
#' diversity scales
#' (rare vs. abundant isoforms) drive the most divergence on average across
#' the dataset,
#' complementing gene-specific divergence profiles.
#'
#' @param divergence_results_se A `SummarizedExperiment` containing
#' pre-computed divergence values.
#' Rows = genes, columns = q-values. Column names should indicate q-values
#' (e.g., 'q_0.5', 'q_1.0').
#' @param gene Optional character. If provided, plot divergence spectrum for
#' this specific gene.
#'   If NULL, plot global divergence curve (aggregated across all genes).
#' @param lm_res Optional data.frame with columns for gene identifiers and
#' p-values. Used to select
#'   top genes when gene = NULL and lm_res is provided. Default: NULL.
#' @param n_genes Integer; number of top genes to plot when showing
#' multi-gene spectra (default: 4).
#' Must be positive. Genes are sorted by p-value significance (lowest
#' p-values first).
#' @param ncol Integer; number of columns in grid layout for multi-gene
#' plots (default: 2).
#' Must be positive. Number of rows is automatically calculated as
#' ceiling(n_genes / ncol).
#' @param metric Character. Summary statistic for global curve: 'median' or
#' 'mean'. Default: 'median'.
#'   Only used when gene = NULL.
#' @param variability_metric Character. Error bar type for global curve:
#' 'sd' or 'iqr'. Default: 'iqr'.
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
#' - Shows which diversity scales (q-values) drive the most divergence on
#' average
#' - Useful for identifying dominant biological mechanisms (rare vs.
#' abundant isoform driven)
#'
#' **Interpretation**: Compare with `plot_tsallis_q_curve` (entropy) to
#' understand
#' the relationship between entropy changes and divergence patterns.
#'
#' @importFrom ggplot2 ggplot aes geom_line geom_point geom_ribbon labs
#' theme_minimal element_text
#' @importFrom SummarizedExperiment assay rowData
#' @importFrom stats IQR var
#'
#' @examples
#' # Create synthetic divergence data
#' set.seed(123)
#' divergence_matrix <- matrix(
#'   rnorm(80, mean = 0.5, sd = 0.1),
#'   nrow = 20, ncol = 4
#' )
#' rownames(divergence_matrix) <- paste0('gene_', 1:20)
#' colnames(divergence_matrix) <- c('q_0.5', 'q_1.0', 'q_1.5', 'q_2.0')
#' divergence_se <- SummarizedExperiment::SummarizedExperiment(
#'   assays = list(divergence = divergence_matrix)
#' )
#' 
#' # Global divergence curve (all genes aggregated)
#' p_global <- .plot_divergence_spectrum(divergence_se)
#' 
#' # Gene-specific divergence spectrum
#' p_gene <- .plot_divergence_spectrum(divergence_se, gene = 'gene_1')
#'

#' @noRd
.spectrum_extract_q_values <- function(col_names) {
    extracted <- gsub("^q[_=]", "", col_names)
    q_vals <- as.numeric(extracted)
    if (any(is.na(q_vals))) {
        stop("Cannot extract numeric q-values from column names. ", "Expected format like 'q_0.5'. Got: ",
            paste(head(col_names, 3), collapse = ", "), call. = FALSE)
    }
    return(q_vals)
}

#' @noRd
.spectrum_get_gene_identifiers <- function(se, div_mat) {
    rd <- SummarizedExperiment::rowData(se)
    if (!is.null(rd) && "gene_name" %in% colnames(rd))
        return(rd$gene_name)
    rn <- rownames(div_mat)
    if (is.null(rn)) {
        stop("'divergence_results_se' has no gene identifiers in rowData or rownames",
            call. = FALSE)
    }
    return(rn)
}

#' @noRd
.spectrum_find_gene_column <- function(df) {
    valid_cols <- c("gene", "gene_name", "gene_id")
    found <- valid_cols[valid_cols %in% colnames(df)]
    if (length(found) == 0) {
        stop("'lm_res' must have a column named 'gene', 'gene_name', or 'gene_id'",
            call. = FALSE)
    }
    return(found[1])
}

#' @noRd
.spectrum_find_pvalue_column <- function(df) {
    valid_cols <- c("adj_p_interaction", "p_interaction", "adj_p_value", "p_value")
    found <- valid_cols[valid_cols %in% colnames(df)]
    if (length(found) == 0) {
        stop("'lm_res' must have a p-value column like adj_p_interaction or p_value",
            call. = FALSE)
    }
    return(found[1])
}

#' @noRd
.spectrum_validate_inputs <- function(divergence_results_se, gene, lm_res, n_genes,
    ncol) {
    if (!inherits(divergence_results_se, "SummarizedExperiment")) {
        stop("'divergence_results_se' must be a SummarizedExperiment", call. = FALSE)
    }
    div_mat <- tryCatch({
        SummarizedExperiment::assay(divergence_results_se, 1)
    }, error = function(e) {
        stop("Failed to extract assay: ", e$message, call. = FALSE)
    })
    if (is.null(div_mat) || nrow(div_mat) == 0 || ncol(div_mat) == 0) {
        stop("'divergence_results_se' assay is empty", call. = FALSE)
    }
    if (!is.null(gene) && (!is.character(gene) || length(gene) != 1)) {
        stop("'gene' must be a single character string or NULL", call. = FALSE)
    }
    if (!is.null(lm_res) && (!is.data.frame(lm_res) || nrow(lm_res) == 0)) {
        stop("'lm_res' must be a non-empty data.frame or NULL", call. = FALSE)
    }
    if (!is.numeric(n_genes) || n_genes < 1)
        stop("'n_genes' must be positive", call. = FALSE)
    if (!is.numeric(ncol) || ncol < 1)
        stop("'ncol' must be positive", call. = FALSE)
    return(div_mat)
}

#' @noRd
.spectrum_plot_single_gene <- function(gene_name, div_mat_sorted, q_vals_sorted,
    gene_names) {
    require_pkgs("ggplot2")
    
    gene_idx <- which(gene_names == gene_name)[1]
    if (is.na(gene_idx)) {
        stop("Gene '", gene_name, "' not found. Available: ", paste(head(gene_names,
            5), collapse = ", "), call. = FALSE)
    }
    
    plot_df <- data.frame(q = q_vals_sorted, divergence = as.numeric(div_mat_sorted[gene_idx, ]),
        stringsAsFactors = FALSE)
    
    # Use helper to create line + point plot (fixed color, no grouping)
    p <- .create_simple_line_plot(plot_df, x_col = "q", y_col = "divergence",
        group_col = NULL, points = TRUE, line_width = 1.2, point_size = 3.5,
        line_color = .palette_blue_red()[1])
    
    p <- .apply_publication_theme(p, base_theme = "theme_base", base_size = 11,
        title = paste("Divergence Spectrum:", gene_name))
    
    p <- p + ggplot2::labs(x = "q value", y = "Divergence D_q")
    
    p
}

#' @noRd
.spectrum_plot_top_genes <- function(lm_res, n_genes_use, ncol, div_mat_sorted, q_vals_sorted,
    gene_names, metric, divergence_results_se) {
    gene_col <- .spectrum_find_gene_column(lm_res)
    p_col <- .spectrum_find_pvalue_column(lm_res)
    lm_sorted <- lm_res[order(lm_res[[p_col]], na.last = TRUE), , drop = FALSE]
    top_genes_vec <- head(lm_sorted[[gene_col]], n_genes_use)

    if (length(top_genes_vec) == 0) {
        stop("No genes found in 'lm_res'", call. = FALSE)
    }

    gene_indices <- match(top_genes_vec, gene_names)
    unmatched <- is.na(gene_indices)
    if (all(unmatched)) {
        stop("None of the top genes found in divergence matrix", call. = FALSE)
    }
    if (any(unmatched)) {
        warning("Some genes not found. Proceeding with ", sum(!unmatched), " matches",
            call. = FALSE)
        gene_indices <- gene_indices[!unmatched]
        top_genes_vec <- top_genes_vec[!unmatched]
    }

    # Check if bootstrap CI assays are available
    has_ci_assays <- FALSE
    ci_lower_mat <- NULL
    ci_upper_mat <- NULL

    if (!is.null(divergence_results_se) && is(divergence_results_se, "SummarizedExperiment")) {
        assay_names <- names(SummarizedExperiment::assays(divergence_results_se))
        if ("ci_lower" %in% assay_names && "ci_upper" %in% assay_names) {
            has_ci_assays <- TRUE
            ci_lower_mat <- SummarizedExperiment::assay(divergence_results_se, "ci_lower")[,
                colnames(div_mat_sorted)]
            ci_upper_mat <- SummarizedExperiment::assay(divergence_results_se, "ci_upper")[,
                colnames(div_mat_sorted)]
        }
    }

    # Build data frame with bootstrap CIs if available
    plot_list <- lapply(seq_along(gene_indices), function(i) {
        df <- data.frame(q = q_vals_sorted, divergence = as.numeric(div_mat_sorted[gene_indices[i],
            ]), gene = gene_names[gene_indices[i]], p_value = lm_sorted[[p_col]][i],
            stringsAsFactors = FALSE)

        # Add bootstrap CI columns if available
        if (has_ci_assays && !all(is.na(ci_lower_mat)) && !all(is.na(ci_upper_mat))) {
            df$ci_lower <- as.numeric(ci_lower_mat[gene_indices[i], ])
            df$ci_upper <- as.numeric(ci_upper_mat[gene_indices[i], ])
        }

        df
    })

    multi_gene_df <- do.call(rbind, plot_list)
    rownames(multi_gene_df) <- NULL
    gene_order <- multi_gene_df[!duplicated(multi_gene_df$gene), ]
    gene_order <- gene_order[order(gene_order$p_value), ]$gene
    multi_gene_df$gene <- factor(multi_gene_df$gene, levels = gene_order)

    # Build plot with optional bootstrap CI ribbon
    p <- ggplot2::ggplot(multi_gene_df, ggplot2::aes(x = q, y = divergence))

    # Add confidence ribbon if bootstrap CIs are available and valid
    if (has_ci_assays && "ci_lower" %in% colnames(multi_gene_df) && !all(is.na(multi_gene_df$ci_lower)) &&
        !all(is.na(multi_gene_df$ci_upper))) {
        p <- p + ggplot2::geom_ribbon(ggplot2::aes(ymin = ci_lower, ymax = ci_upper),
            alpha = 0.15, fill = "#4575B4", color = NA)
    }

    p <- .apply_facet_styling(p, ncol = ncol, facet_var = "gene", scales = "free_y") +
        ggplot2::geom_line(color = "#4575B4", linewidth = 1.2, alpha = 0.8) +
        ggplot2::geom_point(color = "#4575B4", size = 3, alpha = 0.8)
    p <- .apply_publication_theme(p, base_theme = "theme_base", base_size = 11,
        title = "Divergence Spectra: Per-gene Comparisons",
        subtitle = if (has_ci_assays && "ci_lower" %in% colnames(multi_gene_df)) {
            paste0("Ranked by interaction significance (", metric, ") | Bootstrap CI (95%)")
        } else {
            paste0("Ranked by interaction significance (", metric, ")")
        })
    p <- p + ggplot2::labs(x = "q value", y = expression("Divergence D[q]"))

    return(p)
}

#' @noRd
.spectrum_plot_global <- function(div_mat_sorted, q_vals_sorted, metric, variability_metric,
    divergence_results_se = NULL) {
    require_pkgs(c("ggplot2", "SummarizedExperiment"))
    
    # Check if bootstrap CI assays are available
    has_ci_assays <- FALSE
    ci_lower_mat <- NULL
    ci_upper_mat <- NULL

    if (!is.null(divergence_results_se) && methods::is(divergence_results_se, "SummarizedExperiment")) {
        assay_names <- names(SummarizedExperiment::assays(divergence_results_se))
        if ("ci_lower" %in% assay_names && "ci_upper" %in% assay_names) {
            has_ci_assays <- TRUE
            ci_lower_mat <- SummarizedExperiment::assay(divergence_results_se, "ci_lower")[,
                colnames(div_mat_sorted)]
            ci_upper_mat <- SummarizedExperiment::assay(divergence_results_se, "ci_upper")[,
                colnames(div_mat_sorted)]
        }
    }

    # Unified subtitle construction function to avoid duplication
    construct_subtitle <- function(metric_label, ci_source, spread_label = NULL, n_genes) {
        if (!is.null(ci_source)) {
            paste0(metric_label, " with ", ci_source, " CI (", n_genes, " genes)")
        } else {
            paste0(metric_label, " +/- ", spread_label, " (", n_genes, " genes)")
        }
    }

    # Bootstrap CI branch
    if (has_ci_assays && !all(is.na(ci_lower_mat)) && !all(is.na(ci_upper_mat))) {
        summary_stats <- data.frame(q = q_vals_sorted, central = colMeans(div_mat_sorted,
            na.rm = TRUE), ci_lower = colMeans(ci_lower_mat, na.rm = TRUE), ci_upper = colMeans(ci_upper_mat,
            na.rm = TRUE), stringsAsFactors = FALSE)

        # Use CI ribbon helper with publication theme
        p <- .create_ci_ribbon_plot(summary_stats, x_col = "q", y_col = "central",
            ci_lower_col = "ci_lower", ci_upper_col = "ci_upper",
            ribbon_alpha = 0.1, line_width = 1.3, show_points = TRUE)
        p <- .apply_publication_theme(p, base_theme = "theme_base", base_size = 11,
            title = expression("Global Divergence Spectrum: Average " * D[q]),
            subtitle = construct_subtitle("Mean", "Bootstrap (95%)", n_genes = nrow(div_mat_sorted)))
        p <- p + ggplot2::labs(x = "q value", y = expression("Divergence D[q]"))
        return(p)
    }

    # Fallback to IQR/SD computation
    if (variability_metric == "iqr") {
        summary_stats <- data.frame(q = q_vals_sorted, 
            central = apply(div_mat_sorted, 2, function(x) {
                if (metric == "median") median(x, na.rm = TRUE) else mean(x, na.rm = TRUE)
            }), 
            spread = apply(div_mat_sorted, 2, function(x) stats::IQR(x, na.rm = TRUE)),
            stringsAsFactors = FALSE)
        spread_factor <- 0.5
        spread_label <- "IQR"
    } else {
        summary_stats <- data.frame(q = q_vals_sorted, 
            central = apply(div_mat_sorted, 2, function(x) {
                if (metric == "median") median(x, na.rm = TRUE) else mean(x, na.rm = TRUE)
            }), 
            spread = apply(div_mat_sorted, 2, function(x) sqrt(stats::var(x, na.rm = TRUE))),
            stringsAsFactors = FALSE)
        spread_factor <- 1
        spread_label <- "SD"
    }

    if (all(is.na(summary_stats$spread)) || all(is.na(summary_stats$central))) {
        stop("Cannot compute statistics. Check divergence matrix values", call. = FALSE)
    }

    # Create lower/upper bounds for ribbon plot
    summary_stats$ci_lower <- summary_stats$central - summary_stats$spread * spread_factor
    summary_stats$ci_upper <- summary_stats$central + summary_stats$spread * spread_factor

    # Use unified CI ribbon approach
    p <- .create_ci_ribbon_plot(summary_stats, x_col = "q", y_col = "central",
        ci_lower_col = "ci_lower", ci_upper_col = "ci_upper",
        ribbon_alpha = 0.1, line_width = 1.3, show_points = TRUE)
    
    metric_label <- if (metric == "median") "Median" else "Mean"
    p <- .apply_publication_theme(p, base_theme = "theme_base", base_size = 11,
        title = expression("Global Divergence Spectrum: Average " * D[q]),
        subtitle = construct_subtitle(metric_label, NULL, spread_label, nrow(div_mat_sorted)))
    p <- p + ggplot2::labs(x = "q value", y = expression("Divergence D[q]"))
    
    return(p)
}

# PLOT DISPATCH WRAPPER - Route to single gene, top genes, or global plot
#' @noRd
.plot_divergence_spectrum <- function(divergence_results_se, gene = NULL, lm_res = NULL,
    n_genes = 4, ncol = 2, metric = c("median", "mean"), variability_metric = c("iqr",
        "sd")) {
    # Validate inputs and extract matrix
    metric <- match.arg(metric)
    variability_metric <- match.arg(variability_metric)
    n_genes <- as.integer(n_genes)
    ncol <- as.integer(ncol)

    div_mat <- .spectrum_validate_inputs(divergence_results_se, gene, lm_res, n_genes,
        ncol)

    # Prepare data
    q_vals <- .spectrum_extract_q_values(colnames(div_mat))
    gene_names <- .spectrum_get_gene_identifiers(divergence_results_se, div_mat)

    if (length(gene_names) != nrow(div_mat)) {
        stop("Gene identifier count doesn't match matrix rows", call. = FALSE)
    }

    sort_idx <- order(q_vals)
    q_vals_sorted <- q_vals[sort_idx]
    div_mat_sorted <- div_mat[, sort_idx]

    # Dispatch to appropriate case
    if (!is.null(gene)) {
        return(.spectrum_plot_single_gene(gene, div_mat_sorted, q_vals_sorted, gene_names))
    }

    if (!is.null(lm_res)) {
        return(.spectrum_plot_top_genes(lm_res, n_genes, ncol, div_mat_sorted, q_vals_sorted,
            gene_names, metric, divergence_results_se))
    }

    # Pass SummarizedExperiment to global plot so it can use bootstrap CIs if
    # available
    return(.spectrum_plot_global(div_mat_sorted, q_vals_sorted, metric, variability_metric,
        divergence_results_se = divergence_results_se))
}
