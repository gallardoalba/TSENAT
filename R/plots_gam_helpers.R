# ============================================================================
# GAM PLOT HELPERS
# Extracted from plots_helpers.R — July 2026 refactoring (I11)
# ============================================================================


# ============================================================================

#' Build sample-to-group mapping from colData
#'
#' @param cdata SummarizedExperiment colData with sample metadata
#' @param condition_col Column name for group assignments
#' @return Named character vector: sample name -> group value

#' @noRd
.prepare_sample_group_mapping <- function(cdata, condition_col) {
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
#' @return Data frame with columns: sample, group, q, entropy (or NULL if
#' invalid)

#' @noRd
.plot_gam_prepare_gene_data <- function(gene, mat, sample_to_group) {
    if (!(gene %in% rownames(mat))) {
        return(NULL)
    }

    gene_vals <- mat[gene, ]
    col_names_full <- colnames(mat)

    # Parse column names: 'Sample_q=value'
    col_sample_names <- sub("_q=.*", "", col_names_full)
    col_q_values <- as.numeric(sub(".*_q=", "", col_names_full))

    # Look up group for each column
    col_groups <- unname(sample_to_group[col_sample_names])

    if (any(is.na(col_groups))) {
        return(NULL)
    }

    # Build long-format data frame
    plot_df <- data.frame(sample = col_sample_names, group = col_groups, q = col_q_values,
        entropy = as.numeric(gene_vals), stringsAsFactors = FALSE)

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
#' @return List with $plot_data and $pred_data data frames (or NULL if
#' fitting fails)

#' @noRd
.plot_gam_fit_group <- function(plot_df) {

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

        tryCatch({
            k <- min(10, max(2, round(nrow(subset_data)/2)))
            gam_fit <- mgcv::gam(entropy ~ s(q, k = k), data = subset_data)

            pred_data <- data.frame(q = pred_q)
            pred_vals <- stats::predict(gam_fit, newdata = pred_data, se.fit = TRUE)

            pred_list[[as.character(gr)]] <- data.frame(group = gr, q = pred_q, entropy_fit = pred_vals$fit,
                se = pred_vals$se.fit, stringsAsFactors = FALSE)
        }, error = function(e) {
            # Silently skip failed fits
            NULL
        })
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
#' @param sait_res Data frame with gene and p-value columns
#' @param genes Optional character vector of specific genes
#' @param n_top Number of top genes to select
#' @param sig_alpha Significance threshold
#' @return Character vector of gene IDs to plot (or NULL if none selected)

#' @noRd
.plot_select_genes <- function(sait_res, genes = NULL, n_top = 6, sig_alpha = 0.05) {
    if (!is.null(genes)) {
        if (!is.character(genes)) {
            stop("genes must be a character vector of gene names", call. = FALSE)
        }
        return(genes)
    }

    # Identify p-value column
    if ("adj_p_interaction" %in% colnames(sait_res)) {
        p_col <- "adj_p_interaction"
    } else if ("p_interaction" %in% colnames(sait_res)) {
        p_col <- "p_interaction"
    } else {
        stop("sait_res must contain 'adj_p_interaction' or 'p_interaction' column",
            call. = FALSE)
    }

    # Filter to significant genes (p-value < sig_alpha)
    sig_genes <- sait_res[sait_res[[p_col]] < sig_alpha, , drop = FALSE]

    # If no significant genes, return NULL
    if (nrow(sig_genes) == 0) {
        return(NULL)
    }

    # Sort by p-value and select top n_top genes
    top_idx <- order(sig_genes[[p_col]])[seq_len(min(n_top, nrow(sig_genes)))]
    top_genes <- sig_genes[top_idx, , drop = FALSE]

    # Return gene names sorted by p-value
    top_genes$gene
}

#' Prepare gene CI data for plotting
#'
#' Converts CI matrices into long-format data for gene-specific bootstrap CI
#' plotting.
#'
#' @param long_data Long-format data with Gene, group, q, tsallis columns
#' @param ci_lower_mat CI lower bounds matrix (genes x samples*q)
#' @param ci_upper_mat CI upper bounds matrix (genes x samples*q)
#' @param genes Character vector of gene IDs to extract
#'
#' @return Data frame with columns: Gene, group, q, median, ci_lower, ci_upper
#'
#' @noRd
.prepare_gene_ci_data <- function(long_data, ci_lower_mat, ci_upper_mat, genes) {

    # Aggregate to get median per gene, group, q
    stats_df <- dplyr::summarise(dplyr::group_by(long_data, Gene, group, q), median = median(tsallis,
        na.rm = TRUE), .groups = "drop")

    # Extract CI values for each gene, group, q combination
    plot_df <- stats_df
    plot_df$ci_lower <- NA_real_
    plot_df$ci_upper <- NA_real_

    # Map CI assay columns to gene/group/q combinations
    if (nrow(ci_lower_mat) > 0) {
        colnames_ci <- colnames(ci_lower_mat)

        # Parse column names (e.g., 'sample1_q=1.000')
        ci_samples <- sub("_q=.*", "", colnames_ci)
        ci_q_values <- as.numeric(sub(".*_q=", "", colnames_ci))

        for (i in seq_len(nrow(plot_df))) {
            g <- as.character(plot_df$Gene[i])
            gr <- as.character(plot_df$group[i])
            q_val <- as.numeric(plot_df$q[i])

            # Find indices in long_data for this gene/group/q
            matching_rows <- which(as.character(long_data$Gene) == g & as.character(long_data$group) ==
                gr & abs(as.numeric(as.character(long_data$q)) - q_val) < 1e-06)

            if (length(matching_rows) > 0) {
                # Get samples for this group from long_data
                samples_for_group <- unique(as.character(long_data$sample[matching_rows]))

                # Find CI columns for these samples at this q
                ci_col_mask <- (ci_samples %in% samples_for_group) & (abs(ci_q_values -
                  q_val) < 1e-06)
                ci_col_indices <- which(ci_col_mask)

                if (length(ci_col_indices) > 0) {
                  # Get CI bounds for these columns
                  gene_idx <- which(rownames(ci_lower_mat) == g)
                  if (length(gene_idx) > 0) {
                    # Use only first match (shouldn't have duplicates but be
                    # safe)
                    gene_idx <- gene_idx[1]
                    ci_lower_vals <- as.numeric(ci_lower_mat[gene_idx, ci_col_indices])
                    ci_upper_vals <- as.numeric(ci_upper_mat[gene_idx, ci_col_indices])

                    # Use median of CI values across samples in this group
                    plot_df$ci_lower[i] <- median(ci_lower_vals, na.rm = TRUE)
                    plot_df$ci_upper[i] <- median(ci_upper_vals, na.rm = TRUE)
                  }
                }
            }
        }
    }

    plot_df
}

# ============================================================================
# GAM PLOT HELPERS (for main function refactoring)
# ============================================================================

#' Handle and Validate Inputs for GAM Interaction Plot
#'
#' Validates SE object and handles flexible sait_res input formats.
#' Extracts results dataframe from list or validates dataframe directly.
#'
#' @param se A \code{SummarizedExperiment} object
#' @param sait_res Either a data.frame with 'gene' column or list with
#'   $results and $model_data
#'
#' @return List with validated components:
#'   - se: validated SummarizedExperiment
#'   - sait_res: extracted results dataframe
#'   - model_data: extracted model_data (or NULL)
#'
#' @noRd
.plot_gam_handle_inputs <- function(se, sait_res) {
    if (!inherits(se, "SummarizedExperiment")) {
        stop("se must be a SummarizedExperiment", call. = FALSE)
    }

    model_data <- NULL

    # Handle flexible input: sait_res can be either: 1. A data.frame with results
    # (traditional usage) 2. A list with $results and $model_data
    # (return_model_data = TRUE format)
    if (is.list(sait_res) && !is.data.frame(sait_res)) {
        # sait_res is a list with components
        if ("results" %in% names(sait_res) && is.data.frame(sait_res$results)) {
            # Extract results and model_data from the list
            extracted_results <- sait_res$results

            # If model_data provided, extract from sait_res
            if ("model_data" %in% names(sait_res)) {
                model_data <- sait_res$model_data
            }

            sait_res <- extracted_results
        } else {
            stop("sait_res is a list but does not contain 'results' data.frame component",
                call. = FALSE)
        }
    }

    if (!is.data.frame(sait_res) || !("gene" %in% colnames(sait_res))) {
        stop("sait_res must be either:\n  1. A data.frame with 'gene' column from .calculate_sait()\n  2. A list with $results and $model_data from return_model_data = TRUE",
            call. = FALSE)
    }

    if (nrow(sait_res) == 0) {
        stop("sait_res has no rows; .calculate_sait() returned no genes", call. = FALSE)
    }

    list(se = se, sait_res = sait_res, model_data = model_data)
}

#' Validate and Extract Q-Values from Model Data
#'
#' Validates model_data and extracts/normalizes q-values for GAM analysis.
#'
#' @param model_data List from .calculate_sait(...,
#' return_model_data = TRUE)$model_data
#'
#' @return Numeric vector of q-values
#'
#' @noRd
.plot_gam_extract_q_values <- function(model_data) {
    if (is.null(model_data)) {
        stop("model_data is required. Provide it as a parameter or pass full sait_res list with $model_data component",
            call. = FALSE)
    }

    if (!is.list(model_data)) {
        stop("model_data must be a list from .calculate_sait(..., return_model_data = TRUE)",
            call. = FALSE)
    }

    # Extract required metadata - handle both wrapped (from JSON) and unwrapped
    # formats JSON returns arrays: method = [['gam']], need [[1]] Direct list
    # returns: method = 'gam', no [[1]] needed
    q_values <- model_data$q_values
    if (is.null(q_values)) {
        stop("model_data must contain 'q_values' from the original analysis", call. = FALSE)
    }

    # Normalize q_values in case it's wrapped in list
    if (is.list(q_values) && length(q_values) == 1) {
        q_values <- unlist(q_values)
    } else {
        q_values <- unlist(q_values)
    }

    q_values
}

#' Match and Filter Genes Between SE and Results
#'
#' Finds genes present in both SE rownames and sait_res results.
#' Subsets both objects to matching genes only.
#'
#' @param se A \code{SummarizedExperiment}
#' @param sait_res Results dataframe with 'gene' column
#'
#' @return List with:
#'   - se: subset SE
#'   - sait_res: subset results
#'
#' @noRd
.plot_gam_match_genes <- function(se, sait_res) {
    # Match and filter genes between SE and sait_res After calculate_diversity,
    # rownames(SE) are gene names sait_res$gene column also contains gene names
    gene_names_in_results <- sait_res$gene
    gene_names_in_se <- rownames(se)

    # Find genes that exist in both
    available_genes <- gene_names_in_se[gene_names_in_se %in% gene_names_in_results]

    if (length(available_genes) == 0) {
        stop(sprintf("No genes from sait_res found in rownames(se). \n  Examples from sait_res: %s\n  Examples from SE: %s",
            paste(head(gene_names_in_results, 3), collapse = ", "), paste(head(gene_names_in_se,
                3), collapse = ", ")), call. = FALSE)
    }

    # Subset SE to only genes that are in results
    se <- se[available_genes, ]

    # Subset results to only genes that are in SE
    sait_res <- sait_res[sait_res$gene %in% rownames(se), ]

    list(se = se, sait_res = sait_res)
}

#' Create Gene ID to Display Name Mapping
#'
#' Builds mapping for display names from sait_res. Uses gene_name column if
#' available, otherwise uses gene IDs.
#'
#' @param sait_res Results dataframe with 'gene' column and optional
#' 'gene_name' column
#'
#' @return Named character vector mapping gene IDs to display names
#'
#' @noRd
.plot_gam_create_gene_map <- function(sait_res) {
    # Create gene ID to display name mapping from sait_res
    gene_name_map <- setNames(sait_res$gene, sait_res$gene)  # default: use gene ID

    # If sait_res has a gene_name column (e.g., from return_model_data), use it
    if ("gene_name" %in% colnames(sait_res)) {
        gene_name_map <- setNames(sait_res$gene_name, sait_res$gene)
    }

    gene_name_map
}

#' Arrange GAM Plot Grid with Title and Legend
#'
#' Creates final gridded layout of plots with title, plots, and legend.
#'
#' @param plots List of ggplot objects for each gene
#' @param condition_col Column name used for legend title
#' @param font_sizes List with legend font sizes (legend_title, legend_text)
#'
#' @return Single combined ggplot object
#'
#' @noRd
.plot_gam_arrange_grid <- function(plots, condition_col, font_sizes) {

    n_plots <- length(plots)
    n_cols <- 2
    n_rows <- ceiling(n_plots/n_cols)

    # Add margins to plots for spacing, particularly between rows
    plots_with_margins <- lapply(seq_along(plots), function(i) {
        p <- plots[[i]]
        # Add larger bottom margin for plots in the first row to create space
        # before second row
        if (i <= n_cols) {
            p <- p + ggplot2::theme(plot.margin = ggplot2::margin(b = 15, unit = "pt"))
        }
        p
    })

    # Extract legend from first plot
    p_for_legend <- .configure_legend(plots[[1]], position = "bottom", text_size = font_sizes$legend_text,
        title_size = font_sizes$legend_title)
    legend <- cowplot::get_legend(p_for_legend)

    # Create grid without legends
    combined_plot <- cowplot::plot_grid(plotlist = plots_with_margins, nrow = n_rows,
        ncol = n_cols, align = "hv", axis = "lr")

    # Add main title and subtitle above the grid
    title_plot <- .create_title_grob("q-curve: Top genes with group interaction",
        subtitle = "Fitted smooth curves by group", title_size = 20, subtitle_size = 16)

    # Combine title, plots, and single legend at bottom
    final_plot <- cowplot::plot_grid(title_plot, combined_plot, legend, nrow = 3,
        rel_heights = c(0.12, 1, 0.08))

    final_plot
}

#' Save Plot to File
#'
#' Saves a ggplot object to file with optional dimensions and adaptive font sizing.
#'
#' @param plot A ggplot object
#' @param output_file File path for output (if NULL, skips saving)
#' @param width Plot width in inches (default: 12)
#' @param height Plot height in inches (default: 10.3)
#'
#' @return Invisible NULL; plot saved as side effect
#'
#' @details
#' Uses .calculate_plot_dims() to compute consistent aspect ratios and
#' .scale_font_by_area() to adapt font sizes based on output dimensions.
#'
#' @noRd
.plot_gam_save_plot <- function(plot, output_file, width = NULL, height = NULL) {
    if (!is.null(output_file)) {
        # Use provided dimensions or defaults
        save_width <- if (is.null(width))
            12 else width
        save_height <- if (is.null(height))
            10.3 else height

        # Determine aspect type based on provided dimensions tall: height/width
        # ratio > 0.8 (e.g., 10/12 = 0.833) standard: height/width ratio <= 0.8
        # (e.g., 7.2/12 = 0.6)
        aspect_type <- if (save_height/save_width > 0.8)
            "tall" else "standard"

        # Calculate dimensions via .calculate_plot_dims for consistency
        plot_dims <- .calculate_plot_dims(width_inches = save_width, aspect_type = aspect_type,
            dpi_output = 100)

        # Compute adaptive font scaling based on actual area
        font_scale <- .scale_font_by_area(plot_dims$width, plot_dims$height)

        # Apply adaptive font scaling if dimensions deviate significantly from
        # reference (96 sq in) Only scale if deviation is >10% to avoid
        # excessive changes
        if (abs(font_scale - 1) > 0.1) {
            plot <- plot + ggplot2::theme(text = ggplot2::element_text(size = 11 *
                font_scale), plot.title = ggplot2::element_text(size = 14 * font_scale),
                axis.title = ggplot2::element_text(size = 12 * font_scale), axis.text = ggplot2::element_text(size = 10 *
                  font_scale), legend.text = ggplot2::element_text(size = 10 * font_scale),
                legend.title = ggplot2::element_text(size = 11 * font_scale))
        }

        ggplot2::ggsave(output_file, plot = plot, width = plot_dims$width, height = plot_dims$height,
            dpi = plot_dims$dpi, create.dir = TRUE)
    }
    invisible(NULL)
}

#' Create Single GAM Plot for One Gene
#'
#' Generates ggplot object for one gene with fitted GAM curves by group.
#'
#' @param gene Gene ID to plot
#' @param gene_display_name Display name for gene (if NULL, looked up from
#' map)
#' @param gene_name_map Named vector mapping gene IDs to display names
#' @param mat Assay matrix with samples in columns
#' @param sample_to_group Named vector mapping sample names to groups
#' @param condition_col Column name used for legend/facets
#'
#' @return ggplot object or NULL if plot generation fails
#'
#' @noRd
.plot_gam_make_plot <- function(gene, gene_display_name = NULL, gene_name_map, mat,
    sample_to_group, condition_col) {

    # Use provided gene name, or look it up from mapping, or default to gene ID
    if (is.null(gene_display_name)) {
        if (gene %in% names(gene_name_map)) {
            gene_display_name <- gene_name_map[[gene]]
        } else {
            gene_display_name <- gene
        }
    }

    # Prepare plot data for this gene using helper
    plot_df <- .plot_gam_prepare_gene_data(gene, mat, sample_to_group)
    if (is.null(plot_df)) {
        return(NULL)
    }

    # Fit GAM models and generate predictions using helper
    gam_result <- .plot_gam_fit_group(plot_df)
    if (is.null(gam_result)) {
        return(NULL)
    }

    plot_df <- gam_result$plot_data
    pred_df <- gam_result$pred_data
    group_levels <- gam_result$group_levels

    # Create explicit color mapping for all groups
    palette_colors <- .palette_blue_red()
    color_mapping <- c()
    for (i in seq_along(group_levels)) {
        color_idx <- ((i - 1)%%length(palette_colors)) + 1
        color_mapping[group_levels[i]] <- palette_colors[color_idx]
    }

    # Create plot with explicit color scale
    p <- ggplot2::ggplot(plot_df, ggplot2::aes(x = q, y = entropy, color = group)) +
        ggplot2::geom_point(data = plot_df, ggplot2::aes(x = q, y = entropy, color = group),
            alpha = 0.5, size = 2) + ggplot2::geom_line(data = pred_df, ggplot2::aes(x = q,
        y = entropy_fit, color = group, linetype = "Model fit"), linewidth = 1, alpha = 0.9) +
        ggplot2::scale_color_manual(values = color_mapping, name = condition_col,
            breaks = group_levels) + ggplot2::scale_linetype_manual(values = c(`Model fit` = 1),
        name = "") + ggplot2::labs(x = "q parameter", y = "Tsallis entropy", title = ifelse(gene_display_name !=
        gene, sprintf("%s (%s)", gene_display_name, gene), gene_display_name)) +
        .theme_spectrum(base_size = 12)

    p <- .configure_legend(p, position = "none")

    p
}

# ============================================================================
# PLOT COMPOSITION MEGA-HELPERS (Consolidation Phase)
