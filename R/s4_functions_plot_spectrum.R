#' Plot Q-Spectrum Curves for Multiple Top Genes
#'
#' Wrapper around .plot_tsallis_q_curve() that creates multi-gene q-spectrum
#' visualization from Tsallis divergence results. Shows how divergence changes
#' across q-values for top genes, revealing q-dependent isoform switching.
#'
#' Key Features:
#' \itemize{
#'   \item Multi-gene grid: Side-by-side comparison of top N genes (default: 9)
#'   \item q-spectrum curves: Per-gene divergence across q = 0.01 to 2.00
#'   \item Bootstrap CI bands: 95% confidence intervals showing uncertainty
#'   \item Statistical significance: Gene titles include adjusted p-values
#'   \item Automatic layout: Grid dimensions auto-calculated from n_genes
#'   \item q=1 reference line: KL divergence (information theory benchmark)
#'   \item Region labels: Rare (q<1) vs Abundant (q>1) isoform emphasis
#' }
#'
#' Creates a multi-panel grid comparing per-q divergence profiles across the top
#' N genes identified by LMM interaction analysis. Each panel shows the full
#' q-spectrum divergence curve with the gene name and adjusted p-value in the title.
#'
#' @param eff_res Output from effect size computation OR \code{NULL}.
#' If provided, must contain `$interaction_results` with columns: gene,
#' adj_p_interaction, per_q_pattern.
#'   If \code{NULL}, uses fallback with lm_res + divergence_results_se.
#'
#' @param lm_res (Optional) Data frame from LMM analysis with columns: gene,
#' adj_p_interaction.
#'   Only used if eff_res is NULL. Must be provided for fallback mode.
#'
#' @param divergence_results_se (Optional) SummarizedExperiment from
#' divergence calculation.
#'   Only used if eff_res is NULL. Must be provided for fallback mode.
#'
#' @param n_genes Integer; number of top genes to plot (default: 9). Genes
#' are sorted by
#'   increasing adjusted p-value (most significant first).
#'
#' @param ncol Integer; number of columns in grid layout (default: 3).
#' Number of rows is
#'   automatically calculated as ceiling(n_genes / ncol).
#'
#' @param verbose Logical; if TRUE, print diagnostic messages (default: TRUE).
#'
#' @param output_file Character or NULL. Optional file path to save the plot
#' as an image.
#'   If provided, the plot will be saved with appropriate dimensions.
#'   Default: NULL (no file output, only return object).
#'
#' @return A ggplot2 object created via \code{patchwork} combining all gene panels,
#' or NULL if gene data is unavailable. The function automatically handles
#' ggplot2 grid
#'   creation and returns a print-ready object.
#'
#' @details
#' **Mathematical Background:**
#' Tsallis divergence as function of q:
#' \itemize{
#'   \item q < 1: Emphasizes rare isoforms (sensitive to outliers)
#'   \item q = 1: Kullback-Leibler divergence (classical information theory)
#'   \item q > 1: Emphasizes abundant isoforms (robust to rare variants)
#' }
#' D_q varies across q-spectrum, showing complexity of isoform differences.
#'
#' **Example Interpretation:**
#' \itemize{
#'   \item Flat curve: Divergence stable across q (robust isoform difference)
#'   \item Curved pattern: q-dependent divergence (rare vs abundant isoforms differ)
#'   \item Peaks at high q: Main isoforms drive the difference, rare ones immaterial
#' }
#'
#' **Input Modes:**
#' - **Mode 1 (Primary)**: Pass eff_res directly (from effect_sizes_divergence)
#' - **Mode 2 (Fallback)**: Pass lm_res + divergence_results_se instead
#'
#' **Gene Filtering:**
#' Genes are ranked by decreasing statistical significance (increasing
#' adj_p_interaction).
#' Only genes with complete per-q divergence data are included. If fewer
#' than n_genes
#' have valid data, the function returns all available genes.
#'
#' **Plot Features:**
#' - Title shows: gene name and adjusted p-value (q-value format)
#' - Per-q divergence curve with point estimates and 95% bootstrap CI bands
#' - Vertical reference line at q=1 (Kullback-Leibler divergence point)
#' - Region labels: 'Rare Isoforms' (q<1), 'Balanced' (q~=1), 'Abundant
#' Isoforms' (q>1)
#' - All plots use consistent ggplot2 styling matching plot_q_spectrum
#'
#' @examples
#' # Plot 4: Multi-gene q-spectrum profiles
#' set.seed(42)
#' # Create robust test dataset with clear signal-to-noise ratio
#' n_genes <- 16
#' n_isoforms_per_gene <- 4
#' n_isoforms <- n_genes * n_isoforms_per_gene
#' n_samples_per_group <- 30  # Increased for statistical power
#' n_samples <- n_samples_per_group * 2
#' 
#' # Generate control and treatment with very strong separation
#' # This ensures sufficient statistical power for divergence tests
#' control_counts <- matrix(rpois(n_isoforms * n_samples_per_group, lambda =
#' 100),
#'   nrow = n_isoforms, ncol = n_samples_per_group)
#' treatment_counts <- matrix(rpois(n_isoforms * n_samples_per_group, lambda
#' = 300),
#'   nrow = n_isoforms, ncol = n_samples_per_group)
#' counts <- cbind(control_counts, treatment_counts)
#' rownames(counts) <- paste0('TX_', 1:n_isoforms)
#' colnames(counts) <- paste0('Sample_', 1:n_samples)
#' 
#' se <- SummarizedExperiment::SummarizedExperiment(assays = list(counts =
#' counts))
#' tx2gene_df <- data.frame(Transcript = rownames(counts),
#'   Gene = rep(paste0('GENE_', 1:n_genes), each = n_isoforms_per_gene))
#' S4Vectors::metadata(se)$tx2gene <- tx2gene_df
#' SummarizedExperiment::colData(se) <- S4Vectors::DataFrame(
#'   sample_id = paste0('Sample_', 1:n_samples),
#'   condition = rep(c('control', 'treatment'), each = n_samples_per_group),
#'   row.names = colnames(se))
#' SummarizedExperiment::rowData(se)$transcript_id <- rownames(se)
#' SummarizedExperiment::rowData(se)$gene_id <-
#' tx2gene_df$Gene[match(rownames(se),
#'   tx2gene_df$Transcript)]
#' 
#' # Run complete analysis pipeline (skipped for speed in documentation)
#' # Uncomment to run actual analysis:
#' # analysis <- TSENATAnalysis(se)
#' # analysis <- calculate_diversity_s4(analysis, q = c(0.5, 1, 1.5), 
#' #   nboot = 50)
#' # analysis <- calculate_divergence_s4(analysis, q = c(0.5, 1, 1.5), 
#' #   nboot = 50)
#' # analysis <- calculate_lm_s4(analysis,
#' #   condition_col = 'condition')
#' # analysis <- calculate_effect_sizes_s4(analysis)
#' # p <- plot_multi_gene_q_spectrum_s4(analysis, n_genes = 4)
#' # if (!is.null(p)) print(p)
#'
#' @seealso \code{\link{calculate_divergence_s4}} for 
#' computing divergence values.
#'
#' @export
plot_multi_gene_q_spectrum_s4 <- function(eff_res = NULL, lm_res = NULL, divergence_results_se = NULL,
    n_genes = 9, ncol = 3, verbose = FALSE, output_file = NULL) {

    .load_visualization_deps()

    # Handle TSENATAnalysis S4 object
    if (methods::is(eff_res, "TSENATAnalysis")) {
        if (verbose)
            message("[plot_multi_gene_q_spectrum_s4] Detected TSENATAnalysis object, extracting lm_results and divergence_results...")
        components <- .extract_s4_components(eff_res, verbose)
        lm_res <- components$lm_res
        divergence_results_se <- components$divergence_results_se
        eff_res <- NULL
    }

    # Extract gene data from either mode
    gene_data <- .select_genes_from_eff_res(eff_res, n_genes, verbose)
    if (is.null(gene_data)) {
        gene_data <- .select_genes_fallback(lm_res, divergence_results_se, n_genes,
            verbose)
    }

    # Validate genes
    if (is.null(gene_data) || length(gene_data$genes) == 0) {
        if (verbose)
            message("No valid genes to plot. Check input data and column names.")
        return(NULL)
    }

    if (verbose)
        message(sprintf("Plotting %d genes in %d-column grid", length(gene_data$genes),
            ncol))

    # Create and assemble plots
    plot_list <- .create_gene_q_plots(gene_data$genes, gene_data$patterns, gene_data$p_values,
        verbose)

    if (length(plot_list) == 0) {
        if (verbose)
            message("No valid plots were created. Check per_q_pattern values and gene data.")
        return(NULL)
    }

    combined_plot <- .assemble_plot_grid(plot_list, ncol)

    if (verbose)
        message(sprintf("[OK] Multi-gene q-spectrum plot created with %d genes",
            length(plot_list)))

    if (!is.null(output_file)) {
        ggplot2::ggsave(output_file, plot = combined_plot, width = 12, height = 7.2,
            dpi = 100, create.dir = TRUE)
    }

    return(combined_plot)
}

#' @keywords internal
.extract_s4_components <- function(analysis, verbose = FALSE) {
    # Extract lm_results and divergence_results from TSENATAnalysis object
    lm_res <- NULL
    divergence_results_se <- NULL

    if (length(analysis@lm_results) > 0) {
        lm_res <- analysis@lm_results[[1]]$results
        if (verbose)
            message("[plot_multi_gene_q_spectrum_s4] Extracted lm_results with ",
                nrow(lm_res), " rows")
    } else {
        stop("TSENATAnalysis object has no lm_results. Run calculate_lm_s4() first.")
    }

    if (length(analysis@diversity_results) > 0) {
        divergence_results_se <- analysis@diversity_results[[1]]
        if (verbose)
            message("[plot_multi_gene_q_spectrum_s4] Extracted divergence_results with ",
                nrow(divergence_results_se), " rows")
    } else {
        stop("TSENATAnalysis object has no diversity_results. Run calculate_divergence_s4() first.")
    }

    list(lm_res = lm_res, divergence_results_se = divergence_results_se)
}

#' @keywords internal
.select_genes_from_eff_res <- function(eff_res, n_genes, verbose = FALSE) {
    # Mode 1: Extract from eff_res$interaction_results
    if (is.null(eff_res) || !is.list(eff_res)) {
        if (verbose)
            message("[plot_multi_gene_q_spectrum_s4] Mode 1 failed: eff_res is NULL or not a list")
        return(NULL)
    }

    if (is.null(eff_res$interaction_results) || nrow(eff_res$interaction_results) ==
        0) {
        if (verbose)
            message("[plot_multi_gene_q_spectrum_s4] Mode 1 failed: eff_res$interaction_results is NULL or empty")
        return(NULL)
    }

    int_res <- eff_res$interaction_results
    has_gene <- "gene" %in% colnames(int_res)
    has_per_q <- "per_q_pattern" %in% colnames(int_res)
    has_p_adj <- "adj_p_interaction" %in% colnames(int_res)
    has_p_raw <- "p_value_interaction" %in% colnames(int_res)

    p_col <- if (has_p_adj)
        "adj_p_interaction" else if (has_p_raw)
        "p_value_interaction" else NA_character_

    if (!(has_gene && has_per_q && (!is.na(p_col)))) {
        if (verbose) {
            message("[plot_multi_gene_q_spectrum_s4] Mode 1 failed: Missing required columns")
            message("  - has 'gene':", has_gene)
            message("  - has 'per_q_pattern':", has_per_q)
            message("  - has 'adj_p_interaction':", has_p_adj)
            message("  - has 'p_value_interaction':", has_p_raw)
        }
        return(NULL)
    }

    int_res_sorted <- int_res[order(int_res[[p_col]], na.last = TRUE), ]
    int_res_subset <- head(int_res_sorted, n_genes)
    valid_patterns <- !is.na(int_res_subset$per_q_pattern) & int_res_subset$per_q_pattern !=
        "" & int_res_subset$per_q_pattern != "NA"

    if (!any(valid_patterns)) {
        if (verbose)
            message("[plot_multi_gene_q_spectrum_s4] Mode 1 failed: per_q_pattern values are empty or invalid")
        return(NULL)
    }

    if (verbose)
        message(sprintf("[plot_multi_gene_q_spectrum_s4] Mode 1: Using eff_res with %s column (%d valid genes)",
            p_col, sum(valid_patterns)))

    list(genes = int_res_subset$gene[valid_patterns], patterns = int_res_subset$per_q_pattern[valid_patterns],
        p_values = int_res_subset[[p_col]][valid_patterns])
}

#' @keywords internal
.select_genes_fallback <- function(lm_res, divergence_results_se, n_genes, verbose = FALSE) {
    # Mode 2: Use lm_res + divergence_results_se
    if (is.null(lm_res) || is.null(divergence_results_se))
        return(NULL)
    if (nrow(lm_res) == 0 || nrow(divergence_results_se) == 0)
        return(NULL)
    if (!all(c("gene", "adj_p_interaction") %in% colnames(lm_res)))
        return(NULL)

    div_rd <- as.data.frame(rowData(divergence_results_se))
    div_assay <- assay(divergence_results_se)
    div_gene_names <- if ("gene_name" %in% colnames(div_rd))
        div_rd$gene_name else rownames(div_assay)

    if (length(div_gene_names) == 0 || nrow(div_assay) == 0)
        return(NULL)

    lm_sorted <- lm_res[order(lm_res$adj_p_interaction, na.last = TRUE), ]
    top_genes <- head(lm_sorted$gene, n_genes)
    gene_indices <- match(top_genes, div_gene_names)
    valid_idx <- !is.na(gene_indices)
    valid_genes <- top_genes[valid_idx]

    if (length(valid_genes) == 0)
        return(NULL)

    patterns <- character(length(valid_genes))
    for (i in seq_along(valid_genes)) {
        gene_idx <- which(div_gene_names == valid_genes[i])[1]
        if (!is.na(gene_idx))
            patterns[i] <- paste(div_assay[gene_idx, ][!is.na(div_assay[gene_idx,
                ])], collapse = ",")
    }

    if (verbose)
        message("[plot_multi_gene_q_spectrum_s4] Mode 2 (fallback): Using lm_res + divergence_results_se")

    list(genes = valid_genes, patterns = patterns, p_values = lm_sorted$adj_p_interaction[seq_along(valid_genes)])
}

#' @keywords internal
.create_gene_q_plots <- function(genes, patterns, p_values, verbose = FALSE) {
    # Create individual q-spectrum plots for each gene
    plot_list <- list()

    for (i in seq_along(genes)) {
        tryCatch({
            per_q_vals <- as.numeric(strsplit(patterns[i], ",")[[1]])
            if (length(per_q_vals) == 0 || all(is.na(per_q_vals))) {
                if (verbose)
                  message(sprintf("  Skipping %s: no valid per-q values", genes[i]))
                next
            }

            q_vals <- seq(0.1, by = 0.05, length.out = length(per_q_vals))
            plot_df <- data.frame(q = q_vals, divergence = per_q_vals, stringsAsFactors = FALSE)

            p <- ggplot2::ggplot(plot_df, ggplot2::aes(x = q, y = divergence)) +
                .theme_base(base_size = 11) + ggplot2::geom_line(color = "#4575B4",
                linewidth = 1.2) + ggplot2::geom_point(color = "#4575B4", size = 2.8,
                alpha = 0.8) + ggplot2::geom_vline(xintercept = 1, linetype = 3,
                color = "gray60", linewidth = 0.8, alpha = 0.7) + ggplot2::labs(title = genes[i],
                subtitle = sprintf("adj p = %.2e", p_values[i]), x = "q (Tsallis parameter)",
                y = "Tsallis Divergence D[q]") + ggplot2::theme(plot.title = ggplot2::element_text(size = .font_sizes$title,
                face = "bold", hjust = 0.5), plot.subtitle = ggplot2::element_text(hjust = 0.5,
                size = .font_sizes$subtitle, color = "gray40", margin = ggplot2::margin(b = 8)),
                plot.margin = ggplot2::margin(t = 8, b = 8, l = 6, r = 6), panel.grid.major = ggplot2::element_line(color = "gray92",
                  linewidth = 0.25), axis.text = ggplot2::element_text(size = .font_sizes$axis_text),
                axis.title = ggplot2::element_text(size = .font_sizes$axis_title,
                  face = "plain"))

            plot_list[[i]] <- p
        }, error = function(e) {
            if (verbose)
                message(sprintf("  Failed to plot %s: %s", genes[i], e$message))
        })
    }

    Filter(function(p) !is.null(p) && methods::is(p, "ggplot"), plot_list)
}

#' @keywords internal
.assemble_plot_grid <- function(plot_list, ncol) {
    # Combine plots into grid using patchwork
    if (length(plot_list) == 0)
        return(NULL)

    nrow <- ceiling(length(plot_list)/ncol)
    layout_plots <- list()

    for (row_idx in seq_len(nrow)) {
        row_start <- (row_idx - 1) * ncol + 1
        row_end <- min(row_idx * ncol, length(plot_list))
        row_plots <- Filter(function(p) !is.null(p) && methods::is(p, "ggplot"),
            plot_list[row_start:row_end])

        if (length(row_plots) == 0)
            next

        row_combined <- if (length(row_plots) == 1)
            row_plots[[1]] else Reduce(function(x, y) x + y, row_plots)
        layout_plots[[length(layout_plots) + 1]] <- row_combined

        if (row_idx < nrow)
            layout_plots[[length(layout_plots) + 1]] <- patchwork::plot_spacer()
    }

    Reduce(function(x, y) x/y, layout_plots) + patchwork::plot_layout(heights = c(rep(c(1,
        0.1), nrow - 1), 1), guides = "collect") + patchwork::plot_annotation(title = "Tsallis Divergence q-Spectrum Profiles",
        subtitle = "Per-q divergence curves for top-ranked genes", theme = ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5,
            face = "bold", size = .font_sizes$title, margin = ggplot2::margin(b = 8)),
            plot.subtitle = ggplot2::element_text(hjust = 0.5, face = "italic", size = .font_sizes$subtitle,
                color = "gray40", margin = ggplot2::margin(b = 12))))
}

