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
#'   theme element_text scale_color_manual scale_fill_manual geom_hline
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
#' # p <- .plot_tsallis_divergence_profile(
#' #   ts_se, gene = c("gene_1", "gene_2")
#' # )
#'

#' @noRd

.plot_tsallis_divergence_profile <- function(se,
                                            gene = NULL,
                                            lm_res = NULL,
                                            readcounts = NULL,
                                            tx2gene_map = NULL,
                                            group_col = "group",
                                            n_top = 3,
                                            assay_name = "diversity",
                                            arrange_type = c("facet", "list"),
                                            signed = TRUE) {
    # ===== INPUT VALIDATION =====
    
    # Validate SummarizedExperiment
    if (!inherits(se, "SummarizedExperiment")) {
        stop("Argument 'se' must be a SummarizedExperiment object, got: ", 
             class(se)[1])
    }
    
    # Validate assay_name exists
    if (!is.character(assay_name) || length(assay_name) != 1) {
        stop("Argument 'assay_name' must be a single character string, got: ", 
             class(assay_name)[1])
    }
    
    if (!assay_name %in% names(SummarizedExperiment::assays(se))) {
        stop("Assay '", assay_name, "' not found in SummarizedExperiment. ",
             "Available assays: ", paste(names(SummarizedExperiment::assays(se)), collapse = ", "))
    }
    
    # Validate arrange_type
    arrange_type <- match.arg(arrange_type)
    
    # Validate logical parameters
    if (!isTRUE(signed) && !isFALSE(signed)) {
        stop("Argument 'signed' must be TRUE or FALSE, got: ", signed)
    }
    
    # Validate numeric parameters
    if (!is.numeric(n_top) || length(n_top) != 1 || n_top < 1) {
        stop("Argument 'n_top' must be a positive integer, got: ", n_top)
    }
    n_top <- as.integer(n_top)
    
    # Validate group_col exists in colData
    if (!is.character(group_col) || length(group_col) != 1) {
        stop("Argument 'group_col' must be a single character string, got: ", 
             class(group_col)[1])
    }
    
    col_data <- as.data.frame(SummarizedExperiment::colData(se))
    if (!group_col %in% colnames(col_data)) {
        stop("Column '", group_col, "' not found in colData(se). ",
             "Available columns: ", paste(colnames(col_data), collapse = ", "))
    }
    
    # Validate groups
    groups <- unique(col_data[[group_col]])
    if (length(groups) != 2) {
        stop("Exactly 2 groups required in '", group_col, "' column; found ", 
             length(groups), " groups: ", paste(groups, collapse = ", "))
    }
    
    # Determine genes to plot
    genes <- .select_genesselect_genes(gene, lm_res, n_top)
    
    if (length(genes) == 0) {
        stop("No genes selected for plotting after validation")
    }
    
    # Extract q values from column names
    q_info <- .select_genesextract_q_values(se)
    unique_q <- q_info$unique_q
    q_values <- q_info$q_values
    
    if (length(unique_q) < 2) {
        stop("SummarizedExperiment must contain multiple q-values in column names ",
             "(expected format: *_q=0.5). Found ", length(unique_q), " unique q-value(s).")
    }
    
    # ===== COMPUTE DIVERGENCE =====
    
    plot_data <- .select_genescompute_divergence(
        se = se,
        genes = genes,
        q_values = q_values,
        unique_q = unique_q,
        groups = groups,
        group_col = group_col,
        assay_name = assay_name,
        readcounts = readcounts,
        tx2gene_map = tx2gene_map,
        signed = signed
    )
    
    if (nrow(plot_data) == 0) {
        stop("No valid divergence values computed. Check input data structure and parameters.")
    }
    
    # ===== BUILD PLOT =====
    
    if (arrange_type == "list") {
        return(.select_genesbuild_list_plots(plot_data, genes, groups, signed, assay_name))
    } else {
        return(.select_genesbuild_facet_plot(plot_data, groups, signed))
    }
}

# ===== HELPER FUNCTIONS =====

#' Select genes for plotting

#' @noRd

.select_genesselect_genes <- function(gene, lm_res, n_top) {
    if (!is.null(gene)) {
        # User provided specific genes
        if (!is.character(gene)) {
            stop("Argument 'gene' must be character vector, got: ", class(gene)[1])
        }
        return(as.character(unique(gene)))
    }
    
    if (is.null(lm_res)) {
        stop("Either 'gene' or 'lm_res' must be provided")
    }
    
    if (!is.data.frame(lm_res)) {
        stop("Argument 'lm_res' must be a data.frame, got: ", class(lm_res)[1])
    }
    
    if (!("gene" %in% colnames(lm_res))) {
        stop("'lm_res' must have a 'gene' column. Found columns: ",
             paste(colnames(lm_res), collapse = ", "))
    }
    
    # Find p-value column
    p_col <- NULL
    for (col in c("adj_p_lmm", "adj_p_interaction", "p_value_interaction")) {
        if (col %in% colnames(lm_res)) {
            p_col <- col
            break
        }
    }
    
    if (is.null(p_col)) {
        stop("'lm_res' must have one of: adj_p_lmm, adj_p_interaction, or p_value_interaction. ",
             "Found columns: ", paste(colnames(lm_res), collapse = ", "))
    }
    
    # Select top n_top genes by p-value
    genes_ordered <- unique(as.character(lm_res$gene[order(lm_res[[p_col]])]))
    return(head(genes_ordered, n_top))
}

#' Extract q values from SE column names

#' @noRd

.select_genesextract_q_values <- function(se) {
    col_names <- colnames(se)
    
    extract_q <- function(name) {
        if (grepl("_q=", name)) {
            # as.numeric() naturally produces NA for non-numeric strings
            as.numeric(gsub(".*_q=", "", name))
        } else {
            NA_real_
        }
    }
    
    q_values <- vapply(col_names, extract_q, FUN.VALUE = numeric(1))
    unique_q <- sort(unique(q_values[!is.na(q_values)]))
    
    list(q_values = q_values, unique_q = unique_q)
}

#' Compute divergence for all gene-q combinations

#' @noRd

.select_genescompute_divergence <- function(se, genes, q_values, unique_q, groups, group_col,
                                        assay_name, readcounts, tx2gene_map, signed) {
    plot_data_list <- list()
    
    for (gene_name in genes) {
        divergences <- vapply(unique_q, function(q) {
            .select_genescalc_div_for_gene_q(
                se = se,
                gene_name = gene_name,
                q_val = q,
                q_values = q_values,
                groups = groups,
                group_col = group_col,
                assay_name = assay_name,
                readcounts = readcounts,
                tx2gene_map = tx2gene_map,
                signed = signed
            )
        }, FUN.VALUE = numeric(1))
        
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
    
    # Remove NA values
    all_plot_data <- all_plot_data[!is.na(all_plot_data$divergence), ]
    
    # Add direction indicator for signed divergence
    if (nrow(all_plot_data) > 0 && signed) {
        all_plot_data$direction <- ifelse(
            all_plot_data$divergence > 0, 
            paste0("Positive: ", groups[2], " higher"),
            paste0("Negative: ", groups[1], " higher")
        )
    }
    
    all_plot_data
}

#' Calculate divergence for a single gene at a specific q value

#' @noRd

.select_genescalc_div_for_gene_q <- function(se, gene_name, q_val, q_values, groups, group_col,
                                         assay_name, readcounts, tx2gene_map, signed) {
    # Get columns matching this q value
    cols_q <- which(q_values == q_val)
    if (length(cols_q) == 0) {
        return(NA_real_)
    }
    
    # Extract data
    diversity_matrix <- SummarizedExperiment::assay(se, assay_name)
    if (!gene_name %in% rownames(diversity_matrix)) {
        return(NA_real_)
    }
    
    entropy_vals <- diversity_matrix[gene_name, cols_q]
    col_data <- as.data.frame(SummarizedExperiment::colData(se))
    group_vals <- col_data[[group_col]][cols_q]
    
    # TIER 1: True Tsallis divergence if readcounts available
    if (!is.null(readcounts) && !is.null(tx2gene_map)) {
        div <- tryCatch({
            calculate_tsallis_divergence_paired_gene(
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
                entropy_pred = data.frame(
                    entropy_pred = entropy_vals,
                    q = q_val,
                    group = group_vals,
                    stringsAsFactors = FALSE
                ),
                group_levels = groups
            )
        }, error = function(e) {
            NA_real_
        })
        
        if (!is.na(div)) {
            return(if (signed) div else abs(div))
        }
    }
    
    # TIER 2: Entropy-based approximation (fallback)
    group1_vals <- entropy_vals[group_vals == groups[1]]
    group2_vals <- entropy_vals[group_vals == groups[2]]
    
    if (length(group1_vals) == 0 || length(group2_vals) == 0) {
        return(NA_real_)
    }
    
    mean1 <- mean(group1_vals, na.rm = TRUE)
    mean2 <- mean(group2_vals, na.rm = TRUE)
    
    if (signed) {
        mean2 - mean1
    } else {
        max(0, abs(mean1 - mean2))
    }
}

#' Build faceted plot with all genes

#' @noRd

.select_genesbuild_facet_plot <- function(plot_data, groups, signed) {
    if (signed) {
        p <- ggplot2::ggplot(plot_data, ggplot2::aes(x = q, y = divergence, color = direction)) +
            ggplot2::geom_hline(yintercept = 0, linetype = "dashed", color = "gray50", linewidth = 0.8) +
            ggplot2::geom_line(linewidth = 1.1) +
            ggplot2::geom_point(size = 3, alpha = 0.7) +
            ggplot2::scale_color_manual(
                name = "Divergence Direction:",
                values = setNames(
                    c("#E63946", "#1D3557"),
                    c(paste0("Negative: ", groups[1], " higher"),
                      paste0("Positive: ", groups[2], " higher"))
                )
            ) +
            ggplot2::labs(
                title = "Tsallis Divergence Profile: Directional (Signed)",
                x = "q value (diversity scale parameter)",
                y = "Divergence D[q]"
            ) +
            .theme_spectrum(base_size = 11) +
            ggplot2::theme(
                plot.title = ggplot2::element_text(size = .font_sizes$title),
                axis.title = ggplot2::element_text(size = .font_sizes$axis_title),
                axis.text = ggplot2::element_text(size = .font_sizes$axis_text)
            )
    } else {
        p <- ggplot2::ggplot(plot_data, ggplot2::aes(x = q, y = divergence, color = gene)) +
            ggplot2::geom_line(linewidth = 1.1) +
            ggplot2::geom_point(size = 3, alpha = 0.7) +
            ggplot2::labs(
                title = "Tsallis Divergence Profile Across q-Spectrum",
                x = "q value (diversity scale parameter)",
                y = "Divergence D[q] (Absolute)",
                color = "Gene"
            ) +
            .theme_spectrum(base_size = 11) +
            ggplot2::theme(
                plot.title = ggplot2::element_text(size = .font_sizes$title),
                axis.title = ggplot2::element_text(size = .font_sizes$axis_title),
                axis.text = ggplot2::element_text(size = .font_sizes$axis_text)
            )
    }
    
    # Add faceting if multiple genes
    if (length(unique(plot_data$gene)) > 1) {
        p <- p + ggplot2::facet_wrap(~gene, scales = "free_y") +
            ggplot2::theme(legend.position = "top")
    }
    
    p
}

#' Build list of individual plots per gene

#' @noRd

.select_genesbuild_list_plots <- function(plot_data, genes, groups, signed, assay_name) {
    plots <- list()
    
    for (gene_name in genes) {
        df_gene <- plot_data[plot_data$gene == gene_name, ]
        if (nrow(df_gene) == 0) {
            next
        }
        
        if (signed) {
            p_gene <- ggplot2::ggplot(df_gene, ggplot2::aes(x = q, y = divergence, fill = direction, color = direction)) +
                ggplot2::geom_hline(yintercept = 0, linetype = "dashed", color = "gray50", linewidth = 0.8) +
                ggplot2::geom_line(linewidth = 1.2) +
                ggplot2::geom_point(size = 3.5, alpha = 0.8) +
                ggplot2::scale_color_manual(
                    name = "Divergence Direction:",
                    values = setNames(
                        c("#E63946", "#1D3557"),
                        c(paste0("Negative: ", groups[1], " higher"),
                          paste0("Positive: ", groups[2], " higher"))
                    )
                ) +
                ggplot2::scale_fill_manual(
                    name = "Divergence Direction:",
                    values = setNames(
                        c("#E63946", "#1D3557"),
                        c(paste0("Negative: ", groups[1], " higher"),
                          paste0("Positive: ", groups[2], " higher"))
                    )
                ) +
                ggplot2::labs(
                    title = paste("Divergence Profile (Signed):", gene_name),
                    x = "q value",
                    y = expression("Divergence D[q]"),
                    subtitle = paste0("Red = ", groups[1], " higher | Blue = ", groups[2], " higher")
                ) +
                .theme_spectrum(base_size = 11) +
                ggplot2::theme(
                    plot.title = ggplot2::element_text(hjust = 0.5, size = .font_sizes$title, face = "bold"),
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
                .theme_spectrum(base_size = 11) +
                ggplot2::theme(
                    plot.title = ggplot2::element_text(hjust = 0.5, size = .font_sizes$title, face = "bold"),
                    panel.grid.minor = ggplot2::element_blank()
                )
        }
        plots[[gene_name]] <- p_gene
    }
    
    plots
}
