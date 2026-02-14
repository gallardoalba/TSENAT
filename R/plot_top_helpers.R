# Helpers for plot_top_transcripts internals

.ptt_select_genes_from_res <- function(res, top_n) {
    if (is.null(res)) {
        stop("Either 'gene' or 'res' must be provided")
    }
    if (!("genes" %in% colnames(res))) {
        stop("Provided 'res' must contain a 'genes' column")
    }
    if ("adjusted_p_values" %in% colnames(res)) {
        ord <- order(res$adjusted_p_values, na.last = NA)
    } else if ("raw_p_values" %in% colnames(res)) {
        ord <- order(res$raw_p_values, na.last = NA)
    } else {
        ord <- seq_len(nrow(res))
    }
    genes_sel <- as.character(res$genes[ord])
    genes_sel <- unique(genes_sel)
    head(genes_sel, top_n)
}

.ptt_infer_samples_from_coldata <- function(coldata, counts, sample_type_col) {
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
        as.character(cdf[colnames(counts), sample_type_col])
    } else {
        sample_id_cols <- c("sample", "Sample", "sample_id", "id")
        sid <- intersect(sample_id_cols, colnames(cdf))
        if (length(sid) > 0) {
            sid <- sid[1]
            if (!all(colnames(counts) %in% as.character(cdf[[sid]]))) {
                stop("coldata sample id column does not match column names of counts")
            }
            row_ix <- match(colnames(counts), as.character(cdf[[sid]]))
            as.character(cdf[[sample_type_col]][row_ix])
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

.ptt_build_plot_from_summary <- function(df_summary, agg_label_unique, fill_limits = NULL) {
    p <- ggplot2::ggplot(df_summary, ggplot2::aes(x = group, y = tx, fill = log2expr)) +
        ggplot2::geom_tile(color = "white", width = 0.95, height = 0.95) + ggplot2::scale_fill_viridis_c(option = "viridis",
        direction = -1, na.value = "grey80", limits = fill_limits) + ggplot2::theme_minimal(base_size = 14) +
        ggplot2::labs(title = agg_label_unique, x = NULL, y = NULL, fill = "log2(expr)") +
        ggplot2::theme(axis.text.y = ggplot2::element_text(size = 12), axis.text.x = ggplot2::element_text(size = 12),
            plot.title = ggplot2::element_text(size = 16, hjust = 0.6), legend.position = "bottom",
            legend.key.width = ggplot2::unit(1.2, "cm"), plot.margin = ggplot2::margin(4,
                4, 4, 4)) + ggplot2::guides(fill = ggplot2::guide_colorbar(title.position = "top",
        barwidth = 6, barheight = 0.35))
    p
}

.ptt_combine_patchwork <- function(plots, agg_label_unique) {
    combined <- Reduce(`+`, plots) + patchwork::plot_layout(nrow = 1, guides = "collect") &
        ggplot2::theme(legend.position = "bottom")
    combined <- combined + patchwork::plot_annotation(title = agg_label_unique, theme = ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.6,
        size = 16, margin = ggplot2::margin(b = 10))))
    combined
}

.ptt_combine_cowplot <- function(plots, output_file = NULL, agg_label_unique) {
    p_for_legend <- plots[[1]] + ggplot2::theme(legend.position = "bottom")
    legend <- cowplot::get_legend(p_for_legend)
    plots_nolegend <- lapply(plots, function(pp) pp + ggplot2::theme(legend.position = "none"))
    grid <- cowplot::plot_grid(plotlist = plots_nolegend, nrow = 1, align = "hv")
    title_grob <- cowplot::ggdraw() + cowplot::draw_label(agg_label_unique, fontface = "plain",
        x = 0.6, hjust = 0.5, size = 16)
    result_plot <- cowplot::plot_grid(title_grob, grid, legend, ncol = 1, rel_heights = c(0.08,
        1, 0.08))
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
    n <- length(grobs)
    heights <- grid::unit.c(grid::unit(0.6, "cm"), grid::unit(1, "null"), grid::unit(0.7,
        "cm"))
    if (!is.null(output_file)) {
        png(filename = output_file, width = 800 * n, height = 480, res = 150)
        .draw_transcript_grid(grobs, agg_label_unique, legend_grob, n, heights, to_file = output_file)
        invisible(NULL)
    } else {
        .draw_transcript_grid(grobs, agg_label_unique, legend_grob, n, heights)
        invisible(NULL)
    }
}
