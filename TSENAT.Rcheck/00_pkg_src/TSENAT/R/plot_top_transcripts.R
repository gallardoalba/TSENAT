#' Plot top transcripts for a gene
#'
#' For a given gene, find transcripts using a tx->gene mapping, compute per-transcript
#' statistics between two sample groups, select the top N transcripts by p-value and
#' plot their expression across groups.
#'
#' @param counts A numeric matrix of transcript-level expression (rows = transcripts, columns = samples).
#' @param gene Character; gene symbol to inspect.
#' @param samples Character vector of sample group labels (length = ncol(counts)).
#' @param tx2gene Path to a two-column tab-delimited file with columns `Transcript` and `Gen`, or a data.frame with those columns. Required.
#' @param top_n Integer; number of transcripts to show (default = 3). If NULL, all transcripts for the gene are plotted.
#' @param pseudocount Numeric value added when computing log2 fold-change to avoid division by zero (default = 1e-6).
#' @param output_file Optional path to save the plot (ggsave will be used). If NULL the ggplot object is returned.
#' @return A `ggplot` object (or invisibly saved file if `output_file` provided).
#' @importFrom utils read.delim
#' @export
#' @name plot_top_transcripts
if (getRversion() >= "2.15.1") utils::globalVariables(c("tx", "expr", "group", "tx_cond", "sample", "log2expr"))
plot_top_transcripts <- function(counts, gene, samples, tx2gene = NULL, top_n = 3, pseudocount = 1e-6, output_file = NULL) {
  if (!is.matrix(counts) && !is.data.frame(counts)) stop("`counts` must be a matrix or data.frame with transcripts as rownames")
  counts <- as.matrix(counts)
  if (is.null(rownames(counts))) stop("`counts` must have rownames corresponding to transcript identifiers")
  if (length(samples) != ncol(counts)) stop("Length of `samples` must equal number of columns in `counts`")

  # tx2gene must be supplied as a path or data.frame
  if (is.character(tx2gene) && length(tx2gene) == 1) {
    if (!file.exists(tx2gene)) stop("tx2gene file not found: ", tx2gene)
    mapping <- utils::read.delim(tx2gene, stringsAsFactors = FALSE, header = TRUE)
  } else if (is.data.frame(tx2gene)) {
    mapping <- tx2gene
  } else {
    stop("`tx2gene` must be provided as a file path or data.frame")
  }

  if (!all(c("Transcript", "Gen") %in% colnames(mapping))) stop("tx2gene must have columns 'Transcript' and 'Gen'")

  if (!requireNamespace("ggplot2", quietly = TRUE)) stop("ggplot2 required for plotting")

  if (length(samples) != ncol(counts)) stop("Length of `samples` must equal number of columns in `counts`")

  make_plot_for_gene <- function(gene_single, fill_limits = NULL) {
    txs <- mapping$Transcript[mapping$Gen == gene_single]
    txs <- intersect(txs, rownames(counts))
    if (length(txs) == 0) stop("No transcripts found for gene: ", gene_single)
    if (!is.null(top_n)) txs <- head(txs, top_n)

    mat <- counts[txs, , drop = FALSE]
    df_all <- as.data.frame(mat)
    df_all$tx <- rownames(mat)
    df_long <- tidyr::pivot_longer(df_all, -tx, names_to = "sample", values_to = "expr")
    df_long$group <- rep(samples, times = length(txs))

    # summarize median per transcript x group and produce compact heatmap
    df_summary <- aggregate(expr ~ tx + group, data = df_long, FUN = function(x) stats::median(x, na.rm = TRUE))
    df_summary$log2expr <- log2(df_summary$expr + pseudocount)
    df_summary$tx <- factor(df_summary$tx, levels = unique(df_summary$tx))

    # show transcripts on y (readable labels) and groups on x
    p <- ggplot2::ggplot(df_summary, ggplot2::aes(x = group, y = tx, fill = log2expr)) +
      ggplot2::geom_tile(color = "white", width = 0.95, height = 0.95) +
      ggplot2::scale_fill_viridis_c(option = "viridis", direction = -1, na.value = "grey80", limits = fill_limits) +
      ggplot2::theme_minimal(base_size = 10) +
      ggplot2::labs(title = paste0(gene_single), x = NULL, y = NULL, fill = "log2(expr)") +
      ggplot2::theme(axis.text.y = ggplot2::element_text(size = 8),
                     axis.text.x = ggplot2::element_text(size = 8),
                     axis.ticks = ggplot2::element_blank(),
                     panel.grid = ggplot2::element_blank(),
                     plot.title = ggplot2::element_text(size = 10, hjust = 0.5),
                     legend.position = "bottom",
                     legend.key.width = ggplot2::unit(1.2, "cm"),
                     plot.margin = ggplot2::margin(4, 4, 4, 4)) +
      ggplot2::guides(fill = ggplot2::guide_colorbar(title.position = "top", barwidth = 6, barheight = 0.35))

    return(p)
  }

  # Produce plots (single or multiple). Do not save inside helper - save once below.
  if (length(gene) > 1) {
    # compute global fill limits across all genes so color scale is comparable
    mins <- c(); maxs <- c()
    for (g in gene) {
      txs <- mapping$Transcript[mapping$Gen == g]
      txs <- intersect(txs, rownames(counts))
      if (length(txs) == 0) next
      if (!is.null(top_n)) txs <- head(txs, top_n)
      mat <- counts[txs, , drop = FALSE]
      df_all <- as.data.frame(mat)
      df_all$tx <- rownames(mat)
      df_long <- tidyr::pivot_longer(df_all, -tx, names_to = "sample", values_to = "expr")
      df_long$group <- rep(samples, times = length(txs))
      df_summary <- aggregate(expr ~ tx + group, data = df_long, FUN = function(x) stats::median(x, na.rm = TRUE))
      df_summary$log2expr <- log2(df_summary$expr + pseudocount)
      mins <- c(mins, min(df_summary$log2expr, na.rm = TRUE))
      maxs <- c(maxs, max(df_summary$log2expr, na.rm = TRUE))
    }
    if (length(mins) == 0) stop("No transcripts found for provided genes")
    fill_limits <- c(min(mins, na.rm = TRUE), max(maxs, na.rm = TRUE))

    plots <- lapply(gene, function(g) make_plot_for_gene(g, fill_limits = fill_limits))
    # try patchwork first (collect guides), otherwise cowplot with shared legend
    if (requireNamespace("patchwork", quietly = TRUE)) {
      combined <- Reduce(`+`, plots) + patchwork::plot_layout(nrow = 1, guides = "collect") & ggplot2::theme(legend.position = "bottom")
      result_plot <- combined
    } else if (requireNamespace("cowplot", quietly = TRUE)) {
      # extract legend from first plot
      legend <- cowplot::get_legend(plots[[1]] + ggplot2::theme(legend.position = "bottom"))
      plots_nolegend <- lapply(plots, function(pp) pp + ggplot2::theme(legend.position = "none"))
      grid <- cowplot::plot_grid(plotlist = plots_nolegend, nrow = 1, align = "hv")
      result_plot <- cowplot::plot_grid(grid, legend, ncol = 1, rel_heights = c(1, 0.08))
    } else {
      # fallback: arrange grobs horizontally using base grid (no extra packages)
      # remove legends from individual plots and extract a single legend grob
      plots_nolegend <- lapply(plots, function(pp) pp + ggplot2::theme(legend.position = "none"))
      grobs <- lapply(plots_nolegend, ggplot2::ggplotGrob)
      # extract legend from original first plot
      g_full <- ggplot2::ggplotGrob(plots[[1]])
      legend_idx <- which(sapply(g_full$grobs, function(x) x$name) == "guide-box")
      legend_grob <- if (length(legend_idx)) g_full$grobs[[legend_idx[1]]] else NULL
      n <- length(grobs)
      # layout: plots in row 1, legend spanning row 2
      if (!is.null(output_file)) {
        png(filename = output_file, width = 800 * n, height = 420, res = 150)
        grid::grid.newpage()
        grid::pushViewport(grid::viewport(layout = grid::grid.layout(2, n, heights = grid::unit.c(grid::unit(1, "null"), grid::unit(0.7, "cm")))))
        for (i in seq_along(grobs)) {
          vp <- grid::viewport(layout.pos.row = 1, layout.pos.col = i)
          grid::pushViewport(vp)
          grid::grid.draw(grobs[[i]])
          grid::upViewport()
        }
        if (!is.null(legend_grob)) {
          vp_leg <- grid::viewport(layout.pos.row = 2, layout.pos.col = 1:n)
          grid::pushViewport(vp_leg)
          grid::grid.draw(legend_grob)
          grid::upViewport()
        }
        dev.off()
        return(invisible(NULL))
      } else {
        grid::grid.newpage()
        grid::pushViewport(grid::viewport(layout = grid::grid.layout(2, n, heights = grid::unit.c(grid::unit(1, "null"), grid::unit(0.7, "cm")))))
        for (i in seq_along(grobs)) {
          vp <- grid::viewport(layout.pos.row = 1, layout.pos.col = i)
          grid::pushViewport(vp)
          grid::grid.draw(grobs[[i]])
          grid::upViewport()
        }
        if (!is.null(legend_grob)) {
          vp_leg <- grid::viewport(layout.pos.row = 2, layout.pos.col = 1:n)
          grid::pushViewport(vp_leg)
          grid::grid.draw(legend_grob)
          grid::upViewport()
        }
        return(invisible(NULL))
      }
    }
  } else {
    result_plot <- make_plot_for_gene(gene)
  }

  if (!is.null(output_file)) {
    ggplot2::ggsave(output_file, result_plot)
    invisible(NULL)
  } else {
    return(result_plot)
  }
}
