#' Plot Composition & Theme Helpers for TSENAT Visualization
#'
#' This module provides utilities for combining multiple plots,
#' managing themes, and creating consistent ggplot2 scales.
#' Enables modular, reusable visualization composition.
#'
#' @name plot_helpers

#' @noRd
NULL

# ============================================================================
# MULTI-PLOT COMPOSITION: Patchwork
# ============================================================================

#' Combine Transcript Plots with Patchwork
#'
#' Uses patchwork package to compose transcript-level plots
#' in a 2-column grid with proper spacing and legends.
#'
#' @param plots List of ggplot2 objects (one per gene).
#' @param agg_label_unique Character: aggregation metric label for subtitle.
#'
#' @return Combined patchwork plot object.
#'

#' @noRd

.combine_plots_patchwork <- function(plots, agg_label_unique) {
  require_pkgs("patchwork")

  # Use 2 columns (2 genes per row) with controlled spacing
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
      row_plots[[1]] <- row_plots[[1]] +
        ggplot2::theme(plot.margin = ggplot2::margin(r = 1.0, unit = "cm"))
    }

    # Use patchwork composition (| for horizontal)
    if (length(row_plots) == 1) {
      row_combined <- row_plots[[1]]
    } else if (length(row_plots) == 2) {
      row_combined <- row_plots[[1]] | row_plots[[2]]
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
      heights_spec <- c(heights_spec, 0.17)  # Reduced spacing
    }
  }

  # Combine all elements
  combined_plots_section <- Reduce(`/`, combined_elements) +
    patchwork::plot_layout(heights = heights_spec)

  # Create title annotation
  combined <- patchwork::plot_spacer() / combined_plots_section +
    patchwork::plot_annotation(
      title = "Transcript level expression",
      subtitle = paste0("Top genes with metric ", agg_label_unique),
      theme = ggplot2::theme(
        plot.title = ggplot2::element_text(
          hjust = 0.5,
          size = .tsenat_font_sizes$title,
          face = "bold",
          margin = ggplot2::margin(t = 10, b = 10)
        ),
        plot.subtitle = ggplot2::element_text(
          hjust = 0.5,
          size = .tsenat_font_sizes$subtitle,
          face = "italic",
          margin = ggplot2::margin(t = 5, b = 0.4)
        )
      )
    ) +
    patchwork::plot_layout(heights = c(0.045, 1), guides = "collect")

  return(combined)
}

# ============================================================================
# MULTI-PLOT COMPOSITION: Cowplot
# ============================================================================

#' Combine Transcript Plots with Cowplot
#'
#' Uses cowplot package to compose transcript-level plots
#' with proper titles and legend handling.
#'
#' @param plots List of ggplot2 objects (one per gene).
#' @param output_file Character: optional file path to save result.
#' @param agg_label_unique Character: aggregation metric label for subtitle.
#'
#' @return Combined plot object (or invisible NULL if output_file provided).
#'

#' @noRd

.combine_plots_cowplot <- function(plots, output_file = NULL, agg_label_unique) {
  require_pkgs("cowplot")

  # Extract legend from first plot
  p_for_legend <- plots[[1]] + ggplot2::theme(legend.position = "bottom")
  legend <- cowplot::get_legend(p_for_legend)

  # Remove legends from all plots
  plots_nolegend <- lapply(plots, function(pp) {
    pp + ggplot2::theme(legend.position = "none")
  })

  # Compose in 2 columns
  ncol <- 2
  nrow_val <- ceiling(length(plots_nolegend) / ncol)

  grid <- cowplot::plot_grid(
    plotlist = plots_nolegend,
    ncol = ncol,
    nrow = nrow_val,
    align = "hv"
  )

  # Create title and subtitle
  title_grob <- cowplot::ggdraw() +
    cowplot::draw_label("Transcript level expression",
      fontface = "bold",
      x = 0.5, hjust = 0.5,
      size = .tsenat_font_sizes$title
    )

  subtitle_grob <- cowplot::ggdraw() +
    cowplot::draw_label(paste0("Top genes with metric ", agg_label_unique),
      fontface = "italic",
      x = 0.5, hjust = 0.5,
      size = .tsenat_font_sizes$subtitle,
      color = "gray40"
    )

  spacer_grob <- cowplot::ggdraw() + ggplot2::theme_void()

  # Combine all elements
  result_plot <- cowplot::plot_grid(
    title_grob, subtitle_grob, spacer_grob, grid, legend,
    ncol = 1,
    rel_heights = c(0.05, 0.04, 0.0015, 1, 0.08),
    align = "h", axis = "l"
  )

  if (!is.null(output_file)) {
    ggplot2::ggsave(output_file, result_plot)
    return(invisible(NULL))
  }

  return(result_plot)
}

# ============================================================================
# MULTI-PLOT COMPOSITION: Grid (Low-level)
# ============================================================================

#' Combine Transcript Plots with Grid
#'
#' Uses base grid package for low-level composition.
#' Useful for fine-grained control over plot layout.
#'
#' @param plots List of ggplot2 objects (one per gene).
#' @param output_file Character: optional file path to save as PNG.
#' @param agg_label_unique Character: aggregation metric label for subtitle.
#'
#' @return Invisible NULL. Outputs to file or current device.
#'

#' @noRd

.combine_plots_grid <- function(plots, output_file = NULL, agg_label_unique) {
  require_pkgs("grid")

  # Remove legends from all plots
  plots_nolegend <- lapply(plots, function(pp) {
    pp + ggplot2::theme(legend.position = "none")
  })

  # Convert to grobs
  grobs <- lapply(plots_nolegend, ggplot2::ggplotGrob)

  # Extract legend from first plot
  g_full <- ggplot2::ggplotGrob(plots[[1]])
  legend_idx <- which(vapply(g_full$grobs,
    function(x) x$name,
    character(1)
  ) == "guide-box")

  legend_grob <- if (length(legend_idx)) {
    g_full$grobs[[legend_idx[1]]]
  } else {
    NULL
  }

  # Default to 2 columns
  ncol <- min(2, length(grobs))
  nrow <- ceiling(length(grobs) / ncol)

  # Create height specification
  plot_heights <- list()
  for (i in seq_len(nrow)) {
    plot_heights[[length(plot_heights) + 1]] <- grid::unit(1, "null")
    if (i < nrow) {
      plot_heights[[length(plot_heights) + 1]] <- grid::unit(0.17, "cm")
    }
  }

  all_heights <- c(
    list(grid::unit(0.55, "cm")),
    plot_heights,
    list(grid::unit(0.7, "cm"))
  )
  heights <- do.call(grid::unit.c, all_heights)

  if (!is.null(output_file)) {
    png_width <- 800 * ncol
    png_height <- 480 * nrow
    png(filename = output_file, width = png_width, height = png_height, res = 150)
    .draw_transcript_grid(grobs, agg_label_unique, legend_grob, ncol, heights,
      to_file = output_file
    )
    dev.off()
    return(invisible(NULL))
  } else {
    .draw_transcript_grid(grobs, agg_label_unique, legend_grob, ncol, heights)
    return(invisible(NULL))
  }
}

# ============================================================================
# GRID DRAWING UTILITY
# ============================================================================

#' Draw Transcript Grid Layout
#'
#' Internal utility for drawing transcript plots using base grid system.
#'
#' @param grobs List of ggplot grobs (one per plot).
#' @param agg_label_unique Character: metric label.
#' @param legend_grob Grid grob: legend to include.
#' @param ncol Integer: number of columns.
#' @param heights Unit specification: row heights.
#' @param to_file Character: optional file path (internal use).
#'

#' @noRd

.draw_transcript_grid <- function(grobs,
                                 agg_label_unique,
                                 legend_grob,
                                 ncol,
                                 heights,
                                 to_file = NULL) {
  require_pkgs("grid")

  nrow <- ceiling(length(grobs) / ncol)

  # Create viewport structure
  vp_top <- grid::viewport(
    x = 0, y = 0.95, width = 1, height = 0.05,
    just = c("left", "bottom"), name = "title"
  )

  grid::pushViewport(vp_top)
  grid::grid.text("Transcript level expression",
    x = 0.5, y = 0.5,
    gp = grid::gpar(
      fontsize = .tsenat_font_sizes$title,
      fontface = "bold"
    )
  )
  grid::upViewport()

  # Draw plot grid
  grid::grid.newpage()
  grid::pushViewport(grid::viewport(layout = grid::grid.layout(
    nrow = nrow,
    ncol = ncol
  )))

  for (i in seq_along(grobs)) {
    row <- ((i - 1) %/% ncol) + 1
    col <- ((i - 1) %% ncol) + 1
    grid::pushViewport(grid::viewport(
      layout.pos.row = row,
      layout.pos.col = col
    ))
    grid::grid.draw(grobs[[i]])
    grid::upViewport()
  }

  grid::upViewport()
}

# ============================================================================
# SCALE CREATORS
# ============================================================================

#' Create Discrete Color Scale for Heatmap
#'
#' Creates a discrete color scale using TSENAT's standard palette.
#'
#' @param palette Character: "blue_red" (default), "continuous_diverging", or custom.
#' @param direction Integer: 1 (default) or -1 to reverse colors.
#' @param name Character: legend title.
#'
#' @return ggplot2 scale object (ggplot2::scale_color_manual or similar).
#'

#' @noRd

.create_color_scale <- function(palette = "blue_red", direction = 1, name = NULL) {
  require_pkgs("ggplot2")

  if (palette == "blue_red") {
    colors <- .tsenat_palette_blue_red()
  } else if (palette == "continuous_diverging") {
    colors <- .tsenat_palette_continuous_diverging()
  } else if (is.character(palette)) {
    colors <- palette
  } else {
    colors <- .tsenat_palette_blue_red()
  }

  if (direction == -1) {
    colors <- rev(colors)
  }

  if (length(colors) > 1) {
    ggplot2::scale_color_manual(values = colors, name = name)
  } else {
    NULL
  }
}

#' Create Fill Scale for Heatmap
#'
#' Creates a fill scale for heatmap-style plots.
#'
#' @param palette Character: "blue_red" (default) or "continuous_diverging".
#' @param direction Integer: 1 (default) or -1 to reverse.
#' @param name Character: legend title.
#' @param breaks Integer: number of color breaks (default: 50).
#'
#' @return ggplot2 scale object.
#'

#' @noRd

.create_fill_scale <- function(palette = "blue_red",
                              direction = 1,
                              name = NULL,
                              breaks = 50) {
  require_pkgs("ggplot2")

  if (palette == "continuous_diverging") {
    colors <- .tsenat_palette_continuous_diverging(n = breaks)
  } else {
    colors <- .tsenat_palette_blue_red()
  }

  if (direction == -1) {
    colors <- rev(colors)
  }

  ggplot2::scale_fill_gradient(
    low = colors[1],
    high = colors[length(colors)],
    name = name
  )
}

# ============================================================================
# THEME UTILITIES
# ============================================================================

#' Apply TSENAT Base Theme
#'
#' Applies standard TSENAT styling: minimal theme with centered titles.
#'
#' @param base_size Integer: base font size (default: 11).
#' @param color_palette Character: "blue_red" (default) or other.
#'
#' @return ggplot2 theme object.
#'

#' @noRd

.apply_tsenat_theme <- function(base_size = 11, color_palette = "blue_red") {
  require_pkgs("ggplot2")

  theme_result <- .tsenat_theme_base(base_size = base_size) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(
        hjust = 0.5,
        size = .tsenat_font_sizes$title,
        face = "bold"
      ),
      plot.subtitle = ggplot2::element_text(
        hjust = 0.5,
        size = .tsenat_font_sizes$subtitle,
        face = "italic"
      )
    )
  
  return(theme_result)
}

#' Modify Plot Title and Subtitle
#'
#' Convenience function to update title/subtitle in existing plot.
#'
#' @param plot ggplot2 object.
#' @param title Character: new title.
#' @param subtitle Character: new subtitle.
#' @param title_size Integer: title font size (default: from constants).
#' @param subtitle_size Integer: subtitle font size (default: from constants).
#'
#' @return Modified ggplot2 object.
#'

#' @noRd

.set_plot_title <- function(plot,
                           title = NULL,
                           subtitle = NULL,
                           title_size = .tsenat_font_sizes$title,
                           subtitle_size = .tsenat_font_sizes$subtitle) {
  require_pkgs("ggplot2")

  if (!is.null(title)) {
    plot <- plot + ggplot2::labs(title = title)
  }

  if (!is.null(subtitle)) {
    plot <- plot + ggplot2::labs(subtitle = subtitle)
  }

  plot <- plot +
    ggplot2::theme(
      plot.title = ggplot2::element_text(
        hjust = 0.5,
        size = title_size,
        face = "bold"
      ),
      plot.subtitle = ggplot2::element_text(
        hjust = 0.5,
        size = subtitle_size,
        face = "italic"
      )
    )

  return(plot)
}

# ============================================================================
# HEATMAP STYLING
# ============================================================================

#' Create Heatmap with Unified Styling
#'
#' Wraps pheatmap with TSENAT standard settings.
#'
#' @param mat Matrix: data to plot.
#' @param title Character: plot title.
#' @param colors Vector: color specification for heatmap.
#' @param breaks Vector: color break points.
#' @param fontsize_row Integer: row label font size.
#' @param fontsize_col Integer: column label font size.
#' @param ... Additional arguments passed to pheatmap.
#'
#' @return Heatmap object from pheatmap.
#'

#' @noRd

.create_tsenat_heatmap <- function(mat,
                                  title = NULL,
                                  colors = NULL,
                                  breaks = NULL,
                                  fontsize_row = 11,
                                  fontsize_col = 11,
                                  ...) {
  require_pkgs("pheatmap")

  if (is.null(colors)) {
    colors <- .tsenat_palette_continuous_diverging(n = 100)
  }

  pheatmap::pheatmap(mat,
    main = title,
    color = colors,
    breaks = breaks,
    fontsize_row = fontsize_row,
    fontsize_col = fontsize_col,
    ...
  )
}
