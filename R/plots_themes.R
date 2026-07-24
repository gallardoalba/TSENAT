# ============================================================================
# PLOT THEMES, STYLING, FORMATTING & GRID COMPOSITION
# Extracted from plots_helpers.R — July 2026 refactoring (I11)
# ============================================================================


# ============================================================================
# THEME UTILITIES
# ============================================================================

#' Apply TSENAT Base Theme
#'
#' Applies standard TSENAT styling: minimal theme with centered titles.
#'
#' @param base_size Integer: base font size (default: 11).
#' @param color_palette Character: 'blue_red' (default) or other.
#'
#' @return ggplot2 theme object.
#'

#' @noRd

.apply_tsenat_theme <- function(base_size = 11, color_palette = "blue_red") {

    theme_result <- .theme_base(base_size = base_size) + ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5,
        size = .font_sizes$title, face = "bold"), plot.subtitle = ggplot2::element_text(hjust = 0.5,
        size = .font_sizes$subtitle, face = "italic"))

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

.set_plot_title <- function(plot, title = NULL, subtitle = NULL, title_size = .font_sizes$title,
    subtitle_size = .font_sizes$subtitle) {

    if (!is.null(title)) {
        plot <- plot + ggplot2::labs(title = title)
    }

    if (!is.null(subtitle)) {
        plot <- plot + ggplot2::labs(subtitle = subtitle)
    }

    plot <- plot + ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5,
        size = title_size, face = "bold"), plot.subtitle = ggplot2::element_text(hjust = 0.5,
        size = subtitle_size, face = "italic"))

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

.create_tsenat_heatmap <- function(mat, title = NULL, colors = NULL, breaks = NULL,
    fontsize_row = 11, fontsize_col = 11, ...) {

    if (is.null(colors)) {
        colors <- .palette_continuous_diverging(n = 100)
    }

    pheatmap::pheatmap(mat, main = title, color = colors, breaks = breaks, fontsize_row = fontsize_row,
        fontsize_col = fontsize_col, ...)
}


# ============================================================================

#' Apply Publication-Ready Theme with Centered Titles
#'
#' Consolidated helper that applies base theme + centered/bolded title/subtitle.
#' Replaces repeated 35+ line pattern across all plot files.
#'
#' AESTHETIC PRESERVATION: Uses exact same font sizes and styling as originals.
#' - Title: .font_sizes$title, bold, centered
#' - Subtitle: .font_sizes$subtitle, italic, centered  
#' - Base theme: .theme_base (or .theme_spectrum for spectrum plots)
#'
#' @param plot ggplot2 object to style
#' @param title Character: plot title (optional)
#' @param subtitle Character: plot subtitle (optional)
#' @param base_theme Character: 'theme_base' (default) or 'theme_spectrum'
#' @param base_size Integer: base font size (default: 11, matches .theme_base default)
#' @param title_size Integer: title font size (default: from .font_sizes constants)
#' @param subtitle_size Integer: subtitle font size (default: from .font_sizes constants)
#'
#' @return Modified ggplot2 object with applied theme
#'
#' @noRd
.apply_publication_theme <- function(plot, title = NULL, subtitle = NULL, base_theme = "theme_base",
    base_size = 11, title_size = .font_sizes$title, subtitle_size = .font_sizes$subtitle) {

    # Apply base theme (either .theme_base or .theme_spectrum) Add dot prefix
    # if not already present
    theme_name <- if (startsWith(base_theme, "."))
        base_theme else paste0(".", base_theme)
    theme_fn <- get(theme_name)
    result <- plot + theme_fn(base_size = base_size)

    # Apply title/subtitle styling (always centered, bold/italic as per TSENAT
    # convention)
    result <- result + ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5,
        size = title_size, face = "bold"), plot.subtitle = ggplot2::element_text(hjust = 0.5,
        size = subtitle_size, face = "italic"))

    # Add title/subtitle labels if provided
    if (!is.null(title) && !is.null(subtitle)) {
        result <- result + ggplot2::labs(title = title, subtitle = subtitle)
    } else if (!is.null(title)) {
        result <- result + ggplot2::labs(title = title)
    } else if (!is.null(subtitle)) {
        result <- result + ggplot2::labs(subtitle = subtitle)
    }

    result
}

# ============================================================================
# PHASE 5: LEGEND & FORMATTING HELPERS
# ============================================================================

#' Configure Legend Positioning, Sizing, and Styling
#'
#' Consolidated helper that standardizes legend appearance across all plots.
#' Replaces repeated 20+ occurrences of legend.position, legend.key.width, 
#' legend.text, etc. customizations throughout the codebase.
#'
#' @param plot ggplot2 object to modify
#' @param position Character: 'bottom', 'right', 'left', 'top', or 'none' (default: 'bottom')
#' @param width_cm Numeric: width of legend key in cm (default: NULL = don't override)
#' @param height_cm Numeric: height of legend key in cm (default: NULL)
#' @param text_size Numeric: font size for legend text (default: NULL = use plot theme)
#' @param title_size Numeric: font size for legend title (default: NULL)
#' @param justification Character: 'left', 'center', 'right' (default: NULL = auto)
#' @param background_color Character: fill color for legend background (default: NULL)
#' @param border_color Character: border color for legend box (default: NULL)
#' @param spacing_lines Numeric: line spacing in legend (default: 2)
#'
#' @return Modified ggplot2 object with configured legend
#'
#' @examples
#' \dontrun{
#' # Simple usage: position and width
#' p <- ggplot2::ggplot(mtcars, ggplot2::aes(x = wt, y = mpg, color = factor(cyl))) +
#'     ggplot2::geom_point()
#' p_leg <- .configure_legend(p, position = 'bottom', width_cm = 2)
#'
#' # With text sizing
#' p_leg <- .configure_legend(p, position = 'right', text_size = 10, title_size = 11)
#'
#' # No legend
#' p_leg <- .configure_legend(p, position = 'none')
#' }
#'
#' @noRd

.configure_legend <- function(plot, position = "bottom", width_cm = NULL, height_cm = NULL,
    text_size = NULL, title_size = NULL, justification = NULL, background_color = NULL,
    border_color = NULL, spacing_lines = 2) {

    theme_list <- list()

    # Handle legend positioning
    if (!is.null(position) && position != "none") {
        theme_list$legend.position <- position
    } else if (position == "none") {
        theme_list$legend.position <- "none"
    }

    # Handle legend justification
    if (!is.null(justification)) {
        theme_list$legend.justification <- justification
    } else if (!is.null(position) && position == "bottom") {
        theme_list$legend.justification <- "center"
    }

    # Handle legend key dimensions
    if (!is.null(width_cm)) {
        theme_list$legend.key.width <- ggplot2::unit(width_cm, "cm")
    }
    if (!is.null(height_cm)) {
        theme_list$legend.key.height <- ggplot2::unit(height_cm, "cm")
    }

    # Handle text sizes
    if (!is.null(text_size)) {
        theme_list$legend.text <- ggplot2::element_text(size = text_size)
    }
    if (!is.null(title_size)) {
        theme_list$legend.title <- ggplot2::element_text(size = title_size, face = "bold")
    }

    # Handle legend background
    if (!is.null(background_color)) {
        theme_list$legend.background <- ggplot2::element_rect(fill = background_color,
            color = border_color %||% "black")
    } else if (!is.null(border_color)) {
        theme_list$legend.background <- ggplot2::element_rect(fill = NA, color = border_color)
    }

    # Handle spacing
    theme_list$legend.spacing.y <- ggplot2::unit(spacing_lines, "mm")

    if (length(theme_list) > 0) {
        plot <- plot + do.call(ggplot2::theme, theme_list)
    }

    return(plot)
}

#' Add Reference Lines (Horizontal and Vertical)
#'
#' Consolidated helper for adding reference/threshold lines to plots.
#' Replaces repeated 12+ occurrences of geom_hline + geom_vline patterns.
#'
#' @param plot ggplot2 object to modify
#' @param h_intercept Numeric vector: y-coordinates for horizontal lines (default: NULL)
#' @param v_intercept Numeric vector: x-coordinates for vertical lines (default: NULL)
#' @param h_color Character: color for horizontal lines (default: 'gray50')
#' @param v_color Character: color for vertical lines (default: 'gray50')
#' @param h_linetype Character: linetype for horizontal lines (default: 'dashed')
#' @param v_linetype Character: linetype for vertical lines (default: 'dashed')
#' @param h_size Numeric: line width for horizontal lines (default: 0.8)
#' @param v_size Numeric: line width for vertical lines (default: 0.8)
#' @param h_alpha Numeric: transparency for horizontal lines (default: 0.7)
#' @param v_alpha Numeric: transparency for vertical lines (default: 0.7)
#'
#' @return Modified ggplot2 object with reference lines
#'
#' @examples
#' \dontrun{
#' p <- ggplot2::ggplot(mtcars, ggplot2::aes(x = wt, y = mpg)) +
#'     ggplot2::geom_point()
#'
#' # Add threshold lines
#' p_ref <- .add_reference_lines(p, h_intercept = 20, v_intercept = 3)
#'
#' # Customized colors
#' p_ref <- .add_reference_lines(p, h_intercept = 20, h_color = 'red',
#'                               v_intercept = c(2.5, 3.5), v_color = 'blue')
#' }
#'
#' @noRd

.add_reference_lines <- function(plot, h_intercept = NULL, v_intercept = NULL, h_color = "gray50",
    v_color = "gray50", h_linetype = "dashed", v_linetype = "dashed", h_size = 0.8,
    v_size = 0.8, h_alpha = 0.7, v_alpha = 0.7) {

    # Add horizontal reference lines
    if (!is.null(h_intercept)) {
        for (yint in h_intercept) {
            plot <- plot + ggplot2::geom_hline(yintercept = yint, color = h_color,
                linetype = h_linetype, linewidth = h_size, alpha = h_alpha)
        }
    }

    # Add vertical reference lines
    if (!is.null(v_intercept)) {
        for (xint in v_intercept) {
            plot <- plot + ggplot2::geom_vline(xintercept = xint, color = v_color,
                linetype = v_linetype, linewidth = v_size, alpha = v_alpha)
        }
    }

    return(plot)
}

#' Format Axis Labels and Titles
#'
#' Consolidated helper for axis label styling including rotation, sizing, and face.
#' Replaces repeated 8+ occurrences of axis.text.x/y + axis.title customizations.
#'
#' @param plot ggplot2 object to modify
#' @param x_angle Numeric: rotation angle for x-axis labels (default: 0)
#' @param y_angle Numeric: rotation angle for y-axis labels (default: 0)
#' @param x_hjust Numeric: horizontal justification for x-axis (default: NULL = auto)
#' @param y_hjust Numeric: horizontal justification for y-axis (default: NULL = auto)
#' @param x_size Numeric: font size for x-axis labels (default: NULL = no override)
#' @param y_size Numeric: font size for y-axis labels (default: NULL = no override)
#' @param x_face Character: font face ('plain', 'bold', 'italic') for x-axis
#' @param y_face Character: font face for y-axis
#' @param x_color Character: text color for x-axis (default: 'black')
#' @param y_color Character: text color for y-axis (default: 'black')
#' @param bold_title Logical: make axis titles bold (default: TRUE)
#'
#' @return Modified ggplot2 object with formatted axis labels
#'
#' @examples
#' \dontrun{
#' p <- ggplot2::ggplot(mtcars, ggplot2::aes(x = factor(cyl), y = mpg)) +
#'     ggplot2::geom_point()
#'
#' # Horizontal x-axis labels at 90 degrees
#' p_fmt <- .format_axis_labels(p, x_angle = 90, x_size = 12, y_size = 11)
#'
#' # Bold axis titles only
#' p_fmt <- .format_axis_labels(p, bold_title = TRUE)
#' }
#'
#' @noRd

.format_axis_labels <- function(plot, x_angle = 0, y_angle = 0, x_hjust = NULL, y_hjust = NULL,
    x_size = NULL, y_size = NULL, x_face = "plain", y_face = "plain", x_color = "black",
    y_color = "black", bold_title = TRUE) {

    theme_list <- list()

    # X-axis label formatting
    if (x_angle != 0 || !is.null(x_size) || x_face != "plain" || x_color != "black") {
        x_hjust_use <- x_hjust %||% (if (x_angle != 0)
            1 else 0.5)
        x_vjust_use <- if (x_angle != 0)
            0.5 else 1
        theme_list$axis.text.x <- ggplot2::element_text(angle = x_angle, hjust = x_hjust_use,
            vjust = x_vjust_use, size = x_size, face = x_face, color = x_color)
    }

    # Y-axis label formatting
    if (y_angle != 0 || !is.null(y_size) || y_face != "plain" || y_color != "black") {
        y_hjust_use <- y_hjust %||% (if (y_angle != 0)
            1 else 0.5)
        theme_list$axis.text.y <- ggplot2::element_text(angle = y_angle, hjust = y_hjust_use,
            size = y_size, face = y_face, color = y_color)
    }

    # Axis titles
    if (bold_title) {
        theme_list$axis.title <- ggplot2::element_text(face = "bold")
    }

    if (length(theme_list) > 0) {
        plot <- plot + do.call(ggplot2::theme, theme_list)
    }

    return(plot)
}

#' Calculate Scaled Font Sizes
#'
#' Standardized font size calculation for proportional scaling across output dimensions.
#' Used when plots need to scale font sizes relative to output size (e.g., for heatmaps).
#'
#' @param base_size Numeric: base font size (default: 11)
#' @param scale_factor Numeric: overall scaling multiplier (default: 1)
#' @param font_multipliers List: named ratios relative to base size
#'
#' @return List with named elements: base, scaled, axis_text, axis_title, title, legend, subtitle, caption
#'
#' @examples
#' \dontrun{
#' # Standard font scaling
#' fonts <- .calculate_scaled_fonts(base_size = 11, scale_factor = 1.2)
#' # Use: fonts$title, fonts$axis_text, fonts$legend, etc.
#'
#' # Custom multipliers
#' custom_mult <- list(axis_text = 1.0, title = 1.8, legend = 0.8)
#' fonts <- .calculate_scaled_fonts(scale_factor = 2, font_multipliers = custom_mult)
#' }
#'
#' @noRd

.calculate_scaled_fonts <- function(base_size = 11, scale_factor = 1, font_multipliers = list(axis_text = 12/11,
    axis_title = 14/11, title = 16/11, legend = 9/11)) {

    scaling <- base_size * scale_factor

    result <- list(base = base_size, scaled = scaling, axis_text = round(scaling *
        font_multipliers$axis_text), axis_title = round(scaling * font_multipliers$axis_title),
        title = round(scaling * font_multipliers$title), legend = round(scaling *
            font_multipliers$legend), subtitle = round(scaling * 0.9), caption = round(scaling *
            0.8))

    return(result)
}

#' Apply Group Aesthetic Scales (Color + Fill + Legend)
#'
#' Consolidated helper for color palette + manual scales + legend styling.
#' Replaces repeated 21x pattern of palette + scale_color_manual + scale_fill_manual.
#'
#' AESTHETIC PRESERVATION: Uses exact same colors and mappings as originals.
#' - Palette: .palette_blue_red() by default (standard TSENAT convention)
#' - Color/Fill Manual: with name='Group' (standard legend title)
#'
#' @param plot ggplot2 object
#' @param palette Character: palette function name (e.g. 'palette_blue_red') OR 
#'   a vector of colors. If character, will call the .palette_* function.
#' @param legend_name Character: legend title (default: 'Group')
#' @param legend_position Character: legend position (default: 'bottom')
#' @param direction Integer: 1 (normal) or -1 (reversed palette)
#'
#' @return Modified ggplot2 object with applied color scales
#'
#' @noRd
# Get palette colors by name (pure function - cycle-free)
# BREAKING CYCLE: Extracted as independent function to allow testing without .apply_group_aesthetics
.get_palette_colors <- function(palette_name = "palette_blue_red") {
    if (is.character(palette_name) && length(palette_name) == 1) {
        # palette is a function name string, call the function
        tryCatch({
            palette_fn <- get(paste0(".", palette_name))  # e.g., .palette_blue_red
            return(palette_fn())
        }, error = function(e) {
            stop(sprintf("Palette '%s' not found", palette_name), call. = FALSE)
        })
    } else if (is.character(palette_name)) {
        # palette is a vector of color names/codes
        return(palette_name)
    } else {
        # palette is already a vector of colors
        return(palette_name)
    }
}

# Apply aesthetic colors to plot (pure function - cycle-free)
# BREAKING CYCLE: Core logic separated from palette lookup to enable independent testing
.apply_aesthetics_colors <- function(plot, colors, legend_name = "Group", 
                                    legend_position = "bottom", direction = 1) {
    # Validate colors input
    if (length(colors) == 0) {
        stop("colors must be a non-empty vector", call. = FALSE)
    }
    
    # Reverse if needed
    if (direction == -1) {
        colors <- rev(colors)
    }
    
    # Apply color and fill scales
    result <- plot + 
        ggplot2::scale_color_manual(values = colors, name = legend_name, na.value = "gray50") +
        ggplot2::scale_fill_manual(values = colors, name = legend_name, na.value = "gray50") + 
        ggplot2::theme(legend.position = legend_position)
    
    result
}

.apply_group_aesthetics <- function(plot, palette = "palette_blue_red", legend_name = "Group",
    legend_position = "bottom", direction = 1) {
    
    # Get colors (can be palette name or vector)
    colors <- .get_palette_colors(palette)
    
    # Apply aesthetics
    .apply_aesthetics_colors(plot, colors, legend_name, legend_position, direction)
}

#' Create Confidence Interval Line Plot Base
#'
#' Consolidated helper for ribbon + line + point layer pattern.
#' Replaces repeated 18x pattern of geom_ribbon + geom_line + geom_point.
#'
#' AESTHETIC PRESERVATION: Uses exact styling from q-curve plots:
#' - Ribbon: alpha=0.15, color=NA (transparent, no outline)
#' - Line: linewidth=1.2
#' - Point: size=3.5, alpha=0.8
#'
#' @param data Data frame with plot data
#' @param x_col Character: name of x column (default: 'q')
#' @param y_col Character: name of y column (default: 'median')
#' @param group_col Character: optional grouping column (NULL = single series, default: NULL)
#' @param ci_lower_col Character: name of CI lower column (default: 'ci_lower')
#' @param ci_upper_col Character: name of CI upper column (default: 'ci_upper')
#' @param ribbon_alpha Numeric: ribbon transparency (default: 0.15)
#' @param line_width Numeric: line width (default: 1.2)
#' @param point_size Numeric: point size (default: 3.5)
#' @param show_points Logical: include geom_point layer? (default: TRUE)
#'
#' @return Base ggplot2 object with ribbon/line/point layers (unthemed)
#'
#' @noRd
.create_ci_ribbon_plot <- function(data, x_col = "q", y_col = "median", group_col = NULL,
    ci_lower_col = "ci_lower", ci_upper_col = "ci_upper", ribbon_alpha = 0.15, line_width = 1.2,
    point_size = 2.8, show_points = TRUE, default_color = "#4575B4") {

    # Check if we have valid CI data to plot ribbons
    has_valid_ci <- FALSE
    if (ci_lower_col %in% colnames(data) && ci_upper_col %in% colnames(data)) {
        has_valid_ci <- any(!is.na(data[[ci_lower_col]]) & !is.infinite(data[[ci_lower_col]]) &
            !is.na(data[[ci_upper_col]]) & !is.infinite(data[[ci_upper_col]]))
    }

    # Only keep rows with valid CI data for ribbon layer to avoid ggplot
    # warnings But keep data as-is for line/point layers (they use median
    # values, not CI bounds)
    data_ci <- data
    if (has_valid_ci) {
        # Filter to rows with valid CI for ribbon layer only
        valid_ci_rows <- which(!is.na(data[[ci_lower_col]]) & !is.infinite(data[[ci_lower_col]]) &
            !is.na(data[[ci_upper_col]]) & !is.infinite(data[[ci_upper_col]]))
        data_ci <- data[valid_ci_rows, ]
    }

    # Build base aesthetics - include group color/fill only if group_col
    # provided and exists
    if (!is.null(group_col) && group_col %in% colnames(data)) {
        p <- ggplot2::ggplot(data, ggplot2::aes(x = .data[[x_col]], y = .data[[y_col]],
            color = .data[[group_col]], fill = .data[[group_col]], group = .data[[group_col]]))
        has_grouping <- TRUE
    } else {
        # No grouping - simple x/y aesthetics (color applied as fixed
        # aesthetic)
        p <- ggplot2::ggplot(data, ggplot2::aes(x = .data[[x_col]], y = .data[[y_col]]))
        has_grouping <- FALSE
    }

    # Check if we have valid CI data to plot
    has_valid_ci <- FALSE
    if (ci_lower_col %in% colnames(data) && ci_upper_col %in% colnames(data)) {
        has_valid_ci <- any(!is.na(data[[ci_lower_col]]) & !is.infinite(data[[ci_lower_col]]) &
            !is.na(data[[ci_upper_col]]) & !is.infinite(data[[ci_upper_col]]))
    }

    # Add ribbon layer (CI bounds) only if we have valid CI data
    if (has_valid_ci) {
        if (has_grouping) {
            # When grouping, fill aesthetic is inherited from base aes
            p <- p + ggplot2::geom_ribbon(ggplot2::aes(ymin = .data[[ci_lower_col]],
                ymax = .data[[ci_upper_col]]), alpha = ribbon_alpha, color = NA)
        } else {
            # No grouping - apply default fill color
            p <- p + ggplot2::geom_ribbon(ggplot2::aes(ymin = .data[[ci_lower_col]],
                ymax = .data[[ci_upper_col]]), alpha = ribbon_alpha, color = NA,
                fill = default_color)
        }
    }

    # Add line layer If group_col is provided, color aesthetic from aes()
    # applies it Otherwise, apply default color (consistency with
    # .create_simple_line_plot)
    if (has_grouping) {
        # When grouping, color aesthetic is inherited from base aes
        p <- p + ggplot2::geom_line(linewidth = line_width)
    } else {
        # No grouping - apply default color
        p <- p + ggplot2::geom_line(linewidth = line_width, color = default_color)
    }

    # Add point layer if requested If group_col is provided, color aesthetic
    # from aes() applies it Otherwise, apply default color
    if (show_points) {
        if (has_grouping) {
            # When grouping, color aesthetic is inherited from base aes
            p <- p + ggplot2::geom_point(size = point_size, alpha = 0.8)
        } else {
            # No grouping - apply default color
            p <- p + ggplot2::geom_point(size = point_size, alpha = 0.8, color = default_color)
        }
    }

    p
}

#' Assemble Multi-Plot Grid with Shared Legend
#'
#' Consolidated helper for legend extraction + grid composition pattern.
#' Replaces repeated 21x pattern of get_legend + plot_nolegend + plot_grid assembly.
#'
#' AESTHETIC PRESERVATION: Uses exact parameters from original code:
#' - Legend position: 'bottom' (default, customizable)
#' - Legend direction: 'horizontal' (standard for TSENAT)
#' - Grid alignment: 'hv' (both axes aligned)
#' - Title/subtitle: uses .font_sizes constants
#'
#' @param plots List of ggplot2 objects (one per subplot)
#' @param ncol Integer: number of columns (default: 2)
#' @param nrow Integer: number of rows (default: auto-calculated)
#' @param title Character: main title (optional)
#' @param subtitle Character: subtitle under title (optional)
#' @param legend_position Character: position for legend - 'bottom', 'top', 'left', 'right', 'none' 
#'   (default: 'bottom')
#' @param extract_legend Logical: whether to extract and place legend separately 
#'   (default: TRUE). If FALSE, plots retain their individual legends.
#' @param rel_heights Numeric vector: relative heights for title/plots/legend
#'   (default: c(0.08, 1, 0.08) - 8% title, 100% plots, 8% legend)
#'
#' @return Combined ggplot2/cowplot object ready for display/saving
#'
#' @details
#' Automatically:
#' - Extracts legend from first plot (if extract_legend=TRUE)
#' - Removes legends from all individual plots (if extract_legend=TRUE)
#' - Arranges in grid with specified layout
#' - Adds optional title/subtitle at top
#' - Places shared legend at specified position
#'
#' @noRd
.assemble_grid_plot <- function(plots, ncol = 2, nrow = NULL, title = NULL, subtitle = NULL,
    legend_position = "bottom", extract_legend = TRUE, rel_heights = c(0.08, 1, 0.08)) {

    if (length(plots) == 0) {
        stop("plots list cannot be empty", call. = FALSE)
    }

    # Calculate nrow if not provided
    if (is.null(nrow)) {
        nrow <- ceiling(length(plots)/ncol)
    }

    # Extract legend from first plot if requested
    legend_obj <- NULL
    if (extract_legend) {
        legend_obj <- cowplot::get_legend(plots[[1]] + ggplot2::theme(legend.position = legend_position,
            legend.direction = "horizontal"))
    }

    # Remove legends from all plots
    plots_nolegend <- lapply(plots, function(p) {
        .configure_legend(p, position = "none")
    })

    # Compose grid without legend
    grid_plot <- cowplot::plot_grid(plotlist = plots_nolegend, ncol = ncol, nrow = nrow,
        align = "hv")

    # If title/subtitle provided, create title grobs and assemble all 3
    # components
    if (!is.null(title) || !is.null(subtitle)) {
        title_plot <- cowplot::ggdraw()

        if (!is.null(title)) {
            title_plot <- title_plot + cowplot::draw_label(title, fontface = "bold",
                size = .font_sizes$title, x = 0.5, hjust = 0.5)
        }

        if (!is.null(subtitle)) {
            y_pos <- if (is.null(title))
                0.5 else 0.25
            subtitle_plot <- cowplot::ggdraw() + cowplot::draw_label(subtitle, fontface = "italic",
                size = .font_sizes$subtitle, x = 0.5, hjust = 0.5, color = "gray40")
            title_plot <- cowplot::plot_grid(title_plot, subtitle_plot, nrow = 2,
                rel_heights = c(1, 0.6))
        }

        # Assemble title + grid + legend (if extracted)
        if (extract_legend && !is.null(legend_obj)) {
            return(cowplot::plot_grid(title_plot, grid_plot, legend_obj, nrow = 3,
                rel_heights = rel_heights))
        } else {
            # No legend: just title + grid
            return(cowplot::plot_grid(title_plot, grid_plot, nrow = 2, rel_heights = rel_heights[c(1,
                2)]))
        }
    }

    # No title/subtitle: just grid + legend (if extracted)
    if (extract_legend && !is.null(legend_obj)) {
        return(cowplot::plot_grid(grid_plot, legend_obj, nrow = 2, rel_heights = rel_heights[c(2,
            3)]))
    } else {
        return(grid_plot)
    }
}

# ============================================================================
# PHASE 2 HELPERS: MEDIUM-IMPACT PATTERN CONSOLIDATION

# ============================================================================

#' Create Cowplot Title Grob with Optional Subtitle
#'
#' Consolidates cowplot::ggdraw() + draw_label() pattern (7x occurrences).
#' Creates a title-only or title+subtitle grob for use in grid layouts.
#'
#' @param title Character: main title text
#' @param subtitle Character: optional subtitle text
#' @param title_size Numeric: title font size (default: .font_sizes$title)
#' @param subtitle_size Numeric: subtitle font size (default: .font_sizes$subtitle)
#' @param title_face Character: title font face ('bold', 'italic', etc.)
#' @param subtitle_face Character: subtitle font face
#' @param title_color Character: title color (default: 'black')
#' @param subtitle_color Character: subtitle color (default: 'gray40')
#'
#' @return cowplot/ggplot2 grob object ready for plot_grid assembly
#'
#' @noRd
.create_title_grob <- function(title, subtitle = NULL, title_size = .font_sizes$title,
    subtitle_size = .font_sizes$subtitle, title_face = "bold", subtitle_face = "italic",
    title_color = "black", subtitle_color = "gray40") {

    # Start with title grob
    title_grob <- cowplot::ggdraw() + cowplot::draw_label(title, fontface = title_face,
        size = title_size, x = 0.5, hjust = 0.5, color = title_color)

    # Add subtitle if provided
    if (!is.null(subtitle)) {
        subtitle_grob <- cowplot::ggdraw() + cowplot::draw_label(subtitle, fontface = subtitle_face,
            size = subtitle_size, x = 0.5, hjust = 0.5, color = subtitle_color)

        # Combine title + subtitle
        title_grob <- cowplot::plot_grid(title_grob, subtitle_grob, nrow = 2, rel_heights = c(1,
            0.6))
    }

    title_grob
}

#' Apply Facet Styling with Panel Spacing and Strip Text
#'
#' Consolidates facet_wrap() + panel.spacing + strip.text pattern (8x occurrences).
#' Applies consistent faceting and panel styling across all plot types.
#'
#' @param plot ggplot2 object
#' @param ncol Integer: number of columns for facet layout
#' @param nrow Integer: number of rows (optional, usually auto-calculated)
#' @param facet_var Character: variable name to facet by (unquoted expression as string)
#' @param scales Character: 'fixed', 'free_x', 'free_y', or 'free' (default: 'free_y')
#' @param strip_text_size Numeric: font size for strip labels (default: .font_sizes$subtitle)
#' @param panel_spacing_lines Numeric: spacing between panels in lines (default: 1.5)
#'
#' @return Modified ggplot2 object with faceting and styling applied
#'
#' @noRd
.apply_facet_styling <- function(plot, ncol = 2, nrow = NULL, facet_var = NULL, scales = "free_y",
    strip_text_size = .font_sizes$subtitle, panel_spacing_lines = 1.5) {

    # Apply facet wrap if variable specified
    if (!is.null(facet_var)) {
        facet_formula <- stats::as.formula(paste0("~", facet_var))
        plot <- plot + ggplot2::facet_wrap(facet_formula, ncol = ncol, nrow = nrow,
            scales = scales)
    }

    # Apply panel and strip styling
    plot <- plot + ggplot2::theme(panel.spacing = ggplot2::unit(panel_spacing_lines,
        "lines"), strip.text = ggplot2::element_text(face = "bold", size = strip_text_size))

    plot
}

#' Prepare Long Format Data with Metadata Validation
#'
#' Wrapper around .prepare_tsallis_long() that adds validation and 
#' consistent parameter handling. Consolidates 10x data preparation pattern.
#'
#' @param se SummarizedExperiment: diversity/entropy assay object
#' @param assay_name Character: name of assay to use (default: 'diversity')
#' @param condition_col Character: column name for condition/group variable
#'   (optional, uses metadata config if NULL)
#' @param validate Logical: validate output structure (default: TRUE)
#'
#' @return Data frame in long format with columns: q, tsallis, group, Gene
#'
#' @noRd
.prepare_long_format <- function(se, assay_name = "diversity", condition_col = NULL,
    validate = TRUE) {

    # Use condition_col from metadata config if not provided
    if (is.null(condition_col)) {
        if (!is.null(S4Vectors::metadata(se)$condition_col)) {
            condition_col <- S4Vectors::metadata(se)$condition_col
        } else {
            condition_col <- "sample_type"  # TSENAT default
        }
    }

    # Prepare long format
    long <- .prepare_tsallis_long(se, assay_name = assay_name, condition_col = condition_col)

    # Validate structure
    if (validate) {
        required_cols <- c("q", "tsallis", "group", "Gene")
        missing_cols <- setdiff(required_cols, colnames(long))
        if (length(missing_cols) > 0) {
            warning("Prepared data missing columns: ", paste(missing_cols, collapse = ", "))
        }
    }

    long
}

#' Create Centered Theme Element Components
#'
#' Helper for extracting and applying centered/bolded title/subtitle styling.
#' Consolidates 12x centering + bolding theme pattern.
#'
#' @param include_title Logical: include title element (default: TRUE)
#' @param include_subtitle Logical: include subtitle element (default: TRUE)
#' @param title_size Numeric: title size (default: .font_sizes$title)
#' @param subtitle_size Numeric: subtitle size (default: .font_sizes$subtitle)
#' @param hjust Numeric: horizontal justification (default: 0.5 = centered)
#'
#' @return ggplot2::theme() object with centering/styling
#'
#' @noRd
.create_centered_theme <- function(include_title = TRUE, include_subtitle = TRUE,
    title_size = .font_sizes$title, subtitle_size = .font_sizes$subtitle, hjust = 0.5) {

    theme_list <- list()

    if (include_title) {
        theme_list$plot.title <- ggplot2::element_text(hjust = hjust, size = title_size,
            face = "bold")
    }

    if (include_subtitle) {
        theme_list$plot.subtitle <- ggplot2::element_text(hjust = hjust, size = subtitle_size,
            face = "italic")
    }

    do.call(ggplot2::theme, theme_list)
}



#' Create Publication-Ready Line Plot
#'
#' Consolidates simple line + point layer patterns (8x occurrences).
#' Creates line + point layers with optional grouping.
#'
#' @param data Data frame with x and y columns
#' @param x_col Character: x-axis column name
#' @param y_col Character: y-axis column name
#' @param group_col Character: optional grouping column (NULL = single series with fixed color)
#' @param points Logical: add point layer (default: TRUE)
#' @param line_width Numeric: line width (default: 1.2)
#' @param point_size Numeric: point size (default: 2.5)
#' @param alpha Numeric: transparency (default: 0.8)
#' @param line_color Character: fixed line color when no grouping (default: '#4575B4')
#'
#' @return ggplot2 object with line and optional point layers (unthemed)
#'
#' @noRd
.create_simple_line_plot <- function(data, x_col, y_col, group_col = NULL, points = TRUE,
    line_width = 1.2, point_size = 2, alpha = 0.8, line_color = "#4575B4") {

    # NO grouping: simple single-series plot with fixed color
    if (is.null(group_col)) {
        p <- ggplot2::ggplot(data, ggplot2::aes(x = !!rlang::sym(x_col), y = !!rlang::sym(y_col))) +
            ggplot2::geom_line(linewidth = line_width, alpha = alpha, color = line_color)

        if (isTRUE(points)) {
            p <- p + ggplot2::geom_point(size = point_size, alpha = alpha, color = line_color)
        }
        return(p)
    }

    # WITH grouping: map color to group column
    if (!(group_col %in% colnames(data))) {
        stop("Column '", group_col, "' not found in data")
    }

    p <- ggplot2::ggplot(data, ggplot2::aes(x = !!rlang::sym(x_col), y = !!rlang::sym(y_col),
        color = !!rlang::sym(group_col))) + ggplot2::geom_line(linewidth = line_width,
        alpha = alpha)

    if (isTRUE(points)) {
        p <- p + ggplot2::geom_point(size = point_size, alpha = alpha)
    }

    p
}

#' Prepare Grouped Long Format Data with Aggregation
#'
#' Consolidates data aggregation + IQR/SD calculation patterns (6x occurrences).
#' Transforms data from wide format with groups into long format with statistics.
#'
#' @param se SummarizedExperiment: input data
#' @param assay_name Character: assay to transform
#' @param group_by_col Character: column to group by
#' @param stat_funcs List of functions: functions to apply (default: median, IQR)
#'   Named list like list(median = median, iqr = function(x) diff(quantile(x, c(0.25, 0.75))))
#'
#' @return Data frame in long format with grouped statistics
#'   Columns: q (from rownames or metadata), group, value, lower, upper
#'
#' @noRd
.prepare_grouped_long_format <- function(se, assay_name = "diversity", group_by_col = "condition",
    stat_funcs = list(median = median, iqr = function(x) diff(quantile(x, c(0.25,
        0.75), na.rm = TRUE)))) {

    # Extract assay
    assay_mat <- SummarizedExperiment::assay(se, assay_name)

    # Get group info
    coldata <- SummarizedExperiment::colData(se)
    if (!group_by_col %in% colnames(coldata)) {
        stop("Column '", group_by_col, "' not found in colData")
    }
    groups <- coldata[[group_by_col]]

    # Aggregate by group
    long_list <- list()

    for (i in seq_len(nrow(assay_mat))) {
        gene_name <- rownames(assay_mat)[i]
        gene_data <- assay_mat[i, ]

        for (grp in unique(groups)) {
            grp_indices <- which(groups == grp)
            grp_values <- gene_data[grp_indices]
            grp_values <- grp_values[!is.na(grp_values)]

            if (length(grp_values) > 0) {
                median_val <- stat_funcs$median(grp_values)
                iqr_val <- stat_funcs$iqr(grp_values)

                long_list[[paste0(gene_name, "_", grp)]] <- data.frame(Gene = gene_name,
                  group = grp, value = median_val, lower = median_val - (iqr_val/2),
                  upper = median_val + (iqr_val/2), stringsAsFactors = FALSE)
            }
        }
    }

    do.call(rbind, long_list)
}  # ============================================================================
# PHASE 6: PLOT CONSOLIDATION HELPERS

# ============================================================================
# These helpers consolidate remaining high-value patterns identified in Phase
# 6: - Theme + aesthetics merging - Plot saving with unified dimensions -
# Distribution statistics abstraction - Bootstrap CI detection and extraction

#' Apply Publication Theme + Group Aesthetics Combined
#'
#' Consolidates the common pattern of applying both publication theme and
#' group-based color aesthetics in a single call.
#'
#' @param plot ggplot2 object to modify
#' @param title Character: plot title (optional)
#' @param base_size Numeric: base font size (default: 11)
#' @param base_theme Character: theme function name - 'theme_base' or 'theme_spectrum'
#' @param group_col Character: column name for group mapping (optional)
#' @param palette Character: palette name - 'blue_red', custom function name (default: 'blue_red')
#' @param group_levels Character vector: ordered factor levels (optional)
#' @param subtitle Character: plot subtitle (optional)
#'
#' @return Modified ggplot2 object
#'
#' @details
#' This function is a convenient wrapper around `.apply_publication_theme()` and
#' `.apply_group_aesthetics()` for the common case where both are applied sequentially.
#'
#' If `group_col` is NULL, only the publication theme is applied.
#'
#' @examples
#' \dontrun{
#' # Apply theme + group colors in one call
#' p <- ggplot2::ggplot(df, ggplot2::aes(x = q, y = entropy, color = group)) +
#'     ggplot2::geom_point()
#' p <- .apply_publication_aesthetics(p, title = 'Entropy Trend',
#'                                    base_size = 11, base_theme = 'theme_base',
#'                                    group_col = 'group', palette = 'blue_red')
#' }
#'
#' @noRd
.apply_publication_aesthetics <- function(plot, title = NULL, base_size = 11, base_theme = "theme_base",
    group_col = NULL, palette = "blue_red", group_levels = NULL, subtitle = NULL,
    legend_name = "Group", legend_position = "bottom") {

    # Apply publication theme first
    p <- .apply_publication_theme(plot, title = title, base_size = base_size, base_theme = base_theme,
        subtitle = subtitle)

    # Apply group aesthetics if group column specified
    if (!is.null(group_col)) {
        p <- .apply_group_aesthetics(p, palette = palette, legend_name = legend_name,
            legend_position = legend_position)
    }

    return(p)
}


# ============================================================================

#' Save Plot with Standard Dimensions
#'
#' Consolidates the common pattern of calculating plot dimensions and saving
#' with ggsave into a single function call.
#'
#' @param plot ggplot2 object to save
#' @param filename Character: output file path (PNG, PDF, etc.)
#' @param width_inches Numeric: chart width in inches (default: 12)
#' @param aspect_type Character: aspect ratio - 'standard' (16:9), 'square' (1:1),
#'   'wide' (21:9), 'tall' (9:16). Default: 'standard'
#' @param dpi_output Numeric: resolution in DPI (default: 100)
#' @param width_cm Numeric: override width in cm (optional)
#' @param height_cm Numeric: override height in cm (optional)
#'
#' @return Invisible NULL. Saves plot to file as side effect.
#'
#' @details
#' This function eliminates the repetitive pattern:
#' ```
#' plot_dims <- .calculate_plot_dims(width_inches = 12, aspect_type = 'standard')
#' ggplot2::ggsave(file, plot = p, width = plot_dims$width, height = plot_dims$height, dpi = plot_dims$dpi)
#' ```
#'
#' **Consolidation Impact**: 20+ occurrences × 3-4 lines = 60+ LOC saved
#'
#' @examples
#' \dontrun{
#' # Save plot with standard dimensions
#' p <- ggplot2::ggplot(mtcars, ggplot2::aes(x = wt, y = mpg)) +
#'     ggplot2::geom_point()
#' .save_plot_standard(p, 'output/my_plot.png', width_inches = 12)
#' }
#'
#' @noRd
.save_plot_standard <- function(plot, filename, width_inches = 12, aspect_type = "standard",
    dpi_output = 100, width_cm = NULL, height_cm = NULL) {

    # Calculate dimensions
    plot_dims <- .calculate_plot_dims(width_inches = width_inches, aspect_type = aspect_type,
        dpi_output = dpi_output)

    # Override with cm if provided (convert cm to inches)
    if (!is.null(width_cm)) {
        plot_dims$width <- width_cm/2.54
    }
    if (!is.null(height_cm)) {
        plot_dims$height <- height_cm/2.54
    }

    # Save plot
    ggplot2::ggsave(filename, plot = plot, width = plot_dims$width, height = plot_dims$height,
        dpi = plot_dims$dpi)

    invisible(NULL)
}

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

    # Use 2 columns (2 genes per row) with controlled spacing
    n_cols <- 2
    n_rows <- ceiling(length(plots)/n_cols)

    # Build rows of 2 plots each with spacing between columns
    plot_rows <- list()
    for (row in seq_len(n_rows)) {
        start_idx <- (row - 1) * n_cols + 1
        end_idx <- min(row * n_cols, length(plots))
        row_plots <- plots[start_idx:end_idx]

        # Add right margin to first plot to create column spacing
        if (length(row_plots) >= 1) {
            row_plots[[1]] <- row_plots[[1]] + ggplot2::theme(plot.margin = ggplot2::margin(r = 1,
                unit = "cm"))
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
    combined_plots_section <- Reduce(`/`, combined_elements) + patchwork::plot_layout(heights = heights_spec)

    # Create title annotation
    combined <- patchwork::plot_spacer()/combined_plots_section + patchwork::plot_annotation(title = "Transcript level expression",
        subtitle = paste0("Top genes with metric ", agg_label_unique), theme = ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5,
            size = .font_sizes$title, face = "bold", margin = ggplot2::margin(t = 10,
                b = 10)), plot.subtitle = ggplot2::element_text(hjust = 0.5, size = .font_sizes$subtitle,
            face = "italic", margin = ggplot2::margin(t = 5, b = 0.4)))) + patchwork::plot_layout(heights = c(0.045,
        1), guides = "collect")

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

    # Extract legend from first plot
    p_for_legend <- .configure_legend(plots[[1]], position = "bottom")
    legend <- cowplot::get_legend(p_for_legend)

    # Remove legends from all plots
    plots_nolegend <- lapply(plots, function(pp) {
        .configure_legend(pp, position = "none")
    })

    # Compose in 2 columns
    ncol <- 2
    nrow_val <- ceiling(length(plots_nolegend)/ncol)

    grid <- cowplot::plot_grid(plotlist = plots_nolegend, ncol = ncol, nrow = nrow_val,
        align = "hv")

    # Create title and subtitle
    title_grob <- cowplot::ggdraw() + cowplot::draw_label("Transcript level expression",
        fontface = "bold", x = 0.5, hjust = 0.5, size = .font_sizes$title)

    subtitle_grob <- cowplot::ggdraw() + cowplot::draw_label(paste0("Top genes with metric ",
        agg_label_unique), fontface = "italic", x = 0.5, hjust = 0.5, size = .font_sizes$subtitle,
        color = "gray40")

    spacer_grob <- cowplot::ggdraw() + ggplot2::theme_void()

    # Combine all elements
    result_plot <- cowplot::plot_grid(title_grob, subtitle_grob, spacer_grob, grid,
        legend, ncol = 1, rel_heights = c(0.05, 0.04, 0.0015, 1, 0.08), align = "h",
        axis = "l")

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

    # Remove legends from all plots
    plots_nolegend <- lapply(plots, function(pp) {
        .configure_legend(pp, position = "none")
    })

    # Convert to grobs
    grobs <- lapply(plots_nolegend, ggplot2::ggplotGrob)

    # Extract legend from first plot
    g_full <- ggplot2::ggplotGrob(plots[[1]])
    legend_idx <- which(vapply(g_full$grobs, function(x) x$name, character(1)) ==
        "guide-box")

    legend_grob <- if (length(legend_idx)) {
        g_full$grobs[[legend_idx[1]]]
    } else {
        NULL
    }

    # Default to 2 columns
    ncol <- min(2, length(grobs))
    nrow <- ceiling(length(grobs)/ncol)

    # Create height specification
    plot_heights <- list()
    for (i in seq_len(nrow)) {
        plot_heights[[length(plot_heights) + 1]] <- grid::unit(1, "null")
        if (i < nrow) {
            plot_heights[[length(plot_heights) + 1]] <- grid::unit(0.17, "cm")
        }
    }

    all_heights <- c(list(grid::unit(0.55, "cm")), plot_heights, list(grid::unit(0.7,
        "cm")))
    heights <- do.call(grid::unit.c, all_heights)

    if (!is.null(output_file)) {
        png_width <- 800 * ncol
        png_height <- 480 * nrow
        png(filename = output_file, width = png_width, height = png_height, res = 150)
        .draw_transcript_grid(grobs, agg_label_unique, legend_grob, ncol, heights,
            to_file = output_file)
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

.draw_transcript_grid <- function(grobs, agg_label_unique, legend_grob, ncol, heights,
    to_file = NULL) {

    nrow <- ceiling(length(grobs)/ncol)

    # Create viewport structure
    vp_top <- grid::viewport(x = 0, y = 0.95, width = 1, height = 0.05, just = c("left",
        "bottom"), name = "title")

    grid::pushViewport(vp_top)
    grid::grid.text("Transcript level expression", x = 0.5, y = 0.5, gp = grid::gpar(fontsize = .font_sizes$title,
        fontface = "bold"))
    grid::upViewport()

    # Draw plot grid
    grid::grid.newpage()
    grid::pushViewport(grid::viewport(layout = grid::grid.layout(nrow = nrow, ncol = ncol)))

    for (i in seq_along(grobs)) {
        row <- ((i - 1)%/%ncol) + 1
        col <- ((i - 1)%%ncol) + 1
        grid::pushViewport(grid::viewport(layout.pos.row = row, layout.pos.col = col))
        grid::grid.draw(grobs[[i]])
        grid::upViewport()
    }

    grid::upViewport()
}

