#' Import dplyr and stats functions with proper precedence
#'
#' This documentation block explicitly imports specific functions from dplyr
#' and stats packages to avoid namespace conflicts. DPlyr functions like
#' filter, lag, and select need to take precedence over stats equivalents.
#'
#' @name package_imports
#' @noRd
#' @importFrom dplyr arrange filter group_by mutate pull select summarise %>%
#' @importFrom stats IQR aggregate anova aov chisq.test coef cor fitted friedman.test formula kruskal.test lm loess loess.control median model.frame model.matrix na.omit p.adjust pchisq pnorm pt qnorm quantile residuals rmultinom sd setNames t.test var weighted.mean wilcox.test xtabs
#' @importFrom utils capture.output
#' @importFrom rlang .data
#' @importFrom S4Vectors metadata
NULL

# # Declare package global variables to satisfy R CMD check NOTES
if (getRversion() >= "2.15.1") {
    utils::globalVariables(c("xval", "x", "y", "padj", "padj_num", "padj_clean",
        "label_flag", "sample_q", "qnum", "significant", "prcomp", "predict",
        ".data", "divergence", "lower", "upper", "central", "spread", 
        "ci_lower", "ci_upper", "entropy", "entropy_fit", "p_gam", "p_friedman",
        "agreement", "p_value", "method", "direction", "calculate_tsallis_divergence_paired_gene"))
}

# ============================================================================
# METHOD DEPENDENCIES MAP (Gap 2: Explicit dependency declarations)
# ============================================================================

#' Method Dependencies for TSENAT Orchestration
#'
#' Defines explicit dependencies between analysis methods.
#' Used by \code{\link{tsenat}} to validate method combinations and 
#' prevent invalid execution orders.
#'
#' @format Named list mapping method names to their required dependencies:
#' \describe{
#'   \item{\code{diversity}}{No dependencies}
#'   \item{\code{jackknife}}{Requires diversity}
#'   \item{\code{divergence}}{Requires diversity}
#'   \item{\code{q_interactions}}{Requires diversity}
#'   \item{\code{lm_interaction}}{Requires diversity}
#' }
#'
#' @keywords internal
#' @noRd
DEPENDENCIES <- list(
  diversity = character(0),              # No dependencies
  jackknife = "diversity",               # Requires diversity
  divergence = "diversity",              # Requires diversity
  q_interactions = "diversity",          # Requires diversity
  lm_interaction = "diversity"           # Requires diversity
)

#' Method Execution Order
#'
#' Recommended execution order for TSENAT methods.
#' Respects the dependency graph defined in \code{\link{DEPENDENCIES}}.
#'
#' @format Character vector with methods in dependency order.
#'
#' @keywords internal
#' @noRd
METHOD_ORDER <- c("diversity", "jackknife", "lm_interaction", "divergence", "q_interactions")

# ============================================================================
# COLOR UTILITIES - Publication-quality visualization standards
# ============================================================================

#' TSENAT Discrete Color Palette
#'
#' Returns a colorblind-safe discrete palette following RColorBrewer best practices.
#' Uses the "Dark2" palette which is tested for accessibility across color blindness types.
#' Supports up to 8 distinct categories.
#'
#' @param n Numeric; number of colors to return (default: 8, max: 8).
#'   If n > 8, wraps around by repeating palette.
#'
#' @return Character vector of hex color codes.
#'
#' @details
#' The Dark2 palette from RColorBrewer is selected because:
#' - Colorblind-safe (protanopia, deuteranopia, tritanopia tested)
#' - Suitable for publication in Nature, Science, Cell journals
#' - Perceptually distinct colors even when printed grayscale
#' - Works well at screen and print resolutions
#'
#' This palette is applied to all categorical variables (groups, conditions)
#' across TSENAT plots for visual consistency.
#'
#' @examples
#' # Get default 8 colors
#' pal <- .tsenat_palette_discrete()
#' 
#' # Get subset
#' pal_4 <- .tsenat_palette_discrete(4)
#'
#' @keywords internal
#' @noRd
.tsenat_palette_discrete <- function(n = 8) {
  palette_dark2 <- c(
    "#1B9E77", "#D95F02", "#7570B3", "#E7298A",
    "#66A61E", "#E6AB02", "#A6761D", "#666666"
  )
  
  if (n <= 0) {
    return(character(0))
  } else if (n <= 8) {
    palette_dark2[seq_len(n)]
  } else {
    # Wrap around if more than 8 colors requested
    palette_dark2[seq_len(n) %% 8 + 1]
  }
}

#' TSENAT Blue-Red Harmonized Color Palette
#'
#' Returns a discrete palette harmonized with the heatmap's blue-red diverging scheme.
#' Creates visual continuity across all plot types by using colors from the diverging palette.
#'
#' @param n Numeric; number of colors to return (default: 8).
#'   If n > 8, wraps around by repeating palette.
#'
#' @return Character vector of hex color codes harmonized with diverging scheme.
#'
#' @details
#' This palette incorporates:
#' - **Blues (primary)**: #4575B4 (heatmap blue), #74ADD1 (light blue)
#' - **Reds (secondary)**: #D73027 (heatmap red), #F46D43 (light red)
#' - **Support colors**: Dark blue, dark red, pale blue, pale yellow
#'
#' Recommended when you want categorical plots to match the heatmap visual theme.
#' Use when all plots should have cohesive blue-red aesthetic (e.g., comparing
#' upregulated vs downregulated genes across multiple visualizations).
#'
#' @keywords internal
#' @noRd
.tsenat_palette_blue_red <- function(n = 8) {
  palette_blue_red <- c(
    "#4575B4",  # Blue (heatmap primary)
    "#D73027",  # Red (heatmap primary)
    "#74ADD1",  # Light blue
    "#F46D43",  # Light red/orange
    "#1A5490",  # Dark blue
    "#A50026",  # Dark red
    "#ABD9E9",  # Pale blue
    "#FFFFBF"   # Pale yellow
  )
  
  if (n <= 0) {
    return(character(0))
  } else if (n <= 8) {
    palette_blue_red[seq_len(n)]
  } else {
    palette_blue_red[seq_len(n) %% 8 + 1]
  }
}

#' TSENAT Continuous Diverging Color Palette
#'
#' Returns a continuous diverging palette (blue-white-red) matching the heatmap scheme.
#' Generates smooth color gradients for continuous data visualization.
#'
#' @param n Numeric; number of color stops to generate (default: 100).
#'   Higher values produce smoother gradients.
#'
#' @return Character vector of hex color codes for gradient scale.
#'
#' @details
#' Creates a smooth diverging gradient:
#' - **Blue (#4575B4)**: Low/negative values (e.g., downregulated)
#' - **White (#FFFFFF)**: Neutral center (e.g., no change)
#' - **Red (#D73027)**: High/positive values (e.g., upregulated)
#'
#' This is the EXACT palette used in heatmaps for consistency.
#' Use with ggplot2::scale_color_gradientn() or similar.
#'
#' @keywords internal
#' @noRd
.tsenat_palette_continuous_diverging <- function(n = 100) {
  grDevices::colorRampPalette(c("#4575B4", "#FFFFFF", "#D73027"))(n)
}

#' TSENAT Significance Testing Color Scheme
#'
#' Returns a standardized color mapping for significance testing across all plots.
#' Harmonizes with the blue-red heatmap theme.
#'
#' @return Named character vector mapping significance to colors:
#'   - "non-significant": Grey (#CCCCCC) - de-emphasized
#'   - "significant": Red (#D73027) - matches heatmap red
#'
#' @details
#' Design choices:
#' - **Non-significant**: Light grey (#CCCCCC) for neutral, subtle appearance
#' - **Significant**: Red (#D73027) from heatmap palette for consistency
#'
#' Apply in ggplot2 with:
#' \code{ggplot2::scale_color_manual(values = .tsenat_significance_colors())}
#'
#' This replaces ad-hoc color choices like "black" and "red" across different plots,
#' ensuring visual consistency in significance indicator colors.
#'
#' @keywords internal
#' @noRd
.tsenat_significance_colors <- function() {
  c(
    "non-significant" = "#CCCCCC",  # Grey (neutral, understated)
    "significant" = "#D73027"        # Red (heatmap red, emphatic)
  )
}

#' TSENAT Base Plot Theme
#'
#' Returns a publication-ready ggplot2 theme with standardized visual properties.
#' Applies consistent font sizes, gridlines, borders, and color scales across all plots.
#'
#' @param base_size Numeric; base font size in points (default: 11).
#'   All other sizes scale proportionally from this.
#'
#' @return List of ggplot2 theme elements that can be added to plots with `+`.
#'
#' @details
#' The theme includes:
#' - **Base**: theme_minimal() for clean appearance
#' - **Title**: bold, centered, 1.3x base size with bottom margin
#' - **Grid**: major gridlines only (no minor)
#' - **Border**: subtle grey panel border (linewidth=0.3)
#' - **Colors**: RColorBrewer Dark2 palette for categorical data
#' - **Legend**: inherited positioning, customizable per plot if needed
#'
#' Font sizing follows publication standards:
#' - Body text: 11pt
#' - Titles: 14pt (1.3x)
#' - Axis titles: 12pt (1.1x)
#' - Axis labels: 10pt (0.9x)
#'
#' @examples
#' \dontrun{
#' # Apply to a ggplot2 plot
#' ggplot2::ggplot(data) +
#'   .tsenat_theme_base(base_size = 11) +
#'   ggplot2::geom_point()
#' }
#'
#' @keywords internal
#' @noRd
.tsenat_theme_base <- function(base_size = 11) {
  list(
    ggplot2::theme_minimal(base_size = base_size),
    ggplot2::theme(
      plot.title = ggplot2::element_text(
        hjust = 0.5, 
        face = "bold", 
        size = base_size * 1.3,
        margin = ggplot2::margin(b = 8)
      ),
      plot.subtitle = ggplot2::element_text(
        hjust = 0.5,
        face = "italic",
        size = base_size * 0.95
      ),
      axis.title = ggplot2::element_text(size = base_size * 1.1),
      axis.text = ggplot2::element_text(size = base_size * 0.9),
      panel.grid.minor = ggplot2::element_blank(),
      panel.border = ggplot2::element_rect(
        color = "grey85", 
        fill = NA, 
        linewidth = 0.3
      )
    )
  )
}

# ============================================================================
# UNIFIED FONT SIZING CONSTANTS - Applied across all visualizations
# ============================================================================

#' TSENAT Global Font Sizes
#'
#' Unified font sizing constants used across all plots to ensure
#' consistent, readable typography regardless of plot type.
#'
#' These constants define the complete font hierarchy for the package:
#' - larger titles (main plot titles)
#' - axis titles and labels
#' - legend text
#' - heatmap annotations
#'
#' @format List with named numeric elements (all in points):
#' \describe{
#'   \item{\code{title}}{Main plot title: 18pt (large, bold)}
#'   \item{\code{subtitle}}{Plot subtitle: 14pt (medium, italic)}
#'   \item{\code{axis_title}}{Axis labels (X/Y): 14pt}
#'   \item{\code{axis_text}}{Axis tick labels: 12pt}
#'   \item{\code{legend_title}}{Legend heading: 12pt}
#'   \item{\code{legend_text}}{Legend entries: 11pt}
#'   \item{\code{heatmap_main}}{Heatmap panel titles: 12pt}
#'   \item{\code{heatmap_labels}}{Heatmap row/col labels: 10pt}
#' }
#'
#' @details
#' All sizes are absolute points and do NOT scale with base_size.
#' This ensures heatmaps and other grid-based plots maintain
#' readability regardless of panel dimensions.
#'
#' Usage in plots:
#' \code{ggplot2::element_text(size = .tsenat_font_sizes$title)}
#'
#' @keywords internal
#' @noRd
.tsenat_font_sizes <- list(
  title = 19,           # Main plot title (increased from 18)
  subtitle = 15,        # Plot subtitle (increased from 14)
  axis_title = 15,      # X/Y axis titles (increased from 14)
  axis_text = 13,       # Axis tick labels (increased from 12)
  legend_title = 13,    # Legend heading (increased from 12)
  legend_text = 12,     # Legend entries (increased from 11)
  heatmap_main = 13,    # Heatmap panel title (increased from 12)
  heatmap_labels = 11   # Heatmap row/column labels (increased from 10)
)

# ============================================================================
# THEME VARIANTS - Specialized themes for different plot types
# ============================================================================

#' TSENAT Spectrum Theme Variant
#'
#' Returns a specialized variant of the base theme optimized for spectrum/profile plots.
#' Inherits from \code{\link{.tsenat_theme_base}} and adds spectrum-specific overrides.
#'
#' @param base_size Numeric; base font size in points (default: 11).
#'
#' @return List of ggplot2 theme elements that can be added to plots with `+`.
#'
#' @details
#' Spectrum-specific adjustments:
#' - Right-aligned legend for profile plots
#' - Optional major gridlines for q-value axis
#' - Wider plot margins for axis labels
#'
#' Used by: plot_tsallis_q_curve_s4(), plot_lm_interaction_gam(),
#'          plot_tsallis_divergence_profile() and similar spectrum/profile plots.
#'
#' @keywords internal
#' @noRd
.tsenat_theme_spectrum <- function(base_size = 11) {
  list(
    .tsenat_theme_base(base_size = base_size),
    ggplot2::theme(
      legend.position = "right",
      panel.grid.major.y = ggplot2::element_line(color = "gray90", linewidth = 0.25),
      plot.margin = ggplot2::margin(t = 5, r = 8, b = 5, l = 5, unit = "mm")
    )
  )
}

#' TSENAT Heatmap Theme Variant
#'
#' Returns a specialized variant of the base theme optimized for heatmap visualizations.
#' Inherits from \code{\link{.tsenat_theme_base}} and adds heatmap-specific overrides.
#'
#' @param base_size Numeric; base font size in points (default: 11).
#'
#' @return List of ggplot2 theme elements that can be added to plots with `+`.
#'
#' @details
#' Heatmap-specific adjustments:
#' - No gridlines (not applicable to heatmaps)
#' - Compact margins to maximize heatmap area
#' - Bottom legend position for multi-panel layouts
#' - Reduced plot title margins
#'
#' Used by: plot_multiq_delta_influence_heatmaps_s4(), 
#'          plot_top_transcripts_s4() when using ggplot2 heatmap geoms.
#'
#' @keywords internal
#' @noRd
.tsenat_theme_heatmap <- function(base_size = 11) {
  list(
    .tsenat_theme_base(base_size = base_size),
    ggplot2::theme(
      panel.grid.major = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank(),
      legend.position = "bottom",
      plot.margin = ggplot2::margin(t = 3, r = 3, b = 3, l = 3, unit = "mm"),
      plot.title = ggplot2::element_text(
        hjust = 0.5,
        face = "bold",
        size = base_size * 1.3,
        margin = ggplot2::margin(b = 4)
      )
    )
  )
}

#' TSENAT Distribution Theme Variant
#'
#' Returns a specialized variant of the base theme optimized for distribution plots
#' (violin, box, histogram, density).
#'
#' @param base_size Numeric; base font size in points (default: 11).
#'
#' @return List of ggplot2 theme elements that can be added to plots with `+`.
#'
#' @details
#' Distribution-specific adjustments:
#' - Major gridlines on y-axis for easier value reading
#' - Legend on the right for comparison groups
#' - Rotated x-axis labels if many categories
#'
#' Used by: plot_tsallis_violin_density_grid_s4() and similar distribution plots.
#'
#' @keywords internal
#' @noRd
.tsenat_theme_distribution <- function(base_size = 11) {
  list(
    .tsenat_theme_base(base_size = base_size),
    ggplot2::theme(
      legend.position = "right",
      panel.grid.major.y = ggplot2::element_line(color = "gray90", linewidth = 0.25),
      axis.text.x = ggplot2::element_text(
        angle = 45,
        hjust = 1,
        vjust = 1,
        size = base_size * 0.85
      )
    )
  )
}

# ============================================================================
# SPECIALIZED THEME VARIANTS - For specific plot types
# ============================================================================

#' TSENAT Spectrum Plot Theme
#'
#' Theme variant for q-spectrum (diversity/divergence) plots.
#' Adds gridlines for easier value reading.
#'
#' @param base_size Numeric; base font size in points (default: 11).
#'
#' @return List of ggplot2 theme elements.
#'
#' @keywords internal
#' @noRd
.tsenat_theme_spectrum <- function(base_size = 11) {
  list(
    .tsenat_theme_base(base_size = base_size),
    ggplot2::theme(
      panel.grid.major.y = ggplot2::element_line(color = "gray90", linewidth = 0.3),
      legend.position = "right"
    )
  )
}

#' TSENAT Heatmap Theme
#'
#' Theme variant for heatmap plots (via ComplexHeatmap or pheatmap).
#' Minimal borders, adjusted annotation text.
#'
#' @param base_size Numeric; base font size in points (default: 11).
#'
#' @return List of ggplot2 theme elements.
#'
#' @keywords internal
#' @noRd
.tsenat_theme_heatmap <- function(base_size = 11) {
  list(
    .tsenat_theme_base(base_size = base_size),
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1, size = base_size * 0.85),
      axis.text.y = ggplot2::element_text(size = base_size * 0.85),
      panel.grid = ggplot2::element_blank(),
      legend.position = "bottom"
    )
  )
}

#' TSENAT Distribution Plot Theme
#'
#' Theme variant for distribution plots (violin, density, boxplot).
#' Optimized for categorical comparisons.
#'
#' @param base_size Numeric; base font size in points (default: 11).
#'
#' @return List of ggplot2 theme elements.
#'
#' @keywords internal
#' @noRd
.tsenat_theme_distribution <- function(base_size = 11) {
  list(
    .tsenat_theme_base(base_size = base_size),
    ggplot2::theme(
      panel.grid.major.x = ggplot2::element_blank(),
      panel.grid.major.y = ggplot2::element_line(color = "gray90", linewidth = 0.3),
      legend.position = "bottom",
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1)
    )
  )
}

# ============================================================================
# SCALE WRAPPER FUNCTIONS - Consistent color/fill scales
# ============================================================================

#' TSENAT Discrete Color Scale
#'
#' Wrapper for ggplot2 color scale using the TSENAT discrete palette.
#'
#' @param ... Additional arguments passed to \code{ggplot2::scale_color_manual()}.
#'
#' @return ggplot2 scale layer.
#'
#' @keywords internal
#' @noRd
scale_color_tsenat_discrete <- function(...) {
  ggplot2::scale_color_manual(
    values = .tsenat_palette_discrete(),
    ...
  )
}

#' TSENAT Discrete Fill Scale
#'
#' Wrapper for ggplot2 fill scale using the TSENAT discrete palette.
#'
#' @param ... Additional arguments passed to \code{ggplot2::scale_fill_manual()}.
#'
#' @return ggplot2 scale layer.
#'
#' @keywords internal
#' @noRd
scale_fill_tsenat_discrete <- function(...) {
  ggplot2::scale_fill_manual(
    values = .tsenat_palette_discrete(),
    ...
  )
}

#' TSENAT Diverging Color Scale
#'
#' Wrapper for ggplot2 color scale using the TSENAT diverging palette.
#' Use for continuous diverging data (e.g., log fold-change).
#'
#' @param ... Additional arguments passed to \code{ggplot2::scale_color_gradientn()}.
#'
#' @return ggplot2 scale layer.
#'
#' @keywords internal
#' @noRd
scale_color_tsenat_diverging <- function(...) {
  ggplot2::scale_color_gradientn(
    colors = .tsenat_palette_continuous_diverging(n = 100),
    ...
  )
}

#' TSENAT Diverging Fill Scale
#'
#' Wrapper for ggplot2 fill scale using the TSENAT diverging palette.
#' Use for continuous diverging data (e.g., heatmap values).
#'
#' @param ... Additional arguments passed to \code{ggplot2::scale_fill_gradientn()}.
#'
#' @return ggplot2 scale layer.
#'
#' @keywords internal
#' @noRd
scale_fill_tsenat_diverging <- function(...) {
  ggplot2::scale_fill_gradientn(
    colors = .tsenat_palette_continuous_diverging(n = 100),
    ...
  )
}

# ============================================================================
# RESPONSIVE SIZING SYSTEM - Aspect ratio and dimension calculation
# ============================================================================

#' TSENAT Plot Aspect Ratio Constants
#'
#' Defines standardized aspect ratios (width/height) for different plot types.
#' Used to maintain consistent visual proportions across the TSENAT package.
#'
#' @details
#' Aspect ratios defined:
#' - `.aspect_ratio_standard`: 1.67 (16:10) - Default for most plots
#' - `.aspect_ratio_tall`: 1.20 (6:5) - For plots with many rows (heatmaps, multi-gene panels)
#' - `.aspect_ratio_wide`: 2.40 (12:5) - For time-series or wide categorical plots
#'
#' These constants ensure plots maintain readable font sizes and proportions
#' regardless of absolute output dimensions.
#'
#' @keywords internal
#' @noRd
.aspect_ratio_standard <- 12 / 7.2    # 1.67 (16:10, publication standard)
.aspect_ratio_tall <- 12 / 10         # 1.20 (taller, for heatmaps)
.aspect_ratio_wide <- 12 / 5          # 2.40 (wider, for spectral profiles)

#' Calculate Plot Dimensions from Width and Aspect Ratio
#'
#' Helper function to compute height from width and aspect ratio type.
#' Returns both dimensions and DPI as a list for use with ggsave and grid functions.
#'
#' @param width_inches Numeric; plot width in inches (default: 12).
#' @param aspect_type Character; one of "standard", "tall", or "wide" (default: "standard").
#' @param dpi_output Numeric; output DPI for PNG/TIFF (default: 100).
#'
#' @return List with elements: $width (inches), $height (inches), $dpi (integer).
#'
#' @details
#' Example usage in plot functions:
#'
#' ```r
#' # Standard aspect ratio plot
#' dims <- .calculate_plot_dims(width_inches = 12, aspect_type = "standard")
#' ggplot2::ggsave("plot.pdf", plot_obj, width = dims$width, 
#'                  height = dims$height, dpi = dims$dpi)
#'
#' # Heatmap with many rows (use tall aspect)
#' dims <- .calculate_plot_dims(width_inches = 12, aspect_type = "tall")
#' # ... render and save
#' ```
#'
#' This ensures all plots maintain:
#' - Readable font sizes
#' - Consistent proportions
#' - Professional appearance across all outputs
#'
#' @keywords internal
#' @noRd
.calculate_plot_dims <- function(width_inches = 12, 
                                  aspect_type = "standard",
                                  dpi_output = 100) {
  # Validate aspect_type
  valid_aspects <- c("standard", "tall", "wide")
  if (!aspect_type %in% valid_aspects) {
    warning("aspect_type '", aspect_type, "' not recognized. Using 'standard'.")
    aspect_type <- "standard"
  }
  
  # Get aspect ratio constant by direct lookup
  aspect_ratio <- switch(aspect_type,
    standard = .aspect_ratio_standard,
    tall = .aspect_ratio_tall,
    wide = .aspect_ratio_wide,
    .aspect_ratio_standard  # Fallback to standard
  )
  
  # Calculate height from width and aspect ratio
  # aspect_ratio = width / height, so:
  # height = width / aspect_ratio
  height_inches <- width_inches / aspect_ratio
  
  list(
    width = width_inches,
    height = height_inches,
    dpi = as.integer(dpi_output)
  )
}

#' Scale Font Sizes Responsively Based on Plot Area
#'
#' Adjusts base font size scaling factor based on plot output dimensions.
#' Useful for maintaining readability across different output sizes.
#'
#' @param width_inches Numeric; plot width in inches.
#' @param height_inches Numeric; plot height in inches.
#' @param reference_area Numeric; reference area in square inches for baseline scaling.
#'   Default (96 sq in) represents a 12x8 inch plot.
#'
#' @return Numeric scaling factor to multiply with base font sizes.
#'
#' @details
#' The scaling factor is computed as: sqrt(actual_area / reference_area)
#'
#' This ensures font readability is maintained proportionally with plot size.
#' A plot at 6x4 inches (1/4 the area) gets 0.5x font scaling,
#' while a 24x16 inch plot (4x the area) gets 2x scaling.
#'
#' Example:
#' ```r
#' # For a heatmap that's 15 inches wide and 12 inches tall
#' scale <- .scale_font_by_area(width_inches = 15, height_inches = 12)
#' # Results in scale ~ 1.38, increasing fonts for larger output
#'
#' # Use in custom plotting:
#' font_scale <- .scale_font_by_area(12, 10)
#' base_font <- 11 * font_scale
#' ```
#'
#' @keywords internal
#' @noRd
.scale_font_by_area <- function(width_inches, height_inches, 
                                 reference_area = 96) {
  actual_area <- width_inches * height_inches
  sqrt(actual_area / reference_area)
}

