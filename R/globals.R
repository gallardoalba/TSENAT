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

