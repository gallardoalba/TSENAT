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
#' @export
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
#' @export
METHOD_ORDER <- c("diversity", "jackknife", "lm_interaction", "divergence", "q_interactions")

