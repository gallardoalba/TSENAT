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

