#' Linear-model interaction test for Tsallis entropy   For each gene, fit a
#' linear model of the form `entropy ~ q * group` and  extract the p-value for
#' the interaction term (whether the effect of `q`  differs between groups). The
#' function expects a `SummarizedExperiment`  produced by
#' `.calculate_diversity()` when multiple `q` values have been  computed (column
#' names contain `_q=`).
#' @param se A `SummarizedExperiment` containing a `diversity` assay produced
#' by `.calculate_diversity(..., q = <vector>)`.
#' @param condition_col Column name in `colData(se)` that contains
#' a grouping factor for samples (character). If `NULL`, the function will
#' attempt to infer group from column names (suffix `_N` interpreted as
#' 'Normal').
#' @param min_obs Minimum number of non-NA observations required to fit a
#' model for a gene (default: 10).
#' @param method Modeling method to use for interaction testing: one of
#' \code{c('lmm', 'gam', 'fpca', 'gee')} (default: 'lmm').
#' DEPRECATED: 'linear' (OLS) has been removed because it treats q-values as independent,
#' which violates Tsallis entropy properties. Q-values are mathematically dependent
#' (Papers S168-S175: AR(1) covariance structures). Use 'lmm' instead.
#' 'lmm': linear mixed models with AR(1) covariance for q-ordered measurements (requires nlme).
#' 'gam': generalized additive models with flexible smoothing (requires mgcv).
#' 'fpca': functional principal components analysis for curve data. Respects q-value ordering
#' by treating q-values as ordered measurements in a functional data framework (Papers S168-S171).
#' PCA on ordered curves implicitly captures AR(1) correlation structure.
#' 'gee': generalized estimating equations for clustered/paired data (requires geepack).
#' GEE is particularly useful for longitudinal designs with repeated q-measures.
#' @param pvalue Type of p-value to compute: one of
#' \code{c('satterthwaite', 'lrt', 'both')} (default: 'satterthwaite').
#' Note: For method='lmm', only LRT p-values are available (Satterthwaite requires lme4
#' which does not support AR(1) covariance structures). This parameter is ignored for LMM.
#' @param subject_col Optional column name in `colData(se)` that contains
#' subject/individual identifiers for paired or repeated-measures designs
#' (character). If `NULL` and `paired = TRUE`, automatically uses the third 
#' column of `colData(se)`. If provided with `method = 'lmm'`, used as random effect.
#' @param paired Logical; whether samples are paired (default: FALSE).
#' @param nthreads Number of threads (mc.cores) to use for parallel processing
#' (default: 1).
#' @param assay_name Name of the assay in the SummarizedExperiment to use
#' (default: 'diversity').
#' @param pcorr P-value correction method applied to Wilcoxon rank test results 
#' (default: 'BH'). Options: \code{c('BH', 'bonferroni', 'hochberg', 'holm')}.
#' Note: This is distinct from `multicorr` which adjusts for correlation across 
#' multiple q-values in the interaction test. `pcorr` is legacy and may not be 
#' used in all methods. See `multicorr` for the primary multiple testing correction.
#' @param verbose Logical; whether to print progress messages during execution
#' (default: FALSE).
#' @param corstr Correlation structure for GEE method: one of
#' \code{c('ar1', 'exchangeable', 'independence')} (default: 'ar1').
#' 'ar1': First-order autoregressive (recommended for q-ordered measurements).
#' 'exchangeable': Compound symmetry (appropriate when measurements are unordered).
#' 'independence': Independent observations (no correlation structure).
#' Paper S171: Zimmerman & Harville (1991) validates AR(1) for ordered data.
#' This parameter only affects method='gee'.
#' @param bias_correction Logical; whether to apply Kauermann-Carroll (K-C) bias correction
#' for GEE with small number of clusters (default: TRUE). When TRUE and the number of
#' clusters is less than 20, uses t-distribution instead of normal distribution for p-value
#' computation, which maintains Type I error rate for small sample GEE analyses. Reference: 
#' Li & Redden (2015), Statistics in Medicine. This parameter only affects method='gee'.
#' @param regularization Dimensionality reduction method for FPCA analysis: one of
#' \code{c('pca', 'lasso', 'elasticnet')} (default: 'pca'). 'pca' uses classical PCA
#' extracting the first principal component (PC1). 'lasso' uses L1-penalized logistic
#' regression with cross-validation to select important q-values (References: Friedman et al. 2010,
#' Bloch 2020). 'elasticnet' uses elastic net (L1+L2 penalty, alpha=0.5) for improved stability
#' (References: Friedman et al. 2010, Byliskii 2015). Regularization methods can provide
#' higher statistical power than PC1 by automatically selecting informative q-values. This
#' parameter only affects method='fpca'.
#' @param multicorr Method for adjusting p-values across multiple q-values to account for 
#' correlation structure in Tsallis entropy (default: 'hochberg'). The interaction 
#' p-values from linear models naturally exhibit AR(1) correlation for different q-values 
#' of the same gene (Papers S168-S175). This parameter selects the primary multiple testing
#' correction method:
#' 'hochberg': Hochberg stepup procedure (FWER <= alpha under positive regression dependence). 
#' Closed-form, computationally efficient. Recommended for strong signal detection with 
#' family-wise error control.
#' 'westfall-young': True Westfall-Young permutation procedure (FWER <= alpha). Uses resampling 
#' to empirically control FWER by tracking the minima across all tests. More powerful than 
#' Hochberg under dependence but computationally expensive (refits LMM for each permutation).
#' 'benjamini-yekutieli': Benjamini-Yekutieli FDR control (FDR <= alpha under arbitrary dependence). 
#' Valid under any correlation structure. More conservative than Hochberg but makes fewer 
#' power loss assumptions. Reference: Papers S190, S193.
#' @param wy_randomizations Number of permutation randomizations for Westfall-Young 
#' correction (default: 1000). Only used when `multicorr = 'westfall-young'`. 
#' Higher values improve accuracy of empirical null distribution but increase computation time.
#' Minimum: 100. Typical values: 500-2000. Note: Westfall-Young is computationally expensive 
#' as it requires refitting models for each randomization.
#' @param storey Logical; whether to apply Storey's adaptive FDR ?0 estimation after the 
#' selected multicorr method (default: FALSE). When TRUE, adapts the error threshold based 
#' on estimated proportion of true null hypotheses, increasing power when many true signals 
#' are present. Can be applied to any multicorr method. Computationally light enhancement.
#' Requires: .estimate_storey_pi0() and .compute_storey_qvalues() functions. Reference: Storey (2002).
#' @param adaptive_knots Logical; whether to use adaptive spline knot selection for GAM method
#' (default: TRUE). When TRUE, automatically adjusts the number of basis functions (k) per gene
#' based on entropy curve complexity, measured as coefficient of variation of slopes across
#' ordered q-values. Simple curves get fewer knots (min=2), complex non-monotonic curves get
#' more knots (max=10), improving model fit efficiency. When FALSE, uses fixed knot selection
#' based on number of unique q-values. Reference: Wood (2017) Section 4.1.5 Basis dimension.
#' This parameter only affects method='gam'.
#' @param return_model_data Logical; whether to return model metadata alongside results
#' (default: FALSE). When TRUE, returns a list with two elements:
#' \itemize{
#'   \item `$results`: The standard results data.frame (same as returned when FALSE)
#'   \item `$model_data`: A list containing model metadata (method, q-values, sample info, 
#'     test configuration, genes analyzed, etc.) useful for generating diagnostic plots and 
#'     understanding model structure. Can be passed to plotting functions for visualization.
#' }
#' When FALSE, returns only the results data.frame (backward compatible with existing code).
#' This enables users to access comprehensive model information for diagnostics and visualization
#' while maintaining full backward compatibility.
#' @return When `return_model_data = FALSE` (default): A data.frame with columns `gene`, 
#' `p_interaction`, and `adj_p_interaction`, ordered by ascending `p_interaction`.
#' 
#' When `return_model_data = TRUE`: A list with components:
#' \itemize{
#'   \item `$results`: The standard results data.frame
#'   \item `$model_data`: Metadata list containing method, q-values, sample names, test configuration,
#'     and other information useful for downstream visualization and diagnostics
#' }
#'
#' @section Sample Metadata Parameters (Unified Naming Convention):
#' TSENAT functions use consistent parameter names for sample grouping and subject identification:
#' \itemize{
#'   \item{\code{condition_col}: Character string specifying the colData column 
#'         containing sample group/condition labels (e.g., "Normal", "Tumor", "control", "treatment"). 
#'         Default: "condition". Map your grouping variable into this column before calling TSENAT functions.}
#'   \item{\code{subject_col}: For paired/blocked/repeated-measures designs, character string specifying 
#'         the colData column with subject/individual/patient identifiers. Default: NULL. 
#'         Required when \code{paired = TRUE}.}
#' }
#' All functions use \code{SummarizedExperiment::colData()} as the single source of truth 
#' for sample metadata. This eliminates parameter fragmentation and improves API discoverability 
#' across the TSENAT package.
#'
#' @references
#' Kutner, M. H., Nachtsheim, C. J., Neter, J., & Li, W. (2005).
#' \emph{Applied Linear Statistical Models} (5th ed.). McGraw-Hill.
#' Comprehensive treatment of linear regression methodology including interaction
#' models. Chapter 8.2.7 covers testing hypotheses in multiple linear regression.
#' NOTE: Linear method removed (see \code{method} parameter documentation).
#' 
#' Wood, S. N. (2017). \emph{Generalized Additive Models: An Introduction with R}
#' (2nd ed.). Chapman & Hall/CRC. Comprehensive treatment of GAM and GAMM
#' methodology, smooth basis selection, and model comparison for \code{method='gam'}.
#' 
#' Liang, K. Y., & Zeger, S. L. (1986). Longitudinal data analysis using 
#' generalized linear models. \emph{Biometrika}, 73(1), 13-22. Foundational paper
#' introducing Generalized Estimating Equations (GEE) for clustered and repeated 
#' measurement data. Primary reference for \code{method='gee'}.
#' 
#' Song, P. X. K. (2007). \emph{Correlated Data Analysis: Modeling, Analytics, 
#' and Applications}. Springer. Comprehensive treatment of GEE, MEANSURE models,
#' and advanced methods for handling correlated data structures common in longitudinal
#' and spatial studies. Extended methodology reference for \code{method='gee'}.
#' 
#' Benjamini, Y., & Hochberg, Y. (1995). Controlling the false discovery rate:
#' A practical and powerful approach to multiple testing. \emph{Journal of the Royal
#' Statistical Society}, Series B, 57, 289-300. Used for multiple testing correction
#' via \code{p.adjust(..., method='BH')}.
#' 
#' Hochberg, Y. (1988). A sharper Bonferroni procedure for multiple tests of 
#' significance. \emph{Biometrika}, 75(4), 800-802. Stepup procedure for FWER control
#' used in multicorr='westfall-young' option. More powerful than Bonferroni.
#' 
#' Benjamini, Y., & Yekutieli, D. (2001). The control of the false discovery rate 
#' in multiple testing under dependency. \emph{Annals of Statistics}, 29(4), 1165-1188.
#' FDR control under arbitrary dependence (Papers S190, S193). Used in 
#' multicorr='benjamini-yekutieli' option.
#' 
#' Storey, J. D. (2002). A direct approach to false discovery rates. 
#' \emph{Journal of the Royal Statistical Society}, Series B, 64(3), 479-498. 
#' Adaptive FDR estimation via λ0 proportion (used in multicorr='westfall-young-storey').
#' More powerful than Hochberg when substantial proportion of nulls are true.

#' @noRd
#' @examples
#' # Create example data
#' set.seed(123)
#' # Simulate read counts: 5 genes, 3 transcripts each, 4 samples
#' counts <- matrix(
#'   sample(1:100, 60, replace = TRUE),
#'   nrow = 15, ncol = 4
#' )
#' rownames(counts) <- paste0("tx_", 1:15)
#' colnames(counts) <- paste0("sample_", 1:4)
#' genes <- rep(paste0("gene_", 1:5), each = 3)
#' 
#' # Calculate diversity at multiple q values
#' se <- .calculate_diversity(counts, genes = genes, q = c(0.5, 1.0, 1.5), norm = TRUE)
#' 
#' # Add sample metadata
#' SummarizedExperiment::colData(se) <- S4Vectors::DataFrame(
#'   condition = rep(c('Normal', 'Tumor'), length.out = ncol(se)),
#'   row.names = colnames(se)
#' )
#' 
#' # Run linear model interaction analysis
#' results <- .calculate_lm_interaction(se, condition_col = "condition")
.calculate_lm_interaction <- function(
    se,
    condition_col = "condition",
    min_obs = 10,
    method = c("lmm", "gam", "fpca", "gee"),
    pvalue = c("satterthwaite", "lrt", "both"),
    subject_col = NULL,
    paired = FALSE,
    nthreads = 1,
    assay_name = "diversity",
    pcorr = "BH",
    verbose = FALSE,
    bias_correction = TRUE,
    regularization = c("pca", "lasso", "elasticnet", "gamsel",
                       "spline"),
    corstr = c("ar1", "exchangeable", "independence"),
    multicorr = c("hochberg", "westfall-young",
                  "benjamini-yekutieli"),
    storey = FALSE,
    wy_randomizations = 1000,
    adaptive_knots = TRUE,
    return_model_data = FALSE
) {
    # Normalize and validate arguments
    method <- match.arg(method)
    corstr <- match.arg(corstr)
    pvalue <- match.arg(pvalue)
    regularization <- match.arg(regularization)
    pcorr <- match.arg(pcorr, c("BH", "bonferroni", "hochberg",
                                "holm"))
    multicorr <- match.arg(multicorr)

    if (!requireNamespace("SummarizedExperiment", quietly = TRUE)) {
        stop("SummarizedExperiment required")
    }

    if (verbose) {
        message("[calculate_lm_interaction] method=", method)
    }

    # Validate input parameters
    validated <- .validate_lm_interaction_input(
        method = method,
        pvalue = pvalue,
        corstr = corstr,
        regularization = regularization,
        multicorr = multicorr,
        pcorr = pcorr,
        storey = storey,
        wy_randomizations = wy_randomizations,
        paired = paired,
        subject_col = subject_col,
        se = se,
        verbose = verbose
    )

    # Update subject_col from validated params (may be auto-detected)
    subject_col <- validated$subject_col

    # Parse sample metadata and q-values
    metadata <- .parse_sample_metadata(
        se = se,
        condition_col = condition_col,
        assay_name = assay_name,
        verbose = verbose
    )

    mat <- SummarizedExperiment::assay(se, assay_name)

    # Fit models to all genes
    res <- .fit_all_genes(
        mat = mat,
        se = se,
        metadata = metadata,
        method = method,
        pvalue = pvalue,
        subject_col = subject_col,
        paired = paired,
        min_obs = min_obs,
        nthreads = nthreads,
        verbose = verbose,
        bias_correction = bias_correction,
        regularization = regularization,
        corstr = corstr,
        adaptive_knots = adaptive_knots
    )

    if (nrow(res) == 0) {
        return(res)
    }

    # Adjust p-values for multiple q-values
    res$adj_p_interaction <- .adjust_pvalues_multicorr(
        p_values = res$p_interaction,
        multicorr = multicorr,
        wy_randomizations = wy_randomizations,
        metadata = metadata,
        verbose = verbose,
        storey = storey
    )

    # Sort by adjusted p-values, then raw p-values
    res <- res[order(res$adj_p_interaction, res$p_interaction),
               , drop = FALSE]
    rownames(res) <- NULL

    .report_fit_summary(res, verbose = verbose)

    # Map gene identifiers to annotations
    res <- .map_gene_annotations(
        res = res,
        se = se,
        verbose = verbose
    )

    # Optionally return model data alongside results
    if (return_model_data) {
        model_data <- .assemble_model_metadata(
            se = se,
            res = res,
            mat = mat,
            metadata = metadata,
            method = method,
            pvalue = pvalue,
            multicorr = multicorr,
            assay_name = assay_name,
            bias_correction = bias_correction,
            regularization = regularization,
            corstr = corstr,
            adaptive_knots = adaptive_knots
        )

        return(list(
            results = res,
            model_data = model_data
        ))
    }

    return(res)
}

#' Extract Results from calculate_lm_interaction Output
#'
#' Helper function to extract results data.frame from calculate_lm_interaction output,
#' which may be either a data.frame (when return_model_data=FALSE) or a list 
#' (when return_model_data=TRUE). This ensures compatibility with plotting and 
#' analysis functions regardless of return format.
#'
#' @param lm_result Result from .calculate_lm_interaction(), either a data.frame or a list
#'
#' @return The results data.frame with columns gene, p_interaction, adj_p_interaction, etc.
#'

#' @noRd
.extract_lm_results <- function(lm_result) {
    if (is.data.frame(lm_result)) {
        return(lm_result)
    } else if (is.list(lm_result) && "results" %in% names(lm_result)) {
        return(lm_result$results)
    } else {
        stop("lm_result must be either a data.frame or a list with 'results' component from .calculate_lm_interaction()")
    }
}