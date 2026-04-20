#' Regularized/penalized regression interaction testing for Tsallis entropy
#'
#' For each gene, fit a regularized/penalized regression model (GAM, LMM, GEE, or FPCA) 
#' of the form `entropy ~ q * group` with AR(1) correlation structure and extract the 
#' p-value for the interaction term (whether the effect of `q` differs between groups). 
#' The function expects a `SummarizedExperiment` produced by `.calculate_diversity()` 
#' when multiple `q` values have been computed (column names contain `_q=`).
#'
#' @param se A `SummarizedExperiment` containing a `diversity` assay produced
#'   by `.calculate_diversity(..., q = <vector>)`.
#' @param condition_col Column name in `colData(se)` that contains
#'   a grouping factor for samples (character). If `NULL`, the function will
#'   attempt to infer group from column names (suffix `_N` interpreted as
#'   'Normal').
#' @param min_obs Minimum number of non-NA observations required to fit a
#'   model for a gene (default: 10).
#' @param method Modeling method to use for interaction testing: one of
#'   \code{c('lmm', 'gam', 'fpca', 'gee')} (default: 'lmm').
#'   \itemize{
#'     \item `'lmm'`: Linear mixed models with AR(1) covariance for q-ordered
#'       measurements (requires nlme). Recommended for strong signal detection.
#'     \item `'gam'`: Generalized additive models with flexible smoothing
#'       (requires mgcv). Useful for non-monotonic entropy patterns.
#'     \item `'fpca'`: Functional principal components analysis for curve data.
#'       Respects q-value ordering by treating q-values as ordered measurements
#'       in a functional data framework (Papers S168-S171). PCA on ordered curves
#'       implicitly captures AR(1) correlation structure.
#'     \item `'gee'`: Generalized estimating equations for clustered/paired data
#'       (requires geepack). Particularly useful for longitudinal designs with
#'       repeated q-measures.
#'   }
#'   Q-values are mathematically dependent (Papers S168-S175: AR(1) covariance
#'   structures). Regularized regression models without correlation structure should not be used.
#' @param pvalue Type of p-value to compute: one of
#'   \code{c('satterthwaite', 'lrt', 'both')} (default: 'satterthwaite').
#'   Note: For method='lmm', only LRT p-values are available (Satterthwaite
#'   requires lme4 which does not support AR(1) covariance structures).
#'   This parameter is ignored for LMM.
#' @param subject_col Optional column name in `colData(se)` that contains
#'   subject/individual identifiers for paired or repeated-measures designs
#'   (character). If `NULL` and `paired = TRUE`, automatically uses the third
#'   column of `colData(se)`. If provided with `method = 'lmm'`, used as
#'   random effect.
#' @param paired Logical; whether samples are paired (default: FALSE).
#' @param nthreads Number of threads (mc.cores) to use for parallel processing
#'   (default: 1).
#' @param assay_name Name of the assay in the SummarizedExperiment to use
#'   (default: 'diversity').
#' @param pcorr P-value correction method applied to Wilcoxon rank test results
#'   (default: 'BH'). Options: \code{c('BH', 'bonferroni', 'hochberg', 'holm')}.
#'   Note: This is distinct from `multicorr` which adjusts for correlation
#'   across multiple q-values in the interaction test. `pcorr` is legacy and
#'   may not be used in all methods. See `multicorr` for the primary multiple
#'   testing correction.
#' @param verbose Logical; whether to print progress messages during execution
#'   (default: FALSE).
#' @param corstr Correlation structure for GEE method: one of
#'   \code{c('ar1', 'exchangeable', 'independence')} (default: 'ar1').
#'   \itemize{
#'     \item `'ar1'`: First-order autoregressive (recommended for q-ordered
#'       measurements).
#'     \item `'exchangeable'`: Compound symmetry (appropriate when measurements
#'       are unordered).
#'     \item `'independence'`: Independent observations (no correlation).
#'   }
#'   Paper S171: Zimmerman & Harville (1991) validates AR(1) for ordered data.
#'   This parameter only affects `method='gee'`.
#' @param bias_correction Logical; whether to apply Kauermann-Carroll (K-C)
#'   bias correction for GEE with small number of clusters (default: TRUE).
#'   When TRUE and the number of clusters is less than 20, uses t-distribution
#'   instead of normal distribution for p-value computation, which maintains
#'   Type I error rate for small sample GEE analyses. Reference:
#'   Li & Redden (2015), Statistics in Medicine. This parameter only affects
#'   `method='gee'`.
#' @param regularization Dimensionality reduction method for FPCA analysis:
#'   one of \code{c('pca', 'lasso', 'elasticnet')} (default: 'pca').
#'   \itemize{
#'     \item `'pca'`: Classical PCA extracting the first principal component (PC1).
#'     \item `'lasso'`: L1-penalized logistic regression with cross-validation to
#'       select important q-values (References: Friedman et al. 2010, Bloch 2020).
#'     \item `'elasticnet'`: Elastic net (L1+L2 penalty, alpha=0.5) for improved
#'       stability (References: Friedman et al. 2010, Byliskii 2015).
#'   }
#'   Regularization methods can provide higher statistical power than PC1 by
#'   automatically selecting informative q-values. This parameter only affects
#'   `method='fpca'`.
#' @param multicorr Method for adjusting p-values across multiple q-values
#'   to account for correlation structure in Tsallis entropy (default: 'hochberg').
#'   The interaction p-values from Scale-Adaptive Interaction Models naturally exhibit AR(1)
#'   correlation for different q-values of the same gene (Papers S168-S175).
#'   \itemize{
#'     \item `'hochberg'`: Hochberg stepup procedure (FWER <= alpha under positive
#'       regression dependence). Closed-form, computationally efficient.
#'       Recommended for strong signal detection with family-wise error control.
#'     \item `'westfall-young'`: True Westfall-Young permutation procedure
#'       (FWER <= alpha). Uses resampling to empirically control FWER by
#'       tracking the minima across all tests. More powerful than Hochberg under
#'       dependence but computationally expensive (refits LMM for each permutation).
#'     \item `'benjamini-yekutieli'`: Benjamini-Yekutieli FDR control
#'       (FDR <= alpha under arbitrary dependence). Valid under any correlation
#'       structure. More conservative than Hochberg but makes fewer power loss
#'       assumptions. Reference: Papers S190, S193.
#'   }
#' @param wy_randomizations Number of permutation randomizations for
#'   Westfall-Young correction (default: 1000). Only used when
#'   `multicorr = 'westfall-young'`. Higher values improve accuracy of
#'   empirical null distribution but increase computation time.
#'   Minimum: 100. Typical values: 500-2000. Note: Westfall-Young is
#'   computationally expensive as it requires refitting models for each
#'   randomization.
#' @param storey Logical; whether to apply Storey's adaptive FDR π0
#'   estimation after the selected multicorr method (default: FALSE).
#'   When TRUE, adapts the error threshold based on estimated proportion of
#'   true null hypotheses, increasing power when many true signals are present.
#'   Can be applied to any multicorr method. Computationally light enhancement.
#'   Requires: .estimate_storey_pi0() and .compute_storey_qvalues() functions.
#'   Reference: Storey (2002).
#' @param adaptive_knots Logical; whether to use adaptive spline knot
#'   selection for GAM method (default: TRUE). When TRUE, automatically adjusts
#'   the number of basis functions (k) per gene based on entropy curve complexity,
#'   measured as coefficient of variation of slopes across ordered q-values.
#'   Simple curves get fewer knots (min=2), complex non-monotonic curves get
#'   more knots (max=10), improving model fit efficiency. When FALSE, uses
#'   fixed knot selection based on number of unique q-values. Reference:
#'   Wood (2017) Section 4.1.5 Basis dimension. This parameter only affects
#'   `method='gam'`.
#' @param return_model_data Logical; whether to return model metadata
#'   alongside results (default: FALSE). When TRUE, returns a list with two
#'   elements:
#'   \itemize{
#'     \item `$results`: The standard results data.frame (same as returned
#'       when FALSE)
#'     \item `$model_data`: A list containing model metadata (method, q-values,
#'       sample info, test configuration, genes analyzed, etc.) useful for
#'       generating diagnostic plots and understanding model structure. Can be
#'       passed to plotting functions for visualization.
#'   }
#'   When FALSE, returns only the results data.frame (backward compatible with
#'   existing code). This enables users to access comprehensive model information
#'   for diagnostics and visualization while maintaining full backward
#'   compatibility.
#' @return When `return_model_data = FALSE` (default): A data.frame with
#'   columns `gene`, `p_interaction`, and `adj_p_interaction`, ordered by
#'   ascending `p_interaction`.
#'
#'   When `return_model_data = TRUE`: A list with components:
#'   \itemize{
#'     \item `$results`: The standard results data.frame
#'     \item `$model_data`: Metadata list containing method, q-values,
#'       sample names, test configuration, and other information useful for
#'       downstream visualization and diagnostics.
#'   }
#'
#' @section Sample Metadata Parameters (Unified Naming Convention):
#'   TSENAT functions use consistent parameter names for sample grouping and
#'   subject identification:
#'   \itemize{
#'     \item{\code{condition_col}: Character string specifying the colData
#'       column containing sample group/condition labels (e.g., 'Normal',
#'       'Tumor', 'control', 'treatment'). Default: 'condition'. Map your
#'       grouping variable into this column before calling TSENAT functions.}
#'     \item{\code{subject_col}: For paired/blocked/repeated-measures designs,
#'       character string specifying the colData column with
#'       subject/individual/patient identifiers. Default: NULL.
#'       Required when `paired = TRUE`.}
#'   }
#'   All functions use `SummarizedExperiment::colData()` as the single source
#'   of truth for sample metadata. This eliminates parameter fragmentation and
#'   improves API discoverability across the TSENAT package.
#'
#' @references
#'   Kutner, M. H., Nachtsheim, C. J., Neter, J., & Li, W. (2005).
#'   \emph{Applied Linear Statistical Models} (5th ed.). McGraw-Hill.
#'   Comprehensive treatment of linear regression methodology including
#'   interaction models. Chapter 8.2.7 covers testing hypotheses in multiple
#'   linear regression. NOTE: Linear method removed (see `method` parameter
#'   documentation).
#'
#'   Wood, S. N. (2017). \emph{Generalized Additive Models: An Introduction
#'   with R} (2nd ed.). Chapman & Hall/CRC. Comprehensive treatment of GAM and
#'   GAMM methodology, smooth basis selection, and model comparison for
#'   `method='gam'`.
#'
#'   Liang, K. Y., & Zeger, S. L. (1986). Longitudinal data analysis using
#'   generalized linear models. \emph{Biometrika}, 73(1), 13-22. Foundational
#'   paper introducing Generalized Estimating Equations (GEE) for clustered and
#'   repeated measurement data. Primary reference for `method='gee'`.
#'
#'   Song, P. X. K. (2007). \emph{Correlated Data Analysis: Modeling, Analytics,
#'   and Applications}. Springer. Comprehensive treatment of GEE, MEASURE models,
#'   and advanced methods for handling correlated data structures common in
#'   longitudinal and spatial studies. Extended methodology reference for
#'   `method='gee'`.
#'
#'   Benjamini, Y., & Hochberg, Y. (1995). Controlling the false discovery rate:
#'   A practical and powerful approach to multiple testing. \emph{Journal of
#'   the Royal Statistical Society}, Series B, 57, 289-300. Used for multiple
#'   testing correction via `p.adjust(..., method='BH')`.
#'
#'   Hochberg, Y. (1988). A sharper Bonferroni procedure for multiple tests of
#'   significance. \emph{Biometrika}, 75(4), 800-802. Stepup procedure for FWER
#'   control used in `multicorr='westfall-young'` option. More powerful than
#'   Bonferroni.
#'
#'   Benjamini, Y., & Yekutieli, D. (2001). The control of the false discovery
#'   rate in multiple testing under dependency. \emph{Annals of Statistics},
#'   29(4), 1165-1188. FDR control under arbitrary dependence (Papers S190, S193).
#'   Used in `multicorr='benjamini-yekutieli'` option.
#'
#'   Storey, J. D. (2002). A direct approach to false discovery rates.
#'   \emph{Journal of the Royal Statistical Society}, Series B, 64(3),
#'   479-498. Adaptive FDR estimation via π0 proportion (used in
#'   `multicorr='westfall-young-storey'`). More powerful than Hochberg when
#'   substantial proportion of nulls are true.
#'
#' @examples
#'   # Create example data
#'   set.seed(123)
#'
#'   # Simulate read counts: 5 genes, 3 transcripts each, 4 samples
#'   counts <- matrix(
#'     sample(1:100, 60, replace = TRUE),
#'     nrow = 15, ncol = 4
#'   )
#'   rownames(counts) <- paste0('tx_', 1:15)
#'   colnames(counts) <- paste0('sample_', 1:4)
#'   genes <- rep(paste0('gene_', 1:5), each = 3)
#'
#'   # Calculate diversity at multiple q values
#'   se <- .calculate_diversity(
#'     counts, genes = genes,
#'     q = c(0.5, 1.0, 1.5, 2.0, 2.5),
#'     norm = TRUE
#'   )
#'
#'   # Add sample metadata
#'   SummarizedExperiment::colData(se) <- S4Vectors::DataFrame(
#'     condition = rep(c('Normal', 'Tumor'), length.out = ncol(se)),
#'     row.names = colnames(se)
#'   )
#'
#'   # Run regularized regression interaction analysis
#'   results <- TSENAT:::.calculate_sait(se, condition_col = 'condition')
#' @noRd
.calculate_sait <- function(se, condition_col = "condition", min_obs = 5, method = c("lmm",
    "gam", "fpca", "gee"), pvalue = c("satterthwaite", "lrt", "both"), subject_col = NULL,
    paired = FALSE, nthreads = 1, assay_name = "diversity", pcorr = "BH", verbose = FALSE,
    bias_correction = TRUE, regularization = c("pca", "lasso", "elasticnet", "gamsel",
        "spline"), corstr = c("ar1", "exchangeable", "independence"), multicorr = c("hochberg",
        "westfall-young", "benjamini-yekutieli"), storey = FALSE, wy_randomizations = 1000,
    adaptive_knots = TRUE, return_model_data = FALSE) {
    
    # ========================================================================
    # STAGE 1: ARGUMENT NORMALIZATION & BASIC VALIDATION
    # ========================================================================
    
    method <- match.arg(method)
    corstr <- match.arg(corstr)
    pvalue <- match.arg(pvalue)
    regularization <- match.arg(regularization)
    pcorr <- match.arg(pcorr, c("BH", "bonferroni", "hochberg", "holm"))
    multicorr <- match.arg(multicorr)
    
    if (!requireNamespace("SummarizedExperiment", quietly = TRUE)) {
        stop("SummarizedExperiment required")
    }
    
    if (!inherits(se, "SummarizedExperiment")) {
        stop("se must be a SummarizedExperiment object", call. = FALSE)
    }
    
    if (verbose) {
        message("[calculate_sait_interaction] method=", method)
    }
    
    # ========================================================================
    # STAGE 2: DEPENDENCY & STRUCTURE VALIDATION
    # ========================================================================
    
    .validate_sait_method_dependencies(method)
    
    validated <- .validate_sait_interaction_input(method = method, pvalue = pvalue,
        corstr = corstr, regularization = regularization, multicorr = multicorr,
        pcorr = pcorr, storey = storey, wy_randomizations = wy_randomizations, 
        paired = paired, subject_col = subject_col, se = se, verbose = verbose)
    
    subject_col <- validated$subject_col
    
    .validate_sait_data_structure(se, condition_col, assay_name)
    
    # ========================================================================
    # STAGE 3: DATA PREPARATION & FITTING
    # ========================================================================
    
    metadata <- .parse_sample_metadata(se = se, condition_col = condition_col, 
        assay_name = assay_name, verbose = verbose)
    
    mat <- SummarizedExperiment::assay(se, assay_name)
    
    if (verbose)
        message("[.calculate_sait] Starting .fit_all_genes() for ", nrow(mat), " genes")
    
    # Wrap fitting in try-error to catch any errors during fitting
    res <- try(.fit_all_genes(mat = mat, se = se, metadata = metadata, method = method,
        pvalue = pvalue, subject_col = subject_col, paired = paired, min_obs = min_obs,
        nthreads = nthreads, verbose = verbose, bias_correction = bias_correction,
        regularization = regularization, corstr = corstr, adaptive_knots = adaptive_knots),
        silent = FALSE)
    
    if (inherits(res, "try-error")) {
        error_msg <- if (!is.null(attr(res, "condition"))) {
            conditionMessage(attr(res, "condition"))
        } else {
            as.character(res)
        }
        warning("[.calculate_sait] .fit_all_genes() failed with: ", error_msg, 
            "\n[Returning empty results]", call. = FALSE)
        res <- data.frame()
    }
    
    if (verbose && nrow(res) > 0)
        message("[.calculate_sait] .fit_all_genes() completed successfully with ",
            nrow(res), " results")
    
    # ========================================================================
    # STAGE 4: RESULTS PROCESSING & RETURN
    # ========================================================================
    
    if (!is.data.frame(res)) {
        stop(".fit_all_genes() should return a data.frame", call. = FALSE)
    }
    
    if (nrow(res) == 0) {
        return(res)
    }
    
    if (!("p_interaction" %in% colnames(res))) {
        stop("Results data.frame missing required 'p_interaction' column", call. = FALSE)
    }
    
    # Adjust p-values for multiple q-values
    res$adj_p_interaction <- .adjust_pvalues_multicorr(p_values = res$p_interaction,
        multicorr = multicorr, wy_randomizations = wy_randomizations, metadata = metadata,
        verbose = verbose, storey = storey)
    
    # Sort by adjusted p-values, then raw p-values
    res <- res[order(res$adj_p_interaction, res$p_interaction), , drop = FALSE]
    rownames(res) <- NULL
    
    .report_fit_summary(res, verbose = verbose)
    
    # Map gene identifiers to annotations
    res <- .map_gene_annotations(res = res, se = se, verbose = verbose)
    
    # Post-process: add gene column, optionally add model data
    .postprocess_sait_results(res = res, return_model_data = return_model_data, 
        se = se, mat = mat, metadata = metadata, method = method,
        pvalue = pvalue, multicorr = multicorr, assay_name = assay_name,
        bias_correction = bias_correction, regularization = regularization, 
        corstr = corstr, adaptive_knots = adaptive_knots)
}

# Validate method-specific dependencies (cyclomatic complexity reducer)
.validate_sait_method_dependencies <- function(method) {
    if (method == "lmm" && !requireNamespace("nlme", quietly = TRUE)) {
        stop("Package 'nlme' is required for method='lmm'", call. = FALSE)
    }
    if (method == "gam" && !requireNamespace("mgcv", quietly = TRUE)) {
        stop("Package 'mgcv' is required for method='gam'", call. = FALSE)
    }
    if (method == "gee" && !requireNamespace("geepack", quietly = TRUE)) {
        stop("Package 'geepack' is required for method='gee'", call. = FALSE)
    }
    invisible(TRUE)
}

# Validate required columns and assays (cyclomatic complexity reducer)
# REFACTORING: Extract validation into separate function
.validate_sait_data_structure <- function(se, condition_col, assay_name) {
    cd_colnames <- colnames(SummarizedExperiment::colData(se))
    if (!(condition_col %in% cd_colnames)) {
        stop(sprintf("condition_col '%s' not found in colData. Available columns: %s",
            condition_col, paste(cd_colnames, collapse = ", ")), call. = FALSE)
    }
    
    if (!(assay_name %in% SummarizedExperiment::assayNames(se))) {
        stop(sprintf("Assay '%s' not found. Available assays: %s", assay_name, 
            paste(SummarizedExperiment::assayNames(se), collapse = ", ")), call. = FALSE)
    }
    invisible(TRUE)
}

# Post-process SAIT results (cyclomatic complexity reducer)
# REFACTORING: Extract post-processing into separate function
.postprocess_sait_results <- function(res, return_model_data = FALSE, se = NULL, 
                                   mat = NULL, metadata = NULL, method = NULL, 
                                   pvalue = NULL, multicorr = NULL, assay_name = NULL,
                                   bias_correction = NULL, regularization = NULL, 
                                   corstr = NULL, adaptive_knots = NULL) {
    # Ensure 'gene' column exists - required by downstream functions
    if (!("gene" %in% colnames(res))) {
        if ("gene_id" %in% colnames(res)) {
            res$gene <- res$gene_id
        } else if ("gene_name" %in% colnames(res)) {
            res$gene <- res$gene_name
        }
    }
    
    # Optionally return model data alongside results
    if (return_model_data) {
        model_data <- .assemble_model_metadata(se = se, res = res, mat = mat, 
            metadata = metadata, method = method, pvalue = pvalue, 
            multicorr = multicorr, assay_name = assay_name,
            bias_correction = bias_correction, regularization = regularization, 
            corstr = corstr, adaptive_knots = adaptive_knots)
        return(list(results = res, model_data = model_data))
    }
    
    return(res)
}



#' Extract Results from calculate_sait Output
#'
#' Helper function to extract results data.frame from
#' calculate_sait output,
#' which may be either a data.frame (when return_model_data=FALSE) or a list 
#' (when return_model_data=TRUE). This ensures compatibility with plotting and 
#' analysis functions regardless of return format.
#'
#' @param sait_result Result from .calculate_sait(), either a
#' data.frame or a list
#'
#' @return The results data.frame with columns gene, p_interaction,
#' adj_p_interaction, etc.
#'

#' @noRd
.extract_sait_result_df <- function(sait_result) {
    if (is.data.frame(sait_result)) {
        return(sait_result)
    } else if (is.list(sait_result) && "results" %in% names(sait_result)) {
        return(sait_result$results)
    } else {
        stop("sait_result must be either a data.frame or a list with 'results' component from .calculate_sait()")
    }
}
