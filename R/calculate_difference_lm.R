#' Linear-model interaction test for Tsallis entropy   For each gene, fit a
#' linear model of the form `entropy ~ q * group` and  extract the p-value for
#' the interaction term (whether the effect of `q`  differs between groups). The
#' function expects a `SummarizedExperiment`  produced by
#' `calculate_diversity()` when multiple `q` values have been  computed (column
#' names contain `_q=`).
#' @param se A `SummarizedExperiment` containing a `diversity` assay produced
#' by `calculate_diversity(..., q = <vector>)`.
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
#' Requires: estimate_storey_pi0() and compute_storey_qvalues() functions. Reference: Storey (2002).
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
#' @keywords internal
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
#' se <- calculate_diversity(counts, genes = genes, q = c(0.5, 1.0, 1.5), norm = TRUE)
#' 
#' # Add sample metadata
#' SummarizedExperiment::colData(se) <- S4Vectors::DataFrame(
#'   condition = rep(c('Normal', 'Tumor'), length.out = ncol(se)),
#'   row.names = colnames(se)
#' )
#' 
#' # Run linear model interaction analysis
#' results <- calculate_lm_interaction(se, condition_col = "condition")
calculate_lm_interaction <- function(se, condition_col = "condition", min_obs = 10, method = c("lmm",
    "gam", "fpca", "gee"), pvalue = c("satterthwaite", "lrt", "both"), subject_col = NULL,
    paired = FALSE, nthreads = 1, assay_name = "diversity", pcorr = "BH", verbose = FALSE, 
    bias_correction = TRUE, regularization = c("pca", "lasso", "elasticnet", "gamsel", "spline"),
    corstr = c("ar1", "exchangeable", "independence"), multicorr = c("hochberg", "westfall-young", "benjamini-yekutieli"),
    storey = FALSE, wy_randomizations = 1000, adaptive_knots = TRUE, return_model_data = FALSE) {
    method <- match.arg(method)
    corstr <- match.arg(corstr)
    pvalue <- match.arg(pvalue)
    regularization <- match.arg(regularization)
    pcorr <- match.arg(pcorr, c("BH", "bonferroni", "hochberg", "holm"))
    multicorr <- match.arg(multicorr)
    
    # Validate storey parameter
    if (!is.logical(storey)) {
        stop("storey must be TRUE or FALSE", call. = FALSE)
    }
    
    # Validate wy_randomizations
    if (!is.numeric(wy_randomizations) || wy_randomizations < 1) {
        stop("wy_randomizations must be numeric and >= 1", call. = FALSE)
    }
    if (wy_randomizations < 100) {
        warning("wy_randomizations < 100 may give unreliable p-values; recommend >= 100", call. = FALSE)
    }
    
    # Auto-detect subject_col from colData if paired=TRUE and subject_col=NULL
    # Prioritize 'paired_samples' or 'sample_base' columns (created by calculate_diversity or map_metadata)
    if (paired && is.null(subject_col)) {
        cd_colnames <- colnames(SummarizedExperiment::colData(se))
        
        # Check for paired_samples or sample_base columns
        if ("paired_samples" %in% cd_colnames) {
            subject_col <- "paired_samples"
            if (verbose) {
                message("[calculate_lm_interaction] paired=TRUE detected; auto-using subject_col='paired_samples'")
            }
        } else if ("sample_base" %in% cd_colnames) {
            subject_col <- "sample_base"
            if (verbose) {
                message("[calculate_lm_interaction] paired=TRUE detected; auto-using subject_col='sample_base'")
            }
        } else {
            # Error if paired=TRUE but no recognized pairing column found
            stop("paired=TRUE requires either 'paired_samples' or 'sample_base' column in colData. ",
                 "Available columns: ", paste(cd_colnames, collapse = ", "),
                 ". Ensure calculate_diversity() or map_metadata() was called with appropriate metadata.",
                 call. = FALSE)
        }
    }

    
    if (verbose) {
        message("[calculate_lm_interaction] method=", method)
    }
    # internal flags: keep these internal to avoid documenting them in Rd
    suppress_lme4_warnings <- TRUE
    progress <- FALSE
    if (!requireNamespace("SummarizedExperiment", quietly = TRUE)) {
        stop("SummarizedExperiment required")
    }

    mat <- SummarizedExperiment::assay(se, assay_name)
    if (is.null(mat)) {
        stop(sprintf("Assay '%s' not found in SummarizedExperiment", assay_name))
    }

    sample_q <- colnames(mat)
    if (is.null(sample_q) || length(sample_q) == 0) {
        stop("No column names found on diversity assay")
    }

    # parse sample names and q values from column names like 'Sample_q=0.01'
    sample_names <- sub("_q=.*", "", sample_q)
    has_q <- grepl("_q=", sample_q)
    if (!any(has_q)) {
        stop("Could not parse q values; expected '_q=' in column names", call. = FALSE)
    }
    if (!all(has_q)) {
        stop("Some column names are missing '_q='; ensure all diversity columns include a q value",
            call. = FALSE)
    }
    q_vals <- as.numeric(sub(".*_q=", "", sample_q))

    # determine group for each sample
    condition_in_coldata <- !is.null(condition_col) && condition_col %in% colnames(SummarizedExperiment::colData(se))
    if (condition_in_coldata) {
        st <- as.character(SummarizedExperiment::colData(se)[, condition_col])
        names(st) <- rownames(SummarizedExperiment::colData(se))
        # when user provides a condition_col, index by the FULL column names (sample_q)
        # not by sample_names (which lose the q-value information)
        group_vec <- unname(st[sample_q])
    } else {
        stop("No sample grouping found: please supply `condition_col` or map sample",
            "types into `colData(se)` before calling calculate_lm_interaction().",
            call. = FALSE)
    }

    if (verbose && progress) {
        message("[calculate_lm_interaction] parsed samples and groups")
    }
    
    # ================================================================================
    
    all_results <- list()
    fit_one <- function(g, group_vec_override = NULL) {
        # CI weighting removed (March 2026) - not supported by literature
        # See: CI_WEIGHTING_VALIDATION_REPORT.txt, BY020 (Kotzen), BY021 (Kleijn), S232 (Bayarri & Berger)
        gene_weights <- NULL
        
        # Use override group_vec if provided (for permutation testing), otherwise use outer scope
        gv <- if (!is.null(group_vec_override)) group_vec_override else group_vec
        
        .tsenat_fit_one_interaction(g = g, se = se, mat = mat, q_vals = q_vals, sample_names = sample_names,
            group_vec = gv, method = method, pvalue = pvalue, subject_col = subject_col,
            paired = paired, min_obs = min_obs, verbose = verbose, suppress_lme4_warnings = suppress_lme4_warnings,
            progress = progress, bias_correction = bias_correction, regularization = regularization, corstr = corstr,
            adaptive_knots = adaptive_knots, weights = gene_weights)
    }

    if (nthreads > 1) {
        res_list <- .tsenat_bplapply(rownames(mat), fit_one, nthreads = nthreads)
    } else {
        res_list <- lapply(rownames(mat), fit_one)
    }
    all_results <- Filter(Negate(is.null), res_list)

    if (length(all_results) == 0) {
        return(data.frame())
    }
    res <- do.call(rbind, all_results)
    
    # VALIDATION: Ensure critical columns exist after rbind
    if (nrow(res) == 0) {
        warning("[calculate_lm_interaction] No genes analyzed (all filtered out)", call. = FALSE)
        return(res)
    }
    
    critical_cols <- c("p_interaction", "gene")
    missing_cols <- setdiff(critical_cols, colnames(res))
    if (length(missing_cols) > 0) {
        stop("[calculate_lm_interaction] CRITICAL: Missing columns in results for ",
             method, " method: ", paste(missing_cols, collapse=", "),
             "\nAvailable columns: ", paste(colnames(res), collapse=", "),
             call. = FALSE)
    }
    
    # Ensure Shapiro-Wilk columns exist for methods that add them
    # (GAM and GEE should add them; ensure consistency)
    if (method %in% c("gam", "gee")) {
        if (!"shapiro_p_value" %in% colnames(res)) {
            res$shapiro_p_value <- NA_real_
        }
        if (!"residuals_normal" %in% colnames(res)) {
            res$residuals_normal <- NA
        }
        if (!"n_residuals_tested" %in% colnames(res)) {
            res$n_residuals_tested <- NA_integer_
        }
    }
    
    # Ensure ci_weighted column exists for all methods (Phase 1 tracking)
    # This is set by all .tsenat_*_interaction helpers, but validate it's present
    if (!"ci_weighted" %in% colnames(res)) {
        res$ci_weighted <- NA  # Should not happen, but provide fallback
        if (verbose) {
            warning("[calculate_lm_interaction] ci_weighted column was missing; added as NAs. ",
                    "This suggests a method helper did not properly set ci_weighted.", call. = FALSE)
        }
    }
    
    # Apply primary multi-q p-value adjustment method
    if (multicorr == "hochberg") {
        # Hochberg stepup procedure (FWER control under positive regression dependence)
        res$adj_p_interaction <- .tsenat_hochberg_stepup(res$p_interaction)
        if (verbose) {
            message("[calculate_lm_interaction] Applied Hochberg stepup adjustment for multi-q correlation")
        }
    } else if (multicorr == "westfall-young") {
        # True Westfall-Young permutation procedure (FWER control via empirical null distribution)
        # Note: This is computationally expensive as it requires refitting models for permutations.
        # For large datasets, consider using 'hochberg' instead.
        
        if (verbose) {
            message("[calculate_lm_interaction] Computing true Westfall-Young via ", 
                    wy_randomizations, " permutations (may be slow)...")
        }
        
        # Save original group vector for safe restoration
        group_vec_orig <- group_vec
        
        # Use helper function for WY permutation machinery
        # This consolidates the permutation loop and p-value aggregation logic
        # that was previously duplicated in detect_q_gene_interactions()
        # CRITICAL: Permutation must preserve q-level structure because:
        # - Multi-q tests exhibit AR(1) autoregressive correlation WITHIN each q-level
        # - Exchangeability assumption only holds within q-levels, not globally
        # - Global shuffling violates this assumption and inflates Type I error
        # Solution: Shuffle group assignments separately within each q-level
        perm_result <- .tsenat_westfall_young_permutation(
            n_genes = nrow(res),
            wy_randomizations = wy_randomizations,
            permute_fn = function() {
                # Generate permutation assignment by shuffling group labels separately within each q-level
                # This preserves the multi-q correlation structure and maintains valid exchangeability
                q_unique <- unique(q_vals)
                perm_assignment <- group_vec_orig
                for (q_val in q_unique) {
                    q_idx <- which(q_vals == q_val)
                    perm_assignment[q_idx] <- sample(group_vec_orig[q_idx])
                }
                return(perm_assignment)
            },
            refit_fn = function(perm_assignment) {
                # Refit models with permuted group assignment
                # Pass permuted group_vec explicitly to avoid modifying outer scope
                perm_pvalues <- numeric(nrow(res))
                # Refit each gene with permuted group assignment
                for (g_idx in seq_along(rownames(mat))) {
                    gene_name <- rownames(mat)[g_idx]
                    tryCatch({
                        gene_result <- fit_one(gene_name, group_vec_override = perm_assignment)
                        if (!is.null(gene_result) && !is.na(gene_result$p_interaction)) {
                            perm_pvalues[g_idx] <- gene_result$p_interaction
                        }
                    }, error = function(e) { NULL })
                }
                
                return(perm_pvalues)
            },
            nthreads = nthreads,
            verbose = verbose
        )
        
        # Adjust p-values based on permutation distribution
        # For each observed p-value, compute proportion of permutations with min_perm <= p_obs
        res$adj_p_interaction <- vapply(res$p_interaction, function(p_obs) {
            pmin(1.0, (sum(perm_result$perm_minima <= p_obs) + 1) / (wy_randomizations + 1))
        }, FUN.VALUE = numeric(1))
        
        if (verbose) {
            message("[calculate_lm_interaction] Applied true Westfall-Young (permutation) adjustment")
        }
    } else if (multicorr == "benjamini-yekutieli") {
        # Benjamini-Yekutieli FDR control (valid under any dependence structure)
        res$adj_p_interaction <- .tsenat_benjamini_yekutieli(res$p_interaction)
        if (verbose) {
            message("[calculate_lm_interaction] Applied Benjamini-Yekutieli adjustment for dependent tests")
        }
    }
    
    # Apply optional Storey adaptive FDR enhancement layer
    if (storey) {
        if (requireNamespace("fdrtool", quietly = TRUE)) {
            tryCatch({
                res$adj_p_interaction <- compute_storey_qvalues(res$adj_p_interaction)
                if (verbose) {
                    message("[calculate_lm_interaction] Applied Storey adaptive FDR ?0 correction to ", 
                            multicorr, " p-values")
                }
            }, error = function(e) {
                if (verbose) {
                    message("[calculate_lm_interaction] Storey adjustment failed: ", conditionMessage(e))
                }
            })
        } else if (verbose) {
            message("[calculate_lm_interaction] fdrtool package not available for Storey (install with: install.packages('fdrtool'))")
        }
    }
    
    # Sort first by adjusted p-values, then by raw p-values for stable ordering
    res <- res[order(res$adj_p_interaction, res$p_interaction), , drop = FALSE]
    rownames(res) <- NULL

    .tsenat_report_fit_summary(res, verbose = verbose)

    # Map gene identifiers to gene names from rowData for downstream analysis
    # This ensures effect_sizes_divergence() can match genes across different identifier systems
    rd <- SummarizedExperiment::rowData(se)
    
    # Look for gene_name column from calculate_diversity or build_se
    gene_name_col <- if ("gene_name" %in% colnames(rd)) "gene_name" else NULL
    
    if (verbose) {
      message("[calculate_lm_interaction] Gene annotations: ", paste(colnames(rd), collapse=", "))
    }
    
    if (!is.null(gene_name_col)) {
      # Directly extract gene_id and gene_name from rowData using rownames as key
      # rowData rows correspond to matrix rows in same order
      
      # Determine gene_id column if it exists
      id_col <- if ("genes" %in% colnames(rd)) {
        "genes"
      } else if ("gene_id" %in% colnames(rd)) {
        "gene_id"
      } else {
        NA  # rownames will be used as ID
      }
      
      # Build lookup tables: rowname -> gene_id and rowname -> gene_name
      if (is.na(id_col)) {
        # rownames ARE the gene IDs
        rowname_to_id <- setNames(
          as.character(rownames(rd)),
          as.character(rownames(rd))
        )
      } else {
        # gene IDs are in a column
        rowname_to_id <- setNames(
          as.character(rd[[id_col]]),
          as.character(rownames(rd))
        )
      }
      
      rowname_to_name <- setNames(
        as.character(rd[[gene_name_col]]),
        as.character(rownames(rd))
      )
      
      # Vectorized lookup: map res$gene (rownames) to gene_id and gene_name
      res$gene_id <- unname(rowname_to_id[as.character(res$gene)])
      res$gene_name <- unname(rowname_to_name[as.character(res$gene)])
      
      # For any unmapped genes, use gene column as fallback
      unmapped_idx <- is.na(res$gene_name)
      n_mapped <- sum(!unmapped_idx)
      n_unmapped <- sum(unmapped_idx)
      
      if (any(unmapped_idx)) {
        res$gene_name[unmapped_idx] <- res$gene[unmapped_idx]
      }
      
      if (verbose && n_unmapped > 0) {
        message("[calculate_lm_interaction] Gene mapping: ", n_mapped, " mapped, ", 
                n_unmapped, " used ID as fallback")
      }
    } else if (verbose) {
      message("[calculate_lm_interaction] gene_name column not found in rowData - using gene ID as fallback")
    }
    
    # Ensure gene_name column is always present and populated (for downstream functions)
    if (is.null(res$gene_name) || !"gene_name" %in% colnames(res)) {
      # If gene_name wasn't set above, use gene column as fallback
      res$gene_name <- res$gene
    }
    
    # Ensure gene_id column is always present and populated
    if (is.null(res$gene_id) || !"gene_id" %in% colnames(res)) {
      # If gene_id wasn't set above, use gene column (which may be rownames or gene symbols)
      res$gene_id <- res$gene
    }

    # Optionally return model data alongside results
    if (return_model_data) {
        # Extract per-group statistics from SE
        per_group_stats <- list()
        
        mat <- SummarizedExperiment::assay(se, assay_name)
        cdata <- SummarizedExperiment::colData(se)
        
        for (gr in unique(group_vec)) {
            gr_idx <- which(group_vec == gr)
            gr_mat <- mat[, gr_idx, drop = FALSE]
            
            per_group_stats[[gr]] <- list(
                group = gr,
                n_samples = length(unique(sample_names[gr_idx])),
                n_observations = ncol(gr_mat),
                entropy_mean = mean(as.numeric(gr_mat), na.rm = TRUE),
                entropy_sd = sd(as.numeric(gr_mat), na.rm = TRUE),
                entropy_min = min(as.numeric(gr_mat), na.rm = TRUE),
                entropy_max = max(as.numeric(gr_mat), na.rm = TRUE),
                entropy_median = median(as.numeric(gr_mat), na.rm = TRUE),
                n_na = sum(is.na(gr_mat))
            )
        }
        
        model_data <- list(
            method = method,
            n_genes = nrow(res),
            n_q_values = length(unique(q_vals)),
            q_values = sort(unique(q_vals)),
            sample_names = unique(sample_names),
            group_levels = levels(factor(group_vec)),
            per_group_statistics = per_group_stats,
            test_configuration = list(
                method = method,
                pvalue_method = pvalue,
                multicorr = multicorr,
                bias_correction = bias_correction,
                regularization = regularization,
                corstr = corstr,
                adaptive_knots = adaptive_knots
            ),
            genes_analyzed = res$gene,
            call_time = Sys.time(),
            notes = "Use this model_data with plotting functions to visualize model fits and diagnostics. See per_group_statistics for condition-specific entropy summaries."
        )
        
        return(list(
            results = res,
            model_data = model_data
        ))
    }

    # Return the result data.frame (do not attach to or return a
    # SummarizedExperiment)
    return(res)
}
