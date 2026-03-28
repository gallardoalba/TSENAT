# Helper functions for .calculate_lm_interaction()
# These internal functions decompose the main function logic into
# focused, testable components that each handle a single responsibility.

#' @title Validate Input Parameters for LM Interaction Testing
#'
#' @description
#' Internal helper that consolidates parameter validation for
#' \code{.calculate_lm_interaction()}. Checks argument types, values,
#' and inter-dependencies to ensure valid model fitting.
#'
#' @param method Character; modeling method (matched from user input)
#' @param pvalue Character; p-value type specification
#' @param corstr Character; correlation structure
#' @param regularization Character; dimensionality reduction method
#' @param multicorr Character; multi-q correction method
#' @param pcorr Character; legacy p-value correction method
#' @param storey Logical; whether to apply Storey correction
#' @param wy_randomizations Integer; number of permutations
#' @param paired Logical; whether design is paired
#' @param subject_col Character or NULL; subject column name
#' @param se SummarizedExperiment object
#' @param verbose Logical; print diagnostic messages
#'
#' @return List with validated and normalized parameters:
#'   \itemize{
#'     \item method: Validated method name
#'     \item pvalue: Validated p-value type
#'     \item corstr: Validated correlation structure
#'     \item regularization: Validated regularization method
#'     \item multicorr: Validated multicorr method
#'     \item pcorr: Validated legacy pcorr
#'     \item subject_col: Auto-detected or user-provided subject column
#'   }
#'

#' @noRd
.tsenat_validate_lm_interaction_input <- function(
    method,
    pvalue,
    corstr,
    regularization,
    multicorr,
    pcorr,
    storey,
    wy_randomizations,
    paired,
    subject_col,
    se,
    verbose
) {
    # Validate storey parameter
    if (!is.logical(storey)) {
        stop("storey must be TRUE or FALSE", call. = FALSE)
    }

    # Validate wy_randomizations
    if (!is.numeric(wy_randomizations) || wy_randomizations < 1) {
        stop("wy_randomizations must be numeric and >= 1", call. = FALSE)
    }
    if (wy_randomizations < 100) {
        warning(
            "wy_randomizations < 100 may give unreliable p-values; ",
            "recommend >= 100",
            call. = FALSE
        )
    }

    # Auto-detect subject_col from colData if paired=TRUE and subject_col=NULL
    # Prioritize 'paired_samples' or 'sample_base' columns
    if (paired && is.null(subject_col)) {
        cd_colnames <- colnames(SummarizedExperiment::colData(se))

        # Check for paired_samples or sample_base columns
        if ("paired_samples" %in% cd_colnames) {
            subject_col <- "paired_samples"
            if (verbose) {
                message(
                    "[calculate_lm_interaction] paired=TRUE detected; ",
                    "auto-using subject_col='paired_samples'"
                )
            }
        } else if ("sample_base" %in% cd_colnames) {
            subject_col <- "sample_base"
            if (verbose) {
                message(
                    "[calculate_lm_interaction] paired=TRUE detected; ",
                    "auto-using subject_col='sample_base'"
                )
            }
        } else {
            # Error if paired=TRUE but no recognized pairing column found
            stop(
                "paired=TRUE requires either 'paired_samples' or ",
                "'sample_base' column in colData. Available columns: ",
                paste(cd_colnames, collapse = ", "),
                ". Ensure .calculate_diversity() or map_metadata() was ",
                "called with appropriate metadata.",
                call. = FALSE
            )
        }
    }

    return(list(
        method = method,
        pvalue = pvalue,
        corstr = corstr,
        regularization = regularization,
        multicorr = multicorr,
        pcorr = pcorr,
        subject_col = subject_col
    ))
}

#' @title Parse Sample Metadata from SummarizedExperiment
#'
#' @description
#' Internal helper that extracts sample names, q-values, and group
#' assignments from the diversity assay column names and colData.
#'
#' @param se SummarizedExperiment object
#' @param condition_col Character; colData column with group assignments
#' @param assay_name Character; name of diversity assay
#' @param verbose Logical; print diagnostic messages
#'
#' @return List containing:
#'   \itemize{
#'     \item sample_q: Full column names with q= values
#'     \item sample_names: Unique sample identifiers
#'     \item q_vals: Parsed q-value parameters
#'     \item group_vec: Group assignment for each observation
#'     \item has_q: Logical vector indicating cols with q=
#'   }
#'

#' @noRd
.tsenat_parse_sample_metadata <- function(
    se,
    condition_col,
    assay_name,
    verbose
) {
    mat <- SummarizedExperiment::assay(se, assay_name)
    if (is.null(mat)) {
        stop(sprintf("Assay '%s' not found in SummarizedExperiment",
                     assay_name))
    }

    sample_q <- colnames(mat)
    if (is.null(sample_q) || length(sample_q) == 0) {
        stop("No column names found on diversity assay")
    }

    # Parse sample names and q values from column names like
    # 'Sample_q=0.01'
    sample_names <- sub("_q=.*", "", sample_q)
    has_q <- grepl("_q=", sample_q)
    if (!any(has_q)) {
        stop(
            "Could not parse q values; expected '_q=' in column names",
            call. = FALSE
        )
    }
    if (!all(has_q)) {
        stop(
            "Some column names are missing '_q='; ensure all diversity ",
            "columns include a q value",
            call. = FALSE
        )
    }
    q_vals <- as.numeric(sub(".*_q=", "", sample_q))

    # Determine group for each sample
    condition_in_coldata <- !is.null(condition_col) &&
        condition_col %in% colnames(SummarizedExperiment::colData(se))
    if (condition_in_coldata) {
        st <- as.character(
            SummarizedExperiment::colData(se)[, condition_col]
        )
        names(st) <- rownames(SummarizedExperiment::colData(se))
        # Index by the FULL column names (sample_q), not by sample_names
        group_vec <- unname(st[sample_q])
    } else {
        stop(
            "No sample grouping found: please supply `condition_col` ",
            "or map sample types into `colData(se)` before calling ",
            ".calculate_lm_interaction().",
            call. = FALSE
        )
    }

    if (verbose) {
        message(
            "[calculate_lm_interaction] parsed samples and groups: ",
            length(unique(sample_names)), " samples, ",
            length(unique(q_vals)), " q-values"
        )
    }

    return(list(
        sample_q = sample_q,
        sample_names = sample_names,
        q_vals = q_vals,
        group_vec = group_vec,
        has_q = has_q
    ))
}

#' @title Fit Linear Models for All Genes
#'
#' @description
#' Internal helper that orchestrates parallel or sequential fitting
#' of models to all genes in the diversity matrix. Consolidates the
#' fitting loop and result collection logic.
#'
#' @param mat Matrix; diversity assay data
#' @param se SummarizedExperiment object
#' @param metadata List; output from .tsenat_parse_sample_metadata()
#' @param method Character; modeling method
#' @param pvalue Character; p-value type
#' @param subject_col Character or NULL; subject column
#' @param paired Logical; whether design is paired
#' @param min_obs Integer; minimum observations per gene
#' @param nthreads Integer; number of parallel threads
#' @param verbose Logical; print diagnostics
#' @param bias_correction Logical; apply KC bias correction (GEE)
#' @param regularization Character; dimensionality reduction method
#' @param corstr Character; correlation structure
#' @param adaptive_knots Logical; adaptive knot selection (GAM)
#'
#' @return Data.frame with fitted model results for all genes
#'

#' @noRd
.tsenat_fit_all_genes <- function(
    mat,
    se,
    metadata,
    method,
    pvalue,
    subject_col,
    paired,
    min_obs,
    nthreads,
    verbose,
    bias_correction,
    regularization,
    corstr,
    adaptive_knots
) {
    suppress_lme4_warnings <- TRUE
    progress <- FALSE
    gene_weights <- NULL

    fit_one <- function(g, group_vec_override = NULL) {
        # Use override group_vec if provided (for permutation testing),
        # otherwise use outer scope
        gv <- if (!is.null(group_vec_override)) {
            group_vec_override
        } else {
            metadata$group_vec
        }

        .tsenat_fit_one_interaction(
            g = g,
            se = se,
            mat = mat,
            q_vals = metadata$q_vals,
            sample_names = metadata$sample_names,
            group_vec = gv,
            method = method,
            pvalue = pvalue,
            subject_col = subject_col,
            paired = paired,
            min_obs = min_obs,
            verbose = verbose,
            suppress_lme4_warnings = suppress_lme4_warnings,
            progress = progress,
            bias_correction = bias_correction,
            regularization = regularization,
            corstr = corstr,
            adaptive_knots = adaptive_knots,
            weights = gene_weights
        )
    }

    if (nthreads > 1) {
        res_list <- .tsenat_bplapply(rownames(mat), fit_one,
                                      nthreads = nthreads)
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
        warning(
            "[calculate_lm_interaction] No genes analyzed (all filtered out)",
            call. = FALSE
        )
        return(res)
    }

    critical_cols <- c("p_interaction", "gene")
    missing_cols <- setdiff(critical_cols, colnames(res))
    if (length(missing_cols) > 0) {
        stop(
            "[calculate_lm_interaction] CRITICAL: Missing columns in ",
            "results for ", method, " method: ",
            paste(missing_cols, collapse = ", "),
            "\nAvailable columns: ",
            paste(colnames(res), collapse = ", "),
            call. = FALSE
        )
    }

    # Ensure Shapiro-Wilk columns exist for methods that add them
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

    # Ensure ci_weighted column exists (Phase 1 tracking)
    if (!"ci_weighted" %in% colnames(res)) {
        res$ci_weighted <- NA  # Fallback
        if (verbose) {
            warning(
                "[calculate_lm_interaction] ci_weighted column was ",
                "missing; added as NAs. This suggests a method helper ",
                "did not properly set ci_weighted.",
                call. = FALSE
            )
        }
    }

    return(res)
}

#' @title Adjust P-Values for Multiple Q-Values
#'
#' @description
#' Internal router function that applies the specified primary multi-q
#' p-value correction method (Hochberg, Westfall-Young, or
#' Benjamini-Yekutieli), then optionally applies Storey adaptive FDR
#' enhancement.
#'
#' @param p_values Numeric vector; raw p-values to adjust
#' @param multicorr Character; primary correction method
#' @param wy_randomizations Integer; number of permutations (WY only)
#' @param fit_one_fn Function; function to refit models (WY only)
#' @param metadata List; output from .tsenat_parse_sample_metadata()
#' @param mat Matrix; diversity assay data (WY only)
#' @param rownames_mat Character; rownames of matrix (WY only)
#' @param se SummarizedExperiment object (WY only)
#' @param assay_name Character; assay name (WY only)
#' @param method Character; modeling method (WY only)
#' @param pvalue Character; p-value type (WY only)
#' @param subject_col Character or NULL; subject column (WY only)
#' @param paired Logical; paired design (WY only)
#' @param min_obs Integer; min observations (WY only)
#' @param nthreads Integer; parallel threads (WY only)
#' @param verbose Logical; print diagnostics
#' @param bias_correction Logical; KC bias correction (WY/GEE)
#' @param regularization Character; dimensionality reduction (WY/FPCA)
#' @param corstr Character; correlation structure (WY/GEE)
#' @param adaptive_knots Logical; adaptive knots (WY/GAM)
#' @param storey Logical; apply Storey after primary correction
#'
#' @return Adjusted p-values vector
#'

#' @noRd
.tsenat_adjust_pvalues_multicorr <- function(
    p_values,
    multicorr,
    wy_randomizations,
    fit_one_fn = NULL,
    metadata = NULL,
    mat = NULL,
    rownames_mat = NULL,
    se = NULL,
    assay_name = "diversity",
    method = NULL,
    pvalue = NULL,
    subject_col = NULL,
    paired = FALSE,
    min_obs = 10,
    nthreads = 1,
    verbose = FALSE,
    bias_correction = TRUE,
    regularization = NULL,
    corstr = "ar1",
    adaptive_knots = TRUE,
    storey = FALSE
) {
    if (multicorr == "hochberg") {
        adj_p <- .tsenat_hochberg_stepup(p_values)
        if (verbose) {
            message(
                "[calculate_lm_interaction] Applied Hochberg stepup ",
                "adjustment for multi-q correlation"
            )
        }
    } else if (multicorr == "westfall-young") {
        if (verbose) {
            message(
                "[calculate_lm_interaction] Computing true ",
                "Westfall-Young via ", wy_randomizations,
                " permutations (may be slow)..."
            )
        }

        # Save original group vector for safe restoration
        group_vec_orig <- metadata$group_vec

        # Run Westfall-Young permutation
        perm_result <- .tsenat_westfall_young_permutation(
            n_genes = length(p_values),
            wy_randomizations = wy_randomizations,
            permute_fn = function() {
                # Shuffle group labels separately within each q-level
                q_unique <- unique(metadata$q_vals)
                perm_assignment <- group_vec_orig
                for (q_val in q_unique) {
                    q_idx <- which(metadata$q_vals == q_val)
                    perm_assignment[q_idx] <-
                        sample(group_vec_orig[q_idx])
                }
                return(perm_assignment)
            },
            refit_fn = function(perm_assignment) {
                # Refit all genes with permuted group assignment
                perm_pvalues <- numeric(length(p_values))
                for (g_idx in seq_along(rownames_mat)) {
                    gene_name <- rownames_mat[g_idx]
                    tryCatch({
                        gene_result <- fit_one_fn(
                            gene_name,
                            group_vec_override = perm_assignment
                        )
                        if (!is.null(gene_result) &&
                            !is.na(gene_result$p_interaction)) {
                            perm_pvalues[g_idx] <-
                                gene_result$p_interaction
                        }
                    }, error = function(e) { NULL })
                }
                return(perm_pvalues)
            },
            nthreads = nthreads,
            verbose = verbose
        )

        # Adjust p-values based on permutation distribution
        adj_p <- vapply(p_values, function(p_obs) {
            pmin(1.0,
                 (sum(perm_result$perm_minima <= p_obs) + 1) /
                     (wy_randomizations + 1))
        }, FUN.VALUE = numeric(1))

        if (verbose) {
            message(
                "[calculate_lm_interaction] Applied true ",
                "Westfall-Young (permutation) adjustment"
            )
        }
    } else if (multicorr == "benjamini-yekutieli") {
        adj_p <- .tsenat_benjamini_yekutieli(p_values)
        if (verbose) {
            message(
                "[calculate_lm_interaction] Applied ",
                "Benjamini-Yekutieli adjustment for dependent tests"
            )
        }
    } else {
        stop("Unknown multicorr method: ", multicorr, call. = FALSE)
    }

    # Apply optional Storey adaptive FDR enhancement layer
    if (storey) {
        if (requireNamespace("fdrtool", quietly = TRUE)) {
            tryCatch({
                adj_p <- .compute_storey_qvalues(adj_p)
                if (verbose) {
                    message(
                        "[calculate_lm_interaction] Applied Storey ",
                        "adaptive FDR pi0 correction to ",
                        multicorr, " p-values"
                    )
                }
            }, error = function(e) {
                if (verbose) {
                    message(
                        "[calculate_lm_interaction] Storey adjustment ",
                        "failed: ", conditionMessage(e)
                    )
                }
            })
        } else if (verbose) {
            message(
                "[calculate_lm_interaction] fdrtool package not ",
                "available for Storey (install with: ",
                "install.packages('fdrtool'))"
            )
        }
    }

    return(adj_p)
}

#' @title Map Gene Identifiers to Annotations
#'
#' @description
#' Internal helper that maps gene rownames to gene_id and gene_name
#' columns using rowData from the SummarizedExperiment. Ensures
#' consistent gene annotation across downstream analyses.
#'
#' @param res Data.frame; results with gene column (rownames)
#' @param se SummarizedExperiment object
#' @param verbose Logical; print diagnostic messages
#'
#' @return Modified data.frame with gene_id and gene_name columns added
#'

#' @noRd
.tsenat_map_gene_annotations <- function(
    res,
    se,
    verbose
) {
    rd <- SummarizedExperiment::rowData(se)

    # Look for gene_name column from calculate_diversity or build_se
    gene_name_col <- if ("gene_name" %in% colnames(rd)) {
        "gene_name"
    } else {
        NULL
    }

    if (verbose) {
        message(
            "[calculate_lm_interaction] Gene annotations: ",
            paste(colnames(rd), collapse = ", ")
        )
    }

    if (!is.null(gene_name_col)) {
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

        # Vectorized lookup: map res$gene to gene_id and gene_name
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
            message(
                "[calculate_lm_interaction] Gene mapping: ",
                n_mapped, " mapped, ", n_unmapped,
                " used ID as fallback"
            )
        }
    } else if (verbose) {
        message(
            "[calculate_lm_interaction] gene_name column not found ",
            "in rowData - using gene ID as fallback"
        )
    }

    # Ensure gene_name column is always present and populated
    if (is.null(res$gene_name) || !"gene_name" %in% colnames(res)) {
        res$gene_name <- res$gene
    }

    # Ensure gene_id column is always present and populated
    if (is.null(res$gene_id) || !"gene_id" %in% colnames(res)) {
        res$gene_id <- res$gene
    }

    return(res)
}

#' @title Assemble Model Metadata for Diagnostics
#'
#' @description
#' Internal helper that builds the comprehensive metadata list returned
#' when \code{return_model_data = TRUE}. Contains method info, q-values,
#' per-group statistics, and test configuration for downstream visualization.
#'
#' @param se SummarizedExperiment object
#' @param res Data.frame; fitted model results
#' @param mat Matrix; diversity assay data
#' @param metadata List; output from .tsenat_parse_sample_metadata()
#' @param method Character; modeling method name
#' @param pvalue Character; p-value type
#' @param multicorr Character; multi-q correction method
#' @param assay_name Character; assay name
#' @param bias_correction Logical; KC bias correction setting
#' @param regularization Character; dimensionality reduction method
#' @param corstr Character; correlation structure
#' @param adaptive_knots Logical; adaptive knot selection setting
#'
#' @return List with comprehensive model metadata:
#'   \itemize{
#'     \item method: Modeling method used
#'     \item n_genes: Number of genes analyzed
#'     \item n_q_values: Number of q parameters
#'     \item q_values: The specific q-value vector
#'     \item sample_names: Unique sample identifiers
#'     \item group_levels: Group factor levels
#'     \item per_group_statistics: Summary statistics by group
#'     \item test_configuration: Complete test settings
#'     \item genes_analyzed: Vector of gene identifiers
#'     \item call_time: Timestamp of analysis
#'     \item notes: Usage information
#'   }
#'

#' @noRd
.tsenat_assemble_model_metadata <- function(
    se,
    res,
    mat,
    metadata,
    method,
    pvalue,
    multicorr,
    assay_name = "diversity",
    bias_correction = TRUE,
    regularization = "pca",
    corstr = "ar1",
    adaptive_knots = TRUE
) {
    # Extract per-group statistics from SE
    per_group_stats <- list()

    for (gr in unique(metadata$group_vec)) {
        gr_idx <- which(metadata$group_vec == gr)
        gr_mat <- mat[, gr_idx, drop = FALSE]

        per_group_stats[[gr]] <- list(
            group = gr,
            n_samples = length(unique(metadata$sample_names[gr_idx])),
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
        n_q_values = length(unique(metadata$q_vals)),
        q_values = sort(unique(metadata$q_vals)),
        sample_names = unique(metadata$sample_names),
        group_levels = levels(factor(metadata$group_vec)),
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
        notes = paste(
            "Use this model_data with plotting functions to ",
            "visualize model fits and diagnostics. See ",
            "per_group_statistics for condition-specific ",
            "entropy summaries."
        )
    )

    return(model_data)
}
