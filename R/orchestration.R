# Unified TSENAT Analysis Orchestration Main entry point and configuration for
# TSENAT pipeline. Coordinates analysis workflow from raw counts to results and
# visualizations.

# ============================================================================
# MAIN ORCHESTRATION FUNCTION
# ============================================================================

#' Run complete TSENAT analysis pipeline
#'
#' Coordinates the full TSENAT workflow: diversity -> jackknife -> LM
#' interactions -> divergence -> gene interactions -> visualizations.
#'
#' @param se \code{SummarizedExperiment} containing expression counts.
#' @param config \code{list} or \code{TSENATConfig}. Configuration from
#'   \code{\link{tsenat_config}}. If NULL, uses defaults.
#' @param methods \code{character}. Specific methods to run (overrides config).
#' @param q_values \code{numeric}. Specific q-values (overrides config).
#' @param generate_plots \code{logical}. Create visualizations. Default: TRUE.
#' @param verbose \code{logical}. Print progress messages. Default: TRUE.
#' @param parallel \code{logical}. Run independent q-values in parallel.
#'   Default: FALSE.
#' @param ... Additional arguments passed to individual wrapper functions.
#'
#' @return \code{TSENATAnalysis} object containing complete analysis results,
#'   plots, and metadata.
#'
#' @details
#' Pipeline execution order (enforced):
#' \enumerate{
#'   \item \code{calculate_diversity_s4()} - Tsallis entropy per q-value
#'   \item \code{jackknife_entropy_outliers_s4()} - Confidence intervals
#'   \item \code{jackknife_isoform_switching_s4()} - Transcript switching (optional)
#'   \item \code{calculate_lm_interaction_s4()} - Statistical tests
#'   \item \code{calculate_difference_s4()} - Pairwise group differences (optional)
#'   \item \code{calculate_divergence_s4()} - Pairwise divergence metrics
#'   \item \code{effect_sizes_divergence_s4()} - Effect sizes (optional)
#'   \item \code{rank_test_q_condition_s4()} - Q-dependent interactions
#'   \item \code{test_rankbased_assumptions_s4()} - Assumption validation (optional)
#'   \item \code{compute_method_concordance_s4()} - Method comparison (optional)
#'   \item Plot generation (if enabled)
#' }
#'
#' @examples
#' library(SummarizedExperiment)
#' se <- SummarizedExperiment(
#'   assays = list(counts = matrix(rpois(100, 10), nrow = 10, ncol = 10,
#'                                  dimnames = list(sprintf('G%d', 1:10), NULL))),
#'   colData = data.frame(
#'     condition = rep(c('A', 'B'), 5),
#'     sample_id = rep(c('S1', 'S2'), 5)
#'   )
#' )
#' cfg <- tsenat_config(q_values = c(0.5, 1.0), generate_plots = FALSE)
#' analysis <- tsenat(se, config = cfg)
#'
#' @export
tsenat <- function(se, config = NULL, methods = NULL, q_values = NULL, generate_plots = TRUE,
    verbose = TRUE, parallel = FALSE, ...) {
    # Validate input
    if (!is(se, "SummarizedExperiment")) {
        stop("'se' must be a SummarizedExperiment object", call. = FALSE)
    }
    if (nrow(se) == 0) {
        stop("SummarizedExperiment is empty", call. = FALSE)
    }

    # Initialize and setup parameters
    analysis <- TSENATAnalysis(se = se, config = config)
    .validate_analysis_object(analysis)
    
    params <- .setup_tsenat_parameters(analysis, methods, q_values, generate_plots, parallel)
    analysis <- params$analysis
    methods_to_run <- params$methods_to_run
    q_vals <- params$q_vals
    condition_col_name <- params$condition_col_name
    do_plots <- params$do_plots

    # Validate and execute pipeline
    .validate_tsenat_methods(methods_to_run)
    if (verbose) .log_pipeline_start(se, methods_to_run, q_vals)

    analysis <- .execute_diversity_step(analysis, q_vals, methods_to_run, verbose, ...)
    analysis <- .execute_jackknife_step(analysis, q_vals, methods_to_run, verbose, ...)
    analysis <- .execute_jackknife_isoform_switching_step(analysis, q_vals, condition_col_name,
        methods_to_run, verbose, ...)
    analysis <- .execute_lm_interaction_step(analysis, methods_to_run, verbose, ...)
    analysis <- .execute_difference_step(analysis, q_vals, condition_col_name, methods_to_run, verbose, ...)
    analysis <- .execute_divergence_step(analysis, q_vals, methods_to_run, verbose, ...)
    analysis <- .execute_effect_sizes_step(analysis, methods_to_run, verbose, ...)
    analysis <- .execute_q_interactions_step(analysis, q_vals, condition_col_name, 
        methods_to_run, verbose, ...)
    analysis <- .execute_rankbased_assumptions_step(analysis, q_vals, methods_to_run, verbose, ...)
    analysis <- .execute_method_concordance_step(analysis, methods_to_run, verbose, ...)
    analysis <- .execute_plot_generation(analysis, do_plots, verbose)
    
    # Track completion metadata
    analysis <- .track_analysis_metadata(analysis, methods_to_run, analysis@config)
    analysis <- .finalize_tsenat_analysis(analysis, verbose)

    analysis
}

# ============================================================================
# UTILITY FUNCTION
# ============================================================================

#' `%||%` operator for default values
#'
#' Returns left operand if not NULL, otherwise right operand.
#'

#' @noRd
`%||%` <- function(x, y) {
    if (is.null(x))
        y else x
}

# ============================================================================
# RESULT ACCESSOR FUNCTIONS (EXPORTED)
# ============================================================================

#' Extract analysis results from TSENATAnalysis object
#'
#' Provides flexible access to diversity, divergence, and statistical test results.
#'
#' @param analysis \code{TSENATAnalysis} object containing computed results.
#' @param type \code{character}. Type of results to extract:
#'   'diversity', 'divergence', 'lm', 'jackknife', or 'q_interactions'.
#'   Default: 'diversity'.
#' @param q \code{numeric}. For diversity results, optionally filter by q-value.
#'   Default: NULL (return all q-values).
#' @param simplify \code{logical}. If TRUE and q is specified, return as 
#'   vector instead of matrix. Default: TRUE.
#'
#' @return Extracted results as data.frame, matrix, or list depending on type.
#'   Returns NULL if requested result type not computed.
#'
#' @details
#' This function provides a consistent interface to access all computed results
#' from the TSENATAnalysis object, abstracting away internal storage details.
#'
#' @examples
#' # Load example data
#' data(readcounts, package = "TSENAT")
#'
#' # Create SummarizedExperiment from counts
#' se <- SummarizedExperiment::SummarizedExperiment(
#'   assays = list(counts = readcounts),
#'   rowData = data.frame(gene_id = rownames(readcounts)),
#'   colData = data.frame(
#'     sample_id = colnames(readcounts),
#'     group = rep(c("A", "B"), length.out = ncol(readcounts))
#'   )
#' )
#'
#' # Create TSENATAnalysis object with configuration
#' analysis <- TSENATAnalysis(
#'   se = se,
#'   config = list(
#'     q_values = c(0.5, 1.0, 2.0),
#'     condition_col = "group"
#'   )
#' )
#'
#' # Run analysis to generate diversity results
#' analysis <- calculate_diversity_s4(analysis)
#'
#' # Extract diversity results
#' div_results <- getResults(analysis, type = "diversity")
#' if (!is.null(div_results)) {
#'   head(div_results, n = 3)
#' }
#'
#' @export
getResults <- function(analysis, type = "diversity", q = NULL, simplify = TRUE) {
    if (!is(analysis, "TSENATAnalysis")) {
        stop("'analysis' must be a TSENATAnalysis object", call. = FALSE)
    }

    result <- switch(type,
        diversity = if (length(analysis@diversity_results) > 0) 
            analysis@diversity_results else NULL,
        divergence = if (length(analysis@divergence_results) > 0) 
            analysis@divergence_results else NULL,
        lm = if (length(analysis@lm_results) > 0) 
            analysis@lm_results else NULL,
        jackknife = if (length(analysis@jackknife_results) > 0) 
            analysis@jackknife_results else NULL,
        q_interactions = if ("q_interactions" %in% names(analysis@lm_results))
            analysis@lm_results$q_interactions else NULL,
        stop("Unknown result type: '", type, "'. Must be one of: ", 
            "diversity, divergence, lm, jackknife, q_interactions", call. = FALSE)
    )

    if (is.null(result)) {
        return(NULL)
    }

    # Filter by q-value if specified and applicable
    if (!is.null(q) && type == "diversity" && is.matrix(result)) {
        if (simplify) {
            result <- result[as.character(q), , drop = TRUE]
        } else {
            result <- result[as.character(q), , drop = FALSE]
        }
    }

    return(result)
}

# ============================================================================
# CONFIG BUILDER: tsenat_config()
# ============================================================================

#' Create and return TSENAT configuration
#'
#' Builds a configuration list for use with \code{\link{tsenat}}().
#' Allows specifying analysis parameters once and reusing across multiple
#' analyses.
#'
#' @param q_values \code{numeric}. Q-values for Tsallis entropy spectrum.
#'   Default: \code{seq(0.5, 2.0, by = 0.5)}.
#' @param condition_col \code{character}. Name of column in \code{colData(se)}
#'   containing experimental conditions/groups. Default: 'condition'.
#' @param subject_col \code{character}. Name of column in \code{colData(se)}
#'   containing subject/sample identifiers for paired/repeated designs.
#'   If provided, enables paired analysis. Default: NULL (unpaired).
#' @param paired \code{logical}. Whether samples are paired/repeated measures.
#'   Default: FALSE. Used by jackknife and difference analysis.
#' @param control \code{character}. Reference/control group label for difference
#'   analysis (e.g., 'control', 'wt'). Only used if 'difference' in methods.
#'   Default: NULL.
#' @param formula \code{formula}. Optional model formula for LM interactions
#'   (e.g., \code{~ treatment + batch}). Default: NULL (uses condition_col).
#' @param p_threshold \code{numeric}. Raw p-value threshold for significance
#'   in LM interaction testing. Default: 0.05.
#' @param fdr_threshold \code{numeric}. Adjusted p-value (FDR/Benjamini-Hochberg)
#'   threshold. Default: 0.05.
#' @param significance_threshold \code{numeric}. Significance cutoff for effect
#'   sizes, assumptions testing, and result filtering. Default: 0.05.
#' @param n_bootstrap \code{integer}. Number of bootstrap resamples for jackknife
#'   confidence intervals. Default: 1000.
#' @param bootstrap_method \code{character}. Bootstrap CI method: 'percentile'
#'   (fast, assumes symmetric distribution) or 'bca' (bias-corrected, better for
#'   skewed data like bounded entropy). Default: 'percentile'.
#' @param methods \code{character}. Analysis steps to run. Core (auto-required):
#'   'diversity', 'lm_interaction', 'jackknife', 'divergence', 'q_interactions'.
#'   Optional: 'jackknife_isoform_switching' (requires condition/subject info),
#'   'difference' (requires control parameter), 'effect_sizes', 'rankbased_assumptions',
#'   'method_concordance'. Default: core methods.
#' @param generate_plots \code{logical}. Generate visualizations.
#'   Default: TRUE.
#' @param plot_types \code{character}. Specific plot types to generate. Options:
#'   'q_curve', 'lm_interaction', 'divergence_distribution', 'divergence_spectrum',
#'   'influence_heatmap', 'volcano', 'method_concordance', 'multi_gene_q_spectrum',
#'   'top_transcripts', 'tsallis_violin_density'. Default: all available types.
#' @param seed \code{numeric}. Random seed for reproducibility (affects
#'   bootstrap resampling). Default: NULL (no fixed seed).
#' @param ... Additional configuration parameters (stored as-is in @config slot).
#'   Examples: \code{q_diff=1.0} (specific q for differences),
#'   \code{alpha=0.05} (significance for assumptions), etc.
#'
#' @return \code{list} with class \code{TSENATConfig} containing all
#'   specified parameters.
#'
#' @details
#' Configuration is stored in the TSENATAnalysis@config slot and used
#' by wrapper functions to configure analysis behavior.
#'
#' @examples
#' # Default config with standard parameters
#' cfg <- tsenat_config()
#'
#' # Custom with optional features: paired analysis, isoform switching, differences
#' cfg <- tsenat_config(
#'   q_values = c(0.5, 1.0, 1.5, 2.0),
#'   condition_col = "treatment",
#'   subject_col = "subject_id",
#'   paired = TRUE,
#'   control = "untreated",
#'   methods = c("diversity", "lm_interaction", "jackknife", "jackknife_isoform_switching",
#'               "divergence", "q_interactions", "difference", "effect_sizes"),
#'   bootstrap_method = "bca",
#'   n_bootstrap = 5000,
#'   significance_threshold = 0.01
#' )
#'
#' @export
tsenat_config <- function(q_values = NULL, condition_col = "condition", subject_col = NULL,
    paired = FALSE, control = NULL, formula = NULL, p_threshold = 0.05, fdr_threshold = 0.05,
    significance_threshold = 0.05, n_bootstrap = 1000, bootstrap_method = "percentile",
    methods = NULL, generate_plots = TRUE, plot_types = NULL, ...) {
    # Build q_values if range specified
    if (is.null(q_values)) {
        q_values <- seq(0.5, 2, by = 0.5)
    }

    # Default methods
    if (is.null(methods)) {
        methods <- c("diversity", "lm_interaction", "jackknife", "divergence", "q_interactions")
    }

    # Validate methods
    valid_methods <- c("diversity", "lm_interaction", "jackknife", "jackknife_isoform_switching",
        "divergence", "q_interactions", "difference", "rankbased_assumptions",
        "method_concordance", "effect_sizes")
    invalid_methods <- setdiff(methods, valid_methods)
    if (length(invalid_methods) > 0) {
        stop("Invalid methods: ", paste(invalid_methods, collapse = ", "), "\n",
            "Valid: ", paste(valid_methods, collapse = ", "), call. = FALSE)
    }

    # Validate bootstrap_method
    valid_bootstrap_methods <- c("percentile", "bca")
    if (!bootstrap_method %in% valid_bootstrap_methods) {
        stop("'bootstrap_method' must be 'percentile' or 'bca'", call. = FALSE)
    }

    # Build config list with all parameters
    config <- list(
        q_values = q_values,
        condition_col = condition_col,
        subject_col = subject_col,
        paired = paired,
        control = control,
        p_threshold = p_threshold,
        fdr_threshold = fdr_threshold,
        significance_threshold = significance_threshold,
        n_bootstrap = n_bootstrap,
        bootstrap_method = bootstrap_method,
        methods = methods,
        generate_plots = generate_plots
    )

    # Add optional parameters
    if (!is.null(formula))
        config$formula <- formula
    if (!is.null(plot_types))
        config$plot_types <- plot_types

    # Add any additional parameters
    extra_args <- list(...)
    if (length(extra_args) > 0) {
        config <- c(config, extra_args)
    }

    # Mark as TSENATConfig (but keep as list for S4 slot)
    attr(config, "class") <- c("TSENATConfig", "list")
    config
}

# ============================================================================
# HELPER FUNCTIONS (INTERNAL - NOT EXPORTED)
# ============================================================================

#' Setup TSENAT parameters from config
#' @noRd
.setup_tsenat_parameters <- function(analysis, methods, q_values, generate_plots, parallel) {
    `%||%` <- function(x, y) if (is.null(x)) y else x
    
    if (!is.null(methods))
        analysis@config$methods <- methods
    if (!is.null(q_values))
        analysis@config$q_values <- q_values

    list(
        methods_to_run = analysis@config$methods %||% c("diversity", "lm_interaction", 
            "jackknife", "divergence", "q_interactions"),
        q_vals = analysis@config$q_values %||% seq(0.5, 2, by = 0.5),
        condition_col_name = analysis@config$condition_col %||% "condition",
        do_plots = generate_plots && (analysis@config$generate_plots %||% TRUE),
        do_parallel = parallel && ("parallel" %in% rownames(utils::installed.packages())),
        analysis = analysis
    )
}

#' Validate method dependencies
#' @noRd
.validate_tsenat_methods <- function(methods_requested) {
    dependencies <- list(jackknife = "diversity", divergence = "diversity", 
        q_interactions = "diversity", lm_interaction = "diversity")
    
    for (method in methods_requested) {
        if (method %in% names(dependencies)) {
            required <- dependencies[[method]]
            if (!(required %in% methods_requested)) {
                stop("Method '", method, "' requires '", required, 
                    "' to be in methods list.", call. = FALSE)
            }
        }
    }
}

#' Log pipeline start
#' @noRd
.log_pipeline_start <- function(se, methods_to_run, q_vals) {
    message("TSENAT Pipeline")
    message("===============")
    message(sprintf("Genes:   %d", nrow(se)))
    message(sprintf("Samples: %d", ncol(se)))
    message(sprintf("Methods: %s", paste(methods_to_run, collapse = ", ")))
    message(sprintf("Q-values: %s", paste(q_vals, collapse = ", ")))
}

#' Execute diversity step
#' @noRd
.execute_diversity_step <- function(analysis, q_vals, methods_to_run, verbose, ...) {
    if (!("diversity" %in% methods_to_run)) return(analysis)
    
    if (verbose) message("Step 1: Calculating diversity...")
    tryCatch({
        analysis <- calculate_diversity_s4(analysis, q = q_vals, ...)
        if (verbose) message(sprintf("  [OK] Diversity for q = %s", 
            paste(q_vals, collapse = ", ")))
    }, error = function(e) stop("Diversity failed:\n", e$message, call. = FALSE))
    analysis
}

#' Execute jackknife step
#' @noRd
.execute_jackknife_step <- function(analysis, q_vals, methods_to_run, verbose, ...) {
    if (!("jackknife" %in% methods_to_run)) return(analysis)
    
    if (length(analysis@diversity_results) == 0) {
        if (verbose) message("Step 2: Skipping jackknife (requires diversity)")
        return(analysis)
    }
    
    if (verbose) message("Step 2: Running jackknife resampling...")
    tryCatch({
        analysis <- jackknife_entropy_outliers_s4(analysis, q = q_vals, verbose = FALSE, ...)
        if (verbose) message("  [OK] Jackknife CIs computed")
    }, error = function(e) warning("Jackknife failed:\n", e$message, call. = FALSE))
    analysis
}

#' Execute LM interaction step
#' @noRd
.execute_lm_interaction_step <- function(analysis, methods_to_run, verbose, ...) {
    if (!("lm_interaction" %in% methods_to_run)) return(analysis)
    
    if (verbose) message("Step 4: Testing LM interactions...")
    tryCatch({
        fdr <- if (is.null(analysis@config$fdr_threshold)) 0.05 else analysis@config$fdr_threshold
        analysis <- calculate_lm_interaction_s4(analysis, fdr_threshold = fdr, ...)
        if (verbose) message("  [OK] LM analysis complete")
    }, error = function(e) warning("LM failed:\n", e$message, call. = FALSE))
    analysis
}

#' Execute difference step
#' @noRd
.execute_difference_step <- function(analysis, q_vals, condition_col_name, methods_to_run, verbose, ...) {
    if (!("difference" %in% methods_to_run)) return(analysis)
    
    if (length(analysis@diversity_results) == 0) {
        if (verbose) message("Step 5: Skipping difference (requires diversity)")
        return(analysis)
    }
    
    if (verbose) message("Step 5: Testing pairwise differences...")
    tryCatch({
        control <- if (!is.null(analysis@config$control)) {
            analysis@config$control
        } else {
            if (verbose) message("  [SKIP] No 'control' group specified in config")
            return(analysis)
        }
        
        q_to_use <- if (!is.null(analysis@config$q_diff)) {
            analysis@config$q_diff
        } else {
            NULL  # Will use first available from diversity_results
        }
        
        analysis <- calculate_difference_s4(analysis, control = control, q = q_to_use, 
            condition_col = condition_col_name, ...)
        if (verbose) message("  [OK] Difference analysis complete")
    }, error = function(e) warning("Difference failed:\n", e$message, call. = FALSE))
    analysis
}

#' Execute divergence step
#' @noRd
.execute_divergence_step <- function(analysis, q_vals, methods_to_run, verbose, ...) {
    if (!("divergence" %in% methods_to_run)) return(analysis)
    
    if (length(analysis@diversity_results) == 0) {
        if (verbose) message("Step 6: Skipping divergence (requires diversity)")
        return(analysis)
    }
    
    if (verbose) message("Step 6: Calculating divergence metrics...")
    tryCatch({
        analysis <- calculate_divergence_s4(analysis, q = q_vals[1])
        if (verbose) message("  [OK] Divergence computed")
    }, error = function(e) warning("Divergence failed:\n", e$message, call. = FALSE))
    analysis
}

#' Execute Q-interactions step
#' @noRd
.execute_q_interactions_step <- function(analysis, q_vals, condition_col_name, 
                                         methods_to_run, verbose, ...) {
    if (!("q_interactions" %in% methods_to_run)) return(analysis)
    
    if (length(analysis@diversity_results) == 0) {
        if (verbose) message("Step 7: Skipping Q-interactions (requires diversity)")
        return(analysis)
    }
    
    if (verbose) message("Step 7: Detecting Q-dependent interactions...")
    tryCatch({
        analysis <- rank_test_q_condition_s4(analysis, condition_col = condition_col_name,
            q = q_vals, ...)
        if (verbose) message("  [OK] Q-interactions detected")
    }, error = function(e) warning("Q-interactions failed:\n", e$message, call. = FALSE))
    analysis
}

#' Execute jackknife isoform switching step
#' @noRd
.execute_jackknife_isoform_switching_step <- function(analysis, q_vals, condition_col_name, 
                                                       methods_to_run, verbose, ...) {
    if (!("jackknife_isoform_switching" %in% methods_to_run)) return(analysis)
    
    if (length(analysis@diversity_results) == 0) {
        if (verbose) message("Step 3: Skipping isoform switching (requires diversity)")
        return(analysis)
    }
    
    if (verbose) message("Step 3: Computing isoform switching jackknife...")
    tryCatch({
        analysis <- jackknife_isoform_switching_s4(analysis, q = q_vals, 
            condition_col = condition_col_name, verbose = FALSE, ...)
        if (verbose) message("  [OK] Isoform switching jackknife complete")
    }, error = function(e) {
        if (verbose) warning("Isoform switching jackknife skipped: ", e$message, call. = FALSE)
    })
    analysis
}

#' Execute additional LM step helpers
#' @noRd
.execute_effect_sizes_step <- function(analysis, methods_to_run, verbose, ...) {
    if (!("effect_sizes" %in% methods_to_run)) return(analysis)
    
    if (length(analysis@divergence_results) == 0 || length(analysis@lm_results) == 0) {
        if (verbose) message("Step 8: Skipping effect sizes (requires divergence and LM)")
        return(analysis)
    }
    
    if (verbose) message("Step 8: Computing effect sizes...")
    tryCatch({
        analysis <- effect_sizes_divergence_s4(analysis, verbose = FALSE, ...)
        if (verbose) message("  [OK] Effect sizes computed")
    }, error = function(e) {
        if (verbose) warning("Effect size computation skipped: ", e$message, call. = FALSE)
    })
    analysis
}

#' Execute rankbased assumptions testing
#' @noRd
.execute_rankbased_assumptions_step <- function(analysis, q_vals, methods_to_run, verbose, ...) {
    if (!("rankbased_assumptions" %in% methods_to_run)) return(analysis)
    
    if (length(analysis@diversity_results) == 0) {
        if (verbose) message("Step 9: Skipping rankbased assumptions (requires diversity)")
        return(analysis)
    }
    
    if (verbose) message("Step 9: Testing rankbased method assumptions...")
    tryCatch({
        analysis <- test_rankbased_assumptions_s4(analysis, q = q_vals[1], verbose = FALSE, ...)
        if (verbose) message("  [OK] Assumption tests complete")
    }, error = function(e) {
        if (verbose) warning("Assumption testing skipped: ", e$message, call. = FALSE)
    })
    analysis
}

#' Execute method concordance computation
#' @noRd
.execute_method_concordance_step <- function(analysis, methods_to_run, verbose, ...) {
    if (!("method_concordance" %in% methods_to_run)) return(analysis)
    
    if (length(analysis@lm_results) < 2) {
        if (verbose) message("Step 10: Skipping method concordance (requires multiple LM methods)")
        return(analysis)
    }
    
    if (verbose) message("Step 10: Computing method concordance...")
    tryCatch({
        analysis <- compute_method_concordance_s4(analysis, verbose = FALSE, ...)
        if (verbose) message("  [OK] Method concordance computed")
    }, error = function(e) {
        if (verbose) warning("Method concordance skipped: ", e$message, call. = FALSE)
    })
    analysis
}

#' Execute plot generation
#' @noRd
.execute_plot_generation <- function(analysis, do_plots, verbose) {
    if (!do_plots || length(analysis@diversity_results) == 0) return(analysis)
    
    if (verbose) message("Step 11: Generating plots...")
    tryCatch({
        `%||%` <- function(x, y) if (is.null(x)) y else x
        plot_types <- analysis@config$plot_types %||% c("q_curve", "lm_interaction",
            "divergence_distribution", "divergence_spectrum", "influence_heatmap", "volcano",
            "method_concordance", "multi_gene_q_spectrum", "top_transcripts", "tsallis_violin_density")
        
        for (ptype in plot_types) {
            tryCatch({
                plot_obj <- .generate_plot_by_type(ptype, analysis)
                if (!is.null(plot_obj)) {
                    analysis <- addPlot(analysis, type = ptype, plot = plot_obj, replace = TRUE)
                }
            }, error = function(e) {
                if (verbose) warning(sprintf("Plot '%s' failed: %s", ptype, e$message), call. = FALSE)
            })
        }
        if (verbose) message(sprintf("  [OK] %d plot(s) generated", length(analysis@plots)))
    }, error = function(e) warning("Plot generation failed:\n", e$message, call. = FALSE))
    analysis
}

#' Generate plot by type
#' @noRd
.generate_plot_by_type <- function(ptype, analysis) {
    switch(ptype,
        q_curve = if (length(analysis@diversity_results) > 0) {
            plot_tsallis_q_curve_s4(analysis) } else NULL,
        lm_interaction = if (length(analysis@lm_results) > 0 && "lm_interaction" %in% names(analysis@lm_results) &&
                              !is.null(analysis@lm_results$lm_interaction$results) && nrow(analysis@lm_results$lm_interaction$results) > 0) {
            plot_lm_interaction_gam_s4(analysis) } else NULL,
        divergence_distribution = if (!is.null(analysis@metadata$effect_sizes_divergence)) {
            plot_divergence_distribution_s4(analysis) } else NULL,
        divergence_spectrum = if (length(analysis@divergence_results) > 0) {
            plot_divergence_spectrum_s4(analysis) } else NULL,
        influence_heatmap = if (length(analysis@lm_results) > 0 && "q_interactions" %in% names(analysis@lm_results)) {
            plot_multiq_delta_influence_heatmaps_s4(analysis) } else NULL,
        volcano = if (!is.null(analysis@pairwise_results$difference) && length(analysis@divergence_results) > 0) {
            plot_volcano_ma_grid_s4(analysis) } else NULL,
        method_concordance = if (!is.null(analysis@metadata$method_concordance)) {
            plot_method_concordance_s4(analysis) } else NULL,
        multi_gene_q_spectrum = if (length(analysis@lm_results) > 0 && length(analysis@diversity_results) > 0) {
            plot_multi_gene_q_spectrum_s4(analysis) } else NULL,
        top_transcripts = if (length(analysis@diversity_results) > 0 && length(analysis@lm_results) > 0) {
            plot_top_transcripts_s4(analysis) } else NULL,
        tsallis_violin_density = if (length(analysis@diversity_results) > 0) {
            plot_tsallis_violin_density_grid_s4(analysis) } else NULL,
        NULL)
}

#' Track analysis metadata
#' @noRd
.track_analysis_metadata <- function(analysis, methods_run, config) {
    analysis@metadata$workflow <- list(
        steps_completed = methods_run,
        completion_time = Sys.time(),
        tsenat_version = utils::packageVersion("TSENAT")
    )
    analysis@metadata$methods_parameters <- list(
        fdr_threshold = config$fdr_threshold,
        q_values = config$q_values,
        condition_col = config$condition_col %||% "condition"
    )
    analysis
}

#' Validate analysis object structure
#' @noRd
.validate_analysis_object <- function(analysis) {
    checks <- list(
        se_valid = !is.null(analysis@se) && nrow(analysis@se) > 0,
        coldata_valid = all(c("condition") %in% colnames(SummarizedExperiment::colData(analysis@se))),
        min_samples = ncol(analysis@se) >= 2,
        min_genes = nrow(analysis@se) >= 10
    )
    
    if (!all(unlist(checks))) {
        failed <- names(checks)[!unlist(checks)]
        stop("Analysis validation failed: ", paste(failed, collapse = ", "), 
            call. = FALSE)
    }
}

#' Finalize analysis and print summary
#' @noRd
.finalize_tsenat_analysis <- function(analysis, verbose) {
    if (verbose) {
        message("Analysis Complete")
        message("=================")
        message("Results summary:")
        if (length(analysis@diversity_results) > 0) message("  [OK] Diversity")
        if (length(analysis@lm_results) > 0) message("  [OK] LM results")
        if (length(analysis@jackknife_results) > 0) message("  [OK] Jackknife CIs")
        if (length(analysis@divergence_results) > 0) message("  [OK] Divergence")
        if (length(analysis@plots) > 0) message(sprintf("  [OK] Plots (%d)", length(analysis@plots)))
        message("\nUse show(analysis) or summary(analysis) for details")
    }
    analysis@metadata$ended_at <- Sys.time()
    analysis
}

