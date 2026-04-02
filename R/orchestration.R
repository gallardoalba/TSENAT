# Unified TSENAT Analysis Orchestration Main entry point and configuration for
# TSENAT pipeline. Coordinates analysis workflow from raw counts to results and
# visualizations.

# ============================================================================
# CONFIG BUILDER
# ============================================================================

#' Create and return TSENAT configuration
#'
#' Builds a configuration list for use with \code{\link{tsenat}}().
#' Allows specifying analysis parameters once and reusing across multiple
#' analyses.
#'
#' @param q_values \code{numeric}. Q-values for Tsallis entropy spectrum.
#'   Default: \code{seq(0.5, 2.0, by = 0.5)}.
#' @param q_range \code{numeric}.  Alternative to q_values:  lower and 
#' upper bounds.
#' @param filter_genome \code{logical}. Remove zero rows before analysis.
#'   Default: TRUE.
#' @param formula \code{formula}.  Model formula for  LM interactions (e. g. ,
#'  \code{~ treatment}).
#' @param p_threshold \code{numeric}. P-value threshold for significance.
#'   Default: 0.05.
#' @param fdr_threshold \code{numeric}. FDR threshold (Benjamini-Hochberg).
#'   Default: 0.05.
#' @param methods \code{character}. Analysis methods to run. Options:
#'   'diversity', 'lm_interaction', 'jackknife', 'divergence',
#'   'q_interactions', 'difference'. Default: all methods.
#' @param generate_plots \code{logical}. Generate visualizations.
#'   Default: TRUE.
#' @param plot_types \code{character}. Specific plots to generate.
#'   Default: all available types.
#' @param seed \code{numeric}. Random seed for reproducibility.
#' @param condition_col \code{character}. Name of condition column in colData.
#'   Default: 'condition'.
#' @param ... Additional configuration parameters (stored as-is).
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
#' # Custom spectrum and formula
#' cfg <- tsenat_config(
#'   q_values = c(0.5, 1.0, 1.5, 2.0),
#'   formula = ~ treatment + batch,
#'   fdr_threshold = 0.01
#' )
#'
#' @export
tsenat_config <- function(q_values = NULL, q_range = NULL, filter_genome = TRUE,
    formula = NULL, p_threshold = 0.05, fdr_threshold = 0.05, methods = NULL, generate_plots = TRUE,
    plot_types = NULL, seed = NULL, condition_col = "condition", ...) {
    # Build q_values if range specified
    if (!is.null(q_range)) {
        if (length(q_range) != 2) {
            stop("'q_range' must be c(lower, upper)", call. = FALSE)
        }
        if (!is.finite(q_range[1]) || !is.finite(q_range[2])) {
            stop("'q_range' values must be finite", call. = FALSE)
        }
        q_values <- seq(q_range[1], q_range[2], by = 0.5)
    }

    # Default q_values
    if (is.null(q_values)) {
        q_values <- seq(0.5, 2, by = 0.5)
    }

    # Default methods
    if (is.null(methods)) {
        methods <- c("diversity", "lm_interaction", "jackknife", "divergence", "q_interactions")
    }

    # Validate methods
    valid_methods <- c("diversity", "lm_interaction", "jackknife", "divergence",
        "q_interactions", "difference")
    invalid_methods <- setdiff(methods, valid_methods)
    if (length(invalid_methods) > 0) {
        stop("Invalid methods: ", paste(invalid_methods, collapse = ", "), "\n",
            "Valid: ", paste(valid_methods, collapse = ", "), call. = FALSE)
    }

    # Build config list
    config <- list(q_values = q_values, filter_genome = filter_genome, p_threshold = p_threshold,
        fdr_threshold = fdr_threshold, methods = methods, generate_plots = generate_plots,
        condition_col = condition_col)

    # Add optional parameters
    if (!is.null(formula))
        config$formula <- formula
    if (!is.null(plot_types))
        config$plot_types <- plot_types
    if (!is.null(seed))
        config$seed <- seed

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
# MAIN ORCHESTRATION FUNCTION
# ============================================================================

#' Run complete TSENAT analysis pipeline
#'
#' Coordinates the full TSENAT workflow: diversity -> jackknife -> LM
#' interactions ->
#' divergence -> gene interactions -> visualizations.
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
#'   \item \code{calculate_lm_interaction_s4()} - Statistical tests
#'   \item \code{calculate_divergence_s4()} - Pairwise divergence metrics
#'   \item \code{rank_test_q_condition_s4()} - Q-dependent interactions
#'   \item Plot generation (if enabled)
#' }
#'
#' Metadata automatically tracks:
#' - Analysis start/end time
#' - TSENAT version
#' - Function execution sequence
#' - Parameter settings
#'
#' @examples
#' library(SummarizedExperiment)
#' # Create minimal SummarizedExperiment
#' se <- SummarizedExperiment(
#'   assays = list(counts = matrix(rpois(100, 10), nrow = 10, ncol = 10)),
#' colData = data.frame(sample_id = paste0('S', 1:10), condition =
#' rep(c('A', 'B'), 5))
#' )
#' cfg <- tsenat_config(q_values = c(0.5, 1.0), generate_plots = FALSE)
#' analysis <- TSENATAnalysis(se, config = cfg)
#' show(analysis)
#'
#' @export
tsenat <- function(se, config = NULL, methods = NULL, q_values = NULL, generate_plots = TRUE,
    verbose = TRUE, parallel = FALSE, ...) {
    # Validate input
    if (!is(se, "SummarizedExperiment")) {
        stop("'se' must be a SummarizedExperiment object", call. = FALSE)
    }

    if (nrow(se) == 0) {
        stop("SummarizedExperiment is empty (0 genes)", call. = FALSE)
    }

    # Helper for null coalescing
    `%||%` <- function(x, y) if (is.null(x))
        y else x

    # Initialize TSENATAnalysis object
    analysis <- TSENATAnalysis(se = se, config = config)

    # Merge/override config
    if (!is.null(methods))
        analysis@config$methods <- methods
    if (!is.null(q_values))
        analysis@config$q_values <- q_values

    # Extract parameters from config with defaults
    methods_to_run <- analysis@config$methods %||% c("diversity", "lm_interaction",
        "jackknife", "divergence", "q_interactions")
    q_vals <- analysis@config$q_values %||% seq(0.5, 2, by = 0.5)
    condition_col_name <- analysis@config$condition_col %||% "condition"
    do_plots <- generate_plots && (analysis@config$generate_plots %||% TRUE)
    do_parallel <- parallel && ("parallel" %in% rownames(utils::installed.packages()))

    # Validate method dependencies (Gap 10A improvement)
    validate_method_dependencies <- function(methods_requested) {
        dependencies <- list(jackknife = "diversity", divergence = "diversity", q_interactions = "diversity",
            lm_interaction = "diversity")

        for (method in methods_requested) {
            if (method %in% names(dependencies)) {
                required <- dependencies[[method]]
                if (!(required %in% methods_requested)) {
                  stop("Method '", method, "' requires '", required, "' to be in methods list.\n",
                    "Add '", required, "' to methods parameter or remove '", method,
                    "'.", call. = FALSE)
                }
            }
        }
    }

    # Validate requested methods
    validate_method_dependencies(methods_to_run)

    # Log start
    if (verbose) {
        message("TSENAT Pipeline")
        message("===============")
        message(sprintf("Genes:   %d", nrow(se)))
        message(sprintf("Samples: %d", ncol(se)))
        message(sprintf("Methods: %s", paste(methods_to_run, collapse = ", ")))
        message(sprintf("Q-values: %s", paste(q_vals, collapse = ", ")))
    }

    # ========== STEP 1: DIVERSITY ==========
    if ("diversity" %in% methods_to_run) {
        if (verbose)
            message("Step 1: Calculating diversity...")

        tryCatch({
            analysis <- calculate_diversity_s4(analysis, q = q_vals, ...)
            if (verbose)
                message(sprintf("  [OK] Diversity calculated for q = %s", paste(q_vals,
                  collapse = ", ")))
        }, error = function(e) {
            stop("Diversity calculation failed:\n", e$message, call. = FALSE)
        })
    }

    # ========== STEP 2: JACKKNIFE ==========
    if ("jackknife" %in% methods_to_run) {
        if (length(analysis@diversity_results) == 0) {
            if (verbose)
                message("Step 2: Skipping jackknife (requires diversity)")
        } else {
            if (verbose)
                message("Step 2: Running jackknife resampling...")

            tryCatch({
                analysis <- jackknife_entropy_outliers_s4(analysis, q = q_vals, verbose = FALSE,
                  ...)
                if (verbose)
                  message("  [OK] Jackknife CIs computed")
            }, error = function(e) {
                warning("Jackknife failed:\n", e$message, call. = FALSE)
            })
        }
    }

    # ========== STEP 3: LM INTERACTIONS ==========
    if ("lm_interaction" %in% methods_to_run) {
        if (verbose)
            message("Step 3: Testing LM interactions...")

        tryCatch({
            analysis <- calculate_lm_interaction_s4(analysis, fdr_threshold = analysis@config$fdr_threshold %||%
                0.05, ...)
            if (verbose)
                message("  [OK] LM analysis complete")
        }, error = function(e) {
            warning("LM interaction calculation failed:\n", e$message, call. = FALSE)
        })
    }

    # ========== STEP 4: DIVERGENCE ==========
    if ("divergence" %in% methods_to_run) {
        if (length(analysis@diversity_results) == 0) {
            if (verbose)
                message("Step 4: Skipping divergence (requires diversity)")
        } else {
            if (verbose)
                message("Step 4: Calculating divergence metrics...")

            tryCatch({
                analysis <- calculate_divergence_s4(analysis, q = q_vals[1])
                if (verbose)
                  message("  [OK] Divergence metrics computed")
            }, error = function(e) {
                warning("Divergence calculation failed:\n", e$message, call. = FALSE)
            })
        }
    }

    # ========== STEP 5: Q-DEPENDENT INTERACTIONS ==========
    if ("q_interactions" %in% methods_to_run) {
        if (length(analysis@diversity_results) == 0) {
            if (verbose)
                message("Step 5: Skipping Q-interactions (requires diversity)")
        } else {
            if (verbose)
                message("Step 5: Detecting Q-dependent interactions...")

            tryCatch({
                analysis <- rank_test_q_condition_s4(analysis, condition_col = condition_col_name,
                  q = q_vals, ...)
                if (verbose)
                  message("  [OK] Q-interactions detected")
            }, error = function(e) {
                warning("Q-interaction detection failed:\n", e$message, call. = FALSE)
            })
        }
    }

    # ========== STEP 6: PLOT GENERATION ==========
    if (do_plots && length(analysis@diversity_results) > 0) {
        if (verbose)
            message("Step 6: Generating plots...")

        tryCatch({
            # Plot types from config or auto-detect
            plot_types <- analysis@config$plot_types %||% c("q_curve", "lm_interaction",
                "divergence_distribution", "divergence_spectrum", "influence_heatmap",
                "volcano")

            for (ptype in plot_types) {
                tryCatch({
                  # Dispatch to appropriate plot function based on type
                  plot_obj <- switch(ptype, q_curve = plot_tsallis_q_curve_s4(analysis@se,
                    analysis@diversity_results), lm_interaction = if ("lm_interaction" %in%
                    names(analysis@lm_results)) {
                    .plot_lm_interaction_gam(analysis@lm_results$lm_interaction)
                  } else NULL, divergence_distribution = if (length(analysis@divergence_results) >
                    0) {
                    .plot_divergence_distribution(analysis@divergence_results)
                  } else NULL, divergence_spectrum = if (length(analysis@divergence_results) >
                    0) {
                    .plot_divergence_spectrum(analysis@divergence_results)
                  } else NULL, influence_heatmap = if ("q_interactions" %in% names(analysis@lm_results)) {
                    .plot_multiq_delta_influence_heatmaps(analysis@lm_results$q_interactions)
                  } else NULL, volcano = if ("lm_interaction" %in% names(analysis@lm_results)) {
                    .plot_volcano_ma_grid(analysis@lm_results$lm_interaction, analysis@divergence_results)
                  } else NULL, NULL)

                  if (!is.null(plot_obj)) {
                    analysis <- addPlot(analysis, type = ptype, plot = plot_obj,
                      replace = TRUE)
                  }
                }, error = function(e) {
                  if (verbose) {
                    message(sprintf("  [WARNING] Plot '%s' failed: %s", ptype, e$message))
                  }
                })
            }

            if (verbose) {
                message(sprintf("  [OK] %d plot(s) generated", length(analysis@plots)))
            }
        }, error = function(e) {
            warning("Plot generation failed:\n", e$message, call. = FALSE)
        })
    }

    # ========== FINALIZE ==========
    if (verbose) {
        message("Analysis Complete")
        message("=================")
        message("Results summary:")
        if (length(analysis@diversity_results) > 0)
            message("  [OK] Diversity")
        if (length(analysis@lm_results) > 0)
            message("  [OK] LM results")
        if (length(analysis@jackknife_results) > 0)
            message("  [OK] Jackknife CIs")
        if (length(analysis@divergence_results) > 0)
            message("  [OK] Divergence")
        if (length(analysis@plots) > 0)
            message(sprintf("  [OK] Plots (%d)", length(analysis@plots)))
        message("\nUse show(analysis) or summary(analysis) for details")
    }

    # Add final timing
    analysis@metadata$ended_at <- Sys.time()

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
