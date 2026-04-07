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
#' @param analysis \code{TSENATAnalysis} object created by \code{\link{build_analysis_s4}}.
#' @param config \code{list} or \code{TSENATConfig}. Configuration from
#'   \code{\link{tsenat_config}}. If NULL, uses configuration from analysis object.
#' @param output_dir \code{character}. Directory to save results and plots.
#'   Default: "tsenat_outputs". Set to NULL to disable automatic output saving.
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
#' Pipeline execution order (enforced, follows TSENAT.Rmd vignette):
#' \enumerate{
#'   \item \code{filter_analysis_s4()} - Filter low-abundance transcripts
#'   \item \code{calculate_diversity_s4()} - Tsallis entropy per q-value
#'   \item \code{plot_tsallis_q_curve_s4()} - Visualize q-spectrum
#'   \item \code{m_estimate_s4()} - Sample influence QC analysis
#'   \item \code{calculate_lm_interaction_s4()} - LM interaction testing
#'   \item \code{plot_lm_interaction_gam_s4()} - GAM visualization of LM results
#'   \item \code{jackknife_isoform_switching_s4()} - Transcript switching detection
#'   \item \code{prepare_gene_switching_tables_s4()} - Prepare gene switching summary tables
#'   \item \code{plot_multiq_delta_influence_heatmaps_s4()} - Multi-q influence heatmap
#'   \item \code{plot_top_transcripts_s4()} - Top transcript visualization
#'   \item \code{calculate_divergence_s4()} - Pairwise divergence metrics
#'   \item \code{effect_sizes_divergence_s4()} - Effect size computation
#'   \item \code{plot_divergence_distribution_s4()} - Divergence distribution plot
#'   \item \code{plot_divergence_spectrum_s4()} - Divergence spectrum plot
#' }
#'
#' @examples
#' data(readcounts, package = "TSENAT")
#' metadata_df <- read.table(
#'   system.file("extdata", "metadata.tsv", package = "TSENAT"),
#'   header = TRUE, sep = "\t"
#' )
#' gff3_file <- system.file("extdata", "annotation.gff3.gz", package = "TSENAT")
#' 
#' config <- tsenat_config(
#'   sample_col = "sample",
#'   condition_col = "condition",
#'   q_values = c(0.5, 1.0, 1.5, 2.0, 2.5),
#'   generate_plots = FALSE
#' )
#' analysis <- build_analysis_s4(
#'   readcounts = as.matrix(readcounts),
#'   tx2gene = gff3_file,
#'   metadata = metadata_df,
#'   config = config,
#'   tpm = tpm,
#'   effective_length = effective_length
#' )
#' 
#' result <- tsenat(analysis, config = config)
#'
#' @export
tsenat <- function(analysis, config = NULL, output_dir = "tsenat_outputs", verbose = TRUE, ...) {
    # Validate input
    if (!is(analysis, "TSENATAnalysis")) {
        stop("'analysis' must be a TSENATAnalysis object created by build_analysis_s4()",
            call. = FALSE)
    }
    if (nrow(se(analysis)) == 0) {
        stop("TSENATAnalysis contains an empty SummarizedExperiment", call. = FALSE)
    }

    # Create output directory if specified
    if (!is.null(output_dir) && !dir.exists(output_dir)) {
        dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
        if (verbose) message("Created output directory: ", output_dir)
    }

    # Update config if provided
    if (!is.null(config)) {
        analysis <- setConfig(analysis, config)
    }
    .validate_analysis_object(analysis)

    # Extract parameters
    cfg <- getConfig(analysis)
    q_vals <- cfg$q_values %||% seq(0.5, 2, by = 0.5)
    condition_col <- cfg$condition_col %||% "condition"

    # Initialize timing
    workflow_start <- Sys.time()
    step_times <- list()

    # Log pipeline start with configuration
    if (verbose)
        .log_pipeline_start(se(analysis), q_vals, cfg)

    # Execute vignette workflow (in order) with timing
    if (verbose)
        cat("━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n")
    
    if (verbose)
        cat(sprintf("▶ [%2d/14] Filtering low-abundance transcripts\n", 1))
    step_start <- Sys.time()
    tryCatch({
        analysis <- filter_analysis_s4(analysis)
        if (verbose)
            cat("          ✓ Complete\n")
    }, error = function(e) warning("Filtering failed:\n", e$message, call. = FALSE))
    step_times[["filtering"]] <- Sys.time() - step_start

    step_start <- Sys.time()
    analysis <- .execute_diversity_s4(analysis, q_vals, verbose, output_dir)
    step_times[["diversity"]] <- Sys.time() - step_start

    step_start <- Sys.time()
    analysis <- .execute_q_curve_plot(analysis, verbose, output_dir)
    step_times[["q_curve"]] <- Sys.time() - step_start

    step_start <- Sys.time()
    analysis <- .execute_m_estimate_qc(analysis, condition_col, verbose, output_dir)
    step_times[["m_estimate"]] <- Sys.time() - step_start

    step_start <- Sys.time()
    analysis <- .execute_lm_interaction_s4(analysis, verbose, output_dir)
    step_times[["lm_interaction"]] <- Sys.time() - step_start

    step_start <- Sys.time()
    analysis <- .execute_lm_interaction_plot(analysis, verbose, output_dir)
    step_times[["lm_plot"]] <- Sys.time() - step_start

    step_start <- Sys.time()
    analysis <- .execute_jackknife_isoform_switching(analysis, q_vals, condition_col,
        verbose, output_dir)
    step_times[["jackknife"]] <- Sys.time() - step_start

    step_start <- Sys.time()
    analysis <- .execute_prepare_gene_switching_tables(analysis, verbose, output_dir)
    step_times[["switching_tables"]] <- Sys.time() - step_start

    step_start <- Sys.time()
    analysis <- .execute_influence_heatmap_plot(analysis, verbose, output_dir)
    step_times[["influence_heatmap"]] <- Sys.time() - step_start

    step_start <- Sys.time()
    analysis <- .execute_top_transcripts_plot(analysis, verbose, output_dir)
    step_times[["top_transcripts"]] <- Sys.time() - step_start

    step_start <- Sys.time()
    analysis <- .execute_divergence_s4(analysis, q_vals, verbose, output_dir)
    step_times[["divergence"]] <- Sys.time() - step_start

    step_start <- Sys.time()
    analysis <- .execute_effect_sizes_s4(analysis, verbose, output_dir)
    step_times[["effect_sizes"]] <- Sys.time() - step_start

    step_start <- Sys.time()
    analysis <- .execute_divergence_dist_plot(analysis, verbose, output_dir)
    step_times[["div_dist_plot"]] <- Sys.time() - step_start

    step_start <- Sys.time()
    analysis <- .execute_divergence_spectrum_plot(analysis, verbose, output_dir)
    step_times[["div_spectrum_plot"]] <- Sys.time() - step_start
    
    if (verbose)
        cat("━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n")

    # Track completion metadata and timing
    total_time <- Sys.time() - workflow_start
    analysis <- .track_analysis_metadata(analysis, analysis@config)
    analysis <- .finalize_tsenat_analysis(analysis, verbose, step_times, total_time, output_dir)

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
# HELPER FUNCTIONS FOR REPORTING
# ============================================================================

#' Format duration for display
#' @noRd
.format_duration <- function(duration) {
    seconds <- as.numeric(duration, units = "secs")
    if (seconds < 60) {
        return(sprintf("%.1fs", seconds))
    } else if (seconds < 3600) {
        mins <- seconds / 60
        return(sprintf("%.1fm", mins))
    } else {
        hours <- seconds / 3600
        return(sprintf("%.1fh", hours))
    }
}

#' Extract result statistics from analysis
#' @noRd
.extract_analysis_statistics <- function(analysis) {
    stats <- list(
        n_transcripts = 0,
        n_q_values = 0,
        n_lm_significant = 0,
        n_jackknife = 0,
        n_divergence = 0
    )
    
    # Diversity stats
    if (length(analysis@diversity_results) > 0) {
        div_res <- analysis@diversity_results
        if (is.matrix(div_res)) {
            stats$n_transcripts <- ncol(div_res)
            stats$n_q_values <- nrow(div_res)
        } else if (is.list(div_res)) {
            stats$n_q_values <- length(div_res)
        }
    }
    
    # LM results stats
    if (length(analysis@lm_results) > 0 && is.list(analysis@lm_results)) {
        if (!is.null(analysis@lm_results$pvalue_results)) {
            lm_pvals <- analysis@lm_results$pvalue_results
            if (is.data.frame(lm_pvals) && nrow(lm_pvals) > 0) {
                if ("p_value" %in% colnames(lm_pvals)) {
                    stats$n_lm_significant <- sum(lm_pvals$p_value < 0.05, na.rm = TRUE)
                } else if ("padj" %in% colnames(lm_pvals)) {
                    stats$n_lm_significant <- sum(lm_pvals$padj < 0.05, na.rm = TRUE)
                }
            }
        }
    }
    
    # Jackknife stats
    if (length(analysis@jackknife_results) > 0 && is.list(analysis@jackknife_results)) {
        if (!is.null(analysis@jackknife_results$switching_summary)) {
            stats$n_jackknife <- nrow(analysis@jackknife_results$switching_summary)
        }
    }
    
    # Divergence stats
    if (length(analysis@divergence_results) > 0 && is.data.frame(analysis@divergence_results)) {
        stats$n_divergence <- nrow(analysis@divergence_results)
    }
    
    return(stats)
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
#' data(readcounts, package = 'TSENAT')
#'
#' # Create TSENATAnalysis from count matrix
#' # For simple count matrices (no tx2gene mapping), use TSENATAnalysis directly
#' config <- tsenat_config(
#'   q_values = c(0.5, 1.0, 2.0),
#'   condition_col = 'group'
#' )
#' se <- SummarizedExperiment::SummarizedExperiment(
#'   assays = list(counts = readcounts),
#'   colData = data.frame(
#'     group = rep(c('A', 'B'), length.out = ncol(readcounts))
#'   )
#' )
#' analysis <- TSENATAnalysis(se = se, config = config)
#'
#' # Run analysis to generate diversity results
#' analysis <- calculate_diversity_s4(analysis)
#'
#' # Extract diversity results
#' div_results <- getResults(analysis, type = 'diversity')
#' if (!is.null(div_results)) {
#'   head(div_results, n = 3)
#' }
#'
#' @export
getResults <- function(analysis, type = "diversity", q = NULL, simplify = TRUE) {
    if (!is(analysis, "TSENATAnalysis")) {
        stop("'analysis' must be a TSENATAnalysis object", call. = FALSE)
    }

    result <- switch(type, diversity = if (length(analysis@diversity_results) > 0) analysis@diversity_results else NULL,
        divergence = if (length(analysis@divergence_results) > 0) analysis@divergence_results else NULL,
        lm = if (length(analysis@lm_results) > 0) analysis@lm_results else NULL,
        jackknife = if (length(analysis@jackknife_results) > 0) analysis@jackknife_results else NULL,
        q_interactions = if ("q_interactions" %in% names(analysis@lm_results)) analysis@lm_results$q_interactions else NULL,
        stop("Unknown result type: '", type, "'. Must be one of: ", "diversity, divergence, lm, jackknife, q_interactions",
            call. = FALSE))

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
#' @param nthreads \code{integer}. Number of threads for parallel computation
#'   where supported (diversity, divergence, LM fitting). Default: 1 (no parallelization).
#'   Use 2+ for multi-core systems to improve performance.
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
#'   q_values = seq(0, 2, by = 0.05),  # Recommended for paired designs: 41 values
#'   condition_col = 'treatment',
#'   subject_col = 'subject_id',
#'   paired = TRUE,
#'   control = 'untreated',
#'   methods = c('diversity', 'lm_interaction', 'jackknife', 'jackknife_isoform_switching',
#'               'divergence', 'q_interactions', 'difference', 'effect_sizes'),
#'   bootstrap_method = 'bca',
#'   n_bootstrap = 5000,
#'   significance_threshold = 0.01
#' )
#'
#' @export
tsenat_config <- function(q_values = NULL, condition_col = "condition", subject_col = NULL,
    sample_col = "sample", paired = FALSE, control = NULL, formula = NULL, p_threshold = 0.05, fdr_threshold = 0.05,
    significance_threshold = 0.05, n_bootstrap = 1000, bootstrap_method = "percentile",
    methods = NULL, generate_plots = TRUE, plot_types = NULL, nthreads = 1, ...) {
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
        "divergence", "q_interactions", "difference", "rankbased_assumptions", "method_concordance",
        "effect_sizes")
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
    config <- list(q_values = q_values, condition_col = condition_col, subject_col = subject_col,
        sample_col = sample_col, paired = paired, control = control, p_threshold = p_threshold, fdr_threshold = fdr_threshold,
        significance_threshold = significance_threshold, n_bootstrap = n_bootstrap,
        bootstrap_method = bootstrap_method, methods = methods, generate_plots = generate_plots, nthreads = nthreads)

    # Add optional parameters
    if (!is.null(formula))
        config$formula <- formula
    if (!is.null(plot_types))
        config$plot_types <- plot_types

    # Add any additional parameters (except metadata - should be explicit to build_analysis_s4)
    extra_args <- list(...)
    # Reject metadata in config to enforce Bioconductor pattern (explicit data parameters)
    if (!is.null(extra_args$metadata)) {
        warning("[tsenat_config] Parameter 'metadata' should not be in config.\n",
                "  Pass metadata directly to build_analysis_s4() as explicit parameter.\n",
                "  Bioconductor pattern: data files are explicit, config is for analysis choices.",
                call. = FALSE)
        extra_args$metadata <- NULL
    }
    if (length(extra_args) > 0) {
        config <- c(config, extra_args)
    }

    # Validate paired design configuration (fail-fast principle)
    if (paired == TRUE) {
        missing_paired_params <- c()
        
        if (is.null(subject_col)) {
            missing_paired_params <- c(missing_paired_params, "subject_col")
        }
        if (is.null(control)) {
            missing_paired_params <- c(missing_paired_params, "control")
        }
        if (is.null(q_values) || length(q_values) < 5) {
            missing_paired_params <- c(missing_paired_params, "q_values (recommended: >= 5 values)")
        }
        
        if (length(missing_paired_params) > 0) {
            warning("[tsenat_config] Paired design (paired=TRUE) requires complete configuration.\n",
                "  Missing or incomplete parameters: ", paste(missing_paired_params, collapse = ", "), "\n",
                "  This will cause downstream analysis failure or empty results (LM interaction, plotting).\n",
                "  Provide all parameters: \n",
                "    config <- tsenat_config(\n",
                "      q_values = seq(0, 2, by = 0.05),       # At least 5 q-values (41 recommended)\n",
                "      condition_col = 'condition',\n",
                "      subject_col = 'paired_samples',        # Required for paired analysis\n",
                "      paired = TRUE,\n",
                "      control = 'normal'                     # Reference group for comparisons\n",
                "    )",
                call. = FALSE)
        }
    }

    # Mark as TSENATConfig (but keep as list for S4 slot)
    attr(config, "class") <- c("TSENATConfig", "list")
    config
}

# ============================================================================
# HELPER FUNCTIONS (INTERNAL - NOT EXPORTED)
# ============================================================================

#' Validate analysis object structure
#' @noRd
.validate_analysis_object <- function(analysis) {
    checks <- list(se_valid = !is.null(analysis@se) && nrow(analysis@se) > 0, min_samples = ncol(analysis@se) >=
        2, min_genes = nrow(analysis@se) >= 10)

    if (!all(unlist(checks))) {
        failed <- names(checks)[!unlist(checks)]
        stop("Analysis validation failed: ", paste(failed, collapse = ", "), call. = FALSE)
    }
}

#' Log pipeline start
#' @noRd
.log_pipeline_start <- function(se, q_vals, cfg) {
    # Use ANSI color codes for prettier output
    cat("\n")
    cat("╔════════════════════════════════════════════════════════════╗\n")
    cat("║          TSENAT: Tsallis Entropy Analysis Toolbox          ║\n")
    cat("╚════════════════════════════════════════════════════════════╝\n\n")
    
    cat("📊 Data Summary\n")
    cat("  Transcripts ........... ", format(nrow(se), big.mark = ","), "\n", sep = "")
    cat("  Samples .............. ", ncol(se), "\n", sep = "")
    n_conditions <- length(unique(se[[cfg$condition_col %||% "condition"]]))
    cat("  Conditions ........... ", n_conditions, "\n", sep = "")
    cat("  Q-spectrum range ...... ", format(round(min(q_vals), 2), width = 4), " to ", 
        format(round(max(q_vals), 2), width = 4), " (", length(q_vals), " values)\n", sep = "")
    
    cat("\n⚙️  Configuration\n")
    cat("  p-value threshold ..... ", format(cfg$p_threshold %||% 0.05, width = 6, nsmall = 3), "\n", sep = "")
    cat("  FDR threshold ......... ", format(cfg$fdr_threshold %||% 0.05, width = 6, nsmall = 3), "\n", sep = "")
    if (!is.null(cfg$n_bootstrap)) {
        cat("  Bootstrap samples ..... ", format(cfg$n_bootstrap, big.mark = ","), "\n", sep = "")
    }
    cat("\n")
}

#' Step 2: Diversity calculation
#' @noRd
.execute_diversity_s4 <- function(analysis, q_vals, verbose, output_dir) {
    if (verbose)
        cat(sprintf("▶ [%2d/14] Computing Tsallis diversity\n", 2))
    tryCatch({
        output_file <- if (!is.null(output_dir)) file.path(output_dir, "diversity_results.tsv") else NULL
        suppressMessages({analysis <- calculate_diversity_s4(analysis, q = q_vals, output_file = output_file)})
        if (verbose)
            cat(sprintf("          ✓ %d q-values processed\n", length(q_vals)))
    }, error = function(e) stop("Diversity failed:\n", e$message, call. = FALSE))
    analysis
}

#' Step 3: Q-curve plot
#' @noRd
.execute_q_curve_plot <- function(analysis, verbose, output_dir) {
    if (verbose)
        cat(sprintf("▶ [%2d/14] Plotting q-spectrum curve\n", 3))
    tryCatch({
        output_file <- if (!is.null(output_dir)) file.path(output_dir, "q_curve_plot.png") else NULL
        p_qcurve <- plot_tsallis_q_curve_s4(analysis, output_file = output_file)
        if (!is.null(p_qcurve)) {
            analysis <- addPlot(analysis, type = "q_curve", plot = p_qcurve, replace = TRUE)
            if (verbose)
                cat("          ✓ Plot generated\n")
        }
    }, error = function(e) {
        if (verbose)
            warning("Q-curve plot failed: ", e$message, call. = FALSE)
    })
    analysis
}

#' Step 4: M-estimate QC analysis
#' @noRd
.execute_m_estimate_qc <- function(analysis, condition_col, verbose, output_dir) {
    if (verbose)
        message("Step 4: Running sample influence QC analysis (m-estimator)...")
    tryCatch({
        output_file <- if (!is.null(output_dir)) file.path(output_dir, "m_estimate_qc.tsv") else NULL
        analysis <- m_estimate_s4(analysis, condition_col = condition_col, output_file = output_file)
        if (verbose)
            message("  [OK] M-estimate QC complete")
    }, error = function(e) {
        if (verbose)
            warning("M-estimate QC failed: ", e$message, call. = FALSE)
    })
    analysis
}

#' Step 5: LM interaction testing
#' @noRd
.execute_lm_interaction_s4 <- function(analysis, verbose, output_dir) {
    if (verbose)
        message("Step 5: Testing LM interactions with GAM smoother...")
    tryCatch({
        cfg <- getConfig(analysis)
        fdr <- cfg$fdr_threshold %||% 0.05
        output_file <- if (!is.null(output_dir)) file.path(output_dir, "lm_interaction_results.tsv") else NULL
        analysis <- calculate_lm_interaction_s4(analysis, fdr_threshold = fdr, output_file = output_file)
        if (verbose)
            message("  [OK] LM interaction analysis complete")
    }, error = function(e) stop("LM interaction failed:\n", e$message, call. = FALSE))
    analysis
}

#' Step 6: LM interaction GAM plot
#' @noRd
.execute_lm_interaction_plot <- function(analysis, verbose, output_dir) {
    if (verbose)
        message("Step 6: Plotting LM interaction GAM smoother...")
    tryCatch({
        output_file <- if (!is.null(output_dir)) file.path(output_dir, "lm_interaction_gam_plot.png") else NULL
        p_lm <- plot_lm_interaction_gam_s4(analysis, output_file = output_file)
        if (!is.null(p_lm)) {
            analysis <- addPlot(analysis, type = "lm_interaction", plot = p_lm, replace = TRUE)
            if (verbose)
                message("  [OK] LM interaction plot generated")
        }
    }, error = function(e) {
        if (verbose)
            warning("LM interaction plot failed: ", e$message, call. = FALSE)
    })
    analysis
}

#' Step 7: Jackknife isoform switching
#' @noRd
.execute_jackknife_isoform_switching <- function(analysis, q_vals, condition_col,
    verbose, output_dir) {
    if (verbose)
        message("Step 7: Computing jackknife isoform switching analysis...")
    tryCatch({
        cfg <- getConfig(analysis)
        output_file <- if (!is.null(output_dir)) file.path(output_dir, "jackknife_isoform_switching.tsv") else NULL
        analysis <- jackknife_isoform_switching_s4(analysis, condition_col = condition_col,
            output_file = output_file, verbose = FALSE)
        if (verbose)
            message("  [OK] Jackknife isoform switching complete")
    }, error = function(e) {
        if (verbose)
            warning("Jackknife isoform switching failed: ", e$message, call. = FALSE)
    })
    analysis
}

#' Step 8: Prepare gene switching tables
#' @noRd
.execute_prepare_gene_switching_tables <- function(analysis, verbose, output_dir) {
    if (verbose)
        message("Step 8: Preparing gene switching tables...")
    tryCatch({
        output_file <- if (!is.null(output_dir)) file.path(output_dir, "gene_switching_tables.rds") else NULL
        tables_result <- prepare_gene_switching_tables_s4(analysis, output_file = output_file, verbose = FALSE)
        if (!is.null(tables_result)) {
            if (verbose)
                message("  [OK] Gene switching tables prepared")
        }
    }, error = function(e) {
        if (verbose)
            warning("Gene switching tables failed: ", e$message, call. = FALSE)
    })
    analysis
}

#' Step 9: Multi-q influence heatmap
#' @noRd
.execute_influence_heatmap_plot <- function(analysis, verbose, output_dir) {
    if (verbose)
        message("Step 9: Plotting multi-q influence heatmap...")
    tryCatch({
        output_file <- if (!is.null(output_dir)) file.path(output_dir, "influence_heatmap.png") else NULL
        p_heatmap <- plot_multiq_delta_influence_heatmaps_s4(analysis, output_file = output_file)
        if (!is.null(p_heatmap)) {
            analysis <- addPlot(analysis, type = "influence_heatmap", plot = p_heatmap,
                replace = TRUE)
            if (verbose)
                message("  [OK] Influence heatmap generated")
        }
    }, error = function(e) {
        if (verbose)
            warning("Influence heatmap failed: ", e$message, call. = FALSE)
    })
    analysis
}

#' Step 10: Top transcripts plot
#' @noRd
.execute_top_transcripts_plot <- function(analysis, verbose, output_dir) {
    if (verbose)
        message("Step 10: Plotting top transcript counts...")
    tryCatch({
        output_file <- if (!is.null(output_dir)) file.path(output_dir, "top_transcripts.png") else NULL
        p_top_tx <- plot_top_transcripts_s4(analysis, output_file = output_file)
        if (!is.null(p_top_tx)) {
            analysis <- addPlot(analysis, type = "top_transcripts", plot = p_top_tx,
                replace = TRUE)
            if (verbose)
                message("  [OK] Top transcripts plot generated")
        }
    }, error = function(e) {
        if (verbose)
            warning("Top transcripts plot failed: ", e$message, call. = FALSE)
    })
    analysis
}

#' Step 11: Divergence calculation
#' @noRd
.execute_divergence_s4 <- function(analysis, q_vals, verbose, output_dir) {
    if (verbose)
        message("Step 11: Computing divergence metrics...")
    tryCatch({
        output_file <- if (!is.null(output_dir)) file.path(output_dir, "divergence_results.tsv") else NULL
        analysis <- calculate_divergence_s4(analysis, q = q_vals, output_file = output_file)
        if (verbose)
            message("  [OK] Divergence computed")
    }, error = function(e) {
        if (verbose)
            warning("Divergence failed: ", e$message, call. = FALSE)
    })
    analysis
}

#' Step 12: Effect sizes
#' @noRd
.execute_effect_sizes_s4 <- function(analysis, verbose, output_dir) {
    if (verbose)
        message("Step 12: Computing effect sizes for divergence...")
    tryCatch({
        output_file <- if (!is.null(output_dir)) file.path(output_dir, "effect_sizes.tsv") else NULL
        analysis <- effect_sizes_divergence_s4(analysis, verbose = FALSE, output_file = output_file)
        if (verbose)
            message("  [OK] Effect sizes computed")
    }, error = function(e) {
        if (verbose)
            warning("Effect size computation failed: ", e$message, call. = FALSE)
    })
    analysis
}

#' Step 13: Divergence distribution plot
#' @noRd
.execute_divergence_dist_plot <- function(analysis, verbose, output_dir) {
    if (verbose)
        message("Step 13: Plotting divergence distribution...")
    tryCatch({
        output_file <- if (!is.null(output_dir)) file.path(output_dir, "divergence_distribution_plot.png") else NULL
        p_div_dist <- plot_divergence_distribution_s4(analysis, output_file = output_file)
        if (!is.null(p_div_dist)) {
            analysis <- addPlot(analysis, type = "divergence_distribution", plot = p_div_dist,
                replace = TRUE)
            if (verbose)
                message("  [OK] Divergence distribution plot generated")
        }
    }, error = function(e) {
        if (verbose)
            warning("Divergence distribution plot failed: ", e$message, call. = FALSE)
    })
    analysis
}

#' Step 14: Divergence spectrum plot
#' @noRd
.execute_divergence_spectrum_plot <- function(analysis, verbose, output_dir) {
    if (verbose)
        message("Step 14: Plotting divergence spectrum...")
    tryCatch({
        # Plot 1: Global spectrum plot (all genes)
        output_file <- if (!is.null(output_dir)) file.path(output_dir, "divergence_spectrum_plot.png") else NULL
        p_div_spec <- plot_divergence_spectrum_s4(analysis, output_file = output_file)
        if (!is.null(p_div_spec)) {
            analysis <- addPlot(analysis, type = "divergence_spectrum", plot = p_div_spec,
                replace = TRUE)
            if (verbose)
                message("  [OK] Global divergence spectrum plot generated")
        }
        
        # Plot 2: Multi-gene spectrum plot with top 4 genes by p-value
        output_file_multi <- if (!is.null(output_dir)) file.path(output_dir, "divergence_spectrum_plot_top_genes.png") else NULL
        p_multi <- plot_divergence_spectrum_s4(analysis, n_genes = 4, use_pvalue_ranking = TRUE,
            output_file = output_file_multi)
        if (!is.null(p_multi)) {
            analysis <- addPlot(analysis, type = "divergence_spectrum_multi", plot = p_multi,
                replace = TRUE)
            if (verbose)
                message("  [OK] Multi-gene divergence spectrum plot generated")
        }
    }, error = function(e) {
        if (verbose)
            warning("Divergence spectrum plot failed: ", e$message, call. = FALSE)
    })
    analysis
}

#' Track analysis metadata
#' @noRd
.track_analysis_metadata <- function(analysis, config) {
    cfg <- getConfig(analysis)
    analysis@metadata$workflow <- list(workflow_type = "isoform_switching_vignette",
        completion_time = Sys.time(), tsenat_version = utils::packageVersion("TSENAT"))
    analysis@metadata$methods_parameters <- list(fdr_threshold = cfg$fdr_threshold %||%
        0.05, q_values = cfg$q_values %||% seq(0.5, 2, by = 0.5), condition_col = cfg$condition_col %||%
        "condition", filter_stringency = cfg$filter_stringency %||% "medium")
    analysis
}

#' Finalize analysis and print summary
#' @noRd
.finalize_tsenat_analysis <- function(analysis, verbose, step_times = NULL, total_time = NULL, output_dir = NULL) {
    if (verbose) {
        cat("\n")
        cat("╔════════════════════════════════════════════════════════════╗\n")
        cat("║               ✓ ANALYSIS COMPLETE                          ║\n")
        cat("╚════════════════════════════════════════════════════════════╝\n\n")
        
        # Extract statistics
        stats <- .extract_analysis_statistics(analysis)
        
        cat("📈 Results Summary\n")
        
        if (stats$n_transcripts > 0)
            cat(sprintf("  ✓ Diversity ........... %d transcripts × %d q-values\n", 
                stats$n_transcripts, stats$n_q_values))
        if (stats$n_lm_significant > 0)
            cat(sprintf("  ✓ LM interactions ..... %d genes (p < 0.05)\n", 
                stats$n_lm_significant))
        if (stats$n_jackknife > 0)
            cat(sprintf("  ✓ Isoform switching ... %d genes\n", 
                stats$n_jackknife))
        if (stats$n_divergence > 0)
            cat(sprintf("  ✓ Divergence metrics .. %d pairwise comparisons\n", 
                stats$n_divergence))
        if (!is.null(analysis@metadata$effect_sizes_divergence))
            cat("  ✓ Effect sizes ........ computed\n")
        if (length(analysis@plots) > 0)
            cat(sprintf("  ✓ Visualizations ...... %d plots\n", length(analysis@plots)))
        
        # Show timing
        if (!is.null(total_time)) {
            time_str <- .format_duration(total_time)
            cat("\n⏱️  Performance\n")
            cat(sprintf("  Total time ........... %s\n", time_str))
            
            if (!is.null(step_times) && length(step_times) > 3) {
                step_durations <- sapply(step_times, function(x) as.numeric(x, units = "secs"))
                slow_steps <- names(sort(step_durations, decreasing = TRUE))[1:min(3, length(step_durations))]
                cat("  Slowest steps:\n")
                for (i in seq_along(slow_steps)) {
                    sname <- slow_steps[i]
                    stime <- step_times[[sname]]
                    pct <- (as.numeric(stime, units = "secs") / as.numeric(total_time, units = "secs")) * 100
                    cat(sprintf("    %d. %-20s %s (%.1f%%)\n", i, sname, 
                        .format_duration(stime), pct))
                }
            }
        }
        
        # Show output
        if (!is.null(output_dir) && dir.exists(output_dir)) {
            n_files <- length(list.files(output_dir, recursive = TRUE))
            cat("📁 Output\n")
            cat(sprintf("  Directory ........... %s\n", output_dir))
            cat(sprintf("  Files saved ......... %d\n", n_files))
        }
        
        cat("\n")
        cat("💡 Next steps:\n")
        cat("  show(analysis)      - View object structure and slots\n")
        cat("  summary(analysis)   - Print detailed statistics\n")
        cat("  getPlot(analysis)   - Extract visualization results\n")
        cat("\n")
    }
    analysis@metadata$ended_at <- Sys.time()
    analysis
}

