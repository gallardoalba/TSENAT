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
#' @param analysis \code{TSENATAnalysis} object created by \code{\link{build_analysis}}.
#' @param output_dir \code{character}. Directory to save results and plots.
#'   Default: "tsenat_outputs". Set to NULL to disable automatic output saving.
#' @param save_output \code{logical}. Whether to save output files (results tables).
#'   Default: TRUE. If FALSE, no TSV/CSV output files are written to disk.
#' @param output_format \code{character}. Format for output files: 'tsv' (tab-separated),
#'   'csv' (comma-separated), 'txt' (text), or 'rds' (R serialized). Default: 'tsv'.
#' @param verbose \code{logical}. Print progress messages. Default: TRUE.
#'
#' @return \code{TSENATAnalysis} object containing complete analysis results,
#'   plots, and metadata.
#'
#' @details
#' Pipeline execution order (enforced, follows TSENAT.Rmd vignette):
#' \enumerate{
#'   \item \code{filter_analysis()} - Filter low-abundance transcripts
#'   \item \code{calculate_diversity()} - Tsallis entropy per q-value
#'   \item \code{plot_diversity_spectrum()} - Visualize q-spectrum
#'   \item \code{calculate_m_estimator()} - Sample influence QC analysis
#'   \item \code{calculate_lm()} - LM interaction testing
#'   \item \code{plot_lm_gam()} - GAM visualization of LM results
#'   \item \code{calculate_jis()} - Transcript switching detection
#'   \item \code{plot_jis_delta()} - Multi-q influence heatmap (gene switching tables computed lazily via results())
#'   \item \code{plot_expression()} - Top transcript visualization
#'   \item \code{calculate_divergence()} - Pairwise divergence metrics
#'   \item \code{calculate_effect_sizes()} - Effect size computation
#'   \item \code{plot_divergence_distribution()} - Divergence distribution plot
#'   \item \code{plot_divergence_spectrum()} - Divergence spectrum plot
#' }
#'
#' @examples
#' \donttest{
#' data(readcounts, package = "TSENAT")
#' metadata_df <- read.table(
#'   system.file("extdata", "metadata.tsv", package = "TSENAT"),
#'   header = TRUE, sep = "\t"
#' )
#' gff3_file <- system.file("extdata", "annotation.gff3.gz", package = "TSENAT")
#' 
#' config <- TSENAT_config(
#'   sample_col = "sample",
#'   condition_col = "condition",
#'   q = seq(0, 2, length.out = 10),
#'   generate_plots = FALSE
#' )
#' analysis <- build_analysis(
#'   readcounts = as.matrix(readcounts),
#'   tx2gene = gff3_file,
#'   metadata = metadata_df,
#'   config = config,
#'   tpm = tpm,
#'   effective_length = effective_length
#' )
#' 
#' result <- TSENAT(analysis)
#' }
#'
#' @export
TSENAT <- function(analysis, output_dir = "tsenat_outputs", save_output = TRUE, output_format = "tsv", verbose = TRUE) {
    # Validate input
    if (!is(analysis, "TSENATAnalysis")) {
        stop("'analysis' must be a TSENATAnalysis object created by build_analysis()",
            call. = FALSE)
    }
    if (nrow(se(analysis)) == 0) {
        stop("TSENATAnalysis contains an empty SummarizedExperiment", call. = FALSE)
    }
    
    # Validate output_format
    valid_formats <- c("tsv", "csv", "txt", "rds")
    if (!(output_format %in% valid_formats)) {
        stop("'output_format' must be one of: ", paste(valid_formats, collapse = ", "),
            call. = FALSE)
    }
    
    # Disable output if save_output is FALSE
    if (!save_output) {
        output_dir <- NULL
        if (verbose) message("[INFO] save_output = FALSE prevents file output")
    }

    # Create output directory if specified
    if (!is.null(output_dir) && !dir.exists(output_dir)) {
        dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
        if (verbose) message("Created output directory: ", output_dir)
    }

    # Validate input object
    .validate_analysis_object(analysis)

    # Extract parameters from config (already embedded in analysis object from build_analysis)
    cfg <- getConfig(analysis)
    
    # Get q values from config (can be vector or single value)
    q_vals <- cfg$q %||% 1.0
    if (!is.vector(q_vals)) q_vals <- c(q_vals)
    
    condition_col <- cfg$condition_col %||% "condition"

    # Initialize timing
    workflow_start <- Sys.time()
    step_times <- list()

    # Log pipeline start with configuration
    if (verbose)
        .log_pipeline_start(se(analysis), q_vals, cfg)

    # Execute vignette workflow (in order) with timing
    if (verbose)
        message("=============================================================")
    
    if (verbose)
        message(sprintf("[>] [%2d/14] Filtering low-abundance transcripts", 1))
    step_start <- Sys.time()
    tryCatch({
        cfg <- getConfig(analysis)
        stringency_level <- cfg$stringency %||% "medium"
        analysis <- filter_analysis(analysis, stringency = stringency_level)
        if (verbose)
            message("          [OK] Complete")
    }, error = function(e) warning("Filtering failed:\n", e$message, call. = FALSE))
    step_times[["filtering"]] <- Sys.time() - step_start

    step_start <- Sys.time()
    analysis <- .execute_diversity_s4(analysis, q_vals, verbose, output_dir, output_format)
    step_times[["diversity"]] <- Sys.time() - step_start

    step_start <- Sys.time()
    analysis <- .execute_q_curve_plot(analysis, verbose, output_dir)
    step_times[["q_curve"]] <- Sys.time() - step_start

    step_start <- Sys.time()
    analysis <- .execute_m_estimate_qc(analysis, condition_col, verbose, output_dir, output_format)
    step_times[["m_estimate"]] <- Sys.time() - step_start

    step_start <- Sys.time()
    analysis <- .execute_lm_interaction_s4(analysis, verbose, output_dir, output_format)
    step_times[["lm_interaction"]] <- Sys.time() - step_start

    step_start <- Sys.time()
    analysis <- .execute_lm_interaction_plot(analysis, verbose, output_dir)
    step_times[["lm_plot"]] <- Sys.time() - step_start

    step_start <- Sys.time()
    analysis <- .execute_jackknife_isoform_switching(analysis, q_vals, condition_col,
        verbose, output_dir, output_format)
    step_times[["jackknife"]] <- Sys.time() - step_start

    # Step 8 removed: Gene switching tables now computed lazily via results(type='switching_tables')
    # No need for explicit computation - results() automatically computes and caches when needed
    step_start <- Sys.time()
    analysis <- .execute_influence_heatmap_plot(analysis, verbose, output_dir)
    step_times[["influence_heatmap"]] <- Sys.time() - step_start

    step_start <- Sys.time()
    analysis <- .execute_top_transcripts_plot(analysis, verbose, output_dir)
    step_times[["top_transcripts"]] <- Sys.time() - step_start

    step_start <- Sys.time()
    analysis <- .execute_divergence_s4(analysis, q_vals, verbose, output_dir, output_format)
    step_times[["divergence"]] <- Sys.time() - step_start

    step_start <- Sys.time()
    analysis <- .execute_effect_sizes_s4(analysis, verbose, output_dir, output_format)
    step_times[["effect_sizes"]] <- Sys.time() - step_start

    step_start <- Sys.time()
    analysis <- .execute_divergence_dist_plot(analysis, verbose, output_dir)
    step_times[["div_dist_plot"]] <- Sys.time() - step_start

    step_start <- Sys.time()
    analysis <- .execute_divergence_spectrum_plot(analysis, verbose, output_dir)
    step_times[["div_spectrum_plot"]] <- Sys.time() - step_start
    
    if (verbose)
        message("=============================================================")

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


# Helper function to convert result format
.convert_result_format <- function(result, format, type) {
    if (format == "auto") return(result)
    
    if (format == "dataframe" && !is.data.frame(result)) {
        if (is.matrix(result)) {
            result <- as.data.frame(result)
        }
    } else if (format == "matrix" && !is.matrix(result)) {
        if (is.data.frame(result)) {
            result <- as.matrix(result)
        }
    } else if (format == "list" && !is.list(result)) {
        # Convert to list with rows as elements
        if (is.data.frame(result)) {
            result <- as.list(result)
        } else if (is.matrix(result)) {
            result <- asplit(result, 1)  # Split by rows
        }
    }
    
    return(result)
}

# ============================================================================
# CONFIG BUILDER: TSENAT_config()
# ============================================================================

#' Create and return TSENAT configuration
#'
#' Builds a configuration list for use with \code{\link{TSENAT}}().
#' Allows specifying analysis parameters once and reusing across multiple
#' analyses.
#'
#' @param q \code{numeric}. Q-value(s) for Tsallis entropy (single value or vector). 
#'   Default: 1.0 (Shannon entropy).
#'   Usage: \code{calculate_diversity/divergence} use this for spectrum computation 
#'   (if vector) or as default fallback (if single).
#' @param condition_col \code{character}. Column name in colData containing conditions. Default: 'condition'.
#' @param subject_col \code{character}. Column name in colData containing subject IDs (for paired designs). Default: NULL.
#' @param sample_col \code{character}. Column name in colData containing sample IDs. Default: 'sample'.
#' @param paired \code{logical}. Whether samples are paired/repeated measures. Default: FALSE.
#' @param control \code{character}. Reference/control group label. Default: NULL.
#' @param p_threshold \code{numeric}. Raw p-value threshold. Default: 0.05.
#' @param fdr_threshold \code{numeric}. FDR-adjusted p-value threshold. Default: 0.05.
#' @param significance_threshold \code{numeric}. Significance cutoff. Default: 0.05.
#' @param bootstrap \code{logical}. Enable bootstrap CIs. Default: FALSE.
#' @param nboot \code{integer}. Bootstrap resamples for CIs. Default: 1000.
#' @param bootstrap_method \code{character}. Bootstrap method: 'percentile' or 'bca'. Default: 'percentile'.
#' @param bootstrap_ci \code{numeric}. Confidence level (0-1). Default: 0.95.
#' @param bootstrap_include_diagnostics \code{logical}. Include diagnostics. Default: TRUE.
#' @param min_valid_frac \code{numeric}. Min valid replicate fraction. Default: 0.75.
#' @param pseudocount \code{numeric}. Pseudocount for sparse data. Default: 0.
#' @param norm \code{logical}. Enable normalization. Default: TRUE.
#' @param norm_method \code{character}. Normalization: NULL, 'zscore', 'log_odds_ratio', 'relative_reference'. Default: NULL.
#' @param shrinkage \code{character}. Variance reduction: 'none' or 'empirical_bayes'. Default: 'none'.
#' @param stringency \code{character}. Filtering stringency: 'lenient', 'medium', 'severe'. Default: 'medium'.
#' @param lm_method \code{character}. LM method: 'gam', 'lmm', 'fpca', 'gee'. Default: 'gam'.
#' @param lm_pcorr \code{character}. P-value correction: 'BH', 'bonferroni', 'hochberg', 'holm'. Default: 'BH'.
#' @param jis_use_lm_fdr \code{logical}. Filter jackknife genes using LM p-values. Default: TRUE.
#' @param divergence_ci \code{numeric}. Confidence level for divergence CIs. Default: 0.95.
#' @param nthreads \code{integer}. Parallel threads. Default: 1.
#' @param ... Additional configuration parameters (stored as-is).
#'
#' @return \code{list} with class \code{TSENATConfig} containing all
#'   specified parameters.
#'
#' @details
#' Configuration is stored in the TSENATAnalysis@config slot and used
#' by wrapper functions to configure analysis behavior.
#' Note: Statistical tests (Wilcoxon, shuffle) work on a single q-value,
#' so only one q-value is specified in config.
#'
#' @examples
#' # Default config with standard parameters (point estimates only)
#' cfg <- TSENAT_config()
#'
#' # For Wilcoxon/shuffle tests (single q-value required in config)
#' cfg <- TSENAT_config(
#'   q = 1.0,                          # Shannon entropy - for rank tests
#'   condition_col = 'treatment',
#'   control = 'untreated'
#' )
#'
#' # For Scheirer-Ray-Hare rank tests (multiple q-values)
#' cfg <- TSENAT_config(
#'   q = seq(0, 2, by = 0.5),          # Multiple q-values for spectrum or advanced testing
#'   condition_col = 'treatment',
#'   control = 'untreated'
#' )
#'
#' # With bootstrap CIs for uncertainty quantification (recommended)
#' cfg <- TSENAT_config(
#'   bootstrap = TRUE,                # Enable bootstrap confidence intervals
#'   bootstrap_method = "bca",         # Bias-corrected (better for skewed entropy)
#'   nboot = 1000,                     # 1000 resamples
#'   bootstrap_ci = 0.95               # 95% CI
#' )
#'
#' # Custom with paired analysis, strict filtering, and normalization
#' cfg <- TSENAT_config(
#'   q = 1.0,                          # Shannon entropy
#'   condition_col = 'treatment',
#'   subject_col = 'subject_id',
#'   paired = TRUE,
#'   control = 'untreated',
#'   stringency = 'severe',            # High-confidence transcripts only
#'   norm_method = 'zscore',           # Cross-study standardization
#'   shrinkage = 'none',               # Empirical estimates
#'   bootstrap = TRUE,
#'   bootstrap_method = 'bca',
#'   nboot = 5000,                     # Higher precision
#'   pseudocount = 0,                  # Disabled by default; set > 0 to add pseudocount
#'   significance_threshold = 0.01     # Stricter significance level
#' )
#'
#' @export
TSENAT_config <- function(q = 1.0, condition_col = "condition", subject_col = NULL,
    sample_col = "sample", paired = FALSE, control = NULL, p_threshold = 0.05, fdr_threshold = 0.05,
    significance_threshold = 0.05, bootstrap = FALSE, nboot = 1000, bootstrap_method = "percentile",
    stringency = "medium", nthreads = 1, norm = TRUE, 
    bootstrap_ci = 0.95, bootstrap_include_diagnostics = TRUE, min_valid_frac = 0.75,
    norm_method = NULL, pseudocount = 0, shrinkage = "none", lm_method = "gam",
    lm_pcorr = "BH", jis_use_lm_fdr = TRUE, divergence_ci = 0.95, ...) {
    # Validate q parameter (single or multiple q-values)
    if (is.null(q)) {
        stop("'q' must be specified (q-value or q-values for diversity/statistics calculations).", 
            call. = FALSE)
    }
    if (!is.numeric(q) || any(q < 0) || any(q > 2)) {
        stop("'q' must be numeric value(s) between 0 and 2", call. = FALSE)
    }

    # Validate bootstrap_method
    valid_bootstrap_methods <- c("percentile", "bca")
    if (!bootstrap_method %in% valid_bootstrap_methods) {
        stop("'bootstrap_method' must be 'percentile' or 'bca'", call. = FALSE)
    }

    # Validate lm_method
    valid_lm_methods <- c("gam", "lmm", "fpca", "gee")
    if (!lm_method %in% valid_lm_methods) {
        stop("'lm_method' must be one of: ", paste(valid_lm_methods, collapse = ", "),
            call. = FALSE)
    }

    # Validate lm_pcorr
    valid_lm_pcorr <- c("BH", "bonferroni", "hochberg", "holm")
    if (!lm_pcorr %in% valid_lm_pcorr) {
        stop("'lm_pcorr' must be one of: ", paste(valid_lm_pcorr, collapse = ", "),
            call. = FALSE)
    }

    # Validate divergence_ci
    if (!is.numeric(divergence_ci) || divergence_ci <= 0 || divergence_ci >= 1) {
        stop("'divergence_ci' must be a probability in (0, 1)", call. = FALSE)
    }

    # Build config list with all parameters
    config <- list(
        q = q,
        condition_col = condition_col,
        subject_col = subject_col,
        sample_col = sample_col,
        paired = paired,
        control = control,
        p_threshold = p_threshold,
        fdr_threshold = fdr_threshold,
        significance_threshold = significance_threshold,
        nboot = nboot,
        bootstrap_method = bootstrap_method,
        bootstrap = bootstrap,
        bootstrap_ci = bootstrap_ci,
        bootstrap_include_diagnostics = bootstrap_include_diagnostics,
        min_valid_frac = min_valid_frac,
        stringency = stringency,
        nthreads = nthreads,
        norm = norm,
        norm_method = norm_method,
        pseudocount = pseudocount,
        shrinkage = shrinkage,
        lm_method = lm_method,
        lm_pcorr = lm_pcorr,
        jis_use_lm_fdr = jis_use_lm_fdr,
        divergence_ci = divergence_ci
    )

    # Add any additional parameters (except metadata - should be explicit to build_analysis)
    extra_args <- list(...)
    # Reject metadata in config to enforce Bioconductor pattern (explicit data parameters)
    if (!is.null(extra_args$metadata)) {
        warning("[TSENAT_config] Parameter 'metadata' should not be in config.\n",
                "  Pass metadata directly to build_analysis() as explicit parameter.\n",
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
        
        if (length(missing_paired_params) > 0) {
            warning("[TSENAT_config] Paired design (paired=TRUE) requires complete configuration.\n",
                "  Missing or incomplete parameters: ", paste(missing_paired_params, collapse = ", "), "\n",
                "  This will cause downstream analysis failure or empty results (LM interaction, plotting).\n",
                "  Provide all parameters: \n",
                "    config <- TSENAT_config(\n",
                "      q = 1.0,                              # Q-value for Tsallis entropy\n",
                "      condition_col = 'condition',\n",
                "      subject_col = 'paired_samples',      # Required for paired analysis\n",
                "      paired = TRUE,\n",
                "      control = 'normal'                    # Reference group for comparisons\n",
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
    n_conditions <- length(unique(se[[cfg$condition_col %||% "condition"]]))
    output <- paste0(
        "\n",
        "+============================================================+\n",
        "|          TSENAT: Tsallis Entropy Analysis Toolbox          |\n",
        "+============================================================+\n",
        "                                                              \n",
        "      Science is an essentially anarchic enterprise.          \n",
        "                                                              \n",
        "                       -- Paul Feyerabend, Against Method     \n",
        "                                                              \n",
        "[DATA] Data Summary\n",
        "  Transcripts .......... ", format(nrow(se), big.mark = ","), "\n",
        "  Samples .............. ", ncol(se), "\n",
        "  Conditions ........... ", n_conditions, "\n",
        sprintf("  Q-spectrum range ..... %g to %g (%d values)\n", 
                round(min(q_vals), 2), round(max(q_vals), 2), length(q_vals)),
        "\n[CONFIG] Analysis Configuration\n",
        "  Design ............... ", if (cfg$paired) "paired" else "unpaired", "\n",
        "  Filter stringency .... ", cfg$stringency %||% "medium", "\n",
        "  Normalization ........ ", if (cfg$norm) "enabled [0-1]" else "disabled", "\n",
        "  Normalization method . ", if (cfg$norm) toupper(cfg$norm_method %||% "RANGE (default)") else "N/A", "\n",
        "  Pseudocount .......... ", if (cfg$pseudocount == 0) "disabled" else as.character(cfg$pseudocount), "\n",
        "  Shrinkage ............ ", if (tolower(cfg$shrinkage %||% "none") == "none") "disabled" else toupper(cfg$shrinkage), "\n",
        "  Significance ......... p < ", format(cfg$p_threshold %||% 0.05, nsmall = 3), 
        " | FDR < ", format(cfg$fdr_threshold %||% 0.05, nsmall = 3), "\n",
        "  LM method ............ ", toupper(cfg$lm_method %||% "GAM"), "\n",
        "  LM p-corr method ..... ", toupper(cfg$lm_pcorr %||% "BH"), "\n",
        "  Jackknife use_lm_fdr . ", if (isTRUE(cfg$jis_use_lm_fdr)) "TRUE" else "FALSE", "\n"
    )
    if (isTRUE(cfg$bootstrap)) {
        output <- paste0(output,
            "  Divergence CI ........ ", format(cfg$divergence_ci %||% 0.95, nsmall = 2), "\n"
        )
    }
    if (isTRUE(cfg$bootstrap) && !is.null(cfg$nboot)) {
        output <- paste0(output,
            "  Bootstrap ........... ", cfg$nboot, " x ", toupper(cfg$bootstrap_method %||% "PERCENTILE"),
            " (", format(cfg$bootstrap_ci %||% 0.95, nsmall = 2), " CI)\n"
        )
    }
    output <- paste0(output, "\n")
    message(output)
}

# ============================================================================
# HELPER FUNCTION: Build output filenames
# ============================================================================

#' Build output filename with correct extension
#'
#' Creates output filenames with the appropriate extension based on format.
#'
#' @param base_name \code{character}. Base filename without extension.
#' @param output_dir \code{character}. Output directory path or NULL.
#' @param output_format \code{character}. Output format ('tsv', 'csv', 'txt', 'rds').
#'
#' @return \code{character} Full path to output file, or NULL if output_dir is NULL.
#' @noRd
.build_output_file <- function(base_name, output_dir, output_format) {
    if (is.null(output_dir)) return(NULL)
    ext <- switch(output_format,
        "tsv" = "tsv",
        "csv" = "csv",
        "txt" = "txt",
        "rds" = "rds",
        "tsv"  # default
    )
    filename <- paste0(base_name, ".", ext)
    file.path(output_dir, filename)
}

#' Step 2: Diversity calculation
#' @noRd
.execute_diversity_s4 <- function(analysis, q_vals, verbose, output_dir, output_format) {
    if (verbose)
        message(sprintf("[>] [%2d/14] Computing Tsallis diversity", 2))
    tryCatch({
        output_file <- .build_output_file("diversity_results", output_dir, output_format)
        
        # Extract bootstrap parameters from config to ensure consistency with vignette
        cfg <- getConfig(analysis)
        bootstrap_method <- cfg$bootstrap_method %||% "percentile"
        nboot <- cfg$nboot %||% 1000
        
        # Only show messages if verbose is explicitly TRUE (not during orchestration)
        should_show_messages <- verbose && !is.null(cfg$verbose) && cfg$verbose == TRUE
        
        analysis <- calculate_diversity(
            analysis,
            q = q_vals,
            norm = TRUE,
            bootstrap_method = bootstrap_method,
            nboot = nboot,
            output_file = output_file,
            verbose = FALSE,
            show_messages = should_show_messages
        )
        if (verbose)
            message(sprintf("          [OK] %d q-values processed", length(q_vals)))
    }, error = function(e) warning("Diversity calculation failed:\n", e$message, call. = FALSE))
    analysis
}

#' Step 3: Q-curve plot
#' @noRd
.execute_q_curve_plot <- function(analysis, verbose, output_dir) {
    if (verbose)
        message(sprintf("[>] [%2d/14] Plotting q-spectrum curve", 3))
    tryCatch({
        output_file <- if (!is.null(output_dir)) file.path(output_dir, "q_curve_plot.png") else NULL
        p_qcurve <- plot_diversity_spectrum(analysis, output_file = output_file)
        if (!is.null(p_qcurve)) {
            analysis <- addPlot(analysis, type = "q_curve", plot = p_qcurve, replace = TRUE)
            if (verbose)
                message("          [OK] Plot generated")
        }
    }, error = function(e) {
        if (verbose)
            warning("Q-curve plot failed: ", e$message, call. = FALSE)
    })
    analysis
}

#' Step 4: M-estimate QC analysis
#' @noRd
.execute_m_estimate_qc <- function(analysis, condition_col, verbose, output_dir, output_format) {
    if (verbose)
        message(sprintf("[>] [%2d/14] Running sample influence QC analysis (m-estimator)", 4))
    tryCatch({
        output_file <- .build_output_file("m_estimate_qc", output_dir, output_format)
        analysis <- calculate_m_estimator(analysis, condition_col = condition_col, output_file = output_file)
        if (verbose)
            message("          [OK] M-estimate QC complete")
    }, error = function(e) {
        if (verbose)
            warning("M-estimate QC failed: ", e$message, call. = FALSE)
    })
    analysis
}

#' Step 5: LM interaction testing
#' @noRd
.execute_lm_interaction_s4 <- function(analysis, verbose, output_dir, output_format) {
    if (verbose)
        message(sprintf("[>] [%2d/14] Testing LM interactions with GAM smoother", 5))
    tryCatch({
        cfg <- getConfig(analysis)
        fdr <- cfg$fdr_threshold %||% 0.05
        lm_method <- cfg$lm_method %||% "gam"
        lm_pcorr <- cfg$lm_pcorr %||% "BH"
        output_file <- .build_output_file("lm_interaction_results", output_dir, output_format)
        analysis <- calculate_lm(analysis, fdr_threshold = fdr, method = lm_method,
            pcorr = lm_pcorr, output_file = output_file)
        if (verbose)
            message("          [OK] LM interaction analysis complete")
    }, error = function(e) warning("LM interaction analysis failed:\n", e$message, call. = FALSE))
    analysis
}

#' Step 6: LM interaction GAM plot
#' @noRd
.execute_lm_interaction_plot <- function(analysis, verbose, output_dir) {
    if (verbose)
        message(sprintf("[>] [%2d/14] Plotting LM interaction GAM smoother", 6))
    tryCatch({
        output_file <- if (!is.null(output_dir)) file.path(output_dir, "lm_interaction_gam_plot.png") else NULL
        p_lm <- plot_lm_gam(analysis, output_file = output_file)
        if (!is.null(p_lm)) {
            analysis <- addPlot(analysis, type = "lm_interaction", plot = p_lm, replace = TRUE)
            if (verbose)
                message("          [OK] LM interaction plot generated")
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
    verbose, output_dir, output_format) {
    if (verbose)
        message(sprintf("[>] [%2d/14] Computing jackknife isoform switching analysis", 7))
    tryCatch({
        cfg <- getConfig(analysis)
        jis_use_lm_fdr <- cfg$jis_use_lm_fdr %||% TRUE
        output_file <- .build_output_file("jackknife_isoform_switching", output_dir, output_format)
        analysis <- calculate_jis(analysis, condition_col = condition_col,
            use_lm_fdr = jis_use_lm_fdr, output_file = output_file, verbose = FALSE)
        if (verbose)
            message("          [OK] Jackknife isoform switching complete")
    }, error = function(e) {
        if (verbose)
            warning("Jackknife isoform switching failed: ", e$message, call. = FALSE)
    })
    analysis
}

#' Step 9: Multi-q influence heatmap
#' @noRd
.execute_influence_heatmap_plot <- function(analysis, verbose, output_dir) {
    if (verbose)
        message(sprintf("[>] [%2d/14] Plotting multi-q influence heatmap", 9))
    tryCatch({
        output_file <- if (!is.null(output_dir)) file.path(output_dir, "influence_heatmap.png") else NULL
        p_heatmap <- plot_jis_delta(analysis, output_file = output_file)
        if (!is.null(p_heatmap)) {
            analysis <- addPlot(analysis, type = "influence_heatmap", plot = p_heatmap,
                replace = TRUE)
            if (verbose)
                message("          [OK] Influence heatmap generated")
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
        message(sprintf("[>] [%2d/14] Plotting top transcript counts", 10))
    tryCatch({
        output_file <- if (!is.null(output_dir)) file.path(output_dir, "top_transcripts.png") else NULL
        p_top_tx <- plot_expression(analysis, output_file = output_file)
        if (!is.null(p_top_tx)) {
            analysis <- addPlot(analysis, type = "top_transcripts", plot = p_top_tx,
                replace = TRUE)
            if (verbose)
                message("          [OK] Top transcripts plot generated")
        }
    }, error = function(e) {
        if (verbose)
            warning("Top transcripts plot failed: ", e$message, call. = FALSE)
    })
    analysis
}

#' Step 11: Divergence calculation
#' @noRd
.execute_divergence_s4 <- function(analysis, q_vals, verbose, output_dir, output_format) {
    if (verbose)
        message(sprintf("[>] [%2d/14] Computing divergence metrics", 11))
    tryCatch({
        cfg <- getConfig(analysis)
        divergence_ci <- cfg$divergence_ci %||% 0.95
        output_file <- .build_output_file("divergence_results", output_dir, output_format)
        analysis <- calculate_divergence(analysis, q = q_vals, output_file = output_file,
            verbose = FALSE, ci = divergence_ci)
        if (verbose)
            message("          [OK] Divergence computed")
    }, error = function(e) {
        if (verbose)
            warning("Divergence failed: ", e$message, call. = FALSE)
    })
    analysis
}

#' Step 12: Effect sizes
#' @noRd
.execute_effect_sizes_s4 <- function(analysis, verbose, output_dir, output_format) {
    if (verbose)
        message(sprintf("[>] [%2d/14] Computing effect sizes for divergence", 12))
    tryCatch({
        output_file <- .build_output_file("effect_sizes", output_dir, output_format)
        analysis <- calculate_effect_sizes(analysis, verbose = FALSE, output_file = output_file)
        if (verbose)
            message("          [OK] Effect sizes computed")
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
        message(sprintf("[>] [%2d/14] Plotting divergence distribution", 13))
    tryCatch({
        output_file <- if (!is.null(output_dir)) file.path(output_dir, "divergence_distribution_plot.png") else NULL
        p_div_dist <- plot_divergence_distribution(analysis, output_file = output_file)
        if (!is.null(p_div_dist)) {
            analysis <- addPlot(analysis, type = "divergence_distribution", plot = p_div_dist,
                replace = TRUE)
            if (verbose)
                message("          [OK] Divergence distribution plot generated")
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
        message(sprintf("[>] [%2d/14] Plotting divergence spectrum", 14))
    tryCatch({
        # Plot 1: Global spectrum plot (all genes)
        output_file <- if (!is.null(output_dir)) file.path(output_dir, "divergence_spectrum_plot.png") else NULL
        p_div_spec <- plot_divergence_spectrum(analysis, output_file = output_file)
        if (!is.null(p_div_spec)) {
            analysis <- addPlot(analysis, type = "divergence_spectrum", plot = p_div_spec,
                replace = TRUE)
            if (verbose)
                message("          [OK] Global divergence spectrum plot generated")
        }
        
        # Plot 2: Multi-gene spectrum plot with top 4 genes by p-value
        output_file_multi <- if (!is.null(output_dir)) file.path(output_dir, "divergence_spectrum_plot_top_genes.png") else NULL
        p_multi <- plot_divergence_spectrum(analysis, n_genes = 4, use_pvalue_ranking = TRUE,
            output_file = output_file_multi)
        if (!is.null(p_multi)) {
            analysis <- addPlot(analysis, type = "divergence_spectrum_multi", plot = p_multi,
                replace = TRUE)
            if (verbose)
                message("          [OK] Multi-gene divergence spectrum plot generated")
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
        0.05, q = cfg$q %||% 1.0, condition_col = cfg$condition_col %||%
        "condition", filter_stringency = cfg$filter_stringency %||% "medium")
    analysis
}

#' Finalize analysis and print summary
#' @noRd
.finalize_tsenat_analysis <- function(analysis, verbose, step_times = NULL, total_time = NULL, output_dir = NULL) {
    if (verbose) {
        stats <- .extract_analysis_statistics(analysis)
        
        output <- paste0(
            "\n",
            "+============================================================+\n",
            "|               [OK] ANALYSIS COMPLETE                       |\n",
            "+============================================================+\n\n",
            "[RESULTS] Results Summary\n"
        )
        
        if (stats$n_transcripts > 0)
            output <- paste0(output, sprintf("  [OK] Diversity ........... %d transcripts x %d q-values\n",
                stats$n_transcripts, stats$n_q_values))
        if (stats$n_lm_significant > 0)
            output <- paste0(output, sprintf("  [OK] LM interactions ..... %d genes (p < 0.05)\n",
                stats$n_lm_significant))
        if (stats$n_jackknife > 0)
            output <- paste0(output, sprintf("  [OK] Isoform switching ... %d genes with robust switching\n",
                stats$n_jackknife))
        if (stats$n_divergence > 0)
            output <- paste0(output, sprintf("  [OK] Divergence metrics .. %d genes analyzed\n",
                stats$n_divergence))
        
        if (!is.null(total_time)) {
            time_str <- .format_duration(total_time)
            output <- paste0(output, "\n[PERF] Performance\n")
            output <- paste0(output, sprintf("  Total time ........... %s\n", time_str))
            
            if (!is.null(step_times) && length(step_times) > 3) {
                step_durations <- vapply(step_times, function(x) as.numeric(x, units = "secs"), numeric(1))
                slow_steps <- names(sort(step_durations, decreasing = TRUE))[seq_len(min(3, length(step_durations)))]
                output <- paste0(output, "  Slowest steps:\n")
                for (i in seq_along(slow_steps)) {
                    sname <- slow_steps[i]
                    stime <- step_times[[sname]]
                    pct <- (as.numeric(stime, units = "secs") / as.numeric(total_time, units = "secs")) * 100
                    output <- paste0(output, sprintf("    %d. %-20s %s (%.1f%%)\n", i, sname,
                        .format_duration(stime), pct))
                }
            }
        }
        
        if (!is.null(output_dir) && dir.exists(output_dir)) {
            n_files <- length(list.files(output_dir, recursive = TRUE))
            output <- paste0(output, "\n[OUTPUT] Output\n")
            output <- paste0(output, sprintf("  Directory ........... %s\n", output_dir))
            output <- paste0(output, sprintf("  Files saved ......... %d\n", n_files))
        }
        
        output <- paste0(output,
            "\n[TIPS] Next steps:\n",
            "  show(result)             - View object structure and slots\n",
            "  summary(result)          - Print detailed statistics summary\n",
            "  results(result)          - Extract numerical results (diversity, divergence, etc.)\n",
            "  getPlot(result, type)    - Retrieve specific visualization (e.g., 'diversity', 'volcano')\n",
            "  metadata(result)         - Access metadata and workflow parameters\n",
            "  se(result)               - Get SummarizedExperiment object for downstream analysis\n",
            "\n"
        )
        
        message(output)
    }
    analysis@metadata$ended_at <- Sys.time()
    analysis
}

