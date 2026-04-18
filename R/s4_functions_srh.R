#' Detect q-dependent gene interactions
#'
#' Wrapper around [.calculate_srh()] that manages TSENATAnalysis object.
#' Tests for genes with condition-specific q-dependent entropy patterns by
#' testing whether the effect of q-values DIFFERS between experimental conditions.
#' This detects disease-relevant or condition-specific isoform switching patterns.
#'
#' ## Key Features
#'
#' - **Q\eqn{\times} Condition Interaction**: Tests if entropy patterns across q-values differ
#'   by condition (main discovery goal)
#' - **Multi-q Analysis**: Combines diversity results for multiple q-values into
#'   a single SummarizedExperiment for joint hypothesis testing
#' - **Rank-Based Statistics**: Scheirer-Ray-Hare test (two-way ANOVA on ranked data,
#' - **Scheirer-Ray-Hare Test**: Two-way non-parametric ANOVA on ranks
#' - **Multiple Testing Correction**: Hochberg, Benjamini-Yekutieli, or permutation
#'   (Westfall-Young) procedures
#' - **AR(1) Correlation Handling**: Westfall-Young preserves q-value spatial
#'   correlations (important for ordered q measurements)
#' - **Effect Sizes**: Eta-squared (\eqn{\eta^2}) for q\eqn{\times} condition interactions
#'
#' ## Statistical Hypotheses
#'
#' Tests the null hypothesis:
#' - \strong{H\emph{0}} = Gene entropy q-effect does NOT differ between conditions (q-independent)
#' - \strong{H\emph{1}} = Gene entropy q-dependence is CONDITION-SPECIFIC (interaction exists)
#'
#' A significant interaction indicates condition-specific patterns in how entropy
#' varies across the q-value spectrum, revealing biological processes specific to
#' that condition.
#'
#' ## Biological Example
#'
#' Gene shows strong isoform switching (q-dependent entropy) in tumor cells but
#' NOT in healthy cells -> Identified as disease-relevant q-dependent gene.
#'
#' For condition-specific q-dependent genes:
#' - **Condition A**: Strong entropy variation across q (q-dependent isoform usage)
#' - **Condition B**: Flat entropy profile across q (uniform isoform usage)
#' - **Interaction**: Condition-specific q-dependence pattern reveals disease-associated
#'   splicing regulation
#'
#' @param analysis \code{TSENATAnalysis} object.
#'  Must have diversity results from \code{calculate_diversity()}.
#' @param output_file \code{character} or  \code{NULL}.
#'  Optional file path to save results.
#'   Supported formats: .rds (for S4 objects). Default: NULL (no file output).
#' @param paired \code{logical} or  \code{NULL}.  If TRUE,
#'  uses paired/blocked design 
#'   (requires \code{subject_col}). If NULL, reads from \code{@config$paired}.
#' @param subject_col \code{character} or  \code{NULL}.  Column name for 
#' subject/block identifiers
#'   (required when  \code{paired=TRUE}).  If NULL,
#'  reads from \code{@config$subject_col}.
#' @param condition_col \code{character}.  Column name for 
#' sample grouping/condition (REQUIRED).
#' Specifies the condition/treatment variable for testing q\eqn{\times} condition
#' interactions.
#'   Example: 'sample_type', 'treatment', 'disease_status'.
#' @param multicorr \code{character}.  Multiple testing correction:
#'  'hochberg' (default),
#'   'benjamini-yekutieli', 'westfall-young', or 'none'.
#' @param entropy_col \code{character}.
#'  Column name containing entropy/diversity data.
#'   Default: 'diversity'.
#' @param q_col \code{character}. Column name containing q-values. Default: 'q'.
#' @param gene_col \code{character}.  Column name containing gene identifiers.
#'  Default:  'gene'.
#' @param wy_randomizations \code{numeric} or  \code{character}.
#'  Number of permutations for 
#'   Westfall-Young correction. Use 'auto' to estimate from data. Default: 500.
#' @param nperm_mode \code{character}.  Mode for 
#' automatic permutation estimation:
#'   'standard' (default), 'conservative', or 'interactive'.
#' @param nthreads \code{numeric} or  \code{NULL}.
#'  Number of parallel threads for  computation.
#'   If NULL, reads from \code{@config$nthreads}.
#' @param verbose \code{logical}.  If TRUE,  prints progress messages.
#'  Default:  FALSE.
#' @param alpha \code{numeric}. Significance level for p-value correction methods
#'  (default: 0.05). Used by all multiple testing correction methods.
#' @param p_threshold \code{numeric}. P-value threshold for classification of 
#'  interaction significance (default: 0.05).
#' @param eta2_threshold_moderate \code{numeric}. Effect size boundary for 'moderate'
#'  classification (default: 0.01).
#' @param eta2_threshold_strong \code{numeric}. Effect size boundary for 'strong'
#'  classification (default: 0.10).
#' @param min_nperm \code{integer}. Minimum permutations for automatic estimation
#'  when wy_randomizations='auto' (default: 100).
#' @param max_nperm \code{integer}. Maximum permutations for automatic estimation
#'  when wy_randomizations='auto' (default: 10000).
#' @param ... Additional arguments passed to the base \code{.calculate_srh()} function.
#'
#' @return Modified TSENATAnalysis with interaction results in @rrm_results.
#'
#' @details
#' Analyzes how gene interactions change across q-value spectrum using
#' rank-based
#' (Scheirer-Ray-Hare) or parametric (GAM) statistical tests.
#'
#' **Parameter resolution priority** (explicit > @config > default/auto-detect):
#' \itemize{
#'   \item \code{condition_col}: REQUIRED - must be explicitly provided
#'   \item \code{q}: ALWAYS auto-detected from diversity_results (all q-values tested together)
#'   \item \code{paired}: explicit arg > \code{@config$paired} > FALSE (default)
#'   \item \code{subject_col}: explicit arg > \code{@config$subject_col}
#'   \item \code{multicorr}:
#'  explicit arg > \code{@config$multicorr} > 'hochberg'
#'   \item \code{nthreads}: explicit arg > \code{@config$nthreads} > 1 (default)
#'   \item \code{test}:
#'  explicit arg > \code{@config$test} > 'auto' (auto-selection)
#'   \item \code{nperm_mode}:
#'  explicit arg > \code{@config$nperm_mode} > 'standard'
#' }
#'
#' @examples
#' # Load example data (matching TSENAT.Rmd workflow)
#' data(readcounts)
#' readcounts <- as.matrix(readcounts)
#' mode(readcounts) <- 'numeric'
#' metadata_df <- read.table(
#'   system.file('extdata', 'metadata.tsv', package = 'TSENAT'),
#'   header = TRUE, sep = '\t'
#' )
#' gff3_dataset <- system.file('extdata', 'annotation.gff3.gz', package =
#' 'TSENAT')
#' 
#' # Create config first (required when metadata is provided)
#' config <- TSENAT_config(sample_col = 'sample', condition_col = 'condition')
#' 
#' # Build analysis from vignette data and create manageable subset
#' analysis <- build_analysis(
#'   readcounts = readcounts,
#'   tx2gene = gff3_dataset,
#'   metadata = metadata_df,
#'   config = config,
#'   tpm = tpm,
#'   effective_length = effective_length
#' )
#' analysis <- filter_analysis(
#'   analysis,
#'   min_samples = 1,
#'   subset_n_genes = 200
#' )
#' analysis <- calculate_diversity(analysis, q = c(0.5, 1.0, 1.5))
#' 
#' # Test Q\eqn{\times} Condition interaction (condition_col is REQUIRED)
#' analysis <- calculate_srh(
#'   analysis,
#'   condition_col = 'condition',
#'   multicorr = 'hochberg'
#' )
#' # View results using unified accessor
#' rank_test_res <- results(analysis, type = 'rank_test')
#' if (!is.null(rank_test_res)) head(rank_test_res)
#'
#' @export
#' @importFrom utils write.table
calculate_srh <- function(analysis, condition_col, output_file = NULL, paired = NULL,
    subject_col = NULL, multicorr = c("hochberg", "benjamini-yekutieli", "westfall-young",
        "none"), entropy_col = "diversity", q_col = "q", gene_col = "gene", wy_randomizations = 500,
    nperm_mode = c("standard", "conservative", "interactive"), nthreads = NULL, alpha = 0.05,
    p_threshold = 0.05, eta2_threshold_moderate = 0.01, eta2_threshold_strong = 0.1,
    min_nperm = 100, max_nperm = 10000, verbose = FALSE, ...) {

    # PHASE 1: Validate input and prerequisites
    condition_col <- .validate_srh_input(analysis, condition_col, verbose)

    # PHASE 2: Resolve parameters from config + explicit args Note: q-values
    # are ALWAYS auto-detected from diversity_results
    param_result <- .resolve_srh_params(analysis, multicorr, nperm_mode, paired,
        subject_col, nthreads, wy_randomizations, entropy_col, q_col, gene_col)
    dots <- param_result$dots
    dots$condition_col <- condition_col
    dots$verbose <- verbose

    # Add exposed statistical parameters
    dots$alpha <- alpha
    dots$p_threshold <- p_threshold
    dots$eta2_threshold_moderate <- eta2_threshold_moderate
    dots$eta2_threshold_strong <- eta2_threshold_strong
    dots$min_nperm <- min_nperm
    dots$max_nperm <- max_nperm

    # PHASE 3: Prepare multi-Q SummarizedExperiment
    se_multi_q <- .prepare_multi_q_se(analysis)

    # PHASE 4: Run core rank-based testing
    result <- tryCatch({
        do.call(.calculate_srh, c(list(data = se_multi_q), dots))
    }, error = function(e) {
        stop("q-interaction detection failed:\n", e$message, call. = FALSE)
    })

    # PHASE 5: Store results and save if requested
    analysis <- .store_srh_results(analysis, result, output_file, verbose)

    analysis
}

#' Internal: Validate SRH input and prerequisites
#'
#' @noRd
.validate_srh_input <- function(analysis, condition_col, verbose = FALSE) {
    if (!is(analysis, "TSENATAnalysis")) {
        stop("'analysis' must be a TSENATAnalysis object", call. = FALSE)
    }

    # condition_col is REQUIRED for Q x Condition interaction testing (Q main effect support removed March 2026)
    if (missing(condition_col) || is.null(condition_col)) {
        if (!is.null(analysis@config) && "condition_col" %in% names(analysis@config)) {
            condition_col <- analysis@config$condition_col
            if (verbose)
                message("Using condition_col='", condition_col, "' from @config")
        } else {
            stop("'condition_col' is REQUIRED for Q x Condition interaction testing. ",
                 "Specify condition_col argument or set analysis@config$condition_col. ",
                 "This function tests genes with CONDITION-SPECIFIC q-dependent patterns only.",
                 call. = FALSE)
        }
    }

    if (length(analysis@diversity_results) == 0) {
        stop("Diversity results required. Run calculate_diversity() first.", call. = FALSE)
    }

    condition_col
}

#' Internal: Resolve SRH parameters from config
#'
#' @noRd
.resolve_srh_params <- function(analysis, multicorr, nperm_mode, paired, subject_col,
    nthreads, wy_randomizations, entropy_col, q_col, gene_col) {
    dots <- list()

    # Match enums early
    if (!missing(multicorr)) {
        multicorr <- match.arg(multicorr, c("hochberg", "benjamini-yekutieli", "westfall-young",
            "none"))
        dots$multicorr <- multicorr
    } else if ("multicorr" %in% names(analysis@config)) {
        dots$multicorr <- analysis@config$multicorr
    }

    if (!missing(nperm_mode)) {
        nperm_mode <- match.arg(nperm_mode, c("standard", "conservative", "interactive"))
        dots$nperm_mode <- nperm_mode
    } else if ("nperm_mode" %in% names(analysis@config)) {
        dots$nperm_mode <- analysis@config$nperm_mode
    }

    # Add column parameters
    dots$entropy_col <- entropy_col
    dots$q_col <- q_col
    dots$gene_col <- gene_col

    # Use resolve_slot_param for remaining parameters Note: q is NOT resolved
    # here - always auto-detected from diversity_results in
    # .prepare_multi_q_se()
    paired <- resolve_slot_param(paired, analysis@config, "paired", NULL)
    subject_col <- resolve_slot_param(subject_col, analysis@config, "subject_col",
        NULL)
    nthreads <- resolve_slot_param(nthreads, analysis@config, "nthreads", 1)

    if (!is.null(paired))
        dots$paired <- paired
    if (!is.null(subject_col))
        dots$subject_col <- subject_col
    dots$nthreads <- nthreads
    dots$wy_randomizations <- wy_randomizations

    list(dots = dots)
}

#' Internal: Prepare multi-Q SummarizedExperiment for testing
#'
#' @noRd
.prepare_multi_q_se <- function(analysis) {
    # Check cache first
    if (!is.null(analysis@metadata$diversity_combined) && is.list(analysis@metadata$diversity_combined) &&
        !is.null(analysis@metadata$diversity_combined$combined_se)) {

        se_multi_q <- analysis@metadata$diversity_combined$combined_se
        if (is(se_multi_q, "SummarizedExperiment") && ncol(se_multi_q) > 0) {
            return(se_multi_q)
        }
    }

    # Fallback: combine per-Q results
    q_keys <- names(analysis@diversity_results)
    combined_assay_list <- list()
    combined_coldata_list <- list()
    common_rownames <- NULL

    for (key in sort(q_keys)) {
        se <- analysis@diversity_results[[key]]
        q_val <- as.numeric(sub("^q_", "", key))

        # Ensure SE format
        if (!is(se, "SummarizedExperiment")) {
            if (is.matrix(se) || is.data.frame(se)) {
                se <- SummarizedExperiment(assays = list(diversity = as.matrix(se)))
            } else {
                stop("Diversity result for ", key, " is not a SummarizedExperiment or matrix",
                  call. = FALSE)
            }
        }

        assay_data <- SummarizedExperiment::assay(se, 1)
        assay_rownames <- rownames(assay_data)

        if (is.null(assay_rownames)) {
            assay_rownames <- paste0("gene_", seq_len(nrow(assay_data)))
        }
        if (is.null(common_rownames)) {
            common_rownames <- assay_rownames
        } else if (!identical(common_rownames, assay_rownames)) {
            if (length(common_rownames) == length(assay_rownames)) {
                assay_data <- assay_data[common_rownames, , drop = FALSE]
            } else {
                stop("Diversity result for ", key, " has different number of genes",
                  call. = FALSE)
            }
        }
        rownames(assay_data) <- common_rownames

        # Rename columns with q-value suffix
        orig_colnames <- colnames(assay_data)
        if (is.null(orig_colnames))
            orig_colnames <- paste0("sample_", seq_len(ncol(assay_data)))
        unique_colnames <- paste0(orig_colnames, "_q=", q_val)
        colnames(assay_data) <- unique_colnames

        # Get and update colData
        cd <- as.data.frame(SummarizedExperiment::colData(se))
        if (nrow(cd) == 0) {
            cd <- data.frame(q = rep(q_val, ncol(assay_data)))
        } else if (!"q" %in% colnames(cd)) {
            cd$q <- q_val
        }
        rownames(cd) <- unique_colnames

        combined_assay_list[[key]] <- assay_data
        combined_coldata_list[[key]] <- cd
    }

    # Combine horizontally
    combined_assay <- do.call(cbind, combined_assay_list)
    combined_coldata_df <- do.call(rbind, combined_coldata_list)
    colnames(combined_assay) <- rownames(combined_coldata_df)

    # Get rowData from first SE
    first_se <- analysis@diversity_results[[sort(q_keys)[1]]]
    if (!is(first_se, "SummarizedExperiment")) {
        first_se <- SummarizedExperiment(assays = list(diversity = as.matrix(first_se)))
    }
    rd <- tryCatch(SummarizedExperiment::rowData(first_se), error = function(e) NULL)

    se_multi_q <- SummarizedExperiment(assays = list(diversity = combined_assay),
        colData = combined_coldata_df)
    if (!is.null(rd) && nrow(rd) > 0) {
        SummarizedExperiment::rowData(se_multi_q) <- rd
    }

    se_multi_q
}

#' Internal: Store SRH results in analysis object
#'
#' @noRd
.store_srh_results <- function(analysis, result, output_file, verbose) {
    # Store in dedicated rank_test_results slot (not in rrm_results)
    if (is.list(analysis@rank_test_results)) {
        analysis@rank_test_results$rank_test <- result
    } else {
        analysis@rank_test_results <- list(rank_test = result)
    }

    if (!is.null(output_file)) {
        result_df <- as.data.frame(result)
        save_analysis_output(result_df, output_file, object = analysis, verbose = verbose,
            func_name = "calculate_srh")
    }

    analysis
}
