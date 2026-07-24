#' Build a Complete TSENATAnalysis Object
#'
#' Convenience wrapper that  combines \code{. build_se()} and 
#' \code{TSENATAnalysis()}
#' into a single function call. This creates a complete analysis object
#' ready for
#' Tsallis entropy computation and downstream analysis.
#'
#' @param readcounts A matrix or data.frame of transcript-level read counts with
#' transcript IDs as row names and sample names as column names. Typically
#' output
#'   from quantification tools (SALMON,  kallisto,  etc. ).  Optional when 
#' \code{salmon_dir}
#'   is provided (in which case readcounts are auto-loaded from quant.sf files).
#'
#' @param salmon_dir Optional character path to directory containing Salmon
#' quantification
#'   output.  Expected structure:  \code{salmon_dir/sample_name/quant. sf}.
#'  When provided,
#'   automatically discovers and  reads all quant. sf files.
#'  If both \code{readcounts} and
#'   \code{salmon_dir} are provided,  \code{salmon_dir} takes precedence.
#'  Default:  \code{NULL}.
#'
#' @param tx2gene Either:
#'   - A path to a GFF3 or GFF3.gz file containing transcript-to-gene mapping
#'   - A path to a TSV file with columns 'Transcript' and 'Gene'
#'   - A data.frame with transcript-to-gene mapping
#'
#' @param assay_name Character. Name for the assay (default: 'counts').
#'
#' @param metadata Optional data.frame with sample metadata. Should have sample
#' names as row names and metadata columns (e.g., sample_type, condition, etc.).
#' If NULL, will attempt to read from \code{config$metadata}.
#' Priority: explicit \code{metadata} argument > \code{config$metadata} > NULL.
#'
#' @param tpm REQUIRED matrix of transcript-level TPM values (Transcripts Per Million).
#' Must be provided and will be stored in the SummarizedExperiment for use during
#'   diverse filtering. Same dimensions as readcounts required (rows = transcripts, 
#'   columns = samples). Typically from SALMON quantification output.
#'   If not provided to build_analysis(), subsequent filter_analysis() calls 
#'   will fail with an explicit error message.
#'
#' @param effective_length REQUIRED numeric vector of transcript effective lengths
#'   (e.g., from SALMON EffectiveLength column). Length must match nrow(readcounts).
#'   Typically obtained as the median effective length across all samples.
#'   If not provided to build_analysis(), the data will not be stored for later use
#'   in length-normalized calculations.
#'
#' @param config Optional list of configuration parameters to store in the
#'   TSENATAnalysis object. Can also contain \code{config$metadata} which will be
#'   used if the \code{metadata} argument is NULL. Following Bioconductor best practices
#'   (fail-fast principle), create configuration via \code{\link{TSENAT_config}()} FIRST,
#'   then pass to \code{build_analysis()} at object construction time. This ensures
#'   invalid parameters are caught immediately, before analysis proceeds.
#'   See examples below for recommended usage pattern.
#'
#' @param skip Logical. If TRUE, allow unmapped transcripts (transcripts not
#' found
#' in tx2gene mapping) and remove them from analysis. If FALSE (default),
#' stop with
#' an error when unmapped transcripts are detected. Useful for handling data
#' with
#'   transcript IDs that don't match the annotation file provided.
#'
#' @param verbose Logical. If TRUE, print informative messages during execution
#'   (e.g., Salmon sample discovery, progress on data loading). Default: TRUE.
#'
#' @details
#' When using \code{salmon_dir}, the function automatically:
#' \enumerate{
#'   \item Discovers all Salmon sample folders and quant.sf files
#'   \item Reads transcript counts (NumReads), TPM, and effective_length
#'   \item Extracts sample names from directory structure
#'   \item Creates count matrix ready for analysis
#' }
#'
#' The \code{salmon_dir} parameter provides a convenient alternative to manually
#' constructing the \code{readcounts} matrix,
#'  especially useful in Galaxy workflows.
#'
#' @return A \code{TSENATAnalysis} S4 object with:
#'   \item{@se}{The SummarizedExperiment containing transcript counts and 
#' metadata}
#'   \item{@config}{Analysis configuration (empty list or user-provided)}
#'   \item{@diversity_results}{Empty list (populated by calculate_diversity())}
#'   \item{@divergence_results}{Empty list (populated by calculate_divergence())}
#'   \item{@sait_results}{Empty list (populated by calculate_sait())}
#'   \item{@jackknife_results}{Empty list (populated by jackknife functions)}
#'   \item{@plots}{Empty list (populated by plotting functions)}
#'   \item{@metadata}{Metadata with package version and creation timestamp}
#'
#'
#' \strong{CRITICAL: TPM and effective_length Requirements}
#'
#' Both \code{tpm} and \code{effective_length} MUST be provided to ensure correct
#' filtering and normalization in downstream analysis:
#' \itemize{
#'   \item \code{filter_analysis()} requires TPM data (stored in metadata).
#'     If TPM is missing, the function will fail with an explicit error message
#'     that guides you to pass it to \code{build_analysis()}.
#'   \item \code{calculate_diversity()} uses \code{effective_length} for 
#'     length-normalized entropy calculations.
#' }
#'
#' Following Bioconductor best practices (fail-fast principle), these are explicit
#' parameters, not optional. They must be passed at object construction time:
#' \preformatted{
#' analysis <- build_analysis(
#'   readcounts = readcounts,
#'   metadata = metadata_df,
#'   tx2gene = gff3_file,
#'   tpm = tpm,                    # REQUIRED from Salmon output
#'   effective_length = effective_length,  # REQUIRED from Salmon output
#'   config = config
#' )
#' }
#'
#' This wrapper combines two steps into one:
#' \enumerate{
#'   \item Call \code{.build_se()} to create a SummarizedExperiment from transcript counts
#'   \item Wrap the result in \code{TSENATAnalysis()} to create the analysis object
#' }
#'
#' The returned object is ready for diversity analysis via \code{calculate_diversity()}.
#'
#' If you need to inspect or filter the SummarizedExperiment before creating the
#' TSENATAnalysis object, call \code{.build_se()} and \code{TSENATAnalysis()} separately.
#'
#' @seealso
#' \code{\link{TSENATAnalysis}} for the S4 class structure
#' \code{\link{calculate_diversity}} for computing Tsallis entropy
#'
#' @examples
#' # Create example transcript count data
#' set.seed(42)
#' n_genes <- 10
#' n_isoforms_per_gene <- 3
#' n_isoforms <- n_genes * n_isoforms_per_gene
#' n_samples <- 10
#'
#' # Generate count matrix
#' counts <- matrix(rpois(n_isoforms * n_samples, lambda = 20),
#'                  nrow = n_isoforms, ncol = n_samples)
#' rownames(counts) <- paste0('TX_', 1:n_isoforms)
#' colnames(counts) <- paste0('Sample_', 1:n_samples)
#'
#' # Create tx2gene mapping
#' tx2gene <- data.frame(
#'   Transcript = rownames(counts),
#'   Gene = rep(paste0('GENE_', 1:n_genes), each = n_isoforms_per_gene))
#'
#' # Create sample metadata
#' metadata <- data.frame(
#'   sample = colnames(counts),
#'   condition = rep(c('control', 'treatment'), each = 5),
#'   row.names = colnames(counts))
#'
#' # Build analysis object - use NAMED parameters to avoid confusion
#' # Method 1: With explicit tx2gene data.frame (most common)
#' config <- TSENAT_config(sample_col = 'sample', condition_col = 'condition')
#' analysis <- build_analysis(
#'   readcounts = counts,
#'   tx2gene = tx2gene,
#'   metadata = metadata,
#'   config = config)
#'
#' # Verify the analysis object was created
#' analysis
#' print(dim(analysis))
#'
#' # Method 2: From Salmon quantification folder
#' # Requires directory structure like:
#' #   salmon_output/
#' #     sample1/quant.sf
#' #     sample2/quant.sf
#' #     ...
#' # 
#' # salmon_dir <- '/path/to/salmon/directory'
#' # 
#' # First create sample metadata matching Salmon sample names
#' # salmon_metadata <- data.frame(
#' #   condition = c('control', 'control', 'treatment', 'treatment'),
#' #   row.names = c('sample1', 'sample2', 'sample3', 'sample4')
#' # )
#' # 
#' # analysis_salmon <- build_analysis(
#' #   salmon_dir = salmon_dir,
#' #   tx2gene = 'annotation.gff3.gz',  # Auto-parsed from GFF3
#' #   metadata = salmon_metadata
#' # )
#' #
#' # Method 3: Hybrid - Salmon counts with manual tx2gene
#' # analysis_hybrid <- build_analysis(
#' #   salmon_dir = salmon_dir,
#' #   tx2gene = tx2gene,  # data.frame instead of file
#' #   metadata = salmon_metadata
#' # )
#' #
#' # Method 4: Pass metadata via config (parameter resolution pattern)
#' # cfg <- TSENAT_config()
#' # cfg$metadata <- metadata
#' # analysis_with_config <- build_analysis(
#' #   readcounts = counts,
#' #   tx2gene = tx2gene,
#' #   config = cfg
#' #   # Note: metadata argument omitted - will be read from config$metadata
#' # )
#'
#' # Advanced: Assigning metadata to assays after object creation
#' # When adding metadata to SummarizedExperiment assays, always use the
#' # S4Vectors namespace to ensure proper method dispatch:
#' #   
#' #   se <- getSE(analysis)
#' #   assay_with_ci <- SummarizedExperiment::assay(se, 'log2fc_ci')
#' #   S4Vectors::metadata(assay_with_ci)$lower <- ci_lower_bounds
#' #   S4Vectors::metadata(assay_with_ci)$upper <- ci_upper_bounds
#' #
#' # Note: Avoid using metadata(assay) without the namespace - this can
#' # cause silent failures in S4 object metadata assignment.
#'
#' @export
build_analysis <- function(readcounts = NULL, salmon_dir = NULL, tx2gene, assay_name = "counts",
    metadata = NULL, tpm = NULL, effective_length = NULL, config = list(), skip = FALSE,
    verbose = FALSE) {
    # Parameter resolution: explicit argument takes priority, then config
    if (is.null(metadata) && !is.null(config$metadata)) {
        metadata <- config$metadata
        if (verbose)
            message("[build_analysis] Reading metadata from config$metadata")
    }

    # Handle salmon_dir parameter - auto-load Salmon quantification data
    if (!is.null(salmon_dir)) {
        # Detect Salmon samples
        salmon_info <- .detect_salmon_samples(salmon_dir)

        if (verbose)
            message("[build_analysis] Found ", salmon_info$count, " Salmon samples")

        # Validate Salmon sample names match metadata (if metadata provided)
        if (!is.null(metadata)) {
            metadata_samples <- rownames(metadata)
            salmon_samples <- salmon_info$sample_names

            # Check if all Salmon samples have corresponding metadata
            missing_in_metadata <- setdiff(salmon_samples, metadata_samples)
            missing_in_salmon <- setdiff(metadata_samples, salmon_samples)

            if (length(missing_in_metadata) > 0 || length(missing_in_salmon) > 0) {
                error_msg <- "[build_analysis] Sample name mismatch between Salmon folder and metadata:\n"

                if (length(missing_in_metadata) > 0) {
                  error_msg <- paste0(error_msg, "  Salmon samples NOT in metadata (",
                    length(missing_in_metadata), "): ", paste(missing_in_metadata,
                      collapse = ", "), "\n")
                }

                if (length(missing_in_salmon) > 0) {
                  error_msg <- paste0(error_msg, "  Metadata samples NOT in Salmon folder (",
                    length(missing_in_salmon), "): ", paste(missing_in_salmon, collapse = ", "),
                    "\n")
                }

                error_msg <- paste0(error_msg, "\n  Suggestion: Ensure sample folder names in Salmon directory exactly match",
                  "\n  the row names of the metadata data.frame (case-sensitive)")

                stop(error_msg)
            }

            if (verbose)
                message("[build_analysis] [OK] Sample names match metadata")
        }

        # Read Salmon quantification files
        salmon_data <- .read_salmon_samples(file_paths = salmon_info$file_paths,
            sample_names = salmon_info$sample_names, include_tpm = is.null(tpm),
            include_eff_length = is.null(effective_length), verbose = verbose)

        # Assign extracted data to parameters
        readcounts <- salmon_data$counts
        if (is.null(tpm))
            tpm <- salmon_data$tpm
        if (is.null(effective_length))
            effective_length <- salmon_data$effective_length
    }

    # Validate readcounts is now available
    if (is.null(readcounts)) {
        stop("[build_analysis] Either 'readcounts' or 'salmon_dir' must be provided\n",
            "  readcounts: matrix/data.frame of transcript counts\n", "  salmon_dir: path to Salmon quantification output directory")
    }

    # Validate column parameters when metadata is provided
    if (!is.null(metadata)) {
        missing_cols <- c()

        if (is.null(config$sample_col)) {
            missing_cols <- c(missing_cols, "sample_col")
        }
        if (is.null(config$condition_col)) {
            missing_cols <- c(missing_cols, "condition_col")
        }

        if (length(missing_cols) > 0) {
            stop("[build_analysis] Metadata provided but required column parameters missing: ",
                paste(missing_cols, collapse = ", "), "\n", "  These parameters MUST be provided in TSENAT_config():\n",
                "    config <- TSENAT_config(\n", "      sample_col = 'sample',        # Column name with sample identifiers\n",
                "      condition_col = 'condition', # Column name with condition/treatment labels\n",
                "      ...\n", "    )\n", "    analysis <- build_analysis(config = config, metadata = metadata_df, ...)",
                call. = FALSE)
        }
    }

    # Build SummarizedExperiment Extract column names from config
    sample_col_value <- if (!is.null(config$sample_col))
        config$sample_col else "sample"
    condition_col_value <- config$condition_col  # May be NULL for non-paired analysis
    subject_col_value <- config$subject_col  # Optional for paired analysis

    se <- .build_se(readcounts = readcounts, tx2gene = tx2gene, assay_name = assay_name,
        metadata = metadata, sample_col = sample_col_value, condition_col = condition_col_value,
        subject_col = subject_col_value, tpm = tpm, effective_length = effective_length,
        skip = skip, verbose = verbose)

    # Ensure sample_id column exists in colData (required by TSENATAnalysis)
    # OPTIMIZATION: Only add if not already present
    if (!"sample_id" %in% colnames(SummarizedExperiment::colData(se))) {
        SummarizedExperiment::colData(se)$sample_id <- colnames(se)
    }

    # Store metadata in config for later use (e.g., in calculate_sait)
    if (!is.null(metadata)) {
        config$metadata <- metadata
    }

    # Wrap in TSENATAnalysis
    analysis <- TSENATAnalysis(se = se, config = config)

    return(analysis)
}
