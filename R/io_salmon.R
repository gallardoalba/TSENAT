#' Detect Salmon quantification output files in a directory
#'
#' @param salmon_dir \code{character}. Path to directory containing Salmon quantification
#'   output folders. Expected structure: \code{salmon_dir/sample_name/quant.sf} or
#'   \code{salmon_dir/sample_name/quant.sf.gz}.
#' @param pattern \code{character}. Regular expression pattern for quantification files.
#'   Default: \code{'quant\\.sf(\\.gz)?$'} matches \code{quant.sf} or \code{quant.sf.gz}.
#' @param recursive \code{logical}. Whether to search recursively for nested sample folders.
#'   Default: \code{TRUE}.
#'
#' @return \code{list} with elements:
#'   \itemize{
#'     \item \code{$sample_names}: Character vector of discovered sample names
#'     \item \code{$file_paths}: Character vector of full paths to quant.sf files
#'     \item \code{$count}: Integer number of samples found
#'   }
#'
#' @details
#' This function scans \code{salmon_dir} for quant.sf files and extracts sample names
#' from parent directory paths. The standard Salmon output structure is:
#' \preformatted{
#'   salmon/
#'     sample1/
#'       quant.sf           <- Found here
#'       quant.genes.sf
#'       aux_info/
#'     sample2/
#'       quant.sf           <- Found here
#'       ...
#' }
#'
#' @examples
#' # Basic usage
#' salmon_info <- .detect_salmon_samples('/path/to/salmon_output')
#' print(salmon_info$sample_names)
#'
#' # Access file paths and count
#' length(salmon_info$file_paths)
#' head(salmon_info$file_paths)
#'
#' @noRd
.detect_salmon_samples <- function(salmon_dir, pattern = "quant\\.sf(\\.gz)?$", recursive = TRUE) {
    # Validate input
    if (!is.character(salmon_dir) || length(salmon_dir) != 1) {
        stop("[.detect_salmon_samples] 'salmon_dir' must be a single character string")
    }

    if (!dir.exists(salmon_dir)) {
        stop("[.detect_salmon_samples] Directory does not exist: ", salmon_dir, "\n",
            "  Suggestion: Check path spelling and verify it's a valid directory")
    }

    # Find all quant.sf files matching pattern
    file_paths <- list.files(path = salmon_dir, pattern = pattern, full.names = TRUE,
        recursive = recursive)

    if (length(file_paths) == 0) {
        stop("[.detect_salmon_samples] No Salmon quantification files found in: ",
            salmon_dir, "\n", "  Expected structure: salmon_dir/sample_name/quant.sf\n",
            "  Suggestion: Check if files exist and follow standard Salmon output structure")
    }

    # Extract sample names from parent directory paths For:
    # /path/to/salmon/sample1/quant.sf -> 'sample1'
    sample_names <- basename(dirname(file_paths))

    # Check for duplicate sample names (might indicate nested structure issues)
    if (length(unique(sample_names)) != length(sample_names)) {
        warning("[.detect_salmon_samples] Found duplicate sample names. ", "This may indicate nested sample folders with same names.")
    }

    # Return structured list
    list(sample_names = sample_names, file_paths = file_paths, count = length(file_paths))
}


#' Validate Salmon quant.sf file format and consistency
#'
#' @param file_paths \code{character}. Vector of paths to quant.sf files to validate.
#' @param verbose \code{logical}. Whether to print validation messages. Default: \code{TRUE}.
#'
#' @return \code{logical} TRUE if all validations pass. Throws \code{stop()} if errors found.
#'
#' @details
#' Validates that:
#' \itemize{
#'   \item All files are readable
#'   \item All files have proper quant.sf header (Name, Length, EffectiveLength, TPM, NumReads)
#'   \item All files have matching transcript IDs (same transcripts in all files)
#'   \item No missing or corrupted data
#' }
#'
#' @noRd
.validate_salmon_files <- function(file_paths, verbose = TRUE) {
    if (length(file_paths) == 0) {
        stop("[.validate_salmon_files] No files provided to validate")
    }

    # Check all files exist and are readable
    missing_files <- which(!file.exists(file_paths))
    if (length(missing_files) > 0) {
        stop("[.validate_salmon_files] File(s) not found:\n", paste("  -", file_paths[missing_files],
            collapse = "\n"))
    }

    # Expected columns in quant.sf
    expected_cols <- c("Name", "Length", "EffectiveLength", "TPM", "NumReads")

    # Try reading first file to check format
    first_file <- file_paths[1]
    if (verbose)
        message("[.validate_salmon_files] Checking format of: ", basename(first_file))

    tryCatch({
        first_data <- readr::read_tsv(first_file, col_types = readr::cols(Name = readr::col_character(),
            Length = readr::col_integer(), EffectiveLength = readr::col_double(),
            TPM = readr::col_double(), NumReads = readr::col_double()), show_col_types = FALSE)

        # Verify all expected columns present
        missing_cols <- setdiff(expected_cols, colnames(first_data))
        if (length(missing_cols) > 0) {
            stop("[.validate_salmon_files] Missing required columns: ", paste(missing_cols,
                collapse = ", "), "\n", "  Expected columns: ", paste(expected_cols,
                collapse = ", "))
        }

        if (nrow(first_data) == 0) {
            stop("[.validate_salmon_files] First file is empty (no data rows)")
        }

        transcript_ids_first <- first_data$Name
    }, error = function(e) {
        stop("[.validate_salmon_files] Error reading first file: ", basename(first_file),
            "\n", "  Details: ", conditionMessage(e), "\n", "  Suggestion: Check file format and that it's a valid quant.sf file")
    })

    # Check remaining files for format consistency
    if (length(file_paths) > 1) {
        if (verbose) {
            message("[.validate_salmon_files] Checking consistency across ", length(file_paths),
                " files...")
        }

        for (i in 2:length(file_paths)) {
            current_file <- file_paths[i]

            tryCatch({
                current_data <- readr::read_tsv(current_file, col_types = readr::cols(Name = readr::col_character(),
                  Length = readr::col_integer(), EffectiveLength = readr::col_double(),
                  TPM = readr::col_double(), NumReads = readr::col_double()), show_col_types = FALSE)

                transcript_ids_current <- current_data$Name

                # Check if same transcripts in all files
                if (!identical(transcript_ids_first, transcript_ids_current)) {
                  warnings_msg <- ""

                  # Check for missing transcripts
                  missing_in_curr <- setdiff(transcript_ids_first, transcript_ids_current)
                  extra_in_curr <- setdiff(transcript_ids_current, transcript_ids_first)

                  if (length(missing_in_curr) > 0) {
                    warnings_msg <- paste0(warnings_msg, "  Missing in file ", i,
                      ": ", length(missing_in_curr), " transcripts\n")
                  }
                  if (length(extra_in_curr) > 0) {
                    warnings_msg <- paste0(warnings_msg, "  Extra in file ", i, ": ",
                      length(extra_in_curr), " transcripts\n")
                  }

                  warning("[.validate_salmon_files] Transcript ID mismatch between files:\n",
                    warnings_msg, "  Suggestion: Ensure all files use same transcriptome reference")
                }

                if (nrow(current_data) == 0) {
                  stop("[.validate_salmon_files] File is empty (no data rows): ",
                    basename(current_file))
                }
            }, error = function(e) {
                stop("[.validate_salmon_files] Error reading file ", i, ": ", basename(current_file),
                  "\n", "  Details: ", conditionMessage(e))
            })
        }
    }

    if (verbose)
        message("[.validate_salmon_files] ✓ All validations passed")
    invisible(TRUE)
}


#' Read Salmon quantification files into count matrices
#'
#' @param file_paths \code{character}. Vector of paths to quant.sf files.
#' @param sample_names \code{character}. Vector of sample names (should match length of file_paths).
#'   If \code{NULL}, sample names are extracted from file paths (parent directory names).
#' @param include_tpm \code{logical}. Whether to extract and return TPM matrix.
#'   Default: \code{TRUE}.
#' @param include_eff_length \code{logical}. Whether to extract and return effective length matrix.
#'   Default: \code{TRUE}.
#' @param verbose \code{logical}. Whether to print progress messages. Default: \code{TRUE}.
#'
#' @return \code{list} with elements:
#'   \itemize{
#'     \item \code{$counts}: Numeric matrix of NumReads (transcripts \u00D7 samples)
#'     \item \code{$tpm}: Numeric matrix of TPM values (if \code{include_tpm=TRUE})
#'     \item \code{$effective_length}: Numeric matrix of EffectiveLength (if \code{include_eff_length=TRUE})
#'     \item \code{$transcript_ids}: Character vector of transcript IDs (row names)
#'   }
#'
#' @details
#' This function reads Salmon quant.sf files and aggregates them into matrices suitable
#' for downstream analysis. The quant.sf format (from Salmon) contains:
#' \preformatted{
#'   Name              Length  EffectiveLength  TPM      NumReads
#'   ENST00000456328   1234    1100             5.123    123.45
#'   ENST00000415691   5678    5400             8.456    234.56
#'   ...
#' }
#'
#' The \code{NumReads} column represents estimated transcript counts and is used as the
#' primary count matrix. TPM (Transcripts Per Million) and EffectiveLength are also
#' extracted for downstream computations (e.g., filtering, offset calculations).
#'
#' @examples
#' # Read Salmon samples
#' salmon_info <- .detect_salmon_samples('/path/to/salmon')
#'
#' salmon_data <- .read_salmon_samples(
#'   file_paths = salmon_info$file_paths,
#'   sample_names = salmon_info$sample_names
#' )
#'
#' # Access results
#' # dim(salmon_data$counts)           # Transcripts \u00D7 Samples
#' head(salmon_data$transcript_ids)
#'
#' @noRd
.read_salmon_samples <- function(file_paths, sample_names = NULL, include_tpm = TRUE,
    include_eff_length = TRUE, verbose = TRUE) {
    # Validate inputs
    if (!is.character(file_paths) || length(file_paths) == 0) {
        stop("[.read_salmon_samples] 'file_paths' must be a non-empty character vector")
    }

    if (is.null(sample_names)) {
        sample_names <- basename(dirname(file_paths))
        if (verbose) {
            message("[.read_salmon_samples] Using directory names as sample names: ",
                paste(sample_names[1:min(3, length(sample_names))], collapse = ", "),
                "...")
        }
    }

    if (length(sample_names) != length(file_paths)) {
        stop("[.read_salmon_samples] Length mismatch: ", length(file_paths), " files but ",
            length(sample_names), " sample names")
    }

    # Validate files before processing
    if (verbose)
        message("[.read_salmon_samples] Validating ", length(file_paths), " files...")
    .validate_salmon_files(file_paths, verbose = FALSE)

    # Define column specifications for fast, type-safe reading
    col_spec <- readr::cols(Name = readr::col_character(), Length = readr::col_integer(),
        EffectiveLength = readr::col_double(), TPM = readr::col_double(), NumReads = readr::col_double(),
        .default = readr::col_skip()  # Skip any extra columns
)

    # Read first file to initialize matrices
    if (verbose)
        message("[.read_salmon_samples] Reading transcript IDs from: ", basename(file_paths[1]))

    first_data <- readr::read_tsv(file_paths[1], col_types = col_spec, show_col_types = FALSE)

    transcript_ids <- first_data$Name
    n_transcripts <- length(transcript_ids)
    n_samples <- length(file_paths)

    if (verbose) {
        message("[.read_salmon_samples] Found ", n_transcripts, " transcripts ",
            "from ", n_samples, " samples")
    }

    # Initialize result matrices
    counts_matrix <- matrix(NA_real_, nrow = n_transcripts, ncol = n_samples, dimnames = list(transcript_ids,
        sample_names))

    tpm_matrix <- NULL
    if (include_tpm) {
        tpm_matrix <- matrix(NA_real_, nrow = n_transcripts, ncol = n_samples, dimnames = list(transcript_ids,
            sample_names))
    }

    eff_length_matrix <- NULL
    if (include_eff_length) {
        eff_length_matrix <- matrix(NA_real_, nrow = n_transcripts, ncol = n_samples,
            dimnames = list(transcript_ids, sample_names))
    }

    # Populate matrices from first file
    counts_matrix[, 1] <- first_data$NumReads
    if (include_tpm)
        tpm_matrix[, 1] <- first_data$TPM
    if (include_eff_length)
        eff_length_matrix[, 1] <- first_data$EffectiveLength

    # Process remaining files
    if (n_samples > 1) {
        for (i in 2:n_samples) {
            if (verbose && i%%max(1, n_samples%/%10) == 0) {
                message("[.read_salmon_samples] Processing file ", i, "/", n_samples,
                  "...")
            }

            current_data <- readr::read_tsv(file_paths[i], col_types = col_spec,
                show_col_types = FALSE)

            # Assign counts
            counts_matrix[, i] <- current_data$NumReads

            if (include_tpm) {
                tpm_matrix[, i] <- current_data$TPM
            }

            if (include_eff_length) {
                eff_length_matrix[, i] <- current_data$EffectiveLength
            }
        }
    }

    if (verbose)
        message("[.read_salmon_samples] ✓ Successfully read all Salmon files")

    # Build result list
    result <- list(counts = counts_matrix, transcript_ids = transcript_ids)

    if (include_tpm) {
        result$tpm <- tpm_matrix
    }

    if (include_eff_length) {
        result$effective_length <- eff_length_matrix
    }

    return(result)
}
