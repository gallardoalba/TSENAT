## Filter SummarizedExperiment by low-expression transcripts
#' Filter transcripts in a `SummarizedExperiment` by minimum TPM/sample
#'
#' Subset a `SummarizedExperiment` to keep transcripts with more than
#' `min_tpm` TPM (transcripts per million) threshold in strictly more than
#' `min_samples` samples. The function updates assays, `rowData`, and
#' relevant entries in `metadata()` (for example `readcounts` and `tx2gene`)
#' so downstream helpers receive a consistent object.
#'
#' **Filtering Strategy**: Uses TPM-based filtering for results comparable across
#' studies (data from SALMON quantification is already TPM-normalized; Soneson et al. 2015, Law et al. 2014).
#'
#' **TPM Data Source**: Function automatically locates TPM data via:
#' 1. `tpm_assay_name` parameter (specify assay name containing TPM values)
#' 2. Metadata: looks for `salmon_tpm` or `tpm` in SummarizedExperiment metadata
#' 3. Assay names: searches for assay named "tpm" or "abundance"
#' If TPM data NOT found, filtering parameters are compared against raw counts
#' but this is NOT recommended (will produce incorrect results).
#'
#' @param se A `SummarizedExperiment` object containing transcript-level
#'   assays (rows = transcripts, columns = samples). Can contain:
#'   - SALMON TPM values (recommended, already normalized)
#'   - SALMON raw counts (NumReads) if TPM not available
#' @param min_tpm Numeric TPM threshold (default 1.0).
#'   Keeps transcripts with TPM >= `min_tpm` in >= `min_samples` samples.
#'   TPM >= 1 is recommended for typical sequencing depth (Soneson et al. 2015, Law et al. 2014).
#'   Ignored if `stringency` is specified; when stringency is provided, `min_tpm` is
#'   auto-estimated from data using quantile-based approach (Law et al. limma-voom methodology).
#' @param tpm_assay_name Character; name of assay containing TPM data (default: NULL).
#'   If NULL, function searches for TPM data in: metadata$salmon_tpm -> metadata$tpm -> assay named "tpm" -> metadata lookup.
#'   Set explicitly (e.g., `tpm_assay_name = "tpm"`) to use a specific assay by name.
#'   Example from preprocessing: tpm_assay_name = "abundance" for tximport objects.
#' @param min_samples Integer minimum number of samples exceeding
#'   the threshold required to keep a transcript (default 5L). Ignored if
#'   `stringency` is specified.
#' @param stringency Character; auto-calculate `min_samples` and `min_tpm` based on
#'   paired design stringency. Options:
#'   - "soft" (permissive): 25% of samples (min 2), min_tpm = Q1 (25th %ile of mean TPM)
#'   - "medium" (balanced): 50% of samples (min 3), min_tpm = Q2 (median of mean TPM)
#'   - "severe" (stringent): 75% of samples, min_tpm = Q3 (75th %ile of mean TPM)
#'   - NULL (default): use explicit `min_samples` and `min_tpm`
#'   When stringency is specified, both filtering parameters are auto-estimated
#'   from the data distribution (Law et al. limma-voom methodology).
#'   Requires `pair_col` in colData.
#' @param pair_col Character; column name in colData containing pair IDs.
#'   Required when `stringency` is specified. If NULL and `stringency` is set,
#'   function attempts to auto-detect from colData.
#' @param min_tx_per_gene Integer minimum number of transcripts per gene
#'   required to keep a gene (default 2L). Genes with fewer transcripts are
#'   removed since entropy is always 0 for single-transcript genes.
#'   Ignored if `stringency` is specified; when stringency is provided, this is
#'   auto-adjusted (soft=2, medium=2, severe=3) to ensure isoform diversity matches
#'   filtering stringency.
#' @param assay_name Name or index of the assay to use for filtering
#'   (DEPRECATED: use `tpm_assay_name` instead, default: 'counts'). 
#'   **WARNING**: If set to 'counts', filtering uses raw counts with TPM thresholds,
#'   which is incorrect. Use `tpm_assay_name` to specify TPM assay.
#' @param verbose Logical; print before/after counts and filtering parameters when TRUE.
#' @return A filtered `SummarizedExperiment`.
#' @export
#' @examples
#' mat <- matrix(c(0, 6, 7, 2, 8, 9), nrow = 3, dimnames = list(paste0('tx', 1:3), paste0('S', 1:2)))
#' se <- SummarizedExperiment::SummarizedExperiment(assays = list(counts = mat))
#' filt <- filter_se(se, min_samples = 1)
#' class(filt)
filter_se <- function(se, min_samples = 5L, stringency = NULL,
    pair_col = NULL, min_tpm = 1.0, tpm_assay_name = NULL,
    min_tx_per_gene = 2L, assay_name = "counts", verbose = TRUE) {
    if (!is(se, "SummarizedExperiment")) {
        stop("'se' must be a SummarizedExperiment", call. = FALSE)
    }
    
    # ========================================================================
    # LOCATE TPM DATA FOR FILTERING
    # ========================================================================
    # Priority: explicit parameter > metadata > assay names
    assays_list <- SummarizedExperiment::assays(se)
    tpm_assay_mat <- NULL
    tpm_source <- NULL
    
    # Strategy 1: Use explicitly specified TPM assay
    if (!is.null(tpm_assay_name)) {
        if (tpm_assay_name %in% names(assays_list)) {
            tpm_assay_mat <- as.matrix(assays_list[[tpm_assay_name]])
            tpm_source <- sprintf("assay '%s' (user-specified)", tpm_assay_name)
        } else {
            warning(sprintf("TPM assay '%s' not found. Available assays: %s",
                          tpm_assay_name, paste(names(assays_list), collapse=", ")),
                   call. = FALSE)
        }
    }
    
    # Strategy 2: Look in metadata for salmon_tpm or tpm
    if (is.null(tpm_assay_mat)) {
        md <- S4Vectors::metadata(se)
        if (!is.null(md$salmon_tpm) && is.matrix(md$salmon_tpm)) {
            tpm_assay_mat <- as.matrix(md$salmon_tpm)
            tpm_source <- "metadata$salmon_tpm (SALMON preprocessed)"
        } else if (!is.null(md$tpm) && is.matrix(md$tpm)) {
            tpm_assay_mat <- as.matrix(md$tpm)
            tpm_source <- "metadata$tpm"
        }
    }
    
    # Strategy 3: Search for assay named "tpm" or "abundance"
    if (is.null(tpm_assay_mat)) {
        if ("tpm" %in% names(assays_list)) {
            tpm_assay_mat <- as.matrix(assays_list[["tpm"]])
            tpm_source <- "assay 'tpm' (auto-detected)"
        } else if ("abundance" %in% names(assays_list)) {
            tpm_assay_mat <- as.matrix(assays_list[["abundance"]])
            tpm_source <- "assay 'abundance' (tximport format)"
        }
    }
    
    # If no TPM found, warn and fall back to specified/default assay (likely incorrect)
    if (is.null(tpm_assay_mat)) {
        warning("No TPM data found in assays or metadata. Falling back to assay '", assay_name, "'.",
                "\nThis may produce INCORRECT results if '", assay_name, "' contains raw counts.",
                "\nEnsure TPM data is added as an assay or in metadata with salmon_tpm/tpm.",
                call. = FALSE)
        tpm_source <- sprintf("assay '%s' (fallback - NOT TPM!)", assay_name)
        tpm_assay_mat <- NULL  # Will use assay_name below
    }
    
    # Handle stringency-based filtering for paired designs
    # Automatically estimates BOTH min_samples AND min_tpm from data
    if (!is.null(stringency)) {
        if (!(stringency %in% c("soft", "medium", "severe"))) {
            stop("'stringency' must be one of: 'soft', 'medium', 'severe', or NULL", call. = FALSE)
        }
        
        # Detect pair column if not provided
        if (is.null(pair_col)) {
            col_data <- SummarizedExperiment::colData(se)
            # Look for common pair ID column names
            pair_candidates <- c("pair", "pair_id", "paired_samples", "subject", "subject_id", "individual")
            
            # Try colData first
            pair_col <- pair_candidates[pair_candidates %in% colnames(col_data)][1]
            
            # If colData is empty, try metadata
            if (is.na(pair_col)) {
                meta <- S4Vectors::metadata(se)
                if (!is.null(meta$coldata) && is.data.frame(meta$coldata)) {
                    pair_col <- pair_candidates[pair_candidates %in% colnames(meta$coldata)][1]
                    if (!is.na(pair_col)) {
                        # Use metadata coldata for sampling calculation
                        col_data <- meta$coldata
                    }
                }
            }
            
            if (is.na(pair_col)) {
                stop("Could not auto-detect pair column in colData or metadata. Available columns: ",
                     paste(colnames(col_data), collapse = ", "),
                     ". Please specify 'pair_col' parameter.", call. = FALSE)
            }
            if (verbose) {
                message(sprintf("Auto-detected pair column: '%s'", pair_col))
            }
        } else {
            col_data <- SummarizedExperiment::colData(se)
        }
        
        # Verify pair column exists
        if (!(pair_col %in% colnames(col_data))) {
            stop(sprintf("Pair column '%s' not found. Available columns: %s", 
                         pair_col, paste(colnames(col_data), collapse = ", ")), 
                 call. = FALSE)
        }
        
        # Calculate min_samples and min_tpm based on stringency
        n_samples <- ncol(se)
        n_pairs <- length(unique(col_data[[pair_col]]))
        
        if (stringency == "soft") {
            # Low stringency: 25% of samples (permissive but more stringent than previous default)
            min_samples <- max(2L, ceiling(0.25 * n_samples))
            min_tx_per_gene <- 2L  # Allow 2+ isoforms
        } else if (stringency == "medium") {
            # Medium stringency: 50% of samples
            min_samples <- max(3L, ceiling(0.5 * n_samples))
            min_tx_per_gene <- 2L  # Standard: 2+ isoforms
        } else if (stringency == "severe") {
            # High stringency: 75% of samples
            min_samples <- ceiling(0.75 * n_samples)
            min_tx_per_gene <- 3L  # Strict: require 3+ isoforms for high confidence
        }
        
        # Auto-estimate min_tpm based on stringency level using quantiles
        # Calculate mean TPM per gene for quantile estimation
        # USE TPM data if available, otherwise use count data
        if (!is.null(tpm_assay_mat)) {
            assay_mat_temp <- tpm_assay_mat
        } else {
            assays_list_temp <- SummarizedExperiment::assays(se)
            if (is.character(assay_name) && assay_name %in% names(assays_list_temp)) {
                assay_mat_temp <- as.matrix(assays_list_temp[[assay_name]])
            } else if (is.numeric(assay_name) && assay_name >= 1 && assay_name <= length(assays_list_temp)) {
                assay_mat_temp <- as.matrix(assays_list_temp[[assay_name]])
            } else {
                assay_mat_temp <- as.matrix(assays_list_temp[[1]])
            }
        }
        
        mean_tpm <- rowMeans(assay_mat_temp)
        mean_tpm_nonzero <- mean_tpm[mean_tpm > 0]
        
        # Select quantile based on stringency level
        if (stringency == "soft") {
            # Permissive: Q1 (25th percentile) - removes extreme bottom 25%
            quantile_prob <- 0.25
            quant_label <- "Q1"
        } else if (stringency == "medium") {
            # Balanced: Q2/Median (50th percentile) - removes bottom 50%
            quantile_prob <- 0.50
            quant_label <- "Q2/Median"
        } else if (stringency == "severe") {
            # Stringent: Q3 (75th percentile) - removes bottom 75%
            quantile_prob <- 0.75
            quant_label <- "Q3"
        }
        
        min_tpm_estimated <- as.numeric(quantile(mean_tpm_nonzero, probs = quantile_prob, na.rm = TRUE))
        # Bound between 0.1 and 5.0 for practical constraints
        min_tpm_estimated <- max(0.1, min(min_tpm_estimated, 5.0))
        min_tpm <- min_tpm_estimated
        
        if (verbose) {
            message(sprintf("Stringency: '%s' (n_pairs=%d, n_samples=%d) -> min_samples=%d, min_tx_per_gene=%d",
                            stringency, n_pairs, n_samples, min_samples, min_tx_per_gene))
            message(sprintf("Auto-estimated min_tpm from data: %.3f (%s of mean TPM distribution)", 
                           min_tpm, quant_label))
        }
    }

    assays_list <- SummarizedExperiment::assays(se)
    # Use TPM data if found, otherwise use specified assay
    if (!is.null(tpm_assay_mat)) {
        assay_mat <- tpm_assay_mat
    } else {
        if (is.character(assay_name) && assay_name %in% names(assays_list)) {
            assay_mat <- as.matrix(assays_list[[assay_name]])
        } else if (is.numeric(assay_name) && assay_name >= 1 && assay_name <= length(assays_list)) {
            assay_mat <- as.matrix(assays_list[[assay_name]])
        } else {
            # fallback to first assay
            assay_mat <- as.matrix(assays_list[[1]])
            warning("Requested assay not found; using first assay.", call. = FALSE)
        }
    }

    if (!is.numeric(assay_mat)) {
        stop("Assay data must be numeric.", call. = FALSE)
    }

    before <- nrow(assay_mat)
    
    # Cap min_samples to the number of available samples to avoid impossible
    # thresholds on small datasets and use an inclusive threshold so callers
    # requesting e.g. 'min_samples = ncol' keep rows present in all samples.
    # However, if min_samples is explicitly greater than ncol, warn and don't cap
    # to allow filtering to 0 rows when the condition is impossible
    if (min_samples > ncol(assay_mat)) {
        warning(sprintf("min_samples (%d) is greater than number of samples (%d). This will likely result in 0 rows being kept.", 
                       min_samples, ncol(assay_mat)), call. = FALSE)
        eff_min_samples <- as.integer(min_samples)  # Don't cap, let it filter to 0
    } else {
        eff_min_samples <- min(as.integer(min_samples), ncol(assay_mat))
    }
    
    # TPM-based filtering (data is already TPM-normalized from SALMON)
    tokeep <- rowSums(assay_mat >= min_tpm) >= eff_min_samples
    
    if (verbose) {
        msg <- sprintf("Filtering data source: %s", tpm_source)
        message(msg)
        msg <- sprintf("TPM-based filtering: min_tpm = %.3f", min_tpm)
        if (!is.null(stringency)) {
            msg <- paste0(msg, " (data-driven, auto-estimated from ", stringency, " stringency)")
        } else {
            msg <- paste0(msg, " (manual specification)")
        }
        message(msg)
    }
    
    # Filter out genes with fewer than min_tx_per_gene transcripts (entropy = 0 for single transcripts)
    if (min_tx_per_gene > 1) {
        genes_vec <- .get_gene_ids(se)
        if (!is.null(genes_vec)) {
            # Count transcripts per gene among those that pass count filtering
            tx_per_gene <- table(genes_vec[tokeep])
            genes_to_keep <- names(tx_per_gene)[tx_per_gene >= min_tx_per_gene]
            # Only keep transcripts from genes with sufficient transcripts
            tokeep <- tokeep & (genes_vec %in% genes_to_keep)
        }
    }
    
    after <- sum(tokeep)

    if (verbose) {
        message(sprintf("Transcripts: before = %d, after = %d", before, after))
    }

    if (after == 0L) {
        warning("Filtering removed all transcripts; returning empty SummarizedExperiment.",
            call. = FALSE)
    }

    # subset all assays
    new_assays <- S4Vectors::SimpleList(lapply(assays_list, function(a) {
        if (is.matrix(a) || is.data.frame(a)) {
            as.matrix(a)[tokeep, , drop = FALSE]
        } else {
            a
        }
    }))
    names(new_assays) <- names(assays_list)

    # subset rowData if present
    rd <- NULL
    if (nrow(SummarizedExperiment::rowData(se)) > 0) {
        rd <- SummarizedExperiment::rowData(se)[tokeep, , drop = FALSE]
    }

    # subset metadata readcounts and tx2gene if present
    md <- S4Vectors::metadata(se)
    if (!is.null(md$readcounts) && is.matrix(md$readcounts)) {
        md$readcounts <- as.matrix(md$readcounts)[tokeep, , drop = FALSE]
    }
    if (!is.null(md$tx2gene) && is.data.frame(md$tx2gene)) {
        txmap <- md$tx2gene
        txcol <- colnames(txmap)[1]
        md$tx2gene <- txmap[txmap[[txcol]] %in% rownames(assay_mat)[tokeep], , drop = FALSE]
    }

    # construct new SE with same colData
    new_se <- SummarizedExperiment::SummarizedExperiment(assays = new_assays, rowData = if (!is.null(rd)) {
        rd
    } else {
        S4Vectors::DataFrame()
    }, colData = SummarizedExperiment::colData(se), metadata = c(S4Vectors::metadata(se),
        list(filtered = list(min_samples = min_samples, 
                            min_tpm = min_tpm, 
                            min_tx_per_gene = min_tx_per_gene,
                            stringency = stringency))))

    # attach updated metadata pieces
    S4Vectors::metadata(new_se)$readcounts <- if (!is.null(md$readcounts)) {
        md$readcounts
    } else {
        NULL
    }
    S4Vectors::metadata(new_se)$tx2gene <- if (!is.null(md$tx2gene)) {
        md$tx2gene
    } else {
        NULL
    }

    return(new_se)
}
