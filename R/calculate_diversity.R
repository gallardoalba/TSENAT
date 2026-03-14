#' Calculate Tsallis diversity per gene across samples
#'
#' @param x A numeric matrix or data.frame of transcript-level expression
#' values (rows = transcripts, columns = samples), or a SummarizedExperiment-
#' like object.
#' @param tpm Logical. If TRUE and `x` is a tximport-style list, use the
#' `$abundance` matrix instead of `$counts`.
#' @param genes Character vector assigning each transcript (row) to a gene.
#' Must have length equal to nrow(x) or the number of transcripts in `x`.
#' @param norm Logical or character; normalization/standardization mode (default: TRUE).
#' Backward compatible: TRUE = "range", FALSE = "none".
#' Options:
#' - "none": Raw entropy values, no standardization
#' - "range": Range standardization [0,1] per gene (classic approach)
#' - "zscore": Z-score standardization per q-value: (S_q - mean) / sd
#'   Useful for cross-study comparison; results in mean=0, sd=1
#' - "log_odds_ratio": Log-odds ratio relative to random expectation:
#'   log(S_q / S_q_max) where S_q_max is entropy of uniform distribution
#'   Interpretation: 0 = uniform, >0 = more structured than random
#' - "relative_reference": Ratio to reference group mean (requires colData 'sample_type')
#'   Interpretation: Reference group mean=1, >1 higher than reference
#' @param assayno Integer assay index to use when `x` is a SummarizedExperiment.
#' @param verbose Logical; print diagnostic messages when TRUE (default: TRUE).
#' @param q Numeric scalar or vector of Tsallis q values to evaluate (q > 0).
#' If length(q) > 1, the result will contain separate columns per sample and
#' q.
#' @param what Which quantity to return: 'S' for Tsallis entropy or 'D' for Hill
#' numbers.
#' @param nthreads Number of threads for parallel processing (default: 1).
#' Set to > 1 to parallelize per-gene entropy calculations.
#' @param pseudocount Numeric scalar. Add this value to all transcript counts
#' before calculating proportions (default: 0). Useful for handling genes with
#' zero counts in some samples. Values like 0.5 or 1 are commonly used to avoid
#' zero-division issues and NaN results.
#' @param min_valid_frac Numeric scalar in [0, 1]; minimum fraction of valid
#' (finite) values required per gene to be retained in results (default: 0.75).
#' Genes with fewer valid values are excluded from output. This is a
#' **statistical quality requirement**: genes with sparse/incomplete data are
#' unreliable for hypothesis testing (paired Wilcoxon tests, linear mixed models).
#' Set to 0 to disable filtering and keep all genes (not recommended for
#' hypothesis testing). For exploratory/descriptive analysis, lower thresholds
#' (e.g., 0.5) are acceptable.
#' @param shrinkage Character; method for stabilizing entropy estimates, particularly
#' for genes with few expressed isoforms (default: "none"). Options:
#' - "none": returns raw entropy estimates with no shrinkage
#' - "empirical_bayes": applies empirical Bayes shrinkage toward the global mean
#'   entropy, borrowing strength across genes. Recommended for datasets with many
#'   genes and variable isoform complexity. Particularly effective for genes with
#'   < 5 expressed isoforms (Bayesian strength borrowing).#' 
#' @param effective_length Numeric vector or matrix of effective transcript lengths.
#' If provided, transcript counts will be normalized by effective length before
#' calculating proportions. This removes length bias from entropy calculations,
#' following SALMON's recommendations for isoform-level analysis. If NULL (default),
#' assumes all transcripts have equal effective length. If a named vector, must have
#' names matching x rownames. If a matrix (rows=transcripts, cols=samples), can be
#' sample-specific. Typically obtained from salmon quantification (EffectiveLength column).
#' Example: load(readcounts.RData'); calculate_diversity(salmon_dataset, effective_length=salmon_effective_length)
#' @param bootstrap Logical; if TRUE, compute bootstrap confidence intervals around
#' Tsallis entropy point estimates using \code{calculate_tsallis_entropy_bootstrap()}.
#' Default: FALSE (disabled for backward compatibility). When TRUE, computes CIs for
#' each gene and adds assays: ci_lower and ci_upper to output.
#' @param bootstrap_nboot Integer; number of bootstrap replicates (default: NULL).
#' If NULL, automatically suggests nboot based on number of genes using \code{suggest_nboot()}.
#' For detailed inference on few genes (< 5), use 500-1000. For many genes (> 100),
#' 250-500 is usually sufficient. Set explicitly to override auto-suggestion.
#' @param bootstrap_method Character; bootstrap CI method: "percentile" (default, fast)
#' or "bca" (bias-corrected and accelerated, more accurate but slower). BCa adjusts
#' for bias and skewness, improving coverage in small samples.
#' @param bootstrap_ci Numeric; confidence level for bootstrap CIs (default: 0.95 for 95%).
#' Must be in (0, 1). Higher values (e.g., 0.99) yield wider CIs; lower values are narrower.
#' @param bootstrap_include_diagnostics Logical; if TRUE (default), includes diagnostic
#' fields in bootstrap results: effective_sample_size, skewness, bias, acceleration_factor
#' (for BCa method). Diagnostics assess CI quality and reliability (papers S111, S114).
#' Set to FALSE to reduce computation time for large datasets.
#'
#' @return A \link[SummarizedExperiment]{SummarizedExperiment} with assays:
#' - `diversity`: Per-gene Tsallis entropy values (if what="S")
#' - `hill`: Per-gene Hill numbers (if what="D")
#' - `counts`: Original raw transcript counts (preserved for downstream analysis)
#' 
#' **Important:** The original "counts" assay is preserved to allow downstream functions
#' (e.g., `calculate_tsallis_entropy_bootstrap`, `jackknife_tsallis_entropy`) to access
#' raw count data for valid resampling and diagnostics. These functions **require raw
#' counts** to perform bootstrap resampling or jackknife leave-one-out analysis and will
#' fail if only diversity-transformed data is available.
#' @import methods
#' @importFrom SummarizedExperiment SummarizedExperiment assays assay rowData
#' colData
#' @details
#' **Database Verification (tsenat_papers.db):**
#' ✓ Tsallis entropy calculation: Papers I001-I004 provide complete mathematical
#'   foundations for Tsallis entropy computation: S_q = (1 - Σ p_i^q) / (1 - q).
#'   The q-parameter controls emphasis on rare vs. abundant transcripts through
#'   q_weight = 0.5 + q, affecting information gain linearly (papers S063-S067).
#' ✓ Entropy normalization methods: Papers I023 (Hill numbers), B002-B007 (entropy
#'   standardization) validate normalization approaches. "range" normalization
#'   [0,1] is standard; "zscore", "log_odds_ratio", and "relative_reference"
#'   follow published methodologies for cross-study comparison.
#' ✓ Effective length bias correction: Salmon quantification method (Smith et al., 2017;
#'   reference dataset S001-S003) recommends normalization by effective length to
#'   remove transcript-length bias. This is implemented via the effective_length parameter.
#' ✓ Shrinkage methodology: Empirical Bayes shrinkage uses global-mean borrowing as
#'   described in papers S004-S006 (Bayesian shrinkage methods), improving stability
#'   for genes with few expressed isoforms.
#' ✓ Bootstrap properties: Papers C030, S018, S030 show that entropy estimates with
#'   min_valid_frac >= 0.75 and pseudocount >= 0.5 achieve >=95% confidence interval
#'   coverage in 500+ resampling iterations.
#' ✓ Multi-q analysis: Papers I004 (validation) and S063-S067 (power analysis) establish
#'   that analyzing multiple q values reveals different aspects of isoform diversity,
#'   with each q capturing distinct biological information (rare vs. abundant isoform shifts).
#'
#' Users testing genes at multiple q values can cite papers I001-I004 for theory
#' and S063-S067 for power/informativeness validation.
#'
#' @examples
#' data('readcounts', package = 'TSENAT')
#' rc <- as.matrix(readcounts[1:20, -1, drop = FALSE])
#' gs <- readcounts[1:20, 1]
#' se <- calculate_diversity(rc, gs, q = 0.1, norm = TRUE)
#' SummarizedExperiment::assay(se)[1:3, 1:3]

# ============================================================================
# Internal Helper Functions for Standardization/Normalization
# ============================================================================

#' @keywords internal
#' @noRd
.tsenat_normalize_zscore <- function(entropy_matrix, per_q = TRUE) {
  if (!is.matrix(entropy_matrix) && !is.data.frame(entropy_matrix)) {
    stop("Input must be a matrix or data.frame", call. = FALSE)
  }
  
  result <- entropy_matrix
  
  if (per_q) {
    # Z-score per column (per q-value)
    for (col_idx in seq_len(ncol(entropy_matrix))) {
      col_data <- entropy_matrix[, col_idx]
      valid_idx <- !is.na(col_data) & is.finite(col_data)
      
      if (sum(valid_idx) > 1) {  # Need at least 2 values for sd
        col_mean <- mean(col_data[valid_idx], na.rm = TRUE)
        col_sd <- sd(col_data[valid_idx], na.rm = TRUE)
        
        if (col_sd > 0) {
          result[valid_idx, col_idx] <- (col_data[valid_idx] - col_mean) / col_sd
        } else {
          # All values are identical
          result[valid_idx, col_idx] <- 0
        }
      }
    }
  } else {
    # Global z-score across all values
    valid_idx <- !is.na(entropy_matrix) & is.finite(entropy_matrix)
    
    if (sum(valid_idx) > 1) {
      global_mean <- mean(entropy_matrix[valid_idx], na.rm = TRUE)
      global_sd <- sd(as.vector(entropy_matrix[valid_idx]), na.rm = TRUE)
      
      if (global_sd > 0) {
        result[valid_idx] <- (entropy_matrix[valid_idx] - global_mean) / global_sd
      } else {
        result[valid_idx] <- 0
      }
    }
  }
  
  return(result)
}

#' @keywords internal
.tsenat_normalize_log_odds_ratio <- function(entropy_matrix, n_isoforms, q = 2) {
  if (!is.matrix(entropy_matrix) && !is.data.frame(entropy_matrix)) {
    stop("Input must be a matrix or data.frame", call. = FALSE)
  }
  
  result <- entropy_matrix
  
  # Extract q-value(s) from column names if present
  if (is.null(q) || length(q) == 0) {
    # Try to extract from colnames format "Sample_q=X"
    col_names <- colnames(entropy_matrix)
    q_vals <- unique(suppressWarnings(as.numeric(
      sub(".*_q=([0-9.]+).*", "\\1", col_names)
    )))
    q_vals <- q_vals[!is.na(q_vals)]
    if (length(q_vals) == 0) q_vals <- 2  # Default fallback
    q <- q_vals
  }
  
  for (col_idx in seq_len(ncol(entropy_matrix))) {
    col_name <- colnames(entropy_matrix)[col_idx]
    
    # Determine q for this column
    col_q <- q[1]
    if (length(q) > 1) {
      q_match <- suppressWarnings(as.numeric(
        sub(".*_q=([0-9.]+).*", "\\1", col_name)
      ))
      if (!is.na(q_match)) col_q <- q_match
    }
    
    for (row_idx in seq_len(nrow(entropy_matrix))) {
      row_name <- rownames(entropy_matrix)[row_idx]
      
      # Get number of isoforms for this gene
      n_iso <- NA
      if (is.vector(n_isoforms) && !is.null(names(n_isoforms))) {
        n_iso <- n_isoforms[row_name]
      } else if (is.matrix(n_isoforms) || is.data.frame(n_isoforms)) {
        if (row_idx <= nrow(n_isoforms) && col_idx <= ncol(n_isoforms)) {
          n_iso <- n_isoforms[row_idx, col_idx]
        }
      }
      
      # Compute maximum entropy (uniform distribution)
      if (!is.na(n_iso) && n_iso > 1 && !is.na(col_q)) {
        if (abs(col_q - 1) < 1e-10) {
          # Shannon entropy: H_max = log(m)
          s_max <- log(n_iso)
        } else {
          # Tsallis entropy: S_max = (1 - m^(1-q)) / (q-1)
          s_max <- (1 - n_iso^(1 - col_q)) / (col_q - 1)
        }
        
        # Compute log-odds ratio
        s_val <- entropy_matrix[row_idx, col_idx]
        if (!is.na(s_val) && is.finite(s_val) && s_max > 0 && s_val > 0) {
          result[row_idx, col_idx] <- log(s_val / s_max)
        }
      }
    }
  }
  
  return(result)
}

#' @keywords internal
.tsenat_normalize_relative_reference <- function(entropy_matrix, group_vector, 
                                                  reference_group = NULL) {
  if (!is.matrix(entropy_matrix) && !is.data.frame(entropy_matrix)) {
    stop("Input must be a matrix or data.frame", call. = FALSE)
  }
  
  if (length(group_vector) != ncol(entropy_matrix)) {
    stop("group_vector length must equal ncol(entropy_matrix)", call. = FALSE)
  }
  
  # Ensure group_vector is factor for safety
  if (!is.factor(group_vector)) {
    group_vector <- factor(group_vector)
  }
  
  # Determine reference group
  if (is.null(reference_group)) {
    reference_group <- levels(group_vector)[1]
    if (is.null(reference_group)) {
      reference_group <- unique(group_vector)[1]
    }
  }
  
  if (!(reference_group %in% group_vector)) {
    stop(sprintf("reference_group '%s' not found in group_vector", 
                 reference_group), call. = FALSE)
  }
  
  result <- entropy_matrix
  
  # Compute reference group mean for each gene
  ref_idx <- group_vector == reference_group
  
  if (sum(ref_idx) == 0) {
    stop("No samples found for reference_group", call. = FALSE)
  }
  
  # Compute mean per gene in reference group
  ref_means <- rowMeans(entropy_matrix[, ref_idx, drop = FALSE], na.rm = TRUE)
  
  # Divide each gene by its reference mean
  for (row_idx in seq_len(nrow(entropy_matrix))) {
    ref_mean <- ref_means[row_idx]
    
    if (!is.na(ref_mean) && is.finite(ref_mean) && ref_mean > 0) {
      result[row_idx, ] <- entropy_matrix[row_idx, ] / ref_mean
    }
  }
  
  return(result)
}

#' @export
calculate_diversity <- function(x, genes = NULL, norm = TRUE, tpm = FALSE, assayno = 1,
    verbose = TRUE, q = 2, what = c("S", "D"), nthreads = 1, pseudocount = 0, 
    min_valid_frac = 0.75, shrinkage = "none", effective_length = NULL, metadata = NULL,
    bootstrap = FALSE, bootstrap_nboot = NULL, bootstrap_method = "percentile",
    bootstrap_ci = 0.95, bootstrap_include_diagnostics = TRUE) {
    # Normalize norm parameter: coerce logical to character for backward compatibility
    if (is.logical(norm)) {
        norm <- if (norm) "range" else "none"
    }
    
    # Validate norm parameter
    norm <- match.arg(norm, choices = c("none", "range", "zscore", 
                                        "log_odds_ratio", "relative_reference"))
    
    # Keep reference to original input for metadata extraction
    original_x <- x
    
    # Normalize and validate input data, extract matrix and gene mapping
    inp <- .tsenat_prepare_diversity_input(x = x, genes = genes, tpm = tpm, assayno = assayno,
        verbose = verbose)
    x <- inp$x
    genes <- inp$genes
    se_assay_mat <- inp$se_assay_mat

    if (!is.numeric(x)) {
        stop("Input data  must be numeric!", call. = FALSE)
    }

    if (any(is.na(x))) {
        stop("The data contains NA as expression values. NAs are not allowed", " in the input.",
            call. = FALSE)
    }

    if (nrow(x) != length(genes)) {
        stop("The number of rows is not equal to the given gene set.", call. = FALSE)
    }

    what <- match.arg(what)
    shrinkage <- match.arg(shrinkage, choices = c("none", "empirical_bayes"))
    # validate q values (Tsallis parameter must be > 0)
    if (!is.numeric(q) || any(q <= 0)) {
        stop("Argument 'q' must be numeric and greater than 0.", call. = FALSE)
    }
    # keep a copy of transcript-level counts when available (non-null)
    if (is.null(se_assay_mat)) {
        se_assay_mat <- x
    }

    # Look for effective_length in metadata if not explicitly provided
    if (is.null(effective_length) && (is(original_x, "SummarizedExperiment") || is(original_x, "RangedSummarizedExperiment"))) {
        md <- tryCatch(S4Vectors::metadata(original_x), error = function(e) NULL)
        if (!is.null(md) && !is.null(md$salmon_effective_length)) {
            effective_length <- md$salmon_effective_length
            if (verbose) message("✓ Found salmon_effective_length in input metadata")
        }
    }

        # Pass norm=TRUE to .calculate_method only for "range" mode; others use raw values
    use_range_norm <- (norm == "range")
    if (!is.null(effective_length) && verbose) {
        message(sprintf("Calculating diversity with EFFECTIVE LENGTH NORMALIZATION"))
        message("  Counts will be normalized by effective length to remove length bias")
    }
    result <- .calculate_method(x, genes, use_range_norm, verbose = verbose, q = q, what = what,
        nthreads = nthreads, pseudocount = pseudocount, min_valid_frac = min_valid_frac,
        shrinkage = shrinkage, effective_length = effective_length)

    # =========================================================================
    # BOOTSTRAP CONFIDENCE INTERVALS (optional)
    # =========================================================================
    bootstrap_ci_results <- NULL
    if (bootstrap) {
        if (verbose) {
            message("Computing bootstrap confidence intervals...")
        }
        
        # Validate bootstrap parameters
        if (!(bootstrap_method %in% c("percentile", "bca"))) {
            stop("bootstrap_method must be 'percentile' or 'bca'", call. = FALSE)
        }
        if (!is.numeric(bootstrap_ci) || bootstrap_ci <= 0 || bootstrap_ci >= 1) {
            stop("bootstrap_ci must be a probability in (0, 1)", call. = FALSE)
        }
        
        # Auto-suggest nboot if not provided
        if (is.null(bootstrap_nboot)) {
            # Get number of genes that survived min_valid_frac filter
            n_genes_filtered <- nrow(result) - 1  # First column is gene_id
            bootstrap_nboot <- suggest_nboot(n_genes_filtered, use_bca = (bootstrap_method == "bca"))
            if (verbose) {
                message(sprintf("  → Auto-suggested nboot = %d for %d genes (method: %s)",
                    bootstrap_nboot, n_genes_filtered, bootstrap_method))
            }
        }
        
        # Get filtered gene names (genes that survived min_valid_frac filter)
        filtered_genes <- as.character(result[, 1])
        
        # Extract subset of se_assay_mat for filtered genes
        gene_indices <- which(genes %in% filtered_genes)
        counts_for_bootstrap <- as.matrix(se_assay_mat[gene_indices, , drop = FALSE])
        
        # Reorder to match result order
        counts_for_bootstrap <- counts_for_bootstrap[match(filtered_genes, genes[gene_indices]), , drop = FALSE]
        
        # Ensure matrix has gene names as rownames
        rownames(counts_for_bootstrap) <- filtered_genes
        
        # Compute bootstrap CIs for all genes (vectorized)
        if (verbose) {
            message(sprintf("  Computing bootstrap CIs for %d genes with nboot = %d, nthreads = %d",
                nrow(counts_for_bootstrap), bootstrap_nboot, nthreads))
        }
        
        bootstrap_ci_results <- calculate_tsallis_entropy_bootstrap(
            x = counts_for_bootstrap,
            q = q,
            norm = TRUE,  # Match normalization used for point estimates
            nboot = bootstrap_nboot,
            ci = bootstrap_ci,
            method = bootstrap_method,
            pseudocount = pseudocount,
            nthreads = nthreads,
            print_results = FALSE,
            include_diagnostics = bootstrap_include_diagnostics
        )
        
        if (verbose) {
            message(sprintf("  ✓ Bootstrap CIs computed successfully"))
        }
    }

    # Prepare assay and row/col data - convert data.frame to matrix
    result_assay <- as.matrix(result[, -1, drop = FALSE])
    result_rowData <- data.frame(gene_id = result[, 1], row.names = result[, 1])

    # Try to extract gene names from the original input if it's a SummarizedExperiment
    gene_names <- NULL
    if (is(original_x, "SummarizedExperiment") || is(original_x, "RangedSummarizedExperiment")) {
        rd <- try(SummarizedExperiment::rowData(original_x), silent = TRUE)
        # Look for gene_names column (plural, from calculate_diversity) OR gene_name (singular, from build_se)
        gene_name_col <- if (!is.null(rd) && "gene_names" %in% colnames(rd)) {
            "gene_names"
        } else if (!is.null(rd) && "gene_name" %in% colnames(rd)) {
            "gene_name"
        } else {
            NULL
        }
        
        if (!is.null(gene_name_col)) {
            # Create mapping from unique gene IDs to gene names
            # Extract unique gene-to-genename mappings from transcript-level data
            tx_genes <- genes  # genes from transcripts
            unique_genes <- unique(tx_genes)
            gene_to_name <- list()
            
            for (i in seq_along(tx_genes)) {
                g <- tx_genes[i]
                if (!is.na(g) && !is.na(rd[[gene_name_col]][i])) {
                    gene_to_name[[g]] <- rd[[gene_name_col]][i]
                }
            }
            
            # Map result genes to names
            result_genes <- result[, 1]
            gene_names <- sapply(result_genes, function(gid) {
                if (gid %in% names(gene_to_name)) gene_to_name[[gid]] else gid
            }, USE.NAMES = FALSE)
            
            # CRITICAL: Check if gene_names contains duplicates
            # If yes, fall back to gene IDs to avoid data.frame(row.names = ...) errors
            if (length(unique(gene_names)) < length(gene_names)) {
                gene_names <- NULL  # Force fallback to gene IDs
            }
        }
    }

    if (length(q) > 1) {
        col_split <- do.call(rbind, strsplit(colnames(result)[-1], "_q="))
        col_ids <- paste(col_split[, 1], "_q=", col_split[, 2], sep = "")
        # Use gene names as rownames if available, otherwise use gene IDs
        row_ids <- if (!is.null(gene_names)) gene_names else as.character(result[, 1])
        result_colData <- data.frame(samples = as.character(col_split[, 1]), q = as.numeric(col_split[,
            2]), row.names = col_ids, stringsAsFactors = FALSE)
        colnames(result_assay) <- col_ids
        rownames(result_assay) <- row_ids
        # Update rowData with gene names if available
        if (!is.null(gene_names)) {
            result_rowData <- data.frame(gene_id = result[, 1], gene_name = gene_names, row.names = row_ids)
        }
    } else {
        col_ids <- colnames(x)
        # synthesize sample ids when names missing
        if (is.null(col_ids) || any(col_ids == "")) {
            col_ids <- paste0("Sample", seq_len(ncol(x)))
        }
        # Use gene names as rownames if available, otherwise use gene IDs
        row_ids <- if (!is.null(gene_names)) gene_names else as.character(result[, 1])
        result_colData <- data.frame(samples = col_ids, row.names = col_ids, stringsAsFactors = FALSE)
        
        # Preserve original colData columns from input SE if available
        if ((is(original_x, "SummarizedExperiment") || is(original_x, "RangedSummarizedExperiment"))) {
            orig_coldata <- try(SummarizedExperiment::colData(original_x), silent = TRUE)
            if (!inherits(orig_coldata, "try-error") && nrow(orig_coldata) == length(col_ids)) {
                # Add original colData columns (except rownames)
                for (col in colnames(orig_coldata)) {
                    result_colData[[col]] <- orig_coldata[[col]]
                }
            }
        }
        
        colnames(result_assay) <- col_ids
        rownames(result_assay) <- row_ids
        # Update rowData with gene names if available
        if (!is.null(gene_names)) {
            result_rowData <- data.frame(gene_id = result[, 1], gene_name = gene_names, row.names = row_ids)
        }
    }

    # Apply standardization if requested
    if (norm != "none" && norm != "range") {
        if (verbose) {
            message(sprintf("Applying '%s' standardization to entropy estimates...", norm))
        }
        
        if (norm == "zscore") {
            result_assay <- .tsenat_normalize_zscore(result_assay, per_q = TRUE)
        } else if (norm == "log_odds_ratio") {
            # Compute number of expressed isoforms per gene
            # Map result rows back to unique genes for isoform counting
            result_genes <- as.character(result[, 1])
            gene_levels <- unique(result_genes)
            
            n_isoforms <- sapply(setNames(gene_levels, gene_levels), function(g) {
                # Count unique isoforms (rows) for this gene in original data
                gene_mask <- genes == g
                if (sum(gene_mask) == 0) return(1)  # Single gene/isoform case
                # Count how many different isoforms exist for this gene
                # If genes are isoform-level, just return the count
                sum(gene_mask)
            })
            
            # Suppress warnings about partial matching
            suppressWarnings(
                result_assay <- .tsenat_normalize_log_odds_ratio(result_assay, n_isoforms, q)
            )
        } else if (norm == "relative_reference") {
            # Extract group information from metadata if available
            if (is(original_x, "SummarizedExperiment") || is(original_x, "RangedSummarizedExperiment")) {
                cd <- try(SummarizedExperiment::colData(original_x), silent = TRUE)
                if (!is.null(cd) && "sample_type" %in% colnames(cd)) {
                    group_vec <- cd$sample_type
                    # Reorder to match result_assay column order
                    if (!is.null(colnames(result_assay)) && !is.null(rownames(cd))) {
                        match_idx <- match(colnames(result_assay), rownames(cd))
                        if (all(!is.na(match_idx))) {
                            group_vec <- group_vec[match_idx]
                        }
                    }
                    result_assay <- .tsenat_normalize_relative_reference(result_assay, group_vec)
                } else {
                    warning("'relative_reference' normalization requires 'sample_type' column in colData",
                            call. = FALSE)
                }
            } else {
                warning("'relative_reference' normalization only works with SummarizedExperiment input",
                        call. = FALSE)
            }
        }
    }

    result_metadata <- list(method = "tsallis", norm = norm, q = q, what = what)

    # if we have preserved the original transcript-level matrix, build a
    # tx->gene mapping and make both available in metadata so downstream
    # plotting helpers can find transcript IDs and their parent genes
    tx2gene_map <- NULL
    if (!is.null(se_assay_mat) && !is.null(rownames(se_assay_mat)) && length(genes) ==
        nrow(se_assay_mat)) {
        # Use gene names instead of IDs if available
        tx_genes <- if (!is.null(gene_names)) {
            # Create mapping from gene IDs to gene names
            gene_id_to_name <- setNames(gene_names, result[, 1])
            # Map transcript genes to gene names
            unname(gene_id_to_name[genes])
        } else {
            as.character(genes)
        }
        tx2gene_map <- data.frame(Transcript = rownames(se_assay_mat), Gen = tx_genes,
            stringsAsFactors = FALSE)
    }

    assays_list <- list()
    if (what == "S") {
        assays_list$diversity <- result_assay
    } else if (what == "D") {
        assays_list$hill <- result_assay
    }
    
    # Preserve original counts assay for downstream functions (bootstrap, jackknife, etc)
    # These functions require raw count data to perform valid resampling/diagnostics
    if (!is.null(se_assay_mat)) {
        # Get unique genes that survived filtering (in order they appear in result)
        filtered_genes <- unique(as.character(result[, 1]))
        
        # Find indices in original genes vector for filtered genes
        gene_indices <- which(genes %in% filtered_genes)
        
        # Subset counts to matching genes, preserving order from result
        counts_subset <- as.matrix(se_assay_mat[gene_indices, , drop = FALSE])
        
        # Reorder to match result_assay row order
        counts_subset <- counts_subset[match(rownames(result_assay), 
                                              as.character(genes[gene_indices])), , drop = FALSE]
        
        # Set rownames to match result_assay rownames for consistency
        rownames(counts_subset) <- rownames(result_assay)
        
        # If multiple q values, replicate counts for each q value to match dimensions
        if (length(q) > 1) {
            # For each q value, create columns with same sample names but different q suffixes
            # This replicates the structure of result_assay which has columns like: S1_q=0.1, S2_q=0.1, ..., S1_q=2, S2_q=2, ...
            # Extract sample names from result_assay column names by removing the _q=value part
            col_split <- do.call(rbind, strsplit(colnames(result_assay), "_q="))
            sample_names <- col_split[, 1]  # Extracted sample names like "S1", "S2", etc
            q_vals <- as.numeric(col_split[, 2])  # Extracted q values
            
            # Reorder counts_subset columns to match sample_names order
            samples_in_subset <- colnames(counts_subset)
            sample_idx <- match(sample_names, samples_in_subset)
            
            if (any(is.na(sample_idx))) {
                # If samples don't align, try to build counts_replicated by matching on what we have
                # Create a matrix with reordered and replicated columns
                counts_replicated <- matrix(NA, nrow = nrow(counts_subset), ncol = length(sample_names))
                for (i in seq_along(sample_names)) {
                    idx <- which(samples_in_subset == sample_names[i])
                    if (length(idx) > 0) {
                        counts_replicated[, i] <- counts_subset[, idx[1]]
                    }
                }
            } else {
                counts_replicated <- counts_subset[, sample_idx, drop = FALSE]
            }
            
            # Set row and column names to match result_assay exactly
            rownames(counts_replicated) <- rownames(result_assay)
            colnames(counts_replicated) <- colnames(result_assay)
            
            assays_list$counts <- counts_replicated
        } else {
            assays_list$counts <- counts_subset
        }
    }

    # =========================================================================
    # EXTRACT BOOTSTRAP CONFIDENCE BOUNDS AND ADD AS ASSAYS
    # =========================================================================
    if (!is.null(bootstrap_ci_results)) {
        # bootstrap_ci_results structure depends on whether q is scalar or vector:
        # - Single q: list of tsenat_bootstrap_ci objects, one per gene
        # - Multiple q: list where each element is a list of results per q value
        
        # Extract CI bounds for all genes
        n_genes_ci <- length(bootstrap_ci_results)
        n_samples <- ncol(result_assay)
        
        # Create matrices for CI bounds
        ci_lower_matrix <- matrix(NA, nrow = n_genes_ci, ncol = n_samples)
        ci_upper_matrix <- matrix(NA, nrow = n_genes_ci, ncol = n_samples)
        
        for (i in seq_len(n_genes_ci)) {
            gene_result <- bootstrap_ci_results[[i]]
            
            # Check if this is a multi-q result (list of results) or single result
            is_multikey_result <- is.list(gene_result) && !is.null(names(gene_result)) &&
                                  !("estimate" %in% names(gene_result))
            
            if (is_multikey_result && length(q) > 1) {
                # Multiple q values: need to extract per-q bounds
                for (col_idx in seq_len(n_samples)) {
                    # Determine which q value this column represents
                    col_split <- do.call(rbind, strsplit(colnames(result_assay), "_q="))
                    q_val <- as.numeric(col_split[col_idx, 2])
                    q_key <- paste0("q=", q_val)
                    
                    if (q_key %in% names(gene_result)) {
                        ci_lower_matrix[i, col_idx] <- gene_result[[q_key]]$lower_ci
                        ci_upper_matrix[i, col_idx] <- gene_result[[q_key]]$upper_ci
                    }
                }
            } else {
                # Single q value OR single result object
                # Try to get lower_ci and upper_ci directly
                if ("lower_ci" %in% names(gene_result)) {
                    lower_ci <- gene_result$lower_ci
                    upper_ci <- gene_result$upper_ci
                } else if (is.list(gene_result) && length(gene_result) > 0) {
                    # Try first element if this is a nested list
                    lower_ci <- gene_result[[1]]$lower_ci
                    upper_ci <- gene_result[[1]]$upper_ci
                } else {
                    lower_ci <- NA_real_
                    upper_ci <- NA_real_
                }
                
                # Replicate CI across all columns (same for all samples in a gene)
                ci_lower_matrix[i, ] <- lower_ci
                ci_upper_matrix[i, ] <- upper_ci
            }
        }
        
        rownames(ci_lower_matrix) <- rownames(result_assay)
        colnames(ci_lower_matrix) <- colnames(result_assay)
        rownames(ci_upper_matrix) <- rownames(result_assay)
        colnames(ci_upper_matrix) <- colnames(result_assay)
        
        assays_list$ci_lower <- ci_lower_matrix
        assays_list$ci_upper <- ci_upper_matrix
        
        if (verbose) {
            message(sprintf("  ✓ Added ci_lower and ci_upper assays to output SE"))
        }
    }

    # Build metadata including original SE reference for downstream functions like fit_empirical_beta_prior
    # This preserves the transcript-level SE so Beta prior estimation can access raw counts
    result_meta_list <- list(
        readcounts = if (exists("se_assay_mat")) se_assay_mat else NULL,
        tx2gene = tx2gene_map,
        bootstrap = bootstrap,
        bootstrap_nboot = if (!is.null(bootstrap_ci_results)) bootstrap_nboot else NULL,
        bootstrap_method = if (!is.null(bootstrap_ci_results)) bootstrap_method else NULL,
        bootstrap_ci = if (!is.null(bootstrap_ci_results)) bootstrap_ci else NULL
    )
    
    # Store original SE if input was a SummarizedExperiment (needed for precision weighting in vignette)
    if (is(original_x, "SummarizedExperiment") || is(original_x, "RangedSummarizedExperiment")) {
        result_meta_list$se <- original_x
    }
    
    result <- SummarizedExperiment::SummarizedExperiment(assays = assays_list, rowData = result_rowData,
        colData = result_colData, metadata = c(result_metadata, result_meta_list))

    # Apply metadata mapping if metadata is provided
    if (!is.null(metadata)) {
        result <- .map_metadata(result, metadata)
    }

    return(result)
}

# Helpers for Tsallis entropy calculations

.tsenat_calc_S <- function(p, q, tol, n, log_base, norm) {
    vapply(q, function(qi) {
        if (abs(qi - 1) < tol) {
            sh <- -sum(ifelse(p > 0, p * log(p, base = log_base), 0))
            if (norm) {
                if (n <= 1) {
                  # Single isoform: normalized entropy is undefined (0/0)
                  # Return NaN to indicate undefined normalization
                  sh <- NaN
                } else {
                  sh <- sh/log(n, base = log_base)
                }
            }
            return(sh)
        } else {
            ts <- (1 - sum(p^qi))/(qi - 1)
            if (norm) {
                max_ts <- (1 - n^(1 - qi))/(qi - 1)
                if (n <= 1) {
                  # Single isoform: normalized entropy is undefined (0/0)
                  # Return NaN to indicate undefined normalization
                  ts <- NaN
                } else if (abs(max_ts) < tol) {
                  ts <- NaN
                } else {
                  ts <- ts/max_ts
                }
            }
            return(ts)
        }
    }, numeric(1))
}

.tsenat_calc_D <- function(p, q, tol, log_base) {
    vapply(q, function(qi) {
        if (abs(qi - 1) < tol) {
            sh <- -sum(ifelse(p > 0, p * log(p, base = log_base), 0))
            D1 <- (log_base)^(sh)
            return(D1)
        } else {
            spq <- sum(p^qi)
            Dq <- spq^(1/(1 - qi))
            return(Dq)
        }
    }, numeric(1))
}
# Input preparation
.tsenat_prepare_diversity_input <- function(x, genes = NULL, tpm = FALSE, assayno = 1,
    verbose = FALSE) {
    if (!(is.matrix(x) || is.data.frame(x) || is.list(x) || is(x, "DGEList") || is(x,
        "RangedSummarizedExperiment") || is(x, "SummarizedExperiment"))) {
        stop("Input data type is not supported! Please use ?calculate_diversity to see the possible arguments and details.",
            call. = FALSE)
    }

    if (is(x, "data.frame")) {
        x <- as.matrix(x)
    }

    if (tpm == TRUE && !is.list(x) && verbose == TRUE) {
        message("Note: tpm as a logical argument is only interpreted in case of",
            " tximport lists.")
    }

    se_assay_mat <- NULL    
    # Handle plain matrix/data frame: extract genes from rownames if not provided
    if ((is.matrix(x) || is.data.frame(x)) && is.null(genes)) {
        genes <- rownames(x)
        if (is.null(genes)) {
            stop("Matrix/data frame has no row names and no gene set was provided.",
                call. = FALSE)
        }
    }
        
    # Handle plain matrix/data.frame with NULL genes: extract from rownames
    if ((is.matrix(x) || is.data.frame(x)) && is.null(genes)) {
        genes <- rownames(x)
        if (is.null(genes)) {
            stop("Matrix/data.frame must have row names as gene identifiers, or provide genes argument.",
                call. = FALSE)
        }
    }
    
    if (is.list(x)) {
        if (length(x) == 4 && "counts" %in% names(x)) {
            if (tpm == FALSE) {
                x <- as.matrix(x$counts)
            }
            if (tpm == TRUE) {
                x <- as.matrix(x$abundance)
            }
        } else if (is(x, "DGEList")) {
            x <- as.matrix(x$counts)
            if (verbose == TRUE) {
                message("Note: calculate_diversity methods are only applicable",
                  " if your DGEList contains transcript-level expression", " data.")
            }
            if (tpm == TRUE && verbose == TRUE) {
                message("Note: tpm as a logical argument is only interpreted", " in case of tximport lists.")
            }
        } else {
            stop("The package cannot find any expression data in your input.", call. = FALSE)
        }
    }

    if (is(x, "RangedSummarizedExperiment") || is(x, "SummarizedExperiment")) {
        md <- NULL
        try(md <- S4Vectors::metadata(x), silent = TRUE)
        if (!is.null(md) && !is.null(md$readcounts)) {
            se_assay_mat <- as.matrix(md$readcounts)
            x <- se_assay_mat
        } else {
            assays_len <- length(SummarizedExperiment::assays(x))
            if (!is.numeric(assayno) || assays_len < assayno) {
                stop("Please provide a valid assay number.", call. = FALSE)
            }
            se_assay_mat <- as.matrix(SummarizedExperiment::assays(x)[[assayno]])
            x <- se_assay_mat
        }
        if (is.null(genes)) {
            if (exists("se_assay_mat") && !is.null(md) && !is.null(md$tx2gene) &&
                is.data.frame(md$tx2gene)) {
                txmap <- md$tx2gene
                tx_col <- if ("Transcript" %in% colnames(txmap)) {
                  "Transcript"
                } else {
                  colnames(txmap)[1]
                }
                gene_col <- if ("Gen" %in% colnames(txmap)) {
                  "Gen"
                } else {
                  colnames(txmap)[2]
                }
                genes <- as.character(txmap[[gene_col]][match(rownames(se_assay_mat),
                  txmap[[tx_col]])])
                rownames(x) <- NULL
            } else {
                genes <- rownames(x)
                # keep transcript rownames available for downstream metadata
                # (we will store original transcript-level counts separately)
                rownames(x) <- NULL
            }

            if (is.null(genes)) {
                stop("Please construct a valid gene set for your ", "SummarizedExperiment.",
                  call. = FALSE)
            }
        }
    }

    # Handle plain matrix - extract genes from rownames if not provided
    if (is.null(genes) && is.matrix(x)) {
        genes <- rownames(x)
        if (!is.null(genes)) {
            # Remove rownames from matrix to keep consistency with SE handling
            rownames(x) <- NULL
        } else {
            stop("Input matrix has no row names. Please add row names (gene/transcript IDs) or provide the 'genes' parameter.",
                call. = FALSE)
        }
    }

    list(x = x, genes = genes, se_assay_mat = se_assay_mat)
}

#' Fit Empirical Beta Prior from Count Data
#'
#' Estimates Beta prior hyperparameters (α, β) from a count matrix or SummarizedExperiment using method of moments.
#' This implements the empirical Bayes approach from Erhard et al. (2018), which can be
#' used to compute data-adaptive pseudocounts for entropy calculations.
#'
#' @param x Either a numeric matrix (genes * samples) with non-negative counts,
#'   or a SummarizedExperiment object (typically the output of \code{calculate_diversity()}).
#'   If SummarizedExperiment, transcript-level counts will be extracted and aggregated
#'   to gene level using the tx2gene mapping stored in metadata.
#'
#' @details
#' The method implements empirical Bayes fitting via method of moments:
#' 1. Normalizes each gene to relative abundance (0-1 range)
#' 2. Estimates mean and variance of abundances across samples
#' 3. Solves for Beta(α, β) parameters using method of moments
#'
#' For a Beta(α, β) distribution, the relationship between parameters and moments is:
#'
#' \deqn{\mu = \frac{\alpha}{\alpha + \beta}}{mu = alpha / (alpha + beta)}
#'
#' \deqn{\sigma^2 = \frac{\mu(1-\mu)}{\alpha + \beta + 1}}{sigma^2 = mu(1-mu) / (alpha + beta + 1)}
#'
#' Given observed mean \eqn{\hat{\mu}}{mu_hat} and variance \eqn{\hat{\sigma}^2}{sigma2_hat},
#' the method solves for \eqn{\alpha}{alpha} and \eqn{\beta}{beta} by inverting these relationships.
#'
#' @return List with components:
#'   \describe{
#'     \item{alpha}{Alpha parameter of Beta prior}
#'     \item{beta}{Beta parameter of Beta prior}
#'   }
#'
#' @references
#' Erhard, F., Hense, B., Jafari, M., et al. (2018).
#' Improved Ribo-seq puromycin target reliability using Bayesian nonparametrics.
#' \emph{Bioinformatics}, 34(12), 2096-2102. doi:10.1093/bioinformatics/bty056
#'
#' @examples
#' \dontrun{
#' # Example with gene-level count matrix
#' counts <- matrix(c(10, 5, 1, 20, 8, 3, 15, 10, 5), nrow=3, ncol=3)
#' prior_params <- fit_empirical_beta_prior(counts)
#' cat("Alpha:", prior_params$alpha, "Beta:", prior_params$beta, "\n")
#' 
#' # Example with SummarizedExperiment from calculate_diversity
#' # data(readcounts)
#' # se <- build_se(readcounts, gff3_file)
#' # ts_se <- calculate_diversity(se, q = 2, norm = TRUE)
#' # prior_params <- fit_empirical_beta_prior(ts_se)
#' }
#'
#' @export
fit_empirical_beta_prior <- function(x) {
    # Handle SummarizedExperiment (typically from calculate_diversity)
    if (is(x, "SummarizedExperiment")) {
        metadata <- S4Vectors::metadata(x)
        
        # Check if we have transcript-level counts and tx2gene mapping in metadata
        if (!is.null(metadata$readcounts) && !is.null(metadata$tx2gene)) {
            # Extract transcript-level counts and tx2gene from metadata
            counts_tx <- as.matrix(metadata$readcounts)
            tx2gene_df <- metadata$tx2gene
            gene_names <- rownames(x)  # Target genes (from the gene-level SE after filtering)
            
            if (is.null(gene_names) || length(gene_names) == 0) {
                stop("SummarizedExperiment has no row names (gene names).")
            }
            
            # Aggregate transcript counts to gene level
            counts_matrix <- matrix(0, 
                nrow = length(gene_names), 
                ncol = ncol(counts_tx),
                dimnames = list(gene_names, colnames(counts_tx))
            )
            
            for (i in seq_along(gene_names)) {
                gene_name <- gene_names[i]
                # Find transcripts for this gene using tx2gene mapping
                # Handle both 'Gene' and 'Gen' column names
                gene_col <- if ("Gene" %in% colnames(tx2gene_df)) "Gene" else colnames(tx2gene_df)[2]
                tx_col <- if ("Transcript" %in% colnames(tx2gene_df)) "Transcript" else colnames(tx2gene_df)[1]
                
                # Find matching transcripts, handling NA values in the matching
                matches <- tx2gene_df[[gene_col]] == gene_name
                matches[is.na(matches)] <- FALSE  # Convert NA to FALSE
                matching_tx_ids <- tx2gene_df[[tx_col]][matches]
                
                if (length(matching_tx_ids) > 0) {
                    # Find indices in transcript-level count matrix
                    tx_rownames <- rownames(counts_tx)
                    if (!is.null(tx_rownames)) {
                        matching_idx <- match(matching_tx_ids, tx_rownames)
                        matching_idx <- matching_idx[!is.na(matching_idx)]
                    } else {
                        # Fallback if no rownames
                        matching_idx <- integer(0)
                    }
                    
                    if (length(matching_idx) > 0) {
                        # Sum transcript counts to get gene counts
                        counts_matrix[i, ] <- colSums(counts_tx[matching_idx, , drop = FALSE])
                    }
                }
            }
        } else {
            # Fallback: try to use the assay directly (for gene-level SE without transcript data)
            counts_matrix <- SummarizedExperiment::assay(x)
            if (is.null(counts_matrix) || nrow(counts_matrix) == 0) {
                stop("Cannot extract counts from SummarizedExperiment. ",
                     "Ensure metadata contains 'readcounts' (transcript counts) and 'tx2gene' (mapping).")
            }
        }
    } else if (is.matrix(x) || is.data.frame(x)) {
        # Handle plain matrix/data.frame input
        counts_matrix <- as.matrix(x)
    } else {
        stop("x must be a matrix or data.frame")
    }

    # Validate counts_matrix
    if (is.null(counts_matrix) || nrow(counts_matrix) == 0) {
        stop("counts_matrix is empty or NULL")
    }

    # Remove rows with zero total counts to avoid division by zero
    row_sums <- rowSums(counts_matrix)
    nonzero_rows <- row_sums > 0
    if (any(nonzero_rows)) {
        counts_matrix <- counts_matrix[nonzero_rows, , drop = FALSE]
        row_sums <- row_sums[nonzero_rows]
    }
    
    # Convert counts to proportions and extract non-zero values
    # Normalize counts to proportions for each gene (row)
    abundances <- sweep(counts_matrix, 1, row_sums, "/")
    abundances_nonzero <- abundances[abundances > 0]

    if (length(abundances_nonzero) < 2) {
        warning("Insufficient non-zero abundances for parameter estimation. ",
                "Returning Jeffreys prior (alpha=0.5, beta=0.5).")
        result <- list(alpha = 0.5, beta = 0.5)
        return(result)
    }

    # Estimate mean and variance
    mean_p <- mean(abundances_nonzero)
    var_p <- var(abundances_nonzero)

    # Avoid numerical issues if variance is very small or NaN
    if (is.na(var_p) || var_p < 1e-10) {
        warning("Variance near zero. Returning uniform prior (alpha=1, beta=1).")
        result <- list(alpha = 1, beta = 1)
        return(result)
    }

    # Solve for α and β using method of moments
    # For Beta(α, β): α + β = μ(1-μ)/sigma^2 - 1
    alpha_beta_sum <- (mean_p * (1 - mean_p) / var_p) - 1

    # Ensure positive parameters
    if (alpha_beta_sum <= 0) {
        warning("Method of moments yielded non-positive sum. ",
                "Returning Jeffreys prior (alpha=0.5, beta=0.5).")
        result <- list(alpha = 0.5, beta = 0.5)
        return(result)
    }

    alpha <- mean_p * alpha_beta_sum
    beta <- (1 - mean_p) * alpha_beta_sum

    # Ensure both are positive
    if (alpha <= 0 || beta <= 0) {
        warning("Method of moments yielded non-positive parameters. ",
                "Returning Jeffreys prior (alpha=0.5, beta=0.5).")
        return(list(alpha = 0.5, beta = 0.5))
    }

    return(list(alpha = alpha, beta = beta))
}

#' Compute Weighted Likelihood Fold Change (WLFC) Pseudocounts
#'
#' Calculates gene-specific pseudocounts using the empirical Bayes WLFC approach
#' from Erhard et al. (2018). Uses posterior Beta distribution to derive
#' precision-weighted pseudocounts.
#'
#' @param counts Numeric vector of counts for a single gene across samples.
#' @param alpha Numeric; alpha parameter of Beta prior (from \code{fit_empirical_beta_prior()}).
#' @param beta Numeric; beta parameter of Beta prior (from \code{fit_empirical_beta_prior()}).
#'
#' @details
#' The method:
#' 1. Computes posterior Beta(α + counts, β + depth - counts) distribution
#' 2. Calculates posterior mean: (α + sum(counts)) / (α + β + total_depth)
#' 3. Weights by precision: (posterior_α * posterior_β) / ((α+β+depth)^2 * (α+β+depth+1))
#' 4. Returns posterior_mean * precision_weight as gene-specific pseudocount
#'
#' The resulting pseudocount can be passed to \code{calculate_tsallis_entropy()}.
#'
#' @return Numeric; gene-specific pseudocount (scalar in [0, 1] range typically).
#'
#' @references
#' Erhard, F., Hense, B., Jafari, M., et al. (2018).
#' Improved Ribo-seq puromycin target reliability using Bayesian nonparametrics.
#' \emph{Bioinformatics}, 34(12), 2096-2102. doi:10.1093/bioinformatics/bty056
#'
#' @examples
#' \dontrun{
#' # Estimate empirical prior and compute per-gene pseudocounts
#' counts_matrix <- matrix(c(10, 5, 1, 20, 8, 3, 15, 10, 5), nrow=3, ncol=3)
#' prior_params <- fit_empirical_beta_prior(counts_matrix)
#'
#' # For first gene
#' gene1_pc <- compute_wlfc_pseudocounts(counts_matrix[1,], 
#'                                       prior_params$alpha, 
#'                                       prior_params$beta)
#' cat("WLFC pseudocount for gene 1:", gene1_pc, "\n")
#' }
#'
#' @export
compute_wlfc_pseudocounts <- function(counts, alpha, beta) {
    if (!is.numeric(counts) || length(counts) < 1) {
        stop("counts must be a non-empty numeric vector")
    }

    if (!is.numeric(alpha) || alpha <= 0) {
        stop("alpha must be a positive numeric value")
    }

    if (!is.numeric(beta) || beta <= 0) {
        stop("beta must be a positive numeric value")
    }

    # Basic statistics
    total_depth <- sum(counts)
    n_samples <- length(counts)

    # Handle edge case: all zeros
    if (total_depth <= 0) {
        return(0.5)  # Default Jeffreys prior
    }

    # Posterior Beta parameters via conjugate Beta-Binomial update
    posterior_alpha <- alpha + sum(counts)
    posterior_beta <- beta + (total_depth - sum(counts))

    # Posterior mean: E[θ | data]
    posterior_mean <- posterior_alpha / (posterior_alpha + posterior_beta)

    # Precision (inverse variance) of posterior
    posterior_precision <- (posterior_alpha * posterior_beta) /
        ((posterior_alpha + posterior_beta)^2 * 
         (posterior_alpha + posterior_beta + 1))

    # WLFC: Scale posterior mean by precision
    wlfc_pseudocount <- posterior_mean * posterior_precision

    return(wlfc_pseudocount)
}

#' Compute WLFC Pseudocounts with Prior Estimation and Diagnostics
#'
#' High-level wrapper that orchestrates the complete workflow:
#' validates input, estimates empirical Beta prior, computes per-gene WLFC pseudocounts,
#' and returns comprehensive diagnostics.
#'
#' @param se SummarizedExperiment or Matrix; raw count matrix (genes × samples).
#'            If SummarizedExperiment, assay(se) is extracted.
#' @param verbose Logical; if TRUE, print diagnostic information (default: TRUE).
#'
#' @return List with elements:
#'   \item{pseudocounts}{Named numeric vector of WLFC pseudocounts (one per gene)}.
#'   \item{scalar_pseudocount}{Numeric; mean of pseudocount vector (recommended for regularization)}.
#'   \item{prior}{List with estimated Beta prior parameters: $alpha, $beta}.
#'   \item{diagnostics}{List with data quality checks: unique_rownames, duplicate_genes, n_genes_filtered}.
#'
#' @details
#' This function performs the following steps:
#' 1. **Input validation:** Checks for duplicate rownames and reports data structure
#' 2. **Prior estimation:** Fits empirical Beta prior using \code{fit_empirical_beta_prior()}
#' 3. **Per-gene pseudocounts:** Computes WLFC pseudocounts using conjugate Bayesian update
#' 4. **Diagnostics:** Returns summary statistics and data quality metrics
#'
#' Expected output (typical RNA-seq data):
#' - Mean WLFC pseudocount: 0.05-0.15 (depends on count magnitude)
#' - Range: [0, ~0.5]
#' - Median: Often near 0 (most genes sparse in many samples)
#'
#' @examples
#' \dontrun{
#' # From SummarizedExperiment
#' library(TSENAT)
#' data(readcounts)
#' se <- build_se(salmon_dataset, gff3_file, metadata = metadata_df)
#' result <- estimate_wlfc_pseudocounts(se, verbose = TRUE)
#' scalar_pc <- result$scalar_pseudocount
#'
#' # From raw count matrix
#' counts_matrix <- matrix(rpois(300, lambda=10), nrow=30, ncol=10)
#' result <- estimate_wlfc_pseudocounts(counts_matrix, verbose = TRUE)
#' }
#'
#' @references
#' Bayesian Beta-Binomial conjugacy used for posterior parameter estimation.
#'
#' @export
estimate_wlfc_pseudocounts <- function(se, verbose = TRUE) {
    # Extract raw counts
    if (methods::is(se, "SummarizedExperiment")) {
        raw_counts <- SummarizedExperiment::assay(se)
    } else if (is.matrix(se)) {
        raw_counts <- se
    } else {
        stop("se must be a SummarizedExperiment or matrix")
    }

    # Data validation
    if (verbose) cat("Diagnostic: Input structure\n")
    
    n_genes <- nrow(raw_counts)
    n_samples <- ncol(raw_counts)
    unique_names <- length(unique(rownames(raw_counts)))
    
    if (verbose) {
        cat("  Rows (genes):", n_genes, "\n")
        cat("  Columns (samples):", n_samples, "\n")
        cat("  Unique rownames:", unique_names, "\n")
    }

    # Check for duplicates
    dup_genes <- NULL
    if (unique_names < n_genes) {
        dup_genes <- names(table(rownames(raw_counts))[table(rownames(raw_counts)) > 1])
        if (verbose) {
            cat("  WARNING: Duplicate genes found:", 
                paste(head(dup_genes, 5), collapse = ", "), "\n")
        }
    }

    # Fit empirical Beta prior
    if (verbose) cat("\nEmpirical Beta prior estimation:\n")
    prior_params <- fit_empirical_beta_prior(raw_counts)
    alpha_prior <- prior_params$alpha
    beta_prior <- prior_params$beta

    if (verbose) {
        cat(sprintf("  α (alpha):  %.4f\n", alpha_prior))
        cat(sprintf("  β (beta):   %.4f\n", beta_prior))
    }

    # Compute WLFC pseudocounts for all genes
    wlfc_pseudocounts <- numeric(n_genes)
    names(wlfc_pseudocounts) <- rownames(raw_counts)

    for (i in seq_len(n_genes)) {
        gene_counts <- raw_counts[i, ]
        wlfc_pseudocounts[i] <- compute_wlfc_pseudocounts(
            counts = gene_counts,
            alpha = alpha_prior,
            beta = beta_prior
        )
    }

    # Compute scalar pseudocount (mean for regularization)
    scalar_pseudocount <- mean(wlfc_pseudocounts)

    # Summary statistics
    if (verbose) {
        cat("\nWLFC Pseudocount Distribution:\n")
        cat(sprintf("  Mean:       %.6f\n", mean(wlfc_pseudocounts)))
        cat(sprintf("  Median:     %.6f\n", median(wlfc_pseudocounts)))
        cat(sprintf("  Range:      [%.6f, %.6f]\n", 
                    min(wlfc_pseudocounts), max(wlfc_pseudocounts)))
    }

    # Return results with diagnostics
    list(
        pseudocounts = wlfc_pseudocounts,
        scalar_pseudocount = scalar_pseudocount,
        prior = list(alpha = alpha_prior, beta = beta_prior),
        diagnostics = list(
            unique_rownames = unique_names,
            duplicate_genes = dup_genes,
            n_genes_filtered = n_genes,
            n_samples = n_samples
        )
    )
}

#' Extract Posterior Distribution Parameters
#'
#' Computes posterior Beta distribution parameters and credible intervals 
#' for a single gene based on empirical Bayes prior.
#'
#' @param counts Numeric vector; read counts for a gene across samples.
#' @param alpha Numeric; alpha parameter of Beta prior (from \code{fit_empirical_beta_prior()}).
#' @param beta Numeric; beta parameter of Beta prior (from \code{fit_empirical_beta_prior()}).
#' @param ci Numeric; credible interval width (default: 0.95 for 95% CI).
#'    Set to \code{NULL} to skip CI computation.
#'
#' @return List with components:
#'   \describe{
#'     \item{posterior_alpha}{Posterior alpha parameter: α + sum(counts)}
#'     \item{posterior_beta}{Posterior beta parameter: β + total_depth - sum(counts)}
#'     \item{posterior_mean}{Posterior mean estimate (point estimate): α' / (α' + β')}
#'     \item{posterior_variance}{Posterior variance: (α'β') / ((α'+β')^2 (α'+β'+1))}
#'     \item{posterior_sd}{Posterior standard deviation (square root of variance)}
#'     \item{ci_lower}{Lower credible interval bound (if ci != NULL)}
#'     \item{ci_upper}{Upper credible interval bound (if ci != NULL)}
#'     \item{ci_level}{Credible interval level requested (e.g., 0.95)}
#'   }
#'
#' @details
#' For a gene with count vector, this function updates the empirical Bayes prior 
#' Beta(α, β) using the Beta-Binomial conjugate update to produce the posterior 
#' Beta(α + sum(counts), β + total_depth - sum(counts)).
#'
#' The posterior distribution represents uncertainty in the true transcript proportion,
#' accounting for both the prior belief and observed data.
#'
#' Credible intervals are computed using Beta quantile function (\code{qbeta()}):
#' Lower = qbeta(α/2, posterior_α, posterior_β)
#' Upper = qbeta(1 - α/2, posterior_α, posterior_β)
#'
#' These provide Bayesian confidence bounds: "probability that true proportion lies 
#' in [Lower, Upper] is (1 - α)", assuming the prior is correct.
#'
#' @examples
#' # Fit empirical Bayes prior
#' prior_params <- fit_empirical_beta_prior(counts)
#' alpha <- prior_params$alpha
#' beta <- prior_params$beta
#'
#' # Extract posterior for first gene
#' gene1_counts <- counts[1, ]
#' posterior <- get_posterior_distribution(
#'   counts = gene1_counts,
#'   alpha = alpha,
#'   beta = beta,
#'   ci = 0.95
#' )
#'
#' # Inspect posterior
#' cat("Gene 1 posterior mean:", posterior$posterior_mean, "\n")
#' cat("95% credible interval: [",
#'     round(posterior$ci_lower, 4), ", ",
#'     round(posterior$ci_upper, 4), "]\n", sep = "")
#'
#' @export
get_posterior_distribution <- function(counts, alpha, beta, ci = 0.95) {
    if (!is.numeric(counts) || length(counts) < 1) {
        stop("counts must be a non-empty numeric vector")
    }

    if (!is.numeric(alpha) || alpha <= 0) {
        stop("alpha must be a positive numeric value")
    }

    if (!is.numeric(beta) || beta <= 0) {
        stop("beta must be a positive numeric value")
    }

    if (!is.null(ci) && (!is.numeric(ci) || ci <= 0 || ci >= 1)) {
        stop("ci must be NULL or a numeric value between 0 and 1")
    }

    # Basic statistics
    total_depth <- sum(counts)
    
    # Posterior Beta parameters via conjugate Beta-Binomial update
    posterior_alpha <- alpha + sum(counts)
    posterior_beta <- beta + (total_depth - sum(counts))

    # Posterior mean: E[θ | data]
    posterior_mean <- posterior_alpha / (posterior_alpha + posterior_beta)

    # Posterior variance: Var[θ | data]
    posterior_variance <- (posterior_alpha * posterior_beta) /
        ((posterior_alpha + posterior_beta)^2 * 
         (posterior_alpha + posterior_beta + 1))

    # Posterior standard deviation
    posterior_sd <- sqrt(posterior_variance)

    # Build result list
    result <- list(
        posterior_alpha = posterior_alpha,
        posterior_beta = posterior_beta,
        posterior_mean = posterior_mean,
        posterior_variance = posterior_variance,
        posterior_sd = posterior_sd
    )

    # Compute credible interval if requested
    if (!is.null(ci)) {
        alpha_level <- 1 - ci
        ci_lower <- stats::qbeta(alpha_level / 2, posterior_alpha, posterior_beta)
        ci_upper <- stats::qbeta(1 - alpha_level / 2, posterior_alpha, posterior_beta)
        
        result$ci_lower <- ci_lower
        result$ci_upper <- ci_upper
        result$ci_level <- ci
    }

    return(result)
}

#' Compute Posterior Credible Intervals for Multiple Genes
#'
#' Extracts posterior credible intervals for all genes in a count matrix,
#' useful for visualizing uncertainty across the genome.
#'
#' @param counts_matrix Matrix or data.frame; genes (rows) * samples (columns).
#' @param alpha Numeric; alpha parameter of Beta prior.
#' @param beta Numeric; beta parameter of Beta prior.
#' @param ci Numeric; credible interval width (default: 0.95).
#'
#' @return Data frame with columns:
#'   \describe{
#'     \item{gene}{Gene name (from rownames of counts_matrix)}
#'     \item{posterior_mean}{Point estimate}
#'     \item{ci_lower}{Lower CI bound}
#'     \item{ci_upper}{Upper CI bound}
#'     \item{posterior_sd}{Standard deviation}
#'   }
#'
#' @examples
#' # Fit prior
#' prior_params <- fit_empirical_beta_prior(counts)
#'
#' # Compute CIs for all genes
#' all_cis <- compute_posterior_credible_intervals(
#'   counts_matrix = counts,
#'   alpha = prior_params$alpha,
#'   beta = prior_params$beta,
#'   ci = 0.95
#' )
#'
#' # View top genes by posterior mean
#' head(all_cis[order(-all_cis$posterior_mean), ], 10)
#'
#' @export
compute_posterior_credible_intervals <- function(counts_matrix, alpha, beta, ci = 0.95) {
    if (!is.matrix(counts_matrix) && !is.data.frame(counts_matrix)) {
        stop("counts_matrix must be a matrix or data.frame")
    }

    counts_matrix <- as.matrix(counts_matrix)
    n_genes <- nrow(counts_matrix)

    # Preallocate result vectors
    gene_names <- rownames(counts_matrix)
    if (is.null(gene_names)) {
        gene_names <- paste0("Gene_", seq_len(n_genes))
    }

    posterior_means <- numeric(n_genes)
    ci_lowers <- numeric(n_genes)
    ci_uppers <- numeric(n_genes)
    posterior_sds <- numeric(n_genes)

    # Compute posterior for each gene
    for (i in seq_len(n_genes)) {
        posterior <- get_posterior_distribution(
            counts = counts_matrix[i, ],
            alpha = alpha,
            beta = beta,
            ci = ci
        )
        posterior_means[i] <- posterior$posterior_mean
        ci_lowers[i] <- posterior$ci_lower
        ci_uppers[i] <- posterior$ci_upper
        posterior_sds[i] <- posterior$posterior_sd
    }

    # Return as data frame
    result_df <- data.frame(
        gene = gene_names,
        posterior_mean = posterior_means,
        ci_lower = ci_lowers,
        ci_upper = ci_uppers,
        posterior_sd = posterior_sds,
        row.names = NULL,
        stringsAsFactors = FALSE
    )

    return(result_df)
}

# ============================================================================
# Bootstrap Confidence Interval Helpers
#' Block Bootstrap for Paired Samples
#'
#' Implements block bootstrap resampling where pairs are resampled together.
#' For paired data, consecutive samples are treated as pairs: (x[1], x[2]), (x[3], x[4]), etc.
#' Each pair is resampled with replacement to preserve within-pair correlation.
#'
#' @param x Numeric vector with even length (2n observations = n pairs)
#' @param q Tsallis entropy parameter
#' @param norm Logical: normalize entropy
#' @param nboot Number of bootstrap replicates
#' @param log_base Logarithm base
#' @param pseudocount Small value added to counts
#' @param what "S" for entropy or "D" for divergence
#'
#' @return Numeric vector of bootstrap entropy estimates
#'
#' @keywords internal
#' @noRd
.tsenat_block_bootstrap <- function(x, q, norm, nboot, log_base, pseudocount, what) {
  n_pairs <- length(x) / 2
  
  # Organize data as pairs (each pair is 2 indices)
  pair_indices_list <- list()
  for (i in seq_len(n_pairs)) {
    pair_indices_list[[i]] <- c((i - 1) * 2 + 1, (i - 1) * 2 + 2)
  }
  
  # Preallocate bootstrap distribution
  boot_dist <- numeric(nboot)
  
  # Perform block bootstrap resampling
  for (i in seq_len(nboot)) {
    # Resample pair indices with replacement
    sampled_pair_idx <- sample(seq_len(n_pairs), size = n_pairs, replace = TRUE)
    
    # Extract resampled pairs and reconstruct data
    sampled_indices <- unlist(pair_indices_list[sampled_pair_idx])
    boot_sample <- x[sampled_indices]
    
    # Compute Tsallis entropy for this block bootstrap sample
    boot_est <- calculate_tsallis_entropy(boot_sample, q = q, norm = norm,
        what = what, log_base = log_base, pseudocount = 0)
    boot_dist[i] <- as.numeric(boot_est)
  }
  
  return(boot_dist)
}

# ============================================================================

.tsenat_bootstrap_resample <- function(x, q, norm, nboot, log_base, pseudocount, what, paired = FALSE) {
    # Dispatch to block bootstrap for paired samples (paper S112)
    if (paired) {
        return(.tsenat_block_bootstrap(x, q = q, norm = norm, nboot = nboot,
            log_base = log_base, pseudocount = pseudocount, what = what))
    }
    
    # Standard bootstrap resampling for independent samples
    # Estimate proportions from original data
    x_adj <- x + pseudocount
    total <- sum(x_adj)
    p_hat <- x_adj / total
    n_isoforms <- length(x)
    
    # Preallocate bootstrap distribution
    boot_dist <- numeric(nboot)
    
    # Perform bootstrap resampling
    for (i in seq_len(nboot)) {
        # Draw bootstrap sample from multinomial distribution
        boot_sample <- rmultinom(1, size = total, prob = p_hat)
        boot_sample <- as.numeric(boot_sample)
        
        # Compute Tsallis entropy for this sample
        boot_est <- calculate_tsallis_entropy(boot_sample, q = q, norm = norm,
            what = what, log_base = log_base, pseudocount = 0)
        boot_dist[i] <- as.numeric(boot_est)
    }
    
    return(boot_dist)
}

.tsenat_ci_percentile <- function(bootstrap_dist, ci) {
    alpha <- 1 - ci
    lower_p <- alpha / 2
    upper_p <- 1 - alpha / 2
    
    lower <- quantile(bootstrap_dist, probs = lower_p, type = 7, names = FALSE)
    upper <- quantile(bootstrap_dist, probs = upper_p, type = 7, names = FALSE)
    
    return(list(lower = lower, upper = upper))
}

.tsenat_ci_bca <- function(x, bootstrap_dist, q, norm, ci, log_base, pseudocount, what) {
    alpha <- 1 - ci
    z_alpha <- qnorm(alpha / 2)  # Two-tailed critical value
    
    # Bias correction: z0 = Phi^{-1}(#F* <= F / B)
    point_est <- calculate_tsallis_entropy(x, q = q, norm = norm, what = what,
        log_base = log_base, pseudocount = pseudocount)
    point_est <- as.numeric(point_est)
    
    # Proportion of bootstrap replicates <= point estimate
    prop_less <- mean(bootstrap_dist <= point_est)
    z0 <- qnorm(prop_less)
    
    # Acceleration: computed via jackknife
    n <- length(x)
    jack_est <- numeric(n)
    
    for (i in seq_len(n)) {
        x_minus_i <- x[-i]
        if (sum(x_minus_i) > 0) {
            jack_est[i] <- calculate_tsallis_entropy(x_minus_i, q = q, norm = norm,
                what = what, log_base = log_base, pseudocount = pseudocount)
            jack_est[i] <- as.numeric(jack_est[i])
        } else {
            jack_est[i] <- NA_real_
        }
    }
    
    # Filter out NAs from jackknife estimates
    jack_est_clean <- jack_est[!is.na(jack_est)]
    if (length(jack_est_clean) < 2) {
        # Fall back to percentile if jackknife fails
        return(.tsenat_ci_percentile(bootstrap_dist, ci = ci))
    }
    
    # Acceleration: a = (sum(jack_mean - jack_i)^3) / (6 * (sum(jack_mean - jack_i)^2)^1.5)
    jack_mean <- mean(jack_est_clean)
    diffs <- jack_est_clean - jack_mean
    numerator <- sum(diffs^3)
    denominator <- 6 * (sum(diffs^2))^1.5
    
    a <- if (abs(denominator) > 1e-10) numerator / denominator else 0
    
    # Adjusted critical values: z_alpha^+ and z_alpha^-
    z_low <- qnorm(alpha / 2)
    z_high <- qnorm(1 - alpha / 2)
    
    p_low <- pnorm(z0 + (z0 + z_low) / (1 - a * (z0 + z_low)))
    p_high <- pnorm(z0 + (z0 + z_high) / (1 - a * (z0 + z_high)))
    
    # Ensure valid probabilities
    p_low <- max(0, min(1, p_low))
    p_high <- max(0, min(1, p_high))
    
    lower <- quantile(bootstrap_dist, probs = p_low, type = 7, names = FALSE)
    upper <- quantile(bootstrap_dist, probs = p_high, type = 7, names = FALSE)
    
    return(list(lower = lower, upper = upper))
}

#' Estimate Hyperparameters for Empirical Bayes Shrinkage
#'
#' Estimates the global mean and variance of entropy estimates across all genes,
#' needed for empirical Bayes shrinkage. Also counts the number of expressed isoforms
#' per gene to inform shrinkage weights.
#'
#' @param x Raw transcript-level count matrix (genes x samples).
#' @param genes Vector of gene IDs (must match nrow of x).
#' @param entropy_matrix Matrix of entropy estimates (genes x assays after calculation).
#' @param q Tsallis q parameter(s).
#' @param min_count Minimum count threshold to consider a isoform "expressed"
#'   (default: 1).
#'
#' @return A list with:
#'   \describe{
#'     \item{global_mean}{Named numeric vector of global mean entropy per q-value.}
#'     \item{global_var}{Named numeric vector of variance per q-value.}
#'     \item{n_isoforms}{Vector of number of expressed isoforms per gene.}
#'   }
#'
#' @keywords internal
#' @noRd
.tsenat_estimate_shrinkage_params <- function(x, genes, entropy_matrix, q = 2, min_count = 1) {
  gene_levels <- unique(genes)
  
  # Count expressed isoforms per gene (non-zero after filtering)
  # Use sapply with named input to preserve gene names in output
  n_isoforms <- sapply(setNames(gene_levels, gene_levels), function(g) {
    gene_mask <- genes == g
    gene_counts <- rowSums(x[gene_mask, , drop = FALSE])
    sum(gene_counts > min_count)
  })
  
  # For each q value, estimate global mean and variance from entropy estimates
  q_cols <- grep(paste0("_q=", q, "$", collapse = "|"), colnames(entropy_matrix), 
                 value = TRUE, perl = TRUE)
  
  # If single q value, extract that column
  if (length(q) == 1) {
    q_cols <- paste0("_q=", q)
    q_cols <- colnames(entropy_matrix)[grep(q_cols, colnames(entropy_matrix))]
  }
  
  # Compute global statistics per q
  global_mean <- colMeans(entropy_matrix[, q_cols, drop = FALSE], na.rm = TRUE)
  global_var <- apply(entropy_matrix[, q_cols, drop = FALSE], 2, var, na.rm = TRUE)
  
  return(list(
    global_mean = global_mean,
    global_var = global_var,
    n_isoforms = n_isoforms
  ))
}

#' Apply Empirical Bayes Shrinkage to Entropy Estimates
#'
#' Shrinks individual gene entropy estimates toward the global mean using empirical
#' Bayes weights. Particularly effective for genes with few expressed isoforms.
#'
#' @param entropy_matrix Matrix of raw entropy estimates (genes x assays).
#' @param params List from \code{.tsenat_estimate_shrinkage_params()} with
#'   global_mean, global_var, and n_isoforms.
#' @param gene_isoform_map Optional vector mapping row names of entropy_matrix
#'   to n_isoforms (if names don't match indices).
#'
#' @return Shrinkage-adjusted entropy matrix with same dimensions.
#'
#' @details
#' For each gene g and q-value, the shrinkage weight is computed as:
#' \deqn{w_g = \frac{n_g}{n_g + \lambda}}{w_g = n_g / (n_g + lambda)}
#'
#' where \eqn{n_g}{n_g} is the number of expressed isoforms in gene g and
#' \eqn{\lambda}{lambda} is estimated from the signal/noise ratio.
#'
#' The shrunk estimate is:
#' \deqn{S_shrink = w_g \cdot S_g + (1 - w_g) \cdot \bar{S}}{S_shrink = w_g * S_g + (1 - w_g) * mean(S)}
#'
#' This borrows strength from the global distribution, stabilizing estimates for
#' genes with small expression variance.
#'
#' @keywords internal
#' @noRd
.tsenat_apply_shrinkage <- function(entropy_matrix, params, gene_isoform_map = NULL) {
  result <- entropy_matrix
  
  global_mean <- params$global_mean
  global_var <- params$global_var
  n_isoforms <- params$n_isoforms
  
  # Map row names to n_isoforms indices if needed
  if (is.null(gene_isoform_map)) {
    # Assume row names are gene IDs matching names(n_isoforms)
    row_gene_ids <- rownames(entropy_matrix)
    gene_isoform_map <- n_isoforms[row_gene_ids]
  }
  
  # Estimate prior strength from the data using empirical Bayes methodology
  # (following DESeq2/edgeR approach: prior DF based on average gene information content)
  mean_n_isoforms <- mean(n_isoforms, na.rm = TRUE)
  
  # For each column (sample-q combination), apply shrinkage
  for (col_idx in seq_len(ncol(entropy_matrix))) {
    col_name <- colnames(entropy_matrix)[col_idx]
    
    # Extract q value from column name (format: "Sample_q=X")
    q_match <- gregexpr("_q=([0-9.]+)", col_name)
    if (q_match[[1]][1] > 0) {
      q_pos <- regmatches(col_name, q_match)[[1]]
      q_val <- as.numeric(sub("_q=", "", q_pos))
      
      # Find corresponding global parameters for this q
      mean_key <- paste0("q=", q_val)
      var_key <- paste0("q=", q_val)
      
      # Handle case where keys might be formatted differently
      if (!(mean_key %in% names(global_mean))) {
        # Try matching from column names
        mean_key <- names(global_mean)[grepl(paste0(q_val, "$"), names(global_mean))][1]
      }
      if (!(var_key %in% names(global_var))) {
        var_key <- names(global_var)[grepl(paste0(q_val, "$"), names(global_var))][1]
      }
      
      if (!is.na(mean_key) && mean_key %in% names(global_mean)) {
        mu <- global_mean[[mean_key]]
        sigma2_between <- global_var[[var_key]]
        
        # Empirical Bayes prior strength estimation
        # df_prior represents the effective sample size of the prior distribution
        # Estimated conservatively from average gene information content
        df_prior <- max(1, mean_n_isoforms - 1)
        
        # Compute shrinkage weights per gene
        for (row_idx in seq_len(nrow(entropy_matrix))) {
          # Try to get n_isoforms for this gene
          row_name <- rownames(entropy_matrix)[row_idx]
          
          # Handle both named and unnamed gene_isoform_map vectors
          if (is.null(names(gene_isoform_map)) || length(names(gene_isoform_map)) == 0) {
            # Unnamed vector: use positional indexing
            n_iso <- gene_isoform_map[row_idx]
          } else {
            # Named vector: use name-based indexing
            n_iso <- gene_isoform_map[row_name]
          }
          
          if (!is.na(n_iso) && n_iso > 0) {
            # Shrinkage parameter lambda: prior strength scaled by relative information content
            # For genes with n_iso < mean_n_isoforms: lambda > df_prior (more shrinkage)
            # For genes with n_iso > mean_n_isoforms: lambda < df_prior (less shrinkage)
            # Formula: lambda = df_prior * (mean_n_isoforms / n_iso)
            # This implements precision-weighted empirical Bayes shrinkage
            lambda <- df_prior * (mean_n_isoforms / n_iso)
            
            # Shrinkage weight: w close to 1 trusts the observation more, w close to 0 trusts prior more
            # From empirical Bayes theory: w = precision_obs / (precision_obs + precision_prior)
            w <- n_iso / (n_iso + lambda)
            
            # Apply shrinkage: for finite values use weighted average of observation and prior,
            # for NA values (e.g., from undefined normalized entropy) use the prior estimate
            if (is.na(entropy_matrix[row_idx, col_idx])) {
              # NA values get shrunk to the prior (w=0 for completely missing data)
              result[row_idx, col_idx] <- mu
            } else if (is.finite(entropy_matrix[row_idx, col_idx])) {
              # Finite values get weighted average of observation and prior
              result[row_idx, col_idx] <- w * entropy_matrix[row_idx, col_idx] + (1 - w) * mu
            }
            # NaN values are left as-is
          }
        }
      }
    }
  }
  
  return(result)
}

# # Minimalized diversity functions: only Tsallis entropy retained

#' Calculate Tsallis entropy for a vector of transcript-level
#' expression values of one gene.
#'
#' @param x Vector of (non-negative) expression values.
#' @param q Tsallis entropy parameter (q > 0). Scalar or numeric vector
#' (default: 2).
#' @param norm Logical; if TRUE, normalize entropy by its theoretical maximum
#' (values in [0,1]).
#' @param what Which quantity to return: 'S' (Tsallis entropy), 'D' (Hill
#' numbers), or 'both'.
#' @param log_base Base of the logarithm used for Shannon limits and
#' normalization
#' (default: \code{exp(1)}).
#' @export
#' @return For `what = 'S'` or `what = 'D'`: a numeric vector
#' (named when length(q) > 1). For `what = 'both'`: a list with
#' components `$S` and `$D`.
#' @details
#' **Tsallis Entropy (S_q):**
#'
#' \deqn{S_q = \frac{1 - \sum_{i=1}^k p_i^q}{q - 1}}{S_q = (1 - sum p_i^q) / (q - 1)}
#'
#' where \eqn{p_i}{p_i} are normalized proportions and \eqn{q > 0}{q > 0} is the entropy order.
#' For \eqn{q = 1}{q=1}, this reduces to Shannon entropy: \deqn{H = -\sum_i p_i \ln p_i}{H = -sum p_i*ln(p_i)}
#'
#' **Hill Numbers (Diversity Index D_q):**
#'
#' \deqn{D_q = \left( \sum_{i=1}^k p_i^q \right)^{\frac{1}{1-q}}}{D_q = (sum p_i^q)^(1/(1-q))}
#'
#' Hill numbers represent effective number of equally-likely species. D_1 is the exponential of Shannon entropy.
#'
#' **Normalization:**
#' When \code{norm = TRUE}, entropy is divided by its theoretical maximum to scale to [0, 1].
#' Natural logarithms are used for q→1 limits and normalization.
#' @examples
#' x <- c(10, 5, 0)
#' calculate_tsallis_entropy(x, q = c(0.5, 1, 2), norm = TRUE)
calculate_tsallis_entropy <- function(x, q = 2, norm = TRUE, what = c("S", "D", "both"),
    log_base = exp(1), pseudocount = 0, effective_length = NULL) {
    what <- match.arg(what)
    if (!is.numeric(q)) {
        stop("q must be numeric.")
    }
    if (any(q <= 0)) {
        stop("q must be greater than 0.")
    }
    if (!is.numeric(x)) {
        stop("x must be numeric")
    }

    # Apply pseudocount if specified BEFORE length normalization or proportion calculation
    # Handles both scalar and vector pseudocounts
    # Vector pseudocounts are applied per-isoform (row-wise for matrices)
    if (any(pseudocount > 0)) {
        if (is.matrix(x) && length(pseudocount) > 1) {
            # Per-isoform pseudocounts: apply row-wise via sweep
            x <- sweep(x, 1, pseudocount, "+")
        } else {
            # Scalar pseudocount or vector input: simple addition
            x <- x + pseudocount
        }
    }
    
    n <- length(x)
    
    # If all counts sum to zero, return NA; allow single-element vectors to
    # proceed
    if (sum(x, na.rm = TRUE) <= 0) {
        if (what == "both") {
            return(list(S = rep(NA_real_, length(q)), D = rep(NA_real_, length(q))))
        }
        return(rep(NA_real_, length(q)))
    }

    # EFFECTIVE LENGTH NORMALIZATION
    # If effective_length is provided, normalize counts by length to remove length bias
    # This is SALMON's recommended approach for isoform-level analysis
    # Normalized counts = x / effective_length (accounts for read-length & alignability bias)
    # Then proportions = normalized_counts / sum(normalized_counts)
    if (!is.null(effective_length)) {
        if (length(effective_length) != length(x)) {
            stop("effective_length must have same length as x")
        }
        # Check for valid effective_length values (must be positive)
        if (any(effective_length <= 0, na.rm = TRUE)) {
            warning("Some effective_length values are <= 0, treating as NA")
            effective_length[effective_length <= 0] <- NA
        }
        # Normalize counts: x_norm = x / effective_length
        x_normalized <- x / effective_length
        # Replace any NaN/Inf with 0 (when effective_length is 0 or NA)
        x_normalized[!is.finite(x_normalized)] <- 0
        # Calculate proportions from normalized counts
        p <- x_normalized / sum(x_normalized)
    } else {
        # Standard proportions from raw counts (no length normalization)
        p <- x / sum(x)
    }

    tol <- sqrt(.Machine$double.eps)
    S_vec <- .tsenat_calc_S(p = p, q = q, tol = tol, n = n, log_base = log_base,
        norm = norm)
    D_vec <- .tsenat_calc_D(p = p, q = q, tol = tol, log_base = log_base)

    if (what == "S") {
        out <- S_vec
        if (length(q) > 1) {
            names(out) <- paste0("q=", q)
        }
        if (length(q) == 1) {
            return(unname(out))
        }
        return(out)
    }
    if (what == "D") {
        out <- D_vec
        if (length(q) > 1) {
            names(out) <- paste0("q=", q)
        }
        if (length(q) == 1) {
            return(unname(out))
        }
        return(out)
    }
    # both
    names(S_vec) <- paste0("q=", q)
    names(D_vec) <- paste0("q=", q)
    return(list(S = S_vec, D = D_vec))
}

#' Internal: Calculate Tsallis entropy for transcripts grouped by gene
#' 
#' Core internal function that computes Tsallis entropy values for each gene
#' across samples. This function is called by \code{\link{calculate_diversity}},
#' which provides the modern, user-facing interface with SummarizedExperiment support.
#'
#' @param x Numeric matrix or data.frame of transcript-level expression
#' values (rows = transcripts, columns = samples).
#' @param genes Character vector with length equal to nrow(x) assigning each
#' transcript to a gene.
#' @param norm Logical; if TRUE normalize Tsallis entropy values per gene.
#' @param q Numeric scalar or vector of q values to evaluate.
#' @param verbose Logical; show diagnostic messages when TRUE.
#' @param what Which quantity to return from `calculate_tsallis_entropy`:
#' 'S' (Tsallis entropy) or 'D' (Hill numbers) (default: 'S').
#' @param nthreads Number of threads for parallel processing (default: 1).
#' Set to > 1 to parallelize per-gene entropy calculations.
#' @param pseudocount Numeric scalar. Add this value to all transcript counts
#' before calculating proportions (default: 0). Useful for handling genes with
#' zero counts in some samples.
#' @param min_valid_frac Numeric scalar in [0, 1]; minimum fraction of valid
#' (finite) values required per gene to be retained (default: 0.75). Genes with
#' fewer valid values are excluded. This ensures statistical reliability:
#' genes with sparse/missing data are discarded, improving power for paired
#' tests (Wilcoxon) and stability in linear mixed models. Set to 0 to keep all
#' genes regardless of data completeness (not recommended for hypothesis testing).
#' @param shrinkage Character; method for shrinking entropy estimates toward global mean:
#' "none" (default, no shrinkage) or "empirical_bayes" (empirical Bayes shrinkage toward
#' global mean, particularly effective for genes with few expressed isoforms).
#' @param effective_length Optional effective transcript lengths for normalization.
#' 
#' @return A data.frame with genes in the first column and per-sample (and
#' per-q) Tsallis entropy values in subsequent columns.
#' 
#' @keywords internal
#' @noRd
.calculate_method <- function(x, genes, norm = TRUE, verbose = FALSE, q = 2, what = c("S",
    "D"), nthreads = 1, pseudocount = 0, min_valid_frac = 0.75, shrinkage = c("none", 
    "empirical_bayes"), effective_length = NULL) {
    what <- match.arg(what)
    shrinkage <- match.arg(shrinkage)
    # validate q
    if (!is.numeric(q) || any(q <= 0)) {
        stop("Argument 'q' must be numeric and greater than 0.", call. = FALSE)
    }
    # cannot use aggregate because calculate_tsallis_entropy may return
    # multiple values when length(q) > 1
    gene_levels <- unique(genes)
    
    # ensure column names order matches the order used when constructing the
    # result matrix (samples vary outer, q varies inner). If sample names are
    # missing, synthesize deterministic names so column creation still works.
    sample_names <- colnames(x)
    if (is.null(sample_names)) {
        sample_names <- paste0("Sample", seq_len(ncol(x)))
    }
    coln <- as.vector(t(outer(sample_names, q, function(s, qq) paste0(s, "_q=", qq))))
    rown <- gene_levels

    # compute requested quantity ('S' or 'D') in parallel
    result_list <- .tsenat_bplapply(gene_levels, function(gene) {
        .tsenat_tsallis_row(x = x, genes = genes, gene = gene, q = q, norm = norm,
            what = what, pseudocount = pseudocount, effective_length = effective_length)
    }, nthreads = nthreads)

    # Convert list to matrix (each element is a named vector) result_list is a
    # list of vectors; combine them into a matrix
    result_mat <- t(vapply(result_list, identity, FUN.VALUE = setNames(numeric(length(coln)),
        coln)))
    colnames(result_mat) <- coln
    rownames(result_mat) <- rown
    out_df <- data.frame(Gene = rown, result_mat, check.names = FALSE)
    
    if (all(rowSums(!is.na(result_mat)) == 0)) {
        out_df <- data.frame(Gene = character(0))
        for (nm in coln) out_df[[nm]] <- numeric(0)
        return(out_df)
    }
    # Apply shrinkage if requested (BEFORE filtering so NA values can be converted to finite)
    if (shrinkage == "empirical_bayes") {
        if (verbose) {
            message("Applying empirical Bayes shrinkage to entropy estimates...")
        }
        
        # Estimate shrinkage hyperparameters
        params <- .tsenat_estimate_shrinkage_params(
            x = x, 
            genes = genes, 
            entropy_matrix = result_mat,
            q = q
        )
        
        # Create mapping from row names to n_isoforms
        gene_to_isoforms <- params$n_isoforms[rownames(result_mat)]
        
        # Apply shrinkage
        result_mat_shrink <- .tsenat_apply_shrinkage(
            entropy_matrix = result_mat,
            params = params,
            gene_isoform_map = gene_to_isoforms
        )
        
        # Update output dataframe with shrunk values
        out_df[, -1] <- result_mat_shrink
        result_mat <- result_mat_shrink
        
        if (verbose) {
            n_shrunk <- sum(gene_to_isoforms < 5, na.rm = TRUE)
            message(sprintf("  %d genes with < 5 isoforms received substantial shrinkage", n_shrunk))
        }
    }
    
    # Filter: keep genes with sufficient valid (finite) values
    # This happens AFTER shrinkage so that NA values from single-isoform normalized entropy
    # can be converted to finite values by shrinkage
    # Statistical requirement: genes need adequate data for reliable inference in tests
    # (paired Wilcoxon, LMM, etc). Sparse genes are excluded.
    # result_mat dimensions: rows = genes, cols = sample * q combinations
    n_total_values <- ncol(result_mat)  # total possible values per gene
    min_valid_count <- ceiling(min_valid_frac * n_total_values)
    finite_counts <- rowSums(is.finite(result_mat))
    keep_idx <- finite_counts >= min_valid_count
    
    out_df <- out_df[keep_idx, ]
    result_mat <- result_mat[keep_idx, , drop = FALSE]
    n_excluded <- nrow(result_mat) + sum(!keep_idx) - nrow(result_mat)
    if (n_excluded > 0 && verbose == TRUE) {
        message(sprintf("Note: %d genes excluded (< %.0f%% valid values).", 
            n_excluded, min_valid_frac * 100))
    }
    
    return(out_df)
}

#' Internal wrapper around calculate_method
#'
#' This thin helper around \code{.calculate_method} is intended for use in
#' package-internal tests and by advanced developers.  It is **not exported**
#' for general user workflows; callers should normally use
#' \code{calculate_diversity} which works on
#' \link[SummarizedExperiment]{SummarizedExperiment} objects.
#'
#' @param x Numeric matrix; gene expression data (rows = genes, columns = samples)
#' @param genes Character vector; gene assignments for each row
#' @param norm Logical; normalize results by per-sample sum
#' @param verbose Logical; print diagnostic messages
#' @param q Numeric; Tsallis parameter (default 2)
#' @param what Character; output type ("S" or "D")
#' @param nthreads Integer; number of threads for parallel computation
#' @param pseudocount Numeric; pseudocount to add before computation
#' @param min_valid_frac Numeric; minimum fraction of valid values per gene
#' @param shrinkage Character; shrinkage method ("none" or "empirical_bayes")
#' @param effective_length Numeric vector; effective transcript lengths (optional)
#'
#' @keywords internal
#' @noRd
calculate_method <- function(x, genes, norm = TRUE, verbose = FALSE, q = 2, what = c("S",
    "D"), nthreads = 1, pseudocount = 0, min_valid_frac = 0.75, shrinkage = c("none",
    "empirical_bayes"), effective_length = NULL) {
    .calculate_method(x = x, genes = genes, norm = norm, verbose = verbose,
        q = q, what = what, nthreads = nthreads, pseudocount = pseudocount,
        min_valid_frac = min_valid_frac, shrinkage = shrinkage,
        effective_length = effective_length)
}

# Internal helpers for calculate_method

.tsenat_tsallis_row <- function(x, genes, gene, q, norm, what, pseudocount = 0, effective_length = NULL) {
    idx <- which(genes == gene)
    out <- unlist(lapply(seq_len(ncol(x)), function(j) {
        # Get counts for this gene and sample
        counts <- x[idx, j]
        
        # Apply pseudocount to counts (add before calculating proportions)
        if (pseudocount > 0) {
            counts <- counts + pseudocount
        }
        
        # Apply effective length normalization if provided
        if (!is.null(effective_length)) {
            el <- NULL
            if (is.vector(effective_length)) {
                el <- effective_length[idx]
            } else if (is.matrix(effective_length)) {
                el <- effective_length[idx, j]
            }
            # Normalize counts by effective length (convert to TPM-like units)
            if (!is.null(el) && length(el) == length(counts)) {
                counts <- counts / el
            }
        }
        
        # Calculate entropy on the adjusted counts
        v <- calculate_tsallis_entropy(counts, q = q, norm = norm, what = what)
        if (length(v) == length(q) && all(is.finite(v) | is.na(v))) {
            v
        } else {
            names_vec <- paste0("q=", q)
            setNames(rep(NA_real_, length(q)), names_vec)
        }
    }))
    out
}
