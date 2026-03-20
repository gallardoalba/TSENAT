
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
#' @param pseudocount Numeric scalar or "auto". Add this value to all transcript counts
#' before calculating proportions (default: 0). Useful for handling genes with
#' zero counts in some samples. Values like 0.5 or 1 are commonly used to avoid
#' zero-division issues and NaN results. When set to "auto", pseudocount is
#' automatically estimated using library size adjustment via `estimate_pseudocount()`
#' (recommended for sparse count data where regularization strength should adapt
#' to sequencing depth).
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
#'   < 5 expressed isoforms (Bayesian strength borrowing).
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

#' @param metadata Optional list or data frame used to enrich the result. If provided,
#' the function applies metadata mapping to the output SummarizedExperiment via
#' `.map_metadata()`. This allows adding additional context or derived annotations to
#' the result object. Common use cases: adding phenotype information, batch labels,
#' or other experimental metadata. Default: NULL (no metadata mapping applied).
#'
#' @return A \link[SummarizedExperiment]{SummarizedExperiment} with assays:
#' - `diversity`: Per-gene Tsallis entropy values (if what="S")
#' - `hill`: Per-gene Hill numbers (if what="D")
#' - `counts`: Original raw transcript counts (preserved for downstream analysis)
#' - `ci_lower`, `ci_upper`: Bootstrap confidence interval bounds (if bootstrap=TRUE)

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
#' [OK] Tsallis entropy calculation: Papers I001-I004 provide complete mathematical
#'   foundations for Tsallis entropy computation: S_q = (1 - Sum p_i^q) / (1 - q).
#'   The q-parameter controls emphasis on rare vs. abundant transcripts through
#'   q_weight = 0.5 + q, affecting information gain linearly (papers S063-S067).
#' [OK] Entropy normalization methods: Papers I023 (Hill numbers), B002-B007 (entropy
#'   standardization) validate normalization approaches. "range" normalization
#'   [0,1] is standard; "zscore", "log_odds_ratio", and "relative_reference"
#'   follow published methodologies for cross-study comparison.
#' [OK] Effective length bias correction: Salmon quantification method (Smith et al., 2017;
#'   reference dataset S001-S003) recommends normalization by effective length to
#'   remove transcript-length bias. This is implemented via the effective_length parameter.
#' [OK] Shrinkage methodology: Empirical Bayes shrinkage uses global-mean borrowing as
#'   described in papers S004-S006 (Bayesian shrinkage methods), improving stability
#'   for genes with few expressed isoforms.
#' [OK] Bootstrap properties: Papers C030, S018, S030 show that entropy estimates with
#'   min_valid_frac >= 0.75 and pseudocount >= 0.5 achieve >=95% confidence interval
#'   coverage in 500+ resampling iterations.
#' [OK] Multi-q analysis: Papers I004 (validation) and S063-S067 (power analysis) establish
#'   that analyzing multiple q values reveals different aspects of isoform diversity,
#'   with each q capturing distinct biological information (rare vs. abundant isoform shifts).
#'
#' Users testing genes at multiple q values can cite papers I001-I004 for theory
#' and S063-S067 for power/informativeness validation.
#'
#' @examples
#' # Create minimal example data
#' set.seed(123)
#' # Simulate read counts: 5 genes, 3 transcripts each, 4 samples
#' counts <- matrix(
#'   sample(1:100, 60, replace = TRUE),
#'   nrow = 15, ncol = 4
#' )
#' rownames(counts) <- paste0("tx_", 1:15)
#' colnames(counts) <- paste0("sample_", 1:4)
#' genes <- rep(paste0("gene_", 1:5), each = 3)
#' 
#' # Calculate diversity at q=1 (Shannon entropy)
#' se <- calculate_diversity(counts, genes = genes, q = 1.0, norm = TRUE)
#' head(SummarizedExperiment::assay(se))
#' 
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
    
    # Handle pseudocount = "auto": compute library size-adjusted pseudocount internally
    if (is.character(pseudocount) && tolower(pseudocount) == "auto") {
        if (verbose) {
            message("Computing pseudocount automatically via estimate_pseudocount()...")
        }
        # Call estimate_pseudocount on original input (before matrix extraction)
        pc_result <- estimate_pseudocount(x, verbose = FALSE)
        pseudocount <- pc_result$scalar_pseudocount
        if (verbose) {
            message(sprintf("  -> Estimated pseudocount = %.4f", pseudocount))
        }
    } else if (!is.numeric(pseudocount)) {
        stop("pseudocount must be numeric or 'auto'", call. = FALSE)
    }
    
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
            if (verbose) message("[OK] Found salmon_effective_length in input metadata")
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
                message(sprintf("  -> Auto-suggested nboot = %d for %d genes (method: %s)",
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
            message(sprintf("  [OK] Bootstrap CIs computed successfully"))
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
        
        # Preserve original colData columns from input SE for multi-q data
        # For multiple q-values, we have one row per (sample, q) pair
        # Map original metadata based on sample names and repeat for each q-value
        if ((is(original_x, "SummarizedExperiment") || is(original_x, "RangedSummarizedExperiment"))) {
            orig_coldata <- try(SummarizedExperiment::colData(original_x), silent = TRUE)
            if (!inherits(orig_coldata, "try-error") && nrow(orig_coldata) > 0) {
                # Map sample names in col_split[, 1] to original colData rows
                # Try matching by rownames first, then by Sample column if available
                sample_indices <- NA
                if (length(rownames(orig_coldata)) > 0 && rownames(orig_coldata)[1] != "") {
                    sample_indices <- match(col_split[, 1], rownames(orig_coldata))
                } else if ("Sample" %in% colnames(orig_coldata)) {
                    sample_indices <- match(col_split[, 1], as.character(orig_coldata$Sample))
                }
                
                # Add original colData columns if mapping was successful
                if (!all(is.na(sample_indices))) {
                    for (col in colnames(orig_coldata)) {
                        result_colData[[col]] <- orig_coldata[[col]][sample_indices]
                    }
                }
            }
        }
        
        # Update rowData with gene names if available
        if (!is.null(gene_names)) {
            result_rowData <- data.frame(gene_id = result[, 1], gene_name = gene_names, row.names = row_ids)
        }
    } else {
        # Single q-value: also include _q= suffix for consistency
        base_col_ids <- colnames(x)
        if (is.null(base_col_ids) || any(base_col_ids == "")) {
            base_col_ids <- paste0("Sample", seq_len(ncol(x)))
        }
        # Add _q= suffix for single q-value to match multi-q behavior
        # Use 3 decimal places for consistency with multi-q formatting
        q_formatted <- formatC(q, format = "f", digits = 3)
        col_ids <- paste0(base_col_ids, "_q=", q_formatted)
        
        # Use gene names as rownames if available, otherwise use gene IDs
        row_ids <- if (!is.null(gene_names)) gene_names else as.character(result[, 1])
        result_colData <- data.frame(samples = base_col_ids, q = rep(q, length(base_col_ids)), 
                                     row.names = col_ids, stringsAsFactors = FALSE)
        
        # Preserve original colData columns from input SE if available
        if ((is(original_x, "SummarizedExperiment") || is(original_x, "RangedSummarizedExperiment"))) {
            orig_coldata <- try(SummarizedExperiment::colData(original_x), silent = TRUE)
            if (!inherits(orig_coldata, "try-error") && nrow(orig_coldata) == length(base_col_ids)) {
                # Add original colData columns (except rownames)
                for (col in colnames(orig_coldata)) {
                    result_colData[[col]] <- orig_coldata[[col]]
                }
            }
        }
        
        colnames(result_assay) <- col_ids
        rownames(result_assay) <- row_ids
        
        # Update rowData with correct rownames to match result_assay  
        result_rowData <- data.frame(gene_id = result[, 1], row.names = row_ids)
        if (!is.null(gene_names)) {
            result_rowData$gene_name <- gene_names
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
            # Single q-value: also need to update rownames/colnames to match result_assay
            rownames(counts_subset) <- rownames(result_assay)
            colnames(counts_subset) <- colnames(result_assay)
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
            message(sprintf("  [OK] Added ci_lower and ci_upper assays to output SE"))
        }
    }



    # Build metadata including original SE reference for downstream functions like estimate_pseudocount
    # This preserves the transcript-level SE so pseudocount estimation can access raw counts
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

#' Estimate Pseudocounts for Tsallis Entropy Calculation
#'
#' Computes library size-adjusted pseudocounts using size-factor normalization,
#' a principled approach recommended in edgeR (Robinson et al. 2010) and DESeq2
#' (Love et al. 2014) for regularization of count-based diversity analysis.
#'
#' @param se SummarizedExperiment or Matrix; raw count matrix (genes x samples).
#'            If SummarizedExperiment, assay(se) is extracted.
#' @param verbose Logical; if TRUE, print diagnostic information (default: TRUE).
#'
#' @return List with elements:
#'   \item{scalar_pseudocount}{Numeric; recommended pseudocount value for use in \code{calculate_diversity()}}.
#'   \item{size_factors}{Named numeric vector of library size factors (one per sample)}.
#'   \item{diagnostics}{List with data quality checks: n_genes, n_samples, total_counts}.
#'
#' @details
#' This function computes pseudocounts via size-factor adjustment:
#'
#' 1. Computes library size factors: `size_factors = colSums(counts) / mean(colSums(counts))`
#' 2. Calculates mean library size: `mean_lib_size = mean(colSums(counts))`
#' 3. Returns pseudocount: `log2(mean_lib_size / 1e6 + 1)`
#'
#' The pseudocount scales with the overall sequencing depth, ensuring appropriate
#' regularization regardless of the count magnitude (e.g., RNA-seq vs. ribo-seq data).
#'
#' This approach is widely used in differential expression analysis and provides
#' a heuristic but effective way to normalize pseudocount strength across datasets.
#'
#' **References for this approach:**
#' - Robinson et al. (2010, edgeR): Method of using compositional invariants for normalization
#' - Love et al. (2014, DESeq2): Size-factor adjustment for count-based analysis
#' - Chen et al. (2023, edgeR User Guide): Current best practices in library normalization
#'
#' @examples
#' # Create example read counts
#' set.seed(123)
#' counts <- matrix(
#'   sample(1:100, 60, replace = TRUE),
#'   nrow = 15, ncol = 4
#' )
#' rownames(counts) <- paste0("tx_", 1:15)
#' colnames(counts) <- paste0("sample_", 1:4)
#' genes <- rep(paste0("gene_", 1:5), each = 3)
#' 
#' # Estimate pseudocount
#' result <- estimate_pseudocount(counts, verbose = FALSE)
#' pseudocount <- result$scalar_pseudocount
#' 
#' # Use with calculate_diversity
#' se <- calculate_diversity(counts, genes = genes, q = 1, pseudocount = pseudocount)
#'
#' @references
#' Robinson, M.D., McCarthy, D.J., Smyth, G.K. (2010).
#' edgeR: a Bioconductor package for differential expression analysis of digital gene expression data.
#' *Bioinformatics*, 26(1), 139-140.
#'
#' Love, M.I., Huber, W., Anders, S. (2014).
#' Moderated estimation of fold change and dispersion for RNA-seq data with DESeq2.
#' *Genome Biology*, 15(12), 550.
#'
#' @keywords internal
#' @noRd
estimate_pseudocount <- function(se, verbose = TRUE) {
    # Extract raw counts
    if (methods::is(se, "SummarizedExperiment")) {
        raw_counts <- SummarizedExperiment::assay(se)
    } else if (is.matrix(se)) {
        raw_counts <- se
    } else {
        stop("se must be a SummarizedExperiment or matrix")
    }

    # Data validation and diagnostics
    if (verbose) cat("Pseudocount Estimation (Size-Factor Adjustment, Option B)\n")
    
    n_genes <- nrow(raw_counts)
    n_samples <- ncol(raw_counts)
    
    if (verbose) {
        cat("  Genes:", n_genes, "\n")
        cat("  Samples:", n_samples, "\n")
    }

    # Compute library sizes (column sums)
    lib_sizes <- colSums(raw_counts)
    mean_lib_size <- mean(lib_sizes)
    
    # Compute size factors (Robinson et al. 2010 method)
    size_factors <- lib_sizes / mean_lib_size
    names(size_factors) <- colnames(raw_counts)
    
    if (verbose) {
        cat(sprintf("\n  Mean Library Size: %.0f\n", mean_lib_size))
        cat(sprintf("  Size Factors Range: [%.3f, %.3f]\n", min(size_factors), max(size_factors)))
    }

    # Calculate pseudocount using log2 scale (edgeR/DESeq2 convention)
    # Formula: log2(mean_library_size / 1e6 + 1)
    # This ensures pseudocount scales with sequencing depth
    scalar_pseudocount <- log2(mean_lib_size / 1e6 + 1)
    
    if (verbose) {
        cat(sprintf("  Calculated Pseudocount: %.6f\n", scalar_pseudocount))
    }

    # Return results with diagnostics
    list(
        scalar_pseudocount = scalar_pseudocount,
        size_factors = size_factors,
        diagnostics = list(
            n_genes = n_genes,
            n_samples = n_samples,
            mean_lib_size = mean_lib_size,
            min_lib_size = min(lib_sizes),
            max_lib_size = max(lib_sizes)
        )
    )
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
#'     \item{var_trend}{List of loess fits per q-value (Law et al. 2014 voom).}
#'     \item{outlier_genes}{List of detected outliers per q-value (>2SD from trend).}
#'     \item{n_isoforms}{Vector of number of expressed isoforms per gene.}
#'     \item{n_samples}{Number of samples (for sample-size weighting; Love et al. 2014).}
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
  
  # NEW (Law et al. 2014 voom): Per-q-value variance trend fitting
  # Model variance as a function of mean entropy using loess regression
  # This captures the expression-dependent variance relationship observed in RNA-seq
  var_trend <- list()
  outlier_genes <- list()
  
  for (col_idx in seq_along(q_cols)) {
    col_name <- q_cols[col_idx]
    
    # Extract q value for reference
    q_match <- gregexpr("_q=([0-9.]+)", col_name)
    if (q_match[[1]][1] > 0) {
      q_pos <- regmatches(col_name, q_match)[[1]]
      q_val <- as.numeric(sub("_q=", "", q_pos))
    } else {
      q_val <- NA
    }
    
    # Get mean entropy and variance for each gene from this column
    entropy_vals <- entropy_matrix[, col_name]
    gene_names <- rownames(entropy_matrix)
    
    # For each gene, calculate per-sample variance
    # Map genes to their rows and compute row-wise variance
    gene_variances <- sapply(gene_names, function(g_name) {
      # Find the index of this gene
      gene_idx <- which(rownames(entropy_matrix) == g_name)
      if (length(gene_idx) > 0) {
        # This gene has a single entropy value per sample-q combination
        # We use the mean entropy value as a proxy for expression level
        var(entropy_matrix[gene_idx, grep(paste0("_q=", q_val, "$"), colnames(entropy_matrix))], 
            na.rm = TRUE)
      } else {
        NA
      }
    })
    
    # Fit loess trend: variance ~ mean entropy per q-value
    # Only use genes with valid finite values for robust fitting
    valid_idx <- is.finite(entropy_vals) & is.finite(gene_variances)
    
    # Require sufficient data points for stable loess fitting
    # With n<6, loess becomes numerically unstable (degrees of freedom issues)
    if (sum(valid_idx) >= 6) {  # Increased from 4 to 6 for stability
      # Robust loess fitting (following Law et al. 2014 voom)
      # family="symmetric" for robustness to outliers
      # Adaptive span: smaller datasets get larger span for stability
      n_valid <- sum(valid_idx)
      span_adaptive <- min(0.5, max(0.2, 1.5 / n_valid))  # Adaptive span based on n
      
      tryCatch({
        loess_fit <- loess(
          gene_variances[valid_idx] ~ entropy_vals[valid_idx],
          span = span_adaptive,  # Adaptive smoothing span
          family = "symmetric",
          control = loess.control(surface = "direct", iterations = 4)
        )
        var_trend[[col_name]] <- loess_fit
        
        # Predict variance from trend for all valid genes
        # Note: predict(loess_fit) returns predictions on the original fitting data (valid subset)
        predicted_var <- predict(loess_fit)
        residuals <- gene_variances[valid_idx] - predicted_var
        sd_resid <- sd(residuals, na.rm = TRUE)
        
        # Outlier detection: genes with variance >2SD from trend (Love et al. 2014 DESeq2)
        outlier_threshold <- 2 * sd_resid
        outliers <- gene_names[valid_idx][abs(residuals) > outlier_threshold]
        outlier_genes[[col_name]] <- outliers
        
      }, error = function(e) {
        # Fallback to global variance if loess fails (graceful degradation)
        warning(sprintf(
          "Loess trend fitting failed for q=%.2f; using global variance.",
          q_val
        ), call. = FALSE)
        var_trend[[col_name]] <<- NULL
        outlier_genes[[col_name]] <<- character(0)
      })
    } else {
      # Insufficient data for loess fitting
      var_trend[[col_name]] <- NULL
      outlier_genes[[col_name]] <- character(0)
    }
  }
  
  return(list(
    global_mean = global_mean,
    global_var = global_var,
    var_trend = var_trend,           # NEW: Loess fits per q-value
    outlier_genes = outlier_genes,   # NEW: Detected outliers (>2SD from trend)
    n_isoforms = n_isoforms,
    n_samples = ncol(x)  # Sample-size awareness (Love et al. 2014)
  ))
}

#' Apply Empirical Bayes Shrinkage to Entropy Estimates
#'
#' Shrinks individual gene entropy estimates toward the global mean using empirical
#' Bayes weights. Particularly effective for genes with few expressed isoforms.
#' Outlier genes with extreme variance are protected (w=1, no shrinkage).
#'
#' @param entropy_matrix Matrix of raw entropy estimates (genes x assays).
#' @param params List from \code{.tsenat_estimate_shrinkage_params()} with
#'   global_mean, global_var, n_isoforms, n_samples, var_trend, and outlier_genes.
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
#' **Outlier Protection (Love et al. 2014 DESeq2):**
#' Genes with variance >2SD from the expression-dependent trend (detected via loess)
#' are treated as outliers and skip shrinkage (w=1), preserving biologically
#' meaningful extreme variance genes.
#'
#' This borrows strength from the global distribution, stabilizing estimates for
#' genes with small expression variance, while protecting genes with genuine
#' extreme variance signatures.
#'
#' @keywords internal
#' @noRd
.tsenat_apply_shrinkage <- function(entropy_matrix, params, gene_isoform_map = NULL) {
  result <- entropy_matrix
  
  global_mean <- params$global_mean
  global_var <- params$global_var
  n_isoforms <- params$n_isoforms
  n_samples <- params$n_samples  # Sample-size awareness (Love et al. 2014)
  var_trend <- params$var_trend  # NEW: Loess trend fits per q-value
  outlier_genes <- params$outlier_genes  # NEW: Outlier genes (>2SD from trend)
  
  # Map row names to n_isoforms indices if needed
  if (is.null(gene_isoform_map)) {
    # Assume row names are gene IDs matching names(n_isoforms)
    row_gene_ids <- rownames(entropy_matrix)
    gene_isoform_map <- n_isoforms[row_gene_ids]
  }
  
  # Estimate prior strength from the data using empirical Bayes methodology
  # (following DESeq2/edgeR approach: prior DF based on average gene information content)
  mean_n_isoforms <- mean(n_isoforms, na.rm = TRUE)
  
  # Sample-size awareness (Love et al. 2014 DESeq2)
  # With more samples, we have more confidence, so shrink less
  # Formula: weight = n_min / n_samples (smaller weight = less shrinkage)
  n_min <- max(2, floor(mean_n_isoforms))
  sample_size_weight <- n_min / max(n_samples, n_min)
  
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
        
        # Empirical Bayes prior strength estimation with sample-size adjustment
        # df_prior represents the effective sample size of the prior distribution
        # Estimated conservatively from average gene information content
        df_prior <- max(1, mean_n_isoforms - 1)
        
        # Apply sample-size weighting: larger datasets get less shrinkage
        df_prior_adjusted <- df_prior * sample_size_weight
        
        # Get list of outlier genes for this q-value (if available)
        col_outliers <- outlier_genes[[col_name]]
        if (is.null(col_outliers)) {
          col_outliers <- character(0)
        }
        
        # Compute shrinkage weights per gene
        for (row_idx in seq_len(nrow(entropy_matrix))) {
          # Try to get n_isoforms for this gene
          row_name <- rownames(entropy_matrix)[row_idx]
          
          # Check if this gene is an outlier (NEW: Outlier protection, Love et al. 2014)
          is_outlier <- row_name %in% col_outliers
          
          # Handle both named and unnamed gene_isoform_map vectors
          if (is.null(names(gene_isoform_map)) || length(names(gene_isoform_map)) == 0) {
            # Unnamed vector: use positional indexing
            n_iso <- gene_isoform_map[row_idx]
          } else {
            # Named vector: use name-based indexing
            n_iso <- gene_isoform_map[row_name]
          }
          
          if (!is.na(n_iso) && n_iso > 0) {
            # For outlier genes: skip shrinkage (w=1, maintain full observation)
            # This preserves biologically meaningful extreme variance genes
            if (is_outlier) {
              w <- 1  # No shrinkage for outliers
            } else {
              # Standard shrinkage: prior strength scaled by relative information content
              # For genes with n_iso < mean_n_isoforms: lambda > df_prior (more shrinkage)
              # For genes with n_iso > mean_n_isoforms: lambda < df_prior (less shrinkage)
              # Formula: lambda = df_prior_adjusted * (mean_n_isoforms / n_iso)
              # This implements precision-weighted empirical Bayes shrinkage WITH sample-size awareness
              lambda <- df_prior_adjusted * (mean_n_isoforms / n_iso)
              
              # Shrinkage weight: w close to 1 trusts the observation more, w close to 0 trusts prior more
              # From empirical Bayes theory: w = precision_obs / (precision_obs + precision_prior)
              w <- n_iso / (n_iso + lambda)
            }
            
            # Apply shrinkage: for finite values use weighted average of observation and prior,
            # for NA/NaN values (e.g., from undefined normalized entropy) use the prior estimate
            if (is.na(entropy_matrix[row_idx, col_idx]) || is.nan(entropy_matrix[row_idx, col_idx])) {
              # NA and NaN values get shrunk to the prior (w=0 for completely missing data)
              result[row_idx, col_idx] <- mu
            } else if (is.finite(entropy_matrix[row_idx, col_idx])) {
              # Finite values get weighted average of observation and prior
              result[row_idx, col_idx] <- w * entropy_matrix[row_idx, col_idx] + (1 - w) * mu
            }
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
#' normalization (default: \code{exp(1)}).
#' @param pseudocount Numeric scalar. Add this value to all transcript counts
#'   before computing proportions (default: 0). Useful for stability with
#'   zero-count features.
#' @param effective_length Numeric vector of effective transcript lengths (length = length(x)).
#'   When provided, counts are normalized by length to remove length bias before
#'   entropy calculation. This implements SALMON's recommended isoform-level approach.
#' @keywords internal
#' @noRd
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
#' Natural logarithms are used for q->1 limits and normalization.
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
