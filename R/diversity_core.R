

# ============================================================================
# INTERNAL HELPERS FOR REFACTORED .calculate_diversity()
# ============================================================================

#' Validate and normalize diversity parameters

#' @noRd
.validate_diversity_parameters <- function(norm, q, what, shrinkage, pseudocount) {
  # Coerce logical norm to character for backward compatibility
  if (is.logical(norm)) {
    norm <- if (norm) "range" else "none"
  }
  norm <- match.arg(norm, choices = c("none", "range", "zscore", 
                                      "log_odds_ratio", "relative_reference"))
  
  # Validate q values
  if (!is.numeric(q) || any(q < 0)) {
    stop("Argument 'q' must be numeric and >= 0 (q=0 represents species richness).", call. = FALSE)
  }
  
  # Validate other parameters
  what <- match.arg(what, choices = c("S", "D"))
  shrinkage <- match.arg(shrinkage, choices = c("none", "empirical_bayes"))
  
  # Validate pseudocount
  if (!is.numeric(pseudocount) && !(is.character(pseudocount) && tolower(pseudocount) == "auto")) {
    stop("pseudocount must be numeric or 'auto'", call. = FALSE)
  }
  
  list(norm = norm, q = q, what = what, shrinkage = shrinkage, pseudocount = pseudocount)
}

#' Handle pseudocount auto-estimation

#' @noRd
.handle_pseudocount_auto <- function(pseudocount, x, verbose) {
  if (!is.character(pseudocount) || tolower(pseudocount) != "auto") {
    return(pseudocount)  # Return as-is if not "auto"
  }
  
  if (verbose) {
    message("Computing pseudocount automatically via .estimate_pseudocount()...")
  }
  pc_result <- .estimate_pseudocount(x, verbose = FALSE)
  pseudocount <- pc_result$scalar_pseudocount
  
  if (verbose) {
    message(sprintf("  -> Estimated pseudocount = %.4f", pseudocount))
  }
  
  pseudocount
}

#' Extract gene names from SummarizedExperiment rowData

#' @noRd
.extract_gene_names <- function(original_x, genes, result) {
  gene_names <- NULL
  
  if (!is(original_x, "SummarizedExperiment") && !is(original_x, "RangedSummarizedExperiment")) {
    return(NULL)
  }
  
  rd <- try(SummarizedExperiment::rowData(original_x), silent = TRUE)
  if (inherits(rd, "try-error") || is.null(rd)) {
    return(NULL)
  }
  
  # Look for gene_names or gene_name column
  gene_name_col <- if ("gene_names" %in% colnames(rd)) {
    "gene_names"
  } else if ("gene_name" %in% colnames(rd)) {
    "gene_name"
  } else {
    NULL
  }
  
  if (is.null(gene_name_col)) {
    return(NULL)
  }
  
  # Use tapply for vectorized gene-to-name mapping
  tx_genes <- genes
  gene_names_col <- rd[[gene_name_col]]
  
  gene_to_name <- tapply(gene_names_col, tx_genes, function(x) {
    x_valid <- x[!is.na(x)]
    if (length(x_valid) > 0) x_valid[1] else NA
  }, simplify = FALSE)
  
  result_genes <- result[, 1]
  gene_names <- unname(gene_to_name[as.character(result_genes)])
  gene_names[is.na(gene_names)] <- result_genes[is.na(gene_names)]
  
  # Check for duplicates; if found, return NULL to fall back to gene IDs
  if (length(unique(gene_names)) < length(gene_names)) {
    return(NULL)
  }
  
  gene_names
}

#' Prepare column and row data for diversity result#' Prepare and validate diversity input data
#'
#' Internal helper: Prepares input matrix, validates dimensions, looks up effective_length
#' in metadata, and computes initial diversity values via .calculate_method().
#'

#' @noRd
.prepare_diversity_data <- function(x, genes, original_x, effective_length, 
    norm, q, what, nthreads, shrinkage, pseudocount, min_valid_frac, verbose, 
    tpm, assayno) {
    
    # Prepare input data
    inp <- .prepare_diversity_input(x = x, genes = genes, tpm = tpm, assayno = assayno,
        verbose = verbose)
    x <- inp$x
    genes <- inp$genes
    se_assay_mat <- inp$se_assay_mat
    
    # Validate data
    if (!is.numeric(x) || any(is.na(x))) {
        stop("Input data must be numeric and contain no NAs!", call. = FALSE)
    }
    if (nrow(x) != length(genes)) {
        stop("The number of rows is not equal to the given gene set.", call. = FALSE)
    }
    
    if (is.null(se_assay_mat)) {
        se_assay_mat <- x
    }
    
    # Look for effective_length in metadata if not provided
    if (is.null(effective_length) && (is(original_x, "SummarizedExperiment") || is(original_x, "RangedSummarizedExperiment"))) {
        md <- tryCatch(S4Vectors::metadata(original_x), error = function(e) NULL)
        if (!is.null(md) && !is.null(md$salmon_effective_length)) {
            effective_length <- md$salmon_effective_length
            if (verbose) message("[OK] Found salmon_effective_length in input metadata")
        }
    }
    
    # Calculate diversity
    use_range_norm <- (norm == "range")
    if (!is.null(effective_length) && verbose) {
        message("Calculating diversity with EFFECTIVE LENGTH NORMALIZATION")
    }
    
    if (verbose) {
        message("[DEBUG .prepare_diversity_data] About to call .calculate_method:")
        message("  x dimensions: ", nrow(x), " x ", ncol(x))
        message("  genes: ", paste(unique(genes), collapse=", "))
        message("  q values: ", paste(q, collapse=", "))
    }
    
    result <- .calculate_method(x, genes, use_range_norm, verbose = verbose, q = q, what = what,
        nthreads = nthreads, pseudocount = pseudocount, min_valid_frac = min_valid_frac,
        shrinkage = shrinkage, effective_length = effective_length)
    
    if (verbose) {
        message("[DEBUG .prepare_diversity_data] After .calculate_method:")
        message("  result dimensions: ", nrow(result), " x ", ncol(result))
    }
    
    list(result = result, x = x, genes = genes, se_assay_mat = se_assay_mat, effective_length = effective_length)
}

# NOTE (March 2026): .bootstrap_diversity_ci() moved to bootstrap.R for consolidation

#' Build SummarizedExperiment output for diversity
#'
#' Internal helper: Constructs output SE with assays, rowData, colData, and metadata.
#'

#' @noRd
.build_diversity_se_output <- function(result, output_structure, original_x, se_assay_mat,
    bootstrap_ci_results, bootstrap, metadata, verbose, what, q, genes) {
    
    result_assay <- output_structure$result_assay
    # Do NOT set rownames/colnames here - let SE constructor handle alignment with rowData/colData
    
    # Build assays list: include both diversity and counts assays
    # Aggregate transcript-level se_assay_mat to gene-level to match result_assay rows
    # The gene_ids are in the first column of result (before we extracted the assay)
    filtered_gene_ids <- as.character(result[, 1])
    
    # For each filtered gene, aggregate its transcript counts
    if (length(filtered_gene_ids) > 0) {
        counts_assay <- do.call(rbind, lapply(filtered_gene_ids, function(gene_id) {
            # Find transcripts for this gene using the genes vector mapping
            # The genes vector maps transcripts to gene names
            tx_mask <- which(genes == gene_id)
            if (length(tx_mask) > 0) {
                colSums(se_assay_mat[tx_mask, , drop = FALSE])
            } else {
                rep(0, ncol(se_assay_mat))
            }
        }))
    } else {
        # Empty result - no genes passed filtering
        counts_assay <- matrix(nrow = 0, ncol = ncol(se_assay_mat))
    }
    
    # For multi-q case, replicate counts for each q value
    # Check independently of result_assay rows to handle empty result case
    n_q <- length(output_structure$col_ids) / ncol(counts_assay)
    if (is.finite(n_q) && n_q == as.integer(n_q) && n_q > 1) {
        # Need to replicate counts_assay for each q value
        counts_assay <- do.call(cbind, replicate(as.integer(n_q), counts_assay, simplify = FALSE))
    }
    
    # Set assay dimnames to match rowData/colData rownames
    # This ensures compatibility with SummarizedExperiment constructor
    rownames(result_assay) <- rownames(output_structure$rowData)
    colnames(result_assay) <- rownames(output_structure$colData)
    rownames(counts_assay) <- rownames(output_structure$rowData)
    colnames(counts_assay) <- rownames(output_structure$colData)
    
    # Create assays list with appropriate name based on what parameter
    assay_name <- if (what[1] == "D") "hill" else "diversity"
    assays_list <- if (nrow(result_assay) > 0) {
        list(result_assay, counts_assay)
    } else {
        list(result_assay)  # Empty result - just include the empty assay
    }
    names(assays_list)[1] <- assay_name
    if (length(assays_list) > 1) names(assays_list)[2] <- "counts"
    if (!is.null(bootstrap_ci_results) && !is.null(bootstrap_ci_results$bootstrap_ci_results)) {
        # Extract CI matrices from bootstrap results
        bootstrap_out <- bootstrap_ci_results$bootstrap_ci_results
        
        # Create CI matrices: initialize with NA, then populate from bootstrap output
        ci_lower <- result_assay * NA_real_
        ci_upper <- result_assay * NA_real_
        
        # Extract bootstrap CIs if available in output
        # NEW STRUCTURE: each bootstrap result is for ONE (gene, sample) pair
        if (is.list(bootstrap_out) && length(bootstrap_out) > 0) {
            if (verbose) {
                message("[CI EXTRACT] Processing bootstrap_out: ", length(bootstrap_out), " items")
                if (length(bootstrap_out) > 0) {
                    message("[CI EXTRACT] First item structure:")
                    str(bootstrap_out[[1]], max.level=2)
                    message("[CI EXTRACT] First item name: ", names(bootstrap_out)[1])
                }
            }
            
            # Get row names from the result_assay to map (gene, sample) pairs to indices
            result_row_names  <- rownames(result_assay)
            result_col_names <- colnames(result_assay)
            
            # Process each bootstrap result
            # Names should be like "gene_name_sample_1", "gene_name_sample_2", etc.
            for (i in seq_along(bootstrap_out)) {
                boot_item <- bootstrap_out[[i]]
                if (is.null(boot_item)) next
                
                # Parse bootstrap result name: "gene_name_sample_INDEX"
                boot_name <- names(bootstrap_out)[i]
                if (is.null(boot_name) || is.na(boot_name)) next
                
                m <- regexec("^(.+)_sample_([0-9]+)$", boot_name)
                parts <- regmatches(boot_name, m)
                if (length(parts[[1]]) != 3) next
                
                gene_name <- parts[[1]][2]
                sample_idx_str <- parts[[1]][3]
                sample_idx <- as.integer(sample_idx_str)
                
                # Get actual sample name from original se_assay_mat
                if (sample_idx < 1 || sample_idx > ncol(se_assay_mat)) next
                sample_name <- colnames(se_assay_mat)[sample_idx]
                
                # Find gene row in result matrix
                gene_row_idx <- which(result_row_names == gene_name)[1]
                if (is.na(gene_row_idx)) next
                
                # Find columns for this sample in result_assay
                # Columns should have format "SAMPLE_q=NUMBER"
                col_pattern <- paste0("^", gsub("([.^$*+?{}\\(\\)\\[\\]|\\\\])", "\\\\\\1", sample_name), "_q=")
                col_indices <- grep(col_pattern, result_col_names)
                if (length(col_indices) == 0) next
                
                # Extract q values from matching columns
                col_q_values <- sub(
                    paste0("^", gsub("([.^$*+?{}\\(\\)\\[\\]|\\\\])", "\\\\\\1", sample_name), "_q="),
                    "",
                    result_col_names[col_indices])
                
                # Process bootstrap result depending on structure
                if (is.list(boot_item) && !is.null(names(boot_item)) && all(grepl("^q=", names(boot_item)))) {
                    # Multi-q case: bootstrap result is a list with names like "q=1.0"
                    for (j in seq_along(boot_item)) {
                        q_name <- names(boot_item)[j]  # "q=1.0"
                        q_val <- sub("^q=", "", q_name)
                        q_result <- boot_item[[j]]
                        
                        if (!is.null(q_result$lower_ci) && !is.null(q_result$upper_ci)) {
                            # Find which columns match this q value
                            matching_q_idx <- which(col_q_values == q_val)
                            if (length(matching_q_idx) > 0) {
                                target_col_indices <- col_indices[matching_q_idx]
                                
                                if (verbose && i <= 1) {
                                    # Get the point estimate from the result matrix to compare
                                    result_value <- result_assay[gene_row_idx, target_col_indices[1]]
                                    message("[DEBUG CI POPULATE] gene=", gene_name, " sample=", sample_idx, " q=", q_val)
                                    message("  result_value (point_est should match this): ", result_value)
                                    message("  bootstrap point_est:", q_result$estimate)
                                    message("  ci_lower:   ", q_result$lower_ci)
                                    message("  ci_upper:   ", q_result$upper_ci)
                                    message("  target col indices: ", paste(target_col_indices, collapse=", "))
                                }
                                
                                ci_lower[gene_row_idx, target_col_indices] <- as.numeric(q_result$lower_ci)[1]
                                ci_upper[gene_row_idx, target_col_indices] <- as.numeric(q_result$upper_ci)[1]
                            }
                        }
                    }
                } else if (!is.null(boot_item$lower_ci) && !is.null(boot_item$upper_ci)) {
                    # Single-q case: bootstrap result is a tsenat_bootstrap_ci object
                    ci_lower[gene_row_idx, col_indices] <- as.numeric(boot_item$lower_ci)[1]
                    ci_upper[gene_row_idx, col_indices] <- as.numeric(boot_item$upper_ci)[1]
                }
            }
        }
        
        # Set dimnames on CI matrices to match result_assay
        rownames(ci_lower) <- rownames(result_assay)
        colnames(ci_lower) <- colnames(result_assay)
        rownames(ci_upper) <- rownames(result_assay)
        colnames(ci_upper) <- colnames(result_assay)
        
        # Only add CI assays if they have actual values (not all NA)
        if (!all(is.na(ci_lower)) && !all(is.na(ci_upper))) {
            assays_list$ci_lower <- ci_lower
            assays_list$ci_upper <- ci_upper
            if (verbose) message("  [OK] Added ci_lower and ci_upper assays to output SE")
        } else if (verbose) {
            message("  [INFO] Bootstrap CIs extracted but no valid values found; skipping CI assays")
        }
    }
    
    # Build metadata
    result_meta_list <- list(
        q = q,  # Store q values used
        what = what[1],  # Store which metric (S or D)
        readcounts = se_assay_mat,
        bootstrap = bootstrap,
        bootstrap_nboot = if (!is.null(bootstrap_ci_results)) bootstrap_ci_results$bootstrap_nboot else NULL,
        bootstrap_method = if (!is.null(bootstrap_ci_results)) bootstrap_ci_results$bootstrap_method else NULL,
        bootstrap_ci = if (!is.null(bootstrap_ci_results)) bootstrap_ci_results$bootstrap_ci else NULL
    )
    if (is(original_x, "SummarizedExperiment") || is(original_x, "RangedSummarizedExperiment")) {
        result_meta_list$se <- original_x
    }
    
    # Build and return SummarizedExperiment
    result <- SummarizedExperiment::SummarizedExperiment(
        assays = assays_list,
        rowData = output_structure$rowData,
        colData = output_structure$colData,
        metadata = result_meta_list
    )
    
    # Apply metadata mapping if provided
    if (!is.null(metadata)) {
        result <- .map_metadata_se(result, metadata)
    }
    
    result
}

#' @noRd
.prepare_diversity_metadata <- function(x, result, original_x, genes, q, gene_names = NULL) {
  # Get gene IDs from first column of result
  gene_ids <- as.character(result[, 1])
  row_ids <- if (!is.null(gene_names)) gene_names else gene_ids
  
  # Create rowData with row_ids as rownames (not gene_ids)
  result_rowData <- data.frame(gene_id = gene_ids, row.names = row_ids)
  result_rowData$gene_name <- gene_names
  
  if (length(q) > 1) {
    col_split <- do.call(rbind, strsplit(colnames(result)[-1], "_q="))
    col_ids <- paste0(col_split[, 1], "_q=", col_split[, 2])
    
    result_colData <- data.frame(
      samples = as.character(col_split[, 1]),
      q = as.numeric(col_split[, 2]),
      row.names = col_ids,
      stringsAsFactors = FALSE
    )
    
    # Preserve original colData if available
    if (is(original_x, "SummarizedExperiment") || is(original_x, "RangedSummarizedExperiment")) {
      orig_coldata <- try(SummarizedExperiment::colData(original_x), silent = TRUE)
      if (!inherits(orig_coldata, "try-error") && nrow(orig_coldata) > 0) {
        sample_indices <- NA
        if (length(rownames(orig_coldata)) > 0 && rownames(orig_coldata)[1] != "") {
          sample_indices <- match(col_split[, 1], rownames(orig_coldata))
        } else if ("Sample" %in% colnames(orig_coldata)) {
          sample_indices <- match(col_split[, 1], as.character(orig_coldata$Sample))
        }
        
        if (!all(is.na(sample_indices))) {
          for (col in colnames(orig_coldata)) {
            result_colData[[col]] <- orig_coldata[[col]][sample_indices]
          }
        }
      }
    }
  } else {
    # Single q-value
    base_col_ids <- colnames(x)
    if (is.null(base_col_ids) || any(base_col_ids == "")) {
      base_col_ids <- paste0("Sample", seq_len(ncol(x)))
    }
    q_formatted <- formatC(q, format = "f", digits = 3)
    col_ids <- paste0(base_col_ids, "_q=", q_formatted)
    
    result_colData <- data.frame(
      samples = base_col_ids,
      q = rep(q, length(base_col_ids)),
      row.names = col_ids,
      stringsAsFactors = FALSE
    )
    
    # Preserve original colData
    if (is(original_x, "SummarizedExperiment") || is(original_x, "RangedSummarizedExperiment")) {
      orig_coldata <- try(SummarizedExperiment::colData(original_x), silent = TRUE)
      if (!inherits(orig_coldata, "try-error") && nrow(orig_coldata) > 0) {
        for (col in colnames(orig_coldata)) {
          result_colData[[col]] <- orig_coldata[[col]]
        }
      }
    }
  }
  
  list(colData = result_colData, rowData = result_rowData, 
       col_ids = col_ids, row_ids = row_ids, result_assay = as.matrix(result[, -1, drop = FALSE]))
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
#' automatically estimated using library size adjustment via `.estimate_pseudocount()`
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
#' Example: load(readcounts.RData'); .calculate_diversity(salmon_dataset, effective_length=salmon_effective_length)
#' @param bootstrap Logical; if TRUE, compute bootstrap confidence intervals around
#' Tsallis entropy point estimates using \code{.calculate_tsallis_entropy_bootstrap()}.
#' Default: FALSE (disabled for backward compatibility). When TRUE, computes CIs for
#' each gene and adds assays: ci_lower and ci_upper to output.
#' @param bootstrap_nboot Integer; number of bootstrap replicates (default: NULL).
#' If NULL, automatically suggests nboot based on number of genes using \code{.suggest_nboot()}.
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
#' `.map_metadata_se()`. This allows adding additional context or derived annotations to
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
#' se <- .calculate_diversity(counts, genes = genes, q = 1.0, norm = TRUE)
#' head(SummarizedExperiment::assay(se))
#' 

#' @noRd

.calculate_diversity <- function(x, genes = NULL, norm = TRUE, tpm = FALSE, assayno = 1,
    verbose = FALSE, q = 2, what = c("S", "D"), nthreads = 1, pseudocount = 0, 
    min_valid_frac = 0.75, shrinkage = "none", effective_length = NULL, metadata = NULL,
    bootstrap = FALSE, bootstrap_nboot = NULL, bootstrap_method = "percentile",
    bootstrap_ci = 0.95, bootstrap_include_diagnostics = TRUE, seed = NULL) {
    
    # Store original input and validate parameters
    original_x <- x
    validated <- .validate_diversity_parameters(norm, q, what, shrinkage, pseudocount)
    norm <- validated$norm
    q <- validated$q
    what <- validated$what
    shrinkage <- validated$shrinkage
    
    # Handle pseudocount auto-estimation
    pseudocount <- .handle_pseudocount_auto(pseudocount, x, verbose)
    
    # Prepare input and calculate diversity
    prep <- .prepare_diversity_data(x, genes, original_x, effective_length,
        norm, q, what, nthreads, shrinkage, pseudocount, min_valid_frac, verbose,
        tpm, assayno)
    result <- prep$result
    x <- prep$x
    genes <- prep$genes
    se_assay_mat <- prep$se_assay_mat
    effective_length <- prep$effective_length  # Extract effective_length from prep result
    
    if (verbose) {
        message("[DEBUG .calculate_diversity] effective_length after prep:")
        message("  is.null=", is.null(effective_length))
        if (!is.null(effective_length)) {
            message("  class=", class(effective_length))
            message("  length=", length(effective_length))
        }
    }
    
    if (verbose && nrow(result) == 0) {
        message("[WARN] Result from .prepare_diversity_data() is empty (0 rows)")
        message("       This happens when bootstrap=", bootstrap)
    }
    
    # Optional: Compute bootstrap CIs
    bootstrap_ci_results <- .bootstrap_diversity_ci(bootstrap, result, genes, se_assay_mat,
        bootstrap_method, bootstrap_ci, bootstrap_nboot, q, pseudocount, nthreads,
        bootstrap_include_diagnostics, verbose, seed, effective_length)
    
    # Prepare output structure
    gene_names <- .extract_gene_names(original_x, genes, result)
    output_structure <- .prepare_diversity_metadata(x, result, original_x, genes, q, gene_names)
    
    # Build and return SummarizedExperiment
    .build_diversity_se_output(result, output_structure, original_x, se_assay_mat,
        bootstrap_ci_results, bootstrap, metadata, verbose, what, q, genes)
}

