#' Create a Complete TSENATAnalysis with Test Data
#' 
#' Factory function that generates a fully initialized TSENATAnalysis object
#' with realistic data and biological signal. Useful for examples, testing,
#' and documentation. Eliminates ~70% of boilerplate in documentation examples.
#'
#' @param n_genes Number of genes to simulate (default: 8)
#' @param n_samples_per_group Samples per condition (default: 20)
#' @param control_lambda Poisson lambda for control condition (default: 40)
#' @param treatment_lambda Poisson lambda for treatment condition (default: 150)
#' @param q_values Vector of q-values for diversity calculation (default: c(0.5, 1.0, 1.5))
#' @param include_divergence If TRUE, compute divergence results (default: TRUE)
#' @param include_lm_results If TRUE, add placeholder LM results (default: TRUE)
#' @param seed Random seed for reproducibility (default: 42)
#' @param verbose Logical for progress messages (default: FALSE)
#'
#' @return TSENATAnalysis object with:
#'   - SummarizedExperiment with count matrix, rowData, and colData
#'   - Computed diversity results across q-values
#'   - Computed divergence results (optional)
#'   - Placeholder LM results (optional)
#'   - Proper tx2gene metadata mapping
#'   - TPM data in metadata
#'
#' @details
#' The factory ensures:
#' - Biological signal: control (lambda=40) vs treatment (lambda=150) contrast
#' - Sufficient samples: 40 total (20 per group) for stable LM fitting
#' - Multiple q-values: c(0.5, 1.0, 1.5) avoids rank deficiency
#' - Valid S4 object structure: passes all TSENATAnalysis validity checks
#' - Optional divergence and LM results to support testing without warnings
#' - All required metadata (tx2gene, TPM, rowData) pre-configured
#'
#' @examples
#' # Create with defaults (8 genes, 20 samples/group, multi-q, with all results)
#' analysis <- create_test_analysis()
#' 
#' # Create minimal analysis for quick testing
#' analysis <- create_test_analysis(
#'   n_genes = 4,
#'   n_samples_per_group = 10,
#'   q_values = c(0.5, 1.0),
#'   include_divergence = FALSE
#' )
#' 
#' # Create with custom parameters
#' analysis <- create_test_analysis(
#'   n_genes = 16,
#'   n_samples_per_group = 30,
#'   control_lambda = 50,
#'   treatment_lambda = 200,
#'   q_values = c(0.1, 0.5, 1.0, 1.5, 2.0),
#'   include_divergence = TRUE,
#'   include_lm_results = TRUE
#' )
#'
#' @keywords internal
#' @noRd
.create_test_analysis <- function(
    n_genes = 8,
    n_samples_per_group = 20,
    control_lambda = 40,
    treatment_lambda = 150,
    q_values = c(0.5, 1.0, 1.5),
    include_divergence = TRUE,
    include_lm_results = TRUE,
    seed = 42,
    verbose = FALSE) {
  
  set.seed(seed)
  
  # Dimensions
  n_samples <- n_samples_per_group * 2
  n_transcripts <- n_genes * 50
  
  if (verbose) {
    message("[create_test_analysis] Generating ", n_transcripts, " transcripts across ",
            n_genes, " genes with ", n_samples, " samples")
  }
  
  # Generate counts with biological signal
  control_idx <- seq(1, n_samples, by = 2)
  treatment_idx <- seq(2, n_samples, by = 2)
  
  counts <- matrix(0, nrow = n_transcripts, ncol = n_samples)
  for (j in seq_len(n_samples)) {
    if (j %in% control_idx) {
      counts[, j] <- rpois(n_transcripts, lambda = control_lambda)
    } else {
      counts[, j] <- rpois(n_transcripts, lambda = treatment_lambda)
    }
  }
  counts <- pmax(counts, 50)
  
  rownames(counts) <- paste0("TX_", 1:n_transcripts)
  colnames(counts) <- paste0("Sample_", 1:n_samples)
  
  # Create rowData with gene mappings
  rowData <- S4Vectors::DataFrame(
    transcript_id = rownames(counts),
    gene_id = paste0("GENE_", rep(1:n_genes, each = 50, length.out = n_transcripts)),
    row.names = rownames(counts)
  )
  
  # Create colData with experimental design
  colData <- S4Vectors::DataFrame(
    sample_id = colnames(counts),
    condition = rep(c("control", "treatment"), length.out = n_samples),
    sample_type = rep(c("typeA", "typeB"), length.out = n_samples),
    subject = rep(paste0("S", 1:10), length.out = n_samples),
    paired_samples = rep(paste0("pair", 1:10), length.out = n_samples),
    row.names = colnames(counts)
  )
  
  # Create SummarizedExperiment
  se <- SummarizedExperiment::SummarizedExperiment(
    assays = list(counts = counts),
    rowData = rowData,
    colData = colData
  )
  
  # Add tx2gene metadata
  tx2gene_df <- data.frame(
    Transcript = rownames(counts),
    Gene = rowData$gene_id,
    stringsAsFactors = FALSE
  )
  S4Vectors::metadata(se)$tx2gene <- tx2gene_df
  
  # Generate synthetic TPM data (matching counts dimensions)
  tpm <- counts
  for (j in seq_len(ncol(tpm))) {
    lib_size <- colSums(tpm[, j, drop = FALSE])
    if (lib_size > 0) {
      tpm[, j] <- (tpm[, j] / lib_size) * 1e6
    }
  }
  rownames(tpm) <- rownames(counts)
  colnames(tpm) <- colnames(counts)
  S4Vectors::metadata(se)$salmon_tpm <- tpm
  
  # Initialize TSENATAnalysis
  analysis <- TSENATAnalysis(se = se, config = list())
  
  # Calculate diversity
  analysis <- calculate_diversity_s4(
    analysis,
    q = q_values,
    verbose = FALSE,
    min_valid_frac = 0
  )
  
  # Calculate divergence if requested
  if (include_divergence) {
    analysis <- tryCatch({
      calculate_divergence_s4(
        analysis,
        verbose = FALSE
      )
    }, error = function(e) {
      # If divergence fails, continue without it
      if (verbose) {
        message("[create_test_analysis] Warning: divergence calculation failed: ", e$message)
      }
      analysis
    })
  }
  
  # Add placeholder LM results if requested
  if (include_lm_results) {
    # Create a simple placeholder LM result (empty data frame structure)
    # This prevents "No LM results found" warnings in tests
    lm_placeholder <- list(
      overall = data.frame(
        gene = character(0),
        term = character(0),
        estimate = numeric(0),
        std.error = numeric(0),
        statistic = numeric(0),
        p.value = numeric(0)
      )
    )
    analysis@lm_results <- lm_placeholder
  }
  
  if (verbose) {
    message("[create_test_analysis] Analysis created with ",
            length(q_values), " q-values")
    message("[create_test_analysis] Diversity results: ",
            nrow(diversity(analysis)), " genes")
    if (include_divergence) {
      message("[create_test_analysis] Divergence results included")
    }
    if (include_lm_results) {
      message("[create_test_analysis] Placeholder LM results included")
    }
  }
  
  return(analysis)
}


#' Subset a TSENATAnalysis Object for Testing and Examples
#'
#' Create a smaller, representative subset of a TSENATAnalysis object
#' for use in testing, examples, or documentation. Preserves all analysis
#' metadata and computed results while reducing dataset size.
#'
#' @param analysis TSENATAnalysis object to subset
#' @param n_genes Positive integer. Number of genes to retain (default: 10).
#'   Set to NULL to keep all genes.
#' @param n_samples Positive integer. Number of samples to retain (default: NULL,
#'   keep all samples). If specified, samples are selected to balance
#'   conditions when possible.
#' @param genes Character vector of specific gene IDs to retain. If provided,
#'   overrides n_genes argument (default: NULL).
#' @param samples Character vector of specific sample IDs to retain. If provided,
#'   overrides n_samples argument (default: NULL).
#' @param select_by One of "variance" (select genes with highest variance),
#'   "mean" (select genes with highest mean expression), or "random"
#'   (random selection). Default: "variance"
#' @param seed Random seed for reproducible subsetting (default: 42)
#' @param min_count Minimum total transcript count (across all samples) required
#'   for a gene to be included. Filters genes to ensure adequate data density
#'   for statistical operations like jackknife (default: NULL, no filtering).
#'   Useful values: 5-10 for robust estimates.
#' @param verbose Logical. Print progress messages (default: FALSE)
#'
#' @return A new TSENATAnalysis object containing only the specified genes
#'   and samples. All computed results (diversity, LM, jackknife, divergence)
#'   are automatically subsetted to match the new dimensions while preserving
#'   analysis configuration and metadata.
#'
#' @details
#' This function provides a convenient wrapper around the `[` subsetting
#' operator for TSENATAnalysis. It simplifies common use cases:
#'
#' \itemize{
#'   \item \strong{By gene count}: Automatically selects top genes by variance,
#'     mean expression, or at random
#'   \item \strong{By sample count}: Intelligently balances sample selection
#'     across experimental conditions
#'   \item \strong{By specific IDs}: Explicitly specify genes and samples to keep
#'   \item \strong{By statistic}: Select informative genes (high variance = more
#'     informative for testing)
#'   \item \strong{By abundance}: Filter by minimum total count to ensure
#'     adequate data density for operations like jackknife
#' }
#'
#' The function preserves all analysis structure:
#' - Diversity results are subsetted to match gene and sample selection
#' - Jackknife results maintain confidence intervals for selected samples
#' - LM results (gene-level statistics) remain intact
#' - Divergence calculations are updated if applicable
#' - Configuration and metadata are unchanged
#'
#' @examples
#' # Create minimal example TSENATAnalysis
#' library(SummarizedExperiment)
#' se <- SummarizedExperiment(
#'   assays = list(counts = matrix(rpois(1000, 10), nrow = 100, ncol = 10)),
#'   rowData = data.frame(
#'     gene_id = rep(paste0("GENE_", 1:10), 10),
#'     row.names = paste0("TX_", 1:100)
#'   ),
#'   colData = data.frame(
#'     sample_id = paste0("S_", 1:10),
#'     condition = rep(c("control", "treatment"), 5),
#'     row.names = paste0("S_", 1:10)
#'   )
#' )
#' analysis <- TSENATAnalysis(se)
#'
#' # Subset to 5 genes and 4 samples (top genes by variance)
#' small_analysis <- subset_analysis(analysis, n_genes = 30, n_samples = 8)
#'
#' # Subset to genes with minimum 1 total count across all samples
#' # This ensures data adequacy filtering (recommended: min_count = 10-20 for robust estimates)
#' filtered <- subset_analysis(analysis, n_genes = 100, min_count = 1)
#'
#' # Subset to specific genes only
#' subset_genes <- subset_analysis(
#'   analysis,
#'   genes = c("TX_1", "TX_2", "TX_3"),
#'   n_samples = 5
#' )
#'
#' # Random selection of genes (reproducible with seed)
#' random_subset <- subset_analysis(
#'   analysis,
#'   n_genes = 8,
#'   n_samples = 6,
#'   select_by = "random",
#'   seed = 123
#' )
#'
#' # Keep specific samples only
#' control_only <- subset_analysis(
#'   analysis,
#'   samples = colnames(analysis@se)[
#'     colData(analysis@se)$condition == "control"
#'   ]
#' )
#'
#' @export
#' @rdname subset_analysis
subset_analysis <- function(
    analysis,
    n_genes = 10,
    n_samples = NULL,
    genes = NULL,
    samples = NULL,
    select_by = c("variance", "mean", "random"),
    seed = 42,
    min_count = NULL,
    verbose = FALSE) {
  
  # Validate input
  if (!inherits(analysis, "TSENATAnalysis")) {
    stop("analysis must be a TSENATAnalysis object", call. = FALSE)
  }
  
  select_by <- match.arg(select_by)
  set.seed(seed)
  
  se <- analysis@se
  n_genes_total <- nrow(se)
  n_samples_total <- ncol(se)
  
  # ============================================================================
  # Select genes
  # ============================================================================
  
  if (!is.null(genes)) {
    # Use explicitly provided genes
    if (verbose) {
      message("[subset_analysis] Using provided genes: ", 
              paste(head(genes, 3), collapse = ", "), 
              if (length(genes) > 3) "...")
    }
    gene_idx <- match(genes, rownames(se))
    if (any(is.na(gene_idx))) {
      missing_genes <- genes[is.na(gene_idx)]
      stop("Genes not found in analysis: ", 
           paste(head(missing_genes, 3), collapse = ", "),
           call. = FALSE)
    }
  } else if (!is.null(n_genes)) {
    # Select by criterion
    n_genes <- as.integer(n_genes)
    if (n_genes < 1) {
      stop("n_genes must be >= 1", call. = FALSE)
    }
    if (n_genes > n_genes_total) {
      warning("n_genes (", n_genes, ") exceeds available genes (", 
              n_genes_total, "). Using all genes.", call. = FALSE)
      n_genes <- n_genes_total
    }
    
    counts_matrix <- assay(se, "counts")
    
    if (select_by == "variance") {
      # Compute variance per transcript
      tx_vars <- matrixStats::rowVars(counts_matrix)
      gene_idx <- order(tx_vars, decreasing = TRUE)[1:n_genes]
      if (verbose) {
        message("[subset_analysis] Selected ", n_genes, 
                " genes with highest variance")
      }
    } else if (select_by == "mean") {
      # Compute mean per transcript
      tx_means <- rowMeans(counts_matrix)
      gene_idx <- order(tx_means, decreasing = TRUE)[1:n_genes]
      if (verbose) {
        message("[subset_analysis] Selected ", n_genes, 
                " genes with highest mean expression")
      }
    } else if (select_by == "random") {
      # Random selection
      gene_idx <- sample(n_genes_total, n_genes)
      if (verbose) {
        message("[subset_analysis] Randomly selected ", n_genes, 
                " genes (seed = ", seed, ")")
      }
    }
  } else {
    # Keep all genes
    gene_idx <- seq_len(n_genes_total)
    if (verbose) {
      message("[subset_analysis] Keeping all ", n_genes_total, " genes")
    }
  }
  
  # ============================================================================
  # Select samples
  # ============================================================================
  
  if (!is.null(samples)) {
    # Use explicitly provided samples
    if (verbose) {
      message("[subset_analysis] Using provided samples: ", 
              paste(head(samples, 3), collapse = ", "), 
              if (length(samples) > 3) "...")
    }
    sample_idx <- match(samples, colnames(se))
    if (any(is.na(sample_idx))) {
      missing_samples <- samples[is.na(sample_idx)]
      stop("Samples not found in analysis: ", 
           paste(head(missing_samples, 3), collapse = ", "),
           call. = FALSE)
    }
  } else if (!is.null(n_samples)) {
    # Intelligently select samples across conditions
    n_samples <- as.integer(n_samples)
    if (n_samples < 1) {
      stop("n_samples must be >= 1", call. = FALSE)
    }
    if (n_samples > n_samples_total) {
      warning("n_samples (", n_samples, ") exceeds available samples (", 
              n_samples_total, "). Using all samples.", call. = FALSE)
      n_samples <- n_samples_total
    }
    
    coldata <- colData(se)
    
    # Try to balance by condition if available
    if ("condition" %in% colnames(coldata)) {
      conditions <- unique(coldata$condition)
      samples_per_cond <- ceiling(n_samples / length(conditions))
      
      sample_idx <- c()
      for (cond in conditions) {
        cond_idx <- which(coldata$condition == cond)
        selected <- sample(cond_idx, min(samples_per_cond, length(cond_idx)))
        sample_idx <- c(sample_idx, selected)
      }
      sample_idx <- head(sample_idx, n_samples)
      if (verbose) {
        message("[subset_analysis] Selected ", n_samples, 
                " samples balanced across ", length(conditions), " conditions")
      }
    } else if ("sample_type" %in% colnames(coldata)) {
      # Alternative: try sample_type
      types <- unique(coldata$sample_type)
      samples_per_type <- ceiling(n_samples / length(types))
      
      sample_idx <- c()
      for (type in types) {
        type_idx <- which(coldata$sample_type == type)
        selected <- sample(type_idx, min(samples_per_type, length(type_idx)))
        sample_idx <- c(sample_idx, selected)
      }
      sample_idx <- head(sample_idx, n_samples)
      if (verbose) {
        message("[subset_analysis] Selected ", n_samples, 
                " samples balanced across ", length(types), " sample types")
      }
    } else {
      # Random selection of samples
      sample_idx <- sample(n_samples_total, n_samples)
      if (verbose) {
        message("[subset_analysis] Randomly selected ", n_samples, 
                " samples (seed = ", seed, ")")
      }
    }
  } else {
    # Keep all samples
    sample_idx <- seq_len(n_samples_total)
    if (verbose) {
      message("[subset_analysis] Keeping all ", n_samples_total, " samples")
    }
  }
  
  # ============================================================================
  # Filter by minimum count (data adequacy check)
  # ============================================================================
  
  if (!is.null(min_count)) {
    min_count <- as.numeric(min_count)
    if (min_count < 0) {
      stop("min_count must be >= 0", call. = FALSE)
    }
    
    counts_matrix <- assay(se, "counts")[gene_idx, sample_idx, drop = FALSE]
    total_counts <- rowSums(counts_matrix)
    genes_keep <- total_counts >= min_count
    
    n_before_filter <- length(gene_idx)
    gene_idx_filtered <- gene_idx[genes_keep]
    n_after_filter <- length(gene_idx_filtered)
    n_removed <- n_before_filter - n_after_filter
    
    if (n_removed > 0) {
      warning("Filtered out ", n_removed, " gene(s) with total count < ", 
              min_count, " (", n_after_filter, " genes remain)", 
              call. = FALSE)
      if (verbose) {
        message("[subset_analysis] Total count range in selected samples: ", 
                format(min(total_counts), trim = TRUE), " - ", 
                format(max(total_counts), trim = TRUE))
        message("[subset_analysis] Retained genes with count >= ", min_count, 
                ": ", n_after_filter)
      }
    } else if (verbose) {
      message("[subset_analysis] All ", n_after_filter, " genes meet minimum count threshold (", 
              min_count, ")")
    }
    
    gene_idx <- gene_idx_filtered
    
    if (length(gene_idx) == 0) {
      stop("No genes meet the minimum count threshold (min_count = ", 
           min_count, "). Consider lowering min_count or using more/different samples.", 
           call. = FALSE)
    }
  }
  
  # ============================================================================
  # Subset using the [ operator
  # ============================================================================
  
  analysis_subset <- analysis[gene_idx, sample_idx]
  
  # Sync metadata to match the subsetted SE
  # This ensures consistency when calculate_diversity_s4() extracts data
  
  # 1. Sync tx2gene mapping (transcript-to-gene mapping)
  if (!is.null(S4Vectors::metadata(analysis@se)$tx2gene)) {
    tx2gene_full <- S4Vectors::metadata(analysis@se)$tx2gene
    tx2gene_subset <- tx2gene_full[tx2gene_full$Transcript %in% rownames(analysis_subset@se), ]
    S4Vectors::metadata(analysis_subset@se)$tx2gene <- tx2gene_subset
    
    if (verbose) {
      message("[subset_analysis] Filtered tx2gene: ", 
              nrow(tx2gene_full), " → ", nrow(tx2gene_subset), " transcripts")
    }
  }
  
  # 2. Sync readcounts (original count matrix stored in metadata)
  # This is critical because .prepare_diversity_input() checks for md$readcounts
  if (!is.null(S4Vectors::metadata(analysis@se)$readcounts)) {
    readcounts_full <- S4Vectors::metadata(analysis@se)$readcounts
    readcounts_subset <- readcounts_full[gene_idx, sample_idx, drop = FALSE]
    S4Vectors::metadata(analysis_subset@se)$readcounts <- readcounts_subset
    
    if (verbose) {
      message("[subset_analysis] Filtered readcounts: ", 
              nrow(readcounts_full), " → ", nrow(readcounts_subset), " transcripts")
    }
  }
  
  # 3. Sync salmon_tpm (TPM matrix stored in metadata) if present
  if (!is.null(S4Vectors::metadata(analysis@se)$salmon_tpm)) {
    tpm_full <- S4Vectors::metadata(analysis@se)$salmon_tpm
    tpm_subset <- tpm_full[gene_idx, sample_idx, drop = FALSE]
    S4Vectors::metadata(analysis_subset@se)$salmon_tpm <- tpm_subset
    
    if (verbose) {
      message("[subset_analysis] Filtered salmon_tpm: ", 
              nrow(tpm_full), " → ", nrow(tpm_subset), " transcripts")
    }
  }
  
  # 4. Sync salmon_effective_length (effective lengths) if present
  if (!is.null(S4Vectors::metadata(analysis@se)$salmon_effective_length)) {
    eff_len_full <- S4Vectors::metadata(analysis@se)$salmon_effective_length
    
    # Could be vector or matrix, handle both
    if (is.vector(eff_len_full)) {
      # Named vector - subset by matching names to selected genes
      tx_names_subset <- rownames(analysis_subset@se)
      eff_len_subset <- eff_len_full[na.omit(match(tx_names_subset, names(eff_len_full)))]
      if (length(eff_len_subset) > 0) {
        S4Vectors::metadata(analysis_subset@se)$salmon_effective_length <- eff_len_subset
      }
    } else if (is.matrix(eff_len_full)) {
      # Matrix - subset by rows and columns
      eff_len_subset <- eff_len_full[gene_idx, sample_idx, drop = FALSE]
      S4Vectors::metadata(analysis_subset@se)$salmon_effective_length <- eff_len_subset
    }
  }
  
  if (verbose) {
    message("[subset_analysis] Subset complete: ", 
            nrow(analysis_subset@se), " genes × ", 
            ncol(analysis_subset@se), " samples")
  }
  
  return(analysis_subset)
}


#' Load and Subset Vignette Data for Examples
#'
#' Loads the full TSENAT vignette dataset and returns a representative subset
#' for use in roxygen documentation examples. Useful for demonstrating analysis
#' workflows with real biological data while keeping examples computationally
#' efficient.
#'
#' @param n_genes Number of genes to retain in subset (default: 10)
#' @param n_samples Number of samples to retain in subset (default: 8)
#' @param seed Random seed for reproducible subset selection (default: 42)
#' @param select_top_var If TRUE, select genes by highest variance; if FALSE,
#'   select randomly (default: TRUE)
#'
#' @return TSENATAnalysis object created from vignette data subset
#'
#' @details
#' This function loads the built-in readcounts dataset (SALMON transcript-level
#' estimates) along with associated metadata and annotation. It then:
#' 1. Builds a complete TSENATAnalysis object from the full dataset
#' 2. Selects a representative subset of genes (optionally by variance)
#' 3. Selects samples across both conditions
#' 4. Returns the subsetted analysis for use in documentation
#'
#' The vignette data contains 3089 transcripts across 16 samples (8 control,
#' 8 treatment) from a paired experimental design. This function extracts
#' a small, representative subset suitable for fast documentation examples.
#'
#' @examples
#' # Load default subset (10 genes, 8 samples)
#' vignette_subset <- .get_vignette_data_subset()
#' 
#' # Load larger subset (50 genes, 12 samples) with random selection
#' vignette_subset <- .get_vignette_data_subset(
#'   n_genes = 50,
#'   n_samples = 12,
#'   select_top_var = FALSE
#' )
#'
#' @keywords internal
#' @noRd
.get_vignette_data_subset <- function(
    n_genes = 10,
    n_samples = 8,
    seed = 42,
    select_top_var = TRUE) {
  
  set.seed(seed)
  
  # Load vignette data
  data("readcounts", package = "TSENAT", envir = environment())
  
  # Ensure numeric format
  readcounts <- as.matrix(salmon_dataset)
  mode(readcounts) <- "numeric"
  
  # Get associated metadata
  salmon_tpm_data <- salmon_tpm  # noqa: object_name_linter
  
  # Load sample metadata
  metadata_df <- read.table(
    system.file("extdata", "metadata.tsv", package = "TSENAT"),
    header = TRUE, sep = "\t"
  )
  
  # Get GFF3 annotation
  gff3_dataset_path <- system.file("extdata", "annotation.gff3.gz", package = "TSENAT")
  
  # Build full analysis from vignette data
  analysis_full <- build_analysis_s4(
    readcounts,
    gff3_dataset_path,
    metadata = metadata_df,
    tpm = salmon_tpm_data,
    effective_length = salmon_effective_length  # noqa: object_name_linter
  )
  
  # Select genes: by variance or random
  if (select_top_var) {
    # Compute variance per gene (aggregate transcripts to genes)
    se_full <- se(analysis_full)
    counts_matrix <- assay(se_full, "counts")
    gene_ids <- rowData(se_full)$gene_id
    
    # Calculate variance by gene
    gene_vars <- tapply(
      seq_len(nrow(counts_matrix)),
      gene_ids,
      function(idx) {
        var(rowMeans(counts_matrix[idx, , drop = FALSE]))
      }
    )
    
    # Select top n_genes by variance
    top_genes <- names(sort(gene_vars, decreasing = TRUE)[1:min(n_genes, length(gene_vars))])
    gene_idx <- which(gene_ids %in% top_genes)
  } else {
    # Random subset of genes
    all_genes <- unique(rowData(se(analysis_full))$gene_id)
    selected_genes <- sample(all_genes, min(n_genes, length(all_genes)))
    gene_idx <- which(rowData(se(analysis_full))$gene_id %in% selected_genes)
  }
  
  # Select balanced sample set across conditions
  se_full <- se(analysis_full)
  coldata <- colData(se_full)
  
  if ("sample_type" %in% colnames(coldata)) {
    # Try to balance by condition if available
    cond_levels <- unique(coldata$sample_type)
    samples_per_cond <- ceiling(n_samples / length(cond_levels))
    
    sample_idx <- c()
    for (cond in cond_levels) {
      cond_idx <- which(coldata$sample_type == cond)
      selected <- sample(cond_idx, min(samples_per_cond, length(cond_idx)))
      sample_idx <- c(sample_idx, selected)
    }
    sample_idx <- sample_idx[1:min(n_samples, length(sample_idx))]
  } else {
    # Random sample selection
    sample_idx <- sample(ncol(se_full), min(n_samples, ncol(se_full)))
  }
  
  # Return subsetted analysis using the new [ method
  analysis_subset <- analysis_full[gene_idx, sample_idx]
  
  # Update metadata to sync with subsetted genes
  if (!is.null(S4Vectors::metadata(analysis_full@se)$tx2gene)) {
    # Subset tx2gene mapping to only include selected transcripts
    tx2gene_full <- S4Vectors::metadata(analysis_full@se)$tx2gene
    tx2gene_subset <- tx2gene_full[tx2gene_full$Transcript %in% rownames(analysis_subset@se), ]
    S4Vectors::metadata(analysis_subset@se)$tx2gene <- tx2gene_subset
  }
  
  return(analysis_subset)
}
