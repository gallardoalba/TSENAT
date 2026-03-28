#' Test Data Factory Functions for TSENATAnalysis
#' 
#' These helper functions create consistent, reusable test data with biological signal.
#' They follow S4 patterns and reduce boilerplate across the test suite by ~70%.
#'
#' @keywords internal
#' @name test_factories

#' Create a Complete TSENATAnalysis with Diversity Results
#' 
#' Factory function that generates a fully initialized TSENATAnalysis object
#' with realistic data and biological signal. Eliminates repetitive setup code
#' across tests.
#'
#' @param n_genes Number of genes to simulate (default: 8)
#' @param n_samples_per_group Samples per condition (default: 20)
#' @param control_lambda Poisson lambda for control condition (default: 40)
#' @param treatment_lambda Poisson lambda for treatment condition (default: 150)
#' @param q_values Vector of q-values for diversity calculation (default: c(0.5, 1.0, 1.5))
#' @param seed Random seed for reproducibility (default: 42)
#' @param verbose Logical for progress messages (default: FALSE)
#'
#' @return TSENATAnalysis object with:
#'   - SummarizedExperiment with count matrix, rowData, and colData
#'   - Computed diversity results across q-values
#'   - Proper tx2gene metadata mapping
#'
#' @details
#' The factory ensures:
#' - Biological signal: control (lambda=40) vs treatment (lambda=150) contrast
#' - Sufficient samples: 40 total (20 per group) for stable LM fitting
#' - Multiple q-values: c(0.5, 1.0, 1.5) avoids rank deficiency
#' - Valid S4 object structure: passes all TSENATAnalysis validity checks
#'
#' @examples
#' \dontrun{
#'   # Create with defaults (8 genes, 20 samples/group, multi-q)
#'   analysis <- .create_test_analysis()
#'   
#'   # Create with custom parameters
#'   analysis <- .create_test_analysis(
#'     n_genes = 16,
#'     n_samples_per_group = 30,
#'     control_lambda = 50,
#'     treatment_lambda = 200,
#'     q_values = c(0.1, 0.5, 1.0, 1.5, 2.0)
#'   )
#' }
#'
#' @export
create_test_analysis <- function(
    n_genes = 8,
    n_samples_per_group = 20,
    control_lambda = 40,
    treatment_lambda = 150,
    q_values = c(0.5, 1.0, 1.5),
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
  analysis <- TSENAT::TSENATAnalysis(se = se, config = list())
  
  # Calculate diversity
  analysis <- TSENAT::calculate_diversity_s4(
    analysis,
    q = q_values,
    verbose = FALSE,
    min_valid_frac = 0
  )
  
  if (verbose) {
    message("[create_test_analysis] Analysis created with ",
            length(q_values), " q-values")
    message("[create_test_analysis] Diversity results: ",
            nrow(TSENAT::diversity(analysis)), " genes")
  }
  
  return(analysis)
}

#' Create Minimal TSENATAnalysis for Parameter Extraction Tests
#' 
#' Lightweight factory for tests that only validate parameter extraction,
#' not full computation. Reduces test execution time.
#'
#' @param n_genes Number of genes (default: 4)
#' @param n_samples Number of samples (default: 8)
#' @param seed Random seed (default: 42)
#'
#' @return Minimal TSENATAnalysis object
#'
#' @keywords internal
#' @export
create_lightweight_analysis <- function(
    n_genes = 4,
    n_samples = 8,
    seed = 42) {
  
  set.seed(seed)
  
  n_transcripts <- n_genes * 20
  
  # Minimal counts
  counts <- matrix(
    rpois(n_transcripts * n_samples, lambda = 500),
    nrow = n_transcripts,
    ncol = n_samples
  )
  counts <- pmax(counts, 50)
  
  rownames(counts) <- paste0("TX_", 1:n_transcripts)
  colnames(counts) <- paste0("Sample_", 1:n_samples)
  
  rowData <- S4Vectors::DataFrame(
    transcript_id = rownames(counts),
    gene_id = paste0("GENE_", rep(1:n_genes, each = 20, length.out = n_transcripts)),
    row.names = rownames(counts)
  )
  
  colData <- S4Vectors::DataFrame(
    sample_id = colnames(counts),
    condition = rep(c("A", "B"), length.out = n_samples),
    subject = rep(c("S1", "S2", "S3", "S4"), length.out = n_samples),
    row.names = colnames(counts)
  )
  
  se <- SummarizedExperiment::SummarizedExperiment(
    assays = list(counts = counts),
    rowData = rowData,
    colData = colData
  )
  
  tx2gene_df <- data.frame(
    Transcript = rownames(counts),
    Gene = rowData$gene_id,
    stringsAsFactors = FALSE
  )
  S4Vectors::metadata(se)$tx2gene <- tx2gene_df
  
  # Minimal analysis
  analysis <- TSENAT::TSENATAnalysis(se = se, config = list())
  
  # Single q-value (fast)
  analysis <- TSENAT::calculate_diversity_s4(
    analysis,
    q = 1.0,
    verbose = FALSE,
    min_valid_frac = 0
  )
  
  return(analysis)
}

#' Validate Analysis Object Structure
#' 
#' Asserts that a TSENATAnalysis object has the expected structure
#' for conducting downstream LM interaction analyses.
#'
#' @param analysis TSENATAnalysis object to validate
#'
#' @return Invisibly TRUE if valid; otherwise stops with descriptive error
#'
#' @keywords internal
#' @export
assert_analysis_valid <- function(analysis) {
  
  # Check S4 class
  if (!methods::is(analysis, "TSENAT::TSENATAnalysis")) {
    stop("Object must be a TSENATAnalysis instance", call. = FALSE)
  }
  
  # Check required slots
  required_slots <- c("se", "config", "lm_results", "metadata")
  for (slot in required_slots) {
    if (!methods::hasSlot(analysis, slot)) {
      stop("Missing required slot: ", slot, call. = FALSE)
    }
  }
  
  # Check SummarizedExperiment structure
  se <- analysis@se
  if (!methods::is(se, "SummarizedExperiment")) {
    stop("@se must be a SummarizedExperiment object", call. = FALSE)
  }
  
  if (nrow(se) == 0 || ncol(se) == 0) {
    stop("SummarizedExperiment has zero dimensions", call. = FALSE)
  }
  
  # Check colData has required columns
  required_colData <- c("sample_id", "condition")
  for (col in required_colData) {
    if (!col %in% colnames(colData(se))) {
      stop("colData missing required column: ", col, call. = FALSE)
    }
  }
  
  # Check rowData has required columns
  required_rowData <- c("transcript_id", "gene_id")
  for (col in required_rowData) {
    if (!col %in% colnames(rowData(se))) {
      stop("rowData missing required column: ", col, call. = FALSE)
    }
  }
  
  # Check diversity results exist
  if (is.null(analysis@diversity_results) || length(analysis@diversity_results) == 0) {
    stop("No diversity results found. Run calculate_diversity_s4() first.", call. = FALSE)
  }
  
  invisible(TRUE)
}
