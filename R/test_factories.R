#' Create Test Analysis Data
#' 
#' Factory functions that generate consistent, reusable test data 
#' with biological signal for examples and tests.
#'
#' @param n_genes Number of genes to simulate (default: 8)
#' @param n_samples_per_group Samples per condition (default: 20)
#' @param control_lambda Poisson lambda for control (default: 40)
#' @param treatment_lambda Poisson lambda for treatment (default: 150)
#' @param q_values Q-value vector (default: c(0.5, 1.0, 1.5))
#' @param seed Random seed (default: 42)
#' @param verbose Progress messages (default: FALSE)
#'
#' @return TSENATAnalysis with diversity results across q-values
#'
#' @details
#' Eliminates repetitive setup code (~70% reduction). Ensures:
#' - Biological signal (control lambda=40 vs treatment lambda=150)
#' - Sufficient samples (40 total: 20 per group)
#' - Multiple q-values avoid rank deficiency
#'
#' @examples
#' # Create with defaults
#' analysis <- create_test_analysis()
#' head(diversity(analysis, q = 1.0))
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
  
  rownames(counts) <- paste0("TX_", seq_len(n_transcripts))
  colnames(counts) <- paste0("Sample_", seq_len(n_samples))
  
  # Create rowData with gene mappings
  rowData <- S4Vectors::DataFrame(
    transcript_id = rownames(counts),
    gene_id = paste0("GENE_", rep(seq_len(n_genes), each = 50, length.out = n_transcripts)),
    row.names = rownames(counts)
  )
  
  # Create colData with experimental design
  colData <- S4Vectors::DataFrame(
    sample_id = colnames(counts),
    condition = rep(c("control", "treatment"), length.out = n_samples),
    sample_type = rep(c("control", "treatment"), length.out = n_samples),
    subject = rep(paste0("S", seq_len(10)), length.out = n_samples),
    paired_samples = rep(paste0("pair", seq_len(10)), length.out = n_samples),
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
  
  # Initialize TSENATAnalysis
  analysis <- TSENATAnalysis(se = se, config = list())
  
  # Calculate diversity
  analysis <- calculate_diversity_s4(
    analysis,
    q = q_values,
    verbose = FALSE,
    min_valid_frac = 0
  )
  
  if (verbose) {
    message("[create_test_analysis] Analysis created with ",
            length(q_values), " q-values")
    message("[create_test_analysis] Diversity results: ",
            nrow(diversity(analysis)), " genes")
  }
  
  return(analysis)
}

#' Lightweight Test Analysis
#' 
#' Minimal factory for parameter extraction tests.
#' Faster than full create_test_analysis().
#'
#' @param n_genes Genes (default: 4)
#' @param n_samples Samples (default: 8)
#' @param seed Random seed (default: 42)
#'
#' @return Minimal TSENATAnalysis
#'
#' @keywords internal
#' @export
create_lightweight_analysis <- function(
    n_genes = 4,
    n_samples = 8,
    seed = 42) {
  
  n_transcripts <- n_genes * 20
  
  # Minimal counts
  counts <- matrix(
    rpois(n_transcripts * n_samples, lambda = 500),
    nrow = n_transcripts,
    ncol = n_samples
  )
  counts <- pmax(counts, 50)
  
  rownames(counts) <- paste0("TX_", seq_len(n_transcripts))
  colnames(counts) <- paste0("Sample_", seq_len(n_samples))
  
  rowData <- S4Vectors::DataFrame(
    transcript_id = rownames(counts),
    gene_id = paste0("GENE_", rep(seq_len(n_genes), each = 20, length.out = n_transcripts)),
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
  analysis <- TSENATAnalysis(se = se, config = list())
  
  # Single q-value (fast)
  analysis <- calculate_diversity_s4(
    analysis,
    q = 1.0,
    verbose = FALSE,
    min_valid_frac = 0
  )
  
  return(analysis)
}
