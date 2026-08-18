context("S4 Filter Analysis: Numerical Reproducibility and Consistency")

# ===========================================================================
# Helper Functions for Creating Test Data
# ===========================================================================

#' Create a standardized test TSENATAnalysis object
#' @param n_genes Number of genes
#' @param n_samples Number of samples
#' @param n_transcripts_per_gene Transcripts per gene
#' @param seed Random seed for reproducibility
make_test_tsenat_for_filter <- function(
    n_genes = 10,
    n_samples = 6,
    n_transcripts_per_gene = 5,
    seed = 42) {
  
  set.seed(seed)
  
  n_transcripts <- n_genes * n_transcripts_per_gene
  
  # Create count matrix with biological variance
  counts <- matrix(0, nrow = n_transcripts, ncol = n_samples)
  for (j in seq_len(n_samples)) {
    lambda <- if (j > n_samples / 2) 150 else 40
    counts[, j] <- rpois(n_transcripts, lambda = lambda)
  }
  counts <- pmax(counts, 1)  # Ensure minimum non-zero count
  
  rownames(counts) <- paste0("TX_", 1:n_transcripts)
  colnames(counts) <- paste0("Sample_", 1:n_samples)
  
  # Create rowData with gene mappings
  rowData <- S4Vectors::DataFrame(
    transcript_id = rownames(counts),
    gene_id = paste0("GENE_", rep(1:n_genes, each = n_transcripts_per_gene)),
    row.names = rownames(counts)
  )
  
  # Create colData
  colData <- S4Vectors::DataFrame(
    sample_id = colnames(counts),
    condition = rep(c("control", "treatment"), length.out = n_samples),
    pair_id = rep(paste0("pair_", 1:(n_samples / 2)), length.out = n_samples),
    row.names = colnames(counts)
  )
  
  # Create SummarizedExperiment
  se <- SummarizedExperiment::SummarizedExperiment(
    assays = list(counts = counts),
    rowData = rowData,
    colData = colData
  )
  
  # Add TPM data to metadata
  tpm <- counts
  for (j in seq_len(ncol(tpm))) {
    lib_size <- colSums(tpm[, j, drop = FALSE])
    if (lib_size > 0) {
      tpm[, j] <- (tpm[, j] / lib_size) * 1e6
    }
  }
  rownames(tpm) <- rownames(counts)
  colnames(tpm) <- colnames(counts)
  S4Vectors::metadata(se)$tpm <- tpm
  
  # Add tx2gene metadata
  tx2gene_df <- data.frame(
    Transcript = rownames(counts),
    Gene = rowData$gene_id,
    stringsAsFactors = FALSE
  )
  S4Vectors::metadata(se)$tx2gene <- tx2gene_df
  
  # Create TSENATAnalysis
  analysis <- TSENAT::TSENATAnalysis(se = se, config = list())
  
  return(analysis)
}

# ===========================================================================
# TEST GROUP 1: Filtering Reproducibility
# ===========================================================================

test_that("filter_analysis produces identical results with identical parameters", {
  
  analysis1 <- make_test_tsenat_for_filter(n_genes = 8, n_samples = 6, seed = 100)
  analysis2 <- make_test_tsenat_for_filter(n_genes = 8, n_samples = 6, seed = 100)
  
  # Apply identical filtering
  filtered1 <- TSENAT::filter_analysis(
    analysis1,
    min_tpm = 1,
    min_samples = 2,
    verbose = FALSE
  )
  
  filtered2 <- TSENAT::filter_analysis(
    analysis2,
    min_tpm = 1,
    min_samples = 2,
    verbose = FALSE
  )
  
  # Check dimensions match
  expect_equal(nrow(filtered1@se), nrow(filtered2@se))
  expect_equal(ncol(filtered1@se), ncol(filtered2@se))
  
  # Check assay values match exactly
  assay1 <- SummarizedExperiment::assay(filtered1@se, "counts")
  assay2 <- SummarizedExperiment::assay(filtered2@se, "counts")
  expect_equal(assay1, assay2)
  
  # Check rownames and colnames match
  expect_equal(rownames(filtered1@se), rownames(filtered2@se))
  expect_equal(colnames(filtered1@se), colnames(filtered2@se))
})

test_that("filter_analysis is deterministic across multiple runs", {
  
  analysis <- make_test_tsenat_for_filter(n_genes = 12, n_samples = 8, seed = 200)
  
  # Run filtering three times with identical parameters
  filtered_runs <- lapply(1:3, function(i) {
    TSENAT::filter_analysis(
      analysis,
      min_tpm = 2,
      min_samples = 3,
      verbose = FALSE
    )
  })
  
  # All runs should produce identical results
  assays_run1 <- SummarizedExperiment::assay(filtered_runs[[1]]@se, "counts")
  assays_run2 <- SummarizedExperiment::assay(filtered_runs[[2]]@se, "counts")
  assays_run3 <- SummarizedExperiment::assay(filtered_runs[[3]]@se, "counts")
  
  expect_equal(assays_run1, assays_run2)
  expect_equal(assays_run2, assays_run3)
  
  # Dimensions should match
  expect_equal(nrow(filtered_runs[[1]]@se), nrow(filtered_runs[[2]]@se))
  expect_equal(nrow(filtered_runs[[2]]@se), nrow(filtered_runs[[3]]@se))
})

test_that("filter_analysis filtering masks are consistent", {
  
  analysis1 <- make_test_tsenat_for_filter(n_genes = 10, n_samples = 6, seed = 300)
  analysis2 <- make_test_tsenat_for_filter(n_genes = 10, n_samples = 6, seed = 300)
  
  # Filter with specific stringency
  filtered1 <- TSENAT::filter_analysis(
    analysis1,
    min_tpm = 5,
    min_samples = 2,
    verbose = FALSE
  )
  
  filtered2 <- TSENAT::filter_analysis(
    analysis2,
    min_tpm = 5,
    min_samples = 2,
    verbose = FALSE
  )
  
  # The filtered row counts should be identical (same filtering mask applied)
  expect_equal(nrow(filtered1@se), nrow(filtered2@se),
               info = "Number of transcripts after filtering should match")
})

# ===========================================================================
# TEST GROUP 2: Stringency-Based Filtering Reproducibility
# ===========================================================================

test_that("filter_analysis with stringency='soft' is reproducible", {
  
  analysis1 <- make_test_tsenat_for_filter(n_genes = 8, n_samples = 8, seed = 400)
  analysis2 <- make_test_tsenat_for_filter(n_genes = 8, n_samples = 8, seed = 400)
  
  filtered1 <- TSENAT::filter_analysis(
    analysis1,
    stringency = "soft",
    pair_col = "pair_id",
    verbose = FALSE
  )
  
  filtered2 <- TSENAT::filter_analysis(
    analysis2,
    stringency = "soft",
    pair_col = "pair_id",
    verbose = FALSE
  )
  
  # Should produce identical results
  expect_equal(nrow(filtered1@se), nrow(filtered2@se))
  expect_equal(
    SummarizedExperiment::assay(filtered1@se, "counts"),
    SummarizedExperiment::assay(filtered2@se, "counts")
  )
})

test_that("filter_analysis with stringency='medium' is reproducible", {
  
  analysis1 <- make_test_tsenat_for_filter(n_genes = 8, n_samples = 8, seed = 500)
  analysis2 <- make_test_tsenat_for_filter(n_genes = 8, n_samples = 8, seed = 500)
  
  filtered1 <- TSENAT::filter_analysis(
    analysis1,
    stringency = "medium",
    pair_col = "pair_id",
    verbose = FALSE
  )
  
  filtered2 <- TSENAT::filter_analysis(
    analysis2,
    stringency = "medium",
    pair_col = "pair_id",
    verbose = FALSE
  )
  
  expect_equal(nrow(filtered1@se), nrow(filtered2@se))
  expect_equal(
    SummarizedExperiment::assay(filtered1@se, "counts"),
    SummarizedExperiment::assay(filtered2@se, "counts")
  )
})

test_that("filter_analysis with stringency='severe' is reproducible", {
  
  analysis1 <- make_test_tsenat_for_filter(n_genes = 8, n_samples = 8, seed = 600)
  analysis2 <- make_test_tsenat_for_filter(n_genes = 8, n_samples = 8, seed = 600)
  
  filtered1 <- TSENAT::filter_analysis(
    analysis1,
    stringency = "severe",
    pair_col = "pair_id",
    verbose = FALSE
  )
  
  filtered2 <- TSENAT::filter_analysis(
    analysis2,
    stringency = "severe",
    pair_col = "pair_id",
    verbose = FALSE
  )
  
  expect_equal(nrow(filtered1@se), nrow(filtered2@se))
  expect_equal(
    SummarizedExperiment::assay(filtered1@se, "counts"),
    SummarizedExperiment::assay(filtered2@se, "counts")
  )
})

# ===========================================================================
# TEST GROUP 3: Subsetting Reproducibility (Seed-Based)
# ===========================================================================

test_that("filter_analysis with subset_n_genes and seed is reproducible", {
  
  analysis1 <- make_test_tsenat_for_filter(n_genes = 20, n_samples = 6, seed = 700)
  analysis2 <- make_test_tsenat_for_filter(n_genes = 20, n_samples = 6, seed = 700)
  
  # Filter and subset with identical seed
  filtered1 <- TSENAT::filter_analysis(
    analysis1,
    min_tpm = 1,
    min_samples = 2,
    subset_n_genes = 10,
    subset_seed = 888,
    verbose = FALSE
  )
  
  filtered2 <- TSENAT::filter_analysis(
    analysis2,
    min_tpm = 1,
    min_samples = 2,
    subset_n_genes = 10,
    subset_seed = 888,
    verbose = FALSE
  )
  
  # Same genes should be selected
  expect_equal(nrow(filtered1@se), nrow(filtered2@se))
  
  # Assays should be identical
  expect_equal(
    SummarizedExperiment::assay(filtered1@se, "counts"),
    SummarizedExperiment::assay(filtered2@se, "counts")
  )
})

test_that("filter_analysis random subsetting with same seed produces same genes", {
  
  analysis <- make_test_tsenat_for_filter(n_genes = 30, n_samples = 6, seed = 800)
  
  # Run random subsetting twice with same seed
  filtered1 <- TSENAT::filter_analysis(
    analysis,
    min_tpm = 1,
    min_samples = 1,
    subset_n_genes = 12,
    subset_select_by = "random",
    subset_seed = 999,
    verbose = FALSE
  )
  
  # Make a fresh copy for second run
  analysis_copy <- make_test_tsenat_for_filter(n_genes = 30, n_samples = 6, seed = 800)
  
  filtered2 <- TSENAT::filter_analysis(
    analysis_copy,
    min_tpm = 1,
    min_samples = 1,
    subset_n_genes = 12,
    subset_select_by = "random",
    subset_seed = 999,
    verbose = FALSE
  )
  
  # Both runs should select same number of transcripts
  expect_equal(nrow(filtered1@se), nrow(filtered2@se))
  
  # Same transcripts should be selected with identical seed
  expect_equal(rownames(filtered1@se), rownames(filtered2@se))
})

test_that("filter_analysis variance-based subsetting is reproducible", {
  
  analysis1 <- make_test_tsenat_for_filter(n_genes = 24, n_samples = 6, seed = 900)
  analysis2 <- make_test_tsenat_for_filter(n_genes = 24, n_samples = 6, seed = 900)
  
  # Variance-based subsetting doesn't require seed but should be deterministic
  filtered1 <- TSENAT::filter_analysis(
    analysis1,
    min_tpm = 1,
    min_samples = 1,
    subset_n_genes = 15,
    subset_select_by = "variance",
    verbose = FALSE
  )
  
  filtered2 <- TSENAT::filter_analysis(
    analysis2,
    min_tpm = 1,
    min_samples = 1,
    subset_n_genes = 15,
    subset_select_by = "variance",
    verbose = FALSE
  )
  
  # Should select same high-variance genes
  expect_equal(nrow(filtered1@se), nrow(filtered2@se))
  expect_equal(rownames(filtered1@se), rownames(filtered2@se))
})

# ===========================================================================
# TEST GROUP 4: Combined Filtering + Subsetting Reproducibility
# ===========================================================================

test_that("filter_analysis combined filtering and subsetting is reproducible", {
  
  analysis1 <- make_test_tsenat_for_filter(n_genes = 25, n_samples = 8, seed = 1000)
  analysis2 <- make_test_tsenat_for_filter(n_genes = 25, n_samples = 8, seed = 1000)
  
  # Apply both filtering and subsetting
  filtered1 <- TSENAT::filter_analysis(
    analysis1,
    min_tpm = 3,
    min_samples = 2,
    subset_n_genes = 10,
    subset_select_by = "random",
    subset_seed = 777,
    verbose = FALSE
  )
  
  filtered2 <- TSENAT::filter_analysis(
    analysis2,
    min_tpm = 3,
    min_samples = 2,
    subset_n_genes = 10,
    subset_select_by = "random",
    subset_seed = 777,
    verbose = FALSE
  )
  
  # Results should be identical
  expect_equal(nrow(filtered1@se), nrow(filtered2@se))
  expect_equal(ncol(filtered1@se), ncol(filtered2@se))
  expect_equal(
    SummarizedExperiment::assay(filtered1@se, "counts"),
    SummarizedExperiment::assay(filtered2@se, "counts")
  )
})

test_that("filter_analysis different seeds produce different subsets", {
  
  analysis <- make_test_tsenat_for_filter(n_genes = 30, n_samples = 6, seed = 1100)
  
  # Apply subsetting with two different seeds
  filtered_seed1 <- TSENAT::filter_analysis(
    analysis,
    min_tpm = 1,
    min_samples = 1,
    subset_n_genes = 10,
    subset_select_by = "random",
    subset_seed = 111,
    verbose = FALSE
  )
  
  # Make fresh copy
  analysis_copy <- make_test_tsenat_for_filter(n_genes = 30, n_samples = 6, seed = 1100)
  
  filtered_seed2 <- TSENAT::filter_analysis(
    analysis_copy,
    min_tpm = 1,
    min_samples = 1,
    subset_n_genes = 10,
    subset_select_by = "random",
    subset_seed = 222,
    verbose = FALSE
  )
  
  # Different seeds should produce different gene selections (with very high probability)
  expect_false(identical(
    rownames(filtered_seed1@se),
    rownames(filtered_seed2@se)
  ))
})

# ===========================================================================
# TEST GROUP 5: Data Integrity After Filtering
# ===========================================================================

test_that("filter_analysis preserves assay values without modification", {
  
  analysis <- make_test_tsenat_for_filter(n_genes = 12, n_samples = 6, seed = 1200)
  
  # Store original counts
  original_counts <- SummarizedExperiment::assay(analysis@se, "counts")
  
  # Apply filtering
  filtered <- TSENAT::filter_analysis(
    analysis,
    min_tpm = 2,
    min_samples = 2,
    verbose = FALSE
  )
  
  filtered_counts <- SummarizedExperiment::assay(filtered@se, "counts")
  
  # Check that filtered counts are a proper subset of original counts
  # (same values, just subset of rows)
  for (row in rownames(filtered_counts)) {
    expect_true(row %in% rownames(original_counts),
                info = sprintf("Row %s should exist in original", row))
    expect_equal(
      filtered_counts[row, ],
      original_counts[row, ],
      info = sprintf("Row %s values should be unchanged", row)
    )
  }
})

test_that("filter_analysis maintains sample order and metadata", {
  
  analysis <- make_test_tsenat_for_filter(n_genes = 10, n_samples = 8, seed = 1300)
  
  # Store original sample info
  original_samples <- colnames(analysis@se)
  original_coldata <- SummarizedExperiment::colData(analysis@se)
  
  # Apply filtering
  filtered <- TSENAT::filter_analysis(
    analysis,
    min_tpm = 1,
    min_samples = 2,
    verbose = FALSE
  )
  
  filtered_samples <- colnames(filtered@se)
  filtered_coldata <- SummarizedExperiment::colData(filtered@se)
  
  # Samples should be in same order
  expect_equal(filtered_samples, original_samples)
  
  # Column metadata should be identical
  expect_equal(
    as.data.frame(filtered_coldata),
    as.data.frame(original_coldata)
  )
})

# ===========================================================================
# TEST GROUP 6: Numerical Consistency with Extreme Parameters
# ===========================================================================

test_that("filter_analysis with min_tpm=0 is reproducible", {
  
  analysis1 <- make_test_tsenat_for_filter(n_genes = 8, n_samples = 6, seed = 1400)
  analysis2 <- make_test_tsenat_for_filter(n_genes = 8, n_samples = 6, seed = 1400)
  
  filtered1 <- TSENAT::filter_analysis(
    analysis1,
    min_tpm = 0,
    min_samples = 1,
    verbose = FALSE
  )
  
  filtered2 <- TSENAT::filter_analysis(
    analysis2,
    min_tpm = 0,
    min_samples = 1,
    verbose = FALSE
  )
  
  # Should have identical results (minimal filtering)
  expect_equal(nrow(filtered1@se), nrow(filtered2@se))
  expect_equal(
    SummarizedExperiment::assay(filtered1@se, "counts"),
    SummarizedExperiment::assay(filtered2@se, "counts")
  )
})

test_that("filter_analysis with strict filtering is reproducible", {
  
  analysis1 <- make_test_tsenat_for_filter(n_genes = 12, n_samples = 6, seed = 1500)
  analysis2 <- make_test_tsenat_for_filter(n_genes = 12, n_samples = 6, seed = 1500)
  
  # Very strict filtering
  filtered1 <- TSENAT::filter_analysis(
    analysis1,
    min_tpm = 100,
    min_samples = 5,
    verbose = FALSE
  )
  
  filtered2 <- TSENAT::filter_analysis(
    analysis2,
    min_tpm = 100,
    min_samples = 5,
    verbose = FALSE
  )
  
  # Should have identical results
  expect_equal(nrow(filtered1@se), nrow(filtered2@se))
  expect_equal(ncol(filtered1@se), ncol(filtered2@se))
})

# ===========================================================================
# TEST GROUP 7: Subsetting Selection Methods Reproducibility
# ===========================================================================

test_that("filter_analysis mean-based subsetting is reproducible", {
  
  analysis1 <- make_test_tsenat_for_filter(n_genes = 24, n_samples = 6, seed = 1600)
  analysis2 <- make_test_tsenat_for_filter(n_genes = 24, n_samples = 6, seed = 1600)
  
  # Mean-based subsetting is deterministic
  filtered1 <- TSENAT::filter_analysis(
    analysis1,
    min_tpm = 1,
    min_samples = 1,
    subset_n_genes = 14,
    subset_select_by = "mean",
    verbose = FALSE
  )
  
  filtered2 <- TSENAT::filter_analysis(
    analysis2,
    min_tpm = 1,
    min_samples = 1,
    subset_n_genes = 14,
    subset_select_by = "mean",
    verbose = FALSE
  )
  
  # Should select same high-mean genes
  expect_equal(rownames(filtered1@se), rownames(filtered2@se))
  expect_equal(nrow(filtered1@se), nrow(filtered2@se))
})

# ===========================================================================
# TEST GROUP 8: Isoform-Level Filtering Reproducibility
# ===========================================================================

test_that("filter_analysis with min_isoform_abundance is reproducible", {
  
  analysis1 <- make_test_tsenat_for_filter(n_genes = 10, n_samples = 6, seed = 1700)
  analysis2 <- make_test_tsenat_for_filter(n_genes = 10, n_samples = 6, seed = 1700)
  
  filtered1 <- TSENAT::filter_analysis(
    analysis1,
    min_tpm = 1,
    min_samples = 1,
    min_isoform_abundance = 0.05,
    verbose = FALSE
  )
  
  filtered2 <- TSENAT::filter_analysis(
    analysis2,
    min_tpm = 1,
    min_samples = 1,
    min_isoform_abundance = 0.05,
    verbose = FALSE
  )
  
  # Should produce identical results
  expect_equal(nrow(filtered1@se), nrow(filtered2@se))
  expect_equal(
    SummarizedExperiment::assay(filtered1@se, "counts"),
    SummarizedExperiment::assay(filtered2@se, "counts")
  )
})

test_that("filter_analysis with min_tx_per_gene is reproducible", {
  
  analysis1 <- make_test_tsenat_for_filter(n_genes = 8, n_samples = 6, seed = 1800)
  analysis2 <- make_test_tsenat_for_filter(n_genes = 8, n_samples = 6, seed = 1800)
  
  filtered1 <- TSENAT::filter_analysis(
    analysis1,
    min_tpm = 1,
    min_samples = 1,
    min_tx_per_gene = 3,
    verbose = FALSE
  )
  
  filtered2 <- TSENAT::filter_analysis(
    analysis2,
    min_tpm = 1,
    min_samples = 1,
    min_tx_per_gene = 3,
    verbose = FALSE
  )
  
  expect_equal(nrow(filtered1@se), nrow(filtered2@se))
  expect_equal(
    SummarizedExperiment::assay(filtered1@se, "counts"),
    SummarizedExperiment::assay(filtered2@se, "counts")
  )
})

# ===========================================================================
# TEST GROUP 9: Filtering Output Structure Consistency
# ===========================================================================

test_that("filter_analysis output has valid S4 structure after filtering", {
  
  analysis <- make_test_tsenat_for_filter(n_genes = 10, n_samples = 6, seed = 1900)
  
  filtered <- TSENAT::filter_analysis(
    analysis,
    min_tpm = 2,
    min_samples = 2,
    verbose = FALSE
  )
  
  # Check S4 class
  expect_s4_class(filtered, "TSENATAnalysis")
  
  # Check slots exist
  slot_names <- slotNames(filtered)
  expect_true("se" %in% slot_names)
  expect_true("config" %in% slot_names)
  
  # Check SE is valid
  expect_s4_class(filtered@se, "SummarizedExperiment")
  
  # Check assay exists
  expect_true(length(SummarizedExperiment::assays(filtered@se)) > 0)
})

test_that("filter_analysis maintains config across filtering", {
  
  analysis <- make_test_tsenat_for_filter(n_genes = 8, n_samples = 6, seed = 2000)
  analysis@config <- list(param1 = "value1", param2 = 42)
  
  filtered <- TSENAT::filter_analysis(
    analysis,
    min_tpm = 1,
    min_samples = 2,
    verbose = FALSE
  )
  
  # Config should be preserved
  expect_equal(filtered@config$param1, "value1")
  expect_equal(filtered@config$param2, 42)
})

# ===========================================================================
# TEST GROUP 10: Multiple min_tpm/min_samples Combinations
# ===========================================================================

test_that("filter_analysis with various min_tpm thresholds is reproducible", {
  
  test_tpm_values <- c(0.5, 1, 2, 5, 10)
  results_first_run <- list()
  results_second_run <- list()
  
  for (tpm_val in test_tpm_values) {
    analysis1 <- make_test_tsenat_for_filter(n_genes = 12, n_samples = 6, seed = 2100)
    analysis2 <- make_test_tsenat_for_filter(n_genes = 12, n_samples = 6, seed = 2100)
    
    filtered1 <- TSENAT::filter_analysis(
      analysis1,
      min_tpm = tpm_val,
      min_samples = 2,
      verbose = FALSE
    )
    
    filtered2 <- TSENAT::filter_analysis(
      analysis2,
      min_tpm = tpm_val,
      min_samples = 2,
      verbose = FALSE
    )
    
    # Each tpm threshold should produce consistent results
    expect_equal(nrow(filtered1@se), nrow(filtered2@se),
                 info = sprintf("min_tpm=%.1f should be reproducible", tpm_val))
    
    results_first_run[[as.character(tpm_val)]] <- nrow(filtered1@se)
    results_second_run[[as.character(tpm_val)]] <- nrow(filtered2@se)
  }
  
  # Verify results are stored and consistent
  expect_true(length(results_first_run) > 0)
})

test_that("filter_analysis with various min_samples values is reproducible", {
  
  test_samples_values <- c(1L, 2L, 3L, 4L, 5L)
  
  for (min_samp in test_samples_values) {
    analysis1 <- make_test_tsenat_for_filter(n_genes = 10, n_samples = 6, seed = 2200)
    analysis2 <- make_test_tsenat_for_filter(n_genes = 10, n_samples = 6, seed = 2200)
    
    filtered1 <- TSENAT::filter_analysis(
      analysis1,
      min_tpm = 1,
      min_samples = min_samp,
      verbose = FALSE
    )
    
    filtered2 <- TSENAT::filter_analysis(
      analysis2,
      min_tpm = 1,
      min_samples = min_samp,
      verbose = FALSE
    )
    
    # Each min_samples threshold should produce consistent results
    expect_equal(nrow(filtered1@se), nrow(filtered2@se),
                 info = sprintf("min_samples=%d should be reproducible", min_samp))
  }
})

test_that("filter_analysis with combined min_tpm and min_samples variations is reproducible", {
  
  param_combinations <- list(
    list(tpm = 0.5, samples = 1L),
    list(tpm = 1, samples = 2L),
    list(tpm = 2, samples = 3L),
    list(tpm = 5, samples = 4L),
    list(tpm = 10, samples = 5L)
  )
  
  for (combo in param_combinations) {
    analysis1 <- make_test_tsenat_for_filter(n_genes = 15, n_samples = 8, seed = 2300)
    analysis2 <- make_test_tsenat_for_filter(n_genes = 15, n_samples = 8, seed = 2300)
    
    filtered1 <- TSENAT::filter_analysis(
      analysis1,
      min_tpm = combo$tpm,
      min_samples = combo$samples,
      verbose = FALSE
    )
    
    filtered2 <- TSENAT::filter_analysis(
      analysis2,
      min_tpm = combo$tpm,
      min_samples = combo$samples,
      verbose = FALSE
    )
    
    expect_equal(nrow(filtered1@se), nrow(filtered2@se),
                 info = sprintf("tpm=%.1f, samples=%d should be reproducible", combo$tpm, combo$samples))
  }
})

# ===========================================================================
# TEST GROUP 11: min_tx_per_gene Combinations
# ===========================================================================

test_that("filter_analysis with various min_tx_per_gene values is reproducible", {
  
  test_tx_values <- c(1L, 2L, 3L, 4L, 5L)
  
  for (min_tx in test_tx_values) {
    analysis1 <- make_test_tsenat_for_filter(n_genes = 10, n_samples = 6, seed = 2400)
    analysis2 <- make_test_tsenat_for_filter(n_genes = 10, n_samples = 6, seed = 2400)
    
    filtered1 <- TSENAT::filter_analysis(
      analysis1,
      min_tpm = 1,
      min_samples = 1,
      min_tx_per_gene = min_tx,
      verbose = FALSE
    )
    
    filtered2 <- TSENAT::filter_analysis(
      analysis2,
      min_tpm = 1,
      min_samples = 1,
      min_tx_per_gene = min_tx,
      verbose = FALSE
    )
    
    expect_equal(nrow(filtered1@se), nrow(filtered2@se),
                 info = sprintf("min_tx_per_gene=%d should be reproducible", min_tx))
  }
})

test_that("filter_analysis with min_tx_per_gene and other filters combined is reproducible", {
  
  analysis1 <- make_test_tsenat_for_filter(n_genes = 15, n_samples = 8, seed = 2500)
  analysis2 <- make_test_tsenat_for_filter(n_genes = 15, n_samples = 8, seed = 2500)
  
  filtered1 <- TSENAT::filter_analysis(
    analysis1,
    min_tpm = 2,
    min_samples = 3,
    min_tx_per_gene = 2,
    verbose = FALSE
  )
  
  filtered2 <- TSENAT::filter_analysis(
    analysis2,
    min_tpm = 2,
    min_samples = 3,
    min_tx_per_gene = 2,
    verbose = FALSE
  )
  
  expect_equal(nrow(filtered1@se), nrow(filtered2@se))
  expect_equal(
    SummarizedExperiment::assay(filtered1@se, "counts"),
    SummarizedExperiment::assay(filtered2@se, "counts")
  )
})

# ===========================================================================
# TEST GROUP 12: Sample Subsetting Reproducibility
# ===========================================================================

test_that("filter_analysis with subset_n_samples is reproducible", {
  
  analysis1 <- make_test_tsenat_for_filter(n_genes = 10, n_samples = 10, seed = 2600)
  analysis2 <- make_test_tsenat_for_filter(n_genes = 10, n_samples = 10, seed = 2600)
  
  filtered1 <- TSENAT::filter_analysis(
    analysis1,
    min_tpm = 1,
    min_samples = 1,
    subset_n_samples = 6,
    subset_seed = 555,
    verbose = FALSE
  )
  
  filtered2 <- TSENAT::filter_analysis(
    analysis2,
    min_tpm = 1,
    min_samples = 1,
    subset_n_samples = 6,
    subset_seed = 555,
    verbose = FALSE
  )
  
  # Should have same number of samples
  expect_equal(ncol(filtered1@se), ncol(filtered2@se))
  
  # Same samples should be selected
  expect_equal(colnames(filtered1@se), colnames(filtered2@se))
})

test_that("filter_analysis with subset_n_samples different seeds produce different samples", {
  
  analysis <- make_test_tsenat_for_filter(n_genes = 10, n_samples = 10, seed = 2700)
  
  filtered_seed1 <- TSENAT::filter_analysis(
    analysis,
    min_tpm = 1,
    min_samples = 1,
    subset_n_samples = 6,
    subset_seed = 111,
    verbose = FALSE
  )
  
  analysis_copy <- make_test_tsenat_for_filter(n_genes = 10, n_samples = 10, seed = 2700)
  
  filtered_seed2 <- TSENAT::filter_analysis(
    analysis_copy,
    min_tpm = 1,
    min_samples = 1,
    subset_n_samples = 6,
    subset_seed = 222,
    verbose = FALSE
  )
  
  # Same number of samples but different samples selected
  expect_equal(ncol(filtered_seed1@se), ncol(filtered_seed2@se))
  expect_false(identical(colnames(filtered_seed1@se), colnames(filtered_seed2@se)))
})

# ===========================================================================
# TEST GROUP 13: Combining Genes and Samples Subsetting
# ===========================================================================

test_that("filter_analysis with both subset_n_genes and subset_n_samples is reproducible", {
  
  analysis1 <- make_test_tsenat_for_filter(n_genes = 20, n_samples = 10, seed = 2800)
  analysis2 <- make_test_tsenat_for_filter(n_genes = 20, n_samples = 10, seed = 2800)
  
  filtered1 <- TSENAT::filter_analysis(
    analysis1,
    min_tpm = 1,
    min_samples = 1,
    subset_n_genes = 10,
    subset_n_samples = 6,
    subset_seed = 666,
    verbose = FALSE
  )
  
  filtered2 <- TSENAT::filter_analysis(
    analysis2,
    min_tpm = 1,
    min_samples = 1,
    subset_n_genes = 10,
    subset_n_samples = 6,
    subset_seed = 666,
    verbose = FALSE
  )
  
  # Both dimensions should match
  expect_equal(nrow(filtered1@se), nrow(filtered2@se))
  expect_equal(ncol(filtered1@se), ncol(filtered2@se))
  
  # Content should be identical
  expect_equal(
    SummarizedExperiment::assay(filtered1@se, "counts"),
    SummarizedExperiment::assay(filtered2@se, "counts")
  )
})

# ===========================================================================
# TEST GROUP 14: Isoform Abundance with Other Parameters
# ===========================================================================

test_that("filter_analysis with min_isoform_abundance and other parameters is reproducible", {
  
  analysis1 <- make_test_tsenat_for_filter(n_genes = 12, n_samples = 6, seed = 2900)
  analysis2 <- make_test_tsenat_for_filter(n_genes = 12, n_samples = 6, seed = 2900)
  
  filtered1 <- TSENAT::filter_analysis(
    analysis1,
    min_tpm = 1,
    min_samples = 2,
    min_isoform_abundance = 0.05,
    min_tx_per_gene = 2,
    verbose = FALSE
  )
  
  filtered2 <- TSENAT::filter_analysis(
    analysis2,
    min_tpm = 1,
    min_samples = 2,
    min_isoform_abundance = 0.05,
    min_tx_per_gene = 2,
    verbose = FALSE
  )
  
  expect_equal(nrow(filtered1@se), nrow(filtered2@se))
  expect_equal(
    SummarizedExperiment::assay(filtered1@se, "counts"),
    SummarizedExperiment::assay(filtered2@se, "counts")
  )
})

test_that("filter_analysis with different isoform_abundance thresholds is reproducible", {
  
  test_abund_values <- c(0.01, 0.05, 0.1, 0.2)
  
  for (abund_val in test_abund_values) {
    analysis1 <- make_test_tsenat_for_filter(n_genes = 10, n_samples = 6, seed = 3000)
    analysis2 <- make_test_tsenat_for_filter(n_genes = 10, n_samples = 6, seed = 3000)
    
    filtered1 <- TSENAT::filter_analysis(
      analysis1,
      min_tpm = 1,
      min_samples = 1,
      min_isoform_abundance = abund_val,
      verbose = FALSE
    )
    
    filtered2 <- TSENAT::filter_analysis(
      analysis2,
      min_tpm = 1,
      min_samples = 1,
      min_isoform_abundance = abund_val,
      verbose = FALSE
    )
    
    expect_equal(nrow(filtered1@se), nrow(filtered2@se),
                 info = sprintf("min_isoform_abundance=%.2f should be reproducible", abund_val))
  }
})

# ===========================================================================
# TEST GROUP 15: Complex Multi-Parameter Combinations
# ===========================================================================

test_that("filter_analysis with complex parameter combinations is reproducible", {
  
  complex_params <- list(
    list(min_tpm = 1, min_samples = 2, min_tx_per_gene = 2, min_isoform_abundance = 0.05),
    list(min_tpm = 2, min_samples = 3, min_tx_per_gene = 3, min_isoform_abundance = 0.1),
    list(min_tpm = 0.5, min_samples = 1, min_tx_per_gene = 1, min_isoform_abundance = NULL),
    list(min_tpm = 5, min_samples = 4, min_tx_per_gene = 2, min_isoform_abundance = 0.02)
  )
  
  for (i in seq_along(complex_params)) {
    params <- complex_params[[i]]
    
    analysis1 <- make_test_tsenat_for_filter(n_genes = 15, n_samples = 8, seed = 3100 + i)
    analysis2 <- make_test_tsenat_for_filter(n_genes = 15, n_samples = 8, seed = 3100 + i)
    
    filtered1 <- TSENAT::filter_analysis(
      analysis1,
      min_tpm = params$min_tpm,
      min_samples = params$min_samples,
      min_tx_per_gene = params$min_tx_per_gene,
      min_isoform_abundance = params$min_isoform_abundance,
      verbose = FALSE
    )
    
    filtered2 <- TSENAT::filter_analysis(
      analysis2,
      min_tpm = params$min_tpm,
      min_samples = params$min_samples,
      min_tx_per_gene = params$min_tx_per_gene,
      min_isoform_abundance = params$min_isoform_abundance,
      verbose = FALSE
    )
    
    expect_equal(nrow(filtered1@se), nrow(filtered2@se))
    expect_equal(
      SummarizedExperiment::assay(filtered1@se, "counts"),
      SummarizedExperiment::assay(filtered2@se, "counts")
    )
  }
})

# ===========================================================================
# TEST GROUP 16: Edge Cases and Boundary Conditions
# ===========================================================================

test_that("filter_analysis with very small dataset is reproducible", {
  
  analysis1 <- make_test_tsenat_for_filter(n_genes = 2, n_samples = 3, seed = 3200)
  analysis2 <- make_test_tsenat_for_filter(n_genes = 2, n_samples = 3, seed = 3200)
  
  filtered1 <- TSENAT::filter_analysis(
    analysis1,
    min_tpm = 0.1,
    min_samples = 1,
    verbose = FALSE
  )
  
  filtered2 <- TSENAT::filter_analysis(
    analysis2,
    min_tpm = 0.1,
    min_samples = 1,
    verbose = FALSE
  )
  
  expect_equal(nrow(filtered1@se), nrow(filtered2@se))
})

test_that("filter_analysis with filtering that removes many transcripts is reproducible", {
  
  analysis1 <- make_test_tsenat_for_filter(n_genes = 20, n_samples = 6, seed = 3300)
  analysis2 <- make_test_tsenat_for_filter(n_genes = 20, n_samples = 6, seed = 3300)
  
  # Very strict filtering
  filtered1 <- TSENAT::filter_analysis(
    analysis1,
    min_tpm = 50,
    min_samples = 5,
    verbose = FALSE
  )
  
  filtered2 <- TSENAT::filter_analysis(
    analysis2,
    min_tpm = 50,
    min_samples = 5,
    verbose = FALSE
  )
  
  # Should result in same (possibly small or zero) number of rows
  expect_equal(nrow(filtered1@se), nrow(filtered2@se))
})

test_that("filter_analysis chained subsetting operations is reproducible", {
  
  analysis1 <- make_test_tsenat_for_filter(n_genes = 30, n_samples = 10, seed = 3400)
  analysis2 <- make_test_tsenat_for_filter(n_genes = 30, n_samples = 10, seed = 3400)
  
  # Chain with multiple subsetting operations
  filtered1 <- TSENAT::filter_analysis(
    analysis1,
    min_tpm = 1,
    min_samples = 1,
    subset_n_genes = 20,
    subset_n_samples = 8,
    subset_select_by = "variance",
    subset_seed = 444,
    verbose = FALSE
  )
  
  filtered2 <- TSENAT::filter_analysis(
    analysis2,
    min_tpm = 1,
    min_samples = 1,
    subset_n_genes = 20,
    subset_n_samples = 8,
    subset_select_by = "variance",
    subset_seed = 444,
    verbose = FALSE
  )
  
  expect_equal(nrow(filtered1@se), nrow(filtered2@se))
  expect_equal(ncol(filtered1@se), ncol(filtered2@se))
  expect_equal(rownames(filtered1@se), rownames(filtered2@se))
  expect_equal(colnames(filtered1@se), colnames(filtered2@se))
})


# =============================================================================
# CONTEXT: Filter Analysis Workflow
# =============================================================================

context("Analysis Workflow: Filter Analysis")

# Helper functions
make_test_se <- function(n_genes = 100, n_samples = 20) {
  counts <- matrix(rpois(n_genes * n_samples, lambda = 50), nrow = n_genes)
  rownames(counts) <- paste0("TX_", 1:n_genes)
  colnames(counts) <- paste0("Sample_", 1:n_samples)
  
  # Add TPM data to avoid warnings during filtering
  # Simple TPM calculation: scale counts to sum to 1 million per sample
  tpm <- t(t(counts) / colSums(counts) * 1e6)
  
  SummarizedExperiment::SummarizedExperiment(
    assays = list(counts = counts, tpm = tpm)
  )
}

test_that("filter_analysis modifies SE in analysis object", {
  se <- make_test_se(n_genes = 100, n_samples = 20)
  coldata <- S4Vectors::DataFrame(
    pair_id = rep(1:10, each = 2),
    row.names = colnames(se)
  )
  SummarizedExperiment::colData(se) <- coldata
  analysis <- TSENATAnalysis(se)
  
  # Apply filtering
  filtered_analysis <- filter_analysis(analysis, stringency = "severe", verbose = FALSE)
  
  # Check that analysis is returned
  expect_true(inherits(filtered_analysis, "TSENATAnalysis"))
  
  # Check that SE was modified
  expect_true(nrow(filtered_analysis@se) <= nrow(analysis@se))
})

test_that("filter_analysis preserves colData", {
  se <- make_test_se(n_genes = 100, n_samples = 20)
  coldata <- S4Vectors::DataFrame(
    pair_id = rep(1:10, each = 2),
    condition = rep(c("control", "treatment"), 10),
    row.names = colnames(se)
  )
  SummarizedExperiment::colData(se) <- coldata
  
  analysis <- TSENATAnalysis(se)
  filtered_analysis <- filter_analysis(analysis, stringency = "soft", verbose = FALSE)
  
  # colData should preserve original columns (sample_id is added by constructor)
  filtered_coldata <- SummarizedExperiment::colData(filtered_analysis@se)
  expect_true("pair_id" %in% colnames(filtered_coldata))
  expect_true("condition" %in% colnames(filtered_coldata))
  expect_true("sample_id" %in% colnames(filtered_coldata))
  # Check expected column count: pair_id + condition + sample_id (added by constructor)
  expect_equal(ncol(filtered_coldata), 3)
})

test_that("filter_analysis validates input type", {
  bad_input <- "not_an_analysis"
  
  expect_error(
    filter_analysis(bad_input),
    "TSENATAnalysis"
  )
})

test_that("filter_analysis accepts stringency parameter", {
  se <- make_test_se(n_genes = 200, n_samples = 30)
  coldata <- S4Vectors::DataFrame(
    pair_id = rep(1:15, each = 2),
    row.names = colnames(se)
  )
  SummarizedExperiment::colData(se) <- coldata
  analysis <- TSENATAnalysis(se)
  
  # Test different stringency levels
  for (stringency in c("soft", "medium", "severe")) {
    result <- filter_analysis(analysis, stringency = stringency, verbose = FALSE)
    expect_true(inherits(result, "TSENATAnalysis"))
  }
})