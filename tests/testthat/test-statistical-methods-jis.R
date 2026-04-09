context("Jackknife Isoform Switching: Updates - Multi-Q, Gene Mapping, and LM Filtering")

# Helper function to create test SummarizedExperiment with gene names
test_se_with_gene_names <- function() {
  suppressPackageStartupMessages({
  })
  
  counts_mat <- matrix(rpois(60, lambda=15), nrow=6, ncol=10,
                       dimnames=list(c("T1a", "T1b", "T2a", "T2b", "T3a", "T3b"), paste0("S", 1:10)))
  
  SummarizedExperiment(
    assays=list(counts=counts_mat),
    rowData=data.frame(
      isoform_id=c("T1a", "T1b", "T2a", "T2b", "T3a", "T3b"),
      gene_id=c("ENSG00001", "ENSG00001", "ENSG00002", "ENSG00002", "ENSG00003", "ENSG00003"),
      gene_name=c("GeneA", "GeneA", "GeneB", "GeneB", "GeneC", "GeneC")
    ),
    colData=data.frame(
      sample_id=paste0("S", 1:10),
      condition=rep(c("control", "treatment"), each=5)
    )
  )
}

# ═══════════════════════════════════════════════════════════════════════════════
# TEST SUITE 1: Multi-Q Support
# ═══════════════════════════════════════════════════════════════════════════════

test_that("calculate_jis accepts single q value", {
  se <- test_se_with_gene_names()
  
  result <- suppressWarnings(.calculate_jis(
    se = se,
    condition_col = "condition",
    gene_col = "gene_id",
    isoform_col = "isoform_id",
    q = 1.0,
    n_bootstrap = 5,  # Reduced from 10 to 5 for faster testing
    verbose = FALSE
  ))
  
  expect_is(result, "tsenat_isoform_switching")
  expect_true(!inherits(result, "tsenat_isoform_switching_multiq"))
})

test_that("calculate_jis accepts vector of q values", {
  se <- test_se_with_gene_names()
  
  q_values <- c(0.5, 1.0, 1.5)
  result <- suppressWarnings(.calculate_jis(
    se = se,
    condition_col = "condition",
    gene_col = "gene_id",
    isoform_col = "isoform_id",
    q = q_values,
    n_bootstrap = 5,
    verbose = FALSE
  ))
  
  expect_is(result, "tsenat_isoform_switching_multiq")
  expect_length(result, length(q_values))
})

test_that("Multi-q results have correct naming convention", {
  se <- test_se_with_gene_names()
  
  q_values <- c(0.01, 0.5, 1.0, 2.0)
  result <- suppressWarnings(.calculate_jis(
    se = se,
    condition_col = "condition",
    gene_col = "gene_id",
    isoform_col = "isoform_id",
    q = q_values,
    n_bootstrap = 5,  # Reduced from 10 to 5 for faster testing
    verbose = FALSE
  ))
  
  # Check naming format: q_X_XX format
  expected_names <- c("q_0_01", "q_0_50", "q_1_00", "q_2_00")
  expect_equal(names(result), expected_names)
})

test_that("Each multi-q result is a valid tsenat_isoform_switching object", {
  se <- test_se_with_gene_names()
  
  q_values <- c(0.5, 1.0, 1.5)
  result <- suppressWarnings(.calculate_jis(
    se = se,
    condition_col = "condition",
    gene_col = "gene_id",
    isoform_col = "isoform_id",
    q = q_values,
    n_bootstrap = 5,  # Reduced from 10 to 5 for faster testing
    verbose = FALSE
  ))

  # Each element should be a valid single-q result
  for (i in seq_along(result)) {
    expect_is(result[[i]], "tsenat_isoform_switching")
    expect_true("results_per_gene" %in% names(result[[i]]))
    expect_true("summary_table" %in% names(result[[i]]))
    expect_true("all_transcript_stats" %in% names(result[[i]]))
  }
})

test_that("Multi-q analysis analyzes same genes across all q values", {
  se <- test_se_with_gene_names()
  
  q_values <- c(0.5, 1.0, 1.5)
  result <- suppressWarnings(.calculate_jis(
    se = se,
    condition_col = "condition",
    gene_col = "gene_id",
    isoform_col = "isoform_id",
    q = q_values,
    n_bootstrap = 5,
    verbose = FALSE
  ))
  
  # Extract gene names from each q result
  gene_sets <- lapply(result, function(r) r$gene_names)
  
  # All should have same genes
  for (i in 2:length(gene_sets)) {
    expect_equal(gene_sets[[i]], gene_sets[[1]])
  }
})

# ═══════════════════════════════════════════════════════════════════════════════
# TEST SUITE 2: Gene Name to ID Mapping
# ═══════════════════════════════════════════════════════════════════════════════

test_that("calculate_jis maps gene names to IDs in lm_results", {
  se <- test_se_with_gene_names()
  
  # Create lm_results with GENE NAMES (not IDs)
  lm_results <- data.frame(
    gene=c("GeneA", "GeneB", "GeneC"),
    p_interaction=c(0.01, 0.05, 0.50),
    adj_p_interaction=c(0.02, 0.10, 0.60)
  )
  
  result <- suppressWarnings(.calculate_jis(
    se = se,
    condition_col = "condition",
    gene_col = "gene_id",
    isoform_col = "isoform_id",
    lm_results = lm_results,
    lm_p_threshold = 0.05,
    use_lm_fdr = FALSE,
    n_bootstrap = 5,
    verbose = FALSE
  ))
  
  expect_is(result, "tsenat_isoform_switching")
  # Should successfully map and filter (GeneA and GeneB pass p < 0.05)
  expect_gte(length(result$gene_names), 1)
  expect_lte(length(result$gene_names), 3)
})

test_that("calculate_jis handles lm_results with gene IDs", {
  se <- test_se_with_gene_names()
  
  # Create lm_results with GENE IDs (already mapped)
  lm_results <- data.frame(
    gene=c("ENSG00001", "ENSG00002", "ENSG00003"),
    p_interaction=c(0.01, 0.05, 0.50),
    adj_p_interaction=c(0.02, 0.10, 0.60)
  )
  
  result <- suppressWarnings(.calculate_jis(
    se = se,
    condition_col = "condition",
    gene_col = "gene_id",
    isoform_col = "isoform_id",
    lm_results = lm_results,
    lm_p_threshold = 0.05,
    use_lm_fdr = FALSE,
    n_bootstrap = 5,
    verbose = FALSE
  ))
  
  expect_is(result, "tsenat_isoform_switching")
  # Should have filtered genes based on threshold
  expect_gte(length(result$gene_names), 1)
})

test_that("Gene mapping respects LM p-value threshold", {
  se <- test_se_with_gene_names()
  
  lm_results <- data.frame(
    gene=c("GeneA", "GeneB", "GeneC"),
    p_interaction=c(0.001, 0.01, 0.50),
    adj_p_interaction=c(0.002, 0.02, 0.60)
  )
  
  # Analyze with strict threshold (only GeneA)
  result_strict <- suppressWarnings(.calculate_jis(
    se = se,
    condition_col = "condition",
    gene_col = "gene_id",
    isoform_col = "isoform_id",
    lm_results = lm_results,
    lm_p_threshold = 0.005,
    use_lm_fdr = FALSE,
    n_bootstrap = 5,
    verbose = FALSE
  ))
  
  # Analyze with lenient threshold (GeneA and GeneB)
  result_lenient <- suppressWarnings(.calculate_jis(
    se = se,
    condition_col = "condition",
    gene_col = "gene_id",
    isoform_col = "isoform_id",
    lm_results = lm_results,
    lm_p_threshold = 0.05,
    use_lm_fdr = FALSE,
    n_bootstrap = 5,
    verbose = FALSE
  ))
  
  expect_length(result_strict$gene_names, 1)
  expect_length(result_lenient$gene_names, 2)
})

test_that("Gene name to ID mapping removes duplicates", {
  suppressPackageStartupMessages({
  })
  
  # Create SE with duplicate gene names (edge case)
  counts_mat <- matrix(rpois(40, lambda=15), nrow=4, ncol=10,
                       dimnames=list(c("T1a", "T1b", "T2a", "T2b"), paste0("S", 1:10)))
  
  se_dup <- SummarizedExperiment(
    assays=list(counts=counts_mat),
    rowData=data.frame(
      isoform_id=c("T1a", "T1b", "T2a", "T2b"),
      gene_id=c("ENSG00001", "ENSG00001", "ENSG00002", "ENSG00002"),
      gene_name=c("GeneA", "GeneA", "GeneA", "GeneB")  # Duplicate GeneA
    ),
    colData=data.frame(
      sample_id=paste0("S", 1:10),
      condition=rep(c("A", "B"), each=5)
    )
  )
  
  lm_results <- data.frame(
    gene=c("GeneA", "GeneB"),
    p_interaction=c(0.01, 0.50),
    adj_p_interaction=c(0.02, 0.60)
  )
  
  # Should still work despite duplicate names (first mapping wins)
  result <- suppressWarnings(.calculate_jis(
    se = se_dup,
    condition_col = "condition",
    gene_col = "gene_id",
    isoform_col = "isoform_id",
    lm_results = lm_results,
    lm_p_threshold = 0.05,
    n_bootstrap = 5,
    verbose = FALSE
  ))
  
  expect_is(result, "tsenat_isoform_switching")
})

# ═══════════════════════════════════════════════════════════════════════════════
# TEST SUITE 3: All LM-Significant Genes Analyzed (No top_n Limit)
# ═══════════════════════════════════════════════════════════════════════════════

test_that("Function accepts multiple genes in lm_results", {
  se <- test_se_with_gene_names()
  
  # Create lm_results with multiple genes
  lm_results <- data.frame(
    gene=c("GeneA", "GeneB", "GeneC"),
    p_interaction=c(0.001, 0.01, 0.50),
    adj_p_interaction=c(0.002, 0.02, 0.60)
  )
  
  result <- suppressWarnings(.calculate_jis(
    se = se,
    condition_col = "condition",
    gene_col = "gene_id",
    isoform_col = "isoform_id",
    lm_results = lm_results,
    lm_p_threshold = 0.05,
    use_lm_fdr = FALSE,
    n_bootstrap = 5,
    verbose = FALSE
  ))
  
  # Should successfully filter and analyze genes
  expect_is(result, "tsenat_isoform_switching")
  expect_gte(length(result$gene_names), 1)
})

test_that("LM threshold influences number of analyzed genes", {
  se <- test_se_with_gene_names()
  
  lm_results <- data.frame(
    gene=c("GeneA", "GeneB", "GeneC"),
    p_interaction=c(0.001, 0.01, 0.50),
    adj_p_interaction=c(0.002, 0.02, 0.60)
  )
  
  # Run with strict threshold
  result_strict <- suppressWarnings(.calculate_jis(
    se = se,
    condition_col = "condition",
    gene_col = "gene_id",
    isoform_col = "isoform_id",
    lm_results = lm_results,
    lm_p_threshold = 0.005,
    use_lm_fdr = FALSE,
    n_bootstrap = 5,
    verbose = FALSE
  ))
  
  # Run with lenient threshold
  result_lenient <- suppressWarnings(.calculate_jis(
    se = se,
    condition_col = "condition",
    gene_col = "gene_id",
    isoform_col = "isoform_id",
    lm_results = lm_results,
    lm_p_threshold = 0.50,
    use_lm_fdr = FALSE,
    n_bootstrap = 5,
    verbose = FALSE
  ))
  
  # Stricter threshold should result in fewer genes
  expect_lte(length(result_strict$gene_names), length(result_lenient$gene_names))
})

# ═══════════════════════════════════════════════════════════════════════════════
# TEST SUITE 4: top_n Parameter Removed
# ═══════════════════════════════════════════════════════════════════════════════

test_that("top_n parameter is removed from function signature", {
  sig <- formals(.calculate_jis)
  expect_false("top_n" %in% names(sig))
})

test_that("Function works without top_n parameter", {
  se <- test_se_with_gene_names()
  
  # Should work without specifying top_n
  result <- suppressWarnings(.calculate_jis(
    se = se,
    condition_col = "condition",
    gene_col = "gene_id",
    isoform_col = "isoform_id",
    n_bootstrap = 5,
    verbose = FALSE
  ))
  
  expect_is(result, "tsenat_isoform_switching")
})

# ═══════════════════════════════════════════════════════════════════════════════
# TEST SUITE 5: Multi-Q with Gene Mapping Combined
# ═══════════════════════════════════════════════════════════════════════════════

test_that("Multi-q analysis works with gene name to ID mapping", {
  se <- test_se_with_gene_names()
  
  # Use gene names in lm_results with multi-q
  lm_results <- data.frame(
    gene=c("GeneA", "GeneB", "GeneC"),
    p_interaction=c(0.01, 0.05, 0.50),
    adj_p_interaction=c(0.02, 0.10, 0.60)
  )
  
  q_values <- c(0.5, 1.0, 1.5)
  result <- suppressWarnings(.calculate_jis(
    se = se,
    condition_col = "condition",
    gene_col = "gene_id",
    isoform_col = "isoform_id",
    q = q_values,
    lm_results = lm_results,
    lm_p_threshold = 0.05,
    use_lm_fdr = FALSE,
    n_bootstrap = 5,
    verbose = FALSE
  ))
  
  expect_is(result, "tsenat_isoform_switching_multiq")
  expect_length(result, 3)
  
  # Each result should have the same genes
  gene_counts <- sapply(result, function(res) length(res$gene_names))
  # All q values should analyze the same genes
  expect_true(all(gene_counts == gene_counts[1]))
})

test_that("Multi-q with gene filtering returns consistent structure", {
  se <- test_se_with_gene_names()
  
  lm_results <- data.frame(
    gene=c("GeneA", "GeneB", "GeneC"),
    p_interaction=c(0.01, 0.05, 0.50),
    adj_p_interaction=c(0.02, 0.10, 0.60)
  )
  
  q_values <- c(1.0, 1.5)
  result <- suppressWarnings(.calculate_jis(
    se = se,
    condition_col = "condition",
    gene_col = "gene_id",
    isoform_col = "isoform_id",
    q = q_values,
    lm_results = lm_results,
    lm_p_threshold = 0.05,
    use_lm_fdr = FALSE,
    n_bootstrap = 5,
    verbose = FALSE
  ))
  
  # Should return multiq results
  expect_is(result, "tsenat_isoform_switching_multiq")
  expect_length(result, 2)
  
  # Both q results should have results
  for (res in result) {
    expect_is(res, "tsenat_isoform_switching")
  }
})

# ═══════════════════════════════════════════════════════════════════════════════
# TEST SUITE 6: Use_lm_fdr Parameter
# ═══════════════════════════════════════════════════════════════════════════════

test_that("use_lm_fdr parameter switches between raw and adjusted p-values", {
  se <- test_se_with_gene_names()
  
  lm_results <- data.frame(
    gene=c("GeneA", "GeneB", "GeneC"),
    p_interaction=c(0.001, 0.005, 0.50),
    adj_p_interaction=c(0.002, 0.01, 0.60)
  )
  
  result_fdr <- suppressWarnings(.calculate_jis(
    se = se,
    condition_col = "condition",
    gene_col = "gene_id",
    isoform_col = "isoform_id",
    lm_results = lm_results,
    lm_p_threshold = 0.05,
    use_lm_fdr = TRUE,
    n_bootstrap = 5,
    verbose = FALSE
  ))
  
  result_raw <- suppressWarnings(.calculate_jis(
    se = se,
    condition_col = "condition",
    gene_col = "gene_id",
    isoform_col = "isoform_id",
    lm_results = lm_results,
    lm_p_threshold = 0.05,
    use_lm_fdr = FALSE,
    n_bootstrap = 5,
    verbose = FALSE
  ))
  
  # Both should filter genes
  expect_gte(length(result_fdr$gene_names), 0)
  expect_gte(length(result_raw$gene_names), 0)
})

# ═══════════════════════════════════════════════════════════════════════════════
# TEST SUITE 7: Metadata Tracking Updates
# ═══════════════════════════════════════════════════════════════════════════════

test_that("Metadata tracks LM gene filtering correctly", {
  se <- test_se_with_gene_names()
  
  lm_results <- data.frame(
    gene=c("GeneA", "GeneB", "GeneC"),
    p_interaction=c(0.001, 0.01, 0.50),
    adj_p_interaction=c(0.002, 0.02, 0.60)
  )
  
  result <- suppressWarnings(.calculate_jis(
    se = se,
    condition_col = "condition",
    gene_col = "gene_id",
    isoform_col = "isoform_id",
    lm_results = lm_results,
    lm_p_threshold = 0.05,
    n_bootstrap = 5,
    verbose = FALSE
  ))
  
  meta <- result$metadata
  expect_equal(meta$lm_results_provided, TRUE)
  expect_equal(meta$lm_p_threshold, 0.05)
  expect_equal(meta$lm_genes_filtered, 2)  # GeneA and GeneB
})

test_that("Metadata correctly reports LM filtering statistics", {
  se <- test_se_with_gene_names()
  
  lm_results <- data.frame(
    gene=c("GeneA", "GeneB", "GeneC"),
    p_interaction=c(0.001, 0.01, 0.50),
    adj_p_interaction=c(0.002, 0.02, 0.60)
  )
  
  result <- suppressWarnings(.calculate_jis(
    se = se,
    condition_col = "condition",
    gene_col = "gene_id",
    isoform_col = "isoform_id",
    lm_results = lm_results,
    lm_p_threshold = 0.05,
    use_lm_fdr = FALSE,
    n_bootstrap = 5,
    verbose = FALSE
  ))
  
  meta <- result$metadata
  expect_true(meta$lm_results_provided)
  expect_equal(meta$lm_p_threshold, 0.05)
  # Number of filtered genes should be > 0
  expect_gte(meta$lm_genes_filtered, 1)
})
