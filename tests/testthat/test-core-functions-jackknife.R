context("Core Internal Functions: Jackknife Isoform Switching")

# Load helpers
source("helper-jackknife-consolidation.R")

# Setup test data
set.seed(42)
create_test_se <- function() {
  counts_matrix <- matrix(
    c(100, 150, 80, 50, 75, 60,    # Gene g1, transcripts tx1-tx3: condition A (2 samples)
      120, 140, 70, 60, 85, 65,    # Gene g1, transcripts tx1-tx3: condition B (2 samples)
      200, 180, 90, 100, 110, 95,  # Gene g2, transcripts tx4-tx6: condition A
      215, 195, 105, 120, 130, 100),# Gene g2, transcripts tx4-tx6: condition B
    nrow = 6, ncol = 4, byrow = FALSE
  )
  
  SummarizedExperiment::SummarizedExperiment(
    assays = list(counts = counts_matrix),
    rowData = data.frame(
      transcript_id = c("tx1", "tx2", "tx3", "tx4", "tx5", "tx6"),
      gene_id = c("g1", "g1", "g1", "g2", "g2", "g2"),
      gene_name = c("GENE1", "GENE1", "GENE1", "GENE2", "GENE2", "GENE2"),
      stringsAsFactors = FALSE
    ),
    colData = data.frame(
      sample = c("s1", "s2", "s3", "s4"),
      condition = c("A", "A", "B", "B"),
      individual_id = c("ind1", "ind1", "ind1", "ind1"),
      stringsAsFactors = FALSE
    )
  )
}

create_paired_se <- function() {
  counts_matrix <- matrix(
    c(100, 150, 80, 50, 75, 60,    # Gene g1, transcripts tx1-tx3
      120, 140, 70, 60, 85, 65,    # Gene g1, transcripts tx1-tx3
      200, 180, 90, 100, 110, 95,  # Gene g2, transcripts tx4-tx6
      215, 195, 105, 120, 130, 100),# Gene g2, transcripts tx4-tx6
    nrow = 6, ncol = 4, byrow = FALSE
  )
  
  SummarizedExperiment::SummarizedExperiment(
    assays = list(counts = counts_matrix),
    rowData = data.frame(
      transcript_id = c("tx1", "tx2", "tx3", "tx4", "tx5", "tx6"),
      gene_id = c("g1", "g1", "g1", "g2", "g2", "g2"),
      gene_name = c("GENE1", "GENE1", "GENE1", "GENE2", "GENE2", "GENE2"),
      stringsAsFactors = FALSE
    ),
    colData = data.frame(
      sample = c("s1", "s2", "s3", "s4"),
      condition = c("A", "A", "B", "B"),
      individual_id = c("ind1", "ind2", "ind1", "ind2"),  # Paired design
      stringsAsFactors = FALSE
    )
  )
}

# ============================================================================
# TEST: .jis_validate_input()
# ============================================================================

test_that(".jis_validate_input validates SE is a SummarizedExperiment", {
  se <- create_test_se()
  invalid_se <- list(data = "not_a_se")
  
  # Valid input should not throw error
  expect_error(TSENAT:::.jis_validate_input(se, "condition", "gene_id", "transcript_id"), NA)
  
  # Invalid SE should throw error
  expect_error(
    TSENAT:::.jis_validate_input(invalid_se, "condition", "gene_id", "transcript_id"),
    "SummarizedExperiment"
  )
})

test_that(".jis_validate_input detects missing condition_col", {
  se <- create_test_se()
  
  expect_error(
    TSENAT:::.jis_validate_input(se, "missing_col", "gene_id", "transcript_id"),
    "not found in colData"
  )
})

test_that(".jis_validate_input detects exactly 2 conditions", {
  se <- create_test_se()
  
  # Add a third condition
  se_bad <- se
  cd <- SummarizedExperiment::colData(se_bad)
  cd$condition <- factor(c("A", "A", "B", "C"))
  SummarizedExperiment::colData(se_bad) <- cd
  
  expect_error(
    TSENAT:::.jis_validate_input(se_bad, "condition", "gene_id", "transcript_id"),
    "Exactly 2 conditions"
  )
})

test_that(".jis_validate_input detects missing gene_col", {
  se <- create_test_se()
  
  expect_error(
    TSENAT:::.jis_validate_input(se, "condition", "missing_gene", "transcript_id"),
    "not found in rowData"
  )
})

test_that(".jis_validate_input detects missing isoform_col", {
  se <- create_test_se()
  
  expect_error(
    TSENAT:::.jis_validate_input(se, "condition", "gene_id", "missing_iso"),
    "not found in rowData"
  )
})

test_that(".jis_validate_input returns sorted conditions", {
  se <- create_test_se()
  
  conditions <- TSENAT:::.jis_validate_input(se, "condition", "gene_id", "transcript_id")
  
  expect_equal(conditions, c("A", "B"))  # Should be sorted alphabetically
})

# ============================================================================
# TEST: .setup_paired_design_jis()
# ============================================================================

test_that(".setup_paired_design_jis returns FALSE for unpaired design", {
  se <- create_test_se()
  
  result <- TSENAT:::.setup_paired_design_jis(se, NULL, "condition")
  
  expect_false(result$is_paired)
  expect_null(result$pair_info)
  expect_null(result$subject_col)
})

test_that(".setup_paired_design_jis detects paired design", {
  se <- create_paired_se()
  
  result <- TSENAT:::.setup_paired_design_jis(se, "individual_id", "condition")
  
  expect_true(result$is_paired)
  expect_true("n_pairs" %in% names(result$pair_info))
  expect_true("matched_pairs" %in% names(result$pair_info))
  expect_equal(result$subject_col, "individual_id")
})

test_that(".setup_paired_design_jis handles unpaired scenario with high ratio", {
  se <- create_paired_se()
  
  # Even with high pairing ratio, if subject_col is NULL, should be unpaired
  result <- TSENAT:::.setup_paired_design_jis(se, NULL, "condition")
  
  expect_false(result$is_paired)
})

test_that(".setup_paired_design_jis detects missing subject_col", {
  se <- create_paired_se()
  
  expect_error(
    TSENAT:::.setup_paired_design_jis(se, "missing_col", "condition"),
    "not found in colData"
  )
})

# ============================================================================
# TEST: .build_gene_id_mapping()
# ============================================================================

test_that(".build_gene_id_mapping creates gene_id to gene_name mapping", {
  se <- create_test_se()
  
  mapping <- TSENAT:::.build_gene_id_mapping(se, "gene_id")
  
  expect_true(is.character(mapping))
  expect_true(length(mapping) > 0)
  # Check that mapping contains g1 or g2
  expect_true(any(c("g1", "g2") %in% names(mapping)))
})

test_that(".build_gene_id_mapping handles genes with same transpose multiple rows", {
  # Some genes appear in multiple rows, mapping should have unique genes only
  se <- create_test_se()
  
  mapping <- TSENAT:::.build_gene_id_mapping(se, "gene_id")
  
  # Should have only 2 unique genes, not 6 (one per transcript)
  expect_true(length(mapping) <= 2)
})

test_that(".build_gene_id_mapping returns empty for missing gene_name column", {
  se <- create_test_se()
  
  # Remove gene_name column
  rd <- SummarizedExperiment::rowData(se)
  rd$gene_name <- NULL
  SummarizedExperiment::rowData(se) <- rd
  
  mapping <- TSENAT:::.build_gene_id_mapping(se, "gene_id")
  
  expect_true(length(mapping) == 0 || all(is.na(mapping)))
})

# ============================================================================
# TEST: .jis_tsallis_entropy()
# ============================================================================

test_that(".jis_tsallis_entropy_fast calculates entropy via C++", {
  # Convert vector to matrix (rows=transcripts, cols=samples)
  counts <- matrix(c(100, 50, 75, 200, 80, 120), nrow = 6, ncol = 1)
  
  result <- TSENAT:::.jis_tsallis_entropy_fast(counts, q = 1, norm = TRUE, log_base = exp(1), pseudocount = 0)
  
  expect_true(is.numeric(result))
  expect_equal(length(result), 1)
  expect_true(result >= 0)
})

test_that(".jis_tsallis_entropy_fast with q=1 via C++", {
  counts <- matrix(c(100, 100, 100, 100), nrow = 4, ncol = 1)  # Balanced distribution
  
  # Balanced distribution should have high entropy
  result <- TSENAT:::.jis_tsallis_entropy_fast(counts, q = 1, norm = FALSE, log_base = exp(1), pseudocount = 0)
  
  expect_true(result > 1.3)  # Shannon entropy of balanced 4-item distribution
})

test_that(".jis_tsallis_entropy_fast with q=2 via C++", {
  counts <- matrix(c(100, 50, 75), nrow = 3, ncol = 1)
  
  result <- TSENAT:::.jis_tsallis_entropy_fast(counts, q = 2, norm = TRUE, log_base = exp(1), pseudocount = 0)
  
  # Tsallis entropy with q > 1 can be negative - this is mathematically correct
  # The normalized value should be in [-1, 1] approximately
  expect_true(is.finite(result))
  expect_true(result >= -1.1)
  expect_true(result <= 1.1)
})

test_that(".jis_tsallis_entropy_fast handles normalization via C++", {
  counts <- matrix(c(100, 50, 75), nrow = 3, ncol = 1)
  
  result_norm <- TSENAT:::.jis_tsallis_entropy_fast(counts, q = 1, norm = TRUE, log_base = exp(1), pseudocount = 0)
  result_no_norm <- TSENAT:::.jis_tsallis_entropy_fast(counts, q = 1, norm = FALSE, log_base = exp(1), pseudocount = 0)
  
  expect_true(result_norm <= result_no_norm)  # Normalized should be smaller or equal
})

test_that(".jis_tsallis_entropy_fast handles matrix input (per-sample via C++)", {
  counts_matrix <- matrix(c(100, 50, 75, 110, 45, 80), nrow = 3, ncol = 2)
  
  result <- TSENAT:::.jis_tsallis_entropy_fast(counts_matrix, q = 1, norm = TRUE, log_base = exp(1), pseudocount = 0)
  
  expect_true(is.numeric(result))
  expect_equal(length(result), 2)  # One entropy per column (sample)
})

test_that(".jis_tsallis_entropy_fast handles pseudocount via C++", {
  counts <- matrix(c(100, 0, 75), nrow = 3, ncol = 1)  # Has zero count
  
  # Without pseudocount might have numerical issues
  result_with_pc <- TSENAT:::.jis_tsallis_entropy_fast(counts, q = 1, norm = TRUE, log_base = exp(1), pseudocount = 0.5)
  
  expect_true(is.finite(result_with_pc))
})

test_that(".jis_tsallis_entropy_fast respects log base via C++", {
  counts <- matrix(c(100, 50, 75), nrow = 3, ncol = 1)
  
  result_e <- TSENAT:::.jis_tsallis_entropy_fast(counts, q = 1, norm = FALSE, log_base = exp(1), pseudocount = 0)
  result_2 <- TSENAT:::.jis_tsallis_entropy_fast(counts, q = 1, norm = FALSE, log_base = 2, pseudocount = 0)
  
  # Results should differ due to different log base
  expect_false(isTRUE(all.equal(result_e, result_2)))
  # log_2(x) = log_e(x) / log_e(2), so log base 2 produces LARGER entropy values
  # Expected relationship: result_2 > result_e (approximately by factor of ln(2) ≈ 1.443)
  expect_true(result_2 > result_e)
  # Validate the log base conversion ratio
  ratio <- result_2 / result_e
  expect_true(ratio > 1.4 && ratio < 1.5)  # Should be close to 1/ln(2) ≈ 1.443
})

# ============================================================================
# TEST: .jackknife_influences_jis()
# ============================================================================

test_that(".jis_jackknife_influences_fast calculates influences via C++", {
  counts_matrix <- matrix(c(100, 50, 75, 110, 45, 80), nrow = 3, ncol = 2)
  
  influences <- TSENAT:::.jis_jackknife_influences_fast(counts_matrix, q = 1, norm = TRUE, log_base = exp(1), pseudocount = 0)
  
  expect_true(is.numeric(influences))
  expect_equal(length(influences), 3)  # One influence per row (transcript)
  expect_true(all(influences >= 0))
})

test_that(".jis_jackknife_influences_fast identifies outliers via C++", {
  # Create data where first transcript is dominant
  counts_matrix <- matrix(c(1000, 50, 75, 900, 45, 80), nrow = 3, ncol = 2)
  
  influences <- TSENAT:::.jis_jackknife_influences_fast(counts_matrix, q = 1, norm = TRUE, log_base = exp(1), pseudocount = 0)
  
  # First transcript (dominant) should have higher influence
  expect_true(influences[1] > influences[2])
  expect_true(influences[1] > influences[3])
})

test_that(".jis_jackknife_influences_fast with n_tx_fixed via C++", {
  counts_matrix <- matrix(c(100, 50, 75, 110, 45, 80), nrow = 3, ncol = 2)
  
  influences_fixed <- TSENAT:::.jis_jackknife_influences_fast(counts_matrix, q = 1, norm = TRUE, 
                                                               log_base = exp(1), pseudocount = 0, n_tx_fixed = 5)
  influences_unfixed <- TSENAT:::.jis_jackknife_influences_fast(counts_matrix, q = 1, norm = TRUE,
                                                                 log_base = exp(1), pseudocount = 0, n_tx_fixed = NULL)
  
  # Results should differ when n_tx_fixed is specified
  expect_false(isTRUE(all.equal(influences_fixed, influences_unfixed)))
})

# ============================================================================
# TEST: .jis_apply_fdr()
# ============================================================================

test_that(".jis_apply_fdr corrects p-values with Benjamini-Hochberg method", {
  # Create mock gene results with p-values
  results_per_gene <- list(
    g1 = list(
      transcript_ids = c("tx1", "tx2"),
      delta_pvalue = c(0.001, 0.01)
    ),
    g2 = list(
      transcript_ids = c("tx3", "tx4"),
      delta_pvalue = c(0.05, 0.1)
    )
  )
  
  all_pvalues <- list(
    list(gene = "g1", transcript = "tx1", pvalue = 0.001),
    list(gene = "g1", transcript = "tx2", pvalue = 0.01),
    list(gene = "g2", transcript = "tx3", pvalue = 0.05),
    list(gene = "g2", transcript = "tx4", pvalue = 0.1)
  )
  
  results_per_gene <- TSENAT:::.jis_apply_fdr(results_per_gene, all_pvalues)
  
  # Check that FDR values were added to results_per_gene
  expect_true("delta_fdr" %in% names(results_per_gene$g1))
  expect_equal(length(results_per_gene$g1$delta_fdr), 2)
  expect_true(all(results_per_gene$g1$delta_fdr < 1))
})

test_that(".jis_apply_fdr maintains p-value ordering relationships", {
  results_per_gene <- list(
    g1 = list(
      transcript_ids = c("tx1", "tx2", "tx3"),
      delta_pvalue = c(0.001, 0.01, 0.05)
    )
  )
  
  all_pvalues <- list(
    list(gene = "g1", transcript = "tx1", pvalue = 0.001),
    list(gene = "g1", transcript = "tx2", pvalue = 0.01),
    list(gene = "g1", transcript = "tx3", pvalue = 0.05)
  )
  
  results_per_gene <- TSENAT:::.jis_apply_fdr(results_per_gene, all_pvalues)
  
  # FDR values should maintain order (smallest p-value → smallest FDR)
  fdr_vals <- results_per_gene$g1$delta_fdr
  expect_true(fdr_vals[1] <= fdr_vals[2])
  expect_true(fdr_vals[2] <= fdr_vals[3])
})

# ============================================================================
# TEST: .jis_handle_multi_q()
# ============================================================================

test_that(".jis_handle_multi_q returns list of results for multiple q values", {
  skip_on_cran()
  
  se <- create_test_se()
  
  q_values <- c(1.0, 1.5)
  q_params <- list(
    condition_col = "condition",
    gene_col = "gene_name",
    isoform_col = "transcript_id",
    norm = TRUE,
    log_base = exp(1),
    pseudocount = 0,
    n_bootstrap = 10
  )
  
  result <- TSENAT:::.jis_handle_multi_q(se, q_values, q_params, verbose = FALSE)
  
  expect_true(is.list(result))
  expect_equal(names(result), c("q_1_00", "q_1_50"))
})

# ============================================================================
# TEST: Integration tests for helper functions
# ============================================================================

test_that("Validation and pairing helpers work together", {
  se <- create_paired_se()
  
  # Validate input first
  conditions <- TSENAT:::.jis_validate_input(se, "condition", "gene_id", "transcript_id")
  expect_equal(conditions, c("A", "B"))
  
  # Then setup paired design
  paired_info <- TSENAT:::.setup_paired_design_jis(se, "individual_id", "condition")
  expect_true(paired_info$is_paired)
})

test_that("Gene mapping and entropy calculation work together", {
  se <- create_test_se()
  
  # Build gene mapping
  mapping <- TSENAT:::.build_gene_id_mapping(se, "gene_id")
  expect_true(length(mapping) > 0)
  
  # Extract counts for a gene and calculate entropy
  gene_mask <- rowData(se)$gene_id == "g1"
  counts_matrix <- assays(se)$counts[gene_mask, ]
  
  entropy <- TSENAT:::.jis_tsallis_entropy_fast(counts_matrix, q = 1, norm = TRUE, log_base = exp(1), pseudocount = 0)
  expect_true(is.numeric(entropy))
  expect_true(all(entropy >= 0 & entropy <= 1))
})
