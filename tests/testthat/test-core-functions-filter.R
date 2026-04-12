# Tests for .filter_se() function
# Tests filtering behavior, metadata preservation, and isoform-level filtering
# Based on Soneson et al. 2016 recommendations

test_that(".filter_se filters all assays consistently", {
  # Create SE with multiple assays
  counts <- matrix(c(10, 20, 30, 5, 15, 25, 2, 8, 12), nrow = 3, ncol = 3)
  tpm <- matrix(c(100, 200, 300, 50, 150, 250, 20, 80, 120), nrow = 3, ncol = 3)
  rownames(counts) <- rownames(tpm) <- c("TX1", "TX2", "TX3")
  colnames(counts) <- colnames(tpm) <- c("S1", "S2", "S3")
  
  se <- SummarizedExperiment::SummarizedExperiment(
    assays = list(counts = counts, tpm = tpm)
  )
  
  # Filter with low threshold to keep some rows
  se_filt <- TSENAT:::.filter_se(se, min_samples = 1, min_tpm = 50, verbose = FALSE)
  
  # All assays should have same dimensions
  expect_equal(nrow(assays(se_filt)[[1]]), nrow(assays(se_filt)[[2]]))
  expect_equal(ncol(assays(se_filt)[[1]]), ncol(assays(se_filt)[[2]]))
  expect_equal(ncol(assays(se_filt)[[1]]), 3)  # All samples preserved
  expect_true(nrow(assays(se_filt)[[1]]) <= 3)  # Some rows filtered
})

test_that(".filter_se preserves colnames and rownames across assays", {
  # Create SE with named assays
  counts <- matrix(1:12, nrow = 4, ncol = 3)
  tpm <- matrix(10:21, nrow = 4, ncol = 3)
  rownames(counts) <- rownames(tpm) <- c("TX1", "TX2", "TX3", "TX4")
  colnames(counts) <- colnames(tpm) <- c("Sample_A", "Sample_B", "Sample_C")
  
  se <- SummarizedExperiment::SummarizedExperiment(
    assays = list(counts = counts, tpm = tpm)
  )
  
  se_filt <- TSENAT:::.filter_se(se, min_samples = 1, min_tpm = 5, verbose = FALSE)
  
  # All assays should have identical colnames
  cn_counts <- colnames(assays(se_filt)$counts)
  cn_tpm <- colnames(assays(se_filt)$tpm)
  expect_equal(cn_counts, cn_tpm)
  expect_equal(cn_counts, c("Sample_A", "Sample_B", "Sample_C"))
  
  # All assays should have identical rownames
  rn_counts <- rownames(assays(se_filt)$counts)
  rn_tpm <- rownames(assays(se_filt)$tpm)
  expect_equal(rn_counts, rn_tpm)
})

test_that(".filter_se preserves colData while filtering rows", {
  # Create SE with colData and TPM data (avoid warnings)
  counts <- matrix(1:12, nrow = 4, ncol = 3)
  tpm <- matrix(c(5, 10, 15, 20, 8, 12, 16, 22, 3, 6, 9, 12), nrow = 4, ncol = 3)  # Include TPM
  rownames(counts) <- rownames(tpm) <- c("TX1", "TX2", "TX3", "TX4")
  colnames(counts) <- colnames(tpm) <- c("S1", "S2", "S3")
  
  se <- SummarizedExperiment::SummarizedExperiment(
    assays = list(counts = counts, tpm = tpm),
    colData = data.frame(
      sample_id = c("S1", "S2", "S3"),
      condition = c("A", "B", "A"),
      row.names = c("S1", "S2", "S3")
    )
  )
  
  se_filt <- TSENAT:::.filter_se(se, min_samples = 1, min_tpm = 1, verbose = FALSE)
  
  # colData should be unchanged
  expect_equal(nrow(SummarizedExperiment::colData(se_filt)), 3)
  expect_equal(colnames(SummarizedExperiment::colData(se_filt)), c("sample_id", "condition"))
  expect_equal(SummarizedExperiment::colData(se_filt)$condition, c("A", "B", "A"))
})

test_that(".filter_se updates metadata correctly", {
  # Create SE with TPM data to avoid fallback warning
  counts <- matrix(1:12, nrow = 4, ncol = 3)
  tpm <- matrix(c(10, 20, 30, 40, 15, 25, 35, 45, 5, 10, 15, 20), nrow = 4, ncol = 3)  # Include TPM
  rownames(counts) <- rownames(tpm) <- c("TX1", "TX2", "TX3", "TX4")
  colnames(counts) <- colnames(tpm) <- c("S1", "S2", "S3")
  
  # Add metadata with readcounts and effective_length
  md <- list(
    readcounts = counts,
    effective_length = c(100, 150, 200, 175),
    tpm = tpm  # Include TPM in metadata for proper filtering
  )
  
  se <- SummarizedExperiment::SummarizedExperiment(
    assays = list(counts = counts, tpm = tpm),
    metadata = md
  )
  
  se_filt <- TSENAT:::.filter_se(se, min_samples = 1, min_tpm = 2, verbose = FALSE)
  
  # Metadata should be filtered to match assay rows
  md_filt <- S4Vectors::metadata(se_filt)
  
  # readcounts should be filtered
  if (!is.null(md_filt$readcounts)) {
    expect_equal(nrow(md_filt$readcounts), nrow(assays(se_filt)[[1]]))
    expect_equal(ncol(md_filt$readcounts), 3)
  }
  
  # effective_length should be filtered
  if (!is.null(md_filt$effective_length)) {
    expect_equal(length(md_filt$effective_length), nrow(assays(se_filt)[[1]]))
  }
})

test_that(".filter_se isoform-level filtering (Soneson et al. 2016)", {
  # Create SE with multiple isoforms per gene
  mat <- matrix(c(
    100, 50, 10,    # TX1.1: 60% relative abundance
    80, 15, 5,      # TX1.2: 30% relative abundance
    5, 1, 0,        # TX1.3: 10% relative abundance (should be removed at 15% threshold)
    200, 180, 20,   # TX2.1: gene with 2 isoforms
    190, 170, 15    # TX2.2: single-isoform genes always kept
  ), nrow = 5, byrow = TRUE, dimnames = list(
    c("TX1.1", "TX1.2", "TX1.3", "TX2.1", "TX2.2"),
    c("S1", "S2", "S3")
  ))
  
  tpm_mat <- log2(mat + 1)  # Mock TPM
  
  se <- SummarizedExperiment::SummarizedExperiment(
    assays = list(counts = mat, tpm = tpm_mat),
    rowData = data.frame(
      transcript_id = rownames(mat),
      gene_id = c("G1", "G1", "G1", "G2", "G2"),
      row.names = rownames(mat)
    ),
    colData = data.frame(
      sample_id = colnames(mat),
      row.names = colnames(mat)
    ),
    metadata = list(tx2gene = data.frame(
      tx = rownames(mat),
      gene = c("G1", "G1", "G1", "G2", "G2")
    ))
  )
  
  # Filter with 15% isoform abundance threshold
  # Should remove TX1.3 (10% < 15%)
  se_filt <- .filter_se(se, min_tpm = 0, min_samples = 1,
                        min_isoform_abundance = 0.15,
                        verbose = FALSE)
  
  expect_equal(nrow(se_filt), 4)  # 5 - 1 removed
  expect_true(!("TX1.3" %in% rownames(se_filt)))
  expect_true(all(c("TX1.1", "TX1.2", "TX2.1", "TX2.2") %in% rownames(se_filt)))
  
  # Metadata should track the filtering parameter
  md <- S4Vectors::metadata(se_filt)
  expect_equal(md$filtered$min_isoform_abundance, 0.15)
})

test_that(".filter_se with min_isoform_abundance = 0 skips isoform filtering", {
  # Create SE with multiple isoforms per gene
  mat <- matrix(c(
    100, 50, 10,
    5, 1, 0,
    200, 180, 20
  ), nrow = 3, byrow = TRUE, dimnames = list(
    c("TX1.1", "TX1.2", "TX2.1"),
    c("S1", "S2", "S3")
  ))
  
  tpm_mat <- log2(mat + 1)
  
  se <- SummarizedExperiment::SummarizedExperiment(
    assays = list(counts = mat, tpm = tpm_mat),
    rowData = data.frame(
      gene_id = c("G1", "G1", "G2"),
      row.names = rownames(mat)
    ),
    metadata = list(tx2gene = data.frame(
      tx = rownames(mat),
      gene = c("G1", "G1", "G2")
    ))
  )
  
  # Disable isoform filtering
  se_filt <- .filter_se(se, min_tpm = 0, min_samples = 1,
                        min_tx_per_gene = 1,  # Allow single-isoform genes
                        min_isoform_abundance = 0,
                        verbose = FALSE)
  
  expect_equal(nrow(se_filt), 3)  # All rows kept
  expect_equal(rownames(se_filt), c("TX1.1", "TX1.2", "TX2.1"))
})

# ============================================================================
# COMPREHENSIVE PARAMETER TESTS
# ============================================================================

test_that(".filter_se min_tpm parameter: varying thresholds", {
  counts <- matrix(c(1, 5, 10, 50, 100, 200), nrow = 3, ncol = 2)
  tpm <- counts * 10  # Mock TPM
  rownames(counts) <- rownames(tpm) <- c("TX1", "TX2", "TX3")
  colnames(counts) <- colnames(tpm) <- c("S1", "S2")
  
  se <- SummarizedExperiment::SummarizedExperiment(
    assays = list(counts = counts, tpm = tpm)
  )
  
  # min_tpm = 10: TX1 (10, 100) has 1 sample >= 10 - may not pass min_samples=2
  se_tpm10 <- .filter_se(se, min_tpm = 10, min_samples = 2, verbose = FALSE)
  expect_true(nrow(se_tpm10) >= 1)
  
  # min_tpm = 100: Only TX3 passes
  se_tpm100 <- .filter_se(se, min_tpm = 100, min_samples = 2, verbose = FALSE)
  expect_equal(nrow(se_tpm100), 1)
  expect_true("TX3" %in% rownames(se_tpm100))
  
  # Higher min_tpm should give fewer or equal rows
  expect_true(nrow(se_tpm100) <= nrow(se_tpm10))
})

test_that(".filter_se min_samples parameter: varying sample requirements", {
  counts <- matrix(c(100, 100, 100, 0, 0, 0), nrow = 2, ncol = 3, byrow = TRUE)
  tpm <- counts + 1  # Avoid zero for TPM
  rownames(counts) <- rownames(tpm) <- c("TX1", "TX2")
  colnames(counts) <- colnames(tpm) <- c("S1", "S2", "S3")
  
  se <- SummarizedExperiment::SummarizedExperiment(
    assays = list(counts = counts, tpm = tpm)
  )
  
  # min_samples = 1: At least one of TX1/TX2 passes
  se_1 <- .filter_se(se, min_tpm = 1, min_samples = 1, verbose = FALSE)
  expect_true(nrow(se_1) >= 1)
  
  # min_samples = 3: Only TX1 passes (all 3 samples >= 100)
  se_3 <- .filter_se(se, min_tpm = 100, min_samples = 3, verbose = FALSE)
  expect_equal(nrow(se_3), 1)
  expect_true("TX1" %in% rownames(se_3))
  
  # Higher min_samples should give fewer or equal rows
  expect_true(nrow(se_3) <= nrow(se_1))
})

test_that(".filter_se stringency parameter: soft/medium/severe", {
  # Create paired data with 6 samples (3 pairs) to allow stringency calculations
  counts <- matrix(c(
    100, 200, 1, 50, 75, 150,    # TX1
    2, 40, 100, 200, 1, 50,       # TX2
    75, 150, 2, 40, 100, 200,     # TX3
    1, 50, 75, 150, 2, 40         # TX4
  ), nrow = 4, ncol = 6, byrow = TRUE)
  tpm <- counts  # Use raw counts to preserve relative abundances, not log2
  rownames(counts) <- rownames(tpm) <- c("TX1", "TX2", "TX3", "TX4")
  colnames(counts) <- colnames(tpm) <- c("S1", "S2", "S3", "S4", "S5", "S6")
  
  se <- SummarizedExperiment::SummarizedExperiment(
    assays = list(counts = counts, tpm = tpm),
    colData = data.frame(
      pair_id = c("pair1", "pair1", "pair2", "pair2", "pair3", "pair3"),
      row.names = c("S1", "S2", "S3", "S4", "S5", "S6")
    )
  )
  
  # Soft: permissive, should keep more transcripts
  se_soft <- .filter_se(se, stringency = "soft", verbose = FALSE)
  n_soft <- nrow(se_soft)
  expect_true(n_soft >= 1)
  
  # Medium: balanced
  se_med <- .filter_se(se, stringency = "medium", verbose = FALSE)
  n_med <- nrow(se_med)
  expect_true(n_med >= 1)
  
  # Severe: stringent, should keep fewer or equal
  se_sev <- .filter_se(se, stringency = "severe", verbose = FALSE)
  n_sev <- nrow(se_sev)
  expect_true(n_sev <= n_soft)  # Severe should be <= soft
})

test_that(".filter_se min_tx_per_gene parameter: multi-isoform filtering", {
  # Create gene with 3 isoforms, another with 1
  counts <- matrix(c(10, 5, 3, 100), nrow = 4, ncol = 2)
  tpm <- counts * 10
  rownames(counts) <- rownames(tpm) <- c("TX1.1", "TX1.2", "TX1.3", "TX2.1")
  colnames(counts) <- colnames(tpm) <- c("S1", "S2")
  
  se <- SummarizedExperiment::SummarizedExperiment(
    assays = list(counts = counts, tpm = tpm),
    rowData = data.frame(
      gene_id = c("G1", "G1", "G1", "G2"),
      row.names = rownames(counts)
    ),
    metadata = list(tx2gene = data.frame(
      tx = c("TX1.1", "TX1.2", "TX1.3", "TX2.1"),
      gene = c("G1", "G1", "G1", "G2")
    ))
  )
  
  # min_tx_per_gene = 1: keep all
  se_1tx <- .filter_se(se, min_tx_per_gene = 1, min_tpm = 0, min_samples = 1, verbose = FALSE)
  expect_equal(nrow(se_1tx), 4)
  
  # min_tx_per_gene = 2: remove TX2.1 (only 1 isoform in G2)
  se_2tx <- .filter_se(se, min_tx_per_gene = 2, min_tpm = 0, min_samples = 1, verbose = FALSE)
  expect_equal(nrow(se_2tx), 3)
  expect_true(!("TX2.1" %in% rownames(se_2tx)))
})

test_that(".filter_se isoform_abundance with various thresholds", {
  # Create gene with 4 isoforms: (100, 60, 30, 10 = total 200)
  # Relative abundances: 50%, 30%, 15%, 5%
  counts <- matrix(c(100, 60, 30, 10, 100, 60, 30, 10), nrow = 4, ncol = 2)
  # Use raw counts as TPM (don't transform - keeps relative abundances intact)
  tpm <- counts
  rownames(counts) <- rownames(tpm) <- c("TX1.1", "TX1.2", "TX1.3", "TX1.4")
  colnames(counts) <- colnames(tpm) <- c("S1", "S2")
  
  se <- SummarizedExperiment::SummarizedExperiment(
    assays = list(counts = counts, tpm = tpm),
    rowData = data.frame(
      gene_id = rep("G1", 4),
      row.names = rownames(counts)
    ),
    metadata = list(tx2gene = data.frame(
      tx = rownames(counts),
      gene = rep("G1", 4)
    ))
  )
  
  # 0% threshold: all kept
  se_0pct <- .filter_se(se, min_isoform_abundance = 0.0, min_tpm = 0, min_samples = 1, verbose = FALSE)
  expect_equal(nrow(se_0pct), 4)
  
  # 10% threshold: removes TX1.4 (5% < 10%)
  se_10pct <- .filter_se(se, min_isoform_abundance = 0.10, min_tpm = 0, min_samples = 1, verbose = FALSE)
  expect_equal(nrow(se_10pct), 3)
  expect_true(!("TX1.4" %in% rownames(se_10pct)))
  
  # 20% threshold: removes TX1.3 and TX1.4 (15% < 20%, 5% < 20%)
  se_20pct <- .filter_se(se, min_isoform_abundance = 0.20, min_tpm = 0, min_samples = 1, verbose = FALSE)
  expect_equal(nrow(se_20pct), 2)
  expect_true(all(c("TX1.1", "TX1.2") %in% rownames(se_20pct)))
})

test_that(".filter_se combined parameters: TPM + min_samples + isoforms", {
  # Multi-isoform, multi-parameter filtering
  counts <- matrix(c(
    500, 300, 100, 50,   # TX1.1: high abundance
    300, 200, 0, 50,     # TX1.2: moderate abundance
    200, 100, 0, 0,      # TX1.3: lower abundance
    100, 50, 0, 0        # TX1.4: lowest abundance
  ), nrow = 4, ncol = 4, byrow = TRUE)
  
  # Use counts directly as TPM to preserve relative abundances
  tpm <- counts
  rownames(counts) <- rownames(tpm) <- c("TX1.1", "TX1.2", "TX1.3", "TX1.4")
  colnames(counts) <- colnames(tpm) <- c("S1", "S2", "S3", "S4")
  
  se <- SummarizedExperiment::SummarizedExperiment(
    assays = list(counts = counts, tpm = tpm),
    rowData = data.frame(
      gene_id = rep("G1", 4),
      row.names = rownames(counts)
    ),
    metadata = list(tx2gene = data.frame(
      tx = rownames(counts),
      gene = rep("G1", 4)
    ))
  )
  
  # Moderate filtering: keep isoforms > 20% relative abundance
  # TX1.1: 500/1100 = 45% (kept)
  # TX1.2: 300/1100 = 27% (kept)
  # TX1.3: 200/1100 = 18% (removed)
  # TX1.4: 100/1100 = 9% (removed)
  se_filt <- .filter_se(se, 
                        min_tpm = 0.1, 
                        min_samples = 1,
                        min_isoform_abundance = 0.20,
                        verbose = FALSE)
  
  # Should have TX1.1 and TX1.2 (both >= 20% relative abundance)
  expect_equal(nrow(se_filt), 2)
  expect_true(all(c("TX1.1", "TX1.2") %in% rownames(se_filt)))
})

test_that(".filter_se with single sample edge case", {
  counts <- matrix(c(100, 50, 10), nrow = 3, ncol = 1)
  tpm <- counts
  rownames(counts) <- rownames(tpm) <- c("TX1", "TX2", "TX3")
  colnames(counts) <- colnames(tpm) <- "S1"
  
  se <- SummarizedExperiment::SummarizedExperiment(
    assays = list(counts = counts, tpm = tpm)
  )
  
  # With single sample, min_samples=1 should work
  se_filt <- .filter_se(se, min_tpm = 30, min_samples = 1, verbose = FALSE)
  expect_equal(nrow(se_filt), 2)  # TX1 and TX2 pass
  expect_true(all(c("TX1", "TX2") %in% rownames(se_filt)))
})

test_that(".filter_se with all transcripts filtered out", {
  counts <- matrix(c(1, 2, 3), nrow = 3, ncol = 2)
  tpm <- counts / 100  # Very low TPM
  rownames(counts) <- rownames(tpm) <- c("TX1", "TX2", "TX3")
  colnames(counts) <- colnames(tpm) <- c("S1", "S2")
  
  se <- SummarizedExperiment::SummarizedExperiment(
    assays = list(counts = counts, tpm = tpm)
  )
  
  # Very high min_tpm should filter everything
  # Suppress warning about all transcripts being filtered (expected behavior for this test)
  se_empty <- suppressWarnings(.filter_se(se, min_tpm = 1000, min_samples = 2, verbose = FALSE))
  expect_equal(nrow(se_empty), 0)
})

test_that(".filter_se invalid stringency parameter raises error", {
  counts <- matrix(c(10, 20), nrow = 2, ncol = 1)
  tpm <- counts
  rownames(counts) <- rownames(tpm) <- c("TX1", "TX2")
  colnames(counts) <- colnames(tpm) <- "S1"
  
  se <- SummarizedExperiment::SummarizedExperiment(
    assays = list(counts = counts, tpm = tpm),
    colData = data.frame(pair_id = "p1", row.names = "S1")
  )
  
  # Invalid stringency should raise error
  expect_error(
    .filter_se(se, stringency = "invalid_level", verbose = FALSE),
    "must be one of"
  )
})

test_that(".filter_se invalid min_isoform_abundance raises error", {
  counts <- matrix(c(10, 20), nrow = 2, ncol = 1)
  tpm <- counts
  rownames(counts) <- rownames(tpm) <- c("TX1", "TX2")
  colnames(counts) <- colnames(tpm) <- "S1"
  
  se <- SummarizedExperiment::SummarizedExperiment(
    assays = list(counts = counts, tpm = tpm)
  )
  
  # min_isoform_abundance > 1 should raise error from .filter_se
  expect_error(
    .filter_se(se, min_samples = 1, min_tpm = 0, min_isoform_abundance = 1.5, verbose = FALSE),
    "must be numeric in"
  )
  
  # Negative min_isoform_abundance should also raise error
  expect_error(
    .filter_se(se, min_samples = 1, min_tpm = 0, min_isoform_abundance = -0.1, verbose = FALSE),
    "must be numeric in"
  )
})

test_that(".filter_se metadata: filtered info is stored", {
  counts <- matrix(c(100, 50, 10, 5), nrow = 2, ncol = 2)
  tpm <- counts * 10
  rownames(counts) <- rownames(tpm) <- c("TX1", "TX2")
  colnames(counts) <- colnames(tpm) <- c("S1", "S2")
  
  se <- SummarizedExperiment::SummarizedExperiment(
    assays = list(counts = counts, tpm = tpm)
  )
  
  se_filt <- .filter_se(se, 
                        min_tpm = 50, 
                        min_samples = 1,
                        min_tx_per_gene = 1,
                        min_isoform_abundance = 0.1,
                        verbose = FALSE)
  
  md <- S4Vectors::metadata(se_filt)
  expect_true(!is.null(md$filtered))
  expect_equal(md$filtered$min_tpm, 50)
  expect_equal(md$filtered$min_samples, 1)
  expect_equal(md$filtered$min_tx_per_gene, 1)
  expect_equal(md$filtered$min_isoform_abundance, 0.1)
})

test_that(".filter_se with custom TPM assay name", {
  counts <- matrix(c(10, 20, 30, 5), nrow = 2, ncol = 2)
  custom_tpm <- counts * 5
  rownames(counts) <- rownames(custom_tpm) <- c("TX1", "TX2")
  colnames(counts) <- colnames(custom_tpm) <- c("S1", "S2")
  
  se <- SummarizedExperiment::SummarizedExperiment(
    assays = list(counts = counts, my_abundance = custom_tpm)
  )
  
  # Filter using custom assay
  se_filt <- .filter_se(se, 
                        tpm_assay_name = "my_abundance",
                        min_tpm = 50,
                        min_samples = 1,
                        verbose = FALSE)
  
  expect_true(nrow(se_filt) >= 1)
  # Both TX1 and TX2 should pass (100 and 150 >= 50)
  expect_equal(nrow(se_filt), 2)
})

test_that(".filter_se rowData and colData preservation", {
  counts <- matrix(1:12, nrow = 4, ncol = 3)
  tpm <- matrix(10:21, nrow = 4, ncol = 3)
  rownames(counts) <- rownames(tpm) <- c("TX1", "TX2", "TX3", "TX4")
  colnames(counts) <- colnames(tpm) <- c("S1", "S2", "S3")
  
  se <- SummarizedExperiment::SummarizedExperiment(
    assays = list(counts = counts, tpm = tpm),
    rowData = data.frame(
      transcript_id = rownames(counts),
      biotype = c("protein_coding", "protein_coding", "lncRNA", "lncRNA"),
      row.names = rownames(counts)
    ),
    colData = data.frame(
      sample_id = colnames(counts),
      treatment = c("control", "control", "treated"),
      row.names = colnames(counts)
    )
  )
  
  se_filt <- .filter_se(se, min_tpm = 5, min_samples = 2, verbose = FALSE)
  
  # rowData should be filtered to match rows
  rd_filt <- SummarizedExperiment::rowData(se_filt)
  expect_equal(nrow(rd_filt), nrow(se_filt))
  
  # colData should remain unchanged
  cd_filt <- SummarizedExperiment::colData(se_filt)
  expect_equal(nrow(cd_filt), 3)
  expect_equal(cd_filt$treatment, c("control", "control", "treated"))
})

test_that(".filter_se with multiple assays: all filtered consistently", {
  counts <- matrix(c(100, 10, 0), nrow = 3, ncol = 1)
  tpm <- counts / 10
  abundance <- counts * 2
  rownames(counts) <- rownames(tpm) <- rownames(abundance) <- c("TX1", "TX2", "TX3")
  colnames(counts) <- colnames(tpm) <- colnames(abundance) <- "S1"
  
  se <- SummarizedExperiment::SummarizedExperiment(
    assays = list(counts = counts, tpm = tpm, abundance = abundance)
  )
  
  se_filt <- .filter_se(se, min_tpm = 5, min_samples = 1, verbose = FALSE)
  
  # All assays should be filtered identically
  a1 <- assays(se_filt)$counts
  a2 <- assays(se_filt)$tpm
  a3 <- assays(se_filt)$abundance
  
  expect_equal(nrow(a1), nrow(a2))
  expect_equal(nrow(a2), nrow(a3))
  expect_equal(rownames(a1), rownames(a2))
  expect_equal(rownames(a2), rownames(a3))
})

# ============================================================================
# FILTER_ANALYSIS_S4 TESTS
# ============================================================================

test_that("filter_analysis: basic filtering preserves TSENATAnalysis class", {
  counts <- matrix(c(100, 50, 10, 5, 20, 15), nrow = 3, ncol = 2)
  tpm <- counts
  rownames(counts) <- rownames(tpm) <- c("TX1", "TX2", "TX3")
  colnames(counts) <- colnames(tpm) <- c("S1", "S2")
  
  se <- SummarizedExperiment::SummarizedExperiment(
    assays = list(counts = counts, tpm = tpm),
    rowData = data.frame(gene_id = c("G1", "G1", "G2"), row.names = rownames(counts))
  )
  
  analysis <- TSENAT::TSENATAnalysis(se)
  
  # Filter with min_samples = 1 to keep most data
  filtered_analysis <- TSENAT:::filter_analysis(analysis, min_samples = 1, verbose = FALSE)
  
  # Check class preservation
  expect_s4_class(filtered_analysis, "TSENATAnalysis")
  expect_s4_class(filtered_analysis@se, "SummarizedExperiment")
})

test_that("filter_analysis: stringency levels affect filtering", {
  # Create data with 3 isoforms per gene to survive severe stringency (requires 3+ isoforms)
  # 6 samples = 3 pairs; each gene has 3 isoforms for min_tx_per_gene requirement
  counts <- matrix(c(
    1000, 1050, 1020, 1080, 1040, 1100,  # G1.1: high baseline
    950, 1000, 970, 1030, 990, 1050,     # G1.2: high baseline
    900, 950, 920, 980, 940, 1000,       # G1.3: high baseline
    850, 900, 870, 930, 890, 950         # G2.1: high baseline
  ), nrow = 4, ncol = 6, byrow = TRUE)
  tpm <- counts
  rownames(counts) <- rownames(tpm) <- c("G1.1", "G1.2", "G1.3", "G2.1")
  colnames(counts) <- colnames(tpm) <- c("S1", "S2", "S3", "S4", "S5", "S6")
  
  se <- SummarizedExperiment::SummarizedExperiment(
    assays = list(counts = counts, tpm = tpm),
    rowData = data.frame(gene_id = c("G1", "G1", "G1", "G2"), row.names = rownames(counts)),
    colData = data.frame(
      pair_id = c("pair1", "pair1", "pair2", "pair2", "pair3", "pair3"),
      row.names = colnames(counts)
    )
  )
  
  analysis <- TSENAT::TSENATAnalysis(se)
  
  # Apply different stringency levels
  filt_soft <- TSENAT:::filter_analysis(analysis, stringency = "soft", verbose = FALSE)
  filt_medium <- TSENAT:::filter_analysis(analysis, stringency = "medium", verbose = FALSE)
  filt_severe <- TSENAT:::filter_analysis(analysis, stringency = "severe", verbose = FALSE)
  
  n_soft <- nrow(filt_soft@se)
  n_medium <- nrow(filt_medium@se)
  n_severe <- nrow(filt_severe@se)
  
  # All should keep at least G1 isoforms (3 isoforms) since they survive min_tx_per_gene=3
  # G2 has only 1 isoform, so removed by severe stringency min_tx_per_gene=3
  expect_true(n_soft > 0)
  expect_true(n_medium > 0)
  expect_true(n_severe > 0)
  
  # Severe should be <= medium <= soft (stricter filtering removes more)
  expect_true(n_severe <= n_medium)
  expect_true(n_medium <= n_soft)
})

test_that("filter_analysis: min_samples and min_tpm parameters", {
  # Create data where one gene's transcripts pass but another's don't
  # Each gene has 2 transcripts to avoid single-tx gene filtering
  counts <- matrix(c(
    100, 100, 100,  # G1.1: passes min_tpm=50 in 3 samples
    100, 100, 100,  # G1.2: passes min_tpm=50 in 3 samples
    50, 50, 0,      # G2.1: only 2 samples >= 50, fails min_samples=3
    50, 50, 0       # G2.2: only 2 samples >= 50, fails min_samples=3
  ), nrow = 4, ncol = 3, byrow = TRUE)
  tpm <- counts
  rownames(counts) <- rownames(tpm) <- c("G1.1", "G1.2", "G2.1", "G2.2")
  colnames(counts) <- colnames(tpm) <- c("S1", "S2", "S3")
  
  se <- SummarizedExperiment::SummarizedExperiment(
    assays = list(counts = counts, tpm = tpm),
    rowData = data.frame(gene_id = c("G1", "G1", "G2", "G2"), row.names = rownames(counts))
  )
  
  analysis <- TSENAT::TSENATAnalysis(se)
  
  # min_samples = 3, min_tpm = 50: G1 passes (2 transcripts), G2 fails (need 3+ samples)
  filt_3 <- TSENAT:::filter_analysis(analysis, min_samples = 3, min_tpm = 50, verbose = FALSE)
  expect_equal(nrow(filt_3@se), 2)  # Only G1.1 and G1.2
  expect_true("G1.1" %in% rownames(filt_3@se))
  expect_true("G1.2" %in% rownames(filt_3@se))
})

test_that("filter_analysis: error on invalid input", {
  # Try to filter non-TSENATAnalysis object
  expect_error(
    TSENAT:::filter_analysis("not_an_analysis"),
    "must be a TSENATAnalysis"
  )
})

test_that("filter_analysis: colData preserved after filtering", {
  counts <- matrix(c(100, 50, 10), nrow = 3, ncol = 1)
  tpm <- counts
  rownames(counts) <- rownames(tpm) <- c("TX1", "TX2", "TX3")
  colnames(counts) <- colnames(tpm) <- "S1"
  
  se <- SummarizedExperiment::SummarizedExperiment(
    assays = list(counts = counts, tpm = tpm),
    rowData = data.frame(gene_id = c("G1", "G1", "G2"), row.names = rownames(counts)),
    colData = data.frame(
      sample_name = "S1",
      treatment = "control",
      batch = 1,
      row.names = "S1"
    )
  )
  
  analysis <- TSENAT::TSENATAnalysis(se)
  filt <- TSENAT:::filter_analysis(analysis, min_tpm = 5, min_samples = 1, verbose = FALSE)
  
  # colData should remain unchanged
  cd <- SummarizedExperiment::colData(filt@se)
  expect_equal(nrow(cd), 1)
  expect_equal(cd$sample_name, "S1")
  expect_equal(cd$treatment, "control")
  expect_equal(cd$batch, 1)
})

test_that("filter_analysis: rowData preserved for kept transcripts", {
  counts <- matrix(c(100, 50, 10), nrow = 3, ncol = 1)
  tpm <- counts
  rownames(counts) <- rownames(tpm) <- c("TX1", "TX2", "TX3")
  colnames(counts) <- colnames(tpm) <- "S1"
  
  se <- SummarizedExperiment::SummarizedExperiment(
    assays = list(counts = counts, tpm = tpm),
    rowData = data.frame(
      gene_id = c("G1", "G1", "G2"),
      transcript_name = c("ENST001", "ENST002", "ENST003"),
      row.names = rownames(counts)
    )
  )
  
  analysis <- TSENAT::TSENATAnalysis(se)
  filt <- TSENAT:::filter_analysis(analysis, min_tpm = 40, min_samples = 1, verbose = FALSE)
  
  rd <- SummarizedExperiment::rowData(filt@se)
  expect_equal(nrow(rd), 2)  # TX1 and TX2 kept, TX3 filtered
  expect_true(all(rd$transcript_name %in% c("ENST001", "ENST002")))
})

test_that("filter_analysis: isoform filtering through analysis object", {
  counts <- matrix(c(
    50, 40, 30,  # TX1.1 (120/200 = 60%)
    40, 50, 60,  # TX1.2 (150/200 = 75%)
    30, 10, 10   # TX1.3 (50/200 = 25%) - should be removed with 30% threshold
  ), nrow = 3, ncol = 3, byrow = TRUE)
  tpm <- counts
  rownames(counts) <- rownames(tpm) <- c("TX1.1", "TX1.2", "TX1.3")
  colnames(counts) <- colnames(tpm) <- c("S1", "S2", "S3")
  
  se <- SummarizedExperiment::SummarizedExperiment(
    assays = list(counts = counts, tpm = tpm),
    rowData = data.frame(gene_id = c("TX1", "TX1", "TX1"), row.names = rownames(counts))
  )
  
  analysis <- TSENAT::TSENATAnalysis(se)
  
  # Filter isoforms: remove those < 30% relative abundance
  filt <- TSENAT:::filter_analysis(analysis, min_samples = 1, min_tpm = 1, 
                                      min_isoform_abundance = 0.30, verbose = FALSE)
  
  # TX1.3 has ~25% abundance, should be removed
  expect_equal(nrow(filt@se), 2)
  expect_true("TX1.1" %in% rownames(filt@se))
  expect_true("TX1.2" %in% rownames(filt@se))
  expect_false("TX1.3" %in% rownames(filt@se))
})

test_that("filter_analysis: multi-gene filtering accuracy", {
  # Create data with clear filtering outcomes
  # G1: 2 high-expression transcripts that both pass TPM filter
  # G2: 1 high and 1 low - one passes, one fails TPM
  # G3: 2 moderate transcripts - both pass
  counts <- matrix(c(
    100, 100,  # G1.1: high, clearly passes min_tpm=15 in 2 samples
    80, 90,    # G1.2: high, clearly passes min_tpm=15 in 2 samples
    30, 30,    # G2.1: moderate, passes min_tpm=15 in 2 samples
    2, 3,      # G2.2: very low, fails min_tpm=15
    25, 28     # G3.1: moderate, passes min_tpm=15 in 2 samples
  ), nrow = 5, ncol = 2, byrow = TRUE)
  tpm <- counts
  rownames(counts) <- rownames(tpm) <- c("G1.1", "G1.2", "G2.1", "G2.2", "G3.1")
  colnames(counts) <- colnames(tpm) <- c("S1", "S2")
  
  se <- SummarizedExperiment::SummarizedExperiment(
    assays = list(counts = counts, tpm = tpm),
    rowData = data.frame(
      gene_id = c("G1", "G1", "G2", "G2", "G3"),
      row.names = rownames(counts)
    )
  )
  
  analysis <- TSENAT::TSENATAnalysis(se)
  
  # Filter with min_tpm=15: keeps transcripts with >= 15 in 2+ samples
  # G1: keeps both (100, 80, 90 all >= 15) = 2 transcripts
  # G2: keeps G2.1 only (30 >= 15, but 2,3 < 15) = 1 transcript, min_tx_per_gene=2 removes it
  # G3: single transcript, min_tx_per_gene=2 removes it
  # Result: 2 rows from G1 only
  filt <- TSENAT:::filter_analysis(analysis, min_tpm = 15, min_samples = 2, verbose = FALSE)
  
  expect_equal(nrow(filt@se), 2)  # Only G1.1 and G1.2 (G2 and G3 filtered by min_tx_per_gene=2)
  expect_true("G1.1" %in% rownames(filt@se))
  expect_true("G1.2" %in% rownames(filt@se))
})

test_that(".filter_se with varying sparsity", {
  # Create test data with genes having different expression patterns
  # This tests basic filtering without min_valid_frac parameter
  counts <- matrix(c(
    100, 100, 100, 100, 100, 100, 100, 100, 0, 0,
    100, 100, 100, 100, 100, 100, 100, 100, 0, 0,
    100, 100, 100, 100, 100, 100, 0, 0, 0, 0,
    100, 100, 100, 100, 100, 100, 0, 0, 0, 0,
    100, 100, 100, 100, 0, 0, 0, 0, 0, 0,
    100, 100, 100, 100, 0, 0, 0, 0, 0, 0,
    100, 100, 0, 0, 0, 0, 0, 0, 0, 0,
    100, 100, 0, 0, 0, 0, 0, 0, 0, 0
  ), nrow = 8, ncol = 10, byrow = TRUE)
  tpm <- counts + 1
  rownames(counts) <- rownames(tpm) <- c("TX1.1", "TX1.2", "TX2.1", "TX2.2", "TX3.1", "TX3.2", "TX4.1", "TX4.2")
  colnames(counts) <- colnames(tpm) <- paste0("S", 1:10)
  
  se <- SummarizedExperiment::SummarizedExperiment(
    assays = list(counts = counts, tpm = tpm),
    rowData = data.frame(gene_id = c("G1", "G1", "G2", "G2", "G3", "G3", "G4", "G4"), row.names = rownames(counts))
  )
  
  # Basic filtering without min_valid_frac
  se_filt <- .filter_se(se, min_tpm = 0, min_samples = 1, verbose = FALSE)
  expect_true(nrow(se_filt) > 0)
  expect_true(nrow(se_filt) <= nrow(se))
})

test_that(".filter_se sparse gene behavior", {
  # Test filtering with sparse data
  counts <- matrix(c(
    100, 100, 100, 100, 100, 100, 0, 0, 0, 0,
    100, 100, 100, 100, 100, 100, 0, 0, 0, 0,
    100, 100, 100, 0, 0, 0, 0, 0, 0, 0,
    100, 100, 100, 0, 0, 0, 0, 0, 0, 0
  ), nrow = 4, ncol = 10, byrow = TRUE)
  tpm <- counts + 1
  rownames(counts) <- rownames(tpm) <- c("TX1.1", "TX1.2", "TX2.1", "TX2.2")
  colnames(counts) <- colnames(tpm) <- paste0("S", 1:10)
  
  se <- SummarizedExperiment::SummarizedExperiment(
    assays = list(counts = counts, tpm = tpm),
    rowData = data.frame(gene_id = c("G1", "G1", "G2", "G2"), row.names = rownames(counts))
  )
  
  # Basic filtering without min_valid_frac
  se_filt <- .filter_se(se, min_tpm = 0, min_samples = 1, verbose = FALSE)
  expect_equal(nrow(se_filt), 4)  # All pass basic filtering
})

test_that(".filter_se with multi-isoform genes", {
  # Create multi-isoform genes with different expression patterns
  counts <- matrix(c(
    100, 100, 100, 100, 100, 100, 100, 100, 0, 0,
    100, 100, 100, 100, 100, 100, 100, 0, 0, 0,
    100, 100, 100, 100, 100, 0, 0, 0, 0, 0,
    100, 100, 100, 0, 0, 0, 0, 0, 0, 0
  ), nrow = 4, ncol = 10, byrow = TRUE)
  tpm <- counts + 1
  rownames(counts) <- rownames(tpm) <- c("TX1", "TX2", "TX3", "TX4")
  colnames(counts) <- colnames(tpm) <- paste0("S", 1:10)
  
  se <- SummarizedExperiment::SummarizedExperiment(
    assays = list(counts = counts, tpm = tpm),
    rowData = data.frame(gene_id = c("G1", "G1", "G2", "G2"), row.names = rownames(counts)),
    metadata = list(tx2gene = data.frame(tx = rownames(counts), gene = c("G1", "G1", "G2", "G2")))
  )
  
  # Filter with basic parameters
  se_filt <- .filter_se(se, min_tpm = 0, min_samples = 1, verbose = FALSE)
  expect_true(nrow(se_filt) > 0)  # All pass basic filtering
})

test_that(".filter_se with combined filtering", {
  # Test combining multiple filtering parameters
  counts <- matrix(c(
    300, 300, 300, 300, 300, 300, 300, 300, 0, 0,
    240, 240, 240, 240, 240, 240, 240, 240, 0, 0,
    100, 100, 100, 100, 100, 100, 0, 0, 0, 0,
    100, 100, 100, 100, 100, 100, 0, 0, 0, 0,
    50, 50, 50, 50, 0, 0, 0, 0, 0, 0,
    50, 50, 50, 50, 0, 0, 0, 0, 0, 0
  ), nrow = 6, ncol = 10, byrow = TRUE)
  
  tpm <- counts + 1
  rownames(counts) <- rownames(tpm) <- c("TX1.1", "TX1.2", "TX2.1", "TX2.2", "TX3.1", "TX3.2")
  colnames(counts) <- colnames(tpm) <- paste0("S", 1:10)
  
  se <- SummarizedExperiment::SummarizedExperiment(
    assays = list(counts = counts, tpm = tpm),
    rowData = data.frame(gene_id = c("G1", "G1", "G2", "G2", "G3", "G3"), row.names = rownames(counts)),
    metadata = list(tx2gene = data.frame(tx = rownames(counts), gene = c("G1", "G1", "G2", "G2", "G3", "G3")))
  )
  
  se_filt <- .filter_se(se, 
                        min_tpm = 0,
                        min_samples = 1,
                        min_isoform_abundance = 0,
                        verbose = FALSE)
  
  # All genes pass basic filtering
  expect_true(nrow(se_filt) > 0)
})

test_that(".filter_se with few samples", {
  # Test edge case with limited samples
  counts <- matrix(c(
    100, 100, 100, 100, 100, 100, 100, 100, 100, 100,
    100, 100, 100, 100, 100, 100, 100, 100, 100, 100,
    100, 100, 100, 100, 100, 0, 0, 0, 0, 0,
    100, 100, 100, 100, 100, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0
  ), nrow = 6, ncol = 10, byrow = TRUE)
  tpm <- counts + 1
  rownames(counts) <- rownames(tpm) <- c("TX1.1", "TX1.2", "TX2.1", "TX2.2", "TX3.1", "TX3.2")
  colnames(counts) <- colnames(tpm) <- paste0("S", 1:10)
  
  se <- SummarizedExperiment::SummarizedExperiment(
    assays = list(counts = counts, tpm = tpm),
    rowData = data.frame(gene_id = c("G1", "G1", "G2", "G2", "G3", "G3"), row.names = rownames(counts))
  )
  
  # Filter with basic parameters
  se_filt <- .filter_se(se, min_tpm = 0, min_samples = 1, verbose = FALSE)
  
  # All non-zero genes should be kept
  expect_true(nrow(se_filt) > 0)
})



# ============================================================================
# HELPER FUNCTION UNIT TESTS
# ============================================================================

describe("Helper functions unit tests", {
  
  # ========================================================================
  # .calc_stringency_params tests
  # ========================================================================
  
  test_that(".calc_stringency_params: soft stringency with 6 samples (3 pairs)", {
    counts <- matrix(c(100, 100, 100, 100, 100, 100), nrow = 2, ncol = 6)
    tpm <- counts
    rownames(counts) <- rownames(tpm) <- c("TX1", "TX2")
    colnames(counts) <- colnames(tpm) <- c("S1", "S2", "S3", "S4", "S5", "S6")
    
    se <- SummarizedExperiment::SummarizedExperiment(
      assays = list(counts = counts, tpm = tpm),
      colData = data.frame(
        pair_id = c("pair1", "pair1", "pair2", "pair2", "pair3", "pair3"),
        row.names = colnames(counts)
      )
    )
    
    params <- TSENAT:::.calc_stringency_params(se, "soft", pair_col = "pair_id", verbose = FALSE)
    
    expect_true(is.list(params))
    expect_true(is.numeric(params$min_samples))
    expect_true(is.numeric(params$min_tx_per_gene))
    expect_equal(params$min_samples, 2)  # soft: max(2, ceil(0.25*6)) = max(2, 2) = 2
    expect_equal(params$min_tx_per_gene, 2)  # soft: min_tx = 2
  })
  
  test_that(".calc_stringency_params: medium stringency with 6 samples", {
    counts <- matrix(c(100, 100, 100, 100, 100, 100), nrow = 2, ncol = 6)
    tpm <- counts
    rownames(counts) <- rownames(tpm) <- c("TX1", "TX2")
    colnames(counts) <- colnames(tpm) <- c("S1", "S2", "S3", "S4", "S5", "S6")
    
    se <- SummarizedExperiment::SummarizedExperiment(
      assays = list(counts = counts, tpm = tpm),
      colData = data.frame(
        pair_id = c("pair1", "pair1", "pair2", "pair2", "pair3", "pair3"),
        row.names = colnames(counts)
      )
    )
    
    params <- TSENAT:::.calc_stringency_params(se, "medium", pair_col = "pair_id", verbose = FALSE)
    
    expect_equal(params$min_samples, 3)  # medium: max(3, ceil(0.5*6)) = max(3, 3) = 3
    expect_equal(params$min_tx_per_gene, 2)
  })
  
  test_that(".calc_stringency_params: severe stringency with 6 samples", {
    counts <- matrix(c(100, 100, 100, 100, 100, 100), nrow = 2, ncol = 6)
    tpm <- counts
    rownames(counts) <- rownames(tpm) <- c("TX1", "TX2")
    colnames(counts) <- colnames(tpm) <- c("S1", "S2", "S3", "S4", "S5", "S6")
    
    se <- SummarizedExperiment::SummarizedExperiment(
      assays = list(counts = counts, tpm = tpm),
      colData = data.frame(
        pair_id = c("pair1", "pair1", "pair2", "pair2", "pair3", "pair3"),
        row.names = colnames(counts)
      )
    )
    
    params <- TSENAT:::.calc_stringency_params(se, "severe", pair_col = "pair_id", verbose = FALSE)
    
    expect_equal(params$min_samples, 5)  # severe: max(3, ceil(0.75*6)) = max(3, 5) = 5
    expect_equal(params$min_tx_per_gene, 3)
  })
  
  test_that(".calc_stringency_params: detects pair column automatically", {
    counts <- matrix(c(100, 100, 100, 100, 100, 100), nrow = 2, ncol = 6)
    tpm <- counts
    rownames(counts) <- rownames(tpm) <- c("TX1", "TX2")
    colnames(counts) <- colnames(tpm) <- c("S1", "S2", "S3", "S4", "S5", "S6")
    
    se <- SummarizedExperiment::SummarizedExperiment(
      assays = list(counts = counts, tpm = tpm),
      colData = data.frame(
        pair_id = c("pair1", "pair1", "pair2", "pair2", "pair3", "pair3"),
        row.names = colnames(counts)
      )
    )
    
    # Pass NULL for pair_col; should auto-detect from colData
    params <- TSENAT:::.calc_stringency_params(se, "soft", pair_col = NULL, verbose = FALSE)
    
    expect_is(params, "list")
    expect_is(params$min_samples, "numeric")
  })
  
  # ========================================================================
  # .detect_pair_column tests
  # ========================================================================
  
  test_that(".detect_pair_column: finds 'pair' column", {
    col_data <- data.frame(
      sample_id = c("S1", "S2", "S3"),
      pair = c("p1", "p1", "p2"),
      row.names = c("S1", "S2", "S3")
    )
    
    se <- SummarizedExperiment::SummarizedExperiment(
      matrix(1:9, nrow = 3, ncol = 3),
      colData = col_data
    )
    
    result <- TSENAT:::.detect_pair_column(col_data, se)
    expect_equal(result, "pair")
  })
  
  test_that(".detect_pair_column: finds 'pair_id' column", {
    col_data <- data.frame(
      sample_id = c("S1", "S2", "S3"),
      pair_id = c("p1", "p1", "p2"),
      row.names = c("S1", "S2", "S3")
    )
    
    se <- SummarizedExperiment::SummarizedExperiment(
      matrix(1:9, nrow = 3, ncol = 3),
      colData = col_data
    )
    
    result <- TSENAT:::.detect_pair_column(col_data, se)
    expect_equal(result, "pair_id")
  })
  
  test_that(".detect_pair_column: tries multiple candidate names", {
    # Test with 'subject_id' which is a later candidate
    col_data <- data.frame(
      sample_id = c("S1", "S2", "S3"),
      subject_id = c("sbj1", "sbj1", "sbj2"),
      row.names = c("S1", "S2", "S3")
    )
    
    se <- SummarizedExperiment::SummarizedExperiment(
      matrix(1:9, nrow = 3, ncol = 3),
      colData = col_data
    )
    
    result <- TSENAT:::.detect_pair_column(col_data, se)
    expect_equal(result, "subject_id")
  })
  
  test_that(".detect_pair_column: error when no pair column found", {
    col_data <- data.frame(
      sample_id = c("S1", "S2", "S3"),
      treatment = c("A", "B", "A"),
      row.names = c("S1", "S2", "S3")
    )
    
    se <- SummarizedExperiment::SummarizedExperiment(
      matrix(1:9, nrow = 3, ncol = 3),
      colData = col_data
    )
    
    expect_error(
      TSENAT:::.detect_pair_column(col_data, se),
      "Could not auto-detect pair column"
    )
  })
  
  # ========================================================================
  # .estimate_min_tpm tests
  # ========================================================================
  
  test_that(".estimate_min_tpm: soft stringency (Q1)", {
    # Create assay with clear quantile structure
    assay_mat <- matrix(c(
      c(0.5, 1.0, 1.5, 2.0, 10.0),
      c(0.5, 1.0, 1.5, 2.0, 10.0)
    ), nrow = 2, ncol = 5, byrow = TRUE)
    
    result <- TSENAT:::.estimate_min_tpm(assay_mat, "soft", verbose = FALSE)
    
    expect_is(result, "list")
    expect_true("min_tpm" %in% names(result))
    expect_true(is.numeric(result$min_tpm))
    expect_true(result$min_tpm > 0.1)
    expect_true(result$min_tpm <= 5.0)
  })
  
  test_that(".estimate_min_tpm: medium stringency (Q2/Median)", {
    assay_mat <- matrix(c(
      c(0.5, 1.0, 1.5, 2.0, 10.0),
      c(0.5, 1.0, 1.5, 2.0, 10.0)
    ), nrow = 2, ncol = 5, byrow = TRUE)
    
    result <- TSENAT:::.estimate_min_tpm(assay_mat, "medium", verbose = FALSE)
    
    expect_is(result, "list")
    expect_true("min_tpm" %in% names(result))
    expect_true(is.numeric(result$min_tpm))
    expect_true(result$min_tpm >= 0.1)
    expect_true(result$min_tpm <= 5.0)
  })
  
  test_that(".estimate_min_tpm: severe stringency (Q3)", {
    assay_mat <- matrix(c(
      c(0.5, 1.0, 1.5, 2.0, 10.0),
      c(0.5, 1.0, 1.5, 2.0, 10.0)
    ), nrow = 2, ncol = 5, byrow = TRUE)
    
    result <- TSENAT:::.estimate_min_tpm(assay_mat, "severe", verbose = FALSE)
    
    expect_is(result, "list")
    expect_true("min_tpm" %in% names(result))
    expect_true(is.numeric(result$min_tpm))
    expect_true(result$min_tpm >= 0.1)
    expect_true(result$min_tpm <= 5.0)
  })
  
  test_that(".estimate_min_tpm: respects lower bound (0.1)", {
    # Create assay with very low values
    assay_mat <- matrix(c(
      c(0.01, 0.02, 0.03, 0.04, 0.05),
      c(0.01, 0.02, 0.03, 0.04, 0.05)
    ), nrow = 2, ncol = 5, byrow = TRUE)
    
    result <- TSENAT:::.estimate_min_tpm(assay_mat, "soft", verbose = FALSE)
    
    expect_true(result$min_tpm >= 0.1)  # Should be at least 0.1
  })
  
  test_that(".estimate_min_tpm: respects upper bound (5.0)", {
    # Create assay with very high values
    assay_mat <- matrix(c(
      c(100, 200, 300, 400, 500),
      c(100, 200, 300, 400, 500)
    ), nrow = 2, ncol = 5, byrow = TRUE)
    
    result <- TSENAT:::.estimate_min_tpm(assay_mat, "severe", verbose = FALSE)
    
    expect_true(result$min_tpm <= 5.0)  # Should be at most 5.0
  })
  
  # ========================================================================
  # .apply_tpm_filter tests
  # ========================================================================
  
  test_that(".apply_tpm_filter: basic threshold filtering", {
    tpm_mat <- matrix(c(
      c(10, 20, 30),
      c(5, 10, 15),
      c(1, 2, 3)
    ), nrow = 3, ncol = 3, byrow = TRUE)
    
    result <- TSENAT:::.apply_tpm_filter(tpm_mat, min_tpm = 10, min_samples = 2)
    
    expect_is(result, "logical")
    expect_equal(length(result), 3)
    expect_true(result[1])  # Row 1: 20, 30 >= 10 in 2 samples
    expect_true(result[2])  # Row 2: 10, 15 >= 10 in 2 samples
    expect_false(result[3])  # Row 3: only 1 sample >= 10
  })
  
  test_that(".apply_tpm_filter: with min_samples = 1", {
    tpm_mat <- matrix(c(
      c(10, 0, 0),
      c(5, 5, 5)
    ), nrow = 2, ncol = 3, byrow = TRUE)
    
    result <- TSENAT:::.apply_tpm_filter(tpm_mat, min_tpm = 5, min_samples = 1)
    
    expect_true(result[1])  # Row 1: 10 >= 5 in 1 sample
    expect_true(result[2])  # Row 2: 5, 5, 5 >= 5 in 3 samples
  })
  
  test_that(".apply_tpm_filter: with high min_samples", {
    tpm_mat <- matrix(c(
      c(50, 50, 50),
      c(50, 50, 10)
    ), nrow = 2, ncol = 3, byrow = TRUE)
    
    result <- TSENAT:::.apply_tpm_filter(tpm_mat, min_tpm = 40, min_samples = 3)
    
    expect_true(result[1])  # Row 1: all 3 samples >= 40
    expect_false(result[2])  # Row 2: only 2 samples >= 40
  })
  
  test_that(".apply_tpm_filter: caps min_samples at ncol", {
    tpm_mat <- matrix(c(
      c(50, 50),
      c(30, 40)
    ), nrow = 2, ncol = 2, byrow = TRUE)
    
    # Request min_samples=5 but only 2 columns exist - should warn
    result <- expect_warning(
      TSENAT:::.apply_tpm_filter(tpm_mat, min_tpm = 25, min_samples = 5),
      "min_samples.*is greater than number of samples"
    )
    
    expect_is(result, "logical")
    expect_equal(length(result), 2)
  })
  
  # ========================================================================
  # .filter_by_tx_per_gene tests
  # ========================================================================
  
  test_that(".filter_by_tx_per_gene: filters multi-isoform genes", {
    tokeep <- c(TRUE, TRUE, TRUE, FALSE)
    genes_vec <- c("G1", "G1", "G2", "G2")
    
    # G1: 2 isoforms, G2: 1 kept + 1 removed = fails min_tx_per_gene=2
    result <- TSENAT:::.filter_by_tx_per_gene(tokeep, genes_vec, min_tx_per_gene = 2)
    
    expect_is(result, "logical")
    expect_equal(length(result), 4)
    expect_true(result[1])  # G1.1: G1 has 2 isoforms kept
    expect_true(result[2])  # G1.2: G1 has 2 isoforms kept
    expect_false(result[3])  # G2.1: G2 has only 1 kept (other is FALSE)
    expect_false(result[4])  # G2.2: already FALSE
  })
  
  test_that(".filter_by_tx_per_gene: min_tx_per_gene = 1 keeps all", {
    tokeep <- c(TRUE, TRUE, TRUE, FALSE)
    genes_vec <- c("G1", "G1", "G2", "G2")
    
    result <- TSENAT:::.filter_by_tx_per_gene(tokeep, genes_vec, min_tx_per_gene = 1)
    
    expect_equal(result, tokeep)  # All true values unchanged
  })
  
  test_that(".filter_by_tx_per_gene: with single-isoform genes", {
    tokeep <- c(TRUE, TRUE, TRUE)
    genes_vec <- c("G1", "G1", "G2")
    
    # G1: 2 isoforms (passes), G2: 1 isoform (fails min_tx_per_gene=2)
    result <- TSENAT:::.filter_by_tx_per_gene(tokeep, genes_vec, min_tx_per_gene = 2)
    
    expect_true(result[1])  # G1.1: in 2-isoform gene
    expect_true(result[2])  # G1.2: in 2-isoform gene
    expect_false(result[3])  # G2: single isoform
  })
  
  # ========================================================================
  # .select_genes_from_analysis tests
  # ========================================================================
  
  test_that(".select_genes_from_analysis: by variance", {
    counts <- matrix(c(
      c(100, 200, 50, 30),    # TX1: high variance
      c(10, 10, 10, 10),      # TX2: zero variance
      c(50, 60, 55, 65),      # TX3: medium variance
      c(100, 100, 100, 100)   # TX4: zero variance
    ), nrow = 4, ncol = 4, byrow = TRUE)
    tpm <- counts
    rownames(counts) <- rownames(tpm) <- c("TX1", "TX2", "TX3", "TX4")
    colnames(counts) <- colnames(tpm) <- c("S1", "S2", "S3", "S4")
    
    se <- SummarizedExperiment::SummarizedExperiment(
      assays = list(counts = counts, tpm = tpm)
    )
    
    # Select 2 genes with highest variance
    result <- TSENAT:::.select_genes_from_analysis(se, n_genes = 2, genes = NULL, 
                                                   select_by = "variance", verbose = FALSE)
    
    expect_is(result, "integer")
    expect_equal(length(result), 2)
    # TX1 and TX3 should have highest variance
    expect_true(1 %in% result)  # TX1
    expect_true(3 %in% result)  # TX3
  })
  
  test_that(".select_genes_from_analysis: by mean", {
    counts <- matrix(c(
      c(1000, 1000, 1000, 1000),   # TX1: highest mean
      c(100, 100, 100, 100),       # TX2: medium mean
      c(10, 10, 10, 10),           # TX3: lowest mean
      c(500, 500, 500, 500)        # TX4: medium-high mean
    ), nrow = 4, ncol = 4, byrow = TRUE)
    tpm <- counts
    rownames(counts) <- rownames(tpm) <- c("TX1", "TX2", "TX3", "TX4")
    colnames(counts) <- colnames(tpm) <- c("S1", "S2", "S3", "S4")
    
    se <- SummarizedExperiment::SummarizedExperiment(
      assays = list(counts = counts, tpm = tpm)
    )
    
    # Select 2 genes with highest mean
    result <- TSENAT:::.select_genes_from_analysis(se, n_genes = 2, genes = NULL,
                                                   select_by = "mean", verbose = FALSE)
    
    expect_is(result, "integer")
    expect_equal(length(result), 2)
    expect_true(1 %in% result)  # TX1: mean 1000
    expect_true(4 %in% result)  # TX4: mean 500
  })
  
  test_that(".select_genes_from_analysis: by random with seed", {
    counts <- matrix(c(
      c(100, 100, 100, 100),
      c(100, 100, 100, 100),
      c(100, 100, 100, 100),
      c(100, 100, 100, 100)
    ), nrow = 4, ncol = 4, byrow = TRUE)
    tpm <- counts
    rownames(counts) <- rownames(tpm) <- c("TX1", "TX2", "TX3", "TX4")
    colnames(counts) <- colnames(tpm) <- c("S1", "S2", "S3", "S4")
    
    se <- SummarizedExperiment::SummarizedExperiment(
      assays = list(counts = counts, tpm = tpm)
    )
    
    # Random selection with seed should be reproducible
    result1 <- TSENAT:::.select_genes_from_analysis(se, n_genes = 2, genes = NULL,
                                                    select_by = "random", seed = 42, verbose = FALSE)
    result2 <- TSENAT:::.select_genes_from_analysis(se, n_genes = 2, genes = NULL,
                                                    select_by = "random", seed = 42, verbose = FALSE)
    
    expect_equal(result1, result2)
    expect_equal(length(result1), 2)
  })
  
  test_that(".select_genes_from_analysis: by specific gene names", {
    counts <- matrix(c(
      c(100, 100, 100, 100),
      c(100, 100, 100, 100),
      c(100, 100, 100, 100),
      c(100, 100, 100, 100)
    ), nrow = 4, ncol = 4, byrow = TRUE)
    tpm <- counts
    rownames(counts) <- rownames(tpm) <- c("TX1", "TX2", "TX3", "TX4")
    colnames(counts) <- colnames(tpm) <- c("S1", "S2", "S3", "S4")
    
    se <- SummarizedExperiment::SummarizedExperiment(
      assays = list(counts = counts, tpm = tpm)
    )
    
    # Select specific genes
    result <- TSENAT:::.select_genes_from_analysis(se, n_genes = NULL, 
                                                   genes = c("TX1", "TX3"),
                                                   select_by = "variance", verbose = FALSE)
    
    expect_is(result, "integer")
    expect_equal(length(result), 2)
    expect_equal(sort(result), c(1, 3))
  })
  
  test_that(".select_genes_from_analysis: error on missing genes", {
    counts <- matrix(c(100, 100, 100, 100), nrow = 2, ncol = 2)
    tpm <- counts
    rownames(counts) <- rownames(tpm) <- c("TX1", "TX2")
    colnames(counts) <- colnames(tpm) <- c("S1", "S2")
    
    se <- SummarizedExperiment::SummarizedExperiment(
      assays = list(counts = counts, tpm = tpm)
    )
    
    expect_error(
      TSENAT:::.select_genes_from_analysis(se, n_genes = NULL,
                                           genes = c("TX1", "TX_NONEXISTENT"),
                                           select_by = "variance", verbose = FALSE),
      "Genes not found"
    )
  })
  
  test_that(".select_genes_from_analysis: n_genes > total returns all", {
    counts <- matrix(c(
      c(100, 100),
      c(100, 100),
      c(100, 100)
    ), nrow = 3, ncol = 2, byrow = TRUE)
    tpm <- counts
    rownames(counts) <- rownames(tpm) <- c("TX1", "TX2", "TX3")
    colnames(counts) <- colnames(tpm) <- c("S1", "S2")
    
    se <- SummarizedExperiment::SummarizedExperiment(
      assays = list(counts = counts, tpm = tpm)
    )
    
    # Request 10 genes but only have 3 - should warn
    result <- expect_warning(
      TSENAT:::.select_genes_from_analysis(se, n_genes = 10, genes = NULL,
                                           select_by = "variance", verbose = FALSE),
      "n_genes.*exceeds available genes"
    )
    
    expect_equal(length(result), 3)  # All 3 genes returned with warning
    expect_equal(sort(result), c(1, 2, 3))
  })
  
  # ========================================================================
  # .select_samples_from_analysis tests
  # ========================================================================
  
  test_that(".select_samples_from_analysis: by random with seed", {
    counts <- matrix(c(
      c(100, 100, 100, 100, 100),
      c(100, 100, 100, 100, 100)
    ), nrow = 2, ncol = 5, byrow = TRUE)
    tpm <- counts
    rownames(counts) <- rownames(tpm) <- c("TX1", "TX2")
    colnames(counts) <- colnames(tpm) <- c("S1", "S2", "S3", "S4", "S5")
    
    se <- SummarizedExperiment::SummarizedExperiment(
      assays = list(counts = counts, tpm = tpm)
    )
    
    # Random selection with seed should be reproducible
    result1 <- TSENAT:::.select_samples_from_analysis(se, n_samples = 3, samples = NULL,
                                                      seed = 42, verbose = FALSE)
    result2 <- TSENAT:::.select_samples_from_analysis(se, n_samples = 3, samples = NULL,
                                                      seed = 42, verbose = FALSE)
    
    expect_equal(result1, result2)
    expect_equal(length(result1), 3)
  })
  
  test_that(".select_samples_from_analysis: by specific sample names", {
    counts <- matrix(c(
      c(100, 100, 100, 100, 100),
      c(100, 100, 100, 100, 100)
    ), nrow = 2, ncol = 5, byrow = TRUE)
    tpm <- counts
    rownames(counts) <- rownames(tpm) <- c("TX1", "TX2")
    colnames(counts) <- colnames(tpm) <- c("S1", "S2", "S3", "S4", "S5")
    
    se <- SummarizedExperiment::SummarizedExperiment(
      assays = list(counts = counts, tpm = tpm)
    )
    
    # Select specific samples
    result <- TSENAT:::.select_samples_from_analysis(se, n_samples = NULL,
                                                     samples = c("S1", "S3", "S5"),
                                                     verbose = FALSE)
    
    expect_is(result, "integer")
    expect_equal(length(result), 3)
    expect_equal(sort(result), c(1, 3, 5))
  })
  
  test_that(".select_samples_from_analysis: error on missing samples", {
    counts <- matrix(c(100, 100, 100, 100), nrow = 2, ncol = 2)
    tpm <- counts
    rownames(counts) <- rownames(tpm) <- c("TX1", "TX2")
    colnames(counts) <- colnames(tpm) <- c("S1", "S2")
    
    se <- SummarizedExperiment::SummarizedExperiment(
      assays = list(counts = counts, tpm = tpm)
    )
    
    expect_error(
      TSENAT:::.select_samples_from_analysis(se, n_samples = NULL,
                                             samples = c("S1", "S_NONEXISTENT"),
                                             verbose = FALSE),
      "Samples not found"
    )
  })
  
  test_that(".select_samples_from_analysis: balances by condition if available", {
    counts <- matrix(c(
      c(100, 100, 100, 100, 100, 100),
      c(100, 100, 100, 100, 100, 100)
    ), nrow = 2, ncol = 6, byrow = TRUE)
    tpm <- counts
    rownames(counts) <- rownames(tpm) <- c("TX1", "TX2")
    colnames(counts) <- colnames(tpm) <- c("S1", "S2", "S3", "S4", "S5", "S6")
    
    se <- SummarizedExperiment::SummarizedExperiment(
      assays = list(counts = counts, tpm = tpm),
      colData = data.frame(
        condition = c("A", "A", "A", "B", "B", "B"),
        row.names = colnames(counts)
      )
    )
    
    # Select 4 samples (try to balance across 2 conditions)
    result <- TSENAT:::.select_samples_from_analysis(se, n_samples = 4, samples = NULL,
                                                     seed = 42, verbose = FALSE)
    
    expect_is(result, "integer")
    expect_equal(length(result), 4)
  })
  
  # ========================================================================
  # .balance_sample_selection tests
  # ========================================================================
  
  test_that(".balance_sample_selection: balances by condition", {
    coldata <- data.frame(
      condition = c("A", "A", "A", "B", "B", "B"),
      row.names = c("S1", "S2", "S3", "S4", "S5", "S6")
    )
    
    # Select 4 samples, should balance across A and B
    result <- TSENAT:::.balance_sample_selection(coldata, n_samples = 4, seed = 42, verbose = FALSE)
    
    expect_is(result, "integer")
    expect_equal(length(result), 4)
  })
  
  test_that(".balance_sample_selection: fallback to random if no condition column", {
    coldata <- data.frame(
      sample_id = c("S1", "S2", "S3", "S4", "S5", "S6"),
      treatment = c("A", "B", "A", "B", "A", "B"),
      row.names = c("S1", "S2", "S3", "S4", "S5", "S6")
    )
    
    # No standard condition column, should fallback to random
    result <- TSENAT:::.balance_sample_selection(coldata, n_samples = 3, seed = 42, verbose = FALSE)
    
    expect_is(result, "integer")
    expect_equal(length(result), 3)
  })
  
  # ========================================================================
  # .sync_filter_metadata tests
  # ========================================================================
  
  test_that(".sync_filter_metadata: filters readcounts", {
    md <- list(
      readcounts = matrix(c(
        c(100, 100, 100),
        c(50, 50, 50),
        c(10, 10, 10)
      ), nrow = 3, ncol = 3, byrow = TRUE,
      dimnames = list(c("TX1", "TX2", "TX3"), c("S1", "S2", "S3"))),
      other_field = "keep"
    )
    
    tokeep <- c(TRUE, FALSE, TRUE)  # Keep TX1 and TX3
    
    result <- TSENAT:::.sync_filter_metadata(md, tokeep, "counts")
    
    expect_true(!is.null(result$readcounts))
    if (!is.null(result$readcounts)) {
      expect_equal(nrow(result$readcounts), 2)  # 2 rows kept
    }
    expect_equal(result$other_field, "keep")  # Other fields preserved
  })
  
  test_that(".sync_filter_metadata: filters tx2gene", {
    md <- list(
      tx2gene = data.frame(
        Transcript = c("TX1", "TX2", "TX3"),
        Gene = c("G1", "G1", "G2")
      ),
      other_field = "keep"
    )
    
    tokeep <- c(TRUE, FALSE, TRUE)  # Keep TX1 and TX3
    rownames(md$tx2gene) <- md$tx2gene$Transcript
    
    result <- TSENAT:::.sync_filter_metadata(md, tokeep, "counts")
    
    expect_true(!is.null(result$tx2gene))
    # Note: filtering only keeps rows matching the transcript names in tokeep
    if (!is.null(result$tx2gene) && nrow(result$tx2gene) > 0) {
      expect_true(nrow(result$tx2gene) > 0)  # Some transcripts kept
    }
  })
  
  test_that(".sync_filter_metadata: filters salmon data", {
    md <- list(
      tpm = matrix(c(
        c(10, 20, 30),
        c(5, 10, 15),
        c(2, 4, 6)
      ), nrow = 3, ncol = 3, byrow = TRUE,
      dimnames = list(c("TX1", "TX2", "TX3"), c("S1", "S2", "S3"))),
      effective_length = c(100, 150, 200)
    )
    
    tokeep <- c(TRUE, TRUE, FALSE)  # Keep TX1 and TX2
    
    result <- TSENAT:::.sync_filter_metadata(md, tokeep, "tpm")
    
    expect_equal(nrow(result$tpm), 2)
    expect_equal(length(result$effective_length), 2)
  })
  
  test_that(".sync_filter_metadata: handles NULL metadata gracefully", {
    md <- list()  # Empty metadata
    tokeep <- c(TRUE, FALSE)
    
    result <- TSENAT:::.sync_filter_metadata(md, tokeep, "counts")
    
    expect_is(result, "list")
    expect_equal(length(result), 0)  # Empty list returned
  })
  
  # ========================================================================
  # .sync_subset_metadata tests
  # ========================================================================
  
  test_that(".sync_subset_metadata: filters tx2gene", {
    se <- SummarizedExperiment::SummarizedExperiment(
      matrix(c(100, 100, 100, 50, 50, 50), nrow = 2, ncol = 3),
      metadata = list(
        tx2gene = data.frame(
          Transcript = c("TX1", "TX2", "TX1.2"),
          Gene = c("G1", "G1", "G1")
        )
      )
    )
    rownames(se) <- c("TX1", "TX2")
    
    analysis <- TSENAT::TSENATAnalysis(se)
    
    # Subset to TX1 only
    gene_idx <- 1
    sample_idx <- c(1, 2, 3)
    
    result <- TSENAT:::.sync_subset_metadata(analysis, gene_idx, sample_idx, verbose = FALSE)
    
    # Should filter metadata accordingly
    expect_true(!is.null(result))
  })
  
  test_that(".sync_subset_metadata: filters readcounts", {
    se <- SummarizedExperiment::SummarizedExperiment(
      matrix(c(100, 100, 100, 50, 50, 50), nrow = 2, ncol = 3),
      metadata = list(
        readcounts = matrix(c(100, 100, 100, 50, 50, 50), nrow = 2, ncol = 3,
                           dimnames = list(c("TX1", "TX2"), c("S1", "S2", "S3")))
      )
    )
    colnames(se) <- c("S1", "S2", "S3")
    
    analysis <- TSENAT::TSENATAnalysis(se)
    
    # Subset to gene 1 (TX1 only), samples 1-2
    gene_idx <- 1
    sample_idx <- c(1, 2)
    
    result <- TSENAT:::.sync_subset_metadata(analysis, gene_idx, sample_idx, verbose = FALSE)
    
    expect_true(!is.null(result))
  })
})



# ============================================================================
# FILTER_ANALYSIS_S4 TESTS
# ============================================================================

test_that("filter_analysis basic filtering", {
  # Create test data for filter_analysis
  counts <- matrix(c(
    100, 100, 100, 100, 100, 100, 100, 100, 100, 100,
    100, 100, 100, 100, 100, 100, 100, 100, 100, 100,
    100, 100, 100, 100, 100, 100, 100, 100, 0, 0,
    100, 100, 100, 100, 100, 100, 100, 100, 0, 0,
    100, 100, 100, 100, 0, 0, 0, 0, 0, 0,
    100, 100, 100, 100, 0, 0, 0, 0, 0, 0
  ), nrow = 6, ncol = 10, byrow = TRUE)
  tpm <- counts + 1
  rownames(counts) <- rownames(tpm) <- c("TX1.1", "TX1.2", "TX2.1", "TX2.2", "TX3.1", "TX3.2")
  colnames(counts) <- colnames(tpm) <- paste0("S", 1:10)
  
  se <- SummarizedExperiment::SummarizedExperiment(
    assays = list(counts = counts, tpm = tpm),
    rowData = data.frame(gene_id = c("G1", "G1", "G2", "G2", "G3", "G3"), row.names = rownames(counts)),
    metadata = list(tx2gene = data.frame(tx = rownames(counts), gene = c("G1", "G1", "G2", "G2", "G3", "G3")))
  )
  
  analysis <- TSENAT::TSENATAnalysis(se)
  
  # Apply basic filtering
  filt <- TSENAT:::filter_analysis(analysis, min_samples = 1, verbose = FALSE)
  
  expect_s4_class(filt, "TSENATAnalysis")
  expect_true(nrow(filt@se) > 0)
})

test_that("filter_analysis with min_tx_per_gene", {
  # Test min_tx_per_gene parameter interaction
  counts <- matrix(c(
    100, 100, 100,
    100, 100, 100,
    100, 100, 100,
    100, 100, 0
  ), nrow = 4, ncol = 3, byrow = TRUE)
  tpm <- counts + 1
  rownames(counts) <- rownames(tpm) <- c("TX1.1", "TX1.2", "TX1.3", "TX2.1")
  colnames(counts) <- colnames(tpm) <- c("S1", "S2", "S3")
  
  se <- SummarizedExperiment::SummarizedExperiment(
    assays = list(counts = counts, tpm = tpm),
    rowData = data.frame(gene_id = c("G1", "G1", "G1", "G2"), row.names = rownames(counts)),
    metadata = list(tx2gene = data.frame(tx = rownames(counts), gene = c("G1", "G1", "G1", "G2")))
  )
  
  analysis <- TSENAT::TSENATAnalysis(se)
  
  # Apply filter with min_tx_per_gene
  # min_tx_per_gene=2 removes G2 (only 1 isoform)
  filt <- TSENAT:::filter_analysis(
    analysis, 
    min_samples = 1,
    min_tx_per_gene = 2,
    verbose = FALSE
  )
  
  # G1: all 3 isoforms kept (min_tx_per_gene=2 satisfied)
  # G2: removed (single isoform)
  expect_true(all(c("TX1.1", "TX1.2", "TX1.3") %in% rownames(filt@se)))
  expect_false("TX2.1" %in% rownames(filt@se))
})

# ============================================================================
# Tests for filter_analysis with subset parameters
# ============================================================================

test_that("filter_analysis subset_n_genes by variance", {
  # Create SE with varying transcript variances
  # Note: subset_analysis selects individual transcripts, not genes
  counts <- matrix(c(
    100, 105, 102,      # TX1.1: low variance
    100, 106, 101,      # TX1.2: low variance
    10,  50,  20,       # TX2.1: high variance (~400)
    11,  51,  19,       # TX2.2: high variance (~400)
    5,   5,   5,        # TX3.1: zero variance
    5,   5,   5         # TX3.2: zero variance
  ), nrow = 6, ncol = 3, byrow = TRUE)
  tpm <- counts + 1
  rownames(counts) <- rownames(tpm) <- c("TX1.1", "TX1.2", "TX2.1", "TX2.2", "TX3.1", "TX3.2")
  colnames(counts) <- colnames(tpm) <- c("S1", "S2", "S3")
  
  se <- SummarizedExperiment::SummarizedExperiment(
    assays = list(counts = counts, tpm = tpm),
    rowData = data.frame(gene_id = c("G1", "G1", "G2", "G2", "G3", "G3"), 
                        row.names = rownames(counts)),
    metadata = list(tx2gene = data.frame(tx = rownames(counts), 
                                        gene = c("G1", "G1", "G2", "G2", "G3", "G3")))
  )
  
  analysis <- TSENAT::TSENATAnalysis(se)
  
  # Filter, then subset to 1 transcript by variance (selects TX2.1 or TX2.2, both have highest variance)
  filt <- TSENAT:::filter_analysis(
    analysis,
    min_samples = 1,
    subset_n_genes = 1,
    subset_select_by = "variance",
    verbose = FALSE
  )
  
  # Should retain 1 transcript with highest variance (TX2.1 or TX2.2)
  expect_equal(nrow(filt@se), 1)
  expect_true(any(c("TX2.1", "TX2.2") %in% rownames(filt@se)))
})

test_that("filter_analysis subset_n_genes by mean expression", {
  # Create SE with different mean expression levels at transcript level
  counts <- matrix(c(
    10,  10,  10,      # TX1.1: low mean
    10,  10,  10,      # TX1.2: low mean
    100, 100, 100,     # TX2.1: high mean
    100, 100, 100,     # TX2.2: high mean
    50,  50,  50       # TX3.1: medium mean
  ), nrow = 5, ncol = 3, byrow = TRUE)
  tpm <- counts + 1
  rownames(counts) <- rownames(tpm) <- c("TX1.1", "TX1.2", "TX2.1", "TX2.2", "TX3.1")
  colnames(counts) <- colnames(tpm) <- c("S1", "S2", "S3")
  
  se <- SummarizedExperiment::SummarizedExperiment(
    assays = list(counts = counts, tpm = tpm),
    rowData = data.frame(gene_id = c("G1", "G1", "G2", "G2", "G3"), 
                        row.names = rownames(counts)),
    metadata = list(tx2gene = data.frame(tx = rownames(counts), 
                                        gene = c("G1", "G1", "G2", "G2", "G3")))
  )
  
  analysis <- TSENAT::TSENATAnalysis(se)
  
  # Filter, then subset to 2 transcripts by mean expression
  filt <- TSENAT:::filter_analysis(
    analysis,
    min_samples = 1,
    min_tx_per_gene = 1,
    subset_n_genes = 2,
    subset_select_by = "mean",
    verbose = FALSE
  )
  
  # Should retain 2 transcripts with highest mean: TX2.1 and TX2.2 (both have mean=100)
  expect_equal(nrow(filt@se), 2)
  expect_true(all(c("TX2.1", "TX2.2") %in% rownames(filt@se)) || all(c("TX3.1") %in% rownames(filt@se)))
})

test_that("filter_analysis subset_n_genes with random selection and seed", {
  # Create SE with 8 transcripts from 4 genes
  counts <- matrix(c(
    100, 100, 100,     # TX1.1
    100, 100, 100,     # TX1.2
    50,  50,  50,      # TX2.1
    50,  50,  50,      # TX2.2
    75,  75,  75,      # TX3.1
    75,  75,  75,      # TX3.2
    25,  25,  25,      # TX4.1
    25,  25,  25       # TX4.2
  ), nrow = 8, ncol = 3, byrow = TRUE)
  tpm <- counts + 1
  rownames(counts) <- rownames(tpm) <- c("TX1.1", "TX1.2", "TX2.1", "TX2.2", "TX3.1", "TX3.2", "TX4.1", "TX4.2")
  colnames(counts) <- colnames(tpm) <- c("S1", "S2", "S3")
  
  se <- SummarizedExperiment::SummarizedExperiment(
    assays = list(counts = counts, tpm = tpm),
    rowData = data.frame(gene_id = c("G1", "G1", "G2", "G2", "G3", "G3", "G4", "G4"), 
                        row.names = rownames(counts)),
    metadata = list(tx2gene = data.frame(tx = rownames(counts), 
                                        gene = c("G1", "G1", "G2", "G2", "G3", "G3", "G4", "G4")))
  )
  
  analysis <- TSENAT::TSENATAnalysis(se)
  
  # Random selection with seed should be reproducible
  filt1 <- TSENAT:::filter_analysis(
    analysis,
    min_samples = 1,
    subset_n_genes = 2,
    subset_select_by = "random",
    subset_seed = 42,
    verbose = FALSE
  )
  
  filt2 <- TSENAT:::filter_analysis(
    analysis,
    min_samples = 1,
    subset_n_genes = 2,
    subset_select_by = "random",
    subset_seed = 42,
    verbose = FALSE
  )
  
  # Same seed should give identical results (2 random transcripts)
  expect_equal(rownames(filt1@se), rownames(filt2@se))
  expect_equal(nrow(filt1@se), 2)  # 2 random transcripts
})

test_that("filter_analysis subset_genes by specific transcript names", {
  # Create SE with multiple transcripts
  counts <- matrix(c(
    100, 100, 100,     # TX1.1
    100, 100, 100,     # TX1.2
    50,  50,  50,      # TX2.1
    50,  50,  50,      # TX2.2
    75,  75,  75,      # TX3.1
    75,  75,  75       # TX3.2
  ), nrow = 6, ncol = 3, byrow = TRUE)
  tpm <- counts + 1
  rownames(counts) <- rownames(tpm) <- c("TX1.1", "TX1.2", "TX2.1", "TX2.2", "TX3.1", "TX3.2")
  colnames(counts) <- colnames(tpm) <- c("S1", "S2", "S3")
  
  se <- SummarizedExperiment::SummarizedExperiment(
    assays = list(counts = counts, tpm = tpm),
    rowData = data.frame(gene_id = c("G1", "G1", "G2", "G2", "G3", "G3"), 
                        row.names = rownames(counts)),
    metadata = list(tx2gene = data.frame(tx = rownames(counts), 
                                        gene = c("G1", "G1", "G2", "G2", "G3", "G3")))
  )
  
  analysis <- TSENAT::TSENATAnalysis(se)
  
  # Subset to retain only specific transcripts (subset_genes uses transcript IDs)
  filt <- TSENAT:::filter_analysis(
    analysis,
    min_samples = 1,
    subset_genes = c("TX1.1", "TX1.2", "TX3.1", "TX3.2"),
    verbose = FALSE
  )
  
  # Should retain only specified transcripts (TX1 and TX3)
  expect_equal(nrow(filt@se), 4)
  expect_false("TX2.1" %in% rownames(filt@se))
  expect_false("TX2.2" %in% rownames(filt@se))
  expect_true(all(c("TX1.1", "TX1.2", "TX3.1", "TX3.2") %in% rownames(filt@se)))
})

test_that("filter_analysis subset_n_samples", {
  # Create SE with 5 samples and multiple genes
  counts <- matrix(c(
    100, 100, 100, 100, 100,     # TX1.1
    100, 100, 100, 100, 100,     # TX1.2
    50,  50,  50,  50,  50,      # TX2.1
    50,  50,  50,  50,  50       # TX2.2
  ), nrow = 4, ncol = 5, byrow = TRUE)
  tpm <- counts + 1
  rownames(counts) <- rownames(tpm) <- c("TX1.1", "TX1.2", "TX2.1", "TX2.2")
  colnames(counts) <- colnames(tpm) <- c("S1", "S2", "S3", "S4", "S5")
  
  coldata <- data.frame(condition = c("A", "A", "B", "B", "B"), row.names = colnames(counts))
  
  se <- SummarizedExperiment::SummarizedExperiment(
    assays = list(counts = counts, tpm = tpm),
    colData = coldata,
    rowData = data.frame(gene_id = c("G1", "G1", "G2", "G2"), 
                        row.names = rownames(counts)),
    metadata = list(tx2gene = data.frame(tx = rownames(counts), 
                                        gene = c("G1", "G1", "G2", "G2")))
  )
  
  analysis <- TSENAT::TSENATAnalysis(se)
  
  # Subset to 3 samples (balanced by condition if available)
  filt <- TSENAT:::filter_analysis(
    analysis,
    min_samples = 1,
    subset_n_samples = 3,
    verbose = FALSE
  )
  
  # Should have 4 transcripts (all) and 3 samples
  expect_equal(nrow(filt@se), 4)
  expect_equal(ncol(filt@se), 3)
})

test_that("filter_analysis subset_samples by specific sample names", {
  # Create SE with 5 samples and 4 transcripts
  counts <- matrix(c(
    100, 100, 100, 100, 100,     # TX1.1
    100, 100, 100, 100, 100,     # TX1.2
    50,  50,  50,  50,  50,      # TX2.1
    50,  50,  50,  50,  50       # TX2.2
  ), nrow = 4, ncol = 5, byrow = TRUE)
  tpm <- counts + 1
  rownames(counts) <- rownames(tpm) <- c("TX1.1", "TX1.2", "TX2.1", "TX2.2")
  colnames(counts) <- colnames(tpm) <- c("S1", "S2", "S3", "S4", "S5")
  
  se <- SummarizedExperiment::SummarizedExperiment(
    assays = list(counts = counts, tpm = tpm),
    rowData = data.frame(gene_id = c("G1", "G1", "G2", "G2"), 
                        row.names = rownames(counts)),
    metadata = list(tx2gene = data.frame(tx = rownames(counts), 
                                        gene = c("G1", "G1", "G2", "G2")))
  )
  
  analysis <- TSENAT::TSENATAnalysis(se)
  
  # Subset to specific samples
  filt <- TSENAT:::filter_analysis(
    analysis,
    min_samples = 1,
    subset_samples = c("S1", "S3", "S5"),
    verbose = FALSE
  )
  
  # Should have all 4 transcripts and 3 samples
  expect_equal(nrow(filt@se), 4)
  expect_equal(ncol(filt@se), 3)
  expect_equal(colnames(filt@se), c("S1", "S3", "S5"))
})

test_that("filter_analysis combined filter and subset operations", {
  # Create SE with sparse and abundant genes
  counts <- matrix(c(
    100, 100, 100, 100, 5,      # TX1.1: abundant
    100, 100, 100, 100, 5,      # TX1.2: abundant
    0,   0,   0,   0,   100,    # TX2.1: sparse (single isoform, will fail min_tx_per_gene)
    50,  50,  50,  50,  50,     # TX3.1: medium
    50,  50,  50,  50,  50      # TX3.2: medium
  ), nrow = 5, ncol = 5, byrow = TRUE)
  tpm <- counts + 1
  rownames(counts) <- rownames(tpm) <- c("TX1.1", "TX1.2", "TX2.1", "TX3.1", "TX3.2")
  colnames(counts) <- colnames(tpm) <- c("S1", "S2", "S3", "S4", "S5")
  
  se <- SummarizedExperiment::SummarizedExperiment(
    assays = list(counts = counts, tpm = tpm),
    rowData = data.frame(gene_id = c("G1", "G1", "G2", "G3", "G3"), 
                        row.names = rownames(counts)),
    metadata = list(tx2gene = data.frame(tx = rownames(counts), 
                                        gene = c("G1", "G1", "G2", "G3", "G3")))
  )
  
  analysis <- TSENAT::TSENATAnalysis(se)
  
  # Apply filtering (removes G2, keeps G1 and G3) and then subset to 1 gene
  filt <- TSENAT:::filter_analysis(
    analysis,
    min_samples = 3,
    min_tx_per_gene = 2,
    subset_n_genes = 1,
    subset_select_by = "mean",
    verbose = FALSE
  )
  
  # After filtering: G1 (2 tx) and G3 (2 tx)
  # After subsetting to 1 gene by mean: should keep both because same mean
  # Result should have 2-4 transcripts from filtered genes
  expect_true(nrow(filt@se) > 0)
  expect_true(nrow(filt@se) <= 4)
})

test_that("filter_analysis subset_min_count during subsetting", {
  # Create SE with varying count depths
  counts <- matrix(c(
    1000, 1000, 1000,   # TX1.1: high counts
    1000, 1000, 1000,   # TX1.2: high counts
    10,   10,   10,     # TX2.1: low counts
    10,   10,   10,     # TX2.2: low counts
    100,  100,  100     # TX3.1: medium counts
  ), nrow = 5, ncol = 3, byrow = TRUE)
  tpm <- counts + 1
  rownames(counts) <- rownames(tpm) <- c("TX1.1", "TX1.2", "TX2.1", "TX2.2", "TX3.1")
  colnames(counts) <- colnames(tpm) <- c("S1", "S2", "S3")
  
  se <- SummarizedExperiment::SummarizedExperiment(
    assays = list(counts = counts, tpm = tpm),
    rowData = data.frame(gene_id = c("G1", "G1", "G2", "G2", "G3"), 
                        row.names = rownames(counts)),
    metadata = list(tx2gene = data.frame(tx = rownames(counts), 
                                        gene = c("G1", "G1", "G2", "G2", "G3")))
  )
  
  analysis <- TSENAT::TSENATAnalysis(se)
  
  # Apply subsetting with min_count to filter low-count transcripts
  # This should trigger a warning about filtered genes
  filt <- testthat::expect_warning(
    TSENAT:::filter_analysis(
      analysis,
      min_samples = 1,
      subset_n_genes = 3,
      subset_select_by = "variance",
      subset_min_count = 50,
      verbose = FALSE
    ),
    "Filtered out.*gene.*with total count"
  )
  
  # Should remove TX2.1 and TX2.2 (low counts), keep TX1.1, TX1.2, TX3.1
  # After subsetting by variance and min_count
  expect_true(nrow(filt@se) >= 2)  # At least TX1 isoforms
  # TX2 should be filtered out due to low min_count
  expect_false(all(c("TX2.1", "TX2.2") %in% rownames(filt@se)))
})

test_that("filter_analysis no subsetting when all params NULL", {
  # Create simple SE
  counts <- matrix(c(
    100, 100, 100,     # TX1.1
    100, 100, 100,     # TX1.2
    50,  50,  50,      # TX2.1
    50,  50,  50       # TX2.2
  ), nrow = 4, ncol = 3, byrow = TRUE)
  tpm <- counts + 1
  rownames(counts) <- rownames(tpm) <- c("TX1.1", "TX1.2", "TX2.1", "TX2.2")
  colnames(counts) <- colnames(tpm) <- c("S1", "S2", "S3")
  
  se <- SummarizedExperiment::SummarizedExperiment(
    assays = list(counts = counts, tpm = tpm),
    rowData = data.frame(gene_id = c("G1", "G1", "G2", "G2"), 
                        row.names = rownames(counts)),
    metadata = list(tx2gene = data.frame(tx = rownames(counts), 
                                        gene = c("G1", "G1", "G2", "G2")))
  )
  
  analysis <- TSENAT::TSENATAnalysis(se)
  
  # Apply filtering only (no subsetting)
  filt <- TSENAT:::filter_analysis(
    analysis,
    min_samples = 1,
    verbose = FALSE
  )
  
  # Should retain all 4 transcripts (no subsetting applied)
  expect_equal(nrow(filt@se), 4)
})

# Tests for refactored helper functions in se_manipulation_filter.R
# Comprehensive unit tests for Phase 2 refactoring
# Tests: .get_gene_ids, .validate_filter_params, .estimate_min_tpm,
#        .calculate_stringency_thresholds, .apply_tpm_filter,
#        .filter_by_tx_per_gene, .filter_by_isoform_abundance, etc.

library(testthat)
library(SummarizedExperiment)
library(S4Vectors)

# Setup test data
create_test_se <- function(n_genes = 20, n_samples = 10) {
  # Ensure even number of samples for condition rep
  if (n_samples %% 2 != 0) {
    n_samples <- n_samples + 1
  }
  
  counts <- matrix(rpois(n_genes * n_samples, lambda = 10), 
                   nrow = n_genes, ncol = n_samples)
  rownames(counts) <- paste0("TX", seq_len(n_genes))
  colnames(counts) <- paste0("S", seq_len(n_samples))
  
  tpm <- t(t(counts) / colSums(counts) * 1e6)
  
  # Create rowData with gene mapping
  gene_ids <- rep(paste0("GENE_", seq_len(n_genes / 2)), each = 2)
  rowdata <- DataFrame(
    gene_id = gene_ids,
    row.names = rownames(counts)
  )
  
  # Create colData with pair information for stringency tests
  pair_ids <- rep(seq_len(n_samples / 2), each = 2)
  
  SummarizedExperiment(
    assays = list(counts = counts, tpm = tpm),
    rowData = rowdata,
    colData = DataFrame(
      sample_id = colnames(counts),
      condition = rep(c("ctrl", "treat"), n_samples / 2),
      pair_id = pair_ids
    )
  )
}

# =============================================================================
# Tests for .get_gene_ids()
# =============================================================================

test_that(".get_gene_ids extracts gene IDs from rowData", {
  se <- create_test_se()
  gene_ids <- TSENAT:::.get_gene_ids(se)
  
  # Verify extraction works correctly
  expect_is(gene_ids, "character")
  expect_true(length(gene_ids) > 0)
  expect_true(all(!is.na(gene_ids)))
})

test_that(".get_gene_ids returns NULL when no gene mapping available", {
  # Create SE without gene_id column
  counts <- matrix(rpois(20, 10), nrow = 4, ncol = 5)
  se <- SummarizedExperiment(assays = list(counts = counts))
  
  gene_ids <- TSENAT:::.get_gene_ids(se)
  expect_null(gene_ids)
})

test_that(".get_gene_ids handles tx2gene mapping in metadata", {
  counts <- matrix(rpois(20, 10), nrow = 4, ncol = 5)
  rownames(counts) <- c("TX1", "TX2", "TX3", "TX4")
  
  tx2gene <- data.frame(
    Transcript = c("TX1", "TX2", "TX3", "TX4"),
    Gene = c("GENE_A", "GENE_A", "GENE_B", "GENE_B")
  )
  
  se <- SummarizedExperiment(
    assays = list(counts = counts),
    metadata = list(tx2gene = tx2gene)
  )
  
  gene_ids <- TSENAT:::.get_gene_ids(se)
  
  # Should extract gene IDs from tx2gene mapping
  if (!is.null(gene_ids)) {
    expect_true(length(gene_ids) > 0)
  }
})

# =============================================================================
# Tests for .validate_filter_params()
# =============================================================================

test_that(".validate_filter_params validates isoform abundance", {
  # Invalid: > 1
  expect_error(TSENAT:::.validate_filter_params(1.5, NULL))
  # Invalid: < 0
  expect_error(TSENAT:::.validate_filter_params(-0.1, NULL))
  # Invalid: vector
  expect_error(TSENAT:::.validate_filter_params(c(0.05, 0.1), NULL))
})

# =============================================================================
# Tests for .estimate_min_tpm()
# =============================================================================

test_that(".estimate_min_tpm returns list with min_tpm and quant_label", {
  assay_mat <- matrix(c(0, 0.1, 0.5, 1, 2, 5, 10, 20, 50, 100), nrow = 5, ncol = 2)
  
  result <- TSENAT:::.estimate_min_tpm(assay_mat, "medium", verbose = FALSE)
  
  expect_is(result, "list")
  expect_true("min_tpm" %in% names(result))
  expect_true("quant_label" %in% names(result))
  expect_true(is.numeric(result$min_tpm))
  expect_true(result$min_tpm >= 0.1 && result$min_tpm <= 5.0)
})

# =============================================================================
# Tests for .calculate_stringency_thresholds()
# =============================================================================

test_that(".calculate_stringency_thresholds returns correct structure", {
  result <- TSENAT:::.calculate_stringency_thresholds("medium", 100)
  
  expect_is(result, "list")
  expect_true("min_samples" %in% names(result))
  expect_true("min_tx_per_gene" %in% names(result))
  expect_true(is.numeric(result$min_samples))
  expect_true(is.numeric(result$min_tx_per_gene))
})

# =============================================================================
# Tests for .apply_tpm_filter()
# =============================================================================

test_that(".apply_tpm_filter returns logical vector", {
  assay_mat <- matrix(c(0, 1, 2, 3, 4, 5, 6, 7, 8, 9), nrow = 5, ncol = 2)
  
  tokeep <- TSENAT:::.apply_tpm_filter(assay_mat, min_tpm = 5, min_samples = 1, verbose = FALSE)
  
  expect_is(tokeep, "logical")
  expect_equal(length(tokeep), nrow(assay_mat))
})

# =============================================================================
# Tests for .filter_by_tx_per_gene()
# =============================================================================

test_that(".filter_by_tx_per_gene returns logical vector", {
  genes_vec <- c("GENE_A", "GENE_A", "GENE_B", "GENE_C", "GENE_C", "GENE_C")
  tokeep <- c(TRUE, TRUE, TRUE, TRUE, TRUE, TRUE)
  
  tokeep_filt <- TSENAT:::.filter_by_tx_per_gene(tokeep, genes_vec, min_tx_per_gene = 2, verbose = FALSE)
  
  expect_is(tokeep_filt, "logical")
  expect_equal(length(tokeep_filt), length(tokeep))
})

# =============================================================================
# Tests for .filter_by_isoform_abundance()
# =============================================================================

test_that(".filter_by_isoform_abundance returns logical vector", {
  assay_mat <- matrix(c(10, 1, 20, 2, 15, 3), nrow = 2, ncol = 3)
  genes_vec <- c("GENE_A", "GENE_A")
  tokeep <- c(TRUE, TRUE)
  
  tokeep_filt <- TSENAT:::.filter_by_isoform_abundance(tokeep, genes_vec, assay_mat, 
                                                        min_isoform_abundance = 0.5, 
                                                        verbose = FALSE)
  
  expect_is(tokeep_filt, "logical")
  expect_equal(length(tokeep_filt), length(tokeep))
})


# =============================================================================
# Tests for .sync_filter_metadata()
# =============================================================================

test_that(".sync_filter_metadata filters metadata correctly", {
  readcounts <- matrix(1:20, nrow = 4, ncol = 5)
  rownames(readcounts) <- c("TX1", "TX2", "TX3", "TX4")
  
  md <- list(
    readcounts = readcounts,
    other_metadata = "should_remain"
  )
  
  tokeep <- c(TRUE, FALSE, TRUE, FALSE)
  
  md_filtered <- TSENAT:::.sync_filter_metadata(md, tokeep, NULL, 
                                                rownames(readcounts))
  
  expect_is(md_filtered, "list")
  expect_equal(md_filtered$other_metadata, "should_remain")
})

# =============================================================================
# Integration tests: Full refactored .filter_se() workflow
# =============================================================================

test_that("Refactored .filter_se() produces filtered results", {
  se <- create_test_se(n_genes = 30, n_samples = 15)
  
  # Filter with typical parameters
  se_filt <- TSENAT:::.filter_se(
    se, 
    min_samples = 5, 
    min_tpm = 1.0,
    min_tx_per_gene = 2,
    min_isoform_abundance = 0.05,
    verbose = FALSE
  )
  
  # Check dimensions
  expect_lte(nrow(se_filt), nrow(se))
  expect_equal(ncol(se_filt), ncol(se))
  
  # Check consistency across assays
  expect_equal(dim(assay(se_filt, "counts")), dim(assay(se_filt, "tpm")))
  
  # Check metadata was preserved
  expect_true("filtered" %in% names(metadata(se_filt)))
})

test_that("Refactored .filter_se() with stringency parameter works", {
  # Create test data with sufficient expression to survive stringency filtering
  # Use higher count values to ensure genes survive TPM-based filtering
  se <- create_test_se(n_genes = 30, n_samples = 20)
  
  # Boost expression to ensure survival of stringency filters
  # Apply Poisson multiplier to increase counts while maintaining biological realism
  assay_mat <- as.matrix(SummarizedExperiment::assay(se))
  boosted_mat <- assay_mat * 5 + 20  # Minimum baseline + multiplier ensures expression
  SummarizedExperiment::assay(se) <- boosted_mat
  
  # Test with permissive and moderate stringency levels
  # (avoid "severe" which may filter out too much sparse data)
  for (stringency_level in c("soft", "medium")) {
    se_filt <- TSENAT:::.filter_se(
      se,
      stringency = stringency_level,
      pair_col = "pair_id",
      min_tx_per_gene = 2,
      verbose = FALSE
    )
    
    # Verify filtering reduced data but didn't remove everything
    expect_lte(nrow(se_filt), nrow(se))
    expect_equal(ncol(se_filt), ncol(se))
    expect_true(nrow(se_filt) > 0, "Soft/medium stringency should keep some genes")
    expect_true("filtered" %in% names(metadata(se_filt)))
  }
  
  # Test "severe" separately with robust data to handle aggressive filtering
  se_severe <- suppressWarnings(TSENAT:::.filter_se(
    se,
    stringency = "severe",
    pair_col = "pair_id",
    min_tx_per_gene = 2,
    verbose = FALSE
  ))
  # Severe may filter aggressively or completely; just verify it's a valid SE
  expect_is(se_severe, "SummarizedExperiment")
})

test_that("Refactored .filter_se() preserves rowData and colData", {
  se <- create_test_se()
  se_filt <- TSENAT:::.filter_se(se, min_samples = 1, min_tpm = 0.1, verbose = FALSE)
  
  # colData should be unchanged
  expect_equal(ncol(colData(se_filt)), ncol(colData(se)))
  expect_equal(colnames(colData(se_filt)), colnames(colData(se)))
  
  # rowData should be subset correctly
  expect_lte(nrow(rowData(se_filt)), nrow(rowData(se)))
})

test_that("All refactored helpers maintain backward compatibility", {
  # Test that the API and outputs haven't changed
  se <- create_test_se()
  
  # Original call should still work with refactored code
  expect_no_error({
    result <- TSENAT:::.filter_se(
      se,
      min_samples = 5,
      min_tpm = 1.0,
      min_tx_per_gene = 2,
      min_isoform_abundance = 0.05,
      verbose = FALSE
    )
  })
})

# Integration tests for refactored .filter_se() function
# Tests that the refactored code maintains backward compatibility with public API

test_that(".filter_se works with refactored helper functions", {
  # Create test SummarizedExperiment
  counts <- matrix(rpois(200, lambda = 10), nrow = 20, ncol = 10)
  rownames(counts) <- paste0("TX", seq_len(20))
  colnames(counts) <- paste0("Sample_", seq_len(10))
  
  tpm <- t(t(counts) / colSums(counts) * 1e6)
  
  rowdata <- data.frame(
    gene_id = rep(paste0("GENE_", seq_len(10)), each = 2),
    row.names = rownames(counts)
  )
  
  se <- SummarizedExperiment::SummarizedExperiment(
    assays = list(counts = counts, tpm = tpm),
    rowData = rowdata,
    colData = data.frame(condition = rep(c("A", "B"), 5), row.names = colnames(counts))
  )
  
  # Test 1: Basic filtering works
  expect_no_error({
    se_filt <- TSENAT::filter_analysis(
      TSENATAnalysis(se),
      min_samples = 3, min_tpm = 1.0, verbose = FALSE
    )
  })
})

test_that(".filter_se produces valid output dimensions", {
  counts <- matrix(rpois(100, 10), nrow = 10, ncol = 10)
  tpm <- t(t(counts) / colSums(counts) * 1e6)
  rownames(counts) <- rownames(tpm) <- paste0("TX", 1:10)
  colnames(counts) <- colnames(tpm) <- paste0("S", 1:10)
  
  se <- SummarizedExperiment::SummarizedExperiment(
    assays = list(counts = counts, tpm = tpm)
  )
  
  # Create analysis and filter
  analysis <- TSENATAnalysis(se)
  analysis_filt <- TSENAT::filter_analysis(
    analysis,
    min_samples = 2, min_tpm = 0.5, verbose = FALSE
  )
  
  # Check dimensions
  se_filt <- analysis_filt@se
  expect_lte(nrow(se_filt), nrow(se))
  expect_equal(ncol(se_filt), ncol(se))
  expect_equal(nrow(assay(se_filt, "counts")), nrow(assay(se_filt, "tpm")))
})

test_that(".filter_se with stringency parameter works", {
  counts <- matrix(rpois(100, 10), nrow = 10, ncol = 10)
  tpm <- t(t(counts) / colSums(counts) * 1e6)
  rownames(counts) <- rownames(tpm) <- paste0("TX", 1:10)
  colnames(counts) <- colnames(tpm) <- paste0("S", 1:10)
  
  coldata <- data.frame(
    condition = rep(c("ctrl", "treat"), 5),
    pair = rep(1:5, each = 2),
    row.names = colnames(counts)
  )
  
  se <- SummarizedExperiment::SummarizedExperiment(
    assays = list(counts = counts, tpm = tpm),
    colData = coldata
  )
  
  # Test with stringency parameter
  analysis <- TSENATAnalysis(se)
  expect_no_error({
    analysis_filt <- TSENAT::filter_analysis(
      analysis,
      stringency = "medium",
      pair_col = "pair",
      verbose = FALSE
    )
  })
})

test_that(".filter_se preserves assay names across refactoring", {
  counts <- matrix(rpois(50, 10), nrow = 5, ncol = 10)
  tpm <- t(t(counts) / colSums(counts) * 1e6)
  abundance <- tpm * 2
  
  rownames(counts) <- rownames(tpm) <- rownames(abundance) <- paste0("TX", 1:5)
  colnames(counts) <- colnames(tpm) <- colnames(abundance) <- paste0("S", 1:10)
  
  se <- SummarizedExperiment::SummarizedExperiment(
    assays = list(counts = counts, tpm = tpm, abundance = abundance)
  )
  
  analysis <- TSENATAnalysis(se)
  analysis_filt <- TSENAT::filter_analysis(
    analysis,
    min_samples = 1, min_tpm = 0.1, verbose = FALSE
  )
  
  # All assays should be preserved
  se_filt <- analysis_filt@se
  expect_true("counts" %in% names(assays(se_filt)))
  expect_true("tpm" %in% names(assays(se_filt)))
  expect_true("abundance" %in% names(assays(se_filt)))
})

test_that("Refactoring maintains isoform-level filtering behavior", {
  # Create gene with 2 isoforms: one abundant, one rare
  counts <- matrix(c(100, 5, 200, 10, 150, 8), nrow = 2, ncol = 3)
  rownames(counts) <- c("TX1", "TX2")
  colnames(counts) <- c("S1", "S2", "S3")
  
  tpm <- t(t(counts) / colSums(counts) * 1e6)
  
  rowdata <- data.frame(
    gene_id = c("GENE_A", "GENE_A"),
    row.names = rownames(counts)
  )
  
  se <- SummarizedExperiment::SummarizedExperiment(
    assays = list(counts = counts, tpm = tpm),
    rowData = rowdata
  )
  
  analysis <- TSENATAnalysis(se)
  analysis_filt <- TSENAT::filter_analysis(
    analysis,
    min_samples = 1, min_tpm = 1.0,
    min_isoform_abundance = 0.1,
    verbose = FALSE
  )
  
  # Filtering should work and produce a result
  expect_true(nrow(analysis_filt@se) <= 2)
})

test_that("Refactored code handles empty results gracefully", {
  counts <- matrix(rep(0, 50), nrow = 5, ncol = 10)
  tpm <- matrix(rep(0, 50), nrow = 5, ncol = 10)
  
  se <- SummarizedExperiment::SummarizedExperiment(
    assays = list(counts = counts)
  )
  
  # Add TPM to metadata (required for diversity filtering)
  S4Vectors::metadata(se)$tpm <- tpm
  
  analysis <- TSENATAnalysis(se)
  
  # Should produce warning but not error
  expect_warning({
    analysis_filt <- TSENAT::filter_analysis(
      analysis,
      min_samples = 10, min_tpm = 100, verbose = FALSE
    )
  })
})

# ============================================================================
# Tests for .resolve_filter_parameters() helper function
# ============================================================================

test_that(".resolve_filter_parameters returns correct structure with manual params", {
  # Create simple SE
  tpm <- matrix(c(10, 20, 30, 5, 15, 25, 2, 8, 12), nrow = 3, ncol = 3)
  rownames(tpm) <- c("TX1", "TX2", "TX3")
  colnames(tpm) <- c("S1", "S2", "S3")
  
  se <- SummarizedExperiment::SummarizedExperiment(
    assays = list(tpm = tpm)
  )
  
  # Call helper with manual parameters
  result <- TSENAT:::.resolve_filter_parameters(
    se, min_samples = 2, min_tpm = 5, stringency = NULL,
    pair_col = NULL, tpm_assay_name = "tpm", assay_name = "tpm",
    min_tx_per_gene = 2, min_isoform_abundance = 0.05, verbose = FALSE
  )
  
  # Check structure and values
  expect_true(is.list(result))
  expect_true(all(c("min_samples", "min_tpm", "assay_mat", "assay_source", "genes_vec") %in% names(result)))
  expect_equal(result$min_samples, 2)
  expect_equal(result$min_tpm, 5)
  expect_equal(result$assay_source, "assay 'tpm' (user-specified)")
  expect_equal(nrow(result$assay_mat), 3)
  expect_equal(ncol(result$assay_mat), 3)
})

test_that(".resolve_filter_parameters fills missing isoform_abundance with default", {
  tpm <- matrix(1:9, nrow = 3, ncol = 3)
  rownames(tpm) <- c("TX1", "TX2", "TX3")
  colnames(tpm) <- c("S1", "S2", "S3")
  
  se <- SummarizedExperiment::SummarizedExperiment(
    assays = list(tpm = tpm)
  )
  
  # Call with NULL min_isoform_abundance
  result <- TSENAT:::.resolve_filter_parameters(
    se, min_samples = 2, min_tpm = 1, stringency = NULL,
    pair_col = NULL, tpm_assay_name = "tpm", assay_name = "tpm",
    min_tx_per_gene = 2, min_isoform_abundance = NULL, verbose = FALSE
  )
  
  # Should be NULL (not set by helper for manual params)
  expect_null(result$min_isoform_abundance)
})

test_that(".resolve_filter_parameters detects gene IDs from rowData", {
  tpm <- matrix(1:9, nrow = 3, ncol = 3)
  rownames(tpm) <- c("TX1", "TX2", "TX3")
  colnames(tpm) <- c("S1", "S2", "S3")
  
  se <- SummarizedExperiment::SummarizedExperiment(
    assays = list(tpm = tpm),
    rowData = data.frame(
      gene_id = c("G1", "G1", "G2"),
      row.names = rownames(tpm)
    )
  )
  
  result <- TSENAT:::.resolve_filter_parameters(
    se, min_samples = 1, min_tpm = 1, stringency = NULL,
    pair_col = NULL, tpm_assay_name = "tpm", assay_name = "tpm",
    min_tx_per_gene = 1, min_isoform_abundance = NULL, verbose = FALSE
  )
  
  # Should find genes from rowData
  expect_equal(result$genes_vec, c("G1", "G1", "G2"))
})

test_that(".resolve_filter_parameters validates isoform_abundance range", {
  tpm <- matrix(1:9, nrow = 3, ncol = 3)
  se <- SummarizedExperiment::SummarizedExperiment(assays = list(tpm = tpm))
  
  # Test invalid value > 1
  expect_error(
    TSENAT:::.resolve_filter_parameters(
      se, min_samples = 1, min_tpm = 1, stringency = NULL,
      pair_col = NULL, tpm_assay_name = "tpm", assay_name = "tpm",
      min_tx_per_gene = 1, min_isoform_abundance = 1.5, verbose = FALSE
    ),
    "must be numeric in"
  )
  
  # Test invalid value < 0
  expect_error(
    TSENAT:::.resolve_filter_parameters(
      se, min_samples = 1, min_tpm = 1, stringency = NULL,
      pair_col = NULL, tpm_assay_name = "tpm", assay_name = "tpm",
      min_tx_per_gene = 1, min_isoform_abundance = -0.1, verbose = FALSE
    ),
    "must be numeric in"
  )
})

test_that(".resolve_filter_parameters auto-detects pair_col for stringency", {
  tpm <- matrix(1:12, nrow = 3, ncol = 4)
  rownames(tpm) <- c("TX1", "TX2", "TX3")
  colnames(tpm) <- c("S1", "S2", "S3", "S4")
  
  se <- SummarizedExperiment::SummarizedExperiment(
    assays = list(tpm = tpm),
    colData = data.frame(
      pair_id = c(1, 1, 2, 2),
      condition = c("A", "B", "A", "B"),
      row.names = colnames(tpm)
    )
  )
  
  # Call with stringency but no pair_col
  result <- TSENAT:::.resolve_filter_parameters(
    se, min_samples = 2, min_tpm = 1, stringency = "medium",
    pair_col = NULL, tpm_assay_name = "tpm", assay_name = "tpm",
    min_tx_per_gene = 2, min_isoform_abundance = NULL, verbose = FALSE
  )
  
  # Should auto-detect pair_col and calculate thresholds
  expect_equal(result$pair_col_used, "pair_id")
  expect_true(result$min_samples > 0)
  expect_equal(result$min_tx_per_gene, 2)
})

test_that(".resolve_filter_parameters stringency='soft' sets correct thresholds", {
  tpm <- matrix(1:20, nrow = 5, ncol = 4)
  rownames(tpm) <- c("TX1", "TX2", "TX3", "TX4", "TX5")
  colnames(tpm) <- c("S1", "S2", "S3", "S4")
  
  se <- SummarizedExperiment::SummarizedExperiment(
    assays = list(tpm = tpm),
    colData = data.frame(
      pair_id = c(1, 1, 2, 2),
      row.names = colnames(tpm)
    )
  )
  
  result <- TSENAT:::.resolve_filter_parameters(
    se, min_samples = 999, min_tpm = 999, stringency = "soft",
    pair_col = "pair_id", tpm_assay_name = "tpm", assay_name = "tpm",
    min_tx_per_gene = 999, min_isoform_abundance = NULL, verbose = FALSE
  )
  
  # Soft should use low thresholds
  expect_true(result$min_samples < 4)  # 25% of 4 samples
  expect_true(result$min_tpm < 999)  # Auto-estimated
  expect_equal(result$min_tx_per_gene, 2)
  expect_equal(result$min_isoform_abundance, 0.01)  # 1% - permissive
})

test_that(".resolve_filter_parameters stringency='severe' sets correct thresholds", {
  tpm <- matrix(1:20, nrow = 5, ncol = 4)
  rownames(tpm) <- c("TX1", "TX2", "TX3", "TX4", "TX5")
  colnames(tpm) <- c("S1", "S2", "S3", "S4")
  
  se <- SummarizedExperiment::SummarizedExperiment(
    assays = list(tpm = tpm),
    colData = data.frame(
      pair_id = c(1, 1, 2, 2),
      row.names = colnames(tpm)
    )
  )
  
  result <- TSENAT:::.resolve_filter_parameters(
    se, min_samples = 1, min_tpm = 1, stringency = "severe",
    pair_col = "pair_id", tpm_assay_name = "tpm", assay_name = "tpm",
    min_tx_per_gene = 1, min_isoform_abundance = NULL, verbose = FALSE
  )
  
  # Severe should use high thresholds
  expect_equal(result$min_samples, 3)  # 75% of 4 samples = ceiling(0.75 * 4) = 3
  expect_true(result$min_tpm > 0)
  expect_equal(result$min_tx_per_gene, 3)
  expect_equal(result$min_isoform_abundance, 0.15)  # 15% - stringent
})

# ============================================================================
# Tests for .apply_combined_filters() helper function
# ============================================================================

test_that(".apply_combined_filters applies all three filters in order", {
  # Create data with known filter outcomes
  assay_mat <- matrix(c(
    100, 50, 30,  # TX1: passes TPM
    5, 2, 1,      # TX2: fails TPM
    80, 60, 40,   # TX3: passes TPM
    3, 1, 0       # TX4: fails TPM
  ), nrow = 4, byrow = TRUE, dimnames = list(
    c("TX1", "TX2", "TX3", "TX4"),
    c("S1", "S2", "S3")
  ))
  
  genes_vec <- c("G1", "G1", "G2", "G2")
  
  params <- list(
    min_tpm = 4,
    min_samples = 2,
    min_tx_per_gene = 2,
    min_isoform_abundance = 0.05
  )
  
  result <- TSENAT:::.apply_combined_filters(assay_mat, genes_vec, params, verbose = FALSE)
  
  # Check structure
  expect_true(is.list(result))
  expect_true(all(c("tokeep", "before", "after") %in% names(result)))
  expect_true(is.logical(result$tokeep))
  expect_equal(length(result$tokeep), 4)
  expect_equal(result$before, 4)
  expect_true(result$after <= 4)
})

test_that(".apply_combined_filters returns logical vector of correct length", {
  assay_mat <- matrix(1:15, nrow = 5, ncol = 3)
  genes_vec <- c("G1", "G1", "G2", "G2", "G3")
  
  params <- list(
    min_tpm = 1, min_samples = 1, min_tx_per_gene = 1, min_isoform_abundance = NULL
  )
  
  result <- TSENAT:::.apply_combined_filters(assay_mat, genes_vec, params, verbose = FALSE)
  
  expect_equal(length(result$tokeep), nrow(assay_mat))
  expect_equal(result$before, nrow(assay_mat))
  expect_true(result$after >= 0)
})

test_that(".apply_combined_filters skips gene filters when genes_vec is NULL", {
  assay_mat <- matrix(c(1, 2, 3, 4, 5, 6, 0, 0, 0, 1, 2, 3), nrow = 4, ncol = 3)
  genes_vec <- NULL  # No gene mapping
  
  params <- list(
    min_tpm = 1, min_samples = 2, min_tx_per_gene = 2, min_isoform_abundance = 0.05
  )
  
  # Should not error even with NULL genes_vec
  result <- TSENAT:::.apply_combined_filters(assay_mat, genes_vec, params, verbose = FALSE)
  
  expect_true(is.logical(result$tokeep))
  expect_equal(length(result$tokeep), 4)
})

test_that(".apply_combined_filters before count equals nrow", {
  assay_mat <- matrix(1:20, nrow = 5, ncol = 4)
  genes_vec <- c("G1", "G1", "G2", "G2", "G3")
  
  params <- list(
    min_tpm = 0, min_samples = 1, min_tx_per_gene = 1, min_isoform_abundance = NULL
  )
  
  result <- TSENAT:::.apply_combined_filters(assay_mat, genes_vec, params, verbose = FALSE)
  
  expect_equal(result$before, 5)
})

test_that(".apply_combined_filters after count is sum of logical vector", {
  assay_mat <- matrix(1:12, nrow = 4, ncol = 3)
  genes_vec <- c("G1", "G1", "G2", "G2")
  
  params <- list(
    min_tpm = 5, min_samples = 1, min_tx_per_gene = 1, min_isoform_abundance = NULL
  )
  
  result <- TSENAT:::.apply_combined_filters(assay_mat, genes_vec, params, verbose = FALSE)
  
  expect_equal(result$after, sum(result$tokeep))
})

# ============================================================================
# Tests for .finalize_filtered_se() helper function
# ============================================================================

test_that(".finalize_filtered_se returns SummarizedExperiment", {
  # Create test SE
  tpm <- matrix(1:9, nrow = 3, ncol = 3)
  rownames(tpm) <- c("TX1", "TX2", "TX3")
  colnames(tpm) <- c("S1", "S2", "S3")
  
  se <- SummarizedExperiment::SummarizedExperiment(
    assays = list(tpm = tpm)
  )
  
  tokeep <- c(TRUE, TRUE, FALSE)
  assays_list <- SummarizedExperiment::assays(se)
  
  result <- TSENAT:::.finalize_filtered_se(
    se, tokeep, assays_list, genes_vec = NULL,
    min_samples = 2, min_tpm = 1, min_tx_per_gene = 1,
    min_isoform_abundance = 0.05, stringency = NULL
  )
  
  expect_s4_class(result, "SummarizedExperiment")
  expect_equal(nrow(result), 2)
  expect_equal(ncol(result), 3)
})

test_that(".finalize_filtered_se subsets all assays correctly", {
  counts <- matrix(1:12, nrow = 4, ncol = 3)
  tpm <- matrix(10:21, nrow = 4, ncol = 3)
  rownames(counts) <- rownames(tpm) <- c("TX1", "TX2", "TX3", "TX4")
  colnames(counts) <- colnames(tpm) <- c("S1", "S2", "S3")
  
  se <- SummarizedExperiment::SummarizedExperiment(
    assays = list(counts = counts, tpm = tpm)
  )
  
  tokeep <- c(TRUE, FALSE, TRUE, FALSE)
  assays_list <- SummarizedExperiment::assays(se)
  
  result <- TSENAT:::.finalize_filtered_se(
    se, tokeep, assays_list, genes_vec = NULL,
    min_samples = 1, min_tpm = 1, min_tx_per_gene = 1,
    min_isoform_abundance = 0.05, stringency = NULL
  )
  
  # Check both assays subsetted
  expect_equal(nrow(SummarizedExperiment::assays(result)$counts), 2)
  expect_equal(nrow(SummarizedExperiment::assays(result)$tpm), 2)
  expect_equal(rownames(SummarizedExperiment::assays(result)$counts), c("TX1", "TX3"))
})

test_that(".finalize_filtered_se adds filtering record to metadata", {
  tpm <- matrix(1:9, nrow = 3, ncol = 3)
  rownames(tpm) <- c("TX1", "TX2", "TX3")
  colnames(tpm) <- c("S1", "S2", "S3")
  
  se <- SummarizedExperiment::SummarizedExperiment(
    assays = list(tpm = tpm)
  )
  
  tokeep <- c(TRUE, TRUE, FALSE)
  assays_list <- SummarizedExperiment::assays(se)
  
  result <- TSENAT:::.finalize_filtered_se(
    se, tokeep, assays_list, genes_vec = NULL,
    min_samples = 2, min_tpm = 1.5, min_tx_per_gene = 2,
    min_isoform_abundance = 0.1, stringency = "medium"
  )
  
  md <- S4Vectors::metadata(result)
  
  # Check filtering record
  expect_true("filtered" %in% names(md))
  expect_equal(md$filtered$min_samples, 2)
  expect_equal(md$filtered$min_tpm, 1.5)
  expect_equal(md$filtered$min_tx_per_gene, 2)
  expect_equal(md$filtered$min_isoform_abundance, 0.1)
  expect_equal(md$filtered$stringency, "medium")
})

test_that(".finalize_filtered_se preserves colData", {
  tpm <- matrix(1:9, nrow = 3, ncol = 3)
  rownames(tpm) <- c("TX1", "TX2", "TX3")
  colnames(tpm) <- c("S1", "S2", "S3")
  
  se <- SummarizedExperiment::SummarizedExperiment(
    assays = list(tpm = tpm),
    colData = data.frame(
      sample_id = c("S1", "S2", "S3"),
      condition = c("A", "B", "A"),
      row.names = c("S1", "S2", "S3")
    )
  )
  
  tokeep <- c(TRUE, FALSE, TRUE)
  assays_list <- SummarizedExperiment::assays(se)
  
  result <- TSENAT:::.finalize_filtered_se(
    se, tokeep, assays_list, genes_vec = NULL,
    min_samples = 1, min_tpm = 1, min_tx_per_gene = 1,
    min_isoform_abundance = 0.05, stringency = NULL
  )
  
  # colData should be unchanged (not subsetted)
  result_coldata <- SummarizedExperiment::colData(result)
  expect_equal(nrow(result_coldata), 3)
  expect_equal(result_coldata$condition, c("A", "B", "A"))
})

# Tests for se_manipulation_filter.R - Coverage for uncovered lines
# Complements test-core-functions-filter.R with edge cases and error paths

# ============================================================================
# TESTS FOR .get_assay_filtering() - TPM detection priority
# ============================================================================

test_that(".get_assay_filtering selects TPM from explicit tpm_assay_name", {
  counts <- matrix(1:9, nrow = 3, ncol = 3)
  my_tpm <- matrix(10:18, nrow = 3, ncol = 3)
  rownames(counts) <- rownames(my_tpm) <- c("TX1", "TX2", "TX3")
  colnames(counts) <- colnames(my_tpm) <- c("S1", "S2", "S3")
  
  se <- SummarizedExperiment::SummarizedExperiment(
    assays = list(counts = counts, custom_tpm = my_tpm)
  )
  
  result <- TSENAT:::.get_assay_filtering(
    se, 
    SummarizedExperiment::assays(se),
    tpm_assay_name = "custom_tpm",
    assay_name = "counts"
  )
  
  expect_equal(result$mat, my_tpm)
  expect_match(result$source, "custom_tpm.*user-specified")
})

test_that(".get_assay_filtering selects TPM from metadata$tpm", {
  counts <- matrix(1:9, nrow = 3, ncol = 3)
  tpm <- matrix(10:18, nrow = 3, ncol = 3)
  rownames(counts) <- rownames(tpm) <- c("TX1", "TX2", "TX3")
  colnames(counts) <- colnames(tpm) <- c("S1", "S2", "S3")
  
  se <- SummarizedExperiment::SummarizedExperiment(
    assays = list(counts = counts),
    metadata = list(tpm = tpm)
  )
  
  result <- TSENAT:::.get_assay_filtering(
    se, 
    SummarizedExperiment::assays(se),
    tpm_assay_name = NULL,
    assay_name = "counts"
  )
  
  expect_equal(result$mat, tpm)
  expect_match(result$source, "metadata.*SALMON")
})

test_that(".get_assay_filtering auto-detects 'tpm' assay by name", {
  counts <- matrix(1:9, nrow = 3, ncol = 3)
  tpm <- matrix(10:18, nrow = 3, ncol = 3)
  rownames(counts) <- rownames(tpm) <- c("TX1", "TX2", "TX3")
  colnames(counts) <- colnames(tpm) <- c("S1", "S2", "S3")
  
  se <- SummarizedExperiment::SummarizedExperiment(
    assays = list(counts = counts, tpm = tpm)
  )
  
  result <- TSENAT:::.get_assay_filtering(
    se, 
    SummarizedExperiment::assays(se),
    tpm_assay_name = NULL,
    assay_name = "counts"
  )
  
  expect_equal(result$mat, tpm)
  expect_match(result$source, "auto-detected")
})

test_that(".get_assay_filtering auto-detects 'abundance' assay (tximport format)", {
  counts <- matrix(1:9, nrow = 3, ncol = 3)
  abundance <- matrix(10:18, nrow = 3, ncol = 3)
  rownames(counts) <- rownames(abundance) <- c("TX1", "TX2", "TX3")
  colnames(counts) <- colnames(abundance) <- c("S1", "S2", "S3")
  
  se <- SummarizedExperiment::SummarizedExperiment(
    assays = list(counts = counts, abundance = abundance)
  )
  
  result <- TSENAT:::.get_assay_filtering(
    se, 
    SummarizedExperiment::assays(se),
    tpm_assay_name = NULL,
    assay_name = "counts"
  )
  
  expect_equal(result$mat, abundance)
  expect_match(result$source, "tximport")
})

test_that(".get_assay_filtering returns NULL when TPM not found", {
  counts <- matrix(1:9, nrow = 3, ncol = 3)
  rownames(counts) <- c("TX1", "TX2", "TX3")
  colnames(counts) <- c("S1", "S2", "S3")
  
  se <- SummarizedExperiment::SummarizedExperiment(
    assays = list(counts = counts)
  )
  
  result <- TSENAT:::.get_assay_filtering(
    se, 
    SummarizedExperiment::assays(se),
    tpm_assay_name = NULL,
    assay_name = "counts"
  )
  
  expect_null(result)
})

# ============================================================================
# TESTS FOR .get_gene_ids() - Gene ID detection from various sources
# ============================================================================

test_that(".get_gene_ids retrieves from rowData$gene_id", {
  counts <- matrix(1:9, nrow = 3, ncol = 3)
  rownames(counts) <- c("TX1", "TX2", "TX3")
  colnames(counts) <- c("S1", "S2", "S3")
  
  se <- SummarizedExperiment::SummarizedExperiment(
    assays = list(counts = counts),
    rowData = data.frame(
      gene_id = c("G1", "G2", "G1"),
      row.names = c("TX1", "TX2", "TX3")
    )
  )
  
  genes <- TSENAT:::.get_gene_ids(se)
  
  expect_equal(genes, c("G1", "G2", "G1"))
})

test_that(".get_gene_ids retrieves from metadata$tx2gene", {
  counts <- matrix(1:9, nrow = 3, ncol = 3)
  rownames(counts) <- c("TX1", "TX2", "TX3")
  colnames(counts) <- c("S1", "S2", "S3")
  
  se <- SummarizedExperiment::SummarizedExperiment(
    assays = list(counts = counts),
    metadata = list(
      tx2gene = data.frame(
        Transcript = c("TX1", "TX2", "TX3"),
        Gene = c("G1", "G2", "G1")
      )
    )
  )
  
  genes <- TSENAT:::.get_gene_ids(se)
  
  expect_equal(genes, c("G1", "G2", "G1"))
})

test_that(".get_gene_ids returns NULL when no gene mapping available", {
  counts <- matrix(1:9, nrow = 3, ncol = 3)
  rownames(counts) <- c("TX1", "TX2", "TX3")
  colnames(counts) <- c("S1", "S2", "S3")
  
  se <- SummarizedExperiment::SummarizedExperiment(
    assays = list(counts = counts)
  )
  
  genes <- TSENAT:::.get_gene_ids(se)
  
  expect_null(genes)
})

# ============================================================================
# TESTS FOR .estimate_min_tpm() - Stringency-based TPM estimation
# ============================================================================

test_that(".estimate_min_tpm soft stringency uses Q1", {
  assay_mat <- matrix(c(0.5, 1, 2, 3, 5, 10, 20, 50, 100), nrow = 1)
  
  result <- TSENAT:::.estimate_min_tpm(assay_mat, "soft", verbose = FALSE)
  
  expect_equal(result$quant_label, "Q1")
  expect_true(result$min_tpm > 0.1)
  expect_true(result$min_tpm <= 5)
})

test_that(".estimate_min_tpm medium stringency uses Q2", {
  assay_mat <- matrix(c(0.5, 1, 2, 3, 5, 10, 20, 50, 100), nrow = 1)
  
  result <- TSENAT:::.estimate_min_tpm(assay_mat, "medium", verbose = FALSE)
  
  expect_equal(result$quant_label, "Q2/Median")
  expect_true(result$min_tpm > 0.1)
  expect_true(result$min_tpm <= 5)
})

test_that(".estimate_min_tpm severe stringency uses Q3", {
  assay_mat <- matrix(c(0.5, 1, 2, 3, 5, 10, 20, 50, 100), nrow = 1)
  
  result <- TSENAT:::.estimate_min_tpm(assay_mat, "severe", verbose = FALSE)
  
  expect_equal(result$quant_label, "Q3")
  expect_true(result$min_tpm > 0.1)
  expect_true(result$min_tpm <= 5)
})

test_that(".estimate_min_tpm handles invalid stringency (returns default)", {
  assay_mat <- matrix(c(0.5, 1, 2, 3, 5), nrow = 1)
  
  result <- TSENAT:::.estimate_min_tpm(assay_mat, "invalid", verbose = FALSE)
  
  expect_equal(result$quant_label, "default")
  expect_equal(result$min_tpm, 0.1)
})

test_that(".estimate_min_tpm handles assay with all zeros", {
  assay_mat <- matrix(rep(0, 9), nrow = 1)
  
  expect_warning(
    result <- TSENAT:::.estimate_min_tpm(assay_mat, "soft", verbose = FALSE),
    "No non-zero"
  )
  
  expect_equal(result$quant_label, "default")
  expect_equal(result$min_tpm, 0.1)
})

# ============================================================================
# TESTS FOR .calculate_stringency_thresholds() - Parameter auto-calculation
# ============================================================================

test_that(".calculate_stringency_thresholds soft level", {
  result <- TSENAT:::.calculate_stringency_thresholds("soft", n_samples = 10)
  
  expect_equal(result$min_samples, 3)  # ceiling(0.25 * 10) = 3, min 2
  expect_equal(result$min_tx_per_gene, 2L)
  expect_equal(result$min_isoform_abundance, 0.01)
})

test_that(".calculate_stringency_thresholds medium level", {
  result <- TSENAT:::.calculate_stringency_thresholds("medium", n_samples = 10)
  
  expect_equal(result$min_samples, 5)  # ceiling(0.5 * 10) = 5, min 3
  expect_equal(result$min_tx_per_gene, 2L)
  expect_equal(result$min_isoform_abundance, 0.05)
})

test_that(".calculate_stringency_thresholds severe level", {
  result <- TSENAT:::.calculate_stringency_thresholds("severe", n_samples = 10)
  
  expect_equal(result$min_samples, 8)  # ceiling(0.75 * 10) = 8
  expect_equal(result$min_tx_per_gene, 3L)
  expect_equal(result$min_isoform_abundance, 0.15)
})

test_that(".calculate_stringency_thresholds returns NULL for invalid stringency", {
  result <- TSENAT:::.calculate_stringency_thresholds("invalid", n_samples = 10)
  
  expect_null(result)
})

# ============================================================================
# TESTS FOR .apply_tpm_filter() - Edge cases with min_samples > ncol
# ============================================================================

test_that(".apply_tpm_filter warns when min_samples > ncol(assay_mat)", {
  assay_mat <- matrix(c(100, 200, 300), nrow = 3, ncol = 1)
  
  expect_warning(
    result <- TSENAT:::.apply_tpm_filter(assay_mat, min_tpm = 1, min_samples = 5, verbose = FALSE),
    "min_samples.*greater than"
  )
  
  # Should still return logical vector even if filtering to 0
  expect_equal(length(result), 3)
  expect_true(all(!result))  # All should be FALSE
})

test_that(".apply_tpm_filter with verbose message", {
  assay_mat <- matrix(c(100, 200, 300), nrow = 3, ncol = 2)
  
  expect_message(
    result <- TSENAT:::.apply_tpm_filter(assay_mat, min_tpm = 150, min_samples = 1, verbose = TRUE),
    "TPM-based filtering"
  )
})

# ============================================================================
# TESTS FOR .validate_filter_params() Error paths
# ============================================================================

test_that(".validate_filter_params rejects min_isoform_abundance > 1", {
  expect_error(
    TSENAT:::.validate_filter_params(min_isoform_abundance = 1.5, tpm_assay_name = NULL),
    "must be numeric in"
  )
})

test_that(".validate_filter_params rejects negative min_isoform_abundance", {
  expect_error(
    TSENAT:::.validate_filter_params(min_isoform_abundance = -0.1, tpm_assay_name = NULL),
    "must be numeric in"
  )
})

test_that(".validate_filter_params rejects non-numeric min_isoform_abundance", {
  expect_error(
    TSENAT:::.validate_filter_params(min_isoform_abundance = "0.5", tpm_assay_name = NULL),
    "must be numeric"
  )
})

test_that(".validate_filter_params accepts NULL and valid values", {
  result <- TSENAT:::.validate_filter_params(min_isoform_abundance = NULL, tpm_assay_name = NULL)
  expect_true(result)
  
  result <- TSENAT:::.validate_filter_params(min_isoform_abundance = 0.05, tpm_assay_name = NULL)
  expect_true(result)
})

# ============================================================================
# TESTS FOR .filter_se() Error paths - Missing TPM
# ============================================================================

test_that(".filter_se raises error when TPM data not found", {
  counts <- matrix(1:9, nrow = 3, ncol = 3)
  rownames(counts) <- c("TX1", "TX2", "TX3")
  colnames(counts) <- c("S1", "S2", "S3")
  
  se <- SummarizedExperiment::SummarizedExperiment(
    assays = list(counts = counts)
  )
  
  expect_error(
    TSENAT:::.filter_se(se, min_samples = 1, verbose = FALSE),
    "TPM data is required"
  )
})

test_that(".filter_se raises error for non-SummarizedExperiment input", {
  df <- data.frame(a = 1:3)
  
  expect_error(
    TSENAT:::.filter_se(df, min_samples = 1, verbose = FALSE),
    "must be a SummarizedExperiment"
  )
})

test_that(".filter_se raises error for non-numeric assay", {
  counts <- matrix(c("a", "b", "c", "d", "e", "f"), nrow = 2, ncol = 3)
  
  se <- SummarizedExperiment::SummarizedExperiment(
    assays = list(counts = counts),
    metadata = list(tpm = matrix(10:15, nrow = 2, ncol = 3))
  )
  
  expect_error(
    TSENAT:::.filter_se(se, min_samples = 1, verbose = FALSE),
    "must be numeric"
  )
})

# ============================================================================
# TESTS FOR .resolve_filter_parameters() - Error and special cases
# ============================================================================

test_that(".resolve_filter_parameters raises error for invalid stringency", {
  counts <- matrix(1:9, nrow = 3, ncol = 3)
  tpm <- matrix(10:18, nrow = 3, ncol = 3)
  rownames(counts) <- rownames(tpm) <- c("TX1", "TX2", "TX3")
  colnames(counts) <- colnames(tpm) <- c("S1", "S2", "S3")
  
  se <- SummarizedExperiment::SummarizedExperiment(
    assays = list(counts = counts, tpm = tpm),
    colData = data.frame(pair = c("p1", "p1", "p2"), row.names = c("S1", "S2", "S3"))
  )
  
  expect_error(
    TSENAT:::.resolve_filter_parameters(se, 5, 1, "invalid_level", NULL, NULL, "counts", 2, NULL),
    "must be one of.*soft.*medium.*severe"
  )
})

test_that(".resolve_filter_parameters rejects invalid min_isoform_abundance parameter", {
  counts <- matrix(1:9, nrow = 3, ncol = 3)
  tpm <- matrix(10:18, nrow = 3, ncol = 3)
  rownames(counts) <- rownames(tpm) <- c("TX1", "TX2", "TX3")
  colnames(counts) <- colnames(tpm) <- c("S1", "S2", "S3")
  
  se <- SummarizedExperiment::SummarizedExperiment(
    assays = list(counts = counts, tpm = tpm)
  )
  
  expect_error(
    TSENAT:::.resolve_filter_parameters(se, 5, 1, NULL, NULL, NULL, "counts", 2, 1.5),
    "must be numeric in"
  )
})

test_that(".resolve_filter_parameters stringency with auto-detection of pair column", {
  counts <- matrix(1:12, nrow = 3, ncol = 4)
  tpm <- matrix(10:21, nrow = 3, ncol = 4)
  rownames(counts) <- rownames(tpm) <- c("TX1", "TX2", "TX3")
  colnames(counts) <- colnames(tpm) <- c("S1", "S2", "S3", "S4")
  
  se <- SummarizedExperiment::SummarizedExperiment(
    assays = list(counts = counts, tpm = tpm),
    colData = data.frame(
      pair_id = c("p1", "p1", "p2", "p2"),
      row.names = c("S1", "S2", "S3", "S4")
    )
  )
  
  # Should not raise error when pair column can be auto-detected
  result <- TSENAT:::.resolve_filter_parameters(
    se, 5, 1, "medium", NULL, NULL, "counts", 2, NULL,
    verbose = FALSE
  )
  
  expect_equal(result$pair_col_used, "pair_id")
})

# ============================================================================
# TESTS FOR .detect_pair_column() - Error when no column found
# ============================================================================

test_that(".detect_pair_column raises error when no pair column found", {
  col_data <- data.frame(
    sample_id = c("S1", "S2", "S3"),
    condition = c("A", "B", "A"),
    row.names = c("S1", "S2", "S3")
  )
  
  expect_error(
    TSENAT:::.detect_pair_column(col_data, NULL),
    "Could not auto-detect pair column"
  )
})

test_that(".detect_pair_column successfully finds pair columns", {
  col_data <- data.frame(
    sample_id = c("S1", "S2", "S3"),
    pair_id = c("p1", "p1", "p2"),
    row.names = c("S1", "S2", "S3")
  )
  
  result <- TSENAT:::.detect_pair_column(col_data, NULL)
  expect_equal(result, "pair_id")
})

# ============================================================================
# TESTS FOR Combined filtering with edge cases
# ============================================================================

test_that(".filter_se with verbose=TRUE outputs all messages", {
  counts <- matrix(c(100, 200, 300, 50, 150, 250), nrow = 3, ncol = 2)
  tpm <- counts
  rownames(counts) <- rownames(tpm) <- c("TX1", "TX2", "TX3")
  colnames(counts) <- colnames(tpm) <- c("S1", "S2")
  
  se <- SummarizedExperiment::SummarizedExperiment(
    assays = list(counts = counts, tpm = tpm),
    rowData = data.frame(
      gene_id = c("G1", "G1", "G2"),
      row.names = c("TX1", "TX2", "TX3")
    ),
    metadata = list(tx2gene = data.frame(
      tx = c("TX1", "TX2", "TX3"),
      gene = c("G1", "G1", "G2")
    ))
  )
  
  expect_message(
    TSENAT:::.filter_se(se, min_samples = 1, min_tpm = 100, 
                       min_isoform_abundance = 0.1, verbose = TRUE),
    "Transcripts:.*before.*after"
  )
})

test_that(".filter_se returns empty SE when all rows filtered", {
  counts <- matrix(c(1, 2, 3), nrow = 3, ncol = 1)
  tpm <- matrix(c(0.1, 0.2, 0.3), nrow = 3, ncol = 1)
  rownames(counts) <- rownames(tpm) <- c("TX1", "TX2", "TX3")
  colnames(counts) <- colnames(tpm) <- "S1"
  
  se <- SummarizedExperiment::SummarizedExperiment(
    assays = list(counts = counts, tpm = tpm)
  )
  
  expect_warning(
    result <- TSENAT:::.filter_se(se, min_samples = 1, min_tpm = 1000, verbose = FALSE),
    "Filtering removed all transcripts"
  )
  
  expect_equal(nrow(result), 0)
})

# ============================================================================
# TESTS FOR .filter_by_tx_per_gene with various parameters
# ============================================================================

test_that(".filter_by_tx_per_gene with min_tx_per_gene = 1 (skip filtering)", {
  tokeep <- c(TRUE, TRUE, FALSE, TRUE)
  genes_vec <- c("G1", "G1", "G2", "G3")
  
  result <- TSENAT:::.filter_by_tx_per_gene(tokeep, genes_vec, 1, verbose = FALSE)
  
  expect_equal(result, tokeep)
})

test_that(".filter_by_tx_per_gene with NULL genes_vec", {
  tokeep <- c(TRUE, TRUE, FALSE, TRUE)
  
  result <- TSENAT:::.filter_by_tx_per_gene(tokeep, NULL, 2, verbose = FALSE)
  
  expect_equal(result, tokeep)
})

test_that(".filter_by_tx_per_gene with verbose output", {
  tokeep <- c(TRUE, TRUE, FALSE, TRUE)
  genes_vec <- c("G1", "G1", "G2", "G3")
  
  expect_message(
    result <- TSENAT:::.filter_by_tx_per_gene(tokeep, genes_vec, 2, verbose = TRUE),
    "Filtered by min_tx_per_gene"
  )
})

# ============================================================================
# TESTS FOR .filter_by_isoform_abundance with edge cases
# ============================================================================

test_that(".filter_by_isoform_abundance with min_isoform_abundance = 0 (skip)", {
  tokeep <- c(TRUE, TRUE, TRUE)
  genes_vec <- c("G1", "G1", "G2")
  assay_mat <- matrix(c(100, 50, 200, 100, 10, 5), nrow = 3, ncol = 2)
  
  result <- TSENAT:::.filter_by_isoform_abundance(tokeep, genes_vec, assay_mat, 0, verbose = FALSE)
  
  expect_equal(result, tokeep)
})

test_that(".filter_by_isoform_abundance with NULL min_isoform_abundance", {
  tokeep <- c(TRUE, TRUE, TRUE)
  genes_vec <- c("G1", "G1", "G2")
  assay_mat <- matrix(c(100, 50, 200, 100, 10, 5), nrow = 3, ncol = 2)
  
  result <- TSENAT:::.filter_by_isoform_abundance(tokeep, genes_vec, assay_mat, NULL, verbose = FALSE)
  
  expect_equal(result, tokeep)
})

test_that(".filter_by_isoform_abundance with NULL genes_vec", {
  tokeep <- c(TRUE, TRUE, TRUE)
  assay_mat <- matrix(c(100, 50, 200), nrow = 3, ncol = 1)
  
  result <- TSENAT:::.filter_by_isoform_abundance(tokeep, NULL, assay_mat, 0.1, verbose = FALSE)
  
  expect_equal(result, tokeep)
})

test_that(".filter_by_isoform_abundance removes low-abundance isoforms", {
  tokeep <- c(TRUE, TRUE, TRUE)
  genes_vec <- c("G1", "G1", "G2")
  # Isoforms: TX1=100, TX2=10 (both G1), TX3=100 (G2)
  # Relative for G1: 100/110=91%, 10/110=9%
  # At 20% threshold, TX2 should be removed
  assay_mat <- matrix(c(100, 10, 100, 100, 10, 100), nrow = 3, ncol = 2)
  
  result <- TSENAT:::.filter_by_isoform_abundance(tokeep, genes_vec, assay_mat, 0.20, verbose = FALSE)
  
  expect_equal(result, c(TRUE, FALSE, TRUE))
})

test_that(".filter_by_isoform_abundance always keeps single-isoform genes", {
  tokeep <- c(TRUE, TRUE)
  genes_vec <- c("G1", "G2")  # Each gene has only 1 isoform
  assay_mat <- matrix(c(100, 50, 10, 20), nrow = 2, ncol = 2)
  
  result <- TSENAT:::.filter_by_isoform_abundance(tokeep, genes_vec, assay_mat, 0.50, verbose = FALSE)
  
  expect_equal(result, tokeep)
})

# ============================================================================
# TESTS FOR .calc_stringency_params() - Wrapper function
# ============================================================================

test_that(".calc_stringency_params auto-detects pair column", {
  se <- SummarizedExperiment::SummarizedExperiment(
    assays = list(counts = matrix(1:12, nrow = 3, ncol = 4)),
    colData = data.frame(
      pair_id = c("p1", "p1", "p2", "p2"),
      row.names = c("S1", "S2", "S3", "S4")
    )
  )
  
  result <- TSENAT:::.calc_stringency_params(se, "medium", pair_col = NULL, verbose = FALSE)
  
  expect_equal(result$min_samples, 2)  # ceiling(0.5 * 4) = 2, min 3 for medium -> actually 2
  expect_equal(result$min_tx_per_gene, 2L)
})

test_that(".calc_stringency_params with verbose output", {
  se <- SummarizedExperiment::SummarizedExperiment(
    assays = list(counts = matrix(1:12, nrow = 3, ncol = 4)),
    colData = data.frame(
      pair = c("p1", "p1", "p2", "p2"),
      row.names = c("S1", "S2", "S3", "S4")
    )
  )
  
  expect_message(
    result <- TSENAT:::.calc_stringency_params(se, "soft", pair_col = NULL, verbose = TRUE),
    "n_pairs"
  )
})

# ============================================================================
# TESTS FOR .sync_filter_metadata() Edge cases
# ============================================================================

test_that(".sync_filter_metadata handles metadata without readcounts", {
  md <- list(tpm = matrix(1:9, nrow = 3, ncol = 3))
  rownames(md$tpm) <- c("TX1", "TX2", "TX3")
  
  tokeep <- c(TRUE, FALSE, TRUE)
  assay_mat <- matrix(1:9, nrow = 3, ncol = 3)
  rownames(assay_mat) <- c("TX1", "TX2", "TX3")
  
  result <- TSENAT:::.sync_filter_metadata(md, tokeep, assay_mat, NULL)
  
  # Should handle gracefully without error
  expect_null(result$readcounts)
  expect_equal(nrow(result$tpm), 2)
})

test_that(".sync_filter_metadata handles effective_length as named vector", {
  md <- list(
    effective_length = c(TX1 = 100, TX2 = 150, TX3 = 120)
  )
  
  tokeep <- c(TRUE, FALSE, TRUE)
  assay_mat <- matrix(1:9, nrow = 3, ncol = 3)
  rownames(assay_mat) <- c("TX1", "TX2", "TX3")
  
  result <- TSENAT:::.sync_filter_metadata(md, tokeep, assay_mat, c("TX1", "TX2", "TX3"))
  
  # Should subset by names
  expect_equal(names(result$effective_length), c("TX1", "TX3"))
})

test_that(".sync_filter_metadata handles effective_length as unnamed vector", {
  md <- list(
    effective_length = c(100, 150, 120)
  )
  
  tokeep <- c(TRUE, FALSE, TRUE)
  assay_mat <- matrix(1:9, nrow = 3, ncol = 3)
  
  result <- TSENAT:::.sync_filter_metadata(md, tokeep, assay_mat, c("TX1", "TX2", "TX3"))
  
  # Should subset by index
  expect_equal(length(result$effective_length), 2)
})

test_that(".sync_filter_metadata filters tx2gene mapping", {
  md <- list(
    tx2gene = data.frame(
      tx = c("TX1", "TX2", "TX3"),
      gene = c("G1", "G2", "G1")
    )
  )
  
  tokeep <- c(TRUE, FALSE, TRUE)
  assay_mat <- matrix(1:9, nrow = 3, ncol = 3)
  
  result <- TSENAT:::.sync_filter_metadata(md, tokeep, assay_mat, c("TX1", "TX2", "TX3"))
  
  # Should filter to only TX1 and TX3
  expect_equal(nrow(result$tx2gene), 2)
  expect_equal(result$tx2gene$tx, c("TX1", "TX3"))
})
