# Comprehensive tests for Salmon I/O functions and build_analysis_s4 integration
# Tests for: .detect_salmon_samples(), .validate_salmon_files(), .read_salmon_samples()
# And: build_analysis_s4() with salmon_dir and salmon data parameters

# ===========================================================================
# SETUP HELPERS: Create mock Salmon output structure
# ===========================================================================

#' Create mock Salmon quant.sf file for testing
#' @noRd
create_mock_salmon_file <- function(filepath, n_transcripts = 100, seed = NULL, high_tpm = FALSE) {
  if (!is.null(seed)) set.seed(seed)
  
  # Generate realistic Salmon-style data with correlated properties
  # All columns represent realistic sequencing output
  
  # 1. Transcript lengths (500-5000 bp, realistic for mRNA)
  lengths <- sample(500:5000, n_transcripts, replace = TRUE)
  
  # 2. Effective length (typically 95-98% of Length due to k-mer bias and fragment size)
  effective_lengths <- round(lengths * runif(n_transcripts, 0.90, 0.98))
  
  # 3. NumReads (Poisson distribution, realistic sequencing depth)
  num_reads <- rpois(n_transcripts, lambda = 75)  # Mean ~75 reads per transcript
  
  # 4. TPM: derived from counts and effective length (realistic correlation)
  # TPM = (NumReads / EffectiveLength) * (1e6 / sum(NumReads / EffectiveLength))
  # Add realistic variation for expression heterogeneity
  tpm_base <- (num_reads / effective_lengths) * 1e6  # Unnormalized
  tpm_vals <- tpm_base / sum(tpm_base) * 1e6  # Normalized to ~1e6
  tpm_vals <- tpm_vals * runif(n_transcripts, 0.5, 2.0)  # Random fold variation
  tpm_vals <- pmax(tpm_vals, 0.1)  # Ensure all >= 0.1 (minimum detection)
  
  data <- data.frame(
    Name = paste0("ENST", sprintf("%011d", seq_len(n_transcripts))),
    Length = lengths,
    EffectiveLength = effective_lengths,
    TPM = tpm_vals,
    NumReads = num_reads
  )
  
  # Ensure directory exists
  dir.create(dirname(filepath), showWarnings = FALSE, recursive = TRUE)
  
  # Write TSV file
  readr::write_tsv(data, filepath)
  invisible(data)
}

#' Create mock Salmon directory structure
#' @noRd
create_mock_salmon_dir <- function(tmpdir, sample_names = c("sample1", "sample2"), 
                                   n_transcripts = 100, same_transcripts = TRUE, high_tpm = FALSE) {
  # Use unique subdirectory for each test to avoid conflicts
  salmon_base <- paste0("salmon_", gsub(":", "-", Sys.time()))
  salmon_dir <- file.path(tmpdir, salmon_base)
  dir.create(salmon_dir, showWarnings = FALSE, recursive = TRUE)
  
  file_paths <- character(length(sample_names))
  
  for (i in seq_along(sample_names)) {
    sample_dir <- file.path(salmon_dir, sample_names[i])
    dir.create(sample_dir, showWarnings = FALSE, recursive = TRUE)
    
    quant_file <- file.path(sample_dir, "quant.sf")
    
    # Use same seed for first sample to ensure identical transcripts across files
    seed_val <- if (same_transcripts) 42 else (42 + i)
    create_mock_salmon_file(quant_file, n_transcripts = n_transcripts, seed = seed_val, high_tpm = high_tpm)
    
    file_paths[i] <- quant_file
  }
  
  list(
    salmon_dir = salmon_dir,
    file_paths = file_paths,
    sample_names = sample_names
  )
}

#' Create mock Salmon file with guaranteed high TPM values
#' @noRd
create_mock_salmon_file_high_tpm <- function(filepath, n_transcripts = 100, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  
  # Create data with consistently high TPM values for filtering tests
  data <- data.frame(
    Name = paste0("ENST", sprintf("%011d", seq_len(n_transcripts))),
    Length = sample(500:5000, n_transcripts, replace = TRUE),
    EffectiveLength = sample(400:4800, n_transcripts, replace = TRUE),
    TPM = runif(n_transcripts, 5, 50),  # Guarantees min 5 TPM
    NumReads = rpois(n_transcripts, lambda = 100)  # Increase read counts
  )
  
  dir.create(dirname(filepath), showWarnings = FALSE, recursive = TRUE)
  readr::write_tsv(data, filepath)
  invisible(data)
}

# ===========================================================================
# TESTS: .detect_salmon_samples()
# ===========================================================================

test_that(".detect_salmon_samples finds quant.sf files in standard structure", {
  with_mock_salmon <- function(code) {
    tmpdir <- tempdir()
    on.exit(unlink(file.path(tmpdir, "salmon"), recursive = TRUE))
    
    salmon_setup <- create_mock_salmon_dir(tmpdir, sample_names = c("S1", "S2", "S3"))
    
    code(salmon_setup$salmon_dir)
  }
  
  with_mock_salmon(function(salmon_dir) {
    result <- .detect_salmon_samples(salmon_dir)
    
    expect_is(result, "list")
    expect_named(result, c("sample_names", "file_paths", "count"))
    expect_equal(result$count, 3)
    expect_equal(length(result$sample_names), 3)
    expect_equal(length(result$file_paths), 3)
    expect_equal(result$sample_names, c("S1", "S2", "S3"))
  })
})

test_that(".detect_salmon_samples works with gzipped quant.sf.gz files", {
  tmpdir <- tempdir()
  on.exit(unlink(file.path(tmpdir, "salmon_gz"), recursive = TRUE))
  
  salmon_dir <- file.path(tmpdir, "salmon_gz", "sample1")
  dir.create(salmon_dir, recursive = TRUE, showWarnings = FALSE)
  
  # Create and gzip a quant.sf file
  data <- create_mock_salmon_file(file.path(salmon_dir, "quant.sf"))
  system(paste("gzip -f", file.path(salmon_dir, "quant.sf")))
  
  result <- .detect_salmon_samples(file.path(tmpdir, "salmon_gz"))
  
  expect_equal(result$count, 1)
  expect_true(grepl("quant.sf.gz$", result$file_paths[1]))
})

test_that(".detect_salmon_samples rejects non-existent directory", {
  expect_error(
    .detect_salmon_samples("/nonexistent/path/to/salmon"),
    "Directory does not exist"
  )
})

test_that(".detect_salmon_samples rejects directory with no quant.sf files", {
  tmpdir <- tempdir()
  on.exit(unlink(file.path(tmpdir, "empty_salmon"), recursive = TRUE))
  
  salmon_dir <- file.path(tmpdir, "empty_salmon")
  dir.create(file.path(salmon_dir, "sample1"), recursive = TRUE, showWarnings = FALSE)
  
  expect_error(
    .detect_salmon_samples(salmon_dir),
    "No Salmon quantification files found"
  )
})

test_that(".detect_salmon_samples extracts sample names correctly from paths", {
  tmpdir <- tempdir()
  on.exit(unlink(file.path(tmpdir, "salmon_extract"), recursive = TRUE))
  
  salmon_setup <- create_mock_salmon_dir(
    tmpdir, 
    sample_names = c("treatment_rep1", "treatment_rep2", "control_rep1")
  )
  
  result <- .detect_salmon_samples(salmon_setup$salmon_dir)
  
  # Sort both for order-independent comparison (list.files order is not guaranteed)
  expect_equal(
    sort(result$sample_names),
    sort(c("treatment_rep1", "treatment_rep2", "control_rep1"))
  )
})

test_that(".detect_salmon_samples warns on duplicate sample names", {
  tmpdir <- tempdir()
  on.exit(unlink(file.path(tmpdir, "salmon_dup"), recursive = TRUE))
  
  salmon_dir <- file.path(tmpdir, "salmon_dup")
  
  # Create nested structure with duplicate folder names
  dir.create(file.path(salmon_dir, "batch1", "sample1"), recursive = TRUE, showWarnings = FALSE)
  dir.create(file.path(salmon_dir, "batch2", "sample1"), recursive = TRUE, showWarnings = FALSE)
  
  create_mock_salmon_file(file.path(salmon_dir, "batch1", "sample1", "quant.sf"))
  create_mock_salmon_file(file.path(salmon_dir, "batch2", "sample1", "quant.sf"))
  
  expect_warning(
    .detect_salmon_samples(salmon_dir, recursive = TRUE),
    "duplicate sample names"
  )
})

# ===========================================================================
# TESTS: validate_salmon_files()
# ===========================================================================

test_that("validate_salmon_files passes for valid quant.sf files", {
  tmpdir <- tempdir()
  on.exit(unlink(file.path(tmpdir, "validate_test"), recursive = TRUE))
  
  salmon_setup <- create_mock_salmon_dir(tmpdir, sample_names = c("S1", "S2"))
  
  expect_true(.validate_salmon_files(salmon_setup$file_paths, verbose = FALSE))
})

test_that(".validate_salmon_files rejects non-existent files", {
  expect_error(
    .validate_salmon_files(c("/nonexistent/quant.sf", "/another/nonexistent.sf")),
    "File\\(s\\) not found"
  )
})

test_that("validate_salmon_files checks for required columns", {
  tmpdir <- tempdir()
  on.exit(unlink(file.path(tmpdir, "bad_format"), recursive = TRUE))
  
  bad_dir <- file.path(tmpdir, "bad_format", "sample1")
  dir.create(bad_dir, recursive = TRUE, showWarnings = FALSE)
  bad_file <- file.path(bad_dir, "quant.sf")
  
  # Create valid TSV but with missing required columns
  # Write file with only 3 of 5 required columns
  bad_content <- "Name\tLength\tTPM\n"
  bad_content <- paste0(bad_content, "ENST00000000001\t1000\t5.5\n")
  bad_content <- paste0(bad_content, "ENST00000000002\t1001\t3.2\n")
  writeLines(bad_content, bad_file)
  
  # Expect both the warning from readr about missing column specs and our validation error
  expect_warning(
    expect_error(
      .validate_salmon_files(bad_file, verbose = FALSE),
      "Missing required columns"
    ),
    "named parsers don't match"
  )
})

test_that(".validate_salmon_files detects transcript ID mismatches", {
  tmpdir <- tempdir()
  on.exit(unlink(file.path(tmpdir, "mismatch_test"), recursive = TRUE))
  
  # Create two files with different transcripts
  dir1 <- file.path(tmpdir, "mismatch_test", "s1")
  dir2 <- file.path(tmpdir, "mismatch_test", "s2")
  dir.create(dir1, recursive = TRUE, showWarnings = FALSE)
  dir.create(dir2, recursive = TRUE, showWarnings = FALSE)
  
  # First file
  file1_data <- data.frame(
    Name = paste0("ENST", 1:10),
    Length = 1000:1009,
    EffectiveLength = 900:909,
    TPM = runif(10),
    NumReads = rpois(10, 50)
  )
  readr::write_tsv(file1_data, file.path(dir1, "quant.sf"))
  
  # Second file with different transcripts
  file2_data <- data.frame(
    Name = paste0("ENST", 5:14),  # Different range
    Length = 1000:1009,
    EffectiveLength = 900:909,
    TPM = runif(10),
    NumReads = rpois(10, 50)
  )
  readr::write_tsv(file2_data, file.path(dir2, "quant.sf"))
  
  expect_warning(
    .validate_salmon_files(
      c(file.path(dir1, "quant.sf"), file.path(dir2, "quant.sf")),
      verbose = FALSE
    ),
    "Transcript ID mismatch"
  )
})

test_that(".validate_salmon_files rejects empty files", {
  tmpdir <- tempdir()
  on.exit(unlink(file.path(tmpdir, "empty_test"), recursive = TRUE))
  
  empty_dir <- file.path(tmpdir, "empty_test", "sample1")
  dir.create(empty_dir, recursive = TRUE, showWarnings = FALSE)
  empty_file <- file.path(empty_dir, "quant.sf")
  
  # Create file with header but no data
  empty_data <- data.frame(
    Name = character(),
    Length = integer(),
    EffectiveLength = numeric(),
    TPM = numeric(),
    NumReads = numeric()
  )
  readr::write_tsv(empty_data, empty_file)
  
  expect_error(
    .validate_salmon_files(empty_file, verbose = FALSE),
    "empty"
  )
})

# ===========================================================================
# TESTS: read_salmon_samples()
# ===========================================================================

test_that("read_salmon_samples builds count matrix from quant.sf files", {
  tmpdir <- tempdir()
  on.exit(unlink(file.path(tmpdir, "read_test"), recursive = TRUE))
  
  salmon_setup <- create_mock_salmon_dir(tmpdir, sample_names = c("S1", "S2", "S3"))
  
  result <- .read_salmon_samples(
    file_paths = salmon_setup$file_paths,
    sample_names = salmon_setup$sample_names,
    include_tpm = FALSE,
    include_eff_length = FALSE,
    verbose = FALSE
  )
  
  expect_is(result, "list")
  expect_named(result, c("counts", "transcript_ids"))
  expect_equal(nrow(result$counts), 100)
  expect_equal(ncol(result$counts), 3)
  expect_equal(colnames(result$counts), c("S1", "S2", "S3"))
  expect_equal(result$transcript_ids, rownames(result$counts))
})

test_that(".read_salmon_samples includes TPM when requested", {
  tmpdir <- tempdir()
  on.exit(unlink(file.path(tmpdir, "tpm_test"), recursive = TRUE))
  
  salmon_setup <- create_mock_salmon_dir(tmpdir, sample_names = c("S1", "S2"))
  
  result <- .read_salmon_samples(
    file_paths = salmon_setup$file_paths,
    sample_names = salmon_setup$sample_names,
    include_tpm = TRUE,
    include_eff_length = FALSE,
    verbose = FALSE
  )
  
  expect_is(result$tpm, "matrix")
  expect_equal(dim(result$tpm), c(100, 2))
  expect_true(all(result$tpm >= 0))  # TPM should be non-negative
})

test_that(".read_salmon_samples includes effective length when requested", {
  tmpdir <- tempdir()
  on.exit(unlink(file.path(tmpdir, "efflen_test"), recursive = TRUE))
  
  salmon_setup <- create_mock_salmon_dir(tmpdir, sample_names = c("S1", "S2"))
  
  result <- .read_salmon_samples(
    file_paths = salmon_setup$file_paths,
    sample_names = salmon_setup$sample_names,
    include_tpm = FALSE,
    include_eff_length = TRUE,
    verbose = FALSE
  )
  
  expect_is(result$effective_length, "matrix")
  expect_equal(dim(result$effective_length), c(100, 2))
  expect_true(all(result$effective_length > 0))  # Effective length should be positive
})

test_that(".read_salmon_samples infers sample names from file paths", {
  tmpdir <- tempdir()
  on.exit(unlink(file.path(tmpdir, "infer_names"), recursive = TRUE))
  
  salmon_setup <- create_mock_salmon_dir(
    tmpdir,
    sample_names = c("treatment_r1", "control_r1")
  )
  
  # Don't provide sample_names
  result <- .read_salmon_samples(
    file_paths = salmon_setup$file_paths,
    sample_names = NULL,
    verbose = FALSE
  )
  
  expect_equal(colnames(result$counts), c("treatment_r1", "control_r1"))
})

test_that(".read_salmon_samples rejects mismatched file_paths and sample_names", {
  tmpdir <- tempdir()
  on.exit(unlink(file.path(tmpdir, "mismatch"), recursive = TRUE))
  
  salmon_setup <- create_mock_salmon_dir(tmpdir, sample_names = c("S1", "S2"))
  
  expect_error(
    .read_salmon_samples(
      file_paths = salmon_setup$file_paths,
      sample_names = c("S1"),  # Only one name for two files
      verbose = FALSE
    ),
    "Length mismatch"
  )
})

test_that(".read_salmon_samples rejects empty file_paths", {
  expect_error(
    .read_salmon_samples(file_paths = character(), verbose = FALSE),
    "non-empty character vector"
  )
})

# ===========================================================================
# TESTS: build_analysis_s4() with Salmon integration
# ===========================================================================

test_that("build_analysis_s4 works with salmon_dir parameter", {
  skip_if_not_installed("TSENAT")
  
  tmpdir <- tempdir()
  on.exit(unlink(file.path(tmpdir, "analysis_dir"), recursive = TRUE))
  
  salmon_setup <- create_mock_salmon_dir(tmpdir, sample_names = c("S1", "S2", "S3"))
  
  # Create metadata
  metadata <- data.frame(
    sample = salmon_setup$sample_names,
    condition = rep(c("control", "treatment"), length.out = 3),
    row.names = salmon_setup$sample_names
  )
  
  # Create tx2gene mapping file
  tx2gene_file <- file.path(tmpdir, "tx2gene.tsv")
  data <- readr::read_tsv(salmon_setup$file_paths[1], show_col_types = FALSE)
  tx2gene_df <- data.frame(
    Transcript = data$Name,
    Gene = paste0("ENSG", sprintf("%011d", seq_along(data$Name)))
  )
  readr::write_tsv(tx2gene_df, tx2gene_file)
  
  # Build analysis using salmon_dir
  analysis <- build_analysis_s4(
    salmon_dir = salmon_setup$salmon_dir,
    tx2gene = tx2gene_file,
    metadata = metadata,
    skip = TRUE,
    verbose = FALSE
  )
  
  expect_is(analysis, "TSENATAnalysis")
  se <- getSE(analysis)
  expect_equal(ncol(se), 3)  # 3 samples
  expect_gt(nrow(se), 0)  # Some genes/transcripts
})

test_that("build_analysis_s4 with salmon_dir stores TPM and effective_length", {
  skip_if_not_installed("TSENAT")
  
  tmpdir <- tempdir()
  on.exit(unlink(file.path(tmpdir, "salmon_meta"), recursive = TRUE))
  
  salmon_setup <- create_mock_salmon_dir(tmpdir, sample_names = c("S1", "S2"))
  
  metadata <- data.frame(
    sample = salmon_setup$sample_names,
    condition = rep("control", 2),
    row.names = salmon_setup$sample_names
  )
  
  tx2gene_file <- file.path(tmpdir, "tx2gene.tsv")
  data <- readr::read_tsv(salmon_setup$file_paths[1], show_col_types = FALSE)
  tx2gene_df <- data.frame(
    Transcript = data$Name,
    Gene = paste0("ENSG", sprintf("%011d", seq_along(data$Name)))
  )
  readr::write_tsv(tx2gene_df, tx2gene_file)
  
  analysis <- build_analysis_s4(
    salmon_dir = salmon_setup$salmon_dir,
    tx2gene = tx2gene_file,
    metadata = metadata,
    skip = TRUE,
    verbose = FALSE
  )
  
  se <- getSE(analysis)
  
  # Check metadata storage
  expect_true(!is.null(metadata(se)$tpm))
  expect_true(!is.null(metadata(se)$effective_length))
  expect_equal(dim(metadata(se)$tpm), c(100, 2))
  expect_equal(length(metadata(se)$effective_length), 100)
})

test_that("build_analysis_s4 with salmon direct parameters works", {
  skip_if_not_installed("TSENAT")
  
  tmpdir <- tempdir()
  on.exit(unlink(file.path(tmpdir, "direct_salmon*"), recursive = TRUE))
  
  salmon_setup <- create_mock_salmon_dir(tmpdir, sample_names = c("S1", "S2"))
  
  # Read Salmon data directly
  salmon_data <- .read_salmon_samples(
    file_paths = salmon_setup$file_paths,
    sample_names = salmon_setup$sample_names,
    include_tpm = TRUE,
    include_eff_length = TRUE,
    verbose = FALSE
  )
  
  metadata <- data.frame(
    sample = salmon_setup$sample_names,
    condition = rep("control", 2),
    row.names = salmon_setup$sample_names
  )
  
  tx2gene_file <- file.path(tmpdir, "tx2gene.tsv")
  tx2gene_df <- data.frame(
    Transcript = salmon_data$transcript_ids,
    Gene = paste0("ENSG", sprintf("%011d", seq_along(salmon_data$transcript_ids)))
  )
  readr::write_tsv(tx2gene_df, tx2gene_file)
  
  # Build directly with salmon data
  # NOTE: effective_length should be a vector of length n_transcripts, not a matrix
  effective_length_vec <- if (is.matrix(salmon_data$effective_length)) {
    salmon_data$effective_length[, 1]  # Extract first sample column
  } else {
    salmon_data$effective_length
  }
  
  analysis <- build_analysis_s4(
    readcounts = salmon_data$counts,
    tx2gene = tx2gene_file,
    metadata = metadata,
    tpm = salmon_data$tpm,
    effective_length = effective_length_vec,
    skip = TRUE,
    verbose = FALSE
  )
  
  expect_is(analysis, "TSENATAnalysis")
  se <- getSE(analysis)
  expect_equal(dim(se), c(100, 2))
})

test_that("build_analysis_s4 handles skip=TRUE for unmapped transcripts", {
  skip_if_not_installed("TSENAT")
  
  tmpdir <- tempdir()
  on.exit(unlink(file.path(tmpdir, "skip_test"), recursive = TRUE))
  
  # Create Salmon data with 100 transcripts
  salmon_setup <- create_mock_salmon_dir(tmpdir, sample_names = c("S1", "S2"))
  salmon_data <- .read_salmon_samples(
    file_paths = salmon_setup$file_paths,
    sample_names = salmon_setup$sample_names,
    include_tpm = TRUE,
    include_eff_length = TRUE,
    verbose = FALSE
  )
  
  # Create tx2gene with only 80 transcripts (skip 20)
  tx2gene_file <- file.path(tmpdir, "tx2gene_partial.tsv")
  tx2gene_df <- data.frame(
    Transcript = salmon_data$transcript_ids[1:80],  # Only first 80
    Gene = paste0("ENSG", sprintf("%011d", 1:80))
  )
  readr::write_tsv(tx2gene_df, tx2gene_file)
  
  metadata <- data.frame(
    sample = c("S1", "S2"),
    condition = rep("control", 2),
    row.names = c("S1", "S2")
  )
  
  # With skip=TRUE, should use only mapped transcripts
  analysis <- build_analysis_s4(
    readcounts = salmon_data$counts,
    tx2gene = tx2gene_file,
    metadata = metadata,
    tpm = salmon_data$tpm,
    effective_length = salmon_data$effective_length,
    skip = TRUE,
    verbose = FALSE
  )
  
  se <- getSE(analysis)
  expect_equal(nrow(se), 80)  # Only mapped transcripts retained
})

test_that("build_analysis_s4 with verbose=TRUE shows progress", {
  skip_if_not_installed("TSENAT")
  
  tmpdir <- tempdir()
  on.exit(unlink(file.path(tmpdir, "verbose_test*"), recursive = TRUE))
  
  salmon_setup <- create_mock_salmon_dir(tmpdir, sample_names = c("S1", "S2"))
  
  metadata <- data.frame(
    sample = salmon_setup$sample_names,
    condition = rep("control", 2),
    row.names = salmon_setup$sample_names
  )
  
  tx2gene_file <- file.path(tmpdir, "tx2gene.tsv")
  data <- readr::read_tsv(salmon_setup$file_paths[1], show_col_types = FALSE)
  tx2gene_df <- data.frame(
    Transcript = data$Name,
    Gene = paste0("ENSG", sprintf("%011d", seq_along(data$Name)))
  )
  readr::write_tsv(tx2gene_df, tx2gene_file)
  
  # Test that verbose=TRUE produces messages
  # Use expect_message to verify messages are produced during execution
  expect_message(
    analysis <- build_analysis_s4(
      salmon_dir = salmon_setup$salmon_dir,
      tx2gene = tx2gene_file,
      metadata = metadata,
      skip = TRUE,
      verbose = TRUE
    ),
    # Should produce messages about Salmon detection or data reading
    pattern = "Salmon|detected|read|sample"
  )
  
  expect_is(analysis, "TSENATAnalysis")
  expect_gt(nrow(getSE(analysis)), 0)
})

# ===========================================================================
# TESTS: Integration - Full Salmon workflow
# ===========================================================================

test_that("Full Salmon workflow: detect → read → build_analysis creates correct object", {
  skip_if_not_installed("TSENAT")
  
  tmpdir <- tempdir()
  on.exit(unlink(file.path(tmpdir, "full_workflow*"), recursive = TRUE))
  
  # Create mock Salmon data
  salmon_setup <- create_mock_salmon_dir(
    tmpdir,
    sample_names = c("ctrl_1", "ctrl_2", "treat_1", "treat_2")
  )
  
  # Step 1: Detect samples
  salmon_info <- .detect_salmon_samples(salmon_setup$salmon_dir)
  expect_equal(salmon_info$count, 4)
  
  # Step 2: Read Salmon data
  salmon_data <- .read_salmon_samples(
    file_paths = salmon_info$file_paths,
    sample_names = salmon_info$sample_names,
    verbose = FALSE
  )
  expect_equal(ncol(salmon_data$counts), 4)
  n_transcripts <- nrow(salmon_data$counts)
  
  # Step 3: Create metadata matching Salmon samples
  metadata <- data.frame(
    sample = salmon_info$sample_names,
    condition = rep(c("control", "treatment"), each = 2),
    pair_col = rep(c(1, 2), 2),
    row.names = salmon_info$sample_names
  )
  
  # Step 4: Create tx2gene for building analysis
  tx2gene_file <- file.path(tmpdir, "tx2gene.tsv")
  tx2gene_df <- data.frame(
    Transcript = salmon_data$transcript_ids,
    Gene = paste0("ENSG", sprintf("%011d", seq_along(salmon_data$transcript_ids)))
  )
  readr::write_tsv(tx2gene_df, tx2gene_file)
  
  # Extract effective_length as vector (use first sample if matrix)
  effective_length_vec <- if (is.matrix(salmon_data$effective_length)) {
    salmon_data$effective_length[, 1]
  } else {
    salmon_data$effective_length
  }
  
  # Step 5: Build analysis
  analysis <- build_analysis_s4(
    readcounts = salmon_data$counts,
    tx2gene = tx2gene_file,
    metadata = metadata,
    tpm = salmon_data$tpm,
    effective_length = effective_length_vec,
    skip = TRUE,
    verbose = FALSE
  )
  
  # Validate object structure and content
  expect_is(analysis, "TSENATAnalysis")
  
  se <- getSE(analysis)
  
  # Validate dimensions
  expect_equal(ncol(se), 4)  # 4 samples
  expect_equal(nrow(se), n_transcripts)  # Same number of transcripts
  
  # Validate assays are present and correct
  expect_true("counts" %in% assayNames(se))
  
  # Validate assay content has correct dimensions
  expect_equal(nrow(assay(se, "counts")), n_transcripts)
  expect_equal(ncol(assay(se, "counts")), 4)
  
  # Validate counts are non-negative and numeric
  expect_true(all(assay(se, "counts") >= 0, na.rm = TRUE))
  
  # Validate TPM stored in metadata (Salmon data stored there, not as assay)
  obj_meta <- S4Vectors::metadata(se)
  expect_true("tpm" %in% names(obj_meta) || "tpm" %in% names(obj_meta))
  if ("tpm" %in% names(obj_meta)) {
    tpm_mat <- obj_meta$tpm
    expect_equal(nrow(tpm_mat), n_transcripts)
    expect_equal(ncol(tpm_mat), 4)
    # TPM values should be numeric and non-negative
    expect_true(all(tpm_mat >= 0, na.rm = TRUE))
  }
  
  # Validate effective_length stored in metadata
  expect_true("effective_length" %in% names(obj_meta) || "effective_length" %in% names(obj_meta))
  if ("effective_length" %in% names(obj_meta)) {
    eff_length <- obj_meta$effective_length
    expect_equal(length(eff_length), n_transcripts)
    expect_true(all(eff_length > 0, na.rm = TRUE))
  }
  
  # Validate rowData (gene mapping)
  rd <- rowData(se)
  expect_true(!is.null(rd$gene_id))
  expect_equal(length(rd$gene_id), n_transcripts)
  expect_true(all(grepl("^ENSG", rd$gene_id)))  # Gene IDs should have ENSG prefix
  
  # Validate colData (sample metadata)
  cd <- colData(se)
  expect_equal(nrow(cd), 4)
  expect_equal(cd$condition, c("control", "control", "treatment", "treatment"))
  
  # Validate metadata was preserved
  obj_metadata <- S4Vectors::metadata(analysis)
  expect_true(!is.null(obj_metadata))
  
  # Test successfully validates:
  # ✓ Salmon sample detection works
  # ✓ Salmon file reading captures all data columns
  # ✓ Building analysis creates valid S4 object with correct structure
  # ✓ All assays populated with proper dimensions
  # ✓ rowData contains correct gene mapping
  # ✓ colData contains correct sample metadata
})
