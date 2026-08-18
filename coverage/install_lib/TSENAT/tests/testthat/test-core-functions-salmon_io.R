# Comprehensive tests for Salmon I/O functions and build_analysis integration
library(TSENAT)
library(testthat)

# Comprehensive tests for Salmon I/O functions and build_analysis integration
# Tests for: .detect_salmon_samples(), .validate_salmon_files(), .read_salmon_samples()
# And: build_analysis() with salmon_dir and salmon data parameters

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

test_that(".detect_salmon_samples stops on duplicate sample names", {
  tmpdir <- tempdir()
  on.exit(unlink(file.path(tmpdir, "salmon_dup"), recursive = TRUE))
  
  salmon_dir <- file.path(tmpdir, "salmon_dup")
  
  # Create nested structure with duplicate folder names
  dir.create(file.path(salmon_dir, "batch1", "sample1"), recursive = TRUE, showWarnings = FALSE)
  dir.create(file.path(salmon_dir, "batch2", "sample1"), recursive = TRUE, showWarnings = FALSE)
  
  create_mock_salmon_file(file.path(salmon_dir, "batch1", "sample1", "quant.sf"))
  create_mock_salmon_file(file.path(salmon_dir, "batch2", "sample1", "quant.sf"))
  
  # AUDIT S6: duplicate sample names are a hard error (previously a warning
  # that let duplicate matrix columns corrupt downstream analysis)
  expect_error(
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

test_that(".validate_salmon_files stops on transcript ID mismatches", {
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
  
  # AUDIT S5: mismatch is a hard error (previously a warning followed by
  # positional matrix fill that silently corrupted downstream entropy)
  expect_error(
    .validate_salmon_files(
      c(file.path(dir1, "quant.sf"), file.path(dir2, "quant.sf")),
      verbose = FALSE
    ),
    "Transcript ID mismatch"
  )
})

test_that(".read_salmon_samples reorders shuffled transcripts by ID, not position", {
  tmpdir <- tempdir()
  on.exit(unlink(file.path(tmpdir, "reorder_test"), recursive = TRUE))

  dir1 <- file.path(tmpdir, "reorder_test", "s1")
  dir2 <- file.path(tmpdir, "reorder_test", "s2")
  dir.create(dir1, recursive = TRUE, showWarnings = FALSE)
  dir.create(dir2, recursive = TRUE, showWarnings = FALSE)

  n_tx <- 20
  base <- data.frame(
    Name = paste0("TX", seq_len(n_tx)),
    Length = seq_len(n_tx) * 10,
    EffectiveLength = seq_len(n_tx) * 9,
    TPM = rep(1, n_tx),
    NumReads = seq_len(n_tx) * 100
  )
  readr::write_tsv(base, file.path(dir1, "quant.sf"))

  # Same transcripts, DIFFERENT order, different counts per transcript
  shuffled <- base[sample(n_tx), , drop = FALSE]
  shuffled$NumReads <- rev(shuffled$NumReads)
  shuffled$TPM <- rev(shuffled$TPM)
  shuffled$EffectiveLength <- rev(shuffled$EffectiveLength)
  readr::write_tsv(shuffled, file.path(dir2, "quant.sf"))

  res <- .read_salmon_samples(c(file.path(dir1, "quant.sf"), file.path(dir2, "quant.sf")),
    sample_names = c("s1", "s2"), verbose = FALSE)

  expect_identical(unname(res$counts[, 1]), base$NumReads)
  # Counts for s2 must align with the s1 transcript order
  expect_identical(unname(res$counts[, 2]), shuffled$NumReads[match(base$Name, shuffled$Name)])
  expect_identical(res$transcript_ids, base$Name)
})

test_that(".read_salmon_samples stops when transcripts are missing in a file", {
  tmpdir <- tempdir()
  on.exit(unlink(file.path(tmpdir, "missing_tx_test"), recursive = TRUE))

  dir1 <- file.path(tmpdir, "missing_tx_test", "s1")
  dir2 <- file.path(tmpdir, "missing_tx_test", "s2")
  dir.create(dir1, recursive = TRUE, showWarnings = FALSE)
  dir.create(dir2, recursive = TRUE, showWarnings = FALSE)

  base <- data.frame(
    Name = paste0("TX", 1:10), Length = 1:10, EffectiveLength = 1:10,
    TPM = 1, NumReads = 10
  )
  subset_file <- base[1:5, ]
  readr::write_tsv(base, file.path(dir1, "quant.sf"))
  readr::write_tsv(subset_file, file.path(dir2, "quant.sf"))

  expect_error(
    .read_salmon_samples(c(file.path(dir1, "quant.sf"), file.path(dir2, "quant.sf")),
      sample_names = c("s1", "s2"), verbose = FALSE),
    "Missing in file 2"
  )
})

test_that(".validate_salmon_files rejects negative/non-finite quantification values", {
  tmpdir <- tempdir()
  on.exit(unlink(file.path(tmpdir, "numeric_integrity"), recursive = TRUE))

  bad_dir <- file.path(tmpdir, "numeric_integrity", "sample1")
  dir.create(bad_dir, recursive = TRUE, showWarnings = FALSE)
  bad_file <- file.path(bad_dir, "quant.sf")

  bad_content <- "Name\tLength\tEffectiveLength\tTPM\tNumReads\n"
  bad_content <- paste0(bad_content, "ENST00000000001\t1000\t900\t5.5\t-3.2\n")
  bad_content <- paste0(bad_content, "ENST00000000002\t1001\t901\t3.2\t10.0\n")
  writeLines(bad_content, bad_file)

  expect_error(
    .validate_salmon_files(bad_file, verbose = FALSE),
    "Negative values"
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
# TESTS: build_analysis() with Salmon integration
# ===========================================================================

test_that("build_analysis works with salmon_dir parameter", {
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
  config <- TSENAT_config(
    sample_col = "sample",
    condition_col = "condition"
  )
  
  analysis <- build_analysis(
    config = config,
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

test_that("build_analysis with salmon_dir stores TPM and effective_length", {
  skip_if_not_installed("TSENAT")
  
  tmpdir <- tempdir()
  on.exit(unlink(file.path(tmpdir, "salmon_meta"), recursive = TRUE))
  
  salmon_setup <- create_mock_salmon_dir(tmpdir, sample_names = c("S1", "S2"))
  
  metadata <- data.frame(
    sample = salmon_setup$sample_names,
    condition = c("control", "treatment"),
    row.names = salmon_setup$sample_names
  )
  
  tx2gene_file <- file.path(tmpdir, "tx2gene.tsv")
  data <- readr::read_tsv(salmon_setup$file_paths[1], show_col_types = FALSE)
  tx2gene_df <- data.frame(
    Transcript = data$Name,
    Gene = paste0("ENSG", sprintf("%011d", seq_along(data$Name)))
  )
  readr::write_tsv(tx2gene_df, tx2gene_file)
  
  config <- TSENAT_config(
    sample_col = "sample",
    condition_col = "condition"
  )
  
  analysis <- build_analysis(
    config = config,
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

test_that("build_analysis with salmon direct parameters works", {
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
    condition = c("control", "treatment"),
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
  
  config <- TSENAT_config(
    sample_col = "sample",
    condition_col = "condition"
  )
  
  analysis <- build_analysis(
    config = config,
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

test_that("build_analysis handles skip=TRUE for unmapped transcripts", {
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
    condition = c("control", "treatment"),
    row.names = c("S1", "S2")
  )
  
  # With skip=TRUE, should use only mapped transcripts
  config <- TSENAT_config(
    sample_col = "sample",
    condition_col = "condition"
  )
  
  analysis <- build_analysis(
    config = config,
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

test_that("build_analysis with verbose=TRUE shows progress", {
  skip_if_not_installed("TSENAT")
  
  tmpdir <- tempdir()
  on.exit(unlink(file.path(tmpdir, "verbose_test*"), recursive = TRUE))
  
  salmon_setup <- create_mock_salmon_dir(tmpdir, sample_names = c("S1", "S2"))
  
  metadata <- data.frame(
    sample = salmon_setup$sample_names,
    condition = c("control", "treatment"),
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
  config <- TSENAT_config(
    sample_col = "sample",
    condition_col = "condition"
  )
  
  expect_message(
    analysis <- build_analysis(
      config = config,
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
  config <- TSENAT_config(
    sample_col = "sample",
    condition_col = "condition"
  )
  
  analysis <- build_analysis(
    config = config,
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

# ===========================================================================
# TESTS: Salmon Data with Decimal Values (Bug Fix Verification)
# ===========================================================================

#' Helper: Create Salmon data with realistic decimal Length values
#' Mimics actual Salmon output with mixed integer/decimal Length values
#' @noRd
create_realistic_salmon_file <- function(filepath, seed = 42) {
  set.seed(seed)
  
  # Create data with realistic mix of decimal and integer Length values
  # (matching actual Salmon quant.sf format)
  data <- data.frame(
    Name = c(
      "ENST00000162391.8",    # Will get 6289.435
      "ENST00000187762.7",    # Will get 2516.0285
      "ENST00000221130.11",   # Will get 2465.158
      "ENST00000221818.5",    # Will get 1535.5
      "ENST00000261252.4"     # Will get 2351 (integer)
    ),
    Length = c(6289.435, 2516.0285, 2465.158, 1535.5, 2351.0),
    EffectiveLength = c(6289.435, 2516.0285, 2465.158, 1535.5, 2351.0),
    TPM = c(6.062312, 0.736927, 19.426005, 0.102499, 15.234567),
    NumReads = c(757.16, 36.57, 935.366, 3.0, 452.0)
  )
  
  # Add more realistic transcripts to fill out 100 rows
  for (i in 6:100) {
    data <- rbind(data, data.frame(
      Name = paste0("ENST", sprintf("%011d", i)),
      Length = runif(1, 500, 5000) + runif(1, 0, 1),  # Mix of decimal and near-integer
      EffectiveLength = runif(1, 400, 4800) + runif(1, 0, 1),
      TPM = runif(1, 0.1, 50),
      NumReads = rpois(1, 75)
    ))
  }
  
  dir.create(dirname(filepath), showWarnings = FALSE, recursive = TRUE)
  readr::write_tsv(data, filepath)
  invisible(data)
}

test_that("Salmon data with decimals: detect samples", {
  tmpdir <- tempdir()
  on.exit(unlink(file.path(tmpdir, "decimal_test_*"), recursive = TRUE))
  
  salmon_setup <- create_mock_salmon_dir(
    tmpdir,
    sample_names = c("sample1", "sample2", "sample3")
  )
  
  result <- .detect_salmon_samples(salmon_setup$salmon_dir, recursive = TRUE)
  
  expect_is(result, "list")
  expect_named(result, c("sample_names", "file_paths", "count"))
  expect_equal(result$count, 3)
  expect_equal(length(result$sample_names), 3)
  expect_equal(length(result$file_paths), 3)
  
  # Verify file paths are correct
  expect_true(all(file.exists(result$file_paths)))
  expect_true(all(grepl("quant.sf$", result$file_paths)))
})

test_that("Salmon data with decimals: validate files", {
  tmpdir <- tempdir()
  on.exit(unlink(file.path(tmpdir, "validate_dec*"), recursive = TRUE))
  
  salmon_setup <- create_mock_salmon_dir(tmpdir, sample_names = c("S1", "S2"))
  
  # All files should validate successfully
  expect_true(.validate_salmon_files(salmon_setup$file_paths, verbose = FALSE))
})

test_that("Salmon data with decimals: Length column preserves decimal values", {
  tmpdir <- tempdir()
  on.exit(unlink(file.path(tmpdir, "decimal_length*"), recursive = TRUE))
  
  # Create file with realistic decimal values
  test_dir <- file.path(tmpdir, "decimal_length_test", "sample1")
  dir.create(test_dir, recursive = TRUE, showWarnings = FALSE)
  test_file <- file.path(test_dir, "quant.sf")
  
  create_realistic_salmon_file(test_file)
  
  # Read with col_double (correct, after fix)
  data_double <- readr::read_tsv(test_file, 
    col_types = readr::cols(
      Name = readr::col_character(),
      Length = readr::col_double(),
      EffectiveLength = readr::col_double(),
      TPM = readr::col_double(),
      NumReads = readr::col_double()
    ),
    show_col_types = FALSE)
  
  # Verify decimal values are preserved
  enst_6289 <- data_double$Length[data_double$Name == "ENST00000162391.8"]
  expect_equal(enst_6289, 6289.435)  # Exact decimal value preserved
  
  enst_2516 <- data_double$Length[data_double$Name == "ENST00000187762.7"]
  expect_equal(enst_2516, 2516.0285)  # Exact decimal value preserved
  
  enst_2465 <- data_double$Length[data_double$Name == "ENST00000221130.11"]
  expect_equal(enst_2465, 2465.158)  # Exact decimal value preserved
  
  # Verify that col_integer would have truncated these
  expect_false(enst_6289 == 6289)  # Would be truncated to 6289
  expect_false(enst_2516 == 2516)  # Would be truncated to 2516
  expect_false(enst_2465 == 2465)  # Would be truncated to 2465
})

test_that("Salmon data with decimals: read produces correct structures", {
  tmpdir <- tempdir()
  on.exit(unlink(file.path(tmpdir, "read_decimal*"), recursive = TRUE))
  
  salmon_setup <- create_mock_salmon_dir(tmpdir, sample_names = c("S1", "S2", "S3"))
  
  # Read all samples
  salmon_data <- .read_salmon_samples(
    file_paths = salmon_setup$file_paths,
    sample_names = salmon_setup$sample_names,
    include_tpm = TRUE,
    include_eff_length = TRUE,
    verbose = FALSE
  )
  
  # Verify structure
  expect_is(salmon_data, "list")
  expect_true(all(c("counts", "transcript_ids", "tpm", "effective_length") %in% names(salmon_data)))
  
  # Verify dimensions
  expect_equal(ncol(salmon_data$counts), 3)
  expect_equal(ncol(salmon_data$tpm), 3)
  expect_equal(nrow(salmon_data$counts), 100)
  expect_equal(nrow(salmon_data$tpm), 100)
  expect_equal(nrow(salmon_data$effective_length), 100)  # Returns as matrix
  expect_equal(ncol(salmon_data$effective_length), 3)
  
  # Verify all values are numeric and realistic
  expect_true(all(salmon_data$counts >= 0, na.rm = TRUE))
  expect_true(all(salmon_data$tpm >= 0, na.rm = TRUE))
  expect_true(all(salmon_data$effective_length > 0, na.rm = TRUE))
})

test_that("Salmon data with decimals: all samples have matching transcript IDs", {
  tmpdir <- tempdir()
  on.exit(unlink(file.path(tmpdir, "matching_trans*"), recursive = TRUE))
  
  salmon_setup <- create_mock_salmon_dir(
    tmpdir,
    sample_names = c("ctrl_1", "ctrl_2", "treat_1")
  )
  
  # Read all files individually
  all_transcripts <- lapply(salmon_setup$file_paths, function(file) {
    data <- readr::read_tsv(file, 
      col_types = readr::cols(
        Name = readr::col_character(),
        .default = readr::col_skip()
      ),
      show_col_types = FALSE)
    data$Name
  })
  
  # All files should have identical transcript lists
  for (i in 2:length(all_transcripts)) {
    expect_equal(all_transcripts[[1]], all_transcripts[[i]],
      label = paste("Sample", i, "has matching transcripts"))
  }
})

test_that("Salmon data with decimals: realistic distributions", {
  tmpdir <- tempdir()
  on.exit(unlink(file.path(tmpdir, "distrib_test*"), recursive = TRUE))
  
  salmon_setup <- create_mock_salmon_dir(tmpdir, sample_names = c("S1", "S2"))
  
  salmon_data <- .read_salmon_samples(
    file_paths = salmon_setup$file_paths,
    sample_names = salmon_setup$sample_names,
    include_tpm = TRUE,
    include_eff_length = TRUE,
    verbose = FALSE
  )
  
  # Verify realistic distributions
  tpm_matrix <- salmon_data$tpm
  tpm_first <- tpm_matrix[, 1]
  
  # Expect skewed distribution (typical for RNA-seq)
  # Mock data ensures minimum detection of 0.1, so check for variation
  expect_gt(length(unique(tpm_first)), 50)  # Many different values
  expect_gt(max(tpm_first, na.rm = TRUE), min(tpm_first, na.rm = TRUE))  # Has range
  expect_gt(max(tpm_first, na.rm = TRUE), 0)  # Maximum should be positive
  
  # Counts should be non-negative
  counts_matrix <- salmon_data$counts
  expect_true(all(counts_matrix >= 0, na.rm = TRUE))
  
  # Effective length should be in realistic range (200-10000 bp)
  expect_gt(min(salmon_data$effective_length, na.rm = TRUE), 100)
  expect_lt(max(salmon_data$effective_length, na.rm = TRUE), 20000)
})

test_that("Salmon data with decimals: mixed integer and decimal values handled correctly", {
  tmpdir <- tempdir()
  on.exit(unlink(file.path(tmpdir, "mixed_decimal*"), recursive = TRUE))
  
  test_dir <- file.path(tmpdir, "mixed_decimal_test", "sample1")
  dir.create(test_dir, recursive = TRUE, showWarnings = FALSE)
  test_file <- file.path(test_dir, "quant.sf")
  
  create_realistic_salmon_file(test_file)
  
  # Read full data
  data <- readr::read_tsv(test_file, 
    col_types = readr::cols(
      Name = readr::col_character(),
      Length = readr::col_double(),
      EffectiveLength = readr::col_double(),
      TPM = readr::col_double(),
      NumReads = readr::col_double()
    ),
    show_col_types = FALSE)
  
  # Verify both decimal and integer-like values are handled
  lengths <- data$Length
  
  # Find values with decimals and without
  decimal_values <- lengths[lengths != floor(lengths)]
  integer_values <- lengths[lengths == floor(lengths)]
  
  # Should have both types
  expect_gt(length(decimal_values), 0)  # Some decimal values
  expect_gt(length(integer_values), 0)  # Some integer-like values
  
  # All should be valid numbers
  expect_true(all(is.numeric(lengths)))
  expect_true(all(lengths > 0))  # Lengths should always be positive
})

test_that("Salmon data with decimals: column types correct after reading", {
  tmpdir <- tempdir()
  on.exit(unlink(file.path(tmpdir, "col_types_test*"), recursive = TRUE))
  
  salmon_setup <- create_mock_salmon_dir(tmpdir, sample_names = "S1")
  
  # Read first sample
  data <- readr::read_tsv(salmon_setup$file_paths[1], 
    col_types = readr::cols(
      Name = readr::col_character(),
      Length = readr::col_double(),
      EffectiveLength = readr::col_double(),
      TPM = readr::col_double(),
      NumReads = readr::col_double()
    ),
    show_col_types = FALSE)
  
  # Verify column classes
  expect_is(data$Name, "character")
  expect_is(data$Length, "numeric")  # col_double produces numeric
  expect_is(data$EffectiveLength, "numeric")
  expect_is(data$TPM, "numeric")
  expect_is(data$NumReads, "numeric")
})
