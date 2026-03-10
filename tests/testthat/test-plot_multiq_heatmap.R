context("plot_multiq_delta_influence_heatmaps")

test_that("plot_multiq_delta_influence_heatmaps requires valid multi-q results", {
  # Test with wrong class
  expect_error(
    plot_multiq_delta_influence_heatmaps(list(q_0_50 = NULL), n_genes = 2),
    "must be a multi-q result"
  )
})

test_that("plot_multiq_delta_influence_heatmaps rejects results without q-values", {
  # Create a mock object with correct class but no q-value results
  mock_result <- list(
    gene_ids = c("ENSG1", "ENSG2"),
    gene_name_map = c("GENE1", "GENE2"),
    summary_table = data.frame()
  )
  class(mock_result) <- c("tsenat_isoform_switching_multiq", "list")
  
  expect_error(
    plot_multiq_delta_influence_heatmaps(mock_result, n_genes = 2),
    "No multi-q results found"
  )
})

test_that("plot_multiq_delta_influence_heatmaps works with valid multi-q results", {
  # Load test data - salmon_dataset from TSENAT
  data("readcounts", package = "TSENAT", envir = environment())
  
  # Build minimal SE
  se <- SummarizedExperiment::SummarizedExperiment(
    assays = list(counts = as.matrix(salmon_dataset[1:50, 1:10])),
    rowData = data.frame(
      gene_id = rep(paste0("GENE", 1:10), 5),
      transcript_id = paste0("TRANS", 1:50),
      gene_name = rep(paste0("Gene", 1:10), 5),
      stringsAsFactors = FALSE
    ),
    colData = data.frame(
      sample_id = colnames(salmon_dataset)[1:10],
      sample_type = rep(c("A", "B"), 5),
      paired_samples = rep(1:5, 2),
      stringsAsFactors = FALSE
    )
  )
  
  # Run minimal multi-q jackknife (just 2 q-values)
  multi_q_results <- jackknife_isoform_switching(
    se = se,
    condition_col = "sample_type",
    pair_col = "paired_samples",
    gene_col = "gene_id",
    isoform_col = "transcript_id",
    q = c(0.5, 1.0),
    norm = TRUE,
    n_bootstrap = 10,
    print_results = FALSE
  )
  
  # Test that plotting works
  heatmap_file <- plot_multiq_delta_influence_heatmaps(
    switching_results = multi_q_results,
    n_genes = 2
  )
  
  expect_true(is.character(heatmap_file))
  expect_true(file.exists(heatmap_file))
  expect_true(grepl("\\.png$", heatmap_file))
  
  # Clean up
  tryCatch(file.remove(heatmap_file), silent = TRUE)
})

test_that("plot_multiq_delta_influence_heatmaps respects n_genes parameter", {
  # Load test data
  data("readcounts", package = "TSENAT", envir = environment())
  
  # Build minimal SE with more genes
  se <- SummarizedExperiment::SummarizedExperiment(
    assays = list(counts = as.matrix(salmon_dataset[1:100, 1:10])),
    rowData = data.frame(
      gene_id = rep(paste0("GENE", 1:20), 5),
      transcript_id = paste0("TRANS", 1:100),
      gene_name = rep(paste0("Gene", 1:20), 5),
      stringsAsFactors = FALSE
    ),
    colData = data.frame(
      sample_id = colnames(salmon_dataset)[1:10],
      sample_type = rep(c("A", "B"), 5),
      paired_samples = rep(1:5, 2),
      stringsAsFactors = FALSE
    )
  )
  
  # Run multi-q jackknife
  multi_q_results <- jackknife_isoform_switching(
    se = se,
    condition_col = "sample_type",
    pair_col = "paired_samples",
    gene_col = "gene_id",
    isoform_col = "transcript_id",
    q = c(0.5, 1.0),
    norm = TRUE,
    n_bootstrap = 10,
    print_results = FALSE
  )
  
  # Test with different n_genes values
  for (n in c(1, 2, 5)) {
    heatmap_file <- plot_multiq_delta_influence_heatmaps(
      switching_results = multi_q_results,
      n_genes = n
    )
    
    expect_true(file.exists(heatmap_file))
    
    # Clean up
    tryCatch(file.remove(heatmap_file), silent = TRUE)
  }
})

test_that("plot_multiq_delta_influence_heatmaps handles n_genes > available genes", {
  # Load test data
  data("readcounts", package = "TSENAT", envir = environment())
  
  # Build minimal SE
  se <- SummarizedExperiment::SummarizedExperiment(
    assays = list(counts = as.matrix(salmon_dataset[1:50, 1:10])),
    rowData = data.frame(
      gene_id = rep(paste0("GENE", 1:10), 5),
      transcript_id = paste0("TRANS", 1:50),
      gene_name = rep(paste0("Gene", 1:10), 5),
      stringsAsFactors = FALSE
    ),
    colData = data.frame(
      sample_id = colnames(salmon_dataset)[1:10],
      sample_type = rep(c("A", "B"), 5),
      paired_samples = rep(1:5, 2),
      stringsAsFactors = FALSE
    )
  )
  
  # Run multi-q jackknife
  multi_q_results <- jackknife_isoform_switching(
    se = se,
    condition_col = "sample_type",
    pair_col = "paired_samples",
    gene_col = "gene_id",
    isoform_col = "transcript_id",
    q = c(0.5, 1.0),
    norm = TRUE,
    n_bootstrap = 10,
    print_results = FALSE
  )
  
  # Request more genes than available - should gracefully use available genes
  heatmap_file <- plot_multiq_delta_influence_heatmaps(
    switching_results = multi_q_results,
    n_genes = 1000  # More than available
  )
  
  expect_true(file.exists(heatmap_file))
  
  # Clean up
  tryCatch(file.remove(heatmap_file), silent = TRUE)
})

test_that("plot_multiq_delta_influence_heatmaps handles q-values correctly", {
  # Load test data
  data("readcounts", package = "TSENAT", envir = environment())
  
  # Build minimal SE
  se <- SummarizedExperiment::SummarizedExperiment(
    assays = list(counts = as.matrix(salmon_dataset[1:50, 1:10])),
    rowData = data.frame(
      gene_id = rep(paste0("GENE", 1:10), 5),
      transcript_id = paste0("TRANS", 1:50),
      gene_name = rep(paste0("Gene", 1:10), 5),
      stringsAsFactors = FALSE
    ),
    colData = data.frame(
      sample_id = colnames(salmon_dataset)[1:10],
      sample_type = rep(c("A", "B"), 5),
      paired_samples = rep(1:5, 2),
      stringsAsFactors = FALSE
    )
  )
  
  # Run multi-q jackknife with multiple q-values
  multi_q_results <- jackknife_isoform_switching(
    se = se,
    condition_col = "sample_type",
    pair_col = "paired_samples",
    gene_col = "gene_id",
    isoform_col = "transcript_id",
    q = c(0.5, 1.0, 1.5),
    norm = TRUE,
    n_bootstrap = 10,
    print_results = FALSE
  )
  
  # Verify correct number of q-values
  q_keys <- names(multi_q_results)[grepl("^q_", names(multi_q_results))]
  expect_equal(length(q_keys), 3)
  
  # Test plotting works with multiple q-values
  heatmap_file <- plot_multiq_delta_influence_heatmaps(
    switching_results = multi_q_results,
    n_genes = 2
  )
  
  expect_true(file.exists(heatmap_file))
  
  # Clean up
  tryCatch(file.remove(heatmap_file), silent = TRUE)
})

test_that("plot_multiq_delta_influence_heatmaps creates valid PNG file", {
  # Load test data
  data("readcounts", package = "TSENAT", envir = environment())
  
  # Build minimal SE
  se <- SummarizedExperiment::SummarizedExperiment(
    assays = list(counts = as.matrix(salmon_dataset[1:50, 1:10])),
    rowData = data.frame(
      gene_id = rep(paste0("GENE", 1:10), 5),
      transcript_id = paste0("TRANS", 1:50),
      gene_name = rep(paste0("Gene", 1:10), 5),
      stringsAsFactors = FALSE
    ),
    colData = data.frame(
      sample_id = colnames(salmon_dataset)[1:10],
      sample_type = rep(c("A", "B"), 5),
      paired_samples = rep(1:5, 2),
      stringsAsFactors = FALSE
    )
  )
  
  # Run multi-q jackknife
  multi_q_results <- jackknife_isoform_switching(
    se = se,
    condition_col = "sample_type",
    pair_col = "paired_samples",
    gene_col = "gene_id",
    isoform_col = "transcript_id",
    q = c(0.5, 1.0),
    norm = TRUE,
    n_bootstrap = 10,
    print_results = FALSE
  )
  
  heatmap_file <- plot_multiq_delta_influence_heatmaps(
    switching_results = multi_q_results,
    n_genes = 2
  )
  
  # Check file size (PNG should be non-trivial size)
  file_size <- file.size(heatmap_file)
  expect_true(file_size > 1000)  # At least 1KB
  
  # Clean up
  tryCatch(file.remove(heatmap_file), silent = TRUE)
})
