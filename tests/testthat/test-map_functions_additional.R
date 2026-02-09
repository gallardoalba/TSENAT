context("map_functions additional tests")

library(SummarizedExperiment)

test_that("map_tx_to_readcounts assigns rownames when sizes match", {
  rc <- matrix(1:4, nrow = 2)
  txmap <- data.frame(Transcript = c("tx1","tx2"), Gene = c("g1","g2"), stringsAsFactors = FALSE)
  out <- map_tx_to_readcounts(rc, txmap)
  expect_equal(rownames(out), as.character(txmap$Transcript))
})

test_that("map_tx_to_readcounts errors when sizes mismatch and no matching ids", {
  rc <- matrix(1:6, nrow = 3)
  txmap <- data.frame(Transcript = c("a","b"), Gene = c("g1","g2"), stringsAsFactors = FALSE)
  expect_error(map_tx_to_readcounts(rc, txmap), "does not match readcounts rows")
})

test_that("map_tx_to_readcounts matches existing rownames", {
  rc <- matrix(1:4, nrow = 2)
  rownames(rc) <- c("b","a")
  txmap <- data.frame(Transcript = c("a","b","c"), stringsAsFactors = FALSE)
  expect_message(out <- map_tx_to_readcounts(rc, txmap, verbose = TRUE), "Matched and assigned transcript IDs")
  expect_equal(rownames(out), as.character(c("b","a")))
})

test_that("map_samples_to_group uses colData mapping and errors on missing", {
  mat <- matrix(1:4, nrow = 2)
  colnames(mat) <- c("S1_q=1","S2_q=1")
  se <- SummarizedExperiment(assays = SimpleList(diversity = mat))
  SummarizedExperiment::colData(se)$sample_type <- c("A","B")
  # valid mapping
  mapped <- map_samples_to_group(c("S1","S2"), se = se, sample_type_col = "sample_type", mat = NULL)
  expect_equal(mapped, c("A","B"))
  # missing sample should error
  expect_error(map_samples_to_group(c("S1","S3"), se = se, sample_type_col = "sample_type", mat = NULL), "Missing sample_type mapping")
})

test_that("get_assay_long returns long df and respects sample_type_col", {
  mat <- matrix(c(1,2,3,4), nrow = 2)
  rownames(mat) <- c("g1","g2")
  colnames(mat) <- c("S1","S2")
  se <- SummarizedExperiment(assays = SimpleList(diversity = mat))
  rowData(se)$genes <- c("G1","G2")
  SummarizedExperiment::colData(se) <- S4Vectors::DataFrame(sample_type = c("N","T"), row.names = colnames(mat))
  long <- get_assay_long(se, assay_name = "diversity", value_name = "val", sample_type_col = "sample_type")
  expect_true(all(c("Gene","sample","val","sample_type") %in% colnames(long)))
  expect_equal(unique(as.character(long$sample_type)), c("N","T"))
})

test_that("prepare_tsallis_long parses _q= suffixes and maps groups", {
  mat <- matrix(c(1,NA,2,3,4,NA), nrow = 3)
  colnames(mat) <- c("S1_q=1","S1_q=2")
  rownames(mat) <- c("g1","g2","g3")
  se <- SummarizedExperiment(assays = SimpleList(diversity = mat))
  rowData(se)$genes <- c("G1","G2","G3")
  SummarizedExperiment::colData(se) <- S4Vectors::DataFrame(sample_type = c("A","A"), row.names = colnames(mat))
  long <- prepare_tsallis_long(se, assay_name = "diversity", sample_type_col = "sample_type")
  # q should be a factor and group mapping present
  expect_true("q" %in% colnames(long))
  expect_true(is.factor(long$q))
  expect_true(all(long$group == "A"))
  # rows with NA tsallis removed
  expect_false(any(is.na(long$tsallis)))
})

test_that("map_metadata handles NULL/invalid coldata and maps sample types + metadata", {
  mat <- matrix(1:6, nrow = 3)
  colnames(mat) <- c("S1_q=1","S2_q=1")
  se <- SummarizedExperiment(assays = SimpleList(diversity = mat))
  # NULL coldata returns same
  se2 <- map_metadata(se, NULL)
  expect_identical(se2, se)
  # invalid coldata returns same
  bad <- data.frame(A=1:2, B=c("x","y"))
  se3 <- map_metadata(se, bad)
  expect_identical(se3, se)
  # proper mapping: create coldata with Sample and Condition
  coldata <- data.frame(Sample = c("S1","S2"), Condition = c("N","T"), stringsAsFactors = FALSE)
  # attach readcounts and tx2gene into globalenv to test metadata attaching
  readcounts <- matrix(1:6, nrow = 3)
  tx2gene <- data.frame(Transcript = paste0("tx",1:3), Gen = c("g1","g1","g2"), stringsAsFactors = FALSE)
  assign("readcounts", readcounts, envir = globalenv())
  assign("tx2gene", tx2gene, envir = globalenv())
  on.exit({ rm(readcounts, envir = globalenv()); rm(tx2gene, envir = globalenv()) }, add = TRUE)
  se4 <- map_metadata(se, coldata)
  expect_true("sample_type" %in% colnames(SummarizedExperiment::colData(se4)))
  expect_true("sample_base" %in% colnames(SummarizedExperiment::colData(se4)))
  md <- S4Vectors::metadata(se4)
  expect_true(!is.null(md$readcounts))
  expect_true(!is.null(md$tx2gene))
})
