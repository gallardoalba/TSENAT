context("map_functions more tests")

library(SummarizedExperiment)

test_that("map_tx_to_readcounts accepts file path input", {
  rc <- matrix(1:6, nrow = 3)
  txmap <- data.frame(Transcript = paste0("tx", 1:4), Gene = c("g1","g1","g2","g2"), stringsAsFactors = FALSE)
  tf <- tempfile(fileext = ".tsv")
  write.table(txmap, tf, sep = "\t", row.names = FALSE, quote = FALSE)
  expect_error(map_tx_to_readcounts(rc, tf), "does not match readcounts rows")
  unlink(tf)
})

test_that("map_samples_to_group with mat provided returns single group when no mapping", {
  mat <- matrix(1:4, nrow = 2)
  colnames(mat) <- c("A_q=1","B_q=1")
  res <- map_samples_to_group(c("A","B"), se = NULL, sample_type_col = NULL, mat = mat)
  expect_equal(res, c("Group","Group"))
})

test_that("get_assay_long errors when assay missing", {
  se <- SummarizedExperiment(assays = S4Vectors::SimpleList(dummy = matrix(1, nrow = 1)))
  expect_error(get_assay_long(se, assay_name = "nope"))
})

test_that("prepare_tsallis_long handles no _q suffix and default group", {
  mat <- matrix(1:4, nrow = 2)
  colnames(mat) <- c("S1","S2")
  rownames(mat) <- c("g1","g2")
  se <- SummarizedExperiment(assays = S4Vectors::SimpleList(diversity = mat))
  rowData(se)$genes <- c("G1","G2")
  long <- prepare_tsallis_long(se, assay_name = "diversity", sample_type_col = NULL)
  expect_true(all(long$group == "Group"))
  expect_true(all(is.na(long$q)))
})

test_that("map_tx_to_readcounts supports custom tx_col and data.frame input", {
  rc <- data.frame(c1 = 1:3, c2 = 4:6)
  txmap <- data.frame(ID = paste0("t",1:3), stringsAsFactors = FALSE)
  out <- map_tx_to_readcounts(rc, txmap, tx_col = "ID")
  expect_equal(rownames(out), as.character(txmap$ID))
})

test_that("map_metadata paired validation errors on unpaired bases", {
  mat <- matrix(1:8, nrow = 2)
  colnames(mat) <- c("B1_q=1","B2_q=1","B1_q=2","B2_q=2")
  se <- SummarizedExperiment(assays = S4Vectors::SimpleList(diversity = mat))
  # create coldata missing one condition for B2
  coldata <- data.frame(Sample = c("B1","B1","B2"), Condition = c("N","T","N"), stringsAsFactors = FALSE)
  expect_error(map_metadata(se, coldata, paired = TRUE), "Unpaired samples found in coldata")
})
