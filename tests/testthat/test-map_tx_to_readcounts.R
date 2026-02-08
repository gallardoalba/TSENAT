test_that("map_tx_to_readcounts assigns rownames from tx2gene data.frame", {
  rc <- matrix(1:6, nrow = 3, ncol = 2)
  txmap <- data.frame(Transcript = paste0("tx", 1:3), Gene = c("G1","G1","G2"), stringsAsFactors = FALSE)
  out <- map_tx_to_readcounts(rc, txmap)
  expect_equal(rownames(out), txmap$Transcript)
})

test_that("map_tx_to_readcounts reads tx2gene from file and errors on mismatch", {
  rc <- matrix(1:4, nrow = 2)
  txmap <- data.frame(Transcript = paste0("tx", 1:3), stringsAsFactors = FALSE)
  tmp <- tempfile(fileext = ".tsv")
  utils::write.table(txmap, file = tmp, sep = "\t", row.names = FALSE, quote = FALSE)
  expect_error(map_tx_to_readcounts(rc, tmp))
})
