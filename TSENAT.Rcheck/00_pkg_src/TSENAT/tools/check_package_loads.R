pkgs <- c("MatrixGenerics","matrixStats","GenomicRanges","stats4","BiocGenerics","generics","S4Vectors","IRanges","Seqinfo","Biobase","SummarizedExperiment")
for (p in pkgs) {
  cat("Loading package:", p, "\n")
  tryCatch({
    library(p, character.only = TRUE)
    cat("Loaded", p, "OK\n")
  }, error = function(e) {
    cat("ERROR loading", p, "\n")
    message(e)
    quit(status = 1)
  })
}
cat("All packages loaded\n")
