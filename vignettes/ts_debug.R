# Debug script: reproduce map_coldata_to_se call from vignette
library(pkgload)
pkgload::load_all('.')
set.seed(1)
# small synthetic dataset: 30 transcripts, 3 samples
readcounts <- matrix(rpois(30 * 3, lambda = 10), nrow = 30, ncol = 3)
colnames(readcounts) <- c("S1_N", "S2_T", "S3_N")
genes <- rep(paste0("G", 1:10), length.out = nrow(readcounts))
# compute diversity with single q
q <- 0.1
ts_se <- calculate_diversity(readcounts, genes, q = q, norm = TRUE)
# Provide simple coldata matching base sample names
coldata_df <- data.frame(
  Sample = c("S1_N", "S2_T", "S3_N"),
  Condition = c("Normal", "Tumor", "Normal"),
  stringsAsFactors = FALSE
)
# Call map_coldata_to_se with paired = TRUE
cat('Calling map_coldata_to_se with paired=TRUE...\n')
ts_se2 <- map_coldata_to_se(ts_se, coldata_df, paired = TRUE)
cat('Success. colData names:', paste(colnames(SummarizedExperiment::assay(ts_se2)), collapse=', '), '\n')
