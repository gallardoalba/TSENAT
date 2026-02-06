suppressPackageStartupMessages({
  library(pkgload)
  pkgload::load_all(".")
})

tx2gene_tsv <- system.file("extdata","tx2gene.tsv", package = "TSENAT")
coldata_tsv <- system.file("extdata","coldata.tsv", package = "TSENAT")
data("tcga_brca_luma_dataset", package = "TSENAT", envir = globalenv())

genes <- tcga_brca_luma_dataset[, 1]
readcounts <- tcga_brca_luma_dataset[, -1]

# mapping and metadata
txmap <- utils::read.delim(tx2gene_tsv, header = TRUE, stringsAsFactors = FALSE)
rownames(readcounts) <- as.character(txmap$Transcript)
coldata_df <- read.table(coldata_tsv, header = TRUE, sep = "\t", stringsAsFactors = FALSE)

# filter lowly-expressed transcripts
tokeep <- rowSums(readcounts > 5) > 5
readcounts <- readcounts[tokeep, ]
genes <- genes[tokeep]

# calculate diversity and map coldata
ts_se <- calculate_diversity(readcounts, genes, q = 0.1, norm = TRUE)
ts_se <- map_coldata_to_se(ts_se, coldata_df, paired = TRUE)

# prepare inputs for calculate_difference
sample_base_names <- sub("_q=.*", "", colnames(assay(ts_se)))
samples <- as.character(colData(ts_se)$sample_type)
div_df <- as.data.frame(assay(ts_se))
div_df <- cbind(genes = rowData(ts_se)$genes, div_df)

# run differential test
res <- calculate_difference(div_df, samples,
    control = "Normal",
  method = "median", test = "wilcoxon",
  paired = TRUE,
  pseudocount = 1e-6
)

# diagnostics
cat("NA count in adjusted_p_values:", sum(is.na(res$adjusted_p_values)), "\n")
cat("Rows with NA adjusted_p_values (up to 20):\n")
print(head(res[is.na(res$adjusted_p_values), ], 20))
cat("\nTable of is.na(res$p_value):\n")
print(table(is.na(res$p_value)))
cat("\nRows with NA p_value (up to 20):\n")
print(head(res[is.na(res$p_value), ], 20))

## Recalculate adjusted p-values using the correct raw p-value column
if ("raw_p_values" %in% colnames(res)) {
  res$adjusted_p_values <- NA
  idx <- which(!is.na(res$raw_p_values))
  res$adjusted_p_values[idx] <- p.adjust(res$raw_p_values[idx], method = "BH")
} else if ("p_value" %in% colnames(res)) {
  res$adjusted_p_values <- NA
  idx <- which(!is.na(res$p_value))
  res$adjusted_p_values[idx] <- p.adjust(res$p_value[idx], method = "BH")
} else {
  warning("No recognizable raw p-value column found in res; skipping recalculation")
}
res <- res[order(is.na(res$adjusted_p_values), res$adjusted_p_values), , drop = FALSE]
cat("\nAfter recalculation, head(res):\n")
print(head(res))

if (!dir.exists("tmp")) dir.create("tmp")
saveRDS(res, file = "tmp/res_diagnostic.rds")
cat("Saved diagnostic res to tmp/res_diagnostic.rds\n")
