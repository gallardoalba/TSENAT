#!/usr/bin/env Rscript
# Demo: generate tsallis q-curve for specified genes as groups
library(SummarizedExperiment)
if (requireNamespace("pkgload", quietly = TRUE)) pkgload::load_all(".")

genes <- c("GFPT1", "MBD2", "EEF2K", "C1orf86")
qs <- c(0.01, 0.03, 0.05, 0.1)
cols <- unlist(lapply(genes, function(g) paste0(g, "_q=", qs)))

# construct matrix: rows = genes, cols = samples (gene_q)
mat <- matrix(0, nrow = length(genes), ncol = length(cols),
              dimnames = list(genes, cols))
for (i in seq_along(genes)) {
  g <- genes[i]
  # generate deterministic values that vary with q per gene
  vals <- (i) * (1 + qs * 5)
  mat[i, grepl(paste0('^', g, '_q='), colnames(mat))] <- vals
}

# colData: sample_type = base sample name (gene)
base_names <- sub("_q=.*", "", colnames(mat))
coldata <- DataFrame(sample_type = base_names, row.names = colnames(mat))
rowdata <- DataFrame(genes = genes, row.names = genes)
se <- SummarizedExperiment(assays = list(diversity = mat), colData = coldata, rowData = rowdata)

# plot
p <- plot_tsallis_q_curve(se, assay_name = "diversity", sample_type_col = "sample_type")

out_png <- "inst/doc/tsallis_demo.png"
if (!dir.exists(dirname(out_png))) dir.create(dirname(out_png), recursive = TRUE)
if (requireNamespace("ggplot2", quietly = TRUE)) {
  ggplot2::ggsave(out_png, p, width = 7, height = 4)
  message("Saved plot to: ", out_png)
} else {
  message("ggplot2 not available; printing plot object instead")
  print(p)
}
