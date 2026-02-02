# Script to generate a synthetic tcga_brca_luma_dataset.RData for vignette
# Produces inst/extdata/tcga_brca_luma_dataset.RData with object name
# tcga_brca_luma_dataset (data.frame): first column genes, remaining sample
# counts.

dir.create("inst/extdata", recursive = TRUE, showWarnings = FALSE)
set.seed(2026)

n_genes <- 996
n_pairs <- 20 # will create 40 sample columns (20 normal/tumor pairs)

genes <- paste0("G", seq_len(n_genes))

# build sample names: P1_N, P1_T, P2_N, P2_T, ...
sample_names <- unlist(lapply(seq_len(n_pairs), function(i) c(paste0("P", i, "_N"), paste0("P", i, "_T"))))

# generate counts with slight group differences: tumor samples have slightly
# higher lambda
counts_mat <- matrix(nrow = n_genes, ncol = length(sample_names))
for (i in seq_len(n_genes)) {
    # baseline expression per gene
    base <- round(runif(1, 5, 200))
    # normal/tumor variation
    lambdas <- ifelse(grepl("_N$", sample_names), base, pmax(1, round(base * runif(1, 1.05, 1.5))))
    counts_mat[i, ] <- rpois(length(sample_names), lambda = lambdas)
}
colnames(counts_mat) <- sample_names

# assemble data.frame: first column gene names, then samples
tcga_brca_luma_dataset <- data.frame(Gene = genes, counts_mat, stringsAsFactors = FALSE)

# Save to inst/extdata so system.file can find it during vignette build
save(tcga_brca_luma_dataset, file = file.path("inst", "extdata", "tcga_brca_luma_dataset.RData"))
cat("Wrote", file.path("inst", "extdata", "tcga_brca_luma_dataset.RData"), "\n")
