## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----setup--------------------------------------------------------------------
library(abyssEdge)
library(SummarizedExperiment)

# load example dataset (if available in package)
data_exists <- file.exists(system.file("extdata/tcga_brca_luma_dataset.RData", package = "abyssEdge"))
if (data_exists) {
  data(tcga_brca_luma_dataset)
  genes <- tcga_brca_luma_dataset[, 1]
  readcounts <- tcga_brca_luma_dataset[, -1]
} else {
  # fall back to a small synthetic example
  set.seed(1)
  readcounts <- matrix(rpois(300, lambda = 10), nrow = 100, ncol = 3)
  colnames(readcounts) <- c("Sample1_N", "Sample2_T", "Sample3_N")
  genes <- rep(paste0("G", 1:20), length.out = nrow(readcounts))
}

# Ensure we have at least two sample groups for the vignette examples; if not,
# fall back to a balanced synthetic example with two groups.
sample_names <- colnames(readcounts)
groups <- if (!is.null(sample_names)) {
  ifelse(grepl("_N$", sample_names), "Normal",
         ifelse(grepl("_T$", sample_names) | grepl("_Tumor$", sample_names), "Tumor", NA))
} else {
  rep(NA_character_, ncol(readcounts))
}
if (length(unique(na.omit(groups))) < 2) {
  set.seed(42)
  readcounts <- matrix(rpois(100 * 6, lambda = 10), nrow = 100, ncol = 6)
  colnames(readcounts) <- c(paste0("S", 1:3, "_N"), paste0("S", 1:3, "_T"))
  genes <- paste0("G", seq_len(nrow(readcounts)))
}

## ----tsallis-calc-------------------------------------------------------------
# compute Tsallis entropy for q = 0.01 and q = 2 (normalized)
qvec <- c(0.01, 2)
ts_se <- calculate_diversity(readcounts, genes, method = "tsallis", q = qvec, norm = TRUE)
assay(ts_se)[1:5, ]

## ----difference---------------------------------------------------------------
# create a sample grouping vector inferred from sample names
# account for per-q column names like 'Sample_q=0.01' by stripping the '_q=...' suffix
sample_base_names <- sub("_q=.*", "", colnames(assay(ts_se)))
samples <- ifelse(grepl("_N$", sample_base_names), "Normal", "Tumor")

# prepare diversity table as data.frame with gene names in first column
div_df <- as.data.frame(assay(ts_se))
div_df <- cbind(genes = rowData(ts_se)$genes, div_df)

# run difference analysis (mean + wilcoxon)
res <- calculate_difference(div_df, samples, control = "Normal", method = "mean", test = "wilcoxon")
head(res)

## ----lm-interaction-----------------------------------------------------------
# ensure the SummarizedExperiment contains sample names with group suffixes
# (the function infers group from sample name suffix _N -> Normal)
lm_res <- calculate_lm_interaction(ts_se, sample_type_col = NULL, min_obs = 8)
head(lm_res)

## ----plots, fig.width=6, fig.height=4-----------------------------------------
# density plot of diversity values
p1 <- plot_diversity_density(ts_se, assay_name = "diversity")
print(p1)

# violin of per-gene means
p2 <- plot_mean_violin(ts_se, assay_name = "diversity")
print(p2)

# q-curve: median ± IQR across q values by group
p3 <- plot_tsallis_q_curve(readcounts, genes, q_values = c(0.01, 0.1, 0.5, 2))
print(p3)

## -----------------------------------------------------------------------------
sessionInfo()

