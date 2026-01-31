## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  fig.width = 10,
  fig.height = 6
)

## ----setup--------------------------------------------------------------------
# Load packages
suppressPackageStartupMessages({
  library(TSENAT)
  library(ggplot2)
  library(SummarizedExperiment)
  library(mgcv)
})

# Load required files
coldata_tsv <- system.file("extdata", "coldata.tsv", package = "TSENAT")
tx2gene_tsv <- system.file("extdata","tx2gene.tsv", package = "TSENAT")

# Load dataset
data(tcga_brca_luma_dataset)

# Extract gene names and read count data (do not reference ts_se yet)
genes <- tcga_brca_luma_dataset[, 1]
readcounts <- tcga_brca_luma_dataset[, -1]

# If a tx2gene mapping is available and matches the `genes` column, assign
# transcript IDs as rownames of `readcounts` so downstream transcript-level
# plotting functions can use them.
if (nzchar(tx2gene_tsv) && file.exists(tx2gene_tsv)) {
  txmap <- utils::read.delim(tx2gene_tsv, header = TRUE, stringsAsFactors = FALSE)
  if (nrow(txmap) == nrow(readcounts) && all(as.character(txmap$Gen) == as.character(genes))) {
    rownames(readcounts) <- as.character(txmap$Transcript)
  }
}

# Check gene names
head(genes)

# Check read count dataset
dim(readcounts)
head(readcounts[, 1:5])


## ----read-coldata, echo=TRUE--------------------------------------------------
# Read external coldata if available and show a preview
coldata_df <- NULL
if (nzchar(coldata_tsv) && file.exists(coldata_tsv)) {
  coldata_df <- read.table(coldata_tsv, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
  head(coldata_df)
} else {
  message("No external coldata.tsv found via system.file()")
}

## ----readfilter---------------------------------------------------------------
tokeep <- rowSums(readcounts > 5) > 5

readcounts <- readcounts[tokeep, ]
genes      <- genes[tokeep]

## ----pseudocount-example, echo=TRUE-------------------------------------------
# Option A: add a small integer pseudocount (document this choice)
readcounts_pc <- readcounts + 1L

# Option B: keep original counts but filter lowly-observed transcripts (recommended)
# readcounts <- readcounts[rowSums(readcounts > 5) > 5, ]

## ----tsallis-calc-single------------------------------------------------------
# compute Tsallis entropy for q = 0.1
q <- 0.1
ts_se <- calculate_diversity(readcounts, genes, method = "tsallis", q = q, norm = TRUE)
head(assay(ts_se)[1:5,1:5])


## ----apply-coldata-to-se-single-----------------------------------------------
# If we read an external coldata table, map its Condition values to the
# SummarizedExperiment so plotting and downstream functions can use it.
if (!is.null(coldata_df)) {
  sample_base_names <- sub("_q=.*", "", colnames(SummarizedExperiment::assay(ts_se)))
  st_map <- setNames(coldata_df$Condition, coldata_df$Sample)
  sample_types <- unname(st_map[sample_base_names])
  missing_idx <- which(is.na(sample_types))
  if (length(missing_idx) > 0) {
    sample_types[missing_idx] <- infer_sample_group(sample_base_names[missing_idx])
  }
  SummarizedExperiment::colData(ts_se)$sample_type <- sample_types
}

## ----difference---------------------------------------------------------------
# create a sample grouping vector inferred from sample names
# account for per-q column names like 'Sample_q=0.01' by stripping the '_q=...' suffix
sample_base_names <- sub("_q=.*", "", colnames(assay(ts_se)))
if (exists("coldata_df") && !is.null(coldata_df) && "sample_type" %in% colnames(SummarizedExperiment::colData(ts_se))) {
  samples <- as.character(SummarizedExperiment::colData(ts_se)$sample_type)
} else {
  samples <- ifelse(grepl("_N$", sample_base_names), "Normal", "Tumor")
}

# prepare diversity table as data.frame with gene names in first column
div_df <- as.data.frame(assay(ts_se))
div_df <- cbind(genes = rowData(ts_se)$genes, div_df)

# run difference analysis (mean + wilcoxon)
# run difference analysis (mean + paired wilcoxon)
# samples in the example dataset are matched pairs (Normal/Tumor), so use a paired test
res <- calculate_difference(div_df, samples, control = "Normal", method = "mean", test = "wilcoxon", paired = TRUE)
# sort results by adjusted p-value
if ("adjusted_p_values" %in% colnames(res)) {
  res <- res[order(res$adjusted_p_values), , drop = FALSE]
}
head(res)


## ----ma-and-volcano, fig.width=8, fig.height=4--------------------------------
# MA plot using helper
p_ma <- plot_ma(res)
print(p_ma)

# Volcano plot: mean difference vs -log10(adjusted p-value)
if (all(c("mean_difference", "adjusted_p_values") %in% colnames(res))) {
  res$label <- apply(res[, c("mean_difference", "adjusted_p_values")], 1,
                     function(x) ifelse(abs(as.numeric(x[1])) >= 0.1 & as.numeric(x[2]) < 0.05,
                                        "significant", "non-significant"))
  p_volcano <- plot_volcano(res)
  print(p_volcano)
}

## ----top-transcripts-singleq, fig.width=10, fig.height=6----------------------
if (is.data.frame(res) && nrow(res) > 0) {
  # select significant genes (adj p < 0.05 when available) and take top 3 by significance
  if ("adjusted_p_values" %in% colnames(res)) {
    sig_res <- res[res$adjusted_p_values < 0.05, , drop = FALSE]
  } else if ("p_values" %in% colnames(res)) {
    sig_res <- res[res$p_values < 0.05, , drop = FALSE]
  } else {
    sig_res <- res
  }
  if (nrow(sig_res) == 0) {
    message("No significant genes found; falling back to top genes by p-value.")
    top_genes <- head(res$genes, 3)
  } else {
    top_genes <- head(as.character(sig_res$genes), 3)
  }
  # prepare sample grouping vector (from earlier in vignette)
  sample_base_names <- sub("_q=.*", "", colnames(assay(ts_se)))
  samples_vec <- if (exists("coldata_df") && !is.null(coldata_df) && "sample_type" %in% colnames(SummarizedExperiment::colData(ts_se))) {
    as.character(SummarizedExperiment::colData(ts_se)$sample_type)
  } else {
    ifelse(grepl("_N$", sample_base_names), "Normal", "Tumor")
  }
  counts_for_plot <- readcounts
  if (is.null(rownames(counts_for_plot))) {
    message("Skipping transcript plots: `readcounts` does not have transcript rownames.",
            " Provide transcript-level counts with rownames to use `plot_top_transcripts()`.")
  } else {
    p_comb <- plot_top_transcripts(counts_for_plot, gene = top_genes, samples = samples_vec, tx2gene = tx2gene_tsv, top_n = NULL)
    print(p_comb)
  }
}

## ----tsallis-calc-doble-------------------------------------------------------
# compute Tsallis entropy for q = 1 (normalized)
q <- c(0.1, 2)
ts_se <- calculate_diversity(readcounts, genes, method = "tsallis", q = q, norm = TRUE)
head(assay(ts_se)[1:5,1:5])


## ----apply-coldata-to-se-doble------------------------------------------------
# If we read an external coldata table, map its Condition values to the
# SummarizedExperiment so plotting and downstream functions can use it.
if (!is.null(coldata_df)) {
  sample_base_names <- sub("_q=.*", "", colnames(SummarizedExperiment::assay(ts_se)))
  st_map <- setNames(coldata_df$Condition, coldata_df$Sample)
  sample_types <- unname(st_map[sample_base_names])
  missing_idx <- which(is.na(sample_types))
  if (length(missing_idx) > 0) {
    sample_types[missing_idx] <- infer_sample_group(sample_base_names[missing_idx])
  }
  SummarizedExperiment::colData(ts_se)$sample_type <- sample_types
}

## ----plots, fig.width=10, fig.height=6----------------------------------------
## Use package plot helpers for multi-q violin and density
p_violin <- plot_tsallis_violin_multq(ts_se, assay_name = "diversity")
print(p_violin)

p_density <- plot_tsallis_density_multq(ts_se, assay_name = "diversity")
print(p_density)


## ----tsallis-calc-sequence----------------------------------------------------
# compute Tsallis entropy for a sequence of values (normalized)
qvec <- seq(0.01, 2, by = 0.1)
ts_se <- calculate_diversity(readcounts, genes, method = "tsallis", q = qvec, norm = TRUE)
head(assay(ts_se)[1:5,1:5])

## ----apply-coldata-to-se-sequence---------------------------------------------
# If we read an external coldata table, map its Condition values to the
# SummarizedExperiment so plotting and downstream functions can use it.
if (!is.null(coldata_df)) {
  sample_base_names <- sub("_q=.*", "", colnames(SummarizedExperiment::assay(ts_se)))
  st_map <- setNames(coldata_df$Condition, coldata_df$Sample)
  sample_types <- unname(st_map[sample_base_names])
  missing_idx <- which(is.na(sample_types))
  if (length(missing_idx) > 0) {
    sample_types[missing_idx] <- infer_sample_group(sample_base_names[missing_idx])
  }
  SummarizedExperiment::colData(ts_se)$sample_type <- sample_types
}

## ----plot q-curve, fig.width=10, fig.height=6---------------------------------

# q-curve: median ± IQR across q values by group
p3 <- plot_tsallis_q_curve(readcounts, genes, q_values = qvec)
print(p3)

## ----lm-interaction-----------------------------------------------------------
# ensure the SummarizedExperiment contains sample names with group suffixes
# (the function infers group from sample name suffix _N -> Normal)
lm_res <- calculate_lm_interaction(ts_se, sample_type_col = "sample_type", min_obs = 8)
head(lm_res)

## -----------------------------------------------------------------------------
sessionInfo()

