pkgname <- "TSENAT"
source(file.path(R.home("share"), "R", "examples-header.R"))
options(warn = 1)
base::assign(".ExTimings", "TSENAT-Ex.timings", pos = 'CheckExEnv')
base::cat("name\tuser\tsystem\telapsed\n", file=base::get(".ExTimings", pos = 'CheckExEnv'))
base::assign(".format_ptime",
function(x) {
  if(!is.na(x[4L])) x[1L] <- x[1L] + x[4L]
  if(!is.na(x[5L])) x[2L] <- x[2L] + x[5L]
  options(OutDec = '.')
  format(x[1L:3L], digits = 7L)
},
pos = 'CheckExEnv')

### * </HEADER>
library('TSENAT')

base::assign(".oldSearch", base::search(), pos = 'CheckExEnv')
base::assign(".old_wd", base::getwd(), pos = 'CheckExEnv')
cleanEx()
nameEx("calculate_difference")
### * calculate_difference

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: calculate_difference
### Title: Calculate splicing diversity changes between two conditions.
### Aliases: calculate_difference

### ** Examples

# data.frame with splicing diversity values
x <- data.frame(Genes = letters[seq_len(10)], matrix(runif(80), ncol = 8))

# sample categories
samples <- c(rep("Healthy", 4), rep("Pathogenic", 4))

# To calculate the difference of splicing diversity changes between the
# 'Healthy' and 'Pathogenic' condition together with the significance values,
# using mean and Wilcoxon rank sum test, use:
calculate_difference(x, samples, control = "Healthy", method = "mean", test = "wilcoxon")



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("calculate_difference", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("calculate_diversity")
### * calculate_diversity

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: calculate_diversity
### Title: Calculate Tsallis diversity per gene across samples
### Aliases: calculate_diversity

### ** Examples

data("tcga_brca_luma_dataset", package = "TSENAT")
rc <- as.matrix(tcga_brca_luma_dataset[1:20, -1, drop = FALSE])
gs <- tcga_brca_luma_dataset$genes[1:20]
se <- calculate_diversity(rc, gs, q = 0.1, norm = TRUE)
SummarizedExperiment::assay(se)[1:3, 1:3]



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("calculate_diversity", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("calculate_lm_interaction")
### * calculate_lm_interaction

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: calculate_lm_interaction
### Title: Linear-model interaction test for Tsallis entropy
### Aliases: calculate_lm_interaction

### ** Examples

data("tcga_brca_luma_dataset", package = "TSENAT")
rc <- as.matrix(tcga_brca_luma_dataset[1:20, -1, drop = FALSE])
gs <- tcga_brca_luma_dataset$genes[1:20]
se <- calculate_diversity(rc, gs, q = c(0.1, 1), norm = TRUE)
calculate_lm_interaction(se)



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("calculate_lm_interaction", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("calculate_tsallis_entropy")
### * calculate_tsallis_entropy

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: calculate_tsallis_entropy
### Title: Calculate Tsallis entropy for a vector of transcript-level
###   expression values of one gene.
### Aliases: calculate_tsallis_entropy

### ** Examples

# Basic usage with a small numeric vector
x <- c(10, 5, 0)
calculate_tsallis_entropy(x, q = c(0.5, 1, 2), norm = TRUE)



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("calculate_tsallis_entropy", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("infer_sample_group")
### * infer_sample_group

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: infer_sample_group
### Title: Infer sample group from sample names
### Aliases: infer_sample_group

### ** Examples

# Basic usage: returns raw suffix tokens or TCGA two-digit codes when
# no mapping is supplied
infer_sample_group(c("S1_N", "TCGA-XX-01A"))

# Provide a suffix_map to translate tokens
infer_sample_group(c("S1_N", "S2_T"), suffix_map = c(N = "Normal", T = "Tumor"))

# Provide a TCGA mapping and prefer TCGA codes over suffixes
tcga_map <- c("01" = "Tumor")
infer_sample_group(c("TCGA-XX-01A", "Sample_N"), tcga_map = tcga_map,
  prefer_suffix = FALSE)



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("infer_sample_group", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("map_coldata_to_se")
### * map_coldata_to_se

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: map_coldata_to_se
### Title: Map external coldata into a SummarizedExperiment
### Aliases: map_coldata_to_se

### ** Examples

data("tcga_brca_luma_dataset", package = "TSENAT")
rc <- as.matrix(tcga_brca_luma_dataset[1:20, -1, drop = FALSE])
gs <- tcga_brca_luma_dataset$genes[1:20]
se <- calculate_diversity(rc, gs, q = 0.1, norm = TRUE)
sample_names <- sub("_q=.*", "", colnames(SummarizedExperiment::assay(se)))
coldata_df <- data.frame(
  Sample = sample_names,
  Condition = rep(c("A", "B"), length.out = ncol(se))
)
map_coldata_to_se(se, coldata_df)



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("map_coldata_to_se", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("plot_diversity_density")
### * plot_diversity_density

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: plot_diversity_density
### Title: Plot diversity distributions (density) by sample type
### Aliases: plot_diversity_density

### ** Examples

data("tcga_brca_luma_dataset", package = "TSENAT")
rc <- as.matrix(tcga_brca_luma_dataset[1:20, -1, drop = FALSE])
gs <- tcga_brca_luma_dataset$genes[1:20]
se <- calculate_diversity(rc, gs, q = 0.1, norm = TRUE)
plot_diversity_density(se)



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("plot_diversity_density", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("plot_ma")
### * plot_ma

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: plot_ma
### Title: Plot MA plot for difference results
### Aliases: plot_ma

### ** Examples

# Minimal fake diff_df
df <- data.frame(
  gene = paste0("g", seq_len(10)),
  sampleA_mean = runif(10),
  sampleB_mean = runif(10),
  log2_fold_change = rnorm(10),
  adjusted_p_values = runif(10)
)
plot_ma(df)



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("plot_ma", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("plot_mean_violin")
### * plot_mean_violin

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: plot_mean_violin
### Title: Plot violin of per-gene mean diversity by sample type
### Aliases: plot_mean_violin

### ** Examples

data("tcga_brca_luma_dataset", package = "TSENAT")
rc <- as.matrix(tcga_brca_luma_dataset[1:20, -1, drop = FALSE])
gs <- tcga_brca_luma_dataset$genes[1:20]
se <- calculate_diversity(rc, gs, q = 0.1, norm = TRUE)
plot_mean_violin(se)



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("plot_mean_violin", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("plot_top_transcripts")
### * plot_top_transcripts

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: plot_top_transcripts
### Title: Plot top transcripts for a gene
### Aliases: plot_top_transcripts

### ** Examples

# Toy transcript-level example: 6 transcripts across 4 samples
tx_counts <- matrix(sample(1:100, 24, replace = TRUE), nrow = 6)
rownames(tx_counts) <- paste0("tx", seq_len(nrow(tx_counts)))
colnames(tx_counts) <- paste0("S", seq_len(ncol(tx_counts)))
# Map transcripts to genes (3 genes, 2 transcripts each)
tx2gene <- data.frame(
  Transcript = rownames(tx_counts),
  Gen = rep(paste0("G", seq_len(3)), each = 2),
  stringsAsFactors = FALSE
)
# Sample group labels (length = ncol(tx_counts))
samples <- rep(c("Normal", "Tumor"), length.out = ncol(tx_counts))
plot_top_transcripts(
  tx_counts,
  gene = c("G1", "G2"),
  samples = samples,
  tx2gene = tx2gene,
  top_n = 2
)



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("plot_top_transcripts", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("plot_tsallis_density_multq")
### * plot_tsallis_density_multq

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: plot_tsallis_density_multq
### Title: Density plot of Tsallis entropy for multiple q values
### Aliases: plot_tsallis_density_multq

### ** Examples

data("tcga_brca_luma_dataset", package = "TSENAT")
rc <- as.matrix(tcga_brca_luma_dataset[1:20, -1, drop = FALSE])
gs <- tcga_brca_luma_dataset$genes[1:20]
se <- calculate_diversity(rc, gs, q = c(0.1, 1), norm = TRUE)
plot_tsallis_density_multq(se)



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("plot_tsallis_density_multq", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("plot_tsallis_q_curve")
### * plot_tsallis_q_curve

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: plot_tsallis_q_curve
### Title: Plot median +- IQR of Tsallis entropy across q values by group
### Aliases: plot_tsallis_q_curve

### ** Examples

data("tcga_brca_luma_dataset", package = "TSENAT")
rc <- as.matrix(tcga_brca_luma_dataset[1:40, -1, drop = FALSE])
gs <- tcga_brca_luma_dataset$genes[1:40]
p <- plot_tsallis_q_curve(rc, gs, q_values = seq(0.01, 0.1, by = 0.03))
p



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("plot_tsallis_q_curve", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("plot_tsallis_violin_multq")
### * plot_tsallis_violin_multq

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: plot_tsallis_violin_multq
### Title: Violin plot of Tsallis entropy for multiple q values
### Aliases: plot_tsallis_violin_multq

### ** Examples

data("tcga_brca_luma_dataset", package = "TSENAT")
rc <- as.matrix(tcga_brca_luma_dataset[1:20, -1, drop = FALSE])
gs <- tcga_brca_luma_dataset$genes[1:20]
se <- calculate_diversity(rc, gs, q = c(0.1, 1), norm = TRUE)
plot_tsallis_violin_multq(se)



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("plot_tsallis_violin_multq", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("plot_volcano")
### * plot_volcano

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: plot_volcano
### Title: Volcano plot for differential results
### Aliases: plot_volcano

### ** Examples

# Minimal fake diff_df
df <- data.frame(
  gene = paste0("g", seq_len(10)),
  mean_difference = runif(10),
  adjusted_p_values = runif(10)
)
plot_volcano(df, x_col = "mean_difference", padj_col = "adjusted_p_values")



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("plot_volcano", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
### * <FOOTER>
###
cleanEx()
options(digits = 7L)
base::cat("Time elapsed: ", proc.time() - base::get("ptime", pos = 'CheckExEnv'),"\n")
grDevices::dev.off()
###
### Local variables: ***
### mode: outline-minor ***
### outline-regexp: "\\(> \\)?### [*]+" ***
### End: ***
quit('no')
