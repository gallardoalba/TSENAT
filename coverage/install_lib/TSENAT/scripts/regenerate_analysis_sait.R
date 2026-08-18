# =============================================================================
# regenerate_analysis_sait.R — regenera inst/extdata/analysis_sait.rds
# =============================================================================
# Replicates the vignettes/TSENAT.Rmd pipeline (build → filter → diversity →
# m-estimation → calculate_sait(method="gam")) with the CURRENT package and
# reduces the TSENATAnalysis to what is strictly needed for
# vignettes/TSENAT_appendix_B.Rmd (GAMM vs ART concordance):
#   - @se: original SummarizedExperiment (gene annotations)
#   - @config: analysis configuration
#   - @sait_results: results of the q × condition interaction
# The remaining slots are emptied (diversity/rank/jackknife/divergence/plots)
# and metadata keeps created_at + package_version.
#
# Usage: Rscript inst/scripts/regenerate_analysis_sait.R
# =============================================================================

# Resolve the package root: this script lives in <pkg>/inst/scripts/
script_dir <- dirname(sub("--file=", "", grep("--file=", commandArgs(), value = TRUE)))
if (length(script_dir) == 0 || !nzchar(script_dir)) {
    script_dir <- normalizePath(".")
}
pkg_root <- normalizePath(file.path(script_dir, "..", ".."))

pkgload::load_all(pkg_root, quiet = TRUE)
set.seed(42)

# --- Example data (lazy-loaded as 'readcounts') -----------------------------
utils::data(readcounts, package = "TSENAT")
readcounts <- as.matrix(readcounts)

metadata_df <- read.table(
    system.file("extdata", "metadata.tsv", package = "TSENAT"),
    header = TRUE, sep = "\t"
)
gff3_file <- system.file("extdata", "annotation.gff3.gz", package = "TSENAT")

# --- Configuration (identical to TSENAT.Rmd) --------------------------------
config <- TSENAT_config(
    sample_col = "sample",
    condition_col = "condition",
    subject_col = "paired_samples",
    q = seq(0, 2, by = 0.05),
    nthreads = 2,
    paired = TRUE,
    control = "normal"
)

# --- Vignette pipeline -------------------------------------------------------
analysis <- build_analysis(
    config = config,
    readcounts = readcounts,
    metadata = metadata_df,
    tx2gene = gff3_file,
    tpm = tpm,
    effective_length = effective_length
)

analysis <- filter_analysis(analysis, stringency = "medium")

set.seed(12345)
analysis <- calculate_diversity(analysis, norm = TRUE)

analysis <- calculate_m_estimator(
    analysis,
    loss_type = "huber",
    q_combine_method = "mean",
    influence_threshold = 0.75
)

analysis <- calculate_sait(
    analysis,
    method = "gam",
    multicorr = "hochberg"
)

sait_res <- results(analysis, type = "sait")
cat(sprintf("[regenerate_analysis_sait] SAIT: %d genes × %d columns\n",
    nrow(sait_res), ncol(sait_res)))

# --- Reduce to what is needed -----------------------------------------------
analysis@diversity_results <- list()
analysis@rank_test_results <- list()
analysis@jackknife_results <- list()
analysis@divergence_results <- list()
analysis@plots <- list()
analysis@metadata <- list(
    created_at = Sys.time(),
    package_version = as.character(utils::packageVersion("TSENAT"))
)

out <- file.path(pkg_root, "inst", "extdata", "analysis_sait.rds")
saveRDS(analysis, out, version = 2)
cat(sprintf("[regenerate_analysis_sait] saved: %s (%d bytes)\n", out,
    file.info(out)$size))
