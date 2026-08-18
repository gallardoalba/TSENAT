### Name: plot_concordance
### Title: Plot method concordance results from TSENATAnalysis
### Aliases: plot_concordance plot_concordance,TSENATAnalysis-method

### ** Examples

# Load example data (matching TSENAT.Rmd workflow)
data(readcounts)
readcounts <- as.matrix(readcounts)
mode(readcounts) <- 'numeric'
metadata_df <- read.table(
  system.file('extdata', 'metadata.tsv', package = 'TSENAT'),
  header = TRUE, sep = '\t'
)
gff3_dataset <- system.file('extdata', 'annotation.gff3.gz', package =
'TSENAT')

# Build analysis from vignette data and create small subset
config <- TSENAT_config(sample_col = 'sample', condition_col = 'condition')
analysis <- build_analysis(readcounts = readcounts, tx2gene =
gff3_dataset, metadata = metadata_df, config = config,
  tpm = tpm, effective_length = effective_length)
analysis <- filter_analysis(analysis, min_samples = 1, subset_n_genes
= 200)

# Note: calculate_concordance requires additional LM and Conover-Iman Rank Transform
# results computed. For demo purposes, we show that
# plot_concordance needs pre-computed concordance in @metadata




