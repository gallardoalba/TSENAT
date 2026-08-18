### Name: calculate_divergence
### Title: Calculate divergence metrics and store in TSENATAnalysis
### Aliases: calculate_divergence

### ** Examples

# Load example data (matching TSENAT.Rmd workflow)
data(readcounts)
readcounts <- as.matrix(readcounts)
mode(readcounts) <- 'numeric'
metadata_df <- read.table(
  system.file('extdata', 'metadata.tsv', package = 'TSENAT'),
  header = TRUE, sep = '\t'
)
gff3_dataset <- system.file('extdata', 'annotation.gff3.gz', package = 'TSENAT')

# Build analysis from vignette data and create manageable subset
config <- TSENAT_config(sample_col = 'sample', condition_col = 'condition')
analysis <- build_analysis(
  readcounts = readcounts,
  tx2gene = gff3_dataset,
  metadata = metadata_df,
  config = config,
  tpm = tpm,
  effective_length = effective_length
)
# Use 200+ genes to ensure diversity filtering doesn't remove all genes
analysis <- filter_analysis(analysis, min_samples = 1, subset_n_genes = 200)

# Compute diversity first (required for divergence)
analysis <- calculate_diversity(analysis, q = c(0.5, 1.0, 1.5))

# Calculate divergence across q-values
analysis <- calculate_divergence(analysis, q = c(0.5, 1.0, 1.5))

# Check divergence results using unified accessor
head(results(analysis, type = 'divergence'))




