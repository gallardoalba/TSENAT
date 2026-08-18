### Name: calculate_diversity
### Title: Calculate diversity and store in TSENATAnalysis
### Aliases: calculate_diversity

### ** Examples

# Load vignette data and build analysis
data(readcounts)
metadata_df <- read.table(
  system.file('extdata', 'metadata.tsv', package = 'TSENAT'),
  header = TRUE, sep = '\t'
)
gff3_dataset <- system.file('extdata', 'annotation.gff3.gz', package =
'TSENAT')
readcounts <- as.matrix(readcounts)
mode(readcounts) <- 'numeric'

config <- TSENAT_config(sample_col = 'sample', condition_col = 'condition')
analysis <- build_analysis(readcounts = readcounts, tx2gene =
gff3_dataset, metadata = metadata_df, config = config,
  tpm = tpm, effective_length = effective_length)

# Filter to manageable size (use 200+ genes to survive diversity filtering)
analysis <- filter_analysis(analysis, min_samples = 1, subset_n_genes
= 200)

# Compute diversity and access results using unified accessor
analysis <- calculate_diversity(analysis, q = c(0.5, 1.0), verbose =
FALSE)
head(results(analysis, type = 'diversity', q = 1.0))




