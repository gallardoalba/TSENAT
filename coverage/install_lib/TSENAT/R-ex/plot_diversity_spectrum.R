### Name: plot_diversity_spectrum
### Title: Plot Tsallis Entropy q-Curve
### Aliases: plot_diversity_spectrum

### ** Examples

# Plot 7: Tsallis entropy q-curve (combined across all sample diversity)
data(readcounts)
metadata_df <- read.table(
  system.file('extdata', 'metadata.tsv', package = 'TSENAT'),
  header = TRUE, sep = '\t'
)
gff3_dataset <- system.file('extdata', 'annotation.gff3.gz', package =
'TSENAT')
readcounts <- as.matrix(readcounts)
mode(readcounts) <- 'numeric'

# Create configuration (required when metadata is provided)
config <- TSENAT_config(sample_col = 'sample', condition_col = 'condition')
analysis <- build_analysis(readcounts = readcounts, tx2gene =
gff3_dataset, metadata = metadata_df, config = config,
  tpm = tpm, effective_length = effective_length)
analysis <- filter_analysis(analysis, min_samples = 1, subset_n_genes
= 200)
analysis <- calculate_diversity(analysis, q = seq(0, 2, by = 0.5),
)
p <- plot_diversity_spectrum(analysis)
# print(p)




