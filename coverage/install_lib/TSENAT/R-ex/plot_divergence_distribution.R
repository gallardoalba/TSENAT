### Name: plot_divergence_distribution
### Title: Plot Tsallis Divergence Effect Size Distribution (S4 Wrapper)
### Aliases: plot_divergence_distribution

### ** Examples

# Plot 2: Distribution of effect sizes across genes
data(readcounts)
readcounts <- as.matrix(readcounts)
mode(readcounts) <- 'numeric'
metadata_df <- read.table(
  system.file('extdata', 'metadata.tsv', package = 'TSENAT'),
  header = TRUE, sep = '\t'
)
gff3_dataset <- system.file('extdata', 'annotation.gff3.gz', package =
'TSENAT')
# Configure analysis parameters first
config <- TSENAT_config(
  sample_col = 'sample',
  condition_col = 'condition',
  control = 'normal'
)

# Build analysis with configured parameters
analysis <- build_analysis(
  readcounts = readcounts,
  tx2gene = gff3_dataset,
  metadata = metadata_df,
  config = config,
  tpm = tpm,
  effective_length = effective_length
)

analysis <- filter_analysis(analysis, stringency = 'severe')
analysis <- calculate_diversity(
  analysis,
  q = seq(0.2, 2, by = 0.4),
  verbose = FALSE
)
analysis <- calculate_divergence(
  analysis,
  q = seq(0.2, 2, by = 0.4),
  verbose = FALSE
)
analysis <- suppressWarnings(calculate_sait(analysis, method = 'lmm'))
analysis <- calculate_effect_sizes(analysis)
p_dist <- plot_divergence_distribution(analysis)
# print(p_dist)




