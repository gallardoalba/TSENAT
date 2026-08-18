### Name: plot_jis_delta
### Title: Plot Multi-Q Delta Influence Heatmaps from TSENATAnalysis Object
### Aliases: plot_jis_delta

### ** Examples

# Plot 5: Multi-q delta influence (isoform switching) heatmaps
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
  subject_col = 'paired_samples',
  paired = TRUE,
  control = 'normal',
  q = seq(0, 2, by = 0.1)
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
  q = seq(0.2, 2, by = 0.4)
)
analysis <- suppressWarnings(calculate_sait(analysis, method = 'lmm'))
analysis <- calculate_jis(
  analysis,
  q = seq(0.2, 2, by = 0.4),
  n_bootstrap = 20
)
heatmap_file <- plot_jis_delta(analysis, n_genes = 2)




