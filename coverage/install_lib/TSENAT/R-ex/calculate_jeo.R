### Name: calculate_jeo
### Title: Jackknife resampling with confidence intervals
### Aliases: calculate_jeo

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

# TPM and effective_length REQUIRED for filter_analysis()
tpm <- matrix(runif(nrow(readcounts) * ncol(readcounts), 0.1, 10),
              nrow = nrow(readcounts), ncol = ncol(readcounts),
              dimnames = dimnames(readcounts))
effective_length <- matrix(100, nrow = nrow(readcounts), ncol = ncol(readcounts))

# Create config (metadata passed as explicit parameter to build_analysis)
config <- TSENAT_config(
  sample_col = 'sample',
  condition_col = 'condition',
  q = seq(0, 2, by = 0.05),
  paired = FALSE
)

# Build analysis from vignette data - metadata as explicit parameter
analysis <- build_analysis(
  readcounts = readcounts,
  metadata = metadata_df,
  tx2gene = gff3_dataset,
  config = config,
  tpm = tpm,
  effective_length = effective_length
)

# Filter low-abundance genes (required for reliable jackknife estimates)
analysis <- filter_analysis(analysis, stringency = 'severe')

# Compute diversity first (required for jackknife)
analysis <- calculate_diversity(analysis, q = c(0.5, 1.0, 1.5))

# Run jackknife estimation
analysis <- calculate_jeo(analysis, q = c(0.5, 1.0, 1.5))
# Check jackknife results using unified accessor
jackknife_res <- results(analysis, type = 'jackknife')
if (!is.null(jackknife_res)) names(jackknife_res)




