### Name: calculate_effect_sizes
### Title: Compute Effect Sizes from Divergence Results (S4 Wrapper)
### Aliases: calculate_effect_sizes

### ** Examples

# Setup: Create test analysis with divergence and SAIT interaction results
data(readcounts)
readcounts <- as.matrix(readcounts)
mode(readcounts) <- 'numeric'
metadata_df <- read.table(
  system.file('extdata', 'metadata.tsv', package = 'TSENAT'),
  header = TRUE, sep = '\t'
)
gff3_dataset <- system.file('extdata', 'annotation.gff3.gz', package =
'TSENAT')

# Configure analysis parameters (best practice for reproducibility)
## Not run: 
##D config <- TSENAT_config(
##D   sample_col = 'sample',
##D   condition_col = 'condition',
##D   subject_col = 'paired_samples',
##D   paired = TRUE,
##D   control = 'normal',
##D   q = seq(0, 2, by = 0.5)  # Multiple q-values for SAIT (5 unique values)
##D )
##D analysis <- build_analysis(readcounts = readcounts, tx2gene =
##D gff3_dataset, metadata = metadata_df, config = config,
##D   tpm = tpm, effective_length = effective_length)
##D 
##D analysis <- filter_analysis(analysis, stringency = 'severe')
##D analysis <- calculate_diversity(analysis)
##D analysis <- calculate_divergence(analysis)
##D analysis <- suppressWarnings(calculate_sait(analysis, method = 'gam'))
##D 
##D # Compute effect sizes from divergence results
##D analysis <- calculate_effect_sizes(analysis,
##D   significance_threshold = 0.05)
##D 
##D # Access results using unified results accessor
##D effect_size_results <- results(analysis, type = 'effect_sizes_divergence')
##D 
##D # View structure of results
##D str(effect_size_results, max.level = 1)
## End(Not run)




