### Name: plot_sait
### Title: Plot GAM q-curves from TSENATAnalysis object
### Aliases: plot_sait

### ** Examples

## Not run: 
##D # Plot 3: GAM q-curves for genes with q-by-condition interactions
##D data(readcounts)
##D readcounts <- as.matrix(readcounts)
##D mode(readcounts) <- 'numeric'
##D metadata_df <- read.table(
##D   system.file('extdata', 'metadata.tsv', package = 'TSENAT'),
##D   header = TRUE, sep = '\t'
##D )
##D gff3_dataset <- system.file('extdata', 'annotation.gff3.gz', package =
##D 'TSENAT')
##D 
##D # Configure analysis parameters first
##D config <- TSENAT_config(
##D   sample_col = 'sample',
##D   condition_col = 'condition',
##D   subject_col = 'paired_samples',
##D   paired = TRUE,
##D   control = 'normal',
##D   q = seq(0.2, 2, by = 0.4)  # 5 unique q-values: 0.2, 0.6, 1.0, 1.4, 1.8
##D )
##D 
##D # Build analysis with configured parameters
##D analysis <- build_analysis(
##D   readcounts = readcounts,
##D   tx2gene = gff3_dataset,
##D   metadata = metadata_df,
##D   config = config,
##D   tpm = tpm,
##D   effective_length = effective_length
##D )
##D 
##D analysis <- filter_analysis(analysis, stringency = 'severe')
##D analysis <- calculate_diversity(analysis, q = seq(0.2, 2, by = 0.4))
##D analysis <- suppressWarnings(calculate_sait(analysis, method = 'gam'))
##D 
##D p_gam <- plot_sait(analysis, n_top = 2, sig_alpha = 0.15)
##D # print(p_gam)
## End(Not run)




