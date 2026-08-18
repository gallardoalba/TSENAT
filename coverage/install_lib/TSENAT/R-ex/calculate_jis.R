### Name: calculate_jis
### Title: Jackknife isoform switching analysis on TSENATAnalysis object
### Aliases: calculate_jis

### ** Examples

data(readcounts)
metadata_df <- read.table(system.file('extdata', 'metadata.tsv', package
= 'TSENAT'),
                          header = TRUE, sep = '\t')
gff3_file <- system.file('extdata', 'annotation.gff3.gz', package = 'TSENAT')
config <- TSENAT_config(sample_col = 'sample', condition_col = 'condition')
analysis <- build_analysis(readcounts = readcounts, tx2gene =
gff3_file, metadata = metadata_df, config = config,
                             tpm = tpm,
 effective_length = effective_length)
analysis <- filter_analysis(analysis, min_samples = 1, subset_n_genes
= 20, subset_n_samples = 8)
analysis <- calculate_diversity(analysis, q = 1)




