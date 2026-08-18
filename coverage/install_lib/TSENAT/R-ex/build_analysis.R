### Name: build_analysis
### Title: Build a Complete TSENATAnalysis Object
### Aliases: build_analysis

### ** Examples

# Create example transcript count data
set.seed(42)
n_genes <- 10
n_isoforms_per_gene <- 3
n_isoforms <- n_genes * n_isoforms_per_gene
n_samples <- 10

# Generate count matrix
counts <- matrix(rpois(n_isoforms * n_samples, lambda = 20),
                 nrow = n_isoforms, ncol = n_samples)
rownames(counts) <- paste0('TX_', 1:n_isoforms)
colnames(counts) <- paste0('Sample_', 1:n_samples)

# Create tx2gene mapping
tx2gene <- data.frame(
  Transcript = rownames(counts),
  Gene = rep(paste0('GENE_', 1:n_genes), each = n_isoforms_per_gene))

# Create sample metadata
metadata <- data.frame(
  sample = colnames(counts),
  condition = rep(c('control', 'treatment'), each = 5),
  row.names = colnames(counts))

# Build analysis object - use NAMED parameters to avoid confusion
# Method 1: With explicit tx2gene data.frame (most common)
config <- TSENAT_config(sample_col = 'sample', condition_col = 'condition')
analysis <- build_analysis(
  readcounts = counts,
  tx2gene = tx2gene,
  metadata = metadata,
  config = config)

# Verify the analysis object was created
analysis
print(dim(analysis))

# Method 2: From Salmon quantification folder
# Requires directory structure like:
#   salmon_output/
#     sample1/quant.sf
#     sample2/quant.sf
#     ...
# 
# salmon_dir <- '/path/to/salmon/directory'
# 
# First create sample metadata matching Salmon sample names
# salmon_metadata <- data.frame(
#   condition = c('control', 'control', 'treatment', 'treatment'),
#   row.names = c('sample1', 'sample2', 'sample3', 'sample4')
# )
# 
# analysis_salmon <- build_analysis(
#   salmon_dir = salmon_dir,
#   tx2gene = 'annotation.gff3.gz',  # Auto-parsed from GFF3
#   metadata = salmon_metadata
# )
#
# Method 3: Hybrid - Salmon counts with manual tx2gene
# analysis_hybrid <- build_analysis(
#   salmon_dir = salmon_dir,
#   tx2gene = tx2gene,  # data.frame instead of file
#   metadata = salmon_metadata
# )
#
# Method 4: Pass metadata via config (parameter resolution pattern)
# cfg <- TSENAT_config()
# cfg$metadata <- metadata
# analysis_with_config <- build_analysis(
#   readcounts = counts,
#   tx2gene = tx2gene,
#   config = cfg
#   # Note: metadata argument omitted - will be read from config$metadata
# )

# Advanced: Assigning metadata to assays after object creation
# When adding metadata to SummarizedExperiment assays, always use the
# S4Vectors namespace to ensure proper method dispatch:
#   
#   se <- getSE(analysis)
#   assay_with_ci <- SummarizedExperiment::assay(se, 'log2fc_ci')
#   S4Vectors::metadata(assay_with_ci)$lower <- ci_lower_bounds
#   S4Vectors::metadata(assay_with_ci)$upper <- ci_upper_bounds
#
# Note: Avoid using metadata(assay) without the namespace - this can
# cause silent failures in S4 object metadata assignment.




