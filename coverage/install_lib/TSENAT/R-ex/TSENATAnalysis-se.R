### Name: se
### Title: Extract SummarizedExperiment from TSENATAnalysis
### Aliases: se se,TSENATAnalysis-method

### ** Examples

# Create a simple TSENATAnalysis object
library(SummarizedExperiment)
counts <- matrix(rpois(100, lambda = 10), nrow = 10, ncol = 10)
colnames(counts) <- paste0('sample_', 1:10)
rownames(counts) <- paste0('gene_', 1:10)

se <- SummarizedExperiment(assays = list(counts = counts))
analysis <- new('TSENATAnalysis', se = se, config = list())

# Extract the SummarizedExperiment
extracted_se <- se(analysis)
dim(extracted_se)



