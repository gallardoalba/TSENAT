### Name: metadata,TSENATAnalysis-method
### Title: Metadata Accessor Methods
### Aliases: metadata,TSENATAnalysis-method
###   metadata<-,TSENATAnalysis-method

### ** Examples

# Create a TSENATAnalysis object
library(SummarizedExperiment)
se <- SummarizedExperiment(assays = list(counts = matrix(1:100, nrow = 10)))
analysis <- new('TSENATAnalysis', se = se, config = list())

# Get metadata (empty by default)
metadata(analysis)

## No test: 
# Set metadata
metadata(analysis) <- list(processing_date = Sys.Date(), method = 'test')

# Retrieve all metadata
metadata(analysis)

# Retrieve specific metadata key
metadata(analysis, key = 'method')
## End(No test)




