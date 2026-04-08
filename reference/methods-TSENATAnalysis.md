# TSENATAnalysis S4 Methods

Implementation of accessor and utility methods for TSENATAnalysis
objects. These methods follow Bioconductor best practices for accessing
object slots.

## Usage

``` r
# S4 method for class 'TSENATAnalysis'
getSE(object)

# S4 method for class 'TSENATAnalysis'
getMeta(object, key = NULL)

# S4 method for class 'TSENATAnalysis'
getConfig(object)

# S4 method for class 'TSENATAnalysis'
getPlot(object, type = NULL)

# S4 method for class 'TSENATAnalysis'
addPlot(object, type, plot, replace = FALSE)
```

## Value

Methods return different components of TSENATAnalysis: - getSE:
SummarizedExperiment object with count data - lmResults: list of linear
model fitting results - jackKnife: list of jackknife diagnostics -
diversity: list of SummarizedExperiments (one per q-value) or single SE
for specific q - divergence: list of divergence analysis results -
getMeta: list or atomic value of metadata - getConfig: list of analysis
configuration - getPlot: ggplot object or NULL - addPlot:
invisible(object) (adds plot to cache)
