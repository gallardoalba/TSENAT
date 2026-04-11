# Set metadata for TSENATAnalysis object

Access or set the metadata list stored in a TSENATAnalysis object. The
getter function retrieves all metadata or a specific key-value. The
setter function replaces the entire metadata list.

## Usage

``` r
metadata(x) <- value

# S4 method for class 'TSENATAnalysis'
metadata(x, key = NULL)

# S4 method for class 'TSENATAnalysis'
metadata(x) <- value
```

## Arguments

- x:

  A
  [`TSENATAnalysis`](https://gallardoalba.github.io/TSENAT/reference/TSENATAnalysis-class.md)
  object

- value:

  A list of metadata to assign

- key:

  Optional character string specifying a metadata key to retrieve

## Value

TSENATAnalysis object with updated metadata

- `metadata`:

  Returns the full metadata list, or a single value if `key` is
  specified

- `metadata<-`:

  Returns the modified `TSENATAnalysis` object

## Details

Get or Set Metadata

## Examples

``` r
# \donttest{
# For complete examples including the metadata<- setter method, 
# see the TSENATAnalysis-metadata help page
library(SummarizedExperiment)
#> Loading required package: MatrixGenerics
#> Loading required package: matrixStats
#> 
#> Attaching package: ‘MatrixGenerics’
#> The following objects are masked from ‘package:matrixStats’:
#> 
#>     colAlls, colAnyNAs, colAnys, colAvgsPerRowSet, colCollapse,
#>     colCounts, colCummaxs, colCummins, colCumprods, colCumsums,
#>     colDiffs, colIQRDiffs, colIQRs, colLogSumExps, colMadDiffs,
#>     colMads, colMaxs, colMeans2, colMedians, colMins, colOrderStats,
#>     colProds, colQuantiles, colRanges, colRanks, colSdDiffs, colSds,
#>     colSums2, colTabulates, colVarDiffs, colVars, colWeightedMads,
#>     colWeightedMeans, colWeightedMedians, colWeightedSds,
#>     colWeightedVars, rowAlls, rowAnyNAs, rowAnys, rowAvgsPerColSet,
#>     rowCollapse, rowCounts, rowCummaxs, rowCummins, rowCumprods,
#>     rowCumsums, rowDiffs, rowIQRDiffs, rowIQRs, rowLogSumExps,
#>     rowMadDiffs, rowMads, rowMaxs, rowMeans2, rowMedians, rowMins,
#>     rowOrderStats, rowProds, rowQuantiles, rowRanges, rowRanks,
#>     rowSdDiffs, rowSds, rowSums2, rowTabulates, rowVarDiffs, rowVars,
#>     rowWeightedMads, rowWeightedMeans, rowWeightedMedians,
#>     rowWeightedSds, rowWeightedVars
#> Loading required package: GenomicRanges
#> Loading required package: stats4
#> Loading required package: BiocGenerics
#> Loading required package: generics
#> 
#> Attaching package: ‘generics’
#> The following objects are masked from ‘package:base’:
#> 
#>     as.difftime, as.factor, as.ordered, intersect, is.element, setdiff,
#>     setequal, union
#> 
#> Attaching package: ‘BiocGenerics’
#> The following objects are masked from ‘package:stats’:
#> 
#>     IQR, mad, sd, var, xtabs
#> The following objects are masked from ‘package:base’:
#> 
#>     Filter, Find, Map, Position, Reduce, anyDuplicated, aperm, append,
#>     as.data.frame, basename, cbind, colnames, dirname, do.call,
#>     duplicated, eval, evalq, get, grep, grepl, is.unsorted, lapply,
#>     mapply, match, mget, order, paste, pmax, pmax.int, pmin, pmin.int,
#>     rank, rbind, rownames, sapply, saveRDS, table, tapply, unique,
#>     unsplit, which.max, which.min
#> Loading required package: S4Vectors
#> 
#> Attaching package: ‘S4Vectors’
#> The following object is masked from ‘package:TSENAT’:
#> 
#>     metadata<-
#> The following object is masked from ‘package:utils’:
#> 
#>     findMatches
#> The following objects are masked from ‘package:base’:
#> 
#>     I, expand.grid, unname
#> Loading required package: IRanges
#> Loading required package: Seqinfo
#> Loading required package: Biobase
#> Welcome to Bioconductor
#> 
#>     Vignettes contain introductory material; view with
#>     'browseVignettes()'. To cite Bioconductor, see
#>     'citation("Biobase")', and for packages 'citation("pkgname")'.
#> 
#> Attaching package: ‘Biobase’
#> The following object is masked from ‘package:MatrixGenerics’:
#> 
#>     rowMedians
#> The following objects are masked from ‘package:matrixStats’:
#> 
#>     anyMissing, rowMedians
se <- SummarizedExperiment(assays = list(counts = matrix(1:100, nrow = 10)))
analysis <- new('TSENATAnalysis', se = se, config = list())

# Get metadata (this works immediately)
metadata(analysis)
#> $function_calls
#> character(0)
#> 
#> $function_timestamps
#> character(0)
#> 
# }
# Create a TSENATAnalysis object
library(SummarizedExperiment)
se <- SummarizedExperiment(assays = list(counts = matrix(1:100, nrow = 10)))
analysis <- new('TSENATAnalysis', se = se, config = list())

# Get metadata (empty by default)
metadata(analysis)
#> $function_calls
#> character(0)
#> 
#> $function_timestamps
#> character(0)
#> 

# \donttest{
# Set metadata
metadata(analysis) <- list(processing_date = Sys.Date(), method = "test")
#> Error: unable to find an inherited method for function ‘metadata<-’ for signature ‘x = "TSENATAnalysis"’

# Retrieve all metadata
metadata(analysis)
#> $function_calls
#> character(0)
#> 
#> $function_timestamps
#> character(0)
#> 

# Retrieve specific metadata key
metadata(analysis, key = "method")
#> NULL
# }
```
