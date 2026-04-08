# Extract analysis results from TSENATAnalysis object

Provides flexible access to diversity, divergence, and statistical test
results.

## Usage

``` r
getResults(analysis, type = "diversity", q = NULL, simplify = TRUE)
```

## Arguments

- analysis:

  `TSENATAnalysis` object containing computed results.

- type:

  `character`. Type of results to extract: 'diversity', 'divergence',
  'lm', 'jackknife', or 'q_interactions'. Default: 'diversity'.

- q:

  `numeric`. For diversity results, optionally filter by q-value.
  Default: NULL (return all q-values).

- simplify:

  `logical`. If TRUE and q is specified, return as vector instead of
  matrix. Default: TRUE.

## Value

Extracted results as data.frame, matrix, or list depending on type.
Returns NULL if requested result type not computed.

## Details

This function provides a consistent interface to access all computed
results from the TSENATAnalysis object, abstracting away internal
storage details.

## Examples

``` r
# Load example data
data(readcounts, package = 'TSENAT')

# Create TSENATAnalysis from count matrix
# For simple count matrices (no tx2gene mapping), use TSENATAnalysis directly
config <- tsenat_config(
  q_values = c(0.5, 1.0, 2.0),
  condition_col = 'group'
)
se <- SummarizedExperiment::SummarizedExperiment(
  assays = list(counts = readcounts),
  colData = data.frame(
    group = rep(c('A', 'B'), length.out = ncol(readcounts))
  )
)
analysis <- TSENATAnalysis(se = se, config = config)

# Run analysis to generate diversity results
analysis <- calculate_diversity_s4(analysis)
#> Warning: [calculate_diversity_s4] Assay for q=0.5 is empty
#> Warning: [calculate_diversity_s4] Assay for q=1 is empty
#> Warning: [calculate_diversity_s4] Assay for q=2 is empty

# Extract diversity results
div_results <- getResults(analysis, type = 'diversity')
if (!is.null(div_results)) {
  head(div_results, n = 3)
}
#> $q_0.500
#> class: SummarizedExperiment 
#> dim: 0 16 
#> metadata(8): q what ... bootstrap_ci se
#> assays(1): diversity
#> rownames(0):
#> rowData names(1): gene_id
#> colnames(16): SRR14800475 SRR14800476 ... SRR14800489 SRR14800490
#> colData names(2): group sample_id
#> 
#> $q_1.000
#> class: SummarizedExperiment 
#> dim: 0 16 
#> metadata(8): q what ... bootstrap_ci se
#> assays(1): diversity
#> rownames(0):
#> rowData names(1): gene_id
#> colnames(16): SRR14800475 SRR14800476 ... SRR14800489 SRR14800490
#> colData names(2): group sample_id
#> 
#> $q_2.000
#> class: SummarizedExperiment 
#> dim: 0 16 
#> metadata(8): q what ... bootstrap_ci se
#> assays(1): diversity
#> rownames(0):
#> rowData names(1): gene_id
#> colnames(16): SRR14800475 SRR14800476 ... SRR14800489 SRR14800490
#> colData names(2): group sample_id
#> 
```
