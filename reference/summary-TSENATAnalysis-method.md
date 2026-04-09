# Summary method for TSENATAnalysis

Summary method for TSENATAnalysis

Detailed summary of TSENATAnalysis

## Usage

``` r
# S4 method for class 'TSENATAnalysis'
summary(object)

# S4 method for class 'TSENATAnalysis'
summary(object)
```

## Arguments

- object:

  `TSENATAnalysis` object.

## Value

Prints detailed summary (invisibly returns object).

## Details

Provides comprehensive analysis summary including dimensions, results
counts, validation status, and metadata tracking.

## Examples

``` r
# Load real TSENAT data
data(readcounts)
metadata_df <- read.table(system.file('extdata', 'metadata.tsv', package
= 'TSENAT'),
  header = TRUE, sep = '\t')
gff3_file <- system.file('extdata', 'annotation.gff3.gz', package = 'TSENAT')
config <- tsenat_config(sample_col = 'sample', condition_col = 'condition')
analysis <- build_analysis_s4(readcounts = readcounts, tx2gene =
gff3_file, metadata = metadata_df, config = config,
  tpm = tpm, effective_length = effective_length)
analysis <- filter_analysis_s4(analysis, min_samples = 1, subset_n_genes
= 200)
summary(analysis)
#> TSENAT Analysis Summary
#> =======================
#> DATA STRUCTURE:
#>   Genes:        200
#>   Samples:      16
#>   Assays:       counts
#> 
#> CONFIGURATION:
#>   Design & Filtering:
#>     - Design: unpaired
#>     - Stringency: medium
#>     - Normalization: enabled
#>   Diversity Metrics:
#>     - Q-spectrum: 0.00 to 2.00 (5 values)
#>     - Pseudocount: disabled
#>   Statistical Methods:
#>     - LM fitting: GAM
#>     - P-value correction: BH
#>     - Jackknife filtering: LM-based
#> 
#> ANALYSIS RESULTS:
#>   ✗ Diversity: not computed
#>   ✗ LM Interaction: not computed
#>   ✗ Jackknife Switching: not computed
#>   ✗ Divergence Metrics: not computed
#>   ✗ Visualizations: not generated
#> 
#> PROCESSING & METADATA:
#>   Created: 2026-04-09 02:26:14
#>   Package version: 0.99.0
#> 
```
