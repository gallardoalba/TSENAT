# Replace analysis configuration

Replace analysis configuration

## Usage

``` r
setConfig(object, value)

# S4 method for class 'TSENATAnalysis'
setConfig(object, value)
```

## Arguments

- object:

  `TSENATAnalysis` object.

- value:

  List or `TSENATConfig` object containing configuration settings.

## Value

Modified `TSENATAnalysis` object with updated @config.

Updated `TSENATAnalysis` object with replaced configuration.

## Details

Provides type-safe replacement of @config slot. Typically called once at
the start of an analysis via
[`tsenat_config()`](https://gallardoalba.github.io/TSENAT/reference/tsenat_config.md)
rather than directly.

Replace entire configuration in TSENATAnalysis

Replaces the entire configuration of a TSENATAnalysis object. This
method is useful when you need to apply a new set of configuration
parameters to an existing analysis object. All previous configuration
values are replaced with those in the new value object.

## Examples

``` r
# Load real TSENAT data
data(readcounts)
metadata_df <- read.table(system.file('extdata', 'metadata.tsv', package
= 'TSENAT'),
  header = TRUE, sep = '\t')
gff3_file <- system.file('extdata', 'annotation.gff3.gz', package = 'TSENAT')

# Build and subset analysis
config <- tsenat_config(sample_col = 'sample', condition_col = 'condition')
analysis <- build_analysis_s4(readcounts = readcounts, tx2gene =
gff3_file, metadata = metadata_df, config = config,
  tpm = tpm, effective_length = effective_length)
analysis <- filter_analysis_s4(analysis, min_samples = 1, subset_n_genes
= 200)

# Replace configuration with new settings
new_config <- tsenat_config(q_values = c(0.5, 1.0, 1.5), seed = 42)
analysis <- setConfig(analysis, new_config)

# Verify the new configuration was applied
current_config <- getConfig(analysis)
print(current_config$q_values)  # Shows c(0.5, 1.0, 1.5)
#> [1] 0.5 1.0 1.5

# Load real TSENAT data
data(readcounts)
metadata_df <- read.table(system.file('extdata', 'metadata.tsv', package
= 'TSENAT'),
  header = TRUE, sep = '\t')
gff3_file <- system.file('extdata', 'annotation.gff3.gz', package = 'TSENAT')

# Build and subset analysis
config <- tsenat_config(sample_col = 'sample', condition_col = 'condition')
analysis <- build_analysis_s4(readcounts = readcounts, tx2gene =
gff3_file, metadata = metadata_df, config = config,
  tpm = tpm, effective_length = effective_length)
analysis <- filter_analysis_s4(analysis, min_samples = 1, subset_n_genes
= 200)

# Use setConfig via the method (called by setConfig generic)
new_config <- tsenat_config(q_values = c(0.5, 1.0, 1.5, 2.0), seed = 123)
analysis <- setConfig(analysis, new_config)

# Verify and display
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
#>     - Q-spectrum: 0.50 to 2.00 (4 values)
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
#>   Created: 2026-04-09 02:26:12
#>   Package version: 0.99.0
#> 
```
