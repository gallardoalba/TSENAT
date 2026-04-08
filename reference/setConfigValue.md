# Set a single configuration value in TSENATAnalysis

Set a single configuration value in TSENATAnalysis

## Usage

``` r
setConfigValue(object, key, value)

# S4 method for class 'TSENATAnalysis'
setConfigValue(object, key, value)
```

## Arguments

- object:

  `TSENATAnalysis` object.

- key:

  Character. Name of the config item to set.

- value:

  Value to set for the config item.

## Value

Updated `TSENATAnalysis` object with the new config value.

## Details

Convenience method for setting a single configuration value without
needing to retrieve, merge, and set the entire config list. This
preserves all other configuration values while updating only the
specified key.

Unlike
[`setConfig`](https://gallardoalba.github.io/TSENAT/reference/setConfig.md),
which replaces the entire configuration, `setConfigValue` performs a
targeted update. It retrieves the current config, updates one key-value
pair, and stores the modified config back.

## See also

[`setConfig`](https://gallardoalba.github.io/TSENAT/reference/setConfig.md)
for replacing entire configuration,
[`getConfig`](https://gallardoalba.github.io/TSENAT/reference/getConfig.md)
for retrieving configuration,
[`tsenat_config`](https://gallardoalba.github.io/TSENAT/reference/tsenat_config.md)
for creating configuration objects

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

# Update a single configuration value while preserving others
analysis <- setConfigValue(analysis, 'q_values', c(0.5, 1.0, 1.5))

# Verify the update
config <- getConfig(analysis)
print(config$q_values)  # Shows c(0.5, 1.0, 1.5)
#> [1] 0.5 1.0 1.5

# Load real TSENAT data
data(readcounts)
metadata_df <- read.table(system.file('extdata', 'metadata.tsv', package
= 'TSENAT'),
  header = TRUE, sep = '\t')
gff3_file <- system.file('extdata', 'annotation.gff3.gz', package = 'TSENAT')

# Build and subset analysis with initial configuration
config <- tsenat_config(sample_col = 'sample', condition_col = 'condition')
analysis <- build_analysis_s4(readcounts = readcounts, tx2gene =
gff3_file, metadata = metadata_df, config = config,
  tpm = tpm, effective_length = effective_length)
analysis <- filter_analysis_s4(analysis, min_samples = 1, subset_n_genes
= 200)

# Use setConfigValue to update single configuration values
# This preserves all other config values
analysis <- setConfigValue(analysis, 'q_values', c(0.5, 1.5, 2.0))
analysis <- setConfigValue(analysis, 'seed', 456)

# Verify the updates
summary(analysis)
#> === TSENAT Analysis Summary ===
#> DATA:
#>   Genes:       200
#>   Samples:      16
#>   Assays:        1 (counts)
#> 
#> CONFIGURATION:
#>   q_values: 0.5, 1.5, 2
#>   condition_col: condition
#>   subject_col: <NULL>
#>   sample_col: sample
#>   paired: <logical>
#>   control: <NULL>
#>   p_threshold: 0.05
#>   fdr_threshold: 0.05
#>   significance_threshold: 0.05
#>   nboot: 1000
#>   bootstrap_method: percentile
#>   stringency: medium
#>   nthreads: 1
#>   norm: <logical>
#>   bootstrap: <logical>
#>   bootstrap_ci: 0.95
#>   norm_method: <NULL>
#>   pseudocount: <NULL>
#>   metadata: <data.frame>
#>   seed: 456
#> 
#> RESULTS:
#> 
#> METADATA:
#>   Created: 2026-04-08 04:30:23
#>   Package: TSENAT 0.99.0
#> 
```
