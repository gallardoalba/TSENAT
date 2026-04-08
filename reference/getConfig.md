# Get configuration from TSENATAnalysis

Retrieve the configuration parameters stored in a TSENATAnalysis object.
These parameters control analysis behavior including q-values,
normalization, and output settings.

## Usage

``` r
getConfig(object)
```

## Arguments

- object:

  `TSENATAnalysis` object.

## Value

`list` containing configuration parameters (q_values, method, etc.)

## Details

Configuration is stored in the @config slot and controls how downstream
analyses are performed. Use
[`setConfig`](https://gallardoalba.github.io/TSENAT/reference/setConfig.md)
to replace the entire configuration or
[`setConfigValue`](https://gallardoalba.github.io/TSENAT/reference/setConfigValue.md)
for targeted updates.

## See also

[`setConfig`](https://gallardoalba.github.io/TSENAT/reference/setConfig.md)
for replacing configuration,
[`setConfigValue`](https://gallardoalba.github.io/TSENAT/reference/setConfigValue.md)
for single value updates,
[`tsenat_config`](https://gallardoalba.github.io/TSENAT/reference/tsenat_config.md)
for creating configuration objects

## Examples

``` r
data(readcounts)
metadata_df <- read.table(system.file('extdata', 'metadata.tsv', 
  package = 'TSENAT'), header = TRUE, sep = '\t')
gff3_file <- system.file('extdata', 'annotation.gff3.gz', package = 'TSENAT')
config <- tsenat_config(sample_col = 'sample', condition_col = 'condition')
analysis <- build_analysis_s4(readcounts = readcounts, tx2gene = gff3_file,
  metadata = metadata_df, config = config, tpm = tpm, effective_length = effective_length)
config <- getConfig(analysis)
print(config$q_values)
#> [1] 0.0 0.5 1.0 1.5 2.0
```
