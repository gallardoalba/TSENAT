# Constructor for TSENATAnalysis objects

Creates a new TSENATAnalysis object with a SummarizedExperiment base and
optional initial configuration.

## Usage

``` r
TSENATAnalysis(se, config = list())
```

## Arguments

- se:

  `SummarizedExperiment`. The base expression data object.

- config:

  `list`. Optional initial configuration (usually set via
  [`TSENAT_config()`](https://gallardoalba.github.io/TSENAT/reference/TSENAT_config.md)
  instead).

## Value

A new `TSENATAnalysis` object.

## Details

The constructor initializes all slots with empty lists except @se, which
must be provided. The @metadata slot automatically records: - creation
timestamp - TSENAT package version - initial function call

## Examples

``` r
# Load real TSENAT data
data(readcounts)
metadata_df <- read.table(system.file('extdata', 'metadata.tsv', package
= 'TSENAT'),
  header = TRUE, sep = '\t')
gff3_file <- system.file('extdata', 'annotation.gff3.gz', package = 'TSENAT')
config <- TSENAT_config(sample_col = 'sample', condition_col = 'condition')
analysis <- build_analysis(readcounts = readcounts, tx2gene =
gff3_file, metadata = metadata_df, config = config,
  tpm = tpm, effective_length = effective_length)
analysis <- filter_analysis(analysis, min_samples = 1, subset_n_genes
= 200)
```
