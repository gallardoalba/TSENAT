# Extract SummarizedExperiment from TSENATAnalysis

Extract SummarizedExperiment from TSENATAnalysis

## Usage

``` r
se(object)

# S4 method for class 'TSENATAnalysis'
se(object)
```

## Arguments

- object:

  `TSENATAnalysis` object.

## Value

The `SummarizedExperiment` containing transcript/gene counts.

## Details

Provides type-safe accessor for the embedded `SummarizedExperiment`.

## See also

Other TSENATAnalysis accessors:
[`diversity`](https://gallardoalba.github.io/TSENAT/reference/diversity.md),
[`divergence`](https://gallardoalba.github.io/TSENAT/reference/divergence.md),
[`lmResults`](https://gallardoalba.github.io/TSENAT/reference/lmResults.md),
[`metadata`](https://rdrr.io/pkg/S4Vectors/man/Annotated-class.html)

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
```
