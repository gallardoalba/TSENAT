# Setter for divergence results stored in @divergence_results

Setter for divergence results stored in @divergence_results

Extract divergence results

## Usage

``` r
divergence(object) <- value

# S4 method for class 'TSENATAnalysis'
divergence(object) <- value

divergence(object, component = NULL)

# S4 method for class 'TSENATAnalysis'
divergence(object, component = NULL)
```

## Arguments

- object:

  `TSENATAnalysis` object.

- value:

  list. Named list of divergence results.

- component:

  `character`. Component to extract: NULL (all), 'tsallis_divergence',
  'effect_sizes', etc.

## Value

SummarizedExperiment or data.frame with divergence metrics.

## Details

Divergence results are stored in @divergence_results with component
names corresponding to different divergence metrics. Use
`divergence(analysis)` to retrieve all components or specify a component
for targeted extraction.

## See also

Other TSENATAnalysis accessors:
[`diversity`](https://gallardoalba.github.io/TSENAT/reference/diversity.md),
[`lmResults`](https://gallardoalba.github.io/TSENAT/reference/lmResults.md),
[`se`](https://gallardoalba.github.io/TSENAT/reference/se.md),
[`metadata`](https://rdrr.io/pkg/S4Vectors/man/Annotated-class.html)

## Examples

``` r
# Load real TSENAT data and calculate divergence
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
analysis <- calculate_diversity_s4(analysis, q = c(0.5, 1.0, 1.5))
#> Note: 12 genes excluded (< 75% valid values).
analysis <- calculate_divergence_s4(analysis)
div_res <- divergence(analysis)
```
