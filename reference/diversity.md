# Extract diversity results

Extract diversity results

Setter for diversity results stored in @diversity_results

## Usage

``` r
diversity(object, q = NULL)

# S4 method for class 'TSENATAnalysis'
diversity(object, q = NULL)

diversity(object) <- value

# S4 method for class 'TSENATAnalysis'
diversity(object) <- value
```

## Arguments

- object:

  TSENATAnalysis object

- q:

  `numeric`. Q-value to extract (e.g., 1.0, 2.0). If NULL (default),
  returns list of all q-values.

- value:

  list. Named list of diversity results per q-value.

## Value

SummarizedExperiment or list of SummarizedExperiment objects containing
diversity values keyed by q-value.

## Details

Results stored in @diversity_results with names like 'q_0.5', 'q_1.0',
etc. Use `diversity(analysis)` to get all results as a list, or
`diversity(analysis, q=1.0)` for a specific q-value.

## See also

Other TSENATAnalysis accessors:
[`divergence`](https://gallardoalba.github.io/TSENAT/reference/divergence.md),
[`lmResults`](https://gallardoalba.github.io/TSENAT/reference/lmResults.md),
[`se`](https://gallardoalba.github.io/TSENAT/reference/se.md),
[`metadata`](https://rdrr.io/pkg/S4Vectors/man/Annotated-class.html)

## Examples

``` r
# Load real TSENAT data
data(readcounts)
metadata_df <- read.table(system.file('extdata', 'metadata.tsv',
  package = 'TSENAT'), header = TRUE, sep = '\t')
gff3_file <- system.file('extdata', 'annotation.gff3.gz',
  package = 'TSENAT')
config <- tsenat_config(sample_col = 'sample', condition_col = 'condition')
analysis <- build_analysis_s4(readcounts = readcounts, tx2gene =
gff3_file,
  metadata = metadata_df, config = config, tpm = tpm,
  effective_length = effective_length)
analysis <- filter_analysis_s4(analysis, min_samples = 1, subset_n_genes
= 200)
analysis <- calculate_diversity_s4(analysis, q = 1)
#> Note: 12 genes excluded (< 75% valid values).
diversity_results <- diversity(analysis)
```
