# Compute Effect Sizes from Divergence Results (S4 Wrapper)

S4 wrapper for `. calculate_effect_sizes()` that extracts divergence and
LM results directly from a TSENATAnalysis object.

## Usage

``` r
calculate_effect_sizes(
  analysis,
  significance_threshold = NULL,
  enrich_per_q_pattern = NULL,
  verbose = NULL,
  output_file = NULL,
  ...
)
```

## Arguments

- analysis:

  `TSENATAnalysis`. An S4 object containing divergence results (from
  [`calculate_divergence()`](https://gallardoalba.github.io/TSENAT/reference/calculate_divergence.md))
  and SAIT interaction results (from
  [`calculate_sait()`](https://gallardoalba.github.io/TSENAT/reference/calculate_sait.md)).

- significance_threshold:

  `numeric`. Adjusted p-value threshold for filtering significant genes
  (default: 0.05).

- enrich_per_q_pattern:

  `logical`. If TRUE, enriches results with per-q divergence patterns
  (default: TRUE).

- verbose:

  `logical`. If TRUE, print diagnostic messages (default: TRUE).

- output_file:

  `character` or `NULL`. Optional file path to save results. Supported
  formats: .rds (for S4 objects). Default: NULL (no file output).

- ...:

  Additional arguments passed to the base function.

## Value

Modified TSENATAnalysis with effect size results stored via
`metadata(analysis)$effect_sizes_divergence`. Returns the analysis
object visibly to support piping and method chaining.

## Details

\*\*Workflow steps:\*\*

- Validating:

  Input analysis object and required results

- Extracting:

  Divergence SE and SAIT results from analysis slots

- Enriching:

  Divergence SE with gene names via tx2gene or direct mapping

- Computing:

  Effect sizes using base `.calculate_effect_sizes()`

- Storing:

  Results in metadata with function call tracking

\*\*Parameter resolution priority\*\* (explicit \> config \> default):

- `significance_threshold`: Uses explicit arg, else
  `analysis@config$significance_threshold`, else 0.05

- `enrich_per_q_pattern`: Uses explicit arg, else
  `analysis@config$enrich_per_q_pattern`, else TRUE

- `verbose`: Uses explicit arg, else `analysis@config$verbose`, else
  TRUE

Results are accessed via: `metadata(analysis)$effect_sizes_divergence`

## See also

[`calculate_divergence`](https://gallardoalba.github.io/TSENAT/reference/calculate_divergence.md)
for divergence wrapper,
[`calculate_sait`](https://gallardoalba.github.io/TSENAT/reference/calculate_sait.md)
for SAIT interaction wrapper

## Examples

``` r
# Setup: Create test analysis with divergence and SAIT interaction results
data(readcounts)
readcounts <- as.matrix(readcounts)
mode(readcounts) <- 'numeric'
metadata_df <- read.table(
  system.file('extdata', 'metadata.tsv', package = 'TSENAT'),
  header = TRUE, sep = '\t'
)
gff3_dataset <- system.file('extdata', 'annotation.gff3.gz', package =
'TSENAT')

# Configure analysis parameters (best practice for reproducibility)
if (FALSE) { # \dontrun{
config <- TSENAT_config(
  sample_col = 'sample',
  condition_col = 'condition',
  subject_col = 'paired_samples',
  paired = TRUE,
  control = 'normal',
  q = seq(0, 2, by = 0.5)  # Multiple q-values for SAIT (5 unique values)
)
analysis <- build_analysis(readcounts = readcounts, tx2gene =
gff3_dataset, metadata = metadata_df, config = config,
  tpm = tpm, effective_length = effective_length)

analysis <- filter_analysis(analysis, stringency = 'severe')
analysis <- calculate_diversity(analysis)
analysis <- calculate_divergence(analysis)
analysis <- suppressWarnings(calculate_sait(analysis, method = 'gam'))

# Compute effect sizes from divergence results
analysis <- calculate_effect_sizes(analysis,
  significance_threshold = 0.05)

# Access results using unified results accessor
effect_size_results <- results(analysis, type = 'effect_sizes_divergence')

# View structure of results
str(effect_size_results, max.level = 1)
} # }
```
