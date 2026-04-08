# Compute Effect Sizes from Divergence Results (S4 Wrapper)

S4 wrapper for `. effect_sizes_divergence()` that extracts divergence
and LM results directly from a TSENATAnalysis object.

## Usage

``` r
effect_sizes_divergence_s4(
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
  [`calculate_divergence_s4()`](https://gallardoalba.github.io/TSENAT/reference/calculate_divergence_s4.md))
  and LM interaction results (from
  [`calculate_lm_interaction_s4()`](https://gallardoalba.github.io/TSENAT/reference/calculate_lm_interaction_s4.md)).

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

  Divergence SE and LM results from analysis slots

- Enriching:

  Divergence SE with gene names via tx2gene or direct mapping

- Computing:

  Effect sizes using base `.effect_sizes_divergence()`

- Storing:

  Results in metadata with function call tracking

\*\*Parameter resolution priority\*\* (explicit \> metadata \> default):

- `significance_threshold`: Uses explicit arg, else
  `metadata(analysis)$significance_threshold`, else 0.05

- `enrich_per_q_pattern`: Uses explicit arg, else
  `metadata(analysis)$enrich_per_q_pattern`, else TRUE

- `verbose`: Uses explicit arg, else `metadata(analysis)$verbose`, else
  TRUE

Results are accessed via: `metadata(analysis)$effect_sizes_divergence`

## See also

[`calculate_divergence_s4`](https://gallardoalba.github.io/TSENAT/reference/calculate_divergence_s4.md)
for divergence wrapper,
[`calculate_lm_interaction_s4`](https://gallardoalba.github.io/TSENAT/reference/calculate_lm_interaction_s4.md)
for LM interaction wrapper

## Examples

``` r
# Setup: Create test analysis with divergence and LM interaction results
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
config <- tsenat_config(
  sample_col = 'sample',
  condition_col = 'condition',
  subject_col = 'paired_samples',
  paired = TRUE,
  control = 'normal',
  q_values = seq(0, 2, by = 0.05)  # Recommended: 41 q-values for paired designs
)
analysis <- build_analysis_s4(readcounts = readcounts, tx2gene =
gff3_dataset, metadata = metadata_df, config = config,
  tpm = tpm, effective_length = effective_length)

analysis <- filter_analysis_s4(analysis, stringency = 'severe')
analysis <- calculate_diversity_s4(analysis, q = c(0.5, 1.0, 1.5, 2.0, 2.5))
#> Note: 7 genes excluded (< 75% valid values).
analysis <- calculate_divergence_s4(analysis, q = c(0.5, 1.0, 1.5, 2.0, 2.5))
analysis <- calculate_lm_interaction_s4(analysis, method = 'gam')
#> Warning: nlminb problem, convergence error code = 1
#>   message = singular convergence (7)
#> Warning: nlminb problem, convergence error code = 1
#>   message = iteration limit reached without convergence (10)
#> Warning: nlminb problem, convergence error code = 1
#>   message = iteration limit reached without convergence (10)

# Compute effect sizes from divergence results
analysis <- effect_sizes_divergence_s4(analysis,
  significance_threshold = 0.05)

# Access results using metadata accessor
effect_size_results <- getMeta(analysis, 'effect_sizes_divergence')

# View structure of results
str(effect_size_results, max.level = 1)
#> List of 2
#>  $ interaction_results:'data.frame': 5 obs. of  22 variables:
#>  $ validation_stats   :List of 5
```
