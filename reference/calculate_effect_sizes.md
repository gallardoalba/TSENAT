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
  and RRM interaction results (from
  [`calculate_rrm()`](https://gallardoalba.github.io/TSENAT/reference/calculate_rrm.md)).

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

  Divergence SE and RRM results from analysis slots

- Enriching:

  Divergence SE with gene names via tx2gene or direct mapping

- Computing:

  Effect sizes using base `.calculate_effect_sizes()`

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

[`calculate_divergence`](https://gallardoalba.github.io/TSENAT/reference/calculate_divergence.md)
for divergence wrapper,
[`calculate_rrm`](https://gallardoalba.github.io/TSENAT/reference/calculate_rrm.md)
for RRM interaction wrapper

## Examples

``` r
# Setup: Create test analysis with divergence and RRM interaction results
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
config <- TSENAT_config(
  sample_col = 'sample',
  condition_col = 'condition',
  subject_col = 'paired_samples',
  paired = TRUE,
  control = 'normal',
  q = seq(0, 2, by = 0.5)  # Multiple q-values for RRM interaction analysis (5 unique: 0, 0.5, 1, 1.5, 2)
)
analysis <- build_analysis(readcounts = readcounts, tx2gene =
gff3_dataset, metadata = metadata_df, config = config,
  tpm = tpm, effective_length = effective_length)

analysis <- filter_analysis(analysis, stringency = 'severe')
analysis <- calculate_diversity(analysis)
analysis <- calculate_divergence(analysis)
analysis <- suppressWarnings(calculate_rrm(analysis, method = 'gam'))

# Compute effect sizes from divergence results
analysis <- calculate_effect_sizes(analysis,
  significance_threshold = 0.05)

# Access results using unified results accessor
effect_size_results <- results(analysis, type = 'effect_sizes_divergence')

# View structure of results
str(effect_size_results, max.level = 1)
#> 'data.frame':    8 obs. of  12 variables:
#>  $ gene               : chr  "ZNF493" "CENPV" "VRK2" "TLN2" ...
#>  $ p_value_interaction: num  4.69e-08 1.21e-05 2.70e-04 2.40e-03 2.82e-03 ...
#>  $ slope_diff         : num  0.0404 0.0842 0.1841 0.1512 -0.1946 ...
#>  $ effect_size_D_q0_01: num  0.0579 0.5738 0.0842 1 0.2055 ...
#>  $ effect_size_D_q0_5 : num  0.0767 0.7017 0.1084 0.9172 0.2594 ...
#>  $ effect_size_D_q1   : num  0.0787 0.6995 0.1084 0.8251 0.2565 ...
#>  $ effect_size_D_q1_5 : num  0.0686 0.6239 0.0925 0.7475 0.2181 ...
#>  $ effect_size_D_q2   : num  0.0513 0.5025 0.0682 0.6592 0.1615 ...
#>  $ per_q_pattern      : chr  "Balanced" "Balanced" "Balanced" "Rare driven" ...
#>  $ div_rare_median    : num  0.0673 0.6377 0.0963 0.9586 0.2325 ...
#>  $ div_abundant_median: num  0.0599 0.5632 0.0803 0.7034 0.1898 ...
#>  $ q_ratio            : num  1.12 1.13 1.2 1.36 1.22 ...
```
