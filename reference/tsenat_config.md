# Create and return TSENAT configuration

Builds a configuration list for use with
[`tsenat`](https://gallardoalba.github.io/TSENAT/reference/tsenat.md)().
Allows specifying analysis parameters once and reusing across multiple
analyses.

## Usage

``` r
tsenat_config(
  q_values = NULL,
  condition_col = "condition",
  subject_col = NULL,
  sample_col = "sample",
  paired = FALSE,
  control = NULL,
  p_threshold = 0.05,
  fdr_threshold = 0.05,
  significance_threshold = 0.05,
  nboot = 1000,
  bootstrap_method = "percentile",
  stringency = "medium",
  nthreads = 1,
  ...
)
```

## Arguments

- q_values:

  `numeric`. Q-values for Tsallis entropy spectrum. Default:
  `seq(0, 2, by = 0.5)`.

- condition_col:

  `character`. Name of column in `colData(se)` containing experimental
  conditions/groups. Default: 'condition'.

- subject_col:

  `character`. Name of column in `colData(se)` containing subject/sample
  identifiers for paired/repeated designs. If provided, enables paired
  analysis. Default: NULL (unpaired).

- sample_col:

  `character`. Name of column in `colData(se)` containing sample
  identifiers. Default: 'sample'.

- paired:

  `logical`. Whether samples are paired/repeated measures. Default:
  FALSE. Used by jackknife and difference analysis.

- control:

  `character`. Reference/control group label for difference analysis
  (e.g., 'control', 'wt'). Only used if 'difference' in methods.
  Default: NULL.

- p_threshold:

  `numeric`. Raw p-value threshold for significance in LM interaction
  testing. Default: 0.05.

- fdr_threshold:

  `numeric`. Adjusted p-value (FDR/Benjamini-Hochberg) threshold.
  Default: 0.05.

- significance_threshold:

  `numeric`. Significance cutoff for effect sizes, assumptions testing,
  and result filtering. Default: 0.05.

- nboot:

  `integer`. Number of bootstrap resamples for jackknife confidence
  intervals. Default: 1000.

- bootstrap_method:

  `character`. Bootstrap CI method: 'percentile' (fast, assumes
  symmetric distribution) or 'bca' (bias-corrected, better for skewed
  data like bounded entropy). Default: 'percentile'.

- stringency:

  `character`. Transcript filtering stringency level. Options: 'lenient'
  (minimal filtering), 'medium' (default, reasonable filtering),
  'severe' (strict filtering, recommended for high-confidence results).
  Controls which transcripts/genes are retained in initial filtering
  step. Default: 'medium'.

- nthreads:

  `integer`. Number of threads for parallel computation where supported
  (diversity, divergence, LM fitting). Default: 1 (no parallelization).
  Use 2+ for multi-core systems to improve performance.

- ...:

  Additional configuration parameters (stored as-is in @config slot).
  Examples: `q_diff=1.0` (specific q for differences), `alpha=0.05`
  (significance for assumptions), etc.

## Value

`list` with class `TSENATConfig` containing all specified parameters.

## Details

Configuration is stored in the TSENATAnalysis@config slot and used by
wrapper functions to configure analysis behavior.

## Examples

``` r
# Default config with standard parameters
cfg <- tsenat_config()

# Custom with paired analysis and strict filtering for high-confidence results
cfg <- tsenat_config(
  q_values = seq(0, 2, by = 0.05),  # Recommended for paired designs: 41 values
  condition_col = 'treatment',
  subject_col = 'subject_id',
  paired = TRUE,
  control = 'untreated',
  stringency = 'severe',
  bootstrap_method = 'bca',
  nboot = 5000,
  significance_threshold = 0.01
)
```
