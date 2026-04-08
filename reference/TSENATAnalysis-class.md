# TSENATAnalysis S4 Class Definition

Central container for unified analysis workflows in TSENAT.

Central container for unified analysis workflows in TSENAT.

## Value

An S4 object of class `TSENATAnalysis` containing:

- Raw data (SummarizedExperiment)

- Analysis configuration

- Results from diversity, linear model, jackknife, and divergence
  analyses

- Cached visualization objects

- Reproducibility metadata and function history

## Details

The `TSENATAnalysis` class encapsulates all components of a complete
TSENAT analysis: raw data, configuration metadata, and results from each
analytical step. This unified object ensures metadata is never lost
through the analysis pipeline and provides consistent accessor methods
for result retrieval.

Access results via accessor methods (recommended): `diversity(obj, q)`
for diversity, `lmResults(obj)` for models, `jeoResults(obj, q)` for
entropy outlier jackknife, `jisResults(obj, q)` for isoform switching
jackknife, `getMeta(obj)` for metadata.

The `TSENATAnalysis` class encapsulates all components of a complete
TSENAT analysis: raw data, configuration metadata, and results from each
analytical step. This unified object ensures metadata is never lost
through the analysis pipeline and provides consistent accessor methods
for result retrieval.

## Slots

- `se`:

  `SummarizedExperiment`. The base expression data object (genes x
  samples) with assays and colData.

- `config`:

  `list`. Configuration metadata specifying analysis parameters that
  persist through the workflow (q-values, sample grouping columns,
  etc.). Set once via
  [`tsenat_config()`](https://gallardoalba.github.io/TSENAT/reference/tsenat_config.md)
  and used by all downstream wrapper functions.

- `diversity_results`:

  `list`. Named list of diversity calculation results. Each name
  corresponds to a q-value (e.g., 'q_0.5', 'q_1.0'). Values are
  SummarizedExperiment objects or data.frames containing entropy values
  for each gene at that q-value.

- `lm_results`:

  `list`. Complex results from linear model and statistical testing.
  Top-level names identify analysis type:

  `lm_interaction`

  :   LM/GAM/GEE model results (list with `$results` data.frame,
      `$models` list, etc.)

  `q_interactions`

  :   Friedman/rank-based test results

  `divergence_difference`

  :   Differential divergence comparison

- `pairwise_results`:

  `list`. Pairwise group comparison results. Contains differential
  testing output computed by
  [`calculate_difference_s4()`](https://gallardoalba.github.io/TSENAT/reference/calculate_difference_s4.md),
  stored under the `difference` component.

- `jackknife_results`:

  `list`. Resampling-based confidence intervals. Names correspond to
  q-values (e.g., 'q_0.5', 'q_1.0'). Values are jackknife result objects
  containing resamples, CI bounds, and diagnostics.

- `divergence_results`:

  `list`. Divergence metric calculations. Typically contains:

  `tsallis_divergence`

  :   SummarizedExperiment with divergence values

  `effect_sizes`

  :   data.frame with Cohen's d, etc.

- `plots`:

  `list`. Cached visualization objects (ggplot). Names identify plot
  type (e.g., 'q_curve', 'lm_interaction', 'influence'). Populated by
  [`tsenat()`](https://gallardoalba.github.io/TSENAT/reference/tsenat.md)
  if `generate_plots=TRUE`.

- `metadata`:

  `list`. Reproducibility and tracking metadata. Automatically
  maintained by wrapper functions. Includes:

  `created_at`

  :   Timestamp of object creation

  `function_calls`

  :   Vector of wrapper functions called

  `function_timestamps`

  :   Timestamps for each function call

  `package_version`

  :   TSENAT version at creation

- `se`:

  `SummarizedExperiment`. The base expression data object (genes x
  samples) with assays and colData.

- `config`:

  `list`. Configuration metadata specifying analysis parameters that
  persist through the workflow (q-values, sample grouping columns,
  etc.). Set once via
  [`tsenat_config()`](https://gallardoalba.github.io/TSENAT/reference/tsenat_config.md)
  and used by all downstream wrapper functions.

- `diversity_results`:

  `list`. Named list of diversity calculation results. Each name
  corresponds to a q-value (e.g., 'q_0.5', 'q_1.0'). Values are
  SummarizedExperiment objects or data.frames containing entropy values
  for each gene at that q-value.

- `lm_results`:

  `list`. Complex results from linear model and statistical testing.
  Top-level names identify analysis type:

  `lm_interaction`

  :   LM/GAM/GEE model results (list with `$results` data.frame,
      `$models` list, etc.)

  `q_interactions`

  :   Friedman/rank-based test results

  `divergence_difference`

  :   Differential divergence comparison

- `jackknife_results`:

  `list`. Resampling-based confidence intervals. Names correspond to
  q-values (e.g., 'q_0.5', 'q_1.0'). Values are jackknife result objects
  containing resamples, CI bounds, and diagnostics.

- `divergence_results`:

  `list`. Divergence metric calculations. Typically contains:

  `tsallis_divergence`

  :   SummarizedExperiment with divergence values

  `effect_sizes`

  :   data.frame with Cohen's d, etc.

- `plots`:

  `list`. Cached visualization objects (ggplot). Names identify plot
  type (e.g., 'q_curve', 'lm_interaction', 'influence'). Populated by
  [`tsenat()`](https://gallardoalba.github.io/TSENAT/reference/tsenat.md)
  if `generate_plots=TRUE`.

- `metadata`:

  `list`. Reproducibility and tracking metadata. Automatically
  maintained by wrapper functions. Includes:

  `created_at`

  :   Timestamp of object creation

  `function_calls`

  :   Vector of wrapper functions called

  `function_timestamps`

  :   Timestamps for each function call

  `package_version`

  :   TSENAT version at creation

## Accessor Methods

- `diversity(object, q=NULL)`:

  Extract diversity results for q-value

- `lmResults(object, component=NULL)`:

  Extract LM results

- `jeoResults(object, q=NULL)`:

  Extract jackknife entropy outlier results

- `jisResults(object, q=NULL)`:

  Extract jackknife isoform switching results

- `divergence(object)`:

  Extract divergence results

- `getPlot(object, type=NULL)`:

  Retrieve cached plot

- `addPlot(object, type, plot)`:

  Add/cache a new plot

- `show(object)`:

  Display object summary

- `summary(object)`:

  Get detailed analysis summary

## Validation

Validity is checked at object construction. Ensures @se is a
SummarizedExperiment and all slots are correct types.

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
