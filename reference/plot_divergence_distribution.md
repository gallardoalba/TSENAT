# Plot Tsallis Divergence Effect Size Distribution (S4 Wrapper)

S4 wrapper that extracts effect size results from a TSENATAnalysis
object and generates a histogram visualization of Tsallis divergence
effect sizes across genes.

## Usage

``` r
plot_divergence_distribution(
  analysis,
  threshold = 0.1,
  output_file = NULL,
  width = 12,
  height = 6,
  verbose = FALSE,
  ...
)
```

## Arguments

- analysis:

  `TSENATAnalysis` object with effect sizes computed (typically via
  [`calculate_effect_sizes`](https://gallardoalba.github.io/TSENAT/reference/calculate_effect_sizes.md)).

- threshold:

  `numeric`. Effect size threshold for visual marking in the plot.
  Default is 0.1 (information-theoretic significance level).

- output_file:

  `character`. Optional file path to save the plot. If NULL, plot is
  returned but not saved.

- width:

  `numeric`. Plot width in inches. Default is 10.

- height:

  `numeric`. Plot height in inches. Default is 6.

- verbose:

  `logical`. Print status messages. Default is TRUE.

- ...:

  Additional arguments passed to the underlying plotting function.

## Value

Invisibly returns the file path if saved, otherwise the ggplot object.
If the plot cannot be created (missing data, ggplot2 not available),
returns NULL invisibly with an informative message.

## Details

This wrapper extracts the interaction results (with effect size columns)
from `analysis@metadata$effect_sizes_divergence$interaction_results` and
passes them to the base `.plot_divergence_distribution()` function.

The function visualizes the distribution of effect sizes using the
median q-value's divergence (typically around q=1.0, close to Shannon
entropy). A red dashed line marks the information-theoretic significance
threshold.

\*\*Data Requirements:\*\*

- Effect sizes must be computed via
  [`calculate_effect_sizes()`](https://gallardoalba.github.io/TSENAT/reference/calculate_effect_sizes.md)

- `@metadata$effect_sizes_divergence$interaction_results` must contain
  columns matching pattern `effect_size_D_q*`

## See also

[`calculate_effect_sizes`](https://gallardoalba.github.io/TSENAT/reference/calculate_effect_sizes.md)
for computing effect sizes.

## Examples

``` r
# Plot 2: Distribution of effect sizes across genes
data(readcounts)
readcounts <- as.matrix(readcounts)
mode(readcounts) <- 'numeric'
metadata_df <- read.table(
  system.file('extdata', 'metadata.tsv', package = 'TSENAT'),
  header = TRUE, sep = '\t'
)
gff3_dataset <- system.file('extdata', 'annotation.gff3.gz', package =
'TSENAT')
# Configure analysis parameters first
config <- TSENAT_config(
  sample_col = 'sample',
  condition_col = 'condition',
  control = 'normal'
)

# Build analysis with configured parameters
analysis <- build_analysis(
  readcounts = readcounts,
  tx2gene = gff3_dataset,
  metadata = metadata_df,
  config = config,
  tpm = tpm,
  effective_length = effective_length
)

analysis <- filter_analysis(analysis, stringency = 'severe')
analysis <- calculate_diversity(analysis, q = c(0.5, 1.0, 1.5, 2.0, 2.5), verbose
= FALSE)
analysis <- calculate_divergence(analysis, q = c(0.5, 1.0, 1.5, 2.0, 2.5), verbose
= FALSE)
analysis <- calculate_lm(analysis, method = 'gam')
analysis <- calculate_effect_sizes(analysis)
p_dist <- plot_divergence_distribution(analysis)
print(p_dist)

```
