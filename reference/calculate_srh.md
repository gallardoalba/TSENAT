# Detect q-dependent gene interactions

Wrapper around \[.calculate_srh()\] that manages TSENATAnalysis object.
Tests for genes with condition-specific q-dependent entropy patterns by
testing whether the effect of q-values DIFFERS between experimental
conditions. This detects disease-relevant or condition-specific isoform
switching patterns.

## Usage

``` r
calculate_srh(
  analysis,
  condition_col,
  output_file = NULL,
  paired = NULL,
  subject_col = NULL,
  multicorr = c("hochberg", "benjamini-yekutieli", "westfall-young", "none"),
  entropy_col = "diversity",
  q_col = "q",
  gene_col = "gene",
  wy_randomizations = 500,
  nperm_mode = c("standard", "conservative", "interactive"),
  nthreads = NULL,
  alpha = 0.05,
  p_threshold = 0.05,
  eta2_threshold_moderate = 0.01,
  eta2_threshold_strong = 0.1,
  min_nperm = 100,
  max_nperm = 10000,
  verbose = FALSE,
  ...
)
```

## Arguments

- analysis:

  `TSENATAnalysis` object. Must have diversity results from
  [`calculate_diversity()`](https://gallardoalba.github.io/TSENAT/reference/calculate_diversity.md).

- condition_col:

  `character`. Column name for sample grouping/condition (REQUIRED).
  Specifies the condition/treatment variable for testing q\\\times\\
  condition interactions. Example: 'sample_type', 'treatment',
  'disease_status'.

- output_file:

  `character` or `NULL`. Optional file path to save results. Supported
  formats: .rds (for S4 objects). Default: NULL (no file output).

- paired:

  `logical` or `NULL`. If TRUE, uses paired/blocked design (requires
  `subject_col`). If NULL, reads from `@config$paired`.

- subject_col:

  `character` or `NULL`. Column name for subject/block identifiers
  (required when `paired=TRUE`). If NULL, reads from
  `@config$subject_col`.

- multicorr:

  `character`. Multiple testing correction: 'hochberg' (default),
  'benjamini-yekutieli', 'westfall-young', or 'none'.

- entropy_col:

  `character`. Column name containing entropy/diversity data. Default:
  'diversity'.

- q_col:

  `character`. Column name containing q-values. Default: 'q'.

- gene_col:

  `character`. Column name containing gene identifiers. Default: 'gene'.

- wy_randomizations:

  `numeric` or `character`. Number of permutations for Westfall-Young
  correction. Use 'auto' to estimate from data. Default: 500.

- nperm_mode:

  `character`. Mode for automatic permutation estimation: 'standard'
  (default), 'conservative', or 'interactive'.

- nthreads:

  `numeric` or `NULL`. Number of parallel threads for computation. If
  NULL, reads from `@config$nthreads`.

- alpha:

  `numeric`. Significance level for p-value correction methods (default:
  0.05). Used by all multiple testing correction methods.

- p_threshold:

  `numeric`. P-value threshold for classification of interaction
  significance (default: 0.05).

- eta2_threshold_moderate:

  `numeric`. Effect size boundary for 'moderate' classification
  (default: 0.01).

- eta2_threshold_strong:

  `numeric`. Effect size boundary for 'strong' classification (default:
  0.10).

- min_nperm:

  `integer`. Minimum permutations for automatic estimation when
  wy_randomizations='auto' (default: 100).

- max_nperm:

  `integer`. Maximum permutations for automatic estimation when
  wy_randomizations='auto' (default: 10000).

- verbose:

  `logical`. If TRUE, prints progress messages. Default: FALSE.

- ...:

  Additional arguments passed to the base `.calculate_srh()` function.

## Value

Modified TSENATAnalysis with interaction results in @lm_results.

## Details

\## Key Features

\- \*\*Q\\\times\\ Condition Interaction\*\*: Tests if entropy patterns
across q-values differ by condition (main discovery goal) - \*\*Multi-q
Analysis\*\*: Combines diversity results for multiple q-values into a
single SummarizedExperiment for joint hypothesis testing -
\*\*Rank-Based Statistics\*\*: Scheirer-Ray-Hare test (two-way ANOVA on
ranked data, - \*\*Scheirer-Ray-Hare Test\*\*: Two-way non-parametric
ANOVA on ranks - \*\*Multiple Testing Correction\*\*: Hochberg,
Benjamini-Yekutieli, or permutation (Westfall-Young) procedures -
\*\*AR(1) Correlation Handling\*\*: Westfall-Young preserves q-value
spatial correlations (important for ordered q measurements) - \*\*Effect
Sizes\*\*: Eta-squared (\\\eta^2\\) for q\\\times\\ condition
interactions

\## Statistical Hypotheses

Tests the null hypothesis: - **H*0*** = Gene entropy q-effect does NOT
differ between conditions (q-independent) - **H*1*** = Gene entropy
q-dependence is CONDITION-SPECIFIC (interaction exists)

A significant interaction indicates condition-specific patterns in how
entropy varies across the q-value spectrum, revealing biological
processes specific to that condition.

\## Biological Example

Gene shows strong isoform switching (q-dependent entropy) in tumor cells
but NOT in healthy cells -\> Identified as disease-relevant q-dependent
gene.

For condition-specific q-dependent genes: - \*\*Condition A\*\*: Strong
entropy variation across q (q-dependent isoform usage) - \*\*Condition
B\*\*: Flat entropy profile across q (uniform isoform usage) -
\*\*Interaction\*\*: Condition-specific q-dependence pattern reveals
disease-associated splicing regulation

Analyzes how gene interactions change across q-value spectrum using
rank-based (Scheirer-Ray-Hare) or parametric (GAM) statistical tests.

\*\*Parameter resolution priority\*\* (explicit \> @config \>
default/auto-detect):

- `condition_col`: REQUIRED - must be explicitly provided

- `q`: ALWAYS auto-detected from diversity_results (all q-values tested
  together)

- `paired`: explicit arg \> `@config$paired` \> FALSE (default)

- `subject_col`: explicit arg \> `@config$subject_col`

- `multicorr`: explicit arg \> `@config$multicorr` \> 'hochberg'

- `nthreads`: explicit arg \> `@config$nthreads` \> 1 (default)

- `test`: explicit arg \> `@config$test` \> 'auto' (auto-selection)

- `nperm_mode`: explicit arg \> `@config$nperm_mode` \> 'standard'

## Examples

``` r
# Load example data (matching TSENAT.Rmd workflow)
data(readcounts)
readcounts <- as.matrix(readcounts)
mode(readcounts) <- 'numeric'
metadata_df <- read.table(
  system.file('extdata', 'metadata.tsv', package = 'TSENAT'),
  header = TRUE, sep = '\t'
)
gff3_dataset <- system.file('extdata', 'annotation.gff3.gz', package =
'TSENAT')

# Create config first (required when metadata is provided)
config <- TSENAT_config(sample_col = 'sample', condition_col = 'condition')

# Build analysis from vignette data and create manageable subset
analysis <- build_analysis(
  readcounts = readcounts,
  tx2gene = gff3_dataset,
  metadata = metadata_df,
  config = config,
  tpm = tpm,
  effective_length = effective_length
)
analysis <- filter_analysis(
  analysis,
  min_samples = 1,
  subset_n_genes = 200
)
analysis <- calculate_diversity(analysis, q = c(0.5, 1.0, 1.5))

# Test Q\eqn{\times} Condition interaction (condition_col is REQUIRED)
analysis <- calculate_srh(
  analysis,
  condition_col = 'condition',
  multicorr = 'hochberg'
)
# View results using unified accessor
rank_test_res <- results(analysis, type = 'rank_test')
if (!is.null(rank_test_res)) head(rank_test_res)
#>       gene n_q_values_tested f_statistic   p_value adj_p_value ss_interaction
#> 1 SH3PXD2A                 3  0.10283830 0.9024993           1     0.35754547
#> 2     ATG5                 3  0.31874647 0.7288031           1     0.41779676
#> 3   FAXDC2                 3  0.07408088 0.9287176           1     0.53373469
#> 4     CYBB                 3  0.19803026 0.8211066           1     0.33734628
#> 5   PICALM                 3  0.07600418 0.9269395           1     0.04742451
#> 6  CCDC149                 3  0.09984740 0.9051897           1     0.45310777
#>   ss_residual df_interaction df_residual effect_size_eta2 interaction_class
#> 1   0.1850105              2          45        0.6590021   Robust across q
#> 2   0.4288228              2          45        0.4934882   Robust across q
#> 3   0.6214901              2          45        0.4620180   Robust across q
#> 4   0.4670058              2          45        0.4194013   Robust across q
#> 5   0.1261091              2          45        0.2732872   Robust across q
#> 6   1.2613314              2          45        0.2642892   Robust across q
#>    test_method heteroscedastic boundary_clustered highly_skewed
#> 1 srh_unpaired           FALSE              FALSE         FALSE
#> 2 srh_unpaired           FALSE              FALSE         FALSE
#> 3 srh_unpaired           FALSE              FALSE         FALSE
#> 4 srh_unpaired           FALSE              FALSE         FALSE
#> 5 srh_unpaired           FALSE              FALSE         FALSE
#> 6 srh_unpaired           FALSE              FALSE         FALSE
```
