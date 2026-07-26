# Detect q-dependent gene interactions

Wrapper around \[.calculate_rank_transform()\] that manages
TSENATAnalysis object. Tests for genes with condition-specific
q-dependent entropy patterns by testing whether the effect of q-values
DIFFERS between experimental conditions. This detects disease-relevant
or condition-specific isoform switching patterns.

## Usage

``` r
calculate_rank_transform(
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
  method = c("art", "rt"),
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

- method:

  `character`. Non-parametric test method:

  - `'art'` (default): Aligned Rank Transform via ARTool package.
    State-of-the-art for non-parametric interaction testing. Properly
    handles factorial interactions by stripping main effects before
    ranking (Higgins & Tashtoush 1994, Wobbrock et al. 2011).

  - `'rt'`: Conover-Iman Rank Transform. Legacy non-parametric
    procedure. Valid for main effects but has known limitations for
    interaction testing (Conover & Iman 1981). ART is strongly
    recommended.

- ...:

  Additional arguments passed to the base `.calculate_rank_transform()`
  function.

## Value

Modified TSENATAnalysis with interaction results in @rank_test_results.

## Details

\## Key Features

\- \*\*Q\\\times\\ Condition Interaction\*\*: Tests if entropy patterns
across q-values differ by condition (main discovery goal) - \*\*Multi-q
Analysis\*\*: Combines diversity results for multiple q-values into a
single SummarizedExperiment for joint hypothesis testing - \*\*Aligned
Rank Transform (ART)\*\*: State-of-the-art non-parametric interaction
testing via ARTool package (default). Strips main effects before ranking
to preserve interaction structure (Higgins & Tashtoush 1994, Wobbrock et
al. 2011). - \*\*Conover-Iman Rank Transform\*\*: Two-way non-parametric
ANOVA on ranks (fallback via \`method='rt'\`) - \*\*Multiple Testing
Correction\*\*: Hochberg, Benjamini-Yekutieli, or permutation
(Westfall-Young) procedures - \*\*AR(1) Correlation Handling\*\*:
Westfall-Young preserves q-value spatial correlations (important for
ordered q measurements) - \*\*Effect Sizes\*\*: Eta-squared (\\\eta^2\\)
for q\\\times\\ condition interactions

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
rank-based (Conover-Iman Rank Transform) statistical tests.

\*\*Parameter resolution priority\*\* (explicit \> @config \>
default/auto-detect):

- `condition_col`: REQUIRED - must be explicitly provided

- `q`: ALWAYS auto-detected from diversity_results (all q-values tested
  together)

- `paired`: explicit arg \> `@config$paired` \> FALSE (default)

- `subject_col`: explicit arg \> `@config$subject_col`

- `multicorr`: explicit arg \> `@config$multicorr` \> 'hochberg'

- `nthreads`: explicit arg \> `@config$nthreads` \> 1 (default)

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
analysis <- calculate_rank_transform(
  analysis,
  condition_col = 'condition',
  multicorr = 'hochberg'
)
# View results using unified accessor
rank_test_res <- results(analysis, type = 'rank_test')
if (!is.null(rank_test_res)) head(rank_test_res)
#>       gene n_q_values_tested f_statistic   p_value adj_p_value ss_interaction
#> 1    FBXO3                 3  0.24069643 0.7871571   0.9981785     101.791667
#> 2 FAM114A2                 3  0.50372281 0.6078809   0.9981785     215.041667
#> 3 TMEM183A                 3  0.09054441 0.9136116   0.9981785      39.500000
#> 4   ZNF493                 3  0.04868521 0.9525346   0.9981785      21.291667
#> 5    HDAC2                 3  0.51866394 0.5990786   0.9981785     220.166667
#> 6     CYBB                 3  0.01900408 0.9811838   0.9981785       8.041667
#>   ss_residual df_interaction df_residual effect_size_eta2 interaction_class
#> 1     8881.00              2          45      0.012441534   Robust across q
#> 2     8965.00              2          45      0.011374346   Robust across q
#> 3     9161.25              2          45      0.005771208   Robust across q
#> 4     9184.00              2          45      0.004486701   Robust across q
#> 5     8914.25              2          45      0.003787202   Robust across q
#> 6     8886.25              2          45      0.003673146   Robust across q
#>    test_method heteroscedastic boundary_clustered highly_skewed
#> 1 art_unpaired           FALSE              FALSE         FALSE
#> 2 art_unpaired           FALSE              FALSE         FALSE
#> 3 art_unpaired           FALSE              FALSE         FALSE
#> 4 art_unpaired           FALSE              FALSE         FALSE
#> 5 art_unpaired           FALSE              FALSE         FALSE
#> 6 art_unpaired           FALSE              FALSE         FALSE
```
