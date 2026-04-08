# Detect q-dependent gene interactions

Detect q-dependent gene interactions

## Usage

``` r
rank_test_q_condition_s4(
  analysis,
  condition_col,
  q = NULL,
  output_file = NULL,
  paired = NULL,
  subject_col = NULL,
  test = c("auto", "kruskal-wallis", "friedman", "art"),
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

  `TSENATAnalysis` object.

- condition_col:

  `character`. Column name for sample grouping/condition (REQUIRED).
  Specifies the condition/treatment variable for testing q×condition
  interactions. Example: 'sample_type', 'treatment', 'disease_status'.

- q:

  `numeric` or `NULL`. Q-values to test across spectrum. If NULL,
  auto-detects from `@config$q_values` or diversity results.

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

- test:

  `character`. Test method: 'auto' (default), 'kruskal-wallis'
  (unpaired), 'friedman' (paired), or 'art' (aligned rank transform).

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

  Additional arguments passed to the base `.rank_test_q_condition()`
  function.

## Value

Modified TSENATAnalysis with interaction results in @lm_results.

## Details

Analyzes how gene interactions change across q-value spectrum using
rank-based (Friedman/Kruskal-Wallis) or parametric (GAM) statistical
tests.

\*\*Parameter resolution priority\*\* (explicit \> @config \>
default/auto-detect):

- `condition_col`: REQUIRED - must be explicitly provided

- `q`: explicit arg \> `@config$q_values` \> extract from
  diversity_results keys

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
config <- tsenat_config(sample_col = 'sample', condition_col = 'condition')

# Build analysis from vignette data and create manageable subset
analysis <- build_analysis_s4(readcounts = readcounts, tx2gene =
gff3_dataset, metadata = metadata_df, config = config,
  tpm = tpm, effective_length = effective_length)
analysis <- filter_analysis_s4(analysis, min_samples = 1, subset_n_genes
= 200)
analysis <- calculate_diversity_s4(analysis, q = c(0.5, 1.0, 1.5))
#> Note: 12 genes excluded (< 75% valid values).

# Test Q×Condition interaction (condition_col is REQUIRED)
analysis <- rank_test_q_condition_s4(analysis, condition_col = 'condition', 
                                           multicorr = 'hochberg')
# View results
results <- lmResults(analysis)
#> Warning: No LM interaction results found. Run calculate_lm_interaction_s4() first.
if (!is.null(results)) head(results$q_interactions)
```
