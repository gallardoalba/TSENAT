# Calculate SAIT interactions and store in TSENATAnalysis

Statistical interaction testing for entropy diversity using flexible
scale-adaptive models (GAM, LMM, GEE, FPCA) with AR(1) correlation
structure for repeated measures. Tests for significant q-by-condition
interactions across entropic indices.

## Usage

``` r
calculate_sait(
  analysis,
  fdr_threshold = NULL,
  formula = NULL,
  condition_col = NULL,
  method = "gam",
  paired = NULL,
  subject_col = NULL,
  nthreads = NULL,
  multicorr = NULL,
  corstr = NULL,
  pcorr = NULL,
  verbose = NULL,
  return_model_data = NULL,
  output_file = NULL,
  ...
)
```

## Arguments

- analysis:

  `TSENATAnalysis` object.

- fdr_threshold:

  `numeric`. FDR cutoff for significance. Default: 0.05.

- formula:

  `formula` or NULL. Reserved for future use.

- condition_col:

  `character` or `NULL`. Column name in colData identifying sample
  conditions. If NULL, reads from `@config$condition_col` or
  auto-detects common column names.

- method:

  `character`. Statistical method (e. g. , 'lmm', 'gam', 'gee'). If
  NULL, uses method from @config\$method or defaults to 'lmm'. For
  paired designs, `'gam'` dispatches to GAMM:
  [`nlme::lme`](https://rdrr.io/pkg/nlme/man/lme.html) with regression
  splines (`ns(q, df = 3) x condition`), subject random intercept and
  AR(1) within subject x condition (marginal F-test for the
  interaction); [`mgcv::gamm()`](https://rdrr.io/pkg/mgcv/man/gamm.html)
  and standard GAM are fallbacks.

- paired:

  `logical`. Whether to use paired design. Default: FALSE. If not
  specified, reads from `@config$paired` if available.

- subject_col:

  `character` or `NULL`. Column name identifying subject IDs for paired
  designs. If NULL, reads from `@config$subject_col` if available.

- nthreads:

  `numeric` or `NULL`. Number of CPU threads for parallel processing. If
  NULL, reads from `@config$nthreads` (or defaults to NULL, letting base
  function decide).

- multicorr:

  `character` or `NULL`. Multiple comparison correction method. Options:
  'hochberg', 'westfall-young', 'benjamini-yekutieli'. If NULL, uses
  method from @config or base function defaults.

- corstr:

  `character` or `NULL`. Correlation structure for GEE models. Options:
  'ar1', 'exchangeable', 'independence', 'auto'. 'auto' selects the best
  structure via QIC. If NULL, uses method from @config or base function
  defaults.

- pcorr:

  `character` or `NULL`. P-value correction method. Default: 'BH'
  (Benjamini-Hochberg). If NULL, reads from `@config$pcorr` if
  available.

- verbose:

  `logical`. Print progress messages. Default: FALSE.

- return_model_data:

  `logical`. Return model data for visualization. Default: TRUE.

- output_file:

  `character` or `NULL`. Optional file path to save results. Supported
  formats: .rds (for S4 objects), .tsv, .csv, .txt (for tables).
  Default: NULL (no file output).

- ...:

  Additional arguments passed to the base LM function, including:
  pvalue, min_obs, assay_name, bias_correction, regularization, storey,
  wy_randomizations, adaptive_knots, block_col, strata_col,
  permutation_scheme (exchangeability-compatible Westfall-Young
  permutations), etc.

## Value

Modified TSENATAnalysis with results in @sait_results\$sait_interaction.

## Details

Extracts diversity results from @diversity_results (prerequisite),
combines across q-values into single SummarizedExperiment, then runs
`.calculate_sait()`.

\*\*Parameter Priority Resolution:\*\*

- `nthreads`: Priority: explicit \> @config \> NULL

Parameters are resolved in priority order: 1. Explicit arguments passed
to function 2. Values from analysis@config (if present) 3. Function
defaults

## Examples

``` r
if (FALSE) { # \dontrun{
# Create test analysis with appropriate sample structure for paired design
# Note: requires lme4 package for LMM fitting; uses synthetic data
set.seed(42)

# Create transcript-level counts with biological signal
# Note: Use adequate complexity (transcripts/genes, samples, expression) 
# to avoid filtering away all genes during diversity computation
n_genes <- 50
n_transcripts_per_gene <- 30
n_transcripts <- n_genes * n_transcripts_per_gene
n_samples <- 16  # 8 subjects x 2 conditions (paired design)

# Generate counts with clear biological signal
control_idx <- seq(1, n_samples, by = 2)
treatment_idx <- seq(2, n_samples, by = 2)

counts <- matrix(0, nrow = n_transcripts, ncol = n_samples)
for (j in seq_len(n_samples)) {
  lambda <- if (j %in% control_idx) 100 else 180
  counts[, j] <- rpois(n_transcripts, lambda = lambda)
}
counts <- pmax(counts, 50)  # Ensure minimum expression

rownames(counts) <- paste0('TX_', seq_len(n_transcripts))
colnames(counts) <- paste0('Sample_', seq_len(n_samples))

# Create rowData with gene mapping (tx2gene structure)
rowdata <- data.frame(
  transcript_id = rownames(counts),
  gene_id = rep(paste0('GENE_', 1:n_genes), 
                each = n_transcripts_per_gene),
  row.names = rownames(counts)
)

# Create colData with paired design metadata
coldata <- data.frame(
  sample_id = colnames(counts),
  condition = rep(c('control', 'treatment'), 
                  length.out = n_samples),
  subject = rep(paste0('Subject_', 1:8), 
                length.out = n_samples),
  row.names = colnames(counts)
)

# Build SummarizedExperiment
se <- SummarizedExperiment::SummarizedExperiment(
  assays = list(counts = counts),
  rowData = S4Vectors::DataFrame(rowdata),
  colData = S4Vectors::DataFrame(coldata)
)

# Add tx2gene metadata for gene-level aggregation
S4Vectors::metadata(se)$tx2gene <- 
  data.frame(Transcript = rowdata$transcript_id,
             Gene = rowdata$gene_id)

# Initialize TSENATAnalysis
analysis <- TSENATAnalysis(se = se, config = list())

# Compute diversity (prerequisite for SAIT interaction analysis)
analysis <- calculate_diversity(
  analysis, 
  q = c(0.5, 1.0, 1.5, 2.0, 2.5)
)

# Calculate q x condition interactions using GAM
analysis <- suppressWarnings(calculate_sait(
  analysis,
  condition_col = 'condition',
  method = 'gam'
))

# View top interaction results using unified accessor (first 3 genes)
res <- results(analysis, type = "sait")
if (!is.null(res)) head(res, 3)
} # }
```
