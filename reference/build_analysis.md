# Build a Complete TSENATAnalysis Object

Convenience wrapper that combines `. build_se()` and
[`TSENATAnalysis()`](https://gallardoalba.github.io/TSENAT/reference/TSENATAnalysis.md)
into a single function call. This creates a complete analysis object
ready for Tsallis entropy computation and downstream analysis.

## Usage

``` r
build_analysis(
  readcounts = NULL,
  salmon_dir = NULL,
  tx2gene,
  assay_name = "counts",
  metadata = NULL,
  tpm = NULL,
  effective_length = NULL,
  config = list(),
  skip = FALSE,
  verbose = FALSE
)
```

## Arguments

- readcounts:

  A matrix or data.frame of transcript-level read counts with transcript
  IDs as row names and sample names as column names. Typically output
  from quantification tools (SALMON, kallisto, etc. ). Optional when
  `salmon_dir` is provided (in which case readcounts are auto-loaded
  from quant.sf files).

- salmon_dir:

  Optional character path to directory containing Salmon quantification
  output. Expected structure: `salmon_dir/sample_name/quant. sf`. When
  provided, automatically discovers and reads all quant. sf files. If
  both `readcounts` and `salmon_dir` are provided, `salmon_dir` takes
  precedence. Default: `NULL`.

- tx2gene:

  Either: - A path to a GFF3 or GFF3.gz file containing
  transcript-to-gene mapping - A path to a TSV file with columns
  'Transcript' and 'Gene' - A data.frame with transcript-to-gene mapping

- assay_name:

  Character. Name for the assay (default: 'counts').

- metadata:

  Optional data.frame with sample metadata. Should have sample names as
  row names and metadata columns (e.g., sample_type, condition, etc.).
  If NULL, will attempt to read from `config$metadata`. Priority:
  explicit `metadata` argument \> `config$metadata` \> NULL.

- tpm:

  REQUIRED matrix of transcript-level TPM values (Transcripts Per
  Million). Must be provided and will be stored in the
  SummarizedExperiment for use during diverse filtering. Same dimensions
  as readcounts required (rows = transcripts, columns = samples).
  Typically from SALMON quantification output. If not provided to
  build_analysis(), subsequent filter_analysis() calls will fail with an
  explicit error message.

- effective_length:

  REQUIRED numeric vector of transcript effective lengths (e.g., from
  SALMON EffectiveLength column). Length must match nrow(readcounts).
  Typically obtained as the median effective length across all samples.
  If not provided to build_analysis(), the data will not be stored for
  later use in length-normalized calculations.

- config:

  Optional list of configuration parameters to store in the
  TSENATAnalysis object. Can also contain `config$metadata` which will
  be used if the `metadata` argument is NULL. Following Bioconductor
  best practices (fail-fast principle), create configuration via
  [`TSENAT_config()`](https://gallardoalba.github.io/TSENAT/reference/TSENAT_config.md)
  FIRST, then pass to `build_analysis()` at object construction time.
  This ensures invalid parameters are caught immediately, before
  analysis proceeds. See examples below for recommended usage pattern.

- skip:

  Logical. If TRUE, allow unmapped transcripts (transcripts not found in
  tx2gene mapping) and remove them from analysis. If FALSE (default),
  stop with an error when unmapped transcripts are detected. Useful for
  handling data with transcript IDs that don't match the annotation file
  provided.

- verbose:

  Logical. If TRUE, print informative messages during execution (e.g.,
  Salmon sample discovery, progress on data loading). Default: TRUE.

## Value

A `TSENATAnalysis` S4 object with:

- @se:

  The SummarizedExperiment containing transcript counts and metadata

- @config:

  Analysis configuration (empty list or user-provided)

- @diversity_results:

  Empty list (populated by calculate_diversity())

- @divergence_results:

  Empty list (populated by calculate_divergence())

- @lm_results:

  Empty list (populated by calculate_lm())

- @jackknife_results:

  Empty list (populated by jackknife functions)

- @plots:

  Empty list (populated by plotting functions)

- @metadata:

  Metadata with package version and creation timestamp

**CRITICAL: TPM and effective_length Requirements**

Both `tpm` and `effective_length` MUST be provided to ensure correct
filtering and normalization in downstream analysis:

- [`filter_analysis()`](https://gallardoalba.github.io/TSENAT/reference/filter_analysis.md)
  requires TPM data (stored in metadata). If TPM is missing, the
  function will fail with an explicit error message that guides you to
  pass it to `build_analysis()`.

- [`calculate_diversity()`](https://gallardoalba.github.io/TSENAT/reference/calculate_diversity.md)
  uses `effective_length` for length-normalized entropy calculations.

Following Bioconductor best practices (fail-fast principle), these are
explicit parameters, not optional. They must be passed at object
construction time:


    analysis <- build_analysis(
      readcounts = readcounts,
      metadata = metadata_df,
      tx2gene = gff3_file,
      tpm = tpm,                    # REQUIRED from Salmon output
      effective_length = effective_length,  # REQUIRED from Salmon output
      config = config
    )

This wrapper combines two steps into one:

1.  Call `.build_se()` to create a SummarizedExperiment from transcript
    counts

2.  Wrap the result in
    [`TSENATAnalysis()`](https://gallardoalba.github.io/TSENAT/reference/TSENATAnalysis.md)
    to create the analysis object

The returned object is ready for diversity analysis via
[`calculate_diversity()`](https://gallardoalba.github.io/TSENAT/reference/calculate_diversity.md).

If you need to inspect or filter the SummarizedExperiment before
creating the TSENATAnalysis object, call `.build_se()` and
[`TSENATAnalysis()`](https://gallardoalba.github.io/TSENAT/reference/TSENATAnalysis.md)
separately.

## Details

When using `salmon_dir`, the function automatically:

1.  Discovers all Salmon sample folders and quant.sf files

2.  Reads transcript counts (NumReads), TPM, and effective_length

3.  Extracts sample names from directory structure

4.  Creates count matrix ready for analysis

The `salmon_dir` parameter provides a convenient alternative to manually
constructing the `readcounts` matrix, especially useful in Galaxy
workflows.

## See also

[`TSENATAnalysis`](https://gallardoalba.github.io/TSENAT/reference/TSENATAnalysis.md)
for the S4 class structure
[`calculate_diversity`](https://gallardoalba.github.io/TSENAT/reference/calculate_diversity.md)
for computing Tsallis entropy

## Examples

``` r
# Create example transcript count data
set.seed(42)
n_genes <- 10
n_isoforms_per_gene <- 3
n_isoforms <- n_genes * n_isoforms_per_gene
n_samples <- 10

# Generate count matrix
counts <- matrix(rpois(n_isoforms * n_samples, lambda = 20),
                 nrow = n_isoforms, ncol = n_samples)
rownames(counts) <- paste0('TX_', 1:n_isoforms)
colnames(counts) <- paste0('Sample_', 1:n_samples)

# Create tx2gene mapping
tx2gene <- data.frame(
  Transcript = rownames(counts),
  Gene = rep(paste0('GENE_', 1:n_genes), each = n_isoforms_per_gene))

# Create sample metadata
metadata <- data.frame(
  sample = colnames(counts),
  condition = rep(c('control', 'treatment'), each = 5),
  row.names = colnames(counts))

# Build analysis object - use NAMED parameters to avoid confusion
# Method 1: With explicit tx2gene data.frame (most common)
config <- TSENAT_config(sample_col = 'sample', condition_col = 'condition')
analysis <- build_analysis(
  readcounts = counts,
  tx2gene = tx2gene,
  metadata = metadata,
  config = config)

# Verify the analysis object was created
analysis
#> TSENATAnalysis Object
#> ====================
#> Genes:       30
#> Samples:     10
#> Configuration: 27 parameters
#> Analysis status: EMPTY
#> Created: 2026-04-14 02:01:58.221242
#> 
print(dim(analysis))
#> NULL

# \donttest{
# Method 2: From Salmon quantification folder
# Requires directory structure like:
#   salmon_output/
#     sample1/quant.sf
#     sample2/quant.sf
#     ...
# 
# salmon_dir <- '/path/to/salmon/directory'
# 
# First create sample metadata matching Salmon sample names
# salmon_metadata <- data.frame(
#   condition = c('control', 'control', 'treatment', 'treatment'),
#   row.names = c('sample1', 'sample2', 'sample3', 'sample4')
# )
# 
# analysis_salmon <- build_analysis(
#   salmon_dir = salmon_dir,
#   tx2gene = 'annotation.gff3.gz',  # Auto-parsed from GFF3
#   metadata = salmon_metadata
# )
#
# Method 3: Hybrid - Salmon counts with manual tx2gene
# analysis_hybrid <- build_analysis(
#   salmon_dir = salmon_dir,
#   tx2gene = tx2gene,  # data.frame instead of file
#   metadata = salmon_metadata
# )
#
# Method 4: Pass metadata via config (parameter resolution pattern)
# cfg <- TSENAT_config()
# cfg$metadata <- metadata
# analysis_with_config <- build_analysis(
#   readcounts = counts,
#   tx2gene = tx2gene,
#   config = cfg
#   # Note: metadata argument omitted - will be read from config$metadata
# )
# }

# Advanced: Assigning metadata to assays after object creation
# When adding metadata to SummarizedExperiment assays, always use the
# S4Vectors namespace to ensure proper method dispatch:
#   
#   se <- getSE(analysis)
#   assay_with_ci <- SummarizedExperiment::assay(se, 'log2fc_ci')
#   S4Vectors::metadata(assay_with_ci)$lower <- ci_lower_bounds
#   S4Vectors::metadata(assay_with_ci)$upper <- ci_upper_bounds
#
# Note: Avoid using metadata(assay) without the namespace - this can
# cause silent failures in S4 object metadata assignment.
```
