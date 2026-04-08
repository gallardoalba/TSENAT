# Extract SummarizedExperiment from TSENATAnalysis

Retrieve the underlying SummarizedExperiment object containing count
data and sample metadata.

## Usage

``` r
getSE(object, ...)
```

## Arguments

- object:

  TSENATAnalysis object

- ...:

  Additional arguments (for method compatibility)

## Value

SummarizedExperiment containing count matrix and sample metadata

## Details

The SummarizedExperiment object returned by getSE() contains: - assays:
count matrices (transcript-level read counts) - rowData: transcript
information and gene assignments - colData: sample metadata (sample
types, conditions, etc.)

## Examples

``` r
# Load example data and build analysis object
data(readcounts, package = 'TSENAT')
metadata_df <- read.table(system.file('extdata', 'metadata.tsv', package = 'TSENAT'),
  header = TRUE, sep = '\t')
gff3_file <- system.file('extdata', 'annotation.gff3.gz', package = 'TSENAT')

# Build TSENATAnalysis object
config <- tsenat_config(sample_col = 'sample', condition_col = 'condition')
analysis <- build_analysis_s4(readcounts = readcounts, 
                             tx2gene = gff3_file, 
                             metadata = metadata_df,
                             config = config,
                             tpm = tpm,
                             effective_length = effective_length)

# Extract the underlying SummarizedExperiment
se <- getSE(analysis)

# Explore the SummarizedExperiment structure
nrow(se)  # Number of transcripts
#> [1] 3089
ncol(se)  # Number of samples
#> [1] 16
SummarizedExperiment::assayNames(se)  # Available assay matrices
#> [1] "counts"
SummarizedExperiment::colData(se)  # Sample metadata
#> DataFrame with 16 rows and 5 columns
#>               condition sample_type sample_base paired_samples   sample_id
#>             <character> <character> <character>    <character> <character>
#> SRR14800481      normal      normal           A              A SRR14800481
#> SRR14800480      normal      normal           B              B SRR14800480
#> SRR14800477      normal      normal           C              C SRR14800477
#> SRR14800476      normal      normal           D              D SRR14800476
#> SRR14800475      normal      normal           E              E SRR14800475
#> ...                 ...         ...         ...            ...         ...
#> SRR14800486       tumor       tumor           D              D SRR14800486
#> SRR14800485       tumor       tumor           E              E SRR14800485
#> SRR14800484       tumor       tumor           F              F SRR14800484
#> SRR14800483       tumor       tumor           G              G SRR14800483
#> SRR14800482       tumor       tumor           H              H SRR14800482
```
