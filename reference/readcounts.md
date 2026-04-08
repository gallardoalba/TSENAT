# Example transcript-level read counts dataset

A dataset containing transcript-level transcript-level read count
abundances from salmon quantification. The dataset includes samples as
columns (SRR14800475-SRR14800490) and transcript IDs as row names.

## Usage

``` r
data(readcounts)
```

## Format

A data frame with 3514 transcripts (rows) and 16 samples (columns). Row
names are transcript IDs (ENST format) and column names are sample IDs.
Values represent TPM-normalized read counts.

## Value

A matrix (data frame) containing transcript-level TPM-normalized read
counts from salmon quantification. Rows represent transcripts (ENST
format IDs) and columns represent samples. Aliases include `readcounts`
(main read count matrix), `tpm` (TPM values), and `effective_length`
(transcript lengths from salmon).

## Examples

``` r
data(readcounts, package = 'TSENAT')
dim(readcounts)
#> [1] 3089   16
head(readcounts[1:4, 1:3])
#>                    SRR14800475 SRR14800476 SRR14800477
#> ENST00000162391.8      757.160     452.836     329.451
#> ENST00000187762.7       36.570      83.967      56.887
#> ENST00000221130.11     935.366    3356.274    3008.108
#> ENST00000221818.5        3.000       0.000       0.000
```
