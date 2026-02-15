# TCGA Luminal A breast cancer dataset (transcript-level)

Data from The Cancer Genome Atlas, containing transcript-level read
counts of 20 patients with Luminal A type breast cancer (primary tumor
and solid normal samples, 40 samples total). The dataset includes
transcript IDs in the first column and read counts for 40 samples in
subsequent columns.

## Usage

``` r
data(tcga_brca_luma)
```

## Format

A data frame with 1100 rows and 41 columns. The first column contains
transcript IDs, all additional columns contain RNA-sequencing read
counts for samples.

## Source

[TCGA via GDC](https://portal.gdc.cancer.gov/)

## References

The Cancer Genome Atlas Network (2012) Nature 490, 61–70
[doi:10.1038/nature11412](https://doi.org/10.1038/nature11412)

## Examples

``` r
data(tcga_brca_luma)
dim(tcga_brca_luma)
#> [1] 1100   41
head(tcga_brca_luma[1:4, 1:3])
#>   Transcript TCGA-A7-A0CH_N TCGA-A7-A0CH_T
#> 1    MXRA8.1        2858.04         743.56
#> 2    MXRA8.2         127.82          21.28
#> 3    MXRA8.3         370.22          94.38
#> 4    MXRA8.4        7472.00        3564.87
```
