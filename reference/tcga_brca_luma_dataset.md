# TCGA Luminal A breast cancer dataset

Data from The Cancer Genome Atlas, downloaded on 08th September, 2020.
It contains transcript level read counts of 20 patients with Luminal A
type breast cancer (primary tumor and solid normal samples).

## Usage

``` r
data(tcga_brca_luma_dataset)
```

## Format

A data frame with 996 rows and 41 columns. The first column contains
gene names, all additional columns contain RNA-sequencing read counts
for samples.

## Source

<https://portal.gdc.cancer.gov/>

## References

The Cancer Genome Atlas Network (2012) Nature 490, 61–70
([doi:10.1038/nature11412](https://doi.org/10.1038/nature11412) )
