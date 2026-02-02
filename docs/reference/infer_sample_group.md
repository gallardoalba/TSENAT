# Infer sample group from sample names

Infer sample group from sample names

## Usage

``` r
infer_sample_group(
  sample_names,
  suffix_sep = "_",
  suffix_map = NULL,
  tcga_map = NULL,
  coldata = NULL,
  coldata_sample_col = "Sample",
  coldata_condition_col = "Condition",
  prefer_suffix = TRUE,
  default = NA_character_
)
```

## Arguments

- sample_names:

  Character vector of sample names.

- suffix_sep:

  Character separator to detect suffix groups (default "\_").

- suffix_map:

  Named character vector mapping suffix tokens (case- insensitive) to
  group labels, e.g. c(N = "Normal", T = "Tumor").

- tcga_map:

  Named character vector mapping TCGA two-digit codes to group labels,
  e.g. c("01" = "Tumor", "11" = "Normal").

- coldata:

  Optional data.frame or named vector providing mapping from sample
  names to conditions. If a data.frame, specify \`coldata_sample_col\`
  and \`coldata_condition_col\` for the relevant columns.

- coldata_sample_col:

  Column name in \`coldata\` indicating sample IDs (default "Sample").

- coldata_condition_col:

  Column name in \`coldata\` for condition/label (default "Condition").

- prefer_suffix:

  Logical; if TRUE, prefer suffix-based inference when both patterns
  match.

- default:

  Character scalar returned when no mapping applies (default
  NA_character\_).

## Value

Character vector of group labels (or the \`default\` value) with same
length as \`sample_names\`.

## Examples

``` r
infer_sample_group(c("S1_N", "TCGA-XX-01A"))
#> [1] "N"  "01"
infer_sample_group(c("S1_N", "S2_T"), suffix_map = c(
    N = "Normal", T =
        "Tumor"
))
#> [1] "Normal" "Tumor" 
tcga_map <- c("01" = "Tumor")
infer_sample_group(c("TCGA-XX-01A", "Sample_N"),
    tcga_map = tcga_map,
    prefer_suffix = FALSE
)
#> [1] "Tumor" "N"    
```
