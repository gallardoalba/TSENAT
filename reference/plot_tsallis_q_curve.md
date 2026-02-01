# Plot median +- IQR of Tsallis entropy across q values by group

This reproduces the \`tsallis-q-curve-mean-sd\` plot from the vignette:
for each q value, compute per-gene Tsallis entropy per sample, then
summarize across genes by group (median and IQR) and plot median with a
ribbon spanning median +- IQR/2.

## Usage

``` r
plot_tsallis_q_curve(
  readcounts,
  genes,
  q_values = seq(0.01, 2, by = 0.01),
  group_pattern = "_N$",
  group_names = c("Normal", "Tumor")
)
```

## Arguments

- readcounts:

  Numeric matrix or data.frame with transcripts as rows and samples as
  columns.

- genes:

  Character vector assigning a gene id to each row of \`readcounts\`.

- q_values:

  Numeric vector of q values to evaluate (default
  \`seq(0.01,2,by=0.01)\`).

- group_pattern:

  Regular expression used to detect the first group in sample names
  (default \`"\_N\$"\`).

- group_names:

  Character vector of length 2 with names for groups (default
  \`c("Normal","Tumor")\`).

## Value

A \`ggplot\` object showing median +- IQR across q values by group.
