# Changelog

## TSENAT 0.99.0

TSENAT computes Tsallis entropy (S_q) and related Hill numbers (D_q) to
quantify isoform/transcript-level diversity per gene. The package
accepts matrices, tximport-style lists, and SummarizedExperiment
objects, returns per-gene diversity assays, and provides plotting and
simple differential-analysis helpers. A vignette demonstrates a full
workflow on a TCGA Luminal A example.

Highlights:

``` R
Core: calculate_tsallis_entropy() — compute S_q and D_q (supports multiple q values and the q->1 Shannon limit).
Batch processing: calculate_diversity() applies the calculation across transcripts/genes and returns a SummarizedExperiment with a "diversity" assay.
Differential helpers: group comparisons, fold changes and simple tests.
Visualizations: density, violin, q-curve, MA, volcano, and top-transcripts plots (return ggplot objects).
```
