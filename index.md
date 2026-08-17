![TSENAT: Tsallis Entropy Analysis Toolbox. TSENAT is a R package for
quantifying and modeling isoform-usage complexity across RNA-seq samples
using Tsallis entropy - a scale-dependent information-theoretic measure
of transcript
heterogeneity.](https://github.com/user-attachments/assets/b89df329-c697-40ce-a037-0f82fc782a23)

## The Problem

RNA-seq analysis typically addresses one of three complementary
biological questions: (1) which genes change in total abundance?, (2)
which individual transcripts shift their abundance?, (3) which
transcripts shift their *proportions* within genes, independent of total
abundance?

However, a fourth question remains largely invisible: **how does
isoform-usage complexity change?**

## The Solution

TSENAT captures isoform-usage complexity throught Tsallis entropy, whose
entropic index parameter (q) enables to evaluate different data
dimensions:

- **Low q** (e.g., 0.5): Focuses on rare isoforms - detects if diversity
  is maintained or collapsed

- **Mid q** (e.g., 1.0): Balanced view (Shannon entropy) - overall
  isoform complexity

- **High q** (e.g., 2.0): Focuses on dominant isoforms - detects
  dominance shifts

By examining diversity across multiple entropic indices (q-values),
TSENAT allows to identify scale-dependent diversity changes - the
hallmark of coordinate isoform switching.

## The Mathematics Behind Tsallis Entropy

Tsallis entropy is a parametric family of diversity measures that
generalizes Shannon entropy and enables tuning sensitivity to different
scales of isoform organization.

**Mathematical Definition**: For a discrete probability vector
$`p = (p_1, \ldots, p_n)`$ representing isoform proportions within a
gene, Tsallis entropy is:

``` math
S_q = \frac{1 - \sum_{i=1}^{n} p_i^q}{q - 1}
```

This elegant formula unifies diverse diversity concepts at specific
entropic indices (q-values):

- **q = 0**: Richness — Simple count of expressed isoforms; emphasizes
  rare variants most strongly.
- **q = 1**: Shannon entropy — Standard information-theoretic measure;
  balanced weighting across scales.
- **q = 2**: Gini-Simpson index — Probability that two randomly-drawn
  transcripts are different; robust to rare variants.

### Divergence Analysis: Measuring Information-Theoretic Distance Between Conditions

While Tsallis entropy quantifies diversity *within* a single
distribution, **Tsallis divergence** $`D_q`$ measures the
information-theoretic distance *between* two distributions.

**Mathematical Definition** (Furuichi formula): For two probability
distributions $`P`$ and $`Q`$ representing isoform proportions in
control and treatment conditions, Tsallis divergence is:

``` math
D_q(P||Q) = \frac{\sum_i p_i^q \cdot q_i^{1-q} - 1}{q-1}
```

where $`p_i`$ and $`q_i`$ are the probability values at position $`i`$.

Tsallis divergence enables the quantification of how fundamentally
different the isoform complexity patterns are between experimental
conditions.

## Installation

**Requirements:** R \>= 4.5.0

The easiest way to install TSENAT is through Bioconda:

``` bash
conda install -c bioconda r-tsenat
```

To install the latest development version from GitHub:

``` r

remotes::install_github("gallardoalba/TSENAT")
```

## Quick Start

### Load Example Data

Start by loading the built-in example dataset from TSENAT, which
includes transcript-level read counts, TPM values, and effective lengths
from Salmon quantification. Then load the sample metadata and annotation
file that describe your experimental design.

``` r

suppressMessages({
  library(TSENAT)
  library(SummarizedExperiment)})

# Load example dataset (includes readcounts, tpm, and effective_length)
data(readcounts)
readcounts <- as.matrix(readcounts)

# Load sample metadata and annotation
metadata_df <- read.table(
  system.file("extdata", "metadata.tsv", package = "TSENAT"),
  header = TRUE,
  sep = "\t")

gff3_file <- system.file(
  "extdata", "annotation.gff3.gz", package = "TSENAT")
```

### Create configuration and build analysis

Create a configuration object specifying your experimental design
parameters (sample/condition columns from metadata) before building the
analysis object. This fail-fast pattern ensures invalid parameters are
caught immediately before processing begins.

``` r

## Create configuration file
config <- TSENAT_config(
  q = seq(0, 2, by = 0.05),
  sample_col = "sample",
  condition_col = "condition",
  paired = TRUE,
  subject_col = "paired_samples",
  control = "normal",
  nthreads = 2
)

## Build TSENATAnalysis object
analysis <- build_analysis(
  readcounts = readcounts, 
  tx2gene = gff3_file,
  metadata = metadata_df,
  config = config,
  tpm = tpm,
  effective_length = effective_length)
```

### Orchestration Function

The
[`TSENAT()`](https://gallardoalba.github.io/TSENAT/reference/TSENAT.md)
function provides a complete, automated analysis pipeline in a single
call. It takes your configured `TSENATAnalysis` object and executes all
downstream analysis steps: entropy computation, statistical testing for
entropic index (q-value) by condition interactions, and rich
visualization. This is the recommended entry point for most users—it
orchestrates the full workflow while respecting your configuration
parameters (entropic indices, design, bootstrap settings, etc.) and
handles output management seamlessly.

``` r

# Returns: Fully configured TSENATAnalysis object
result <- TSENAT(analysis)
```

### Accessing Results

After running the analysis, retrieve results using the unified
[`results()`](https://gallardoalba.github.io/TSENAT/reference/results.md)
interface:

``` r

# Unified interface: request results by type
diversity_results <- results(result, type = "diversity")
sait_results <- results(result, type = "sait")
jackknife_results <- results(result, type = "jackknife", q = 1)

# Filter and rank results flexibly
ranked_sait <- results(result, type = "sait", rankBy = "padj", n = 50)
filtered_diversity <- results(result, type = "diversity", q = 1.0)

# Retrieve plots
diversity_plot <- results(result, type = "diversity", plot=TRUE)
print(diversity_plot)
```

![Isoform diversity profiles across
q-values](https://raw.githubusercontent.com/gallardoalba/TSENAT/gh-pages/articles/TSENAT_files/figure-html/fig-1-isoform-diversity-profiles-1.png)

Isoform diversity profiles across q-values

``` r

# Extract the scale-dependent genes plot
sait_plot <- results(result, type = "sait", plot = TRUE)
print(sait_plot)
```

![Scale-dependent
genes](https://raw.githubusercontent.com/gallardoalba/TSENAT/gh-pages/articles/TSENAT_files/figure-html/fig-2-scale-dependent-genes-1.png)

Scale-dependent genes

## Documentation

For a complete walkthrough of the analysis pipeline with real biological
examples, see the [main package
vignette](https://gallardoalba.github.io/TSENAT/articles/TSENAT.html).
This includes theory background, step-by-step explanations of each
analysis function, and interpretation guidance for understanding your
results. It covers the entire workflow from loading Salmon-quantified
data through entropy computation, statistical testing, and
visualization.

## Statistical Inference Methods

TSENAT provides a flexible statistical framework optimized for
entropy-based diversity analysis. The recommended main workflow (paired
design) relies on **Generalized Additive Mixed Models (GAMM)**: paired
designs are fit with
[`nlme::lme`](https://rdrr.io/pkg/nlme/man/lme.html) using regression
splines (`ns(q, df = 3) × condition`), a subject random intercept, and
AR(1) correlation within each subject × condition block, with the
interaction tested via a marginal F-test.

The statistical methods available in TSENAT include:

- **Scale-adaptive interaction tests (SAIT)**: Multiple modeling
  approaches for repeated q-ordered measurements. GAM/GAMM, LMM, GEE and
  FPCA model within-subject correlation with a working AR(1) structure
  within each subject × condition block, and are parametrized to handle
  the non-normality and heteroscedasticity characteristic of entropy
  data.
- **Aligned Rank Transform (ART)**: State-of-the-art non-parametric
  interaction testing via the ARTool package (Kay et al. 2021). Strips
  main effects before ranking (“alignment”) to properly preserve
  interaction structure — addressing the known limitation of classical
  rank-transform methods for factorial designs. The Conover-Iman Rank
  Transform remains available as a fallback via `method='rt'`.
- **M-estimation**: Robust location estimation for group comparison
  using iteratively re-weighted least squares, resistant to outliers.
- **Jackknife isoform switching (JIS)**: Leave-one-out resampling to
  identify transcripts with condition-specific switching patterns and
  quantify their influence on entropy differences.

## Related Packages

Below is how TSENAT complements other Bioconductor tools:

| Tool | Answers | TSENAT Difference |
|----|----|----|
| **edgeR, DESeq2** | Do individual transcripts increase/decrease in expression? Do genes change in total abundance? | TSENAT measures isoform-usage complexity changes independent of transcript or gene abundance |
| **DEXSeq, DRIMSeq** | Which transcripts shift their *proportions* within genes, independent of abundance changes? | TSENAT detects whether the isoform landscape consolidates or fragments |
| **IsoformSwitchAnalyzeR** | Which *individual isoforms* switch; what are the *functional consequences*? | TSENAT measures overall isoform diversity and diversity shifts rather than cataloging individual transcript switches or predicting functional consequences; complements switch identification with diversity patterns |
| **SplicingFactory** | What is the *overall isoform diversity*? | TSENAT extends with scale-dependent diversity (q-spectrum) vs fixed measures |
| **Kallisto, Salmon** | How many reads per transcript? | TSENAT uses their quantification as input; adds diversity analysis layer |

## Native Salmon Integration

TSENAT is specifically engineered to work seamlessly with Salmon
quantification output. Rather than requiring manual parsing or format
conversion, TSENAT automatically discovers transcript-level
quantification files across your Salmon output directory and integrates
them directly into the analysis pipeline. This tight integration means
you can move from Salmon quantification to entropy analysis without
intermediate data manipulation—the raw `quant.sf` files are all you
need. TSENAT discovers these files automatically, validates their
compatibility with your experimental design, and handles
length-correction and normalization as part of the diversity computation
workflow.

To get started with Salmon-quantified data:

``` r

suppressMessages(library(TSENAT))

# Prepare configuration FIRST
config <- TSENAT_config(
  q = seq(0, 2, by = 0.1),
  condition_col = "condition",
  sample_col = "sample",
  subject_col = "paired_samples"
)

# Build analysis directly from Salmon output directory
analysis <- build_analysis(
  salmon_dir = "path/to/salmon_output",  # Auto-discovers all quant.sf files
  tx2gene = "path/to/annotation.gff3.gz",  # Transcript-to-gene mapping
  metadata = metadata_df,  # Sample metadata (see structure below)
  config = config
)

# Run analysis pipeline
analysis <- TSENAT(analysis)
```

### Salmon Directory Structure

Your Salmon output must be organized with one folder per sample, each
containing a quant.sf file with transcript-level quantification. Sample
folder names must exactly match the row names in your metadata file
(case-sensitive) to ensure correct sample attribution.

TSENAT expects Salmon output organized with one subdirectory per sample:

    salmon_output/
    ├── sample_1/
    │   └── quant.sf
    ├── sample_2/
    │   └── quant.sf
    ├── sample_3/
    │   └── quant.sf
    └── sample_4/
        └── quant.sf

### Metadata File Structure

Expected TSV file format:

    sample      condition    paired_samples
    sample_1    normal            A
    sample_2    normal            B
    sample_3    tumor             A
    sample_4    tumor             B

Key requirements:

- `sample` column: must match Salmon folder names exactly.
- `condition` column: experimental groups.
- `paired_samples` column: required if using paired designs; identifier
  for matched samples.

## Tests Coverage

Testing is vital in research as it ensures the validity and reliability
of results, which is essential for accurately interpreting findings. The
report about the current testing coverage can be found
[here](https://app.codecov.io/gh/gallardoalba/TSENAT).

## Citation

If you use TSENAT in your research, please cite:

``` bibtex
@software{gallardo2026tsenat,
  title={TSENAT: Tsallis Entropy Analysis Toolbox},
  author={Gallardo Alba, Cristóbal},
  url={https://github.com/gallardoalba/TSENAT},
  year={2026}
}
```

## License and Attribution

This project is licensed under the GNU General Public License v3.0
(GPL-3). See [LICENSE](https://gallardoalba.github.io/TSENAT/LICENSE)
for details.

TSENAT builds upon the [SplicingFactory
package](https://github.com/esebesty/SplicingFactory), extending it with
specialized focus on Tsallis entropy analysis.

> **“If I ever come back from the past, it’s to create a cyclone.”**
>
> – Juan José Lozano
