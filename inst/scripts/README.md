## Reproduction of processed example dataset

This package ships a small example dataset for vignette and tests: `data/tcga_brca_luma_dataset.RData` and a few small support files under `inst/extdata/` (for example `tx2gene.tsv` and `coldata.tsv`). The raw TCGA isoform expression files used to build the example dataset are NOT included in the package due to their size and licensing.

If you need to reproduce the example dataset locally or re-create it using the original raw files, follow these steps:

1. Obtain raw files
   - The raw isoform expression quantification files were downloaded from the GDC via `TCGAbiolinks` (TCGA-BRCA, legacy isoform expression quantification). See the `data-raw/tcga_brca_luma_dataset.R` script for details on which files and query parameters were used.

   Example R commands (run interactively or in a script):

   ```r
   # Install if necessary
   if (!requireNamespace("TCGAbiolinks", quietly = TRUE)) install.packages("TCGAbiolinks")
   library(TCGAbiolinks)

   query <- GDCquery(
     project = "TCGA-BRCA",
     legacy = TRUE,
     data.category = "Gene expression",
     data.type = "Isoform expression quantification",
     experimental.strategy = "RNA-Seq",
     sample.type = c("Primary Tumor", "Solid Tissue Normal")
   )

   # Download raw files to a directory of your choice
   GDCdownload(query, directory = "path/to/download_dir")
   ```

   Or run the packaged script which implements the above workflow (adjust paths inside the script as needed):

   ```sh
   Rscript data-raw/tcga_brca_luma_dataset.R
   ```

2. Recreate the example dataset
   - Run the R script in `data-raw/tcga_brca_luma_dataset.R`. It contains detailed steps for query and preprocessing and will write `tcga_brca_luma_dataset.RData` to the current working directory.
   - Optionally, copy the resulting `.RData` into `data/` or use `usethis::use_data()` to include in the package during development.

3. What is provided in `inst/extdata`
   - Small auxiliary lookup files used by examples and vignettes are included here (for example `tx2gene.tsv`, `coldata.tsv`, and `tcga_sample_ids.tsv`). These are small text files safe for inclusion in a package and intended for demonstration and reproducibility of vignette examples.

4. Automating the small processed files
   - A helper script is provided at `inst/scripts/generate_extdata.R`. This script documents the expected outputs, optionally runs `data-raw/tcga_brca_luma_dataset.R` via `Rscript` (when invoked with `--run-data-raw`), and writes an MD5 checksum into `inst/extdata/` for verification.

5. Licensing and access notes
   - Raw TCGA files are subject to GDC access and licensing terms and are not distributed with this package. The `data-raw/` script documents the source (GDC via `TCGAbiolinks`) and the commands to obtain these files; consult the GDC documentation for any controlled-access or token requirements.