#!/usr/bin/env Rscript
# Helper to (re)generate the small processed example dataset and extdata files.
# This script documents and automates the steps used to create files in
# `inst/extdata/` from the TCGA raw files. It does NOT download raw files on its
# own; see `data-raw/tcga_brca_luma_dataset.R` for download instructions.

# Usage:
#   Rscript inst/scripts/generate_extdata.R [--run-data-raw]
# Options:
#   --run-data-raw  : Run `data-raw/tcga_brca_luma_dataset.R` (may require
#                     adjusting path variables in that script and availability
#                     of raw files downloaded from GDC).

args <- commandArgs(trailingOnly = TRUE)
run_data_raw <- any(grepl("--run-data-raw", args))

ensure <- function(x) stopifnot(file.exists(x))

if (run_data_raw) {
    message("Running data-raw/tcga_brca_luma_dataset.R (it may download files and/or process downloaded files)")
    # This will execute the data-raw script. Ensure you have adjusted `location`
    # and other path variables inside that script to point to your downloaded
    # files before running this helper with --run-data-raw.
    system2("Rscript", c("data-raw/tcga_brca_luma_dataset.R"), wait = TRUE)
}

# Check for expected outputs
expected_rdata <- "data/tcga_brca_luma_dataset.RData"
if (!file.exists(expected_rdata)) {
    message("Expected dataset not found: ", expected_rdata, "\nRun the data-raw script to generate it, or copy the produced file into data/.")
} else {
    message("Found dataset: ", expected_rdata)
}

# Ensure inst/extdata contains the small lookup files
ext_files <- c("inst/extdata/tx2gene.tsv", "inst/extdata/coldata.tsv")
for (f in ext_files) {
    if (!file.exists(f)) {
        message("Warning: expected extdata file missing: ", f)
    } else {
        message("Found extdata file: ", f)
    }
}

# Compute and write checksum for the generated dataset (helpful for verification)
if (file.exists(expected_rdata)) {
    md5 <- tools::md5sum(expected_rdata)
    out <- file.path("inst/extdata", "tcga_brca_luma_dataset.md5")
    writeLines(paste(md5, basename(expected_rdata)), con = out)
    message("Wrote checksum to: ", out)
}

message("Done. See data-raw/tcga_brca_luma_dataset.R and inst/scripts/README.md for details and licensing notes.")
