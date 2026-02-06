#!/usr/bin/env Rscript
# Generate synthetic TCGA BRCA dataset for vignette examples
#
# Usage: Rscript scripts/generate_tcga_rdata.R

tryCatch(
    {
        cat("Generating synthetic TCGA BRCA LUMA dataset...\n")
        
        dir.create("inst/extdata", recursive = TRUE, showWarnings = FALSE)
        set.seed(2026)
        
        n_genes <- 996
        n_pairs <- 20  # creates 40 sample columns (20 normal/tumor pairs)
        
        genes <- paste0("G", seq_len(n_genes))
        
        # Build sample names: P1_N, P1_T, P2_N, P2_T, ...
        sample_names <- unlist(
            lapply(seq_len(n_pairs),
                   function(i) c(paste0("P", i, "_N"), paste0("P", i, "_T")))
        )
        
        # Generate counts with slight group differences
        # Tumor samples have slightly higher expression lambda
        counts_mat <- matrix(nrow = n_genes, ncol = length(sample_names))
        
        for (i in seq_len(n_genes)) {
            is_tumor <- grepl("_T$", sample_names)
            lambda_tumor <- stats::rgamma(sum(is_tumor), shape = 2, rate = 0.05)
            lambda_normal <- stats::rgamma(sum(!is_tumor), shape = 2, rate = 0.08)
            
            counts_tumor <- stats::rpois(sum(is_tumor), lambda = lambda_tumor)
            counts_normal <- stats::rpois(sum(!is_tumor), lambda = lambda_normal)
            
            counts_mat[i, is_tumor] <- counts_tumor
            counts_mat[i, !is_tumor] <- counts_normal
        }
        
        # Create data frame with genes in first column
        tcga_brca_luma_dataset <- as.data.frame(counts_mat)
        tcga_brca_luma_dataset <- cbind(genes = genes, tcga_brca_luma_dataset)
        colnames(tcga_brca_luma_dataset)[2:ncol(tcga_brca_luma_dataset)] <- sample_names
        
        # Save as RData
        output_file <- "inst/extdata/tcga_brca_luma_dataset.RData"
        save(tcga_brca_luma_dataset, file = output_file)
        
        cat("✓ Dataset generated:", output_file, "\n")
        cat("  Dimensions:", nrow(tcga_brca_luma_dataset), "genes ×", 
            ncol(tcga_brca_luma_dataset)-1, "samples\n")
    },
    error = function(e) {
        cat("✗ ERROR:", conditionMessage(e), "\n")
        quit(status = 1)
    }
)
