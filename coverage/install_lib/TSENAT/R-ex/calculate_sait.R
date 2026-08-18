### Name: calculate_sait
### Title: Calculate SAIT interactions and store in TSENATAnalysis
### Aliases: calculate_sait

### ** Examples

## Not run: 
##D # Create test analysis with appropriate sample structure for paired design
##D # Note: requires lme4 package for LMM fitting; uses synthetic data
##D set.seed(42)
##D 
##D # Create transcript-level counts with biological signal
##D # Note: Use adequate complexity (transcripts/genes, samples, expression) 
##D # to avoid filtering away all genes during diversity computation
##D n_genes <- 50
##D n_transcripts_per_gene <- 30
##D n_transcripts <- n_genes * n_transcripts_per_gene
##D n_samples <- 16  # 8 subjects x 2 conditions (paired design)
##D 
##D # Generate counts with clear biological signal
##D control_idx <- seq(1, n_samples, by = 2)
##D treatment_idx <- seq(2, n_samples, by = 2)
##D 
##D counts <- matrix(0, nrow = n_transcripts, ncol = n_samples)
##D for (j in seq_len(n_samples)) {
##D   lambda <- if (j %in% control_idx) 100 else 180
##D   counts[, j] <- rpois(n_transcripts, lambda = lambda)
##D }
##D counts <- pmax(counts, 50)  # Ensure minimum expression
##D 
##D rownames(counts) <- paste0('TX_', seq_len(n_transcripts))
##D colnames(counts) <- paste0('Sample_', seq_len(n_samples))
##D 
##D # Create rowData with gene mapping (tx2gene structure)
##D rowdata <- data.frame(
##D   transcript_id = rownames(counts),
##D   gene_id = rep(paste0('GENE_', 1:n_genes), 
##D                 each = n_transcripts_per_gene),
##D   row.names = rownames(counts)
##D )
##D 
##D # Create colData with paired design metadata
##D coldata <- data.frame(
##D   sample_id = colnames(counts),
##D   condition = rep(c('control', 'treatment'), 
##D                   length.out = n_samples),
##D   subject = rep(paste0('Subject_', 1:8), 
##D                 length.out = n_samples),
##D   row.names = colnames(counts)
##D )
##D 
##D # Build SummarizedExperiment
##D se <- SummarizedExperiment::SummarizedExperiment(
##D   assays = list(counts = counts),
##D   rowData = S4Vectors::DataFrame(rowdata),
##D   colData = S4Vectors::DataFrame(coldata)
##D )
##D 
##D # Add tx2gene metadata for gene-level aggregation
##D S4Vectors::metadata(se)$tx2gene <- 
##D   data.frame(Transcript = rowdata$transcript_id,
##D              Gene = rowdata$gene_id)
##D 
##D # Initialize TSENATAnalysis
##D analysis <- TSENATAnalysis(se = se, config = list())
##D 
##D # Compute diversity (prerequisite for SAIT interaction analysis)
##D analysis <- calculate_diversity(
##D   analysis, 
##D   q = c(0.5, 1.0, 1.5, 2.0, 2.5)
##D )
##D 
##D # Calculate q x condition interactions using GAM
##D analysis <- suppressWarnings(calculate_sait(
##D   analysis,
##D   condition_col = 'condition',
##D   method = 'gam'
##D ))
##D 
##D # View top interaction results using unified accessor (first 3 genes)
##D res <- results(analysis, type = "sait")
##D if (!is.null(res)) head(res, 3)
## End(Not run)




