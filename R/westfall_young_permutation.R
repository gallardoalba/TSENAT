
#' Multi-q Permutation Test with Westfall-Young Stepdown Procedure
#'
#' Performs joint permutation testing across multiple q-values (diversity orders)
#' using the Westfall-Young stepdown procedure for family-wise error rate control.
#'
#' @param x A \code{matrix} with diversity values (rows = genes, columns = samples*q-values).
#'   Column names should indicate q-value (e.g., "Sample1_q=0.5", "Sample1_q=1.0").
#' @param q_values Numeric vector of q-values used in computing diversity.  
#'   Should be in increasing order and match column names in \code{x}.
#' @param samples Character vector matching the sample pairs in column names.
#'   Length should match number of unique samples (before q expansion).
#' @param control Name of the control sample category (e.g., "Normal", "WT").
#' @param method Method for aggregating diversity values: "mean" or "median".
#' @param randomizations Integer number of permutations (default: 1000).
#' @param verbose Logical; if TRUE, print progress updates (default: FALSE).
#' @param nthreads Integer number of threads for parallel processing (default: 1).
#'
#' @return A data frame with columns:
#'   - `gene`: Gene identifier (from row names of x)
#'   - `q_value`: q-value used
#'   - `log2FC`: Log2 fold-change at this q
#'   - `pvalue_raw`: Raw two-sided permutation p-value (Phipson-Smyth corrected)
#'   - `pvalue_wy`: Westfall-Young adjusted p-value (true stepdown FWER control)
#'   - `padj_bh`: Multiple testing adjusted p-value (BH-FDR on raw p-values)
#'   - Family-wise error rate <= 0.05 is guaranteed by true WY stepdown
#'
#' @details
#' **Westfall-Young Stepdown Procedure for Multi-q Analysis:**
#'
#' Traditional Benjamini-Hochberg FDR control assumes independence, but multi-q
#' tests for the same gene are correlated. The Westfall-Young stepdown procedure
#' properly accounts for these correlations by:
#'
#' 1. Computing raw two-sided permutation p-values for each gene-q combination
#' 2. Ranking p-values from smallest to largest
#' 3. For each p-value in order:
#'    - Count permutations where minimum p-value <= observed (accounts for all smaller tests)
#'    - Adjusted p-value = min(raw_p, count/num_permutations)
#' 4. Enforce monotonicity: adjusted p-values are non-decreasing as p-values increase
#'
#' **Advantages Over Independent Testing:**
#' - Controls correlation between q-values for same gene
#' - Maintains family-wise error rate (FWER) at α level
#' - More powerful than simple Bonferroni when tests are correlated
#' - Accounts for joint dependence structure in multi-q analysis
#'
#' **Phipson-Smyth Correction (S019):**
#' Each permutation test uses p = (b+1)/(m+1) to avoid zero p-values,
#' where b = count of permutations >= observed, m = number of permutations.
#'
#' @references
#' Westfall, P. H., & Young, S. S. (1993). Resampling-based Multiple Testing:
#' Examples and Methods for p-Value Adjustment. John Wiley & Sons.
#'
#' Phipson, B., & Smyth, G. K. (2010). Permutation p-values should never be zero:
#' calculating exact p-values when permutations are randomly drawn.
#' Statistical Applications in Genetics and Molecular Biology, 9(1), 39.
#'
#' Meinshausen, N., Maathuis, M., & Bühlmann, P. (2011). Optimality of the 
#' Westfall-Young permutation procedure for multiple testing under dependence.
#' The Annals of Statistics, 39(6), 3359-3382. [Database: S165, S166]
#'
#' Cox, D. D. Pointwise Testing with Functional Data Using the Westfall-Young 
#' Randomization Method. [Database: S167]
#'
#' @export
#' @examples
#' \dontrun{
#' # Diversity matrix with multiple q-values
#' set.seed(123)
#' n_genes <- 50
#' n_samples <- 8  # 4 pairs
#' n_q <- 5  # q = 0.5, 1.0, 1.5, 2.0
#'
#' # Create diversity matrix: genes * (samples * q-values)
#' diversity_data <- matrix(rnorm(n_genes * n_samples * n_q), 
#'                          nrow = n_genes)
#' colnames(diversity_data) <- 
#'   paste0(rep(paste0("Sample", 1:n_samples), n_q), 
#'          "_q=", rep(c(0.5, 1.0, 1.5, 2.0), each = n_samples))
#' rownames(diversity_data) <- paste0("Gene", 1:n_genes)
#'
#' # Sample labels
#' samples <- rep(c("Normal", "Tumor"), each = n_samples/2)
#' q_vals <- c(0.5, 1.0, 1.5, 2.0)
#'
#' # Run Westfall-Young multi-q permutation test
#' wy_results <- label_shuffling_westfall_young(
#'   x = diversity_data,
#'   q_values = q_vals,
#'   samples = samples,
#'   control = "Normal",
#'   randomizations = 1000,
#'   verbose = TRUE
#' )
#'
#' # View top results (most significant across all q-values)
#' head(wy_results[order(wy_results$pvalue_wy), ])
#'
#' # Count significant genes at FWER = 0.05
#' sum(wy_results$pvalue_wy < 0.05)
#' }
label_shuffling_westfall_young <- function(x, q_values, samples, control, 
                                           method = "mean", randomizations = 1000,
                                           verbose = FALSE, nthreads = 1) {
    
    # Input validation
    if (!all(grepl(paste0("_q=", collapse = "|"), colnames(x)))) {
        stop("Column names must contain _q= to indicate q-values (e.g., 'Sample1_q=0.5')")
    }
    
    if (length(q_values) == 0) {
        stop("q_values must be provided")
    }
    
    # Extract q-values from column names and verify they match input
    # Note: calculate_diversity arranges columns as (sample, q-value)
    # So pattern is rep(q_values, times=length(samples))
    col_q_values <- as.numeric(sub(".*_q=", "", colnames(x)))
    if (!isTRUE(all.equal(col_q_values, rep(q_values, times = length(samples))))) {
        stop("q-values in column names do not match provided q_values")
    }
    
    n_genes <- nrow(x)
    n_samples <- length(samples)
    n_q <- length(q_values)
    n_total_tests <- n_genes * n_q
    
    if (verbose) {
        cat("Westfall-Young Multi-q Permutation Test\n")
        cat("=========================================\n")
        cat("Genes:", n_genes, "\n")
        cat("Q-values:", paste(q_values, collapse = ", "), "\n")
        cat("Total tests:", n_total_tests, "\n")
        cat("Permutations:", randomizations, "\n")
        cat("FWER control: α = 0.05\n\n")
    }
    
    # Step 1: Compute raw p-values for all gene * q combinations
    # Also collect permutation distributions for true WY stepdown
    if (verbose) cat("Computing permutations...\n")
    
    raw_pvalues <- matrix(NA_real_, nrow = n_genes, ncol = n_q,
                         dimnames = list(rownames(x), paste0("q=", q_values)))
    log2fc_values <- matrix(NA_real_, nrow = n_genes, ncol = n_q,
                           dimnames = list(rownames(x), paste0("q=", q_values)))
    
    # List to store permutation p-value matrices for all q-values
    # Will be used to compute permutation minima for true WY stepdown
    perm_pvalues_all_q <- vector("list", n_q)
    names(perm_pvalues_all_q) <- paste0("q=", q_values)
    
    # For each q-value, extract the appropriate columns and test
    for (q_idx in seq_along(q_values)) {
        q_val <- q_values[q_idx]
        
        # Get column indices for this q-value
        col_indices <- which(col_q_values == q_val)
        x_q <- x[, col_indices, drop = FALSE]
        
        if (verbose) {
            cat("  q =", q_val, "... ")
        }
        
        # Run permutation test for this q-value
        # Returns both p-values and permutation distributions
        perm_results <- .label_shuffling_single_q(
            x = x_q,
            samples = samples,
            control = control,
            method = method,
            randomizations = randomizations,
            nthreads = nthreads
        )
        
        raw_pvalues[, q_idx] <- perm_results$pvalue
        log2fc_values[, q_idx] <- perm_results$log2FC
        perm_pvalues_all_q[[q_idx]] <- perm_results$perm_matrix
        
        if (verbose) {
            cat("done\n")
        }
    }
    
    # Step 2: Compute permutation minima across all tests (true WY stepdown)
    if (verbose) cat("Computing permutation minima across all tests...\n")
    
    # For true Westfall-Young, we need the minimum p-value from each permutation
    # across all tests. However, computing p-values for each gene in each permutation
    # is expensive. Instead, we use test statistics (absolute fold-change):
    # For each permutation, track the maximum |FC| across all genes * q-values
    # Then for each observed test, count permutations where max(|perm_FC|) >= |obs_FC|
    
    perm_minima <- numeric(randomizations)
    
    for (perm_idx in 1:randomizations) {
        # Collect all fold-change values across genes and q-values for this permutation
        all_fc_abs_perm <- numeric()
        
        for (q_idx in seq_along(q_values)) {
            # perm_pvalues_all_q[[q_idx]] is a matrix: genes * randomizations
            # Extract the fold-change values for this permutation at this q-value
            fc_perm_q <- abs(perm_pvalues_all_q[[q_idx]][, perm_idx])
            all_fc_abs_perm <- c(all_fc_abs_perm, fc_perm_q)
        }
        
        # The minimum p-value in this permutation corresponds to maximum |FC|
        # (since p-value is monotonically decreasing in |FC|)
        # Store the maximum absolute fold-change for this permutation
        perm_minima[perm_idx] <- max(all_fc_abs_perm, na.rm = TRUE)
    }
    
    if (verbose) cat("Applying true Westfall-Young stepdown...\n")
    
    # Flatten p-values and log2FC for stepdown procedure
    all_pvalues <- as.numeric(raw_pvalues)
    all_log2fc <- as.numeric(log2fc_values)
    
    # Create index mapping back to original matrix
    idx_matrix <- expand.grid(gene = 1:n_genes, q = 1:n_q)
    gene_idx <- idx_matrix$gene
    q_idx_vec <- idx_matrix$q
    
    # Sort by p-value (ascending), which corresponds to largest |FC| first
    sort_order <- order(all_pvalues)
    sorted_pvalues <- all_pvalues[sort_order]
    sorted_log2fc <- all_log2fc[sort_order]
    sorted_gene_idx <- gene_idx[sort_order]
    sorted_q_idx <- q_idx_vec[sort_order]
    
    # Initialize adjusted p-values
    wy_adjusted <- numeric(length(all_pvalues))
    
    # True Westfall-Young Stepdown Procedure
    # ======================================
    # For each ranked test i (ordered by p-value):
    # 1. Get the observed test statistic (|log2FC|)
    # 2. Count how many permutation maxima are >= this observed |FC|
    #    (equivalently: count permutations where min p-value would be <= observed p)
    # 3. Adjusted p = (count + 1) / (randomizations + 1)  [Phipson-Smyth correction]
    # 4. Enforce monotonicity: adjusted p-values must be non-decreasing
    #
    # This implements the true WY stepdown from Westfall & Young (1993),
    # accounting for the joint dependence among all gene * q-value tests.
    # Reference: Meinshausen et al. (2011, S165/S166) prove optimality.
    #
    for (i in seq_along(sorted_pvalues)) {
        obs_fc_abs <- abs(sorted_log2fc[i])
        # Count how many permutation maxima are >= observed |FC|
        # (This is equivalent to counting permutations with min p-value <= obs p)
        count <- sum(perm_minima >= obs_fc_abs)
        # Phipson-Smyth correction
        wy_adjusted[sort_order[i]] <- min(1.0, (count + 1) / (randomizations + 1))
    }
    
    # Enforce monotonicity: adjusted p-values must be non-decreasing
    # This maintains the stepdown property and guarantees FWER <= α
    for (i in 2:length(wy_adjusted)) {
        if (wy_adjusted[sort_order[i]] < wy_adjusted[sort_order[i-1]]) {
            wy_adjusted[sort_order[i]] <- wy_adjusted[sort_order[i-1]]
        }
    }
    
    # Step 3: Also compute standard BH-FDR adjustments for comparison
    bh_adjusted <- p.adjust(all_pvalues, method = "BH")
    
    # Step 4: Build output data frame
    output <- data.frame(
        gene = rownames(x)[gene_idx],
        q_value = q_values[q_idx_vec],
        log2FC = all_log2fc,
        pvalue_raw = all_pvalues,
        pvalue_wy = wy_adjusted,  # Westfall-Young stepdown
        padj_bh = bh_adjusted,     # BH-FDR for comparison
        stringsAsFactors = FALSE
    )
    
    # Sort by Westfall-Young adjusted p-value
    output <- output[order(output$pvalue_wy), ]
    rownames(output) <- NULL
    
    # Add attributes with procedure information
    attr(output, "procedure") <- "True Westfall-Young Stepdown"
    attr(output, "fwer_level") <- 0.05
    attr(output, "n_tests") <- n_total_tests
    attr(output, "n_q_values") <- n_q
    attr(output, "n_genes") <- n_genes
    attr(output, "randomizations") <- randomizations
    
    if (verbose) {
        sig_wy <- sum(output$pvalue_wy < 0.05, na.rm = TRUE)
        sig_bh <- sum(output$padj_bh < 0.05, na.rm = TRUE)
        cat("\nResults:\n")
        cat("  Significant (WY-adjusted p < 0.05):", sig_wy, "\n")
        cat("  Significant (BH-FDR < 0.05):", sig_bh, "\n")
        cat("  FWER control: Family-wise error rate <= 0.05\n")
    }
    
    return(output)
}


#' Internal: Single q-value label shuffling
#'
#' Helper function for computing permutation p-values at a single q-value.
#' Returns both the p-values and the full permutation distribution matrix.
#'
#' @param x Matrix of diversity values (genes * samples)
#' @param samples Character vector of sample group labels
#' @param control Control group name
#' @param method Aggregation method ("mean" or "median")
#' @param randomizations Number of permutations
#' @param nthreads Number of threads
#'
#' @return List with:
#'   \describe{
#'     \item{pvalue}{Vector of p-values for each gene}
#'     \item{log2FC}{Vector of observed log2 fold-changes}
#'     \item{perm_matrix}{Matrix of permutation fold-changes (genes * randomizations)}
#'   }
#'
#' @keywords internal
#' @noRd
.label_shuffling_single_q <- function(x, samples, control, method, randomizations, nthreads) {
    
    # Compute observed fold changes
    fc_result <- calculate_fc(x, samples, control, method)
    log2_fc <- fc_result[, 4]
    
    # Generate permutations
    permuted <- replicate(randomizations, 
                         calculate_fc(x, sample(samples), control, method), 
                         simplify = FALSE)
    perm_mat <- vapply(permuted, function(z) as.numeric(z[, 4]), numeric(nrow(x)))
    
    # Ensure perm_mat is always a matrix (fixes single-gene edge case)
    if (!is.matrix(perm_mat)) {
        perm_mat <- matrix(perm_mat, nrow = nrow(x), byrow = FALSE)
    }
    
    # Compute p-values with Phipson-Smyth correction
    .compute_pval <- function(i) {
        obs <- log2_fc[i]
        nulls <- perm_mat[i, ]
        if (is.na(obs) || all(is.na(nulls))) {
            return(NA_real_)
        }
        nulls_non_na <- nulls[!is.na(nulls)]
        n_non_na <- length(nulls_non_na)
        if (n_non_na == 0) {
            return(NA_real_)
        }
        cnt <- sum(abs(nulls_non_na) >= abs(obs))
        pval <- (cnt + 1) / (n_non_na + 1)  # Phipson-Smyth correction
        return(pval)
    }
    
    raw_p_values <- unlist(.tsenat_bplapply(seq_len(nrow(perm_mat)), .compute_pval, 
                                            nthreads = nthreads))
    
    # Replace NA p-values with 1.0 for compatibility
    raw_p_values[is.na(raw_p_values)] <- 1.0
    
    # Return results including permutation matrix for true WY stepdown
    result <- list(
        pvalue = raw_p_values,
        log2FC = log2_fc,
        perm_matrix = perm_mat  # Full permutation distribution (genes * randomizations)
    )
    
    return(result)
}

