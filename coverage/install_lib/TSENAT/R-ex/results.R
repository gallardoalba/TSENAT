### Name: results
### Title: Extract analysis results from TSENATAnalysis object
### Aliases: results

### ** Examples

# Load example data
data(readcounts, package = 'TSENAT')

# Create TSENATAnalysis from count matrix
config <- TSENAT_config(
  q = 1.0,
  condition_col = 'group'
)
se <- SummarizedExperiment::SummarizedExperiment(
  assays = list(counts = readcounts),
  colData = data.frame(
    group = rep(c('A', 'B'), length.out = ncol(readcounts))
  )
)
analysis <- TSENATAnalysis(se = se, config = config)
analysis <- calculate_diversity(analysis)

# Get all diversity results (list of SummarizedExperiment objects, one per q)
div_all <- results(analysis, type = 'diversity')

# Get diversity for specific q-value (single SummarizedExperiment table for display)
div_q1 <- results(analysis, type = 'diversity', q = 1.0)

# Get diversity as SummarizedExperiment for downstream processing
div_q1_se <- results(analysis, type = 'diversity', q = 1.0, format = 'se')

# Get results ranked by p-value, top 20 genes
# Using accessor function instead of @ slot access
top_sait <- results(analysis, type = "sait", rankBy = 'pvalue', n = 20)

# Get effect size results, top 6 genes by p-value (most significant first)
top_effect_sizes <- results(analysis, type = 'effect_sizes_divergence', 
                             top_n = 6, sort_by = 'adj_p_interaction')

# Get effect sizes sorted by mean divergence (largest effect sizes first)
large_effects <- results(analysis, type = 'effect_sizes_divergence',
                         top_n = 10, sort_by = 'Mean_Divergence')

# Get switching tables - automatically computed if prerequisites exist
# (no need to call prepare_gene_switching_tables_s4 separately)
# Default format='text' returns structured list for vignette rendering
switching <- results(analysis, type = 'switching_tables')

# Get switching tables in raw format (named list of data frames) for direct manipulation
switching_raw <- results(analysis, type = 'switching_tables', format = 'raw')

# ========================================================================
# RETRIEVE CACHED PLOTS using the plot parameter
# ========================================================================

# Get the diversity spectrum plot
diversity_plot <- results(analysis, type = 'diversity', plot = TRUE)

# Get the LM/regularized regression (GAM, LMM, GEE, FPCA) interaction plot
sait_plot <- results(analysis, type = "sait", plot = TRUE)

# Get the divergence distribution plot
div_dist_plot <- results(analysis, type = 'divergence', plot = TRUE)

# Get the influence/m-estimator plot
influence_plot <- results(analysis, type = 'influence', plot = TRUE)

# Available plot types correspond to analysis types:
# - type = 'diversity': Returns Tsallis entropy q-spectrum visualization
# - type = "sait": Returns regularized/penalized regression (GAM, LMM, GEE, FPCA) interaction plot
# - type = 'divergence': Returns distribution of divergence metrics across genes
# - type = 'influence': Returns m-estimator sample influence analysis
# - type = 'rank_test': Returns Conover-Iman Rank Transform interaction results
# - type = 'concordance': Returns method concordance comparison (LM vs rank test)




