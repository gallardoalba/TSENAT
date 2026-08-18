### Name: TSENAT_config
### Title: Create and return TSENAT configuration
### Aliases: TSENAT_config

### ** Examples

# Default config with standard parameters (point estimates only)
cfg <- TSENAT_config()

# For Wilcoxon/shuffle tests (single q-value required in config)
cfg <- TSENAT_config(
  q = 1.0,                          # Shannon entropy - for rank tests
  condition_col = 'treatment',
  control = 'untreated'
)

# For Conover-Iman Rank Transform tests (multiple q-values)
cfg <- TSENAT_config(
  q = seq(0, 2, by = 0.5),          # Multiple q-values for spectrum or advanced testing
  condition_col = 'treatment',
  control = 'untreated'
)

# With bootstrap CIs for uncertainty quantification (recommended)
cfg <- TSENAT_config(
  bootstrap = TRUE,                # Enable bootstrap confidence intervals
  bootstrap_method = 'bca',         # Bias-corrected (better for skewed entropy)
  nboot = 1000,                     # 1000 resamples
  bootstrap_ci = 0.95               # 95% CI
)

# Custom with paired analysis, strict filtering, and normalization
cfg <- TSENAT_config(
  q = 1.0,                          # Shannon entropy
  condition_col = 'treatment',
  subject_col = 'subject_id',
  paired = TRUE,
  control = 'untreated',
  stringency = 'severe',            # High-confidence transcripts only
  norm_method = 'zscore',           # Cross-study standardization
  shrinkage = 'none',               # Empirical estimates
  bootstrap = TRUE,
  bootstrap_method = 'bca',
  nboot = 5000,                     # Higher precision
  pseudocount = 0,                  # Disabled by default; set > 0 to add pseudocount
  significance_threshold = 0.01     # Stricter significance level
)




