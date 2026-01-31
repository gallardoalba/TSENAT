# Helper to install package dependencies for development/CI
if (!requireNamespace("remotes", quietly = TRUE)) install.packages("remotes")
# Install dependencies (Imports, Depends, Suggests)
remotes::install_deps(dependencies = TRUE)
