#!/usr/bin/env bash
# Wrapper to run the R coverage script.
set -euo pipefail
Rscript tools/run_coverage.R "$@"
