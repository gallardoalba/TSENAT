#!/usr/bin/env bash
# Wrapper to run the R coverage script.
set -euo pipefail
Rscript scripts/run_coverage.R "$@"
