# PR Summary: Comprehensive Statistical Bug Audit & Fix (July 2026)

## Description

Comprehensive audit and fix of statistical implementation bugs across 15
R source files, identified via systematic code review with numerical
verification (`bugs_audit.txt`). Covers core statistical methods
(Hochberg correction, GEE bias-correction, bootstrap CIs, AR(1)
estimation, LMM/GEE model validation) plus jackknife diagnostics,
M-estimation, paired bootstrap, Storey π₀, Westfall-Young permutation,
and entropy/divergence helpers.

Also includes supporting code refactoring: pipeline stage extraction,
input validation helpers, and modular assumption checks.

**Related Issues:** N/A (proactive audit)

------------------------------------------------------------------------

## Checklist

Run `R CMD check` and ensure no errors

Add tests for new behavior (12 test files enhanced, ~800+ new test
lines)

Update `NEWS.md` if relevant

------------------------------------------------------------------------

## Bugs Fixed

### Critical

| ID | File | Issue | Fix |
|----|----|----|----|
| **C1** | `rank_transform_helpers.R` | `.hochberg_stepup` used Holm (`cummax`) instead of Hochberg step-up | `rev(cummin(rev(...)))` — now matches `p.adjust(method="hochberg")` |
| **C2** | `sait_gee.R` | GEE bias-correction: `<20` vs `<30` threshold mismatch between code and documentation | Unified to `<30` clusters threshold |
| **C3** | `sait_gee.R` | KC/HC1 multiplier computed but never applied (`vcov=NULL` path dropped result) | Now passes `vcov(fit_alt)` + `coef_value` through correction pipeline |
| **C4** | `bootstrap.R` | Bootstrap CIs destroyed by cross-gene normalization applied before CI computation | CIs remain on raw (per-gene) scale; normalization deferred |
| **C5** | `jackknife_diagnostics.R` | Tsallis (q≠1) jackknife estimates divided by `log(log_base)`, putting them on a different scale than the point estimate from `.entropy_single()` — spurious influence values when `log_base ≠ exp(1)` | Removed `/log(log_base)` from `.jackknife_compute_estimates()` and `.jackknife_process_vector_core()`. Tsallis entropy is scale-invariant. |

### High

| ID | File | Issue | Fix |
|----|----|----|----|
| **H1** | `sait_helpers.R` | `regularization="gamsel"` + `method="lmm"` → `match.arg` error (invalid combo silently accepted) | Validate `(method, regularization)` combination before dispatch |
| **H2** | `sait_gee.R` | Two competing AR(1) formulas; GEE path used wrong asymptotic form | Unified to finite-sample corrected formula; removed dead-code asymptotic path |
| **H3** | `westfall_young_permutation.R` | WY permutation broke paired structure (permuted across subjects) | Permute within subjects for paired designs; preserve subject blocking |
| **H4** | `sait_gee.R` | ARIMA `group[1]` collapsed to single factor level, losing contrast | `group[-1]` preserves both factor levels for valid design matrix |
| **H5** | `diversity_helpers.R` | `log_base` dropped in diversity matrix path (5 intermediate functions) | Threaded `log_base` parameter through all intermediate functions |
| **H6** | `divergence_helpers.R` | `method="bca"` returned percentile CIs labeled “bca” (BCa not implemented for divergence) | Warns + reports `"percentile"` in method field |
| **H7** | `divergence_helpers.R` | `d_max = log(2)` / `1` invalid for Tsallis divergence normalization | Data-driven column maximum replaces hardcoded constant |
| **H8** | `m_estimation.R` | Paired design used `df <- n_obs - 2` (two-parameter model) instead of `df <- n_obs - 1` (one-parameter location model), producing anti-conservative p-values | `df <- if (use_intercept) n_obs - 2 else n_obs - 1` |
| **H9** | `m_estimation.R` | Huber weight floored at 0.01 (`1/pmax(abs_resid, 0.01)`) instead of standard `1/\|u\|`, breaking monotonic decrease of influence (Huber 1981 §4.2) | Changed to `1/abs_resid` in both IRLS locations; `abs_resid > 1` condition already guards against division by zero |
| **H10** | `jackknife_isoform_switching.R` | Paired bootstrap silently resampled unrelated subjects in lockstep when upstream matching failed (no assertion that `ncol(A) == ncol(B)`) | Added `stopifnot(ncol(counts_A) == ncol(counts_B))` before paired resampling |
| **H11** | `jackknife_isoform_switching.R` | Per-transcript p-values pooled across all transcripts for global BH correction, while switching_status uses per-transcript CIs — two criteria can disagree with no user signal | Documented inconsistency; significance axis and switching direction axis use different tests by design |

### Medium

| ID | File | Issue | Fix |
|----|----|----|----|
| **M1** | `sait_helpers.R` | ~280 lines dead stationarity testing code (never called, unreachable) | Removed |
| **M2** | `sait_helpers.R` | LMM `min_obs=3` (internal) vs public documentation `min_obs=5` | Standardized to 5 across all paths |
| **M3** | `sait_helpers.R` | `corstr="auto"` unreachable from public API (not in `match.arg`) | Exposed in public `match.arg` choices |
| **M4** | `sait_fpca.R` | FPCA MANOVA fallback: `min(pc_pvals_valid)` anti-conservative for multiple PCs | Returns `NA_real_` when no single PC is significant |
| **M5** | `sait_gee.R` | AR(1) ρ estimated from concatenated residuals across clusters (contaminated by between-cluster covariance) | Within-cluster estimation + Fisher z-transform pooling |
| **M6** | `divergence_helpers.R` | q=0 divergence hardcoded to 0 (mathematically incorrect) | Jaccard support-difference: D₀ = \|supp(p) Δ supp(r)\| / \|supp(p) ∪ supp(r)\| |
| **M7** | `divergence_helpers.R` | Paired bootstrap silent fallback to unpaired when pairing fails | [`warning()`](https://rdrr.io/r/base/warning.html) instead of silent [`message()`](https://rdrr.io/r/base/message.html) |
| **M8** | `westfall_young_permutation.R` | Storey bootstrap π₀ selects λ by minimum CV (stability), not the canonical Storey MSE criterion — label mismatch, may bias π₀ upward | Added documentation comment; users should prefer `pi0_method = "smoother"` for canonical Storey |
| **M9** | `westfall_young_permutation.R` | Westfall-Young maxT fallback converts p-values to `-log(p)` — different scale than F-statistics, creating a latent landmine if any permutation mixes return types | Added warning comments in both parallel and serial paths |
| **M10** | `divergence_helpers.R` | `.jis_resample_paired_data()` called [`stop()`](https://rdrr.io/r/base/stop.html) inside a bootstrap loop, aborting the entire computation on one bad replicate | Changed to [`warning()`](https://rdrr.io/r/base/warning.html) + graceful return with `failed = TRUE` flag |
| **M11** | `m_estimation.R` | S-estimator scale uses fixed ±10% multipliers without bracketing — may oscillate and fail to converge within tolerance | Documented; proper bisection would require non-trivial refactoring |
| **M12** | `divergence_helpers.R` | `abs(div)` masks small-negative numerical noise but also flips genuinely negative q\<1 divergence (real signal of “p more uniform than r”) | Documented trade-off; q\<1 is uncommon in practice |

### Low

| ID | File | Issue | Fix |
|----|----|----|----|
| **L1** | `entropy_core.R` | Pseudocount contract across `.entropy_single()` → `.entropy_vectorized()` → `.entropy_core()` stack was implicit | Added explicit documentation in `.entropy_single()` |
| **L2** | `divergence_effect_sizes.R` | `.createResultsDataFrame()` used `sprintf("%.1f", q)` while `.formatMultiQResult()` used `as.character(q)`, producing different column names for integer q values | Unified to `as.character(q_values)` in both functions |
| **L3** | `jackknife_diagnostics.R` | `@exportS3method` in `@noRd` block may be silently dropped by roxygen2 | Flagged for verification with `methods::registeredS3method` |
| **L4** | `rank_transform_core.R` | eta² computed via [`lm()`](https://rdrr.io/r/stats/lm.html) without subject/block term — for paired designs this eta² is not comparable to the (paired) F-test | Added comment noting limitation |

------------------------------------------------------------------------

## Code Refactoring

### Pipeline Orchestration (`orchestration.R`, +178 lines)

Extracted monolithic pipeline into named stages: - **STAGE 1:
PREPROCESSING** — filter low-abundance transcripts - **STAGE 2:
DIVERSITY** — Tsallis entropy computation + q-spectrum plots - **STAGE
3: QC** — M-estimator influence analysis - **STAGE 4: INTERACTION** —
SAIT model fitting + visualization

### Assumptions Module (`assumptions.R`, 510-line restructure)

- Extracted helper functions: `.expand_assumptions_checks()`, rank-based
  checks, method-specific checks
- Reduced orchestrator from ~300 to ~80 lines
- Added guard against empty data inputs

### Bootstrap Validation (`bootstrap.R`, +141 lines)

- Added `.validate_bootstrap_input()` — centralized validation for NULL,
  non-numeric, zero-length, NA/Inf, negative values, all-zero edges

### GEE Module (`sait_gee.R`, 326-line restructure)

- Extracted Kauermann-Carroll bias correction into separate helper
  (Phase 7)
- Removed duplicated KC code paths
- Cleaned up phase numbering

### New Helpers

- `plots_themes.R` (new, 74 lines): Publication-ready ggplot2 theme
  helpers
- `diversity_helpers.R`, `sait_core.R`, `sait_helpers.R`: Modular helper
  additions
