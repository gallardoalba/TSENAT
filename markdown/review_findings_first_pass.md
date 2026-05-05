# TSENAT R Folder Review – First Pass Findings

Date: 2026-07-21

## Scope
This note summarizes the first-pass static review of the R scripts in the TSENAT package. The goal was to identify likely bugs, dangerous fallbacks, and maintainability problems visible from source inspection.

## High-confidence findings

### 1. Duplicate helper definition: `.get_gene_ids()` is defined twice

Evidence:
- [R/se_manipulation_build.R](R/se_manipulation_build.R#L220-L247)
- [R/se_manipulation_filter.R](R/se_manipulation_filter.R#L197-L223)

Issue:
The same internal helper is implemented in two different files with different behavior. In R, the last loaded definition wins. This can make behavior dependent on source order and is a common cause of hidden regressions.

Impact:
- Non-deterministic behavior across environments
- Hard-to-debug downstream inconsistencies
- Increased maintenance burden

Recommendation:
Keep one canonical implementation and remove the duplicate.

### 2. Silent fallback in assay selection can mask incorrect inputs

Evidence:
- [R/se_manipulation_filter.R](R/se_manipulation_filter.R#L149-L160)

Issue:
`.resolve_assay_index()` silently falls back to assay 1 whenever an assay name is unknown or an out-of-range numeric index is provided.

Impact:
- Bad input can quietly target the wrong assay
- The user may get a result without realizing the configuration is wrong
- Scientific conclusions can be silently altered

Recommendation:
Only allow fallback when the parameter is omitted; otherwise, raise an error.

### 3. Constructor mutates the caller's `SummarizedExperiment` in place

Evidence:
- [R/s4_class.R](R/s4_class.R#L33-L76)

Issue:
`TSENATAnalysis()` appends `sample_id` to `colData(se)` when missing. This mutates the input object directly.

Impact:
- Side effects on user-managed objects
- Harder debugging and reproducibility
- Unexpected changes before analysis logic begins

Recommendation:
Avoid mutating the incoming object. Build a modified copy or create the new metadata in a local object.

### 4. Auto-detection heuristics may silently choose the wrong group/control

Evidence:
- [R/divergence_helpers.R](R/divergence_helpers.R#L1-L97)

Issue:
The group and control auto-detection logic relies on hard-coded column priority and a heuristic “smallest group count” fallback.

Impact:
- A wrong group column or reference group can be selected silently
- This directly affects divergence/contrast interpretation

Recommendation:
Use these heuristics only when no ambiguity exists. Otherwise, stop and request explicit parameters.

## Medium-confidence findings

### 5. Too much best-effort fallback across the API

Evidence:
- [R/se_manipulation_filter.R](R/se_manipulation_filter.R#L149-L160)
- [R/divergence_helpers.R](R/divergence_helpers.R#L1-L97)

Issue:
Several helpers silently recover from invalid or ambiguous inputs instead of failing fast.

Impact:
- Lower discoverability of configuration mistakes
- Greater risk of incorrect but non-crashing results

Recommendation:
Prefer explicit validation and early errors for user-facing analysis parameters.

## Second pass findings

### 6. The orchestration layer converts major step failures into warnings and keeps going

Evidence:
- [R/orchestration.R](R/orchestration.R#L654-L673)
- [R/orchestration.R](R/orchestration.R#L700-L708)
- [R/orchestration.R](R/orchestration.R#L716-L727)
- [R/orchestration.R](R/orchestration.R#L755-L766)

Issue:
The pipeline helper functions wrap many analytical steps in `tryCatch()` and then emit a warning instead of stopping. The workflow can continue with missing intermediate results, so the run may complete with a partially populated object and no hard failure.

Impact:
- Partial analysis artifacts can be invisible to the user
- Later steps may fail unpredictably or compare results from different stages
- This is a high-risk design for a scientific workflow package

Recommendation:
Use a strict mode for essential pipeline stages, and either fail fast or record a clearly marked step status that prevents downstream steps from running on incomplete state.

### 7. The package has a wide pattern of “best-effort” downstream recovery that can hide data problems

Evidence:
- [R/orchestration.R](R/orchestration.R#L654-L892)
- [R/sait_helpers.R](R/sait_helpers.R#L1500-L1537)

Issue:
There are several fallback paths that keep the pipeline alive even when data mapping or preprocessing is incomplete. In some places, the code falls back to using gene IDs or unpaired tests without a clear signal to the caller.

Impact:
- The package may return a result that looks valid but is built on a degraded preprocessing path
- Reproducibility and auditing become harder

Recommendation:
Keep fallbacks for visualization metadata only, and use hard errors for statistical or mapping failures that would change the biological interpretation.

## Hardening applied

The highest-risk silent-recovery behaviors are now hardened in source:

1. Fail-fast assay selection
- [R/se_manipulation_filter.R](R/se_manipulation_filter.R#L149-L167)
- `.resolve_assay_index()` now raises an explicit error for unknown assay names and out-of-range or malformed assay references instead of silently selecting assay 1.

2. Non-mutating constructor behavior
- [R/s4_class.R](R/s4_class.R#L54-L76)
- `TSENATAnalysis()` now copies the incoming `SummarizedExperiment` before adding `sample_id`, so the constructor no longer modifies the caller’s `colData` object in place.

3. Strict orchestration failure handling
- [R/orchestration.R](R/orchestration.R#L654-L823)
- Critical downstream steps now stop on error rather than continuing with a partial pipeline state and only issuing a warning.

4. SAIT fail-fast propagation
- [R/s4_functions_sait.R](R/s4_functions_sait.R#L183-L223)
- [R/sait_core.R](R/sait_core.R#L300-L320)
- The SAIT wrapper now propagates real computation errors directly instead of converting them into a warning-only empty result object.

## Verification

The hardening was validated with targeted test runs:

- `Rscript -e "library(testthat); devtools::load_all('.'); test_file('tests/testthat/test-core-functions-filter.R', reporter='summary')"`
  - Result: `══ DONE ════════════════════════════════════════════════════════════════════════` with no failed tests.
- `Rscript -e "library(testthat); devtools::load_all('.'); test_file('tests/testthat/test-core-functions-s4.R', reporter='summary')"`
  - Result: `══ DONE ════════════════════════════════════════════════════════════════════════` with no failed tests.
- `Rscript -e "library(testthat); devtools::load_all('.'); test_file('tests/testthat/test-core-functions-s4-sait.R', reporter='summary')"`
  - Result: `══ DONE ════════════════════════════════════════════════════════════════════════` with no failed tests.

### 8. SAIT result handling can silently degrade a real failure into an empty analysis result

Evidence:
- [R/s4_functions_sait.R](R/s4_functions_sait.R#L183-L223)
- [R/s4_functions_sait.R](R/s4_functions_sait.R#L403-L426)
- [R/sait_core.R](R/sait_core.R#L300-L320)

Issue:
The public SAIT wrapper (`calculate_sait()`) wraps the core computation in `tryCatch()` and converts any error into a warning plus an empty `data.frame()`. The downstream validation helper then emits another warning and returns an empty `results` slot instead of stopping the workflow. The object can therefore remain “successful” from the caller’s point of view while carrying no scientific signal.

Impact:
- A genuine model/setup failure can be misinterpreted as a valid no-signal result
- Users may continue downstream steps on an object that is structurally empty
- This is a high-risk silent recovery pattern for a scientific analysis API

Recommendation:
Prefer a fail-fast policy for user- or config-triggered problems, and reserve empty-result behavior only for well-defined, explicitly documented no-signal situations.

## Statistical Analysis Findings (Third Pass - Deep Evaluation)

### 9. ADF test p-value uses t-distribution instead of Dickey-Fuller distribution [CONFIRMED BUG]

Evidence:
- [R/sait_helpers.R](R/sait_helpers.R#L536-L537)

Issue:
The Augmented Dickey-Fuller test calculates p-values using `2 * pt(t_stat, df = ...)`, which assumes a standard t-distribution. The ADF test statistic follows the Dickey-Fuller distribution, which is non-standard and depends on sample size and regression specification. Using the t-distribution produces incorrect p-values.

The code comments acknowledge this: "Approximate p-value based on t-statistic comparison to critical value. This is a simplified approximation; exact p-values require special distribution."

Impact:
- Unit root tests may give wrong conclusions about stationarity
- ARIMA(1,1,0) differencing decisions may be based on incorrect statistical evidence
- **Severity: Medium** - The function is used for diagnostic reporting, not for automated decision-making in the main pipeline. The `stationary` boolean decision uses critical values (not p-values), which are correctly approximated.

Recommendation:
Use a proper Dickey-Fuller distribution table or implement the Phillips-Perron test via the `tseries` package. At minimum, the existing documentation should explicitly state that p-values are approximate and should not be used for formal hypothesis testing.

### 10. ADF test named as "augmented" but implements simple Dickey-Fuller [CONFIRMED - NAMING ISSUE]

Evidence:
- [R/sait_helpers.R](R/sait_helpers.R#L453-L545)

Issue:
The function is named `.adf_test` and the comments describe an "augmented" ADF test with lagged differences (Δy_{t-i}), but the implementation only uses y_{t-1} as a predictor without any lagged difference terms. Line 494 explicitly states: "Simple regression: just use y_{t-1} without augmentation for stability."

Impact:
- **Severity: Low** - The simple DF test is appropriate for short time series (typical q-value sequences have 5-10 observations). Augmentation is mainly needed to handle serial correlation in longer series. The function works correctly for its intended use case.

Recommendation:
Either rename to `.df_test()` to reflect the implementation, or add optional lag augmentation. The current naming is misleading but the function is functionally correct for its use case.

### 11. KPSS test uses crude p-value approximation [CONFIRMED - DOCUMENTED LIMITATION]

Evidence:
- [R/sait_helpers.R](R/sait_helpers.R#L617-L629)

Issue:
The KPSS p-value is computed using fixed critical values without interpolation, resulting in a step function: p ∈ {0.1, 0.05, 0.025, 0.01, 0.001} for the reject case, and 1-p for the non-reject case.

Impact:
- **Severity: Low** - This is a diagnostic function for exploratory analysis. The boolean `stationary` decision is based on the critical value comparison (which is correct), not the p-value. The p-value is used only for reporting.

Recommendation:
This is an acceptable trade-off for avoiding external dependencies. The function documentation should note the crude approximation. For formal testing, users should use the `tseries` package.

### 12. Breusch-Pagan statistic formula [NOT A BUG - VERIFIED CORRECT]

Evidence:
- [R/sait_helpers.R](R/sait_helpers.R#L886-L893)

Analysis:
The code computes:
- `fitted_sq_sum = sum((fitted(fit_aux) - mean(fitted(fit_aux)))^2)` = Regression Sum of Squares (SSR)
- `rss_aux = sum(residuals(fit_aux)^2)` = Residual Sum of Squares (RSS)
- `tss_aux = fitted_sq_sum + rss_aux` = Total Sum of Squares (TSS)
- `bp_stat = (fitted_sq_sum/tss_aux) * nrow(df)` = `n * (SSR/TSS)` = `n * R²`

**This IS the correct Breusch-Pagan formula**: BP = n × R² from the auxiliary regression of squared residuals on predictors. The chi-squared p-value calculation with appropriate degrees of freedom is also correct.

Conclusion: **This is NOT a bug.** The implementation follows the standard BP test correctly.

### 13. AR(1) design effect uses loop instead of closed-form [PERFORMANCE ISSUE, NOT A BUG]

Evidence:
- [R/sait_helpers.R](R/sait_helpers.R#L107-L119)

Issue:
The design effect is computed using a loop that accumulates `summed += lambda_k * (rho^k)`. The comments mention closed-form formulas but don't use them.

Analysis:
- For the typical use case (cluster_size = 5-10 q-values, rho < 0.9), the loop is numerically stable
- The formula used `D_eff = 1 + 2 * Σ (1 - k/m) * ρ^k` is mathematically correct for AR(1) with edge effects
- The closed-form `D_eff = (1 + phi) / (1 - phi)` is only valid for m → ∞, not for small m

Impact:
- **Severity: Very Low** - Numerical precision is adequate for typical use cases. The loop is O(m) which is fine for m ≤ 50.

Recommendation:
This is acceptable as-is. For very large cluster sizes (>100), consider switching to the asymptotic formula.

### 14. Concurvity function computes model complexity, not concurvity [CONFIRMED - MISLEADING NAME]

Evidence:
- [R/srh_helpers.R](R/srh_helpers.R#L457-L589)

Issue:
The function `.compute_concurvity_index()` is documented as "Detects collinearity among smooth terms" and the return value is named `overall_concurvity`, but the implementation computes "model complexity" via effective degrees of freedom (EDF) from GAM smooths. True concurvity measures how well one smooth term can be approximated by other smooth terms.

Impact:
- **Severity: Medium** - Users expecting concurvity diagnostics will get complexity metrics instead. This could lead to missed model diagnostics.

Recommendation:
Either:
1. Implement true concurvity using `mgcv::concurvity()` (preferred), or
2. Rename the function to `.compute_model_complexity()` and update documentation to clarify what is measured.

### 15. Westfall-Young permutation shuffles within q-levels [DESIGN CHOICE - NOT NECESSARILY A BUG]

Evidence:
- [R/sait_helpers.R](R/sait_helpers.R#L1402-L1412)

Issue:
The Westfall-Young permutation shuffles group labels separately within each q-level rather than permuting subject labels.

Analysis:
- The permutation is applied to test the q × condition interaction, not to preserve AR(1) structure
- By shuffling within each q-level, the null hypothesis is that group labels are independent of the response at each q-level
- This is a valid permutation strategy for testing interaction effects
- The AR(1) correlation is between q-values, not between subjects, so permuting subject labels would not be appropriate

Impact:
- **Severity: Low** - The permutation strategy is valid for the hypothesis being tested. The documentation in `westfall_young_permutation.R` correctly describes the approach.

Conclusion: **This is NOT a bug.** The permutation strategy is appropriate for testing q × condition interactions.

## Summary of Validated Findings

| # | Finding | Status | Severity |
|---|---------|--------|----------|
| 9 | ADF p-value uses t-distribution | **CONFIRMED BUG** | Medium |
| 10 | ADF named as "augmented" but implements DF | **CONFIRMED - NAMING ISSUE** | Low |
| 11 | KPSS crude p-value approximation | **CONFIRMED - DOCUMENTED LIMITATION** | Low |
| 12 | Breusch-Pagan formula | **NOT A BUG** | N/A |
| 13 | AR(1) design effect loop | **PERFORMANCE ISSUE** | Very Low |
| 14 | Concurvity misnamed | **CONFIRMED - MISLEADING NAME** | Medium |
| 15 | Westfall-Young permutation | **NOT A BUG** | N/A |

## Notes
The third pass review identified 2 confirmed bugs (findings 9, 14), 2 naming/documentation issues (findings 10, 11), and 2 non-issues (findings 12, 15). The findings with Medium severity (9, 14) should be addressed to improve statistical rigor and user trust. The Low severity items are acceptable trade-offs for avoiding external dependencies and maintaining simplicity.
