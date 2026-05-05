# Upgrade Conover-Iman Rank Transform → Aligned Rank Transform (ART)

## Summary

Replace the current Conover-Iman Rank Transform implementation in `calculate_rank_transform()` with the **Aligned Rank Transform (ART)** using the `ARTool` R package for proper non-parametric interaction testing in factorial designs.

## Background

### Current Implementation

The package currently uses the **Conover-Iman Rank Transform** (Conover & Iman, 1981) for non-parametric Q×Condition interaction testing. The implementation in `.test_q_condition_interaction()` (`R/rank_transform_helpers.R`) works as follows:

```r
# Rank the data
data$ranks <- rank(data$entropy)

# Fit standard linear model on ranks
sait_model <- lm(ranks ~ q * condition, data = data)
anova_result <- anova(sait_model)

# Extract interaction F-statistic and p-value
f_stat <- anova_result$`F value`[interaction_row]
p_val  <- anova_result$`Pr(>F)`[interaction_row]
```

### Why This Needs Upgrading

While the Conover-Iman rank transform is a valid non-parametric procedure for **main effects**, it has a well-documented limitation for **interaction testing**:

- **Rank transformation does not preserve interaction structure.** When data is ranked and then a linear model with interactions is fitted, the interaction term can be distorted because ranking is a nonlinear transformation that affects the additive decomposition of effects.
- The F-test on rank-transformed data for interactions can be **anti-conservative** (inflated Type I error) or **conservative** (reduced power), depending on data structure and design balance.
- Sawilowsky (1990) and Headrick & Sawilowsky (2000) demonstrated these failures in multiple scenarios.

### The Solution: Aligned Rank Transform (ART)

The **Aligned Rank Transform (ART)** (Higgins & Tashtoush, 1994; Wobbrock et al., 2011) is the state-of-the-art non-parametric method for factorial designs with interactions:

1. **Align**: Remove the effect of all other factors before ranking — strips out main effects so the interaction is isolated
2. **Rank**: Rank the aligned residuals
3. **Test**: Apply standard ANOVA to the aligned ranks

The ART properly handles interactions by:
- Ensuring the interaction term is valid after alignment
- Maintaining Type I error control across balanced and unbalanced designs
- Providing a unified framework for main effects AND interactions

The R package **`ARTool`** (Kay et al., 2021, CRAN: https://cran.r-project.org/package=ARTool) provides a well-tested, peer-reviewed implementation with:
- Support for fixed and random effects
- Post-hoc contrast tests
- Effect size computation
- Proper handling of repeated measures

## Proposed Changes

### Files to Modify

| File | Change |
|------|--------|
| `R/rank_transform_helpers.R` | Replace `.test_q_condition_interaction()` implementation to use ART via `ARTool` |
| `R/rank_transform_core.R` | Update `.detect_q_analyze_gene()` to work with ART output format |
| `R/s4_functions_rank_transform.R` | Update wrapper documentation |
| `DESCRIPTION` | Add `ARTool` to `Imports` or `Suggests` |
| `tests/testthat/test-statistical-methods-rank.R` | Update test expectations for ART output |
| `man/calculate_rank_transform.Rd` | Update documentation |

### Key Implementation Details

1. **ART formula**: `art(entropy ~ q * condition + Error(subject/q))` for paired designs, `art(entropy ~ q * condition)` for unpaired
2. **ANOVA extraction**: Use `anova.art()` to get the ART ANOVA table with proper F-values
3. **Effect sizes**: ART-provided partial η² should replace the current manually-computed η²
4. **Backward compatibility**: Consider keeping the existing function name `.test_q_condition_interaction()` but replacing the internals, or add a new method parameter (e.g., `method = c("art", "conover-iman")`)
5. **Optional dependency**: Consider `ARTool` in `Suggests` with a fallback to the current Conover-Iman implementation if not installed

### Pseudo-code for New Implementation

```r
.test_q_condition_interaction <- function(data, value_col, q_col, condition_col,
                                           paired, subject_col, method = "art") {
    if (method == "art") {
        if (!requireNamespace("ARTool", quietly = TRUE)) {
            warning("ARTool not installed, falling back to Conover-Iman")
            method <- "conover-iman"
        }
    }
    
    if (method == "art") {
        # Ensure factors
        data[[q_col]] <- factor(data[[q_col]])
        data[[condition_col]] <- factor(data[[condition_col]])
        
        # Fit ART model
        if (paired && !is.null(subject_col)) {
            data[[subject_col]] <- factor(data[[subject_col]])
            art_formula <- as.formula(
                paste(value_col, "~", q_col, "*", condition_col, 
                      "+ Error(", subject_col, "/", q_col, ")"))
        } else {
            art_formula <- as.formula(
                paste(value_col, "~", q_col, "*", condition_col))
        }
        
        art_model <- ARTool::art(art_formula, data = data)
        art_anova <- anova(art_model)
        
        # Extract interaction row
        interaction_term <- paste0(q_col, ":", condition_col)
        int_row <- art_anova[interaction_term, ]
        
        return(list(
            statistic = int_row[["F value"]],
            p_value   = int_row[["Pr(>F)"]],
            method    = "Aligned Rank Transform (ART)",
            test_type = if (paired) "art_paired" else "art_unpaired",
            ...
        ))
    }
    # ... existing Conover-Iman fallback ...
}
```

## References

- Conover, W. J., & Iman, R. L. (1981). Rank transformations as a bridge between parametric and nonparametric statistics. *The American Statistician*, 35(3), 124–129.
- Higgins, J. J., & Tashtoush, S. (1994). An aligned rank transform test for interaction. *Nonlinear World*, 1(2), 201–211.
- Wobbrock, J. O., Findlater, L., Gergle, D., & Higgins, J. J. (2011). The aligned rank transform for nonparametric factorial analyses using only ANOVA procedures. *CHI 2011*, 143–146.
- Kay, M., Elkin, L. A., Higgins, J. J., & Wobbrock, J. O. (2021). ARTool: Aligned Rank Transform for Nonparametric Factorial ANOVAs. R package version 0.11.1. https://cran.r-project.org/package=ARTool
- Sawilowsky, S. S. (1990). Nonparametric tests of interaction in experimental design. *Review of Educational Research*, 60(1), 91–126.

## Priority & Effort Estimate

| Aspect | Detail |
|--------|--------|
| **Priority** | Medium — current Conover-Iman is valid non-parametric, ART is an improvement |
| **Effort** | ~1–2 days for implementation + testing |
| **Risk** | Low if ARTool is in `Suggests` with fallback; Medium if it becomes a hard dependency |
| **Breaking change** | Minimal if `method` parameter is added; tests will need updating |

## Related Issues

- I3: Conover-Iman mislabeled as Scheirer-Ray-Hare (already fixed — correct labeling in place)
- I5: η² computed on ranks (already fixed — now on original entropy scale)
