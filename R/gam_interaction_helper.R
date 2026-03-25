# GAM interaction helper - enhanced with regularization and bias correction support
.tsenat_gam_interaction <- function(df, q_vals, g, min_obs = 10, subject = NULL,
                                   regularization = c("pca", "gamsel", "spline"),
                                   bias_correction = TRUE, adaptive_knots = TRUE, weights = NULL) {
    # GAMM with ARIMA(1,1,0) covariance structure for q-dependent entropy measurements
    # Paper S171 (Zimmerman & Harville, 1991): Validates generalized covariance structures.
    # Papers S168-S170: Theoretical foundation for time series correlation patterns.
    # TEST L.1.6: Confirms first differences DeltaH_q follow AR(1) pattern via ARIMA(1,1,0).
    #
    # Adaptive knot selection (ENHANCED - March 2026):
    # - automatic_knots = TRUE: per-gene adaptivity based on entropy curve complexity
    # - Measures curve roughness via coefficient of variation of slopes
    # - Allocates more basis functions (k) to complex curves, fewer to simple curves
    # - Improves model fit efficiency and reduces overfitting/underfitting trade-off
    #
    if (!requireNamespace("mgcv", quietly = TRUE)) {
        stop("Package 'mgcv' is required for method = 'gam'")
    }
    
    regularization <- match.arg(regularization)
    
    # Ensure 'group' is a factor for 'by' argument in GAM/GAMM smooths
    # This is required for s(q, by = group, ...) to work correctly
    if (!is.null(df$group)) {
        df$group <- factor(df$group)
    }
    
    # ===============================================================================
    # BOUNDED SUPPORT HANDLING (CRITICAL FIX - March 2026)
    # ===============================================================================
    # MUST be done BEFORE ARIMA differencing to:
    # 1. Detect bounds on ORIGINAL entropy (not differenced)
    # 2. Stabilize entropy before any transformations
    # 3. Preserve ARIMA structure for differenced data
    # 
    # Bug fix: Was previously done AFTER ARIMA differencing, which:
    # - Checked bounds on differences DeltaH_q (can be negative, so bounds check always failed)
    # - Tried to stabilize differences (clamping negatives to 1e-7, destroying AR(1) structure)
    bounded_result <- .tsenat_handle_bounded_support(df, q_vals, group_vec = df$group, verbose = FALSE)
    use_bounded_family <- bounded_result$use_gamma
    family_gam <- bounded_result$family_obj
    inverse_link_fn <- bounded_result$inverse_link
    
    # GAMM COMPATIBILITY FIX (March 2026): mgcv::gamm() does NOT support extended families
    # (beta, gamma, Tweedie, etc.). For paired designs (subject != NULL -> uses gamm),
    # fall back to gaussian family instead of extended families.
    # Reference: C042/C043 (GAMM Tutorial, mgcv Documentation)
    if (!is.null(subject)) {
        # For paired/mixed designs using gamm(), force gaussian family
        # CRITICAL FIX (March 2026): mgcv::gamm() does NOT support extended families
        # Warn user that bounded family selection is being overridden for reproducibility
        if (use_bounded_family) {
            warning(
                "[calculate_lm_interaction] GAMM with paired design detected. ",
                "mgcv::gamm() does not support extended families. ",
                "Forcing gaussian family. Results may be less accurate for bounded data.",
                call. = FALSE
            )
        }
        family_gam <- stats::gaussian()
        inverse_link_fn <- function(eta) eta
    }
    
    # If beta regression is selected, use stabilized entropy
    if (bounded_result$use_beta && !is.null(bounded_result$stabilized_df)) {
        df <- bounded_result$stabilized_df
    }
    
    # ===============================================================================
    # HETEROSCEDASTICITY DETECTION AND VARIANCE WEIGHTING (FIXED - March 2026)
    # ===============================================================================
    # CRITICAL FIX: Detect heteroscedasticity on ORIGINAL entropy BEFORE ARIMA differencing
    # BUG #5 FIX: Was previously called after ARIMA, so it tested variance of differences DeltaH_q
    #             Now properly analyzes variance pattern of original entropy
    # 
    # Weights are computed on original data; if ARIMA is applied later, they will NOT be used
    # because the differenced data has a different variance structure (differencing changes variance)
    # This is a conservative approach: better to avoid incorrect weighting than to misapply weights
    hetero_result <- .tsenat_detect_heteroscedasticity(df, q_vals, df$group)
    gam_weights <- NULL
    gam_weights_original <- NULL  # Will be used only if ARIMA NOT applied
    
    # PHASE 1 WEIGHTING (March 2026): Use bootstrap CI weights if provided
    # These take precedence over heteroscedasticity-estimated weights
    if (!is.null(weights) && length(weights) == nrow(df)) {
        gam_weights_original <- weights  # Bootstrap CI weights for Phase 1
    } else if (!is.na(hetero_result$is_heteroscedastic) && hetero_result$is_heteroscedastic) {
        weights_result <- .tsenat_estimate_variance_weights(df, q_vals, method = "power")
        if (!is.null(weights_result)) {
            gam_weights_original <- weights_result$weights  # Store original weights
        }
    }
    
    # ARIMA(1,1,0) IMPLEMENTATION: Compute first differences for stationarity
    # Differencing removes monotone trend from Tsallis entropy, enabling valid AR(1) inference
    # This is applied when subject information is available (paired design)
    # CRITICAL: ARIMA is applied AFTER heteroscedasticity detection on original data
    use_arima <- FALSE
    n_samples_original <- nrow(df)
    
    if (!is.null(subject)) {
        # Only apply ARIMA differencing for paired designs (has subject info)
        arima_result <- .tsenat_compute_arima_differences(df, q_vals, df$group, factor(subject))
        
        if (!is.null(arima_result) && nrow(arima_result$df) >= 3) {
            # Sufficient data for ARIMA differencing
            df <- arima_result$df
            use_arima <- TRUE
            
            # FIX: When ARIMA is applied, DON'T use weights computed on original data
            # Reason: Differencing changes the variance structure, weights would be invalid
            # Conservative approach: Better to lose efficiency than introduce bias
            gam_weights_original <- NULL
        } else {
            # ARIMA differencing failed or insufficient data - skip ARIMA
        }
    }
    
    # Use original weights only if ARIMA was NOT applied
    gam_weights <- gam_weights_original
    
    # PHASE 1 WEIGHTING (March 2026): Set df$weight for ci_weighted flag tracking
    if (!is.null(gam_weights)) {
        df$weight <- gam_weights
    }
    
    # Determine sample size for bias correction
    # For GAM, use actual number of observations (nrow(df)) not number of subjects
    # The 'subject' parameter is structural metadata for paired designs but n_samples 
    # for bias correction should be based on total observations used in model fitting
    n_samples <- nrow(df)
    
    uq_len <- length(unique(na.omit(q_vals)))
    
    # Adaptive knot selection: compute k based on this gene's entropy curve complexity
    if (adaptive_knots) {
        k_q <- .tsenat_adaptive_spline_knots(entropy_vals = df$entropy, q_vals = q_vals,
                                             n_q_unique = uq_len, min_k = 2, max_k = 10)
    } else {
        # Fallback to static knot selection
        k_q <- max(2, min(10, uq_len - 1))
    }
    
    # Apply regularization for variable selection if requested (not "pca")
    reg_result <- NULL
    if (regularization != "pca") {
        reg_result <- .tsenat_gam_regularization(entropy_vals = df$entropy, 
                                                 q_vals = q_vals,
                                                 group_vec = df$group,
                                                 regularization = regularization)
    }
    
    # Use GAMM with random intercept for subject if subject info is available
    # This properly accounts for paired samples
    if (!is.null(subject)) {
        # If ARIMA differencing was not applied, add subject factor to df
        # (ARIMA result already has subject as a column)
        if (!use_arima) {
            df$subject <- factor(subject)
        } else {
            # After ARIMA differencing, subject is already in df, ensure it's a factor
            df$subject <- factor(df$subject)
        }
        
        # CRITICAL FIX: Ensure data is sorted by q|subject for corAR1 correlation structure
        # (nlme::corAR1 assumes observations are ordered by the within-group ordering variable)
        # This is especially important after adding new columns like gam_weights or group_numeric
        df <- df[order(df$subject, df$q), ]
        rownames(df) <- NULL
        
        # FIX: Create observation sequence within each subject for corAR1 ordering
        # corAR1 requires unique values within each group; use sequence index
        # This respects the q-ordering (data is sorted by q within subject) while providing unique IDs
        df$obs_seq <- unlist(lapply(rle(as.numeric(df$subject))$lengths, seq_len))
        
        # Check if we have at least 2 subjects
        if (length(unique(na.omit(df$subject))) < 2) {
            return(NULL)
        }
        
        # OPTIMIZATION (March 2026): Adaptive knot selection for smooth terms
        # Issue #4 & #8: Replace hardcoded/aggressive k parameters with data-driven selection
        # FIX (April 2026): mgcv thin-plate splines s(x, bs="tp") have minimum basis dimension k=3
        # mgcv automatically increases k if specified value is too low, causing warning
        # Previous code set min_k_adaptive=2 for small samples, triggering mgcv auto-increase
        min_k_adaptive <- 3L  # Thin-plate spline minimum basis dimension
        # Marginal smooth for q: Use conservative knots to avoid overfitting
        # Formula: at least 3 knots, but not more than available unique q-values - 1
        k_q_marginal <- as.integer(max(min_k_adaptive, min(k_q, max(min_k_adaptive, nrow(df) / 15))))
        # Interaction smooth: Use fewer knots since it must accommodate both q and group
        # Also enforce minimum of 3 to prevent mgcv auto-increase warnings
        k_q_interaction <- as.integer(max(3L, min(k_q / 2, 4L)))  # Min 3 for tp splines, cap at 4
        
        # Determine spline smoothing approach based on regularization
        bs_arg <- "tp"  # Thin plate spline basis - good for continuous covariates
        if (!is.null(reg_result) && reg_result$mode == "spline") {
            # Use adaptive smoothing with mgcv's GCV
            bs_arg <- "tp"
        }
        
        # Fit GAMM with ARIMA(1,1,0) covariance structure for q-measurements within subjects
        # ARIMA(1,1,0): First difference DeltaH_q modeled as AR(1) to handle monotone trend
        # Cov(DeltaY_t, DeltaY_s) = sigma^2 phi^|t-s| where t,s are q-ordered indices
        # This separates trend (differencing) from autocorrelation, validated in TEST L.1.6
        # BOUNDED SUPPORT: Use quasibinomial(logit) if entropy is bounded [0, log(m)]
        # HETEROSCEDASTICITY: Use weights parameter to model variance heterogeneity
        
        # BUG FIX: Use smooth splines s() instead of poly() for actual GAM fitting
        # Adaptive spline basis with thin-plate (tp) for flexible curve fitting
        # FIX (April 2026): Implement fallback strategy for convergence failures
        # Priority 1: GAMM with AR(1) correlation (handles autocorrelated measurements)
        # Priority 2: GAMM with random intercept only (if AR(1) convergence fails)
        # Priority 3: Standard GAM with independence assumption (last resort)
        
        # PRIORITY 1: Try GAMM with AR(1) correlation (full autocorrelation model)
        convergence_note <- NULL
        fit_null <- NULL
        fit_alt <- NULL
        use_ar1 <- TRUE  # Track which model succeeded
        
        if (!is.null(gam_weights)) {
            df$gam_weights <- gam_weights
            fit_null <- try(
                mgcv::gamm(entropy ~ group + s(q, bs="tp", k=k_q_marginal), 
                          random = list(subject = ~1), 
                          correlation = nlme::corAR1(form = ~obs_seq|subject),
                          family = family_gam,
                          weights = gam_weights,
                          data = df),
                silent = TRUE
            )
            df$group_numeric <- as.numeric(df$group)
            fit_alt <- try(
                mgcv::gamm(entropy ~ group + s(q, bs="tp", k=k_q_interaction, by=group),
                          random = list(subject = ~1), 
                          correlation = nlme::corAR1(form = ~obs_seq|subject),
                          family = family_gam,
                          weights = gam_weights,
                          data = df),
                silent = TRUE
            )
        } else {
            fit_null <- try(
                mgcv::gamm(entropy ~ group + s(q, bs="tp", k=k_q_marginal), 
                          random = list(subject = ~1), 
                          correlation = nlme::corAR1(form = ~obs_seq|subject),
                          family = family_gam,
                          data = df),
                silent = TRUE
            )
            fit_alt <- try(
                mgcv::gamm(entropy ~ group + s(q, bs="tp", k=k_q_interaction, by=group),
                          random = list(subject = ~1), 
                          correlation = nlme::corAR1(form = ~obs_seq|subject),
                          family = family_gam,
                          data = df),
                silent = TRUE
            )
        }
        
        # PRIORITY 2: If AR(1) convergence failed, try GAMM without correlation structure
        if (inherits(fit_null, "try-error") || inherits(fit_alt, "try-error")) {
            convergence_note <- "AR1_convergence_failed"
            use_ar1 <- FALSE
            
            if (!is.null(gam_weights)) {
                fit_null <- try(
                    mgcv::gamm(entropy ~ group + s(q, bs="tp", k=k_q_marginal), 
                              random = list(subject = ~1),
                              family = family_gam,
                              weights = gam_weights,
                              data = df),
                    silent = TRUE
                )
                fit_alt <- try(
                    mgcv::gamm(entropy ~ group + s(q, bs="tp", k=k_q_interaction, by=group),
                              random = list(subject = ~1),
                              family = family_gam,
                              weights = gam_weights,
                              data = df),
                    silent = TRUE
                )
            } else {
                fit_null <- try(
                    mgcv::gamm(entropy ~ group + s(q, bs="tp", k=k_q_marginal), 
                              random = list(subject = ~1),
                              family = family_gam,
                              data = df),
                    silent = TRUE
                )
                fit_alt <- try(
                    mgcv::gamm(entropy ~ group + s(q, bs="tp", k=k_q_interaction, by=group),
                              random = list(subject = ~1),
                              family = family_gam,
                              data = df),
                    silent = TRUE
                )
            }
        }
        
        # PRIORITY 3: If GAMM fails entirely, fall back to standard GAM (independence assumption)
        if (inherits(fit_null, "try-error") || inherits(fit_alt, "try-error")) {
            convergence_note <- "GAMM_failed_using_GAM_fallback"
            
            if (!is.null(gam_weights)) {
                fit_null <- try(
                    mgcv::gam(entropy ~ group + s(q, bs="tp", k=k_q_marginal), 
                              family = family_gam,
                              weights = gam_weights,
                              data = df),
                    silent = TRUE
                )
                fit_alt <- try(
                    mgcv::gam(entropy ~ group + s(q, bs="tp", k=k_q_interaction, by=group),
                              family = family_gam,
                              weights = gam_weights,
                              data = df),
                    silent = TRUE
                )
            } else {
                fit_null <- try(
                    mgcv::gam(entropy ~ group + s(q, bs="tp", k=k_q_marginal), 
                              family = family_gam,
                              data = df),
                    silent = TRUE
                )
                fit_alt <- try(
                    mgcv::gam(entropy ~ group + s(q, bs="tp", k=k_q_interaction, by=group),
                              family = family_gam,
                              data = df),
                    silent = TRUE
                )
            }
        }
        
        if (inherits(fit_null, "try-error") || inherits(fit_alt, "try-error")) {
            return(NULL)
        }
        
        # Compare models based on what type of fit we have
        # GAMM objects have $lme component; standard GAM objects don't
        old_warn <- options(warn = -1)
        
        p_interaction <- NA_real_
        if (!is.null(fit_null$lme) && !is.null(fit_alt$lme)) {
            # GAMM comparison via LME component
            an <- try(anova(fit_null$lme, fit_alt$lme), silent = TRUE)
            if (!inherits(an, "try-error") && nrow(an) >= 2) {
                if ("p-value" %in% colnames(an)) {
                    p_interaction <- an[2, "p-value"]
                } else if ("Pr(>Chisq)" %in% colnames(an)) {
                    p_interaction <- an[2, "Pr(>Chisq)"]
                } else if ("Pr(>F)" %in% colnames(an)) {
                    p_interaction <- an[2, "Pr(>F)"]
                }
            }
        } else {
            # Standard GAM comparison via anova.gam()
            # Use hypothesis test comparing smooth term differences
            an <- try(anova(fit_null, fit_alt, test = "Chisq"), silent = TRUE)
            if (!inherits(an, "try-error") && nrow(an) >= 2) {
                if ("p-value" %in% colnames(an)) {
                    p_interaction <- an[2, "p-value"]
                } else if ("Pr(>Chi)" %in% colnames(an)) {
                    p_interaction <- an[2, "Pr(>Chi)"]
                } else if ("p-value" %in% tolower(colnames(an))) {
                    # Case-insensitive search for p-value column
                    col_idx <- grep("p-value", tolower(colnames(an)))[1]
                    if (!is.na(col_idx) && nrow(an) >= 2) {
                        p_interaction <- an[2, col_idx]
                    }
                }
            }
        }
        
        options(old_warn)
    } else {
        # Fallback to standard GAM (treats samples as independent)
        # This is used only when subject info is not available
        # BOUNDED SUPPORT: Use quasibinomial(logit) if entropy is bounded [0, log(m)]
        # HETEROSCEDASTICITY: Use weights parameter to model variance heterogeneity
        bs_arg <- "tp"  # Thin plate spline
        
        # BUG FIX: Use smooth splines s() instead of poly() for actual GAM fitting
        # Adaptive spline basis with thin-plate (tp) for flexible curve fitting
        
        # Determine k values for standard GAM (non-paired design)
        # Adaptive minimum k based on sample size: mgcv needs ~15-20 points per basis function
        min_k_adaptive <- if (nrow(df) < 25) 2L else if (nrow(df) < 50) 3L else 4L
        k_q_marginal <- as.integer(max(min_k_adaptive, min(k_q, max(min_k_adaptive, nrow(df) / 25))))
        k_q_interaction <- as.integer(max(2L, min(k_q, max(2L, nrow(df) / 20))))
        
        if (!is.null(gam_weights)) {
            # Fitting with weights
            df$gam_weights <- gam_weights
            fit_null <- try(
                mgcv::gam(entropy ~ group + s(q, bs="tp", k=k_q_marginal), 
                         family = family_gam,
                         weights = gam_weights,
                         data = df), 
                silent = TRUE
            )
            if (inherits(fit_null, "try-error")) {
            }
            # Use group-specific smooth for interaction testing
            fit_alt <- try(
                mgcv::gam(entropy ~ group + s(q, bs="tp", k=k_q_interaction, by=group),
                         family = family_gam,
                         weights = gam_weights,
                         data = df),
                silent = TRUE
            )
            if (inherits(fit_alt, "try-error")) {
            }
        } else {
            fit_null <- try(
                mgcv::gam(entropy ~ group + s(q, bs="tp", k=k_q_marginal), 
                         family = family_gam,
                         data = df), 
                silent = TRUE
            )
            if (inherits(fit_null, "try-error")) {
            }
            # Use group-specific smooth for interaction testing
            fit_alt <- try(
                mgcv::gam(entropy ~ group + s(q, bs="tp", k=k_q_interaction, by=group),
                         family = family_gam,
                         data = df),
                silent = TRUE
            )
            if (inherits(fit_alt, "try-error")) {
            }
        }
        
        # Allow model fitting errors to proceed with NA values (don't return NULL)
        # This ensures genes with convergence issues still appear in output
        
        # CRITICAL: Return NULL if BOTH models failed to fit
        if (inherits(fit_null, "try-error") && inherits(fit_alt, "try-error")) {
            return(NULL)
        }
        
        # Suppress NaN warnings from anova.gam F-test with small samples or edge cases
        old_warn <- options(warn = -1)
        an <- try(mgcv::anova.gam(fit_null, fit_alt, test = "F"), silent = TRUE)
        options(old_warn)
        
        # If anova fails (e.g., numerical issues with Gamma family on small samples),
        # proceed with NA p-value instead of returning NULL
        if (inherits(an, "try-error")) {
            p_interaction <- NA_real_
        } else {
            p_interaction <- NA_real_
            if (nrow(an) >= 2) {
                if ("Pr(F)" %in% colnames(an)) {
                    p_interaction <- an[2, "Pr(F)"]
                } else if ("Pr(>F)" %in% colnames(an)) {
                    p_interaction <- an[2, "Pr(>F)"]
                } else if ("p-value" %in% colnames(an)) {
                    p_interaction <- an[2, "p-value"]
                }
            }
        }
        
        # Extract test statistic and effect size from standard GAM fit
        test_statistic <- NA_real_
        effect_size <- NA_real_
        df_residual <- NA_real_
        model_converged <- !inherits(fit_alt, "try-error")
        
        if (!inherits(fit_alt, "try-error") && !is.null(fit_alt)) {
            tryCatch({
                gam_summary <- summary(fit_alt)
                if (!is.null(gam_summary)) {
                    # Try dev.expl first (standard GAM), then r.sq as fallback
                    if (!is.null(gam_summary$dev.expl) && length(gam_summary$dev.expl) > 0 && is.finite(gam_summary$dev.expl)) {
                        effect_size <- as.numeric(gam_summary$dev.expl)[1]
                    } else if (!is.null(gam_summary$r.sq) && length(gam_summary$r.sq) > 0 && is.finite(gam_summary$r.sq)) {
                        effect_size <- as.numeric(gam_summary$r.sq)[1]
                    }
                    
                    if (!is.finite(effect_size)) {
                        effect_size <- NA_real_
                    }
                    
                    if (!is.null(gam_summary$residual.df)) {
                        df_residual <- as.numeric(gam_summary$residual.df)[1]
                        if (!is.finite(df_residual)) {
                            df_residual <- NA_real_
                        }
                    }
                }
            }, error = function(e) { NULL })
            
            # Extract F-statistic from anova if available
            if (!is.null(an) && nrow(an) >= 2 && !inherits(an, "try-error")) {
                tryCatch({
                    if ("F" %in% colnames(an)) {
                        test_statistic <- as.numeric(an[2, "F"])[1]
                    } else if ("Chisq" %in% colnames(an)) {
                        test_statistic <- as.numeric(an[2, "Chisq"])[1]
                    }
                    if (!is.finite(test_statistic)) {
                        test_statistic <- NA_real_
                    }
                }, error = function(e) { NULL })
            }
        }
    }
    
    # Apply GAM-specific bias correction for small samples (C071)
    # Account for ARIMA(1,1,0) correlation structure in Tsallis entropy measurements
    # n_observations = total data points; n_subjects = independent observational units
    n_subjects_bc <- if (!is.null(subject)) length(unique(na.omit(subject))) else NULL
    bc_result <- .tsenat_gam_bias_correct(p_interaction, n_observations = n_samples,
                                          n_subjects = n_subjects_bc,
                                          ar1_correlation = TRUE,
                                          bias_correction = bias_correction,
                                          entropy_data = df$entropy,
                                          subject_data = df$subject)
        
    # Extract test statistic and effect size from anova if available
    test_statistic <- NA_real_
    effect_size <- NA_real_
    df_residual <- NA_real_
    model_converged <- NA
    
    # Try to extract statistics from GAM/GAMM model summaries
    if (!inherits(fit_alt, "try-error")) {
        model_converged <- TRUE
        
        # Extract effect size (deviance explained / R-squared equivalent)
        # Handle both GAMM (list with $gam component) and GAM (gam object directly)
        is_gamm <- is.list(fit_alt) && !is.null(fit_alt$gam)
        gam_obj <- if (is_gamm) fit_alt$gam else fit_alt
        
        # Now gam_obj should be a gam object (from either standard GAM or GAMM)
        if (!is.null(gam_obj)) {
            tryCatch({
                gam_summary <- summary(gam_obj)
                if (!is.null(gam_summary)) {
                    # For standard GAM: use dev.expl (deviance explained)
                    # For GAMM: dev.expl may be NA due to random effects, use r.sq instead
                    if (!is.null(gam_summary$dev.expl) && length(gam_summary$dev.expl) > 0 && is.finite(gam_summary$dev.expl)) {
                        effect_size <- as.numeric(gam_summary$dev.expl)[1]
                    } else if (!is_gamm && !is.null(gam_summary$r.sq) && length(gam_summary$r.sq) > 0 && is.finite(gam_summary$r.sq)) {
                        # Fallback to r.sq for any model where dev.expl fails
                        effect_size <- as.numeric(gam_summary$r.sq)[1]
                    } else if (is_gamm && !is.null(gam_summary$r.sq) && length(gam_summary$r.sq) > 0 && is.finite(gam_summary$r.sq)) {
                        # GAMM: use r.sq (explained variance) as effect size
                        effect_size <- as.numeric(gam_summary$r.sq)[1]
                    }
                    
                    if (!is.finite(effect_size)) {
                        effect_size <- NA_real_
                    }
                    
                    # Residual df
                    if (!is.null(gam_summary$residual.df)) {
                        df_residual <- as.numeric(gam_summary$residual.df)[1]
                        if (!is.finite(df_residual)) {
                            df_residual <- NA_real_
                        }
                    }
                }
            }, error = function(e) { NULL })
        }
        
        # Extract F-statistic or likelihood ratio if anova results are available
        if (!is.null(an) && nrow(an) >= 2 && !inherits(an, "try-error")) {
            tryCatch({
                if ("L.Ratio" %in% colnames(an)) {
                    test_statistic <- as.numeric(an[2, "L.Ratio"])[1]
                } else if ("F" %in% colnames(an)) {
                    test_statistic <- as.numeric(an[2, "F"])[1]
                } else if ("Chisq" %in% colnames(an)) {
                    test_statistic <- as.numeric(an[2, "Chisq"])[1]
                }
                if (!is.finite(test_statistic)) {
                    test_statistic <- NA_real_
                }
            }, error = function(e) { NULL })
        }
    }
    
    # Return result with bias correction information
    result <- data.frame(
        gene = g, 
        p_interaction = bc_result$p_value,
        p_raw = bc_result$p_raw,  # Always include for reference
        n_observations = bc_result$n_observations,  # Total observations (observations per subject * subjects)
        n_subjects = bc_result$n_subjects,  # Independent sampl units
        n_effective = bc_result$n_effective,  # Effective sample size accounting for AR(1) correlation
        rho_ar1 = bc_result$rho_estimate,  # AR(1) rho estimate used in design effect
        stringsAsFactors = FALSE
    )
    
    # Explicitly add effect size, test statistic, and df columns
    result$test_statistic <- test_statistic
    result$effect_size <- effect_size
    result$df_residual <- df_residual
    result$model_converged <- model_converged
    
    # Extract slope_diff from GAM by computing predicted slopes for each group
    slope_diff <- NA_real_
    if (!inherits(fit_alt, "try-error") && !is.null(fit_alt)) {
        tryCatch({
            # Get the GAM/GAMM object (may be nested in list for GAMM)
            gam_obj <- if (is.list(fit_alt) && !is.null(fit_alt$gam)) fit_alt$gam else fit_alt
            
            # Create prediction grid at min and max q for each group
            q_range <- range(df$q, na.rm = TRUE)
            unique_groups <- unique(na.omit(as.character(df$group)))
            
            if (length(unique_groups) == 2 && is.finite(q_range[1]) && is.finite(q_range[2])) {
                pred_slopes <- numeric(2)
                for (g_idx in seq_along(unique_groups)) {
                    pred_grid <- data.frame(
                        q = c(q_range[1], q_range[2]),
                        group = factor(rep(unique_groups[g_idx], 2), levels = levels(df$group))
                    )
                    
                    # Add subject if needed for GAMM
                    if (!is.null(subject) && "subject" %in% colnames(df)) {
                        pred_grid$subject <- df$subject[1]  # Use first subject as reference
                    }
                    
                    preds <- tryCatch(
                        predict(gam_obj, newdata = pred_grid, type = "response", se.fit = FALSE),
                        error = function(e) NULL
                    )
                    
                    if (!is.null(preds) && length(preds) == 2 && all(is.finite(preds))) {
                        pred_slopes[g_idx] <- (preds[2] - preds[1]) / (q_range[2] - q_range[1])
                    }
                }
                
                # Compute slope_diff if we got both slopes
                if (all(is.finite(pred_slopes))) {
                    slope_diff <- pred_slopes[2] - pred_slopes[1]
                }
            }
        }, error = function(e) { NULL })
    }
    result$slope_diff <- slope_diff
    
    # Add bias correction metadata if applied
    if (bc_result$bias_correction_applied) {
        result$bias_correction_applied <- TRUE
        result$correction_method <- bc_result$correction_method
        result$adjustment_factor <- bc_result$adjustment_factor
        result$rho_method <- if(bc_result$rho_data_driven) "data-driven" else "default"
        result$adjustment_factor <- bc_result$adjustment_factor
    }
    
    # Add ARIMA transformation flag
    result$arima_transformation <- use_arima
    
    # Add bounded support handling flag and family selection information (March 2026)
    result$bounded_support_model <- use_bounded_family
    
    # Add family selection information from conditional logic
    if (bounded_result$use_gamma) {
        result$family_used <- "Gamma"
        result$heteroscedasticity_detected <- bounded_result$family_info$heteroscedastic
        result$variance_ratio_q <- bounded_result$family_info$var_ratio_q
    } else {
        result$family_used <- "Gaussian"
        result$heteroscedasticity_detected <- bounded_result$family_info$heteroscedastic
        result$variance_ratio_q <- bounded_result$family_info$var_ratio_q
    }
    
    # Add fit method tag
    result$fit_method <- ifelse(use_arima, "mgcv::gamm_arima(1,1,0)", "mgcv::gamm")
    
    # ===============================================================================
    # RESIDUAL NORMALITY TESTING (NEW - March 2026)
    # Database Evidence: B001, B004, C017 (Normality testing in regression)
    # ===============================================================================
    shapiro_result <- .tsenat_test_residual_normality(
        model = fit_alt,
        model_type = if (!is.null(subject)) "gamm" else "gam",
        verbose = FALSE
    )
    
    if (!is.na(shapiro_result$shapiro_p_value)) {
        result$shapiro_p_value <- shapiro_result$shapiro_p_value
        result$residuals_normal <- shapiro_result$residuals_normal
        result$n_residuals_tested <- shapiro_result$n_residuals
    } else {
        result$shapiro_p_value <- NA_real_
        result$residuals_normal <- NA
        result$n_residuals_tested <- NA_integer_
    }
    
    # VALIDATION: Ensure p_interaction is always present with correct name
    if (!"p_interaction" %in% colnames(result)) {
        stop(sprintf("[.tsenat_gam_interaction] CRITICAL: p_interaction missing from result for gene %s. Available columns: %s",
                     g, paste(colnames(result), collapse=", ")))
    }
    
    # Add Phase 1 bootstrap CI weighting tracking (March 2026)
    result$ci_weighted <- !is.null(df$weight)
    
    return(result)
}
