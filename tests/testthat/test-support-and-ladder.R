# assess_support(), estimand_feasibility(), cleanTMLE:::who_is_unsupported(),
# declare_estimand_ladder(), run_estimand_ladder(), and the estimators
# cleanTMLE:::run_att_tmle(), cleanTMLE:::estimate_ato(), cleanTMLE:::run_trimmed_tmle(), cleanTMLE:::implausibility_check().

# A cohort with GOOD overlap and a known constant additive effect.
.good_cohort <- function(n = 400, seed = 31, effect = 0.15) {
  set.seed(seed)
  x1 <- stats::rnorm(n); x2 <- stats::rbinom(n, 1, 0.5)
  g  <- stats::plogis(0.4 * x1 + 0.3 * x2)
  A  <- stats::rbinom(n, 1, g)
  p0 <- stats::plogis(-1.2 + 0.4 * x1 + 0.3 * x2)
  y  <- stats::rbinom(n, 1, pmin(pmax(p0 + effect * A, 0.001), 0.999))
  data.frame(x1 = x1, x2 = x2, treatment = A, outcome = y)
}

# A cohort with a SEVERE practical positivity violation: a strong covariate
# nearly separates the arms, so 1/g weights explode.
.severe_cohort <- function(n = 600, seed = 41, effect = 0.10) {
  set.seed(seed)
  x1 <- stats::rnorm(n)
  x2 <- stats::rbinom(n, 1, 0.3)
  g  <- stats::plogis(-2.5 + 3.5 * x1 + 0.5 * x2)
  A  <- stats::rbinom(n, 1, g)
  p0 <- stats::plogis(-1.5 + 0.5 * x1)
  y  <- stats::rbinom(n, 1, pmin(pmax(p0 + effect * A, 0.001), 0.999))
  data.frame(x1 = x1, x2 = x2, treatment = A, outcome = y)
}

test_that("assess_support grades a known score vector exactly", {
  n <- 500
  set.seed(5)
  dat <- data.frame(x1 = rnorm(n), treatment = rbinom(n, 1, 0.5),
                    outcome = rbinom(n, 1, 0.2))
  # 15 scores below 0.05 (3% outside) and a max of 0.90: FLAG by the
  # percent rule, not the weight rule.
  scores <- c(rep(0.02, 15), runif(n - 16, 0.10, 0.85), 0.90)
  lock <- create_analysis_lock(dat, "treatment", "outcome", "x1")
  psf  <- fit_ps(lock, "external", scores = scores, truncate = NULL)
  sup  <- assess_support(psf, tree_search = FALSE)
  expect_equal(sup$summary$pct_outside_band, 3)
  expect_equal(sup$summary$n_below_band, 15)
  expect_identical(sup$verdict, "FLAG")
  expect_match(sup$caveat, "Mild")
  expect_output(print(sup), "Verdict: FLAG")
})

test_that("assess_support: PASS on good overlap, SEVERE on separation", {
  good <- .good_cohort()
  lock_g <- create_analysis_lock(good, "treatment", "outcome", c("x1", "x2"))
  sup_g <- assess_support(fit_ps(lock_g, "glm"), tree_search = FALSE)
  expect_identical(sup_g$verdict, "PASS")
  expect_gt(sup_g$c_statistic, 0.5)
  expect_lt(sup_g$c_statistic, 0.75)

  sev <- .severe_cohort()
  lock_s <- create_analysis_lock(sev, "treatment", "outcome", c("x1", "x2"))
  sup_s <- assess_support(fit_ps(lock_s, "glm"), tree_search = FALSE)
  expect_true(sup_s$verdict %in% c("SEVERE", "FAIL"))
  expect_gt(sup_s$summary$max_iptw_weight,
            sup_g$summary$max_iptw_weight)
  p <- plot(sup_s)
  expect_s3_class(p, "ggplot")
})

test_that("near-deterministic strata escalate PASS to SEVERE, never to FAIL", {
  n <- 500
  set.seed(7)
  x1 <- stats::rnorm(n)
  # A binary covariate whose 'on' stratum (60 patients) holds 3 treated:
  # globally invisible, conditionally near-deterministic.
  x_flag <- c(rep(1, 60), rep(0, n - 60))
  A <- stats::rbinom(n, 1, stats::plogis(0.3 * x1))
  A[1:60] <- c(rep(1, 3), rep(0, 57))
  y <- stats::rbinom(n, 1, 0.2)
  dat <- data.frame(x1 = x1, x_flag = x_flag, treatment = A, outcome = y)
  lock <- create_analysis_lock(dat, "treatment", "outcome",
                             c("x1", "x_flag"))
  sup <- assess_support(fit_ps(lock, "glm"), tree_search = FALSE)
  expect_true(isTRUE(sup$escalated))
  expect_identical(sup$verdict, "SEVERE")
  expect_true(sup$verdict_global %in% c("PASS", "FLAG"))
  expect_true("x_flag" %in% sup$near_deterministic$covariate)
})

test_that("the tree search names a multivariate violation region", {
  skip_if_not_installed("rpart")
  n <- 600
  set.seed(9)
  x1 <- stats::rbinom(n, 1, 0.5)
  x2 <- stats::rbinom(n, 1, 0.5)
  g <- ifelse(x1 == 1 & x2 == 1, 0.02, 0.5)   # the joint stratum is starved
  A <- stats::rbinom(n, 1, g)
  dat <- data.frame(x1 = x1, x2 = x2, treatment = A,
                    outcome = stats::rbinom(n, 1, 0.2))
  lock <- create_analysis_lock(dat, "treatment", "outcome", c("x1", "x2"))
  sup <- assess_support(fit_ps(lock, "glm"), tree_min_n = 40L)
  expect_false(is.null(sup$violation_regions))
  expect_true(any(sup$violation_regions$p_treated < 0.10))
})

test_that("estimand_feasibility separates the ATE from ATT and ATO on a severe cohort", {
  sev <- .severe_cohort()
  lock <- create_analysis_lock(sev, "treatment", "outcome", c("x1", "x2"))
  psf <- fit_ps(lock, "glm")
  fea <- estimand_feasibility(psf)
  tab <- fea$table
  expect_true(all(c("ATE", "ATT", "ATO") %in% tab$estimand))
  v_ate <- tab$verdict[tab$estimand == "ATE"]
  expect_true(v_ate %in% c("SEVERE", "FAIL"))
  expect_false(tab$feasible[tab$estimand == "ATE"])
  # ATO weights are bounded by one, so its verdict is PASS by construction.
  expect_identical(tab$verdict[tab$estimand == "ATO"], "PASS")
  expect_lte(tab$max_weight[tab$estimand == "ATO"], 1)
  # The ATT's control weights are bounded by max g/(1-g), far below the
  # ATE's max weight on the same fit.
  expect_lt(tab$max_weight[tab$estimand == "ATT"],
            tab$max_weight[tab$estimand == "ATE"])
  # Trimmed rows record who they remove.
  tr <- tab[startsWith(tab$estimand, "trimmed_ATE"), ]
  expect_true(all(tr$n_removed_treated + tr$n_removed_control > 0))
  expect_output(print(fea), "Target populations")
})

test_that("who_is_unsupported profiles the removed patients", {
  n <- 500
  set.seed(3)
  x1 <- stats::rnorm(n)
  dat <- data.frame(x1 = x1, treatment = stats::rbinom(n, 1, 0.5),
                    outcome = stats::rbinom(n, 1, 0.2))
  # Push exactly the high-x1 patients outside the band.
  scores <- ifelse(x1 > stats::quantile(x1, 0.9), 0.97, 0.5)
  lock <- create_analysis_lock(dat, "treatment", "outcome", "x1")
  psf  <- fit_ps(lock, "external", scores = scores, truncate = NULL)
  prof <- cleanTMLE:::who_is_unsupported(psf, vars = "x1")
  expect_true(all(c("all") %in% prof$population))
  sm <- prof$smd[prof$population == "all" & prof$variable == "x1"]
  expect_gt(sm, 1)   # the removed are far higher on x1 by construction
})

test_that("declare_estimand_ladder records the rule on the lock", {
  dat <- .good_cohort(200)
  lock <- create_analysis_lock(dat, "treatment", "outcome", c("x1", "x2"))
  lock <- declare_estimand_ladder(lock, primary = "ATE",
                                  fallbacks = c("trimmed_ATE", "ATT", "ATO"),
                                  trigger = "SEVERE",
                                  bias_to_null_floor = 0.0035)
  expect_identical(lock$estimand_ladder$primary, "ATE")
  expect_identical(lock$estimand_ladder$trigger, "SEVERE")
  expect_true(any(lock$design_log$type == "estimand_ladder"))
  expect_error(declare_estimand_ladder(lock, primary = "banana"))
})

test_that("implausibility_check flags the documented failure modes", {
  set.seed(2)
  y <- stats::rbinom(400, 1, 0.10)
  a <- rep(c(0, 1), 200)
  y[a == 1] <- stats::rbinom(200, 1, 0.05)   # crude difference negative
  crude <- mean(y[a == 1]) - mean(y[a == 0])
  expect_lt(crude, 0)
  # Sign flip.
  g1 <- cleanTMLE:::implausibility_check(0.20, y, a, family = "binomial")
  expect_true(g1$implausible)
  expect_match(g1$implausible_reason, "sign differs")
  # Exceeds the largest arm rate.
  g2 <- cleanTMLE:::implausibility_check(-0.90, y, a, family = "binomial")
  expect_true(g2$implausible)
  expect_match(g2$implausible_reason, "largest arm rate")
  # A sane estimate passes.
  g3 <- cleanTMLE:::implausibility_check(crude, y, a, family = "binomial")
  expect_false(g3$implausible)
  # Continuous: exceeds observed range.
  yc <- stats::rnorm(400)
  g4 <- cleanTMLE:::implausibility_check(50, yc, a, family = "gaussian")
  expect_true(g4$implausible)
})

test_that("run_att_tmle recovers a known effect and reuses the ATE's g spec", {
  skip_if_not_installed("tmle")
  dat <- .good_cohort(500, seed = 51, effect = 0.15)
  lock <- create_analysis_lock(dat, "treatment", "outcome", c("x1", "x2"),
                             sl_library = "SL.glm", seed = 17L)
  att <- cleanTMLE:::run_att_tmle(lock, sl_library = "SL.glm")
  # Constant additive effect: ATT truth = ATE truth = 0.15. A single seed
  # carries binomial noise of about 0.045 SE, so recovery is judged against
  # the fit's own uncertainty.
  expect_lt(abs(att$estimate - 0.15), max(0.10, 2.6 * att$se))
  expect_match(att$estimand, "outcome observed")
  # The same fit's ATE agrees with the ATT under good overlap.
  expect_lt(abs(att$estimate - att$ate_same_fit$estimate), 0.05)
  # The g specification is the shared builder's ATE-path specification.
  ate_args <- cleanTMLE:::.tmle_delegate_args(lock, family = "binomial",
                                              use_delta = FALSE,
                                              sl_library = "SL.glm")
  g_fields <- intersect(c("g.SL.library", "V.g", "prescreenW.g", "gbound"),
                        names(att$spec))
  for (f in g_fields)
    expect_identical(att$spec[[f]], ate_args[[f]])
  expect_output(print(att), "Complete-case ATT")
})

test_that("estimate_ato is unbiased under randomization and balances exactly", {
  dat <- .good_cohort(800, seed = 61, effect = 0.15)
  lock <- create_analysis_lock(dat, "treatment", "outcome", c("x1", "x2"),
                             sl_library = "SL.glm", seed = 19L)
  # Under a constant g = 0.5 the ATO equals the ATE and the augmented
  # estimator reduces to a covariate-adjusted difference in means.
  psf <- fit_ps(lock, "external", truncate = NULL, scores = rep(0.5, nrow(dat)))
  ato <- cleanTMLE:::estimate_ato(lock, psf, sl_library = "SL.glm")
  expect_lt(abs(ato$estimate - 0.15), 0.07)
  expect_lt(ato$ci_lower, 0.15); expect_gt(ato$ci_upper, 0.15)
  # With a fitted logistic g the overlap weights balance the covariates.
  psf2 <- fit_ps(lock, "glm")
  ato2 <- cleanTMLE:::estimate_ato(lock, psf2, sl_library = "SL.glm")
  expect_lt(ato2$max_abs_weighted_smd, 0.05)
  expect_lt(abs(ato2$estimate - ato2$hajek_estimate), 0.05)
  expect_output(print(ato2), "overlap")
})

test_that("run_trimmed_tmle trims, refits, and records the dropped", {
  skip_if_not_installed("tmle")
  sev <- .severe_cohort(700, seed = 71, effect = 0.10)
  lock <- create_analysis_lock(sev, "treatment", "outcome", c("x1", "x2"),
                             sl_library = "SL.glm", seed = 23L)
  psf <- fit_ps(lock, "glm")
  trm <- cleanTMLE:::run_trimmed_tmle(lock, psf, levels = c(0.05, 0.10),
                          sl_library = "SL.glm", verbose = FALSE)
  expect_s3_class(trm, "trimmed_tmle_fit")
  expect_gt(trm$n_dropped_treated + trm$n_dropped_control, 0)
  expect_lt(trm$n, nrow(sev))
  expect_false(trm$verdict %in% c("SEVERE", "FAIL"))
  expect_match(trm$population, "trimmed g in")
  # The Crump rule returns a sensible level: higher under separation than
  # under good overlap.
  g_sev <- psf$ps_raw
  g_good <- fit_ps(create_analysis_lock(.good_cohort(400), "treatment",
                                      "outcome", c("x1", "x2")), "glm")$ps_raw
  expect_gte(cleanTMLE:::.crump_alpha(g_sev),
             cleanTMLE:::.crump_alpha(g_good))
})

test_that("run_estimand_ladder switches off an infeasible primary and logs it", {
  skip_if_not_installed("tmle")
  sev <- .severe_cohort(700, seed = 81, effect = 0.10)
  lock <- create_analysis_lock(sev, "treatment", "outcome", c("x1", "x2"),
                             sl_library = "SL.glm", seed = 29L)
  lock <- declare_estimand_ladder(lock, primary = "ATE",
                                  fallbacks = c("trimmed_ATE", "ATT", "ATO"))
  psf <- fit_ps(lock, "glm")
  res <- run_estimand_ladder(lock, psf, sl_library = "SL.glm",
                             verbose = FALSE)
  expect_s3_class(res, "estimand_ladder_result")
  expect_false(res$primary_feasible)
  expect_false("ATE" %in% res$table$estimand)
  expect_true(all(c("ATT", "ATO") %in% res$table$estimand))
  expect_true(any(res$design_log$type == "estimand_switch"))
  expect_true(all(res$table$support_verdict[res$table$estimand == "ATO"] ==
                    "PASS"))
  expect_output(print(res), "Estimand ladder")
  p <- plot(res)
  expect_s3_class(p, "ggplot")

  # With an override the infeasible primary is estimated and labelled.
  res2 <- run_estimand_ladder(lock, psf, sl_library = "SL.glm",
                              override_reason = "diagnostic comparison",
                              verbose = FALSE)
  expect_true("ATE" %in% res2$table$estimand)
  expect_true(any(res2$design_log$type == "override"))
})

test_that("run_estimand_ladder keeps a feasible primary first", {
  skip_if_not_installed("tmle")
  good <- .good_cohort(400, seed = 91)
  lock <- create_analysis_lock(good, "treatment", "outcome", c("x1", "x2"),
                             sl_library = "SL.glm", seed = 37L)
  lock <- declare_estimand_ladder(lock, primary = "ATE",
                                  fallbacks = c("ATT", "ATO"))
  psf <- fit_ps(lock, "glm")
  res <- run_estimand_ladder(lock, psf, sl_library = "SL.glm",
                             verbose = FALSE)
  expect_true(res$primary_feasible)
  expect_true("ATE" %in% res$table$estimand)
  expect_true(res$table$is_primary[res$table$estimand == "ATE"])
  # Under good overlap and a constant effect the rungs agree.
  ests <- res$table$estimate
  expect_lt(diff(range(ests)), 0.1)
})

test_that("estimate_effect is the one front door and matches its workers", {
  skip_if_not_installed("tmle")
  dat <- .good_cohort(400, seed = 101, effect = 0.15)
  lock <- create_analysis_lock(dat, "treatment", "outcome", c("x1", "x2"),
                               sl_library = "SL.glm", seed = 43L)
  psf <- fit_ps(lock, "glm")
  att1 <- estimate_effect(lock, psf, estimand = "ATT",
                          sl_library = "SL.glm")
  att2 <- cleanTMLE:::run_att_tmle(lock, sl_library = "SL.glm")
  expect_equal(att1$estimate, att2$estimate, tolerance = 1e-10)
  ato <- estimate_effect(lock, psf, estimand = "ATO", sl_library = "SL.glm")
  expect_match(ato$estimand, "overlap")
  crude <- estimate_effect(lock, estimand = "ATE", estimator = "crude")
  expect_true(is.finite(crude$estimate))
  # The ATT ignores an ipcw request with a warning: complete case by design.
  dat2 <- dat; dat2$outcome[1:60] <- NA
  lock2 <- create_analysis_lock(dat2, "treatment", "outcome",
                                c("x1", "x2"), sl_library = "SL.glm",
                                seed = 43L)
  expect_warning(estimate_effect(lock2, psf, estimand = "ATT",
                                 missing = "ipcw", sl_library = "SL.glm"),
                 "complete case")
  # Missing ps_fit is refused where it is needed.
  expect_error(estimate_effect(lock, estimand = "ATO"), "needs a ps_fit")
  # A steps-carrying fit supports both diagnostic plots (the clever plot
  # resolves the targeting step from $steps; regression for the vignette).
  fit_steps <- estimate_effect(lock, psf, estimand = "ATE",
                               estimator = "tmle", sl_library = "SL.glm",
                               return_steps = TRUE)
  expect_s3_class(plot(fit_steps, type = "clever"), "ggplot")
  expect_s3_class(plot(fit_steps, type = "ic"), "ggplot")
  bare <- structure(list(type = "tmle"), class = "tmle_fit")
  expect_error(plot(bare, type = "clever"), "targeting step")
})

test_that("merged arguments work: thresholds list, surface list, profile_vars", {
  sev <- .severe_cohort(500, seed = 51)
  lock <- create_analysis_lock(sev, "treatment", "outcome", c("x1", "x2"))
  psf <- fit_ps(lock, "glm")
  sup <- assess_support(psf, thresholds = list(flag_max_weight = 5),
                        tree_search = FALSE)
  expect_true(sup$summary$max_iptw_weight > 5)   # stricter threshold binds
  expect_false(is.null(sup$balance))
  expect_true(all(c("smd_unweighted", "smd_weighted") %in%
                    names(sup$balance)))
  p <- love_plot(sup)
  expect_s3_class(p, "ggplot")
  expect_error(assess_support(psf, thresholds = list(nope = 1)), "Unknown")
  fea <- estimand_feasibility(psf, profile_vars = c("x1", "x2"))
  expect_false(is.null(fea$unsupported))
  sim <- simulate_support(lock, psf,
                          surface = list(confounding = 0, modification = 0,
                                         effect = 0.1, base_rate = 0.15),
                          reps = 5, verbose = FALSE)
  expect_s3_class(sim, "support_simulation")
  expect_error(simulate_support(lock, psf, surface = list(bogus = 1),
                                reps = 5, verbose = FALSE), "Unknown")
  pb <- simulate_support(lock, psf, design = "parametric_bootstrap",
                         reps = 10)
  expect_true(isTRUE(pb$optimistic))
})
