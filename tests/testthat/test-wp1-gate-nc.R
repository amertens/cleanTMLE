# WP1 corrections: the pure locked DQ verdict and tipping points, the
# locked negative-control criteria, the content fingerprint, the design
# log export, and the design report hygiene.

# ── Fixtures: the recorded Scenario A and Scenario C stress rows ─────────
# Provenance: results_new/dq_stress_good_overlap.rds (selected candidate
# "aggressive") and results_new/dq_stress_unmeasured_conf.rds (selected
# candidate "robust"), unmeasured-confounding rows at U odds ratios
# 2, 3, 4, 6, 8 with the clean baseline, as recorded by the June/July
# 2026 runs and documented in revision/STAGE0_AUDIT.md section 2.3.
# These are test DATA (a fixed fixture), not results: the test asserts
# what the corrected pure rule returns on them.
.dq_fixture <- function(candidate, baseline, rows) {
  lv <- sprintf("OR_trt%.1f_out%.1f", c(2, 3, 4, 6, 8), c(2, 3, 4, 6, 8))
  metrics <- rbind(
    data.frame(scenario = "none", level = "0", effect_size = 0.05,
               candidate = candidate, bias = baseline[1],
               rmse = baseline[2], coverage = baseline[3],
               stringsAsFactors = FALSE),
    data.frame(scenario = "unmeasured_U", level = lv, effect_size = 0.05,
               candidate = candidate, bias = rows$bias, rmse = rows$rmse,
               coverage = rows$coverage, stringsAsFactors = FALSE))
  structure(list(metrics = metrics, reps = 30L,
                 baseline = metrics[metrics$scenario == "none", ]),
            class = "plasmode_dq_results")
}
fx_A <- .dq_fixture("aggressive", c(0.0028, 0.019, 0.933),
  list(bias = c(0.01032, 0.02640, 0.04165, 0.07562, 0.09920),
       rmse = c(0.02028, 0.03351, 0.04394, 0.07672, 0.10037),
       coverage = c(0.900, 0.600, 0.300, 0.000, 0.000)))
fx_C <- .dq_fixture("robust", c(-0.00252, 0.01582, 0.967),
  list(bias = c(0.01348, 0.03254, 0.05339, 0.08502, 0.10326),
       rmse = c(0.02599, 0.03940, 0.05727, 0.08726, 0.10484),
       coverage = c(0.833, 0.433, 0.267, 0.000, 0.000)))
.th <- list(max_abs_bias = 0.02, min_coverage = 0.90, max_rmse_ratio = 1.5)

test_that("the same locked grid reads STOP in Scenario A and Scenario C", {
  # Decision D2: the verdict is a pure function of the locked thresholds
  # and the declared grid; applied to the recorded fixtures it returns
  # STOP for both scenarios, and no grid tuning is performed to separate
  # them.
  vA <- dq_locked_verdict(fx_A, .th)
  vC <- dq_locked_verdict(fx_C, .th)
  expect_identical(vA$verdict, "STOP")
  expect_identical(vC$verdict, "STOP")
  expect_true(grepl("bias", vA$breached))
  expect_true(grepl("bias", vC$breached))
})

test_that("tipping points are the primary reading and differ by scenario", {
  tA <- dq_tipping_points(fx_A, .th)
  tC <- dq_tipping_points(fx_C, .th)
  # Scenario A first crosses a locked bound at OR 3.0 (bias 0.0264 over
  # the 0.02 bound; OR 2.0 coverage is exactly 0.90, not under it).
  expect_identical(tA$tipping_level, "OR_trt3.0_out3.0")
  expect_identical(tA$tipping_or, 3)
  expect_identical(tA$tipping_level_bias, "OR_trt3.0_out3.0")
  # Scenario C first crosses at OR 2.0, through the coverage floor
  # (0.833 < 0.90), while its bias tipping point is also OR 3.0.
  expect_identical(tC$tipping_level, "OR_trt2.0_out2.0")
  expect_identical(tC$tipping_or, 2)
  expect_identical(tC$tipping_level_bias, "OR_trt3.0_out3.0")
  expect_identical(tC$tipping_level_coverage, "OR_trt2.0_out2.0")
})

test_that("the verdict is pure: nothing outside metrics and thresholds enters", {
  v1 <- dq_locked_verdict(fx_A, .th)
  fx_A2 <- fx_A
  fx_A2$anything_else <- "changed"
  v2 <- dq_locked_verdict(fx_A2, .th)
  expect_identical(v1, v2)
  # On a mild grid (the OR 2.0 cell only) the same thresholds read GO:
  # bias 0.0103 is under 0.02 and coverage 0.90 meets the floor.
  fx_mild <- fx_A
  fx_mild$metrics <- fx_mild$metrics[fx_mild$metrics$level %in%
                                       c("0", "OR_trt2.0_out2.0"), ]
  v3 <- dq_locked_verdict(fx_mild, list(max_abs_bias = 0.02,
                                        min_coverage = 0.90,
                                        max_rmse_ratio = 1.5))
  expect_identical(v3$verdict, "GO")
})

# ── Locked negative-control criteria ─────────────────────────────────────

.nc_row <- function(cohort, nc, domain, est, lo, hi,
                    status = "estimated", method = "tmle") {
  data.frame(cohort = cohort, negative_control = nc, domain = domain,
             n = 100L, estimate = est, ci_lower = lo, ci_upper = hi,
             p_value = 0.5, flagged = FALSE, status = status,
             method = method, stringsAsFactors = FALSE)
}

test_that("nc_ladder_verdict reads only the locked criteria", {
  tab <- rbind(
    .nc_row("full cohort", "nc1", "ses",       0.010, -0.005, 0.025),
    .nc_row("full cohort", "nc2", "ses",       0.015,  0.001, 0.030),
    .nc_row("full cohort", "nc3", "indication", 0.050,  0.030, 0.070),
    .nc_row("full cohort", "nc4", "rare", NA, NA, NA,
            status = "inestimable (a cell below 5 events)"))
  cr <- list(null_band = 0.02)
  v <- nc_ladder_verdict(tab, cr)
  # ses passes on the point rule (both under 0.02), indication fails,
  # rare is insufficient: the rung verdict is STOP (a domain failed).
  expect_identical(
    v$by_domain$reading[v$by_domain$domain == "ses"], "pass")
  expect_identical(
    v$by_domain$reading[v$by_domain$domain == "indication"], "fail")
  expect_identical(
    v$by_domain$reading[v$by_domain$domain == "rare"], "insufficient")
  expect_identical(v$by_rung$verdict, "STOP")

  # The CI rule is stricter: nc1's interval leaves the band, so ses
  # fails under ci_in_band while passing under the point rule.
  v_ci <- nc_ladder_verdict(tab, list(null_band = 0.02,
                                      rule = "ci_in_band"))
  expect_identical(
    v_ci$by_domain$reading[v_ci$by_domain$domain == "ses"], "fail")

  # Majority consistency: one of two in band fails "all" but not
  # "majority" needs a strict majority, so exactly half still fails.
  tab2 <- rbind(
    .nc_row("full cohort", "nc1", "ses", 0.010, -0.005, 0.018),
    .nc_row("full cohort", "nc2", "ses", 0.050,  0.030, 0.070),
    .nc_row("full cohort", "nc5", "ses", 0.005, -0.010, 0.018))
  v_all <- nc_ladder_verdict(tab2, list(null_band = 0.02))
  v_maj <- nc_ladder_verdict(tab2, list(null_band = 0.02,
                                        consistency = "majority_within_band"))
  expect_identical(v_all$by_rung$verdict, "STOP")
  expect_identical(v_maj$by_rung$verdict, "GO")

  # min_per_domain above the estimable count reads insufficient (FLAG).
  v_min <- nc_ladder_verdict(tab2, list(null_band = 0.10,
                                        min_per_domain = 4L))
  expect_identical(v_min$by_rung$verdict, "FLAG")
})

test_that("the ladder carries method and domain and grades from the lock", {
  dat  <- sim_func1(n = 400, seed = 42)
  lock <- create_analysis_lock(dat, "treatment", "event_24",
                               c("age", "sex", "biomarker"), seed = 7,
                               negative_controls = "nc_outcome",
                               nc_criteria = list(null_band = 0.05))
  rungs <- list("older half" = dat$age >= stats::median(dat$age))
  ncl <- run_negative_control_ladder(lock, rungs, method = "unadjusted",
                                     verbose = FALSE)
  expect_true(all(c("method", "domain") %in% names(ncl$table)))
  expect_true(all(ncl$table$method == "unadjusted"))
  expect_false(is.null(ncl$verdict))
  expect_true(all(ncl$verdict$by_rung$verdict %in% c("GO", "FLAG", "STOP")))
  expect_output(print(ncl), "Check Point 3")

  # Without locked criteria there is no verdict, and the print says so.
  lock2 <- create_analysis_lock(dat, "treatment", "event_24",
                                c("age", "sex", "biomarker"), seed = 7,
                                negative_controls = "nc_outcome")
  ncl2 <- run_negative_control_ladder(lock2, rungs, method = "unadjusted",
                                      verbose = FALSE)
  expect_null(ncl2$verdict)
})

test_that("run_negative_control_tmle never falls back silently", {
  # Force a Q-fit failure with an sl_library of nonexistent learners.
  dat  <- sim_func1(n = 200, seed = 42)
  lock <- create_analysis_lock(dat, "treatment", "event_24",
                               c("age", "sex", "biomarker"), seed = 7,
                               negative_controls = "nc_outcome")
  ps <- fit_ps(lock, method = "glm")
  expect_error(
    suppressWarnings(run_negative_control_tmle(
      lock, "nc_outcome", ps, sl_library = "SL.does.not.exist")),
    "no fallback")
})

# ── Content fingerprint ──────────────────────────────────────────────────

test_that("the lock fingerprint commits to design-data content, not the outcome", {
  dat  <- sim_func1(n = 200, seed = 42)
  covs <- c("age", "sex", "biomarker")
  lock <- create_analysis_lock(dat, "treatment", "event_24", covs, seed = 7)
  expect_message(validate_analysis_lock(lock), "validated")

  # Tampering with a covariate value invalidates the fingerprint.
  tampered <- lock
  tampered$data$age[1] <- tampered$data$age[1] + 100
  expect_error(validate_analysis_lock(tampered), "hash mismatch")

  # Permuting or masking the outcome does not: the fingerprint commits
  # to nothing about outcome values.
  perm <- dat; perm$event_24 <- rev(perm$event_24)
  lock_perm <- create_analysis_lock(perm, "treatment", "event_24", covs,
                                    seed = 7)
  expect_identical(lock$lock_hash, lock_perm$lock_hash)
  expect_message(validate_analysis_lock(mask_outcome(lock)), "validated")
})

test_that("lock declarations are validated and fingerprinted", {
  dat  <- sim_func1(n = 200, seed = 42)
  covs <- c("age", "sex", "biomarker")
  expect_error(create_analysis_lock(dat, "treatment", "event_24", covs,
                                    nc_criteria = list(bad = 1)),
               "null_band")
  expect_error(create_analysis_lock(dat, "treatment", "event_24", covs,
                                    dq_thresholds = list(max_abs_bias = 0.02)),
               "dq_thresholds")
  expect_error(create_analysis_lock(dat, "treatment", "event_24", covs,
                                    nc_criteria = list(null_band = 0.02,
                                                       scale = "ratio")),
               "risk-difference")
  l1 <- create_analysis_lock(dat, "treatment", "event_24", covs, seed = 7)
  l2 <- create_analysis_lock(dat, "treatment", "event_24", covs, seed = 7,
                             dq_thresholds = list(max_abs_bias = 0.02,
                                                  min_coverage = 0.9,
                                                  max_rmse_ratio = 1.5))
  expect_false(identical(l1$lock_hash, l2$lock_hash))
})

# ── Design log export ────────────────────────────────────────────────────

test_that("export_design_log writes the Muntner Table S1 columns", {
  dat  <- sim_func1(n = 200, seed = 42)
  lock <- create_analysis_lock(dat, "treatment", "event_24",
                               c("age", "sex", "biomarker"), seed = 7)
  lock <- declare_estimand_ladder(lock, primary = "ATE",
                                  fallbacks = c("trimmed_ATE", "ATT"),
                                  trigger = "SEVERE")
  m <- export_design_log(lock, format = "muntner")
  expect_identical(names(m), c("date", "stage", "issue", "decision",
                               "rationale", "decided_by"))
  expect_true(nrow(m) >= 1)
  tidy <- export_design_log(lock, format = "tidy")
  expect_true(all(c("timestamp", "type", "note", "stage", "decision",
                    "rationale", "decided_by") %in% names(tidy)))
})

# ── Design report hygiene and DQ wiring ──────────────────────────────────

test_that("design_report stores summaries only and reads the locked verdicts", {
  dat  <- sim_func1(n = 300, seed = 42)
  lock <- create_analysis_lock(dat, "treatment", "event_24",
                               c("age", "sex", "biomarker"), seed = 7,
                               roles = list(analyst = "A. Mertens",
                                            review_team = "TBD"))
  lock <- declare_estimand_ladder(lock, primary = "ATE",
                                  fallbacks = c("trimmed_ATE", "ATT"),
                                  trigger = "SEVERE")
  ps  <- fit_ps(lock, method = "glm")
  sup <- assess_support(ps)
  fea <- estimand_feasibility(ps)

  dq <- fx_A
  dq$verdict <- dq_locked_verdict(dq, .th)
  dq$tipping <- dq_tipping_points(dq, .th)

  dr <- design_report(lock, sup, fea, dq = dq)
  expect_null(dr$support$g)
  expect_null(dr$support$A)
  expect_true(isTRUE(dr$support$individual_vectors_removed))
  expect_true(grepl("STOP", dr$recommendation))
  expect_true(grepl("tipping point: OR 3.0", dr$recommendation))
  expect_output(print(dr), "Roles:")
  expect_output(print(dr), "DQ stress:")
  # The original support object keeps its vectors for plotting.
  expect_identical(length(sup$g), nrow(dat))
})