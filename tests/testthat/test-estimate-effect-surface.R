# Workflow-verb surface coverage: every estimand and estimator branch
# of estimate_effect(), the external propensity path of fit_ps(), the
# optional components of design_report(), and the file/empty branches
# of export_design_log().

.surface_fixture <- function(n = 400L, seed = 31L) {
  dat  <- sim_func1(n = n, seed = seed)
  lock <- create_analysis_lock(dat, "treatment", "event_24",
                               c("age", "sex", "biomarker"), seed = seed,
                               negative_controls = "nc_outcome")
  ps <- fit_ps(lock, method = "glm")
  list(dat = dat, lock = lock, ps = ps)
}

test_that("estimate_effect covers the comparator estimators", {
  fx <- .surface_fixture()
  crude <- estimate_effect(fx$lock, estimand = "ATE",
                           estimator = "crude")
  iptw  <- estimate_effect(fx$lock, fx$ps, estimand = "ATE",
                           estimator = "iptw")
  match <- estimate_effect(fx$lock, fx$ps, estimand = "ATE",
                           estimator = "match")
  for (r in list(crude, iptw, match))
    expect_true(is.numeric(r$estimate) && is.finite(r$estimate))
})

test_that("estimate_effect covers the non-ATE estimands", {
  skip_on_cran()
  fx <- .surface_fixture()
  lib <- c("SL.glm")
  trimmed <- estimate_effect(fx$lock, fx$ps, estimand = "trimmed_ATE",
                             sl_library = lib, verbose = FALSE)
  expect_true(is.data.frame(trimmed$results) || is.list(trimmed))
  att <- estimate_effect(fx$lock, fx$ps, estimand = "ATT",
                         sl_library = lib)
  expect_true(is.numeric(att$estimate))
  ato <- estimate_effect(fx$lock, fx$ps, estimand = "ATO",
                         sl_library = lib)
  expect_true(is.numeric(ato$estimate))
  matt <- estimate_effect(fx$lock, fx$ps, estimand = "matched_ATT",
                          sl_library = lib)
  expect_true(is.numeric(matt$estimates$ATE$estimate))
  expect_gt(matt$n_matched, 10L)
})

test_that("estimate_effect requires a ps_fit where one is needed and covers IPCW", {
  skip_on_cran()
  fx <- .surface_fixture()
  expect_error(estimate_effect(fx$lock, estimand = "ATT"),
               "needs a ps_fit")
  dat_na <- fx$dat
  dat_na$event_24[seq_len(60)] <- NA
  lock_na <- create_analysis_lock(dat_na, "treatment", "event_24",
                                  c("age", "sex", "biomarker"),
                                  seed = 31)
  ps_na <- fit_ps(lock_na, method = "glm")
  ipcw <- estimate_effect(lock_na, ps_na, estimand = "ATE",
                          estimator = "tmle", missing = "ipcw",
                          sl_library = c("SL.glm"))
  expect_true(is.numeric(ipcw$estimate))
})

test_that("fit_ps external method wraps supplied scores and errors without them", {
  fx <- .surface_fixture()
  expect_error(fit_ps(fx$lock, method = "external"), "scores")
  ext <- fit_ps(fx$lock, method = "external", scores = fx$ps$ps)
  expect_s3_class(ext, "ps_fit")
  sup <- assess_support(ext)
  expect_true(sup$verdict %in% c("PASS", "FLAG", "SEVERE", "FAIL"))
})

test_that("design_report carries the optional simulation, NC, and DQ components", {
  skip_on_cran()
  fx <- .surface_fixture()
  lock <- create_analysis_lock(fx$dat, "treatment", "event_24",
                               c("age", "sex", "biomarker"), seed = 31,
                               negative_controls = "nc_outcome",
                               dq_thresholds = list(max_abs_bias = 0.02,
                                                    min_coverage = 0.90,
                                                    max_rmse_ratio = 1.5),
                               nc_criteria = list(null_band = 0.02))
  lock <- declare_estimand_ladder(lock, primary = "ATE",
                                  fallbacks = c("trimmed_ATE", "ATT"),
                                  trigger = "SEVERE")
  ps   <- fit_ps(lock, method = "glm")
  sup  <- assess_support(ps)
  feas <- estimand_feasibility(ps)
  sim  <- simulate_support(lock, ps, reps = 5L, verbose = FALSE)
  ncl  <- suppressWarnings(negative_control_ladder(
    lock, restrictions = list("younger" = fx$dat$age < median(fx$dat$age)),
    method = "unadjusted", verbose = FALSE))
  cand <- define_candidates("glm", g_library = "SL.glm",
                            truncation = 0.01)
  dq <- suppressMessages(stress_test(
    lock, cand, reps = 2L, verbose = FALSE,
    threats = list(covariate_missingness = list(fractions = 0.10))))
  rep <- design_report(lock, sup, feas, simulation = sim,
                       nc_ladder = ncl, dq = dq)
  expect_s3_class(rep, "design_report")
  expect_true(nzchar(rep$recommendation))
  expect_false(is.null(rep$nc_reading))
  expect_false(is.null(rep$dq_reading))
  expect_output(print(rep), "Support map passes")
  expect_output(print(rep), "DQ stress")
})

test_that("export_design_log handles the empty log and the file argument", {
  dat  <- sim_func1(n = 100, seed = 32)
  lock <- create_analysis_lock(dat, "treatment", "event_24",
                               c("age", "sex"), seed = 32)
  empty <- export_design_log(lock)
  expect_true(is.data.frame(empty))
  expect_equal(nrow(empty), 0L)

  lock <- declare_negative_controls(lock, "nc_outcome")
  f <- tempfile(fileext = ".csv")
  out <- export_design_log(lock, format = "muntner", file = f)
  expect_true(file.exists(f))
  back <- utils::read.csv(f, stringsAsFactors = FALSE)
  expect_equal(nrow(back), nrow(out))
  expect_true(all(c("date", "stage", "issue", "decision", "rationale",
                    "decided_by") %in% names(back)))
  unlink(f)
})