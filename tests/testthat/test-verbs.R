# The 0.3.0 verb wrappers: stress_test -> select_candidate round trip,
# and the negative_control_ladder verb. The engines behind them keep
# their own tests (test-plasmode_dq.R, test-design-tools.R); this file
# tests the wrapper contracts.

test_that("stress_test with no threats is the clean-data baseline and select_candidate works", {
  dat  <- sim_func1(n = 300, seed = 21)
  lock <- create_analysis_lock(dat, "treatment", "event_24",
                               c("age", "sex", "biomarker"), seed = 21)
  cands <- define_candidates(grid = list(
    truncations = c(0.01, 0.10),
    libraries   = list(glm = "SL.glm")))
  st <- stress_test(lock, cands, threats = NULL, reps = 3L,
                    verbose = FALSE)
  expect_s3_class(st, "ct_stress")
  expect_true(all(st$metrics$scenario == "none"))
  expect_equal(sort(unique(st$metrics$candidate)),
               sort(names(cands)))
  # The result carries the candidate specifications, so selection can
  # run on it directly.
  expect_length(st$tmle_candidates, 2L)

  best <- select_candidate(st, rule = "min_rmse")
  expect_s3_class(best, "tmle_selected_spec")
  expect_true(best$candidate_id %in% names(cands))
})

test_that("a baseline-only run on a thresholds-carrying lock has no verdict", {
  dat  <- sim_func1(n = 300, seed = 25)
  lock <- create_analysis_lock(dat, "treatment", "event_24",
                               c("age", "sex", "biomarker"), seed = 25,
                               dq_thresholds = list(max_abs_bias = 0.02,
                                                    min_coverage = 0.90,
                                                    max_rmse_ratio = 1.5))
  cand <- define_candidates("glm", g_library = "SL.glm",
                            truncation = 0.01)
  st <- stress_test(lock, cand, threats = NULL, reps = 2L,
                    verbose = FALSE)
  expect_null(st$verdict)
  expect_null(st$tipping)
  expect_false(is.null(st$thresholds))
  expect_output(print(st), "Baseline-only run")
})

test_that("select_candidate min_max_rmse uses the degraded rows", {
  dat  <- sim_func1(n = 300, seed = 22)
  lock <- create_analysis_lock(dat, "treatment", "event_24",
                               c("age", "sex", "biomarker"), seed = 22)
  cands <- define_candidates(grid = list(
    truncations = c(0.01, 0.10),
    libraries   = list(glm = "SL.glm")))
  st <- suppressMessages(stress_test(
    lock, cands, reps = 2L, verbose = FALSE,
    threats = list(covariate_missingness = list(fractions = 0.20))))
  expect_true(any(st$metrics$scenario != "none"))
  best <- select_candidate(st, rule = "min_max_rmse")
  expect_s3_class(best, "tmle_selected_spec")
})

test_that("select_candidate rejects foreign objects and empty baselines", {
  expect_error(select_candidate(list()), "stress_test")
  fake <- structure(list(
    metrics = data.frame(scenario = "cov_miss", candidate = "a",
                         stringsAsFactors = FALSE),
    tmle_candidates = list()), class = "plasmode_dq_results")
  expect_error(select_candidate(fake), "no baseline")
})

test_that("stress_test on an external_pilot lock runs masked", {
  dat  <- sim_func1(n = 300, seed = 23)
  lock <- create_analysis_lock(dat, "treatment", "event_24",
                               c("age", "sex", "biomarker"), seed = 23,
                               dgp_mode = "external_pilot")
  masked <- mask_outcome(lock)
  cand <- define_candidates("glm", g_library = "SL.glm",
                            truncation = 0.01)
  st <- suppressMessages(stress_test(
    masked, cand, reps = 2L, verbose = FALSE,
    pilot_q0 = function(W) stats::plogis(-1 + 0.01 * W$age)))
  expect_s3_class(st, "ct_stress")
  expect_identical(st$dgp_mode, "external_pilot")
})

test_that("negative_control_ladder wraps the engine and classes the result", {
  dat  <- sim_func1(n = 400, seed = 24)
  lock <- create_analysis_lock(dat, "treatment", "event_24",
                               c("age", "sex", "biomarker"), seed = 24,
                               negative_controls = "nc_outcome")
  ncl <- suppressWarnings(negative_control_ladder(
    lock,
    restrictions = list("biomarker below median" =
                          dat$biomarker < stats::median(dat$biomarker)),
    method = "unadjusted", verbose = FALSE))
  expect_s3_class(ncl, "ct_nc_ladder")
  expect_s3_class(ncl, "nc_ladder")
  expect_true(is.data.frame(as.data.frame(ncl)))
  expect_true(all(c("full cohort", "biomarker below median") %in%
                    ncl$table$cohort))
})

test_that("a threat list naming only the MNAR mechanism runs only MNAR", {
  # `$` on the scenario list partially matched covariate_missingness to
  # covariate_missingness_mnar, so declaring MNAR alone also ran MCAR.
  dat  <- sim_func1(n = 300, seed = 26)
  lock <- create_analysis_lock(dat, "treatment", "event_24",
                               c("age", "sex", "biomarker"), seed = 26)
  cand <- define_candidates("glm", g_library = "SL.glm",
                            truncation = 0.01)
  st <- suppressMessages(stress_test(
    lock, cand, reps = 2L, verbose = FALSE,
    threats = list(covariate_missingness_mnar = list(fractions = 0.10))))
  expect_setequal(unique(st$metrics$scenario), c("none", "cov_miss_mnar"))
})
