test_that("a binary plan expands 'default' and fills defaults", {
  p <- plan_fx(outcome = "y", estimands = c("ATE", "ATO"))
  expect_s3_class(p, "cr_plan")
  expect_identical(p$outcome_type, "binary")
  expect_identical(p$ps_library, c("glmnet", "earth", "xgboost", "nnet"))
  expect_identical(p$candidates$library$sl, c("glmnet", "earth", "xgboost", "nnet"))
  expect_identical(p$candidates$library$glm, "glm")
  expect_equal(p$tolerance, list(bias = 0.02, coverage = 0.90))
  expect_equal(p$surfaces$baseline_risk, 0.10)
  expect_identical(p$K, 5L)
  expect_identical(p$max_reps, 1000L)
})

test_that("plan validation rejects what the spec rules out", {
  expect_error(plan_fx("y", "ATX"), "estimands")
  expect_error(plan_fx("y", c("ATE", "trimmed_ATE")), "trim_band")
  expect_error(plan_fx("y", "ATE", ps_library = "SL.glm"), "Unknown learner")
  expect_error(plan_fx(c(time = "t", event = "s"), "ATE", target_time = 365),
               "hazards")
  expect_error(plan_fx(c(time = "t", event = "s"), c("ATE", "ATT"),
                             target_time = 365, hazards = list(`1` = "coxnet")),
               "not available for a time-to-event")
  expect_error(plan_fx("y", "ATE", candidates = list(
    library = list(glm = "glm"), truncation = 0.6)), "truncation")
})

test_that("negative controls get default criteria", {
  p <- plan_fx("y", "ATE", negative_controls = c(nc_visit = "care use"))
  expect_equal(p$nc_criteria$null_band, c(-0.02, 0.02))
  expect_identical(p$nc_criteria$rule, "interval")
})

test_that(".candidates crosses libraries with truncation levels", {
  p <- plan_fx("y", "ATE", candidates = list(
    library = list(glm = "glm", sl = c("glm", "glmnet")), truncation = c(0.05, 0.01)))
  ids <- vapply(.candidates(p), `[[`, "", "id")
  expect_identical(ids, c("glm_t0.01", "glm_t0.05", "sl_t0.01", "sl_t0.05"))
})

test_that("missing or non-whole numbers are refused, not passed through", {
  expect_error(plan_fx("y", "ATE", K = NA), "`K`")
  expect_error(plan_fx("y", "ATE", seed = NA), "`seed`")
  expect_error(plan_fx("y", "ATE", K = 2.7), "whole number")
  expect_error(plan_fx("y", "ATE", V = NA_integer_), "`V`")
  expect_error(plan_fx("y", "ATE", reps = c(10, 20)), "`reps`")
  expect_error(plan_fx("y", "ATE", K = 0), "K >= 1")
  expect_error(plan_fx("y", "ATE", tolerance = list(bias = NA)), "tolerance")
  expect_error(plan_fx("y", "ATE", surfaces = list(log_or = NA_real_)), "log_or")
})

test_that("max_reps is a whole number no smaller than reps", {
  expect_error(plan_fx("y", "ATE", reps = 50, max_reps = 40), "max_reps >= reps")
  expect_error(plan_fx("y", "ATE", max_reps = 300.5), "`max_reps`")
  expect_error(plan_fx("y", "ATE", max_reps = NA), "`max_reps`")
  expect_identical(plan_fx("y", "ATE", reps = 50, max_reps = 50)$max_reps, 50L)
})

test_that("surface vectors are validated and the grid size is announced", {
  expect_error(plan_fx("y", "ATE", surfaces = list(prognostic = 0)), "prognostic")
  expect_error(plan_fx("y", "ATE", surfaces = list(heterogeneity = c(0.5, -1))),
               "heterogeneity")
  p <- plan_fx("y", "ATE", surfaces = list(log_or = c(0.5, 0.5, 1)))
  expect_identical(p$surfaces$log_or, c(0.5, 1))
  expect_identical(.n_surfaces(p$surfaces), 8L)
  expect_message(plan_fx("y", "ATE", surfaces = list(log_or = c(0, 0.5, 1),
                                                         prognostic = c(0.5, 1, 2))),
                 "36 outcome surfaces")
  expect_silent(plan_fx("y", "ATE"))
})

test_that("a ci_halfwidth tolerance must be one positive number", {
  expect_error(plan_fx("y", "ATE", tolerance = list(ci_halfwidth = -1)), "ci_halfwidth")
  expect_error(plan_fx("y", "ATE", tolerance = list(ci_halfwidth = c(0.1, 0.2))),
               "ci_halfwidth")
  expect_equal(plan_fx("y", "ATE", tolerance = list(ci_halfwidth = 0.05))$tolerance$ci_halfwidth,
               0.05)
  expect_null(plan_fx("y", "ATE")$tolerance$ci_halfwidth)
})

test_that("the negative-control rule accepts point, ci and interval, and defaults to interval", {
  nc <- c(nc_visit = "care use")
  for (r in c("point", "ci", "interval"))
    expect_identical(plan_fx("y", "ATE", negative_controls = nc,
                                   nc_criteria = list(rule = r))$nc_criteria$rule, r)
  expect_error(plan_fx("y", "ATE", negative_controls = nc,
                             nc_criteria = list(rule = "bayes")), "point.*ci.*interval")
})

test_that("nc_exclude is validated against the declared controls", {
  nc <- c(nc_visit = "care use", nc_bad = "other")
  p <- plan_fx("y", "ATE", negative_controls = nc,
                     nc_exclude = list(nc_visit = c("w1", "w1", "w2")))
  expect_identical(p$nc_exclude, list(nc_visit = c("w1", "w2")))
  expect_null(plan_fx("y", "ATE", negative_controls = nc)$nc_exclude)
  expect_null(plan_fx("y", "ATE", negative_controls = nc, nc_exclude = list())$nc_exclude)
  expect_error(plan_fx("y", "ATE", negative_controls = nc,
                             nc_exclude = list(nc_other = "w1")), "not declared negative controls: `nc_other`")
  expect_error(plan_fx("y", "ATE", negative_controls = nc,
                             nc_exclude = list(a = "w1", b = "w1")), "`a`, `b`")
  expect_error(plan_fx("y", "ATE", nc_exclude = list(nc_visit = "w1")),
               "not declared negative controls")
  expect_error(plan_fx("y", "ATE", negative_controls = nc, nc_exclude = "w1"),
               "named list")
  expect_error(plan_fx("y", "ATE", negative_controls = nc,
                             nc_exclude = list(nc_visit = character())), "character vector")
  expect_error(plan_fx("y", "ATE", negative_controls = nc,
                             nc_exclude = list(nc_visit = 1)), "character vector")
  expect_error(plan_fx("y", "ATE", negative_controls = nc,
                             nc_exclude = list(nc_visit = "w1", nc_visit = "w2")), "named list")
})

test_that("the plan print shows nc_exclude when it is declared", {
  nc <- c(nc_visit = "care use")
  shown <- capture.output(print(plan_fx("y", "ATE", negative_controls = nc,
                                              nc_exclude = list(nc_visit = "w1"))))
  expect_true(any(grepl("^  NC excluded: +nc_visit: w1$", shown)))
  expect_false(any(grepl("left out", capture.output(
    print(plan_fx("y", "ATE", negative_controls = nc))))))
})

# time zero and eligibility
test_that("analysis_plan requires and normalizes time_zero", {
  expect_error(analysis_plan("y", "ATE"), "`time_zero` is required")
  p <- analysis_plan("y", "ATE", time_zero = "treatment assignment")
  expect_identical(p$time_zero$event, "treatment assignment")
  expect_identical(p$time_zero$follow_up_from, "treatment assignment")
  p2 <- analysis_plan("y", "ATE", time_zero = list(event = "transport decision",
                      measured_as = "time of injury"))
  expect_identical(p2$time_zero$measured_as, "time of injury")
})
test_that("follow-up from another event needs a response", {
  tz <- list(event = "transport decision", follow_up_from = "arrival")
  expect_error(analysis_plan("y", "ATE", time_zero = tz), "Follow-up starts at `arrival`, not at time zero")
  expect_error(analysis_plan("y", "ATE", time_zero = c(tz, response = "  ")), "not at time zero")
  expect_silent(analysis_plan("y", "ATE", time_zero = c(tz, response = "Arrival is minutes after the decision.")))
})
test_that("eligibility assessed after time zero needs a response", {
  el <- list(arrived_alive = list(description = "Arrived alive", assessed = "after"))
  expect_error(analysis_plan("y", "ATE", time_zero = "t0", eligibility = el),
               "Eligibility criterion `arrived_alive` is assessed after time zero")
  el$arrived_alive$response <- "Bias analysis registered before unblinding."
  expect_identical(analysis_plan("y", "ATE", time_zero = "t0", eligibility = el)$eligibility$arrived_alive$assessed, "after")
  expect_error(analysis_plan("y", "ATE", time_zero = "t0",
               eligibility = list(x = list(description = "d", assessed = "later"))), "before, at, after")
  expect_error(analysis_plan("y", "ATE", time_zero = "t0", eligibility = list(list(description = "d", assessed = "at"))),
               "named list")
})
test_that("the plan table shows time zero and eligibility", {
  p <- plan_fx(outcome = "y", estimands = "ATE", eligibility = list(adult = list(description = "Age 18+", assessed = "at")))
  tab <- cleanTMLE:::.plan_table(cleanTMLE:::.plan_summary(p))
  expect_true(all(c("Time zero", "Eligibility (when assessed)") %in% tab$setting))
})
