test_that("a binary plan expands 'default' and fills defaults", {
  p <- analysis_plan(outcome = "y", estimands = c("ATE", "ATO"))
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
  expect_error(analysis_plan("y", "ATX"), "estimands")
  expect_error(analysis_plan("y", c("ATE", "trimmed_ATE")), "trim_band")
  expect_error(analysis_plan("y", "ATE", ps_library = "SL.glm"), "Unknown learner")
  expect_error(analysis_plan(c(time = "t", event = "s"), "ATE", target_time = 365),
               "hazards")
  expect_error(analysis_plan(c(time = "t", event = "s"), c("ATE", "ATT"),
                             target_time = 365, hazards = list(`1` = "coxnet")),
               "not available for a time-to-event")
  expect_error(analysis_plan("y", "ATE", candidates = list(
    library = list(glm = "glm"), truncation = 0.6)), "truncation")
})

test_that("negative controls get default criteria", {
  p <- analysis_plan("y", "ATE", negative_controls = c(nc_visit = "care use"))
  expect_equal(p$nc_criteria$null_band, c(-0.02, 0.02))
  expect_identical(p$nc_criteria$rule, "interval")
})

test_that(".candidates crosses libraries with truncation levels", {
  p <- analysis_plan("y", "ATE", candidates = list(
    library = list(glm = "glm", sl = c("glm", "glmnet")), truncation = c(0.05, 0.01)))
  ids <- vapply(.candidates(p), `[[`, "", "id")
  expect_identical(ids, c("glm_t0.01", "glm_t0.05", "sl_t0.01", "sl_t0.05"))
})

test_that("missing or non-whole numbers are refused, not passed through", {
  expect_error(analysis_plan("y", "ATE", K = NA), "`K`")
  expect_error(analysis_plan("y", "ATE", seed = NA), "`seed`")
  expect_error(analysis_plan("y", "ATE", K = 2.7), "whole number")
  expect_error(analysis_plan("y", "ATE", V = NA_integer_), "`V`")
  expect_error(analysis_plan("y", "ATE", reps = c(10, 20)), "`reps`")
  expect_error(analysis_plan("y", "ATE", K = 0), "K >= 1")
  expect_error(analysis_plan("y", "ATE", tolerance = list(bias = NA)), "tolerance")
  expect_error(analysis_plan("y", "ATE", surfaces = list(log_or = NA_real_)), "log_or")
})

test_that("max_reps is a whole number no smaller than reps", {
  expect_error(analysis_plan("y", "ATE", reps = 50, max_reps = 40), "max_reps >= reps")
  expect_error(analysis_plan("y", "ATE", max_reps = 300.5), "`max_reps`")
  expect_error(analysis_plan("y", "ATE", max_reps = NA), "`max_reps`")
  expect_identical(analysis_plan("y", "ATE", reps = 50, max_reps = 50)$max_reps, 50L)
})

test_that("surface vectors are validated and the grid size is announced", {
  expect_error(analysis_plan("y", "ATE", surfaces = list(prognostic = 0)), "prognostic")
  expect_error(analysis_plan("y", "ATE", surfaces = list(heterogeneity = c(0.5, -1))),
               "heterogeneity")
  p <- analysis_plan("y", "ATE", surfaces = list(log_or = c(0.5, 0.5, 1)))
  expect_identical(p$surfaces$log_or, c(0.5, 1))
  expect_identical(.n_surfaces(p$surfaces), 8L)
  expect_message(analysis_plan("y", "ATE", surfaces = list(log_or = c(0, 0.5, 1),
                                                         prognostic = c(0.5, 1, 2))),
                 "36 outcome surfaces")
  expect_silent(analysis_plan("y", "ATE"))
})

test_that("a ci_halfwidth tolerance must be one positive number", {
  expect_error(analysis_plan("y", "ATE", tolerance = list(ci_halfwidth = -1)), "ci_halfwidth")
  expect_error(analysis_plan("y", "ATE", tolerance = list(ci_halfwidth = c(0.1, 0.2))),
               "ci_halfwidth")
  expect_equal(analysis_plan("y", "ATE", tolerance = list(ci_halfwidth = 0.05))$tolerance$ci_halfwidth,
               0.05)
  expect_null(analysis_plan("y", "ATE")$tolerance$ci_halfwidth)
})

test_that("the negative-control rule accepts point, ci and interval, and defaults to interval", {
  nc <- c(nc_visit = "care use")
  for (r in c("point", "ci", "interval"))
    expect_identical(analysis_plan("y", "ATE", negative_controls = nc,
                                   nc_criteria = list(rule = r))$nc_criteria$rule, r)
  expect_error(analysis_plan("y", "ATE", negative_controls = nc,
                             nc_criteria = list(rule = "bayes")), "point.*ci.*interval")
})

test_that("nc_exclude is validated against the declared controls", {
  nc <- c(nc_visit = "care use", nc_bad = "other")
  p <- analysis_plan("y", "ATE", negative_controls = nc,
                     nc_exclude = list(nc_visit = c("w1", "w1", "w2")))
  expect_identical(p$nc_exclude, list(nc_visit = c("w1", "w2")))
  expect_null(analysis_plan("y", "ATE", negative_controls = nc)$nc_exclude)
  expect_null(analysis_plan("y", "ATE", negative_controls = nc, nc_exclude = list())$nc_exclude)
  expect_error(analysis_plan("y", "ATE", negative_controls = nc,
                             nc_exclude = list(nc_other = "w1")), "not declared negative controls: `nc_other`")
  expect_error(analysis_plan("y", "ATE", negative_controls = nc,
                             nc_exclude = list(a = "w1", b = "w1")), "`a`, `b`")
  expect_error(analysis_plan("y", "ATE", nc_exclude = list(nc_visit = "w1")),
               "not declared negative controls")
  expect_error(analysis_plan("y", "ATE", negative_controls = nc, nc_exclude = "w1"),
               "named list")
  expect_error(analysis_plan("y", "ATE", negative_controls = nc,
                             nc_exclude = list(nc_visit = character())), "character vector")
  expect_error(analysis_plan("y", "ATE", negative_controls = nc,
                             nc_exclude = list(nc_visit = 1)), "character vector")
  expect_error(analysis_plan("y", "ATE", negative_controls = nc,
                             nc_exclude = list(nc_visit = "w1", nc_visit = "w2")), "named list")
})

test_that("the plan print shows nc_exclude when it is declared", {
  nc <- c(nc_visit = "care use")
  shown <- capture.output(print(analysis_plan("y", "ATE", negative_controls = nc,
                                              nc_exclude = list(nc_visit = "w1"))))
  expect_true(any(grepl("^  NC excluded: +nc_visit: w1$", shown)))
  expect_false(any(grepl("left out", capture.output(
    print(analysis_plan("y", "ATE", negative_controls = nc))))))
})
