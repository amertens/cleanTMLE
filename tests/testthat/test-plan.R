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
  expect_identical(p$nc_criteria$rule, "point")
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
