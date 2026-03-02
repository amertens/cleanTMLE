test_that("make_table1 works on data.frame", {
  dat <- sim_func1(n = 100, seed = 30)
  tbl <- make_table1(dat, vars = c("age", "sex", "biomarker"),
                     by = "treatment")

  expect_true(is.data.frame(tbl))
  expect_true("variable" %in% names(tbl))
  expect_true("mean" %in% names(tbl))
  expect_true("smd" %in% names(tbl))
})

test_that("make_table1 without grouping", {
  dat <- sim_func1(n = 50, seed = 31)
  tbl <- make_table1(dat, vars = c("age", "sex"))

  expect_true(is.data.frame(tbl))
  expect_equal(unique(tbl$group), "Overall")
})

test_that("make_table1 works on cumrisk object", {
  dat <- sim_func1(n = 100, seed = 32)
  spec <- specify_models(data = dat) |>
    identify_outcome(event, type = "time_to_event") |>
    identify_treatment(treatment, formula = ~ age + sex)

  fit <- estimate_ipwrisk(spec, risk_time = c(12))

  tbl <- make_table1(fit, vars = c("age", "sex"))
  expect_true(is.data.frame(tbl))
})

test_that("make_table2 works on cumrisk object", {
  dat <- sim_func1(n = 200, seed = 33)
  spec <- specify_models(data = dat) |>
    identify_outcome(event, type = "time_to_event") |>
    identify_treatment(treatment, formula = ~ age)

  fit <- estimate_ipwrisk(spec, risk_time = c(12, 24))
  tbl <- make_table2(fit, risk_time = 24)

  expect_true(is.data.frame(tbl))
  expect_true("cumulative_risk" %in% names(tbl))
  expect_equal(nrow(tbl), 2)
})

test_that("make_wt_summary_table works", {
  dat <- sim_func1(n = 200, seed = 34)
  spec <- specify_models(data = dat) |>
    identify_outcome(event, type = "time_to_event") |>
    identify_treatment(treatment, formula = ~ age + sex)

  fit <- estimate_ipwrisk(spec, risk_time = c(12))
  tbl <- make_wt_summary_table(fit)

  expect_true(is.data.frame(tbl))
  expect_true("mean" %in% names(tbl))
  expect_true("p50" %in% names(tbl))
  expect_equal(nrow(tbl), 2)
})

test_that("extreme_weights returns top-k", {
  dat <- sim_func1(n = 200, seed = 35)
  spec <- specify_models(data = dat) |>
    identify_outcome(event, type = "time_to_event") |>
    identify_treatment(treatment, formula = ~ age + sex + biomarker)

  fit <- estimate_ipwrisk(spec, risk_time = c(12))
  ew <- extreme_weights(fit, k = 5)

  expect_true(is.data.frame(ew))
  expect_equal(nrow(ew), 5)
  expect_true("weight" %in% names(ew))
})

test_that("inspect_ipw_weights extracts weights", {
  dat <- sim_func1(n = 100, seed = 36)
  spec <- specify_models(data = dat) |>
    identify_outcome(event, type = "time_to_event") |>
    identify_treatment(treatment, formula = ~ age)

  fit <- estimate_ipwrisk(spec, risk_time = c(12))

  w <- inspect_ipw_weights(fit, type = "iptw")
  expect_true(is.numeric(w))
  expect_length(w, 100)
})
