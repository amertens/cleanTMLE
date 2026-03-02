test_that("specify_models creates a cr_spec object", {
  dat <- data.frame(x = 1:10, y = rnorm(10))
  spec <- specify_models(data = dat)

  expect_s3_class(spec, "cr_spec")
  expect_equal(nrow(spec$data), 10)
  expect_null(spec$outcome)
  expect_null(spec$treatment)
})

test_that("specify_models rejects non-data.frame", {
  expect_error(specify_models(data = "not_a_df"), "data.frame")
})

test_that("identify_outcome adds outcome info", {
  dat <- data.frame(event = c(0, 1, 0, 1), time = c(1, 2, 3, 4))
  spec <- specify_models(data = dat) |>
    identify_outcome(event, type = "time_to_event")

  expect_equal(spec$outcome$name, "event")
  expect_equal(spec$outcome$type, "time_to_event")
})

test_that("identify_treatment adds treatment info", {
  dat <- data.frame(trt = c(0, 1, 0, 1), x = rnorm(4))
  spec <- specify_models(data = dat) |>
    identify_treatment(trt, formula = ~ x, model = "glm")

  expect_equal(spec$treatment$name, "trt")
  expect_equal(spec$treatment$model, "glm")
  expect_true(inherits(spec$treatment$formula, "formula"))
})

test_that("identify_censoring appends to censoring list", {
  dat <- data.frame(cens = c(0, 1, 0, 0), x = rnorm(4))
  spec <- specify_models(data = dat) |>
    identify_censoring(cens, formula = ~ x, model = "glm")

  expect_length(spec$censoring, 1)
  expect_equal(spec$censoring[[1]]$name, "cens")
})

test_that("identify_competing_risk stores event_value", {
  dat <- data.frame(etype = c(0, 1, 2, 1))
  spec <- specify_models(data = dat) |>
    identify_competing_risk(etype, event_value = 1)

  expect_equal(spec$competing_risk$name, "etype")
  expect_equal(spec$competing_risk$event_value, 1)
})

test_that("identify_subject stores subject ID", {
  dat <- data.frame(id = 1:5, x = rnorm(5))
  spec <- specify_models(data = dat) |>
    identify_subject(id)

  expect_equal(spec$subject$name, "id")
})

test_that("identify_interval stores start and stop", {
  dat <- data.frame(t0 = 0:4, t1 = 1:5)
  spec <- specify_models(data = dat) |>
    identify_interval(t0, t1)

  expect_equal(spec$interval$start, "t0")
  expect_equal(spec$interval$stop, "t1")
})

test_that("piping multiple identify_ calls works", {
  dat <- sim_func1(n = 50, seed = 1)
  spec <- specify_models(data = dat) |>
    identify_outcome(event, type = "time_to_event") |>
    identify_treatment(treatment, formula = ~ age + sex) |>
    identify_censoring(censored, formula = ~ age, model = "coxph") |>
    identify_competing_risk(event_type, event_value = 1) |>
    identify_subject(id)

  expect_s3_class(spec, "cr_spec")
  expect_equal(spec$outcome$name, "event")
  expect_equal(spec$treatment$name, "treatment")
  expect_length(spec$censoring, 1)
  expect_equal(spec$competing_risk$event_value, 1)
  expect_equal(spec$subject$name, "id")
})

test_that("print.cr_spec outputs without error", {
  dat <- sim_func1(n = 20, seed = 1)
  spec <- specify_models(data = dat) |>
    identify_outcome(event, type = "time_to_event") |>
    identify_treatment(treatment, formula = ~ age)

  expect_output(print(spec), "cleanTMLE Model Specification")
})
