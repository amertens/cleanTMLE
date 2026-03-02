test_that("estimate_ipwrisk works with basic spec (no treatment)", {
  dat <- sim_func1(n = 200, seed = 1)
  spec <- specify_models(data = dat) |>
    identify_outcome(event, type = "time_to_event")

  fit <- estimate_ipwrisk(spec, risk_time = c(6, 12, 24))

  expect_s3_class(fit, "cumrisk")
  expect_true("risk" %in% names(fit))
  expect_equal(nrow(fit$risk), 3)
  expect_true(all(fit$risk$risk >= 0 & fit$risk$risk <= 1))
  expect_equal(fit$risk$group, rep("overall", 3))
})

test_that("estimate_ipwrisk works with treatment (no PS model)", {
 dat <- sim_func1(n = 200, seed = 2)
  spec <- specify_models(data = dat) |>
    identify_outcome(event, type = "time_to_event") |>
    identify_treatment(treatment)

  fit <- estimate_ipwrisk(spec, risk_time = c(12, 24))

  expect_s3_class(fit, "cumrisk")
  groups <- unique(fit$risk$group)
  expect_length(groups, 2)
})

test_that("estimate_ipwrisk works with IPTW", {
  dat <- sim_func1(n = 300, seed = 3)
  spec <- specify_models(data = dat) |>
    identify_outcome(event, type = "time_to_event") |>
    identify_treatment(treatment, formula = ~ age + sex + biomarker)

  fit <- estimate_ipwrisk(spec, risk_time = c(12, 24))

  expect_s3_class(fit, "cumrisk")
  # Weights should not all be 1
  expect_false(all(fit$weights$iptw == 1))
})

test_that("estimate_ipwrisk with IPCW", {
  dat <- sim_func1(n = 200, seed = 4)
  spec <- specify_models(data = dat) |>
    identify_outcome(event, type = "time_to_event") |>
    identify_censoring(censored, formula = ~ age, model = "coxph")

  fit <- estimate_ipwrisk(spec, risk_time = c(12, 24))

  expect_s3_class(fit, "cumrisk")
  expect_false(all(fit$weights$ipcw == 1))
})

test_that("estimate_ipwrisk with bootstrap CIs", {
  dat <- sim_func1(n = 100, seed = 5)
  spec <- specify_models(data = dat) |>
    identify_outcome(event, type = "time_to_event") |>
    identify_treatment(treatment, formula = ~ age)

  fit <- estimate_ipwrisk(spec, risk_time = c(12), nboot = 10, seed = 99)

  expect_s3_class(fit, "cumrisk")
  expect_true("ci_lower" %in% names(fit$risk))
  expect_true("ci_upper" %in% names(fit$risk))
})

test_that("estimate_ipwrisk with trimming", {
  dat <- sim_func1(n = 200, seed = 6)
  spec <- specify_models(data = dat) |>
    identify_outcome(event, type = "time_to_event") |>
    identify_treatment(treatment, formula = ~ age + sex + biomarker)

  fit <- estimate_ipwrisk(spec, risk_time = c(12), trim = 0.05)

  expect_s3_class(fit, "cumrisk")
  # Some weights should be zero (trimmed)
  expect_true(any(fit$weights$iptw == 0))
})

test_that("estimate_ipwrisk with truncation", {
  dat <- sim_func1(n = 200, seed = 7)
  spec <- specify_models(data = dat) |>
    identify_outcome(event, type = "time_to_event") |>
    identify_treatment(treatment, formula = ~ age + sex + biomarker)

  fit <- estimate_ipwrisk(spec, risk_time = c(12), trunc = 0.99)

  expect_s3_class(fit, "cumrisk")
})

test_that("estimate_ipwrisk with SMR weights", {
  dat <- sim_func1(n = 200, seed = 8)
  spec <- specify_models(data = dat) |>
    identify_outcome(event, type = "time_to_event") |>
    identify_treatment(treatment, formula = ~ age + sex)

  fit <- estimate_ipwrisk(spec, risk_time = c(12), weight_type = "smr")

  expect_s3_class(fit, "cumrisk")
  expect_equal(fit$weight_type, "smr")
})

test_that("estimate_ipwrisk with competing risks", {
  dat <- sim_func1(n = 300, seed = 9)
  spec <- specify_models(data = dat) |>
    identify_outcome(event, type = "time_to_event") |>
    identify_treatment(treatment, formula = ~ age + sex) |>
    identify_competing_risk(event_type, event_value = 1)

  fit <- estimate_ipwrisk(spec, risk_time = c(12, 24))

  expect_s3_class(fit, "cumrisk")
})

test_that("print.cumrisk outputs without error", {
  dat <- sim_func1(n = 50, seed = 10)
  spec <- specify_models(data = dat) |>
    identify_outcome(event, type = "time_to_event")

  fit <- estimate_ipwrisk(spec, risk_time = c(12))
  expect_output(print(fit), "IPW Cumulative Risk")
})

test_that("plot.cumrisk returns ggplot", {
  dat <- sim_func1(n = 100, seed = 11)
  spec <- specify_models(data = dat) |>
    identify_outcome(event, type = "time_to_event") |>
    identify_treatment(treatment, formula = ~ age)

  fit <- estimate_ipwrisk(spec, risk_time = c(6, 12, 18, 24))

  p <- plot(fit)
  expect_s3_class(p, "ggplot")

  p_rd <- plot(fit, effect = "RD")
  expect_s3_class(p_rd, "ggplot")
})
