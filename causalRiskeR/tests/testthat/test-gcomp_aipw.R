test_that("estimate_gcomprisk works", {
  dat <- sim_func1(n = 200, seed = 40)
  spec <- specify_models(data = dat) |>
    identify_outcome(event, formula = ~ treatment + age + sex,
                     type = "time_to_event") |>
    identify_treatment(treatment)

  fit <- estimate_gcomprisk(spec, risk_time = c(12, 24))

  expect_s3_class(fit, "gcomp")
  expect_equal(nrow(fit$risk), 4)  # 2 groups x 2 times
  expect_true(all(fit$risk$risk >= 0))
})

test_that("estimate_gcomprisk print and plot", {
  dat <- sim_func1(n = 100, seed = 41)
  spec <- specify_models(data = dat) |>
    identify_outcome(event, formula = ~ treatment + age,
                     type = "time_to_event") |>
    identify_treatment(treatment)

  fit <- estimate_gcomprisk(spec, risk_time = c(12, 24))

  expect_output(print(fit), "G-Computation")
  p <- plot(fit)
  expect_s3_class(p, "ggplot")
})

test_that("estimate_aipwrisk works", {
  dat <- sim_func1(n = 200, seed = 42)
  spec <- specify_models(data = dat) |>
    identify_outcome(event, formula = ~ treatment + age + sex,
                     type = "time_to_event") |>
    identify_treatment(treatment, formula = ~ age + sex)

  fit <- estimate_aipwrisk(spec, risk_time = c(12, 24))

  expect_s3_class(fit, "aipw")
  expect_equal(nrow(fit$risk), 4)
  expect_true(all(fit$risk$risk >= 0 & fit$risk$risk <= 1))
})

test_that("estimate_aipwrisk requires treatment", {
  dat <- sim_func1(n = 50, seed = 43)
  spec <- specify_models(data = dat) |>
    identify_outcome(event, type = "time_to_event")

  expect_error(estimate_aipwrisk(spec), "Treatment must be specified")
})

test_that("estimate_aipwrisk print and plot", {
  dat <- sim_func1(n = 100, seed = 44)
  spec <- specify_models(data = dat) |>
    identify_outcome(event, formula = ~ treatment + age,
                     type = "time_to_event") |>
    identify_treatment(treatment, formula = ~ age)

  fit <- estimate_aipwrisk(spec, risk_time = c(12))

  expect_output(print(fit), "AIPW")
  p <- plot(fit)
  expect_s3_class(p, "ggplot")
})
