test_that("update_treatment modifies spec", {
  dat <- sim_func1(n = 50, seed = 50)
  spec <- specify_models(data = dat) |>
    identify_outcome(event, type = "time_to_event") |>
    identify_treatment(treatment, formula = ~ age)

  spec2 <- update_treatment(spec, treatment, formula = ~ age + sex)
  expect_true(grepl("sex", deparse(spec2$treatment$formula)))
})

test_that("compare_fits works", {
  dat <- sim_func1(n = 200, seed = 51)

  spec1 <- specify_models(data = dat) |>
    identify_outcome(event, type = "time_to_event") |>
    identify_treatment(treatment, formula = ~ age)

  spec2 <- specify_models(data = dat) |>
    identify_outcome(event, formula = ~ treatment + age,
                     type = "time_to_event") |>
    identify_treatment(treatment)

  fit1 <- estimate_ipwrisk(spec1, risk_time = c(24))
  fit2 <- estimate_gcomprisk(spec2, risk_time = c(24))

  comp <- compare_fits(IPW = fit1, GComp = fit2, risk_time = 24)

  expect_true(is.data.frame(comp))
  expect_true("fit_name" %in% names(comp))
  expect_true(nrow(comp) >= 2)
})

test_that("re_estimate works on cumrisk", {
  dat <- sim_func1(n = 100, seed = 52)
  spec <- specify_models(data = dat) |>
    identify_outcome(event, type = "time_to_event") |>
    identify_treatment(treatment, formula = ~ age)

  fit1 <- estimate_ipwrisk(spec, risk_time = c(12))
  fit2 <- re_estimate(fit1, risk_time = c(12, 24))

  expect_s3_class(fit2, "cumrisk")
  expect_equal(nrow(fit2$risk), 4)  # 2 groups x 2 times
})
