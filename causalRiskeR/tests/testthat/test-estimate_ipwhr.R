test_that("estimate_ipwhr produces hr object", {
  dat <- sim_func1(n = 300, seed = 20)
  spec <- specify_models(data = dat) |>
    identify_outcome(event, type = "time_to_event") |>
    identify_treatment(treatment, formula = ~ age + sex + biomarker)

  fit <- estimate_ipwhr(spec)

  expect_s3_class(fit, "hr")
  expect_true("hr_table" %in% names(fit))
  expect_true(nrow(fit$hr_table) >= 1)
  expect_true(all(c("term", "hr", "ci_lower", "ci_upper") %in%
                    names(fit$hr_table)))
})

test_that("hr_data extracts table", {
  dat <- sim_func1(n = 200, seed = 21)
  spec <- specify_models(data = dat) |>
    identify_outcome(event, type = "time_to_event") |>
    identify_treatment(treatment, formula = ~ age + sex)

  fit <- estimate_ipwhr(spec)
  hrd <- hr_data(fit)

  expect_true(is.data.frame(hrd))
  expect_true("hr" %in% names(hrd))
})

test_that("forest_plot returns ggplot", {
  dat <- sim_func1(n = 200, seed = 22)
  spec <- specify_models(data = dat) |>
    identify_outcome(event, type = "time_to_event") |>
    identify_treatment(treatment, formula = ~ age)

  fit <- estimate_ipwhr(spec)
  p <- forest_plot(fit)

  expect_s3_class(p, "ggplot")
})

test_that("estimate_ipwhr with covariates", {
  dat <- sim_func1(n = 200, seed = 23)
  spec <- specify_models(data = dat) |>
    identify_outcome(event, type = "time_to_event") |>
    identify_treatment(treatment, formula = ~ age + sex)

  fit <- estimate_ipwhr(spec, covariates = c("comorbidity"))

  expect_true(nrow(fit$hr_table) >= 2)
})

test_that("print.hr outputs without error", {
  dat <- sim_func1(n = 100, seed = 24)
  spec <- specify_models(data = dat) |>
    identify_outcome(event, type = "time_to_event") |>
    identify_treatment(treatment, formula = ~ age)

  fit <- estimate_ipwhr(spec)
  expect_output(print(fit), "IPW Hazard Ratio")
})

test_that("make_table2.hr produces results table", {
  dat <- sim_func1(n = 200, seed = 25)
  spec <- specify_models(data = dat) |>
    identify_outcome(event, type = "time_to_event") |>
    identify_treatment(treatment, formula = ~ age + sex)

  fit <- estimate_ipwhr(spec)
  tbl <- make_table2(fit)

  expect_true(is.data.frame(tbl))
  expect_true("group" %in% names(tbl))
  expect_true("n" %in% names(tbl))
})
