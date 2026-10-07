test_that("on cr_example the vignette plan moves the ladder off the ATE", {
  skip_on_cran()
  lock <- create_analysis_lock(cr_example$design, "A",
                               c("age", "female", "severity", "prior_visits"), vignette_plan())
  design <- assess_design(lock)
  sim <- suppressMessages(simulate_design(lock, design))
  st <- stats::setNames(sim$verdict$status, sim$verdict$estimand)
  expect_identical(st[["ATE"]], "infeasible")
  expect_identical(st[["trimmed_ATE"]], "feasible")
  expect_identical(st[["ATO"]], "feasible")
})
