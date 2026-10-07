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

test_that("on cr_example the negative controls are read on the primary estimand and do not stop it", {
  skip_on_cran()
  lock <- create_analysis_lock(cr_example$design, "A",
                               c("age", "female", "severity", "prior_visits"), vignette_plan())
  expect_identical(lock$plan$nc_criteria$rule, "interval")
  design <- assess_design(lock)
  sim <- suppressMessages(simulate_design(lock, design))
  nc <- negative_control_ladder(lock, design)
  expect_identical(names(nc$verdict), lock$plan$estimands)
  dossier <- design_report(lock, design, sim, nc)
  # The ATE is infeasible here, so its control estimates (extrapolation) must not gate.
  expect_identical(dossier$decision$primary, "trimmed_ATE")
  expect_identical(dossier$decision$nc_estimand, "trimmed_ATE")
  expect_true(is.character(nc$verdict[["ATE"]]))
  expect_identical(dossier$decision$nc_verdict, nc$verdict[["trimmed_ATE"]])
  # The interval of the control is wider than the null band, so it is too imprecise
  # to call: FLAG, which does not block unblinding.
  expect_false(identical(dossier$decision$nc_verdict, "STOP"))
  expect_identical(dossier$decision$nc_verdict, "FLAG")
  ub <- unblind(lock, dossier, cr_example$outcomes, approved_by = "Review team")
  expect_s3_class(ub, "cr_unblinded")
  expect_true(is.na(ub$approval$override))
})
