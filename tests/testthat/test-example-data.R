# The design, the simulation and the ladder on cr_example under the vignette plan,
# computed once and shared by the tests in this file (the simulation is the slow part).
cr_vignette <- local({
  cache <- NULL
  function() {
    if (is.null(cache)) {
      lock <- create_analysis_lock(cr_example$design, "A",
                                   c("age", "female", "severity", "prior_visits"),
                                   vignette_plan())
      design <- assess_design(lock)
      sim <- suppressMessages(simulate_design(lock, design))
      cache <<- list(lock = lock, design = design, sim = sim,
                     nc = negative_control_ladder(lock, design))
    }
    cache
  }
})

test_that("on cr_example the vignette plan moves the ladder off the ATE", {
  skip_on_cran()
  sim <- cr_vignette()$sim
  st <- stats::setNames(sim$verdict$status, sim$verdict$estimand)
  expect_identical(st[["ATE"]], "infeasible")
  expect_identical(st[["trimmed_ATE"]], "feasible")
  expect_identical(st[["ATO"]], "feasible")
})

test_that("on cr_example the negative controls are read on the primary estimand and do not stop it", {
  skip_on_cran()
  fx <- cr_vignette()
  lock <- fx$lock; design <- fx$design; sim <- fx$sim; nc <- fx$nc
  expect_identical(lock$plan$nc_criteria$rule, "interval")
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
