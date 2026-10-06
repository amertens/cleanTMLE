test_that("assess_design grades overlap and reports per-estimand weights", {
  good <- make_design(400, strength = 0.3)
  poor <- make_design(400, strength = 3)
  covs <- c("w1", "w2", "w3")
  dg <- assess_design(create_analysis_lock(good$design, "A", covs, fast_plan()))
  dp <- assess_design(create_analysis_lock(poor$design, "A", covs, fast_plan()))
  expect_s3_class(dg, "cr_design")
  expect_length(dg$g, 400L)
  expect_true(all(dg$g > 0 & dg$g < 1))
  expect_identical(dg$overlap_grade, "good")
  expect_identical(dp$overlap_grade, "poor")
  expect_setequal(dp$ess$estimand, c("ATE", "trimmed_ATE", "ATT", "ATO"))
  ess <- function(d, e) d$ess$ess_control[d$ess$estimand == e]
  expect_lt(ess(dp, "ATE"), ess(dp, "ATO"))
  expect_true(all(is.finite(dg$balance$max_abs_smd)))
  expect_lt(dg$balance$max_abs_smd[dg$balance$estimand == "ATO"], 0.1)
})

test_that("the world propensity score is a function of the covariates alone", {
  d <- make_design(300)
  p <- fast_plan(K = 2L)
  folds <- .make_folds(300, p$K, p$seed)
  j <- which(folds != folds[1])[1]
  d$design[j, c("w1", "w2", "w3")] <- d$design[1, c("w1", "w2", "w3")]
  ds <- assess_design(create_analysis_lock(d$design, "A", c("w1", "w2", "w3"), p))
  # The cross-fitted score depends on the fold, which is not a covariate ...
  expect_false(isTRUE(all.equal(ds$g[1], ds$g[j])))
  # ... but the single fit that drives the simulated world does not.
  expect_length(ds$g_world, 300L)
  expect_identical(ds$g_world[1], ds$g_world[j])
  expect_true(all(ds$g_world > 0 & ds$g_world < 1))
})

test_that("the design is stamped and refuses a modified lock", {
  d <- make_design(300)
  lk <- create_analysis_lock(d$design, "A", c("w1", "w2", "w3"), fast_plan())
  ds <- assess_design(lk)
  expect_true(.check_stamp(lk, ds, "design"))
  lk$plan$K <- 9L
  expect_error(assess_design(lk), "hash mismatch")
})

test_that("a balance failure warns and returns NA instead of passing silently", {
  d <- make_design(300)
  lk <- create_analysis_lock(d$design, "A", c("w1", "w2", "w3"), fast_plan())
  testthat::local_mocked_bindings(bal.tab = function(...) stop("boom"), .package = "cobalt")
  warns <- character()
  ds <- withCallingHandlers(assess_design(lk), warning = function(w) {
    warns <<- c(warns, conditionMessage(w))
    invokeRestart("muffleWarning")
  })
  expect_length(warns, 4L)
  expect_true(all(grepl("could not be computed", warns)))
  expect_match(warns[1], "boom")
  expect_true(all(is.na(ds$balance$max_abs_smd)))
  expect_setequal(ds$balance$estimand, c("ATE", "trimmed_ATE", "ATT", "ATO"))
})
