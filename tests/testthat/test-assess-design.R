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

test_that("the design is stamped and refuses a modified lock", {
  d <- make_design(300)
  lk <- create_analysis_lock(d$design, "A", c("w1", "w2", "w3"), fast_plan())
  ds <- assess_design(lk)
  expect_true(.check_stamp(lk, ds, "design"))
  lk$plan$K <- 9L
  expect_error(assess_design(lk), "hash mismatch")
})
