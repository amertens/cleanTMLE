covs <- c("w1", "w2", "w3")
ladder_for <- function(nc, ...) {
  d <- make_design(500, seed = 5)
  p <- fast_plan(negative_controls = nc, ...)
  lk <- create_analysis_lock(d$design, "A", covs, p)
  negative_control_ladder(lk, assess_design(lk))
}

test_that("no declared controls gives an empty ladder with no verdict", {
  d <- make_design(300)
  lk <- create_analysis_lock(d$design, "A", covs, fast_plan())
  nc <- negative_control_ladder(lk, assess_design(lk))
  expect_false(nc$declared)
  expect_true(is.na(nc$verdict))
})

test_that("a null control passes and a treatment-driven control stops", {
  ok <- ladder_for(c(nc_visit = "care use"), nc_criteria = list(null_band = c(-0.1, 0.1)),
                   restrictions = list(no_transfer = ~ transfer == 0))
  expect_identical(ok$by_rung$rung, c("full cohort", "no_transfer"))
  expect_identical(ok$verdict, "GO")
  bad <- ladder_for(c(nc_bad = "care use"), nc_criteria = list(null_band = c(-0.1, 0.1)))
  expect_identical(bad$verdict, "STOP")
})

test_that("an inestimable control makes the domain insufficient", {
  d <- make_design(500, seed = 5)
  d$design$nc_rare <- 0L
  d$design$nc_rare[1:2] <- 1L
  p <- fast_plan(negative_controls = c(nc_rare = "rare"))
  lk <- create_analysis_lock(d$design, "A", covs, p)
  nc <- negative_control_ladder(lk, assess_design(lk))
  expect_match(nc$table$status, "inestimable")
  expect_identical(nc$verdict, "FLAG")
})

test_that("a rung with one arm is recorded, not fatal", {
  d <- make_design(500, seed = 5)
  p <- fast_plan(negative_controls = c(nc_visit = "care use"),
                 restrictions = list(treated_only = ~ A == 1))
  lk <- create_analysis_lock(d$design, "A", covs, p)
  nc <- negative_control_ladder(lk, assess_design(lk))
  expect_match(nc$table$status[nc$table$rung == "treated_only"], "inestimable|failed")
})
