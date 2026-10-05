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

test_that(".nc_grade applies the point and ci rules and the precedence of readings", {
  tab <- function(...) {
    x <- data.frame(...)
    x$negative_control <- paste0("nc", seq_len(nrow(x)))
    x$n <- 100L
    x
  }
  crit <- function(rule, min_per_domain = 1) {
    list(null_band = c(-0.1, 0.1), rule = rule, min_per_domain = min_per_domain)
  }
  est <- "estimated"
  # point rule: judged on the estimate, so an interval that leaves the band still passes
  t1 <- tab(rung = "r", domain = "d", estimate = 0.05, ci_lower = -0.2, ci_upper = 0.3,
            status = est)
  expect_identical(.nc_grade(t1, crit("point"))$by_domain$reading, "pass")
  # ci rule: the whole interval must sit in the band
  expect_identical(.nc_grade(t1, crit("ci"))$by_domain$reading, "fail")
  t2 <- tab(rung = "r", domain = "d", estimate = 0.05, ci_lower = -0.05, ci_upper = 0.08,
            status = est)
  expect_identical(.nc_grade(t2, crit("ci"))$by_domain$reading, "pass")
  # point rule fail
  t3 <- tab(rung = "r", domain = "d", estimate = 0.2, ci_lower = 0.1, ci_upper = 0.3,
            status = est)
  expect_identical(.nc_grade(t3, crit("point"))$by_domain$reading, "fail")
  # insufficient: too few estimable controls, none out of band
  t4 <- tab(rung = "r", domain = "d", estimate = c(0.01, NA), ci_lower = c(-0.1, NA),
            ci_upper = c(0.1, NA), status = c(est, "inestimable (a cell below 5 events)"))
  g4 <- .nc_grade(t4, crit("point", min_per_domain = 2))
  expect_identical(g4$by_domain$reading, "insufficient")
  expect_identical(g4$by_domain$n_estimable, 1L)
  expect_identical(g4$by_rung$verdict, "FLAG")
  # fail takes precedence over insufficient: one estimable, out of band, two required
  t5 <- tab(rung = "r", domain = "d", estimate = c(0.4, NA), ci_lower = c(0.3, NA),
            ci_upper = c(0.5, NA), status = c(est, "failed: x"))
  g5 <- .nc_grade(t5, crit("point", min_per_domain = 2))
  expect_identical(g5$by_domain$reading, "fail")
  expect_identical(g5$by_rung$verdict, "STOP")
  # rung verdict: STOP > FLAG > GO, over domains within a rung
  t6 <- tab(rung = c("r1", "r2", "r2", "r3", "r3", "r3"),
            domain = c("a", "a", "b", "a", "b", "c"),
            estimate = c(0, 0, 0.5, 0, NA, 0.5),
            ci_lower = c(0, 0, 0.4, 0, NA, 0.4), ci_upper = c(0, 0, 0.6, 0, NA, 0.6),
            status = c(est, est, est, est, "failed: x", est))
  g6 <- .nc_grade(t6, crit("point", min_per_domain = 1))
  expect_identical(g6$by_rung$rung, c("r1", "r2", "r3"))
  expect_identical(g6$by_rung$verdict, c("GO", "STOP", "STOP"))
  t7 <- tab(rung = c("r1", "r2", "r2"), domain = c("a", "a", "b"),
            estimate = c(0, 0, NA), ci_lower = c(0, 0, NA), ci_upper = c(0, 0, NA),
            status = c(est, est, "failed: x"))
  g7 <- .nc_grade(t7, crit("point", min_per_domain = 1))
  expect_identical(g7$by_rung$verdict, c("GO", "FLAG"))
})

test_that("the verdict is the last rung's, which can differ from the first rung's", {
  # nc_bad is driven by treatment, so the full cohort stops; the second rung keeps only the
  # untreated, so the treated arm is empty and the control cannot be estimated (FLAG).
  nc <- ladder_for(c(nc_bad = "care use"), nc_criteria = list(null_band = c(-0.1, 0.1)),
                   restrictions = list(untreated_only = ~ A == 0))
  expect_identical(nc$by_rung$rung, c("full cohort", "untreated_only"))
  expect_identical(nc$by_rung$verdict, c("STOP", "FLAG"))
  expect_identical(nc$verdict, "FLAG")
  expect_identical(nc$verdict, nc$by_rung$verdict[nrow(nc$by_rung)])
})

test_that("a design stamped by a different lock is rejected", {
  d <- make_design(300, seed = 5)
  p <- fast_plan(negative_controls = c(nc_visit = "care use"))
  lk <- create_analysis_lock(d$design, "A", covs, p)
  other <- create_analysis_lock(d$design, "A", covs,
                                fast_plan(negative_controls = c(nc_visit = "care use"),
                                          seed = 99L))
  expect_error(negative_control_ladder(lk, assess_design(other)), "not computed from this lock")
})
