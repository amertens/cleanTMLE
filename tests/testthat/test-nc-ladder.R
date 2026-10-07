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
  ests <- c("ATE", "trimmed_ATE", "ATT", "ATO")
  expect_identical(ok$by_rung$rung, rep(c("full cohort", "no_transfer"), each = 4))
  expect_identical(ok$by_rung$estimand, rep(ests, 2))
  expect_identical(ok$verdict, stats::setNames(rep("GO", 4), ests))
  bad <- ladder_for(c(nc_bad = "care use"), nc_criteria = list(null_band = c(-0.1, 0.1)))
  expect_identical(names(bad$verdict), ests)
  expect_identical(bad$verdict[["ATE"]], "STOP")
})

test_that("each control is estimated for every declared estimand, one row each", {
  nc <- ladder_for(c(nc_visit = "care use"), nc_criteria = list(null_band = c(-0.1, 0.1)),
                   restrictions = list(no_transfer = ~ transfer == 0))
  ests <- c("ATE", "trimmed_ATE", "ATT", "ATO")
  expect_true("estimand" %in% names(nc$table))
  expect_identical(nrow(nc$table), 2L * 4L)
  expect_identical(nc$table$estimand, rep(ests, 2))
  expect_identical(nc$table$rung, rep(c("full cohort", "no_transfer"), each = 4))
  expect_identical(nc$by_domain$estimand, rep(ests, 2))
  expect_true(all(nc$table$status == "estimated"))
  # The estimates differ by estimand: they are not the ATE repeated.
  expect_gt(length(unique(round(nc$table$estimate[nc$table$rung == "full cohort"], 8))), 1L)
  # A plan that declares one estimand gets one verdict, named by it.
  one <- ladder_for(c(nc_visit = "care use"), estimands = "ATO",
                    nc_criteria = list(null_band = c(-0.1, 0.1)))
  expect_identical(names(one$verdict), "ATO")
  expect_identical(one$table$estimand, "ATO")
})

test_that("the verdict vector is the full-cohort rung of by_rung, named by estimand", {
  d <- make_design(500, seed = 5)
  withr::with_seed(7, {
    d$design$grp <- stats::rbinom(500, 1, 0.5)
    d$design$nc_sub <- stats::rbinom(500, 1, stats::plogis(-1 + 3 * d$design$A * d$design$grp))
  })
  p <- fast_plan(negative_controls = c(nc_sub = "care use"),
                 nc_criteria = list(null_band = c(-0.1, 0.1)))
  lk <- create_analysis_lock(d$design, "A", covs, p)
  nc <- negative_control_ladder(lk, assess_design(lk))
  expect_identical(nc$by_rung$estimand, p$estimands)
  expect_identical(nc$verdict, stats::setNames(nc$by_rung$verdict, nc$by_rung$estimand))
  expect_identical(nc$verdict[["ATE"]], "STOP")
})

test_that("an inestimable control makes the domain insufficient", {
  d <- make_design(500, seed = 5)
  d$design$nc_rare <- 0L
  d$design$nc_rare[1:2] <- 1L
  p <- fast_plan(negative_controls = c(nc_rare = "rare"))
  lk <- create_analysis_lock(d$design, "A", covs, p)
  nc <- negative_control_ladder(lk, assess_design(lk))
  expect_match(nc$table$status, "inestimable")
  expect_identical(nrow(nc$table), 4L)
  expect_identical(nc$verdict, stats::setNames(rep("FLAG", 4),
                                               c("ATE", "trimmed_ATE", "ATT", "ATO")))
})

test_that("a rung with one arm is recorded, not fatal", {
  d <- make_design(500, seed = 5)
  p <- fast_plan(negative_controls = c(nc_visit = "care use"),
                 restrictions = list(treated_only = ~ A == 1))
  lk <- create_analysis_lock(d$design, "A", covs, p)
  nc <- negative_control_ladder(lk, assess_design(lk))
  expect_match(nc$table$status[nc$table$rung == "treated_only"], "inestimable|failed")
  expect_identical(sum(nc$table$rung == "treated_only"), 4L)
})

test_that(".nc_grade applies the point and ci rules and the precedence of readings", {
  tab <- function(...) {
    x <- data.frame(...)
    x$negative_control <- paste0("nc", seq_len(nrow(x)))
    x$n <- 100L
    if (is.null(x$estimand)) x$estimand <- "ATE"
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
  # grading is per rung and estimand: the same control can pass on one and fail on another
  t8 <- tab(rung = "r", estimand = c("ATE", "ATO"), domain = "d",
            estimate = c(0.5, 0.01), ci_lower = c(0.4, -0.02), ci_upper = c(0.6, 0.04),
            status = est)
  g8 <- .nc_grade(t8, crit("point"))
  expect_identical(g8$by_domain$estimand, c("ATE", "ATO"))
  expect_identical(g8$by_domain$reading, c("fail", "pass"))
  expect_identical(g8$by_rung$estimand, c("ATE", "ATO"))
  expect_identical(g8$by_rung$verdict, c("STOP", "GO"))
  # domains are separate within an estimand
  t9 <- tab(rung = "r", estimand = c("ATE", "ATE", "ATO", "ATO"), domain = c("a", "b", "a", "b"),
            estimate = c(0, 0.5, 0, 0), ci_lower = c(0, 0.4, 0, 0), ci_upper = c(0, 0.6, 0, 0),
            status = est)
  g9 <- .nc_grade(t9, crit("point"))
  expect_identical(g9$by_domain$reading, c("pass", "fail", "pass", "pass"))
  expect_identical(g9$by_rung$verdict, c("STOP", "GO"))
})

test_that("the verdict is the full cohort's for each estimand, whatever the restricted rungs read", {
  # nc_sub is driven by treatment only where grp == 1, so the full cohort stops, the
  # grp == 0 rung reads the control as null (GO), and the untreated-only rung empties the
  # treated arm (FLAG). estimate_effect() analyses the full cohort, so the gate is STOP.
  d <- make_design(500, seed = 5)
  withr::with_seed(7, {
    d$design$grp <- stats::rbinom(500, 1, 0.5)
    d$design$nc_sub <- stats::rbinom(500, 1, stats::plogis(-1 + 3 * d$design$A * d$design$grp))
  })
  p <- fast_plan(negative_controls = c(nc_sub = "care use"),
                 nc_criteria = list(null_band = c(-0.1, 0.1)),
                 restrictions = list(grp0 = ~ grp == 0, untreated_only = ~ A == 0),
                 tolerance = list(bias = 0.1, coverage = 0.5), reps = 10L)
  lk <- create_analysis_lock(d$design, "A", covs, p)
  ds <- assess_design(lk)
  nc <- negative_control_ladder(lk, ds)
  ests <- p$estimands
  expect_identical(nc$by_rung$rung, rep(c("full cohort", "grp0", "untreated_only"), each = 4))
  expect_identical(nc$by_rung$estimand, rep(ests, 3))
  rg <- function(r) nc$by_rung$verdict[nc$by_rung$rung == r]
  expect_identical(rg("full cohort"), unname(nc$verdict))
  expect_identical(rg("full cohort")[ests == "ATE"], "STOP")
  expect_identical(rg("grp0")[ests == "ATE"], "GO")
  expect_identical(rg("untreated_only")[ests == "ATE"], "FLAG")
  expect_identical(nc$verdict[["ATE"]], "STOP")
  sim <- suppressMessages(simulate_design(lk, ds))
  dr <- design_report(lk, ds, sim, nc)
  expect_identical(dr$decision$nc_verdict, nc$verdict[[dr$decision$primary]])
  expect_error(unblind(lk, dr, d$outcomes, approved_by = "Review team"),
               "negative-control verdict is STOP")
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

test_that(".nc_grade relabels an estimated row with a non-finite estimate as failed", {
  tab <- data.frame(rung = "full cohort", estimand = "ATE", negative_control = c("nc1", "nc2"),
                    domain = "d",
                    n = 100L, estimate = c(NaN, 0.01), ci_lower = c(-0.1, -0.05),
                    ci_upper = c(0.1, NA), status = "estimated", stringsAsFactors = FALSE)
  g <- .nc_grade(tab, list(null_band = c(-0.1, 0.1), rule = "point", min_per_domain = 1))
  expect_identical(g$table$status, rep("failed: non-finite estimate", 2))
  expect_identical(g$by_domain$n_estimable, 0L)
  expect_identical(g$by_domain$reading, "insufficient")
  expect_identical(g$by_rung$verdict, "FLAG")
})
