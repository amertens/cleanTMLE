covs <- c("w1", "w2", "w3")
# The simulation is only needed to say which candidate each control uses, so a
# few repetitions will do.
sim_of <- function(lk, ds) suppressMessages(simulate_design(lk, ds))
ladder_for <- function(nc, ...) {
  d <- make_design(500, seed = 5)
  args <- utils::modifyList(list(reps = 6L), list(...))
  p <- do.call(fast_plan, c(list(negative_controls = nc), args))
  lk <- create_analysis_lock(d$design, "A", covs, p)
  ds <- assess_design(lk)
  negative_control_ladder(lk, ds, sim_of(lk, ds))
}

test_that("no declared controls gives an empty ladder with no verdict", {
  d <- make_design(300)
  lk <- create_analysis_lock(d$design, "A", covs, fast_plan(reps = 6L))
  ds <- assess_design(lk)
  nc <- negative_control_ladder(lk, ds, sim_of(lk, ds))
  expect_false(nc$declared)
  expect_true(is.na(nc$verdict))
})

test_that("a null control passes and a treatment-driven control stops", {
  ok <- ladder_for(c(nc_visit = "care use"), nc_criteria = list(null_band = c(-0.1, 0.1), rule = "point"),
                   restrictions = list(no_transfer = ~ transfer == 0))
  ests <- c("ATE", "trimmed_ATE", "ATT", "ATO")
  expect_identical(ok$by_rung$rung, rep(c("full cohort", "no_transfer"), each = 4))
  expect_identical(ok$by_rung$estimand, rep(ests, 2))
  expect_identical(ok$verdict, stats::setNames(rep("GO", 4), ests))
  bad <- ladder_for(c(nc_bad = "care use"), nc_criteria = list(null_band = c(-0.1, 0.1), rule = "point"))
  expect_identical(names(bad$verdict), ests)
  expect_identical(bad$verdict[["ATE"]], "STOP")
})

test_that("each control is estimated for every declared estimand, one row each", {
  nc <- ladder_for(c(nc_visit = "care use"), nc_criteria = list(null_band = c(-0.1, 0.1), rule = "point"),
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
                    nc_criteria = list(null_band = c(-0.1, 0.1), rule = "point"))
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
                 nc_criteria = list(null_band = c(-0.1, 0.1), rule = "point"), reps = 6L)
  lk <- create_analysis_lock(d$design, "A", covs, p)
  ds <- assess_design(lk)
  nc <- negative_control_ladder(lk, ds, sim_of(lk, ds))
  expect_identical(nc$by_rung$estimand, p$estimands)
  expect_identical(nc$verdict, stats::setNames(nc$by_rung$verdict, nc$by_rung$estimand))
  expect_identical(nc$verdict[["ATE"]], "STOP")
})

test_that("an inestimable control makes the domain insufficient", {
  d <- make_design(500, seed = 5)
  d$design$nc_rare <- 0L
  d$design$nc_rare[1:2] <- 1L
  p <- fast_plan(negative_controls = c(nc_rare = "rare"), reps = 6L)
  lk <- create_analysis_lock(d$design, "A", covs, p)
  ds <- assess_design(lk)
  nc <- negative_control_ladder(lk, ds, sim_of(lk, ds))
  expect_match(nc$table$status, "inestimable")
  expect_identical(nrow(nc$table), 4L)
  expect_identical(nc$verdict, stats::setNames(rep("FLAG", 4),
                                               c("ATE", "trimmed_ATE", "ATT", "ATO")))
})

test_that("a rung with one arm is recorded, not fatal", {
  d <- make_design(500, seed = 5)
  p <- fast_plan(negative_controls = c(nc_visit = "care use"),
                 restrictions = list(treated_only = ~ A == 1), reps = 6L)
  lk <- create_analysis_lock(d$design, "A", covs, p)
  ds <- assess_design(lk)
  nc <- negative_control_ladder(lk, ds, sim_of(lk, ds))
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
                 nc_criteria = list(null_band = c(-0.1, 0.1), rule = "point"),
                 restrictions = list(grp0 = ~ grp == 0, untreated_only = ~ A == 0),
                 tolerance = list(bias = 0.1, coverage = 0.5), reps = 10L)
  lk <- create_analysis_lock(d$design, "A", covs, p)
  ds <- assess_design(lk)
  sim <- sim_of(lk, ds)
  nc <- negative_control_ladder(lk, ds, sim)
  ests <- p$estimands
  expect_identical(nc$by_rung$rung, rep(c("full cohort", "grp0", "untreated_only"), each = 4))
  expect_identical(nc$by_rung$estimand, rep(ests, 3))
  rg <- function(r) nc$by_rung$verdict[nc$by_rung$rung == r]
  expect_identical(rg("full cohort"), unname(nc$verdict))
  expect_identical(rg("full cohort")[ests == "ATE"], "STOP")
  expect_identical(rg("grp0")[ests == "ATE"], "GO")
  expect_identical(rg("untreated_only")[ests == "ATE"], "FLAG")
  expect_identical(nc$verdict[["ATE"]], "STOP")
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
  expect_error(negative_control_ladder(lk, assess_design(other), sim_of(lk, assess_design(lk))),
               "not computed from this lock")
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

test_that("the interval rule reads each control from its interval against the null band", {
  tab <- function(est, lo, hi, status = "estimated", domain = "d", estimand = "ATE",
                  rung = "r") {
    data.frame(rung = rung, estimand = estimand, negative_control = paste0("nc", seq_along(est)),
               domain = domain, n = 100L, estimate = est, ci_lower = lo, ci_upper = hi,
               status = status, stringsAsFactors = FALSE)
  }
  crit <- list(null_band = c(-0.1, 0.1), rule = "interval", min_per_domain = 1L)
  # inside, outside (either side), inconclusive (overlapping either edge or spanning the band)
  g <- .nc_grade(tab(est = c(0, 0.5, -0.5, 0.12, -0.12, 0),
                     lo = c(-0.05, 0.4, -0.6, 0.05, -0.2, -0.3),
                     hi = c(0.05, 0.6, -0.4, 0.2, -0.05, 0.3)), crit)
  expect_identical(g$table$reading,
                   c("inside", "outside", "outside", "inconclusive", "inconclusive", "inconclusive"))
  # an interval exactly on the band edge is inside; touching the band from outside is not outside
  e <- .nc_grade(tab(est = c(0, 0.1, -0.1), lo = c(-0.1, 0.1, -0.2), hi = c(0.1, 0.2, -0.1)), crit)
  expect_identical(e$table$reading, c("inside", "inconclusive", "inconclusive"))
  e2 <- .nc_grade(tab(est = c(0.2, -0.2), lo = c(0.1000001, -0.3), hi = c(0.3, -0.1000001)), crit)
  expect_identical(e2$table$reading, c("outside", "outside"))
  # a missing or infinite interval is a failed control, never "inside"
  nf <- .nc_grade(tab(est = c(0, 0, 0), lo = c(NA, -Inf, -0.05), hi = c(0.05, 0.05, Inf)), crit)
  expect_identical(nf$table$status, rep("failed: non-finite estimate", 3))
  expect_true(all(is.na(nf$table$reading)))
  expect_false(any(nf$table$reading %in% "inside"))
  expect_identical(nf$by_domain$reading, "insufficient")
  expect_identical(nf$by_rung$verdict, "FLAG")
  # in_band is kept: TRUE inside, FALSE outside, NA when inconclusive
  expect_identical(g$table$in_band, c(TRUE, FALSE, FALSE, NA, NA, NA))
  # domain and rung mapping: one control each
  one <- function(lo, hi, ...) .nc_grade(tab(est = (lo + hi) / 2, lo = lo, hi = hi, ...), crit)
  expect_identical(one(-0.05, 0.05)$by_domain$reading, "pass")
  expect_identical(one(-0.05, 0.05)$by_rung$verdict, "GO")
  expect_identical(one(0.2, 0.4)$by_domain$reading, "fail")
  expect_identical(one(0.2, 0.4)$by_rung$verdict, "STOP")
  expect_identical(one(-0.05, 0.3)$by_domain$reading, "inconclusive")
  expect_identical(one(-0.05, 0.3)$by_rung$verdict, "FLAG")
  # fail beats inconclusive; inconclusive beats insufficient; insufficient is FLAG
  two <- .nc_grade(tab(est = c(0.3, 0.1), lo = c(0.2, -0.05), hi = c(0.4, 0.3)), crit)
  expect_identical(two$by_domain$reading, "fail")
  three <- .nc_grade(tab(est = c(0.1, 0), lo = c(-0.05, NA), hi = c(0.3, NA),
                         status = c("estimated", "inestimable (a cell below 5 events)")),
                     utils::modifyList(crit, list(min_per_domain = 2L)))
  expect_identical(three$by_domain$reading, "inconclusive")
  four <- .nc_grade(tab(est = c(0, 0), lo = c(-0.05, NA), hi = c(0.05, NA),
                        status = c("estimated", "failed: x")),
                    utils::modifyList(crit, list(min_per_domain = 2L)))
  expect_identical(four$by_domain$reading, "insufficient")
  expect_identical(four$by_rung$verdict, "FLAG")
  # graded per estimand and per domain
  m <- .nc_grade(rbind(tab(0, -0.05, 0.05, estimand = "ATE"),
                       tab(0.3, 0.2, 0.4, estimand = "ATO"),
                       tab(0.1, -0.05, 0.3, estimand = "ATT")), crit)
  expect_identical(m$by_rung$estimand, c("ATE", "ATO", "ATT"))
  expect_identical(m$by_rung$verdict, c("GO", "STOP", "FLAG"))
})

test_that("point and ci rules map in_band to a reading and are never inconclusive", {
  est <- "estimated"
  tab <- data.frame(rung = "r", estimand = "ATE", negative_control = c("a", "b", "c"),
                    domain = "d", n = 100L, estimate = c(0.05, 0.2, 0.05),
                    ci_lower = c(-0.2, 0.1, -0.05), ci_upper = c(0.3, 0.3, 0.08),
                    status = c(est, est, "failed: x"), stringsAsFactors = FALSE)
  for (rule in c("point", "ci")) {
    g <- .nc_grade(tab, list(null_band = c(-0.1, 0.1), rule = rule, min_per_domain = 1))
    expect_false(any(g$table$reading %in% "inconclusive"))
    expect_identical(g$table$reading[2], "outside")
    expect_true(is.na(g$table$reading[3]))
    expect_identical(g$table$reading[1], if (rule == "point") "inside" else "outside")
    expect_identical(g$table$in_band[1:2], g$table$reading[1:2] == "inside")
  }
})

test_that("a ladder under the interval rule can read FLAG where the point rule reads STOP", {
  # nc_visit is null; the band is narrower than the sampling error, so a point rule
  # fails it and the interval rule calls it inconclusive.
  nc <- ladder_for(c(nc_visit = "care use"),
                   nc_criteria = list(null_band = c(-0.001, 0.001), rule = "interval"))
  expect_true("reading" %in% names(nc$table))
  expect_false(any(nc$verdict == "STOP"))
  expect_true(all(nc$verdict == "FLAG"))
  pt <- ladder_for(c(nc_visit = "care use"),
                   nc_criteria = list(null_band = c(-0.001, 0.001), rule = "point"))
  expect_true(all(pt$verdict == "STOP"))
})

# A simulation whose verdict and metrics are edited, so each estimand selects a
# visibly different candidate. fast_plan() has two candidates (truncation 0.01
# and 0.05). ATE and trimmed_ATE are feasible and select 0.05 and 0.01. ATT and
# ATO are infeasible; the smallest worst-case RMSE is 0.05 for ATT and 0.01 for
# ATO.
forged_sim <- function(fx) {
  sim <- fx$sim
  v <- sim$verdict
  set <- function(e, status, sel) {
    i <- v$estimand == e
    v$status[i] <<- status
    v$selected[i] <<- sel
    v$library[i] <<- if (is.na(sel)) NA_character_ else "glm"
    v$truncation[i] <<- if (is.na(sel)) NA_real_ else as.numeric(sub("glm_t", "", sel))
  }
  set("ATE", "feasible", "glm_t0.05")
  set("trimmed_ATE", "feasible", "glm_t0.01")
  set("ATT", "infeasible", NA_character_)
  set("ATO", "infeasible", NA_character_)
  sim$verdict <- v
  m <- sim$metrics
  m$rmse <- 0.2
  m$rmse[m$estimand == "ATT" & m$candidate == "glm_t0.05"] <- 0.01
  m$rmse[m$estimand == "ATO" & m$candidate == "glm_t0.01"] <- 0.01
  sim$metrics <- m
  sim$stamp <- .stamp(fx$lock, sim)
  sim
}

nc_fx <- local({
  cache <- NULL
  function() {
    if (is.null(cache)) cache <<- fixture_pipeline()
    cache
  }
})

test_that("each estimand's controls use the candidate the simulation selected for it", {
  fx <- nc_fx()
  sim <- forged_sim(fx)
  nc <- negative_control_ladder(fx$lock, fx$design, sim)
  tb <- nc$table
  expect_true(all(c("candidate", "candidate_feasible") %in% names(tb)))
  pick <- function(col) vapply(c("ATE", "trimmed_ATE", "ATT", "ATO"), function(e) {
    x <- unique(tb[[col]][tb$estimand == e])
    expect_length(x, 1L)
    x
  }, FUN.VALUE = tb[[col]][1])
  expect_identical(pick("candidate"), c(ATE = "glm_t0.05", trimmed_ATE = "glm_t0.01",
                                        ATT = "glm_t0.05", ATO = "glm_t0.01"))
  expect_identical(pick("candidate_feasible"),
                   c(ATE = TRUE, trimmed_ATE = TRUE, ATT = FALSE, ATO = FALSE))
  expect_identical(nc$simulation_hash, sim$stamp$hash)
  # The grading columns are unchanged.
  expect_true(all(c("rung", "estimand", "negative_control", "domain", "n", "estimate",
                    "ci_lower", "ci_upper", "status", "reading", "in_band") %in% names(tb)))
  expect_identical(names(nc$by_rung), c("rung", "estimand", "verdict"))
  expect_identical(names(nc$by_domain),
                   c("rung", "estimand", "domain", "n_estimable", "reading"))
})

test_that("fit_candidate receives the whole selected candidate, including a non-glm library", {
  fx <- nc_fx()
  # Two libraries, so the learners differ between candidates. The simulation is
  # the one-library fixture simulation re-stamped for this plan, with alt rows
  # added to its metrics; the ladder only reads its verdict and metrics.
  d <- make_design(400, seed = 2, strength = 0.8)
  p2 <- fast_plan(negative_controls = c(nc_visit = "care use"), reps = 6L,
                  candidates = list(library = list(glm = "glm", alt = c("glm", "glmnet")),
                                    truncation = c(0.01, 0.05)))
  lk2 <- create_analysis_lock(d$design, "A", covs, p2)
  ds2 <- assess_design(lk2)
  sim <- fx$sim
  sim$design_hash <- ds2$stamp$hash
  m <- sim$metrics
  alt <- m
  alt$library <- "alt"
  alt$candidate <- sub("^glm", "alt", alt$candidate)
  m <- rbind(m, alt)
  m$rmse <- 0.2
  m$rmse[m$estimand == "ATT" & m$candidate == "alt_t0.01"] <- 0.01
  m$rmse[m$estimand == "ATO" & m$candidate == "glm_t0.05"] <- 0.01
  sim$metrics <- m
  v <- sim$verdict
  pick <- function(e, status, lib, tr) {
    i <- v$estimand == e
    v$status[i] <<- status
    v$selected[i] <<- if (is.na(lib)) NA_character_ else paste0(lib, "_t", tr)
    v$library[i] <<- lib
    v$truncation[i] <<- tr
  }
  pick("ATE", "feasible", "alt", 0.05)
  pick("trimmed_ATE", "feasible", "glm", 0.01)
  pick("ATT", "infeasible", NA_character_, NA_real_)
  pick("ATO", "infeasible", NA_character_, NA_real_)
  sim$verdict <- v
  sim$stamp <- .stamp(lk2, sim)

  seen <- list()
  local_mocked_bindings(fit_candidate = function(X, A, Y, estimand, candidate, plan, folds,
                                                 seed, g_fit = NULL) {
    seen[[length(seen) + 1L]] <<- list(estimand = estimand, candidate = candidate)
    data.frame(estimand = estimand, truncation = candidate$truncation, estimate = 0,
               se = 0.01, ci_lower = -0.02, ci_upper = 0.02, n_population = length(Y),
               failed = FALSE, message = NA_character_, stringsAsFactors = FALSE)
  })
  nc <- negative_control_ladder(lk2, ds2, sim)
  cand <- function(id, lib, learners, tr)
    list(id = id, library = lib, learners = learners, truncation = tr)
  expected <- list(
    ATE = cand("alt_t0.05", "alt", c("glm", "glmnet"), 0.05),
    trimmed_ATE = cand("glm_t0.01", "glm", "glm", 0.01),
    ATT = cand("alt_t0.01", "alt", c("glm", "glmnet"), 0.01),
    ATO = cand("glm_t0.05", "glm", "glm", 0.05))
  expect_length(seen, 4L)
  for (s in seen) {
    expect_length(s$estimand, 1L)
    expect_identical(s$candidate, expected[[s$estimand]])
  }
  expect_identical(unique(nc$table[c("estimand", "candidate")])$candidate,
                   unname(vapply(expected, function(x) x$id, "")))
})

test_that(".nc_candidates takes the learners from the selected or best library", {
  p <- analysis_plan(outcome = "y", estimands = c("ATE", "ATT", "ATO"),
                     ps_library = "glm",
                     candidates = list(library = list(glm = "glm", pen = c("glmnet", "glm")),
                                       truncation = c(0.01, 0.05)),
                     negative_controls = c(nc1 = "d"),
                     nc_criteria = list(null_band = c(-0.1, 0.1), rule = "point"))
  cands <- vapply(.candidates(p), function(x) x$id, "")
  # One row per candidate and surface; `a` and `b` are the RMSE on two surfaces,
  # in the order glm_t0.01, glm_t0.05, pen_t0.01, pen_t0.05.
  mk <- function(e, a, b) do.call(rbind, lapply(1:2, function(s) data.frame(
    estimand = e, candidate = cands, surface = s,
    library = vapply(.candidates(p), function(x) x$library, ""),
    truncation = vapply(.candidates(p), function(x) x$truncation, 0),
    pass = FALSE, unresolved = FALSE, pass_hw = TRUE, rmse = if (s == 1) a else b,
    stringsAsFactors = FALSE)))
  # ATT: by worst case pen_t0.01 (0.2, 0.2) wins; by mean pen_t0.05 (0.05, 0.40)
  # would. glm_t0.05 has an NA cell and the smallest finite cell, and is skipped.
  metrics <- rbind(mk("ATE", rep(0.3, 4), rep(0.3, 4)),
                   mk("ATT", c(0.3, NA, 0.2, 0.05), c(0.3, 0.01, 0.2, 0.40)),
                   mk("ATO", rep(NA_real_, 4), rep(NA_real_, 4)))
  verdict <- data.frame(estimand = c("ATE", "ATT", "ATO"),
                        status = c("feasible", "infeasible", "infeasible"),
                        selected = c("pen_t0.05", NA, NA), library = c("pen", NA, NA),
                        truncation = c(0.05, NA, NA), stringsAsFactors = FALSE)
  got <- .nc_candidates(list(metrics = metrics, verdict = verdict), p)
  expect_identical(got$ATE$id, "pen_t0.05")
  expect_identical(got$ATE$learners, c("glmnet", "glm"))
  expect_identical(got$ATE$truncation, 0.05)
  expect_true(got$ATE$feasible)
  expect_false(got$ATT$feasible)
  expect_identical(got$ATT$id, "pen_t0.01")
  expect_identical(got$ATT$learners, c("glmnet", "glm"))
  expect_match(got$ATT$basis, "smallest worst-case RMSE")
  # No usable RMSE: the first candidate, and the row says so.
  expect_false(got$ATO$feasible)
  expect_identical(got$ATO$id, .candidates(p)[[1]]$id)
  expect_match(got$ATO$basis, "first candidate")
})

test_that("a simulation from another design or lock is refused", {
  fx <- nc_fx()
  other_design <- fx$sim
  other_design$design_hash <- "not-this-design"
  other_design$stamp <- .stamp(fx$lock, other_design)
  expect_error(negative_control_ladder(fx$lock, fx$design, other_design),
               "`simulation` was not run on this `design`.", fixed = TRUE)
  d <- make_design(400, seed = 2, strength = 0.8)
  lk2 <- create_analysis_lock(d$design, "A", c("w1", "w2", "w3"),
                              fast_plan(negative_controls = c(nc_visit = "care use"),
                                        seed = 99L, reps = 6L))
  sim2 <- sim_of(lk2, assess_design(lk2))
  expect_error(negative_control_ladder(fx$lock, fx$design, sim2),
               "not computed from this lock")
  tampered <- fx$sim
  tampered$reps <- tampered$reps + 1L
  expect_error(negative_control_ladder(fx$lock, fx$design, tampered), "modified")
})

test_that("design_report refuses controls run on a different simulation or by an earlier version", {
  fx <- nc_fx()
  sim2 <- forged_sim(fx)
  nc2 <- negative_control_ladder(fx$lock, fx$design, sim2)
  expect_error(design_report(fx$lock, fx$design, fx$sim, nc2),
               "`nc` was not run on this `simulation`.", fixed = TRUE)
  expect_s3_class(design_report(fx$lock, fx$design, sim2, nc2), "cr_dossier")
  old <- fx$nc
  old$table$candidate <- NULL
  old$stamp <- .stamp(fx$lock, old)
  msg <- "This negative-control ladder was made by an earlier cleanTMLE; rerun negative_control_ladder()."
  expect_error(design_report(fx$lock, fx$design, fx$sim, old), msg, fixed = TRUE)
  # A ladder with no simulation hash is refused whether or not it declares controls.
  old3 <- fx$nc
  old3$simulation_hash <- NULL
  old3$stamp <- .stamp(fx$lock, old3)
  expect_error(design_report(fx$lock, fx$design, fx$sim, old3), msg, fixed = TRUE)
  nn <- fixture_pipeline(nc = NULL)
  expect_false(nn$nc$declared)
  expect_s3_class(nn$dossier, "cr_dossier")
  old4 <- nn$nc
  old4$simulation_hash <- NULL
  old4$stamp <- .stamp(nn$lock, old4)
  expect_error(design_report(nn$lock, nn$design, nn$sim, old4), msg, fixed = TRUE)
})

test_that("negative_control_ladder refuses a simulation made by an earlier cleanTMLE", {
  fx <- nc_fx()
  old <- fx$sim
  old$expected_events <- NULL
  old$stamp <- .stamp(fx$lock, old)
  expect_error(negative_control_ladder(fx$lock, fx$design, old),
               "This simulation was made by an earlier cleanTMLE; rerun simulate_design().",
               fixed = TRUE)
})

# A design with a factor covariate, and two controls, so exclusions can be
# checked on the expanded design-matrix columns.
nc_excl_fx <- function(exclude) {
  d <- make_design(400, seed = 3)
  d$design$site <- factor(withr::with_seed(9, sample(c("A", "B", "C"), 400, TRUE)))
  cv <- c("w1", "w2", "w3", "site")
  p <- fast_plan(negative_controls = c(nc_visit = "care use", nc_bad = "other"),
                 nc_criteria = list(null_band = c(-0.1, 0.1), rule = "point"),
                 nc_exclude = exclude, reps = 6L)
  lk <- create_analysis_lock(d$design, "A", cv, p)
  ds <- assess_design(lk)
  list(d = d, lock = lk, design = ds, sim = sim_of(lk, ds))
}

mock_record <- function(env, d) {
  function(X, A, Y, estimand, candidate, plan, folds, seed, g_fit = NULL) {
    ctl <- if (identical(as.numeric(Y), as.numeric(d$design$nc_visit))) "nc_visit" else "nc_bad"
    env$seen[[length(env$seen) + 1L]] <- list(control = ctl, cols = colnames(X), estimand = estimand)
    data.frame(estimand = estimand, truncation = candidate$truncation, estimate = 0,
               se = 0.01, ci_lower = -0.02, ci_upper = 0.02, n_population = length(Y),
               failed = FALSE, message = NA_character_, stringsAsFactors = FALSE)
  }
}

test_that("each control is fitted without the covariates it excludes, and only that control", {
  fx <- nc_excl_fx(list(nc_visit = c("site", "w1")))
  env <- new.env(); env$seen <- list()
  local_mocked_bindings(fit_candidate = mock_record(env, fx$d))
  nc <- negative_control_ladder(fx$lock, fx$design, fx$sim)
  cols <- function(ctl) unique(lapply(Filter(function(s) s$control == ctl, env$seen),
                                      `[[`, "cols"))
  expect_identical(cols("nc_visit"), list(c("w2", "w3")))
  expect_identical(cols("nc_bad"), list(c("w1", "w2", "w3", "siteB", "siteC")))
  tb <- nc$table
  expect_identical(unique(tb$excluded_covariates[tb$negative_control == "nc_visit"]), "w1, site")
  expect_identical(unique(tb$excluded_covariates[tb$negative_control == "nc_bad"]), "")
  expect_false(anyNA(tb$excluded_covariates))
})

test_that("excluding a factor covariate removes all of its columns", {
  fx <- nc_excl_fx(list(nc_bad = "site"))
  env <- new.env(); env$seen <- list()
  local_mocked_bindings(fit_candidate = mock_record(env, fx$d))
  negative_control_ladder(fx$lock, fx$design, fx$sim)
  bad <- unique(lapply(Filter(function(s) s$control == "nc_bad", env$seen), `[[`, "cols"))
  vis <- unique(lapply(Filter(function(s) s$control == "nc_visit", env$seen), `[[`, "cols"))
  expect_identical(bad, list(c("w1", "w2", "w3")))
  expect_identical(vis, list(c("w1", "w2", "w3", "siteB", "siteC")))
})

test_that("without nc_exclude every control uses every covariate and the column is empty", {
  fx <- nc_excl_fx(NULL)
  env <- new.env(); env$seen <- list()
  local_mocked_bindings(fit_candidate = mock_record(env, fx$d))
  nc <- negative_control_ladder(fx$lock, fx$design, fx$sim)
  expect_identical(unique(lapply(env$seen, `[[`, "cols")),
                   list(c("w1", "w2", "w3", "siteB", "siteC")))
  expect_identical(unique(nc$table$excluded_covariates), "")
})

test_that("controls with the same exclusions share one propensity-score fit per rung", {
  count_g <- function(exclude) {
    fx <- nc_excl_fx(exclude)
    env <- new.env(); env$seen <- list(); env$g <- character()
    real_fit_g <- .fit_g
    local_mocked_bindings(
      fit_candidate = mock_record(env, fx$d),
      .fit_g = function(X, A, ...) { env$g <- c(env$g, paste(colnames(X), collapse = ",")); real_fit_g(X, A, ...) })
    negative_control_ladder(fx$lock, fx$design, fx$sim)
    env$g
  }
  expect_length(count_g(NULL), 1L)
  expect_length(count_g(list(nc_visit = "w1", nc_bad = "w1")), 1L)
  g <- count_g(list(nc_visit = "w1"))
  expect_length(g, 2L)
  expect_false(identical(g[1], g[2]))
})

test_that("a ladder with exclusions and the real fitter returns finite rows", {
  fx <- nc_excl_fx(list(nc_visit = "w1"))
  nc <- negative_control_ladder(fx$lock, fx$design, fx$sim)
  expect_true(all(nc$table$status == "estimated"))
  expect_true(all(is.finite(nc$table$estimate)))
})
