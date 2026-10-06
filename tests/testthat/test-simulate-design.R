covs <- c("w1", "w2", "w3")

test_that("surfaces hit the baseline risk and the true values match a Monte Carlo draw", {
  d <- make_design(500, strength = 1)
  X <- .design_matrix(d$design, covs)
  g <- stats::plogis(-0.2 + d$design$w1 + 0.3 * d$design$w2)
  p <- fast_plan(surfaces = list(baseline_risk = 0.2, heterogeneity = 1))
  s <- .surfaces(p, X, g, p$seed)
  expect_length(s, 4L)
  for (sf in s) expect_equal(mean(sf$q0), 0.2, tolerance = 1e-6)
  sf <- s[[4]]
  tv <- .true_values(sf$q1, sf$q0, g, c(0.05, 0.95))
  withr::local_seed(1)
  idx <- sample.int(500, 1e6, replace = TRUE)
  y1 <- stats::rbinom(1e6, 1, sf$q1[idx]); y0 <- stats::rbinom(1e6, 1, sf$q0[idx])
  a <- stats::rbinom(1e6, 1, g[idx])
  h <- g[idx] * (1 - g[idx])
  expect_equal(unname(tv["ATE"]), mean(y1 - y0), tolerance = 0.005)
  expect_equal(unname(tv["ATT"]), mean((y1 - y0)[a == 1]), tolerance = 0.005)
  expect_lt(abs(unname(tv["ATO"]) - sum(h * (y1 - y0)) / sum(h)), 0.005)
})

# One metrics cell: R ok repetitions with estimates truth + bias +/- spread,
# the first `covered` of them with an interval that covers the truth.
sim_cell <- function(library, R, bias = 0, spread = 0.01, covered = R, failed = 0L,
                     truth = 0.10, estimand = "ATE") {
  n <- R + failed
  est <- truth + bias + rep(c(-spread, spread), length.out = n)
  cov <- seq_len(n) <= covered
  data.frame(rep = seq_len(n), surface = "s", library = library, estimand = estimand,
             truncation = 0.01, failed = seq_len(n) > R, estimate = est,
             ci_lower = ifelse(cov, truth - 0.1, truth + 0.5),
             ci_upper = ifelse(cov, truth + 0.1, truth + 0.6), stringsAsFactors = FALSE)
}

wilson <- function(x, n, z = stats::qnorm(0.975)) {
  p <- x / n
  mid <- (p + z^2 / (2 * n)) / (1 + z^2 / n)
  half <- z / (1 + z^2 / n) * sqrt(p * (1 - p) / n + z^2 / (4 * n^2))
  c(mid - half, mid + half)
}

test_that("metrics decide on point estimates and flag unresolved cells", {
  truths <- data.frame(surface = "s", estimand = "ATE", truth = 0.10)
  res <- rbind(
    sim_cell("clear", 40),                         # bias 0, coverage 1: clearly passes
    sim_cell("cov90", 200, covered = 180),         # coverage exactly at the tolerance
    sim_cell("cov89", 200, covered = 178),         # just below: fails, but unresolved
    sim_cell("biasedge", 40, bias = 0.02),         # bias interval straddles 0.021
    sim_cell("far", 40, bias = 0.2),               # bias far outside: fails, resolved
    sim_cell("fails", 36, bias = 0.02, failed = 4L))  # 10% failed fits
  p <- fast_plan(estimands = "ATE", tolerance = list(bias = 0.021, coverage = 0.9))
  m <- .sim_metrics(res, truths, p)
  m <- m[match(c("clear", "cov90", "cov89", "biasedge", "far", "fails"), m$library), ]
  expect_false("borderline" %in% names(m))
  expect_identical(m$pass, c(TRUE, TRUE, FALSE, TRUE, FALSE, FALSE))
  expect_identical(m$unresolved, c(FALSE, TRUE, TRUE, TRUE, FALSE, FALSE))
  # Wilson 95% interval against the closed form: 180 of 200 covered.
  expect_equal(m$coverage[2], 0.9)
  expect_equal(c(m$coverage_lower[2], m$coverage_upper[2]), wilson(180, 200))
  expect_equal(c(m$coverage_lower[1], m$coverage_upper[1]), wilson(40, 40))
  expect_true(m$coverage_lower[1] > 0.9)
  expect_true(is.finite(m$coverage_mcse[2]))
  # The failure rule is unchanged: more than 5% failed fits fails the cell.
  expect_equal(m$fail_rate[6], 0.1)
  expect_identical(m$reps, c(40L, 200L, 200L, 40L, 40L, 40L))
})

test_that("the verdict is feasible or infeasible and flags a decision resting on unresolved cells", {
  truths <- data.frame(surface = "s", estimand = "ATE", truth = 0.10)
  p <- fast_plan(estimands = "ATE", tolerance = list(bias = 0.02, coverage = 0.9))
  # Feasible; the selected candidate (smallest RMSE) has an unresolved cell.
  m <- .sim_metrics(rbind(sim_cell("clear", 40),
                          sim_cell("cov90", 200, spread = 0.005, covered = 180)), truths, p)
  v <- .sim_verdict(m, p)
  expect_identical(v$status, "feasible")
  expect_identical(v$selected, "cov90_t0.01")
  expect_true(v$unresolved)
  # The same metrics with `select` choosing the resolved candidate.
  p_sel <- fast_plan(estimands = "ATE", tolerance = list(bias = 0.02, coverage = 0.9),
                     select = function(m) "clear_t0.01")
  v <- .sim_verdict(m, p_sel)
  expect_identical(v$status, "feasible")
  expect_false(v$unresolved)
  # Infeasible, but one candidate fails only on an unresolved cell.
  m <- .sim_metrics(rbind(sim_cell("cov89", 200, covered = 178),
                          sim_cell("far", 40, bias = 0.2)), truths, p)
  v <- .sim_verdict(m, p)
  expect_identical(v$status, "infeasible")
  expect_true(is.na(v$selected))
  expect_true(v$unresolved)
  # Infeasible and settled: every failure is clear.
  m <- .sim_metrics(rbind(sim_cell("far", 40, bias = 0.2),
                          sim_cell("fails", 36, failed = 4L)), truths, p)
  v <- .sim_verdict(m, p)
  expect_identical(v$status, "infeasible")
  expect_false(v$unresolved)
})

test_that("a constant simulated outcome is a failed repetition, not an error", {
  d <- make_design(300)
  X <- .design_matrix(d$design, covs)
  g <- rep(0.5, 300)
  s <- list(flat = list(q1 = rep(1e-9, 300), q0 = rep(1e-9, 300)))
  out <- .one_rep(1L, X, g, s, fast_plan())
  expect_true(all(out$results$failed))
})

test_that("good overlap does not make the ATE infeasible", {
  d <- make_design(600, seed = 3, strength = 0.3)
  p <- fast_plan(tolerance = list(bias = 0.05, coverage = 0.75), reps = 30L)
  lk <- create_analysis_lock(d$design, "A", covs, p)
  sim <- suppressMessages(simulate_design(lk, assess_design(lk)))
  expect_s3_class(sim, "cr_simulation")
  expect_identical(sim$verdict$status[sim$verdict$estimand == "ATE"], "feasible")
  expect_true(.check_stamp(lk, sim, "simulation"))
})

test_that("severe positivity makes the ATE infeasible while the ATO survives", {
  d <- make_design(600, seed = 4, strength = 4)
  p <- fast_plan(tolerance = list(bias = 0.02, coverage = 0.9), reps = 30L,
                 surfaces = list(heterogeneity = 2))
  lk <- create_analysis_lock(d$design, "A", covs, p)
  sim <- suppressMessages(simulate_design(lk, assess_design(lk)))
  st <- stats::setNames(sim$verdict$status, sim$verdict$estimand)
  expect_identical(unname(st["ATE"]), "infeasible")
  expect_identical(unname(st["ATO"]), "feasible")
})

test_that("a repetition's draws do not depend on which surface-library pairs run", {
  d <- make_design(300)
  X <- .design_matrix(d$design, covs)
  g <- stats::plogis(-0.2 + 0.8 * d$design$w1)
  p <- fast_plan(estimands = c("ATE", "ATO"), surfaces = list(forms = "linear"),
                 candidates = list(library = list(glm = "glm", glm2 = "glm"),
                                   truncation = 0.05))
  s <- .surfaces(p, X, g, p$seed)
  full <- .one_rep(3L, X, g, s, p)$results
  sub <- .one_rep(3L, X, g, s, p, surfaces_run = names(s)[2], libraries_run = "glm2")$results
  expect_identical(unique(paste(sub$surface, sub$library)), paste(names(s)[2], "glm2"))
  keep <- full$surface == names(s)[2] & full$library == "glm2"
  ref <- full[keep, , drop = FALSE]
  rownames(ref) <- NULL; rownames(sub) <- NULL
  expect_identical(sub, ref)
  # An explicit, non-crossed set of pairs gives the same rows as well.
  pr <- data.frame(surface = names(s), library = c("glm2", "glm"))
  sub2 <- .one_rep(3L, X, g, s, p, pairs = pr)$results
  ref2 <- full[paste(full$surface, full$library) %in% paste(pr$surface, pr$library), ]
  rownames(ref2) <- NULL; rownames(sub2) <- NULL
  expect_identical(sub2[order(sub2$surface, sub2$library), ],
                   ref2[order(ref2$surface, ref2$library), ], ignore_attr = TRUE)
})

test_that("unresolved cells at or above the primary get extra repetitions up to max_reps", {
  d <- make_design(300, seed = 5, strength = 0.5)
  # Four repetitions cannot resolve a coverage tolerance of 0.75 (the Wilson
  # interval for 4 of 4 starts near 0.51), so the ATE cells are unresolved.
  p <- fast_plan(estimands = c("ATE", "ATO"), reps = 4L, max_reps = 12L,
                 tolerance = list(bias = 0.5, coverage = 0.75),
                 surfaces = list(forms = "linear"),
                 candidates = list(library = list(glm = "glm"), truncation = 0.05))
  lk <- create_analysis_lock(d$design, "A", covs, p)
  msgs <- character()
  sim <- withCallingHandlers(simulate_design(lk, assess_design(lk)),
                             message = function(m) {
                               msgs <<- c(msgs, conditionMessage(m))
                               invokeRestart("muffleMessage")
                             })
  m <- sim$metrics
  v <- sim$verdict
  pos <- match("feasible", v$status)
  expect_identical(pos, 1L)  # the ATE is primary, so the ATO sits below it
  expect_true(any(m$reps > p$reps))
  expect_true(all(m$reps <= p$max_reps))
  expect_true(all(m$reps[m$estimand == "ATO"] == p$reps))
  expect_true(nrow(sim$extension) > 0)
  expect_named(sim$extension, c("batch", "reps_from", "reps_to", "pairs",
                                "cells_unresolved_before"))
  expect_identical(sim$extension$reps_from[1], p$reps + 1L)
  expect_identical(sim$max_reps, 12L)
  expect_match(sim$decision_note, "max_reps")
  expect_true(any(grepl("extending .* surface-library pair", msgs)))
  # No extension when max_reps equals reps.
  p0 <- fast_plan(estimands = c("ATE", "ATO"), reps = 4L,
                  tolerance = list(bias = 0.5, coverage = 0.75),
                  surfaces = list(forms = "linear"),
                  candidates = list(library = list(glm = "glm"), truncation = 0.05))
  lk0 <- create_analysis_lock(d$design, "A", covs, p0)
  sim0 <- suppressMessages(simulate_design(lk0, assess_design(lk0)))
  expect_true(all(sim0$metrics$reps == 4L))
  expect_identical(nrow(sim0$extension), 0L)
})

test_that("the simulation uses the plan's K for its folds", {
  d <- make_design(300)
  p <- fast_plan(K = 3L, reps = 2L, estimands = "ATE")
  lk <- create_analysis_lock(d$design, "A", covs, p)
  ds <- assess_design(lk)
  seen <- integer()
  orig <- .make_folds
  local_mocked_bindings(.make_folds = function(n, K, seed) {
    seen <<- c(seen, K); orig(n, K, seed)
  })
  suppressMessages(simulate_design(lk, ds))
  expect_true(length(seen) > 0)
  expect_true(all(seen == 3L))
})

test_that("the simulated world uses the single-fit propensity score, not the cross-fitted one", {
  d <- make_design(300)
  p <- fast_plan(reps = 2L, estimands = "ATE")
  lk <- create_analysis_lock(d$design, "A", covs, p)
  ds <- assess_design(lk)
  seen <- list()
  orig_s <- .surfaces; orig_t <- .true_values; orig_r <- .one_rep
  local_mocked_bindings(
    .surfaces = function(plan, X, g, seed) { seen$surfaces <<- g; orig_s(plan, X, g, seed) },
    .true_values = function(q1, q0, g, band) { seen$truth <<- g; orig_t(q1, q0, g, band) },
    .one_rep = function(r, X, g, surfaces, plan, ...) {
      seen$draw <<- g; orig_r(r, X, g, surfaces, plan, ...)
    })
  suppressMessages(simulate_design(lk, ds))
  for (nm in c("surfaces", "truth", "draw")) expect_identical(seen[[nm]], ds$g_world, label = nm)
})

test_that("bootstrap copies of one row never straddle outer folds", {
  n <- 600
  idx <- withr::with_seed(7, sample.int(n, n, replace = TRUE))
  folds <- .make_folds(n, 5L, 11L)[idx]
  expect_true(all(tapply(folds, idx, function(f) length(unique(f))) == 1L))
  expect_identical(sort(unique(folds)), 1:5)
})

test_that("a repetition with a NaN interval counts as failed, not as a crash", {
  res <- data.frame(rep = 1:40, surface = "s", library = "glm", estimand = "ATE",
                    truncation = 0.01, failed = FALSE,
                    estimate = rep(c(0.09, 0.11), 20), ci_lower = 0.0, ci_upper = 0.2)
  res$ci_lower[3] <- NaN
  res$ci_upper[5] <- NA
  truths <- data.frame(surface = "s", estimand = "ATE", truth = 0.10)
  p <- fast_plan(estimands = "ATE", tolerance = list(bias = 0.02, coverage = 0.9))
  m <- .sim_metrics(res, truths, p)
  expect_identical(m$failed, 2L)
  expect_equal(m$coverage, 1)
  expect_false(is.na(m$pass))
  expect_no_error(.sim_verdict(m, p))
})
