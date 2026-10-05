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

test_that("metrics apply the failure and borderline rules", {
  res <- data.frame(
    rep = rep(1:40, 2), surface = "s", library = "glm", estimand = "ATE",
    truncation = 0.01, failed = c(rep(FALSE, 40), rep(FALSE, 36), rep(TRUE, 4)),
    estimate = c(rep(c(0.09, 0.11), 20), rep(0.10, 40)),
    ci_lower = 0.0, ci_upper = 0.2)
  res$library[41:80] <- "sl"
  truths <- data.frame(surface = "s", estimand = "ATE", truth = 0.10)
  m <- .sim_metrics(res, truths, fast_plan(tolerance = list(bias = 0.02, coverage = 0.9)))
  expect_true(m$pass[m$library == "glm"])
  expect_false(m$borderline[m$library == "glm"])
  expect_equal(m$fail_rate[m$library == "sl"], 0.1)
  expect_false(m$pass[m$library == "sl"])
  v <- .sim_verdict(m, fast_plan(estimands = "ATE"))
  expect_identical(v$status, "feasible")
  expect_identical(v$selected, "glm_t0.01")
})

test_that("a constant simulated outcome is a failed repetition, not an error", {
  d <- make_design(300)
  X <- .design_matrix(d$design, covs)
  g <- rep(0.5, 300)
  s <- list(flat = list(q1 = rep(1e-9, 300), q0 = rep(1e-9, 300)))
  out <- .one_rep(1L, X, g, s, fast_plan())
  expect_true(all(out$results$failed))
})

test_that("good overlap gives a feasible ATE", {
  d <- make_design(600, seed = 3, strength = 0.3)
  p <- fast_plan(tolerance = list(bias = 0.05, coverage = 0.6), reps = 30L)
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
  expect_true(st["ATO"] %in% c("feasible", "borderline"))
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

test_that("bootstrap copies of one row never straddle outer folds", {
  n <- 600
  idx <- withr::with_seed(7, sample.int(n, n, replace = TRUE))
  folds <- .make_folds(n, 5L, 11L)[idx]
  expect_true(all(tapply(folds, idx, function(f) length(unique(f))) == 1L))
  expect_identical(sort(unique(folds)), 1:5)
})
