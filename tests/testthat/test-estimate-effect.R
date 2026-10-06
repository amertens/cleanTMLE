fx <- fixture_pipeline()
ub <- unblind(fx$lock, fx$dossier, fx$d$outcomes, approved_by = "Review team")

test_that("estimate_effect passes the dossier's candidate, K and fold seed to fit_candidate", {
  seen <- NULL
  orig <- fit_candidate
  local_mocked_bindings(fit_candidate = function(X, A, Y, estimand, candidate, plan, folds,
                                                 seed, g_fit = NULL) {
    seen <<- list(estimand = estimand, candidate = candidate, folds = folds, seed = seed)
    orig(X, A, Y, estimand, candidate, plan, folds, seed, g_fit)
  })
  fit <- estimate_effect(ub)
  dec <- fx$dossier$decision
  expect_identical(seen$estimand, dec$primary)
  expect_identical(seen$candidate, dec$candidates[[dec$primary]])
  expect_identical(seen$seed, dec$fold_seed)
  expect_identical(seen$folds, .make_folds(length(seen$folds), dec$K, dec$fold_seed))
  expect_s3_class(fit, "cr_estimate")
  expect_true(fit$prespecified)
  expect_true(fit$stress_tested)
})

test_that("a non-primary estimand is allowed and labeled", {
  dec <- fx$dossier$decision
  alt <- setdiff(names(dec$candidates), dec$primary)[1]
  skip_if(is.na(alt), "no second feasible estimand in this fixture")
  fit <- estimate_effect(ub, estimand = alt)
  expect_false(fit$prespecified)
  expect_match(fit$note, "not the prespecified primary")
})

test_that("weighting matches a direct WeightIt call", {
  fit <- estimate_effect(ub, method = "weighting", estimand = "ATO")
  X <- .design_matrix(fx$lock$data, fx$lock$covariates)
  df <- data.frame(.Y = ub$outcomes$y, .A = fx$lock$data$A, X)
  w <- WeightIt::weightit(.A ~ w1 + w2 + w3, data = df, estimand = "ATO", method = "glm")
  ref <- WeightIt::lm_weightit(.Y ~ .A, data = df, weightit = w)
  expect_equal(fit$estimate, unname(stats::coef(ref)[".A"]))
  expect_false(fit$stress_tested)
})

test_that("matching estimates the ATT only", {
  skip_if_not_installed("MatchIt")
  expect_error(estimate_effect(ub, method = "matching", estimand = "ATE"), "ATT only")
  fit <- estimate_effect(ub, method = "matching", estimand = "ATT")
  expect_true(is.finite(fit$estimate))
})

test_that("factor and oddly named covariates run through every backend", {
  d <- make_design(400, seed = 8)
  withr::with_seed(8, {
    d$design$site <- factor(sample(c("a", "b", "c"), 400, replace = TRUE))
    d$design[["age group"]] <- stats::rnorm(400)
  })
  p <- fast_plan(tolerance = list(bias = 0.1, coverage = 0.5), reps = 20L)
  lk <- create_analysis_lock(d$design, "A", c("w1", "site", "age group"), p)
  ds <- assess_design(lk)
  dr <- design_report(lk, ds, suppressMessages(simulate_design(lk, ds)))
  expect_false(is.na(dr$decision$primary))
  u <- unblind(lk, dr, d$outcomes, approved_by = "Review team")
  expect_true(is.finite(estimate_effect(u, method = "weighting", estimand = "ATE")$estimate))
  expect_true(is.finite(estimate_effect(u)$estimate))
})

test_that("tmle without a candidate names the status the dossier gave the estimand", {
  d <- make_design(300, seed = 3)
  # A near-zero bias tolerance with few repetitions leaves no estimand feasible.
  p <- fast_plan(tolerance = list(bias = 1e-6, coverage = 0.99), reps = 4L)
  lk <- create_analysis_lock(d$design, "A", c("w1", "w2", "w3"), p)
  ds <- assess_design(lk)
  dr <- design_report(lk, ds, suppressMessages(simulate_design(lk, ds)))
  expect_true(is.na(dr$decision$primary))
  u <- unblind(lk, dr, d$outcomes, approved_by = "Review team", override = "test fixture")
  st <- dr$decision$status[["ATE"]]
  expect_false(identical(st, "feasible"))
  expect_error(estimate_effect(u, estimand = "ATE"),
               sprintf("rated ATE %s; only feasible estimands", st))
})

test_that("summary adds an E-value and tidy returns one row", {
  skip_if_not_installed("EValue")
  fit <- estimate_effect(ub)
  s <- summary(fit)
  expect_true(is.numeric(s$evalue["point"]))
  expect_equal(nrow(as.data.frame(fit)), 1L)
})

test_that("a time-to-event plan runs through concrete", {
  skip_on_cran()
  skip_if_not_installed("concrete")
  skip_if_not_installed("data.table")
  skip_if_not_installed("SuperLearner")
  d <- make_design(300, seed = 9)
  withr::with_seed(9, {
    t_event <- stats::rexp(300, 0.002 * exp(0.4 * d$design$A + 0.3 * d$design$w1))
    t_cens <- stats::rexp(300, 0.001)
  })
  oc <- data.frame(id = d$design$id, t = pmin(t_event, t_cens),
                   s = as.integer(t_event <= t_cens))
  p <- fast_plan(outcome = c(time = "t", event = "s"), estimands = "ATE",
                 target_time = 300, hazards = list(`0` = list(Surv(t, s == 0) ~ .),
                                                   `1` = list(Surv(t, s == 1) ~ .)),
                 tolerance = list(bias = 0.1, coverage = 0.5), reps = 10L)
  lk <- create_analysis_lock(d$design, "A", c("w1", "w2", "w3"), p)
  ds <- assess_design(lk)
  dr <- design_report(lk, ds, suppressMessages(simulate_design(lk, ds)))
  u <- unblind(lk, dr, oc, approved_by = "Review team", override = "fixture")
  fit <- estimate_effect(u, estimand = "ATE")
  expect_true(is.finite(fit$estimate))
  expect_identical(fit$method, "concrete")
})
