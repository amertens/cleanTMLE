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

test_that("populations gives each estimand's covariate means and shift", {
  X <- cbind(x = c(0, 0, 1, 1, 2, 2), z = c(1, 0, 1, 0, 1, 0))
  A <- c(0, 0, 0, 1, 1, 1)
  g <- c(0.1, 0.2, 0.5, 0.5, 0.8, 0.9)
  # Each patient five times, so the two removed rows become ten removed patients.
  rows <- rep(1:6, each = 5)
  X <- X[rows, ]; A <- A[rows]; g <- g[rows]
  p <- .populations(X, A, g, c("ATE", "trimmed_ATE", "ATT", "ATO"), c(0.15, 0.85))
  expect_identical(p$covariate, c("x", "z"))
  expect_equal(p$cohort, unname(colMeans(X)))
  expect_equal(p$ATE, p$cohort)
  expect_equal(p$ATE_shift, c(0, 0))
  expect_equal(p$trimmed_ATE, unname(colMeans(X[g >= 0.15 & g <= 0.85, ])))
  expect_equal(p$ATT, unname(colMeans(X[A == 1, ])))
  h <- g * (1 - g)
  expect_equal(p$ATO, unname(colSums(X * h) / sum(h)))
  expect_equal(p$ATT_shift, (p$ATT - p$cohort) / unname(apply(X, 2, stats::sd)))
})

test_that("trimmed_profile counts removed patients by arm and suppresses small groups", {
  set.seed(3)
  X <- cbind(x = stats::rnorm(40))
  A <- rep(0:1, 20)
  g <- c(rep(0.01, 6), rep(0.5, 34))           # six rows below the band
  tp <- .trimmed_profile(X, A, g, c(0.05, 0.95))
  expect_identical(tp$counts$group, c("all", "treated", "control"))
  expect_identical(tp$counts$removed, c(6L, 3L, 3L))
  expect_identical(tp$counts$kept, c(34L, 17L, 17L))
  expect_equal(tp$covariates$removed, mean(X[1:6, 1]))
  expect_equal(tp$covariates$std_diff,
               (mean(X[1:6, 1]) - mean(X[7:40, 1])) / stats::sd(X[, 1]))
  few <- .trimmed_profile(X, A, c(0.01, rep(0.5, 39)), c(0.05, 0.95))
  expect_identical(few$counts$removed[1], 1L)
  expect_true(is.na(few$covariates$removed))
  expect_true(is.na(few$covariates$std_diff))
  two <- .trimmed_profile(X, A, c(0.01, 0.99, rep(0.5, 38)), c(0.05, 0.95))
  expect_identical(two$counts$removed[1], 2L)
  expect_true(is.na(two$covariates$removed))
  expect_true(is.na(two$covariates$kept))
})

test_that("with one to four patients removed, no table lets the reader difference them out", {
  set.seed(4)
  X <- cbind(x = stats::rnorm(40), z = stats::rbinom(40, 1, 0.5))
  A <- rep(0:1, 20)
  ests <- c("ATE", "trimmed_ATE", "ATT", "ATO")
  for (k in c(1L, 4L)) {
    g <- c(rep(0.01, k), rep(0.5, 40 - k))
    p <- .populations(X, A, g, ests, c(0.05, 0.95))
    tp <- .trimmed_profile(X, A, g, c(0.05, 0.95))
    expect_true(all(is.na(p$trimmed_ATE)))
    expect_true(all(is.na(p$trimmed_ATE_shift)))
    expect_true(all(is.na(tp$covariates$removed)))
    expect_true(all(is.na(tp$covariates$kept)))
    expect_true(all(is.na(tp$covariates$std_diff)))
    expect_identical(tp$counts$removed[1], k)     # the counts stay
    expect_false(anyNA(p$cohort))
    expect_false(anyNA(p$ATE)); expect_false(anyNA(p$ATT)); expect_false(anyNA(p$ATO))
  }
  # Five removed patients is the smallest group that is reported.
  g5 <- c(rep(0.01, 5), rep(0.5, 35))
  p5 <- .populations(X, A, g5, ests, c(0.05, 0.95))
  tp5 <- .trimmed_profile(X, A, g5, c(0.05, 0.95))
  expect_equal(p5$trimmed_ATE, unname(colMeans(X[6:40, ])))
  expect_equal(tp5$covariates$kept, unname(colMeans(X[6:40, ])))
  expect_false(anyNA(tp5$covariates$std_diff))
  # Nobody removed: the trimmed population is the cohort, and nothing is hidden.
  p0 <- .populations(X, A, rep(0.5, 40), ests, c(0.05, 0.95))
  expect_equal(p0$trimmed_ATE, p0$cohort)
  expect_equal(p0$trimmed_ATE_shift, c(0, 0))
})

test_that("assess_design carries the population tables, with and without a trim band", {
  d <- make_design(300)
  d$design$site <- factor(rep(c("A", "B", "C"), 100))
  covs <- c("w1", "w2", "site")
  ds <- assess_design(create_analysis_lock(d$design, "A", covs, fast_plan()))
  expect_setequal(ds$populations$covariate, c("w1", "w2", "siteB", "siteC"))
  expect_equal(ds$populations$cohort[ds$populations$covariate == "siteB"], 1 / 3)
  expect_true(all(c("trimmed_ATE", "ATT", "ATO") %in% names(ds$populations)))
  expect_type(ds$trimmed_profile, "list")
  nt <- assess_design(create_analysis_lock(d$design, "A", covs,
                                           fast_plan(estimands = c("ATE", "ATO"), trim_band = NULL)))
  expect_false("trimmed_ATE" %in% names(nt$populations))
  expect_null(nt$trimmed_profile)
})

test_that("the dossier carries the population tables", {
  fx <- fixture_pipeline()
  expect_true(is.data.frame(fx$dossier$design$populations))
  expect_type(fx$dossier$design$trimmed_profile, "list")
})
