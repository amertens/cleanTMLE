# design_report(), the rwe_wide library preset, and the new SL learners.

test_that("SL.glm.pca fits, predicts, and round-trips through SuperLearner", {
  skip_if_not_installed("SuperLearner")
  set.seed(41)
  n <- 300
  X <- as.data.frame(matrix(stats::rnorm(n * 8), n, 8))
  y <- stats::rbinom(n, 1, stats::plogis(0.8 * X$V1 - 0.5 * X$V2))
  fit <- SL.glm.pca(Y = y, X = X, newX = X, family = stats::binomial(),
                    obsWeights = rep(1, n), k = 3)
  expect_length(fit$pred, n)
  expect_true(all(fit$pred >= 0 & fit$pred <= 1))
  p2 <- predict.SL.glm.pca(fit$fit, newdata = X)
  expect_equal(p2, fit$pred, tolerance = 1e-8)
  # Inside SuperLearner, alongside the glmnet variants.
  sl <- SuperLearner::SuperLearner(
    Y = y, X = X, family = stats::binomial(),
    SL.library = c("SL.mean", "SL.glm.pca5", "SL.glmnet.ridge"),
    env = cleanTMLE:::.cleantmle_sl_env())
  expect_true(all(is.finite(sl$SL.predict)))
})

test_that("the rwe_wide preset returns the smooth wide-design library", {
  lib <- build_sl_library(role = "g", n_eff = 1200, preset = "rwe_wide")
  flat <- unlist(lib$library)
  expect_true(all(c("SL.mean", "SL.glmnet", "SL.glmnet.ridge",
                    "SL.glmnet.enet", "SL.bayesglm", "SL.glm.pca5",
                    "SL.glm.pca10") %in% flat))
  expect_false(any(grepl("ranger|xgboost", flat)))
  expect_true(any(vapply(lib$library, function(x)
    length(x) == 2 && x[2] == "screen.corP", logical(1))))
  expect_true(lib$stratify_cv)
})

test_that("design_report reads the ladder against support and feasibility", {
  set.seed(41)
  n <- 600
  x1 <- stats::rnorm(n)
  g <- stats::plogis(-2.5 + 3.5 * x1)
  A <- stats::rbinom(n, 1, g)
  y <- stats::rbinom(n, 1, 0.2)
  dat <- data.frame(x1 = x1, treatment = A, outcome = y)
  lock <- create_analysis_lock(dat, "treatment", "outcome", "x1")
  lock <- declare_estimand_ladder(lock, primary = "ATE",
                                  fallbacks = c("trimmed_ATE", "ATT", "ATO"))
  psf <- fit_ps(lock, "glm")
  sup <- assess_support(psf, tree_search = FALSE)
  fea <- estimand_feasibility(psf)
  rep_ <- design_report(lock, sup, fea)
  expect_s3_class(rep_, "design_report")
  expect_false(rep_$primary_supported)
  expect_match(rep_$recommendation, "NOT supported")
  expect_match(rep_$recommendation, "logged design decision")
  expect_true("ATO" %in% rep_$feasible_estimands)
  expect_output(print(rep_), "Design report")

  # A well-supported design recommends proceeding.
  set.seed(42)
  g2 <- stats::plogis(0.3 * x1)
  A2 <- stats::rbinom(n, 1, g2)
  dat2 <- data.frame(x1 = x1, treatment = A2, outcome = y)
  lock2 <- create_analysis_lock(dat2, "treatment", "outcome", "x1")
  lock2 <- declare_estimand_ladder(lock2, primary = "ATE")
  psf2 <- fit_ps(lock2, "glm")
  rep2 <- design_report(lock2, assess_support(psf2, tree_search = FALSE),
                        estimand_feasibility(psf2))
  expect_true(rep2$primary_supported)
  expect_match(rep2$recommendation, "Proceed to estimation")
})

test_that("emulation_table renders the TARGET two-column table from the lock", {
  dat <- data.frame(x1 = rnorm(200), treatment = rbinom(200, 1, 0.5),
                    outcome = rbinom(200, 1, 0.2))
  lock <- create_analysis_lock(dat, "treatment", "outcome", "x1")
  lock <- declare_estimand_ladder(lock, primary = "ATE",
                                  fallbacks = c("ATT", "ATO"))
  tab <- cleanTMLE:::emulation_table(lock,
                         protocol = c(eligibility = "adults with trauma"))
  expect_s3_class(tab, "data.frame")
  expect_identical(names(tab), c("component", "target_trial", "emulation"))
  expect_equal(nrow(tab), 7L)
  expect_identical(tab$target_trial[1], "adults with trauma")
  expect_match(tab$emulation[tab$component == "Causal contrast (estimand)"],
               "primary ATE")
  expect_match(tab$emulation[tab$component == "Outcome"], "outcome")
})
