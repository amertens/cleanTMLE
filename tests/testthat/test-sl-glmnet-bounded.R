test_that("SL.glmnet.bounded fits and predicts on a near-separable design", {
  skip_if_not_installed("glmnet")
  set.seed(1)
  n <- 400L
  X <- data.frame(x1 = rnorm(n), x2 = rnorm(n), x3 = rnorm(n))
  # Near-separable treatment: a strong linear predictor pushes the propensity
  # toward 0/1 for many units, the regime where unbounded glmnet runs away.
  lp <- 6 * X$x1 + 4 * X$x2
  Y  <- rbinom(n, 1L, plogis(lp))

  elapsed <- system.time({
    fit <- cleanTMLE::SL.glmnet.bounded(
      Y = Y, X = X, newX = X, family = binomial())
  })[["elapsed"]]

  expect_type(fit, "list")
  expect_length(fit$pred, n)
  expect_true(all(fit$pred >= 0 & fit$pred <= 1))
  expect_s3_class(fit$fit, "SL.glmnet.bounded")
  # The bound is the point: the fit must return quickly even when separable.
  expect_lt(elapsed, 20)

  pred2 <- predict(fit$fit, newdata = X)
  expect_length(pred2, n)
  expect_true(all(pred2 >= 0 & pred2 <= 1))
})

test_that("run_plasmode_dq_stress accepts the near_positivity threat", {
  set.seed(2)
  n <- 600L
  W <- data.frame(age = rnorm(n, 55, 10), sex = rbinom(n, 1, 0.5),
                  bm = rnorm(n))
  lp_a <- 1.2 * (0.03 * (W$age - 55) + 0.8 * W$sex + 0.6 * W$bm)
  A <- rbinom(n, 1L, plogis(lp_a))
  Y <- rbinom(n, 1L, plogis(-1 + 0.3 * W$sex + 0.2 * W$bm - 0.3 * A))
  dat <- cbind(W, treatment = A, event_24 = Y)

  lock <- create_analysis_lock(
    dat, "treatment", "event_24", c("age", "sex", "bm"),
    sl_library = c("SL.glm", "SL.mean"), plasmode_reps = 5L, seed = 7L)

  cands <- list(
    tmle_candidate("light", g_library = "SL.glm", truncation = 0.001),
    tmle_candidate("heavy", g_library = "SL.glm", truncation = 0.20))

  dq <- run_plasmode_dq_stress(
    lock, tmle_candidates = cands, effect_sizes = c(0.05), reps = 5L,
    data_quality_scenarios = list(near_positivity = list(slopes = c(2, 3))),
    fit_timeout = 30, verbose = FALSE)

  expect_s3_class(dq, "plasmode_dq_results")
  expect_true("near_positivity" %in% dq$metrics$scenario)
  lvls <- unique(dq$metrics$level[dq$metrics$scenario == "near_positivity"])
  expect_setequal(lvls, c("slope_x2.0", "slope_x3.0"))
  # Under amplified positivity the light-truncation candidate should carry a
  # heavier weight tail and at least as much RMSE as the heavy candidate.
  worst <- dq$metrics[dq$metrics$scenario == "near_positivity", ]
  rmse_light <- max(worst$rmse[worst$candidate == "light"], na.rm = TRUE)
  rmse_heavy <- max(worst$rmse[worst$candidate == "heavy"], na.rm = TRUE)
  expect_gte(rmse_light, rmse_heavy * 0.99)
})
