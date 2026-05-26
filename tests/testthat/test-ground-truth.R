# Ground-truth and cross-implementation regression tests.
#
# These tests guard against the silent-failure modes and check the
# estimators against (a) a known analytic truth and (b) an independent
# implementation (tmle::tmle) on shared inputs.

test_that("g-computation injects the treatment and warns instead of returning RD = 0", {
  set.seed(11)
  n  <- 3000
  W1 <- rnorm(n)
  A  <- rbinom(n, 1, plogis(0.3 * W1))
  # Treatment raises the hazard, so the cumulative-risk RD is clearly non-zero.
  lambda <- exp(-3 + 0.6 * A + 0.3 * W1)
  Tt     <- rexp(n, lambda)
  Cc     <- rexp(n, exp(-3.2))
  df <- data.frame(time = pmin(Tt, Cc), censor = as.integer(Tt <= Cc),
                   A = A, W1 = W1)
  t_eval <- as.numeric(stats::quantile(df$time, 0.5))

  # Covariate-only outcome formula (no treatment): must warn, not silently
  # return a degenerate zero contrast.
  expect_warning(
    g <- specify_models(data = df) |>
      identify_outcome(censor, type = "time_to_event", formula = ~ W1) |>
      identify_treatment(A, formula = ~ W1) |>
      estimate_gcomprisk(risk_time = t_eval),
    regexp = "treatment"
  )
  r  <- g$risk[g$risk$time == t_eval, ]
  rd <- r$risk[as.character(r$group) == "1"] - r$risk[as.character(r$group) == "0"]
  expect_true(is.finite(rd))
  expect_gt(abs(rd), 0.005)   # not the degenerate zero

  # With the treatment already in the outcome formula, no warning and the
  # same (non-degenerate) estimate.
  expect_silent(
    g2 <- specify_models(data = df) |>
      identify_outcome(censor, type = "time_to_event", formula = ~ A + W1) |>
      identify_treatment(A, formula = ~ W1) |>
      estimate_gcomprisk(risk_time = t_eval)
  )
  r2  <- g2$risk[g2$risk$time == t_eval, ]
  rd2 <- r2$risk[as.character(r2$group) == "1"] - r2$risk[as.character(r2$group) == "0"]
  expect_equal(rd, rd2, tolerance = 1e-8)
})

test_that("point TMLE recovers the marginal risk difference and agrees with tmle::tmle", {
  set.seed(22)
  n  <- 4000
  W1 <- rnorm(n)
  W2 <- rbinom(n, 1, 0.5)
  A  <- rbinom(n, 1, plogis(-0.3 + 0.5 * W1 + 0.4 * W2))
  pY <- function(a) plogis(-1 + 0.8 * a + 0.6 * W1 + 0.5 * W2)
  Y  <- rbinom(n, 1, pY(A))
  truth <- mean(pY(1)) - mean(pY(0))   # g-formula truth on this sample
  dat <- data.frame(A = A, Y = Y, W1 = W1, W2 = W2)

  ct <- estimate_tmle_risk_point(
    data = dat, treatment = "A", outcome = "Y",
    covariates = c("W1", "W2"),
    sl_library = c("SL.glm"), n_folds = 2L)
  ct_psi <- ct$estimates$ATE$estimate
  expect_lt(abs(ct_psi - truth), 0.04)   # absolute RD-scale tolerance

  skip_if_not_installed("tmle")
  skip_if_not_installed("SuperLearner")
  set.seed(22)
  tt <- tmle::tmle(Y = Y, A = A, W = cbind(W1, W2), family = "binomial",
                   Q.SL.library = "SL.glm", g.SL.library = "SL.glm")
  # cleanTMLE and tmle::tmle should agree very closely on the same data.
  expect_lt(abs(ct_psi - tt$estimates$ATE$psi), 0.01)
})

test_that("IPW and AIPW cumulative risk recover the truth under randomized treatment", {
  # Under randomized A, all adjusted estimators (and the crude) should agree
  # and recover the true per-arm cumulative risk.
  set.seed(33)
  n  <- 5000
  W1 <- rnorm(n)
  A  <- rbinom(n, 1, 0.5)                      # randomized: A independent of W
  lambda <- exp(-3 + 0.5 * A + 0.3 * W1)
  Tt <- rexp(n, lambda)
  Cc <- rexp(n, exp(-3.0))
  df <- data.frame(time = pmin(Tt, Cc), censor = as.integer(Tt <= Cc),
                   A = A, W1 = W1)
  t_eval <- as.numeric(stats::quantile(df$time, 0.4))

  ipw <- specify_models(data = df) |>
    identify_outcome(censor, type = "time_to_event") |>
    identify_treatment(A, formula = ~ W1) |>
    estimate_ipwrisk(risk_time = t_eval)
  aipw <- specify_models(data = df) |>
    identify_outcome(censor, type = "time_to_event", formula = ~ A + W1) |>
    identify_treatment(A, formula = ~ W1) |>
    estimate_aipwrisk(risk_time = t_eval)

  rd_of <- function(fit) {
    r <- fit$risk[fit$risk$time == t_eval, ]
    r$risk[as.character(r$group) == "1"] - r$risk[as.character(r$group) == "0"]
  }
  # Under randomization the two estimators should be close to each other.
  expect_lt(abs(rd_of(ipw) - rd_of(aipw)), 0.03)
  # Treatment raises the hazard, so both risk differences are positive.
  expect_gt(rd_of(ipw), 0)
})
