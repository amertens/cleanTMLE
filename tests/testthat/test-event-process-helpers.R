# Tests for the weight diagnostics. The pre-outcome weight checkpoint was
# removed with the checkpoint layer in 0.3.0; the governance formatters
# folded back and are tested in test-governance-notes.R.

# clean_weight_diagnostics ---------------------------------------------

test_that("weight diagnostics compute ESS correctly", {
  set.seed(7)
  n <- 500
  A <- rbinom(n, 1, 0.5)
  ps <- plogis(0.1 * rnorm(n))
  w  <- ifelse(A == 1, 1 / ps, 1 / (1 - ps))
  out <- clean_weight_diagnostics(w, treatment = A)
  manual_ess <- sum(w)^2 / sum(w^2)
  expect_equal(out$ess$overall, manual_ess, tolerance = 1e-9)
  expect_s3_class(out, "cleantmle_weight_diag")
  expect_true(all(c("min", "median", "max", "mean") %in%
                    names(out$percentiles$overall)))
})

test_that("SMDs are computed and finite", {
  set.seed(11)
  n <- 600
  A <- rbinom(n, 1, 0.4)
  X <- data.frame(age = rnorm(n) + A * 0.4,
                  sex = rbinom(n, 1, 0.5))
  ps <- plogis(0.5 * X$age + 0.3 * X$sex)
  w  <- ifelse(A == 1, 1 / ps, 1 / (1 - ps))
  out <- clean_weight_diagnostics(w, treatment = A, covariates = X)
  expect_true(all(c("smd_unweight", "smd_weighted") %in% names(out$smd)))
  expect_true(all(is.finite(out$smd$smd_unweight)))
  expect_true(all(is.finite(out$smd$smd_weighted)))
})

test_that("weight diagnostics flag low ESS", {
  w <- c(rep(1, 50), rep(100, 1))
  out <- clean_weight_diagnostics(w, ess_floor = 100)
  expect_true(out$flags$low_ess)
})
