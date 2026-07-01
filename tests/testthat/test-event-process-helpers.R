# Tests for weight diagnostics and the weight checkpoint. The event-process,
# target-population, missing-data, coherence, and risk-report formatters were
# split out into the cleanroomGov package (tested there).

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


# checkpoint_weights ---------------------------------------------------

test_that("checkpoint_weights returns GO when ESS and weights are clean", {
  set.seed(31)
  n <- 600
  A <- rbinom(n, 1, 0.5)
  ps <- plogis(0.1 * rnorm(n))
  w  <- ifelse(A == 1, 1 / ps, 1 / (1 - ps))
  cp <- checkpoint_weights(w, treatment = A,
                           max_weight_threshold = 50,
                           ess_floor = 0.2 * n)
  expect_s3_class(cp, "cleantmle_checkpoint")
  expect_equal(cp$decision, "GO")
})

test_that("checkpoint_weights flags a single failing criterion", {
  # Construct a weight vector with one extreme weight but acceptable ESS
  w <- c(rep(1, 990), rep(50, 10))
  cp <- checkpoint_weights(w, max_weight_threshold = 10,
                           ess_floor = 100,
                           extreme_prop_threshold = 0.001)
  expect_s3_class(cp, "cleantmle_checkpoint")
  # extreme prop 0.01 > 0.001 threshold -> ext_fail = TRUE; ESS still high
  expect_true(cp$decision %in% c("FLAG", "STOP"))
})

test_that("checkpoint_weights returns STOP when both criteria fail", {
  # Heavy concentration -> low ESS and many extreme weights
  w <- c(rep(0.01, 50), rep(100, 50))
  cp <- checkpoint_weights(w, max_weight_threshold = 10,
                           ess_floor = 90,
                           extreme_prop_threshold = 0.05)
  expect_equal(cp$decision, "STOP")
  expect_true(grepl("Weight checkpoint", cp$rationale))
})

test_that("checkpoint_weights metrics include weight_type label", {
  w <- runif(100, 0.5, 1.5)
  cp <- checkpoint_weights(w, weight_type = "censoring")
  expect_equal(cp$metrics$weight_type, "censoring")
})
