# Validation tests against simple reference implementations and against
# theoretical bias direction. These do not certify the estimator; they
# guard against silent regressions.

# 1. clean_weight_diagnostics ESS matches manual ESS ----------------

test_that("ESS matches the standard sum(w)^2 / sum(w^2) formula", {
  set.seed(101)
  for (n in c(50, 200, 1000)) {
    w <- runif(n, 0.5, 2.5)
    ref <- sum(w)^2 / sum(w^2)
    out <- clean_weight_diagnostics(w)
    expect_equal(out$ess$overall, ref, tolerance = 1e-10,
                 info = paste("n =", n))
  }
})

# 2. SMD reference: complete-case manual computation ---------------

test_that("Unweighted SMD matches direct computation on simulated data", {
  set.seed(202)
  n <- 500
  A <- rbinom(n, 1, 0.4)
  X <- data.frame(age = rnorm(n) + A * 0.5,
                  bmi = rnorm(n) + A * 0.2)
  ref_smd <- function(x, A) {
    m1 <- mean(x[A == 1]); m0 <- mean(x[A == 0])
    s1 <- var(x[A == 1]);  s0 <- var(x[A == 0])
    abs(m1 - m0) / sqrt((s1 + s0) / 2)
  }
  ref_age <- ref_smd(X$age, A)
  ref_bmi <- ref_smd(X$bmi, A)

  w <- rep(1, n)
  out <- clean_weight_diagnostics(w, treatment = A, covariates = X)
  expect_equal(out$smd$smd_unweight[out$smd$covariate == "age"],
               ref_age, tolerance = 1e-10)
  expect_equal(out$smd$smd_unweight[out$smd$covariate == "bmi"],
               ref_bmi, tolerance = 1e-10)
})

# 3. (Coherence-check validation lives in test-governance-notes.R.)

# 4. DQ stress directional check: bias under outcome misclassification
#    should attenuate a positive risk difference toward 0 -----------------

test_that("Outcome misclassification attenuates the simulated effect toward 0", {
  set.seed(303)
  # A two-arm Bernoulli population with a known additive RD = 0.10.
  n  <- 4000
  A  <- rbinom(n, 1, 0.5)
  Y0 <- rbinom(n, 1, 0.30)
  Y1 <- rbinom(n, 1, 0.40)
  Y  <- ifelse(A == 1, Y1, Y0)

  rd_observed <- mean(Y[A == 1]) - mean(Y[A == 0])

  # Apply non-differential outcome misclassification with sens = 0.80,
  # spec = 0.95. Theory says the observed RD shrinks toward 0.
  flip <- function(y, sens, spec) {
    obs <- y
    pos <- which(y == 1)
    neg <- which(y == 0)
    obs[pos] <- rbinom(length(pos), 1, sens)
    obs[neg] <- 1 - rbinom(length(neg), 1, spec)
    obs
  }
  Y_mis <- flip(Y, sens = 0.80, spec = 0.95)
  rd_misclassified <- mean(Y_mis[A == 1]) - mean(Y_mis[A == 0])

  # Direction: |rd_misclassified| should be smaller than |rd_observed|.
  expect_lt(abs(rd_misclassified), abs(rd_observed))
})

# 5. (The checkpoint_weights validation died with the checkpoint layer
#    in 0.3.0; clean_weight_diagnostics is validated in sections 1-2.)

# 6. (Risk-report-table validation lives in test-governance-notes.R.)
