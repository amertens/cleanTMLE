test_that("the design matrix expands factors and keeps names", {
  df <- data.frame(w1 = c(1, 2, 3), site = factor(c("a", "b", "c")), `age group` = 1:3,
                   check.names = FALSE)
  X <- .design_matrix(df, c("w1", "site", "age group"))
  expect_equal(ncol(X), 4L)
  expect_true(any(grepl("age group", colnames(X), fixed = TRUE)))
  expect_true(all(c("siteb", "sitec") %in% colnames(X)))
  expect_type(X, "double")
})

test_that("cross-fitted g for fold 1 comes from a fit on fold 2 only", {
  d <- make_design(300)
  X <- .design_matrix(d$design, c("w1", "w2", "w3"))
  A <- d$design$A
  folds <- .make_folds(300, 2, 5)
  gf <- .fit_g(X, A, "glm", folds, 2L, 5L)
  m <- stats::glm(A[folds == 2] ~ X[folds == 2, ], family = stats::binomial())
  expect_equal(gf$g[folds == 1],
               .bound(stats::plogis(drop(cbind(1, X[folds == 1, ]) %*% stats::coef(m))), 1e-6),
               tolerance = 1e-6)
})

test_that("truncation reaches both the ATE and the ATT", {
  withr::local_seed(3)
  n <- 500
  X <- matrix(stats::rnorm(n), n, 1, dimnames = list(NULL, "w1"))
  g <- stats::plogis(3 * X[, 1])
  A <- stats::rbinom(n, 1, g)
  Y <- stats::rbinom(n, 1, stats::plogis(-1 + A + X[, 1]))
  Q1 <- stats::plogis(0.5 * X[, 1]); Q0 <- stats::plogis(-1 + 0.5 * X[, 1])
  a1 <- .target(Y, A, X, g, Q1, Q0, "ATE", 0.01, NULL)
  a2 <- .target(Y, A, X, g, Q1, Q0, "ATE", 0.10, NULL)
  t1 <- .target(Y, A, X, g, Q1, Q0, "ATT", 0.01, NULL)
  t2 <- .target(Y, A, X, g, Q1, Q0, "ATT", 0.10, NULL)
  expect_false(isTRUE(all.equal(a1$estimate, a2$estimate)))
  expect_false(isTRUE(all.equal(t1$estimate, t2$estimate)))
})

test_that("the ATO TMLE recovers the truth and agrees with the augmented estimator", {
  withr::local_seed(7)
  n <- 4000
  w <- stats::rnorm(n)
  g <- stats::plogis(1.2 * w)
  A <- stats::rbinom(n, 1, g)
  q <- function(a, w) stats::plogis(-0.5 + 0.8 * a + 0.6 * w + 0.4 * a * w)
  Y <- stats::rbinom(n, 1, q(A, w))
  r <- .target_ato(Y, A, g, q(1, w), q(0, w))
  wb <- stats::rnorm(1e6); gb <- stats::plogis(1.2 * wb); hb <- gb * (1 - gb)
  truth <- sum(hb * (q(1, wb) - q(0, wb))) / sum(hb)
  expect_lt(abs(r$estimate - truth), 3 * r$se)
  h <- g * (1 - g)
  num <- h * (q(1, w) - q(0, w)) + A * (1 - g) * (Y - q(1, w)) - (1 - A) * g * (Y - q(0, w))
  expect_lt(abs(r$estimate - mean(num) / mean(h)), 0.01)
  expect_equal(r$risk1 - r$risk0, r$estimate, tolerance = 1e-8)
})

test_that(".fit_and_target returns one row per estimand and truncation", {
  d <- make_design(400)
  X <- .design_matrix(d$design, c("w1", "w2", "w3"))
  p <- fast_plan()
  folds <- .make_folds(400, p$K, p$seed)
  out <- .fit_and_target(X, d$design$A, d$outcomes$y, "glm", p$estimands,
                         p$candidates$truncation, p, folds, p$seed)
  expect_equal(nrow(out), 8L)
  expect_false(any(out$failed))
  ato <- out[out$estimand == "ATO", ]
  expect_equal(ato$estimate[1], ato$estimate[2])
  expect_setequal(unique(attr(out, "risks")$nuisance), c("g", "Q"))
  one <- fit_candidate(X, d$design$A, d$outcomes$y, "ATE",
                       .candidates(p)[[1]], p, folds, p$seed)
  expect_equal(one$estimate, out$estimate[out$estimand == "ATE" & out$truncation == 0.01])
})

test_that("a target population with one arm fails as a row, not an error", {
  d <- make_design(300)
  X <- .design_matrix(d$design, c("w1", "w2", "w3"))
  p <- fast_plan(trim_band = c(0.001, 0.002))
  folds <- .make_folds(300, p$K, p$seed)
  out <- .fit_and_target(X, d$design$A, d$outcomes$y, "glm", "trimmed_ATE", 0.01,
                         p, folds, p$seed)
  expect_true(out$failed)
  expect_true(is.na(out$estimate))
})

test_that("the implausibility guard flags a sign flip", {
  y <- c(1, 1, 1, 0, 0, 0, 0, 0); a <- c(1, 1, 1, 1, 0, 0, 0, 0)
  expect_true(.implausibility_check(-0.4, y, a)$implausible)
  expect_false(.implausibility_check(0.6, y, a)$implausible)
})
