sim_xy <- function(n = 300, seed = 1) {
  withr::with_seed(seed, {
    X <- cbind(w1 = stats::rnorm(n), w2 = stats::rbinom(n, 1, 0.5), w3 = stats::rnorm(n))
    y <- stats::rbinom(n, 1, stats::plogis(-0.3 + X[, 1] + 0.5 * X[, 3]^2))
    list(X = X, y = y)
  })
}

test_that("the default library has 17 settings, one learner type per family", {
  s <- .settings(.default_library)
  expect_length(s, 17L)
  expect_identical(sort(unique(vapply(s, `[[`, "", "type"))),
                   c("earth", "glmnet", "nnet", "xgboost"))
})

test_that("xgboost prediction at round k equals a k-round fit", {
  d <- sim_xy()
  st <- .settings("xgboost")[[1]]
  m50 <- .learner_fit(st, d$X, d$y, path_index = 50L)
  m20 <- .learner_fit(st, d$X, d$y, path_index = 20L)
  expect_equal(.learner_predict(m50, st, d$X, 20L), .learner_predict(m20, st, d$X, 20L))
})

test_that("the super learner returns bounded predictions, convex weights and risks", {
  d <- sim_xy()
  sl <- .super_learner(d$X, d$y, c("glm", "glmnet", "xgboost"), list(d$X[1:10, ]),
                       V = 3L, seed = 2L)
  expect_length(sl$pred[[1]], 10L)
  expect_true(all(sl$pred[[1]] > 0 & sl$pred[[1]] < 1))
  expect_equal(sum(sl$weights), 1)
  expect_true(all(sl$weights >= 0))
  expect_equal(nrow(sl$risks), 1L + 3L + 10L)
  expect_true(all(sl$risks$k[sl$risks$type == "xgboost"] <= 500L))
})

test_that("a failing learner is dropped and recorded, not fatal", {
  d <- sim_xy()
  orig <- .learner_fit
  local_mocked_bindings(.learner_fit = function(s, ...) {
    if (s$type == "earth") stop("boom") else orig(s, ...)
  })
  sl <- .super_learner(d$X, d$y, c("glm", "earth"), list(d$X), V = 2L, seed = 1L)
  expect_true(sl$risks$failed[sl$risks$type == "earth"])
  expect_identical(names(sl$weights), "glm")
})

test_that("a single glm skips cross-validation", {
  d <- sim_xy()
  sl <- .super_learner(d$X, d$y, "glm", list(d$X), V = 5L, seed = 1L)
  expect_identical(unname(sl$weights), 1)
})

test_that("the edge check flags a setting that keeps winning on a grid edge", {
  mk <- function(fold, depth, k) data.frame(
    fold = fold, setting = c("a", "b"), type = "xgboost",
    max_depth = c(depth, 1), eta = 0.1, alpha = NA, decay = NA,
    k = c(k, 10), P = 500L, cv_risk = c(0.4, 0.6), weight = NA,
    failed = FALSE, message = NA)
  risks <- rbind(mk(1, 5, 500L), mk(2, 5, 500L), mk(3, 3, 120L))
  e <- .edge_check(risks)
  row <- e[e$edge == "max_depth at the top of its grid", ]
  expect_equal(row$share, 2 / 3)
  expect_true(row$consistent)
  expect_true(e$consistent[e$edge == "number of trees at the cap"])
})

test_that(".nnloglik returns the simplex minimizer of the cross-validated log-loss", {
  withr::with_seed(11, {
    n <- 5000L
    eta <- stats::rnorm(n, 0, 1.5)
    y <- stats::rbinom(n, 1, stats::plogis(eta))
    Z <- cbind(shrunk = stats::plogis(0.5 * eta),
               noisy = stats::plogis(eta + stats::rnorm(n, 0, 1)))
  })
  w <- .nnloglik(Z, y)
  L <- stats::qlogis(.bound(Z, 1e-6))
  ll <- function(b) .logloss(y, stats::plogis(drop(L %*% b)))
  expect_true(all(w >= 0))
  expect_equal(sum(w), 1)
  expect_lte(ll(w), ll(c(1, 0)) + 1e-6)
  expect_lte(ll(w), ll(c(0, 1)) + 1e-6)
  expect_lte(ll(w), ll(c(0.5, 0.5)) + 1e-6)
})

test_that("grouped inner folds keep copies of one row together", {
  withr::local_seed(5)
  idx <- sample.int(200, 200, replace = TRUE)
  inner <- .inner_folds(200L, 5L, 3L, groups = idx)
  expect_length(inner, 200L)
  expect_true(all(tapply(inner, idx, function(f) length(unique(f))) == 1L))
  expect_identical(.inner_folds(50L, 3L, 3L), .make_folds(50L, 3L, 3L))
})

test_that("the super learner accepts groups and fails cleanly when a split has no data", {
  d <- sim_xy(200)
  idx <- withr::with_seed(2, sample.int(200, 200, replace = TRUE))
  sl <- .super_learner(d$X[idx, ], d$y[idx], c("glm", "glmnet"), list(d$X[1:5, ]),
                       V = 3L, seed = 1L, groups = idx)
  expect_true(all(sl$pred[[1]] > 0 & sl$pred[[1]] < 1))
  expect_error(.super_learner(d$X[rep(1, 50), ], rep(0:1, 25), c("glm", "glmnet"),
                              list(d$X[1:5, ]), V = 3L, seed = 1L, groups = rep(1, 50)),
               "Every learner failed")
})
