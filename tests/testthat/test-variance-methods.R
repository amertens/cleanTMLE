# Tests for the native nonparametric bootstrap variance and the FIORD
# second-stage variance-method selector.

make_conf_data <- function(n = 1500, seed = 1) {
  set.seed(seed)
  W1 <- rnorm(n); W2 <- rbinom(n, 1, 0.5)
  A  <- rbinom(n, 1, plogis(-0.2 + 0.8 * W1 + 0.5 * W2))   # confounded
  pY <- function(a) plogis(-1 + 0.6 * a + 0.7 * W1 + 0.5 * W2)
  Y  <- rbinom(n, 1, pY(A))
  list(data = data.frame(A = A, Y = Y, W1 = W1, W2 = W2),
       truth = mean(pY(1)) - mean(pY(0)))
}

test_that("bootstrap_rd_variance returns a valid SE and interval (IPTW)", {
  d <- make_conf_data()
  bs <- bootstrap_rd_variance(d$data, "A", "Y", c("W1", "W2"),
                              estimator = "iptw", R = 300L, seed = 7L)
  expect_true(is.finite(bs$se) && bs$se > 0)
  expect_length(bs$ci, 2)
  expect_lt(bs$ci[1], bs$ci[2])
  expect_gte(bs$R_effective, 250L)
  # Point estimate recovers the truth within sampling error.
  expect_lt(abs(bs$estimate - d$truth), 0.05)
  # The percentile interval covers the truth.
  expect_lte(bs$ci[1], d$truth)
  expect_gte(bs$ci[2], d$truth)
})

test_that("bootstrap_rd_variance works for TMLE", {
  skip_if_not_installed("tmle")
  skip_if_not_installed("SuperLearner")
  d <- make_conf_data(n = 1000, seed = 3)
  bs <- bootstrap_rd_variance(d$data, "A", "Y", c("W1", "W2"),
                              estimator = "tmle", R = 40L,
                              sl_library = "SL.glm", seed = 5L)
  expect_true(is.finite(bs$se) && bs$se > 0)
  expect_lt(abs(bs$estimate - d$truth), 0.06)
})

test_that("bootstrap_rd_variance works for match_tmle", {
  skip_if_not_installed("tmle")
  skip_if_not_installed("SuperLearner")
  d <- make_conf_data(n = 1200, seed = 7)
  bs <- bootstrap_rd_variance(d$data, "A", "Y", c("W1", "W2"),
                              estimator = "match_tmle", R = 30L,
                              sl_library = "SL.glm", seed = 11L)
  expect_true(is.finite(bs$se) && bs$se > 0)
  expect_length(bs$ci, 2)
  expect_lt(bs$ci[1], bs$ci[2])
  expect_gte(bs$R_effective, 20L)
  # Point estimate should be close to truth
  expect_lt(abs(bs$estimate - d$truth), 0.08)
})

test_that("select_variance_method returns a method and per-method coverage", {
  # Build a small set of synthetic datasets with a known truth.
  set.seed(99)
  truth <- 0.10
  dsets <- lapply(1:15, function(i) {
    n <- 1200
    W1 <- rnorm(n); W2 <- rbinom(n, 1, 0.5)
    A  <- rbinom(n, 1, plogis(0.3 * W1 + 0.2 * W2))
    # Construct Y so the marginal RD is approximately `truth`.
    p0 <- plogis(-0.8 + 0.5 * W1 + 0.4 * W2)
    p1 <- pmin(p0 + truth, 0.99)
    Y  <- rbinom(n, 1, ifelse(A == 1, p1, p0))
    data.frame(A = A, Y = Y, W1 = W1, W2 = W2)
  })
  sel <- select_variance_method(dsets, truth = truth,
                                treatment = "A", outcome = "Y",
                                covariates = c("W1", "W2"),
                                estimator = "iptw",
                                methods = c("influence", "bootstrap"),
                                R_boot = 100L)
  expect_true(sel$selected %in% c("influence", "bootstrap"))
  expect_true(all(c("influence", "bootstrap") %in% names(sel$coverage)))
  expect_true(all(sel$coverage >= 0 & sel$coverage <= 1))
})
