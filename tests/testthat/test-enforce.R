# The enforce = TRUE contract: since 0.3.0 the single institutional
# switch is create_analysis_lock(enforce = TRUE), and the single act of
# authorisation is unmask_outcome(lock, original_lock, approved_by =)
# with a named approver written into the design log. There is no token,
# no gate object, and no audit fingerprint; this file replaces
# test-clean-tmle-authorization.R.

.mk_enforce_fixture <- function(n = 300L, seed = 11L) {
  d <- sim_func1(n = n, seed = seed)
  lock <- create_analysis_lock(d, "treatment", "event_24",
                               c("age", "sex", "biomarker"),
                               seed = seed, enforce = TRUE)
  list(data = d, lock = lock, masked = mask_outcome(lock))
}

test_that("an enforce lock refuses Stage 4 estimation until authorised", {
  fx <- .mk_enforce_fixture()
  expect_error(
    run_crude_workflow(fx$lock),
    "enforce = TRUE", fixed = TRUE)
})

test_that("unmask_outcome on an enforce lock requires a named approver", {
  fx <- .mk_enforce_fixture()
  expect_error(unmask_outcome(fx$masked, fx$lock),
               "approved_by", fixed = TRUE)
  expect_error(unmask_outcome(fx$masked, fx$lock, approved_by = "  "),
               "approved_by", fixed = TRUE)
})

test_that("a named unmasking authorises estimation and lands in the design log", {
  fx <- .mk_enforce_fixture()
  un <- unmask_outcome(fx$masked, fx$lock, approved_by = "review team")
  expect_false(isTRUE(un$.outcome_masked))
  expect_true(isTRUE(un$.outcome_authorized))
  expect_identical(cleanTMLE:::.outcome_vector(un),
                   cleanTMLE:::.outcome_vector(fx$lock))

  log <- un$design_log
  row <- log[log$type == "outcome_unmasked", , drop = FALSE]
  expect_equal(nrow(row), 1L)
  expect_equal(row$decided_by, "review team")
  expect_equal(row$stage, "Stage 4 (unmasking)")

  res <- run_crude_workflow(un)
  expect_true(is.numeric(res$estimate))
})

test_that("the deprecated allow_unauthorized maps to a forced approval", {
  fx <- .mk_enforce_fixture()
  expect_warning(
    un <- unmask_outcome(fx$masked, fx$lock, allow_unauthorized = TRUE),
    "deprecated")
  expect_true(isTRUE(un$.outcome_authorized))
  row <- un$design_log[un$design_log$type == "outcome_unmasked", ,
                       drop = FALSE]
  expect_match(row$decided_by, "forced")
})

test_that("a default (enforce = FALSE) lock estimates without unmasking", {
  d <- sim_func1(n = 300, seed = 12)
  lock <- create_analysis_lock(d, "treatment", "event_24",
                               c("age", "sex", "biomarker"), seed = 12)
  res <- run_crude_workflow(lock)
  expect_true(is.numeric(res$estimate))
})

test_that("allow_outcome_access overrides the enforce guard per call", {
  fx <- .mk_enforce_fixture()
  res <- run_crude_workflow(fx$lock, allow_outcome_access = TRUE)
  expect_true(is.numeric(res$estimate))
})

test_that("run_clean_tmle labels its internal lock as unguarded (cleanroom_enabled = FALSE)", {
  d <- sim_func1(n = 300, seed = 5)
  res <- suppressMessages(run_clean_tmle(
    data = d, Avar = "treatment", Yvar = "event_24",
    covariates = c("age", "sex", "biomarker"),
    learner_lib = c("SL.glm", "SL.mean"), truncation = c(0.01, 0.99),
    seed = 5, verbose = FALSE))
  expect_false(isTRUE(res$lock$cleanroom_enabled))
})