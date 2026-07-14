# Tests for the split, enforced run_clean_tmle workflow:
#   authorize_outcome_analysis() -> hash-and-fingerprint-bound token,
#   run_clean_tmle_preoutcome()  -> pre-outcome dossier (Y not authorised),
#   run_clean_tmle_primary()     -> refuses to read Y without a valid token.
#
# Guard tests build the `pre` fixture with fit_ps_glm and no plasmode so they
# stay fast; the orchestration test exercises run_clean_tmle_preoutcome().

# Fast fixture: a `pre`-shaped object with a recorded audit and a gate token,
# built without SuperLearner or plasmode. When stop_gate = TRUE the gate is a
# STOP (Check Point 3 required but never recorded).
.mk_pre <- function(stop_gate = FALSE, seed = 7L, n = 400L) {
  d   <- sim_func1(n = n, seed = seed)
  cov <- c("age", "sex", "biomarker")
  lock  <- create_analysis_lock(d, "treatment", "event_24", cov,
                                sl_library = c("SL.glm", "SL.mean"), seed = seed)
  audit <- create_audit_log(lock)
  audit <- record_checkpoint(audit, checkpoint_cohort_adequacy(lock))
  ps    <- fit_ps_glm(lock)
  audit <- record_checkpoint(
    audit, checkpoint_balance(compute_ps_diagnostics(ps),
                              lock_hash = lock$lock_hash))
  req <- c("Check Point 1", "Check Point 2", "Check Point 3")
  if (!stop_gate) {
    lk2 <- define_negative_control(lock, "nc_outcome")
    nc  <- run_negative_control(lk2, "nc_outcome", ps)
    audit <- record_checkpoint(
      audit, checkpoint_residual_bias(nc, lock_hash = lock$lock_hash))
  }
  gate <- authorize_outcome_analysis(audit, required_stages = req,
                                     allow_flag = TRUE)
  structure(list(lock = lock, audit = audit, gate = gate, ps_fit = ps,
                 required_stages = req),
            class = "clean_tmle_preoutcome")
}

test_that("the gate token carries an audit fingerprint, lock hash, and verdict", {
  pre <- .mk_pre()
  g   <- pre$gate
  expect_s3_class(g, "pre_outcome_gate")
  expect_false(is.null(g$audit_fingerprint))
  expect_false(is.na(g$audit_fingerprint))
  expect_false(is.na(g$lock_hash))
  expect_equal(g$decision, "GO")
  expect_true(isTRUE(g$authorized))
})

test_that("the audit fingerprint changes when the audit is modified (tamper-evident)", {
  pre <- .mk_pre()
  fp1 <- pre$gate$audit_fingerprint
  a2  <- record_checkpoint(pre$audit, checkpoint_cohort_adequacy(pre$lock))
  g2  <- authorize_outcome_analysis(a2, required_stages = pre$required_stages,
                                    allow_flag = TRUE)
  expect_false(identical(fp1, g2$audit_fingerprint))
})

test_that("run_clean_tmle_primary refuses when no authorization is available", {
  pre <- .mk_pre()
  pre$gate <- NULL
  expect_error(run_clean_tmle_primary(pre, authorization = NULL),
               "authoriz", ignore.case = TRUE)
})

test_that("run_clean_tmle_primary refuses a STOP gate", {
  pre <- .mk_pre(stop_gate = TRUE)
  expect_error(run_clean_tmle_primary(pre, verbose = FALSE),
               "STOP|not authoris", ignore.case = TRUE)
})

test_that("run_clean_tmle_primary refuses an authorization minted for a different lock", {
  pre <- .mk_pre()
  bad <- pre$gate
  bad$lock_hash <- "deadbeefdeadbeef"
  expect_error(run_clean_tmle_primary(pre, authorization = bad, verbose = FALSE),
               "match", ignore.case = TRUE)
})

test_that("run_clean_tmle_primary refuses when the audit changed after authorization", {
  pre <- .mk_pre()
  # Tamper: append a checkpoint AFTER the token was minted from the audit.
  pre$audit <- record_checkpoint(pre$audit, checkpoint_cohort_adequacy(pre$lock))
  expect_error(run_clean_tmle_primary(pre, verbose = FALSE),
               "audit|fingerprint|changed", ignore.case = TRUE)
})

test_that("run_clean_tmle_primary runs the primary TMLE under a valid GO token", {
  pre <- .mk_pre()
  fit <- run_clean_tmle_primary(pre, verbose = FALSE)
  expect_true(is.numeric(fit$risk_difference))
  expect_true(isTRUE(fit$lock$.outcome_authorized))
})

test_that("run_clean_tmle_primary honors the allow_outcome_access escape hatch", {
  pre <- .mk_pre(stop_gate = TRUE)
  fit <- suppressWarnings(
    run_clean_tmle_primary(pre, allow_outcome_access = TRUE, verbose = FALSE))
  expect_true(is.numeric(fit$risk_difference))
})

test_that("a token minted from an empty audit (NA fingerprint) is refused (fail-closed)", {
  pre <- .mk_pre()
  # Forge a GO token from an EMPTY audit whose lock_hash still matches pre$lock.
  forged <- authorize_outcome_analysis(create_audit_log(pre$lock),
                                       required_stages = character(0))
  expect_true(isTRUE(forged$authorized))         # looks authorised
  expect_true(is.na(forged$audit_fingerprint))   # but carries no audit fingerprint
  expect_error(
    run_clean_tmle_primary(pre, authorization = forged, verbose = FALSE),
    "fingerprint|not bound|audit", ignore.case = TRUE)
})

test_that("run_clean_tmle_preoutcome builds a dossier without authorizing the outcome", {
  d <- sim_func1(n = 400, seed = 3)
  pre <- run_clean_tmle_preoutcome(
    data = d, Avar = "treatment", Yvar = "event_24",
    covariates = c("age", "sex", "biomarker"),
    learner_lib = c("SL.glm", "SL.mean"),
    ps_method = "glm", run_selection = FALSE, seed = 3, verbose = FALSE)
  expect_s3_class(pre, "clean_tmle_preoutcome")
  expect_s3_class(pre$gate, "pre_outcome_gate")
  expect_false(isTRUE(pre$lock$.outcome_authorized))
  expect_s3_class(pre$dossier, "clean_tmle_dossier")
  expect_true(all(c("estimand", "balance", "decision") %in% names(pre$dossier)))
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

test_that("preoutcome -> authorize -> primary is a working two-pass flow", {
  d   <- sim_func1(n = 400, seed = 3)
  pre <- run_clean_tmle_preoutcome(
    data = d, Avar = "treatment", Yvar = "event_24",
    covariates = c("age", "sex", "biomarker"),
    learner_lib = c("SL.glm", "SL.mean"),
    ps_method = "glm", run_selection = FALSE, seed = 3, verbose = FALSE)
  skip_if(pre$gate$decision == "STOP", "pre-outcome gate did not authorise")
  auth <- authorize_outcome_analysis(pre$audit,
                                     required_stages = pre$required_stages,
                                     allow_flag = TRUE)
  fit <- run_clean_tmle_primary(pre, authorization = auth, verbose = FALSE)
  expect_true(is.numeric(fit$risk_difference))
})

test_that("run_clean_tmle_preoutcome binds the selected candidate to the lock", {
  skip_on_cran()
  d     <- sim_func1(n = 400, seed = 3)
  cands <- list(tmle_candidate("t01", truncation = 0.01),
                tmle_candidate("t20", truncation = 0.20))
  pre <- run_clean_tmle_preoutcome(
    data = d, Avar = "treatment", Yvar = "event_24",
    covariates = c("age", "sex", "biomarker"),
    learner_lib = c("SL.glm", "SL.mean"), tmle_candidates = cands,
    run_selection = TRUE, plasmode_reps = 5L, ps_method = "glm",
    seed = 3, verbose = FALSE)
  expect_false(is.null(pre$lock$primary_tmle_spec))
  expect_equal(pre$lock$primary_tmle_spec$candidate_id, pre$selected$candidate_id)
})

test_that("run_clean_tmle_preoutcome DQ branch accepts dq_scenarios (arg name matches)", {
  skip_on_cran()
  d     <- sim_func1(n = 400, seed = 4)
  cands <- list(tmle_candidate("t01", truncation = 0.01))
  pre <- run_clean_tmle_preoutcome(
    data = d, Avar = "treatment", Yvar = "event_24",
    covariates = c("age", "sex", "biomarker"),
    learner_lib = c("SL.glm", "SL.mean"), tmle_candidates = cands,
    run_selection = TRUE, plasmode_reps = 5L,
    dq_scenarios = list(covariate_missingness = list(fractions = c(0.10))),
    ps_method = "glm", seed = 4, verbose = FALSE)
  expect_false(is.null(pre$dq_results))
})
