# Regression tests for cleanTMLE 0.1.3:
#   decision_thresholds() + dt_* extractors + attach_decision_thresholds()
#   gate_dq() GO/FLAG/STOP
#   authorize_outcome_analysis() masking-bug fix + block_on_any_stop
#   checkpoint_residual_bias(rule = "equivalence")
#   checkpoint_balance(stop_smd) / checkpoint_cohort_adequacy(stop_*)
#   run_plasmode_dq_stress() MAR scenario
#   assess_dgp_fidelity()

# ── decision_thresholds + extractors ─────────────────────────────────────
test_that("decision_thresholds builds a classed object and extractors work", {
  dt <- decision_thresholds(dq_max_abs_bias = 0.03, nco_rule = "equivalence",
                            nco_null_band = 0.025)
  expect_s3_class(dt, "cleantmle_thresholds")
  expect_equal(dt_dq(dt)$max_abs_bias, 0.03)
  expect_equal(dt_nco(dt)$rule, "equivalence")
  expect_equal(dt_nco(dt)$null_band, 0.025)
  expect_true(all(c("max_abs_bias", "coverage_threshold", "se_sd_window") %in%
                    names(dt_plasmode(dt))))
  expect_equal(dt_balance(dt)$stop_smd, 0.20)
  expect_equal(dt_cohort(dt)$stop_n_per_arm, 20L)
})

test_that("attach_decision_thresholds fingerprints the rule and preserves lock_hash", {
  dat  <- sim_func1(n = 120, seed = 1)
  lock <- create_analysis_lock(dat, "treatment", "event_24",
                               c("age", "sex"), seed = 1)
  base_hash <- lock$lock_hash
  lock2 <- attach_decision_thresholds(lock, decision_thresholds())
  expect_false(is.null(lock2$thresholds_hash))
  expect_identical(lock2$lock_hash, base_hash)         # base hash unchanged
  expect_no_error(validate_analysis_lock(lock2))
  # Different thresholds => different fingerprint.
  lock3 <- attach_decision_thresholds(lock,
                                      decision_thresholds(dq_max_abs_bias = 0.05))
  expect_false(identical(lock2$thresholds_hash, lock3$thresholds_hash))
})

# ── gate_dq ──────────────────────────────────────────────────────────────
make_dq <- function(bias_deg = 0.03, rmse_deg = 0.020, cov_deg = 0.80) {
  metrics <- data.frame(
    scenario    = c("none", "unmeasured_U"),
    level       = c("0", "or2"),
    effect_size = c(0.05, 0.05),
    candidate   = c("c1", "c1"),
    bias        = c(0.001, bias_deg),
    rmse        = c(0.010, rmse_deg),
    coverage    = c(0.95, cov_deg),
    emp_sd      = c(0.010, 0.020),
    mean_se     = c(0.010, 0.020),
    se_cal      = c(1, 1),
    n_converged = c(30L, 30L),
    stringsAsFactors = FALSE)
  structure(list(metrics = metrics), class = "plasmode_dq_results")
}

test_that("gate_dq returns STOP when a degraded row breaches the bias envelope", {
  cp <- gate_dq(make_dq(bias_deg = 0.03), candidate = "c1",
                max_abs_bias = 0.02, min_coverage = 0.50, max_rmse_ratio = 5)
  expect_s3_class(cp, "cleantmle_checkpoint")
  expect_equal(cp$decision, "STOP")
  expect_equal(cp$stage, "Check Point 2c: DQ Stress")
})

test_that("gate_dq returns GO when all degraded rows are within the envelope", {
  cp <- gate_dq(make_dq(bias_deg = 0.005, rmse_deg = 0.011, cov_deg = 0.93),
                candidate = "c1", max_abs_bias = 0.02, min_coverage = 0.85,
                max_rmse_ratio = 1.5)
  expect_equal(cp$decision, "GO")
})

test_that("gate_dq reads thresholds from a decision_thresholds object", {
  dt <- decision_thresholds(dq_max_abs_bias = 0.02, dq_min_coverage = 0.85,
                            dq_max_rmse_ratio = 1.5)
  cp <- gate_dq(make_dq(bias_deg = 0.03), candidate = "c1", thresholds = dt)
  expect_equal(cp$decision, "STOP")
})

# ── authorize_outcome_analysis: masking bug + block_on_any_stop ───────────
test_that("a balance STOP is not masked by a later same-prefix DQ checkpoint", {
  dat  <- sim_func1(n = 120, seed = 2)
  lock <- create_analysis_lock(dat, "treatment", "event_24",
                               c("age", "sex"), seed = 2)
  audit <- create_audit_log(lock)
  cp1 <- new_checkpoint("Check Point 1: Cohort Adequacy", "GO",
                        data.frame(x = 1), list(), "ok")
  cp2 <- new_checkpoint("Check Point 2: Treatment Comparability", "STOP",
                        data.frame(x = 1), list(), "bad balance")
  cpdq <- new_checkpoint("Check Point 2c: DQ Stress", "GO",
                         data.frame(x = 1), list(), "dq ok")
  cp3 <- new_checkpoint("Check Point 3: Residual Bias", "GO",
                        data.frame(x = 1), list(), "nc ok")
  audit <- record_checkpoint(audit, cp1)
  audit <- record_checkpoint(audit, cp2)
  audit <- record_checkpoint(audit, cpdq)   # later, shares "Check Point 2" prefix
  audit <- record_checkpoint(audit, cp3)
  g <- authorize_outcome_analysis(audit)
  expect_equal(g$decision, "STOP")           # balance STOP must win
  expect_false(isTRUE(g$authorized))
})

test_that("a recorded DQ STOP blocks authorization via block_on_any_stop", {
  dat  <- sim_func1(n = 120, seed = 3)
  lock <- create_analysis_lock(dat, "treatment", "event_24",
                               c("age", "sex"), seed = 3)
  audit <- create_audit_log(lock)
  audit <- record_checkpoint(audit, new_checkpoint(
    "Check Point 1: Cohort Adequacy", "GO", data.frame(x=1), list(), ""))
  audit <- record_checkpoint(audit, new_checkpoint(
    "Check Point 2: Treatment Comparability", "GO", data.frame(x=1), list(), ""))
  audit <- record_checkpoint(audit, new_checkpoint(
    "Check Point 3: Residual Bias", "GO", data.frame(x=1), list(), ""))
  audit <- record_checkpoint(audit, new_checkpoint(
    "Check Point 2c: DQ Stress", "STOP", data.frame(x=1), list(), "dq breach"))
  g_block <- authorize_outcome_analysis(audit, block_on_any_stop = TRUE)
  expect_equal(g_block$decision, "STOP")
  # With the safety net off, the non-required CP2c STOP no longer blocks.
  g_noblock <- authorize_outcome_analysis(audit, block_on_any_stop = FALSE)
  expect_equal(g_noblock$decision, "GO")
})

# ── checkpoint_residual_bias equivalence rule ────────────────────────────
nc_obj <- function(est, se) structure(
  list(variable = "nc", estimate = est, se = se,
       ci_lower = est - 1.96 * se, ci_upper = est + 1.96 * se,
       p_value = 2 * stats::pnorm(-abs(est / se))),
  class = "cleantmle_nc_result")

test_that("equivalence NCO rule passes a precise near-null and flags an imprecise one", {
  cp_ok  <- checkpoint_residual_bias(nc_obj(0.001, 0.005),
                                     rule = "equivalence", null_band = 0.02)
  expect_equal(cp_ok$decision, "GO")
  cp_bad <- checkpoint_residual_bias(nc_obj(0.001, 0.020),
                                     rule = "equivalence", null_band = 0.02)
  expect_equal(cp_bad$decision, "STOP")          # wide CI cannot prove near-null
  expect_error(checkpoint_residual_bias(nc_obj(0, 0.01), rule = "equivalence"))
})

# ── checkpoint_balance stop_smd / cohort stop floors ─────────────────────
test_that("checkpoint_balance STOP floor is configurable via stop_smd", {
  ps_diag <- structure(list(
    smds = data.frame(variable = c("a", "b"),
                      smd_unweighted = c(0.30, 0.10),
                      smd_weighted   = c(0.25, 0.05)),
    ess  = data.frame(group = c("Treated", "Control", "Total"),
                      ess_pct = c(80, 80, 80))),
    class = "ps_diagnostics")
  expect_equal(checkpoint_balance(ps_diag, max_smd = 0.10,
                                  stop_smd = 0.20)$decision, "STOP")
  expect_equal(checkpoint_balance(ps_diag, max_smd = 0.10,
                                  stop_smd = 0.30)$decision, "FLAG")
})

test_that("checkpoint_cohort_adequacy STOP floor is configurable", {
  set.seed(10)
  df <- data.frame(
    treatment = c(rep(1L, 15), rep(0L, 105)),
    event_24  = rbinom(120, 1, 0.3),
    age = rnorm(120), sex = rbinom(120, 1, 0.5))
  lock <- create_analysis_lock(df, "treatment", "event_24",
                               c("age", "sex"), seed = 1)
  cp <- checkpoint_cohort_adequacy(lock, min_n_per_arm = 50,
                                   stop_n_per_arm = 20)
  expect_equal(cp$decision, "STOP")              # 15 treated < 20
  expect_equal(cp$thresholds$stop_n_per_arm, 20)
})

# ── run_plasmode_dq_stress MAR scenario ──────────────────────────────────
test_that("run_plasmode_dq_stress includes the MAR covariate-missingness scenario", {
  skip_on_cran()
  dat  <- sim_func1(n = 300, seed = 5)
  lock <- create_analysis_lock(dat, "treatment", "event_24",
                               c("age", "sex", "biomarker"), seed = 5)
  cand <- list(tmle_candidate("glm", "glm", g_library = "SL.glm",
                              q_library = "SL.glm", truncation = 0.01))
  dq <- run_plasmode_dq_stress(
    lock, tmle_candidates = cand, reps = 2, verbose = FALSE,
    data_quality_scenarios = list(
      covariate_missingness     = list(fractions = 0.10),
      covariate_missingness_mar = list(fractions = 0.10, treatment_OR = 3)))
  expect_true("cov_miss_mar" %in% unique(dq$metrics$scenario))
  expect_true("cov_miss"     %in% unique(dq$metrics$scenario))
})

# ── assess_dgp_fidelity ──────────────────────────────────────────────────
test_that("assess_dgp_fidelity flags a shifted synthetic covariate", {
  set.seed(11)
  real  <- data.frame(age = rnorm(500), sex = rbinom(500, 1, 0.5),
                      treatment = rbinom(500, 1, 0.5))
  synth_same  <- real[sample.int(500, replace = TRUE), ]
  synth_shift <- real; synth_shift$age <- synth_shift$age + 2
  f_ok <- assess_dgp_fidelity(real, synth_same, c("age", "sex"), "treatment")
  f_bad <- assess_dgp_fidelity(real, synth_shift, c("age", "sex"), "treatment")
  expect_equal(f_ok$decision, "GO")
  expect_equal(f_bad$decision, "FLAG")
  expect_s3_class(f_bad, "cleantmle_dgp_fidelity")
})
