# What remains of the 0.1.x decision-threshold and gate layer after its
# removal in 0.3.0: the DQ stress scenarios themselves (now reached
# through stress_test()) and the DGP fidelity check. The
# decision_thresholds object, gate_dq, authorize_outcome_analysis, and
# the checkpoint constructors were deleted; prespecified thresholds now
# live on the lock (dq_thresholds, nc_criteria) and the locked verdict
# is computed inside the stress test.

# ── stress_test carries the MAR covariate-missingness scenario ───────────
test_that("stress_test includes the MAR covariate-missingness scenario", {
  skip_on_cran()
  dat  <- sim_func1(n = 300, seed = 5)
  lock <- create_analysis_lock(dat, "treatment", "event_24",
                               c("age", "sex", "biomarker"), seed = 5)
  cand <- define_candidates("glm", g_library = "SL.glm",
                            q_library = "SL.glm", truncation = 0.01)
  dq <- stress_test(
    lock, cand, reps = 2, verbose = FALSE,
    threats = list(
      covariate_missingness     = list(fractions = 0.10),
      covariate_missingness_mar = list(fractions = 0.10, treatment_OR = 3)))
  expect_s3_class(dq, "ct_stress")
  expect_true("cov_miss_mar" %in% unique(dq$metrics$scenario))
  expect_true("cov_miss"     %in% unique(dq$metrics$scenario))
  expect_identical(as.data.frame(dq), dq$metrics)
})

# ── assess_dgp_fidelity ──────────────────────────────────────────────────
test_that("assess_dgp_fidelity flags a shifted simulated covariate", {
  withr::local_seed(11)
  real  <- data.frame(age = rnorm(500), sex = rbinom(500, 1, 0.5),
                      treatment = rbinom(500, 1, 0.5))
  synth_same  <- real[sample.int(500, replace = TRUE), ]
  synth_shift <- real; synth_shift$age <- synth_shift$age + 2
  f_ok <- assess_dgp_fidelity(real, synth_same, c("age", "sex"),
                              "treatment")
  f_bad <- assess_dgp_fidelity(real, synth_shift, c("age", "sex"),
                               "treatment")
  expect_equal(f_ok$decision, "GO")
  expect_equal(f_bad$decision, "FLAG")
  expect_s3_class(f_bad, "cleantmle_dgp_fidelity")
})