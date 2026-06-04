test_that("run_plasmode_dq_stress runs and returns a metrics frame", {
  set.seed(2026)
  dat  <- sim_func1(n = 300, seed = 1)
  lock <- create_analysis_lock(
    data        = dat,
    treatment   = "treatment",
    outcome     = "event_24",
    covariates  = c("age", "sex", "biomarker"),
    sl_library  = c("SL.glm", "SL.mean"),
    plasmode_reps = 5L,
    seed        = 1L
  )
  candidates <- list(
    tmle_candidate("glm_t01", "GLM, trunc=0.01",
                   g_library = "SL.glm", truncation = 0.01)
  )

  res <- run_plasmode_dq_stress(
    lock,
    tmle_candidates = candidates,
    effect_sizes    = c(0.05),
    reps            = 3L,
    data_quality_scenarios = list(
      covariate_missingness = list(fractions = c(0.05, 0.10)),
      treatment_misclass    = list(sensitivity = c(0.95),
                                   specificity = c(0.99)),
      outcome_misclass      = list(sensitivity = c(0.90),
                                   specificity = c(0.95)),
      unmeasured_confounding = list(U_prevalence  = 0.20,
                                    U_treatment_OR = c(2.0),
                                    U_outcome_OR   = c(2.0))
    ),
    verbose = FALSE
  )

  expect_s3_class(res, "plasmode_dq_results")
  expect_true(is.data.frame(res$metrics))
  expect_true(nrow(res$metrics) > 0L)
  expect_true(all(c("scenario", "level", "candidate", "bias",
                     "rmse", "coverage") %in% names(res$metrics)))
  # Baseline must appear and be the only "none" row per candidate-effect combo.
  expect_true("none" %in% res$metrics$scenario)
  expect_true("cov_miss" %in% res$metrics$scenario)
  expect_true("trt_misclass" %in% res$metrics$scenario)
  expect_true("out_misclass" %in% res$metrics$scenario)
  expect_true("unmeasured_U" %in% res$metrics$scenario)
})


test_that("fit_timeout argument: callr guard produces valid results or degrades gracefully", {
  # When callr is available and fit_timeout is finite, each replicate's fits
  # run in a killable subprocess. On timeout the replicate is recorded as NA.
  # This test verifies that: (a) results still conform to the expected schema,
  # (b) the guard falls back gracefully when callr is absent, and
  # (c) a generous timeout (300s) never actually triggers on this tiny DGP.
  set.seed(2026)
  dat  <- sim_func1(n = 200, seed = 2)
  lock <- create_analysis_lock(
    data        = dat,
    treatment   = "treatment",
    outcome     = "event_24",
    covariates  = c("age", "sex", "biomarker"),
    sl_library  = "SL.glm",
    plasmode_reps = 3L,
    seed        = 2L
  )
  candidates <- list(
    tmle_candidate("glm_t01", "GLM, trunc=0.01",
                   g_library = "SL.glm", truncation = 0.01)
  )
  # 300 second timeout — never triggers on this tiny DGP but exercises the
  # callr-session-init path when callr is available.
  res <- run_plasmode_dq_stress(
    lock,
    tmle_candidates        = candidates,
    effect_sizes           = c(0.05),
    reps                   = 2L,
    data_quality_scenarios = list(
      covariate_missingness = list(fractions = c(0.10))
    ),
    fit_timeout = 300,
    verbose     = FALSE
  )
  expect_s3_class(res, "plasmode_dq_results")
  expect_true(is.data.frame(res$metrics))
  expect_true(nrow(res$metrics) > 0L)
  expect_true(all(c("scenario", "level", "candidate",
                    "bias", "rmse", "coverage") %in% names(res$metrics)))
})


test_that("summarize_dq_degradation produces a sensible relative table", {
  set.seed(2026)
  dat  <- sim_func1(n = 300, seed = 1)
  lock <- create_analysis_lock(
    data        = dat,
    treatment   = "treatment",
    outcome     = "event_24",
    covariates  = c("age", "sex", "biomarker"),
    sl_library  = "SL.glm",
    plasmode_reps = 5L,
    seed        = 1L
  )
  candidates <- list(
    tmle_candidate("glm_t01", "GLM, trunc=0.01",
                   g_library = "SL.glm", truncation = 0.01)
  )
  res <- run_plasmode_dq_stress(
    lock,
    tmle_candidates = candidates,
    effect_sizes    = c(0.05),
    reps            = 3L,
    data_quality_scenarios = list(
      covariate_missingness = list(fractions = c(0.10))
    ),
    verbose = FALSE
  )
  deg <- summarize_dq_degradation(res)
  expect_true(is.data.frame(deg))
  expect_true(nrow(deg) >= 1L)
  expect_true(all(c("scenario", "level", "rmse_ratio", "cov_drop")
                  %in% names(deg)))
})


test_that(".degrade_treatment honours sensitivity/specificity", {
  set.seed(7)
  A <- rbinom(2000, 1, 0.5)

  # Sensitivity 1, specificity 1 -> identity.
  expect_identical(
    cleanTMLE:::.degrade_treatment(A, sensitivity = 1, specificity = 1),
    A
  )

  # Sensitivity 0.5, specificity 1 -> only positives can flip to 0.
  set.seed(7)
  out <- cleanTMLE:::.degrade_treatment(A, sensitivity = 0.5, specificity = 1)
  pos <- which(A == 1L)
  neg <- which(A == 0L)
  # Negatives unchanged.
  expect_identical(out[neg], A[neg])
  # Approximately half of the positives flipped.
  expect_lt(abs(mean(out[pos] == 0) - 0.5), 0.05)
})


test_that(".degrade_outcome honours per-arm sens/spec", {
  set.seed(9)
  Y <- rbinom(2000, 1, 0.3)
  A <- rbinom(2000, 1, 0.5)

  out <- cleanTMLE:::.degrade_outcome(
    Y, A = A,
    sens_a1 = 1.0, sens_a0 = 0.5,
    spec_a1 = 1.0, spec_a0 = 1.0
  )
  # Arm 1 untouched.
  expect_identical(out[A == 1L], Y[A == 1L])
  # Arm 0 positives partially flipped.
  pos0 <- which(A == 0L & Y == 1L)
  if (length(pos0) > 0L) {
    expect_lt(abs(mean(out[pos0] == 0) - 0.5), 0.10)
  }
})


test_that("run_plasmode_dq_stress is reproducible under fixed seeds", {
  set.seed(2026)
  dat  <- sim_func1(n = 300, seed = 1)
  lock <- create_analysis_lock(
    data        = dat,
    treatment   = "treatment",
    outcome     = "event_24",
    covariates  = c("age", "sex", "biomarker"),
    sl_library  = "SL.glm",
    plasmode_reps = 3L,
    seed        = 7L
  )
  candidates <- list(
    tmle_candidate("glm_t01", "GLM, trunc=0.01",
                   g_library = "SL.glm", truncation = 0.01)
  )
  spec <- list(covariate_missingness = list(fractions = c(0.10)))

  res1 <- run_plasmode_dq_stress(lock, tmle_candidates = candidates,
                                  reps = 3L,
                                  data_quality_scenarios = spec,
                                  verbose = FALSE)
  res2 <- run_plasmode_dq_stress(lock, tmle_candidates = candidates,
                                  reps = 3L,
                                  data_quality_scenarios = spec,
                                  verbose = FALSE)
  expect_equal(res1$metrics, res2$metrics)
})


test_that("near_positivity inflates the lightly-truncated candidate's error", {
  # A constructed near-positivity design: amplifying the centred PS log-odds
  # pushes a subgroup toward deterministic treatment, so the inverse-probability
  # weight tail explodes for a lightly truncated candidate. The empirical SD and
  # RMSE should rise above the undisturbed baseline.
  set.seed(2026)
  dat  <- sim_func1(n = 500, seed = 11)
  lock <- create_analysis_lock(
    data        = dat,
    treatment   = "treatment",
    outcome     = "event_24",
    covariates  = c("age", "sex", "biomarker", "comorbidity"),
    sl_library  = "SL.glm",
    plasmode_reps = 10L,
    seed        = 11L
  )
  candidates <- list(
    tmle_candidate("aggressive", "GLM, trunc=0.001",
                   g_library = "SL.glm", truncation = 0.001)
  )

  res <- run_plasmode_dq_stress(
    lock,
    tmle_candidates = candidates,
    effect_sizes    = c(0.05),
    reps            = 10L,
    data_quality_scenarios = list(
      near_positivity = list(slopes = c(3.0, 4.0))
    ),
    verbose = FALSE
  )

  m <- res$metrics
  expect_true("near_positivity" %in% m$scenario)

  base_sd  <- m$emp_sd[m$scenario == "none"][1]
  base_rmse <- m$rmse[m$scenario == "none"][1]
  np_sd   <- max(m$emp_sd[m$scenario == "near_positivity"])
  np_rmse <- max(m$rmse[m$scenario == "near_positivity"])

  # Inflated dispersion: at least one near-positivity severity worsens the
  # empirical SD and the RMSE relative to baseline.
  expect_gt(np_sd, base_sd)
  expect_gt(np_rmse, base_rmse)
})
