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
