# Tests for cleanroom workflow functions

# ── sim_func1 now includes event_24 ───────────────────────────────────────

test_that("sim_func1 includes event_24 column", {
  dat <- sim_func1(n = 100, seed = 80)
  expect_true("event_24" %in% names(dat))
  expect_true(all(dat$event_24 %in% c(0L, 1L)))
  # event_24 = 1 iff primary event by 24 months
  expected <- as.integer(dat$event == 1 & dat$time <= 24)
  expect_identical(dat$event_24, expected)
})


# ── create_analysis_lock ──────────────────────────────────────────────────

test_that("create_analysis_lock returns a cleanroom_lock", {
  dat  <- sim_func1(n = 200, seed = 81)
  lock <- create_analysis_lock(
    data       = dat,
    treatment  = "treatment",
    outcome    = "event_24",
    covariates = c("age", "sex", "biomarker"),
    seed       = 81L
  )
  expect_s3_class(lock, "cleanroom_lock")
  expect_equal(lock$treatment, "treatment")
  expect_equal(lock$outcome, "event_24")
  expect_equal(lock$covariates, c("age", "sex", "biomarker"))
  expect_equal(lock$seed, 81L)
  expect_true(is.character(lock$lock_hash))
  expect_equal(nrow(lock$data), 200L)
})

test_that("create_analysis_lock errors on bad inputs", {
  dat <- sim_func1(n = 50, seed = 82)

  expect_error(create_analysis_lock(
    data = "not_a_df", treatment = "treatment",
    outcome = "event_24", covariates = "age"
  ), "data.frame")

  expect_error(create_analysis_lock(
    data = dat, treatment = "no_such_var",
    outcome = "event_24", covariates = "age"
  ), "not found")

  expect_error(create_analysis_lock(
    data = dat, treatment = "treatment",
    outcome = "event_24", covariates = c("age", "no_such_col")
  ), "not found")
})


# ── validate_analysis_lock ────────────────────────────────────────────────

test_that("validate_analysis_lock passes for valid lock", {
  dat  <- sim_func1(n = 100, seed = 83)
  lock <- create_analysis_lock(
    data = dat, treatment = "treatment",
    outcome = "event_24", covariates = c("age", "sex"),
    seed = 83L
  )
  expect_message(validate_analysis_lock(lock), "validated successfully")
})

test_that("validate_analysis_lock detects tampered hash", {
  dat  <- sim_func1(n = 100, seed = 84)
  lock <- create_analysis_lock(
    data = dat, treatment = "treatment",
    outcome = "event_24", covariates = c("age", "sex"),
    seed = 84L
  )
  lock$lock_hash <- "000000000"   # tamper
  expect_error(validate_analysis_lock(lock), "hash mismatch")
})

test_that("print.cleanroom_lock works", {
  dat  <- sim_func1(n = 50, seed = 85)
  lock <- create_analysis_lock(
    data = dat, treatment = "treatment",
    outcome = "event_24", covariates = c("age", "sex"),
    seed = 85L
  )
  expect_output(print(lock), "cleanTMLE Analysis Lock")
  expect_output(print(lock), "treatment")
  expect_output(print(lock), "event_24")
})


# ── fit_ps_superlearner ───────────────────────────────────────────────────

test_that("fit_ps_superlearner errors without SuperLearner package", {
  skip_if(requireNamespace("SuperLearner", quietly = TRUE),
          "SuperLearner installed, skipping absence test")

  dat  <- sim_func1(n = 50, seed = 86)
  lock <- create_analysis_lock(
    data = dat, treatment = "treatment",
    outcome = "event_24", covariates = c("age", "sex"),
    seed = 86L
  )
  expect_error(fit_ps_superlearner(lock), "SuperLearner")
})

test_that("fit_ps_superlearner works when SuperLearner is available", {
  skip_if_not_installed("SuperLearner")

  dat  <- sim_func1(n = 200, seed = 87)
  lock <- create_analysis_lock(
    data = dat, treatment = "treatment",
    outcome = "event_24", covariates = c("age", "sex", "biomarker"),
    sl_library = c("SL.glm", "SL.mean"),
    seed = 87L
  )
  ps_fit <- fit_ps_superlearner(lock)

  expect_s3_class(ps_fit, "ps_fit")
  expect_equal(length(ps_fit$ps), 200L)
  expect_true(all(ps_fit$ps > 0 & ps_fit$ps < 1))
  expect_output(print(ps_fit), "Propensity Score")
})


# ── compute_ps_diagnostics ────────────────────────────────────────────────

test_that("compute_ps_diagnostics works with mock ps_fit", {
  dat <- sim_func1(n = 200, seed = 88)

  # Create a mock ps_fit with GLM-estimated PS
  ps_fml <- stats::reformulate(c("age", "sex", "biomarker"),
                                response = "treatment")
  ps_mod <- stats::glm(ps_fml, data = dat, family = binomial())
  ps_hat <- as.numeric(predict(ps_mod, type = "response"))

  mock_ps_fit <- list(
    ps         = ps_hat,
    sl_fit     = ps_mod,
    treatment  = "treatment",
    covariates = c("age", "sex", "biomarker"),
    data       = dat
  )
  class(mock_ps_fit) <- c("ps_fit", "cr_result")

  diag <- compute_ps_diagnostics(mock_ps_fit)

  expect_s3_class(diag, "ps_diagnostics")
  expect_true(is.data.frame(diag$ess))
  expect_true(is.data.frame(diag$smds))
  expect_equal(nrow(diag$smds), 3L)    # 3 covariates
  expect_s3_class(diag$overlap_plot, "ggplot")

  expect_output(print(diag), "Diagnostics")
  p <- plot(diag)
  expect_s3_class(p, "ggplot")
})


# ── run_plasmode_feasibility ──────────────────────────────────────────────

test_that("run_plasmode_feasibility returns plasmode_results with candidates", {
  dat  <- sim_func1(n = 300, seed = 89)
  lock <- create_analysis_lock(
    data = dat, treatment = "treatment",
    outcome = "event_24", covariates = c("age", "sex", "biomarker"),
    plasmode_reps = 5L, seed = 89L
  )

  candidates <- list(
    tmle_candidate("glm_t01", "GLM 0.01", g_library = "SL.glm", truncation = 0.01),
    tmle_candidate("glm_t05", "GLM 0.05", g_library = "SL.glm", truncation = 0.05)
  )

  sim_res <- run_plasmode_feasibility(
    lock            = lock,
    tmle_candidates = candidates,
    effect_sizes    = c(0.05),
    reps            = 5L
  )

  expect_s3_class(sim_res, "plasmode_results")
  expect_true(is.data.frame(sim_res$metrics))
  expect_true(all(c("effect_size", "candidate", "bias", "rmse", "coverage") %in%
                    names(sim_res$metrics)))
  # 2 candidates x 1 effect size
  expect_equal(nrow(sim_res$metrics), 2L)
  expect_output(print(sim_res), "Plasmode")
})

test_that("run_plasmode_feasibility with multiple effect sizes", {
  dat  <- sim_func1(n = 200, seed = 90)
  lock <- create_analysis_lock(
    data = dat, treatment = "treatment",
    outcome = "event_24", covariates = c("age", "sex"),
    plasmode_reps = 3L, seed = 90L
  )

  candidates <- list(
    tmle_candidate("glm_t01", g_library = "SL.glm", truncation = 0.01),
    tmle_candidate("glm_t05", g_library = "SL.glm", truncation = 0.05)
  )

  sim_res <- run_plasmode_feasibility(
    lock            = lock,
    tmle_candidates = candidates,
    effect_sizes    = c(0.03, 0.07),
    reps            = 3L
  )

  expect_equal(nrow(sim_res$metrics), 4L)   # 2 effects x 2 candidates
})


# ── select_tmle_candidate ─────────────────────────────────────────────────

test_that("select_tmle_candidate works with all rules", {
  # Build a mock plasmode_results with candidate specs
  cand_a <- tmle_candidate("cand_a", g_library = "SL.glm", truncation = 0.01)
  cand_b <- tmle_candidate("cand_b", g_library = "SL.glm", truncation = 0.05)

  mock_metrics <- data.frame(
    effect_size = c(0.05, 0.05, 0.10, 0.10),
    candidate   = c("cand_a", "cand_b", "cand_a", "cand_b"),
    bias        = c(0.02, 0.01, 0.03, 0.01),
    rmse        = c(0.05, 0.04, 0.06, 0.03),
    coverage    = c(0.90, 0.93, 0.88, 0.94),
    emp_sd      = c(0.04, 0.03, 0.05, 0.03),
    mean_se     = c(0.04, 0.031, 0.05, 0.031),
    stringsAsFactors = FALSE
  )
  mock_res <- list(metrics = mock_metrics, effect_sizes = c(0.05, 0.10),
                   reps = 5L, tmle_candidates = list(cand_a, cand_b))
  class(mock_res) <- "plasmode_results"

  best_rmse     <- select_tmle_candidate(mock_res, rule = "min_rmse")
  best_bias     <- select_tmle_candidate(mock_res, rule = "min_bias")
  best_coverage <- select_tmle_candidate(mock_res, rule = "max_coverage")

  expect_s3_class(best_rmse, "tmle_selected_spec")
  expect_s3_class(best_bias, "tmle_selected_spec")
  expect_s3_class(best_coverage, "tmle_selected_spec")

  # cand_b has lower rmse, lower bias, higher coverage
 expect_equal(best_rmse$candidate_id, "cand_b")
  expect_equal(best_bias$candidate_id, "cand_b")
  expect_equal(best_coverage$candidate_id, "cand_b")
})


test_that("select_tmle_candidate min_max_rmse picks robust candidate", {
  cand_a <- tmle_candidate("cand_a", g_library = "SL.glm", truncation = 0.01)
  cand_b <- tmle_candidate("cand_b", g_library = "SL.glm", truncation = 0.05)

  # Baseline: cand_a slightly better RMSE
  mock_metrics <- data.frame(
    effect_size = c(0.05, 0.05),
    candidate   = c("cand_a", "cand_b"),
    bias        = c(0.01, 0.012),
    rmse        = c(0.040, 0.042),
    coverage    = c(0.94, 0.94),
    emp_sd      = c(0.04, 0.04),
    mean_se     = c(0.04, 0.04),
    stringsAsFactors = FALSE
  )
  mock_res <- list(metrics = mock_metrics, effect_sizes = 0.05,
                   reps = 5L, tmle_candidates = list(cand_a, cand_b))
  class(mock_res) <- "plasmode_results"

  # DQ: cand_a is fragile (RMSE blows up under one scenario), cand_b is robust
  dq_metrics <- data.frame(
    scenario    = c("none", "none",
                    "covariate_missingness", "covariate_missingness",
                    "outcome_misclass", "outcome_misclass"),
    level       = c("baseline", "baseline", "frac=0.20", "frac=0.20",
                    "sens=0.90,spec=0.95", "sens=0.90,spec=0.95"),
    effect_size = rep(0.05, 6),
    candidate   = rep(c("cand_a", "cand_b"), 3),
    bias        = c(0.01, 0.012, 0.03, 0.013, 0.04, 0.014),
    rmse        = c(0.04, 0.042, 0.20, 0.05, 0.25, 0.06),
    coverage    = rep(0.94, 6),
    emp_sd      = rep(0.04, 6),
    mean_se     = rep(0.04, 6),
    se_cal      = rep(1, 6),
    n_converged = rep(5, 6),
    stringsAsFactors = FALSE
  )
  dq_res <- list(metrics = dq_metrics,
                 scenarios = list(),
                 baseline = dq_metrics[dq_metrics$scenario == "none", ],
                 lock = NULL, reps = 5L, call = sys.call())
  class(dq_res) <- "plasmode_dq_results"

  expect_error(
    select_tmle_candidate(mock_res, rule = "min_max_rmse"),
    "dq_results.*required"
  )

  best <- select_tmle_candidate(mock_res, rule = "min_max_rmse",
                                dq_results = dq_res)
  expect_s3_class(best, "tmle_selected_spec")
  expect_equal(best$candidate_id, "cand_b")
})


# ── run_match_workflow ────────────────────────────────────────────────────

test_that("run_match_workflow returns match_result", {
  dat  <- sim_func1(n = 300, seed = 91)
  lock <- create_analysis_lock(
    data = dat, treatment = "treatment",
    outcome = "event_24", covariates = c("age", "sex", "biomarker"),
    seed = 91L
  )

  ps_fml <- stats::reformulate(c("age", "sex", "biomarker"),
                                response = "treatment")
  ps_mod <- stats::glm(ps_fml, data = dat, family = binomial())
  ps_hat <- as.numeric(predict(ps_mod, type = "response"))
  mock_ps <- list(ps = ps_hat, sl_fit = ps_mod, treatment = "treatment",
                  covariates = c("age", "sex", "biomarker"), data = dat)
  class(mock_ps) <- c("ps_fit", "cr_result")

  m_fit <- run_match_workflow(lock, mock_ps)

  expect_s3_class(m_fit, "match_result")
  expect_true(!is.na(m_fit$estimate))
  expect_true(m_fit$n_matched > 0L)
  expect_true(!is.na(m_fit$se))
  expect_output(print(m_fit), "Matching")
})


# ── run_iptw_workflow ─────────────────────────────────────────────────────

test_that("run_iptw_workflow returns iptw_result", {
  dat  <- sim_func1(n = 300, seed = 92)
  lock <- create_analysis_lock(
    data = dat, treatment = "treatment",
    outcome = "event_24", covariates = c("age", "sex", "biomarker"),
    seed = 92L
  )

  ps_fml <- stats::reformulate(c("age", "sex", "biomarker"),
                                response = "treatment")
  ps_mod <- stats::glm(ps_fml, data = dat, family = binomial())
  ps_hat <- as.numeric(predict(ps_mod, type = "response"))
  mock_ps <- list(ps = ps_hat, sl_fit = ps_mod, treatment = "treatment",
                  covariates = c("age", "sex", "biomarker"), data = dat)
  class(mock_ps) <- c("ps_fit", "cr_result")

  i_fit <- run_iptw_workflow(lock, mock_ps)

  expect_s3_class(i_fit, "iptw_result")
  expect_true(!is.na(i_fit$estimate))
  expect_true(!is.na(i_fit$se))
  expect_true(!is.na(i_fit$ci_lower))
  expect_output(print(i_fit), "IPTW")
})

test_that("run_iptw_workflow with trimming", {
  dat  <- sim_func1(n = 200, seed = 93)
  lock <- create_analysis_lock(
    data = dat, treatment = "treatment",
    outcome = "event_24", covariates = c("age", "sex"),
    seed = 93L
  )
  ps_fml <- stats::reformulate(c("age", "sex"), response = "treatment")
  ps_mod <- stats::glm(ps_fml, data = dat, family = binomial())
  ps_hat <- as.numeric(predict(ps_mod, type = "response"))
  mock_ps <- list(ps = ps_hat, sl_fit = ps_mod, treatment = "treatment",
                  covariates = c("age", "sex"), data = dat)
  class(mock_ps) <- c("ps_fit", "cr_result")

  i_fit <- run_iptw_workflow(lock, mock_ps, trim = 0.01)
  expect_s3_class(i_fit, "iptw_result")
})


# ── Modular TMLE ──────────────────────────────────────────────────────────

test_that("fit_tmle_treatment_mechanism uses existing ps_fit", {
  dat  <- sim_func1(n = 200, seed = 94)
  lock <- create_analysis_lock(
    data = dat, treatment = "treatment",
    outcome = "event_24", covariates = c("age", "sex"),
    seed = 94L
  )

  ps_fml <- stats::reformulate(c("age", "sex"), response = "treatment")
  ps_mod <- stats::glm(ps_fml, data = dat, family = binomial())
  ps_hat <- as.numeric(predict(ps_mod, type = "response"))
  mock_ps <- list(ps = ps_hat, sl_fit = ps_mod, treatment = "treatment",
                  covariates = c("age", "sex"), data = dat)
  class(mock_ps) <- c("ps_fit", "cr_result")

  g_fit <- fit_tmle_treatment_mechanism(lock, ps_fit = mock_ps)

  expect_s3_class(g_fit, "tmle_mechanism")
  expect_equal(g_fit$type, "treatment")
  expect_equal(length(g_fit$ps), 200L)
})

test_that("fit_tmle_outcome_mechanism returns outcome mechanism", {
  dat  <- sim_func1(n = 200, seed = 95)
  lock <- create_analysis_lock(
    data = dat, treatment = "treatment",
    outcome = "event_24", covariates = c("age", "sex"),
    seed = 95L
  )

  ps_fml <- stats::reformulate(c("age", "sex"), response = "treatment")
  ps_mod <- stats::glm(ps_fml, data = dat, family = binomial())
  ps_hat <- as.numeric(predict(ps_mod, type = "response"))
  mock_ps <- list(ps = ps_hat, sl_fit = ps_mod, treatment = "treatment",
                  covariates = c("age", "sex"), data = dat)
  class(mock_ps) <- c("ps_fit", "cr_result")

  g_fit <- fit_tmle_treatment_mechanism(lock, ps_fit = mock_ps)
  Q_fit <- fit_tmle_outcome_mechanism(lock, g_fit)

  expect_s3_class(Q_fit, "tmle_mechanism")
  expect_equal(Q_fit$type, "outcome")
  expect_equal(length(Q_fit$Q_a1), 200L)
  expect_equal(length(Q_fit$Q_a0), 200L)
  expect_true(all(Q_fit$Q_a1 > 0 & Q_fit$Q_a1 < 1))
})

test_that("run_tmle_targeting_step and extract_tmle_estimate work", {
  dat  <- sim_func1(n = 300, seed = 96)
  lock <- create_analysis_lock(
    data = dat, treatment = "treatment",
    outcome = "event_24", covariates = c("age", "sex", "biomarker"),
    seed = 96L
  )

  ps_fml <- stats::reformulate(c("age", "sex", "biomarker"),
                                response = "treatment")
  ps_mod <- stats::glm(ps_fml, data = dat, family = binomial())
  ps_hat <- as.numeric(predict(ps_mod, type = "response"))
  mock_ps <- list(ps = ps_hat, sl_fit = ps_mod, treatment = "treatment",
                  covariates = c("age", "sex", "biomarker"), data = dat)
  class(mock_ps) <- c("ps_fit", "cr_result")

  g_fit    <- fit_tmle_treatment_mechanism(lock, ps_fit = mock_ps)
  Q_fit    <- fit_tmle_outcome_mechanism(lock, g_fit)
  tmle_upd <- run_tmle_targeting_step(g_fit, Q_fit)

  expect_s3_class(tmle_upd, "tmle_update")
  expect_true(!is.na(tmle_upd$psi))
  expect_equal(length(tmle_upd$eic), 300L)

  tmle_est <- extract_tmle_estimate(tmle_upd)

  expect_s3_class(tmle_est, "tmle_fit")
  expect_true("ATE" %in% names(tmle_est$estimates))
  expect_true(!is.na(tmle_est$estimates$ATE$estimate))
  expect_true(!is.na(tmle_est$estimates$ATE$se))
  expect_output(print(tmle_est), "TMLE")
})


# ── summarize_cleanroom_results ───────────────────────────────────────────

test_that("summarize_cleanroom_results creates comparison table", {
  # Mock objects
  mock_match <- list(
    estimate = 0.04, se = 0.02, ci_lower = 0.00, ci_upper = 0.08,
    p_value = 0.05
  )
  class(mock_match) <- c("match_result", "cr_result")

  mock_iptw <- list(
    estimate = 0.05, se = 0.02, ci_lower = 0.01, ci_upper = 0.09,
    p_value = 0.01
  )
  class(mock_iptw) <- c("iptw_result", "cr_result")

  mock_tmle <- list(
    estimates = list(
      ATE = list(estimate = 0.06, se = 0.02, ci_lower = 0.02,
                 ci_upper = 0.10, p_value = 0.003)
    ),
    type = "modular_tmle"
  )
  class(mock_tmle) <- c("tmle_fit", "cr_result")

  tbl <- summarize_cleanroom_results(list(mock_match, mock_iptw, mock_tmle))

  expect_true(is.data.frame(tbl))
  expect_equal(nrow(tbl), 3L)
  expect_true(all(c("method", "estimate", "se", "ci_lower", "ci_upper",
                     "p_value") %in% names(tbl)))
  expect_equal(tbl$method, c("PS Matching", "IPTW", "TMLE"))
})

test_that("summarize_cleanroom_results uses list names when provided", {
  mock_match <- list(
    estimate = 0.04, se = 0.02, ci_lower = 0.00, ci_upper = 0.08,
    p_value = 0.05
  )
  class(mock_match) <- c("match_result", "cr_result")

  mock_iptw <- list(
    estimate = 0.05, se = 0.02, ci_lower = 0.01, ci_upper = 0.09,
    p_value = 0.01
  )
  class(mock_iptw) <- c("iptw_result", "cr_result")

  tbl <- summarize_cleanroom_results(list(
    "My Match" = mock_match, "My IPTW" = mock_iptw
  ))
  expect_equal(tbl$method, c("My Match", "My IPTW"))
})

test_that("summarize_cleanroom_results warns on unsupported class", {
  mock_match <- list(
    estimate = 0.04, se = 0.02, ci_lower = 0.00, ci_upper = 0.08,
    p_value = 0.05
  )
  class(mock_match) <- c("match_result", "cr_result")

  bad_obj <- list(x = 1)
  class(bad_obj) <- "unknown_class"

  expect_warning(
    tbl <- summarize_cleanroom_results(list(mock_match, bad_obj)),
    "Unsupported"
  )
  expect_equal(nrow(tbl), 1L)
})


# ── fit_ps_glm ────────────────────────────────────────────────────────────

test_that("fit_ps_glm returns a ps_fit with method = 'glm'", {
  dat  <- sim_func1(n = 200, seed = 91)
  lock <- create_analysis_lock(
    data = dat, treatment = "treatment",
    outcome = "event_24", covariates = c("age", "sex", "biomarker"),
    seed = 91L
  )
  ps <- fit_ps_glm(lock)

  expect_s3_class(ps, "ps_fit")
  expect_equal(ps$method, "glm")
  expect_equal(length(ps$ps), 200L)
  expect_true(all(ps$ps > 0 & ps$ps < 1))
  expect_true(all(ps$ps >= 0.01 & ps$ps <= 0.99))
})


# ── summarize_plasmode_results ────────────────────────────────────────────

test_that("summarize_plasmode_results prints and returns invisibly", {
  dat  <- sim_func1(n = 200, seed = 92)
  lock <- create_analysis_lock(
    data = dat, treatment = "treatment",
    outcome = "event_24", covariates = c("age", "sex"),
    plasmode_reps = 5L, seed = 92L
  )
  sim_res <- run_plasmode_feasibility(lock, effect_sizes = c(0.05), reps = 5L)
  expect_output(
    result <- summarize_plasmode_results(sim_res),
    "Plasmode"
  )
  expect_s3_class(result, "plasmode_results")
})


# ── fit_final_workflows ───────────────────────────────────────────────────

test_that("fit_final_workflows runs all three workflows", {
  dat  <- sim_func1(n = 300, seed = 93)
  lock <- create_analysis_lock(
    data = dat, treatment = "treatment",
    outcome = "event_24", covariates = c("age", "sex", "biomarker"),
    seed = 93L
  )
  ps_fit <- fit_ps_glm(lock)
  result <- fit_final_workflows(lock, ps_fit)

  expect_true(is.list(result))
  expect_true("match" %in% names(result))
  expect_true("iptw"  %in% names(result))
  expect_true("tmle"  %in% names(result))
  expect_s3_class(result$match, "match_result")
  expect_s3_class(result$iptw,  "iptw_result")
  expect_s3_class(result$tmle,  "tmle_fit")
})


# ── fit_tmle_candidate_set ────────────────────────────────────────────────

test_that("fit_tmle_candidate_set returns results with candidate specs", {
  dat  <- sim_func1(n = 300, seed = 94)
  lock <- create_analysis_lock(
    data = dat, treatment = "treatment",
    outcome = "event_24", covariates = c("age", "sex", "biomarker"),
    seed = 94L
  )

  candidates <- list(
    tmle_candidate("glm_t01", g_library = "SL.glm", truncation = 0.01),
    tmle_candidate("glm_t05", g_library = "SL.glm", truncation = 0.05)
  )
  cands <- fit_tmle_candidate_set(lock, candidates = candidates)

  expect_true(is.list(cands))
  expect_true(length(cands) > 0L)
  expect_true(all(vapply(cands, inherits, logical(1L), "tmle_fit")))
  expect_true(all(c("glm_t01", "glm_t05") %in% names(cands)))
})
