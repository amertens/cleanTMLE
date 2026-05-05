# Regression tests for cleanTMLE 0.1.1 — outcome NA handling, preset
# scenarios, attrition_table polymorphism, make_table1 lock method,
# print_locked_spec, and authorize_outcome_analysis(checkpoints = ...).

test_that("run_plasmode_feasibility tolerates 25% outcome NA", {
  set.seed(1)
  dat <- sim_func1(n = 400, seed = 1)
  dat$event_24[sample(400, 100)] <- NA
  lock <- create_analysis_lock(dat, "treatment", "event_24",
                               c("age", "sex", "biomarker", "comorbidity"),
                               seed = 1L)
  res <- run_plasmode_feasibility(lock, reps = 2L, effect_sizes = 0.05,
                                   verbose = FALSE)
  expect_s3_class(res, "plasmode_results")
  expect_gt(nrow(res$metrics), 0L)
})

test_that("run_plasmode_dq_stress tolerates 25% outcome NA", {
  set.seed(2)
  dat <- sim_func1(n = 400, seed = 2)
  dat$event_24[sample(400, 100)] <- NA
  lock <- create_analysis_lock(dat, "treatment", "event_24",
                               c("age", "sex", "biomarker", "comorbidity"),
                               seed = 2L)
  res <- run_plasmode_dq_stress(lock, reps = 2L, effect_sizes = 0.05,
                                 data_quality_scenarios = default_dq_scenarios("exploratory"),
                                 verbose = FALSE)
  expect_s3_class(res, "plasmode_dq_results")
  expect_gt(nrow(res$metrics), 0L)
})

test_that("plasmode functions error clearly on a fully-masked lock", {
  set.seed(3)
  dat <- sim_func1(n = 200, seed = 3)
  lock <- create_analysis_lock(dat, "treatment", "event_24",
                               c("age", "sex", "biomarker"), seed = 3L)
  masked <- mask_outcome(lock)
  expect_error(run_plasmode_feasibility(masked, reps = 2L,
                                         effect_sizes = 0.05,
                                         verbose = FALSE),
               "outcome|masked|outcome-access|cannot be fit")
})

test_that("default_dq_scenarios returns valid configs for each preset", {
  for (p in c("regulatory_standard", "exploratory", "stress")) {
    cfg <- default_dq_scenarios(p)
    expect_named(cfg, c("covariate_missingness", "treatment_misclass",
                        "outcome_misclass", "unmeasured_confounding"))
    expect_true(!is.null(cfg$covariate_missingness$fractions))
    expect_true(!is.null(cfg$unmeasured_confounding$U_treatment_OR))
  }
})

test_that("attrition_table accepts named-list and named-numeric inputs", {
  out_list <- attrition_table(list(All = 1000, Eligible = 800, Final = 600))
  expect_s3_class(out_list, "cleantmle_attrition")
  expect_equal(out_list$n_remaining, c(1000, 800, 600))
  expect_equal(out_list$n_excluded,  c(0L,   200L, 200L))

  out_num <- attrition_table(c(All = 1000, Eligible = 800, Final = 600))
  expect_equal(out_num$n_remaining, out_list$n_remaining)
})

test_that("make_table1 accepts a cleanroom_lock", {
  dat  <- sim_func1(n = 200, seed = 5)
  lock <- create_analysis_lock(dat, "treatment", "event_24",
                               c("age", "sex", "biomarker"), seed = 5L)
  tbl <- make_table1(lock)
  expect_true(is.data.frame(tbl) || is.list(tbl))
})

test_that("print_locked_spec returns NULL when nothing locked", {
  dat  <- sim_func1(n = 200, seed = 6)
  lock <- create_analysis_lock(dat, "treatment", "event_24",
                               c("age", "sex", "biomarker"), seed = 6L)
  out <- capture.output(spec <- print_locked_spec(lock))
  expect_null(spec)
})

test_that("authorize_outcome_analysis accepts checkpoints =", {
  mk <- function(s, d) new_checkpoint(stage = s, decision = d,
                                       metrics = data.frame(),
                                       thresholds = list())
  cp1 <- mk("Check Point 1", "GO")
  cp2 <- mk("Check Point 2", "FLAG")
  cp3 <- mk("Check Point 3", "GO")
  gate <- authorize_outcome_analysis(checkpoints = list(cp1, cp2, cp3),
                                      allow_flag = TRUE)
  expect_true(isTRUE(gate$authorized))
})

test_that("gate_check accepts ergonomic short form on plasmode_results", {
  cand_a <- tmle_candidate("a", g_library = "SL.glm", truncation = 0.01)
  cand_b <- tmle_candidate("b", g_library = "SL.glm", truncation = 0.05)
  m <- data.frame(
    effect_size = c(0.05, 0.05),
    candidate   = c("a", "b"),
    bias        = c(0.001, 0.002),
    rmse        = c(0.04,  0.05),
    coverage    = c(0.94,  0.93),
    emp_sd      = c(0.04,  0.04),
    mean_se     = c(0.04,  0.04),
    stringsAsFactors = FALSE)
  res <- list(metrics = m, effect_sizes = 0.05, reps = 5L,
              tmle_candidates = list(cand_a, cand_b))
  class(res) <- "plasmode_results"
  out <- gate_check(res, rmse_threshold = 0.06,
                    coverage_threshold = 0.90, method = "a")
  expect_true(out$decision %in% c("GO", "FLAG", "STOP"))
})

test_that("tmle_candidate accepts deprecated Q_library alias", {
  expect_warning(c <- tmle_candidate("x", g_library = "SL.glm",
                                      Q_library = c("SL.glm", "SL.mean")),
                 "deprecated")
  expect_equal(c$q_library, c("SL.glm", "SL.mean"))
})
