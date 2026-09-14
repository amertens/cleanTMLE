# Regression tests for outcome NA handling, preset scenarios,
# attrition_table polymorphism, and the make_table1 lock method.
# The authorize/gate_check/print_locked_spec tests died with that
# layer in 0.3.0.
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
                                 data_quality_scenarios = cleanTMLE:::default_dq_scenarios("exploratory"),
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
    cfg <- cleanTMLE:::default_dq_scenarios(p)
    expect_named(cfg, c("covariate_missingness", "treatment_misclass",
                        "outcome_misclass", "unmeasured_confounding",
                        "near_positivity"))
    expect_true(!is.null(cfg$covariate_missingness$fractions))
    expect_true(!is.null(cfg$unmeasured_confounding$U_treatment_OR))
  }
})

test_that("default_dq_scenarios carries five threats including near_positivity", {
  for (p in c("regulatory_standard", "exploratory", "stress")) {
    cfg <- cleanTMLE:::default_dq_scenarios(p)
    expect_length(cfg, 5L)
    expect_true("near_positivity" %in% names(cfg))
    expect_true(!is.null(cfg$near_positivity$slopes))
    expect_true(all(cfg$near_positivity$slopes > 1))
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

test_that("tmle_candidate accepts deprecated Q_library alias", {
  expect_warning(c <- tmle_candidate("x", g_library = "SL.glm",
                                      Q_library = c("SL.glm", "SL.mean")),
                 "deprecated")
  expect_equal(c$q_library, c("SL.glm", "SL.mean"))
})
