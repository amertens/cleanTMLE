# Tests for the governance-note formatters (event-process
# classification, coherence check, target population, missing-data
# plan, cumulative-risk reporting), folded back from the short-lived
# cleanroomGov companion package in 0.3.0. The estimand and
# sensitivity-plan lock arguments are tested in test-stages.R.

# A. clean_event_process_table -----------------------------------------

test_that("clean_event_process_table returns required columns", {
  out <- clean_event_process_table(list(
    list(event_name = "Y",     event_variable = "y",
         role_in_estimand = "outcome event",
         primary_handling = "fit Q on Y"),
    list(event_name = "Drop",  event_variable = "lost",
         role_in_estimand = "administrative censoring",
         primary_handling = "IPCW")))
  required <- c("event_name", "event_variable", "event_timing_variable",
                "role_in_estimand", "role_in_estimation",
                "ICH_E9R1_strategy", "justification", "primary_handling",
                "sensitivity_handling",
                "affects_identification_assumption", "notes")
  expect_true(all(required %in% names(out)))
  expect_s3_class(out, "cleantmle_event_process")
  expect_equal(nrow(out), 2L)
})

test_that("death-as-censoring without justification triggers a warning", {
  expect_warning(clean_event_process_table(list(
    list(event_name = "Y",       event_variable = "y",
         role_in_estimand = "outcome event"),
    list(event_name = "Death",   event_variable = "dod",
         role_in_estimand = "informative censoring",
         primary_handling = "IPCW"))),
    "death is classified as")
})

test_that("censoring without handling triggers a warning", {
  expect_warning(clean_event_process_table(list(
    list(event_name = "Y",         event_variable = "y",
         role_in_estimand = "outcome event"),
    list(event_name = "Loss",      event_variable = "lost",
         role_in_estimand = "administrative censoring"))),
    "no censoring handling is documented")
})

# B. clean_check_event_processes ---------------------------------------

test_that("coherent component risks pass without flagging", {
  df <- data.frame(time_point = c(180, 365),
                   event_of_interest_risk = c(0.05, 0.09),
                   competing_event_risk   = c(0.03, 0.05),
                   composite_risk         = c(0.08, 0.14))
  out <- clean_check_event_processes(df)
  expect_true(all(out$flag == "OK"))
  expect_equal(out$component_sum, c(0.08, 0.14))
})

test_that("component risks exceeding composite trigger a warning", {
  df <- data.frame(time_point = 365,
                   event_of_interest_risk = 0.10,
                   competing_event_risk   = 0.06,
                   composite_risk         = 0.13)
  expect_warning(out <- clean_check_event_processes(df),
                 "exceeds the composite risk")
  expect_equal(out$flag, "COMPONENT_EXCEEDS_COMPOSITE")
})

test_that("check_event_processes matches hand-computed sums", {
  df <- data.frame(time_point = 365,
                   event_of_interest_risk = 0.09,
                   competing_event_risk   = 0.04,
                   composite_risk         = 0.13)
  out <- clean_check_event_processes(df, tolerance = 1e-10)
  expect_equal(out$component_sum, 0.13)
  expect_equal(out$flag, "OK")
})

# C. clean_target_population -------------------------------------------

test_that("positivity-restricted population without rule warns", {
  expect_warning(clean_target_population(
    "positivity-supported restricted population",
    reference_group = "Control"),
    "no `restriction_rule`")
})

test_that("ordinary call returns target_population object", {
  tp <- clean_target_population("full eligible population",
                                 reference_group = "Control",
                                 rationale = "Marginal ATE")
  expect_s3_class(tp, "cleantmle_target_population")
  expect_equal(tp$target_population, "full eligible population")
})

# D. clean_missing_data_plan -------------------------------------------

test_that("missing-data plan distinguishes baseline / outcome / censoring", {
  plan <- clean_missing_data_plan(list(
    "baseline covariate missingness" = list(
      timing = "baseline", presumed_mechanism = "MCAR",
      handling_primary = "complete-case"),
    "outcome missingness" = list(
      timing = "follow-up", presumed_mechanism = "MAR",
      handling_primary = "IPCW"),
    "censoring or loss to follow-up" = list(
      timing = "follow-up", presumed_mechanism = "administrative",
      handling_primary = "IPCW")))
  expect_s3_class(plan, "cleantmle_missing_plan")
  expect_equal(nrow(plan), 3L)
  expect_setequal(plan$variable_or_process,
                  c("baseline covariate missingness",
                    "outcome missingness",
                    "censoring or loss to follow-up"))
})

# E. clean_risk_report_table -------------------------------------------

test_that("risk report table has required columns and tolerates missing optionals", {
  rows <- list(
    list(time_point = 180, treatment_strategy = "Treated",
         risk = 0.05),
    list(time_point = 180, treatment_strategy = "Control",
         risk = 0.07, risk_ci_lower = 0.05, risk_ci_upper = 0.09,
         estimator = "TMLE"))
  out <- clean_risk_report_table(rows)
  required <- c("time_point", "treatment_strategy", "risk",
                "risk_ci_lower", "risk_ci_upper", "risk_difference",
                "risk_ratio", "estimator", "nuisance_specification",
                "notes", "n_observed", "events", "censoring_events",
                "competing_events", "person_time",
                "effective_sample_size")
  expect_true(all(required %in% names(out)))
  expect_s3_class(out, "cleantmle_risk_report")
  expect_equal(out$risk[1], 0.05)
  expect_true(is.na(out$risk_ci_lower[1]))
})

test_that("risk report table errors on missing required field", {
  expect_error(clean_risk_report_table(list(
    list(time_point = 180, treatment_strategy = "Treated"))),
    "Each row needs")
})

test_that("risk report table preserves numeric types", {
  rows <- list(
    list(time_point = 365, treatment_strategy = "Treated", risk = 0.12,
         events = 34, n_observed = 300, estimator = "TMLE"),
    list(time_point = 365, treatment_strategy = "Control", risk = 0.17,
         events = 51, n_observed = 300, risk_difference = -0.05,
         estimator = "TMLE"))
  out <- clean_risk_report_table(rows)
  expect_type(out$risk, "double")
  expect_type(out$events, "double")
  expect_equal(out$risk, c(0.12, 0.17))
})