# Tests for staged workflow infrastructure (R/stages.R and the lock
# declarations). The checkpoint, audit-log, gate, and old decision-log
# layer was removed in 0.3.0; what remains here are the declarations,
# the design-precision and event-support readers, the candidate
# specifications, and outcome masking under the outcome-store split.

test_that("define_negative_control registers variable", {
  dat  <- sim_func1(n = 100, seed = 1)
  lock <- create_analysis_lock(dat, "treatment", "event_24",
                               c("age", "sex"), seed = 1)
  lock <- cleanTMLE:::define_negative_control(lock, "nc_outcome", description = "test")
  expect_true("nc_outcome" %in% names(lock$negative_controls))
  expect_equal(lock$negative_controls$nc_outcome$type, "outcome")
})

test_that("define_negative_control errors on missing variable", {
  dat  <- sim_func1(n = 100, seed = 1)
  lock <- create_analysis_lock(dat, "treatment", "event_24",
                               c("age", "sex"), seed = 1)
  expect_error(
    cleanTMLE:::define_negative_control(lock, "nonexistent"),
    "not found"
  )
})

# ── declare_negative_controls (the exported verb) ────────────────────────

test_that("declare_negative_controls registers controls with domains and logs", {
  dat  <- sim_func1(n = 100, seed = 1)
  lock <- create_analysis_lock(dat, "treatment", "event_24",
                               c("age", "sex"), seed = 1)
  lock <- declare_negative_controls(lock, "nc_outcome",
                                    domains = "health_seeking_behavior")
  expect_true("nc_outcome" %in% names(lock$negative_controls))
  expect_equal(lock$negative_controls$nc_outcome$domain,
               "health_seeking_behavior")
  expect_true(any(lock$design_log$type == "negative_controls"))
})

test_that("nc_criteria declared after creation validate and are logged", {
  dat  <- sim_func1(n = 100, seed = 1)
  lock <- create_analysis_lock(dat, "treatment", "event_24",
                               c("age", "sex"), seed = 1)
  expect_false(isTRUE(lock$declared_at_creation$nc_criteria))
  lock <- declare_negative_controls(lock, "nc_outcome",
                                    nc_criteria = list(null_band = 0.02))
  expect_equal(lock$nc_criteria$null_band, 0.02)
  # Post-creation criteria live in the design log, not the fingerprint:
  # validation must still pass.
  expect_message(validate_analysis_lock(lock), "validated")
  expect_true(any(lock$design_log$type == "nc_criteria"))
})

test_that("nc_criteria cannot be redeclared", {
  dat  <- sim_func1(n = 100, seed = 1)
  lock <- create_analysis_lock(dat, "treatment", "event_24",
                               c("age", "sex"), seed = 1,
                               nc_criteria = list(null_band = 0.02))
  expect_true(isTRUE(lock$declared_at_creation$nc_criteria))
  expect_message(validate_analysis_lock(lock), "validated")
  expect_error(
    declare_negative_controls(lock, "nc_outcome",
                              nc_criteria = list(null_band = 0.05)),
    "already declared")
})

# ── estimand and sensitivity plans as lock arguments ─────────────────────

test_that("create_analysis_lock accepts estimand and sensitivity_plans", {
  dat  <- sim_func1(n = 100, seed = 1)
  lock <- create_analysis_lock(
    dat, "treatment", "event_24", c("age", "sex"), seed = 1,
    estimand = list(description = "24-week event risk difference",
                    population  = "treated-eligible adults"),
    sensitivity_plans = list(
      truncation_sweep = list(
        description = "Vary the PS truncation bound",
        settings    = list(truncation = c(0.01, 0.05, 0.10)))))
  expect_equal(lock$estimand$description,
               "24-week event risk difference")
  expect_equal(lock$estimand$contrast, "risk_difference")
  expect_equal(lock$sensitivity_plans$truncation_sweep$label,
               "truncation_sweep")
  # Metadata is outside the fingerprint: an identical lock without it
  # hashes the same, and validation passes with it.
  plain <- create_analysis_lock(dat, "treatment", "event_24",
                                c("age", "sex"), seed = 1)
  expect_identical(lock$lock_hash, plain$lock_hash)
  expect_message(validate_analysis_lock(lock), "validated")
})

test_that("unknown estimand fields and unnamed plans are rejected", {
  dat <- sim_func1(n = 100, seed = 1)
  expect_error(
    create_analysis_lock(dat, "treatment", "event_24", c("age", "sex"),
                         seed = 1, estimand = list(nonsense = 1)),
    "Unknown estimand field")
  expect_error(
    create_analysis_lock(dat, "treatment", "event_24", c("age", "sex"),
                         seed = 1,
                         sensitivity_plans = list(list(description = "x"))),
    "named list")
})

test_that("compute_evalue returns correct structure", {
  ev <- compute_evalue(2.0)
  expect_true("e_value" %in% names(ev))
  expect_true(ev["e_value"] > 2)

  ev2 <- compute_evalue(2.0, ci_bound = 1.5)
  expect_true("e_value_ci" %in% names(ev2))

  # RR = 1 should give e_value = 1
  ev3 <- compute_evalue(1.0)
  expect_equal(unname(ev3["e_value"]), 1)
})

test_that("sim_func1 includes nc_outcome column", {
  dat <- sim_func1(n = 100, seed = 1)
  expect_true("nc_outcome" %in% names(dat))
  expect_true(all(dat$nc_outcome %in% c(0L, 1L)))
})

test_that("print methods do not error", {
  dat  <- sim_func1(n = 200, seed = 1)
  lock <- create_analysis_lock(dat, "treatment", "event_24",
                               c("age", "sex", "biomarker"), seed = 1)
  lock <- cleanTMLE:::define_negative_control(lock, "nc_outcome")
  expect_output(print(lock), "Negative Controls")
  expect_output(print(lock), "separate store")
})

# ── tmle_candidate infrastructure ────────────────────────────────────────

test_that("tmle_candidate creates a tmle_candidate_spec", {
  cand <- tmle_candidate("test_cand", "Test Candidate",
                         g_library = c("SL.glm"), truncation = 0.05)
  expect_s3_class(cand, "tmle_candidate_spec")
  expect_equal(cand$candidate_id, "test_cand")
  expect_equal(cand$label, "Test Candidate")
  expect_equal(cand$truncation, 0.05)
  expect_equal(cand$g_library, c("SL.glm"))
})

test_that("validate_tmle_candidates rejects non-specs", {
  cands <- list(
    tmle_candidate("a", g_library = "SL.glm"),
    list(not_a_spec = TRUE)
  )
  expect_error(cleanTMLE:::validate_tmle_candidates(cands), "not a tmle_candidate_spec")
})

test_that("validate_tmle_candidates rejects duplicates", {
  cands <- list(
    tmle_candidate("dup_id", g_library = "SL.glm"),
    tmle_candidate("dup_id", g_library = "SL.glm", truncation = 0.05)
  )
  expect_error(cleanTMLE:::validate_tmle_candidates(cands), "Duplicate")
})

test_that("define_candidates wraps the single and grid constructors", {
  one <- define_candidates("glm_t01", g_library = "SL.glm",
                           truncation = 0.01)
  expect_s3_class(one, "ct_candidates")
  expect_equal(length(one), 1L)
  expect_s3_class(one[["glm_t01"]], "tmle_candidate_spec")

  grid <- define_candidates(grid = list(
    truncations = c(0.01, 0.05),
    libraries   = list(glm = "SL.glm")))
  expect_s3_class(grid, "ct_candidates")
  expect_equal(length(grid), 2L)
  expect_output(print(grid), "Candidate set")
  df <- as.data.frame(grid)
  expect_equal(nrow(df), 2L)
  expect_true(all(c("candidate_id", "truncation") %in% names(df)))

  expect_error(define_candidates("x", grid = list(truncations = 0.01)),
               "not both")
  expect_error(define_candidates(grid = list(bogus_axis = 1)),
               "Unknown grid axis")
  expect_error(define_candidates(), "candidate_id")
})

test_that("lock_primary_tmle_spec and get_primary_tmle_spec work", {
  dat  <- sim_func1(n = 100, seed = 1)
  lock <- create_analysis_lock(dat, "treatment", "event_24",
                               c("age", "sex"), seed = 1)

  expect_null(cleanTMLE:::get_primary_tmle_spec(lock))

  cand <- tmle_candidate("test_lock", g_library = "SL.glm", truncation = 0.02)
  lock <- cleanTMLE:::lock_primary_tmle_spec(lock, cand)

  spec <- cleanTMLE:::get_primary_tmle_spec(lock)
  expect_s3_class(spec, "tmle_candidate_spec")
  expect_equal(spec$candidate_id, "test_lock")
  expect_equal(spec$truncation, 0.02)
})

test_that("print.cleanroom_lock shows primary TMLE spec", {
  dat  <- sim_func1(n = 100, seed = 1)
  lock <- create_analysis_lock(dat, "treatment", "event_24",
                               c("age", "sex"), seed = 1)
  cand <- tmle_candidate("show_spec", "Display Spec",
                         g_library = "SL.glm", truncation = 0.03)
  lock <- cleanTMLE:::lock_primary_tmle_spec(lock, cand)

  expect_output(print(lock), "Primary TMLE Specification")
  expect_output(print(lock), "show_spec")
})

# ── Design precision and event support ───────────────────────────────────

test_that("estimate_design_precision returns correct structure", {
  dat  <- sim_func1(n = 300, seed = 1)
  lock <- create_analysis_lock(dat, "treatment", "event_24",
                               c("age", "sex", "biomarker"), seed = 1)
  dp <- estimate_design_precision(lock)
  expect_s3_class(dp, "design_precision")
  expect_true(dp$n_total == 300)
  expect_true(dp$n_treated + dp$n_control == dp$n_total)
  expect_true(is.numeric(dp$se_proxy))
  expect_true(is.numeric(dp$mdd_80))
  expect_null(dp$target_mdd)
  expect_output(print(dp), "Design-Stage Precision")
})

test_that("estimate_design_precision compares against target MDD", {
  dat  <- sim_func1(n = 300, seed = 1)
  lock <- create_analysis_lock(dat, "treatment", "event_24",
                               c("age", "sex", "biomarker"), seed = 1)
  dp <- estimate_design_precision(lock, target_mdd = 0.50)
  expect_true(isTRUE(dp$mdd_feasible))
  dp2 <- estimate_design_precision(lock, target_mdd = 0.001)
  expect_false(dp2$mdd_feasible)
})

test_that("summarize_event_support is marginal only", {
  dat  <- sim_func1(n = 300, seed = 1)
  lock <- create_analysis_lock(dat, "treatment", "event_24",
                               c("age", "sex", "biomarker"), seed = 1)
  es <- cleanTMLE:::summarize_event_support(lock)
  expect_true(is.data.frame(es))
  expect_equal(nrow(es), 1)
  expect_true(all(c("n", "events", "event_rate") %in% names(es)))
  expect_false("arm" %in% names(es))
  expect_equal(es$events, sum(dat$event_24))
})

test_that("event_support_by_arm requires a reason, warns, and logs", {
  dat  <- sim_func1(n = 300, seed = 1)
  lock <- create_analysis_lock(dat, "treatment", "event_24",
                               c("age", "sex", "biomarker"), seed = 1)
  expect_error(event_support_by_arm(lock), "reason")
  expect_warning(
    esa <- event_support_by_arm(lock, reason = "sparse-cell check"),
    "crude treatment-outcome association")
  expect_s3_class(esa, "event_support_by_arm")
  expect_equal(nrow(esa$table), 3)
  expect_true(all(c("arm", "n", "events", "event_rate") %in%
                    names(esa$table)))
  # The access is on the returned lock's design log.
  expect_true(any(esa$lock$design_log$type == "event_support_by_arm"))
  expect_true(any(grepl("sparse-cell check", esa$lock$design_log$note)))
  # And the per-arm counts reconcile with the marginal summary.
  es <- cleanTMLE:::summarize_event_support(lock)
  expect_equal(sum(esa$table$events[esa$table$arm != "Total"]), es$events)
})

test_that("estimate_design_precision output carries no per-arm outcome split", {
  dat  <- sim_func1(n = 300, seed = 1)
  lock <- create_analysis_lock(dat, "treatment", "event_24",
                               c("age", "sex", "biomarker"), seed = 1)
  dp <- estimate_design_precision(lock)
  expect_null(dp$events_per_arm)
  expect_null(dp$crude_rates)
  expect_false(any(grepl("arm", names(dp$event_support))))
  expect_output(print(dp), "marginal")
})

# ── Outcome masking under the store split ────────────────────────────────

test_that("mask_outcome removes the outcome store", {
  dat  <- sim_func1(n = 100, seed = 1)
  lock <- create_analysis_lock(dat, "treatment", "event_24",
                               c("age", "sex"), seed = 1)
  masked <- mask_outcome(lock)
  expect_null(masked$outcome_store)
  expect_false("event_24" %in% names(masked$data))
  expect_true(isTRUE(masked$.outcome_masked))
})

test_that("unmask_outcome restores the outcome from the original lock", {
  dat  <- sim_func1(n = 100, seed = 1)
  lock <- create_analysis_lock(dat, "treatment", "event_24",
                               c("age", "sex"), seed = 1)
  masked   <- mask_outcome(lock)
  unmasked <- unmask_outcome(masked, lock)
  expect_identical(cleanTMLE:::.outcome_vector(unmasked),
                   cleanTMLE:::.outcome_vector(lock))
  expect_false(isTRUE(unmasked$.outcome_masked))
  expect_true(any(unmasked$design_log$type == "outcome_unmasked"))
})

test_that("Stage 4 functions reject masked outcome", {
  dat  <- sim_func1(n = 200, seed = 1)
  lock <- create_analysis_lock(dat, "treatment", "event_24",
                               c("age", "sex"), seed = 1)
  masked <- mask_outcome(lock)
  expect_error(cleanTMLE:::run_crude_workflow(masked), "Outcome is masked")
})

test_that("Stage 4 functions allow the per-call override", {
  dat  <- sim_func1(n = 200, seed = 1)
  lock <- create_analysis_lock(dat, "treatment", "event_24",
                               c("age", "sex"), seed = 1)
  result <- cleanTMLE:::run_crude_workflow(lock, allow_outcome_access = TRUE)
  expect_true(is.numeric(result$estimate))
})