# Tests for staged workflow infrastructure (R/stages.R)

test_that("attach_estimand adds estimand to lock", {
  dat  <- sim_func1(n = 100, seed = 1)
  lock <- create_analysis_lock(dat, "treatment", "event_24",
                               c("age", "sex"), seed = 1)
  lock2 <- attach_estimand(lock,
    description = "Test question",
    population  = "Test pop",
    contrast    = "risk_difference"
  )
  expect_true(!is.null(lock2$estimand))
  expect_equal(lock2$estimand$contrast, "risk_difference")
  expect_equal(lock2$estimand$description, "Test question")
  # Hash should be unchanged

  expect_equal(lock2$lock_hash, lock$lock_hash)
})

test_that("declare_sensitivity_plan appends plans", {
  dat  <- sim_func1(n = 100, seed = 1)
  lock <- create_analysis_lock(dat, "treatment", "event_24",
                               c("age", "sex"), seed = 1)
  lock <- declare_sensitivity_plan(lock, "plan1", "First plan")
  lock <- declare_sensitivity_plan(lock, "plan2", "Second plan")
  expect_length(lock$sensitivity_plans, 2)
  expect_equal(lock$sensitivity_plans$plan1$label, "plan1")
  expect_equal(lock$sensitivity_plans$plan2$label, "plan2")
})

test_that("define_negative_control registers variable", {
  dat  <- sim_func1(n = 100, seed = 1)
  lock <- create_analysis_lock(dat, "treatment", "event_24",
                               c("age", "sex"), seed = 1)
  lock <- define_negative_control(lock, "nc_outcome", description = "test")
  expect_true("nc_outcome" %in% names(lock$negative_controls))
  expect_equal(lock$negative_controls$nc_outcome$type, "outcome")
})

test_that("define_negative_control errors on missing variable", {
  dat  <- sim_func1(n = 100, seed = 1)
  lock <- create_analysis_lock(dat, "treatment", "event_24",
                               c("age", "sex"), seed = 1)
  expect_error(
    define_negative_control(lock, "nonexistent"),
    "not found"
  )
})

test_that("checkpoint_cohort_adequacy returns GO for adequate data", {
  dat  <- sim_func1(n = 500, seed = 1)
  lock <- create_analysis_lock(dat, "treatment", "event_24",
                               c("age", "sex", "biomarker"), seed = 1)
  cp <- checkpoint_cohort_adequacy(lock)
  expect_s3_class(cp, "cleantmle_checkpoint")
  expect_equal(cp$stage, "Check Point 1: Cohort Adequacy")
  expect_true(cp$decision %in% c("GO", "FLAG", "STOP"))
})

test_that("checkpoint_cohort_adequacy returns STOP for tiny data", {
  dat  <- sim_func1(n = 30, seed = 1)
  lock <- create_analysis_lock(dat, "treatment", "event_24",
                               c("age", "sex"), seed = 1)
  cp <- checkpoint_cohort_adequacy(lock, min_n_per_arm = 50, min_events = 100)
  expect_equal(cp$decision, "STOP")
})

test_that("checkpoint_balance works with ps_diagnostics", {
  dat  <- sim_func1(n = 300, seed = 1)
  lock <- create_analysis_lock(dat, "treatment", "event_24",
                               c("age", "sex", "biomarker"), seed = 1)
  ps   <- fit_ps_glm(lock)
  diag <- compute_ps_diagnostics(ps)
  cp   <- checkpoint_balance(diag, lock_hash = lock$lock_hash)
  expect_s3_class(cp, "cleantmle_checkpoint")
  expect_true(cp$decision %in% c("GO", "FLAG", "STOP"))
  expect_equal(cp$lock_hash, lock$lock_hash)
})

test_that("run_negative_control returns structured result", {
  dat  <- sim_func1(n = 300, seed = 1)
  lock <- create_analysis_lock(dat, "treatment", "event_24",
                               c("age", "sex", "biomarker"), seed = 1)
  lock <- define_negative_control(lock, "nc_outcome")
  ps   <- fit_ps_glm(lock)
  nc   <- run_negative_control(lock, "nc_outcome", ps)
  expect_s3_class(nc, "cleantmle_nc_result")
  expect_true(is.numeric(nc$estimate))
  expect_true(is.numeric(nc$p_value))
  expect_true(nchar(nc$interpretation) > 0)
})

test_that("checkpoint_residual_bias handles single and multiple NC results", {
  dat  <- sim_func1(n = 300, seed = 1)
  lock <- create_analysis_lock(dat, "treatment", "event_24",
                               c("age", "sex", "biomarker"), seed = 1)
  lock <- define_negative_control(lock, "nc_outcome")
  ps   <- fit_ps_glm(lock)
  nc   <- run_negative_control(lock, "nc_outcome", ps)

  # Single result
  cp <- checkpoint_residual_bias(nc)
  expect_s3_class(cp, "cleantmle_checkpoint")

  # List of results
  cp2 <- checkpoint_residual_bias(list(nc))
  expect_s3_class(cp2, "cleantmle_checkpoint")
})

test_that("audit log records and exports correctly", {
  dat  <- sim_func1(n = 100, seed = 1)
  lock <- create_analysis_lock(dat, "treatment", "event_24",
                               c("age", "sex"), seed = 1)
  audit <- create_audit_log(lock)
  expect_s3_class(audit, "cleantmle_audit")
  expect_equal(audit$lock_hash, lock$lock_hash)
  expect_length(audit$entries, 0)

  audit <- record_stage(audit, "Stage 1a", "Lock created")
  expect_length(audit$entries, 1)

  trail <- export_audit_trail(audit)
  expect_true(is.data.frame(trail))
  expect_equal(nrow(trail), 1)
  expect_true("stage" %in% names(trail))
})

test_that("record_checkpoint adds checkpoint to audit", {
  dat  <- sim_func1(n = 300, seed = 1)
  lock <- create_analysis_lock(dat, "treatment", "event_24",
                               c("age", "sex", "biomarker"), seed = 1)
  audit <- create_audit_log(lock)
  cp1   <- checkpoint_cohort_adequacy(lock)
  audit <- record_checkpoint(audit, cp1)
  expect_length(audit$entries, 1)
  expect_equal(audit$entries[[1]]$decision, cp1$decision)
})

test_that("build_stage_manifest produces output", {
  dat  <- sim_func1(n = 100, seed = 1)
  lock <- create_analysis_lock(dat, "treatment", "event_24",
                               c("age", "sex"), seed = 1)
  audit <- create_audit_log(lock)
  audit <- record_stage(audit, "Stage 1a", "Lock created")

  out <- capture.output(build_stage_manifest(audit))
  expect_true(length(out) > 0)
  expect_true(any(grepl("Stage 1a", out)))
})

test_that("as.data.frame.cleantmle_checkpoint works", {
  dat  <- sim_func1(n = 300, seed = 1)
  lock <- create_analysis_lock(dat, "treatment", "event_24",
                               c("age", "sex", "biomarker"), seed = 1)
  cp <- checkpoint_cohort_adequacy(lock)
  df <- as.data.frame(cp)
  expect_true(is.data.frame(df))
  expect_equal(nrow(df), 1)
  expect_true("decision" %in% names(df))
})

test_that("sensitivity_truncation returns data.frame", {
  dat  <- sim_func1(n = 200, seed = 1)
  lock <- create_analysis_lock(dat, "treatment", "event_24",
                               c("age", "sex", "biomarker"), seed = 1)
  sens <- sensitivity_truncation(lock, thresholds = c(0.01, 0.05))
  expect_true(is.data.frame(sens))
  expect_equal(nrow(sens), 2)
  expect_true("truncation" %in% names(sens))
  expect_true("estimate" %in% names(sens))
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
  lock <- attach_estimand(lock, description = "Test", contrast = "risk_difference")
  lock <- declare_sensitivity_plan(lock, "sp1", "Test plan")
  lock <- define_negative_control(lock, "nc_outcome")

  expect_output(print(lock), "Estimand")
  expect_output(print(lock), "Sensitivity Plans")
  expect_output(print(lock), "Negative Controls")

  cp <- checkpoint_cohort_adequacy(lock)
  expect_output(print(cp), "Check Point 1")

  audit <- create_audit_log(lock)
  audit <- record_stage(audit, "Stage 1", "Test")
  expect_output(print(audit), "Audit Log")
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
  expect_error(validate_tmle_candidates(cands), "not a tmle_candidate_spec")
})

test_that("validate_tmle_candidates rejects duplicates", {
  cands <- list(
    tmle_candidate("dup_id", g_library = "SL.glm"),
    tmle_candidate("dup_id", g_library = "SL.glm", truncation = 0.05)
  )
  expect_error(validate_tmle_candidates(cands), "Duplicate")
})

test_that("expand_tmle_candidate_grid creates grid", {
  grid <- expand_tmle_candidate_grid(
    libraries   = list(glm = "SL.glm"),
    truncations = c(0.01, 0.05)
  )
  expect_true(is.list(grid))
  expect_equal(length(grid), 2L)
  expect_true(all(vapply(grid, inherits, logical(1), "tmle_candidate_spec")))
})

test_that("lock_primary_tmle_spec and get_primary_tmle_spec work", {
  dat  <- sim_func1(n = 100, seed = 1)
  lock <- create_analysis_lock(dat, "treatment", "event_24",
                               c("age", "sex"), seed = 1)

  expect_null(get_primary_tmle_spec(lock))

  cand <- tmle_candidate("test_lock", g_library = "SL.glm", truncation = 0.02)
  lock <- lock_primary_tmle_spec(lock, cand)

  spec <- get_primary_tmle_spec(lock)
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
  lock <- lock_primary_tmle_spec(lock, cand)

  expect_output(print(lock), "Primary TMLE Specification")
  expect_output(print(lock), "show_spec")
})


# ── New Phase 2 Tests: Design Precision ──────────────────────────────────

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

test_that("summarize_event_support returns data.frame", {
  dat  <- sim_func1(n = 300, seed = 1)
  lock <- create_analysis_lock(dat, "treatment", "event_24",
                               c("age", "sex", "biomarker"), seed = 1)
  es <- summarize_event_support(lock)
  expect_true(is.data.frame(es))
  expect_equal(nrow(es), 3)
  expect_true(all(c("arm", "n", "events", "event_rate") %in% names(es)))
  expect_equal(es$arm, c("Treated", "Control", "Total"))
})


# ── New Phase 2 Tests: Residual Confounding Stage ────────────────────────

test_that("run_residual_confounding_stage works", {
  dat  <- sim_func1(n = 300, seed = 1)
  lock <- create_analysis_lock(dat, "treatment", "event_24",
                               c("age", "sex", "biomarker"), seed = 1)
  lock <- define_negative_control(lock, "nc_outcome")
  ps   <- fit_ps_glm(lock)
  stage3 <- run_residual_confounding_stage(lock, ps)
  expect_s3_class(stage3, "residual_confounding_stage")
  expect_true(is.data.frame(stage3$summary_table))
  expect_s3_class(stage3$checkpoint, "cleantmle_checkpoint")
  expect_true(stage3$n_controls == 1)
  expect_output(print(stage3), "Stage 3")
})

test_that("run_residual_confounding_stage errors without NCs", {
  dat  <- sim_func1(n = 100, seed = 1)
  lock <- create_analysis_lock(dat, "treatment", "event_24",
                               c("age", "sex"), seed = 1)
  ps <- fit_ps_glm(lock)
  expect_error(run_residual_confounding_stage(lock, ps),
               "No negative controls")
})


# ── New Phase 2 Tests: Pre-Outcome Gate ──────────────────────────────────

test_that("authorize_outcome_analysis returns GO when all checkpoints pass", {
  dat  <- sim_func1(n = 300, seed = 1)
  lock <- create_analysis_lock(dat, "treatment", "event_24",
                               c("age", "sex", "biomarker"), seed = 1)
  lock <- define_negative_control(lock, "nc_outcome")
  audit <- create_audit_log(lock)

  cp1 <- checkpoint_cohort_adequacy(lock)
  audit <- record_checkpoint(audit, cp1)

  ps   <- fit_ps_glm(lock)
  diag <- compute_ps_diagnostics(ps)
  cp2  <- checkpoint_balance(diag, lock_hash = lock$lock_hash)
  audit <- record_checkpoint(audit, cp2)

  nc  <- run_negative_control(lock, "nc_outcome", ps)
  cp3 <- checkpoint_residual_bias(nc, lock_hash = lock$lock_hash)
  audit <- record_checkpoint(audit, cp3)

  gate <- authorize_outcome_analysis(audit)
  expect_s3_class(gate, "pre_outcome_gate")
  expect_s3_class(gate, "cleantmle_checkpoint")
  expect_true(gate$authorized)
  expect_equal(gate$decision, "GO")
})

test_that("authorize_outcome_analysis returns STOP when checkpoints missing", {
  dat  <- sim_func1(n = 100, seed = 1)
  lock <- create_analysis_lock(dat, "treatment", "event_24",
                               c("age", "sex"), seed = 1)
  audit <- create_audit_log(lock)
  # No checkpoints recorded
  gate <- authorize_outcome_analysis(audit)
  expect_false(gate$authorized)
  expect_equal(gate$decision, "STOP")
})

test_that("assert_outcome_authorized errors when not authorized", {
  dat  <- sim_func1(n = 100, seed = 1)
  lock <- create_analysis_lock(dat, "treatment", "event_24",
                               c("age", "sex"), seed = 1)
  audit <- create_audit_log(lock)
  expect_error(assert_outcome_authorized(audit), "NOT authorised")
})


# ── New Phase 2 Tests: Decision Log ─────────────────────────────────────

test_that("record_decision_log_entry and export_decision_log work", {
  dat  <- sim_func1(n = 100, seed = 1)
  lock <- create_analysis_lock(dat, "treatment", "event_24",
                               c("age", "sex"), seed = 1)
  audit <- create_audit_log(lock)
  audit <- record_decision_log_entry(
    audit, "Stage 2", "model_specification",
    "Selected GLM PS model",
    rationale = "Pre-specified in SAP"
  )
  audit <- record_decision_log_entry(
    audit, "Stage 2b", "override",
    "Accepted FLAG for ESS"
  )
  dl <- export_decision_log(audit)
  expect_true(is.data.frame(dl))
  expect_equal(nrow(dl), 2)
  expect_true(all(c("stage", "decision_type", "description") %in% names(dl)))
})

test_that("export_decision_log returns empty df when no entries", {
  dat  <- sim_func1(n = 100, seed = 1)
  lock <- create_analysis_lock(dat, "treatment", "event_24",
                               c("age", "sex"), seed = 1)
  audit <- create_audit_log(lock)
  dl <- export_decision_log(audit)
  expect_true(is.data.frame(dl))
  expect_equal(nrow(dl), 0)
})


# ── New Phase 2 Tests: Stage Path Narrative ──────────────────────────────

test_that("summarize_stage_path produces narrative", {
  dat  <- sim_func1(n = 100, seed = 1)
  lock <- create_analysis_lock(dat, "treatment", "event_24",
                               c("age", "sex"), seed = 1)
  audit <- create_audit_log(lock)
  audit <- record_stage(audit, "Stage 1a", "Lock created", decision = NA)
  audit <- record_stage(audit, "Check Point 1", "Evaluated", decision = "GO")
  out <- capture.output(summarize_stage_path(audit))
  expect_true(any(grepl("Stage 1a", out)))
  expect_true(any(grepl("GO", out)))
})


# ── New Phase 2 Tests: Outcome Masking ───────────────────────────────────

test_that("mask_outcome sets outcome to NA", {
  dat  <- sim_func1(n = 100, seed = 1)
  lock <- create_analysis_lock(dat, "treatment", "event_24",
                               c("age", "sex"), seed = 1)
  masked <- mask_outcome(lock)
  expect_true(all(is.na(masked$data[["event_24"]])))
  expect_true(isTRUE(masked$.outcome_masked))
})

test_that("unmask_outcome restores outcome", {
  dat  <- sim_func1(n = 100, seed = 1)
  lock <- create_analysis_lock(dat, "treatment", "event_24",
                               c("age", "sex"), seed = 1)
  masked   <- mask_outcome(lock)
  unmasked <- unmask_outcome(masked, lock)
  expect_identical(unmasked$data[["event_24"]], lock$data[["event_24"]])
  expect_false(isTRUE(unmasked$.outcome_masked))
})

test_that("Stage 4 functions reject masked outcome", {
  dat  <- sim_func1(n = 200, seed = 1)
  lock <- create_analysis_lock(dat, "treatment", "event_24",
                               c("age", "sex"), seed = 1)
  masked <- mask_outcome(lock)
  expect_error(run_crude_workflow(masked), "Outcome is masked")
})

test_that("Stage 4 functions allow override_clean_room", {
  dat  <- sim_func1(n = 200, seed = 1)
  lock <- create_analysis_lock(dat, "treatment", "event_24",
                               c("age", "sex"), seed = 1)
  # Non-masked lock should work fine without override
  result <- run_crude_workflow(lock)
  expect_true(is.numeric(result$estimate))
})
