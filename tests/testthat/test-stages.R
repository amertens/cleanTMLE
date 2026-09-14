# Tests for staged workflow infrastructure (R/stages.R)
# (attach_estimand and declare_sensitivity_plan moved to cleanroomGov, tested
# there.)

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

test_that("checkpoint_cohort_adequacy returns GO for adequate data", {
  dat  <- sim_func1(n = 500, seed = 1)
  lock <- create_analysis_lock(dat, "treatment", "event_24",
                               c("age", "sex", "biomarker"), seed = 1)
  cp <- cleanTMLE:::checkpoint_cohort_adequacy(lock)
  expect_s3_class(cp, "cleantmle_checkpoint")
  expect_equal(cp$stage, "Check Point 1: Cohort Adequacy")
  expect_true(cp$decision %in% c("GO", "FLAG", "STOP"))
})

test_that("checkpoint_cohort_adequacy returns STOP for tiny data", {
  dat  <- sim_func1(n = 30, seed = 1)
  lock <- create_analysis_lock(dat, "treatment", "event_24",
                               c("age", "sex"), seed = 1)
  cp <- cleanTMLE:::checkpoint_cohort_adequacy(lock, min_n_per_arm = 50, min_events = 100)
  expect_equal(cp$decision, "STOP")
})

test_that("checkpoint_balance works with ps_diagnostics", {
  dat  <- sim_func1(n = 300, seed = 1)
  lock <- create_analysis_lock(dat, "treatment", "event_24",
                               c("age", "sex", "biomarker"), seed = 1)
  ps   <- cleanTMLE:::fit_ps_glm(lock)
  diag <- cleanTMLE:::compute_ps_diagnostics(ps)
  cp   <- cleanTMLE:::checkpoint_balance(diag, lock_hash = lock$lock_hash)
  expect_s3_class(cp, "cleantmle_checkpoint")
  expect_true(cp$decision %in% c("GO", "FLAG", "STOP"))
  expect_equal(cp$lock_hash, lock$lock_hash)
})

test_that("run_negative_control returns structured result", {
  dat  <- sim_func1(n = 300, seed = 1)
  lock <- create_analysis_lock(dat, "treatment", "event_24",
                               c("age", "sex", "biomarker"), seed = 1)
  lock <- cleanTMLE:::define_negative_control(lock, "nc_outcome")
  ps   <- cleanTMLE:::fit_ps_glm(lock)
  nc   <- cleanTMLE:::run_negative_control(lock, "nc_outcome", ps)
  expect_s3_class(nc, "cleantmle_nc_result")
  expect_true(is.numeric(nc$estimate))
  expect_true(is.numeric(nc$p_value))
  expect_true(nchar(nc$interpretation) > 0)
})

test_that("checkpoint_residual_bias handles single and multiple NC results", {
  dat  <- sim_func1(n = 300, seed = 1)
  lock <- create_analysis_lock(dat, "treatment", "event_24",
                               c("age", "sex", "biomarker"), seed = 1)
  lock <- cleanTMLE:::define_negative_control(lock, "nc_outcome")
  ps   <- cleanTMLE:::fit_ps_glm(lock)
  nc   <- cleanTMLE:::run_negative_control(lock, "nc_outcome", ps)

  # Single result
  cp <- cleanTMLE:::checkpoint_residual_bias(nc)
  expect_s3_class(cp, "cleantmle_checkpoint")

  # List of results
  cp2 <- cleanTMLE:::checkpoint_residual_bias(list(nc))
  expect_s3_class(cp2, "cleantmle_checkpoint")
})

test_that("audit log records and exports correctly", {
  dat  <- sim_func1(n = 100, seed = 1)
  lock <- create_analysis_lock(dat, "treatment", "event_24",
                               c("age", "sex"), seed = 1)
  audit <- cleanTMLE:::create_audit_log(lock)
  expect_s3_class(audit, "cleantmle_audit")
  expect_equal(audit$lock_hash, lock$lock_hash)
  expect_length(audit$entries, 0)

  audit <- cleanTMLE:::record_stage(audit, "Stage 1a", "Lock created")
  expect_length(audit$entries, 1)

  trail <- cleanTMLE:::export_audit_trail(audit)
  expect_true(is.data.frame(trail))
  expect_equal(nrow(trail), 1)
  expect_true("stage" %in% names(trail))
})

test_that("record_checkpoint adds checkpoint to audit", {
  dat  <- sim_func1(n = 300, seed = 1)
  lock <- create_analysis_lock(dat, "treatment", "event_24",
                               c("age", "sex", "biomarker"), seed = 1)
  audit <- cleanTMLE:::create_audit_log(lock)
  cp1   <- cleanTMLE:::checkpoint_cohort_adequacy(lock)
  audit <- cleanTMLE:::record_checkpoint(audit, cp1)
  expect_length(audit$entries, 1)
  expect_equal(audit$entries[[1]]$decision, cp1$decision)
})

# (build_stage_manifest moved to cleanroomGov, tested there.)

test_that("as.data.frame.cleantmle_checkpoint works", {
  dat  <- sim_func1(n = 300, seed = 1)
  lock <- create_analysis_lock(dat, "treatment", "event_24",
                               c("age", "sex", "biomarker"), seed = 1)
  cp <- cleanTMLE:::checkpoint_cohort_adequacy(lock)
  df <- as.data.frame(cp)
  expect_true(is.data.frame(df))
  expect_equal(nrow(df), 1)
  expect_true("decision" %in% names(df))
})

test_that("sensitivity_truncation returns data.frame", {
  dat  <- sim_func1(n = 200, seed = 1)
  lock <- create_analysis_lock(dat, "treatment", "event_24",
                               c("age", "sex", "biomarker"), seed = 1)
  lock$.outcome_authorized <- TRUE  # authorised lock; sensitivity re-estimates
  sens <- cleanTMLE:::sensitivity_truncation(lock, thresholds = c(0.01, 0.05))
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
  lock <- cleanTMLE:::define_negative_control(lock, "nc_outcome")

  expect_output(print(lock), "Negative Controls")

  cp <- cleanTMLE:::checkpoint_cohort_adequacy(lock)
  expect_output(print(cp), "Check Point 1")

  audit <- cleanTMLE:::create_audit_log(lock)
  audit <- cleanTMLE:::record_stage(audit, "Stage 1", "Test")
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
  expect_error(cleanTMLE:::validate_tmle_candidates(cands), "not a tmle_candidate_spec")
})

test_that("validate_tmle_candidates rejects duplicates", {
  cands <- list(
    tmle_candidate("dup_id", g_library = "SL.glm"),
    tmle_candidate("dup_id", g_library = "SL.glm", truncation = 0.05)
  )
  expect_error(cleanTMLE:::validate_tmle_candidates(cands), "Duplicate")
})

test_that("expand_tmle_candidate_grid creates grid", {
  grid <- cleanTMLE:::expand_tmle_candidate_grid(
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


# ── New Phase 2 Tests: Residual Confounding Stage ────────────────────────

test_that("run_residual_confounding_stage works", {
  dat  <- sim_func1(n = 300, seed = 1)
  lock <- create_analysis_lock(dat, "treatment", "event_24",
                               c("age", "sex", "biomarker"), seed = 1)
  lock <- cleanTMLE:::define_negative_control(lock, "nc_outcome")
  ps   <- cleanTMLE:::fit_ps_glm(lock)
  stage3 <- cleanTMLE:::run_residual_confounding_stage(lock, ps)
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
  ps <- cleanTMLE:::fit_ps_glm(lock)
  expect_error(cleanTMLE:::run_residual_confounding_stage(lock, ps),
               "No negative controls")
})


# ── New Phase 2 Tests: Pre-Outcome Gate ──────────────────────────────────

test_that("authorize_outcome_analysis returns GO when all checkpoints pass", {
  dat  <- sim_func1(n = 300, seed = 1)
  lock <- create_analysis_lock(dat, "treatment", "event_24",
                               c("age", "sex", "biomarker"), seed = 1)
  lock <- cleanTMLE:::define_negative_control(lock, "nc_outcome")
  audit <- cleanTMLE:::create_audit_log(lock)

  cp1 <- cleanTMLE:::checkpoint_cohort_adequacy(lock)
  audit <- cleanTMLE:::record_checkpoint(audit, cp1)

  ps   <- cleanTMLE:::fit_ps_glm(lock)
  diag <- cleanTMLE:::compute_ps_diagnostics(ps)
  cp2  <- cleanTMLE:::checkpoint_balance(diag, lock_hash = lock$lock_hash)
  audit <- cleanTMLE:::record_checkpoint(audit, cp2)

  nc  <- cleanTMLE:::run_negative_control(lock, "nc_outcome", ps)
  cp3 <- cleanTMLE:::checkpoint_residual_bias(nc, lock_hash = lock$lock_hash)
  audit <- cleanTMLE:::record_checkpoint(audit, cp3)

  gate <- cleanTMLE:::authorize_outcome_analysis(audit)
  expect_s3_class(gate, "pre_outcome_gate")
  expect_s3_class(gate, "cleantmle_checkpoint")
  expect_true(gate$authorized)
  expect_equal(gate$decision, "GO")
})

test_that("authorize_outcome_analysis returns STOP when checkpoints missing", {
  dat  <- sim_func1(n = 100, seed = 1)
  lock <- create_analysis_lock(dat, "treatment", "event_24",
                               c("age", "sex"), seed = 1)
  audit <- cleanTMLE:::create_audit_log(lock)
  # No checkpoints recorded
  gate <- cleanTMLE:::authorize_outcome_analysis(audit)
  expect_false(gate$authorized)
  expect_equal(gate$decision, "STOP")
})

test_that("assert_outcome_authorized errors when not authorized", {
  dat  <- sim_func1(n = 100, seed = 1)
  lock <- create_analysis_lock(dat, "treatment", "event_24",
                               c("age", "sex"), seed = 1)
  audit <- cleanTMLE:::create_audit_log(lock)
  expect_error(cleanTMLE:::assert_outcome_authorized(audit), "NOT authorised")
})

test_that("unmask_outcome is a plain reversal since 0.2.0; the two-pass path stays strict", {
  dat  <- sim_func1(n = 150, seed = 4)
  lock <- create_analysis_lock(dat, "treatment", "event_24",
                               c("age", "sex", "biomarker"), seed = 4)
  masked <- mask_outcome(lock)
  # A plain cleanroom lock unmasks without an audit: masking is the honest
  # blinding device, and no token is demanded (0.2.0 contract).
  un <- unmask_outcome(masked, lock)
  expect_true(isTRUE(un$.outcome_authorized))
  expect_false(isTRUE(un$.outcome_masked))
  expect_identical(un$data[["event_24"]], lock$data[["event_24"]])
  # A lock that opted into enforcement (the two-pass path) still errors
  # without an audit, and forcing still warns.
  strict <- masked
  strict$require_authorization <- TRUE
  expect_error(unmask_outcome(strict, lock), "requires an .audit.")
  expect_warning(
    un2 <- unmask_outcome(strict, lock, allow_unauthorized = TRUE),
    "forced without an audit")
  expect_true(isTRUE(un2$.outcome_authorized))
})

test_that("unmask_outcome errors on a non-authorising gate unless forced (hard)", {
  dat  <- sim_func1(n = 100, seed = 7)
  lock <- create_analysis_lock(dat, "treatment", "event_24",
                               c("age", "sex"), seed = 7)
  masked <- mask_outcome(lock)
  audit  <- cleanTMLE:::create_audit_log(lock)          # no checkpoints -> gate STOP
  expect_error(unmask_outcome(masked, lock, audit = audit), "did NOT authorise")
  expect_warning(
    un <- unmask_outcome(masked, lock, audit = audit, allow_unauthorized = TRUE),
    "forced via allow_unauthorized")
  expect_true(isTRUE(un$.outcome_authorized))
})

test_that("unmask_outcome authorises through a passing gate (hard)", {
  dat  <- sim_func1(n = 300, seed = 9)
  lock <- create_analysis_lock(dat, "treatment", "event_24",
                               c("age", "sex", "biomarker"), seed = 9)
  lock <- cleanTMLE:::define_negative_control(lock, "nc_outcome")
  masked <- mask_outcome(lock)
  audit  <- cleanTMLE:::create_audit_log(lock)
  audit  <- cleanTMLE:::record_checkpoint(audit, cleanTMLE:::checkpoint_cohort_adequacy(lock))
  ps     <- cleanTMLE:::fit_ps_glm(lock)
  diag   <- cleanTMLE:::compute_ps_diagnostics(ps)
  audit  <- cleanTMLE:::record_checkpoint(audit,
              cleanTMLE:::checkpoint_balance(diag, lock_hash = lock$lock_hash))
  nc     <- cleanTMLE:::run_negative_control(lock, "nc_outcome", ps)
  audit  <- cleanTMLE:::record_checkpoint(audit,
              cleanTMLE:::checkpoint_residual_bias(nc, lock_hash = lock$lock_hash))

  un <- unmask_outcome(masked, lock, audit = audit)   # gate GO -> no error
  expect_true(isTRUE(un$.outcome_authorized))
  # A Stage 4 estimator now runs without an override.
  expect_silent(.check_outcome_access(un, caller = "test"))
})

test_that(".check_outcome_access: masking always enforced, authorization opt-in (0.2.0)", {
  dat  <- sim_func1(n = 100, seed = 8)
  lock <- create_analysis_lock(dat, "treatment", "event_24",
                               c("age", "sex"), seed = 8)
  # An unmasked cleanroom lock passes: no token is demanded by default.
  expect_silent(.check_outcome_access(lock, caller = "test"))
  # A masked lock is still refused.
  expect_error(.check_outcome_access(mask_outcome(lock), caller = "test"),
               "masked")
  # A lock that opted into enforcement is refused until authorised.
  strict <- lock
  strict$require_authorization <- TRUE
  expect_error(.check_outcome_access(strict, caller = "test"),
               "requires a recorded pre-outcome authorisation")
  expect_silent(.check_outcome_access(strict, allow_outcome_access = TRUE,
                                      caller = "test"))
  # A plain (cleanroom_enabled = FALSE) lock is exempt from everything.
  simple <- cleanTMLE:::create_simple_lock(dat, "treatment", "event_24",
                               c("age", "sex"), seed = 8)
  expect_silent(.check_outcome_access(simple, caller = "test"))
})


# ── New Phase 2 Tests: Decision Log ─────────────────────────────────────

test_that("record_decision_log_entry and export_decision_log work", {
  dat  <- sim_func1(n = 100, seed = 1)
  lock <- create_analysis_lock(dat, "treatment", "event_24",
                               c("age", "sex"), seed = 1)
  audit <- cleanTMLE:::create_audit_log(lock)
  audit <- cleanTMLE:::record_decision_log_entry(
    audit, "Stage 2", "model_specification",
    "Selected GLM PS model",
    rationale = "Pre-specified in SAP"
  )
  audit <- cleanTMLE:::record_decision_log_entry(
    audit, "Stage 2b", "override",
    "Accepted FLAG for ESS"
  )
  dl <- cleanTMLE:::export_decision_log(audit)
  expect_true(is.data.frame(dl))
  expect_equal(nrow(dl), 2)
  expect_true(all(c("stage", "decision_type", "description") %in% names(dl)))
})

test_that("export_decision_log returns empty df when no entries", {
  dat  <- sim_func1(n = 100, seed = 1)
  lock <- create_analysis_lock(dat, "treatment", "event_24",
                               c("age", "sex"), seed = 1)
  audit <- cleanTMLE:::create_audit_log(lock)
  dl <- cleanTMLE:::export_decision_log(audit)
  expect_true(is.data.frame(dl))
  expect_equal(nrow(dl), 0)
})


# ── New Phase 2 Tests: Stage Path Narrative ──────────────────────────────

# (summarize_stage_path moved to cleanroomGov, tested there.)


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
  unmasked <- unmask_outcome(masked, lock, allow_unauthorized = TRUE)
  expect_identical(unmasked$data[["event_24"]], lock$data[["event_24"]])
  expect_false(isTRUE(unmasked$.outcome_masked))
})

test_that("Stage 4 functions reject masked outcome", {
  dat  <- sim_func1(n = 200, seed = 1)
  lock <- create_analysis_lock(dat, "treatment", "event_24",
                               c("age", "sex"), seed = 1)
  masked <- mask_outcome(lock)
  expect_error(cleanTMLE:::run_crude_workflow(masked), "Outcome is masked")
})

test_that("Stage 4 functions allow override_clean_room", {
  dat  <- sim_func1(n = 200, seed = 1)
  lock <- create_analysis_lock(dat, "treatment", "event_24",
                               c("age", "sex"), seed = 1)
  # A non-masked but unauthorised cleanroom lock now requires the override.
  result <- cleanTMLE:::run_crude_workflow(lock, allow_outcome_access = TRUE)
  expect_true(is.numeric(result$estimate))
})
