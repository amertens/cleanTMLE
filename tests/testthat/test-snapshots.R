# Snapshot tests for the two documents a review team reads: the printed
# design report and the exported Muntner-format decision log. The
# snapshots pin the exact wording and layout, so a drive-by change to
# either release document shows up as a reviewable diff under
# tests/testthat/_snaps/. Timestamps are scrubbed before comparison.

.scrub_times <- function(lines) {
  lines <- gsub("\\d{4}-\\d{2}-\\d{2} \\d{2}:\\d{2}:\\d{2}",
                "<timestamp>", lines)
  gsub("\\d{4}-\\d{2}-\\d{2}", "<date>", lines)
}

.snapshot_fixture <- function() {
  dat  <- sim_func1(n = 400, seed = 42)
  lock <- create_analysis_lock(
    dat, "treatment", "event_24", c("age", "sex", "biomarker"),
    seed = 7, negative_controls = "nc_outcome",
    roles = list(analyst = "programmer A", review_team = "panel B"))
  lock <- declare_estimand_ladder(lock, primary = "ATE",
                                  fallbacks = c("trimmed_ATE", "ATT"),
                                  trigger = "SEVERE")
  ps   <- fit_ps(lock, method = "glm")
  list(lock = lock, ps = ps)
}

test_that("the printed design report is stable (snapshot)", {
  fx   <- .snapshot_fixture()
  sup  <- assess_support(fx$ps)
  feas <- estimand_feasibility(fx$ps)
  rep  <- design_report(fx$lock, sup, feas)
  expect_snapshot(print(rep), transform = .scrub_times)
})

test_that("the Muntner-format decision log export is stable (snapshot)", {
  fx <- .snapshot_fixture()
  lock <- fx$lock
  lock <- declare_negative_controls(lock, "nc_outcome",
                                    domains = "health_seeking_behavior")
  masked <- mask_outcome(lock)
  lock2  <- unmask_outcome(masked, lock, approved_by = "review team")
  log <- export_design_log(lock2, format = "muntner")
  expect_named(log, c("date", "stage", "issue", "decision",
                      "rationale", "decided_by"))
  expect_snapshot(print(log, row.names = FALSE),
                  transform = .scrub_times)
})