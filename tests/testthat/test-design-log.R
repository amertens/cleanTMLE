test_that("the design log carries hashes, decision, approval and estimates", {
  fx <- fixture_pipeline()
  ub <- unblind(fx$lock, fx$dossier, fx$d$outcomes, approved_by = "Review team")
  fit <- estimate_effect(ub)
  f <- withr::local_tempfile(fileext = ".json")
  export_design_log(ub, fit, file = f)
  log <- jsonlite::read_json(f)
  expect_identical(log$lock_hash, fx$lock$lock_hash)
  expect_identical(log$dossier_hash, fx$dossier$dossier_hash)
  expect_identical(log$decision$primary, "ATE")
  expect_identical(log$approval$approved_by, "Review team")
  expect_equal(log$estimates[[1]]$estimate, fit$estimate, tolerance = 1e-10)
})

test_that("a dossier alone exports the pre-outcome log", {
  fx <- fixture_pipeline()
  log <- export_design_log(fx$dossier)
  expect_null(log$approval)
  expect_identical(log$decision$primary, "ATE")
})

test_that("the example dataset keeps the outcome separate and carries the truth", {
  expect_setequal(names(cr_example), c("design", "outcomes", "truth", "dictionary"))
  expect_false("death_30d" %in% names(cr_example$design))
  expect_identical(cr_example$design$id, cr_example$outcomes$id)
  expect_setequal(cr_example$truth$estimand, c("ATE", "trimmed_ATE", "ATT", "ATO"))
  expect_true(all(is.finite(cr_example$truth$truth)))
})

test_that("the package exports exactly the twelve verbs", {
  ns <- readLines(system.file("NAMESPACE", package = "cleanTMLE"))
  ex <- sub("^export\\((.*)\\)$", "\\1", grep("^export\\(", ns, value = TRUE))
  expect_setequal(ex, c("analysis_plan", "assess_design", "create_analysis_lock",
                        "design_report", "estimate_effect", "export_design_log",
                        "load_lock", "negative_control_ladder", "save_lock",
                        "simulate_design", "unblind", "verify_lock"))
})

test_that("the Muntner export has the Table S1 columns and one row per step", {
  fx <- fixture_pipeline()
  ub <- unblind(fx$lock, fx$dossier, fx$d$outcomes, approved_by = "Review team")
  fit <- estimate_effect(ub)
  f <- withr::local_tempfile(fileext = ".csv")
  out <- export_design_log(ub, fit, file = f, format = "muntner")
  expect_identical(names(out), c("date", "stage", "issue", "decision", "rationale", "decided_by"))
  expect_identical(names(utils::read.csv(f)), names(out))
  expect_setequal(unique(out$stage),
                  c("Lock", "Design assessment", "Outcome-free simulation", "Negative controls",
                    "Pre-outcome decision", "Approval", "Estimation"))
  expect_equal(sum(out$stage == "Outcome-free simulation"), length(fx$plan$estimands))
  # Every simulation row states the criterion, so an infeasible row does not read as met.
  expect_match(out$rationale[out$stage == "Outcome-free simulation"], "^Feasible requires ")
  expect_identical(out$decided_by[out$stage == "Approval"], "Review team")
  expect_true(all(is.na(out$decided_by[out$stage == "Lock"])))
  pre <- export_design_log(fx$dossier, format = "muntner")
  expect_false(any(pre$stage %in% c("Approval", "Estimation")))
})

test_that("the Muntner export names roles and records an override", {
  d <- make_design(400, seed = 2, strength = 0.8)
  p <- fast_plan(negative_controls = c(nc_bad = "care use"),
                 nc_criteria = list(null_band = c(-0.01, 0.01), rule = "point"),
                 tolerance = list(bias = 0.1, coverage = 0.5), reps = 10L,
                 roles = list(analyst = "A. Analyst"))
  lk <- lock_fx(d$design, "A", c("w1", "w2", "w3"), p)
  ds <- assess_design(lk)
  sim <- suppressMessages(simulate_design(lk, ds))
  ncl <- negative_control_ladder(lk, ds, sim)
  dos <- design_report(lk, ds, sim, ncl)
  skip_if_not(identical(dos$decision$nc_verdict, "STOP"))
  ub <- unblind(lk, dos, d$outcomes, approved_by = "Review team",
                override = "Proceeding to report the bias the negative control shows.")
  out <- export_design_log(ub, format = "muntner")
  expect_identical(out$decided_by[out$stage == "Lock"], "A. Analyst")
  ov <- out[out$issue == "Override", ]
  expect_equal(nrow(ov), 1L)
  expect_match(ov$rationale, "Proceeding to report")
})

test_that("the Muntner row names the estimand the negative controls were read on", {
  fx <- fixture_pipeline()
  prim <- fx$dossier$decision$primary
  out <- export_design_log(fx$dossier, format = "muntner")
  nc_row <- out[out$stage == "Negative controls", ]
  expect_equal(nrow(nc_row), 1L)
  expect_identical(nc_row$issue,
                   sprintf("Residual confounding check on the %s (the primary estimand)", prim))
  expect_identical(nc_row$decision, fx$dossier$decision$nc_verdict)
  # No feasible estimand: the row says so and does not report a verdict.
  sim <- fx$sim
  sim$verdict$status <- "infeasible"
  sim$stamp <- .stamp(fx$lock, sim)
  none <- design_report(fx$lock, fx$design, sim,
                        negative_control_ladder(fx$lock, fx$design, sim))
  out2 <- export_design_log(none, format = "muntner")
  nc2 <- out2[out2$stage == "Negative controls", ]
  expect_equal(nrow(nc2), 1L)
  expect_identical(nc2$issue, "Residual confounding check (no feasible estimand)")
  expect_false(nc2$decision %in% c("GO", "FLAG", "STOP"))
  # No declared controls: no row.
  nn <- fixture_pipeline(nc = NULL)
  expect_false("Negative controls" %in%
                 export_design_log(nn$dossier, format = "muntner")$stage)
})

test_that("the JSON log carries the by-rung verdicts by estimand and the estimand that gates", {
  fx <- fixture_pipeline()
  f <- withr::local_tempfile(fileext = ".json")
  export_design_log(fx$dossier, file = f)
  log <- jsonlite::read_json(f, simplifyVector = TRUE)
  expect_identical(log$decision$nc_estimand, fx$dossier$decision$primary)
  expect_identical(log$decision$nc_verdict, fx$dossier$decision$nc_verdict)
  br <- log$negative_controls$by_rung
  expect_true("estimand" %in% names(br))
  expect_identical(br$estimand, fx$nc$by_rung$estimand)
  expect_identical(br$verdict, fx$nc$by_rung$verdict)
  expect_identical(names(log$negative_controls$verdict), fx$plan$estimands)
  expect_identical(unlist(log$negative_controls$verdict), fx$nc$verdict)
})

test_that("the Muntner negative-control rationale names the rule", {
  fx <- fixture_pipeline()
  out <- export_design_log(fx$dossier, format = "muntner")
  expect_match(out$rationale[out$stage == "Negative controls"], "point-estimate rule", fixed = TRUE)
  expect_match(out$rationale[out$stage == "Negative controls"], "for the primary estimand",
               fixed = TRUE)
  expect_match(out$rationale[out$stage == "Negative controls"],
               "when any lies outside it", fixed = TRUE)
  fi <- fixture_pipeline(nc_rule = "interval")
  o2 <- export_design_log(fi$dossier, format = "muntner")
  expect_match(o2$rationale[o2$stage == "Negative controls"], "interval rule", fixed = TRUE)
  expect_match(o2$rationale[o2$stage == "Negative controls"],
               "GO when every control's 95% interval lies inside the null band", fixed = TRUE)
  expect_match(o2$rationale[o2$stage == "Negative controls"],
               "FLAG otherwise", fixed = TRUE)
  expect_match(o2$rationale[o2$stage == "Negative controls"], "-0.1", fixed = TRUE)
})

test_that("a dossier built before the estimand was recorded still exports the Muntner log", {
  fx <- fixture_pipeline()
  old <- fx$dossier
  old$decision$nc_estimand <- NULL
  old$dossier_hash <- .dossier_hash(old)
  out <- export_design_log(old, format = "muntner")
  nc_row <- out[out$stage == "Negative controls", ]
  expect_equal(nrow(nc_row), 1L)
  expect_identical(nc_row$issue, "Residual confounding check")
  expect_identical(nc_row$decision, old$decision$nc_verdict)
  expect_false(grepl("primary estimand", nc_row$rationale))
})

test_that("the rule text says how the controls combine into a verdict, for each rule", {
  band <- c(-0.02, 0.02)
  iv <- .nc_rule_text(list(rule = "interval", null_band = band))
  expect_match(iv, "GO when every control's 95% interval lies inside the null band", fixed = TRUE)
  expect_match(iv, "STOP when any lies entirely outside it", fixed = TRUE)
  expect_match(iv, "FLAG otherwise", fixed = TRUE)
  expect_match(iv, "overlaps the band edge", fixed = TRUE)
  expect_match(iv, "fewer estimable controls than the plan requires", fixed = TRUE)
  expect_match(iv, "FLAG does not block unblinding", fixed = TRUE)
  ci <- .nc_rule_text(list(rule = "ci", null_band = band))
  expect_match(ci, "GO when every control's whole 95% interval lies inside the null band", fixed = TRUE)
  expect_match(ci, "STOP when any control's interval is not entirely inside it", fixed = TRUE)
  expect_match(ci, "FLAG when no control is outside but fewer controls can be estimated", fixed = TRUE)
  pt <- .nc_rule_text(list(rule = "point", null_band = band))
  expect_match(pt, "GO when every control's point estimate lies inside the null band", fixed = TRUE)
  expect_match(pt, "STOP when any lies outside it", fixed = TRUE)
  expect_match(pt, "FLAG when none is outside but fewer controls can be estimated", fixed = TRUE)
})

test_that("each identification flag is a Muntner row", {
  ip <- id_pipeline()
  d <- ip$dossier
  log <- export_design_log(d, format = "muntner")
  expect_identical(sum(log$stage == "Identification"), 3L)
  k <- which(log$stage == "Lock")
  expect_identical(which(log$stage == "Identification"), k + 1:3)
  idr <- log[log$stage == "Identification", ]
  expect_true(any(grepl("r1", idr$rationale, fixed = TRUE)))
  fl <- d$identification[d$identification$flag == "flag", ]
  expect_identical(idr$issue, paste0(fl$component, ": ", fl$item))
  expect_identical(idr$decision, fl$declared)
  expect_identical(idr$rationale, fl$response)
  expect_true(all(idr$date == d$lock_created))
  expect_true(all(idr$decided_by == d$plan$roles$analyst))
  expect_identical(export_design_log(d)$identification, d$identification)
})

test_that("a dossier without an identification table has no Identification rows", {
  d <- fixture_pipeline()$dossier
  d$identification <- NULL
  d$dossier_hash <- cleanTMLE:::.dossier_hash(d)
  expect_false("Identification" %in% export_design_log(d, format = "muntner")$stage)
  expect_null(export_design_log(d)$identification)
})
