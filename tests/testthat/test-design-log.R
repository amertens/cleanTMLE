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

test_that("the example dataset keeps the outcome separate", {
  expect_setequal(names(cr_example), c("design", "outcomes"))
  expect_false("death_30d" %in% names(cr_example$design))
  expect_identical(cr_example$design$id, cr_example$outcomes$id)
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
  expect_identical(out$decided_by[out$stage == "Approval"], "Review team")
  expect_true(all(is.na(out$decided_by[out$stage == "Lock"])))
  pre <- export_design_log(fx$dossier, format = "muntner")
  expect_false(any(pre$stage %in% c("Approval", "Estimation")))
})

test_that("the Muntner export names roles and records an override", {
  d <- make_design(400, seed = 2, strength = 0.8)
  p <- fast_plan(negative_controls = c(nc_bad = "care use"),
                 nc_criteria = list(null_band = c(-0.01, 0.01)),
                 tolerance = list(bias = 0.1, coverage = 0.5), reps = 10L,
                 roles = list(analyst = "A. Analyst"))
  lk <- create_analysis_lock(d$design, "A", c("w1", "w2", "w3"), p)
  ds <- assess_design(lk)
  sim <- suppressMessages(simulate_design(lk, ds))
  ncl <- negative_control_ladder(lk, ds)
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
