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
