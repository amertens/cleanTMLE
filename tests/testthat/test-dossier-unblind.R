fx <- fixture_pipeline()
st <- fixture_pipeline(nc = c(nc_bad = "care use"))

test_that("the dossier walks the ladder and records the candidate", {
  d <- fx$dossier
  expect_s3_class(d, "cr_dossier")
  expect_identical(d$decision$primary, "ATE")
  expect_identical(d$decision$candidates$ATE$learners, "glm")
  expect_identical(d$decision$K, fx$plan$K)
  expect_identical(d$decision$fold_seed, fx$plan$seed)
  expect_identical(d$decision$nc_verdict, "GO")
  expect_true(all(unlist(d$decision$status) %in% c("feasible", "infeasible")))
  expect_identical(names(d$decision$unresolved), fx$plan$estimands)
  expect_true(all(vapply(d$decision$unresolved, is.logical, logical(1))))
  expect_identical(unname(unlist(d$decision$unresolved)), fx$sim$verdict$unresolved)
  expect_identical(d$simulation$decision_note, fx$sim$decision_note)
  expect_s3_class(d$simulation$extension, "data.frame")
  expect_identical(d$simulation$max_reps, fx$plan$max_reps)
  expect_null(d$design$g)
  expect_null(d$design$g_world)
  has_g_world <- function(x) is.list(x) &&
    ("g_world" %in% names(x) || any(vapply(x, has_g_world, logical(1))))
  expect_false(has_g_world(unclass(d)))
})

test_that("design_report refuses objects from another lock or design", {
  other <- fixture_pipeline(seed = 99L)
  expect_error(design_report(fx$lock, other$design, fx$sim, fx$nc), "not computed from this lock")
  expect_error(design_report(fx$lock, fx$design, other$sim, fx$nc), "not computed from this lock")
})

test_that("unblind checks the dossier and aligns outcomes by id", {
  ub <- unblind(fx$lock, fx$dossier, fx$d$outcomes, approved_by = "Review team")
  expect_s3_class(ub, "cr_unblinded")
  expect_identical(ub$outcomes$y, fx$d$outcomes$y)
  shuffled <- fx$d$outcomes[sample(nrow(fx$d$outcomes)), ]
  ub2 <- unblind(fx$lock, fx$dossier, shuffled, approved_by = "Review team")
  expect_identical(ub2$outcomes$y, fx$d$outcomes$y)
  partial <- fx$d$outcomes[-(1:5), ]
  ub3 <- unblind(fx$lock, fx$dossier, rbind(partial, data.frame(id = 9999, y = 1)),
                 approved_by = "Review team")
  expect_equal(ub3$approval$n_design_without_outcome, 5)
  expect_equal(ub3$approval$n_outcome_ids_not_in_design, 1)
})

test_that("a dossier from another lock, or an edited one, is rejected", {
  other <- fixture_pipeline(seed = 99L)
  expect_error(unblind(other$lock, fx$dossier, fx$d$outcomes, "Review team"),
               "different lock")
  ed <- fx$dossier
  ed$decision$primary <- "ATO"
  expect_error(unblind(fx$lock, ed, fx$d$outcomes, "Review team"), "modified")
  expect_error(unblind(fx$lock, fx$dossier, fx$d$outcomes, ""), "approved_by")
})

test_that("STOP blocks unblinding unless an override gives a reason", {
  expect_identical(st$dossier$decision$nc_verdict, "STOP")
  expect_error(unblind(st$lock, st$dossier, st$d$outcomes, "Review team"),
               "negative-control verdict is STOP")
  ub <- unblind(st$lock, st$dossier, st$d$outcomes, "Review team",
                override = "Control judged invalid by the clinical lead")
  expect_identical(ub$approval$override, "Control judged invalid by the clinical lead")
  expect_identical(ub$approval$blocked_by, "the negative-control verdict is STOP")
})

test_that("the gate is the primary estimand's verdict, and by_rung has a row per rung and estimand", {
  for (x in list(fx, st)) {
    d <- x$dossier
    expect_identical(d$decision$nc_estimand, d$decision$primary)
    expect_false(is.na(d$decision$primary))
    expect_identical(d$decision$nc_verdict, x$nc$verdict[[d$decision$primary]])
    expect_identical(names(x$nc$verdict), x$plan$estimands)
    rungs <- c("full cohort", names(x$plan$restrictions))
    expect_identical(nrow(x$nc$by_rung), length(rungs) * length(x$plan$estimands))
    expect_false(anyDuplicated(x$nc$by_rung[c("rung", "estimand")]) > 0)
    expect_identical(d$nc$by_rung, x$nc$by_rung)
  }
})

test_that("the gate follows the primary: a different primary reads a different verdict", {
  # Verdicts set by hand to differ by estimand; the dossier takes the primary's.
  prim <- fx$dossier$decision$primary
  for (v in c("STOP", "FLAG")) {
    nc <- fx$nc
    nc$verdict[] <- "GO"
    nc$verdict[[prim]] <- v
    nc$by_rung$verdict[nc$by_rung$rung == "full cohort"] <- unname(nc$verdict)
    nc$stamp <- .stamp(fx$lock, nc)
    d <- design_report(fx$lock, fx$design, fx$sim, nc)
    expect_identical(d$decision$nc_verdict, v)
    if (v == "STOP") {
      expect_identical(.blocked_reasons(d), "the negative-control verdict is STOP")
    } else {
      expect_length(.blocked_reasons(d), 0L)
    }
  }
  # A STOP on another estimand does not block the primary.
  other <- setdiff(fx$plan$estimands, prim)[1]
  nc <- fx$nc
  nc$verdict[] <- "GO"
  nc$verdict[[other]] <- "STOP"
  nc$by_rung$verdict[nc$by_rung$rung == "full cohort"] <- unname(nc$verdict)
  nc$stamp <- .stamp(fx$lock, nc)
  d <- design_report(fx$lock, fx$design, fx$sim, nc)
  expect_identical(d$decision$nc_verdict, "GO")
  expect_length(.blocked_reasons(d), 0L)
})

test_that("with no feasible estimand the nc verdict is NA and the plan is blocked by infeasibility", {
  sim <- fx$sim
  sim$verdict$status <- "infeasible"
  sim$stamp <- .stamp(fx$lock, sim)
  d <- design_report(fx$lock, fx$design, sim, negative_control_ladder(fx$lock, fx$design, sim))
  expect_true(is.na(d$decision$primary))
  expect_true(is.na(d$decision$nc_verdict))
  expect_true(is.na(d$decision$nc_estimand))
  expect_identical(.blocked_reasons(d), "no estimand is feasible")
})

test_that("a ladder made before the estimand column is refused", {
  old <- fx$nc
  old$by_rung$estimand <- NULL
  old$stamp <- .stamp(fx$lock, old)
  msg <- "This negative-control ladder was made by an earlier cleanTMLE; rerun negative_control_ladder()."
  expect_error(design_report(fx$lock, fx$design, fx$sim, old), msg, fixed = TRUE)
  old2 <- fx$nc
  old2$by_rung$estimand <- NULL
  expect_error(design_report(fx$lock, fx$design, fx$sim, old2), msg, fixed = TRUE)
})

test_that("a plan with negative controls cannot produce a dossier without nc", {
  expect_error(design_report(fx$lock, fx$design, fx$sim), "negative controls")
})

test_that("NA and blank override or approved_by do not count", {
  for (bad in list(NA_character_, "   ", ""))
    expect_error(unblind(st$lock, st$dossier, st$d$outcomes, "Review team", override = bad),
                 "negative-control verdict is STOP")
  for (bad in list(NA_character_, "   "))
    expect_error(unblind(fx$lock, fx$dossier, fx$d$outcomes, bad), "approved_by")
})

test_that("unblind errors when no outcome id matches a design id", {
  o <- fx$d$outcomes
  o$id <- sprintf("%04d", o$id)
  expect_error(unblind(fx$lock, fx$dossier, o, "Review team"), "No outcome id matches")
})

test_that("the dossier renders to HTML", {
  skip_on_cran()
  skip_if_not_installed("quarto")
  skip_if(is.null(quarto::quarto_path()))
  f <- withr::local_tempfile(fileext = ".html")
  design_report(fx$lock, fx$design, fx$sim, fx$nc, file = f)
  expect_true(file.exists(f))
  html <- paste(readLines(f, warn = FALSE), collapse = "")
  expect_match(html, "Decision")
  expect_match(html, "verdict")
  expect_match(html, "coverage_lower")
  # The negative controls are read on the primary estimand, and the section says so.
  expect_match(html, "(the primary estimand) gates unblinding", fixed = TRUE)
  expect_match(html, "rung</th>[^/]*estimand</th>[^/]*verdict</th>")
  # The section states the grading rule in one sentence (fx uses the point rule).
  expect_match(html, "point estimate against the null band", fixed = TRUE)
  # Quarto may turn the apostrophe into a typographic one.
  expect_match(html, "every control.s point estimate lies inside the null band")
  # The primary estimand's full-cohort controls are listed with their intervals and readings.
  expect_match(html, "Full-cohort controls on the ATE", fixed = TRUE)
  expect_match(html, "ci_upper</th>[^/]*reading</th>")
  expect_match(html, "nc_visit")
  expect_match(html, "judged on point estimates against the declared tolerances")
  # The trimmed-patient counts table appears, with its lead-in and the
  # covariate table that follows it.
  expect_match(html, "Patients outside the trim band, by arm", fixed = TRUE)
  expect_match(html, "group</th>.*removed</th>.*kept</th>.*treated</td>")
  expect_match(html, "covariate</th>.*std_diff</th>")
})

test_that("a dossier for a plan with no trim band builds and renders", {
  d <- make_design(300, seed = 2)
  plan <- fast_plan(estimands = c("ATE", "ATO"), trim_band = NULL, reps = 10L)
  lock <- lock_fx(d$design, "A", c("w1", "w2", "w3"), plan)
  design <- assess_design(lock)
  sim <- suppressMessages(simulate_design(lock, design))
  dos <- design_report(lock, design, sim, NULL)
  expect_null(dos$design$trimmed_profile)
  expect_false("trimmed_ATE" %in% names(dos$design$populations))
  expect_true(is.data.frame(dos$design$populations))
  skip_on_cran()
  skip_if_not_installed("quarto")
  skip_if(is.null(quarto::quarto_path()))
  f <- withr::local_tempfile(fileext = ".html")
  design_report(lock, design, sim, NULL, file = f)
  html <- paste(readLines(f, warn = FALSE), collapse = "")
  expect_match(html, "Who each estimand describes")
  expect_match(html, "No trim band was declared.")
  expect_no_match(html, "Means of fewer than 5 patients")
})

test_that("the dossier states the interval rule in its negative-controls section", {
  skip_on_cran()
  skip_if_not_installed("quarto")
  skip_if(is.null(quarto::quarto_path()))
  iv <- fixture_pipeline(nc_rule = "interval")
  f <- withr::local_tempfile(fileext = ".html")
  design_report(iv$lock, iv$design, iv$sim, iv$nc, file = f)
  html <- paste(readLines(f, warn = FALSE), collapse = "")
  expect_match(html, "GO when every control.s interval lies inside the null band")
  expect_match(html, "STOP when any lies entirely outside it", fixed = TRUE)
  expect_match(html, "FLAG does not block unblinding", fixed = TRUE)
  expect_match(html, "mixes extrapolation error with confounding", fixed = TRUE)
  expect_match(html, "ci_upper</th>[^/]*reading</th>")
})

test_that("a dossier carries a non-empty extension table and renders it", {
  # Ten repetitions cannot resolve a 0.8 coverage tolerance (the Wilson
  # interval for 10 of 10 starts near 0.72), so the selected cells are extended.
  ex <- fixture_pipeline(tolerance = list(bias = 0.1, coverage = 0.8), max_reps = 20L)
  d <- ex$dossier
  expect_true(nrow(d$simulation$extension) > 0)
  expect_identical(d$simulation$extension, ex$sim$extension)
  expect_true(any(d$simulation$metrics$reps > 10L))
  expect_true(all(d$simulation$metrics$reps <= 20L))
  expect_match(d$simulation$decision_note, "received extra repetitions up to max_reps = 20")
  skip_on_cran()
  skip_if_not_installed("quarto")
  skip_if(is.null(quarto::quarto_path()))
  f <- withr::local_tempfile(fileext = ".html")
  design_report(ex$lock, ex$design, ex$sim, ex$nc, file = f)
  html <- paste(readLines(f, warn = FALSE), collapse = "")
  expect_match(html, "reps_from")
  expect_match(html, "cells_unresolved_before")
})

test_that("unblind fingerprints the dossier, the approval and the outcomes", {
  ub <- unblind(fx$lock, fx$dossier, fx$d$outcomes, approved_by = "Review team")
  expect_identical(ub$unblind_hash,
                   .hash(list(fx$dossier$dossier_hash, ub$approval, ub$outcomes)))
  expect_true(.verify_unblinded(ub))
})

test_that("a hand-built unblinded object past STOP is refused", {
  ub <- unblind(st$lock, st$dossier, st$d$outcomes, "Review team",
                override = "Control judged invalid by the clinical lead")
  expect_true(.verify_unblinded(ub))
  # Built by hand with a consistent hash but no override: the gate is recomputed.
  hb <- structure(list(lock = st$lock, dossier = st$dossier, outcomes = ub$outcomes,
                       approval = utils::modifyList(ub$approval,
                                                    list(override = NA_character_))),
                  class = "cr_unblinded")
  hb$unblind_hash <- .unblind_hash(hb)
  expect_error(estimate_effect(hb, method = "weighting", estimand = "ATE"),
               "negative-control verdict is STOP")
  expect_error(export_design_log(hb), "negative-control verdict is STOP")
  # Built by hand with no hash at all.
  hb$unblind_hash <- NULL
  expect_error(estimate_effect(hb, method = "weighting", estimand = "ATE"), "modified")
})

test_that("editing the override or an outcome after unblind is detected", {
  ub <- unblind(st$lock, st$dossier, st$d$outcomes, "Review team",
                override = "Control judged invalid by the clinical lead")
  e1 <- ub
  e1$approval$override <- "A different reason"
  expect_error(estimate_effect(e1, method = "weighting", estimand = "ATE"), "modified")
  expect_error(export_design_log(e1), "modified")
  e2 <- ub
  e2$outcomes$y[1] <- 1 - e2$outcomes$y[1]
  expect_error(estimate_effect(e2, method = "weighting", estimand = "ATE"), "modified")
  expect_error(export_design_log(e2), "modified")
})

test_that("unblind checks the coding of a binary outcome", {
  o <- fx$d$outcomes
  cont <- o
  cont$y <- withr::with_seed(1, stats::runif(nrow(o)))
  expect_error(unblind(fx$lock, fx$dossier, cont, "Review team"), "coded 0/1")
  two <- o
  two$y <- o$y + 1L
  expect_error(unblind(fx$lock, fx$dossier, two, "Review team"), "coded 0/1")
  chr <- o
  chr$y <- as.character(o$y)
  expect_error(unblind(fx$lock, fx$dossier, chr, "Review team"), "coded 0/1")
  miss <- o
  miss$y[1:3] <- NA
  expect_s3_class(unblind(fx$lock, fx$dossier, miss, "Review team"), "cr_unblinded")
})

test_that("unblind checks the coding of a time-to-event outcome", {
  d <- make_design(300, seed = 9)
  p <- fast_plan(outcome = c(time = "t", event = "s"), estimands = "ATE",
                 target_time = 300, hazards = list(`0` = list(Surv(t, s == 0) ~ .),
                                                   `1` = list(Surv(t, s == 1) ~ .)),
                 candidates = list(library = list(glm = "glm"), truncation = 0.01),
                 surfaces = list(forms = "linear", heterogeneity = 0), reps = 2L)
  lk <- lock_fx(d$design, "A", c("w1", "w2", "w3"), p)
  ds <- assess_design(lk)
  dr <- design_report(lk, ds, suppressMessages(simulate_design(lk, ds)))
  oc <- data.frame(id = d$design$id, t = seq(10, 400, length.out = 300),
                   s = rep(0:1, 150))
  ok <- function(o) unblind(lk, dr, o, "Review team", override = "fixture")
  expect_s3_class(ok(oc), "cr_unblinded")
  neg <- oc
  neg$t[1] <- -5
  expect_error(ok(neg), "greater than 0")
  zero <- oc
  zero$t[2] <- 0
  expect_error(ok(zero), "greater than 0")
  frac <- oc
  frac$s[3] <- 0.5
  expect_error(ok(frac), "whole numbers")
  negev <- oc
  negev$s[3] <- -1
  expect_error(ok(negev), "whole numbers")
})

test_that("the dossier carries the lock's creation time and the plan table", {
  fx <- fixture_pipeline()
  expect_identical(fx$dossier$lock_created, fx$lock$created)
  expect_true(all(c("Estimand ladder (in order)", "Tolerances", "Outcome surfaces") %in%
                    fx$dossier$plan_table$setting))
  bad <- fx$dossier
  bad$lock_created <- "2000-01-01T00:00:00Z"
  expect_error(.verify_dossier(bad), "modified")
})

test_that("a lock without a creation time gives a dossier with NA", {
  fx <- fixture_pipeline()
  old <- fx$lock
  old$created <- NULL
  expect_true(verify_lock(old))
  d <- design_report(old, fx$design, fx$sim, fx$nc)
  expect_true(is.na(d$lock_created))
})

test_that("objects from cleanTMLE 0.4.0 stop with a clear message, not wrong output", {
  fx <- fixture_pipeline()
  msg_design <- "made by an earlier cleanTMLE; rerun assess_design"
  strip <- function(x, nm) { x[nm] <- NULL; x }
  for (nm in c("A", "c_statistic")) {
    old <- strip(fx$design, nm)
    expect_error(simulate_design(fx$lock, old), msg_design)
    expect_error(design_report(fx$lock, old, fx$sim, fx$nc), msg_design)
    expect_error(print(old), msg_design)
    expect_error(plot(old), msg_design)
  }
  expect_error(design_report(fx$lock, strip(fx$design, "populations"), fx$sim, fx$nc),
               msg_design)
  old_sim <- strip(fx$sim, "expected_events")
  expect_error(design_report(fx$lock, fx$design, old_sim, fx$nc), "rerun simulate_design")
  # A current design and simulation still work.
  expect_s3_class(design_report(fx$lock, fx$design, fx$sim, fx$nc), "cr_dossier")
})

test_that("the plan table prints values without padding", {
  expect_identical(.flat(c("ATE", "trimmed_ATE")), "ATE, trimmed_ATE")
  expect_identical(.flat(c(-0.02, 0.02)), "-0.02, 0.02")
  expect_identical(.flat(c(1, 0.5)), "1, 0.5")
  expect_identical(.flat(NULL), "none")
  expect_identical(.flat(list(a = c("x", "y"), b = 2)), "a: x, y; b: 2")
})

test_that("every recorded time is UTC with a trailing Z", {
  utc <- "^\\d{4}-\\d{2}-\\d{2}T\\d{2}:\\d{2}:\\d{2}Z$"
  fx <- fixture_pipeline()
  expect_match(fx$lock$created, utc)
  expect_match(fx$dossier$created, utc)
  ub <- unblind(fx$lock, fx$dossier, fx$d$outcomes, approved_by = "Review team")
  expect_match(ub$approval$approved_at, utc)
})

test_that("the dossier records the covariates left out of a negative control", {
  fx2 <- fixture_pipeline(nc_exclude = list(nc_visit = "w1"))
  expect_identical(unique(fx2$dossier$nc$table$excluded_covariates), "w1")
  row_label <- "Covariates left out of each negative control's adjustment"
  pt <- fx2$dossier$plan_table
  expect_identical(pt$value[pt$setting == row_label], "nc_visit: w1")
  none <- fixture_pipeline()$dossier$plan_table
  expect_identical(none$value[none$setting == row_label], "none")
  skip_on_cran()
  skip_if_not_installed("quarto")
  skip_if(is.null(quarto::quarto_path()))
  f <- withr::local_tempfile(fileext = ".html")
  design_report(fx2$lock, fx2$design, fx2$sim, fx2$nc, file = f)
  html <- paste(readLines(f, warn = FALSE), collapse = "")
  expect_match(html, "near zero by construction", fixed = TRUE)
  expect_match(html, "excluded_covariates</th>")
  # A plan with no exclusions does not carry the sentence.
  g <- withr::local_tempfile(fileext = ".html")
  design_report(fx$lock, fx$design, fx$sim, fx$nc, file = g)
  expect_no_match(paste(readLines(g, warn = FALSE), collapse = ""),
                  "near zero by construction")
})

test_that("design_report refuses a control that is also a covariate, even under a valid hash", {
  lk <- fx$lock
  lk$covariates <- c(lk$covariates, "nc_visit")
  lk$lock_hash <- .lock_hash(lk)
  expect_error(design_report(lk, fx$design, fx$sim, fx$nc),
               "Negative control `nc_visit` is also an adjustment covariate")
})

test_that("the dossier carries the identification table", {
  ip <- id_pipeline()
  d <- ip$dossier
  idt <- d$identification
  expect_identical(names(idt), c("component", "item", "declared", "flag", "response"))
  expect_true(all(vapply(idt, is.character, logical(1))))
  expect_true(all(idt$flag %in% c("flag", "note", "")))
  expect_identical(
    unique(idt$component),
    c("Time zero", "Start of follow-up", "Eligibility", "Treatment",
      "Adjustment set: fixed before time zero", "Adjustment set: stand-in",
      "Node after treatment", "Negative control"))
  expect_identical(sum(idt$flag == "flag"), 3L)
  tz <- idt[idt$component == "Time zero", ]
  expect_identical(tz$declared, "transport decision (measured as time of injury)")
  fu <- idt[idt$component == "Start of follow-up", ]
  expect_identical(c(fu$flag, fu$response), c("flag", "r0"))
  el <- idt[idt$component == "Eligibility", ]
  expect_identical(c(el$flag, el$response), c("flag", "r1"))
  si <- idt[idt$component == "Adjustment set: stand-in", ]
  expect_identical(si$response, "r2; sensitivity: s2")
  adj <- idt[idt$component == "Adjustment set: fixed before time zero", ]
  expect_identical(adj$item, "2 variables")
  expect_match(adj$declared, "w1", fixed = TRUE)
  expect_match(adj$response, "w1: x", fixed = TRUE)
  expect_identical(adj$flag, "note")
  expect_identical(idt$flag[idt$component == "Node after treatment"], "note")
  expect_identical(idt$flag[idt$component == "Treatment"], "")
  expect_identical(d$dossier_hash, cleanTMLE:::.dossier_hash(d))
  # The table is part of the hash: editing a response is detected.
  ed <- d
  ed$identification$response[1] <- "changed"
  expect_error(unblind(ip$lock, ed, ip$d$outcomes, "Review team"), "modified")
})

test_that("a plan with no flags has an identification table without flags", {
  idt <- fx$dossier$identification
  expect_identical(sum(idt$flag == "flag"), 0L)
  expect_identical(idt$declared[idt$component == "Time zero"], "treatment assignment")
  expect_identical(idt$declared[idt$component == "Adjustment set: fixed before time zero"],
                   "all recorded at or before time zero")
  expect_false("Eligibility" %in% idt$component)
})

test_that("the rendered dossier shows the identification section", {
  skip_on_cran()
  skip_if_not_installed("quarto")
  skip_if(is.null(quarto::quarto_path()))
  ip <- id_pipeline()
  f <- withr::local_tempfile(fileext = ".html")
  design_report(ip$lock, ip$design, ip$sim, ip$nc, file = f)
  expect_match(paste(readLines(f, warn = FALSE), collapse = " "),
               "Identification: what the plan declares", fixed = TRUE)
})

test_that("a plan restriction on a variable fixed after time zero is flagged, whatever its role", {
  # The table is built from the lock alone, so no pipeline is needed.
  idt_for <- function(var, edit) {
    d <- make_design()$design
    plan <- fast_plan(negative_controls = c(nc_visit = "care use"),
                      restrictions = list(direct = stats::as.formula(paste("~", var, ">= -100"))))
    covs <- c("w1", "w2", "w3")
    lk <- lock_fx(d, "A", covs, plan,
                  dictionary = edit(fx_dictionary(d, "A", covs, plan)))
    cleanTMLE:::.identification_table(lk)
  }
  # A stand-in that is also a restriction variable.
  idt <- idt_for("w3", function(dict) {
    i <- dict$column == "w3"
    dict$role[i] <- "stand_in"; dict$fixed[i] <- "after"
    dict$reason[i] <- "r2"; dict$sensitivity[i] <- "s2"
    dict
  })
  r <- idt[idt$component == "Restriction variable", ]
  expect_identical(nrow(r), 1L)
  expect_identical(c(r$item, r$flag, r$response), c("w3", "flag", "r2"))
  expect_identical(sum(idt$component == "Adjustment set: stand-in"), 1L)
  # A negative control that is also a restriction variable.
  idt <- idt_for("nc_visit", function(dict) {
    i <- dict$column == "nc_visit"
    dict$fixed[i] <- "after"; dict$reason[i] <- "r3"
    dict
  })
  r <- idt[idt$component == "Restriction variable", ]
  expect_identical(nrow(r), 1L)
  expect_identical(c(r$flag, r$response), c("flag", "r3"))
  expect_identical(idt$flag[idt$component == "Negative control"], "note")
  # A restriction variable fixed before time zero is not flagged.
  idt <- idt_for("w3", identity)
  expect_false("Restriction variable" %in% idt$component)
  # A role-restriction variable fixed after gets one row, not two.
  idt <- idt_for("transfer", function(dict) {
    i <- dict$column == "transfer"
    dict$fixed[i] <- "after"; dict$reason[i] <- "r4"
    dict
  })
  r <- idt[idt$component == "Restriction variable", ]
  expect_identical(c(nrow(r), r$response), c("1", "r4"))
})
