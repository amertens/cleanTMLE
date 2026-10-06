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
  expect_match(html, "judged on point estimates against the declared tolerances")
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
  lk <- create_analysis_lock(d$design, "A", c("w1", "w2", "w3"), p)
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
