fx <- fixture_pipeline()

test_that("the dossier walks the ladder and records the candidate", {
  d <- fx$dossier
  expect_s3_class(d, "cr_dossier")
  expect_identical(d$decision$primary, "ATE")
  expect_identical(d$decision$candidates$ATE$learners, "glm")
  expect_identical(d$decision$K, fx$plan$K)
  expect_identical(d$decision$fold_seed, fx$plan$seed)
  expect_identical(d$decision$nc_verdict, "GO")
  expect_null(d$design$g)
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
  st <- fixture_pipeline(nc = c(nc_bad = "care use"))
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
  st <- fixture_pipeline(nc = c(nc_bad = "care use"))
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
})
