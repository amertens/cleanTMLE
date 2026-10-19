test_that("the lock refuses outcome columns", {
  d <- make_design()
  bad <- merge(d$design, d$outcomes, by = "id")
  expect_error(lock_fx(bad, "A", c("w1", "w2", "w3"), fast_plan()),
               "contains the outcome")
})

test_that("the lock keeps only the columns it needs and verifies", {
  d <- make_design()
  lk <- lock_fx(d$design, "A", c("w1", "w2", "w3"),
                             fast_plan(negative_controls = c(nc_visit = "care"),
                                       restrictions = list(no_transfer = ~ transfer == 0)))
  expect_s3_class(lk, "cr_lock")
  expect_setequal(names(lk$data), c("id", "A", "w1", "w2", "w3", "nc_visit", "transfer"))
  expect_true(verify_lock(lk))
})

test_that("any change to the plan or the design data breaks the hash", {
  d <- make_design()
  lk <- lock_fx(d$design, "A", c("w1", "w2", "w3"), fast_plan())
  lk2 <- lock_fx(d$design, "A", c("w1", "w2", "w3"),
                              fast_plan(tolerance = list(bias = 0.03)))
  expect_false(identical(lk$lock_hash, lk2$lock_hash))
  t1 <- lk; t1$plan$tolerance$bias <- 0.5
  expect_error(verify_lock(t1), "hash mismatch")
  t2 <- lk; t2$data$w1[1] <- 99
  expect_error(verify_lock(t2), "hash mismatch")
})

test_that("bad inputs are refused", {
  d <- make_design()
  x <- d$design; x$w1[2] <- NA
  expect_error(lock_fx(x, "A", c("w1", "w2"), fast_plan()), "complete")
  x <- d$design; x$id[2] <- 1L
  expect_error(lock_fx(x, "A", c("w1", "w2"), fast_plan()), "unique")
  x <- d$design; x$A[1] <- 2
  expect_error(lock_fx(x, "A", c("w1", "w2"), fast_plan()), "0/1")
})

test_that("save and load round-trip, and a version change only warns", {
  d <- make_design()
  lk <- lock_fx(d$design, "A", c("w1", "w2", "w3"), fast_plan())
  path <- withr::local_tempfile(fileext = ".rds")
  save_lock(lk, path)
  expect_identical(load_lock(path)$lock_hash, lk$lock_hash)
  local_mocked_bindings(.pkg_version = function() "9.9.9")
  expect_warning(verify_lock(lk), "created with cleanTMLE")
})

test_that("stamps detect a foreign lock and a modified object", {
  d <- make_design()
  lk <- lock_fx(d$design, "A", c("w1", "w2", "w3"), fast_plan())
  lk2 <- lock_fx(d$design, "A", c("w1", "w2", "w3"), fast_plan(seed = 99L))
  obj <- list(a = 1:3); obj$stamp <- .stamp(lk, obj)
  expect_true(.check_stamp(lk, obj, "obj"))
  expect_error(.check_stamp(lk2, obj, "obj"), "not computed from this lock")
  obj$a[1] <- 9L
  expect_error(.check_stamp(lk, obj, "obj"), "modified")
})

test_that("a captured variable or a function body in the plan changes the lock hash", {
  d <- make_design()
  lock_with <- function(plan) lock_fx(d$design, "A", c("w1", "w2", "w3"), plan)
  mk <- function(delta) function(a, X) stats::plogis(delta * a)
  h <- function(f) lock_with(fast_plan(surfaces = list(custom = list(s = f))))$lock_hash
  expect_identical(h(mk(0.5)), h(mk(0.5)))
  expect_false(identical(h(mk(0.5)), h(mk(0.6))))
  expect_false(identical(
    lock_with(fast_plan(select = function(m) "a"))$lock_hash,
    lock_with(fast_plan(select = function(m) "b"))$lock_hash))
})

test_that(".check_stamp recomputes the lock hash, so an edited lock fails", {
  d <- make_design()
  lk <- lock_fx(d$design, "A", c("w1", "w2", "w3"), fast_plan())
  obj <- list(a = 1:3); obj$stamp <- .stamp(lk, obj)
  lk$plan$tolerance$bias <- 0.5
  expect_error(.check_stamp(lk, obj, "obj"), "hash mismatch")
})

test_that("negative-control columns must be coded 0/1", {
  d <- make_design()
  covs <- c("w1", "w2", "w3")
  p <- fast_plan(negative_controls = c(nc_visit = "care use"))
  f <- d$design; f$nc_visit <- factor(ifelse(f$nc_visit == 1, "yes", "no"))
  expect_error(lock_fx(f, "A", covs, p), "must be coded 0/1: nc_visit")
  o <- d$design; o$nc_visit <- o$nc_visit + 1L
  expect_error(lock_fx(o, "A", covs, p), "must be coded 0/1: nc_visit")
  m <- d$design; m$nc_visit[1:5] <- NA
  expect_s3_class(lock_fx(m, "A", covs, p), "cr_lock")
})

test_that("roles are validated and hashed, and the creation time is not", {
  d <- make_design(100)
  expect_error(fast_plan(roles = list("someone")), "roles")
  p1 <- fast_plan(roles = list(analyst = "A. Analyst"))
  p2 <- fast_plan(roles = list(analyst = "B. Analyst"))
  l1 <- lock_fx(d$design, "A", c("w1", "w2", "w3"), p1)
  l2 <- lock_fx(d$design, "A", c("w1", "w2", "w3"), p2)
  expect_false(identical(l1$lock_hash, l2$lock_hash))
  expect_match(l1$created, "^\\d{4}-\\d{2}-\\d{2}T\\d{2}:\\d{2}:\\d{2}Z$")
  l1$created <- "2000-01-01T00:00:00Z"
  expect_true(verify_lock(l1))
})

test_that("a negative control that is also a covariate is refused", {
  d <- make_design()
  p <- fast_plan(negative_controls = c(nc_visit = "care use"))
  expect_error(lock_fx(d$design, "A", c("w1", "w2", "nc_visit"), p),
               "Negative control `nc_visit` is also an adjustment covariate.*nc_exclude")
  expect_s3_class(lock_fx(d$design, "A", c("w1", "w2"), p), "cr_lock")
})

test_that("nc_exclude values must be covariates and cannot remove every covariate", {
  d <- make_design()
  nc <- c(nc_visit = "care use")
  mk <- function(ex, covs = c("w1", "w2", "w3"))
    lock_fx(d$design, "A", covs,
                         fast_plan(negative_controls = nc, nc_exclude = list(nc_visit = ex)))
  expect_error(mk("nope"), "`nope`.*not an adjustment covariate")
  expect_error(mk(c("w1", "w2", "w3")), "every covariate")
  expect_error(mk("w1", covs = "w1"), "every covariate")
  expect_s3_class(mk(c("w1", "w2")), "cr_lock")
})

test_that("nc_exclude is part of the lock hash", {
  d <- make_design()
  mk <- function(...) lock_fx(d$design, "A", c("w1", "w2", "w3"),
                                           fast_plan(negative_controls = c(nc_visit = "care use"), ...))
  base <- mk()
  a <- mk(nc_exclude = list(nc_visit = "w1"))
  b <- mk(nc_exclude = list(nc_visit = "w2"))
  expect_false(identical(base$lock_hash, a$lock_hash))
  expect_false(identical(a$lock_hash, b$lock_hash))
  expect_identical(a$lock_hash, mk(nc_exclude = list(nc_visit = "w1"))$lock_hash)
  t1 <- a; t1$plan$nc_exclude$nc_visit <- "w2"
  expect_error(verify_lock(t1), "hash mismatch")
})

test_that("nc_exclude does not excuse a control that is also a covariate", {
  d <- make_design()
  p <- fast_plan(negative_controls = c(nc_visit = "care use"),
                 nc_exclude = list(nc_visit = "nc_visit"))
  expect_error(lock_fx(d$design, "A", c("w1", "w2", "nc_visit"), p),
               "also an adjustment covariate.*Remove it from `covariates` or from `negative_controls`")
})

# ---- variable dictionary ----------------------------------------------------

# Edit one row of a dictionary by column name.
dct_set <- function(dct, col, ...) {
  a <- list(...)
  for (n in names(a)) dct[[n]][dct$column == col] <- a[[n]]
  dct
}
# Add a row with before-time-zero defaults.
dct_add <- function(dct, column, ...) {
  row <- data.frame(column = column, variable = column, source = "test fixture",
                    follow_up = FALSE, recorded = "before", fixed = "before",
                    role = "covariate", reason = "", sensitivity = "",
                    stringsAsFactors = FALSE)
  a <- list(...)
  for (n in names(a)) row[[n]] <- a[[n]]
  rbind(dct, row)
}
dct_fx <- function(plan = fast_plan(negative_controls = c(nc_visit = "care use")),
                   design = make_design()$design, covs = c("w1", "w2", "w3")) {
  list(d = design, covs = covs, plan = plan,
       dct = fx_dictionary(design, "A", covs, plan))
}

test_that("the lock requires a dictionary", {
  fx <- dct_fx()
  expect_error(create_analysis_lock(fx$d, "A", fx$covs, fx$plan),
               "needs a variable `dictionary`")
  expect_error(create_analysis_lock(fx$d, "A", fx$covs, fx$plan, dictionary = NULL),
               "needs a variable `dictionary`")
  expect_s3_class(create_analysis_lock(fx$d, "A", fx$covs, fx$plan, dictionary = fx$dct),
                  "cr_lock")
})

test_that("rows must match the kept columns", {
  fx <- dct_fx()
  lk <- function(dct) create_analysis_lock(fx$d, "A", fx$covs, fx$plan, dictionary = dct)
  expect_error(lk(fx$dct[fx$dct$column != "w3", ]), "No dictionary row for column\\(s\\): w3")
  expect_error(lk(dct_add(fx$dct, "zz")), "not in `design_data`: zz")
  expect_error(lk(rbind(fx$dct, fx$dct[fx$dct$column == "w1", ])),
               "Duplicate dictionary rows for: w1")
  expect_error(lk(dct_set(fx$dct, "w1", fixed = "later")),
               "`fixed` must be one of before, at, after")
  expect_error(lk(dct_set(fx$dct, "w1", recorded = "never")),
               "`recorded` must be one of before, at, after")
  expect_error(lk(dct_set(fx$dct, "w1", role = "boss")), "`role` must be one of")
  s <- fx$dct; s$source[2] <- NA
  expect_error(lk(s), "`source`")
  expect_error(lk(fx$dct[, setdiff(names(fx$dct), "fixed")]), "`fixed`")
  expect_error(lk(as.list(fx$dct)), "data frame")
})

test_that("roles must match their use", {
  fx <- dct_fx()
  lk <- function(dct) create_analysis_lock(fx$d, "A", fx$covs, fx$plan, dictionary = dct)
  expect_error(lk(dct_set(fx$dct, "w2", role = "negative_control")),
               "`w2` is in `covariates` but its dictionary role is negative_control")
  expect_error(lk(dct_set(fx$dct, "id", role = "covariate")),
               "`id` is the id column but its dictionary role is covariate")
  expect_error(lk(dct_set(fx$dct, "A", role = "covariate")),
               "`A` is the treatment column but its dictionary role is covariate")
  # A covariate row for a column that is not a covariate.
  expect_error(lk(dct_set(fx$dct, "nc_visit", role = "covariate")),
               "`nc_visit` has dictionary role covariate but is not in `covariates`")
  expect_error(lk(dct_add(fx$dct, "nc_bad", role = "stand_in")),
               "`nc_bad` has dictionary role stand_in but is not in `covariates`")
})

test_that("a covariate fixed after time zero is refused", {
  fx <- dct_fx()
  lk <- function(dct) create_analysis_lock(fx$d, "A", fx$covs, fx$plan, dictionary = dct)
  after <- dct_set(fx$dct, "w1", fixed = "after")
  expect_error(lk(after), "Covariate `w1` is fixed after time zero")
  expect_error(lk(after), "after_treatment.*stand_in")
  si <- dct_set(after, "w1", role = "stand_in")
  expect_error(lk(si), "Stand-in `w1` needs a `reason` and a `sensitivity`")
  expect_error(lk(dct_set(si, "w1", reason = "Closest record of severity at the scene")),
               "Stand-in `w1` needs a `reason` and a `sensitivity`")
  ok <- dct_set(si, "w1", reason = "Closest record of severity at the scene",
                sensitivity = "Lock without w1")
  expect_s3_class(lk(ok), "cr_lock")
})

test_that("follow-up columns are refused unless they observe the outcome or are controls", {
  d <- make_design()$design
  d$fu_item <- withr::with_seed(3, stats::rbinom(nrow(d), 1, 0.5))
  plan <- fast_plan(negative_controls = c(nc_visit = "care use"))
  W <- c("w1", "w2", "w3")
  # As a covariate.
  dct <- dct_set(fx_dictionary(d, "A", c(W, "fu_item"), plan), "fu_item", follow_up = TRUE)
  expect_error(create_analysis_lock(d, "A", c(W, "fu_item"), plan, dictionary = dct),
               "`fu_item` comes from the follow-up source")
  # As an outcome-observation indicator.
  base <- fx_dictionary(d, "A", W, plan)
  obs <- dct_add(base, "fu_item", role = "observation", follow_up = TRUE, fixed = "after")
  expect_error(create_analysis_lock(d, "A", W, plan, dictionary = obs),
               "Outcome-observation indicator `fu_item` needs a `reason`")
  obs <- dct_set(obs, "fu_item", reason = "D3, 2026-10-18")
  lk <- create_analysis_lock(d, "A", W, plan, dictionary = obs)
  expect_true("fu_item" %in% names(lk$data))
  expect_false("fu_item" %in% lk$covariates)
  # A negative control from the follow-up source.
  nc <- dct_set(base, "nc_visit", follow_up = TRUE)
  expect_error(create_analysis_lock(d, "A", W, plan, dictionary = nc),
               "Negative control `nc_visit` comes from the follow-up source")
  nc <- dct_set(nc, "nc_visit", reason = "Pre-injury item from the interview")
  expect_s3_class(create_analysis_lock(d, "A", W, plan, dictionary = nc), "cr_lock")
})

test_that("after_treatment columns are kept but never adjusted for", {
  d <- make_design()$design
  d$arr_hr <- withr::with_seed(4, stats::rnorm(nrow(d)))
  d$arr_hr[c(2, 9)] <- NA
  plan <- fast_plan(negative_controls = c(nc_visit = "care use"))
  W <- c("w1", "w2", "w3")
  dct <- dct_add(fx_dictionary(d, "A", W, plan), "arr_hr",
                 role = "after_treatment", fixed = "after", recorded = "after")
  lk <- create_analysis_lock(d, "A", W, plan, dictionary = dct)
  expect_true("arr_hr" %in% names(lk$data))
  expect_false("arr_hr" %in% lk$covariates)
  d2 <- d; d2$arr_hr[is.na(d2$arr_hr)] <- 0
  expect_error(create_analysis_lock(d2, "A", c(W, "arr_hr"), plan, dictionary = dct),
               "`arr_hr` is a node after treatment")
  obs <- dct_set(dct, "arr_hr", role = "observation", reason = "approved")
  expect_error(create_analysis_lock(d2, "A", c(W, "arr_hr"), plan, dictionary = obs),
               "`arr_hr` is a node after treatment")
})

test_that("the treatment must be fixed at time zero", {
  fx <- dct_fx()
  expect_error(create_analysis_lock(fx$d, "A", fx$covs, fx$plan,
                                    dictionary = dct_set(fx$dct, "A", fixed = "after")),
               "The treatment `A` must be fixed at time zero")
  expect_error(create_analysis_lock(fx$d, "A", fx$covs, fx$plan,
                                    dictionary = dct_set(fx$dct, "A", fixed = "before")),
               "The treatment `A` must be fixed at time zero")
})

test_that("a restriction variable fixed after time zero needs a reason", {
  plan <- fast_plan(negative_controls = c(nc_visit = "care use"),
                    restrictions = list(direct = ~ transfer == 0))
  fx <- dct_fx(plan)
  expect_equal(fx$dct$role[fx$dct$column == "transfer"], "restriction")
  late <- dct_set(fx$dct, "transfer", fixed = "after")
  expect_error(create_analysis_lock(fx$d, "A", fx$covs, plan, dictionary = late),
               "Restriction variable `transfer` is fixed after time zero")
  late <- dct_set(late, "transfer", reason = "Depends on the route; see the bias analysis")
  expect_s3_class(create_analysis_lock(fx$d, "A", fx$covs, plan, dictionary = late),
                  "cr_lock")
})

test_that("a restriction variable that is also a covariate takes its covariate role", {
  plan <- fast_plan(negative_controls = c(nc_visit = "care use"),
                    restrictions = list(direct = ~ transfer == 0))
  d <- make_design()$design
  covs <- c("w1", "w2", "w3", "transfer")
  dct <- fx_dictionary(d, "A", covs, plan)
  expect_equal(dct$role[dct$column == "transfer"], "covariate")
  expect_s3_class(create_analysis_lock(d, "A", covs, plan, dictionary = dct), "cr_lock")
  # The restriction role would then name a column that is a covariate.
  expect_error(create_analysis_lock(d, "A", covs, plan,
                                    dictionary = dct_set(dct, "transfer", role = "restriction")),
               "`transfer` is in `covariates` but its dictionary role is restriction")
})

test_that("negative-control and restriction roles follow the plan", {
  plan <- fast_plan(negative_controls = c(nc_visit = "care use"),
                    restrictions = list(direct = ~ transfer == 0))
  fx <- dct_fx(plan)
  lk <- function(dct) create_analysis_lock(fx$d, "A", fx$covs, plan, dictionary = dct)
  expect_error(lk(dct_set(fx$dct, "nc_visit", role = "restriction")),
               "`nc_visit` is a negative control in the plan but its dictionary role is restriction")
  expect_error(lk(dct_set(fx$dct, "transfer", role = "negative_control")),
               "`transfer` has dictionary role negative_control but is not a negative control")
  expect_error(lk(dct_set(fx$dct, "transfer", role = "after_treatment", fixed = "after")),
               "`transfer` is a restriction variable in the plan but its dictionary role is after_treatment")
})

test_that("dummies of one variable must agree", {
  fx <- dct_fx()
  dct <- dct_set(fx$dct, "w2", variable = "v")
  dct <- dct_set(dct, "w3", variable = "v", fixed = "at")
  expect_error(create_analysis_lock(fx$d, "A", fx$covs, fx$plan, dictionary = dct),
               "Columns of variable `v` disagree on fixed")
  dct <- dct_set(dct, "w3", fixed = "before", recorded = "after")
  expect_error(create_analysis_lock(fx$d, "A", fx$covs, fx$plan, dictionary = dct),
               "Columns of variable `v` disagree on recorded")
  dct <- dct_set(dct, "w3", recorded = "before")
  expect_s3_class(create_analysis_lock(fx$d, "A", fx$covs, fx$plan, dictionary = dct),
                  "cr_lock")
  # One dummy a covariate and another a stand-in: the variable would split.
  mixed <- dct_set(dct, "w3", role = "stand_in", reason = "Closest record",
                   sensitivity = "Lock without w3")
  expect_error(create_analysis_lock(fx$d, "A", fx$covs, fx$plan, dictionary = mixed),
               "Columns of variable `v` disagree on role")
  # follow_up disagreement is named too.
  fu <- dct_set(dct, "w3", follow_up = TRUE)
  expect_error(create_analysis_lock(fx$d, "A", fx$covs, fx$plan, dictionary = fu),
               "Columns of variable `v` disagree on follow_up")
})

test_that("blank responses count as missing", {
  fx <- dct_fx()
  si <- dct_set(fx$dct, "w1", fixed = "after", role = "stand_in",
                reason = "   ", sensitivity = "Lock without w1")
  expect_error(create_analysis_lock(fx$d, "A", fx$covs, fx$plan, dictionary = si),
               "Stand-in `w1` needs a `reason` and a `sensitivity`")
  si <- dct_set(si, "w1", reason = "A real reason", sensitivity = " \t")
  expect_error(create_analysis_lock(fx$d, "A", fx$covs, fx$plan, dictionary = si),
               "Stand-in `w1` needs a `reason` and a `sensitivity`")
  na <- fx$dct; na$reason <- NA_character_
  lk <- create_analysis_lock(fx$d, "A", fx$covs, fx$plan, dictionary = na)
  expect_identical(lk$dictionary$reason, rep("", nrow(lk$dictionary)))
})

test_that("dictionary type and order do not change the lock", {
  fx <- dct_fx()
  mk <- function(dct) create_analysis_lock(fx$d, "A", fx$covs, fx$plan, dictionary = dct)$lock_hash
  base <- mk(fx$dct)
  expect_identical(mk(fx$dct[rev(seq_len(nrow(fx$dct))), ]), base)
  fac <- fx$dct
  fac[] <- lapply(fac, function(x) if (is.character(x)) factor(x) else x)
  expect_identical(mk(fac), base)
  chr <- fx$dct; chr$follow_up <- as.character(chr$follow_up)
  expect_identical(mk(chr), base)
  tbl <- structure(fx$dct, class = c("tbl_df", "tbl", "data.frame"))
  expect_identical(mk(tbl), base)
})

test_that("the lock stores the normalized dictionary", {
  fx <- dct_fx()
  shuffled <- fx$dct[rev(seq_len(nrow(fx$dct))), ]
  lk <- create_analysis_lock(fx$d, "A", fx$covs, fx$plan, dictionary = shuffled)
  expect_identical(names(lk$dictionary),
                   c("column", "variable", "source", "follow_up", "recorded", "fixed",
                     "role", "reason", "sensitivity"))
  expect_identical(lk$dictionary$column, sort(fx$dct$column, method = "radix"))
  expect_type(lk$dictionary$follow_up, "logical")
  expect_identical(class(lk$dictionary), "data.frame")
  expect_identical(attr(lk$dictionary, "row.names"), seq_len(nrow(fx$dct)))
})

test_that("the dictionary is hashed", {
  fx <- dct_fx()
  lk <- create_analysis_lock(fx$d, "A", fx$covs, fx$plan, dictionary = fx$dct)
  expect_true(verify_lock(lk))
  lk$dictionary$fixed[1] <- "after"
  expect_error(verify_lock(lk), "Lock hash mismatch")
})
