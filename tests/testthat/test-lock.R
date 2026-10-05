test_that("the lock refuses outcome columns", {
  d <- make_design()
  bad <- merge(d$design, d$outcomes, by = "id")
  expect_error(create_analysis_lock(bad, "A", c("w1", "w2", "w3"), fast_plan()),
               "contains the outcome")
})

test_that("the lock keeps only the columns it needs and verifies", {
  d <- make_design()
  lk <- create_analysis_lock(d$design, "A", c("w1", "w2", "w3"),
                             fast_plan(negative_controls = c(nc_visit = "care"),
                                       restrictions = list(no_transfer = ~ transfer == 0)))
  expect_s3_class(lk, "cr_lock")
  expect_setequal(names(lk$data), c("id", "A", "w1", "w2", "w3", "nc_visit", "transfer"))
  expect_true(verify_lock(lk))
})

test_that("any change to the plan or the design data breaks the hash", {
  d <- make_design()
  lk <- create_analysis_lock(d$design, "A", c("w1", "w2", "w3"), fast_plan())
  lk2 <- create_analysis_lock(d$design, "A", c("w1", "w2", "w3"),
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
  expect_error(create_analysis_lock(x, "A", c("w1", "w2"), fast_plan()), "complete")
  x <- d$design; x$id[2] <- 1L
  expect_error(create_analysis_lock(x, "A", c("w1", "w2"), fast_plan()), "unique")
  x <- d$design; x$A[1] <- 2
  expect_error(create_analysis_lock(x, "A", c("w1", "w2"), fast_plan()), "0/1")
})

test_that("save and load round-trip, and a version change only warns", {
  d <- make_design()
  lk <- create_analysis_lock(d$design, "A", c("w1", "w2", "w3"), fast_plan())
  path <- withr::local_tempfile(fileext = ".rds")
  save_lock(lk, path)
  expect_identical(load_lock(path)$lock_hash, lk$lock_hash)
  local_mocked_bindings(.pkg_version = function() "9.9.9")
  expect_warning(verify_lock(lk), "created with cleanTMLE")
})

test_that("stamps detect a foreign lock and a modified object", {
  d <- make_design()
  lk <- create_analysis_lock(d$design, "A", c("w1", "w2", "w3"), fast_plan())
  lk2 <- create_analysis_lock(d$design, "A", c("w1", "w2", "w3"), fast_plan(seed = 99L))
  obj <- list(a = 1:3); obj$stamp <- .stamp(lk, obj)
  expect_true(.check_stamp(lk, obj, "obj"))
  expect_error(.check_stamp(lk2, obj, "obj"), "not computed from this lock")
  obj$a[1] <- 9L
  expect_error(.check_stamp(lk, obj, "obj"), "modified")
})
