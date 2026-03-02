test_that("sim_func1 generates correct structure", {
  dat <- sim_func1(n = 100, seed = 60)

  expect_true(is.data.frame(dat))
  expect_equal(nrow(dat), 100)
  expect_true(all(c("id", "age", "sex", "biomarker", "comorbidity",
                     "treatment", "time", "event", "event_type",
                     "censored") %in% names(dat)))
})

test_that("sim_func1 values are valid", {
  dat <- sim_func1(n = 500, seed = 61)

  # Treatment is binary
 expect_true(all(dat$treatment %in% c(0, 1)))

  # Event is binary
  expect_true(all(dat$event %in% c(0, 1)))

  # Event type is 0, 1, or 2
  expect_true(all(dat$event_type %in% c(0, 1, 2)))

  # Censored is binary
  expect_true(all(dat$censored %in% c(0, 1)))

  # Time is positive
  expect_true(all(dat$time > 0))

  # Censored + event should cover all
  # (censored = 1 when event_occurred = 0, i.e., event_type = 0)
  expect_true(all((dat$censored == 1 & dat$event_type == 0) |
                    (dat$censored == 0 & dat$event_type > 0)))
})

test_that("sim_func1 is reproducible with seed", {
  dat1 <- sim_func1(n = 50, seed = 62)
  dat2 <- sim_func1(n = 50, seed = 62)

  expect_identical(dat1, dat2)
})

test_that("sim_func1 respects n parameter", {
  dat <- sim_func1(n = 37, seed = 63)
  expect_equal(nrow(dat), 37)
})
