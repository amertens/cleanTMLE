# create_contrast_locks(), run_negative_control_ladder(),
# check_process_indicators().

.multiarm_data <- function(n = 900, seed = 101) {
  set.seed(seed)
  x1 <- stats::rnorm(n)
  transfer <- stats::rbinom(n, 1, 0.4)
  # Three transport modes; mode probabilities depend on x1.
  p_rco <- stats::plogis(-0.5 + 0.8 * x1)
  p_amb <- stats::plogis(-1 + 0.3 * x1)
  u <- stats::runif(n)
  mode <- ifelse(u < p_rco / 2, "rescueco",
                 ifelse(u < p_rco / 2 + p_amb / 2, "other_amb", "none"))
  # A negative-control covariate associated with mode ONLY through the
  # transfer subgroup: transfers are mostly urban and mostly other_amb.
  urban <- stats::rbinom(n, 1, ifelse(transfer == 1 & mode == "other_amb",
                                      0.95, 0.45))
  y <- stats::rbinom(n, 1, 0.15)
  data.frame(x1 = x1, transfer = transfer, urban = urban,
             mode = mode, outcome = y)
}

test_that("create_contrast_locks builds one lock per contrast with shared NCs", {
  dat <- .multiarm_data()
  locks <- create_contrast_locks(
    dat, "mode",
    contrasts = list(
      C1 = list(treated = "rescueco", control = c("other_amb", "none"),
                label = "rescueco vs everyone"),
      C4 = list(treated = "rescueco", control = "other_amb",
                label = "rescueco vs other ambulance")),
    outcome = "outcome", covariates = c("x1", "transfer"),
    negative_controls = "urban", cleanroom_enabled = FALSE)
  expect_s3_class(locks, "contrast_locks")
  expect_named(locks, c("C1", "C4"))
  n_rco <- sum(dat$mode == "rescueco")
  expect_equal(nrow(locks$C1$data), nrow(dat))
  expect_equal(sum(locks$C1$data$.A_contrast), n_rco)
  expect_equal(nrow(locks$C4$data),
               sum(dat$mode %in% c("rescueco", "other_amb")))
  expect_true("urban" %in% names(locks$C1$negative_controls))
  expect_true("urban" %in% names(locks$C4$negative_controls))
  expect_identical(locks$C4$contrast$label, "rescueco vs other ambulance")
  expect_true(any(locks$C1$design_log$type == "contrast"))
  expect_output(print(locks), "rescueco vs everyone")
  expect_error(create_contrast_locks(
    dat, "mode", contrasts = list(bad = list(treated = "nope",
                                             control = "none")),
    outcome = "outcome", covariates = "x1"), "not in mode")
})

test_that("the negative-control ladder finds where restriction removes the association", {
  dat <- .multiarm_data(2400, seed = 111)
  sub <- dat[dat$mode %in% c("rescueco", "other_amb"), ]
  sub$A <- as.integer(sub$mode == "rescueco")
  lock <- create_analysis_lock(sub, "A", "outcome", c("x1", "transfer"))
  lock <- cleanTMLE:::define_negative_control(lock, "urban")
  lad <- run_negative_control_ladder(
    lock,
    restrictions = list(`transfers excluded` = sub$transfer == 0),
    method = "unadjusted", verbose = FALSE)
  expect_s3_class(lad, "nc_ladder")
  tab <- lad$table
  full <- tab[tab$cohort == "full cohort" & tab$negative_control == "urban", ]
  rest <- tab[tab$cohort == "transfers excluded" &
                tab$negative_control == "urban", ]
  # Urban is associated with treatment through the transfer subgroup only.
  expect_true(full$flagged)
  expect_false(rest$flagged)
  expect_length(lad$turned_null, 1)
  expect_match(lad$turned_null, "restriction, not adjustment")
  expect_output(print(lad), "restriction ladder")
})

test_that("the negative-control ladder's TMLE method runs per rung", {
  dat <- .multiarm_data(600, seed = 121)
  sub <- dat[dat$mode %in% c("rescueco", "other_amb"), ]
  sub$A <- as.integer(sub$mode == "rescueco")
  lock <- create_analysis_lock(sub, "A", "outcome", c("x1",  "transfer"),
                             sl_library = "SL.glm")
  lock <- cleanTMLE:::define_negative_control(lock, "urban")
  lad <- run_negative_control_ladder(
    lock, restrictions = list(`transfers excluded` = sub$transfer == 0),
    method = "tmle", ps_method = "glm", verbose = FALSE)
  expect_true(all(lad$table$status == "estimated"))
  expect_true(all(is.finite(lad$table$estimate)))
})

test_that("check_process_indicators flags what treatment predicts", {
  set.seed(131)
  n <- 700
  x1 <- stats::rnorm(n)
  A <- stats::rbinom(n, 1, stats::plogis(0.5 * x1))
  vitals_missing <- stats::rbinom(n, 1, stats::plogis(-1 + 1.2 * A))
  noise_ind <- stats::rbinom(n, 1, 0.3)
  y <- stats::rbinom(n, 1, 0.2)
  dat <- data.frame(x1 = x1, A = A, vitals_missing = vitals_missing,
                    noise_ind = noise_ind, outcome = y)
  lock <- create_analysis_lock(dat, "A", "outcome",
                             c("x1", "vitals_missing", "noise_ind"))
  chk <- check_process_indicators(lock,
                                  indicators = c("vitals_missing",
                                                 "noise_ind"),
                                  condition_on = "x1")
  expect_true(chk$predicts[chk$indicator == "vitals_missing"])
  expect_false(chk$predicts[chk$indicator == "noise_ind"])
  expect_identical(chk$indicator[1], "vitals_missing")
})
