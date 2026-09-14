# Blinding invariants for the design stage (revision WP0).
#
# Contract under test: every design-stage function must return the same
# output whether or not the primary outcome column is readable. Two
# sentinels are used against a reference lock:
#   (i)  a masked lock (mask_outcome(); outcome column all NA), and
#   (ii) a lock whose outcome column is randomly permuted (this preserves
#        the marginal outcome distribution, so functions whose contract
#        allows marginal outcome summaries must also be invariant).
# Outputs are compared after stripping volatile or data-carrying fields
# (calls, timestamps, embedded data frames, fitted model objects); the
# lock fingerprint is deliberately NOT stripped, so a future fingerprint
# that commits to outcome values will surface here.
#
# Two known Stage 0 violations are marked with skip() and an explicit
# WP1 pointer; WP1 removes the skips when it fixes the functions.

strip_volatile <- function(x) {
  if (is.environment(x) || is.function(x)) return(NULL)
  if (inherits(x, "formula")) return(paste(deparse(x), collapse = " "))
  if (is.data.frame(x)) {
    x$timestamp <- NULL
    return(x)
  }
  if (is.list(x)) {
    drop <- c("call", "data", "glm_fit", "sl_fit", "fit", "fits",
              "locked_at", "created_at", "declared_at", "saved_at")
    nm <- names(x)
    if (!is.null(nm)) x <- x[!(nm %in% drop)]
    if (length(x)) x[] <- lapply(x, strip_volatile)
    return(x)
  }
  x
}

expect_blind_identical <- function(a, b) {
  expect_equal(strip_volatile(a), strip_volatile(b), tolerance = 1e-12)
}

blinding_fixtures <- function(n = 400L) {
  dat <- sim_func1(n = n, seed = 42)
  covs <- c("age", "sex", "biomarker")
  set.seed(90125)
  dat_perm <- dat
  dat_perm$event_24 <- sample(dat_perm$event_24)
  lock <- create_analysis_lock(dat, "treatment", "event_24", covs,
                               seed = 7, negative_controls = "nc_outcome")
  lock_perm <- create_analysis_lock(dat_perm, "treatment", "event_24", covs,
                                    seed = 7, negative_controls = "nc_outcome")
  lock_masked <- mask_outcome(lock)
  list(dat = dat, lock = lock, perm = lock_perm, masked = lock_masked)
}

fx <- blinding_fixtures()

test_that("mask_outcome blanks the column on a copy and flags the lock", {
  expect_true("event_24" %in% names(fx$masked$data))
  expect_true(all(is.na(fx$masked$data$event_24)))
  expect_true(isTRUE(fx$masked$.outcome_masked))
  expect_false(anyNA(fx$lock$data$event_24))
})

test_that("fit_ps(method = 'glm') never reads the outcome", {
  ps  <- fit_ps(fx$lock,   method = "glm")
  psM <- fit_ps(fx$masked, method = "glm")
  psP <- fit_ps(fx$perm,   method = "glm")
  expect_blind_identical(ps, psM)
  expect_blind_identical(ps, psP)
})

test_that("assess_support never reads the outcome", {
  sup  <- assess_support(fit_ps(fx$lock,   method = "glm"))
  supM <- assess_support(fit_ps(fx$masked, method = "glm"))
  supP <- assess_support(fit_ps(fx$perm,   method = "glm"))
  expect_blind_identical(sup, supM)
  expect_blind_identical(sup, supP)
})

test_that("estimand_feasibility (with profile) never reads the outcome", {
  fea  <- estimand_feasibility(fit_ps(fx$lock,   method = "glm"),
                               profile_vars = c("age", "sex"))
  feaM <- estimand_feasibility(fit_ps(fx$masked, method = "glm"),
                               profile_vars = c("age", "sex"))
  feaP <- estimand_feasibility(fit_ps(fx$perm,   method = "glm"),
                               profile_vars = c("age", "sex"))
  expect_blind_identical(fea, feaM)
  expect_blind_identical(fea, feaP)
})

test_that("run_negative_control_ladder reads NC columns, never the outcome", {
  rungs <- list("older half" = fx$dat$age >= stats::median(fx$dat$age))
  ncl  <- run_negative_control_ladder(fx$lock,   rungs,
                                      method = "unadjusted", verbose = FALSE)
  nclM <- run_negative_control_ladder(fx$masked, rungs,
                                      method = "unadjusted", verbose = FALSE)
  nclP <- run_negative_control_ladder(fx$perm,   rungs,
                                      method = "unadjusted", verbose = FALSE)
  expect_blind_identical(ncl, nclM)
  expect_blind_identical(ncl, nclP)
})

test_that("check_process_indicators never reads the outcome", {
  cpi  <- check_process_indicators(fx$lock,   indicators = "censored")
  cpiM <- check_process_indicators(fx$masked, indicators = "censored")
  cpiP <- check_process_indicators(fx$perm,   indicators = "censored")
  expect_blind_identical(cpi, cpiM)
  expect_blind_identical(cpi, cpiP)
})

test_that("declare_estimand_ladder writes the same declaration blind", {
  dl <- function(l) {
    l2 <- declare_estimand_ladder(l, primary = "ATE",
                                  fallbacks = c("trimmed_ATE", "ATT", "ATO"),
                                  trigger = "SEVERE")
    list(ladder = l2$estimand_ladder, log = l2$design_log)
  }
  expect_blind_identical(dl(fx$lock), dl(fx$masked))
  expect_blind_identical(dl(fx$lock), dl(fx$perm))
})

test_that("simulate_support with an explicit base rate never reads the outcome", {
  surf <- list(confounding = c(0, 1), modification = 0, base_rate = 0.10)
  ss  <- simulate_support(fx$lock,   surface = surf, reps = 20L,
                          verbose = FALSE)
  ssM <- simulate_support(fx$masked, surface = surf, reps = 20L,
                          verbose = FALSE)
  ssP <- simulate_support(fx$perm,   surface = surf, reps = 20L,
                          verbose = FALSE)
  expect_blind_identical(ss, ssM)
  expect_blind_identical(ss, ssP)
})

test_that("simulate_support default base rate uses only the marginal outcome", {
  # A permutation preserves the marginal distribution, so the anchored
  # base rate (and everything downstream of it) must not change.
  surf <- list(confounding = c(0, 1), modification = 0)
  ss  <- simulate_support(fx$lock, surface = surf, reps = 10L,
                          verbose = FALSE)
  ssP <- simulate_support(fx$perm, surface = surf, reps = 10L,
                          verbose = FALSE)
  expect_blind_identical(ss, ssP)
})

test_that("design_report is blind and carries no raw data", {
  rep_of <- function(l) {
    l2 <- declare_estimand_ladder(l, primary = "ATE",
                                  fallbacks = c("trimmed_ATE", "ATT", "ATO"),
                                  trigger = "SEVERE")
    ps <- fit_ps(l2, method = "glm")
    design_report(l2, assess_support(ps), estimand_feasibility(ps))
  }
  dr  <- rep_of(fx$lock)
  drM <- rep_of(fx$masked)
  drP <- rep_of(fx$perm)
  expect_blind_identical(dr, drM)
  expect_blind_identical(dr, drP)
  # The report must not embed the analytic data frame anywhere.
  has_df_of_n <- function(x, n) {
    if (is.data.frame(x)) return(nrow(x) == n && "event_24" %in% names(x))
    if (is.list(x)) return(any(vapply(x, has_df_of_n, logical(1), n = n)))
    FALSE
  }
  expect_false(has_df_of_n(dr, nrow(fx$dat)))
})

test_that("estimate_design_precision reads only marginal outcome summaries", {
  skip(paste("WP1 item 1 pending (Stage 0 finding F3):",
             "estimate_design_precision() and summarize_event_support()",
             "return per-arm event counts and crude per-arm rates, the",
             "crude treatment-outcome association. Once WP1 restricts the",
             "default output to totals and the marginal rate, remove this",
             "skip: the assertions below must then pass."))
  dp  <- estimate_design_precision(fx$lock)
  dpP <- estimate_design_precision(fx$perm)
  # A permutation preserves total events and the marginal rate, so a
  # marginal-only summary is permutation invariant.
  expect_blind_identical(dp, dpP)
})

test_that("plasmode results do not carry the primary outcome", {
  skip(paste("WP1 item 2 pending (Stage 0 finding F4):",
             "run_plasmode_feasibility() and run_plasmode_dq_stress()",
             "return the entire lock, including the real outcome column,",
             "as result$lock. Once WP1 stops storing the outcome on the",
             "returned objects, remove this skip."))
  cand <- list(tmle_candidate("c1", g_library = "SL.glm", truncation = 0.01))
  plas <- run_plasmode_feasibility(fx$lock, tmle_candidates = cand,
                                   effect_sizes = 0.05, reps = 2L)
  carries_y <- function(x) {
    if (is.data.frame(x)) return("event_24" %in% names(x) &&
                                   !all(is.na(x$event_24)))
    if (is.list(x)) return(any(vapply(x, carries_y, logical(1))))
    FALSE
  }
  expect_false(carries_y(plas))
})

test_that("plasmode external-pilot mode is fully outcome blind", {
  skip(paste("WP1 item 2 / decision D8 pending: dgp_mode",
             "('hybrid', 'external_pilot') is not implemented yet.",
             "Once implemented, external_pilot runs must be identical on",
             "the reference, masked, and permuted locks."))
})
