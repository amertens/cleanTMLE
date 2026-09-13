# simulate_support(): the outcome-blind support simulation, and the
# generate-treatment fix to the plasmode design (Shaw et al. 2025).

# A small cohort with real confounding and moderate overlap. The treatment
# depends on two covariates, so the propensity direction is meaningful.
.make_support_data <- function(n = 300, seed = 11) {
  set.seed(seed)
  x1 <- stats::rnorm(n)
  x2 <- stats::rbinom(n, 1, 0.4)
  g  <- stats::plogis(-0.3 + 1.2 * x1 + 0.8 * x2)
  A  <- stats::rbinom(n, 1, g)
  y  <- stats::rbinom(n, 1, stats::plogis(-2 + 0.5 * x1 + 0.4 * x2 + 0.3 * A))
  data.frame(x1 = x1, x2 = x2, treatment = A, outcome = y)
}

test_that("support_surfaces validates and prints", {
  s <- support_surfaces()
  expect_s3_class(s, "support_surfaces")
  expect_error(support_surfaces(base_rate = 2), "base_rate")
  expect_output(print(s), "confounding strengths")
})

test_that("simulate_support returns per-estimand metrics with known truths", {
  dat  <- .make_support_data(300)
  lock <- create_simple_lock(dat, "treatment", "outcome", c("x1", "x2"),
                             seed = 7L)
  sim <- simulate_support(
    lock,
    surface = support_surfaces(confounding = c(0, 2), modification = c(0, 1),
                               effect = 0.10, base_rate = 0.15),
    reps = 40, verbose = FALSE)

  expect_s3_class(sim, "support_simulation")
  expect_true(isTRUE(sim$demo))          # < 50 reps is demonstration only
  expect_identical(sim$design, "generate_treatment")
  m <- sim$metrics
  expect_true(all(c("ATE", "ATT", "ATO", "trimmed_ATE", "matched_ATT") %in%
                    m$estimand))
  expect_true(all(m$n_converged > 0))

  # At confounding 0 and modification 0 every estimand has the same truth
  # (a constant additive effect), and the ATE estimator recovers it.
  base <- m[m$confounding == 0 & m$modification == 0, ]
  expect_lt(diff(range(base$truth_mean[base$estimand %in%
                                         c("ATE", "ATT", "ATO")])), 0.005)
  expect_lt(abs(base$bias[base$estimand == "ATE"]), 0.03)

  # With effect modification along the propensity direction the ATT and ATE
  # truths genuinely differ: treated units sit higher on the direction.
  modded <- m[m$confounding == 0 & m$modification == 1, ]
  t_ate <- modded$truth_mean[modded$estimand == "ATE"]
  t_att <- modded$truth_mean[modded$estimand == "ATT"]
  expect_gt(t_att - t_ate, 0.02)

  # Each estimator tracks its own estimand even where the truths differ.
  expect_lt(abs(modded$bias[modded$estimand == "ATT"]), 0.035)
  expect_lt(abs(modded$bias[modded$estimand == "ATO"]), 0.035)

  # Printing summarises feasibility per estimand.
  expect_output(print(sim), "FLAG")
})

test_that("the sample-treatment design warns and is recorded", {
  dat  <- .make_support_data(200)
  lock <- create_simple_lock(dat, "treatment", "outcome", c("x1", "x2"),
                             seed = 3L)
  expect_warning(
    sim <- simulate_support(
      lock, surface = support_surfaces(confounding = 0, modification = 0),
      reps = 5, design = "sample_treatment", verbose = FALSE),
    "Shaw")
  expect_identical(sim$design, "sample_treatment")
})

test_that("q0 sources are enforced: primary outcome needs explicit authorisation", {
  dat  <- .make_support_data(200)
  lock <- create_simple_lock(dat, "treatment", "outcome", c("x1", "x2"),
                             seed = 3L)
  expect_error(
    simulate_support(lock, reps = 5, q0_source = "primary_outcome",
                     verbose = FALSE),
    "not outcome-blind")
  sim <- simulate_support(
    lock, surface = support_surfaces(confounding = 0, modification = 0),
    reps = 5, q0_source = "heldout", verbose = FALSE)
  expect_gt(length(sim$q0$heldout_rows), 0)
})

test_that("the generate-treatment design beats the Shaw et al. artifact", {
  # Under the sample-treatment design the synthetic assignment mechanism is
  # degenerate (P(A = a | W) = 1 at the observed a), so the ATE estimator
  # shows bias and undercoverage that the generate-treatment design does not,
  # on the same cohort, at the same replicate count. Deterministic given the
  # lock seed.
  dat  <- .make_support_data(250, seed = 21)
  lock <- create_simple_lock(dat, "treatment", "outcome", c("x1", "x2"),
                             seed = 13L)
  surf <- support_surfaces(confounding = 2, modification = 0,
                           effect = 0.10, base_rate = 0.15)
  gen <- simulate_support(lock, surface = surf, reps = 80, verbose = FALSE)
  smp <- suppressWarnings(simulate_support(lock, surface = surf, reps = 80,
                                           design = "sample_treatment",
                                           verbose = FALSE))
  g_ate <- gen$metrics[gen$metrics$estimand == "ATE", ]
  s_ate <- smp$metrics[smp$metrics$estimand == "ATE", ]
  # The artifact shows as inflated bias: with treatment fixed, the
  # finite-sample error of the fitted propensity never averages out across
  # replicates. Coverage is distorted too, but its direction depends on how
  # the missing A-variation deflates the empirical SD, so the deterministic
  # claims are the bias ordering and the generate design's calibration.
  expect_lt(abs(g_ate$bias), abs(s_ate$bias))
  expect_gte(g_ate$coverage, 0.85)
})

test_that("run_plasmode_feasibility exposes and records the design", {
  dat  <- .make_support_data(200)
  lock <- create_simple_lock(dat, "treatment", "outcome", c("x1", "x2"),
                             plasmode_reps = 5L, seed = 5L)
  cands <- list(tmle_candidate("glm_t01", g_library = "SL.glm",
                               q_library = "SL.glm", truncation = 0.01))
  plas <- run_plasmode_feasibility(lock, tmle_candidates = cands,
                                   effect_sizes = 0.05, reps = 5L)
  expect_identical(plas$design, "generate_treatment")
  expect_true(all(is.finite(plas$metrics$bias)))
  expect_warning(
    run_plasmode_feasibility(lock, tmle_candidates = cands,
                             effect_sizes = 0.05, reps = 3L,
                             design = "sample_treatment"),
    "Shaw")
})

test_that("run_plasmode_dq_stress exposes and records the design", {
  dat  <- .make_support_data(200)
  lock <- create_simple_lock(dat, "treatment", "outcome", c("x1", "x2"),
                             seed = 5L)
  cands <- list(tmle_candidate("glm_t01", g_library = "SL.glm",
                               q_library = "SL.glm", truncation = 0.01))
  dq <- run_plasmode_dq_stress(
    lock, tmle_candidates = cands, effect_sizes = 0.05, reps = 3L,
    data_quality_scenarios = list(
      covariate_missingness = list(fractions = 0.1)),
    verbose = FALSE)
  expect_identical(dq$design, "generate_treatment")
  expect_true("cov_miss" %in% dq$metrics$scenario)
})

test_that("check_locked_estimator runs and labels itself optimistic", {
  dat  <- .make_support_data(250)
  lock <- create_simple_lock(dat, "treatment", "outcome", c("x1", "x2"),
                             seed = 9L)
  chk <- check_locked_estimator(lock, reps = 30)
  expect_true(isTRUE(chk$optimistic))
  expect_lt(abs(chk$bias), 0.05)
  expect_gt(chk$coverage, 0.8)
})
