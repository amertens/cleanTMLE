# Outcome-blind support simulation.
#
# simulate_support() is the design team's answer to "which estimands can this
# design deliver?". It generates synthetic data under the generate-treatment
# plasmode design (resample W, draw A from the fitted propensity, draw Y from
# a prespecified synthetic outcome family), computes the truth for every
# estimand on every replicate from the generating model, and reports bias,
# RMSE, coverage and SE calibration per estimator per estimand per point on a
# confounding-by-effect-modification grid. The support map summarises where
# each estimand's estimator stays inside a prespecified bias tolerance with
# near-nominal coverage. An estimand whose estimator fails inside the plausible
# region of the grid is infeasible for this design.

#' Prespecify the Synthetic Outcome-Surface Family for simulate_support()
#'
#' The synthetic outcome family spans three axes. The confounding axis scales
#' how strongly the outcome depends on the covariate direction that most
#' predicts treatment (the standardised logit of the fitted propensity), so
#' extrapolation error under non-overlap grows along it. The modification axis
#' makes the treatment effect vary along the same direction, so the ATE, ATT
#' and ATO genuinely differ and each estimator is judged against its own
#' estimand. The complexity axis optionally adds curvature (a squared term in
#' the propensity direction) that a main-terms outcome model cannot represent,
#' so Q misspecification is part of the stress. Positivity bias is
#' extrapolation error: it depends on how the outcome surface behaves where
#' only one arm is observed, which is exactly what this family varies.
#'
#' @param confounding Numeric vector of confounding strengths (log-odds of the
#'   baseline outcome per SD of the propensity direction). Default
#'   `c(0, 1, 2)`.
#' @param modification Numeric vector of effect-modification strengths. At
#'   `m`, the additive effect at propensity-direction value `s` is
#'   `effect * (1 + m * s)`, so `m = 0` is a constant effect and at `m = 1`
#'   the effect doubles one SD into the treated-typical region and vanishes
#'   one SD into the control-typical region. Default `c(0, 1)`.
#' @param complexity Character vector, subset of `c("linear", "interaction")`.
#'   `"interaction"` adds `0.5 * confounding * (s^2 - 1)` to the baseline
#'   log-odds, curvature a main-terms Q model misses. Default `"linear"`.
#' @param effect Numeric; the additive risk difference at `s = 0`.
#'   Default 0.05.
#' @param base_rate Numeric in (0, 1) or `NULL`; the marginal baseline
#'   outcome rate the synthetic surface is anchored to. `NULL` uses the
#'   observed marginal outcome rate when the lock's outcome is readable
#'   (a marginal summary, not the treatment-outcome association), else 0.10.
#'
#' @return An object of class `support_surfaces`: the grid specification.
#' @seealso [simulate_support()]
#' @examples
#' support_surfaces()
#' support_surfaces(confounding = c(0, 0.5, 1, 2), modification = c(0, 0.5, 1))
#' @keywords internal
support_surfaces <- function(confounding = c(0, 1, 2),
                             modification = c(0, 1),
                             complexity = "linear",
                             effect = 0.05,
                             base_rate = NULL) {
  complexity <- match.arg(complexity, c("linear", "interaction"),
                          several.ok = TRUE)
  stopifnot(is.numeric(confounding), all(confounding >= 0),
            is.numeric(modification), all(modification >= 0),
            is.numeric(effect), length(effect) == 1L)
  if (!is.null(base_rate) &&
      (!is.numeric(base_rate) || base_rate <= 0 || base_rate >= 1))
    stop("`base_rate` must be in (0, 1) or NULL.", call. = FALSE)
  out <- list(confounding = sort(unique(confounding)),
              modification = sort(unique(modification)),
              complexity = complexity,
              effect = effect,
              base_rate = base_rate)
  class(out) <- "support_surfaces"
  out
}

#' @export
print.support_surfaces <- function(x, ...) {
  cat("Synthetic outcome-surface family for simulate_support()\n")
  cat(sprintf("  confounding strengths:  %s\n",
              paste(x$confounding, collapse = ", ")))
  cat(sprintf("  effect modification:    %s\n",
              paste(x$modification, collapse = ", ")))
  cat(sprintf("  complexity:             %s\n",
              paste(x$complexity, collapse = ", ")))
  cat(sprintf("  additive effect at s=0: %g\n", x$effect))
  cat(sprintf("  baseline rate:          %s\n",
              if (is.null(x$base_rate)) "observed marginal (or 0.10)"
              else format(x$base_rate)))
  invisible(x)
}


# ── Internal: lean per-replicate estimators ─────────────────────────────────
# These are deliberately cheap (GLM nuisances) because they run
# reps x grid times. The question the simulation answers is about the design
# and the estimand, not about the nuisance library; the locked candidate's
# library is exercised by run_plasmode_feasibility() and the DQ stress test.

.ss_fit_g <- function(A, W) {
  df <- data.frame(A = A, W)
  fit <- stats::glm(A ~ ., data = df, family = stats::binomial())
  as.numeric(stats::predict(fit, type = "response"))
}

.ss_fit_q <- function(Y, A, W) {
  df <- data.frame(Y = Y, A = A, W)
  fit <- stats::glm(Y ~ ., data = df, family = stats::binomial())
  d1 <- df; d1$A <- 1L
  d0 <- df; d0$A <- 0L
  list(Q1 = as.numeric(stats::predict(fit, newdata = d1, type = "response")),
       Q0 = as.numeric(stats::predict(fit, newdata = d0, type = "response")),
       Qa = as.numeric(stats::predict(fit, type = "response")))
}

.ss_est_tmle_ate <- function(Y, A, W, truncation = 0.01) {
  n <- length(Y)
  g <- pmin(pmax(.ss_fit_g(A, W), truncation), 1 - truncation)
  q <- .ss_fit_q(Y, A, W)
  H1 <- 1 / g; H0 <- -1 / (1 - g)
  Ha <- ifelse(A == 1, H1, H0)
  lq <- stats::qlogis(pmax(pmin(q$Qa, 0.999), 0.001))
  eps <- tryCatch(unname(stats::coef(stats::glm(
    Y ~ -1 + Ha + offset(lq), family = stats::binomial()))),
    error = function(e) 0)
  Q1u <- stats::plogis(stats::qlogis(pmax(pmin(q$Q1, .999), .001)) + eps * H1)
  Q0u <- stats::plogis(stats::qlogis(pmax(pmin(q$Q0, .999), .001)) + eps * H0)
  Qau <- stats::plogis(lq + eps * Ha)
  psi <- mean(Q1u) - mean(Q0u)
  eic <- Ha * (Y - Qau) + (Q1u - Q0u) - psi
  se  <- sqrt(stats::var(eic) / n)
  list(est = psi, se = se)
}

.ss_est_att_dr <- function(Y, A, W, truncation = 0.01) {
  n <- length(Y)
  g <- pmin(pmax(.ss_fit_g(A, W), truncation), 1 - truncation)
  q <- .ss_fit_q(Y, A, W)
  p1 <- mean(A)
  if (p1 <= 0) return(list(est = NA_real_, se = NA_real_))
  num <- A * (Y - q$Q0) - (1 - A) * (g / (1 - g)) * (Y - q$Q0)
  psi <- mean(num) / p1
  phi <- (num - A * psi) / p1
  list(est = psi, se = sqrt(stats::var(phi) / n))
}

.ss_est_ato_aug <- function(Y, A, W, truncation = 1e-6) {
  n <- length(Y)
  g <- pmin(pmax(.ss_fit_g(A, W), truncation), 1 - truncation)
  q <- .ss_fit_q(Y, A, W)
  h <- g * (1 - g)
  num <- h * (q$Q1 - q$Q0) +
    A * (1 - g) * (Y - q$Q1) - (1 - A) * g * (Y - q$Q0)
  hbar <- mean(h)
  psi  <- mean(num) / hbar
  phi  <- (num - h * psi) / hbar
  list(est = psi, se = sqrt(stats::var(phi) / n))
}

.ss_est_trimmed_tmle <- function(Y, A, W, band = c(0.05, 0.95),
                                 truncation = 0.01) {
  g0 <- .ss_fit_g(A, W)
  keep <- g0 >= band[1] & g0 <= band[2]
  if (sum(keep) < 50 || sum(A[keep] == 1) < 15 || sum(A[keep] == 0) < 15)
    return(list(est = NA_real_, se = NA_real_, keep = keep))
  # Refit g on the trimmed subset: reusing the full-cohort g would describe
  # the overlap of a population that no longer exists.
  out <- .ss_est_tmle_ate(Y[keep], A[keep], W[keep, , drop = FALSE],
                          truncation = truncation)
  out$keep <- keep
  out
}

.ss_est_matched_att <- function(Y, A, W, caliper_sd = 0.2,
                                truncation = 0.01) {
  g <- pmin(pmax(.ss_fit_g(A, W), truncation), 1 - truncation)
  m <- .greedy_caliper_match(g, A, caliper_sd = caliper_sd)
  if (length(m$treated) < 10)
    return(list(est = NA_real_, se = NA_real_, matched_treated = m$treated))
  diffs <- Y[m$treated] - Y[m$control]
  list(est = mean(diffs), se = sqrt(stats::var(diffs) / length(diffs)),
       matched_treated = m$treated)
}


#' Outcome-Blind Support Simulation and Support Map
#'
#' Decides, before outcome access, which estimands this design can deliver.
#' Synthetic replicates are generated under the generate-treatment plasmode
#' design: covariate rows are resampled with replacement, treatment is drawn
#' from the propensity model fitted on the real data, and the outcome is drawn
#' from the prespecified synthetic family of [support_surfaces()]. On every
#' replicate the truth for each estimand (ATE, ATT among the treated as
#' realised, ATO, trimmed ATE on the prespecified band, matched ATT on the
#' matched set) is computed from the generating model, so the output is bias,
#' RMSE, coverage and SE calibration per estimand per grid point. Because the
#' matched, trimmed, ATT and ATO estimators run on the same replicates as the
#' full-cohort ATE, differences between them decompose into estimand
#' difference (their truths differ), efficiency loss (their empirical SDs
#' differ), and extrapolation error (only the ATE's bias grows along the
#' confounding axis when overlap is poor).
#'
#' @section What the map means:
#' The gate reads the map. A cell passes when absolute bias is at most
#' `bias_tolerance` and coverage is at least `coverage_floor`. An estimand
#' whose estimator fails inside the plausible region of the grid is
#' infeasible for this design; which region is plausible is a subject-matter
#' judgement the design team records. The simulation never uses the real
#' treatment-outcome association unless `q0_source = "primary_outcome"` is
#' explicitly authorised.
#'
#' @param lock A `cleanroom_lock`.
#' @param ps_fit Optional `ps_fit`; when supplied its propensity scores are
#'   the generating propensity, otherwise a logistic model is fitted on the
#'   lock's covariates.
#' @param surface A [support_surfaces()] specification.
#' @param reps Integer; replicates per grid point. Default 200. Runs below
#'   50 are tagged demonstration-only and the printed verdict is FLAG.
#' @param design `"generate_treatment"` (default) or the deprecated
#'   `"sample_treatment"` (Shaw et al. 2025 artifact; warns).
#' @param band Numeric length 2; the trimming band for the trimmed-ATE
#'   estimand and the region used by the support verdicts. Default
#'   `c(0.05, 0.95)`.
#' @param q0_source Where the baseline outcome surface comes from:
#'   `"synthetic"` (default; the [support_surfaces()] family anchored at
#'   `base_rate`), `"negative_control"` or `"auxiliary"` (a covariate-only
#'   GLM fitted to `q0_variable`, with the confounding axis adding outcome
#'   dependence on the propensity direction on top of it), `"heldout"` (fit
#'   on a design split whose rows are recorded in `heldout_rows` and must be
#'   excluded from Stage 4), or `"primary_outcome"` (requires
#'   `allow_outcome_q0 = TRUE`; not outcome-blind).
#' @param q0_variable Column name for `"negative_control"` / `"auxiliary"`.
#' @param heldout_fraction Fraction for `"heldout"`. Default 0.3.
#' @param allow_outcome_q0 Logical; must be `TRUE` to fit Q0 on the primary
#'   outcome. Default `FALSE`.
#' @param bias_tolerance Absolute bias a cell may show and still pass.
#'   Default `0.5 * surface$effect`.
#' @param coverage_floor Minimum coverage for a cell to pass. Default 0.90.
#' @param caliper_sd Caliper for the matched-ATT estimator, in SDs of the
#'   logit propensity. Default 0.2.
#' @param truncation Truncation for the lean estimators' fitted g.
#'   Default 0.01.
#' @param checkpoint_file Optional path; results are saved after every grid
#'   point and a partial run resumes from it.
#' @param verbose Print progress. Default TRUE.
#'
#' @return An object of class `support_simulation`: `metrics` (one row per
#'   grid point per estimand), `map` (the pass/fail matrix), `feasible`
#'   (per-estimand verdict over the whole grid), `surface`, `design`, `q0`
#'   (source, variable, `heldout_rows`), `band`, `reps`, `demo` flag.
#'
#' @references Shaw PA et al. (2025) arXiv:2504.11740. Petersen ML et al.
#'   (2012) Stat Methods Med Res 21:31-54. Li F, Morgan KL, Zaslavsky AM
#'   (2018) JASA 113:390-400. Mao H, Li L, Greene T (2019) SMMR 28:2439-2454.
#'
#' @examples
#' \dontrun{
#' lock <- create_simple_lock(sim_func1(400), "treatment", "event_24",
#'                            c("age", "sex", "biomarker"))
#' sim <- simulate_support(lock, reps = 100)
#' print(sim)
#' plot(sim)
#' }
#' @export
simulate_support <- function(lock,
                             ps_fit = NULL,
                             surface = NULL,
                             reps = 200L,
                             design = c("generate_treatment",
                                        "sample_treatment",
                                        "parametric_bootstrap"),
                             band = c(0.05, 0.95),
                             q0_source = c("synthetic", "negative_control",
                                           "auxiliary", "heldout",
                                           "primary_outcome"),
                             q0_variable = NULL,
                             heldout_fraction = 0.3,
                             allow_outcome_q0 = FALSE,
                             bias_tolerance = NULL,
                             coverage_floor = 0.90,
                             caliper_sd = 0.2,
                             truncation = 0.01,
                             checkpoint_file = NULL,
                             verbose = TRUE) {
  if (!inherits(lock, "cleanroom_lock"))
    stop("`lock` must be a cleanroom_lock object.", call. = FALSE)
  # `surface` is a plain named list of overrides for the prespecified
  # outcome-surface family (confounding, modification, complexity, effect,
  # base_rate); NULL takes the documented defaults.
  if (is.null(surface)) {
    surface <- support_surfaces()
  } else if (!inherits(surface, "support_surfaces")) {
    if (!is.list(surface) || is.null(names(surface)))
      stop("`surface` must be a named list (confounding, modification, ",
           "complexity, effect, base_rate).", call. = FALSE)
    bad <- setdiff(names(surface), names(formals(support_surfaces)))
    if (length(bad))
      stop("Unknown surface field(s): ", paste(bad, collapse = ", "),
           call. = FALSE)
    surface <- do.call(support_surfaces, surface)
  }
  design    <- match.arg(design)
  q0_source <- match.arg(q0_source)
  if (design == "sample_treatment")
    rlang::warn(paste(
      "design = 'sample_treatment' keeps the observed treatment; Shaw et al.",
      "(2025, arXiv:2504.11740) show this induces a positivity violation by",
      "construction. Use the default generate_treatment design."),
      .frequency = "once", .frequency_id = "simulate_support_sample_trt")
  if (design == "parametric_bootstrap") {
    # Petersen et al.'s (2012) check on the locked estimator: the fitted
    # nuisance models generate the truth, so the result is optimistic by
    # construction and labelled as such.
    return(check_locked_estimator(lock, ps_fit = ps_fit, reps = reps,
                                  truncation = truncation))
  }
  if (q0_source == "primary_outcome" && !isTRUE(allow_outcome_q0))
    stop("Fitting Q0 on the primary outcome is not outcome-blind. ",
         "Pass allow_outcome_q0 = TRUE to authorise it explicitly, or use ",
         "q0_source = 'synthetic', 'negative_control', 'auxiliary' or ",
         "'heldout'.", call. = FALSE)

  data       <- lock$data
  covariates <- lock$covariates
  treatment  <- lock$treatment
  A_obs      <- as.integer(data[[treatment]])
  n          <- nrow(data)
  W          <- data[, covariates, drop = FALSE]
  W <- as.data.frame(lapply(W, function(x) {
    x <- if (is.numeric(x)) x else as.numeric(as.factor(x))
    if (anyNA(x)) x[is.na(x)] <- stats::median(x, na.rm = TRUE)
    x
  }))

  # Generating propensity: the supplied fit or a logistic model.
  if (!is.null(ps_fit)) {
    if (!inherits(ps_fit, "ps_fit"))
      stop("`ps_fit` must be a ps_fit object.", call. = FALSE)
    g_gen <- as.numeric(ps_fit$ps)
    if (length(g_gen) != n)
      stop("ps_fit does not match the lock's rows.", call. = FALSE)
  } else {
    g_gen <- .ss_fit_g(A_obs, W)
  }
  g_gen <- pmin(pmax(g_gen, 0.005), 0.995)

  # The propensity direction: standardised logit of the generating g.
  s_all <- stats::qlogis(g_gen)
  s_all <- (s_all - mean(s_all)) / max(stats::sd(s_all), 1e-8)

  # Baseline outcome surface.
  heldout_rows <- integer(0)
  q0_lp_extra  <- rep(0, n)   # log-odds contribution of a fitted Q0 source
  base_rate <- surface$base_rate
  if (q0_source %in% c("negative_control", "auxiliary")) {
    if (is.null(q0_variable) || !q0_variable %in% names(data))
      stop("q0_source = '", q0_source, "' needs `q0_variable`, a column of ",
           "the lock data.", call. = FALSE)
    yv <- as.numeric(data[[q0_variable]])
    fit <- stats::glm(yv ~ ., data = cbind(data.frame(yv = yv), W),
                      family = stats::binomial(), na.action = stats::na.exclude)
    pr  <- as.numeric(stats::predict(fit, newdata = cbind(data.frame(yv = 0), W),
                                     type = "response"))
    pr[!is.finite(pr)] <- mean(pr[is.finite(pr)])
    q0_lp_extra <- stats::qlogis(pmin(pmax(pr, 0.005), 0.995))
    q0_lp_extra <- q0_lp_extra - mean(q0_lp_extra)
    if (is.null(base_rate)) base_rate <- mean(yv, na.rm = TRUE)
  } else if (q0_source %in% c("heldout", "primary_outcome")) {
    yv <- as.numeric(data[[lock$outcome]])
    if (all(is.na(yv)))
      stop("The lock outcome is fully NA (masked); unmask it or choose a ",
           "different q0_source.", call. = FALSE)
    rows <- which(!is.na(yv))
    if (q0_source == "heldout") {
      set.seed(lock$seed + 1717L)
      heldout_rows <- sort(sample(rows, ceiling(length(rows) *
                                                  heldout_fraction)))
      rows <- heldout_rows
    }
    fit <- stats::glm(yv ~ ., data = cbind(data.frame(yv = yv), W),
                      subset = rows, family = stats::binomial())
    pr  <- as.numeric(stats::predict(fit, newdata = cbind(data.frame(yv = 0), W),
                                     type = "response"))
    q0_lp_extra <- stats::qlogis(pmin(pmax(pr, 0.005), 0.995))
    q0_lp_extra <- q0_lp_extra - mean(q0_lp_extra)
    if (is.null(base_rate)) base_rate <- mean(yv[rows], na.rm = TRUE)
  } else if (is.null(base_rate)) {
    yv <- suppressWarnings(as.numeric(data[[lock$outcome]]))
    base_rate <- if (!all(is.na(yv))) mean(yv, na.rm = TRUE) else 0.10
    if (!is.finite(base_rate) || base_rate <= 0 || base_rate >= 1)
      base_rate <- 0.10
  }

  grid <- expand.grid(confounding = surface$confounding,
                      modification = surface$modification,
                      complexity = surface$complexity,
                      stringsAsFactors = FALSE)

  estimands <- c("ATE", "ATT", "ATO", "trimmed_ATE", "matched_ATT")
  eff <- surface$effect

  # Resume support: reload finished grid cells.
  done <- list()
  if (!is.null(checkpoint_file) && file.exists(checkpoint_file)) {
    done <- readRDS(checkpoint_file)
    if (verbose) message("simulate_support: resuming; ", length(done),
                         " grid cells found in ", checkpoint_file)
  }

  metric_rows <- list()
  for (gi in seq_len(nrow(grid))) {
    cf <- grid$confounding[gi]; md <- grid$modification[gi]
    cx <- grid$complexity[gi]
    cell_key <- sprintf("c%.3f_m%.3f_%s", cf, md, cx)
    if (!is.null(done[[cell_key]])) {
      metric_rows[[cell_key]] <- done[[cell_key]]
      next
    }
    if (verbose)
      message(sprintf("simulate_support: confounding=%.2f modification=%.2f %s",
                      cf, md, cx))

    # Per-replicate estimates and truths.
    est_arr <- array(NA_real_, dim = c(reps, length(estimands), 2),
                     dimnames = list(NULL, estimands, c("est", "se")))
    tru_arr <- matrix(NA_real_, nrow = reps, ncol = length(estimands),
                      dimnames = list(NULL, estimands))

    for (r in seq_len(reps)) {
      set.seed(lock$seed + 100000L * gi + r)
      if (design == "generate_treatment") {
        idx <- sample.int(n, n, replace = TRUE)
        A_r <- stats::rbinom(n, 1L, g_gen[idx])
      } else {
        idx <- seq_len(n)
        A_r <- A_obs
      }
      s   <- s_all[idx]
      g_t <- g_gen[idx]
      W_r <- W[idx, , drop = FALSE]

      # Baseline log-odds: anchor + fitted-source contribution + confounding
      # along the propensity direction (+ curvature under "interaction").
      lp0 <- q0_lp_extra[idx] + cf * s
      if (cx == "interaction") lp0 <- lp0 + 0.5 * cf * (s^2 - 1)
      # Solve the intercept so the marginal baseline rate is base_rate.
      a0 <- tryCatch(stats::uniroot(function(a)
        mean(stats::plogis(a + lp0)) - base_rate,
        lower = -12, upper = 12)$root,
        error = function(e) stats::qlogis(base_rate))
      p0 <- stats::plogis(a0 + lp0)
      tau <- eff * (1 + md * s)
      p1  <- pmin(pmax(p0 + tau, 0.001), 0.999)
      p0  <- pmin(pmax(p0, 0.001), 0.999)
      tau_real <- p1 - p0

      Y_r <- stats::rbinom(n, 1L, ifelse(A_r == 1, p1, p0))

      # Truths on this replicate, from the generating model.
      tru_arr[r, "ATE"] <- mean(tau_real)
      tru_arr[r, "ATT"] <- if (any(A_r == 1))
        mean(tau_real[A_r == 1]) else NA_real_
      h <- g_t * (1 - g_t)
      tru_arr[r, "ATO"] <- sum(h * tau_real) / sum(h)
      in_band <- g_t >= band[1] & g_t <= band[2]
      tru_arr[r, "trimmed_ATE"] <- if (any(in_band))
        mean(tau_real[in_band]) else NA_real_

      # Estimators.
      e_ate <- tryCatch(.ss_est_tmle_ate(Y_r, A_r, W_r, truncation),
                        error = function(e) list(est = NA_real_, se = NA_real_))
      e_att <- tryCatch(.ss_est_att_dr(Y_r, A_r, W_r, truncation),
                        error = function(e) list(est = NA_real_, se = NA_real_))
      e_ato <- tryCatch(.ss_est_ato_aug(Y_r, A_r, W_r),
                        error = function(e) list(est = NA_real_, se = NA_real_))
      e_trm <- tryCatch(.ss_est_trimmed_tmle(Y_r, A_r, W_r, band, truncation),
                        error = function(e) list(est = NA_real_, se = NA_real_))
      e_mat <- tryCatch(.ss_est_matched_att(Y_r, A_r, W_r, caliper_sd,
                                            truncation),
                        error = function(e)
                          list(est = NA_real_, se = NA_real_,
                               matched_treated = integer(0)))
      tru_arr[r, "matched_ATT"] <- if (length(e_mat$matched_treated) >= 10)
        mean(tau_real[e_mat$matched_treated]) else NA_real_

      est_arr[r, "ATE", ]         <- c(e_ate$est, e_ate$se)
      est_arr[r, "ATT", ]         <- c(e_att$est, e_att$se)
      est_arr[r, "ATO", ]         <- c(e_ato$est, e_ato$se)
      est_arr[r, "trimmed_ATE", ] <- c(e_trm$est, e_trm$se)
      est_arr[r, "matched_ATT", ] <- c(e_mat$est, e_mat$se)
    }

    cell_rows <- lapply(estimands, function(ed) {
      es <- est_arr[, ed, "est"]; ss <- est_arr[, ed, "se"]
      tr <- tru_arr[, ed]
      ok <- is.finite(es) & is.finite(tr)
      if (!any(ok)) return(NULL)
      es <- es[ok]; ss <- ss[ok]; tr <- tr[ok]
      lo <- es - 1.96 * ss; hi <- es + 1.96 * ss
      emp <- stats::sd(es - tr)
      data.frame(confounding = cf, modification = md, complexity = cx,
                 estimand = ed, estimator = .ss_estimator_label(ed),
                 truth_mean = round(mean(tr), 5),
                 bias = round(mean(es - tr), 5),
                 rmse = round(sqrt(mean((es - tr)^2)), 5),
                 coverage = round(mean(lo <= tr & tr <= hi), 3),
                 emp_sd = round(emp, 5),
                 mean_se = round(mean(ss, na.rm = TRUE), 5),
                 se_cal = round(if (is.finite(emp) && emp > 0)
                   mean(ss, na.rm = TRUE) / emp else NA_real_, 3),
                 n_converged = sum(ok),
                 stringsAsFactors = FALSE)
    })
    cell_df <- do.call(rbind, Filter(Negate(is.null), cell_rows))
    metric_rows[[cell_key]] <- cell_df
    if (!is.null(checkpoint_file)) {
      done[[cell_key]] <- cell_df
      saveRDS(done, checkpoint_file)
    }
  }

  metrics <- do.call(rbind, metric_rows)
  rownames(metrics) <- NULL
  if (is.null(bias_tolerance)) bias_tolerance <- 0.5 * abs(eff)
  metrics$pass <- abs(metrics$bias) <= bias_tolerance &
    metrics$coverage >= coverage_floor

  feasible <- vapply(estimands, function(ed) {
    m <- metrics[metrics$estimand == ed, ]
    nrow(m) > 0 && all(m$pass)
  }, logical(1))

  out <- list(
    metrics = metrics,
    map = metrics[, c("confounding", "modification", "complexity",
                      "estimand", "pass", "bias", "coverage")],
    feasible = feasible,
    bias_tolerance = bias_tolerance,
    coverage_floor = coverage_floor,
    surface = surface,
    design = design,
    band = band,
    reps = reps,
    demo = reps < 50L,
    q0 = list(source = q0_source, variable = q0_variable,
              base_rate = base_rate, heldout_rows = heldout_rows),
    lock_hash = lock$lock_hash,
    call = match.call()
  )
  class(out) <- "support_simulation"
  out
}

.ss_estimator_label <- function(ed) {
  switch(ed,
    ATE = "TMLE (GLM nuisances)",
    ATT = "DR one-step ATT",
    ATO = "augmented overlap weights",
    trimmed_ATE = "TMLE after trim + g refit",
    matched_ATT = "1:1 caliper match, paired difference",
    ed)
}

#' @export
print.support_simulation <- function(x, ...) {
  cat("Outcome-blind support simulation (", x$design, " design)\n", sep = "")
  cat(sprintf("  %d replicates per grid point; Q0 source: %s\n",
              x$reps, x$q0$source))
  if (isTRUE(x$demo))
    cat("  FLAG: fewer than 50 replicates; demonstration only, not",
        "an inferential run.\n")
  cat(sprintf("  pass rule: |bias| <= %.4f and coverage >= %.2f\n\n",
              x$bias_tolerance, x$coverage_floor))
  for (ed in unique(x$metrics$estimand)) {
    m <- x$metrics[x$metrics$estimand == ed, ]
    fails <- m[!m$pass, , drop = FALSE]
    if (nrow(fails) == 0) {
      cat(sprintf("  %-12s feasible across the grid (max |bias| %.4f, min coverage %.2f)\n",
                  ed, max(abs(m$bias)), min(m$coverage)))
    } else {
      worst <- fails[which.max(abs(fails$bias)), ]
      cat(sprintf("  %-12s FAILS at %d of %d grid points (worst: confounding %.1f, modification %.1f, bias %.4f, coverage %.2f)\n",
                  ed, nrow(fails), nrow(m), worst$confounding,
                  worst$modification, worst$bias, worst$coverage))
    }
  }
  invisible(x)
}

#' @export
plot.support_simulation <- function(x, metric = c("pass", "bias", "coverage"),
                                    ...) {
  metric <- match.arg(metric)
  if (!requireNamespace("ggplot2", quietly = TRUE))
    stop("ggplot2 is required.", call. = FALSE)
  m <- x$metrics
  m$conf_f <- factor(m$confounding)
  m$mod_f  <- factor(m$modification)
  fill_lab <- switch(metric, pass = "passes", bias = "bias",
                     coverage = "coverage")
  m$fill <- switch(metric, pass = m$pass, bias = m$bias,
                   coverage = m$coverage)
  p <- ggplot2::ggplot(m, ggplot2::aes(x = .data$conf_f, y = .data$mod_f,
                                       fill = .data$fill)) +
    ggplot2::geom_tile(colour = "white") +
    ggplot2::facet_wrap(~ estimand + complexity) +
    ggplot2::labs(x = "confounding strength (outcome dependence on the propensity direction)",
                  y = "effect-modification strength",
                  fill = fill_lab,
                  title = "Support map: where each estimand's estimator holds") +
    ggplot2::theme_minimal()
  if (metric == "pass")
    p <- p + ggplot2::scale_fill_manual(values = c(`TRUE` = "#2e7d32",
                                                   `FALSE` = "#c62828"))
  p
}


#' Parametric-Bootstrap Check of the Locked Estimator
#'
#' Petersen et al.'s (2012) parametric bootstrap for one fixed estimator: the
#' fitted nuisance models are treated as the truth, data are regenerated from
#' them, and the estimator is refit. Its known limitation is optimism, because
#' the fitted Q generates the truth the estimator is judged against; a clean
#' result here is necessary, not sufficient. Run it after candidate selection
#' as a final check, not instead of [simulate_support()].
#'
#' @param lock A `cleanroom_lock` (the outcome must be readable; this check
#'   is post-lock, pre-report).
#' @param ps_fit Optional `ps_fit` for the generating propensity.
#' @param reps Integer bootstrap replicates. Default 200.
#' @param truncation Truncation applied to the fitted g. Default 0.01.
#' @return A list with the per-replicate estimates, bias against the plug-in
#'   truth, coverage, and the `optimistic = TRUE` label.
#' @references Petersen ML, Porter KE, Gruber S, Wang Y, van der Laan MJ
#'   (2012). Diagnosing and responding to violations in the positivity
#'   assumption. Stat Methods Med Res 21:31-54.
#' @keywords internal
check_locked_estimator <- function(lock, ps_fit = NULL, reps = 200L,
                                   truncation = 0.01) {
  if (!inherits(lock, "cleanroom_lock"))
    stop("`lock` must be a cleanroom_lock object.", call. = FALSE)
  data <- lock$data
  Y <- as.numeric(data[[lock$outcome]])
  A <- as.integer(data[[lock$treatment]])
  if (all(is.na(Y)))
    stop("check_locked_estimator needs a readable outcome.", call. = FALSE)
  W <- data[, lock$covariates, drop = FALSE]
  W <- as.data.frame(lapply(W, function(x) {
    x <- if (is.numeric(x)) x else as.numeric(as.factor(x))
    if (anyNA(x)) x[is.na(x)] <- stats::median(x, na.rm = TRUE)
    x
  }))
  cc <- !is.na(Y)
  g_gen <- if (!is.null(ps_fit)) as.numeric(ps_fit$ps) else .ss_fit_g(A, W)
  g_gen <- pmin(pmax(g_gen, 0.005), 0.995)
  q <- .ss_fit_q(Y[cc], A[cc], W[cc, , drop = FALSE])
  # Plug-in truth from the fitted models (the source of the optimism).
  truth <- mean(q$Q1) - mean(q$Q0)
  n <- sum(cc)
  Wc <- W[cc, , drop = FALSE]; gc_ <- g_gen[cc]
  ests <- ses <- rep(NA_real_, reps)
  for (r in seq_len(reps)) {
    set.seed(lock$seed + 5000L + r)
    idx <- sample.int(n, n, replace = TRUE)
    A_r <- stats::rbinom(n, 1L, gc_[idx])
    p_r <- ifelse(A_r == 1, q$Q1[idx], q$Q0[idx])
    Y_r <- stats::rbinom(n, 1L, pmin(pmax(p_r, 0.001), 0.999))
    e <- tryCatch(.ss_est_tmle_ate(Y_r, A_r, Wc[idx, , drop = FALSE],
                                   truncation),
                  error = function(e) list(est = NA_real_, se = NA_real_))
    ests[r] <- e$est; ses[r] <- e$se
  }
  ok <- is.finite(ests)
  lo <- ests - 1.96 * ses; hi <- ests + 1.96 * ses
  list(truth_plugin = truth,
       bias = mean(ests[ok]) - truth,
       rmse = sqrt(mean((ests[ok] - truth)^2)),
       coverage = mean(lo[ok] <= truth & truth <= hi[ok]),
       n_converged = sum(ok),
       reps = reps,
       optimistic = TRUE,
       note = paste("Parametric bootstrap: the fitted models generate the",
                    "truth, so a clean result is necessary, not sufficient",
                    "(Petersen et al. 2012)."))
}
