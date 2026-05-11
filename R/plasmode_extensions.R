# Plasmode candidate-grid extensions (cleanTMLE 0.2 batch).
#
# Helpers and post-processors motivated by the March-workshop materials
# (Phillips et al. 2023; Gruber et al. 2023 SAP + supplement; Gruber et al.
# 2023 RWE; Gruber & van der Laan 2009). Implements the low-effort and
# medium-effort items in TODO.md A.21-A.25. The high-effort items
# (C-TMLE wrapper, nonparametric bootstrap variance, drtmle) are tracked
# but not implemented in this file.

# ------------------------------------------------------------------
# Effective sample size for super-learner tuning
# ------------------------------------------------------------------

#' Effective Sample Size for SuperLearner Tuning
#'
#' Implements `n_eff = min(n, 5 * n_rare)` for binary outcomes
#' (Phillips et al. 2023, *Stat Med*). For non-binary outcomes the
#' raw sample size is returned.
#'
#' The effective sample size drives the default V-fold CV choice and
#' the SuperLearner library composition used by [build_sl_library()].
#'
#' @param y Outcome vector. Binary outcomes can be 0/1, TRUE/FALSE,
#'   or factor.
#' @param family Either `"binomial"` (default; uses the rare-event
#'   rule) or `"gaussian"` (returns `n`).
#'
#' @return Integer effective sample size.
#'
#' @examples
#' y <- rbinom(500, 1, 0.1)
#' compute_n_eff(y, family = "binomial")
#'
#' @references Phillips RV, et al. (2023) Practical considerations
#'   for specifying a super learner. *Statistics in Medicine*.
#'
#' @export
compute_n_eff <- function(y, family = "binomial") {
  n <- length(y)
  if (!identical(family, "binomial")) return(as.integer(n))
  y_num <- as.integer(as.logical(y))
  n_rare <- min(sum(y_num == 1L, na.rm = TRUE),
                sum(y_num == 0L, na.rm = TRUE))
  as.integer(min(n, 5L * n_rare))
}


# ------------------------------------------------------------------
# Default V from Phillips Step 3
# ------------------------------------------------------------------

#' Recommend a Default V for V-fold Cross-Validation
#'
#' Implements the Phillips et al. 2023 rule of thumb:
#' `V = n - 1` for n_eff < 30; `V = 20` for 30 <= n_eff < 500;
#' `V = 10` for 500 <= n_eff < 5,000; `V = 5` for
#' 5,000 <= n_eff < 10,000; `V = 2` above 10,000.
#'
#' @param n_eff Effective sample size from [compute_n_eff()].
#'
#' @return Integer V.
#'
#' @examples
#' recommend_cv_V(compute_n_eff(rbinom(200, 1, 0.1)))
#'
#' @export
recommend_cv_V <- function(n_eff) {
  if (n_eff < 30L)    return(max(n_eff - 1L, 2L))
  if (n_eff < 500L)   return(20L)
  if (n_eff < 5000L)  return(10L)
  if (n_eff < 10000L) return(5L)
  2L
}


# ------------------------------------------------------------------
# SuperLearner library presets
# ------------------------------------------------------------------

#' Build a SuperLearner Library by Role and Sample Size
#'
#' Returns a role-aware SuperLearner library keyed by effective
#' sample size and dimensionality. Implements the published
#' workshop-default specifications from Gruber et al. 2023 SAP
#' supplement (Checklists C-E) and Phillips et al. 2023 Steps 2-5.
#'
#' Roles:
#' - `"Q"`: initial outcome regression. Family matches the outcome
#'   (binomial here).
#' - `"g"`: propensity-score model. Family is binomial; CV is
#'   stratified.
#' - `"Delta"`: missingness / censoring response model. Same as `g`.
#'
#' Presets (selected by `preset`; `"auto"` picks by `n_eff`):
#' - `"small_n"`: GLM + GLM with main effects + ranger-shallow.
#' - `"default"`: GLM + glmnet + GAM (Phillips workshop default).
#' - `"rich"`: GLM + glmnet + GAM + ranger + xgboost.
#' - `"very_rich"`: above + BART (via `tmle.SL.dbarts2`) + HAL
#'   (if installed).
#'
#' @param role One of `"Q"`, `"g"`, `"Delta"`.
#' @param n_eff Effective sample size (from [compute_n_eff()]).
#' @param p Number of candidate predictors.
#' @param preset One of `"auto"`, `"small_n"`, `"default"`, `"rich"`,
#'   `"very_rich"`.
#' @param include_screeners Logical; if TRUE wrap learners with
#'   `c(learner, "screen.corP")` to pair them with a correlation
#'   screener (Phillips p.1283).
#'
#' @return A list with `library` (character vector or list of
#'   `c(learner, screener)` pairs), `cv_V` (recommended V), and
#'   `stratify_cv` (logical, TRUE for `g` / `Delta`, FALSE for `Q`).
#'
#' @examples
#' build_sl_library(role = "g", n_eff = 300, p = 10, preset = "default")
#'
#' @references Phillips RV, et al. (2023); Gruber S, et al. (2023)
#'   SAP supplement Checklists C-E.
#'
#' @export
build_sl_library <- function(role = c("Q", "g", "Delta"),
                             n_eff,
                             p = NULL,
                             preset = c("auto", "small_n", "default",
                                        "rich", "very_rich"),
                             include_screeners = FALSE) {
  role   <- match.arg(role)
  preset <- match.arg(preset)
  if (preset == "auto") {
    preset <- if (n_eff < 200L) "small_n"
              else if (n_eff < 2000L) "default" else "rich"
  }

  base <- switch(preset,
    small_n   = c("SL.glm", "SL.mean", "SL.ranger"),
    default   = c("SL.glm", "SL.glmnet", "SL.gam", "SL.mean"),
    rich      = c("SL.glm", "SL.glmnet", "SL.gam", "SL.ranger",
                  "SL.xgboost", "SL.mean"),
    very_rich = c("SL.glm", "SL.glmnet", "SL.gam", "SL.ranger",
                  "SL.xgboost", "tmle.SL.dbarts2", "SL.hal9001",
                  "SL.mean")
  )

  if (isTRUE(include_screeners)) {
    base <- lapply(base, function(L) c(L, "screen.corP"))
  }

  list(
    library     = base,
    cv_V        = recommend_cv_V(n_eff),
    stratify_cv = role %in% c("g", "Delta")
  )
}


# ------------------------------------------------------------------
# Truncation rules
# ------------------------------------------------------------------

#' Resolve a Named Truncation Rule
#'
#' Accepts either a numeric truncation bound (returned as-is) or a
#' named rule. The SAP-supplement default is `"sqrt_n_ln_n"` =
#' `5 / (sqrt(n) * ln(n))`.
#'
#' @param rule Either a numeric value in `(0, 0.5)` or one of
#'   `"sqrt_n_ln_n"`, `"fixed_001"`, `"fixed_025"`, `"fixed_05"`.
#' @param n Sample size used by `"sqrt_n_ln_n"`.
#'
#' @return Numeric truncation bound.
#'
#' @examples
#' resolve_truncation_rule("sqrt_n_ln_n", n = 1000)
#' resolve_truncation_rule(0.01)
#'
#' @references Gruber S, et al. (2022) *Am J Epi*; SAP supplement
#'   §8.3.2.
#'
#' @export
resolve_truncation_rule <- function(rule, n = NULL) {
  if (is.numeric(rule)) {
    if (any(rule <= 0 | rule >= 0.5))
      stop("Numeric truncation must be in (0, 0.5).", call. = FALSE)
    return(rule)
  }
  switch(rule,
    sqrt_n_ln_n = {
      if (is.null(n) || n < 2L)
        stop("`n` is required for the sqrt_n_ln_n rule.", call. = FALSE)
      max(5 / (sqrt(n) * log(n)), 1e-4)
    },
    fixed_001 = 0.01,
    fixed_025 = 0.025,
    fixed_05  = 0.05,
    stop("Unknown truncation rule: ", rule, call. = FALSE)
  )
}


# ------------------------------------------------------------------
# Positivity diagnostics consolidation
# ------------------------------------------------------------------

#' Pre-Outcome Positivity Diagnostics
#'
#' Compiles propensity-score positivity diagnostics into a single
#' tidy object: PS C-statistic, per-arm distribution summary, the
#' proportion of subjects with truncated PS values, the bounded-g
#' quantiles, and a near-positivity-violation flag using
#' `eps(n) = max(0.01, 5 / (sqrt(n) * ln(n)))` (Gruber 2009 §3.1).
#'
#' @param ps_fit A `cleanroom_ps_fit` produced by
#'   [fit_ps_superlearner()] or [fit_ps_glm()].
#' @param truncation Numeric truncation bound; defaults to the
#'   `sqrt_n_ln_n` rule applied to the fitted sample size.
#'
#' @return A list with `c_statistic`, `summary_by_arm`,
#'   `pct_truncated`, `g_quantiles`, `eps_threshold`,
#'   `near_violation`, and `truncation`.
#'
#' @examples
#' \dontrun{
#' lock <- create_simple_lock(sim_func1(200), "treatment", "event_24",
#'                            c("age","sex","biomarker","comorbidity"))
#' ps   <- fit_ps_superlearner(lock)
#' diag <- run_positivity_diagnostics(ps)
#' diag$c_statistic
#' diag$near_violation
#' }
#'
#' @references Gruber S, van der Laan MJ (2009); Gruber S et al. (2023) RWE.
#'
#' @export
run_positivity_diagnostics <- function(ps_fit, truncation = NULL) {
  if (!inherits(ps_fit, "cleanroom_ps_fit"))
    stop("`ps_fit` must be a cleanroom_ps_fit.", call. = FALSE)
  ps <- as.numeric(ps_fit$ps)
  a  <- as.integer(ps_fit$treatment)
  n  <- length(ps)
  if (is.null(truncation))
    truncation <- resolve_truncation_rule("sqrt_n_ln_n", n)
  eps <- max(0.01, 5 / (sqrt(n) * log(n)))
  # C-statistic (AUC) for the PS as classifier of A.
  c_stat <- tryCatch({
    ord <- order(ps)
    a_ord <- a[ord]
    n1 <- sum(a_ord); n0 <- sum(1L - a_ord)
    if (n1 == 0 || n0 == 0) NA_real_
    else (sum((cumsum(1L - a_ord) * a_ord)) - n1 * (n1 + 1) / 2) /
         (n0 * n1)
  }, error = function(e) NA_real_)
  summary_by_arm <- vapply(c(0L, 1L), function(level) {
    s <- ps[a == level]
    c(n = length(s), mean = mean(s), sd = stats::sd(s),
      q05 = stats::quantile(s, 0.05, names = FALSE),
      q50 = stats::quantile(s, 0.50, names = FALSE),
      q95 = stats::quantile(s, 0.95, names = FALSE),
      min = min(s), max = max(s))
  }, numeric(8))
  colnames(summary_by_arm) <- c("control", "treated")
  bounded <- pmin(pmax(ps, truncation), 1 - truncation)
  pct_truncated <- mean(ps < truncation | ps > 1 - truncation)
  g_quantiles <- stats::quantile(bounded,
                                 c(0, 0.01, 0.05, 0.5, 0.95, 0.99, 1),
                                 names = FALSE)
  list(
    c_statistic   = c_stat,
    summary_by_arm = summary_by_arm,
    pct_truncated = pct_truncated,
    g_quantiles   = stats::setNames(g_quantiles,
                                     c("min","q01","q05","q50",
                                       "q95","q99","max")),
    eps_threshold = eps,
    near_violation = any(ps < eps) || any(ps > 1 - eps),
    truncation    = truncation
  )
}


# ------------------------------------------------------------------
# G-value: the causal-gap that flips the conclusion
# ------------------------------------------------------------------

#' Compute the G-Value of a Causal Effect Estimate
#'
#' Returns the additive bias that would push the confidence interval
#' across the null. Defined as
#' `G = min(|psi - 1.96 * se - null|, |psi + 1.96 * se - null|)`.
#' Equivalent to the smaller of the two CI-to-null distances; a
#' simple pre-QBA tipping-point readout (Gruber et al. 2023 RWE p.5).
#'
#' @param estimate Numeric point estimate.
#' @param se Numeric standard error. Either `se` or both `ci_lower`
#'   and `ci_upper` must be supplied.
#' @param ci_lower Numeric lower 95% CI bound (optional).
#' @param ci_upper Numeric upper 95% CI bound (optional).
#' @param null Numeric null value (default 0 for risk difference;
#'   use 1 for risk ratio or hazard ratio).
#'
#' @return A list with `g_value`, `direction`, and the CI bounds used.
#'
#' @examples
#' compute_G_value(estimate = 0.031, se = 0.016)
#' compute_G_value(estimate = 0.031, ci_lower = -0.001, ci_upper = 0.063)
#'
#' @references Gruber S, et al. (2023) *Evaluating and improving RWE
#'   with Targeted Learning*.
#'
#' @export
compute_G_value <- function(estimate, se = NULL,
                             ci_lower = NULL, ci_upper = NULL,
                             null = 0) {
  if (is.null(ci_lower) || is.null(ci_upper)) {
    if (is.null(se))
      stop("Provide either `se` or both `ci_lower` and `ci_upper`.",
           call. = FALSE)
    ci_lower <- estimate - 1.96 * se
    ci_upper <- estimate + 1.96 * se
  }
  d_lower <- abs(ci_lower - null)
  d_upper <- abs(ci_upper - null)
  g       <- min(d_lower, d_upper)
  list(
    g_value   = g,
    direction = if (estimate > null) "toward null"
                else if (estimate < null) "away from null"
                else "at null",
    ci_lower  = ci_lower,
    ci_upper  = ci_upper,
    null      = null
  )
}


# ------------------------------------------------------------------
# Delta-sensitivity (causal-gap)
# ------------------------------------------------------------------

#' Causal-Gap Sensitivity Curve
#'
#' Sweeps a causal-gap delta over a grid and returns the corresponding
#' shifted point estimate and 95% CI. Re-implements the SAP §8.3.4
#' delta-sensitivity reporting on top of an existing TMLE fit, with
#' no new estimation: the delta grid is applied as an additive shift
#' to the point estimate, holding the SE fixed.
#'
#' @param estimate Numeric point estimate.
#' @param se Numeric standard error.
#' @param delta_grid Numeric vector of bias values to subtract from
#'   the estimate. Defaults to a symmetric grid centred on zero that
#'   spans up to twice the half-width of the CI.
#'
#' @return A data frame with columns `delta`, `estimate_shifted`,
#'   `ci_lower`, `ci_upper`, `crosses_null`.
#'
#' @examples
#' run_delta_sensitivity(estimate = 0.031, se = 0.016)
#'
#' @references Gruber S, et al. (2023) SAP §8.3.4.
#'
#' @export
run_delta_sensitivity <- function(estimate, se, delta_grid = NULL) {
  if (is.null(delta_grid)) {
    halfw <- 1.96 * se
    delta_grid <- seq(-2 * halfw, 2 * halfw, length.out = 11L)
  }
  out <- data.frame(
    delta            = delta_grid,
    estimate_shifted = estimate - delta_grid,
    ci_lower         = (estimate - delta_grid) - 1.96 * se,
    ci_upper         = (estimate - delta_grid) + 1.96 * se,
    stringsAsFactors = FALSE
  )
  out$crosses_null <- out$ci_lower < 0 & out$ci_upper > 0
  class(out) <- c("cleantmle_delta_sensitivity", class(out))
  out
}


# ------------------------------------------------------------------
# AIPW and one-step estimators (light wrappers around TMLE nuisances)
# ------------------------------------------------------------------

#' Compute an AIPW Estimate from TMLE Nuisance Fits
#'
#' Given the propensity-score fit (`g`) and the initial outcome model
#' fit (`Q`) used inside a TMLE, returns the augmented
#' inverse-probability-weighting estimate of the average treatment
#' effect along with its IF-based standard error. Useful as a
#' candidate-grid comparator in the plasmode loop.
#'
#' @param g_fit A `tmle_mechanism` from
#'   [fit_tmle_treatment_mechanism()].
#' @param Q_fit A `tmle_mechanism` from
#'   [fit_tmle_outcome_mechanism()].
#'
#' @return A list with `estimate`, `se`, `ci_lower`, `ci_upper`.
#'
#' @export
compute_aipw <- function(g_fit, Q_fit) {
  stopifnot(inherits(g_fit, "tmle_mechanism"),
            inherits(Q_fit, "tmle_mechanism"))
  g    <- as.numeric(g_fit$ps)
  Q_a1 <- as.numeric(Q_fit$Q_a1)
  Q_a0 <- as.numeric(Q_fit$Q_a0)
  Q_aw <- as.numeric(Q_fit$Q_aw)
  A    <- as.integer(g_fit$treatment)
  Y    <- as.numeric(Q_fit$outcome)
  # AIPW: psi = mean( (A/g) * (Y - Q_aw) + Q_a1 ) -
  #             mean( ((1-A)/(1-g)) * (Y - Q_aw) + Q_a0 )
  ic_1 <- (A / g) * (Y - Q_aw) + Q_a1
  ic_0 <- ((1 - A) / (1 - g)) * (Y - Q_aw) + Q_a0
  psi  <- mean(ic_1 - ic_0, na.rm = TRUE)
  ic   <- (ic_1 - ic_0) - psi
  ic_obs <- ic[!is.na(ic)]
  n_eff  <- length(ic_obs)
  se <- if (n_eff < 2L) NA_real_ else
        sqrt(stats::var(ic_obs) / n_eff)
  list(
    estimate = psi,
    se       = se,
    ci_lower = psi - 1.96 * se,
    ci_upper = psi + 1.96 * se,
    n_eff    = n_eff
  )
}
