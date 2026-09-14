# Estimators for the estimand ladder: the complete-case ATT, the augmented
# overlap-weighted ATO, trimming with a propensity refit, and the
# implausibility guard attached to every estimate.
#
# All tmle::tmle delegation goes through one argument builder,
# .tmle_delegate_args(), so the ATE and ATT paths cannot silently diverge in
# their treatment-model specification. The Rescue.Co main analysis lost twelve
# hours of compute to an ATT refit that dropped prescreenW.g = FALSE and so
# fitted a different propensity model from the rest of the pipeline; the
# builder plus a regression test make that class of bug structural rather
# than behavioural.

# ── Shared tmle::tmle argument builder ──────────────────────────────────────

#' @keywords internal
.tmle_delegate_args <- function(lock, family = "binomial",
                                use_delta = FALSE,
                                sl_library = NULL,
                                gbound = NULL,
                                cv_folds = 10L,
                                prescreen_g = FALSE) {
  data       <- lock$data
  Y <- as.numeric(data[[lock$outcome]])
  A <- as.numeric(data[[lock$treatment]])
  W <- data[, lock$covariates, drop = FALSE]
  W <- as.data.frame(lapply(W, function(x) {
    x <- if (is.numeric(x)) x else as.numeric(as.factor(x))
    if (anyNA(x)) x[is.na(x)] <- stats::median(x, na.rm = TRUE)
    x
  }))
  if (is.null(sl_library)) {
    spec <- lock$primary_tmle_spec
    sl_library <- if (!is.null(spec) && !is.null(spec$q_library))
      spec$q_library else lock$sl_library
  }
  Delta <- as.integer(!is.na(Y))
  if (is.character(gbound))
    gbound <- resolve_truncation_rule(gbound, n = nrow(data))

  if (use_delta) {
    Yi <- Y
    Yi[Delta == 0L] <- if (identical(family, "binomial")) 0 else
      mean(Y, na.rm = TRUE)
    args <- list(Y = Yi, A = A, W = W, Delta = Delta, family = family,
                 Q.SL.library = sl_library, g.SL.library = sl_library,
                 g.Delta.SL.library = sl_library,
                 V.Q = cv_folds, V.g = cv_folds, V.Delta = cv_folds,
                 cvQinit = TRUE, verbose = FALSE)
  } else {
    k <- Delta == 1L
    args <- list(Y = Y[k], A = A[k], W = W[k, , drop = FALSE],
                 family = family,
                 Q.SL.library = sl_library, g.SL.library = sl_library,
                 V.Q = cv_folds, V.g = cv_folds,
                 cvQinit = TRUE, verbose = FALSE)
  }
  # prescreenW.g = FALSE deviates from the tmle package default: the default
  # LASSO-screens the propensity covariates on the outcome, and in a cohort
  # with few events it can drop a true confounder (the Rescue.Co screen
  # dropped months_since_start). Screening still happens inside the
  # SuperLearner library where its uncertainty reaches the influence curve.
  if ("prescreenW.g" %in% names(formals(tmle::tmle)))
    args$prescreenW.g <- isTRUE(prescreen_g)
  if (!is.null(gbound)) args$gbound <- gbound
  args
}

#' @keywords internal
.rubin_combine <- function(psi, vars) {
  K  <- length(psi)
  est <- mean(psi)
  se  <- sqrt(mean(vars) + (1 + 1 / K) * (if (K > 1) stats::var(psi) else 0))
  c(est = est, se = se, lo = est - 1.96 * se, hi = est + 1.96 * se,
    p = 2 * stats::pnorm(-abs(est / se)))
}

# Greedy 1:1 nearest-neighbour caliper match on the logit propensity.
# Deterministic given the RNG state (treated units are visited in a random
# order); returns index vectors into the original rows.
#' @keywords internal
.greedy_caliper_match <- function(g, A, caliper_sd = 0.2,
                                  order = c("random", "sorted")) {
  # `order` controls the greedy pass over treated units: "random" (the
  # estimator paths, under a caller-controlled seed) or "sorted"
  # (deterministic, hardest-to-match first by descending logit; used by
  # the design diagnostics so they need no RNG state at all).
  order <- match.arg(order)
  lp <- stats::qlogis(pmin(pmax(g, 1e-6), 1 - 1e-6))
  cal <- caliper_sd * stats::sd(lp)
  tr <- which(A == 1); ct <- which(A == 0)
  if (!length(tr) || !length(ct))
    return(list(treated = integer(0), control = integer(0)))
  used <- rep(FALSE, length(ct))
  m_t <- integer(0); m_c <- integer(0)
  tr_seq <- if (order == "sorted") tr[base::order(-lp[tr])] else sample(tr)
  for (i in tr_seq) {
    d <- abs(lp[ct] - lp[i]); d[used] <- Inf
    j <- which.min(d)
    if (is.finite(d[j]) && d[j] <= cal) {
      used[j] <- TRUE
      m_t <- c(m_t, i); m_c <- c(m_c, ct[j])
    }
  }
  list(treated = m_t, control = m_c)
}


#' Implausibility Guard for a Treatment-Effect Estimate
#'
#' A mechanical bounds check that flags, and never suppresses, an estimate
#' the observed data cannot support: a sign that differs from the crude
#' difference, a magnitude above five times the crude difference, a risk
#' difference beyond the largest arm rate for a binary outcome, or an effect
#' beyond the observed range for a continuous one. Noise floors keep a
#' near-zero crude difference from manufacturing a flag. Extreme propensity
#' weights wreck a point estimate while leaving its influence-curve interval
#' narrow, so nothing in an estimator's own output separates such artifacts
#' from real findings; this guard and the support verdict catch different
#' failures and both travel with every estimate the package reports.
#'
#' @param est The estimated difference (ATE, ATT, ATO).
#' @param y_obs Outcome vector, rows with an observed outcome.
#' @param a_obs Treatment indicator for the same rows.
#' @param family `"binomial"` or anything else (treated as continuous).
#' @return A list with `crude_diff`, `implausible`, `implausible_reason`.
#' @keywords internal
implausibility_check <- function(est, y_obs, a_obs, family = "gaussian") {
  none <- list(crude_diff = NA_real_, implausible = FALSE,
               implausible_reason = NA_character_)
  ok <- !is.na(y_obs) & !is.na(a_obs)
  y_obs <- y_obs[ok]; a_obs <- a_obs[ok]
  if (!length(y_obs) || length(unique(a_obs)) < 2 || !is.finite(est))
    return(none)
  m1 <- mean(y_obs[a_obs == 1]); m0 <- mean(y_obs[a_obs == 0])
  cd <- m1 - m0
  sdy  <- stats::sd(y_obs)
  span <- suppressWarnings(diff(range(y_obs)))
  noise <- if (identical(family, "binomial")) 0.005 else
    max(0.05 * sdy, 1e-12, na.rm = TRUE)
  f <- character()
  if (is.finite(cd)) {
    if (abs(est) > noise && abs(cd) > noise && sign(est) != sign(cd))
      f <- c(f, sprintf("sign differs from crude (%.4g)", cd))
    if (abs(cd) > noise && abs(est) > 5 * abs(cd))
      f <- c(f, sprintf("%.1fx crude difference", abs(est / cd)))
  }
  cap <- suppressWarnings(max(m1, m0, na.rm = TRUE))
  if (identical(family, "binomial") && is.finite(cap) && abs(est) > cap)
    f <- c(f, sprintf("risk difference exceeds largest arm rate (%.4g)", cap))
  if (!identical(family, "binomial") && is.finite(span) && span > 0 &&
      abs(est) > span)
    f <- c(f, "exceeds observed outcome range")
  list(crude_diff = cd, implausible = length(f) > 0L,
       implausible_reason = if (length(f)) paste(f, collapse = "; ") else
         NA_character_)
}


#' Complete-Case ATT via TMLE, Never Through the Censoring Mechanism
#'
#' Estimates the effect among the treated on the rows with an observed
#' outcome, named as its own estimand:
#' `E[Y(1) - Y(0) | A = 1, outcome observed]`. The delegation to
#' [tmle::tmle()] never passes `Delta`: under the censoring mechanism the
#' package ATT loses double robustness (with the outcome model misspecified
#' its bias is large and coverage collapses, while the ATE path is
#' unaffected), so the complete-case fit is consistent for its own estimand
#' with no missingness model at all. The ATE from the same fit is returned
#' beside the ATT; with good overlap the two agree, and their gap under poor
#' overlap is itself the diagnostic. The treatment-model specification is
#' built by the same internal builder as every other delegation in the
#' package, so the ATT's propensity model cannot silently differ from the
#' ATE's.
#'
#' @param lock A `cleanroom_lock`.
#' @param family `"binomial"` or `"gaussian"`.
#' @param sl_library SuperLearner library for Q and g; defaults to the locked
#'   candidate's library, else the lock library.
#' @param gbound Numeric bound, named rule (see
#'   [resolve_truncation_rule()]), or `NULL` for the tmle package default.
#' @param cv_folds V for Q and g cross-validation. Default 10.
#' @param prescreen_g Passed to `tmle::tmle(prescreenW.g = )` when the
#'   installed version supports it. Default FALSE; see the builder's note.
#' @param reps Seed-averaging replicates; distinct seeds are drawn from
#'   `seed` and combined with the multiple-imputation rule
#'   `se^2 = mean(var) + (1 + 1/K) var(psi)`. Default 1.
#' @param seed Base seed. Defaults to `lock$seed`.
#' @param allow_outcome_access Bypass the outcome guard. Default FALSE.
#' @return An object of class `att_fit` (inherits `cr_result`): the ATT with
#'   SE, CI and p-value, the same-fit ATE, the named estimand, n, and the
#'   implausibility flags.
#' @keywords internal
run_att_tmle <- function(lock,
                         family = "binomial",
                         sl_library = NULL,
                         gbound = NULL,
                         cv_folds = 10L,
                         prescreen_g = FALSE,
                         reps = 1L,
                         seed = NULL,
                         allow_outcome_access = FALSE) {
  if (!inherits(lock, "cleanroom_lock"))
    stop("`lock` must be a cleanroom_lock object.", call. = FALSE)
  .check_outcome_access(lock, allow_outcome_access, caller = "run_att_tmle")
  if (!requireNamespace("tmle", quietly = TRUE))
    stop("Package 'tmle' is required for run_att_tmle().", call. = FALSE)
  if (is.null(seed)) seed <- lock$seed

  args <- .tmle_delegate_args(lock, family = family, use_delta = FALSE,
                              sl_library = sl_library, gbound = gbound,
                              cv_folds = cv_folds, prescreen_g = prescreen_g)
  n_cc <- length(args$Y)
  if (n_cc < 10L)
    stop("run_att_tmle: fewer than 10 observed outcomes.", call. = FALSE)

  fits <- list()
  for (r in seq_len(reps)) {
    set.seed(seed + r - 1L)
    f <- tryCatch(do.call(tmle::tmle, args), error = function(e) {
      warning("run_att_tmle: tmle fit failed (", conditionMessage(e), ").",
              call. = FALSE); NULL })
    if (!is.null(f)) fits[[length(fits) + 1L]] <- f
  }
  if (!length(fits))
    stop("run_att_tmle: no tmle fit converged.", call. = FALSE)

  att_ok <- all(vapply(fits, function(f)
    !is.null(f$estimates$ATT$psi), logical(1)))
  if (!att_ok)
    stop("run_att_tmle: the installed tmle version did not return an ATT.",
         call. = FALSE)
  att <- .rubin_combine(
    vapply(fits, function(f) f$estimates$ATT$psi, numeric(1)),
    vapply(fits, function(f) f$estimates$ATT$var.psi, numeric(1)))
  ate <- .rubin_combine(
    vapply(fits, function(f) f$estimates$ATE$psi, numeric(1)),
    vapply(fits, function(f) f$estimates$ATE$var.psi, numeric(1)))
  conv <- all(vapply(fits, function(f)
    isTRUE(f$estimates$ATT$converged) || is.null(f$estimates$ATT$converged),
    logical(1)))
  if (!conv)
    warning("run_att_tmle: the ATT targeting step did not converge on every ",
            "replicate; treat with caution.", call. = FALSE)

  guard <- implausibility_check(unname(att["est"]), args$Y, args$A, family)

  out <- list(
    estimate = unname(att["est"]), se = unname(att["se"]),
    ci_lower = unname(att["lo"]), ci_upper = unname(att["hi"]),
    p_value = unname(att["p"]),
    estimand = "E[Y(1)-Y(0) | A=1, outcome observed]",
    ate_same_fit = list(estimate = unname(ate["est"]), se = unname(ate["se"]),
                        ci_lower = unname(ate["lo"]),
                        ci_upper = unname(ate["hi"])),
    n = n_cc,
    crude_diff = guard$crude_diff,
    implausible = guard$implausible,
    implausible_reason = guard$implausible_reason,
    spec = args[setdiff(names(args), c("Y", "A", "W", "Delta"))],
    reps = length(fits),
    tmle_fit = fits[[1]],
    treatment = lock$treatment, outcome = lock$outcome,
    type = "att_tmle", call = match.call())
  class(out) <- c("att_fit", "cr_result")
  out
}

#' @export
print.att_fit <- function(x, ...) {
  cat("Complete-case ATT (TMLE, never through Delta)\n")
  cat(sprintf("  Estimand: %s\n", x$estimand))
  cat(sprintf("  ATT: %.4f (95%% CI %.4f, %.4f), SE %.4f, p = %.4f, n = %d\n",
              x$estimate, x$ci_lower, x$ci_upper, x$se, x$p_value, x$n))
  cat(sprintf("  ATE from the same fit: %.4f (%.4f, %.4f)\n",
              x$ate_same_fit$estimate, x$ate_same_fit$ci_lower,
              x$ate_same_fit$ci_upper))
  if (isTRUE(x$implausible))
    cat("  IMPLAUSIBLE:", x$implausible_reason, "\n")
  invisible(x)
}


#' Augmented Overlap-Weighted Estimator (ATO)
#'
#' The overlap-weighted average treatment effect of Li, Morgan and Zaslavsky
#' (2018): treated units weighted by 1 - g, controls by g, so the target
#' population is the region of clinical equipoise, weights are bounded by
#' one, no one is discarded, and under a logistic propensity model the
#' weighted covariate means balance exactly. This implementation is the
#' augmented (doubly robust) estimator of Mao, Li and Greene (2019), with an
#' influence-function variance and the exact-balance check reported. The
#' unaugmented Hajek point estimate is returned beside it. Estimation is on
#' rows with an observed outcome.
#'
#' @param lock A `cleanroom_lock`.
#' @param ps_fit A `ps_fit`; its untruncated scores are used.
#' @param sl_library SuperLearner library for the outcome model; defaults to
#'   the locked candidate's library, else the lock library. GLM fallback
#'   when SuperLearner is unavailable.
#' @param family `"binomial"` or `"gaussian"`.
#' @param allow_outcome_access Bypass the outcome guard. Default FALSE.
#' @return An object of class `ato_fit` (inherits `cr_result`): estimate,
#'   IF-based SE, CI, p-value, the Hajek point estimate, the maximum
#'   absolute overlap-weighted SMD across the lock covariates, n, and the
#'   implausibility flags.
#' @references Li F, Morgan KL, Zaslavsky AM (2018) JASA 113:390-400.
#'   Mao H, Li L, Greene T (2019) Stat Methods Med Res 28:2439-2454.
#' @keywords internal
estimate_ato <- function(lock, ps_fit,
                         sl_library = NULL,
                         family = "binomial",
                         allow_outcome_access = FALSE) {
  if (!inherits(lock, "cleanroom_lock"))
    stop("`lock` must be a cleanroom_lock object.", call. = FALSE)
  if (!inherits(ps_fit, "ps_fit"))
    stop("`ps_fit` must be a ps_fit object.", call. = FALSE)
  .check_outcome_access(lock, allow_outcome_access, caller = "estimate_ato")

  data <- lock$data
  Y <- as.numeric(data[[lock$outcome]])
  A <- as.integer(data[[lock$treatment]])
  g <- pmin(pmax(as.numeric(ps_fit$ps_raw %||% ps_fit$ps), 1e-6), 1 - 1e-6)
  W <- data[, lock$covariates, drop = FALSE]
  W <- as.data.frame(lapply(W, function(x) {
    x <- if (is.numeric(x)) x else as.numeric(as.factor(x))
    if (anyNA(x)) x[is.na(x)] <- stats::median(x, na.rm = TRUE)
    x
  }))

  cc <- !is.na(Y)
  if (sum(cc) < 10L)
    stop("estimate_ato: fewer than 10 observed outcomes.", call. = FALSE)
  Yc <- Y[cc]; Ac <- A[cc]; gc_ <- g[cc]; Wc <- W[cc, , drop = FALSE]
  n <- length(Yc)

  if (is.null(sl_library)) {
    spec <- lock$primary_tmle_spec
    sl_library <- if (!is.null(spec) && !is.null(spec$q_library))
      spec$q_library else lock$sl_library
  }
  fam_obj <- if (identical(family, "binomial")) stats::binomial() else
    stats::gaussian()
  q <- tryCatch({
    if (requireNamespace("SuperLearner", quietly = TRUE)) {
      set.seed(lock$seed + 31L)
      AW <- cbind(data.frame(.A = Ac), Wc)
      sl <- SuperLearner::SuperLearner(Y = Yc, X = AW, family = fam_obj,
                                       SL.library = sl_library,
                                       env = .cleantmle_sl_env())
      A1 <- AW; A1$.A <- 1L
      A0 <- AW; A0$.A <- 0L
      list(Q1 = as.numeric(predict(sl, newdata = A1)$pred),
           Q0 = as.numeric(predict(sl, newdata = A0)$pred))
    } else stop("no SuperLearner")
  }, error = function(e) {
    df <- cbind(data.frame(.Y = Yc, .A = Ac), Wc)
    fit <- stats::glm(.Y ~ ., data = df, family = fam_obj)
    d1 <- df; d1$.A <- 1L
    d0 <- df; d0$.A <- 0L
    list(Q1 = as.numeric(stats::predict(fit, newdata = d1,
                                        type = "response")),
         Q0 = as.numeric(stats::predict(fit, newdata = d0,
                                        type = "response")))
  })

  h <- gc_ * (1 - gc_)
  num <- h * (q$Q1 - q$Q0) +
    Ac * (1 - gc_) * (Yc - q$Q1) - (1 - Ac) * gc_ * (Yc - q$Q0)
  hbar <- mean(h)
  psi  <- mean(num) / hbar
  phi  <- (num - h * psi) / hbar
  se   <- sqrt(stats::var(phi) / n)

  # Unaugmented Hajek point estimate and the overlap-weighted arm means
  # (the adjusted arm risks sensitivity analyses should read).
  w1 <- (1 - gc_)[Ac == 1]; w0 <- gc_[Ac == 0]
  arm1 <- stats::weighted.mean(Yc[Ac == 1], w1)
  arm0 <- stats::weighted.mean(Yc[Ac == 0], w0)
  hajek <- arm1 - arm0
  risk1_aug <- sum(h * q$Q1 + Ac * (1 - gc_) * (Yc - q$Q1)) / sum(h)
  risk0_aug <- sum(h * q$Q0 + (1 - Ac) * gc_ * (Yc - q$Q0)) / sum(h)

  # Exact-balance check: overlap-weighted SMDs across the lock covariates
  # (exact under a logistic propensity model; approximate under an ensemble).
  wts <- ifelse(Ac == 1, 1 - gc_, gc_)
  wsmd <- vapply(names(Wc), function(v) {
    xv <- Wc[[v]]
    m1 <- stats::weighted.mean(xv[Ac == 1], wts[Ac == 1])
    m0 <- stats::weighted.mean(xv[Ac == 0], wts[Ac == 0])
    s  <- sqrt((stats::var(xv[Ac == 1]) + stats::var(xv[Ac == 0])) / 2)
    if (is.finite(s) && s > 0) (m1 - m0) / s else 0
  }, numeric(1))

  guard <- implausibility_check(psi, Yc, Ac, family)

  out <- list(
    estimate = psi, se = se,
    ci_lower = psi - 1.96 * se, ci_upper = psi + 1.96 * se,
    p_value = 2 * stats::pnorm(-abs(psi / se)),
    estimand = "ATO (overlap-weighted average treatment effect), outcome-observed rows",
    hajek_estimate = hajek,
    risk_treated = risk1_aug, risk_control = risk0_aug,
    max_abs_weighted_smd = max(abs(wsmd)),
    weighted_smds = wsmd,
    n = n,
    crude_diff = guard$crude_diff,
    implausible = guard$implausible,
    implausible_reason = guard$implausible_reason,
    treatment = lock$treatment, outcome = lock$outcome,
    type = "ato", call = match.call())
  class(out) <- c("ato_fit", "cr_result")
  out
}

#' @export
print.ato_fit <- function(x, ...) {
  cat("Augmented overlap-weighted estimator (ATO)\n")
  cat(sprintf("  Estimand: %s\n", x$estimand))
  cat(sprintf("  ATO: %.4f (95%% CI %.4f, %.4f), SE %.4f, p = %.4f, n = %d\n",
              x$estimate, x$ci_lower, x$ci_upper, x$se, x$p_value, x$n))
  cat(sprintf("  Hajek (unaugmented) point estimate: %.4f\n",
              x$hajek_estimate))
  cat(sprintf("  Max |overlap-weighted SMD|: %.3f %s\n",
              x$max_abs_weighted_smd,
              if (x$max_abs_weighted_smd <= 0.1) "(balanced)" else
                "(check balance)"))
  if (isTRUE(x$implausible))
    cat("  IMPLAUSIBLE:", x$implausible_reason, "\n")
  invisible(x)
}


# Crump et al. (2009) data-adaptive trimming threshold: the smallest alpha
# such that 1 / (alpha (1 - alpha)) >= 2 E[ h | h <= 1 / (alpha (1 - alpha)) ]
# with h = 1 / (g (1 - g)).
#' @keywords internal
.crump_alpha <- function(g, grid = seq(0.001, 0.2, by = 0.001)) {
  h <- 1 / (g * (1 - g))
  for (a in grid) {
    lam <- 1 / (a * (1 - a))
    sub <- h[h <= lam]
    if (length(sub) && lam >= 2 * mean(sub)) return(a)
  }
  max(grid)
}


#' Trimmed ATE with a Propensity Refit
#'
#' Trims to the common-support band, refits the propensity model on the
#' trimmed subset (reusing the full-cohort fit would describe the overlap of
#' a population that no longer exists), reassesses support, and estimates
#' the ATE on the subset only when the refit clears the gate; otherwise the
#' next, tighter level is tried. Trimming changes the estimand: the result
#' is the effect among patients whose covariates could plausibly have
#' produced either arm, and every row is labelled with that population.
#'
#' @param lock A `cleanroom_lock`.
#' @param ps_fit A `ps_fit` on the full cohort; its untruncated scores define
#'   the trim.
#' @param levels Trim levels tried in order (g in `[l, 1 - l]`). Default
#'   `c(0.05, 0.10)`.
#' @param rule `"fixed"` uses `levels`; `"crump"` derives the level from the
#'   data (Crump et al. 2009) and prepends it.
#' @param family,sl_library,gbound,cv_folds,prescreen_g,seed As in
#'   [run_att_tmle()].
#' @param use_ipcw Estimate with the censoring mechanism (`Delta`) on the
#'   trimmed subset. Default FALSE (complete case).
#' @param thresholds [support_thresholds()] used to judge the refit.
#' @param min_n,min_per_arm Feasibility floor after trimming. Defaults 100
#'   and 30.
#' @param allow_outcome_access Bypass the outcome guard. Default FALSE.
#' @param verbose Print the per-level path. Default TRUE.
#' @return An object of class `trimmed_tmle_fit`: the estimate at the first
#'   level that cleared the gate, dropped counts by arm, the refit's support
#'   assessment, the population label, the implausibility flags, and the
#'   record of levels tried.
#' @references Crump RK, Hotz VJ, Imbens GW, Mitnik OA (2009) Biometrika
#'   96:187-199.
#' @keywords internal
run_trimmed_tmle <- function(lock, ps_fit,
                             levels = c(0.05, 0.10),
                             rule = c("fixed", "crump"),
                             family = "binomial",
                             sl_library = NULL,
                             gbound = NULL,
                             cv_folds = 10L,
                             prescreen_g = FALSE,
                             use_ipcw = FALSE,
                             thresholds = support_thresholds(),
                             min_n = 100L,
                             min_per_arm = 30L,
                             seed = NULL,
                             allow_outcome_access = FALSE,
                             verbose = TRUE) {
  if (!inherits(lock, "cleanroom_lock"))
    stop("`lock` must be a cleanroom_lock object.", call. = FALSE)
  if (!inherits(ps_fit, "ps_fit"))
    stop("`ps_fit` must be a ps_fit object.", call. = FALSE)
  .check_outcome_access(lock, allow_outcome_access,
                        caller = "run_trimmed_tmle")
  rule <- match.arg(rule)
  if (is.null(seed)) seed <- lock$seed

  g0 <- pmin(pmax(as.numeric(ps_fit$ps_raw %||% ps_fit$ps), 1e-6), 1 - 1e-6)
  A  <- as.integer(lock$data[[lock$treatment]])
  if (rule == "crump") {
    a_star <- .crump_alpha(g0)
    levels <- unique(c(a_star, levels))
    if (verbose) message("run_trimmed_tmle: Crump threshold alpha = ",
                         format(a_star))
  }

  path <- list()
  for (lo in levels) {
    keep <- g0 >= lo & g0 <= (1 - lo)
    dropped_t <- sum(!keep & A == 1); dropped_c <- sum(!keep & A == 0)
    lab <- sprintf("trimmed g in [%.2f, %.2f]", lo, 1 - lo)
    if (sum(keep) < min_n || sum(A[keep] == 1) < min_per_arm ||
        sum(A[keep] == 0) < min_per_arm) {
      path[[length(path) + 1L]] <- list(level = lo, status = "too little left")
      if (verbose) message("  [", lab, "] too little left after trimming; ",
                           "skipped")
      next
    }
    sub_lock <- lock
    sub_lock$data <- lock$data[keep, , drop = FALSE]

    # Refit the propensity on the trimmed subset with the same method.
    refit <- tryCatch({
      if (identical(ps_fit$method, "glm")) fit_ps_glm(sub_lock)
      else fit_ps_superlearner(sub_lock,
                               cv_folds = ps_fit$cv_folds %||% 10L)
    }, error = function(e) NULL)
    if (is.null(refit)) {
      path[[length(path) + 1L]] <- list(level = lo, status = "refit failed")
      if (verbose) message("  [", lab, "] propensity refit failed; skipped")
      next
    }
    sup2 <- assess_support(refit, thresholds = thresholds,
                           tree_search = FALSE)
    if (verbose)
      message(sprintf("  [%s] dropped %d treated / %d control; refit -> %s",
                      lab, dropped_t, dropped_c, sup2$verdict))
    if (sup2$verdict %in% c("SEVERE", "FAIL")) {
      path[[length(path) + 1L]] <- list(level = lo, status = sup2$verdict)
      next
    }

    args <- .tmle_delegate_args(sub_lock, family = family,
                                use_delta = use_ipcw,
                                sl_library = sl_library, gbound = gbound,
                                cv_folds = cv_folds,
                                prescreen_g = prescreen_g)
    set.seed(seed)
    fit <- tryCatch(do.call(tmle::tmle, args), error = function(e) {
      warning("run_trimmed_tmle: tmle fit failed at level ", lo, " (",
              conditionMessage(e), ").", call. = FALSE); NULL })
    if (is.null(fit)) {
      path[[length(path) + 1L]] <- list(level = lo, status = "tmle failed")
      next
    }
    ate <- fit$estimates$ATE
    guard <- implausibility_check(unname(ate$psi), args$Y, args$A, family)
    out <- list(
      estimate = unname(ate$psi),
      se = unname(sqrt(ate$var.psi)),
      ci_lower = unname(ate$CI[1]), ci_upper = unname(ate$CI[2]),
      p_value = unname(ate$pvalue),
      estimand = paste0("ATE on the common-support population, ", lab),
      population = lab,
      level = lo,
      n = sum(keep),
      n_dropped_treated = dropped_t,
      n_dropped_control = dropped_c,
      support = sup2,
      verdict = sup2$verdict,
      caveat = sup2$caveat,
      used_ipcw = isTRUE(use_ipcw),
      crude_diff = guard$crude_diff,
      implausible = guard$implausible,
      implausible_reason = guard$implausible_reason,
      levels_tried = path,
      tmle_fit = fit,
      treatment = lock$treatment, outcome = lock$outcome,
      type = "trimmed_tmle", call = match.call())
    class(out) <- c("trimmed_tmle_fit", "cr_result")
    return(out)
  }
  stop("run_trimmed_tmle: no trim level cleared the gate (",
       paste(vapply(path, function(p)
         paste0(p$level, ": ", p$status), character(1)), collapse = "; "),
       ").", call. = FALSE)
}

#' @export
print.trimmed_tmle_fit <- function(x, ...) {
  cat("Trimmed ATE with propensity refit\n")
  cat(sprintf("  Population: %s (dropped %d treated / %d control)\n",
              x$population, x$n_dropped_treated, x$n_dropped_control))
  cat(sprintf("  ATE: %.4f (95%% CI %.4f, %.4f), SE %.4f, p = %.4f, n = %d\n",
              x$estimate, x$ci_lower, x$ci_upper, x$se, x$p_value, x$n))
  cat(sprintf("  Refit support verdict: %s\n", x$verdict))
  if (isTRUE(x$implausible))
    cat("  IMPLAUSIBLE:", x$implausible_reason, "\n")
  invisible(x)
}
