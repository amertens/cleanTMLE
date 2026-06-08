#' Data-Quality Stress Testing via Plasmode Simulation
#'
#' Extends the plasmode feasibility framework to evaluate TMLE
#' candidate robustness under controlled data-quality degradations:
#' covariate missingness (missing completely at random, missing at
#' random, or missing not at random, with a choice of imputation
#' handler), treatment misclassification (with asymmetric sensitivity
#' and specificity), outcome misclassification (optionally differential
#' by treatment arm), unmeasured confounding, and near-positivity
#' (amplification of the covariate-to-treatment association toward
#' deterministic treatment).
#'
#' @name plasmode_dq
NULL


#' Default DQ Stress Scenario Configurations
#'
#' Returns a recommended scenario list for [run_plasmode_dq_stress()],
#' tunable by preset.  Saves the user from constructing the verbose
#' nested list by hand for routine cases.
#'
#' @param preset Character; one of `"regulatory_standard"` (five threats,
#'   moderate severities; the default), `"exploratory"` (lighter, faster),
#'   or `"stress"` (heavier severities).
#' @return A named list suitable for the `data_quality_scenarios`
#'   argument of [run_plasmode_dq_stress()]. Every preset returns the
#'   same five threats: `covariate_missingness`, `treatment_misclass`,
#'   `outcome_misclass`, `unmeasured_confounding`, and `near_positivity`.
#'   The `near_positivity` element carries a `slopes` vector of
#'   multipliers (greater than 1) applied to the centred propensity-score
#'   log-odds, so a subgroup approaches deterministic treatment and the
#'   estimated propensity scores reach the boundary. This is the same
#'   covariate-to-treatment amplification used by the divergence study.
#' @examples
#' default_dq_scenarios()
#' default_dq_scenarios("exploratory")
#' @export
default_dq_scenarios <- function(preset = c("regulatory_standard",
                                             "exploratory", "stress")) {
  preset <- match.arg(preset)
  switch(preset,
    regulatory_standard = list(
      covariate_missingness = list(fractions = c(0.05, 0.10, 0.20)),
      treatment_misclass    = list(sensitivity = c(0.95, 0.90),
                                   specificity = c(0.99, 0.95)),
      outcome_misclass      = list(sensitivity = c(0.95, 0.90),
                                   specificity = c(0.99, 0.95)),
      unmeasured_confounding = list(U_prevalence  = 0.20,
                                    U_treatment_OR = c(1.5, 2.0),
                                    U_outcome_OR   = c(1.5, 2.0)),
      near_positivity       = list(slopes = c(1.5, 2.0))
    ),
    exploratory = list(
      covariate_missingness = list(fractions = c(0.10)),
      treatment_misclass    = list(sensitivity = 0.90, specificity = 0.95),
      outcome_misclass      = list(sensitivity = 0.90, specificity = 0.95),
      unmeasured_confounding = list(U_prevalence  = 0.20,
                                    U_treatment_OR = 1.5,
                                    U_outcome_OR   = 1.5),
      near_positivity       = list(slopes = 1.5)
    ),
    stress = list(
      covariate_missingness = list(fractions = c(0.10, 0.20, 0.40)),
      treatment_misclass    = list(sensitivity = c(0.90, 0.80),
                                   specificity = c(0.95, 0.90)),
      outcome_misclass      = list(sensitivity = c(0.90, 0.80),
                                   specificity = c(0.95, 0.90)),
      unmeasured_confounding = list(U_prevalence  = 0.20,
                                    U_treatment_OR = c(2.0, 3.0),
                                    U_outcome_OR   = c(2.0, 3.0)),
      near_positivity       = list(slopes = c(2.0, 3.0, 4.0))
    )
  )
}


#' Print the Active Locked TMLE Specification
#'
#' Prints (or returns invisibly) the locked primary TMLE specification
#' carried on a `cleanroom_lock`: candidate id, label, g/Q libraries,
#' truncation, plasmode RMSE that won the selection, and the lock hash.
#' Useful for verifying that downstream estimator calls will pick up the
#' right candidate.
#'
#' @param lock A `cleanroom_lock`.
#' @return Invisibly returns the locked spec (or `NULL` if none).
#' @examples
#' \dontrun{
#' print_locked_spec(lock)
#' }
#' @export
print_locked_spec <- function(lock) {
  if (!inherits(lock, "cleanroom_lock"))
    stop("`lock` must be a cleanroom_lock object.", call. = FALSE)
  spec <- lock$primary_tmle_spec
  if (is.null(spec)) {
    cat("No primary TMLE specification locked. Use ",
        "lock_primary_tmle_spec() after select_tmle_candidate().\n",
        sep = "")
    return(invisible(NULL))
  }
  cat("Locked Primary TMLE Specification\n")
  cat("=================================\n")
  cat(sprintf("Candidate ID:    %s\n", spec$candidate_id %||% NA))
  cat(sprintf("Label:           %s\n", spec$label        %||% NA))
  cat(sprintf("g-library:       %s\n",
              paste(spec$g_library, collapse = ", ")))
  cat(sprintf("Q-library:       %s\n",
              paste(spec$q_library, collapse = ", ")))
  cat(sprintf("Truncation:      %s\n", as.character(spec$truncation)))
  if (!is.null(spec$selection_rule))
    cat(sprintf("Selection rule:  %s\n", spec$selection_rule))
  if (!is.null(spec$metrics) && is.data.frame(spec$metrics)) {
    if ("rmse" %in% names(spec$metrics))
      cat(sprintf("Plasmode RMSE:   %.5f\n", spec$metrics$rmse[1]))
    if ("coverage" %in% names(spec$metrics))
      cat(sprintf("Plasmode cov:    %.3f\n", spec$metrics$coverage[1]))
  }
  cat(sprintf("Lock hash:       %s\n", lock$lock_hash %||% NA))
  invisible(spec)
}


# ── Internal: Single-Rep TMLE Fit ────────────────────────────────────────

#' Fit one plasmode replicate for a single TMLE candidate
#'
#' @param Y_sim Simulated outcome vector.
#' @param A Treatment vector.
#' @param W_data Covariate data.frame (may be degraded).
#' @param treatment Character; treatment column name.
#' @param covariates Character vector of covariate names.
#' @param cand A \code{tmle_candidate_spec}.
#' @param n Sample size.
#'
#' @return List with est, se, ci_lower, ci_upper.
#' @keywords internal
.plasmode_fit_one_candidate <- function(Y_sim, A, W_data, treatment,
                                         covariates, cand, n) {
  g_lib <- cand$g_library
  q_lib <- cand$q_library

  # ── g-model ──
  use_sl <- requireNamespace("SuperLearner", quietly = TRUE) &&
    !identical(g_lib, "SL.glm")

  if (use_sl) {
    g_sl   <- SuperLearner::SuperLearner(
      Y = A, X = W_data[, covariates, drop = FALSE],
      family = binomial(), SL.library = g_lib,
      env = .cleantmle_sl_env()
    )
    ps_hat <- as.numeric(g_sl$SL.predict)
  } else {
    ps_fml <- stats::reformulate(covariates, response = treatment)
    ds_g   <- W_data; ds_g[[treatment]] <- A
    ps_mod <- stats::glm(ps_fml, data = ds_g, family = stats::binomial())
    ps_hat <- as.numeric(stats::predict(ps_mod, type = "response"))
  }
  ps_hat <- pmax(pmin(ps_hat, 1 - cand$truncation), cand$truncation)

  # ── Q-model ──
  ds <- W_data
  ds[[treatment]] <- A
  ds[[".Y_sim."]] <- Y_sim
  AW <- ds[, c(treatment, covariates), drop = FALSE]

  use_sl_q <- requireNamespace("SuperLearner", quietly = TRUE) &&
    !identical(q_lib, "SL.glm")

  if (use_sl_q) {
    q_sl <- SuperLearner::SuperLearner(
      Y = Y_sim, X = AW, family = binomial(), SL.library = q_lib,
      env = .cleantmle_sl_env()
    )
    AW_a1 <- AW; AW_a1[[treatment]] <- 1L
    AW_a0 <- AW; AW_a0[[treatment]] <- 0L
    Q_a1  <- as.numeric(predict(q_sl, newdata = AW_a1)$pred)
    Q_a0  <- as.numeric(predict(q_sl, newdata = AW_a0)$pred)
    Q_aw  <- as.numeric(q_sl$SL.predict)
  } else {
    Q_fml   <- stats::reformulate(c(treatment, covariates),
                                   response = ".Y_sim.")
    Q_fit_s <- stats::glm(Q_fml, data = ds, family = stats::binomial())
    ds_a1 <- ds; ds_a1[[treatment]] <- 1L
    ds_a0 <- ds; ds_a0[[treatment]] <- 0L
    Q_a1 <- as.numeric(stats::predict(Q_fit_s, newdata = ds_a1,
                                       type = "response"))
    Q_a0 <- as.numeric(stats::predict(Q_fit_s, newdata = ds_a0,
                                       type = "response"))
    Q_aw <- as.numeric(stats::predict(Q_fit_s, type = "response"))
  }

  # ── Targeting ──
  H_a1 <-  1 / ps_hat
  H_a0 <- -1 / (1 - ps_hat)
  H_aw <- ifelse(A == 1, H_a1, H_a0)

  Q_aw_logit <- stats::qlogis(pmax(pmin(Q_aw, 0.999), 0.001))
  epsilon <- tryCatch({
    fluc <- stats::glm(Y_sim ~ -1 + H_aw + offset(Q_aw_logit),
                        family = stats::binomial())
    unname(stats::coef(fluc))
  }, error = function(e) 0)

  Q_a1_u <- stats::plogis(
    stats::qlogis(pmax(pmin(Q_a1, 0.999), 0.001)) + epsilon * H_a1)
  Q_a0_u <- stats::plogis(
    stats::qlogis(pmax(pmin(Q_a0, 0.999), 0.001)) + epsilon * H_a0)
  Q_aw_u <- stats::plogis(Q_aw_logit + epsilon * H_aw)

  est <- mean(Q_a1_u) - mean(Q_a0_u)
  eic <- H_aw * (Y_sim - Q_aw_u) + (Q_a1_u - Q_a0_u) - est
  se  <- sqrt(var(eic) / n)

  list(est = est, se = se,
       ci_lower = est - 1.96 * se,
       ci_upper = est + 1.96 * se)
}


#' Fit all candidates for one synthetic replicate, bounded by a wall clock.
#'
#' A degraded synthetic design can send the SuperLearner/glmnet candidate fit
#' into a CPU/memory runaway. That is not an R error (so a plain tryCatch does
#' not help) and in-process elapsed-time limits (setTimeLimit / withTimeout)
#' abort the R front end unsafely on some platforms. We therefore run the
#' replicate's fits inside a persistent, killable `callr` subprocess held in
#' `sess$rs`. If the fits exceed `timeout` seconds the subprocess is killed,
#' a fresh one is started for subsequent replicates, and every candidate for
#' this replicate is recorded as NA -- exactly as a fit error would be. When
#' `sess` is NULL (callr unavailable or no finite timeout) fits run in-process.
#' @keywords internal
#' @noRd
.fit_candidates_bounded <- function(sess, Y_sim, A, W_data, treatment,
                                     covariates, n, tmle_candidates,
                                     timeout, label = "") {
  na1 <- list(est = NA_real_, se = NA_real_,
              ci_lower = NA_real_, ci_upper = NA_real_)
  ids <- vapply(tmle_candidates, function(x) x$candidate_id, character(1))
  na_all <- stats::setNames(rep(list(na1), length(ids)), ids)

  in_process <- function() {
    res <- lapply(tmle_candidates, function(cand)
      tryCatch(.plasmode_fit_one_candidate(
                 Y_sim = Y_sim, A = A, W_data = W_data, treatment = treatment,
                 covariates = covariates, cand = cand, n = n),
               error = function(e) na1))
    stats::setNames(res, ids)
  }

  if (is.null(sess) || is.null(sess$rs) || !is.finite(timeout))
    return(in_process())

  # Self-contained worker: references cleanTMLE by namespace so it does not
  # drag the package environment through serialisation.
  child_fn <- function(Y_sim, A, W_data, treatment, covariates, n, cands) {
    one <- function(cand) tryCatch(
      cleanTMLE:::.plasmode_fit_one_candidate(
        Y_sim = Y_sim, A = A, W_data = W_data, treatment = treatment,
        covariates = covariates, cand = cand, n = n),
      error = function(e) list(est = NA_real_, se = NA_real_,
                               ci_lower = NA_real_, ci_upper = NA_real_))
    stats::setNames(lapply(cands, one),
                    vapply(cands, function(x) x$candidate_id, character(1)))
  }
  environment(child_fn) <- globalenv()

  started <- tryCatch({
    sess$rs$call(child_fn, args = list(Y_sim, A, W_data, treatment,
                                       covariates, n, tmle_candidates))
    TRUE
  }, error = function(e) FALSE)
  if (!started) {                       # session wedged: rebuild and degrade
    try(sess$rs$kill(), silent = TRUE)
    sess$rs <- tryCatch(callr::r_session$new(), error = function(e) NULL)
    return(in_process())
  }

  state <- tryCatch(sess$rs$poll_process(timeout * 1000),
                    error = function(e) "error")
  if (!identical(state, "ready")) {     # timeout or crash: kill the runaway
    message(sprintf(
      "    DQ fit exceeded %.0fs%s -> NA (out-of-process fit killed)",
      timeout, if (nzchar(label)) paste0(" at ", label) else ""))
    try(sess$rs$kill(), silent = TRUE)
    sess$rs <- tryCatch(callr::r_session$new(), error = function(e) NULL)
    return(na_all)
  }
  msg <- tryCatch(sess$rs$read(), error = function(e) NULL)
  if (is.null(msg) || !is.null(msg$error) || is.null(msg$result))
    return(na_all)
  msg$result
}


# ── Internal: DQ Degradation Functions ───────────────────────────────────

#' Introduce MCAR missingness into covariates and impute with column medians
#'
#' Operates on the data.frame in-place (returns a new copy). The outer
#' \code{run_plasmode_dq_stress()} loop sets the seed deterministically
#' before each rep, so this function does not need its own seed argument.
#'
#' @param data data.frame whose covariate columns will be degraded.
#' @param covariates character vector of column names to degrade.
#' @param fraction MCAR fraction in the unit interval.
#' @keywords internal
.degrade_missingness <- function(data, covariates, fraction) {
  for (v in covariates) {
    n <- nrow(data)
    n_miss <- round(n * fraction)
    if (n_miss > 0) {
      idx <- sample.int(n, n_miss)
      data[[v]][idx] <- NA
    }
  }
  # Impute with column medians (matches sanitize_covariates behaviour)
  for (v in covariates) {
    bad <- is.na(data[[v]])
    if (any(bad)) {
      med <- median(data[[v]][!bad], na.rm = TRUE)
      data[[v]][bad] <- if (is.na(med)) 0 else med
    }
  }
  data
}


#' Introduce MAR (treatment-dependent) covariate missingness, then impute
#'
#' Unlike \code{.degrade_missingness} (MCAR), the missingness probability
#' differs by treatment arm, so median imputation is biased rather than
#' merely inefficient. This exercises the missingness *handling* (not just
#' the imputation) and is a stronger test of robustness than MCAR. The
#' control-arm missingness probability is \code{base_fraction}; the treated
#' arm is shifted on the odds scale by \code{treatment_OR}.
#'
#' @param data data.frame whose covariate columns will be degraded.
#' @param covariates character vector of column names to degrade.
#' @param base_fraction control-arm missingness fraction in (0, 1).
#' @param A treatment vector (same length as \code{nrow(data)}).
#' @param treatment_OR odds-ratio multiplier for the treated-arm
#'   missingness probability relative to control. Default 3.
#' @keywords internal
.degrade_missingness_mar <- function(data, covariates, base_fraction, A,
                                     treatment_OR = 3) {
  n  <- nrow(data)
  p0 <- min(max(base_fraction, 1e-6), 1 - 1e-6)
  odds0 <- p0 / (1 - p0)
  p1 <- (odds0 * treatment_OR) / (1 + odds0 * treatment_OR)
  pmiss <- ifelse(A == 1, p1, p0)
  for (v in covariates) {
    miss <- stats::rbinom(n, 1L, pmiss) == 1L
    data[[v]][miss] <- NA
  }
  for (v in covariates) {
    bad <- is.na(data[[v]])
    if (any(bad)) {
      med <- stats::median(data[[v]][!bad], na.rm = TRUE)
      data[[v]][bad] <- if (is.na(med)) 0 else med
    }
  }
  data
}


#' Introduce MNAR covariate missingness (value-dependent), then impute
#'
#' Missingness not at random: the probability that a covariate is missing
#' depends on its own (unobserved) value, so median imputation is biased in a
#' way that cannot be diagnosed from the observed data. Larger standardised
#' values are more likely to be set missing; the column is then median-imputed,
#' which systematically pulls the high tail toward the centre. This is the
#' strongest of the three missingness mechanisms.
#'
#' @param data data.frame whose covariate columns will be degraded.
#' @param covariates character vector of column names to degrade.
#' @param base_fraction baseline missingness fraction in (0, 1) at the column
#'   mean (the logit intercept).
#' @param strength slope (on the logit scale) of the dependence of the
#'   missingness probability on the standardised covariate value. Default 1.5.
#' @keywords internal
.degrade_missingness_mnar <- function(data, covariates, base_fraction,
                                      strength = 1.5) {
  n  <- nrow(data)
  p0 <- min(max(base_fraction, 1e-6), 1 - 1e-6)
  for (v in covariates) {
    x  <- data[[v]]
    sx <- stats::sd(x, na.rm = TRUE)
    z  <- if (is.finite(sx) && sx > 0) (x - mean(x, na.rm = TRUE)) / sx else rep(0, n)
    pmiss <- stats::plogis(stats::qlogis(p0) + strength * z)
    miss  <- stats::rbinom(n, 1L, pmiss) == 1L
    data[[v]][miss] <- NA
  }
  for (v in covariates) {
    bad <- is.na(data[[v]])
    if (any(bad)) {
      med <- stats::median(data[[v]][!bad], na.rm = TRUE)
      data[[v]][bad] <- if (is.na(med)) 0 else med
    }
  }
  data
}


#' Apply asymmetric (sensitivity, specificity) misclassification to A.
#'
#' Interpretation:
#'   sens = P(A* = 1 | A = 1)
#'   spec = P(A* = 0 | A = 0)
#' For backwards compatibility, callers may instead pass a single
#' \code{rate}, in which case the function applies a symmetric
#' non-differential flip (sens = spec = 1 - rate).
#'
#' @param A integer/numeric treatment vector.
#' @param sensitivity numeric in the unit interval or NULL.
#' @param specificity numeric in the unit interval or NULL.
#' @param rate optional symmetric flip rate; used if both
#'   \code{sensitivity} and \code{specificity} are NULL.
#' @keywords internal
.degrade_treatment <- function(A, sensitivity = NULL, specificity = NULL,
                                rate = NULL) {
  n <- length(A)

  if (is.null(sensitivity) && is.null(specificity)) {
    if (is.null(rate)) return(A)
    sensitivity <- 1 - rate
    specificity <- 1 - rate
  } else {
    if (is.null(sensitivity)) sensitivity <- 1
    if (is.null(specificity)) specificity <- 1
  }

  A_out <- A
  pos_idx <- which(A == 1)
  neg_idx <- which(A == 0)

  if (length(pos_idx) > 0L) {
    flip <- stats::rbinom(length(pos_idx), 1L, 1 - sensitivity)
    A_out[pos_idx[flip == 1L]] <- 0L
  }
  if (length(neg_idx) > 0L) {
    flip <- stats::rbinom(length(neg_idx), 1L, 1 - specificity)
    A_out[neg_idx[flip == 1L]] <- 1L
  }
  A_out
}


#' Apply outcome misclassification, optionally differential by treatment arm.
#'
#' If \code{sens_a1}/\code{spec_a1} and \code{sens_a0}/\code{spec_a0} are
#' supplied, sensitivity and specificity may differ between treatment arms,
#' modelling differential outcome ascertainment. Otherwise the marginal
#' (sens, spec) pair is applied to both arms.
#'
#' @param Y integer outcome vector.
#' @param sensitivity marginal sensitivity (used if per-arm values are NULL).
#' @param specificity marginal specificity.
#' @param A optional treatment vector for per-arm differential
#'   misclassification.
#' @param sens_a1,sens_a0,spec_a1,spec_a0 optional per-arm overrides.
#' @keywords internal
.degrade_outcome <- function(Y, sensitivity = NULL, specificity = NULL,
                              A = NULL,
                              sens_a1 = NULL, spec_a1 = NULL,
                              sens_a0 = NULL, spec_a0 = NULL) {
  Y_out <- Y

  has_per_arm <- !is.null(A) &&
    (!is.null(sens_a1) || !is.null(spec_a1) ||
     !is.null(sens_a0) || !is.null(spec_a0))

  if (!has_per_arm) {
    sens <- if (is.null(sensitivity)) 1 else sensitivity
    spec <- if (is.null(specificity)) 1 else specificity
    pos_idx <- which(Y == 1)
    if (length(pos_idx) > 0L) {
      flip <- stats::rbinom(length(pos_idx), 1L, 1 - sens)
      Y_out[pos_idx[flip == 1L]] <- 0L
    }
    neg_idx <- which(Y == 0)
    if (length(neg_idx) > 0L) {
      flip <- stats::rbinom(length(neg_idx), 1L, 1 - spec)
      Y_out[neg_idx[flip == 1L]] <- 1L
    }
    return(Y_out)
  }

  # Per-arm differential misclassification.
  sens_a1 <- if (is.null(sens_a1)) (sensitivity %||% 1) else sens_a1
  sens_a0 <- if (is.null(sens_a0)) (sensitivity %||% 1) else sens_a0
  spec_a1 <- if (is.null(spec_a1)) (specificity %||% 1) else spec_a1
  spec_a0 <- if (is.null(spec_a0)) (specificity %||% 1) else spec_a0

  for (a in c(0L, 1L)) {
    arm_idx <- which(A == a)
    if (length(arm_idx) == 0L) next
    sens_a <- if (a == 1L) sens_a1 else sens_a0
    spec_a <- if (a == 1L) spec_a1 else spec_a0

    pos_idx <- arm_idx[Y[arm_idx] == 1L]
    if (length(pos_idx) > 0L) {
      flip <- stats::rbinom(length(pos_idx), 1L, 1 - sens_a)
      Y_out[pos_idx[flip == 1L]] <- 0L
    }
    neg_idx <- arm_idx[Y[arm_idx] == 0L]
    if (length(neg_idx) > 0L) {
      flip <- stats::rbinom(length(neg_idx), 1L, 1 - spec_a)
      Y_out[neg_idx[flip == 1L]] <- 1L
    }
  }

  Y_out
}


# Local null-coalescing operator.
`%||%` <- function(a, b) if (is.null(a)) b else a


# ── Main: DQ Stress Test ────────────────────────────────────────────────

#' Run Data-Quality Stress Tests via Plasmode Simulation
#'
#' Evaluates TMLE candidate robustness under controlled data-quality
#' degradations. For each DQ scenario and severity level, plasmode
#' outcomes are generated from the real covariate distribution and
#' the fitted Q0 model, then degraded before each prespecified TMLE
#' candidate is fit. The procedure is outcome-blind: only synthetic
#' outcomes are used.
#'
#' @section Clean-room stage: Stage 2b (pre-outcome).
#'
#' @param lock A \code{cleanroom_lock}.
#' @param tmle_candidates A list of \code{tmle_candidate_spec} objects.
#' @param effect_sizes Numeric vector of true risk differences.
#'   Default: \code{c(0.05)}.
#' @param reps Integer; plasmode replicates per scenario. Default 50,
#'   typically smaller than the baseline plasmode for speed.
#' @param data_quality_scenarios A named list of DQ scenarios; see Details.
#' @param q0_library Optional SuperLearner library (character vector) for
#'   the plasmode outcome generator Q0. Default \code{NULL} uses a
#'   covariate-only logistic GLM. Supplying a richer library reduces the
#'   linear-in-logit bias the GLM Q0 imposes on candidate selection when
#'   the true outcome surface is nonlinear.
#' @param verbose Logical; print progress messages.
#'
#' @details
#' The \code{data_quality_scenarios} list may contain any combination of:
#'
#' \describe{
#'   \item{\code{covariate_missingness}}{MCAR. List with \code{fractions}
#'     (numeric vector, e.g. \code{c(0.05, 0.10, 0.20)}) and optional
#'     \code{variables} (character vector; default: all covariates).}
#'   \item{\code{covariate_missingness_mar}}{MAR (treatment-dependent).
#'     List with \code{fractions} (control-arm missingness fractions),
#'     optional \code{variables}, and optional \code{treatment_OR}
#'     (odds-ratio multiplier for treated-arm missingness; default 3).
#'     Because missingness depends on the arm and the column is
#'     median-imputed, this scenario is a stronger test than MCAR: the
#'     imputation is biased rather than merely inefficient.}
#'   \item{\code{covariate_missingness_mnar}}{MNAR (value-dependent).
#'     List with \code{fractions} (baseline missingness fractions at the
#'     column mean), optional \code{variables}, and optional \code{strength}
#'     (logit-scale slope of the dependence of the missingness probability on
#'     the standardised covariate value; default 1.5). The probability that a
#'     covariate is missing depends on its own unobserved value, so median
#'     imputation is biased in a way that cannot be diagnosed from the observed
#'     data. This is the strongest of the three missingness mechanisms.}
#'   \item{\code{treatment_misclass}}{List with optional \code{rates}
#'     (numeric vector for symmetric flips) and / or asymmetric
#'     \code{sensitivity} and \code{specificity} (numeric vectors of
#'     equal length).}
#'   \item{\code{outcome_misclass}}{List with \code{sensitivity} and
#'     \code{specificity} (numeric vectors of equal length). For
#'     differential misclassification by treatment arm, supply
#'     \code{sens_a1}, \code{sens_a0}, \code{spec_a1}, \code{spec_a0}
#'     of equal length instead of (or alongside) the marginal pair.}
#'   \item{\code{unmeasured_confounding}}{List with \code{U_prevalence},
#'     \code{U_treatment_OR}, and \code{U_outcome_OR} (the latter two
#'     of equal length).}
#'   \item{\code{near_positivity}}{List with \code{slopes} (numeric vector of
#'     multipliers greater than 1). Each multiplier scales the centred
#'     log-odds of the lock-data propensity model, amplifying the
#'     covariate-to-treatment association so a subgroup approaches
#'     deterministic treatment and the estimated propensity score reaches the
#'     boundary. This stresses positivity (the inverse-probability weight tail)
#'     rather than the outcome surface, and it is the threat under which a
#'     lightly truncated candidate degrades while a heavily truncated candidate
#'     is protected.}
#' }
#'
#' @return An object of class \code{plasmode_dq_results} containing
#'   \describe{
#'     \item{metrics}{data.frame with columns scenario, level,
#'       effect_size, candidate, bias, rmse, coverage, emp_sd,
#'       mean_se, se_cal, n_converged.}
#'     \item{scenarios}{the input \code{data_quality_scenarios}.}
#'     \item{baseline}{rows of metrics where scenario == "none".}
#'   }
#'
#' @export
run_plasmode_dq_stress <- function(lock,
                                    tmle_candidates = NULL,
                                    effect_sizes    = c(0.05),
                                    reps            = 50L,
                                    data_quality_scenarios = list(),
                                    q0_library      = NULL,
                                    fit_timeout     = Inf,
                                    verbose         = TRUE) {
  if (!inherits(lock, "cleanroom_lock"))
    stop("`lock` must be a cleanroom_lock object.", call. = FALSE)

  if (is.null(tmle_candidates))
    tmle_candidates <- expand_tmle_candidate_grid()
  validate_tmle_candidates(tmle_candidates)

  data       <- lock$data
  treatment  <- lock$treatment
  outcome    <- lock$outcome
  covariates <- lock$covariates
  A          <- data[[treatment]]
  n          <- nrow(data)

  # Guard against locks where the outcome column is fully NA (e.g. a lock
  # that has been mask_outcome()'d). The Q0 fit needs the real outcome.
  y_all <- data[[outcome]]
  n_obs_y <- sum(!is.na(y_all))
  if (n_obs_y == 0L) {
    stop("Q0 model cannot be fit: lock$data[[lock$outcome]] has zero ",
         "non-NA observations. If the outcome is masked, call ",
         "unmask_outcome() before run_plasmode_dq_stress(); plasmode ",
         "uses the marginal Y|W distribution but does not look at the ",
         "treatment-outcome association.", call. = FALSE)
  }
  if (n_obs_y < length(covariates) + 1L) {
    stop("Q0 model cannot be fit: only ", n_obs_y, " non-NA outcome ",
         "rows available for ", length(covariates), " covariates. ",
         "Drop covariates or check the cohort filter.", call. = FALSE)
  }

  # Baseline outcome model for plasmode DGP (Q0; outcome-blind in the sense
  # that the treatment-outcome association is not used here -- this is the
  # covariate-only mean of Y). Fit on complete cases and predict for ALL
  # rows so the resulting probability vector has length n. Without
  # newdata = data, predict.glm() returns predictions only for the rows
  # used in fitting (which fails downstream when Y has NAs).
  # Default Q0 is a covariate-only logistic GLM. When `q0_library` is
  # supplied and SuperLearner is available, fit Q0 with that library
  # instead; this relaxes the linear-in-logit bias that a GLM Q0 imposes
  # on candidate selection when the true outcome surface is nonlinear.
  use_sl_q0 <- !is.null(q0_library) &&
    requireNamespace("SuperLearner", quietly = TRUE)
  if (use_sl_q0) {
    cc <- !is.na(data[[outcome]])
    Q0_sl <- SuperLearner::SuperLearner(
      Y = data[[outcome]][cc],
      X = data[cc, covariates, drop = FALSE],
      family = stats::binomial(), SL.library = q0_library,
      env = .cleantmle_sl_env())
    p_base <- as.numeric(predict(
      Q0_sl, newdata = data[, covariates, drop = FALSE])$pred)
  } else {
    if (!is.null(q0_library))
      warning("q0_library supplied but SuperLearner is not available; ",
              "falling back to logistic-GLM Q0.", call. = FALSE)
    Q0_fml <- stats::reformulate(covariates, response = outcome)
    Q0_fit <- stats::glm(Q0_fml, data = data, family = stats::binomial(),
                         na.action = stats::na.exclude)
    p_base <- as.numeric(stats::predict(Q0_fit, newdata = data,
                                         type = "response"))
  }
  if (length(p_base) != n) {
    stop("Internal error in plasmode Q0 fit: predicted vector length (",
         length(p_base), ") does not match data rows (", n, ").",
         call. = FALSE)
  }
  if (any(is.na(p_base))) p_base[is.na(p_base)] <- mean(p_base, na.rm = TRUE)
  p_base <- pmin(pmax(p_base, 0.001), 0.999)

  # Real propensity (treatment-mechanism) model fitted on the lock data.
  # Used as the *baseline* propensity that the unmeasured-confounding
  # scenario shifts. Earlier prototype versions of this function reused
  # the outcome model for this purpose, which biased the U scenario.
  ps_fml  <- stats::reformulate(covariates, response = treatment)
  ps_mod  <- stats::glm(ps_fml, data = data, family = stats::binomial())
  ps_base <- as.numeric(stats::predict(ps_mod, type = "response"))
  ps_base <- pmin(pmax(ps_base, 0.001), 0.999)

  cand_ids <- vapply(tmle_candidates, function(x) x$candidate_id, character(1))

  # ── Build scenario grid ──────────────────────────────────────────────
  scenario_grid <- data.frame(
    scenario = "none", level = "0", stringsAsFactors = FALSE
  )

  dqs <- data_quality_scenarios

  if (!is.null(dqs$covariate_missingness)) {
    for (f in dqs$covariate_missingness$fractions) {
      scenario_grid <- rbind(scenario_grid, data.frame(
        scenario = "cov_miss", level = as.character(f),
        stringsAsFactors = FALSE))
    }
  }

  if (!is.null(dqs$covariate_missingness_mar)) {
    for (f in dqs$covariate_missingness_mar$fractions) {
      scenario_grid <- rbind(scenario_grid, data.frame(
        scenario = "cov_miss_mar", level = as.character(f),
        stringsAsFactors = FALSE))
    }
  }

  if (!is.null(dqs$covariate_missingness_mnar)) {
    for (f in dqs$covariate_missingness_mnar$fractions) {
      scenario_grid <- rbind(scenario_grid, data.frame(
        scenario = "cov_miss_mnar", level = as.character(f),
        stringsAsFactors = FALSE))
    }
  }

  if (!is.null(dqs$treatment_misclass)) {
    sens_t <- dqs$treatment_misclass$sensitivity
    spec_t <- dqs$treatment_misclass$specificity
    rates_t <- dqs$treatment_misclass$rates
    if (!is.null(sens_t) || !is.null(spec_t)) {
      ns <- max(length(sens_t), length(spec_t))
      if (is.null(sens_t)) sens_t <- rep(1, ns)
      if (is.null(spec_t)) spec_t <- rep(1, ns)
      for (i in seq_len(ns)) {
        scenario_grid <- rbind(scenario_grid, data.frame(
          scenario = "trt_misclass",
          level = sprintf("sens%.2f_spec%.2f", sens_t[i], spec_t[i]),
          stringsAsFactors = FALSE))
      }
    }
    if (!is.null(rates_t)) {
      for (r in rates_t) {
        scenario_grid <- rbind(scenario_grid, data.frame(
          scenario = "trt_misclass",
          level = sprintf("rate%.2f", r),
          stringsAsFactors = FALSE))
      }
    }
  }

  if (!is.null(dqs$outcome_misclass)) {
    sens <- dqs$outcome_misclass$sensitivity
    spec <- dqs$outcome_misclass$specificity
    sens_a1 <- dqs$outcome_misclass$sens_a1
    sens_a0 <- dqs$outcome_misclass$sens_a0
    spec_a1 <- dqs$outcome_misclass$spec_a1
    spec_a0 <- dqs$outcome_misclass$spec_a0
    has_marginal <- !is.null(sens) || !is.null(spec)
    has_perarm   <- !is.null(sens_a1) || !is.null(spec_a1) ||
                    !is.null(sens_a0) || !is.null(spec_a0)
    if (has_marginal) {
      for (i in seq_along(sens)) {
        scenario_grid <- rbind(scenario_grid, data.frame(
          scenario = "out_misclass",
          level = sprintf("sens%.2f_spec%.2f", sens[i], spec[i]),
          stringsAsFactors = FALSE))
      }
    }
    if (has_perarm) {
      ns <- max(length(sens_a1) %||% 0L, length(sens_a0) %||% 0L,
                length(spec_a1) %||% 0L, length(spec_a0) %||% 0L)
      sens_a1 <- if (is.null(sens_a1)) rep(1, ns) else sens_a1
      sens_a0 <- if (is.null(sens_a0)) rep(1, ns) else sens_a0
      spec_a1 <- if (is.null(spec_a1)) rep(1, ns) else spec_a1
      spec_a0 <- if (is.null(spec_a0)) rep(1, ns) else spec_a0
      for (i in seq_len(ns)) {
        scenario_grid <- rbind(scenario_grid, data.frame(
          scenario = "out_misclass_diff",
          level = sprintf("a1[%.2f,%.2f]_a0[%.2f,%.2f]",
                          sens_a1[i], spec_a1[i],
                          sens_a0[i], spec_a0[i]),
          stringsAsFactors = FALSE))
      }
    }
  }

  if (!is.null(dqs$unmeasured_confounding)) {
    u_trt <- dqs$unmeasured_confounding$U_treatment_OR
    u_out <- dqs$unmeasured_confounding$U_outcome_OR
    for (i in seq_along(u_trt)) {
      scenario_grid <- rbind(scenario_grid, data.frame(
        scenario = "unmeasured_U",
        level = sprintf("OR_trt%.1f_out%.1f", u_trt[i], u_out[i]),
        stringsAsFactors = FALSE))
    }
  }

  # Practical positivity violation: the covariate-to-treatment association is
  # amplified so a subgroup approaches deterministic treatment and the estimated
  # propensity score reaches the boundary. The `slopes` are multipliers (> 1) on
  # the centred log-odds of the lock-data propensity model.
  if (!is.null(dqs$near_positivity)) {
    for (s in dqs$near_positivity$slopes) {
      scenario_grid <- rbind(scenario_grid, data.frame(
        scenario = "near_positivity",
        level = sprintf("slope_x%.1f", s),
        stringsAsFactors = FALSE))
    }
  }

  if (verbose) {
    cat(sprintf("DQ Stress Test: %d scenarios x %d effect sizes x %d reps x %d candidates\n",
                nrow(scenario_grid), length(effect_sizes), reps, length(cand_ids)))
  }

  # ── Out-of-process fit guard ──────────────────────────────────────────
  # Each replicate's candidate fits run in a persistent, killable callr
  # subprocess so a degenerate synthetic design that sends glmnet into a
  # runaway cannot wedge the whole study. See .fit_candidates_bounded().
  # Falls back to in-process fitting when callr is unavailable.
  sess <- NULL
  if (is.finite(fit_timeout) && requireNamespace("callr", quietly = TRUE)) {
    sess <- new.env(parent = emptyenv())
    sess$rs <- tryCatch(callr::r_session$new(), error = function(e) NULL)
    if (is.null(sess$rs)) sess <- NULL
  }
  on.exit({
    if (!is.null(sess) && !is.null(sess$rs)) try(sess$rs$close(), silent = TRUE)
  }, add = TRUE)

  # ── Run plasmode across scenario grid ─────────────────────────────────
  all_metrics <- list()

  for (sg_i in seq_len(nrow(scenario_grid))) {
    sc_name <- scenario_grid$scenario[sg_i]
    sc_level <- scenario_grid$level[sg_i]

    if (verbose) cat(sprintf("  Scenario: %s [%s]\n", sc_name, sc_level))

    # Resolve scenario-specific parameters from the level string.
    sens_lvl <- spec_lvl <- rate_lvl <- NA_real_
    sens_a1_lvl <- sens_a0_lvl <- spec_a1_lvl <- spec_a0_lvl <- NA_real_
    u_trt_or_lvl <- u_out_or_lvl <- NA_real_

    if (sc_name == "trt_misclass") {
      if (grepl("^rate", sc_level)) {
        rate_lvl <- as.numeric(sub("rate", "", sc_level))
      } else {
        parts <- strsplit(sc_level, "_")[[1]]
        sens_lvl <- as.numeric(sub("sens", "", parts[1]))
        spec_lvl <- as.numeric(sub("spec", "", parts[2]))
      }
    }
    if (sc_name == "out_misclass") {
      parts <- strsplit(sc_level, "_")[[1]]
      sens_lvl <- as.numeric(sub("sens", "", parts[1]))
      spec_lvl <- as.numeric(sub("spec", "", parts[2]))
    }
    if (sc_name == "out_misclass_diff") {
      m <- regmatches(sc_level,
        regexec("a1\\[([0-9.]+),([0-9.]+)\\]_a0\\[([0-9.]+),([0-9.]+)\\]",
                sc_level))[[1]]
      sens_a1_lvl <- as.numeric(m[2]); spec_a1_lvl <- as.numeric(m[3])
      sens_a0_lvl <- as.numeric(m[4]); spec_a0_lvl <- as.numeric(m[5])
    }
    if (sc_name == "unmeasured_U") {
      m <- regmatches(sc_level,
        regexec("OR_trt([0-9.]+)_out([0-9.]+)", sc_level))[[1]]
      u_trt_or_lvl <- as.numeric(m[2])
      u_out_or_lvl <- as.numeric(m[3])
    }
    pos_slope_lvl <- NA_real_
    if (sc_name == "near_positivity") {
      pos_slope_lvl <- as.numeric(sub("slope_x", "", sc_level))
    }

    for (es in effect_sizes) {
      rep_results <- vector("list", reps)

      for (rep_i in seq_len(reps)) {
        set.seed(lock$seed + rep_i + sg_i * 10000L)

        # Generate synthetic outcome and possibly modify A under U.
        # Clamp BOTH bounds so negative `es` with small p_base does
        # not produce negative probabilities (rbinom -> NA).
        A_rep <- A
        p1_sim <- pmin(pmax(p_base + es, 0.001), 0.999)
        p0_sim <- p_base

        if (sc_name == "unmeasured_U") {
          u_prev <- dqs$unmeasured_confounding$U_prevalence
          U <- stats::rbinom(n, 1, u_prev)

          # Use the *real* propensity model fitted on lock data as the
          # baseline; shift its log-odds by log(OR_A) * U.
          lp_ps   <- stats::qlogis(ps_base)
          lp_ps_u <- lp_ps + log(u_trt_or_lvl) * U
          ps_u    <- stats::plogis(lp_ps_u)
          A_rep   <- stats::rbinom(n, 1L, ps_u)

          # Modify outcome probabilities given U.
          lp1 <- stats::qlogis(pmax(pmin(p1_sim, 0.999), 0.001)) +
            log(u_out_or_lvl) * U
          lp0 <- stats::qlogis(pmax(pmin(p0_sim, 0.999), 0.001)) +
            log(u_out_or_lvl) * U
          p1_sim <- stats::plogis(lp1)
          p0_sim <- stats::plogis(lp0)
        }

        if (sc_name == "near_positivity") {
          # Amplify the covariate-to-treatment slopes about their mean so a
          # subgroup approaches deterministic treatment (PS -> 0/1). The outcome
          # model is left unchanged; this threat stresses positivity, not the
          # outcome surface.
          lp_ps  <- stats::qlogis(ps_base)
          lp_bar <- mean(lp_ps)
          ps_pos <- stats::plogis(lp_bar + pos_slope_lvl * (lp_ps - lp_bar))
          A_rep  <- stats::rbinom(n, 1L, ps_pos)
        }

        p_obs  <- ifelse(A_rep == 1, p1_sim, p0_sim)
        Y_sim  <- stats::rbinom(n, 1L, p_obs)
        truth  <- mean(p1_sim) - mean(p0_sim)

        # Treatment misclassification (asymmetric or symmetric).
        A_fit <- A_rep
        if (sc_name == "trt_misclass") {
          if (!is.na(rate_lvl)) {
            A_fit <- .degrade_treatment(A_rep, rate = rate_lvl)
          } else {
            A_fit <- .degrade_treatment(A_rep,
                                         sensitivity = sens_lvl,
                                         specificity = spec_lvl)
          }
        }

        # Outcome misclassification (marginal or per-arm).
        Y_fit <- Y_sim
        if (sc_name == "out_misclass") {
          Y_fit <- .degrade_outcome(Y_sim,
                                     sensitivity = sens_lvl,
                                     specificity = spec_lvl)
        }
        if (sc_name == "out_misclass_diff") {
          Y_fit <- .degrade_outcome(Y_sim,
                                     A = A_rep,
                                     sens_a1 = sens_a1_lvl,
                                     sens_a0 = sens_a0_lvl,
                                     spec_a1 = spec_a1_lvl,
                                     spec_a0 = spec_a0_lvl)
        }

        # Covariate missingness (MCAR).
        W_fit <- data
        if (sc_name == "cov_miss") {
          frac <- as.numeric(sc_level)
          miss_vars <- dqs$covariate_missingness$variables
          if (is.null(miss_vars)) miss_vars <- covariates
          W_fit <- .degrade_missingness(data, miss_vars, frac)
        }
        # Covariate missingness (MAR; treatment-dependent + median impute).
        if (sc_name == "cov_miss_mar") {
          frac <- as.numeric(sc_level)
          spec <- dqs$covariate_missingness_mar
          miss_vars <- spec$variables
          if (is.null(miss_vars)) miss_vars <- covariates
          or_a <- spec$treatment_OR %||% 3
          W_fit <- .degrade_missingness_mar(data, miss_vars, frac,
                                            A = A_rep, treatment_OR = or_a)
        }
        # Covariate missingness (MNAR; value-dependent + median impute).
        if (sc_name == "cov_miss_mnar") {
          frac <- as.numeric(sc_level)
          spec <- dqs$covariate_missingness_mnar
          miss_vars <- spec$variables
          if (is.null(miss_vars)) miss_vars <- covariates
          strength <- spec$strength %||% 1.5
          W_fit <- .degrade_missingness_mnar(data, miss_vars, frac,
                                             strength = strength)
        }

        # Fit all candidates for this replicate inside a killable subprocess
        # bounded by `fit_timeout` seconds (see .fit_candidates_bounded): a
        # degraded synthetic design can send the SuperLearner/glmnet fit into
        # a CPU/memory runaway that is not an R error, so one pathological draw
        # would otherwise wedge the entire stress test. On timeout the worker
        # is killed and this replicate's candidates are recorded as NA.
        cand_results <- .fit_candidates_bounded(
          sess, Y_sim = Y_fit, A = A_fit, W_data = W_fit,
          treatment = treatment, covariates = covariates, n = n,
          tmle_candidates = tmle_candidates, timeout = fit_timeout,
          label = sprintf("%s[%s] rep %d", sc_name, sc_level, rep_i))

        rep_results[[rep_i]] <- c(cand_results, list(.truth = truth))
      }

      # Aggregate per-candidate metrics.
      truth_v <- vapply(rep_results, function(x) x$.truth, numeric(1))

      for (cid in cand_ids) {
        ests   <- vapply(rep_results, function(x) x[[cid]]$est, numeric(1))
        ses    <- vapply(rep_results, function(x) x[[cid]]$se, numeric(1))
        ci_los <- vapply(rep_results, function(x) x[[cid]]$ci_lower, numeric(1))
        ci_his <- vapply(rep_results, function(x) x[[cid]]$ci_upper, numeric(1))

        valid <- !is.na(ests)
        if (sum(valid) == 0) next

        ests_v   <- ests[valid]
        ses_v    <- ses[valid]
        truth_vv <- truth_v[valid]

        emp_sd  <- sd(ests_v)
        mean_se <- mean(ses_v)

        row <- data.frame(
          scenario    = sc_name,
          level       = sc_level,
          effect_size = es,
          candidate   = cid,
          bias        = round(mean(ests_v - truth_vv), 5),
          rmse        = round(sqrt(mean((ests_v - truth_vv)^2)), 5),
          coverage    = round(mean(ci_los[valid] <= truth_vv &
                                     truth_vv <= ci_his[valid]), 3),
          emp_sd      = round(emp_sd, 5),
          mean_se     = round(mean_se, 5),
          se_cal      = round(if (emp_sd > 0) mean_se / emp_sd else NA, 3),
          n_converged = sum(valid),
          stringsAsFactors = FALSE
        )
        all_metrics <- c(all_metrics, list(row))
      }
    }
  }

  metrics <- do.call(rbind, all_metrics)

  result <- list(
    metrics   = metrics,
    scenarios = data_quality_scenarios,
    baseline  = metrics[metrics$scenario == "none", ],
    lock      = lock,
    reps      = reps,
    call      = match.call()
  )
  class(result) <- "plasmode_dq_results"
  result
}


#' @export
print.plasmode_dq_results <- function(x, ...) {
  cat("Plasmode Data-Quality Stress Test\n")
  cat("==================================\n")
  cat(sprintf("Replicates per scenario: %d\n", x$reps))
  cat(sprintf("Total scenario-candidate rows: %d\n\n", nrow(x$metrics)))

  scenarios <- unique(x$metrics$scenario)
  for (sc in scenarios) {
    sub <- x$metrics[x$metrics$scenario == sc, ]
    cat(sprintf("--- %s ---\n", sc))
    print(sub[, c("level", "candidate", "bias", "rmse", "coverage", "se_cal")],
          row.names = FALSE)
    cat("\n")
  }
  invisible(x)
}


#' @export
plot.plasmode_dq_results <- function(x, metric = c("rmse", "bias", "coverage"),
                                     ...) {
  metric <- match.arg(metric)
  if (!requireNamespace("ggplot2", quietly = TRUE))
    stop("ggplot2 is required for plot.plasmode_dq_results().", call. = FALSE)

  m <- x$metrics
  m <- m[m$scenario != "none", , drop = FALSE]
  if (nrow(m) == 0L) {
    message("No degraded-scenario rows to plot.")
    return(invisible(NULL))
  }

  m$y <- switch(metric,
                rmse     = m$rmse,
                bias     = abs(m$bias),
                coverage = m$coverage)

  ylab_text <- switch(metric,
                      rmse     = "RMSE",
                      bias     = "|Bias|",
                      coverage = "Coverage")

  p <- ggplot2::ggplot(m, ggplot2::aes_string(x = "level", y = "y",
                                                colour = "candidate",
                                                group = "candidate")) +
    ggplot2::geom_point(size = 2.2) +
    ggplot2::geom_line(linewidth = 0.7) +
    ggplot2::facet_wrap(~ scenario, scales = "free_x") +
    ggplot2::labs(x = "DQ severity level", y = ylab_text,
                  title = sprintf("Degradation gradient (%s)", metric),
                  colour = "Candidate") +
    ggplot2::theme_minimal() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 30, hjust = 1,
                                                         size = 8))

  if (metric == "coverage") {
    p <- p + ggplot2::geom_hline(yintercept = 0.95, linetype = "dashed",
                                  colour = "red")
  }

  p
}


#' Pre-Outcome Checkpoint Derived from a DQ Stress Test
#'
#' Converts the degraded-scenario rows of a \code{plasmode_dq_results}
#' object into a \code{cleantmle_checkpoint} so the DQ output can be
#' fed into \code{\link{gate_all}()} alongside the cohort-adequacy,
#' balance, and residual-bias checkpoints. The checkpoint flips to
#' STOP if any degraded row for the supplied candidate exceeds the
#' configured \code{max_abs_bias}, falls below \code{min_coverage},
#' or has an RMSE ratio above \code{max_rmse_ratio}. FLAG triggers
#' when the bias is within tolerance but coverage or RMSE-ratio
#' breach a softer envelope (\code{flag_coverage}, \code{flag_rmse_ratio}).
#'
#' @param dq_results A \code{plasmode_dq_results} object.
#' @param candidate Character; the \code{candidate_id} that the gate
#'   should evaluate. Defaults to the first candidate appearing in
#'   the degraded rows.
#' @param max_abs_bias Numeric; any degraded row with
#'   \code{|bias| > max_abs_bias} returns STOP. Default 0.02.
#' @param min_coverage Numeric; any degraded row with
#'   \code{coverage < min_coverage} returns STOP. Default 0.85.
#' @param max_rmse_ratio Numeric; any degraded row with
#'   \code{rmse / rmse_baseline > max_rmse_ratio} returns STOP.
#'   Default 1.5.
#' @param flag_coverage,flag_rmse_ratio Softer thresholds that demote
#'   a GO to a FLAG when no STOP condition is met. Defaults 0.90 and
#'   1.20.
#' @param scenarios Optional character vector restricting which DQ
#'   scenarios to evaluate (e.g. \code{c("unmeasured_U")}). Default
#'   evaluates every degraded row.
#' @param thresholds Optional \code{\link{decision_thresholds}} object (or
#'   its \code{$dq} sublist). When supplied, the DQ thresholds are drawn
#'   from it, overriding the \code{max_abs_bias} / \code{min_coverage} /
#'   \code{max_rmse_ratio} / \code{flag_*} arguments, so the gate reads
#'   from a single prespecified, fingerprinted source.
#' @param lock_hash Optional lock hash to attach to the checkpoint.
#'
#' @return A \code{cleantmle_checkpoint} of stage
#'   \code{"Check Point 2c: DQ Stress"}.
#'
#' @examples
#' \dontrun{
#' cp_dq <- gate_dq(dq_results, candidate = best$candidate_id)
#' audit <- record_checkpoint(audit, cp_dq)
#' gate  <- gate_all(cp1, cp2, cp_dq, cp3, allow_flag = TRUE)
#' }
#' @export
gate_dq <- function(dq_results,
                    candidate      = NULL,
                    max_abs_bias   = 0.02,
                    min_coverage   = 0.85,
                    max_rmse_ratio = 1.5,
                    flag_coverage  = 0.90,
                    flag_rmse_ratio = 1.20,
                    scenarios      = NULL,
                    thresholds     = NULL,
                    lock_hash      = NA_character_) {
  if (!inherits(dq_results, "plasmode_dq_results"))
    stop("`dq_results` must be a plasmode_dq_results object.", call. = FALSE)

  # If a prespecified decision_thresholds() object (or its $dq sublist) is
  # supplied, draw the DQ thresholds from it so the gate reads from a
  # single fingerprinted source rather than ad hoc arguments.
  if (!is.null(thresholds)) {
    dq <- if (inherits(thresholds, "cleantmle_thresholds")) thresholds$dq
          else thresholds
    if (!is.null(dq$max_abs_bias))    max_abs_bias    <- dq$max_abs_bias
    if (!is.null(dq$min_coverage))    min_coverage    <- dq$min_coverage
    if (!is.null(dq$max_rmse_ratio))  max_rmse_ratio  <- dq$max_rmse_ratio
    if (!is.null(dq$flag_coverage))   flag_coverage   <- dq$flag_coverage
    if (!is.null(dq$flag_rmse_ratio)) flag_rmse_ratio <- dq$flag_rmse_ratio
  }

  deg <- summarize_dq_degradation(dq_results)
  if (nrow(deg) == 0L) {
    return(new_checkpoint(
      stage      = "Check Point 2c: DQ Stress",
      decision   = "FLAG",
      metrics    = data.frame(note = "no degraded rows to evaluate"),
      thresholds = list(max_abs_bias = max_abs_bias,
                        min_coverage = min_coverage,
                        max_rmse_ratio = max_rmse_ratio),
      rationale  = "DQ stress test produced no degraded-scenario rows.",
      lock_hash  = lock_hash
    ))
  }

  if (is.null(candidate)) candidate <- deg$candidate[1]
  if (!candidate %in% deg$candidate)
    stop("candidate '", candidate, "' not found in DQ output.", call. = FALSE)

  sub <- deg[deg$candidate == candidate, , drop = FALSE]
  if (!is.null(scenarios))
    sub <- sub[sub$scenario %in% scenarios, , drop = FALSE]
  if (nrow(sub) == 0L)
    stop("no DQ rows match the requested scenarios.", call. = FALSE)

  bias_v <- abs(sub$bias_degraded)
  cov_v  <- sub$cov_degraded
  rr_v   <- sub$rmse_ratio

  stop_bias <- any(bias_v > max_abs_bias, na.rm = TRUE)
  stop_cov  <- any(cov_v  < min_coverage,  na.rm = TRUE)
  stop_rr   <- any(rr_v   > max_rmse_ratio, na.rm = TRUE)

  flag_cov  <- any(cov_v < flag_coverage,    na.rm = TRUE)
  flag_rr   <- any(rr_v  > flag_rmse_ratio,  na.rm = TRUE)

  decision <- if (stop_bias || stop_cov || stop_rr) "STOP"
              else if (flag_cov || flag_rr)         "FLAG"
              else                                   "GO"

  worst <- sub[which.max(bias_v), , drop = FALSE]
  metrics <- data.frame(
    candidate     = candidate,
    n_rows        = nrow(sub),
    max_abs_bias  = round(max(bias_v, na.rm = TRUE), 5),
    min_coverage  = round(min(cov_v,  na.rm = TRUE), 3),
    max_rmse_ratio = round(max(rr_v,  na.rm = TRUE), 3),
    worst_scenario = worst$scenario,
    worst_level    = worst$level,
    stringsAsFactors = FALSE
  )

  rationale <- switch(decision,
    STOP = sprintf("DQ stress exceeds locked thresholds (worst: %s @ %s).",
                   worst$scenario, worst$level),
    FLAG = "DQ stress within STOP thresholds but past FLAG envelope.",
    GO   = "DQ stress within locked thresholds across scenarios.")

  new_checkpoint(
    stage      = "Check Point 2c: DQ Stress",
    decision   = decision,
    metrics    = metrics,
    thresholds = list(max_abs_bias   = max_abs_bias,
                      min_coverage   = min_coverage,
                      max_rmse_ratio = max_rmse_ratio,
                      flag_coverage  = flag_coverage,
                      flag_rmse_ratio = flag_rmse_ratio,
                      scenarios      = scenarios),
    rationale  = rationale,
    lock_hash  = lock_hash
  )
}


#' Summarise DQ Stress Test as a Degradation Table
#'
#' Computes the relative change in bias, RMSE, and coverage versus the
#' baseline (no degradation) for each DQ scenario, level, and candidate.
#'
#' @param dq_results A \code{plasmode_dq_results} object.
#'
#' @return A data.frame with relative degradation metrics.
#'
#' @export
summarize_dq_degradation <- function(dq_results) {
  if (!inherits(dq_results, "plasmode_dq_results"))
    stop("`dq_results` must be a plasmode_dq_results object.", call. = FALSE)

  m <- dq_results$metrics
  base <- m[m$scenario == "none", ]
  degraded <- m[m$scenario != "none", ]

  if (nrow(degraded) == 0) return(data.frame())

  rows <- lapply(seq_len(nrow(degraded)), function(i) {
    d <- degraded[i, ]
    b <- base[base$candidate == d$candidate &
                base$effect_size == d$effect_size, ]
    if (nrow(b) == 0) return(NULL)

    data.frame(
      scenario      = d$scenario,
      level         = d$level,
      candidate     = d$candidate,
      bias_baseline = b$bias,
      bias_degraded = d$bias,
      bias_change   = round(abs(d$bias) - abs(b$bias), 5),
      rmse_baseline = b$rmse,
      rmse_degraded = d$rmse,
      rmse_ratio    = round(d$rmse / max(b$rmse, 1e-6), 3),
      cov_baseline  = b$coverage,
      cov_degraded  = d$coverage,
      cov_drop      = round(b$coverage - d$coverage, 3),
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, Filter(Negate(is.null), rows))
}


#' Assess Synthetic-vs-Real Data Fidelity for a Plasmode Generator
#'
#' Compares the covariate and treatment distributions of a synthetic
#' (plasmode) sample against the real lock data, so the analyst can defend
#' that the plasmode generator is faithful enough to base candidate
#' selection on. This addresses the central caveat of the DQ stress test:
#' estimator rankings under synthetic outcomes only transfer to the real
#' study when the synthetic-data-generating process is faithful.
#'
#' For each covariate the function reports the absolute standardised mean
#' difference (SMD) between real and synthetic samples and, for continuous
#' covariates, a two-sample Kolmogorov-Smirnov statistic. It also reports
#' the absolute difference in treatment prevalence. A FLAG is raised when
#' any covariate SMD exceeds \code{smd_threshold} or any KS statistic
#' exceeds \code{ks_threshold}.
#'
#' @param real_data,synth_data data.frames sharing the covariate and
#'   treatment columns.
#' @param covariates Character vector of covariate column names.
#' @param treatment Character; treatment column name (optional).
#' @param smd_threshold,ks_threshold Fidelity-FLAG thresholds.
#'
#' @return An object of class \code{cleantmle_dgp_fidelity}: a list with a
#'   per-covariate \code{table}, the treatment-prevalence difference, and a
#'   \code{decision} ("GO" or "FLAG").
#'
#' @examples
#' \dontrun{
#' fid <- assess_dgp_fidelity(real, synth,
#'   covariates = c("age", "sex", "biomarker"), treatment = "treatment")
#' print(fid)
#' }
#' @export
assess_dgp_fidelity <- function(real_data, synth_data, covariates,
                                treatment = NULL,
                                smd_threshold = 0.10,
                                ks_threshold  = 0.10) {
  if (!is.data.frame(real_data) || !is.data.frame(synth_data))
    stop("`real_data` and `synth_data` must be data.frames.", call. = FALSE)
  miss <- covariates[!covariates %in% intersect(names(real_data),
                                                names(synth_data))]
  if (length(miss) > 0L)
    stop("covariates absent from one of the samples: ",
         paste(miss, collapse = ", "), call. = FALSE)

  rows <- lapply(covariates, function(v) {
    xr <- real_data[[v]]; xs <- synth_data[[v]]
    if (is.numeric(xr)) {
      mr <- mean(xr, na.rm = TRUE); ms <- mean(xs, na.rm = TRUE)
      sp <- sqrt((stats::var(xr, na.rm = TRUE) +
                  stats::var(xs, na.rm = TRUE)) / 2)
      smd <- if (is.finite(sp) && sp > 0) abs(mr - ms) / sp else NA_real_
      ks  <- tryCatch(unname(stats::ks.test(xr, xs)$statistic),
                      error = function(e) NA_real_, warning = function(w) {
                        suppressWarnings(unname(stats::ks.test(xr, xs)$statistic))
                      })
    } else {
      # Categorical: SMD on the proportion of the most common real level.
      lev <- names(sort(table(xr), decreasing = TRUE))[1]
      pr  <- mean(xr == lev, na.rm = TRUE); ps <- mean(xs == lev, na.rm = TRUE)
      sp  <- sqrt((pr * (1 - pr) + ps * (1 - ps)) / 2)
      smd <- if (is.finite(sp) && sp > 0) abs(pr - ps) / sp else NA_real_
      ks  <- NA_real_
    }
    data.frame(variable = v, smd = round(smd, 4), ks = round(ks, 4),
               stringsAsFactors = FALSE)
  })
  tbl <- do.call(rbind, rows)

  trt_prev_diff <- NA_real_
  if (!is.null(treatment) && treatment %in% names(real_data) &&
      treatment %in% names(synth_data)) {
    trt_prev_diff <- abs(mean(real_data[[treatment]] == 1, na.rm = TRUE) -
                         mean(synth_data[[treatment]] == 1, na.rm = TRUE))
  }

  flag <- any(tbl$smd > smd_threshold, na.rm = TRUE) ||
          any(tbl$ks  > ks_threshold,  na.rm = TRUE)
  decision <- if (flag) "FLAG" else "GO"

  out <- list(
    table         = tbl,
    trt_prev_diff = round(trt_prev_diff, 4),
    decision      = decision,
    thresholds    = list(smd_threshold = smd_threshold,
                         ks_threshold = ks_threshold)
  )
  class(out) <- "cleantmle_dgp_fidelity"
  out
}

#' @export
print.cleantmle_dgp_fidelity <- function(x, ...) {
  cat("Plasmode DGP fidelity (real vs synthetic)\n")
  cat("=========================================\n")
  print(x$table, row.names = FALSE)
  if (!is.na(x$trt_prev_diff))
    cat(sprintf("\nTreatment-prevalence |diff|: %.4f\n", x$trt_prev_diff))
  cat(sprintf("Decision: %s (SMD>%.2f or KS>%.2f flags)\n",
              x$decision, x$thresholds$smd_threshold,
              x$thresholds$ks_threshold))
  invisible(x)
}
