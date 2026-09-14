# estimate_effect(): the one estimation front door.
#
# Every point-treatment effect the package can estimate is reached from here
# by two arguments: which estimand (ATE, trimmed ATE, ATT, ATO, matched ATT)
# and, for the ATE, which estimator family (TMLE, IPTW, matching, crude).
# The former sibling entry points (run_att_tmle, estimate_ato,
# run_trimmed_tmle, run_ipcw_tmle, run_iptw_workflow, run_match_workflow,
# run_crude_workflow, run_matched_tmle) are internal workers behind it, and
# the implausibility guard is attached to every result.

#' Estimate a Treatment Effect for a Declared Estimand
#'
#' The single estimation entry point. Choose the estimand; the function
#' chooses the right machinery: TMLE (optionally IPCW via the censoring
#' mechanism) for the ATE, TMLE on the common-support subset with a
#' propensity refit for the trimmed ATE, the complete-case ATT that never
#' goes through the censoring mechanism, the augmented overlap-weighted
#' estimator for the ATO, and TMLE on a caliper-matched cohort for the
#' matched ATT. For the ATE, `estimator` selects simpler comparison
#' estimators (stabilised IPTW, matched difference, crude difference) for
#' side-by-side reporting. Every result carries its estimand label and the
#' implausibility flags; support verdicts travel with results produced
#' through [run_estimand_ladder()].
#'
#' @param lock A `cleanroom_lock`.
#' @param ps_fit A `ps_fit` from [fit_ps()]; required for every estimand
#'   except a plain `estimator = "crude"`, and optional for `"ATE"` with
#'   `estimator = "tmle"` (tmle::tmle fits its own g there).
#' @param estimand One of `"ATE"`, `"trimmed_ATE"`, `"ATT"`, `"ATO"`,
#'   `"matched_ATT"`.
#' @param estimator For `estimand = "ATE"` only: `"tmle"` (default),
#'   `"iptw"`, `"match"`, or `"crude"`.
#' @param missing How missing outcomes are handled for the ATE:
#'   `"auto"` (IPCW when the outcome has missing values, else complete
#'   case), `"ipcw"`, or `"complete_case"`. The ATT is always complete
#'   case by construction; see [run_estimand_ladder()] for why.
#' @param family `"binomial"` or `"gaussian"`.
#' @param sl_library SuperLearner library for the nuisance models; defaults
#'   to the ladder candidate's library, else the lock library.
#' @param trim_levels Trim levels for `estimand = "trimmed_ATE"`, tried in
#'   order. Default `c(0.05, 0.10)`.
#' @param trim_rule `"fixed"` or `"crump"` for the trimmed ATE.
#' @param gbound,cv_folds,prescreen_g,seed,caliper_sd Passed to the
#'   underlying machinery; see the package vignette.
#' @param return_steps For `estimand = "ATE", estimator = "tmle",
#'   missing = "complete_case"`: also return the modular four-step pieces
#'   (treatment mechanism, outcome mechanism, targeting, extraction) for
#'   teaching and diagnostics. Default FALSE.
#' @param allow_outcome_access Bypass the outcome guard. Default FALSE.
#' @param verbose Print progress. Default FALSE.
#' @return An effect-estimate object (class depends on the machinery;
#'   all inherit `cr_result`) with `estimate`, `se`, `ci_lower`,
#'   `ci_upper`, `p_value`, `n`, the named `estimand`, and the
#'   implausibility flags.
#' @examples
#' \dontrun{
#' lock <- create_analysis_lock(dat, "treatment", "outcome", covs)
#' ps   <- fit_ps(lock, "glm")
#' estimate_effect(lock, ps, estimand = "ATT")
#' estimate_effect(lock, ps, estimand = "ATO")
#' estimate_effect(lock, ps, estimand = "trimmed_ATE", trim_rule = "crump")
#' estimate_effect(lock, estimand = "ATE", missing = "ipcw")
#' }
#' @export
estimate_effect <- function(lock, ps_fit = NULL,
                            estimand = c("ATE", "trimmed_ATE", "ATT",
                                         "ATO", "matched_ATT"),
                            estimator = c("tmle", "iptw", "match", "crude"),
                            missing = c("auto", "ipcw", "complete_case"),
                            family = "binomial",
                            sl_library = NULL,
                            trim_levels = c(0.05, 0.10),
                            trim_rule = c("fixed", "crump"),
                            gbound = NULL,
                            cv_folds = 10L,
                            prescreen_g = FALSE,
                            seed = NULL,
                            caliper_sd = 0.2,
                            return_steps = FALSE,
                            allow_outcome_access = FALSE,
                            verbose = FALSE) {
  if (!inherits(lock, "cleanroom_lock"))
    stop("`lock` must be a cleanroom_lock object.", call. = FALSE)
  estimand  <- match.arg(estimand)
  estimator <- match.arg(estimator)
  missing   <- match.arg(missing)
  trim_rule <- match.arg(trim_rule)
  if (is.null(seed)) seed <- lock$seed
  Y <- lock$data[[lock$outcome]]
  use_ipcw <- switch(missing, auto = anyNA(Y), ipcw = TRUE,
                     complete_case = FALSE)
  need_ps <- estimand != "ATE" || estimator %in% c("iptw", "match")
  if (need_ps && !inherits(ps_fit, "ps_fit"))
    stop("estimand '", estimand, "' (or estimator '", estimator,
         "') needs a ps_fit from fit_ps().", call. = FALSE)

  if (estimand == "ATT") {
    if (isTRUE(use_ipcw) && identical(missing, "ipcw"))
      warning("The ATT is estimated complete case by construction; under ",
              "the censoring mechanism it loses double robustness. ",
              "Ignoring missing = 'ipcw'.", call. = FALSE)
    return(run_att_tmle(lock, family = family, sl_library = sl_library,
                        gbound = gbound, cv_folds = cv_folds,
                        prescreen_g = prescreen_g, seed = seed,
                        allow_outcome_access = allow_outcome_access))
  }
  if (estimand == "ATO")
    return(estimate_ato(lock, ps_fit, sl_library = sl_library,
                        family = family,
                        allow_outcome_access = allow_outcome_access))
  if (estimand == "trimmed_ATE")
    return(run_trimmed_tmle(lock, ps_fit, levels = trim_levels,
                            rule = trim_rule, family = family,
                            sl_library = sl_library, gbound = gbound,
                            cv_folds = cv_folds, prescreen_g = prescreen_g,
                            use_ipcw = use_ipcw, seed = seed,
                            allow_outcome_access = allow_outcome_access,
                            verbose = verbose))
  if (estimand == "matched_ATT") {
    g <- pmin(pmax(as.numeric(ps_fit$ps_raw %||% ps_fit$ps), 1e-6), 1 - 1e-6)
    A <- as.integer(lock$data[[lock$treatment]])
    set.seed(seed)
    mm <- .greedy_caliper_match(g, A, caliper_sd = caliper_sd)
    if (length(mm$treated) < 10L)
      stop("matched_ATT: fewer than 10 matched pairs.", call. = FALSE)
    return(run_matched_tmle(lock, ps_fit,
                            subset_idx = sort(c(mm$treated, mm$control)),
                            sl_library = sl_library,
                            override_clean_room = allow_outcome_access))
  }

  # estimand == "ATE"
  if (estimator == "crude")
    return(run_crude_workflow(lock,
                              allow_outcome_access = allow_outcome_access))
  if (estimator == "iptw")
    return(run_iptw_workflow(lock, ps_fit,
                             allow_outcome_access = allow_outcome_access))
  if (estimator == "match")
    return(run_match_workflow(lock, ps_fit,
                              allow_outcome_access = allow_outcome_access))
  # TMLE.
  if (use_ipcw)
    return(run_ipcw_tmle(lock, ps_fit = ps_fit,
                         allow_outcome_access = allow_outcome_access))
  if (isTRUE(return_steps)) {
    g_fit <- fit_tmle_treatment_mechanism(lock, ps_fit = ps_fit)
    q_fit <- fit_tmle_outcome_mechanism(lock, g_fit,
                                        sl_library = sl_library,
                                        allow_outcome_access =
                                          allow_outcome_access)
    upd <- run_tmle_targeting_step(g_fit, q_fit)
    est <- extract_tmle_estimate(upd)
    est$steps <- list(treatment = g_fit, outcome = q_fit, targeting = upd)
    guard <- implausibility_check(est$estimates$ATE$estimate,
                                  as.numeric(Y),
                                  as.numeric(lock$data[[lock$treatment]]),
                                  family)
    est$implausible <- guard$implausible
    est$implausible_reason <- guard$implausible_reason
    est$crude_diff <- guard$crude_diff
    est$estimand <- "ATE (complete case)"
    return(est)
  }
  args <- .tmle_delegate_args(lock, family = family, use_delta = FALSE,
                              sl_library = sl_library, gbound = gbound,
                              cv_folds = cv_folds,
                              prescreen_g = prescreen_g)
  .check_outcome_access(lock, allow_outcome_access,
                        caller = "estimate_effect")
  if (!requireNamespace("tmle", quietly = TRUE))
    stop("Package 'tmle' is required.", call. = FALSE)
  set.seed(seed)
  f <- do.call(tmle::tmle, args)
  est <- f$estimates$ATE
  guard <- implausibility_check(unname(est$psi), args$Y, args$A, family)
  out <- list(estimate = unname(est$psi), se = unname(sqrt(est$var.psi)),
              ci_lower = unname(est$CI[1]), ci_upper = unname(est$CI[2]),
              p_value = unname(est$pvalue),
              estimand = "ATE (complete case)",
              n = length(args$Y),
              risk_treated = tryCatch(unname(f$estimates$EY1$psi),
                                      error = function(e) NA_real_),
              risk_control = tryCatch(unname(f$estimates$EY0$psi),
                                      error = function(e) NA_real_),
              crude_diff = guard$crude_diff,
              implausible = guard$implausible,
              implausible_reason = guard$implausible_reason,
              tmle_fit = f, treatment = lock$treatment,
              outcome = lock$outcome, type = "ate_tmle",
              call = match.call())
  class(out) <- c("tmle_fit", "cr_result")
  out
}
