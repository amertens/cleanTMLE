# ============================================================================
# Split, enforced clean-room entry point
# ----------------------------------------------------------------------------
# run_clean_tmle() builds its lock internally and reads the outcome
# unconditionally, so it cannot enforce a pre-outcome gate (the lock hash is
# not known until after the call). These functions split the pipeline at the
# gate:
#
#   pre  <- run_clean_tmle_preoutcome(...)          # Y not authorised
#   auth <- authorize_outcome_analysis(pre$audit)   # the token
#   fit  <- run_clean_tmle_primary(pre, auth)       # refuses without a valid token
#
# The token (a pre_outcome_gate) carries the lock hash and an audit
# fingerprint; run_clean_tmle_primary() refuses to read the outcome unless the
# token authorises, matches the lock, and matches the current audit.
# ============================================================================

#' Pre-Outcome Pass of the Staged Clean-Room Workflow
#'
#' Runs the outcome-blind stages of the workflow (cohort adequacy,
#' propensity-score balance, optional outcome-blind candidate selection and
#' data-quality stress test, optional negative controls) and assembles a
#' reviewer-facing dossier and a pre-outcome gate. The primary
#' treatment-outcome association is **not** estimated and the returned lock is
#' left unauthorised: pass the result to [authorize_outcome_analysis()] and
#' then [run_clean_tmle_primary()] to read the outcome.
#'
#' @param data A data frame.
#' @param Avar,Yvar Treatment and outcome column names.
#' @param covariates Character vector of covariates; if `NULL`, all numeric
#'   columns other than the treatment, outcome, and common id/time columns.
#' @param learner_lib Pre-specified SuperLearner library.
#' @param tmle_candidates Optional list of TMLE candidates for the plasmode
#'   loop; if `NULL`, the default candidate grid is used.
#' @param selection_rule Candidate-selection rule (see
#'   [select_tmle_candidate()]).
#' @param dq_scenarios Optional DQ-stress scenarios (see
#'   [run_plasmode_dq_stress()]); when supplied, a DQ gate is recorded.
#' @param neg_control_outcomes Optional character vector of negative-control
#'   outcome columns.
#' @param plasmode_reps Replicates for the plasmode loop.
#' @param run_selection Logical; run the outcome-blind candidate selection /
#'   DQ stage. Default `TRUE`.
#' @param ps_method `"superlearner"` (default) or `"glm"` for the
#'   propensity-score fit.
#' @param cohort_thresholds,balance_thresholds Named lists of thresholds passed
#'   to [checkpoint_cohort_adequacy()] and [checkpoint_balance()].
#' @param seed Integer seed recorded on the lock.
#' @param allow_flag Logical; whether a FLAG gate still authorises (default
#'   `TRUE`).
#' @param verbose Logical; emit stage messages.
#'
#' @return A `clean_tmle_preoutcome` object with the (unauthorised) `lock`,
#'   the `audit`, the `gate` token, the `dossier`, and the design diagnostics.
#'
#' @seealso [run_clean_tmle_primary()], [authorize_outcome_analysis()]
#' @export
run_clean_tmle_preoutcome <- function(data, Avar, Yvar, covariates = NULL,
                                      learner_lib,
                                      tmle_candidates      = NULL,
                                      selection_rule       = "min_rmse",
                                      dq_scenarios         = NULL,
                                      neg_control_outcomes = NULL,
                                      plasmode_reps        = 50L,
                                      run_selection        = TRUE,
                                      ps_method            = c("superlearner", "glm"),
                                      cohort_thresholds    = list(),
                                      balance_thresholds   = list(),
                                      seed                 = 42L,
                                      allow_flag           = TRUE,
                                      verbose              = TRUE) {
  if (missing(learner_lib))
    stop("A pre-specified learner library (`learner_lib`) is required.",
         call. = FALSE)
  ps_method <- match.arg(ps_method)

  if (is.null(covariates)) {
    exclude <- c(Avar, Yvar, "id", "subject", "time", "event",
                 "event_type", "censored")
    covariates <- setdiff(
      names(data)[vapply(data, is.numeric, logical(1))], exclude)
  }

  lock <- create_analysis_lock(
    data = data, treatment = Avar, outcome = Yvar, covariates = covariates,
    sl_library = learner_lib, plasmode_reps = as.integer(plasmode_reps),
    seed = as.integer(seed))
  audit <- create_audit_log(lock)

  # Check Point 1: cohort adequacy (marginal Y only).
  cp1 <- do.call(checkpoint_cohort_adequacy, c(list(lock), cohort_thresholds))
  audit <- record_checkpoint(audit, cp1)
  if (verbose) message("Check Point 1 (cohort adequacy): ", cp1$decision)

  # Check Point 2: propensity-score overlap and balance (no Y).
  ps_fit <- if (ps_method == "glm") fit_ps_glm(lock) else fit_ps_superlearner(lock)
  psd    <- compute_ps_diagnostics(ps_fit)
  cp2    <- do.call(checkpoint_balance,
                    c(list(psd), balance_thresholds,
                      list(lock_hash = lock$lock_hash)))
  audit  <- record_checkpoint(audit, cp2)
  if (verbose) message("Check Point 2 (balance): ", cp2$decision)

  required_stages <- c("Check Point 1", "Check Point 2")

  # Stage 2b/2c: outcome-blind candidate selection and DQ stress (synthetic Y).
  plasmode_results <- NULL; dq_results <- NULL; selected <- NULL
  if (isTRUE(run_selection)) {
    plasmode_results <- run_plasmode_feasibility(
      lock, tmle_candidates = tmle_candidates,
      reps = as.integer(plasmode_reps), verbose = verbose)
    if (!is.null(dq_scenarios)) {
      dq_results <- run_plasmode_dq_stress(
        lock, tmle_candidates = tmle_candidates, scenarios = dq_scenarios,
        verbose = verbose)
    }
    selected <- select_tmle_candidate(plasmode_results, rule = selection_rule,
                                      dq_results = dq_results)
    # Bind the selected candidate to the lock so the primary analysis honours it
    # (truncation and nuisance library), rather than silently falling back to
    # the default locked library.
    lock <- lock_primary_tmle_spec(lock, selected)
    audit <- record_stage(
      audit, "Stage 2b",
      paste("Outcome-blind candidate selection; selected:",
            selected$candidate_id), decision = "GO")
    if (!is.null(dq_results)) {
      cp_dq <- gate_dq(dq_results, candidate = selected$candidate_id)
      audit <- record_checkpoint(audit, cp_dq)
      if (verbose) message("Check Point 2c (DQ stress): ", cp_dq$decision)
    }
  }

  # Check Point 3: negative-control residual-bias probe (NC Y only).
  nc_results <- NULL
  if (!is.null(neg_control_outcomes)) {
    nc_results <- list()
    for (nc_var in neg_control_outcomes) {
      if (!nc_var %in% names(data)) {
        warning("Negative-control variable '", nc_var,
                "' not found in data; skipping.", call. = FALSE)
        next
      }
      lock_nc <- define_negative_control(lock, nc_var)
      nc_results[[nc_var]] <- run_negative_control(lock_nc, nc_var, ps_fit)
    }
    if (length(nc_results) > 0L) {
      cp3 <- checkpoint_residual_bias(nc_results, lock_hash = lock$lock_hash)
      audit <- record_checkpoint(audit, cp3)
      required_stages <- c(required_stages, "Check Point 3")
      if (verbose) message("Check Point 3 (residual bias): ", cp3$decision)
    }
  }

  gate <- authorize_outcome_analysis(audit, required_stages = required_stages,
                                     allow_flag = allow_flag)
  if (verbose) message("Pre-outcome gate: ", gate$decision)

  dossier <- build_dossier(
    lock = lock, audit = audit, ps_diagnostics = psd,
    plasmode_results = plasmode_results, dq_results = dq_results,
    selected = selected, nc_results = nc_results, gate = gate)

  structure(list(
    lock = lock, audit = audit, gate = gate, dossier = dossier,
    ps_fit = ps_fit, ps_diagnostics = psd,
    plasmode_results = plasmode_results, dq_results = dq_results,
    selected = selected, nc_results = nc_results,
    required_stages = required_stages),
    class = "clean_tmle_preoutcome")
}


#' Primary Pass of the Staged Clean-Room Workflow
#'
#' Reads the observed outcome and runs the locked primary TMLE, but only after
#' verifying the pre-outcome authorisation. The authorisation (a
#' `pre_outcome_gate` from [authorize_outcome_analysis()]) must authorise the
#' analysis, match the lock hash, and match the current audit fingerprint;
#' otherwise the outcome is not read.
#'
#' @param pre A `clean_tmle_preoutcome` object.
#' @param authorization A `pre_outcome_gate` token. If `NULL`, `pre$gate` is
#'   used. The token must have `authorized = TRUE`, a `lock_hash` equal to
#'   `pre$lock$lock_hash`, and an `audit_fingerprint` equal to the current
#'   `pre$audit`.
#' @param allow_outcome_access Logical escape hatch; when `TRUE`, the token
#'   checks are skipped and the outcome is read (a governance override). Default
#'   `FALSE`.
#' @param truncate Propensity-score truncation for the primary TMLE.
#' @param verbose Logical; emit the estimate.
#'
#' @return A `clean_tmle_primary` object with the `tmle_result`, the
#'   `risk_difference`, the CI, the (now authorised) `lock`, and the `dossier`.
#'
#' @seealso [run_clean_tmle_preoutcome()], [authorize_outcome_analysis()]
#' @export
run_clean_tmle_primary <- function(pre, authorization = NULL,
                                   allow_outcome_access = FALSE,
                                   truncate = 0.01, verbose = TRUE) {
  if (!inherits(pre, "clean_tmle_preoutcome"))
    stop("`pre` must be a clean_tmle_preoutcome object from ",
         "run_clean_tmle_preoutcome().", call. = FALSE)
  lock <- pre$lock

  if (!isTRUE(allow_outcome_access)) {
    auth <- if (!is.null(authorization)) authorization else pre$gate
    if (is.null(auth))
      stop("run_clean_tmle_primary(): the primary analysis requires an ",
           "authorization from authorize_outcome_analysis(). None was supplied ",
           "and `pre` carries no gate. Run the pre-outcome pass and authorise ",
           "the gate, or pass allow_outcome_access = TRUE to override.",
           call. = FALSE)
    if (!isTRUE(auth$authorized) ||
        identical(as.character(auth$decision), "STOP"))
      stop("run_clean_tmle_primary(): the pre-outcome gate did not authorise ",
           "outcome analysis (decision: ", auth$decision,
           "). Primary analysis withheld.", call. = FALSE)
    if (!identical(as.character(auth$lock_hash), as.character(lock$lock_hash)))
      stop("run_clean_tmle_primary(): the authorization does not match this ",
           "lock (hash mismatch). The spec or data changed after authorization ",
           "was granted.", call. = FALSE)
    # Fail closed: the token must carry a non-NA audit fingerprint that matches
    # the current audit. A token minted from an empty audit has an NA
    # fingerprint (and a GO decision), so accepting NA would let a forged
    # empty-audit token bypass the gate. Refuse when either side is NA or they
    # differ.
    fp_now <- .audit_fingerprint(pre$audit)
    fp_tok <- if (is.null(auth$audit_fingerprint)) NA_character_
              else as.character(auth$audit_fingerprint)
    if (is.na(fp_tok) || is.na(fp_now) ||
        !identical(fp_tok, as.character(fp_now)))
      stop("run_clean_tmle_primary(): the authorization is not bound to the ",
           "current audit (audit fingerprint missing or mismatched). Re-authorise ",
           "on the recorded checkpoints before reading the outcome.", call. = FALSE)
    # Token verified; stamp the lock so the outcome guard permits Stage 4.
    lock$.outcome_authorized <- TRUE
  } else {
    lock$.outcome_authorized <- TRUE
  }
  .check_outcome_access(lock, allow_outcome_access = allow_outcome_access,
                        caller = "run_clean_tmle_primary")

  # Honour the locked primary candidate (its outcome-model library and
  # truncation) when the pre-outcome pass selected one; otherwise fall back to
  # the lock's default library and the `truncate` argument.
  spec   <- lock$primary_tmle_spec
  sl_lib <- if (!is.null(spec) && !is.null(spec$q_library)) spec$q_library
            else lock$sl_library
  trunc  <- if (!is.null(spec) && !is.null(spec$truncation)) spec$truncation
            else truncate

  tmle_result <- estimate_tmle_risk_point(
    data = lock$data, treatment = lock$treatment, outcome = lock$outcome,
    covariates = lock$covariates, sl_library = sl_lib,
    truncate = trunc)

  ate <- tmle_result$estimates$ATE
  if (verbose && !is.null(ate))
    message(sprintf("Primary TMLE risk difference %.4f (95%% CI %.4f, %.4f)",
                    ate$estimate, ate$ci_lower, ate$ci_upper))

  structure(list(
    tmle_result     = tmle_result,
    risk_difference = if (!is.null(ate)) ate$estimate else NA_real_,
    ci              = if (!is.null(ate)) c(ate$ci_lower, ate$ci_upper) else c(NA, NA),
    lock            = lock,
    authorization   = if (!is.null(authorization)) authorization else pre$gate,
    dossier         = pre$dossier),
    class = "clean_tmle_primary")
}


#' Assemble a Pre-Outcome Study Dossier
#'
#' Bundles the outcome-blind design artefacts a review team reads before
#' authorising the primary analysis: the estimand and lock fingerprint,
#' propensity-score balance, outcome-blind candidate selection and DQ
#' degradation, negative-control estimates, and the pre-outcome decision with
#' its audit trail. No comparative treatment-outcome estimate is included.
#'
#' @param lock A `cleanroom_lock`.
#' @param audit A `cleantmle_audit`.
#' @param ps_diagnostics Optional `ps_diagnostics`.
#' @param plasmode_results,dq_results,selected Optional candidate-selection
#'   artefacts.
#' @param nc_results Optional named list of negative-control results.
#' @param gate Optional `pre_outcome_gate`.
#'
#' @return A `clean_tmle_dossier` object.
#' @export
build_dossier <- function(lock, audit, ps_diagnostics = NULL,
                          plasmode_results = NULL, dq_results = NULL,
                          selected = NULL, nc_results = NULL, gate = NULL) {
  estimand <- list(
    treatment = lock$treatment, outcome = lock$outcome,
    covariates = lock$covariates, estimand = "marginal risk difference",
    sl_library = lock$sl_library, lock_hash = lock$lock_hash)

  balance <- if (!is.null(ps_diagnostics)) list(
    max_weighted_smd = tryCatch(max(abs(ps_diagnostics$smds$smd_weighted)),
                                error = function(e) NA_real_),
    smds = ps_diagnostics$smds, ess = ps_diagnostics$ess) else NULL

  selection <- if (!is.null(selected)) list(
    candidate_id = selected$candidate_id, metrics = selected$metrics,
    dq_degradation = if (!is.null(dq_results))
      tryCatch(summarize_dq_degradation(dq_results),
               error = function(e) NULL) else NULL) else NULL

  negative_controls <- if (!is.null(nc_results) && length(nc_results) > 0L)
    do.call(rbind, lapply(names(nc_results), function(v) data.frame(
      variable = v,
      estimate = nc_results[[v]]$estimate %||% NA_real_,
      p_value  = nc_results[[v]]$p_value %||% NA_real_,
      stringsAsFactors = FALSE))) else NULL

  decision <- list(
    decision   = if (!is.null(gate)) gate$decision else NA_character_,
    authorized = if (!is.null(gate)) isTRUE(gate$authorized) else NA,
    rationale  = if (!is.null(gate)) gate$rationale else NA_character_,
    audit      = tryCatch(export_decision_log(audit), error = function(e) NULL))

  structure(list(
    estimand = estimand, balance = balance, selection = selection,
    negative_controls = negative_controls, decision = decision),
    class = "clean_tmle_dossier")
}


#' @export
print.clean_tmle_dossier <- function(x, ...) {
  cat("== Pre-outcome study dossier ==\n")
  cat(sprintf("Estimand : marginal risk difference of '%s' on '%s'\n",
              x$estimand$treatment, x$estimand$outcome))
  cat(sprintf("Lock hash: %s\n", x$estimand$lock_hash))
  if (!is.null(x$balance) && !is.na(x$balance$max_weighted_smd))
    cat(sprintf("Balance  : max weighted SMD %.3f\n", x$balance$max_weighted_smd))
  if (!is.null(x$selection))
    cat(sprintf("Selected : candidate '%s'\n", x$selection$candidate_id))
  if (!is.null(x$negative_controls))
    cat(sprintf("Neg ctrls: %d evaluated\n", nrow(x$negative_controls)))
  cat(sprintf("Decision : %s (authorised: %s)\n",
              x$decision$decision, x$decision$authorized))
  cat("\n(No comparative treatment-outcome estimate is included in the dossier.)\n")
  invisible(x)
}


#' @export
print.clean_tmle_preoutcome <- function(x, ...) {
  cat("== clean_tmle pre-outcome pass ==\n")
  cat(sprintf("Lock hash: %s\n", x$lock$lock_hash))
  cat(sprintf("Gate     : %s (authorised: %s)\n",
              x$gate$decision, isTRUE(x$gate$authorized)))
  cat("Outcome  : NOT authorised; run authorize_outcome_analysis() then ",
      "run_clean_tmle_primary().\n", sep = "")
  invisible(x)
}


#' @export
print.clean_tmle_primary <- function(x, ...) {
  cat("== clean_tmle primary analysis ==\n")
  if (!is.na(x$risk_difference))
    cat(sprintf("Risk difference %.4f (95%% CI %.4f, %.4f)\n",
                x$risk_difference, x$ci[1], x$ci[2]))
  cat(sprintf("Lock hash: %s\n", x$lock$lock_hash))
  invisible(x)
}
