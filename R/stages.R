# ── Staged Clean-Room Workflow Infrastructure ──────────────────────────────
#
# This module provides S3 classes and helper functions for Muntner-style
# staged clean-room analyses.  Every major stage function returns a typed
# object that carries metadata (stage name, timestamp, lock hash) so that
# an audit trail can be assembled automatically.
#
# Stage mapping (Muntner et al.):
#   Stage 1a  – analysis specification and lock
#   Stage 1b  – cohort adequacy / Check Point 1
#   Stage 2   – treatment-arm comparability / Check Point 2
#   Stage 3   – residual-bias assessment / Check Point 3
#   Stage 4   – comparative analysis
#
# All functions in this file are pre-outcome unless explicitly documented
# otherwise.


# ── Estimand attachment ───────────────────────────────────────────────────

#' Attach an Estimand to an Analysis Lock
#'
#' Enriches a \code{cleanroom_lock} with a structured description of the
#' causal and statistical estimand.
#'
#' @section Clean-room stage: Stage 1a (pre-outcome).
#'
#' @param lock A \code{cleanroom_lock} from \code{\link{create_analysis_lock}}.
#' @param description Character; free-text description of the causal question.
#' @param population Character; description of the target population.
#' @param treatment_strategies Character vector of length 2 describing
#'   the treatment contrast (e.g. \code{c("Drug A", "Placebo")}).
#' @param outcome_label Character; description of the primary outcome.
#' @param followup Character; follow-up horizon description (e.g. "24 months").
#' @param contrast Character; estimand type, e.g. \code{"risk_difference"},
#'   \code{"risk_ratio"}, \code{"hazard_ratio"}.
#' @param statistical_estimand Character; description of the statistical
#'   (observed-data) estimand under identification assumptions.
#'
#' @return A modified \code{cleanroom_lock} with element \code{estimand}.
#'   The lock hash is unchanged (estimand metadata does not invalidate the
#'   analytic specification).
#'
#' @examples
#' dat  <- sim_func1(n = 200, seed = 1)
#' lock <- create_analysis_lock(dat, "treatment", "event_24",
#'                              c("age", "sex", "biomarker"), seed = 1)
#' lock <- attach_estimand(lock,
#'   description           = "Effect of treatment on 24-month event risk",
#'   population            = "Adults with simulated disease",
#'   treatment_strategies  = c("Treatment", "Control"),
#'   outcome_label         = "Primary event by 24 months",
#'   followup              = "24 months",
#'   contrast              = "risk_difference",
#'   statistical_estimand  = "E[Y(1)] - E[Y(0)] under SUTVA + positivity"
#' )
#' print(lock)
#'
#' @export
attach_estimand <- function(lock,
                            description          = NULL,
                            population           = NULL,
                            treatment_strategies = NULL,
                            outcome_label        = NULL,
                            followup             = NULL,
                            contrast             = "risk_difference",
                            statistical_estimand = NULL) {
  if (!inherits(lock, "cleanroom_lock"))
    stop("`lock` must be a cleanroom_lock object.", call. = FALSE)

  lock$estimand <- list(
    description          = description,
    population           = population,
    treatment_strategies = treatment_strategies,
    outcome_label        = outcome_label,
    followup             = followup,
    contrast             = contrast,
    statistical_estimand = statistical_estimand
  )
  lock
}


# ── Sensitivity plan declaration ──────────────────────────────────────────

#' Declare a Sensitivity Analysis Plan
#'
#' Attaches a named sensitivity analysis plan to the analysis lock.  The
#' plan is stored as metadata and does not alter the analytic specification
#' hash.  Multiple calls append additional plans.
#'
#' @section Clean-room stage: Stage 1a (pre-outcome).
#'
#' @param lock A \code{cleanroom_lock}.
#' @param label Character; short name for this sensitivity analysis.
#' @param description Character; what the sensitivity analysis evaluates.
#' @param settings A named list of parameters (e.g.
#'   \code{list(truncation = c(0.01, 0.05, 0.10))}).
#'
#' @return Modified \code{cleanroom_lock} with appended sensitivity plan.
#'
#' @examples
#' dat  <- sim_func1(n = 200, seed = 1)
#' lock <- create_analysis_lock(dat, "treatment", "event_24",
#'                              c("age", "sex", "biomarker"), seed = 1)
#' lock <- declare_sensitivity_plan(lock,
#'   label       = "truncation_sensitivity",
#'   description = "Re-estimate under alternate PS truncation thresholds",
#'   settings    = list(truncation = c(0.01, 0.05, 0.10))
#' )
#' lock$sensitivity_plans
#'
#' @export
declare_sensitivity_plan <- function(lock, label, description = NULL,
                                     settings = list()) {
  if (!inherits(lock, "cleanroom_lock"))
    stop("`lock` must be a cleanroom_lock object.", call. = FALSE)
  if (!is.character(label) || length(label) != 1L)
    stop("`label` must be a single character string.", call. = FALSE)

  plan <- list(
    label       = label,
    description = description,
    settings    = settings,
    declared_at = Sys.time()
  )

  if (is.null(lock$sensitivity_plans)) lock$sensitivity_plans <- list()
  lock$sensitivity_plans[[label]] <- plan
  lock
}


# ── Checkpoint constructor ────────────────────────────────────────────────

#' Create a Stage Checkpoint Object
#'
#' Low-level constructor for \code{cleantmle_checkpoint} objects.  Higher-
#' level helpers (\code{checkpoint_cohort_adequacy},
#' \code{checkpoint_balance}, \code{checkpoint_residual_bias}) call this
#' internally.
#'
#' @param stage Character; stage label (e.g. \code{"Check Point 1"}).
#' @param decision Character; one of \code{"GO"}, \code{"FLAG"}, \code{"STOP"}.
#' @param metrics A data.frame of diagnostic metrics.
#' @param thresholds A named list of thresholds used.
#' @param rationale Character; human-readable justification.
#' @param lock_hash Character; hash from the analysis lock, for traceability.
#'
#' @return An object of class \code{cleantmle_checkpoint}.
#'
#' @keywords internal
#' @export
new_checkpoint <- function(stage, decision, metrics, thresholds,
                           rationale = "", lock_hash = NA_character_) {
  decision <- match.arg(decision, c("GO", "FLAG", "STOP"))
  obj <- list(
    stage      = stage,
    decision   = decision,
    metrics    = metrics,
    thresholds = thresholds,
    rationale  = rationale,
    lock_hash  = lock_hash,
    timestamp  = Sys.time()
  )
  class(obj) <- "cleantmle_checkpoint"
  obj
}


#' @export
print.cleantmle_checkpoint <- function(x, ...) {
  cat(sprintf("=== %s ===\n", x$stage))
  cat(sprintf("Decision:  %s\n", x$decision))
  cat(sprintf("Rationale: %s\n", x$rationale))
  if (!is.na(x$lock_hash))
    cat(sprintf("Lock hash: %s\n", x$lock_hash))
  cat(sprintf("Timestamp: %s\n", format(x$timestamp, "%Y-%m-%d %H:%M:%S")))
  cat("\nMetrics:\n")
  print(x$metrics, row.names = FALSE)
  invisible(x)
}


#' @export
as.data.frame.cleantmle_checkpoint <- function(x, ...) {
  data.frame(
    stage     = x$stage,
    decision  = x$decision,
    rationale = x$rationale,
    lock_hash = x$lock_hash,
    timestamp = format(x$timestamp, "%Y-%m-%d %H:%M:%S"),
    stringsAsFactors = FALSE
  )
}


# ── Check Point 1: Cohort Adequacy ───────────────────────────────────────

#' Check Point 1: Cohort Adequacy and Precision
#'
#' Generates a structured report assessing whether the analytic cohort has
#' adequate sample size, event counts, and positivity for the planned
#' analysis.  Returns a \code{cleantmle_checkpoint} with a GO / FLAG /
#' STOP decision.
#'
#' @section Clean-room stage: Stage 1b (pre-outcome for design-stage
#'   summaries; event counts use the outcome variable).
#'
#' @param lock A \code{cleanroom_lock}.
#' @param min_n_per_arm Integer; minimum acceptable sample size per
#'   treatment arm. Default: 50.
#' @param min_events Integer; minimum acceptable total outcome events.
#'   Default: 20.
#' @param min_prevalence Numeric; minimum outcome prevalence. Default: 0.01.
#'
#' @return A \code{cleantmle_checkpoint} for Check Point 1.
#'
#' @examples
#' dat  <- sim_func1(n = 500, seed = 1)
#' lock <- create_analysis_lock(dat, "treatment", "event_24",
#'                              c("age", "sex", "biomarker"), seed = 1)
#' cp1 <- checkpoint_cohort_adequacy(lock)
#' print(cp1)
#'
#' @export
checkpoint_cohort_adequacy <- function(lock,
                                       min_n_per_arm  = 50L,
                                       min_events     = 20L,
                                       min_prevalence = 0.01) {
  if (!inherits(lock, "cleanroom_lock"))
    stop("`lock` must be a cleanroom_lock object.", call. = FALSE)

  data <- lock$data
  A    <- data[[lock$treatment]]
  Y    <- data[[lock$outcome]]
  n    <- nrow(data)
  n1   <- sum(A == 1, na.rm = TRUE)
  n0   <- sum(A == 0, na.rm = TRUE)
  events_total <- sum(Y == 1, na.rm = TRUE)
  events_trt   <- sum(Y[A == 1] == 1, na.rm = TRUE)
  events_ctrl  <- sum(Y[A == 0] == 1, na.rm = TRUE)
  prevalence   <- mean(Y == 1, na.rm = TRUE)
  if (is.na(prevalence)) prevalence <- 0

  # Simple MDD proxy: risk difference detectable with 80% power
  p_bar <- prevalence
  mdd <- if (n1 > 0 && n0 > 0 && p_bar > 0 && p_bar < 1) {
    2.8 * sqrt(p_bar * (1 - p_bar) * (1/n1 + 1/n0))
  } else NA_real_

  metrics <- data.frame(
    metric = c("N total", "N treated", "N control",
               "Events total", "Events treated", "Events control",
               "Outcome prevalence", "MDD (approx)"),
    value  = c(n, n1, n0, events_total, events_trt, events_ctrl,
               round(prevalence, 4), round(mdd, 4)),
    stringsAsFactors = FALSE
  )

  # Decision logic
  flags <- character(0)
  if (n1 < min_n_per_arm) flags <- c(flags, "treated arm below min N")
  if (n0 < min_n_per_arm) flags <- c(flags, "control arm below min N")
  if (events_total < min_events) flags <- c(flags, "total events below minimum")
  if (is.na(prevalence) || prevalence < min_prevalence)
    flags <- c(flags, "outcome prevalence too low or NA")

  # Positivity red flag: any covariate with zero variance in an arm
  pos_flags <- vapply(lock$covariates, function(v) {
    x <- data[[v]]
    if (is.numeric(x)) {
      v1 <- var(x[A == 1], na.rm = TRUE)
      v0 <- var(x[A == 0], na.rm = TRUE)
      isTRUE(v1 == 0) || isTRUE(v0 == 0)
    } else {
      a1_vals <- x[!is.na(A) & A == 1]
      a0_vals <- x[!is.na(A) & A == 0]
      length(unique(a1_vals[!is.na(a1_vals)])) <= 1 ||
        length(unique(a0_vals[!is.na(a0_vals)])) <= 1
    }
  }, logical(1))
  if (any(pos_flags))
    flags <- c(flags,
               paste("positivity concern:",
                     paste(lock$covariates[pos_flags], collapse = ", ")))

  decision <- if (length(flags) == 0L) "GO"
              else if (events_total < min_events || n1 < 20 || n0 < 20) "STOP"
              else "FLAG"

  rationale <- if (length(flags) == 0L) {
    "Cohort size, event counts, and positivity adequate."
  } else {
    paste("Issues:", paste(flags, collapse = "; "))
  }

  thresholds <- list(
    min_n_per_arm  = min_n_per_arm,
    min_events     = min_events,
    min_prevalence = min_prevalence
  )

  new_checkpoint(
    stage      = "Check Point 1: Cohort Adequacy",
    decision   = decision,
    metrics    = metrics,
    thresholds = thresholds,
    rationale  = rationale,
    lock_hash  = lock$lock_hash
  )
}


# ── Check Point 2: Balance / Comparability ───────────────────────────────

#' Check Point 2: Treatment Arm Comparability
#'
#' Converts propensity-score diagnostics into a structured checkpoint
#' decision.  Evaluates standardised mean differences, effective sample
#' size, and overlap against user-specified thresholds.
#'
#' @section Clean-room stage: Stage 2 (pre-outcome).
#'
#' @param ps_diag A \code{ps_diagnostics} object from
#'   \code{\link{compute_ps_diagnostics}}.
#' @param max_smd Numeric; maximum tolerable absolute SMD after weighting.
#'   Default: 0.10.
#' @param min_ess_pct Numeric; minimum ESS as a percent of the original N.
#'   Default: 50.
#' @param lock_hash Character; optional lock hash for traceability.
#'
#' @return A \code{cleantmle_checkpoint} for Check Point 2.
#'
#' @examples
#' dat  <- sim_func1(n = 500, seed = 1)
#' lock <- create_analysis_lock(dat, "treatment", "event_24",
#'                              c("age", "sex", "biomarker"), seed = 1)
#' ps   <- fit_ps_glm(lock)
#' diag <- compute_ps_diagnostics(ps)
#' cp2  <- checkpoint_balance(diag, lock_hash = lock$lock_hash)
#' print(cp2)
#'
#' @export
checkpoint_balance <- function(ps_diag,
                               max_smd     = 0.10,
                               min_ess_pct = 50,
                               lock_hash   = NA_character_) {
  if (!inherits(ps_diag, "ps_diagnostics"))
    stop("`ps_diag` must be a ps_diagnostics object.", call. = FALSE)

  smds     <- ps_diag$smds
  ess_tbl  <- ps_diag$ess

  max_weighted_smd <- max(abs(smds$smd_weighted))
  total_ess_pct    <- ess_tbl$ess_pct[ess_tbl$group == "Total"]

  flags <- character(0)
  if (max_weighted_smd > max_smd)
    flags <- c(flags,
               sprintf("max weighted SMD = %.3f (> %.2f)",
                       max_weighted_smd, max_smd))
  if (total_ess_pct < min_ess_pct)
    flags <- c(flags,
               sprintf("total ESS%% = %.1f (< %.0f%%)",
                       total_ess_pct, min_ess_pct))

  metrics <- data.frame(
    variable    = c(smds$variable, "Total ESS %"),
    smd_before  = c(smds$smd_unweighted, NA),
    smd_after   = c(smds$smd_weighted, NA),
    ess_pct     = c(rep(NA, nrow(smds)), total_ess_pct),
    stringsAsFactors = FALSE
  )

  decision <- if (length(flags) == 0L) "GO"
              else if (max_weighted_smd > 0.20) "STOP"
              else "FLAG"

  rationale <- if (length(flags) == 0L) {
    "All covariates balanced and ESS adequate."
  } else {
    paste("Issues:", paste(flags, collapse = "; "))
  }

  new_checkpoint(
    stage      = "Check Point 2: Treatment Comparability",
    decision   = decision,
    metrics    = metrics,
    thresholds = list(max_smd = max_smd, min_ess_pct = min_ess_pct),
    rationale  = rationale,
    lock_hash  = lock_hash
  )
}


# ── Negative Control Framework ────────────────────────────────────────────

#' Define a Negative Control Outcome
#'
#' Registers a variable as a negative control outcome in the analysis lock.
#' A negative control outcome is one for which the treatment is believed to
#' have no causal effect; any estimated association would indicate residual
#' confounding.
#'
#' @section Clean-room stage: Stage 1a (pre-outcome).
#'
#' @param lock A \code{cleanroom_lock}.
#' @param variable Character; name of the negative control outcome column
#'   in \code{lock$data}.
#' @param type Character; type of negative control. One of
#'   \code{"outcome"} (negative control outcome) or \code{"exposure"}
#'   (negative control exposure). Default: \code{"outcome"}.
#' @param description Character; optional description of why this variable
#'   serves as a negative control.
#'
#' @return Modified \code{cleanroom_lock} with element
#'   \code{negative_controls}.
#'
#' @examples
#' dat  <- sim_func1(n = 200, seed = 1)
#' lock <- create_analysis_lock(dat, "treatment", "event_24",
#'                              c("age", "sex", "biomarker"), seed = 1)
#' lock <- define_negative_control(lock, "nc_outcome",
#'   description = "Outcome known to be unrelated to treatment")
#'
#' @export
define_negative_control <- function(lock, variable, type = "outcome",
                                    description = NULL) {
  if (!inherits(lock, "cleanroom_lock"))
    stop("`lock` must be a cleanroom_lock object.", call. = FALSE)
  type <- match.arg(type, c("outcome", "exposure"))
  if (!variable %in% names(lock$data))
    stop("Variable '", variable, "' not found in lock data.", call. = FALSE)

  nc <- list(
    variable    = variable,
    type        = type,
    description = description
  )

  if (is.null(lock$negative_controls)) lock$negative_controls <- list()
  lock$negative_controls[[variable]] <- nc
  lock
}


#' Run a Negative Control Analysis
#'
#' Estimates the association between treatment and a negative control
#' outcome using the same analytic design (covariates, PS model) as the
#' primary analysis.  Any non-null association suggests residual
#' confounding.
#'
#' @section Clean-room stage: Stage 3 (accesses the negative control
#'   outcome variable).
#'
#' @param lock A \code{cleanroom_lock} with registered negative controls.
#' @param variable Character; the negative control variable name.
#' @param ps_fit A \code{ps_fit} object from the Stage 2 PS estimation.
#'
#' @return A list of class \code{cleantmle_nc_result} with elements
#'   \code{variable}, \code{estimate}, \code{se}, \code{ci_lower},
#'   \code{ci_upper}, \code{p_value}, and \code{interpretation}.
#'
#' @examples
#' dat  <- sim_func1(n = 500, seed = 1)
#' lock <- create_analysis_lock(dat, "treatment", "event_24",
#'                              c("age", "sex", "biomarker"), seed = 1)
#' lock <- define_negative_control(lock, "nc_outcome")
#' ps   <- fit_ps_glm(lock)
#' nc_result <- run_negative_control(lock, "nc_outcome", ps)
#' print(nc_result)
#'
#' @export
run_negative_control <- function(lock, variable, ps_fit) {
  if (!inherits(lock, "cleanroom_lock"))
    stop("`lock` must be a cleanroom_lock object.", call. = FALSE)
  if (!inherits(ps_fit, "ps_fit"))
    stop("`ps_fit` must be a ps_fit object.", call. = FALSE)
  if (!variable %in% names(lock$data))
    stop("Variable '", variable, "' not found in lock data.", call. = FALSE)

  data <- lock$data
  A    <- data[[lock$treatment]]
  Y_nc <- data[[variable]]
  ps   <- ps_fit$ps
  n    <- nrow(data)

  # Stabilised IPTW
  p_trt <- mean(A)
  w     <- ifelse(A == 1, p_trt / ps, (1 - p_trt) / (1 - ps))

  r1 <- weighted.mean(Y_nc[A == 1], w[A == 1])
  r0 <- weighted.mean(Y_nc[A == 0], w[A == 0])
  rd <- r1 - r0

  se_sq_1 <- sum(w[A == 1]^2 * (Y_nc[A == 1] - r1)^2) / sum(w[A == 1])^2
  se_sq_0 <- sum(w[A == 0]^2 * (Y_nc[A == 0] - r0)^2) / sum(w[A == 0])^2
  se <- sqrt(se_sq_1 + se_sq_0)

  ci_lo   <- rd - 1.96 * se
  ci_hi   <- rd + 1.96 * se
  p_value <- 2 * pnorm(-abs(rd / se))

  # Interpretation
  interpretation <- if (p_value < 0.05) {
    "Association detected: potential residual confounding."
  } else {
    "No significant association: no evidence of residual confounding."
  }

  result <- list(
    variable       = variable,
    estimate       = rd,
    se             = se,
    ci_lower       = ci_lo,
    ci_upper       = ci_hi,
    p_value        = p_value,
    interpretation = interpretation
  )
  class(result) <- "cleantmle_nc_result"
  result
}


#' @export
print.cleantmle_nc_result <- function(x, ...) {
  cat("Negative Control Analysis\n")
  cat("=========================\n")
  cat(sprintf("Variable:    %s\n", x$variable))
  cat(sprintf("Estimate:    %.5f  (95%% CI: %.5f, %.5f)\n",
              x$estimate, x$ci_lower, x$ci_upper))
  cat(sprintf("SE:          %.5f   p-value: %.4f\n", x$se, x$p_value))
  cat(sprintf("Assessment:  %s\n", x$interpretation))
  invisible(x)
}


# ── Check Point 3: Residual Bias ─────────────────────────────────────────

#' Check Point 3: Residual Bias Assessment
#'
#' Evaluates negative control results to assess whether residual
#' confounding is likely present.  Returns a structured checkpoint.
#'
#' @section Clean-room stage: Stage 3 (after negative control analysis).
#'
#' @param nc_results A single \code{cleantmle_nc_result} or a list of them.
#' @param alpha Numeric; significance threshold. Default: 0.05.
#' @param max_nc_estimate Numeric; maximum tolerable absolute NC estimate.
#'   Default: \code{Inf} (rely on p-value only).
#' @param lock_hash Character; optional lock hash for traceability.
#'
#' @return A \code{cleantmle_checkpoint} for Check Point 3.
#'
#' @examples
#' dat  <- sim_func1(n = 500, seed = 1)
#' lock <- create_analysis_lock(dat, "treatment", "event_24",
#'                              c("age", "sex", "biomarker"), seed = 1)
#' lock <- define_negative_control(lock, "nc_outcome")
#' ps   <- fit_ps_glm(lock)
#' nc   <- run_negative_control(lock, "nc_outcome", ps)
#' cp3  <- checkpoint_residual_bias(nc, lock_hash = lock$lock_hash)
#' print(cp3)
#'
#' @export
checkpoint_residual_bias <- function(nc_results,
                                     alpha           = 0.05,
                                     max_nc_estimate = Inf,
                                     lock_hash       = NA_character_) {
  # Accept single result or list

  if (inherits(nc_results, "cleantmle_nc_result"))
    nc_results <- list(nc_results)

  rows <- lapply(nc_results, function(nc) {
    data.frame(
      variable = nc$variable,
      estimate = round(nc$estimate, 5),
      se       = round(nc$se, 5),
      p_value  = round(nc$p_value, 4),
      flagged  = nc$p_value < alpha || abs(nc$estimate) > max_nc_estimate,
      stringsAsFactors = FALSE
    )
  })
  metrics <- do.call(rbind, rows)
  n_flagged <- sum(metrics$flagged)

  decision <- if (n_flagged == 0L) "GO"
              else if (n_flagged <= length(nc_results) / 2) "FLAG"
              else "STOP"

  rationale <- if (n_flagged == 0L) {
    "No negative controls flagged; no evidence of residual confounding."
  } else {
    sprintf("%d of %d negative control(s) flagged.",
            n_flagged, length(nc_results))
  }

  new_checkpoint(
    stage      = "Check Point 3: Residual Bias",
    decision   = decision,
    metrics    = metrics,
    thresholds = list(alpha = alpha, max_nc_estimate = max_nc_estimate),
    rationale  = rationale,
    lock_hash  = lock_hash
  )
}


# ── Audit Trail ──────────────────────────────────────────────────────────

#' Create an Audit Log
#'
#' Initialises an empty audit log that accumulates entries as the
#' analysis progresses through stages.
#'
#' @param lock A \code{cleanroom_lock}.
#'
#' @return An object of class \code{cleantmle_audit}.
#'
#' @export
create_audit_log <- function(lock) {
  if (!inherits(lock, "cleanroom_lock"))
    stop("`lock` must be a cleanroom_lock object.", call. = FALSE)

  obj <- list(
    lock_hash  = lock$lock_hash,
    locked_at  = lock$locked_at,
    entries    = list(),
    created_at = Sys.time()
  )
  class(obj) <- "cleantmle_audit"
  obj
}


#' Record a Stage Entry in the Audit Log
#'
#' Appends a timestamped entry to the audit log.
#'
#' @param audit A \code{cleantmle_audit}.
#' @param stage Character; stage label.
#' @param action Character; description of what was done.
#' @param decision Character; checkpoint decision if applicable
#'   (\code{"GO"}, \code{"FLAG"}, \code{"STOP"}, or \code{NA}).
#' @param details Character; optional additional detail.
#'
#' @return Modified \code{cleantmle_audit}.
#'
#' @export
record_stage <- function(audit, stage, action, decision = NA_character_,
                         details = "") {
  if (!inherits(audit, "cleantmle_audit"))
    stop("`audit` must be a cleantmle_audit object.", call. = FALSE)

  entry <- list(
    stage     = stage,
    action    = action,
    decision  = decision,
    details   = details,
    timestamp = Sys.time()
  )

  audit$entries <- c(audit$entries, list(entry))
  audit
}


#' Record a Checkpoint in the Audit Log
#'
#' Convenience wrapper that extracts stage, decision, and rationale from
#' a \code{cleantmle_checkpoint} and appends them to the audit log.
#'
#' @param audit A \code{cleantmle_audit}.
#' @param checkpoint A \code{cleantmle_checkpoint}.
#'
#' @return Modified \code{cleantmle_audit}.
#'
#' @export
record_checkpoint <- function(audit, checkpoint) {
  if (!inherits(checkpoint, "cleantmle_checkpoint"))
    stop("`checkpoint` must be a cleantmle_checkpoint object.", call. = FALSE)

  record_stage(audit,
    stage    = checkpoint$stage,
    action   = paste("Checkpoint evaluated:", checkpoint$rationale),
    decision = checkpoint$decision
  )
}


#' Export an Audit Trail as a Data Frame
#'
#' Converts the audit log to a tidy data frame suitable for printing
#' or inclusion in a vignette.
#'
#' @param audit A \code{cleantmle_audit}.
#'
#' @return A data.frame with columns \code{stage}, \code{action},
#'   \code{decision}, \code{details}, and \code{timestamp}.
#'
#' @examples
#' dat  <- sim_func1(n = 200, seed = 1)
#' lock <- create_analysis_lock(dat, "treatment", "event_24",
#'                              c("age", "sex", "biomarker"), seed = 1)
#' audit <- create_audit_log(lock)
#' audit <- record_stage(audit, "Stage 1a", "Analysis lock created")
#' export_audit_trail(audit)
#'
#' @export
export_audit_trail <- function(audit) {
  if (!inherits(audit, "cleantmle_audit"))
    stop("`audit` must be a cleantmle_audit object.", call. = FALSE)

  if (length(audit$entries) == 0L) {
    return(data.frame(
      stage     = character(0),
      action    = character(0),
      decision  = character(0),
      details   = character(0),
      timestamp = character(0),
      stringsAsFactors = FALSE
    ))
  }

  rows <- lapply(audit$entries, function(e) {
    data.frame(
      stage     = e$stage,
      action    = e$action,
      decision  = if (is.na(e$decision)) "" else e$decision,
      details   = e$details,
      timestamp = format(e$timestamp, "%Y-%m-%d %H:%M:%S"),
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, rows)
}


#' @export
print.cleantmle_audit <- function(x, ...) {
  cat("cleanTMLE Audit Log\n")
  cat("====================\n")
  cat(sprintf("Lock hash:  %s\n", x$lock_hash))
  cat(sprintf("Entries:    %d\n", length(x$entries)))
  if (length(x$entries) > 0L) {
    cat("\n")
    trail <- export_audit_trail(x)
    print(trail, row.names = FALSE, right = FALSE)
  }
  invisible(x)
}


# ── Stage Manifest / Diagram ─────────────────────────────────────────────

#' Build a Stage Manifest
#'
#' Produces a compact textual summary of the staged analysis path taken,
#' based on the audit log.
#'
#' @param audit A \code{cleantmle_audit}.
#'
#' @return A character string (printed invisibly) showing the stage path.
#'
#' @export
build_stage_manifest <- function(audit) {
  if (!inherits(audit, "cleantmle_audit"))
    stop("`audit` must be a cleantmle_audit object.", call. = FALSE)

  if (length(audit$entries) == 0L) {
    msg <- "No stage entries recorded."
    cat(msg, "\n")
    return(invisible(msg))
  }

  lines <- vapply(seq_along(audit$entries), function(i) {
    e <- audit$entries[[i]]
    dec <- if (!is.na(e$decision)) sprintf(" [%s]", e$decision) else ""
    sprintf("  %d. %s: %s%s", i, e$stage, e$action, dec)
  }, character(1))

  header <- sprintf("Clean-Room Stage Path (lock: %s)", audit$lock_hash)
  sep    <- paste(rep("-", nchar(header)), collapse = "")
  msg    <- paste(c(header, sep, lines), collapse = "\n")
  cat(msg, "\n")
  invisible(msg)
}


# ── Sensitivity Analysis Helper ──────────────────────────────────────────

#' Run Truncation Sensitivity Analysis
#'
#' Re-estimates the primary IPTW-based risk difference under multiple
#' propensity score truncation thresholds to assess sensitivity of the
#' estimate to extreme weights.
#'
#' @section Clean-room stage: Stage 4 (post-outcome).
#'
#' @param lock A \code{cleanroom_lock}.
#' @param thresholds Numeric vector of truncation thresholds to evaluate.
#'   Default: \code{c(0.01, 0.025, 0.05, 0.10)}.
#' @param override_clean_room Logical; if \code{TRUE}, skips the outcome-access check.  Default \code{FALSE}.
#'
#' @return A data.frame with columns \code{truncation}, \code{estimate},
#'   \code{se}, \code{ci_lower}, \code{ci_upper}.
#'
#' @export
sensitivity_truncation <- function(lock,
                                   thresholds = c(0.01, 0.025, 0.05, 0.10),
                                   override_clean_room = FALSE) {
  if (!inherits(lock, "cleanroom_lock"))
    stop("`lock` must be a cleanroom_lock object.", call. = FALSE)
  .check_outcome_access(lock, override_clean_room,
                        caller = "sensitivity_truncation")

  rows <- lapply(thresholds, function(trunc) {
    ps_fit <- fit_ps_glm(lock, truncate = trunc)
    iptw   <- run_iptw_workflow(lock, ps_fit)
    data.frame(
      truncation = trunc,
      estimate   = round(iptw$estimate, 5),
      se         = round(iptw$se, 5),
      ci_lower   = round(iptw$ci_lower, 5),
      ci_upper   = round(iptw$ci_upper, 5),
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, rows)
}


#' Compute E-value for a Risk Ratio
#'
#' Computes the E-value (VanderWeele and Ding, 2017) which quantifies the
#' minimum strength of association that an unmeasured confounder would
#' need with both treatment and outcome to explain away the observed
#' risk ratio.
#'
#' @param rr Numeric; the observed risk ratio (must be > 0).
#' @param ci_bound Numeric; optional confidence interval bound (lower
#'   bound for RR > 1, upper bound for RR < 1) for the E-value of the
#'   confidence limit.
#'
#' @return A named numeric vector with \code{e_value} and optionally
#'   \code{e_value_ci}.
#'
#' @references
#' VanderWeele TJ, Ding P. Sensitivity analysis in observational
#' research: introducing the E-value. \emph{Ann Intern Med.}
#' 2017;167(4):268--274.
#'
#' @export
compute_evalue <- function(rr, ci_bound = NULL) {
  if (rr < 0)
    stop("`rr` must be positive.", call. = FALSE)

  .evalue_calc <- function(r) {
    if (r < 1) r <- 1 / r
    r + sqrt(r * (r - 1))
  }

  ev <- .evalue_calc(rr)
  result <- c(e_value = ev)

  if (!is.null(ci_bound)) {
    if (ci_bound < 0) stop("`ci_bound` must be positive.", call. = FALSE)
    if ((rr >= 1 && ci_bound <= 1) || (rr < 1 && ci_bound >= 1)) {
      result["e_value_ci"] <- 1
    } else {
      result["e_value_ci"] <- .evalue_calc(ci_bound)
    }
  }

  result
}


# ── Lock / Retrieve Primary TMLE Specification ───────────────────────────

#' Lock the Primary TMLE Specification
#'
#' Stores the selected TMLE candidate specification in the analysis lock
#' so that downstream Stage 4 functions automatically use it.
#'
#' @section Clean-room stage: Stage 2b (after candidate selection).
#'
#' @param lock A \code{cleanroom_lock}.
#' @param selected A \code{tmle_selected_spec} from
#'   \code{\link{select_tmle_candidate}}, or any
#'   \code{tmle_candidate_spec}.
#'
#' @return Modified \code{cleanroom_lock} with element
#'   \code{primary_tmle_spec}.
#'
#' @examples
#' dat  <- sim_func1(n = 200, seed = 1)
#' lock <- create_analysis_lock(dat, "treatment", "event_24",
#'                              c("age", "sex", "biomarker"), seed = 1)
#' spec <- tmle_candidate("glm_t01", "GLM trunc=0.01",
#'                        g_library = "SL.glm", truncation = 0.01)
#' lock <- lock_primary_tmle_spec(lock, spec)
#' get_primary_tmle_spec(lock)
#'
#' @export
lock_primary_tmle_spec <- function(lock, selected) {
  if (!inherits(lock, "cleanroom_lock"))
    stop("`lock` must be a cleanroom_lock object.", call. = FALSE)
  if (!inherits(selected, "tmle_candidate_spec"))
    stop("`selected` must be a tmle_candidate_spec or tmle_selected_spec.",
         call. = FALSE)

  lock$primary_tmle_spec <- selected
  lock
}


#' Retrieve the Locked Primary TMLE Specification
#'
#' Extracts the primary TMLE candidate specification from the lock.
#' Returns \code{NULL} if no specification has been locked.
#'
#' @param lock A \code{cleanroom_lock}.
#'
#' @return A \code{tmle_candidate_spec} or \code{tmle_selected_spec},
#'   or \code{NULL}.
#'
#' @export
get_primary_tmle_spec <- function(lock) {
  if (!inherits(lock, "cleanroom_lock"))
    stop("`lock` must be a cleanroom_lock object.", call. = FALSE)
  lock$primary_tmle_spec
}


# ── Design Precision ──────────────────────────────────────────────────────

#' Estimate Design-Stage Precision Summaries
#'
#' Computes pre-outcome precision summaries using only sample size and
#' event information available at the design stage.  No outcome modelling
#' is performed.
#'
#' @section Clean-room stage: Stage 1b (design-stage; uses outcome counts
#'   for power/precision summaries only, not for estimation).
#'
#' @param lock A \code{cleanroom_lock}.
#' @param target_mdd Numeric; an optional target minimum detectable
#'   difference on the risk-difference scale to compare against the
#'   computed MDD.  Default \code{NULL} (no comparison).
#'
#' @return A list of class \code{design_precision} containing:
#'   \describe{
#'     \item{n_total}{Total sample size.}
#'     \item{n_treated}{Number of treated subjects.}
#'     \item{n_control}{Number of control subjects.}
#'     \item{events_total}{Total outcome events.}
#'     \item{events_per_arm}{Named vector with events per arm.}
#'     \item{crude_rates}{Named vector of crude event rates per arm.}
#'     \item{prevalence}{Overall outcome prevalence.}
#'     \item{se_proxy}{SE proxy for the risk difference.}
#'     \item{ci_halfwidth}{95\% CI half-width proxy (1.96 * se_proxy).}
#'     \item{mdd_80}{MDD at 80\% power (2.8 * se_proxy).}
#'     \item{target_mdd}{The user-supplied target MDD (or \code{NULL}).}
#'     \item{mdd_feasible}{Logical; whether computed MDD <= target_mdd
#'       (only present when \code{target_mdd} is not \code{NULL}).}
#'   }
#'
#' @examples
#' dat  <- sim_func1(n = 500, seed = 1)
#' lock <- create_analysis_lock(dat, "treatment", "event_24",
#'                              c("age", "sex", "biomarker"), seed = 1)
#' dp <- estimate_design_precision(lock)
#' print(dp)
#'
#' @export
estimate_design_precision <- function(lock, target_mdd = NULL) {
  if (!inherits(lock, "cleanroom_lock"))
    stop("`lock` must be a cleanroom_lock object.", call. = FALSE)

  data      <- lock$data
  A         <- data[[lock$treatment]]
  Y         <- data[[lock$outcome]]
  n_total   <- nrow(data)
  n_treated <- sum(A == 1L)
  n_control <- sum(A == 0L)

  events_treated <- sum(Y[A == 1L])
  events_control <- sum(Y[A == 0L])
  events_total   <- events_treated + events_control

  rate_treated <- events_treated / n_treated
  rate_control <- events_control / n_control
  crude_rates  <- c(treated = rate_treated, control = rate_control)

  events_per_arm <- c(treated = events_treated, control = events_control)

  prevalence  <- events_total / n_total
  p           <- prevalence
  se_proxy    <- sqrt(p * (1 - p) * (1 / n_treated + 1 / n_control))
  ci_halfwidth <- 1.96 * se_proxy
  mdd_80       <- 2.8  * se_proxy

  result <- list(
    n_total        = n_total,
    n_treated      = n_treated,
    n_control      = n_control,
    events_total   = events_total,
    events_per_arm = events_per_arm,
    crude_rates    = crude_rates,
    prevalence     = prevalence,
    se_proxy       = se_proxy,
    ci_halfwidth   = ci_halfwidth,
    mdd_80         = mdd_80,
    target_mdd     = target_mdd
  )

  if (!is.null(target_mdd)) {
    result$mdd_feasible <- mdd_80 <= target_mdd
  }

  class(result) <- "design_precision"
  result
}


#' @export
print.design_precision <- function(x, ...) {
  cat("=== Design-Stage Precision Summary ===\n")
  tbl <- data.frame(
    Metric = c(
      "N total", "N treated", "N control",
      "Events total", "Events (treated)", "Events (control)",
      "Overall prevalence",
      "SE proxy (RD)", "95% CI half-width", "MDD (80% power)"
    ),
    Value = c(
      x$n_total, x$n_treated, x$n_control,
      x$events_total, x$events_per_arm[["treated"]], x$events_per_arm[["control"]],
      round(x$prevalence,    4),
      round(x$se_proxy,      5),
      round(x$ci_halfwidth,  5),
      round(x$mdd_80,        5)
    ),
    stringsAsFactors = FALSE
  )
  print(tbl, row.names = FALSE, right = FALSE)
  if (!is.null(x$target_mdd)) {
    cat(sprintf(
      "\nTarget MDD: %.5f  |  Feasible: %s\n",
      x$target_mdd,
      if (isTRUE(x$mdd_feasible)) "YES" else "NO"
    ))
  }
  invisible(x)
}


# ── Event Support Summary ─────────────────────────────────────────────────

#' Summarise Event Support by Treatment Arm
#'
#' Returns a data frame with event counts and rates per arm, and prints
#' a warning if any arm has fewer than 10 events (sparse-event flag).
#'
#' @section Clean-room stage: Stage 1b (design diagnostics).
#'
#' @param lock A \code{cleanroom_lock}.
#'
#' @return A data.frame with columns \code{arm}, \code{n},
#'   \code{events}, and \code{event_rate}.
#'
#' @examples
#' dat  <- sim_func1(n = 500, seed = 1)
#' lock <- create_analysis_lock(dat, "treatment", "event_24",
#'                              c("age", "sex", "biomarker"), seed = 1)
#' summarize_event_support(lock)
#'
#' @export
summarize_event_support <- function(lock) {
  if (!inherits(lock, "cleanroom_lock"))
    stop("`lock` must be a cleanroom_lock object.", call. = FALSE)

  data <- lock$data
  A    <- data[[lock$treatment]]
  Y    <- data[[lock$outcome]]

  n1 <- sum(A == 1L);  e1 <- sum(Y[A == 1L])
  n0 <- sum(A == 0L);  e0 <- sum(Y[A == 0L])
  nt <- length(A);     et <- e1 + e0

  out <- data.frame(
    arm        = c("Treated", "Control", "Total"),
    n          = c(n1, n0, nt),
    events     = c(e1, e0, et),
    event_rate = round(c(e1 / n1, e0 / n0, et / nt), 4),
    stringsAsFactors = FALSE
  )

  sparse_arms <- out$arm[out$arm != "Total" & out$events < 10L]
  if (length(sparse_arms) > 0L) {
    message(
      "Sparse-event warning: fewer than 10 events in arm(s): ",
      paste(sparse_arms, collapse = ", "),
      ". Estimates may be unstable."
    )
  }

  out
}


# ── Stage 3 Residual Confounding Wrapper ─────────────────────────────────

#' Stage 3: Residual Confounding Assessment
#'
#' Convenience wrapper that runs \code{\link{run_negative_control}} on
#' every registered negative control variable in the lock and returns a
#' structured stage result with an embedded checkpoint.
#'
#' @section Clean-room stage: Stage 3 (residual confounding).
#'
#' @param lock A \code{cleanroom_lock} with negative controls registered
#'   via \code{\link{define_negative_control}}.
#' @param ps_fit A \code{ps_fit} object from \code{\link{fit_ps_glm}}.
#' @param variables Character vector of negative-control variable names
#'   to evaluate.  If \code{NULL} (default), all variables registered in
#'   \code{lock$negative_controls} are used.
#' @param alpha Numeric; significance threshold for flagging a negative
#'   control.  Default \code{0.05}.
#' @param max_abs_estimate Numeric; absolute estimate threshold for
#'   flagging a negative control.  Default \code{0.05}.
#'
#' @return A list of class \code{residual_confounding_stage} with:
#'   \describe{
#'     \item{nc_results}{Named list of \code{cleantmle_nc_result} objects.}
#'     \item{summary_table}{data.frame summarising each NC variable.}
#'     \item{checkpoint}{A \code{cleantmle_checkpoint} (Check Point 3).}
#'     \item{n_controls}{Number of negative controls evaluated.}
#'     \item{n_flagged}{Number of flagged negative controls.}
#'     \item{proportion_flagged}{Proportion flagged.}
#'   }
#'
#' @examples
#' dat  <- sim_func1(n = 500, seed = 1)
#' lock <- create_analysis_lock(dat, "treatment", "event_24",
#'                              c("age", "sex", "biomarker"), seed = 1)
#' lock <- define_negative_control(lock, "nc_outcome")
#' ps   <- fit_ps_glm(lock)
#' stage3 <- run_residual_confounding_stage(lock, ps)
#' print(stage3)
#'
#' @export
run_residual_confounding_stage <- function(lock,
                                           ps_fit,
                                           variables        = NULL,
                                           alpha            = 0.05,
                                           max_abs_estimate = 0.05) {
  if (!inherits(lock, "cleanroom_lock"))
    stop("`lock` must be a cleanroom_lock object.", call. = FALSE)
  if (!inherits(ps_fit, "ps_fit"))
    stop("`ps_fit` must be a ps_fit object.", call. = FALSE)

  if (is.null(variables)) {
    variables <- names(lock$negative_controls)
    if (length(variables) == 0L)
      stop("No negative controls registered. Use define_negative_control() first.",
           call. = FALSE)
  }

  nc_results <- lapply(variables, function(v) {
    run_negative_control(lock, v, ps_fit)
  })
  names(nc_results) <- variables

  # Build summary table
  summary_table <- do.call(rbind, lapply(nc_results, function(nc) {
    flagged <- nc$p_value < alpha || abs(nc$estimate) > max_abs_estimate
    data.frame(
      variable         = nc$variable,
      estimate         = round(nc$estimate, 5),
      se               = round(nc$se,       5),
      p_value          = round(nc$p_value,  4),
      flagged          = flagged,
      stringsAsFactors = FALSE
    )
  }))
  rownames(summary_table) <- NULL

  n_controls         <- length(nc_results)
  n_flagged          <- sum(summary_table$flagged)
  proportion_flagged <- n_flagged / n_controls

  checkpoint <- checkpoint_residual_bias(
    nc_results,
    alpha           = alpha,
    max_nc_estimate = max_abs_estimate,
    lock_hash       = lock$lock_hash
  )

  result <- list(
    nc_results         = nc_results,
    summary_table      = summary_table,
    checkpoint         = checkpoint,
    n_controls         = n_controls,
    n_flagged          = n_flagged,
    proportion_flagged = proportion_flagged
  )
  class(result) <- "residual_confounding_stage"
  result
}


#' @export
print.residual_confounding_stage <- function(x, ...) {
  cat("=== Stage 3: Residual Confounding Assessment ===\n")
  cat(sprintf("Negative controls evaluated: %d\n", x$n_controls))
  cat(sprintf("Flagged:                     %d (%.1f%%)\n",
              x$n_flagged, 100 * x$proportion_flagged))
  cat("\nSummary table:\n")
  print(x$summary_table, row.names = FALSE)
  cat("\n")
  print(x$checkpoint)
  invisible(x)
}


# ── Pre-Outcome Gate ──────────────────────────────────────────────────────

#' Authorise Outcome Analysis
#'
#' Scans the audit log to verify that all required stage checkpoints have
#' been recorded and none resulted in an unconditional STOP.  Returns a
#' \code{pre_outcome_gate} checkpoint summarising the go/no-go decision.
#'
#' @section Clean-room stage: Pre-Outcome Gate (between Stage 3 and
#'   Stage 4).
#'
#' @param audit A \code{cleantmle_audit}.
#' @param required_stages Character vector of stage labels that must be
#'   present in the audit.  Partial matching is used.  Default:
#'   \code{c("Check Point 1", "Check Point 2", "Check Point 3")}.
#' @param allow_flag Logical; if \code{FALSE} any FLAG decision is treated
#'   as a STOP.  Default \code{TRUE}.
#' @param lock_hash Character; optional lock hash for traceability.
#'
#' @return A \code{cleantmle_checkpoint} of subclass
#'   \code{pre_outcome_gate} with an additional logical field
#'   \code{authorized}.
#'
#' @examples
#' dat   <- sim_func1(n = 500, seed = 1)
#' lock  <- create_analysis_lock(dat, "treatment", "event_24",
#'                               c("age", "sex", "biomarker"), seed = 1)
#' audit <- create_audit_log(lock)
#' ps    <- fit_ps_glm(lock)
#' cp1   <- checkpoint_cohort_adequacy(lock)
#' cp2   <- checkpoint_balance(lock, ps)
#' lock  <- define_negative_control(lock, "nc_outcome")
#' nc    <- run_negative_control(lock, "nc_outcome", ps)
#' cp3   <- checkpoint_residual_bias(nc, lock_hash = lock$lock_hash)
#' audit <- record_checkpoint(audit, cp1)
#' audit <- record_checkpoint(audit, cp2)
#' audit <- record_checkpoint(audit, cp3)
#' gate  <- authorize_outcome_analysis(audit)
#' print(gate)
#'
#' @export
authorize_outcome_analysis <- function(audit,
                                       required_stages = NULL,
                                       allow_flag      = TRUE,
                                       lock_hash       = NULL) {
  if (!inherits(audit, "cleantmle_audit"))
    stop("`audit` must be a cleantmle_audit object.", call. = FALSE)

  if (is.null(required_stages))
    required_stages <- c("Check Point 1", "Check Point 2", "Check Point 3")

  if (is.null(lock_hash))
    lock_hash <- if (!is.null(audit$lock_hash)) audit$lock_hash else NA_character_

  # For each required stage, find matching audit entries (partial match)
  find_decision <- function(stage_pattern) {
    matched <- Filter(
      function(e) grepl(stage_pattern, e$stage, fixed = TRUE),
      audit$entries
    )
    if (length(matched) == 0L) return(NULL)
    # Return the decision from the last matching entry
    matched[[length(matched)]]$decision
  }

  decisions      <- lapply(required_stages, find_decision)
  names(decisions) <- required_stages

  missing_stages <- required_stages[vapply(decisions, is.null, logical(1L))]
  found_decisions <- unlist(decisions[!vapply(decisions, is.null, logical(1L))])

  overall_decision <- "GO"
  rationale_parts  <- character(0)

  if (length(missing_stages) > 0L) {
    overall_decision <- "STOP"
    rationale_parts  <- c(rationale_parts,
      paste("Missing required stages:", paste(missing_stages, collapse = "; ")))
  }

  if (any(found_decisions == "STOP")) {
    overall_decision <- "STOP"
    stop_stages <- names(found_decisions)[found_decisions == "STOP"]
    rationale_parts <- c(rationale_parts,
      paste("STOP decision in:", paste(stop_stages, collapse = "; ")))
  }

  if (!allow_flag && any(found_decisions == "FLAG")) {
    overall_decision <- "STOP"
    flag_stages <- names(found_decisions)[found_decisions == "FLAG"]
    rationale_parts <- c(rationale_parts,
      paste("FLAG decision not permitted in:", paste(flag_stages, collapse = "; ")))
  }

  if (overall_decision == "GO" && any(found_decisions == "FLAG")) {
    overall_decision <- "FLAG"
    flag_stages <- names(found_decisions)[found_decisions == "FLAG"]
    rationale_parts <- c(rationale_parts,
      paste("FLAG in:", paste(flag_stages, collapse = "; "), "(allowed)"))
  }

  rationale <- if (length(rationale_parts) == 0L) {
    "All required stages passed; outcome analysis authorised."
  } else {
    paste(rationale_parts, collapse = " | ")
  }

  metrics <- data.frame(
    stage             = required_stages,
    decision_in_audit = vapply(decisions, function(d) if (is.null(d)) "MISSING" else d,
                               character(1L)),
    stringsAsFactors  = FALSE
  )

  cp <- new_checkpoint(
    stage      = "Pre-Outcome Gate",
    decision   = overall_decision,
    metrics    = metrics,
    thresholds = list(allow_flag = allow_flag,
                      required_stages = required_stages),
    rationale  = rationale,
    lock_hash  = lock_hash
  )

  class(cp)        <- c("pre_outcome_gate", "cleantmle_checkpoint")
  cp$authorized    <- (cp$decision != "STOP")
  cp
}


#' Assert That Outcome Analysis Is Authorised
#'
#' Calls \code{\link{authorize_outcome_analysis}} and raises an error if
#' the outcome analysis has not been authorised.
#'
#' @section Clean-room stage: Pre-Outcome Gate assertion.
#'
#' @param audit A \code{cleantmle_audit}.
#' @param ... Additional arguments passed to
#'   \code{\link{authorize_outcome_analysis}}.
#'
#' @return \code{invisible(TRUE)} if authorised.
#'
#' @examples
#' dat   <- sim_func1(n = 500, seed = 1)
#' lock  <- create_analysis_lock(dat, "treatment", "event_24",
#'                               c("age", "sex", "biomarker"), seed = 1)
#' audit <- create_audit_log(lock)
#' ps    <- fit_ps_glm(lock)
#' cp1   <- checkpoint_cohort_adequacy(lock)
#' cp2   <- checkpoint_balance(lock, ps)
#' lock  <- define_negative_control(lock, "nc_outcome")
#' nc    <- run_negative_control(lock, "nc_outcome", ps)
#' cp3   <- checkpoint_residual_bias(nc, lock_hash = lock$lock_hash)
#' audit <- record_checkpoint(audit, cp1)
#' audit <- record_checkpoint(audit, cp2)
#' audit <- record_checkpoint(audit, cp3)
#' assert_outcome_authorized(audit)
#'
#' @export
assert_outcome_authorized <- function(audit, ...) {
  gate <- authorize_outcome_analysis(audit, ...)
  if (!isTRUE(gate$authorized)) {
    stop(
      "Outcome analysis NOT authorised. Pre-Outcome Gate decision: ",
      gate$decision, ". Rationale: ", gate$rationale,
      call. = FALSE
    )
  }
  invisible(TRUE)
}


# ── Decision Log ──────────────────────────────────────────────────────────

#' Record a Decision Log Entry
#'
#' Appends a structured decision log entry to the audit.  Use this to
#' capture analyst decisions, protocol deviations, or model-specification
#' choices that are not automatically captured by checkpoint functions.
#'
#' @param audit A \code{cleantmle_audit}.
#' @param stage Character; stage label associated with the decision.
#' @param decision_type Character; free-text category such as
#'   \code{"model_specification"}, \code{"protocol_deviation"}, or
#'   \code{"override"}.
#' @param description Character; brief description of the decision.
#' @param rationale Character; optional justification.  Default \code{""}.
#' @param metrics A named list of optional numeric or character metrics.
#'   Default \code{NULL}.
#'
#' @return Modified \code{cleantmle_audit} with the new entry appended to
#'   \code{audit$decision_log}.
#'
#' @examples
#' dat   <- sim_func1(n = 200, seed = 1)
#' lock  <- create_analysis_lock(dat, "treatment", "event_24",
#'                               c("age", "sex", "biomarker"), seed = 1)
#' audit <- create_audit_log(lock)
#' audit <- record_decision_log_entry(
#'   audit, "Stage 2", "model_specification",
#'   "Selected logistic PS model with 3 covariates",
#'   rationale = "Pre-specified in SAP"
#' )
#'
#' @export
record_decision_log_entry <- function(audit,
                                      stage,
                                      decision_type,
                                      description,
                                      rationale = "",
                                      metrics   = NULL) {
  if (!inherits(audit, "cleantmle_audit"))
    stop("`audit` must be a cleantmle_audit object.", call. = FALSE)

  entry <- list(
    stage         = stage,
    decision_type = decision_type,
    description   = description,
    rationale     = rationale,
    metrics       = metrics,
    timestamp     = Sys.time()
  )

  if (is.null(audit$decision_log))
    audit$decision_log <- list()

  audit$decision_log <- c(audit$decision_log, list(entry))
  audit
}


#' Export the Decision Log as a Data Frame
#'
#' Converts all entries in \code{audit$decision_log} to a tidy data frame.
#'
#' @param audit A \code{cleantmle_audit}.
#'
#' @return A data.frame with columns \code{stage}, \code{decision_type},
#'   \code{description}, \code{rationale}, and \code{timestamp}.  Returns
#'   an empty data frame if no decision log entries exist.
#'
#' @examples
#' dat   <- sim_func1(n = 200, seed = 1)
#' lock  <- create_analysis_lock(dat, "treatment", "event_24",
#'                               c("age", "sex", "biomarker"), seed = 1)
#' audit <- create_audit_log(lock)
#' audit <- record_decision_log_entry(audit, "Stage 2", "model_specification",
#'                                    "Logistic PS model selected")
#' export_decision_log(audit)
#'
#' @export
export_decision_log <- function(audit) {
  if (!inherits(audit, "cleantmle_audit"))
    stop("`audit` must be a cleantmle_audit object.", call. = FALSE)

  empty <- data.frame(
    stage         = character(0),
    decision_type = character(0),
    description   = character(0),
    rationale     = character(0),
    timestamp     = character(0),
    stringsAsFactors = FALSE
  )

  if (is.null(audit$decision_log) || length(audit$decision_log) == 0L)
    return(empty)

  rows <- lapply(audit$decision_log, function(e) {
    data.frame(
      stage         = e$stage,
      decision_type = e$decision_type,
      description   = e$description,
      rationale     = e$rationale,
      timestamp     = format(e$timestamp, "%Y-%m-%d %H:%M:%S"),
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, rows)
}


# ── Stage Path Narrative ──────────────────────────────────────────────────

#' Summarise the Stage Path as a Compact Narrative
#'
#' Generates a compact, human-readable narrative string showing the
#' sequence of stages and checkpoint decisions recorded in the audit log.
#'
#' @param audit A \code{cleantmle_audit}.
#'
#' @return The narrative string, returned \code{invisible}ly after
#'   printing.
#'
#' @examples
#' dat   <- sim_func1(n = 200, seed = 1)
#' lock  <- create_analysis_lock(dat, "treatment", "event_24",
#'                               c("age", "sex", "biomarker"), seed = 1)
#' audit <- create_audit_log(lock)
#' audit <- record_stage(audit, "Stage 1a", "Lock created", decision = NA)
#' audit <- record_stage(audit, "Check Point 1", "CP1 evaluated", decision = "GO")
#' summarize_stage_path(audit)
#'
#' @export
summarize_stage_path <- function(audit) {
  if (!inherits(audit, "cleantmle_audit"))
    stop("`audit` must be a cleantmle_audit object.", call. = FALSE)

  if (length(audit$entries) == 0L) {
    msg <- "(No stages recorded)"
    cat(msg, "\n")
    return(invisible(msg))
  }

  parts <- vapply(audit$entries, function(e) {
    if (!is.na(e$decision) && nzchar(e$decision)) {
      sprintf("%s: %s", e$stage, e$decision)
    } else {
      e$stage
    }
  }, character(1L))

  narrative <- paste(parts, collapse = " -> ")
  cat(narrative, "\n")
  invisible(narrative)
}


# ── Outcome Masking ───────────────────────────────────────────────────────

#' Mask the Outcome Column in the Lock
#'
#' Returns a modified copy of the lock in which the outcome column is
#' replaced with \code{NA}.  This enforces the clean-room principle of
#' keeping analysts blind to outcome data during the design stage.
#'
#' @section Clean-room stage: Design stage (pre-outcome blinding).
#'
#' @param lock A \code{cleanroom_lock}.
#'
#' @return A modified \code{cleanroom_lock} with \code{lock$data[[outcome]]}
#'   set to \code{NA} and \code{lock$.outcome_masked} set to \code{TRUE}.
#'
#' @examples
#' dat    <- sim_func1(n = 200, seed = 1)
#' lock   <- create_analysis_lock(dat, "treatment", "event_24",
#'                                c("age", "sex", "biomarker"), seed = 1)
#' masked <- mask_outcome(lock)
#' all(is.na(masked$data[[masked$outcome]]))  # TRUE
#'
#' @export
mask_outcome <- function(lock) {
  if (!inherits(lock, "cleanroom_lock"))
    stop("`lock` must be a cleanroom_lock object.", call. = FALSE)

  lock$data[[lock$outcome]] <- NA
  lock$.outcome_masked      <- TRUE
  lock
}


#' Restore the Outcome Column from an Original Lock
#'
#' Copies the outcome column from \code{original_lock} back into
#' \code{lock}, reversing the effect of \code{\link{mask_outcome}}.
#'
#' @section Clean-room stage: Pre-estimation (outcome unblinding).
#'
#' @param lock A \code{cleanroom_lock} whose outcome column is masked.
#' @param original_lock The original \code{cleanroom_lock} holding the
#'   true outcome values.
#'
#' @return A modified \code{cleanroom_lock} with the outcome column
#'   restored and \code{lock$.outcome_masked} set to \code{FALSE}.
#'
#' @examples
#' dat          <- sim_func1(n = 200, seed = 1)
#' lock         <- create_analysis_lock(dat, "treatment", "event_24",
#'                                      c("age", "sex", "biomarker"), seed = 1)
#' masked       <- mask_outcome(lock)
#' unmasked     <- unmask_outcome(masked, lock)
#' identical(unmasked$data[["event_24"]], lock$data[["event_24"]])  # TRUE
#'
#' @export
unmask_outcome <- function(lock, original_lock) {
  if (!inherits(lock, "cleanroom_lock"))
    stop("`lock` must be a cleanroom_lock object.", call. = FALSE)
  if (!inherits(original_lock, "cleanroom_lock"))
    stop("`original_lock` must be a cleanroom_lock object.", call. = FALSE)

  outcome_col <- lock$outcome
  if (!outcome_col %in% names(original_lock$data))
    stop("Outcome column '", outcome_col,
         "' not found in original_lock$data.", call. = FALSE)

  lock$data[[outcome_col]] <- original_lock$data[[outcome_col]]
  lock$.outcome_masked     <- FALSE
  lock
}
