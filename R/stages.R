

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
#' @param domain Character; optional Muntner et al. (2024) negative-control
#'   domain. One of \code{"confounding_by_indication"},
#'   \code{"functional_status"}, \code{"health_seeking_behavior"},
#'   \code{"access_to_healthcare"}, or \code{"other"}. Recording the
#'   domain helps reviewers see that the registered NCs span the four
#'   confounding domains recommended in the staging-and-clean-room paper.
#'
#' @return Modified \code{cleanroom_lock} with element
#'   \code{negative_controls}.
#'
#' @examples
#' \dontrun{
#' dat  <- sim_func1(n = 200, seed = 1)
#' lock <- create_analysis_lock(dat, "treatment", "event_24",
#'                              c("age", "sex", "biomarker"), seed = 1)
#' lock <- define_negative_control(lock, "nc_outcome",
#'   description = "Outcome known to be unrelated to treatment",
#'   domain = "confounding_by_indication")
#'
#' }
#' @keywords internal
define_negative_control <- function(lock, variable, type = "outcome",
                                    description = NULL, domain = NULL) {
  if (!inherits(lock, "cleanroom_lock"))
    stop("`lock` must be a cleanroom_lock object.", call. = FALSE)
  type <- match.arg(type, c("outcome", "exposure"))
  if (!variable %in% names(lock$data))
    stop("Variable '", variable, "' not found in lock data.", call. = FALSE)
  if (!is.null(domain)) {
    valid_domains <- c("confounding_by_indication", "functional_status",
                       "health_seeking_behavior", "access_to_healthcare",
                       "other")
    domain <- match.arg(domain, valid_domains)
  }

  nc <- list(
    variable    = variable,
    type        = type,
    description = description,
    domain      = domain
  )

  if (is.null(lock$negative_controls)) lock$negative_controls <- list()
  lock$negative_controls[[variable]] <- nc
  lock
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
#' @section Lock-hash impact:
#' Attaching the primary TMLE spec extends \code{lock$primary_tmle_spec}
#' but does not recompute the original \code{lock_hash} captured by
#' \code{\link{create_analysis_lock}}. The hash continues to fingerprint
#' the prespecified analytic plan (treatment, outcome, covariates,
#' SuperLearner library, seed, data shape). The selected candidate is
#' recorded as a downstream choice made under the locked candidate grid
#' and selection rule; switching the spec is logged in the audit trail
#' rather than by re-hashing the lock.
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
#' \dontrun{
#' dat  <- sim_func1(n = 200, seed = 1)
#' lock <- create_analysis_lock(dat, "treatment", "event_24",
#'                              c("age", "sex", "biomarker"), seed = 1)
#' spec <- tmle_candidate("glm_t01", "GLM trunc=0.01",
#'                        g_library = "SL.glm", truncation = 0.01)
#' lock <- lock_primary_tmle_spec(lock, spec)
#' get_primary_tmle_spec(lock)
#'
#' }
#' @keywords internal
lock_primary_tmle_spec <- function(lock, selected) {
  if (!inherits(lock, "cleanroom_lock"))
    stop("`lock` must be a cleanroom_lock object.", call. = FALSE)
  # A single-candidate set from define_candidates() stands in for its
  # one specification.
  if (inherits(selected, "ct_candidates") && length(selected) == 1L)
    selected <- selected[[1L]]
  if (!inherits(selected, "tmle_candidate_spec"))
    stop("`selected` must be a tmle_candidate_spec (or a length-one ",
         "candidate set from define_candidates()).", call. = FALSE)

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
#' @keywords internal
get_primary_tmle_spec <- function(lock) {
  if (!inherits(lock, "cleanroom_lock"))
    stop("`lock` must be a cleanroom_lock object.", call. = FALSE)
  lock$primary_tmle_spec
}


# ── Design Precision ──────────────────────────────────────────────────────

#' Estimate Design-Stage Precision Summaries
#'
#' Computes pre-outcome precision summaries from arm sizes and marginal
#' outcome information only. No outcome modelling is performed, and
#' nothing in the output stratifies the outcome by treatment arm: total
#' events, the marginal event rate, and the precision proxies derived
#' from them are the only outcome quantities read. Arm-specific event
#' counts are the crude treatment-outcome association; they are
#' available only through [event_support_by_arm()], which requires a
#' recorded reason.
#'
#' @section Clean-room stage: Stage 1b (design-stage; marginal outcome
#'   counts only, never the treatment-outcome association).
#'
#' @param lock A \code{cleanroom_lock}.
#' @param target_mdd Numeric; an optional target minimum detectable
#'   difference on the risk-difference scale to compare against the
#'   computed MDD.  Default \code{NULL} (no comparison).
#'
#' @return A list of class `design_precision` with elements: `n_total`
#'   (total sample size); `n_treated`; `n_control` (arm sizes are
#'   treatment-marginal design quantities); `events_total` (marginal
#'   event count); `prevalence` (marginal outcome prevalence among rows
#'   with an observed outcome); `se_proxy` (SE proxy for the risk
#'   difference under the marginal rate); `ci_halfwidth` (95% CI
#'   half-width proxy, 1.96 * se_proxy); `mdd_80` (minimum detectable
#'   difference at 80% power, 2.8 * se_proxy); `target_mdd` (the
#'   user-supplied target, or `NULL`); `mdd_feasible` (logical; present
#'   only when `target_mdd` is supplied); and `event_support` (the
#'   marginal table from [summarize_event_support()]).
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
  Y         <- .outcome_vector(lock)

  if (isTRUE(lock$.outcome_masked) || !.outcome_readable(lock)) {
    stop("estimate_design_precision() requires a readable outcome. ",
         "The lock is masked or the outcome is all NA. Call this on the ",
         "un-masked lock (before mask_outcome()) or call unmask_outcome() ",
         "first.", call. = FALSE)
  }

  n_total   <- nrow(data)
  n_treated <- sum(A == 1L, na.rm = TRUE)
  n_control <- sum(A == 0L, na.rm = TRUE)

  # Marginal outcome quantities only: total events and the marginal rate.
  # No quantity below conditions the outcome on the treatment arm.
  events_total <- sum(Y, na.rm = TRUE)
  n_obs_y      <- sum(!is.na(Y))

  prevalence  <- if (n_obs_y > 0) events_total / n_obs_y else NA_real_
  p           <- prevalence
  se_proxy    <- if (!is.na(p) && n_treated > 0 && n_control > 0)
                   sqrt(p * (1 - p) * (1 / n_treated + 1 / n_control))
                 else NA_real_
  ci_halfwidth <- 1.96 * se_proxy
  mdd_80       <- 2.8  * se_proxy

  result <- list(
    n_total        = n_total,
    n_treated      = n_treated,
    n_control      = n_control,
    events_total   = events_total,
    prevalence     = prevalence,
    se_proxy       = se_proxy,
    ci_halfwidth   = ci_halfwidth,
    mdd_80         = mdd_80,
    target_mdd     = target_mdd
  )

  if (!is.null(target_mdd)) {
    result$mdd_feasible <- mdd_80 <= target_mdd
  }

  # The marginal event-support table rides on the precision object.
  result$event_support <- tryCatch(summarize_event_support(lock),
                                   error = function(e) NULL)

  class(result) <- "design_precision"
  result
}


#' @export
print.design_precision <- function(x, ...) {
  cat("=== Design-Stage Precision Summary (marginal outcome only) ===\n")
  tbl <- data.frame(
    Metric = c(
      "N total", "N treated", "N control",
      "Events total (marginal)",
      "Marginal prevalence",
      "SE proxy (RD)", "95% CI half-width", "MDD (80% power)"
    ),
    Value = c(
      x$n_total, x$n_treated, x$n_control,
      x$events_total,
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
  cat("Arm-specific event counts require event_support_by_arm(lock, reason = ...).\n")
  invisible(x)
}


# ── Event Support Summary ─────────────────────────────────────────────────

#' Summarise Marginal Event Support
#'
#' Returns a one-row data frame with the total sample size, the marginal
#' event count, and the marginal event rate. Nothing here conditions the
#' outcome on the treatment arm: arm-specific counts are the crude
#' treatment-outcome association and are available only through
#' [event_support_by_arm()], which requires a recorded reason.
#'
#' @section Clean-room stage: Stage 1b (design diagnostics; marginal
#'   outcome only).
#'
#' @param lock A \code{cleanroom_lock}.
#'
#' @return A one-row data.frame with columns \code{n}, \code{events},
#'   and \code{event_rate} (all marginal).
#'
#' @examples
#' \dontrun{
#' dat  <- sim_func1(n = 500, seed = 1)
#' lock <- create_analysis_lock(dat, "treatment", "event_24",
#'                              c("age", "sex", "biomarker"), seed = 1)
#' summarize_event_support(lock)
#'
#' }
#' @keywords internal
summarize_event_support <- function(lock) {
  if (!inherits(lock, "cleanroom_lock"))
    stop("`lock` must be a cleanroom_lock object.", call. = FALSE)

  data <- lock$data
  Y    <- .outcome_vector(lock)

  if (isTRUE(lock$.outcome_masked) || !.outcome_readable(lock)) {
    stop("summarize_event_support() requires a readable outcome. ",
         "The lock is masked or the outcome is all NA.",
         call. = FALSE)
  }

  n_obs_y <- sum(!is.na(Y))
  events  <- sum(Y, na.rm = TRUE)

  out <- data.frame(
    n          = nrow(data),
    events     = events,
    event_rate = round(if (n_obs_y > 0) events / n_obs_y else NA_real_, 4),
    stringsAsFactors = FALSE
  )

  if (events < 10L) {
    message("Sparse-event warning: fewer than 10 events in total. ",
            "Estimates may be unstable.")
  }

  out
}


#' Event Counts by Treatment Arm (Logged Access)
#'
#' Tabulates event counts and rates per treatment arm. Arm-specific
#' event counts are the crude treatment-outcome association, so this
#' access path is deliberately separate from the marginal design
#' diagnostics: it requires a stated `reason`, emits a warning, and
#' writes a design-log entry recording that arm-specific counts were
#' viewed. The returned object carries both the table and the updated
#' lock; keep the lock so the access stays on the record:
#' `esa <- event_support_by_arm(lock, reason = "..."); lock <- esa$lock`.
#'
#' @section Clean-room stage: outside the outcome-blind design path; the
#'   access itself is what the design log records.
#'
#' @param lock A \code{cleanroom_lock} with a readable outcome column.
#' @param reason Character; why arm-specific counts are needed (for
#'   example a data-manager review of sparse cells). Mandatory.
#'
#' @return An object of class `event_support_by_arm`: a list with
#'   `table` (columns `arm`, `n`, `events`, `event_rate`), `reason`, and
#'   `lock` (the lock with the access recorded in its design log).
#'
#' @examples
#' dat  <- sim_func1(n = 500, seed = 1)
#' lock <- create_analysis_lock(dat, "treatment", "event_24",
#'                              c("age", "sex", "biomarker"), seed = 1)
#' esa  <- event_support_by_arm(lock,
#'   reason = "data-manager check of sparse cells before locking")
#' esa$table
#' lock <- esa$lock  # keep the logged access on the lock
#'
#' @export
event_support_by_arm <- function(lock, reason) {
  if (!inherits(lock, "cleanroom_lock"))
    stop("`lock` must be a cleanroom_lock object.", call. = FALSE)
  if (missing(reason) || !is.character(reason) || length(reason) != 1L ||
      !nzchar(trimws(reason)))
    stop("event_support_by_arm() requires a non-empty `reason`: ",
         "arm-specific event counts reveal the crude treatment-outcome ",
         "association, and the access is recorded in the design log.",
         call. = FALSE)

  data <- lock$data
  A    <- data[[lock$treatment]]
  Y    <- .outcome_vector(lock)
  if (isTRUE(lock$.outcome_masked) || !.outcome_readable(lock))
    stop("event_support_by_arm() requires a readable outcome. ",
         "The lock is masked or the outcome is all NA.",
         call. = FALSE)

  strat <- lock$estimand$treatment_strategies
  treated_lbl <- if (!is.null(strat) && length(strat) >= 1) strat[1] else "Treated"
  control_lbl <- if (!is.null(strat) && length(strat) >= 2) strat[2] else "Control"

  n1 <- sum(A == 1L, na.rm = TRUE)
  n0 <- sum(A == 0L, na.rm = TRUE)
  e1 <- sum(Y[A == 1L], na.rm = TRUE)
  e0 <- sum(Y[A == 0L], na.rm = TRUE)

  tab <- data.frame(
    arm        = c(treated_lbl, control_lbl, "Total"),
    n          = c(n1, n0, n1 + n0),
    events     = c(e1, e0, e1 + e0),
    event_rate = round(c(if (n1 > 0) e1 / n1 else NA_real_,
                         if (n0 > 0) e0 / n0 else NA_real_,
                         if (n1 + n0 > 0) (e1 + e0) / (n1 + n0) else NA_real_),
                       4),
    stringsAsFactors = FALSE
  )

  warning("Arm-specific event counts reveal the crude treatment-outcome ",
          "association; this access has been recorded in the design log.",
          call. = FALSE)

  sparse_arms <- tab$arm[tab$arm != "Total" & tab$events < 10L]
  if (length(sparse_arms) > 0L)
    message("Sparse-event warning: fewer than 10 events in arm(s): ",
            paste(sparse_arms, collapse = ", "), ".")

  lock <- .log_design_decision(
    lock, "event_support_by_arm",
    sprintf("Arm-specific event counts viewed. Reason: %s", reason))

  out <- list(table = tab, reason = reason, lock = lock)
  class(out) <- "event_support_by_arm"
  out
}

#' @export
print.event_support_by_arm <- function(x, ...) {
  cat("Event support by arm (logged access)\n")
  print(x$table, row.names = FALSE)
  cat("Reason recorded in the design log:", x$reason, "\n")
  invisible(x)
}


# ── Outcome Masking ───────────────────────────────────────────────────────

#' Mask the Outcome
#'
#' Returns a modified copy of the lock with the outcome physically
#' absent: on a store-format lock (0.3.0 and later) the outcome store is
#' removed, so no outcome values exist anywhere on the returned object;
#' on a legacy lock the outcome column is set to \code{NA}. Keep an
#' unmasked lock (or the raw data) to unmask later with
#' \code{\link{unmask_outcome}}.
#'
#' @section Clean-room stage: Design stage (pre-outcome blinding).
#'
#' @param lock A \code{cleanroom_lock}.
#'
#' @return A modified \code{cleanroom_lock} with the outcome removed and
#'   \code{lock$.outcome_masked} set to \code{TRUE}.
#'
#' @examples
#' dat    <- sim_func1(n = 200, seed = 1)
#' lock   <- create_analysis_lock(dat, "treatment", "event_24",
#'                                c("age", "sex", "biomarker"), seed = 1)
#' masked <- mask_outcome(lock)
#' is.null(masked$outcome_store)  # TRUE: physically outcome-free
#'
#' @export
mask_outcome <- function(lock) {
  if (!inherits(lock, "cleanroom_lock"))
    stop("`lock` must be a cleanroom_lock object.", call. = FALSE)

  if (.lock_has_store_format(lock)) {
    lock$outcome_store <- NULL
  } else {
    lock$data[[lock$outcome]] <- NA
  }
  lock$.outcome_masked <- TRUE
  lock
}


#' Restore the Outcome from an Unmasked Original
#'
#' Restores the outcome store (or, on a legacy lock, the outcome column)
#' from \code{original_lock}, reversing \code{\link{mask_outcome}}, and
#' records the unmasking in the design log. On a lock created with
#' \code{enforce = TRUE}, unmasking requires a named \code{approved_by};
#' the name is written into the design log's \code{decided_by} column
#' and the lock is stamped as authorised for estimation. This explicit,
#' logged step is the single institutional switch: there is no token and
#' no gate object, only a person on the record.
#'
#' @section Clean-room stage: Pre-estimation (outcome unblinding).
#'
#' @param lock A \code{cleanroom_lock}, normally masked.
#' @param original_lock An unmasked \code{cleanroom_lock} over the same
#'   data, holding the true outcome values.
#' @param approved_by Character; the person or role authorising the
#'   unmasking. Required when the lock was created with
#'   \code{enforce = TRUE}; recorded in the design log whenever given.
#' @param allow_unauthorized Deprecated (the audit gate it bypassed was
#'   removed in 0.3.0). \code{TRUE} is accepted with a warning and
#'   behaves as an unnamed forced approval, so old scripts keep running.
#'
#' @return A modified \code{cleanroom_lock} with the outcome restored,
#'   \code{lock$.outcome_masked} set to \code{FALSE},
#'   \code{lock$.outcome_authorized} set to \code{TRUE}, and the
#'   unmasking recorded in \code{lock$design_log}.
#'
#' @examples
#' dat    <- sim_func1(n = 200, seed = 1)
#' lock   <- create_analysis_lock(dat, "treatment", "event_24",
#'                                c("age", "sex", "biomarker"), seed = 1)
#' masked <- mask_outcome(lock)
#' unmasked <- unmask_outcome(masked, lock, approved_by = "review team")
#' identical(cleanTMLE:::.outcome_vector(unmasked),
#'           cleanTMLE:::.outcome_vector(lock))  # TRUE
#'
#' @export
unmask_outcome <- function(lock, original_lock, approved_by = NULL,
                           allow_unauthorized = NULL) {
  if (!inherits(lock, "cleanroom_lock"))
    stop("`lock` must be a cleanroom_lock object.", call. = FALSE)
  if (!inherits(original_lock, "cleanroom_lock"))
    stop("`original_lock` must be a cleanroom_lock object.", call. = FALSE)

  if (!is.null(allow_unauthorized)) {
    rlang::warn(paste(
      "unmask_outcome(allow_unauthorized = ) is deprecated: the audit gate",
      "it bypassed was removed in 0.3.0. Unmasking is now recorded in the",
      "design log; pass approved_by = to name the approver."),
      .frequency = "once", .frequency_id = "unmask_allow_unauthorized")
    if (isTRUE(allow_unauthorized) && is.null(approved_by))
      approved_by <- "(forced; legacy allow_unauthorized argument)"
  }

  y <- .outcome_vector(original_lock)
  if (is.null(y) || all(is.na(y)))
    stop("`original_lock` holds no readable outcome '", lock$outcome,
         "'; unmask from the lock created over the full data.",
         call. = FALSE)
  if (length(y) != nrow(lock$data))
    stop("`original_lock` outcome has ", length(y), " rows; the lock ",
         "data has ", nrow(lock$data), ".", call. = FALSE)

  if (isTRUE(lock$cleanroom_enabled) &&
      isTRUE(lock$require_authorization) &&
      (is.null(approved_by) || !nzchar(trimws(approved_by)))) {
    stop("unmask_outcome(): this lock was created with enforce = TRUE; ",
         "unmasking requires a named `approved_by`, which is written ",
         "into the design log.", call. = FALSE)
  }

  if (.lock_has_store_format(lock)) {
    lock$outcome_store <- .new_outcome_store(y, lock$outcome,
                                             lock$lock_hash)
  } else {
    lock$data[[lock$outcome]] <- y
  }
  lock$.outcome_masked     <- FALSE
  lock$.outcome_authorized <- TRUE
  lock <- .log_design_decision(
    lock, "outcome_unmasked",
    sprintf("Outcome '%s' unmasked for estimation.", lock$outcome),
    stage = "Stage 4 (unmasking)",
    decision = "unmask and authorise estimation",
    decided_by = approved_by %||% NA_character_)
  lock
}
