#' Analysis Lock and Staged Workflow Functions
#'
#' These functions implement the staged clean-room workflow: locking the
#' analytic specification, estimating propensity scores with SuperLearner,
#' running plasmode feasibility evaluation, and executing the final
#' conventional (matching, IPTW) and modular TMLE workflows.
#'
#' @name cleanroom_workflow
NULL


# ── Internal: Clean-Room Outcome Guard ──────────────────────────────────

#' Check whether the lock allows outcome access
#'
#' @description Internal helper used by Stage 4 functions to operationalise
#'   the clean-room gate.  Checks for outcome masking and optionally
#'   for an \code{audit} attribute carrying authorisation status.
#'
#' @param lock A \code{cleanroom_lock}.
#' @param allow_outcome_access Logical; if \code{TRUE} the outcome-guard check
#'   is skipped. Use this argument in place of the deprecated
#'   \code{override_clean_room}.
#' @param caller Character; name of the calling function (for error messages).
#'
#' @return \code{invisible(TRUE)} if access is permitted.
#' @keywords internal
.check_outcome_access <- function(lock, allow_outcome_access = FALSE,
                                  caller = "Stage 4 function") {
  if (isTRUE(allow_outcome_access)) return(invisible(TRUE))

  # Lock-level opt-out: if the user created the lock with
  # cleanroom_enabled = FALSE, the package is being used as a plain
  # outcome-blind RWE pipeline without the audit/hash/gate machinery.
  # Stage 4 functions still run; outcome-blindness is the analyst's
  # responsibility rather than software-enforced.
  if (isFALSE(lock$cleanroom_enabled)) return(invisible(TRUE))

  # Check outcome masking
  if (isTRUE(lock$.outcome_masked)) {
    stop(
      caller, ": Outcome is masked. ",
      "Use unmask_outcome() first, set allow_outcome_access = TRUE, ",
      "or create the lock with cleanroom_enabled = FALSE to skip ",
      "the software-enforced outcome guard entirely.",
      call. = FALSE
    )
  }

  # Authorisation check, opt-in via create_analysis_lock(enforce = TRUE).
  # The default lock enforces masking only: the one honest blinding
  # device is a physically absent outcome, and what protects an analysis
  # beyond that is statistical (the support verdict, the estimand
  # ladder, the implausibility guard), not a token. The enforce switch
  # adds one requirement: estimation refuses to run until the outcome
  # was unmasked through unmask_outcome() with a named approver, which
  # writes the authorisation into the design log.
  if (isTRUE(lock$require_authorization) &&
      !isTRUE(lock$.outcome_authorized)) {
    stop(
      caller, ": this lock was created with enforce = TRUE and has not ",
      "been authorised for estimation. Unmask it through ",
      "unmask_outcome(lock, original_lock, approved_by = <name>) so the ",
      "authorisation is on the design log; or set ",
      "allow_outcome_access = TRUE to override this call.",
      call. = FALSE
    )
  }

  invisible(TRUE)
}


# ── Stage 1: Analysis Lock ────────────────────────────────────────────────

#' Create an Analysis Lock
#'
#' Captures the full analytic specification (data, variable names, SuperLearner
#' library, and plasmode simulation settings) as a locked object. The lock
#' is created before the outcome is examined to operationalise the
#' pre-specified analytic plan across all subsequent stages.
#'
#' @section What is fingerprinted in the lock hash:
#' The \code{lock_hash} is a sha256 digest (via \pkg{digest}; with a
#' clearly-labelled non-cryptographic checksum fallback) computed over the
#' treatment name, outcome name, sorted covariate vector, SuperLearner
#' library, integer seed, the row count, column count, and sorted column
#' names of \code{data}, a content digest of every column of \code{data}
#' except the outcome column (so swapping design-data values is
#' detectable, while masking or unmasking the outcome leaves the
#' fingerprint unchanged and the fingerprint commits to nothing about
#' the outcome), the plasmode generator mode (\code{dgp_mode}), and the
#' declared \code{nc_criteria} and \code{dq_thresholds}. Changing any of
#' these invalidates the hash.
#'
#' @section What does NOT invalidate the lock:
#' Design-log entries and user-supplied notes accumulate on the lock for
#' traceability and do not feed into the lock hash; the fingerprint
#' covers the specification, not the running record.
#'
#' @section The outcome store:
#' Since 0.3.0 the lock's data frame holds design data only (covariates,
#' treatment, missingness indicators, negative controls); the primary
#' outcome column never sits in \code{lock$data}. The outcome lives in
#' \code{lock$outcome_store}, a sealed object that only the Stage 4
#' estimators join back, after the outcome guard.
#' \code{\link{mask_outcome}} removes the store from its copy of the
#' lock, so a masked lock is physically outcome-free, and
#' \code{\link{unmask_outcome}} restores it from an unmasked original
#' with a design-log entry recording the unmasking (on an
#' \code{enforce = TRUE} lock, with a named approver).
#'
#' @section Verdicts, not tokens:
#' The design stage records graded verdicts on its result objects: the
#' support verdict of \code{\link{assess_support}}, the per-estimand
#' verdicts of \code{\link{estimand_feasibility}}, the locked
#' data-quality verdict and tipping points of \code{\link{stress_test}},
#' and the Check Point 3 reading of \code{\link{negative_control_ladder}}
#' under the locked \code{nc_criteria}. \code{\link{design_report}}
#' assembles them for the review team; estimand switches and overrides
#' are logged design decisions in \code{lock$design_log}, exported with
#' \code{\link{export_design_log}}.
#'
#' @param data A data.frame containing covariates, treatment, and outcome.
#' @param treatment Character; name of the binary treatment column.
#' @param outcome Character; name of the outcome column.
#' @param covariates Character vector of baseline covariate column names.
#' @param sl_library Character vector of SuperLearner algorithm names.
#'   Default: `c("SL.glm", "SL.mean")`.
#' @param plasmode_reps Integer; number of plasmode replicates for Stage 2b
#'   feasibility evaluation. Default: 100.
#' @param seed Integer; random seed for reproducibility. Default: 42.
#' @param cleanroom_enabled Logical; if `TRUE` (default) the outcome-masking
#'   guard applies to this lock. Since 0.2.0 this is the only default
#'   enforcement; see `enforce`.
#' @param negative_controls Optional character vector of negative-control
#'   columns to register on the lock at creation (before any downstream
#'   filtering can touch them).
#' @param mask Logical; if `TRUE`, the outcome store is removed from the
#'   returned lock, so it is physically outcome-free. Keep an unmasked
#'   lock (or create one later over the same data) to unmask with
#'   [unmask_outcome()]. Default `FALSE`.
#' @param enforce Logical; if `TRUE`, Stage 4 estimators refuse to run
#'   until the outcome has been unmasked through [unmask_outcome()] with
#'   a named `approved_by`, which writes the authorisation into the
#'   design log. This is the single institutional switch. Default
#'   `FALSE`.
#' @param dgp_mode Character; the plasmode outcome-generator mode, a
#'   locked field read by [run_plasmode_feasibility()] and
#'   [run_plasmode_dq_stress()]. `"hybrid"` (default) fits the
#'   generator's baseline outcome model Q0(W) on the lock's real outcome
#'   (covariate-only; the treatment-outcome association is never used);
#'   `"external_pilot"` requires the caller to supply Q0 predictions
#'   fitted on external pilot data, so the plasmode stage reads no real
#'   outcome at all and runs on a masked lock.
#' @param nc_criteria Optional prespecified negative-control decision
#'   criteria, declared before outcome access and fingerprinted: a list
#'   with `null_band` (positive half-width on the `scale`), and
#'   optionally `scale` (`"rd"`, the default, or `"ratio"`), `rule`
#'   (`"point_in_band"`, the default, or `"ci_in_band"`),
#'   `min_per_domain` (default 1), and `consistency`
#'   (`"all_within_band"`, the default, or `"majority_within_band"`).
#'   Read by [run_negative_control_ladder()], which refuses to grade
#'   controls without them.
#' @param dq_thresholds Optional prespecified data-quality decision
#'   thresholds, declared before outcome access and fingerprinted: a
#'   list with `max_abs_bias`, `min_coverage`, `max_rmse_ratio`, and
#'   optionally `flag_coverage` and `flag_rmse_ratio`. Read by
#'   [run_plasmode_dq_stress()], whose verdict and tipping points are a
#'   pure function of these thresholds and the declared threat grid.
#' @param roles Optional named list or character vector recording the
#'   personnel structure (for example `programmer`, `analyst`,
#'   `analytic_advisor`, `review_team`). Recorded and printed in the
#'   design report header; nothing is enforced.
#' @param estimand Optional structured estimand description: a named
#'   list with any of `description`, `population`,
#'   `treatment_strategies` (length-2 character contrast),
#'   `outcome_label`, `followup`, `contrast` (default
#'   `"risk_difference"`), and `statistical_estimand`; or a single
#'   character description. Descriptive metadata; not part of the
#'   fingerprint.
#' @param sensitivity_plans Optional named list of sensitivity-analysis
#'   plans; each element name is the plan label and each element may
#'   carry `description` and `settings` (a named list, for example
#'   `list(truncation = c(0.01, 0.05, 0.10))`). Descriptive metadata;
#'   not part of the fingerprint.
#'
#' @return An object of class `cleanroom_lock` containing all specified
#'   analysis parameters plus a reproducibility fingerprint (`lock_hash`).
#'
#' @examples
#' dat <- sim_func1(n = 200, seed = 1)
#' lock <- create_analysis_lock(
#'   data       = dat,
#'   treatment  = "treatment",
#'   outcome    = "event_24",
#'   covariates = c("age", "sex", "biomarker"),
#'   seed       = 1
#' )
#' print(lock)
#'
#' @export
create_analysis_lock <- function(data, treatment, outcome, covariates,
                                  sl_library        = c("SL.glm", "SL.mean"),
                                  plasmode_reps     = 100L,
                                  seed              = 42L,
                                  cleanroom_enabled = TRUE,
                                  negative_controls = NULL,
                                  mask              = FALSE,
                                  enforce           = FALSE,
                                  dgp_mode          = c("hybrid",
                                                        "external_pilot"),
                                  nc_criteria       = NULL,
                                  dq_thresholds     = NULL,
                                  roles             = NULL,
                                  estimand          = NULL,
                                  sensitivity_plans = NULL) {
  if (!is.data.frame(data))
    stop("`data` must be a data.frame.", call. = FALSE)
  if (!is.character(treatment) || length(treatment) != 1L)
    stop("`treatment` must be a single character string.", call. = FALSE)
  if (!is.character(outcome) || length(outcome) != 1L)
    stop("`outcome` must be a single character string.", call. = FALSE)
  if (!treatment %in% names(data))
    stop("treatment variable '", treatment, "' not found in data.",
         call. = FALSE)
  if (!outcome %in% names(data))
    stop("outcome variable '", outcome, "' not found in data.",
         call. = FALSE)
  missing_cov <- covariates[!covariates %in% names(data)]
  if (length(missing_cov) > 0L)
    stop("covariates not found in data: ",
         paste(missing_cov, collapse = ", "), call. = FALSE)
  dgp_mode      <- match.arg(dgp_mode)
  nc_criteria   <- .normalize_nc_criteria(nc_criteria)
  dq_thresholds <- .normalize_dq_thresholds(dq_thresholds)
  if (!is.null(roles) && is.null(names(roles)))
    stop("`roles` must be a named list or character vector.", call. = FALSE)

  lock_hash <- .compute_lock_hash(list(
    treatment     = treatment,
    outcome       = outcome,
    covariates    = covariates,
    sl_library    = sl_library,
    seed          = as.integer(seed),
    data_nrow     = nrow(data),
    data_ncol     = ncol(data),
    data_names    = paste(sort(names(data)), collapse = "|"),
    data_content  = .design_data_digest(data, outcome),
    dgp_mode      = dgp_mode,
    nc_criteria   = nc_criteria,
    dq_thresholds = dq_thresholds
  ))

  # The outcome-store split: the lock's data frame carries design data
  # only, and the primary outcome moves into a sealed store that Stage 4
  # joins back after the outcome guard.
  design_data <- data
  design_data[[outcome]] <- NULL

  lock <- list(
    data              = design_data,
    outcome_store     = .new_outcome_store(data[[outcome]], outcome,
                                           lock_hash),
    treatment         = treatment,
    outcome           = outcome,
    covariates        = covariates,
    sl_library        = sl_library,
    plasmode_reps     = as.integer(plasmode_reps),
    seed              = as.integer(seed),
    cleanroom_enabled = isTRUE(cleanroom_enabled),
    dgp_mode          = dgp_mode,
    nc_criteria       = nc_criteria,
    dq_thresholds     = dq_thresholds,
    roles             = roles,
    # Which prespecified rules the fingerprint commits to. A rule declared
    # later (declare_negative_controls after creation) is recorded in the
    # design log instead of the hash, and validation hashes NULL for it.
    declared_at_creation = list(nc_criteria   = !is.null(nc_criteria),
                                dq_thresholds = !is.null(dq_thresholds)),
    lock_format       = 2L,
    locked_at         = Sys.time(),
    lock_hash         = lock_hash
  )
  class(lock) <- c("ct_lock", "cleanroom_lock")
  if (isTRUE(enforce)) lock$require_authorization <- TRUE
  if (!is.null(estimand)) lock <- attach_estimand(lock, estimand)
  if (!is.null(sensitivity_plans))
    lock <- declare_sensitivity_plan(lock, sensitivity_plans)
  if (!is.null(negative_controls))
    for (nc in negative_controls) lock <- define_negative_control(lock, nc)
  if (isTRUE(mask)) lock <- mask_outcome(lock)
  lock
}


# Content digest of the design data: every column of `data` except the
# outcome column, in sorted column order. The outcome is excluded so the
# fingerprint commits to nothing about outcome values: masking, unmasking,
# or permuting the outcome leaves it unchanged, while any edit to
# covariates, treatment, missingness indicators, or negative-control
# columns invalidates it.
#' @keywords internal
.design_data_digest <- function(data, outcome) {
  design_cols <- sort(setdiff(names(data), outcome))
  .compute_lock_hash(data[design_cols])
}


# ── The outcome store ─────────────────────────────────────────────────────
# Since 0.3.0 the lock's data frame holds design data only (covariates,
# treatment, missingness indicators, negative controls); the primary
# outcome lives in a separate outcome store that Stage 4 joins back.
# mask_outcome() removes the store from its copy of the lock, so a masked
# lock is physically outcome-free. Locks created by earlier versions
# carry the outcome inside $data; every accessor below branches on that.

#' @keywords internal
.new_outcome_store <- function(y, outcome, lock_hash) {
  store <- list(outcome = outcome, y = y, n = length(y),
                lock_hash = lock_hash)
  class(store) <- "ct_outcome_store"
  store
}

#' @export
print.ct_outcome_store <- function(x, ...) {
  cat(sprintf("Outcome store: '%s' (%d rows; joined at Stage 4; lock %s)\n",
              x$outcome, x$n,
              if (is.null(x$lock_hash)) "?" else substr(x$lock_hash, 1, 12)))
  invisible(x)
}

# TRUE for locks created since the outcome-store split.
#' @keywords internal
.lock_has_store_format <- function(lock) isTRUE(lock$lock_format >= 2L)

# The primary outcome vector, or NULL when the lock does not hold one
# (a masked store-format lock). Legacy locks read the data column.
#' @keywords internal
.outcome_vector <- function(lock) {
  if (.lock_has_store_format(lock)) {
    if (is.null(lock$outcome_store)) return(NULL)
    lock$outcome_store$y
  } else {
    lock$data[[lock$outcome]]
  }
}

# Whether the lock's outcome is readable (present and not all NA).
#' @keywords internal
.outcome_readable <- function(lock) {
  y <- .outcome_vector(lock)
  !is.null(y) && !all(is.na(y))
}

# The full analytic data frame: design data plus the outcome column.
# Stage 4 functions call this once at entry, after the outcome guard.
#' @keywords internal
.join_outcome <- function(lock) {
  if (!.lock_has_store_format(lock)) return(lock$data)
  d <- lock$data
  y <- .outcome_vector(lock)
  if (is.null(y))
    stop("The lock holds no outcome (masked store-format lock). ",
         "Unmask with unmask_outcome() before estimation.", call. = FALSE)
  d[[lock$outcome]] <- y
  d
}

# A transient Stage 4 sub-lock over selected rows, with the outcome
# joined into $data (legacy layout). Used by the matched and trimmed
# paths after the outcome guard has passed; never a design-stage object.
#' @keywords internal
.sublock_with_outcome <- function(lock, idx) {
  sub <- lock
  sub$data <- .join_outcome(lock)[idx, , drop = FALSE]
  sub$outcome_store <- NULL
  sub$lock_format <- 1L
  sub
}

# Normalise (and validate) the prespecified negative-control criteria.
#' @keywords internal
.normalize_nc_criteria <- function(x) {
  if (is.null(x)) return(NULL)
  if (!is.list(x) || is.null(x$null_band) || !is.numeric(x$null_band) ||
      length(x$null_band) != 1L || x$null_band <= 0)
    stop("`nc_criteria` must be a list with a positive `null_band` ",
         "half-width (plus optional scale, rule, min_per_domain, ",
         "consistency).", call. = FALSE)
  if (!is.null(x$scale) && !identical(x$scale, "rd"))
    stop("`nc_criteria$scale` supports only \"rd\" for now: the ",
         "negative-control ladder reports risk differences. Declare the ",
         "band on the risk-difference scale.", call. = FALSE)
  out <- list(
    null_band      = as.numeric(x$null_band),
    scale          = "rd",
    rule           = match.arg(x$rule %||% "point_in_band",
                               c("point_in_band", "ci_in_band")),
    min_per_domain = as.integer(x$min_per_domain %||% 1L),
    consistency    = match.arg(x$consistency %||% "all_within_band",
                               c("all_within_band", "majority_within_band"))
  )
  if (out$min_per_domain < 1L)
    stop("`nc_criteria$min_per_domain` must be at least 1.", call. = FALSE)
  extra <- setdiff(names(x), names(out))
  if (length(extra))
    stop("Unknown nc_criteria field(s): ", paste(extra, collapse = ", "),
         call. = FALSE)
  out
}

# Normalise (and validate) the prespecified data-quality thresholds.
#' @keywords internal
.normalize_dq_thresholds <- function(x) {
  if (is.null(x)) return(NULL)
  need <- c("max_abs_bias", "min_coverage", "max_rmse_ratio")
  if (!is.list(x) || !all(need %in% names(x)))
    stop("`dq_thresholds` must be a list with max_abs_bias, min_coverage, ",
         "and max_rmse_ratio (plus optional flag_coverage, ",
         "flag_rmse_ratio).", call. = FALSE)
  out <- list(
    max_abs_bias    = as.numeric(x$max_abs_bias),
    min_coverage    = as.numeric(x$min_coverage),
    max_rmse_ratio  = as.numeric(x$max_rmse_ratio),
    flag_coverage   = as.numeric(x$flag_coverage %||% NA_real_),
    flag_rmse_ratio = as.numeric(x$flag_rmse_ratio %||% NA_real_)
  )
  if (out$max_abs_bias <= 0 || out$min_coverage <= 0 ||
      out$min_coverage >= 1 || out$max_rmse_ratio <= 1)
    stop("dq_thresholds out of range: need max_abs_bias > 0, ",
         "0 < min_coverage < 1, max_rmse_ratio > 1.", call. = FALSE)
  extra <- setdiff(names(x), names(out))
  if (length(extra))
    stop("Unknown dq_thresholds field(s): ", paste(extra, collapse = ", "),
         call. = FALSE)
  out
}

#' Create a Simple (Non-Clean-Room) Analysis Lock
#'
#' A convenience wrapper around [create_analysis_lock()] that disables the
#' software-enforced clean-room machinery (outcome guard, gate authorisation,
#' audit-trail requirement). The resulting lock can be passed to the
#' estimation functions ([run_clean_tmle()], [estimate_effect()],
#' [estimate_ipwrisk()], [estimate_gcomprisk()], [estimate_aipwrisk()],
#' and the rest) with no unmasking ritual. The lock fingerprint and the design
#' diagnostics (`compute_ps_diagnostics()`, `clean_weight_diagnostics()`,
#' `love_plot()`) remain available.
#'
#' Use this constructor when you want cleanTMLE as a plain outcome-blind
#' RWE pipeline (TMLE, IPTW, g-computation, AIPW, IPCW-TMLE) without the
#' staged-analysis governance layer. Outcome blindness is then the
#' analyst's responsibility, not software-enforced.
#'
#' For high-stakes confirmatory studies, use [create_analysis_lock()]
#' (the default) and the full staged workflow.
#'
#' @inheritParams create_analysis_lock
#'
#' @return A `cleanroom_lock` object with `cleanroom_enabled = FALSE`.
#'
#' @examples
#' \dontrun{
#' dat  <- sim_func1(n = 200, seed = 1)
#' lock <- create_simple_lock(
#'   data       = dat,
#'   treatment  = "treatment",
#'   outcome    = "event_24",
#'   covariates = c("age", "sex", "biomarker"),
#'   seed       = 1
#' )
#' isFALSE(lock$cleanroom_enabled)
#'
#' }
#' @keywords internal
create_simple_lock <- function(data, treatment, outcome, covariates,
                                sl_library    = c("SL.glm", "SL.mean"),
                                plasmode_reps = 100L,
                                seed          = 42L) {
  create_analysis_lock(
    data              = data,
    treatment         = treatment,
    outcome           = outcome,
    covariates        = covariates,
    sl_library        = sl_library,
    plasmode_reps     = as.integer(plasmode_reps),
    seed              = as.integer(seed),
    cleanroom_enabled = FALSE
  )
}


#' @keywords internal
.compute_lock_hash <- function(params) {
  if (requireNamespace("digest", quietly = TRUE)) {
    return(digest::digest(params, algo = "sha256", serialize = TRUE))
  }
  # Fallback (digest unavailable): emit a clearly-labelled non-cryptographic
  # checksum so callers can detect that integrity guarantees are weakened.
  s     <- paste(unlist(lapply(params, as.character)), collapse = "|")
  chars <- utf8ToInt(s)
  val   <- sum(chars * seq_along(chars)) %% 1e9
  paste0("checksum-", sprintf("%09.0f", val))
}


#' Validate an Analysis Lock
#'
#' Checks that a `cleanroom_lock` object is complete, internally consistent,
#' and has not been modified after creation (hash check).
#'
#' @param lock A `cleanroom_lock` object from [create_analysis_lock()].
#'
#' @return Invisibly returns `lock` if valid; otherwise throws an error.
#'
#' @examples
#' \dontrun{
#' dat <- sim_func1(n = 200, seed = 1)
#' lock <- create_analysis_lock(
#'   data       = dat,
#'   treatment  = "treatment",
#'   outcome    = "event_24",
#'   covariates = c("age", "sex", "biomarker"),
#'   seed       = 1
#' )
#' validate_analysis_lock(lock)
#'
#' }
#' @keywords internal
validate_analysis_lock <- function(lock) {
  if (!inherits(lock, "cleanroom_lock"))
    stop("`lock` must be a cleanroom_lock object.", call. = FALSE)

  required <- c("data", "treatment", "outcome", "covariates",
                "sl_library", "plasmode_reps", "seed", "lock_hash")
  missing_fields <- setdiff(required, names(lock))
  if (length(missing_fields) > 0L)
    stop("Lock is missing required fields: ",
         paste(missing_fields, collapse = ", "), call. = FALSE)

  # Recompute and compare hash. The fingerprint was computed over the
  # full input data at creation; a store-format lock reconstructs the
  # input's dimensions and names from the design data plus the outcome
  # name, so masking (removing the store) does not disturb validation.
  # The content digest covers every column except the outcome. Locks
  # from earlier versions fall back to the fingerprint scheme of their
  # era so archived locks stay loadable.
  if (is.null(lock$dgp_mode)) {
    computed_hash <- .compute_lock_hash(list(
      treatment  = lock$treatment,
      outcome    = lock$outcome,
      covariates = lock$covariates,
      sl_library = lock$sl_library,
      seed       = as.integer(lock$seed),
      data_nrow  = nrow(lock$data),
      data_ncol  = ncol(lock$data),
      data_names = paste(sort(names(lock$data)), collapse = "|")
    ))
    message("validate_analysis_lock: pre-content-digest lock; the legacy ",
            "fingerprint (names and dimensions only) was checked.")
  } else {
    if (.lock_has_store_format(lock)) {
      full_names <- sort(c(names(lock$data), lock$outcome))
      full_ncol  <- ncol(lock$data) + 1L
    } else {
      full_names <- sort(names(lock$data))
      full_ncol  <- ncol(lock$data)
    }
    # Hash only the rules the fingerprint committed to at creation; a
    # rule declared afterwards lives in the design log, not the hash.
    # Locks predating the flag treated any present rule as
    # creation-declared, so that is the fallback.
    dac <- lock$declared_at_creation
    nc_hashed <- if (is.null(dac)) lock$nc_criteria
                 else if (isTRUE(dac$nc_criteria)) lock$nc_criteria
    dq_hashed <- if (is.null(dac)) lock$dq_thresholds
                 else if (isTRUE(dac$dq_thresholds)) lock$dq_thresholds
    computed_hash <- .compute_lock_hash(list(
      treatment     = lock$treatment,
      outcome       = lock$outcome,
      covariates    = lock$covariates,
      sl_library    = lock$sl_library,
      seed          = as.integer(lock$seed),
      data_nrow     = nrow(lock$data),
      data_ncol     = full_ncol,
      data_names    = paste(full_names, collapse = "|"),
      data_content  = .design_data_digest(lock$data, lock$outcome),
      dgp_mode      = lock$dgp_mode,
      nc_criteria   = nc_hashed,
      dq_thresholds = dq_hashed
    ))
  }
  if (!identical(lock$lock_hash, computed_hash))
    stop("Lock hash mismatch: the lock object may have been modified.",
         call. = FALSE)

  # Variable presence
  if (!lock$treatment %in% names(lock$data))
    stop("treatment variable '", lock$treatment,
         "' not found in locked data.", call. = FALSE)
  if (.lock_has_store_format(lock)) {
    if (!is.null(lock$outcome_store) &&
        !identical(lock$outcome_store$outcome, lock$outcome))
      stop("outcome store names '", lock$outcome_store$outcome,
           "' but the lock declares '", lock$outcome, "'.", call. = FALSE)
    if (!is.null(lock$outcome_store) &&
        lock$outcome_store$n != nrow(lock$data))
      stop("outcome store has ", lock$outcome_store$n, " rows; the lock ",
           "data has ", nrow(lock$data), ".", call. = FALSE)
  } else if (!lock$outcome %in% names(lock$data)) {
    stop("outcome variable '", lock$outcome,
         "' not found in locked data.", call. = FALSE)
  }
  missing_cov <- lock$covariates[!lock$covariates %in% names(lock$data)]
  if (length(missing_cov) > 0L)
    stop("covariates not found in locked data: ",
         paste(missing_cov, collapse = ", "), call. = FALSE)

  message("Analysis lock validated successfully.")
  invisible(lock)
}


#' @export
print.cleanroom_lock <- function(x, ...) {

  cat("cleanTMLE Analysis Lock\n")
  cat("=======================\n")
  if (.lock_has_store_format(x)) {
    cat("Data:       ", nrow(x$data), "observations,", ncol(x$data),
        "design variables (outcome kept in a separate store)\n")
    cat("Outcome:    ", x$outcome,
        if (is.null(x$outcome_store)) " [store removed: masked]"
        else " [in outcome store]", "\n", sep = "")
  } else {
    cat("Data:       ", nrow(x$data), "observations,", ncol(x$data),
        "variables\n")
    cat("Outcome:    ", x$outcome, "\n")
  }
  cat("Treatment:  ", x$treatment, "\n")
  cat("Covariates: ", paste(x$covariates, collapse = ", "), "\n")
  cat("SL library: ", paste(x$sl_library, collapse = ", "), "\n")
  cat("Plasmode:   ", x$plasmode_reps, "replicates\n")
  cat("Seed:       ", x$seed, "\n")
  if (!is.null(x$locked_at))
    cat("Locked at:  ", format(x$locked_at, "%Y-%m-%d %H:%M:%S"), "\n")
  cat("Hash:       ", x$lock_hash, "\n")
  if (!is.null(x$thresholds_hash))
    cat("Thresholds: ", x$thresholds_hash, "(decision rule fingerprint)\n")

  # Estimand
  if (!is.null(x$estimand)) {
    cat("\nEstimand\n")
    cat("--------\n")
    if (!is.null(x$estimand$description))
      cat("Question:    ", x$estimand$description, "\n")
    if (!is.null(x$estimand$population))
      cat("Population:  ", x$estimand$population, "\n")
    if (!is.null(x$estimand$treatment_strategies))
      cat("Contrast:    ", paste(x$estimand$treatment_strategies,
                                 collapse = " vs. "), "\n")
    if (!is.null(x$estimand$outcome_label))
      cat("Outcome:     ", x$estimand$outcome_label, "\n")
    if (!is.null(x$estimand$followup))
      cat("Follow-up:   ", x$estimand$followup, "\n")
    cat("Estimand:    ", x$estimand$contrast, "\n")
    if (!is.null(x$estimand$statistical_estimand))
      cat("Statistical: ", x$estimand$statistical_estimand, "\n")
  }

  # Sensitivity plans
  if (!is.null(x$sensitivity_plans) && length(x$sensitivity_plans) > 0L) {
    cat("\nSensitivity Plans\n")
    cat("-----------------\n")
    for (sp in x$sensitivity_plans) {
      cat(sprintf("  - %s: %s\n", sp$label,
                  if (!is.null(sp$description)) sp$description else ""))
    }
  }

  # Negative controls
  if (!is.null(x$negative_controls) && length(x$negative_controls) > 0L) {
    cat("\nNegative Controls\n")
    cat("-----------------\n")
    for (nc in x$negative_controls) {
      cat(sprintf("  - %s (%s)\n", nc$variable, nc$type))
    }
  }

  # Primary TMLE specification
  if (!is.null(x$primary_tmle_spec)) {
    spec <- x$primary_tmle_spec
    cat("\nPrimary TMLE Specification\n")
    cat("--------------------------\n")
    cat("Candidate:  ", spec$candidate_id, "\n")
    cat("Label:      ", spec$label, "\n")
    cat("Truncation: ", spec$truncation, "\n")
    cat("G-library:  ", paste(spec$g_library, collapse = ", "), "\n")
    if (!is.null(spec$q_library))
      cat("Q-library:  ", paste(spec$q_library, collapse = ", "), "\n")
    if (!is.null(spec$selection_rule))
      cat("Selected by:", spec$selection_rule, "\n")
  }

  invisible(x)
}


# ── Stage 2a: Propensity Score ────────────────────────────────────────────

#' Fit Propensity Score Using SuperLearner
#'
#' Estimates the treatment propensity score P(A=1|W) using SuperLearner
#' ensemble learning. Only covariates (W) and treatment (A) are used;
#' the outcome is never accessed, so this step is safe at Stage 2a.
#'
#' @param lock A `cleanroom_lock` from [create_analysis_lock()].
#' @param truncate Numeric; propensity-score truncation bound in (0, 0.5).
#'   Predicted scores are clamped to `[truncate, 1 - truncate]`. Default 0.01.
#' @param cv_folds Integer; number of SuperLearner cross-validation folds.
#'   Default 10.
#' @param cluster Optional parallel cluster (from \pkg{parallel}) for the
#'   SuperLearner fit. Default `NULL` (sequential).
#' @param ... Additional arguments passed to `SuperLearner::SuperLearner()`.
#'
#' @return An object of class `ps_fit` containing propensity scores (`ps`),
#'   the fitted SuperLearner object (`sl_fit`), and metadata.
#'
#' @details
#' Requires the \pkg{SuperLearner} package.
#'
#' @examples
#' \dontrun{
#' dat  <- sim_func1(n = 500, seed = 1)
#' lock <- create_analysis_lock(dat, "treatment", "event_24",
#'                              c("age", "sex", "biomarker"), seed = 1)
#' ps_fit <- fit_ps_superlearner(lock)
#' }
#'
#' @keywords internal
fit_ps_superlearner <- function(lock, truncate = 0.01,
                                 cv_folds = 10L,
                                 cluster = NULL,
                                 ...) {
  if (!inherits(lock, "cleanroom_lock"))
    stop("`lock` must be a cleanroom_lock object.", call. = FALSE)
  if (!requireNamespace("SuperLearner", quietly = TRUE))
    stop("Package 'SuperLearner' is required. ",
         "Install with: install.packages('SuperLearner')", call. = FALSE)

  data <- lock$data
  A    <- data[[lock$treatment]]
  W    <- data[, lock$covariates, drop = FALSE]

  # Resolve the SuperLearner library. If a primary TMLE spec has been
  # locked with `g_library`, that takes precedence over the lock-level
  # `sl_library` -- otherwise switching candidates after locking would
  # silently leave the PS fit using the original lock library.
  primary_spec <- lock$primary_tmle_spec
  sl_library   <- lock$sl_library
  if (!is.null(primary_spec) && !is.null(primary_spec$g_library))
    sl_library <- primary_spec$g_library

  # Default to sequential evaluation. The "invalid connection" error
  # observed in 0.1.0 came from a stale snow/foreach backend registered in
  # the user's session, not from SuperLearner itself; we therefore use the
  # plain SuperLearner() entry point unless the caller passes an explicit
  # cluster. With cluster supplied, dispatch to snowSuperLearner().
  sl_args <- list(
    Y          = A,
    X          = W,
    family     = binomial(),
    SL.library = sl_library,
    env        = .cleantmle_sl_env(),
    cvControl  = list(V = cv_folds)
  )

  extra_args <- list(...)
  sl_args[names(extra_args)] <- extra_args

  withr::local_seed(lock$seed)
  sl_fit <- tryCatch({
    if (!is.null(cluster)) {
      sl_args$cluster <- cluster
      do.call(SuperLearner::snowSuperLearner, sl_args)
    } else {
      do.call(SuperLearner::SuperLearner, sl_args)
    }
  }, error = function(e) {
    stop("fit_ps_superlearner: SuperLearner failed (", e$message,
         "). If the error mentions 'invalid connection', a stale ",
         "parallel cluster is registered in the session; restart R.",
         call. = FALSE)
  })

  ps <- as.numeric(sl_fit$SL.predict)
  ps_raw <- ps

  # Truncate to [truncate, 1 - truncate] to bound subsequent IPW weights.
  # The untruncated scores are kept as `ps_raw`: support assessment must see
  # the propensity the model actually produced, because truncation caps the
  # very weights whose size is the diagnostic.
  if (!is.null(truncate)) {
    if (!is.numeric(truncate) || length(truncate) != 1L ||
        truncate <= 0 || truncate >= 0.5)
      stop("`truncate` must be a single numeric in (0, 0.5).", call. = FALSE)
    ps <- pmin(pmax(ps, truncate), 1 - truncate)
  }

  result <- list(
    ps         = ps,
    ps_raw     = ps_raw,
    sl_fit     = sl_fit,
    treatment  = lock$treatment,
    covariates = lock$covariates,
    data       = data,
    sl_library = sl_library,
    truncate   = truncate,
    cv_folds   = cv_folds,
    method     = "superlearner",
    call       = match.call()
  )
  class(result) <- c("ps_fit", "cr_result")
  result
}


#' @export
print.ps_fit <- function(x, ...) {
  method_label <- if (isTRUE(x$method == "glm")) "GLM" else "SuperLearner"
  cat(sprintf("Propensity Score Fit (%s)\n", method_label))
  cat("====================================\n")
  cat("Treatment:  ", x$treatment, "\n")
  cat("Covariates: ", paste(x$covariates, collapse = ", "), "\n")
  cat(sprintf("PS range:    [%.4f, %.4f]\n", min(x$ps), max(x$ps)))
  cat(sprintf("PS mean:      %.4f\n", mean(x$ps)))
  invisible(x)
}


#' Fit a Logistic Regression Propensity Score Model
#'
#' A conventional logistic-regression alternative to [fit_ps_superlearner()]
#' for use when `SuperLearner` is not installed or for rapid sensitivity
#' analyses.
#'
#' @param lock A `cleanroom_lock` from [create_analysis_lock()].
#' @param truncate Numeric in (0, 0.5); propensity scores are bounded to
#'   `[truncate, 1 - truncate]` to avoid extreme weights. Default: `0.01`.
#'
#' @return An object of class `ps_fit` (same structure as
#'   [fit_ps_superlearner()], but with `method = "glm"` and `glm_fit`
#'   instead of `sl_fit`).
#'
#' @examples
#' \dontrun{
#' dat  <- sim_func1(n = 200, seed = 1)
#' lock <- create_analysis_lock(
#'   data = dat, treatment = "treatment", outcome = "event_24",
#'   covariates = c("age", "sex", "biomarker"), seed = 1L
#' )
#' ps_fit <- fit_ps_glm(lock)
#' print(ps_fit)
#'
#' }
#' @keywords internal
fit_ps_glm <- function(lock, truncate = 0.01) {
  if (!inherits(lock, "cleanroom_lock"))
    stop("`lock` must be a cleanroom_lock object.", call. = FALSE)

  data       <- lock$data
  treatment  <- lock$treatment
  covariates <- lock$covariates
  A          <- data[[treatment]]

  fml     <- stats::reformulate(covariates, response = treatment)
  glm_fit <- stats::glm(fml, data = data, family = stats::binomial())
  ps_raw  <- as.numeric(stats::predict(glm_fit, type = "response"))
  ps      <- pmax(pmin(ps_raw, 1 - truncate), truncate)

  result <- list(
    ps         = ps,
    ps_raw     = ps_raw,
    glm_fit    = glm_fit,
    sl_fit     = NULL,
    treatment  = treatment,
    covariates = covariates,
    data       = data,
    sl_library = NULL,
    method     = "glm",
    truncate   = truncate,
    call       = match.call()
  )
  class(result) <- c("ps_fit", "cr_result")
  result
}


#' Compute Propensity Score Diagnostics
#'
#' Computes overlap diagnostics, effective sample size (ESS), and
#' standardized mean differences (SMDs) from a fitted propensity score object.
#'
#' @param ps_fit A `ps_fit` object from [fit_ps_superlearner()].
#' @param ... Currently unused.
#'
#' @return An object of class `ps_diagnostics` containing:
#'   * `ess` - effective sample size table
#'   * `smds` - standardized mean differences before and after weighting
#'   * `overlap_plot` - a ggplot2 histogram of PS by treatment group
#'
#' @examples
#' \dontrun{
#' diag <- compute_ps_diagnostics(ps_fit)
#' print(diag)
#' plot(diag)
#' }
#'
#' @keywords internal
compute_ps_diagnostics <- function(ps_fit, ...) {
  if (!inherits(ps_fit, "ps_fit"))
    stop("`ps_fit` must be a ps_fit object from fit_ps_superlearner().",
         call. = FALSE)

  data <- ps_fit$data
  A    <- data[[ps_fit$treatment]]
  ps   <- ps_fit$ps

  # IPTW weights (unstabilized) for ESS and weighted SMDs
  w <- ifelse(A == 1, 1 / ps, 1 / (1 - ps))

  # Effective sample size (Kish)
  ess_treated <- sum(w[A == 1])^2 / sum(w[A == 1]^2)
  ess_control <- sum(w[A == 0])^2 / sum(w[A == 0]^2)
  n_treated   <- sum(A == 1)
  n_control   <- sum(A == 0)

  # Standardized mean differences
  smd_rows <- lapply(ps_fit$covariates, function(v) {
    x_val <- data[[v]]
    if (!is.numeric(x_val)) x_val <- as.numeric(as.factor(x_val)) - 1L

    m1 <- mean(x_val[A == 1])
    m0 <- mean(x_val[A == 0])
    s1 <- var(x_val[A == 1])
    s0 <- var(x_val[A == 0])
    pooled_sd <- sqrt((s1 + s0) / 2)
    smd_unw <- if (pooled_sd > 0) (m1 - m0) / pooled_sd else 0

    m1w <- weighted.mean(x_val[A == 1], w[A == 1])
    m0w <- weighted.mean(x_val[A == 0], w[A == 0])
    wv1 <- .weighted_var_cl(x_val[A == 1], w[A == 1])
    wv0 <- .weighted_var_cl(x_val[A == 0], w[A == 0])
    pooled_w <- sqrt((wv1 + wv0) / 2)
    smd_w <- if (pooled_w > 0) (m1w - m0w) / pooled_w else 0

    data.frame(
      variable       = v,
      smd_unweighted = round(smd_unw, 4),
      smd_weighted   = round(smd_w,   4),
      stringsAsFactors = FALSE
    )
  })
  smd_table <- do.call(rbind, smd_rows)

  # Overlap histogram
  plot_data <- data.frame(
    ps    = c(ps[A == 1],            ps[A == 0]),
    group = c(rep("Treated (A=1)", sum(A == 1)),
              rep("Control (A=0)", sum(A == 0))),
    stringsAsFactors = FALSE
  )
  overlap_plot <- ggplot2::ggplot(
    plot_data,
    ggplot2::aes(x = .data$ps, fill = .data$group)
  ) +
    ggplot2::geom_histogram(position = "identity", alpha = 0.5, bins = 30L) +
    ggplot2::labs(
      x     = "Propensity Score",
      y     = "Count",
      title = "Propensity Score Overlap",
      fill  = "Group"
    ) +
    ggplot2::theme_minimal()

  result <- list(
    ess = data.frame(
      group   = c("Treated", "Control", "Total"),
      n       = c(n_treated, n_control, n_treated + n_control),
      ess     = round(c(ess_treated, ess_control,
                        ess_treated + ess_control), 1L),
      ess_pct = round(100 * c(ess_treated / n_treated,
                               ess_control / n_control,
                               (ess_treated + ess_control) /
                                 (n_treated + n_control)), 1L),
      stringsAsFactors = FALSE
    ),
    smds         = smd_table,
    overlap_plot = overlap_plot,
    ps_summary   = summary(ps),
    call         = match.call()
  )
  class(result) <- "ps_diagnostics"
  result
}


#' @keywords internal
.weighted_var_cl <- function(x, w) {
  w   <- w / sum(w)
  xbar <- sum(w * x)
  sum(w * (x - xbar)^2)
}


#' @export
print.ps_diagnostics <- function(x, ...) {
  cat("Propensity Score Diagnostics\n")
  cat("============================\n\n")
  cat("Effective Sample Size (Kish ESS):\n")
  print(x$ess, row.names = FALSE)
  cat("\nStandardized Mean Differences:\n")
  print(x$smds, row.names = FALSE)
  cat("\nOverlap plot available via plot(diag).\n")
  invisible(x)
}


#' @export
plot.ps_diagnostics <- function(x, ...) {
  x$overlap_plot
}


# ── TMLE Candidate Specifications ─────────────────────────────────────────

#' Create a TMLE Candidate Specification
#'
#' Defines a single fully-specified TMLE candidate for use in Stage 2b
#' plasmode evaluation. Each candidate specifies the nuisance-learner
#' libraries and PS truncation threshold.
#'
#' @section Clean-room stage: Stage 1a / 2b (pre-outcome).
#'
#' @param candidate_id Character; unique identifier for this candidate.
#' @param label Character; human-readable description.
#' @param g_library Character vector; SuperLearner library for the
#'   treatment mechanism (g-model). Falls back to GLM internally.
#' @param q_library Character vector; SuperLearner library for the
#'   outcome mechanism (Q-model). If \code{NULL}, defaults to
#'   \code{g_library}.
#' @param truncation Numeric; PS truncation threshold in (0, 0.5).
#'   Default: 0.01.
#' @param variance_method Character; variance estimator for the candidate,
#'   one of \code{"IF"} (influence-function; default), \code{"cv_IF"}
#'   (cross-validated IF), or \code{"robust"}.
#' @param cv_scheme Character; cross-fitting scheme, one of \code{"none"}
#'   (default), \code{"cv_tmle"}, or \code{"sample_split"}.
#' @param cv_V Integer or \code{NULL}; number of cross-fitting folds when
#'   \code{cv_scheme} is not \code{"none"}. \code{NULL} uses a default.
#' @param estimator Character; the estimator family, one of \code{"tmle"}
#'   (default), \code{"aipw"}, or \code{"onestep"}.
#' @param discrete_sl Logical; if \code{TRUE}, use the discrete (single
#'   best-learner) SuperLearner rather than the convex ensemble.
#' @param screener Character; variable-screening algorithm prepended to the
#'   learner library, one of \code{"All"} (default), \code{"corP"},
#'   \code{"corRank"}, \code{"glmnet"}, or \code{"randomForest"}.
#' @param tmle_control List or \code{NULL}; optional control parameters
#'   passed through to the targeting step.
#' @param match_spec List or \code{NULL}; optional matched-cohort
#'   specification when the candidate operates on a matched subset.
#' @param ... Reserved for future candidate options.
#'
#' @return A list of class \code{tmle_candidate_spec}.
#'
#' @examples
#' \dontrun{
#' tmle_candidate("glm_t01", "GLM, trunc=0.01",
#'                g_library = "SL.glm", truncation = 0.01)
#'
#' }
#' @keywords internal
tmle_candidate <- function(candidate_id, label = candidate_id,
                           g_library       = c("SL.glm"),
                           q_library       = NULL,
                           truncation      = 0.01,
                           variance_method = c("IF", "cv_IF", "robust"),
                           cv_scheme       = c("none", "cv_tmle", "sample_split"),
                           cv_V            = NULL,
                           estimator       = c("tmle", "aipw", "onestep"),
                           discrete_sl     = FALSE,
                           screener        = c("All", "corP", "corRank",
                                                "glmnet", "randomForest"),
                           tmle_control    = NULL,
                           match_spec      = NULL,
                           ...) {
  if (!is.character(candidate_id) || length(candidate_id) != 1L)
    stop("`candidate_id` must be a single character string.", call. = FALSE)

  # Accept Q_library / Q_libraries as deprecated capital-Q aliases for
  # q_library to match the older vignette spelling.
  dots <- list(...)
  alt_q <- intersect(names(dots), c("Q_library", "Q_libraries", "q_libraries"))
  if (length(alt_q) > 0L) {
    if (is.null(q_library)) q_library <- dots[[alt_q[1]]]
    warning("tmle_candidate(): argument `", alt_q[1],
            "` is deprecated; use `q_library` instead.", call. = FALSE)
    dots[alt_q] <- NULL
  }
  if (length(dots) > 0L)
    stop("Unused arguments: ", paste(names(dots), collapse = ", "),
         call. = FALSE)

  if (is.null(q_library)) q_library <- g_library

  variance_method <- match.arg(variance_method)
  cv_scheme       <- match.arg(cv_scheme)
  estimator       <- match.arg(estimator)
  screener        <- match.arg(screener)

  # Resolve named truncation rule lazily (sample size unknown here)
  truncation_rule <- if (is.character(truncation)) truncation else NULL

  if (is.null(tmle_control)) {
    tmle_control <- list(
      fluctuation = "logistic",
      alpha       = 0.9995,
      target.gwt  = TRUE,
      automate    = FALSE,
      min.retain  = 2L,
      cv_Qinit    = TRUE
    )
  }

  obj <- list(
    candidate_id    = candidate_id,
    label           = label,
    g_library       = g_library,
    q_library       = q_library,
    truncation      = truncation,
    truncation_rule = truncation_rule,
    variance_method = variance_method,
    cv_scheme       = cv_scheme,
    cv_V            = cv_V,
    estimator       = estimator,
    discrete_sl     = isTRUE(discrete_sl),
    screener        = screener,
    tmle_control    = tmle_control,
    match_spec      = match_spec
  )
  class(obj) <- "tmle_candidate_spec"
  obj
}


#' @export
print.tmle_candidate_spec <- function(x, ...) {
  cat(sprintf("TMLE Candidate: %s\n", x$candidate_id))
  cat(sprintf("  Label:           %s\n", x$label))
  cat(sprintf("  Estimator:       %s\n", x$estimator %||% "tmle"))
  cat(sprintf("  CV scheme:       %s\n", x$cv_scheme %||% "none"))
  cat(sprintf("  Variance method: %s\n", x$variance_method %||% "IF"))
  cat(sprintf("  g-library:       %s\n", paste(x$g_library, collapse = ", ")))
  cat(sprintf("  Q-library:       %s\n", paste(x$q_library, collapse = ", ")))
  cat(sprintf("  Discrete SL:     %s\n", isTRUE(x$discrete_sl)))
  cat(sprintf("  Screener:        %s\n", x$screener %||% "All"))
  cat(sprintf("  Truncation:      %s\n", x$truncation))
  if (!is.null(x$match_spec))
    cat("  Matched cohort: yes\n")
  invisible(x)
}

# NOTE: `%||%` is already defined in tmle_clean_room_wrapper.R and
# plasmode_dq.R; relying on those package-internal definitions.


#' Check a List of TMLE Candidate Specifications
#'
#' Checks that all candidates are \code{tmle_candidate_spec} objects with
#' unique IDs and valid settings.
#'
#' @param candidates A list of \code{tmle_candidate_spec} objects.
#'
#' @return Invisibly returns \code{candidates} if valid; errors otherwise.
#'
#' @keywords internal
validate_tmle_candidates <- function(candidates) {
  if (!is.list(candidates) || length(candidates) == 0L)
    stop("`candidates` must be a non-empty list.", call. = FALSE)

  for (i in seq_along(candidates)) {
    if (!inherits(candidates[[i]], "tmle_candidate_spec"))
      stop(sprintf("Element %d is not a tmle_candidate_spec.", i),
           call. = FALSE)
  }

  ids <- vapply(candidates, function(x) x$candidate_id, character(1))
  if (anyDuplicated(ids))
    stop("Duplicate candidate IDs: ",
         paste(ids[duplicated(ids)], collapse = ", "), call. = FALSE)

  invisible(candidates)
}


#' Generate a Default TMLE Candidate Grid
#'
#' Creates a compact grid of TMLE candidates varying the SuperLearner
#' library and PS truncation level. Useful as a default when the user
#' does not supply a custom candidate list.
#'
#' @param truncations Numeric vector of truncation thresholds.
#'   Default: \code{c(0.01, 0.05)}.
#' @param libraries A named list of SuperLearner library vectors.
#'   Default includes GLM-only and GLM+mean.
#' @param estimators Character vector of estimator families to cross with
#'   the library/truncation grid (e.g. \code{"tmle"}). Default \code{"tmle"}.
#' @param cv_schemes Character vector of cross-fitting schemes. Default
#'   \code{"none"}.
#' @param variance_methods Character vector of variance estimators. Default
#'   \code{"IF"}.
#' @param discrete_sl Logical; passed to each candidate (discrete vs convex
#'   SuperLearner). Default \code{FALSE}.
#' @param screeners Character vector of variable screeners to cross into the
#'   grid. Default \code{"All"}.
#' @param max_candidates Integer; cap on the number of generated candidates
#'   (guards against combinatorial blow-up). Default 64.
#' @param ... Reserved for future grid dimensions / back-compatible aliases.
#'
#' @return A list of \code{tmle_candidate_spec} objects.
#'
#' @examples
#' \dontrun{
#' grid <- expand_tmle_candidate_grid()
#' length(grid)
#' grid[[1]]
#'
#' }
#' @keywords internal
expand_tmle_candidate_grid <- function(
    truncations = c(0.01, 0.05),
    libraries   = list(
      glm      = c("SL.glm"),
      glm_mean = c("SL.glm", "SL.mean")
    ),
    estimators       = "tmle",
    cv_schemes       = "none",
    variance_methods = "IF",
    discrete_sl      = FALSE,
    screeners        = "All",
    max_candidates   = 64L,
    ...) {

  # Back-compat: accept g_libraries / Q_libraries from the older docs and
  # treat them as a paired enumeration when both are supplied.
  dots <- list(...)
  if (length(dots) > 0L) {
    if (!is.null(dots$g_libraries) && !is.null(dots$Q_libraries)) {
      warning("expand_tmle_candidate_grid(): `g_libraries`/`Q_libraries` ",
              "are deprecated; use a single `libraries` list (g and Q ",
              "share the library by default).", call. = FALSE)
      g_libs <- dots$g_libraries; q_libs <- dots$Q_libraries
      candidates <- list()
      for (i in seq_along(g_libs)) {
        for (trunc in truncations) {
          cand_id <- sprintf("tmle_lib%d_t%g", i, trunc)
          candidates[[cand_id]] <- tmle_candidate(
            candidate_id = cand_id, label = cand_id,
            g_library = g_libs[[i]], q_library = q_libs[[i]],
            truncation = trunc)
        }
      }
      return(candidates)
    }
    extra <- setdiff(names(dots), c("g_libraries", "Q_libraries"))
    if (length(extra) > 0L)
      stop("Unused arguments: ", paste(extra, collapse = ", "),
           call. = FALSE)
  }

  # Build the Cartesian product across the workshop-driven axes.
  # When all extra axes are at their single-value defaults, the
  # behaviour is identical to the pre-0.2 grid.
  candidates <- list()
  grid <- expand.grid(
    lib_name = names(libraries),
    trunc    = truncations,
    est      = estimators,
    cvs      = cv_schemes,
    var_m    = variance_methods,
    dsl      = discrete_sl,
    scr      = screeners,
    stringsAsFactors = FALSE,
    KEEP.OUT.ATTRS   = FALSE
  )
  if (nrow(grid) > max_candidates)
    stop("Candidate grid has ", nrow(grid), " entries (max_candidates = ",
         max_candidates, "). Narrow the axes or raise max_candidates.",
         call. = FALSE)

  for (r in seq_len(nrow(grid))) {
    g       <- grid[r, ]
    trunc_label <- sub("\\.", "", sprintf("t%s", g$trunc))
    parts   <- c(g$lib_name, trunc_label,
                 if (g$est != "tmle")   g$est,
                 if (g$cvs != "none")   g$cvs,
                 if (g$var_m != "IF")   g$var_m,
                 if (isTRUE(g$dsl))     "dSL",
                 if (g$scr != "All")    g$scr)
    cand_id <- paste(c("tmle", parts), collapse = "_")
    label   <- sprintf(
      "%s: lib=%s, trunc=%.3g, cv=%s, var=%s%s%s",
      toupper(g$est),
      paste(libraries[[g$lib_name]], collapse = "+"),
      g$trunc, g$cvs, g$var_m,
      if (isTRUE(g$dsl)) ", dSL" else "",
      if (g$scr != "All") paste0(", screener=", g$scr) else ""
    )
    candidates[[cand_id]] <- tmle_candidate(
      candidate_id    = cand_id,
      label           = label,
      g_library       = libraries[[g$lib_name]],
      q_library       = libraries[[g$lib_name]],
      truncation      = g$trunc,
      estimator       = g$est,
      cv_scheme       = g$cvs,
      variance_method = g$var_m,
      discrete_sl     = isTRUE(g$dsl),
      screener        = g$scr
    )
  }
  candidates
}


# ── Stage 2b: Plasmode Feasibility ────────────────────────────────────────

# Shared contract checks for the plasmode generator mode. `pilot_q0` is
# required in external-pilot mode and refused in hybrid mode (the mode is
# a locked field; mixing signals a specification error, not a preference).
#' @keywords internal
.check_pilot_q0_contract <- function(dgp_mode, pilot_q0, q0_library) {
  if (dgp_mode == "external_pilot") {
    if (is.null(pilot_q0))
      stop("The lock's dgp_mode is 'external_pilot': supply `pilot_q0` ",
           "(baseline event probabilities fitted on external pilot data, ",
           "as a numeric vector over the lock rows or a function of the ",
           "covariate data frame).", call. = FALSE)
    if (!is.null(q0_library))
      warning("`q0_library` is ignored in external_pilot mode: the ",
              "baseline surface comes from `pilot_q0`.", call. = FALSE)
  } else if (!is.null(pilot_q0)) {
    stop("`pilot_q0` was supplied but the lock's dgp_mode is 'hybrid'. ",
         "The generator mode is a locked field: create the lock with ",
         "dgp_mode = 'external_pilot' to use pilot-based generation.",
         call. = FALSE)
  }
  invisible(TRUE)
}

# Resolve external-pilot baseline probabilities to a length-n vector.
#' @keywords internal
.resolve_pilot_q0 <- function(pilot_q0, data, covariates, n) {
  p <- if (is.function(pilot_q0))
    pilot_q0(data[, covariates, drop = FALSE]) else pilot_q0
  p <- as.numeric(p)
  if (length(p) != n)
    stop("`pilot_q0` must yield one probability per lock row (need ", n,
         ", got ", length(p), ").", call. = FALSE)
  if (anyNA(p) || any(p <= 0) || any(p >= 1))
    stop("`pilot_q0` values must lie strictly inside (0, 1) with no NA.",
         call. = FALSE)
  p
}

# Mann-Whitney c-statistic of a propensity against treatment.
#' @keywords internal
.ps_cstat <- function(g, A) {
  n1 <- sum(A == 1); n0 <- sum(A == 0)
  if (n1 == 0 || n0 == 0) return(NA_real_)
  (sum(rank(g)[A == 1]) - n1 * (n1 + 1) / 2) / (n1 * n0)
}

# Hybrid-mode leakage warning under strong separation: E[Y|W] read
# alongside a well-separated propensity approximates the outcome rate by
# propensity stratum, which approaches the crude treatment-outcome
# association (see the generator-modes section of the plasmode docs).
#' @keywords internal
.warn_hybrid_separation <- function(cstat, caller) {
  if (is.finite(cstat) && cstat > 0.80)
    warning(caller, ": the propensity c-statistic is ",
            sprintf("%.3f", cstat), " (> 0.80). In hybrid mode the ",
            "covariate-only Q0 is fitted on the real outcome, and under ",
            "separation this strong E[Y|W] approximates the outcome rate ",
            "by propensity stratum, i.e. the crude treatment-outcome ",
            "association. Consider dgp_mode = 'external_pilot'.",
            call. = FALSE)
  invisible(cstat)
}


#' Run Plasmode-Simulation Feasibility Evaluation
#'
#' Evaluates the performance of prespecified TMLE candidate specifications
#' using plasmode simulation.  Simulated binary outcomes are generated from
#' a parametric baseline-risk model fit on the real covariates, augmented
#' with a specified additive treatment effect.  Each TMLE candidate is fit
#' on every replicate and performance metrics (bias, RMSE, coverage,
#' empirical SD, mean SE) are computed against the known true effect.
#'
#' @section Clean-room stage: Stage 2b (pre-outcome).  The real
#'   treatment--outcome association is never used; only simulated
#'   outcomes are generated.
#'
#' @section Generator modes and what each can leak:
#' The baseline outcome surface Q0(W) comes from the lock's locked
#' `dgp_mode`. In `"hybrid"` mode (the default) Q0 is fitted on the
#' lock's real outcome against the covariates only, so the analyst's
#' code holds an estimate of E\[Y|W\]. When the covariates strongly
#' separate the arms, E\[Y|W\] read alongside the fitted propensity
#' approximates the outcome rate within high- and low-propensity strata,
#' which approaches the crude treatment-outcome association as the
#' c-statistic grows; the function therefore warns when the propensity
#' c-statistic exceeds 0.80, and neither the fitted Q0 object nor its
#' predictions nor the lock data appear anywhere in the returned object
#' (only per-replicate candidate metrics do). In `"external_pilot"` mode
#' the caller supplies `pilot_q0` fitted on external pilot data; the
#' primary outcome column is never read, and the function runs on a
#' masked lock.
#'
#' @param lock A `cleanroom_lock` from [create_analysis_lock()]; its
#'   locked `dgp_mode` selects the generator mode.
#' @param tmle_candidates A list of \code{\link{tmle_candidate}} objects.
#'   If \code{NULL}, a default grid is generated via
#'   \code{\link{expand_tmle_candidate_grid}}.
#' @param effect_sizes Numeric vector of true risk differences to simulate.
#'   Default: `c(0.05, 0.10)`.
#' @param reps Integer; number of plasmode replicates per effect size.
#'   Defaults to `lock$plasmode_reps`.
#' @param q0_library Optional SuperLearner library used to fit the
#'   baseline outcome model Q0 (covariates only) that generates the
#'   simulated outcomes for candidate selection (`"hybrid"` mode only).
#'   If `NULL` (default), Q0 is a logistic GLM in the covariates, which
#'   biases candidate selection toward learners that handle
#'   linear-in-logit structure well. Set this to a SuperLearner library
#'   (e.g. `c("SL.glm", "SL.glmnet", "SL.gam", "SL.mean")`) to generate
#'   simulated outcomes from a richer Q0 surface.
#' @param pilot_q0 Required when the lock's `dgp_mode` is
#'   `"external_pilot"` (and refused otherwise): baseline event
#'   probabilities for the lock rows, fitted on external pilot data.
#'   Either a numeric vector of length `nrow(lock$data)` with values in
#'   (0, 1), or a function that receives the covariate data frame and
#'   returns such a vector.
#' @param design Character; how each simulated replicate is generated.
#'   `"generate_treatment"` (default) resamples the covariate rows with
#'   replacement and draws treatment from a propensity model fitted on the
#'   real data, so the simulated data satisfy positivity to the same extent
#'   the cohort does. `"sample_treatment"` keeps every subject's observed
#'   treatment and simulates only the outcome, which Shaw et al. (2025,
#'   arXiv:2504.11740) show induces a positivity violation by construction
#'   (propensity-based estimators appear biased and undercover even when the
#'   source population has no violation); it is retained only for
#'   reproducing legacy runs and warns when used.
#' @param verbose Logical; if `TRUE`, print progress messages and emit
#'   warnings when any candidate converges on < 50% of inner reps.
#'   Default: `FALSE`.
#'
#' @references Shaw PA, Gruber S, Williamson BD, Desai R, Shortreed SM,
#'   Krakauer C, Nelson JC, van der Laan MJ (2025). A cautionary note for
#'   plasmode simulation studies in the setting of causal inference.
#'   arXiv:2504.11740.
#'
#' @return An object of class `plasmode_results` containing:
#'   * `metrics` - data.frame with one row per candidate per effect size
#'   * `results` - raw per-replicate estimates (nested list)
#'   * `tmle_candidates` - the candidate list used
#'
#' @examples
#' \dontrun{
#' dat  <- sim_func1(n = 500, seed = 1)
#' lock <- create_analysis_lock(dat, "treatment", "event_24",
#'                              c("age", "sex", "biomarker"),
#'                              plasmode_reps = 20L, seed = 1)
#' plas <- run_plasmode_feasibility(lock, reps = 20L)
#' print(plas)
#' }
#'
#' @keywords internal
run_plasmode_feasibility <- function(lock,
                                      tmle_candidates = NULL,
                                      effect_sizes    = c(0.05, 0.10),
                                      reps            = lock$plasmode_reps,
                                      q0_library      = NULL,
                                      design          = c("generate_treatment",
                                                          "sample_treatment"),
                                      pilot_q0        = NULL,
                                      verbose         = FALSE) {
  if (!inherits(lock, "cleanroom_lock"))
    stop("`lock` must be a cleanroom_lock object.", call. = FALSE)
  design   <- match.arg(design)
  dgp_mode <- lock$dgp_mode %||% "hybrid"
  .check_pilot_q0_contract(dgp_mode, pilot_q0, q0_library)
  if (design == "sample_treatment") {
    rlang::warn(paste(
      "design = 'sample_treatment' keeps each subject's observed treatment and",
      "simulates only the outcome. Shaw et al. (2025, arXiv:2504.11740) show",
      "this induces a positivity violation by construction, so",
      "propensity-based estimators look biased and undercover even when the",
      "source population has no violation. Use the default",
      "design = 'generate_treatment' unless you are reproducing legacy runs."),
      .frequency = "once", .frequency_id = "plasmode_sample_treatment")
  }

  # Default candidate grid if not supplied

  if (is.null(tmle_candidates)) {
    tmle_candidates <- expand_tmle_candidate_grid()
  }
  validate_tmle_candidates(tmle_candidates)

  data       <- lock$data
  treatment  <- lock$treatment
  outcome    <- lock$outcome
  covariates <- lock$covariates

  A <- data[[treatment]]
  n <- nrow(data)

  if (dgp_mode == "external_pilot") {
    # External-pilot mode: the baseline surface comes from pilot data
    # supplied by the caller; the primary outcome column is never read,
    # so this path runs on a masked lock.
    p_base <- .resolve_pilot_q0(pilot_q0, data, covariates, n)
  } else {
    Y <- .outcome_vector(lock)
    n_obs_y <- if (is.null(Y)) 0L else sum(!is.na(Y))
    if (n_obs_y == 0L) {
      stop("Q0 model cannot be fit: the lock holds no readable outcome. ",
           "If the outcome is masked, either call unmask_outcome() before ",
           "run_plasmode_feasibility() or lock dgp_mode = 'external_pilot' ",
           "and supply `pilot_q0`.", call. = FALSE)
    }
    if (n_obs_y < length(covariates) + 1L) {
      stop("Q0 model cannot be fit: only ", n_obs_y, " non-NA outcome ",
           "rows available for ", length(covariates), " covariates.",
           call. = FALSE)
    }

    # Hybrid mode: fit the baseline outcome model on the real outcome
    # (covariates only; the treatment--outcome association is not used).
    # By default Q0 is a logistic GLM; pass `q0_library` (a SuperLearner
    # library) to use a richer Q0 -- important when the true outcome
    # surface is nonlinear, so that the simulated outcomes used to
    # select among candidates are not biased toward linear-in-logit
    # learners.
    if (is.null(q0_library)) {
      q_data <- data[, covariates, drop = FALSE]
      q_data[[outcome]] <- Y
      Q0_fml <- stats::reformulate(covariates, response = outcome)
      Q0_fit <- stats::glm(Q0_fml, data = q_data, family = stats::binomial(),
                           na.action = stats::na.exclude)
      p_base <- as.numeric(stats::predict(Q0_fit, type = "response",
                                           newdata = q_data))
    } else {
      if (!requireNamespace("SuperLearner", quietly = TRUE))
        stop("`q0_library` requested but 'SuperLearner' package is not ",
             "available. Install it or leave q0_library = NULL.",
             call. = FALSE)
      cc <- !is.na(Y) & stats::complete.cases(data[, covariates,
                                                   drop = FALSE])
      p_base <- rep(NA_real_, nrow(data))
      Q0_sl <- SuperLearner::SuperLearner(
        Y          = Y[cc],
        X          = data[cc, covariates, drop = FALSE],
        family     = stats::binomial(),
        SL.library = q0_library,
        env        = .cleantmle_sl_env()
      )
      p_base[cc] <- as.numeric(Q0_sl$SL.predict)
      if (any(!cc)) {
        p_base[!cc] <- as.numeric(
          stats::predict(Q0_sl,
                         newdata = data[!cc, covariates, drop = FALSE])$pred)
      }
    }
  }
  # Handle NAs from quasi-separation or missing covariate values
  if (any(is.na(p_base))) {
    p_base[is.na(p_base)] <- mean(p_base, na.rm = TRUE)
  }
  # Clamp to valid probability range
  p_base <- pmin(pmax(p_base, 0.001), 0.999)

  # Generating propensity for the generate-treatment design: fitted once on
  # the real data, then treatment is redrawn from it on every replicate so the
  # simulated data satisfy positivity exactly to the extent the cohort does.
  # Keeping the observed treatment instead (design = "sample_treatment") makes
  # P(A = a | W) degenerate at the observed a, the artifact Shaw et al. (2025)
  # warn about.
  ps_fml_gen <- stats::reformulate(covariates, response = treatment)
  ps_mod_gen <- stats::glm(ps_fml_gen, data = data,
                           family = stats::binomial())
  ps_base <- pmin(pmax(as.numeric(
    stats::predict(ps_mod_gen, type = "response")), 0.001), 0.999)
  if (dgp_mode == "hybrid")
    .warn_hybrid_separation(.ps_cstat(ps_base, A),
                            "run_plasmode_feasibility")

  cand_ids <- vapply(tmle_candidates, function(x) x$candidate_id, character(1))

  all_results <- vector("list", length(effect_sizes))
  names(all_results) <- as.character(effect_sizes)

  for (es_idx in seq_along(effect_sizes)) {
    es <- effect_sizes[es_idx]
    if (verbose) message("Effect size: ", es)

    rep_results <- vector("list", reps)

    for (rep_i in seq_len(reps)) {
      withr::local_seed(lock$seed + rep_i)

      # Generate-treatment design (default): resample the covariate rows with
      # replacement and draw treatment from the fitted generating propensity.
      # Sample-treatment design (deprecated): keep the observed rows and the
      # observed treatment, simulating only the outcome.
      if (design == "generate_treatment") {
        idx   <- sample.int(n, n, replace = TRUE)
        A_rep <- stats::rbinom(n, 1L, ps_base[idx])
      } else {
        idx   <- seq_len(n)
        A_rep <- A
      }
      data_rep <- data[idx, , drop = FALSE]
      data_rep[[treatment]] <- A_rep

      # Simulated outcomes: additive risk difference es for treated.
      # Clamp BOTH bounds: with negative es and small p_base, an
      # unclamped p_base + es can go negative and stats::rbinom()
      # returns NA, causing every candidate fit to fail downstream.
      p1_sim <- pmin(pmax(p_base[idx] + es, 0.001), 0.999)
      p0_sim <- p_base[idx]
      p_obs  <- ifelse(A_rep == 1, p1_sim, p0_sim)
      Y_sim  <- stats::rbinom(n, 1L, p_obs)
      truth  <- mean(p1_sim) - mean(p0_sim)

      cand_results <- list()

      for (cand in tmle_candidates) {
        cand_result <- tryCatch({
          # Fit PS using the candidate's g-library
          g_lib <- cand$g_library
          q_lib <- cand$q_library

          use_sl <- requireNamespace("SuperLearner", quietly = TRUE) &&
            !identical(g_lib, "SL.glm")

          if (use_sl) {
            W_mat <- data_rep[, covariates, drop = FALSE]
            g_sl  <- SuperLearner::SuperLearner(
              Y = A_rep, X = W_mat, family = binomial(),
              SL.library = g_lib,
              env = .cleantmle_sl_env()
            )
            ps_hat <- as.numeric(g_sl$SL.predict)
          } else {
            ps_fml <- stats::reformulate(covariates, response = treatment)
            ps_mod <- stats::glm(ps_fml, data = data_rep,
                                 family = stats::binomial())
            ps_hat <- as.numeric(stats::predict(ps_mod, type = "response"))
          }
          ps_hat <- pmax(pmin(ps_hat, 1 - cand$truncation), cand$truncation)

          # Fit Q-model using the candidate's q-library
          ds <- data_rep
          ds[[".Y_sim."]] <- Y_sim
          AW <- ds[, c(treatment, covariates), drop = FALSE]

          use_sl_q <- requireNamespace("SuperLearner", quietly = TRUE) &&
            !identical(q_lib, "SL.glm")

          if (use_sl_q) {
            q_sl <- SuperLearner::SuperLearner(
              Y = Y_sim, X = AW, family = binomial(),
              SL.library = q_lib,
              env = .cleantmle_sl_env()
            )
            ds_a1 <- AW; ds_a1[[treatment]] <- 1L
            ds_a0 <- AW; ds_a0[[treatment]] <- 0L
            Q_a1 <- as.numeric(predict(q_sl, newdata = ds_a1)$pred)
            Q_a0 <- as.numeric(predict(q_sl, newdata = ds_a0)$pred)
            Q_aw <- as.numeric(q_sl$SL.predict)
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

          # TMLE targeting step
          H_a1 <- 1 / ps_hat
          H_a0 <- -1 / (1 - ps_hat)
          H_aw <- ifelse(A_rep == 1, H_a1, H_a0)

          epsilon <- tryCatch({
            Q_logit <- stats::qlogis(pmax(pmin(Q_aw, 0.999), 0.001))
            fluc <- stats::glm(
              Y_sim ~ -1 + H_aw + offset(Q_logit),
              family = stats::binomial()
            )
            unname(stats::coef(fluc))
          }, error = function(e) 0)

          Q_a1_u <- stats::plogis(
            stats::qlogis(pmax(pmin(Q_a1, 0.999), 0.001)) + epsilon * H_a1)
          Q_a0_u <- stats::plogis(
            stats::qlogis(pmax(pmin(Q_a0, 0.999), 0.001)) + epsilon * H_a0)
          Q_aw_u <- stats::plogis(
            stats::qlogis(pmax(pmin(Q_aw, 0.999), 0.001)) + epsilon * H_aw)

          est  <- mean(Q_a1_u) - mean(Q_a0_u)
          eic  <- H_aw * (Y_sim - Q_aw_u) + (Q_a1_u - Q_a0_u) - est
          se   <- sqrt(var(eic) / n)

          list(est      = est,
               se       = se,
               ci_lower = est - 1.96 * se,
               ci_upper = est + 1.96 * se)
        }, error = function(e) {
          list(est = NA_real_, se = NA_real_,
               ci_lower = NA_real_, ci_upper = NA_real_)
        })

        cand_results[[cand$candidate_id]] <- cand_result
      }

      rep_results[[rep_i]] <- c(cand_results, list(.truth = truth))
    }

    all_results[[as.character(es)]] <- rep_results
  }

  # Aggregate performance metrics per candidate
  metrics_rows <- lapply(effect_sizes, function(es) {
    rr <- all_results[[as.character(es)]]
    truth_v <- vapply(rr, function(x) x$.truth, numeric(1L))

    lapply(cand_ids, function(cid) {
      ests   <- vapply(rr, function(x) x[[cid]]$est,      numeric(1L))
      ses    <- vapply(rr, function(x) x[[cid]]$se,        numeric(1L))
      ci_los <- vapply(rr, function(x) x[[cid]]$ci_lower, numeric(1L))
      ci_his <- vapply(rr, function(x) x[[cid]]$ci_upper, numeric(1L))

      valid  <- !is.na(ests)
      ests_v <- ests[valid]; ses_v <- ses[valid]
      ci_lo_v <- ci_los[valid]; ci_hi_v <- ci_his[valid]
      truth_vv <- truth_v[valid]

      data.frame(
        effect_size = es,
        candidate   = cid,
        bias        = round(mean(ests_v - truth_vv), 5L),
        rmse        = round(sqrt(mean((ests_v - truth_vv)^2)), 5L),
        coverage    = round(mean(ci_lo_v <= truth_vv & truth_vv <= ci_hi_v), 3L),
        emp_sd      = round(sd(ests_v), 5L),
        mean_se     = round(mean(ses_v), 5L),
        n_converged = sum(valid),
        stringsAsFactors = FALSE
      )
    })
  })
  metrics <- do.call(rbind, unlist(metrics_rows, recursive = FALSE))

  # Add SE calibration ratio
  metrics$se_cal <- round(
    ifelse(metrics$emp_sd > 0, metrics$mean_se / metrics$emp_sd, NA_real_), 3L)

  # Surface poorly-converging candidates so the caller can see when
  # selection is operating on too few successful fits to be reliable.
  # A candidate that converged on < 50% of inner reps cannot be trusted
  # to inform candidate selection; warn the user explicitly.
  poor <- metrics[metrics$n_converged < ceiling(reps / 2), , drop = FALSE]
  if (nrow(poor) > 0L && isTRUE(verbose)) {
    msg <- paste0(
      "run_plasmode_feasibility: ", nrow(poor),
      " (effect_size x candidate) cell(s) had < 50% convergence; ",
      "selection on these candidates is unreliable. ",
      "Worst: ",
      paste(sprintf("%s@es=%g (n_converged=%d/%d)",
                     poor$candidate, poor$effect_size,
                     poor$n_converged, reps),
            collapse = "; "))
    warning(msg, call. = FALSE)
  }

  # The returned object deliberately carries no lock, no data, and no
  # fitted Q0: only per-replicate candidate metrics plus the lock
  # fingerprint and the generator mode (see the generator-modes section).
  result <- list(
    results         = all_results,
    metrics         = metrics,
    tmle_candidates = tmle_candidates,
    lock_hash       = lock$lock_hash,
    dgp_mode        = dgp_mode,
    effect_sizes    = effect_sizes,
    reps            = reps,
    design          = design,
    call            = match.call()
  )
  class(result) <- "plasmode_results"
  result
}


#' @export
print.plasmode_results <- function(x, ...) {
  cat("Plasmode-Simulation Feasibility Evaluation\n")
  cat("============================================\n")
  n_cands <- length(x$tmle_candidates)
  cat(sprintf("TMLE candidates: %d\n", n_cands))
  cat("Effect sizes evaluated:", paste(x$effect_sizes, collapse = ", "), "\n")
  cat("Replicates per effect size:", x$reps, "\n\n")
  cat("Performance Metrics:\n")
  print(x$metrics, row.names = FALSE)
  invisible(x)
}


#' Select Best TMLE Candidate Specification
#'
#' Applies a prespecified selection rule to plasmode-simulation performance
#' metrics to choose the best TMLE candidate specification.  Candidates
#' are prespecified TMLE implementations varying in truncation threshold
#' and/or nuisance-learner library.
#'
#' @section Clean-room stage: Stage 2b (pre-outcome).
#'
#' @param sim_results A `plasmode_results` object from
#'   [run_plasmode_feasibility()].
#' @param rule Character; selection criterion.  Default: \code{"min_rmse"}.
#'   Also supports \code{"min_bias"}, \code{"max_coverage"}, and
#'   \code{"min_max_rmse"} (minimax RMSE across DQ stress scenarios; requires
#'   \code{dq_results}).
#' @param thresholds Optional named list with \code{max_abs_bias},
#'   \code{min_coverage} for pre-filtering.  Candidates failing these
#'   thresholds are excluded before the selection rule is applied.
#'   If all candidates fail, the least-bad candidate is returned with
#'   a warning.
#' @param dq_results Optional \code{plasmode_dq_results} object from
#'   [run_plasmode_dq_stress()].  Required when \code{rule = "min_max_rmse"}:
#'   the rule then minimises the worst-case RMSE that each candidate exhibits
#'   across all DQ degraded scenarios (excluding the \code{"none"} baseline).
#'   Ignored for the other rules.
#' @param composite_weights Named numeric vector of weights for
#'   \code{rule = "composite"}: \code{rmse}, \code{coverage_penalty}, and
#'   \code{se_penalty}. Default \code{c(rmse = 1, coverage_penalty = 2,
#'   se_penalty = 0.5)}.
#' @param fiord_target_coverage Numeric; nominal coverage used by the
#'   \code{"ci_coverage"} and \code{"fiord_two_stage"} rules. Default 0.95.
#' @param fiord_coverage_tol Numeric; tolerance around
#'   \code{fiord_target_coverage} for the \code{"fiord_two_stage"} Stage-1
#'   oracle-coverage screen. Default 0.02.
#'
#' @return An object of class \code{tmle_selected_spec} containing the
#'   full candidate specification, the selection rule, and the metrics
#'   at selection.  Also inherits from \code{tmle_candidate_spec}.
#'   The \code{candidate_id} field can be extracted with
#'   \code{as.character()}.
#'
#' @examples
#' \dontrun{
#' plas <- run_plasmode_feasibility(lock, reps = 20L)
#' best <- select_tmle_candidate(plas, rule = "min_rmse")
#' print(best)
#' as.character(best)
#' }
#'
#' @keywords internal
select_tmle_candidate <- function(sim_results,
                                   rule = c("min_rmse", "min_bias",
                                            "max_coverage", "min_max_rmse",
                                            "ci_coverage", "min_se",
                                            "composite",
                                            "fiord_two_stage"),
                                   thresholds = NULL,
                                   dq_results = NULL,
                                   composite_weights = c(rmse = 1.0,
                                                         coverage_penalty = 2.0,
                                                         se_penalty = 0.5),
                                   fiord_target_coverage = 0.95,
                                   fiord_coverage_tol = 0.02) {
  if (!inherits(sim_results, "plasmode_results"))
    stop("`sim_results` must be a plasmode_results object.", call. = FALSE)
  rule <- match.arg(rule)
  if (rule == "min_max_rmse" && is.null(dq_results))
    stop("`dq_results` is required when rule = 'min_max_rmse'.",
         call. = FALSE)
  if (!is.null(dq_results) && !inherits(dq_results, "plasmode_dq_results"))
    stop("`dq_results` must be a plasmode_dq_results object.", call. = FALSE)

  m <- sim_results$metrics
  cands <- sim_results$tmle_candidates

  # Worst-case RMSE per candidate across DQ degraded scenarios
  worst_rmse <- NULL
  if (!is.null(dq_results)) {
    dqm <- dq_results$metrics
    dqm <- dqm[dqm$scenario != "none", ]
    if (nrow(dqm) > 0) {
      ids <- unique(dqm$candidate)
      worst_rmse <- vapply(ids, function(cid)
        max(dqm$rmse[dqm$candidate == cid], na.rm = TRUE),
        numeric(1))
      names(worst_rmse) <- ids
    }
  }

  # Average metrics across effect sizes per candidate
  cand_ids <- unique(m$candidate)
  summary_m <- do.call(rbind, lapply(cand_ids, function(cid) {
    sub <- m[m$candidate == cid, ]
    data.frame(
      candidate = cid,
      bias      = mean(abs(sub$bias)),
      rmse      = mean(sub$rmse),
      coverage  = mean(sub$coverage),
      emp_sd    = mean(sub$emp_sd),
      mean_se   = mean(sub$mean_se),
      max_rmse  = if (!is.null(worst_rmse) && cid %in% names(worst_rmse))
                    unname(worst_rmse[cid]) else NA_real_,
      stringsAsFactors = FALSE
    )
  }))

  # Pre-filter by thresholds if supplied
  eligible <- rep(TRUE, nrow(summary_m))
  if (!is.null(thresholds)) {
    if (!is.null(thresholds$max_abs_bias))
      eligible <- eligible & summary_m$bias <= thresholds$max_abs_bias
    if (!is.null(thresholds$min_coverage))
      eligible <- eligible & summary_m$coverage >= thresholds$min_coverage
  }

  all_failed <- !any(eligible)
  if (all_failed) {
    warning("No candidates passed thresholds; selecting least-bad candidate.",
            call. = FALSE)
    eligible <- rep(TRUE, nrow(summary_m))
  }

  pool <- summary_m[eligible, , drop = FALSE]

  # Apply selection rule
  best_idx <- switch(rule,
    min_rmse     = which.min(pool$rmse),
    min_bias     = which.min(pool$bias),
    max_coverage = which.max(pool$coverage),
    min_max_rmse = {
      if (all(is.na(pool$max_rmse)))
        stop("dq_results contains no candidates matching sim_results.",
             call. = FALSE)
      which.min(pool$max_rmse)
    },
    # New (workshop-driven; TODO A.23):
    ci_coverage = {
      # Pick the candidate whose coverage is closest to nominal.
      which.min(abs(pool$coverage - fiord_target_coverage))
    },
    min_se = which.min(pool$mean_se),
    composite = {
      # Weighted composite: standardised RMSE + coverage penalty +
      # SE-vs-empirical-SD penalty. Smaller is better.
      cov_pen <- (pool$coverage - fiord_target_coverage)^2
      se_pen  <- (pool$mean_se - pool$emp_sd)^2 /
                 pmax(pool$emp_sd, .Machine$double.eps)^2
      score <- composite_weights["rmse"] * scale(pool$rmse)[, 1L] +
               composite_weights["coverage_penalty"] * scale(cov_pen)[, 1L] +
               composite_weights["se_penalty"] * scale(se_pen)[, 1L]
      which.min(score)
    },
    fiord_two_stage = {
      # Stage 1: keep candidates whose oracle coverage is within
      # fiord_coverage_tol of the nominal target.
      keep <- abs(pool$coverage - fiord_target_coverage) <= fiord_coverage_tol
      if (!any(keep)) {
        warning("FIORD stage 1: no candidate achieves nominal coverage; ",
                "falling back to closest-coverage candidate.",
                call. = FALSE)
        return_idx <- which.min(abs(pool$coverage - fiord_target_coverage))
        return_idx
      } else {
        # Stage 2: among those, choose the candidate with the
        # smallest mean SE (most precise inference at the right
        # coverage). NOTE: variance-method selection per FIORD is a
        # separate step that needs candidates run under each variance
        # method; this rule selects the point-estimator candidate.
        stage2 <- pool[keep, , drop = FALSE]
        stage2_idx <- which.min(stage2$mean_se)
        # Map back to the row of pool
        which(keep)[stage2_idx]
      }
    }
  )
  best_id <- pool$candidate[best_idx]

  # Build selected-spec object from the matching candidate
  cand_match <- NULL
  for (c in cands) {
    if (c$candidate_id == best_id) { cand_match <- c; break }
  }
  if (is.null(cand_match))
    stop("Internal error: selected candidate not found in list.", call. = FALSE)

  best_metrics <- pool[best_idx, , drop = FALSE]

  result <- c(cand_match,
    list(
      selection_rule   = rule,
      metrics          = best_metrics,
      thresholds_used  = thresholds,
      all_failed       = all_failed,
      lock_hash        = sim_results$lock_hash %||%
                           (if (!is.null(sim_results$lock))
                              sim_results$lock$lock_hash else NA_character_)
    )
  )
  class(result) <- c("tmle_selected_spec", "tmle_candidate_spec")
  message("Selected TMLE candidate: '", best_id, "' (rule = '", rule, "')")
  result
}


#' @export
print.tmle_selected_spec <- function(x, ...) {
  cat(sprintf("Selected TMLE Candidate: %s\n", x$candidate_id))
  cat(sprintf("  Label:      %s\n", x$label))
  cat(sprintf("  g-library:  %s\n", paste(x$g_library, collapse = ", ")))
  cat(sprintf("  Q-library:  %s\n", paste(x$q_library, collapse = ", ")))
  cat(sprintf("  Truncation: %s\n", x$truncation))
  cat(sprintf("  Rule:       %s\n", x$selection_rule))
  if (!is.null(x$metrics)) {
    cat(sprintf("  RMSE:       %.5f\n", x$metrics$rmse))
    cat(sprintf("  Bias:       %.5f\n", x$metrics$bias))
    cat(sprintf("  Coverage:   %.3f\n", x$metrics$coverage))
  }
  invisible(x)
}


#' @export
as.character.tmle_selected_spec <- function(x, ...) {
  x$candidate_id
}


# ── Stage 3: Conventional Workflows ──────────────────────────────────────

#' Run Propensity-Score Matching Workflow
#'
#' Performs greedy 1:1 nearest-neighbor matching on the logit of the
#' propensity score, then estimates the average treatment effect in the
#' matched sample as the simple risk difference.
#'
#' @param lock A `cleanroom_lock` from [create_analysis_lock()].
#' @param ps_fit A `ps_fit` object from [fit_ps_superlearner()].
#' @param caliper Numeric; maximum PS distance (logit scale) allowed for a
#'   match. Default: `0.2 * sd(logit(ps))`.
#' @param allow_outcome_access Logical; if \code{TRUE}, skips the
#'   outcome-access check. Default \code{FALSE}.
#' @param override_clean_room Deprecated. Use \code{allow_outcome_access}.
#'
#' @return An object of class `match_result` containing the causal risk
#'   difference estimate, SE, 95% CI, p-value, and the matched dataset.
#'
#' @section Variance:
#' The paired-difference SE in the returned object conditions on the fixed
#' matched dataset and does not propagate matching-draw variance -- the
#' variability from which controls are selected when re-sampling. Abadie and
#' Imbens (2008) show that the standard paired SE is inconsistent for this
#' reason; Abadie and Imbens (2016) show that matching on the estimated
#' propensity score does not restore consistency of the naive bootstrap.
#' The principled solution is a full-pipeline nonparametric bootstrap that
#' re-runs PS estimation, matching, and TMLE on each resample.
#' Use \code{bootstrap_rd_variance(estimator = "match_tmle")} to obtain a
#' bootstrap SE that captures matching-draw variance.
#' \code{select_variance_method(estimator = "match_tmle")} can confirm
#' which variance method achieves nominal oracle coverage for the study's
#' data-generating process and sample size.
#'
#' @references
#' Abadie, A. and Imbens, G. W. (2008). On the failure of the bootstrap for
#' matching estimators. \emph{Econometrica}, 76(6), 1537--1557.
#'
#' Abadie, A. and Imbens, G. W. (2016). Matching on the estimated propensity
#' score. \emph{Econometrica}, 84(2), 781--807.
#'
#' @keywords internal
run_match_workflow <- function(lock, ps_fit, caliper = NULL,
                               allow_outcome_access = FALSE,
                               override_clean_room = NULL) {
  if (!is.null(override_clean_room)) {
    rlang::warn(
      "override_clean_room is deprecated in run_match_workflow(); use allow_outcome_access.",
      .frequency = "once", .frequency_id = "run_match_workflow_deprecated"
    )
    allow_outcome_access <- override_clean_room
  }
  if (!inherits(lock, "cleanroom_lock"))
    stop("`lock` must be a cleanroom_lock object.", call. = FALSE)
  if (!inherits(ps_fit, "ps_fit"))
    stop("`ps_fit` must be a ps_fit object.", call. = FALSE)
  .check_outcome_access(lock, allow_outcome_access,
                        caller = "run_match_workflow")

  data      <- .join_outcome(lock)
  treatment <- lock$treatment
  outcome   <- lock$outcome
  A         <- data[[treatment]]
  Y         <- data[[outcome]]
  ps        <- ps_fit$ps

  logit_ps <- log(ps / (1 - ps))
  if (is.null(caliper))
    caliper <- 0.2 * sd(logit_ps)

  treated_idx  <- which(A == 1)
  control_idx  <- which(A == 0)
  n_treated    <- length(treated_idx)
  used_control <- logical(length(control_idx))
  matched_ctrl <- integer(n_treated)

  for (i in seq_len(n_treated)) {
    dists             <- abs(logit_ps[treated_idx[i]] - logit_ps[control_idx])
    dists[used_control] <- Inf
    best              <- which.min(dists)
    if (dists[best] <= caliper) {
      matched_ctrl[i]    <- control_idx[best]
      used_control[best] <- TRUE
    } else {
      matched_ctrl[i] <- NA_integer_
    }
  }

  valid        <- !is.na(matched_ctrl)
  match_idx    <- c(treated_idx[valid], matched_ctrl[valid])
  matched_data <- data[match_idx, , drop = FALSE]

  A_m  <- matched_data[[treatment]]
  Y_m  <- matched_data[[outcome]]
  n_m  <- sum(valid)

  Y_m <- as.numeric(Y_m)

  # Complete-case handling (NA outcomes drop pairs; warn).
  if (any(is.na(Y_m))) {
    warning("run_match_workflow: ", sum(is.na(Y_m)), " of ", length(Y_m),
            " matched-cohort outcome rows are NA; computing matched RD on ",
            "complete cases.", call. = FALSE)
  }

  n1 <- sum(A_m == 1 & !is.na(Y_m))
  n0 <- sum(A_m == 0 & !is.na(Y_m))
  r1   <- mean(Y_m[A_m == 1], na.rm = TRUE)
  r0   <- mean(Y_m[A_m == 0], na.rm = TRUE)
  rd   <- r1 - r0

  # Paired SE on complete pairs (drop pair if either side is NA)
  Y1_full <- Y_m[A_m == 1]
  Y0_full <- Y_m[A_m == 0][seq_len(n_m)]
  pair_ok <- !is.na(Y1_full) & !is.na(Y0_full)
  if (sum(pair_ok) > 1L) {
    se <- sqrt(stats::var(Y1_full[pair_ok] - Y0_full[pair_ok]) / sum(pair_ok))
  } else {
    se <- NA_real_
  }
  ci_lo <- rd - 1.96 * se
  ci_hi <- rd + 1.96 * se
  p_val <- 2 * stats::pnorm(-abs(rd / se))

  result <- list(
    estimate     = rd,
    se           = se,
    ci_lower     = ci_lo,
    ci_upper     = ci_hi,
    p_value      = p_val,
    r1           = r1,
    r0           = r0,
    n_matched    = n_m,
    n_unmatched  = sum(!valid),
    matched_data = matched_data,
    treatment    = treatment,
    outcome      = outcome,
    call         = match.call()
  )
  class(result) <- c("match_result", "cr_result")
  result
}


#' @export
print.match_result <- function(x, ...) {
  cat("Propensity Score Matching Workflow\n")
  cat("===================================\n")
  cat(sprintf("Matched pairs:   %d\n", x$n_matched))
  if (x$n_unmatched > 0L)
    cat(sprintf("Unmatched (caliper): %d\n", x$n_unmatched))
  cat(sprintf("Risk (treated):  %.4f\n", x$r1))
  cat(sprintf("Risk (control):  %.4f\n", x$r0))
  cat(sprintf("Risk Difference: %.4f  (95%% CI: %.4f, %.4f)\n",
              x$estimate, x$ci_lower, x$ci_upper))
  cat(sprintf("SE:              %.5f   p-value: %.4f\n", x$se, x$p_value))
  invisible(x)
}


#' Run IPTW Workflow
#'
#' Computes stabilized inverse probability of treatment weights from the
#' fitted propensity score and estimates the causal risk difference using
#' the Hajek (normalized) estimator with influence-curve-based variance.
#'
#' @param lock A `cleanroom_lock` from [create_analysis_lock()].
#' @param ps_fit A `ps_fit` object from [fit_ps_superlearner()].
#' @param trim Quantile for weight trimming (e.g., `0.01` trims at 1st and
#'   99th percentile). Default: `NULL` (no trimming).
#' @param allow_outcome_access Logical; if \code{TRUE}, skips the outcome-access check. Default \code{FALSE}.
#' @param override_clean_room Deprecated. Use \code{allow_outcome_access}.
#'
#' @return An object of class `iptw_result` containing the estimated risk
#'   difference, SE, 95% CI, p-value, and IPTW weights.
#'
#' @keywords internal
run_iptw_workflow <- function(lock, ps_fit, trim = NULL,
                              allow_outcome_access = FALSE,
                              override_clean_room = NULL) {
  if (!is.null(override_clean_room)) {
    rlang::warn(
      "override_clean_room is deprecated in run_iptw_workflow(); use allow_outcome_access.",
      .frequency = "once", .frequency_id = "run_iptw_workflow_deprecated"
    )
    allow_outcome_access <- override_clean_room
  }
  if (!inherits(lock, "cleanroom_lock"))
    stop("`lock` must be a cleanroom_lock object.", call. = FALSE)
  if (!inherits(ps_fit, "ps_fit"))
    stop("`ps_fit` must be a ps_fit object.", call. = FALSE)
  .check_outcome_access(lock, allow_outcome_access,
                        caller = "run_iptw_workflow")

  data      <- .join_outcome(lock)
  treatment <- lock$treatment
  outcome   <- lock$outcome
  A         <- data[[treatment]]
  Y         <- data[[outcome]]
  ps        <- ps_fit$ps
  n         <- nrow(data)

  # Outcome NA handling (complete-case IPTW; warn; recommend IPCW for MAR)
  na_y <- is.na(Y)
  if (any(na_y)) {
    warning("run_iptw_workflow: ", sum(na_y), " of ", n,
            " outcome rows are NA; computing complete-case IPTW. ",
            "Inference is valid under MCAR only -- consider IPCW for MAR.",
            call. = FALSE)
  }
  keep <- !na_y
  A_k  <- A[keep]; Y_k <- Y[keep]; ps_k <- ps[keep]
  n_eff <- sum(keep)

  # Stabilized IPTW weights
  p_trt <- mean(A_k)
  w     <- ifelse(A_k == 1, p_trt / ps_k, (1 - p_trt) / (1 - ps_k))

  if (!is.null(trim)) {
    lo <- stats::quantile(w, trim)
    hi <- stats::quantile(w, 1 - trim)
    w  <- pmin(pmax(w, lo), hi)
  }

  # Hajek estimator
  r1 <- stats::weighted.mean(Y_k[A_k == 1], w[A_k == 1])
  r0 <- stats::weighted.mean(Y_k[A_k == 0], w[A_k == 0])
  rd <- r1 - r0

  # Linearized variance for Hajek estimator
  se_sq_1 <- sum(w[A_k == 1]^2 * (Y_k[A_k == 1] - r1)^2) / sum(w[A_k == 1])^2
  se_sq_0 <- sum(w[A_k == 0]^2 * (Y_k[A_k == 0] - r0)^2) / sum(w[A_k == 0])^2
  se      <- sqrt(se_sq_1 + se_sq_0)

  ci_lo <- rd - 1.96 * se
  ci_hi <- rd + 1.96 * se
  p_val <- 2 * stats::pnorm(-abs(rd / se))

  # Re-pad weights to length n with NA at NA-Y rows for downstream alignment
  w_full <- rep(NA_real_, n)
  w_full[keep] <- w

  result <- list(
    estimate  = rd,
    se        = se,
    ci_lower  = ci_lo,
    ci_upper  = ci_hi,
    p_value   = p_val,
    r1        = r1,
    r0        = r0,
    weights   = w_full,
    ps        = ps,
    n         = n,
    n_effective = n_eff,
    treatment = treatment,
    outcome   = outcome,
    call      = match.call()
  )
  class(result) <- c("iptw_result", "cr_result")
  result
}


#' Run IPCW-Weighted TMLE for a Binary Outcome with Missing-At-Random Y
#'
#' Targets the marginal risk difference on the full cohort when the
#' outcome `Y` has missing values, using inverse-probability-of-
#' censoring weighting to correct for missing-at-random follow-up.
#' Internally, a response model
#' \eqn{P(\text{observed} \mid A, W)} is fit by SuperLearner (or GLM if
#' SuperLearner is unavailable), stabilised inverse-probability weights
#' are constructed, and the cohort is passed to `tmle::tmle()` with the
#' \code{Delta} argument so that the package's targeting step uses the
#' censoring weights internally rather than dropping incomplete rows.
#'
#' @section Clean-room stage: Stage 4 (accesses the outcome).
#'
#' @param lock A `cleanroom_lock` from [create_analysis_lock()] with the
#'   outcome unmasked.
#' @param ps_fit Optional `ps_fit` from [fit_ps_superlearner()]. If
#'   \code{NULL} (default), the censoring model still uses the locked
#'   SuperLearner library; the treatment mechanism is re-fit by
#'   `tmle::tmle()`.
#' @param censoring_library SuperLearner library used for the censoring
#'   model `g(Delta = 1 | A, W)`. Default \code{NULL} uses the lock's
#'   SuperLearner library.
#' @param weight_truncation Upper-quantile cap for the IPCW weights to
#'   bound extreme values. Default \code{0.99}.
#' @param allow_outcome_access Logical; if \code{TRUE}, skips the outcome-access check. Default \code{FALSE}.
#' @param override_clean_room Deprecated. Use \code{allow_outcome_access}.
#'
#' @return A list of class `tmle_fit` (and `cr_result`) containing the
#'   ATE estimate, IPCW summary, and the underlying `tmle::tmle()` fit.
#'   The returned object is compatible with [summarize_cleanroom_results()]
#'   and [forest_plot()].
#'
#' @examples
#' \dontrun{
#' lock_unmasked <- unmask_outcome(lock, lock_pre_mask)
#' ps_fit        <- fit_ps_superlearner(lock_unmasked)
#' ipcw_fit      <- run_ipcw_tmle(lock_unmasked, ps_fit)
#' print(ipcw_fit)
#' }
#'
#' @keywords internal
run_ipcw_tmle <- function(lock, ps_fit = NULL,
                          censoring_library = NULL,
                          weight_truncation = 0.99,
                          allow_outcome_access = FALSE,
                          override_clean_room = NULL) {
  if (!is.null(override_clean_room)) {
    rlang::warn(
      "override_clean_room is deprecated in run_ipcw_tmle(); use allow_outcome_access.",
      .frequency = "once", .frequency_id = "run_ipcw_tmle_deprecated"
    )
    allow_outcome_access <- override_clean_room
  }
  if (!inherits(lock, "cleanroom_lock"))
    stop("`lock` must be a cleanroom_lock object.", call. = FALSE)
  .check_outcome_access(lock, allow_outcome_access,
                        caller = "run_ipcw_tmle")
  if (!requireNamespace("tmle", quietly = TRUE))
    stop("Package 'tmle' is required for run_ipcw_tmle(). ",
         "Install with: install.packages('tmle')", call. = FALSE)

  data       <- .join_outcome(lock)
  treatment  <- lock$treatment
  outcome    <- lock$outcome
  covariates <- lock$covariates
  A          <- as.numeric(data[[treatment]])
  Y          <- as.numeric(data[[outcome]])
  W          <- data[, covariates, drop = FALSE]
  n          <- nrow(data)

  R <- as.integer(!is.na(Y))
  if (sum(R) < length(covariates) + 2L)
    stop("run_ipcw_tmle: too few non-NA outcomes (", sum(R),
         ") to fit Q with ", length(covariates), " covariates.",
         call. = FALSE)
  if (mean(R) == 1)
    message("run_ipcw_tmle: outcome has no missing values; ",
            "IPCW weights collapse to 1 -- consider plain TMLE instead.")

  cens_lib <- if (!is.null(censoring_library)) censoring_library
              else lock$sl_library
  q_lib    <- if (!is.null(lock$primary_tmle_spec))
                lock$primary_tmle_spec$q_library else lock$sl_library
  g_lib    <- if (!is.null(lock$primary_tmle_spec))
                lock$primary_tmle_spec$g_library else lock$sl_library

  # Stabilised IPCW weights using a SuperLearner response model.
  resp_X <- data.frame(A = A, W)
  resp_pred <- tryCatch({
    if (requireNamespace("SuperLearner", quietly = TRUE)) {
      withr::local_seed(lock$seed + 99L)
      sl <- SuperLearner::SuperLearner(
        Y = R, X = resp_X, family = stats::binomial(),
        SL.library = cens_lib,
        env = .cleantmle_sl_env())
      as.numeric(sl$SL.predict)
    } else {
      glm_resp <- stats::glm(R ~ ., data = resp_X, family = stats::binomial())
      as.numeric(stats::predict(glm_resp, type = "response"))
    }
  }, error = function(e) {
    glm_resp <- stats::glm(R ~ ., data = resp_X, family = stats::binomial())
    as.numeric(stats::predict(glm_resp, type = "response"))
  })
  pr_marginal <- tapply(R, A, mean)[as.character(A)]
  ipcw <- as.numeric(pr_marginal) / pmax(resp_pred, 0.01)
  ipcw_cap <- stats::quantile(ipcw, weight_truncation, names = FALSE,
                               na.rm = TRUE)
  ipcw <- pmin(ipcw, ipcw_cap)

  # Run tmle::tmle with Delta = R so the package handles the censoring
  # adjustment internally on the full cohort (rather than us dropping NA
  # rows). When tmle()'s Delta path fails, fall back to a complete-case
  # TMLE on the SuperLearner-Q-fit weighted by the IPCW.
  # tmle::tmle's Delta argument tells it to fit a censoring model and use
  # the resulting weights internally. The censoring SL library defaults
  # to g.SL.library on tmle versions that accept it via that path. Older
  # tmle versions reject Delta entirely; in that case we drop to a
  # complete-case TMLE on the SuperLearner Q-fit weighted by the IPCW
  # we computed above.
  method_used <- "tmle_delta_full_cohort"
  fit <- tryCatch({
    tmle::tmle(
      Y = ifelse(is.na(Y), 0L, Y),
      A = A,
      W = as.data.frame(W),
      family = "binomial",
      Q.SL.library = q_lib,
      g.SL.library = g_lib,
      Delta = R,
      verbose = FALSE
    )
  }, error = function(e) {
    method_used <<- "ipcw_weighted_complete_case"
    warning("run_ipcw_tmle: tmle::tmle Delta-path failed (",
            conditionMessage(e), "). Falling back to complete-case ",
            "TMLE weighted by IPCW. The analysis no longer uses the ",
            "full cohort; n_effective = ", sum(R), " (of ", n, "). ",
            "Note that this fallback gives up the joint double-robustness ",
            "in the censoring + outcome models that the Delta-path ",
            "provides: consistency now requires the IPCW SuperLearner ",
            "response model to be correctly specified. Inspect ",
            "method_used and the IPCW weight summary on the returned ",
            "object before reporting.",
            call. = FALSE)
    cc      <- which(R == 1L)
    Y_cc    <- Y[cc]; A_cc <- A[cc]
    W_cc    <- as.data.frame(W[cc, , drop = FALSE])
    w_cc    <- ipcw[cc]
    tmle::tmle(
      Y = Y_cc, A = A_cc, W = W_cc, family = "binomial",
      Q.SL.library = q_lib, g.SL.library = g_lib,
      obsWeights = w_cc, verbose = FALSE)
  })

  ate <- fit$estimates$ATE
  result <- list(
    estimate    = unname(ate$psi),
    se          = unname(sqrt(ate$var.psi)),
    ci_lower    = unname(ate$CI[1]),
    ci_upper    = unname(ate$CI[2]),
    p_value     = unname(ate$pvalue),
    estimates   = list(ATE = list(
      estimate = unname(ate$psi),
      se       = unname(sqrt(ate$var.psi)),
      ci_lower = unname(ate$CI[1]),
      ci_upper = unname(ate$CI[2]),
      p_value  = unname(ate$pvalue))),
    risk_treated = unname(fit$estimates$EY1$psi),
    risk_control = unname(fit$estimates$EY0$psi),
    n            = n,
    n_observed   = sum(R),
    n_missing    = sum(R == 0L),
    pct_missing  = round(100 * mean(R == 0L), 2),
    ipcw         = ipcw,
    ipcw_summary = list(mean = mean(ipcw), sd = stats::sd(ipcw),
                         min = min(ipcw),  max = max(ipcw)),
    tmle_fit     = fit,
    treatment    = treatment,
    outcome      = outcome,
    type         = "ipcw_tmle",
    method_used  = method_used,
    call         = match.call()
  )
  class(result) <- c("tmle_fit", "cr_result")
  result
}


#' @export
print.iptw_result <- function(x, ...) {
  cat("IPTW Workflow\n")
  cat("=============\n")
  cat(sprintf("N:               %d\n", x$n))
  cat(sprintf("Risk (treated):  %.4f\n", x$r1))
  cat(sprintf("Risk (control):  %.4f\n", x$r0))
  cat(sprintf("Risk Difference: %.4f  (95%% CI: %.4f, %.4f)\n",
              x$estimate, x$ci_lower, x$ci_upper))
  cat(sprintf("SE:              %.5f   p-value: %.4f\n", x$se, x$p_value))
  invisible(x)
}


# ── Stage 3: Modular TMLE ────────────────────────────────────────────────

#' Fit TMLE Treatment Mechanism
#'
#' Fits the treatment mechanism g(A|W) = P(A=1|W) for use in the modular
#' TMLE workflow. Uses propensity scores from an existing `ps_fit` if
#' provided; otherwise re-estimates via SuperLearner. The outcome is never
#' accessed, making this step safe at any stage.
#'
#' If the lock contains a locked primary TMLE specification (from
#' [lock_primary_tmle_spec()]), the PS truncation threshold from that
#' specification is applied automatically.
#'
#' @section Clean-room stage: Stage 2 / 4 (pre-outcome; does not access
#'   the outcome).
#'
#' @param lock A `cleanroom_lock` from [create_analysis_lock()].
#' @param ps_fit Optional `ps_fit` from [fit_ps_superlearner()]. If provided,
#'   the already-estimated propensity scores are reused.
#' @param truncation Numeric; PS truncation threshold.  If \code{NULL}
#'   (default), uses the locked primary TMLE spec truncation if available,
#'   otherwise 0.01.
#'
#' @return An object of class `tmle_mechanism` with `type = "treatment"`.
#'
#' @param n_folds Integer; number of cross-fitting folds. Default: 1
#'   (no cross-fitting). When > 1, PS is estimated via out-of-fold
#'   SuperLearner predictions.
#' @param fold_vec Optional integer vector assigning each observation
#'   to a fold.  Overrides \code{n_folds}.
#'
#' @keywords internal
fit_tmle_treatment_mechanism <- function(lock, ps_fit = NULL,
                                          truncation = NULL,
                                          n_folds = 1L,
                                          fold_vec = NULL) {
  if (!inherits(lock, "cleanroom_lock"))
    stop("`lock` must be a cleanroom_lock object.", call. = FALSE)

  # Resolve truncation from locked spec
  primary_spec <- lock$primary_tmle_spec
  if (is.null(truncation)) {
    truncation <- if (!is.null(primary_spec)) primary_spec$truncation else 0.01
  }

  data       <- lock$data
  A          <- data[[lock$treatment]]
  W          <- data[, lock$covariates, drop = FALSE]
  n          <- nrow(data)

  # Named truncation rule support (TODO A.21.6): if `truncation` is a
  # string like "sqrt_n_ln_n", resolve it now that n is known.
  if (is.character(truncation))
    truncation <- resolve_truncation_rule(truncation, n = n)

  # CV scheme: if the candidate spec sets cv_scheme = "cv_tmle" and the
  # caller did not override n_folds, use a sensible default V from the
  # Phillips rule.
  if (!is.null(primary_spec) &&
      identical(primary_spec$cv_scheme, "cv_tmle") &&
      n_folds == 1L && is.null(fold_vec)) {
    n_folds <- if (!is.null(primary_spec$cv_V)) primary_spec$cv_V
               else recommend_cv_V(compute_n_eff(.outcome_vector(lock),
                                                 family = "binomial"))
  }
  use_cv     <- !is.null(fold_vec) || n_folds > 1L

  if (use_cv && is.null(ps_fit)) {
    # Cross-fitted PS estimation
    if (is.null(fold_vec))
      fold_vec <- sample(rep(seq_len(n_folds), length.out = n))
    K <- max(fold_vec)

    if (!requireNamespace("SuperLearner", quietly = TRUE))
      stop("Package 'SuperLearner' is required for cross-fitting.",
           call. = FALSE)

    ps <- numeric(n)
    for (k in seq_len(K)) {
      val_idx   <- which(fold_vec == k)
      train_idx <- which(fold_vec != k)

      withr::local_seed(lock$seed + k)
      g_sl <- SuperLearner::SuperLearner(
        Y = A[train_idx], X = W[train_idx, , drop = FALSE],
        family = binomial(), SL.library = lock$sl_library,
        env = .cleantmle_sl_env()
      )
      ps[val_idx] <- as.numeric(
        predict(g_sl, newdata = W[val_idx, , drop = FALSE])$pred)
    }
    ps <- pmax(pmin(ps, 1 - truncation), truncation)
    g_mod <- NULL

  } else if (!is.null(ps_fit)) {
    if (!inherits(ps_fit, "ps_fit"))
      stop("`ps_fit` must be a ps_fit object.", call. = FALSE)
    ps    <- ps_fit$ps
    g_mod <- ps_fit$sl_fit
    ps <- pmax(pmin(ps, 1 - truncation), truncation)
  } else {
    temp  <- fit_ps_superlearner(lock)
    ps    <- temp$ps
    g_mod <- temp$sl_fit
    ps <- pmax(pmin(ps, 1 - truncation), truncation)
  }

  result <- list(
    type       = "treatment",
    ps         = ps,
    g_fit      = g_mod,
    treatment  = lock$treatment,
    covariates = lock$covariates,
    data       = lock$data,
    lock       = lock,
    truncation = truncation,
    cross_fitted = use_cv,
    n_folds    = if (use_cv) max(fold_vec) else 1L,
    fold_vec   = fold_vec,
    call       = match.call()
  )
  class(result) <- "tmle_mechanism"
  result
}


#' Fit TMLE Outcome Mechanism
#'
#' Fits the outcome mechanism Q(A,W) = E\[Y|A,W\] using SuperLearner (or
#' logistic regression as a fallback). This step accesses the outcome and
#' **must only be called in Stage 4** (after outcome unblinding).
#'
#' If the lock contains a locked primary TMLE specification (from
#' [lock_primary_tmle_spec()]), the Q-library from that specification is
#' used by default.
#'
#' @section Clean-room stage: Stage 4 (accesses the real outcome).
#'
#' @param lock A `cleanroom_lock` from [create_analysis_lock()].
#' @param g_fit A `tmle_mechanism` object with `type = "treatment"` from
#'   [fit_tmle_treatment_mechanism()].
#' @param sl_library Optional SuperLearner library override. If \code{NULL},
#'   uses the locked primary TMLE spec Q-library if available, then falls
#'   back to `lock$sl_library`.
#' @param allow_outcome_access Logical; if \code{TRUE}, skips the outcome-access check. Default \code{FALSE}.
#' @param override_clean_room Deprecated. Use \code{allow_outcome_access}.
#'
#' @return An object of class `tmle_mechanism` with `type = "outcome"`
#'   containing initial outcome predictions `Q_a1`, `Q_a0`, and `Q_aw`.
#'
#' @keywords internal
fit_tmle_outcome_mechanism <- function(lock, g_fit, sl_library = NULL,
                                       allow_outcome_access = FALSE,
                                       override_clean_room = NULL) {
  if (!is.null(override_clean_room)) {
    rlang::warn(
      "override_clean_room is deprecated in fit_tmle_outcome_mechanism(); use allow_outcome_access.",
      .frequency = "once", .frequency_id = "fit_tmle_outcome_mechanism_deprecated"
    )
    allow_outcome_access <- override_clean_room
  }
  if (!inherits(lock, "cleanroom_lock"))
    stop("`lock` must be a cleanroom_lock object.", call. = FALSE)
  if (!inherits(g_fit, "tmle_mechanism") || g_fit$type != "treatment")
    stop("`g_fit` must be a tmle_mechanism of type 'treatment'.",
         call. = FALSE)
  .check_outcome_access(lock, allow_outcome_access,
                        caller = "fit_tmle_outcome_mechanism")

  # Resolve Q-library from locked spec
  primary_spec <- lock$primary_tmle_spec
  if (is.null(sl_library)) {
    if (!is.null(primary_spec)) {
      sl_library <- primary_spec$q_library
    } else {
      sl_library <- lock$sl_library
    }
  }

  data       <- .join_outcome(lock)
  treatment  <- lock$treatment
  outcome    <- lock$outcome
  covariates <- lock$covariates
  AW         <- data[, c(treatment, covariates), drop = FALSE]
  Y          <- data[[outcome]]
  n          <- nrow(data)

  # SuperLearner requires Y to be a numeric vector. Locks built from real
  # data often store the outcome as integer or labelled-integer (haven /
  # data.table conventions); coerce explicitly.
  if (is.factor(Y) || is.character(Y)) Y <- as.numeric(as.character(Y))
  Y <- as.numeric(Y)

  # Outcome NA handling: SuperLearner refuses Y with NAs, but real-world
  # observational data routinely has missing outcomes (loss to follow-up).
  # We fit Q on complete-outcome rows and predict for ALL rows so the
  # downstream targeting step sees a length-n Q vector. Targeting itself
  # uses Y which may have NAs; the fluctuation glm drops them via na.action.
  # This is a complete-case Q fit and is unbiased only under MCAR; users
  # who suspect MAR/MNAR should weight by an inverse-probability-of-
  # censoring model (see vignette).
  na_y <- is.na(Y)
  n_na_y <- sum(na_y)
  if (n_na_y > 0L) {
    warning("fit_tmle_outcome_mechanism: ", n_na_y, " of ", n,
            " outcome rows are NA; fitting Q on complete cases. ",
            "Inference is valid under MCAR only -- consider IPCW for MAR.",
            call. = FALSE)
  }
  fit_idx <- which(!na_y)

  # Inherit cross-fitting from g_fit if present
  use_cv   <- isTRUE(g_fit$cross_fitted)
  fold_vec <- g_fit$fold_vec

  if (use_cv && !is.null(fold_vec) &&
      requireNamespace("SuperLearner", quietly = TRUE)) {
    # Cross-fitted Q estimation
    K <- max(fold_vec)
    Q_a1 <- numeric(n)
    Q_a0 <- numeric(n)
    Q_aw <- numeric(n)

    for (k in seq_len(K)) {
      val_idx   <- which(fold_vec == k)
      train_idx <- intersect(which(fold_vec != k), fit_idx)
      if (length(train_idx) == 0L) next

      withr::local_seed(lock$seed + 1L + k)
      Q_sl_k <- SuperLearner::SuperLearner(
        Y          = Y[train_idx],
        X          = AW[train_idx, , drop = FALSE],
        family     = binomial(),
        SL.library = sl_library,
        env        = .cleantmle_sl_env()
      )

      AW_val <- AW[val_idx, , drop = FALSE]
      Q_aw[val_idx] <- as.numeric(predict(Q_sl_k, newdata = AW_val)$pred)

      AW_a1 <- AW_val; AW_a1[[treatment]] <- 1L
      AW_a0 <- AW_val; AW_a0[[treatment]] <- 0L
      Q_a1[val_idx] <- as.numeric(predict(Q_sl_k, newdata = AW_a1)$pred)
      Q_a0[val_idx] <- as.numeric(predict(Q_sl_k, newdata = AW_a0)$pred)
    }
    Q_fit_o <- NULL

  } else if (requireNamespace("SuperLearner", quietly = TRUE)) {
    withr::local_seed(lock$seed + 1L)
    Q_sl <- SuperLearner::SuperLearner(
      Y          = Y[fit_idx],
      X          = AW[fit_idx, , drop = FALSE],
      family     = binomial(),
      SL.library = sl_library,
      env        = .cleantmle_sl_env()
    )
    AW_a1   <- AW; AW_a1[[treatment]] <- 1L
    AW_a0   <- AW; AW_a0[[treatment]] <- 0L
    Q_a1    <- as.numeric(predict(Q_sl, newdata = AW_a1)$pred)
    Q_a0    <- as.numeric(predict(Q_sl, newdata = AW_a0)$pred)
    Q_aw    <- as.numeric(predict(Q_sl, newdata = AW)$pred)
    Q_fit_o <- Q_sl
  } else {
    Q_fml   <- stats::reformulate(c(treatment, covariates), response = outcome)
    Q_glm   <- stats::glm(Q_fml, data = data, family = stats::binomial(),
                          na.action = stats::na.exclude)
    da1     <- data; da1[[treatment]] <- 1L
    da0     <- data; da0[[treatment]] <- 0L
    Q_a1    <- as.numeric(stats::predict(Q_glm, newdata = da1, type = "response"))
    Q_a0    <- as.numeric(stats::predict(Q_glm, newdata = da0, type = "response"))
    Q_aw    <- as.numeric(stats::predict(Q_glm, newdata = data, type = "response"))
    Q_fit_o <- Q_glm
  }

  result <- list(
    type         = "outcome",
    Q_a1         = Q_a1,
    Q_a0         = Q_a0,
    Q_aw         = Q_aw,
    Q_fit        = Q_fit_o,
    outcome      = outcome,
    treatment    = treatment,
    covariates   = covariates,
    data         = data,
    lock         = lock,
    g_fit        = g_fit,
    cross_fitted = use_cv,
    fold_vec     = fold_vec,
    call         = match.call()
  )
  class(result) <- "tmle_mechanism"
  result
}


#' Run TMLE Targeting Step
#'
#' Performs the TMLE fluctuation (targeting) step by fitting a logistic
#' submodel through the initial outcome estimates using the clever covariate
#' H(A,W) = A/g(W) - (1-A)/(1-g(W)). Returns updated outcome predictions,
#' the fluctuation parameter epsilon, and the efficient influence curve.
#'
#' @param g_fit A `tmle_mechanism` of type `"treatment"` from
#'   [fit_tmle_treatment_mechanism()].
#' @param Q_fit A `tmle_mechanism` of type `"outcome"` from
#'   [fit_tmle_outcome_mechanism()].
#'
#' @return An object of class `tmle_update` containing updated predictions,
#'   the fluctuation parameter, the efficient influence curve, and the
#'   initial point estimate.
#'
#' @keywords internal
run_tmle_targeting_step <- function(g_fit, Q_fit) {
  if (!inherits(g_fit, "tmle_mechanism") || g_fit$type != "treatment")
    stop("`g_fit` must be a tmle_mechanism of type 'treatment'.",
         call. = FALSE)
  if (!inherits(Q_fit, "tmle_mechanism") || Q_fit$type != "outcome")
    stop("`Q_fit` must be a tmle_mechanism of type 'outcome'.",
         call. = FALSE)

  data      <- Q_fit$data
  treatment <- Q_fit$treatment
  outcome   <- Q_fit$outcome
  A         <- data[[treatment]]
  Y         <- data[[outcome]]
  ps        <- g_fit$ps
  n         <- nrow(data)

  Q_a1 <- Q_fit$Q_a1
  Q_a0 <- Q_fit$Q_a0
  Q_aw <- Q_fit$Q_aw

  # Clever covariate
  H_a1 <-  1 / ps
  H_a0 <- -1 / (1 - ps)
  H_aw <- ifelse(A == 1, H_a1, H_a0)

  # Fluctuation via logistic submodel
  Q_aw_logit <- stats::qlogis(pmax(pmin(Q_aw, 0.999), 0.001))

  epsilon <- tryCatch({
    fluc <- stats::glm(
      Y ~ -1 + H_aw + offset(Q_aw_logit),
      family = stats::binomial()
    )
    unname(stats::coef(fluc))
  }, error = function(e) {
    # Surface the failure so callers know the TMLE targeting step did
    # not converge and the returned estimate is the untargeted Q (i.e.
    # behaves like g-computation rather than TMLE).
    warning(
      "run_tmle_targeting_step: targeting-step GLM failed (",
      conditionMessage(e), "). Falling back to epsilon = 0 ",
      "(untargeted plug-in). Inspect the clever-covariate distribution ",
      "and propensity-score truncation before treating the result as a ",
      "TMLE estimate.",
      call. = FALSE
    )
    0
  })

  # Updated predictions
  Q_a1_upd <- stats::plogis(
    stats::qlogis(pmax(pmin(Q_a1, 0.999), 0.001)) + epsilon * H_a1
  )
  Q_a0_upd <- stats::plogis(
    stats::qlogis(pmax(pmin(Q_a0, 0.999), 0.001)) + epsilon * H_a0
  )
  Q_aw_upd <- stats::plogis(
    stats::qlogis(pmax(pmin(Q_aw, 0.999), 0.001)) + epsilon * H_aw
  )

  psi <- mean(Q_a1_upd) - mean(Q_a0_upd)

  # Efficient influence curve
  eic <- H_aw * (Y - Q_aw_upd) + (Q_a1_upd - Q_a0_upd) - psi

  result <- list(
    psi      = psi,
    epsilon  = epsilon,
    Q_a1_upd = Q_a1_upd,
    Q_a0_upd = Q_a0_upd,
    Q_aw_upd = Q_aw_upd,
    eic      = eic,
    g_fit    = g_fit,
    Q_fit    = Q_fit,
    treatment = treatment,
    outcome   = outcome,
    n         = n,
    call      = match.call()
  )
  class(result) <- "tmle_update"
  result
}


#' Extract Final TMLE Estimate
#'
#' Extracts the final TMLE point estimate, influence-curve standard error,
#' 95% confidence interval, p-value, and one-step diagnostics from a
#' `tmle_update` object.
#'
#' @param tmle_upd A `tmle_update` object from [run_tmle_targeting_step()].
#'
#' @return An object of class `tmle_fit` (inherits from `cr_result`). Its
#'   `estimates$ATE` holds the estimate, SE, CI, and p-value, and `risks` holds
#'   the model-based treatment-specific risks `risks$treated` and
#'   `risks$control` (the TMLE plug-in arm risks; their difference equals the
#'   ATE). Report arm risks from these fields rather than reconstructing them
#'   from crude means and the ATE.
#'
#' @keywords internal
extract_tmle_estimate <- function(tmle_upd) {
  if (!inherits(tmle_upd, "tmle_update"))
    stop("`tmle_upd` must be a tmle_update object from ",
         "run_tmle_targeting_step().", call. = FALSE)

  psi <- tmle_upd$psi
  eic <- tmle_upd$eic
  n   <- tmle_upd$n

  # Model-based (TMLE plug-in) treatment-specific risks. These are the
  # confounding-adjusted arm risks the estimator actually targeted; the ATE
  # is exactly their difference (psi = mean(Q_a1_upd) - mean(Q_a0_upd)), so
  # reporting code should use these rather than reconstructing arm risks from
  # crude means and the ATE.
  risk_treated <- mean(tmle_upd$Q_a1_upd)
  risk_control <- mean(tmle_upd$Q_a0_upd)

  # eic may contain NAs at rows with missing outcome (complete-case Q fit
  # plus Y NA in the (Y - Q_aw_upd) term). Use the n_eff = sum(!is.na(eic))
  # for the IF-based SE; this is the complete-case TMLE variance estimator.
  eic_obs <- eic[!is.na(eic)]
  n_eff   <- length(eic_obs)
  if (n_eff < 2L) {
    se      <- NA_real_
  } else {
    se      <- sqrt(stats::var(eic_obs) / n_eff)
  }
  ci_lo   <- psi - 1.96 * se
  ci_hi   <- psi + 1.96 * se
  p_value <- 2 * stats::pnorm(-abs(psi / se))

  result <- list(
    estimates = list(
      ATE = list(
        estimate = psi,
        se       = se,
        ci_lower = ci_lo,
        ci_upper = ci_hi,
        p_value  = p_value
      )
    ),
    risks = list(
      treated = risk_treated,
      control = risk_control
    ),
    tmle_obj        = tmle_upd,
    influence_curve = eic,
    diagnostics = list(
      mean_eic = mean(eic),
      epsilon  = tmle_upd$epsilon,
      n        = n
    ),
    treatment = tmle_upd$treatment,
    outcome   = tmle_upd$outcome,
    type      = "modular_tmle",
    call      = match.call()
  )
  class(result) <- c("tmle_fit", "cr_result")
  result
}


# ── Summary ───────────────────────────────────────────────────────────────

#' Summarize Cleanroom Results
#'
#' Creates a side-by-side comparison table of estimates from multiple
#' cleanroom workflow objects (matching, IPTW, TMLE).
#'
#' @param fits A named or unnamed list of fitted workflow objects.
#'   Supported classes: `match_result`, `iptw_result`, `tmle_fit`.
#' @param ... Currently unused.
#'
#' @return A data.frame with one row per workflow containing columns
#'   `method`, `estimate`, `se`, `ci_lower`, `ci_upper`, and `p_value`.
#'
#' @keywords internal
summarize_cleanroom_results <- function(fits, ...) {
  if (!is.list(fits))
    stop("`fits` must be a list of fitted workflow objects.", call. = FALSE)

  rows <- lapply(seq_along(fits), function(i) {
    fit <- fits[[i]]
    nm  <- if (!is.null(names(fits)) && nzchar(names(fits)[i]))
      names(fits)[i] else ""

    if (inherits(fit, "match_result")) {
      method <- if (nzchar(nm)) nm else "PS Matching"
      data.frame(method   = method,
                 estimate = fit$estimate,
                 se       = fit$se,
                 ci_lower = fit$ci_lower,
                 ci_upper = fit$ci_upper,
                 p_value  = fit$p_value,
                 stringsAsFactors = FALSE)

    } else if (inherits(fit, "iptw_result")) {
      method <- if (nzchar(nm)) nm else "IPTW"
      data.frame(method   = method,
                 estimate = fit$estimate,
                 se       = fit$se,
                 ci_lower = fit$ci_lower,
                 ci_upper = fit$ci_upper,
                 p_value  = fit$p_value,
                 stringsAsFactors = FALSE)

    } else if (inherits(fit, "tmle_fit")) {
      method <- if (nzchar(nm)) nm else "TMLE"
      ate    <- fit$estimates$ATE
      data.frame(method   = method,
                 estimate = ate$estimate,
                 se       = ate$se,
                 ci_lower = ate$ci_lower,
                 ci_upper = ate$ci_upper,
                 p_value  = ate$p_value,
                 stringsAsFactors = FALSE)

    } else {
      warning("Unsupported fit class: ",
              paste(class(fit), collapse = ", "), ". Skipping.",
              call. = FALSE)
      NULL
    }
  })

  rows <- rows[!vapply(rows, is.null, logical(1L))]
  if (length(rows) == 0L)
    stop("No supported workflow objects found in `fits`.", call. = FALSE)

  result <- do.call(rbind, rows)
  rownames(result) <- NULL
  result
}


#' Summarize Plasmode Results
#'
#' Convenience function that prints and invisibly returns the performance
#' metrics from [run_plasmode_feasibility()].
#'
#' @param x A `plasmode_results` object from [run_plasmode_feasibility()].
#' @param ... Currently unused.
#'
#' @return Invisibly returns `x`.
#'
#' @keywords internal
summarize_plasmode_results <- function(x, ...) {
  if (!inherits(x, "plasmode_results"))
    stop("`x` must be a plasmode_results object.", call. = FALSE)
  print(x)
  invisible(x)
}


#' Fit All Stage 3 Workflows
#'
#' Convenience wrapper that runs matching, IPTW, and modular TMLE in a
#' single call, returning a named list of all fitted workflow objects.
#'
#' @param lock A `cleanroom_lock` from [create_analysis_lock()].
#' @param ps_fit A `ps_fit` object from [fit_ps_superlearner()] or
#'   [fit_ps_glm()].
#' @param workflows Character vector specifying which workflows to run.
#'   Default: `c("match", "iptw", "tmle")`.  Matching and IPTW serve as
#'   secondary comparators; the primary TMLE uses the locked specification
#'   from [lock_primary_tmle_spec()] when available.
#' @param allow_outcome_access Logical; if \code{TRUE}, skips the outcome-access check. Default \code{FALSE}.
#' @param override_clean_room Deprecated. Use \code{allow_outcome_access}.
#'
#' @return A named list with elements named by the requested workflows.
#'
#' @keywords internal
fit_final_workflows <- function(lock, ps_fit,
                                 workflows = c("match", "iptw", "tmle"),
                                 allow_outcome_access = FALSE,
                                 override_clean_room = NULL) {
  if (!is.null(override_clean_room)) {
    rlang::warn(
      "override_clean_room is deprecated in fit_final_workflows(); use allow_outcome_access.",
      .frequency = "once", .frequency_id = "fit_final_workflows_deprecated"
    )
    allow_outcome_access <- override_clean_room
  }
  .check_outcome_access(lock, allow_outcome_access,
                        caller = "fit_final_workflows")
  workflows <- match.arg(workflows, choices = c("match", "iptw", "tmle"),
                          several.ok = TRUE)
  results <- list()

  if ("match" %in% workflows) {
    results$match <- tryCatch(
      run_match_workflow(lock, ps_fit),
      error = function(e) {
        message("Matching workflow failed: ", e$message)
        NULL
      }
    )
  }

  if ("iptw" %in% workflows) {
    results$iptw <- tryCatch(
      run_iptw_workflow(lock, ps_fit),
      error = function(e) {
        message("IPTW workflow failed: ", e$message)
        NULL
      }
    )
  }

  if ("tmle" %in% workflows) {
    results$tmle <- tryCatch({
      # fit_tmle_treatment_mechanism and fit_tmle_outcome_mechanism
      # automatically resolve truncation and Q-library from the locked
      # primary TMLE spec when present.
      g_fit    <- fit_tmle_treatment_mechanism(lock, ps_fit)
      Q_fit    <- fit_tmle_outcome_mechanism(lock, g_fit)
      tmle_upd <- run_tmle_targeting_step(g_fit, Q_fit)
      extract_tmle_estimate(tmle_upd)
    }, error = function(e) {
      message("TMLE workflow failed: ", e$message)
      NULL
    })
  }

  results[!vapply(results, is.null, logical(1L))]
}


#' Fit a Set of TMLE Candidate Specifications on Real Data
#'
#' Fits multiple TMLE candidate specifications on the real outcome data
#' (Stage 4) and returns all results.
#'
#' @param lock A `cleanroom_lock` from [create_analysis_lock()].
#' @param candidates A list of [tmle_candidate()] objects, or a named list
#'   of SuperLearner library vectors (legacy API). If `NULL`, uses
#'   [expand_tmle_candidate_grid()] defaults.
#' @param ps_fit Optional `ps_fit` object; reused across candidates to
#'   avoid redundant PS estimation.
#' @param allow_outcome_access Logical; if \code{TRUE}, skips the outcome-access check. Default \code{FALSE}.
#' @param override_clean_room Deprecated. Use \code{allow_outcome_access}.
#'
#' @return A named list of TMLE estimate objects (failed candidates dropped).
#'
#' @section Clean-room stage: Stage 4 (accesses the real outcome).
#'
#' @keywords internal
fit_tmle_candidate_set <- function(lock, candidates = NULL, ps_fit = NULL,
                                    allow_outcome_access = FALSE,
                                    override_clean_room = NULL) {
  if (!is.null(override_clean_room)) {
    rlang::warn(
      "override_clean_room is deprecated in fit_tmle_candidate_set(); use allow_outcome_access.",
      .frequency = "once", .frequency_id = "fit_tmle_candidate_set_deprecated"
    )
    allow_outcome_access <- override_clean_room
  }
  if (!inherits(lock, "cleanroom_lock"))
    stop("`lock` must be a cleanroom_lock object.", call. = FALSE)
  .check_outcome_access(lock, allow_outcome_access,
                        caller = "fit_tmle_candidate_set")

  # Default: use grid expansion

  if (is.null(candidates)) {
    candidates <- expand_tmle_candidate_grid()
  }

  # Accept either tmle_candidate_spec list or legacy named-library list
  use_specs <- all(vapply(candidates, function(x)
    inherits(x, "tmle_candidate_spec"), logical(1)))

  if (use_specs) {
    cand_ids <- vapply(candidates, function(x) x$candidate_id, character(1))
    results <- lapply(seq_along(candidates), function(i) {
      cand <- candidates[[i]]
      lock_i <- lock
      lock_i$primary_tmle_spec <- cand
      class(lock_i) <- "cleanroom_lock"

      tryCatch({
        ps_i <- if (!is.null(ps_fit)) ps_fit else fit_ps_glm(lock_i)
        g_i     <- fit_tmle_treatment_mechanism(lock_i, ps_i)
        Q_i     <- fit_tmle_outcome_mechanism(lock_i, g_i)
        upd_i   <- run_tmle_targeting_step(g_i, Q_i)
        extract_tmle_estimate(upd_i)
      }, error = function(e) {
        message("Candidate '", cand$candidate_id, "' failed: ", e$message)
        NULL
      })
    })
    names(results) <- cand_ids
  } else {
    # Legacy: named list of SL library vectors
    results <- lapply(names(candidates), function(nm) {
      lib    <- candidates[[nm]]
      lock_i <- lock
      lock_i$sl_library <- lib
      class(lock_i) <- "cleanroom_lock"

      tryCatch({
        ps_i <- if (!is.null(ps_fit)) ps_fit else fit_ps_glm(lock_i)
        g_i     <- fit_tmle_treatment_mechanism(lock_i, ps_i)
        Q_i     <- fit_tmle_outcome_mechanism(lock_i, g_i)
        upd_i   <- run_tmle_targeting_step(g_i, Q_i)
        extract_tmle_estimate(upd_i)
      }, error = function(e) {
        message("Candidate '", nm, "' failed: ", e$message)
        NULL
      })
    })
    names(results) <- names(candidates)
  }

  results[!vapply(results, is.null, logical(1L))]
}


# ── Crude Workflow ────────────────────────────────────────────────────────

#' Run Crude (Unadjusted) Risk Difference Workflow
#'
#' Computes the unadjusted risk difference between treatment groups
#' without any propensity-score or covariate adjustment. Intended as
#' a benchmark comparator for the adjusted workflows.
#'
#' @param lock A `cleanroom_lock` from [create_analysis_lock()].
#' @param allow_outcome_access Logical; if \code{TRUE}, skips the outcome-access check. Default \code{FALSE}.
#' @param override_clean_room Deprecated. Use \code{allow_outcome_access}.
#'
#' @return A list with elements `estimate`, `se`, `ci_lower`, `ci_upper`,
#'   `p_value`, `r1` (risk in treated), and `r0` (risk in control).
#'
#' @examples
#' \dontrun{
#' dat  <- sim_func1(n = 500, seed = 1)
#' lock <- create_analysis_lock(
#'   data = dat, treatment = "treatment", outcome = "event_24",
#'   covariates = c("age", "sex", "biomarker"), seed = 1
#' )
#' run_crude_workflow(lock)
#'
#' }
#' @keywords internal
run_crude_workflow <- function(lock, allow_outcome_access = FALSE,
                               override_clean_room = NULL) {
  if (!is.null(override_clean_room)) {
    rlang::warn(
      "override_clean_room is deprecated in run_crude_workflow(); use allow_outcome_access.",
      .frequency = "once", .frequency_id = "run_crude_workflow_deprecated"
    )
    allow_outcome_access <- override_clean_room
  }
  if (!inherits(lock, "cleanroom_lock"))
    stop("`lock` must be a cleanroom_lock object.", call. = FALSE)
  .check_outcome_access(lock, allow_outcome_access,
                        caller = "run_crude_workflow")

  data <- .join_outcome(lock)
  A    <- data[[lock$treatment]]
  Y    <- data[[lock$outcome]]
  Y    <- as.numeric(Y)

  if (any(is.na(Y))) {
    warning("run_crude_workflow: ", sum(is.na(Y)), " of ", length(Y),
            " outcome rows are NA; computing crude RD on complete cases.",
            call. = FALSE)
    keep <- !is.na(Y)
    A <- A[keep]; Y <- Y[keep]
  }

  n1 <- sum(A == 1); n0 <- sum(A == 0)
  r1 <- mean(Y[A == 1]); r0 <- mean(Y[A == 0])
  rd <- r1 - r0

  se <- sqrt(stats::var(Y[A == 1]) / n1 + stats::var(Y[A == 0]) / n0)
  ci_lo   <- rd - 1.96 * se
  ci_hi   <- rd + 1.96 * se
  p_value <- 2 * stats::pnorm(-abs(rd / se))

  list(
    estimate = rd,
    se       = se,
    ci_lower = ci_lo,
    ci_upper = ci_hi,
    p_value  = p_value,
    r1       = r1,
    r0       = r0,
    n        = length(Y),
    treatment = lock$treatment,
    outcome   = lock$outcome
  )
}
