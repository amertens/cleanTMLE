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
#' @description Internal helper used by Stage 4 functions to enforce
#'   the clean-room gate.  Checks for outcome masking and optionally
#'   for an \code{audit} attribute carrying authorization status.
#'
#' @param lock A \code{cleanroom_lock}.
#' @param override_clean_room Logical; if \code{TRUE} the check is skipped.
#' @param caller Character; name of the calling function (for error messages).
#'
#' @return \code{invisible(TRUE)} if access is permitted.
#' @keywords internal
.check_outcome_access <- function(lock, override_clean_room = FALSE,
                                  caller = "Stage 4 function") {
  if (isTRUE(override_clean_room)) return(invisible(TRUE))

  # Check outcome masking

  if (isTRUE(lock$.outcome_masked)) {
    stop(
      caller, ": Outcome is masked. ",
      "Use unmask_outcome() first, or set override_clean_room = TRUE.",
      call. = FALSE
    )
  }

  invisible(TRUE)
}


# ── Stage 1: Analysis Lock ────────────────────────────────────────────────

#' Create an Analysis Lock
#'
#' Captures the full analytic specification—data, variable names, SuperLearner
#' library, and plasmode simulation settings—as a locked object. The lock
#' is created before the outcome is examined to enforce the pre-specified
#' analytic plan across all subsequent stages.
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
                                  sl_library    = c("SL.glm", "SL.mean"),
                                  plasmode_reps = 100L,
                                  seed          = 42L) {
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

  lock_hash <- .compute_lock_hash(list(
    treatment  = treatment,
    outcome    = outcome,
    covariates = covariates,
    sl_library = sl_library,
    seed       = as.integer(seed),
    data_nrow  = nrow(data),
    data_ncol  = ncol(data),
    data_names = paste(sort(names(data)), collapse = "|")
  ))

  lock <- list(
    data          = data,
    treatment     = treatment,
    outcome       = outcome,
    covariates    = covariates,
    sl_library    = sl_library,
    plasmode_reps = as.integer(plasmode_reps),
    seed          = as.integer(seed),
    locked_at     = Sys.time(),
    lock_hash     = lock_hash
  )
  class(lock) <- "cleanroom_lock"
  lock
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
#' @export
validate_analysis_lock <- function(lock) {
  if (!inherits(lock, "cleanroom_lock"))
    stop("`lock` must be a cleanroom_lock object.", call. = FALSE)

  required <- c("data", "treatment", "outcome", "covariates",
                "sl_library", "plasmode_reps", "seed", "lock_hash")
  missing_fields <- setdiff(required, names(lock))
  if (length(missing_fields) > 0L)
    stop("Lock is missing required fields: ",
         paste(missing_fields, collapse = ", "), call. = FALSE)

  # Recompute and compare hash
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
  if (!identical(lock$lock_hash, computed_hash))
    stop("Lock hash mismatch: the lock object may have been modified.",
         call. = FALSE)

  # Variable presence
  if (!lock$treatment %in% names(lock$data))
    stop("treatment variable '", lock$treatment,
         "' not found in locked data.", call. = FALSE)
  if (!lock$outcome %in% names(lock$data))
    stop("outcome variable '", lock$outcome,
         "' not found in locked data.", call. = FALSE)
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
  cat("Data:       ", nrow(x$data), "observations,", ncol(x$data),
      "variables\n")
  cat("Treatment:  ", x$treatment, "\n")
  cat("Outcome:    ", x$outcome, "\n")
  cat("Covariates: ", paste(x$covariates, collapse = ", "), "\n")
  cat("SL library: ", paste(x$sl_library, collapse = ", "), "\n")
  cat("Plasmode:   ", x$plasmode_reps, "replicates\n")
  cat("Seed:       ", x$seed, "\n")
  if (!is.null(x$locked_at))
    cat("Locked at:  ", format(x$locked_at, "%Y-%m-%d %H:%M:%S"), "\n")
  cat("Hash:       ", x$lock_hash, "\n")

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
#' @export
fit_ps_superlearner <- function(lock, ...) {
  if (!inherits(lock, "cleanroom_lock"))
    stop("`lock` must be a cleanroom_lock object.", call. = FALSE)
  if (!requireNamespace("SuperLearner", quietly = TRUE))
    stop("Package 'SuperLearner' is required. ",
         "Install with: install.packages('SuperLearner')", call. = FALSE)

  data <- lock$data
  A    <- data[[lock$treatment]]
  W    <- data[, lock$covariates, drop = FALSE]

  set.seed(lock$seed)
  sl_fit <- SuperLearner::SuperLearner(
    Y          = A,
    X          = W,
    family     = binomial(),
    SL.library = lock$sl_library,
    env        = asNamespace("SuperLearner"),
    ...
  )

  ps <- as.numeric(sl_fit$SL.predict)

  result <- list(
    ps         = ps,
    sl_fit     = sl_fit,
    treatment  = lock$treatment,
    covariates = lock$covariates,
    data       = data,
    sl_library = lock$sl_library,
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
#' dat  <- sim_func1(n = 200, seed = 1)
#' lock <- create_analysis_lock(
#'   data = dat, treatment = "treatment", outcome = "event_24",
#'   covariates = c("age", "sex", "biomarker"), seed = 1L
#' )
#' ps_fit <- fit_ps_glm(lock)
#' print(ps_fit)
#'
#' @export
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
#' @export
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
#'
#' @return A list of class \code{tmle_candidate_spec}.
#'
#' @examples
#' tmle_candidate("glm_t01", "GLM, trunc=0.01",
#'                g_library = "SL.glm", truncation = 0.01)
#'
#' @export
tmle_candidate <- function(candidate_id, label = candidate_id,
                           g_library  = c("SL.glm"),
                           q_library  = NULL,
                           truncation = 0.01) {
  if (!is.character(candidate_id) || length(candidate_id) != 1L)
    stop("`candidate_id` must be a single character string.", call. = FALSE)
  if (is.null(q_library)) q_library <- g_library

  obj <- list(
    candidate_id = candidate_id,
    label        = label,
    g_library    = g_library,
    q_library    = q_library,
    truncation   = truncation
  )
  class(obj) <- "tmle_candidate_spec"
  obj
}


#' @export
print.tmle_candidate_spec <- function(x, ...) {
  cat(sprintf("TMLE Candidate: %s\n", x$candidate_id))
  cat(sprintf("  Label:      %s\n", x$label))
  cat(sprintf("  g-library:  %s\n", paste(x$g_library, collapse = ", ")))
  cat(sprintf("  Q-library:  %s\n", paste(x$q_library, collapse = ", ")))
  cat(sprintf("  Truncation: %s\n", x$truncation))
  invisible(x)
}


#' Validate a List of TMLE Candidate Specifications
#'
#' Checks that all candidates are \code{tmle_candidate_spec} objects with
#' unique IDs and valid settings.
#'
#' @param candidates A list of \code{tmle_candidate_spec} objects.
#'
#' @return Invisibly returns \code{candidates} if valid; errors otherwise.
#'
#' @export
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
#'
#' @return A list of \code{tmle_candidate_spec} objects.
#'
#' @examples
#' grid <- expand_tmle_candidate_grid()
#' length(grid)
#' grid[[1]]
#'
#' @export
expand_tmle_candidate_grid <- function(
    truncations = c(0.01, 0.05),
    libraries   = list(
      glm      = c("SL.glm"),
      glm_mean = c("SL.glm", "SL.mean")
    )) {

  candidates <- list()
  for (lib_name in names(libraries)) {
    for (trunc in truncations) {
      trunc_label <- sub("\\.", "", sprintf("t%s", trunc))
      cand_id <- paste0("tmle_", lib_name, "_", trunc_label)
      label   <- sprintf("TMLE: %s, trunc=%.3g",
                          paste(libraries[[lib_name]], collapse = "+"),
                          trunc)
      candidates[[cand_id]] <- tmle_candidate(
        candidate_id = cand_id,
        label        = label,
        g_library    = libraries[[lib_name]],
        q_library    = libraries[[lib_name]],
        truncation   = trunc
      )
    }
  }
  candidates
}


# ── Stage 2b: Plasmode Feasibility ────────────────────────────────────────

#' Run Plasmode-Simulation Feasibility Evaluation
#'
#' Evaluates the performance of prespecified TMLE candidate specifications
#' using plasmode simulation.  Synthetic binary outcomes are generated from
#' a parametric baseline-risk model fit on the real covariates, augmented
#' with a specified additive treatment effect.  Each TMLE candidate is fit
#' on every replicate and performance metrics (bias, RMSE, coverage,
#' empirical SD, mean SE) are computed against the known true effect.
#'
#' @section Clean-room stage: Stage 2b (pre-outcome).  The real
#'   treatment--outcome association is never used; only simulated
#'   outcomes are generated.
#'
#' @param lock A `cleanroom_lock` from [create_analysis_lock()].
#' @param tmle_candidates A list of \code{\link{tmle_candidate}} objects.
#'   If \code{NULL}, a default grid is generated via
#'   \code{\link{expand_tmle_candidate_grid}}.
#' @param effect_sizes Numeric vector of true risk differences to simulate.
#'   Default: `c(0.05, 0.10)`.
#' @param reps Integer; number of plasmode replicates per effect size.
#'   Defaults to `lock$plasmode_reps`.
#' @param verbose Logical; if `TRUE`, print progress messages. Default: `FALSE`.
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
#' @export
run_plasmode_feasibility <- function(lock,
                                      tmle_candidates = NULL,
                                      effect_sizes    = c(0.05, 0.10),
                                      reps            = lock$plasmode_reps,
                                      verbose         = FALSE) {
  if (!inherits(lock, "cleanroom_lock"))
    stop("`lock` must be a cleanroom_lock object.", call. = FALSE)

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
  Y <- data[[outcome]]
  n <- nrow(data)

  # Fit baseline outcome model (covariates only; outcome-blind in the
  # sense that the treatment--outcome association is not used)
  # Use complete cases for GLM fitting, then predict for all rows
  Q0_fml <- stats::reformulate(covariates, response = outcome)
  Q0_fit <- stats::glm(Q0_fml, data = data, family = stats::binomial())
  p_base <- as.numeric(stats::predict(Q0_fit, type = "response",
                                       newdata = data))
  # Handle NAs from quasi-separation or missing covariate values
  if (any(is.na(p_base))) {
    p_base[is.na(p_base)] <- mean(p_base, na.rm = TRUE)
  }
  # Clamp to valid probability range
  p_base <- pmin(pmax(p_base, 0.001), 0.999)

  cand_ids <- vapply(tmle_candidates, function(x) x$candidate_id, character(1))

  all_results <- vector("list", length(effect_sizes))
  names(all_results) <- as.character(effect_sizes)

  for (es_idx in seq_along(effect_sizes)) {
    es <- effect_sizes[es_idx]
    if (verbose) message("Effect size: ", es)

    rep_results <- vector("list", reps)

    for (rep_i in seq_len(reps)) {
      set.seed(lock$seed + rep_i)

      # Synthetic outcomes: additive risk difference es for treated
      p1_sim <- pmin(p_base + es, 0.999)
      p0_sim <- p_base
      p_obs  <- ifelse(A == 1, p1_sim, p0_sim)
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
            W_mat <- data[, covariates, drop = FALSE]
            g_sl  <- SuperLearner::SuperLearner(
              Y = A, X = W_mat, family = binomial(),
              SL.library = g_lib,
              env = asNamespace("SuperLearner")
            )
            ps_hat <- as.numeric(g_sl$SL.predict)
          } else {
            ps_fml <- stats::reformulate(covariates, response = treatment)
            ps_mod <- stats::glm(ps_fml, data = data, family = stats::binomial())
            ps_hat <- as.numeric(stats::predict(ps_mod, type = "response"))
          }
          ps_hat <- pmax(pmin(ps_hat, 1 - cand$truncation), cand$truncation)

          # Fit Q-model using the candidate's q-library
          ds <- data
          ds[[".Y_sim."]] <- Y_sim
          AW <- ds[, c(treatment, covariates), drop = FALSE]

          use_sl_q <- requireNamespace("SuperLearner", quietly = TRUE) &&
            !identical(q_lib, "SL.glm")

          if (use_sl_q) {
            q_sl <- SuperLearner::SuperLearner(
              Y = Y_sim, X = AW, family = binomial(),
              SL.library = q_lib,
              env = asNamespace("SuperLearner")
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
          H_aw <- ifelse(A == 1, H_a1, H_a0)

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

  result <- list(
    results         = all_results,
    metrics         = metrics,
    tmle_candidates = tmle_candidates,
    lock            = lock,
    effect_sizes    = effect_sizes,
    reps            = reps,
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
#' @export
select_tmle_candidate <- function(sim_results,
                                   rule = c("min_rmse", "min_bias",
                                            "max_coverage", "min_max_rmse"),
                                   thresholds = NULL,
                                   dq_results = NULL) {
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
      lock_hash        = if (!is.null(sim_results$lock))
                           sim_results$lock$lock_hash else NA_character_
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
#' @param override_clean_room Logical; if \code{TRUE}, skips the outcome-access check.  Default \code{FALSE}.
#'
#' @return An object of class `match_result` containing the causal risk
#'   difference estimate, SE, 95% CI, p-value, and the matched dataset.
#'
#' @export
run_match_workflow <- function(lock, ps_fit, caliper = NULL,
                               override_clean_room = FALSE) {
  if (!inherits(lock, "cleanroom_lock"))
    stop("`lock` must be a cleanroom_lock object.", call. = FALSE)
  if (!inherits(ps_fit, "ps_fit"))
    stop("`ps_fit` must be a ps_fit object.", call. = FALSE)
  .check_outcome_access(lock, override_clean_room,
                        caller = "run_match_workflow")

  data      <- lock$data
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

  r1   <- mean(Y_m[A_m == 1])
  r0   <- mean(Y_m[A_m == 0])
  rd   <- r1 - r0

  # Paired SE
  Y1   <- Y_m[A_m == 1]
  Y0   <- Y_m[A_m == 0][seq_len(n_m)]
  se   <- sqrt(var(Y1 - Y0) / n_m)
  ci_lo <- rd - 1.96 * se
  ci_hi <- rd + 1.96 * se
  p_val <- 2 * pnorm(-abs(rd / se))

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
#' @param override_clean_room Logical; if \code{TRUE}, skips the outcome-access check.  Default \code{FALSE}.
#'
#' @return An object of class `iptw_result` containing the estimated risk
#'   difference, SE, 95% CI, p-value, and IPTW weights.
#'
#' @export
run_iptw_workflow <- function(lock, ps_fit, trim = NULL,
                              override_clean_room = FALSE) {
  if (!inherits(lock, "cleanroom_lock"))
    stop("`lock` must be a cleanroom_lock object.", call. = FALSE)
  if (!inherits(ps_fit, "ps_fit"))
    stop("`ps_fit` must be a ps_fit object.", call. = FALSE)
  .check_outcome_access(lock, override_clean_room,
                        caller = "run_iptw_workflow")

  data      <- lock$data
  treatment <- lock$treatment
  outcome   <- lock$outcome
  A         <- data[[treatment]]
  Y         <- data[[outcome]]
  ps        <- ps_fit$ps
  n         <- nrow(data)

  # Stabilized IPTW weights
  p_trt <- mean(A)
  w     <- ifelse(A == 1, p_trt / ps, (1 - p_trt) / (1 - ps))

  if (!is.null(trim)) {
    lo <- quantile(w, trim)
    hi <- quantile(w, 1 - trim)
    w  <- pmin(pmax(w, lo), hi)
  }

  # Hajek estimator
  r1 <- weighted.mean(Y[A == 1], w[A == 1])
  r0 <- weighted.mean(Y[A == 0], w[A == 0])
  rd <- r1 - r0

  # Linearized variance for Hajek estimator
  se_sq_1 <- sum(w[A == 1]^2 * (Y[A == 1] - r1)^2) / sum(w[A == 1])^2
  se_sq_0 <- sum(w[A == 0]^2 * (Y[A == 0] - r0)^2) / sum(w[A == 0])^2
  se      <- sqrt(se_sq_1 + se_sq_0)

  ci_lo <- rd - 1.96 * se
  ci_hi <- rd + 1.96 * se
  p_val <- 2 * pnorm(-abs(rd / se))

  result <- list(
    estimate  = rd,
    se        = se,
    ci_lower  = ci_lo,
    ci_upper  = ci_hi,
    p_value   = p_val,
    r1        = r1,
    r0        = r0,
    weights   = w,
    ps        = ps,
    n         = n,
    treatment = treatment,
    outcome   = outcome,
    call      = match.call()
  )
  class(result) <- c("iptw_result", "cr_result")
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
#' @export
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

      set.seed(lock$seed + k)
      g_sl <- SuperLearner::SuperLearner(
        Y = A[train_idx], X = W[train_idx, , drop = FALSE],
        family = binomial(), SL.library = lock$sl_library,
        env = asNamespace("SuperLearner")
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
#' @param override_clean_room Logical; if \code{TRUE}, skips the outcome-
#'   access check.  Default \code{FALSE}.
#'
#' @return An object of class `tmle_mechanism` with `type = "outcome"`
#'   containing initial outcome predictions `Q_a1`, `Q_a0`, and `Q_aw`.
#'
#' @export
fit_tmle_outcome_mechanism <- function(lock, g_fit, sl_library = NULL,
                                       override_clean_room = FALSE) {
  if (!inherits(lock, "cleanroom_lock"))
    stop("`lock` must be a cleanroom_lock object.", call. = FALSE)
  if (!inherits(g_fit, "tmle_mechanism") || g_fit$type != "treatment")
    stop("`g_fit` must be a tmle_mechanism of type 'treatment'.",
         call. = FALSE)
  .check_outcome_access(lock, override_clean_room,
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

  data       <- lock$data
  treatment  <- lock$treatment
  outcome    <- lock$outcome
  covariates <- lock$covariates
  AW         <- data[, c(treatment, covariates), drop = FALSE]
  Y          <- data[[outcome]]
  n          <- nrow(data)

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
      train_idx <- which(fold_vec != k)

      set.seed(lock$seed + 1L + k)
      Q_sl_k <- SuperLearner::SuperLearner(
        Y          = Y[train_idx],
        X          = AW[train_idx, , drop = FALSE],
        family     = binomial(),
        SL.library = sl_library,
        env        = asNamespace("SuperLearner")
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
    set.seed(lock$seed + 1L)
    Q_sl <- SuperLearner::SuperLearner(
      Y          = Y,
      X          = AW,
      family     = binomial(),
      SL.library = sl_library,
      env        = asNamespace("SuperLearner")
    )
    AW_a1   <- AW; AW_a1[[treatment]] <- 1L
    AW_a0   <- AW; AW_a0[[treatment]] <- 0L
    Q_a1    <- as.numeric(predict(Q_sl, newdata = AW_a1)$pred)
    Q_a0    <- as.numeric(predict(Q_sl, newdata = AW_a0)$pred)
    Q_aw    <- as.numeric(Q_sl$SL.predict)
    Q_fit_o <- Q_sl
  } else {
    Q_fml   <- stats::reformulate(c(treatment, covariates), response = outcome)
    Q_glm   <- stats::glm(Q_fml, data = data, family = stats::binomial())
    da1     <- data; da1[[treatment]] <- 1L
    da0     <- data; da0[[treatment]] <- 0L
    Q_a1    <- as.numeric(stats::predict(Q_glm, newdata = da1, type = "response"))
    Q_a0    <- as.numeric(stats::predict(Q_glm, newdata = da0, type = "response"))
    Q_aw    <- as.numeric(stats::predict(Q_glm, type = "response"))
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
#' @export
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
  }, error = function(e) 0)

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
#' @return An object of class `tmle_fit` (inherits from `cr_result`).
#'
#' @export
extract_tmle_estimate <- function(tmle_upd) {
  if (!inherits(tmle_upd, "tmle_update"))
    stop("`tmle_upd` must be a tmle_update object from ",
         "run_tmle_targeting_step().", call. = FALSE)

  psi <- tmle_upd$psi
  eic <- tmle_upd$eic
  n   <- tmle_upd$n

  se      <- sqrt(var(eic) / n)
  ci_lo   <- psi - 1.96 * se
  ci_hi   <- psi + 1.96 * se
  p_value <- 2 * pnorm(-abs(psi / se))

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
#' @export
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
#' @export
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
#' @param override_clean_room Logical; if \code{TRUE}, skips the outcome-
#'   access check.  Default \code{FALSE}.
#'
#' @return A named list with elements named by the requested workflows.
#'
#' @export
fit_final_workflows <- function(lock, ps_fit,
                                 workflows = c("match", "iptw", "tmle"),
                                 override_clean_room = FALSE) {
  .check_outcome_access(lock, override_clean_room,
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
#' @param override_clean_room Logical; if \code{TRUE}, skips the outcome-access check.  Default \code{FALSE}.
#'
#' @return A named list of TMLE estimate objects (failed candidates dropped).
#'
#' @section Clean-room stage: Stage 4 (accesses the real outcome).
#'
#' @export
fit_tmle_candidate_set <- function(lock, candidates = NULL, ps_fit = NULL,
                                    override_clean_room = FALSE) {
  if (!inherits(lock, "cleanroom_lock"))
    stop("`lock` must be a cleanroom_lock object.", call. = FALSE)
  .check_outcome_access(lock, override_clean_room,
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


# ── Gate Decision ─────────────────────────────────────────────────────────

#' Evaluate a Plasmode-Simulation Gate (GO / FLAG / STOP)
#'
#' Compares performance metrics from a plasmode-simulation replicate set
#' against pre-specified targets and returns a structured GO / FLAG / STOP
#' decision. This implements the checkpoint logic described in the
#' clean-room staged workflow: if bias, coverage, and SE-calibration all
#' pass, the decision is **GO**; if bias and coverage pass but
#' SE-calibration fails, the decision is **FLAG** (proceed with caution);
#' otherwise the decision is **STOP**.
#'
#' @param metrics A data.frame with at least columns `bias`, `coverage`,
#'   `emp_sd`, and `mean_se`, plus a `candidate` column (from the updated
#'   plasmode API) or a `method` column (legacy API).
#' @param scenario_name Character label for the scenario (used in output).
#' @param targets A list with elements `max_abs_bias` (maximum tolerable
#'   absolute bias), `min_coverage` (minimum acceptable 95 percent CI
#'   coverage), `se_sd_low` (lower bound for SE / empirical-SD ratio), and
#'   `se_sd_high` (upper bound for SE / empirical-SD ratio).
#' @param method Character; which row to evaluate.
#'   Matched against the `candidate` column if it exists, otherwise the
#'   `method` column. Default: `"TMLE"`.
#'
#' @return A list with elements `decision` (one of `"GO"`, `"FLAG"`, or
#'   `"STOP"`), `table` (a one-row data.frame summarising the evaluation),
#'   and `scenario` (the `scenario_name` used).
#'
#' @examples
#' metrics <- data.frame(
#'   candidate = c("glm_t01", "glm_t05"),
#'   bias      = c(0.002, 0.005),
#'   coverage  = c(0.95,  0.93),
#'   emp_sd    = c(0.03,  0.04),
#'   mean_se   = c(0.031, 0.041),
#'   stringsAsFactors = FALSE
#' )
#' targets <- list(
#'   max_abs_bias = 0.01,
#'   min_coverage = 0.90,
#'   se_sd_low    = 0.8,
#'   se_sd_high   = 1.2
#' )
#' gate_check(metrics, "Example", targets, method = "glm_t01")
#'
#' @export
gate_check <- function(metrics, scenario_name, targets,
                       method = "TMLE") {
  if (!is.data.frame(metrics))
    stop("`metrics` must be a data.frame.", call. = FALSE)

  # Support both new "candidate" column and legacy "method" column
  id_col <- if ("candidate" %in% names(metrics)) "candidate" else "method"

  if (!method %in% metrics[[id_col]])
    stop("'", method, "' not found in metrics$", id_col, ".", call. = FALSE)

  row <- metrics[metrics[[id_col]] == method, , drop = FALSE]
  if (nrow(row) == 0L)
    stop("No rows matched for method '", method, "'.", call. = FALSE)
  row <- row[1L, , drop = FALSE]

  bias_ok <- abs(row$bias) < targets$max_abs_bias
  cov_ok  <- row$coverage >= targets$min_coverage

  se_sd_ratio <- if ("mean_se" %in% names(row) && "emp_sd" %in% names(row) &&
                      !is.na(row$emp_sd) && row$emp_sd > 0) {
    row$mean_se / row$emp_sd
  } else {
    NA_real_
  }
  se_cal <- if (!is.na(se_sd_ratio)) {
    se_sd_ratio >= targets$se_sd_low && se_sd_ratio <= targets$se_sd_high
  } else {
    FALSE
  }

  decision <- if (bias_ok && cov_ok && se_cal) {
    "GO"
  } else if (bias_ok && cov_ok) {
    "FLAG"
  } else {
    "STOP"
  }

  tbl <- data.frame(
    scenario    = scenario_name,
    method      = method,
    bias        = round(row$bias, 5),
    coverage    = round(row$coverage, 3),
    se_sd_ratio = if (!is.na(se_sd_ratio)) round(se_sd_ratio, 3) else NA_real_,
    bias_ok     = bias_ok,
    cov_ok      = cov_ok,
    se_cal      = se_cal,
    decision    = decision,
    stringsAsFactors = FALSE
  )

  list(decision = decision, table = tbl, scenario = scenario_name)
}


# ── Crude Workflow ────────────────────────────────────────────────────────

#' Run Crude (Unadjusted) Risk Difference Workflow
#'
#' Computes the unadjusted risk difference between treatment groups
#' without any propensity-score or covariate adjustment. Intended as
#' a benchmark comparator for the adjusted workflows.
#'
#' @param lock A `cleanroom_lock` from [create_analysis_lock()].
#' @param override_clean_room Logical; if \code{TRUE}, skips the outcome-access check.  Default \code{FALSE}.
#'
#' @return A list with elements `estimate`, `se`, `ci_lower`, `ci_upper`,
#'   `p_value`, `r1` (risk in treated), and `r0` (risk in control).
#'
#' @examples
#' dat  <- sim_func1(n = 500, seed = 1)
#' lock <- create_analysis_lock(
#'   data = dat, treatment = "treatment", outcome = "event_24",
#'   covariates = c("age", "sex", "biomarker"), seed = 1
#' )
#' run_crude_workflow(lock)
#'
#' @export
run_crude_workflow <- function(lock, override_clean_room = FALSE) {
  if (!inherits(lock, "cleanroom_lock"))
    stop("`lock` must be a cleanroom_lock object.", call. = FALSE)
  .check_outcome_access(lock, override_clean_room,
                        caller = "run_crude_workflow")

  data <- lock$data
  A    <- data[[lock$treatment]]
  Y    <- data[[lock$outcome]]

  r1 <- mean(Y[A == 1])
  r0 <- mean(Y[A == 0])
  rd <- r1 - r0

  se <- sqrt(var(Y[A == 1]) / sum(A == 1) +
               var(Y[A == 0]) / sum(A == 0))
  ci_lo   <- rd - 1.96 * se
  ci_hi   <- rd + 1.96 * se
  p_value <- 2 * pnorm(-abs(rd / se))

  list(
    estimate = rd,
    se       = se,
    ci_lower = ci_lo,
    ci_upper = ci_hi,
    p_value  = p_value,
    r1       = r1,
    r0       = r0,
    n        = nrow(data),
    treatment = lock$treatment,
    outcome   = lock$outcome
  )
}
