#' Integration Utilities for Real-World Pipelines
#'
#' Functions that address common integration friction points when
#' deploying cleanTMLE in file-based, multi-stage analysis pipelines.
#' Includes external PS wrapping, parallel SuperLearner, lock
#' serialisation, covariate sanitisation, YAML config support, and
#' matched-cohort TMLE.
#'
#' @name integration
NULL


# ── P0: Accept External PS Scores ────────────────────────────────────────

#' Wrap External Propensity Scores into a ps_fit Object
#'
#' Creates a \code{ps_fit} object from user-supplied propensity scores,
#' allowing cleanTMLE diagnostics, checkpoints, and downstream TMLE to
#' use PS computed outside the package (e.g., via parallel SuperLearner,
#' Bayesian methods, or external software).
#'
#' @param lock A \code{cleanroom_lock} from \code{\link{create_analysis_lock}}.
#' @param ps_scores Numeric vector of propensity scores, one per
#'   observation in \code{lock$data}.  Must be in (0, 1).
#' @param method Character; label for the estimation method.
#'   Default: \code{"external"}.
#' @param truncate Numeric in (0, 0.5) or \code{NULL}; if provided,
#'   scores are bounded to \code{[truncate, 1 - truncate]}.
#'   Default: \code{NULL} (no truncation; scores used as-is).
#'
#' @return An object of class \code{ps_fit} compatible with all
#'   downstream cleanTMLE functions (\code{\link{compute_ps_diagnostics}},
#'   \code{\link{checkpoint_balance}}, \code{\link{run_iptw_workflow}},
#'   etc.).
#'
#' @examples
#' \dontrun{
#' # Estimate PS outside cleanTMLE (e.g., parallel SL)
#' my_ps <- predict(my_external_model, newdata = lock$data)
#' ps_fit <- wrap_ps_fit(lock, ps_scores = my_ps, method = "parallel_sl")
#' diag   <- compute_ps_diagnostics(ps_fit)
#' }
#'
#' @export
wrap_ps_fit <- function(lock, ps_scores, method = "external",
                        truncate = NULL) {
  if (!inherits(lock, "cleanroom_lock"))
    stop("`lock` must be a cleanroom_lock object.", call. = FALSE)
  if (!is.numeric(ps_scores) || length(ps_scores) != nrow(lock$data))
    stop("`ps_scores` must be a numeric vector with one value per observation.",
         call. = FALSE)
  if (any(is.na(ps_scores)))
    stop("`ps_scores` must not contain NA values.", call. = FALSE)

  ps <- ps_scores
  if (!is.null(truncate)) {
    ps <- pmax(pmin(ps, 1 - truncate), truncate)
  }

  result <- list(
    ps         = ps,
    ps_raw     = ps_scores,
    sl_fit     = NULL,
    glm_fit    = NULL,
    treatment  = lock$treatment,
    covariates = lock$covariates,
    data       = lock$data,
    sl_library = NULL,
    method     = method,
    truncate   = truncate,
    call       = match.call()
  )
  class(result) <- c("ps_fit", "cr_result")
  result
}


# ── P0: Parallel SuperLearner ────────────────────────────────────────────

#' Fit Propensity Score with Parallel SuperLearner
#'
#' Like \code{\link{fit_ps_superlearner}} but accepts a PSOCK or FORK
#' cluster for parallel computation via
#' \code{SuperLearner::snowSuperLearner()}.  Essential for reasonable
#' runtime on Windows with non-trivial learner libraries.
#'
#' @inheritParams fit_ps_superlearner
#' @param cluster A \code{parallel} cluster object (e.g., from
#'   \code{parallel::makeCluster()}).  If \code{NULL}, falls back to
#'   sequential \code{SuperLearner::SuperLearner()}.
#' @param truncate Numeric in (0, 0.5) or \code{NULL}; PS truncation.
#'
#' @return A \code{ps_fit} object.
#'
#' @examples
#' \dontrun{
#' cl <- parallel::makeCluster(2, type = "PSOCK")
#' ps_fit <- fit_ps_parallel(lock, cluster = cl)
#' parallel::stopCluster(cl)
#' }
#'
#' @export
fit_ps_parallel <- function(lock, cluster = NULL, truncate = NULL, ...) {
  if (!inherits(lock, "cleanroom_lock"))
    stop("`lock` must be a cleanroom_lock object.", call. = FALSE)
  if (!requireNamespace("SuperLearner", quietly = TRUE))
    stop("Package 'SuperLearner' is required.", call. = FALSE)

  data <- lock$data
  A    <- data[[lock$treatment]]
  W    <- data[, lock$covariates, drop = FALSE]

  set.seed(lock$seed)

  if (!is.null(cluster)) {
    sl_fit <- SuperLearner::snowSuperLearner(
      Y          = A,
      X          = W,
      family     = binomial(),
      SL.library = lock$sl_library,
      cluster    = cluster,
      env        = asNamespace("SuperLearner"),
      ...
    )
  } else {
    sl_fit <- SuperLearner::SuperLearner(
      Y          = A,
      X          = W,
      family     = binomial(),
      SL.library = lock$sl_library,
      env        = asNamespace("SuperLearner"),
      ...
    )
  }

  ps <- as.numeric(sl_fit$SL.predict)
  ps_raw <- ps
  if (!is.null(truncate)) {
    ps <- pmax(pmin(ps, 1 - truncate), truncate)
  }

  result <- list(
    ps         = ps,
    ps_raw     = ps_raw,
    sl_fit     = sl_fit,
    treatment  = lock$treatment,
    covariates = lock$covariates,
    data       = data,
    sl_library = lock$sl_library,
    method     = if (!is.null(cluster)) "parallel_sl" else "sl",
    truncate   = truncate,
    call       = match.call()
  )
  class(result) <- c("ps_fit", "cr_result")
  result
}


# ── P1: Lock Serialisation ───────────────────────────────────────────────

#' Save a cleanroom_lock to Disk
#'
#' Serialises the lock as an RDS file, including the data, hash, and
#' all metadata.  Designed for file-based multi-stage pipelines where
#' each stage loads the lock from the previous stage.
#'
#' @param lock A \code{cleanroom_lock}.
#' @param path Character; file path (must end in \code{.rds}).
#' @param validate Logical; if \code{TRUE} (default), validates the
#'   lock before saving to ensure integrity.
#'
#' @return Invisibly returns \code{path}.
#'
#' @export
save_lock <- function(lock, path, validate = TRUE) {
  if (!inherits(lock, "cleanroom_lock"))
    stop("`lock` must be a cleanroom_lock object.", call. = FALSE)
  if (validate) validate_analysis_lock(lock)

  dir_path <- dirname(path)
  if (!dir.exists(dir_path)) dir.create(dir_path, recursive = TRUE)

  saveRDS(lock, file = path)
  invisible(path)
}


#' Load a cleanroom_lock from Disk
#'
#' Reads a lock from an RDS file and re-validates its hash to confirm
#' the object was not modified outside the clean-room workflow.
#'
#' @param path Character; path to the RDS file.
#' @param validate Logical; if \code{TRUE} (default), re-validates the
#'   lock hash after loading.
#'
#' @return A \code{cleanroom_lock} object.
#'
#' @export
load_lock <- function(path, validate = TRUE) {
  if (!file.exists(path))
    stop("Lock file not found: ", path, call. = FALSE)

  lock <- readRDS(path)

  if (!inherits(lock, "cleanroom_lock"))
    stop("File does not contain a cleanroom_lock object.", call. = FALSE)

  if (validate) validate_analysis_lock(lock)

  lock
}


#' Save an Audit Log to Disk
#'
#' @param audit A \code{cleantmle_audit}.
#' @param path Character; file path (must end in \code{.rds}).
#'
#' @return Invisibly returns \code{path}.
#'
#' @export
save_audit <- function(audit, path) {
  if (!inherits(audit, "cleantmle_audit"))
    stop("`audit` must be a cleantmle_audit object.", call. = FALSE)

  dir_path <- dirname(path)
  if (!dir.exists(dir_path)) dir.create(dir_path, recursive = TRUE)

  saveRDS(audit, file = path)
  invisible(path)
}


#' Load an Audit Log from Disk
#'
#' @param path Character; path to the RDS file.
#'
#' @return A \code{cleantmle_audit} object.
#'
#' @export
load_audit <- function(path) {
  if (!file.exists(path))
    stop("Audit file not found: ", path, call. = FALSE)

  audit <- readRDS(path)

  if (!inherits(audit, "cleantmle_audit"))
    stop("File does not contain a cleantmle_audit object.", call. = FALSE)

  audit
}


# ── P1: TMLE-Based Negative Controls ────────────────────────────────────

#' Run TMLE-Based Negative Control Analysis
#'
#' Estimates the association between treatment and a negative-control
#' outcome using the modular TMLE pipeline (doubly robust), rather than
#' IPTW.  Uses the same g-model (propensity score) as the primary
#' analysis but fits a separate Q-model for the negative-control outcome.
#'
#' @param lock A \code{cleanroom_lock} with registered negative controls.
#' @param variable Character; the negative-control outcome variable name.
#' @param ps_fit A \code{ps_fit} object from Stage 2.
#' @param sl_library Optional SuperLearner library for the Q-model.
#'   Defaults to \code{lock$sl_library}.
#'
#' @return A list of class \code{cleantmle_nc_result} with \code{method = "tmle"}.
#'
#' @export
run_negative_control_tmle <- function(lock, variable, ps_fit,
                                      sl_library = NULL) {
  if (!inherits(lock, "cleanroom_lock"))
    stop("`lock` must be a cleanroom_lock object.", call. = FALSE)
  if (!inherits(ps_fit, "ps_fit"))
    stop("`ps_fit` must be a ps_fit object.", call. = FALSE)
  if (!variable %in% names(lock$data))
    stop("Variable '", variable, "' not found in lock data.", call. = FALSE)

  if (is.null(sl_library)) sl_library <- lock$sl_library

  data       <- lock$data
  A          <- data[[lock$treatment]]
  Y_nc       <- data[[variable]]
  covariates <- lock$covariates
  ps         <- ps_fit$ps
  n          <- nrow(data)

  # Resolve truncation
  trunc <- if (!is.null(lock$primary_tmle_spec))
    lock$primary_tmle_spec$truncation else 0.01
  ps_trunc <- pmax(pmin(ps, 1 - trunc), trunc)

  # Q-model for NC outcome
  AW <- data[, c(lock$treatment, covariates), drop = FALSE]

  Q_fit <- tryCatch({
    if (requireNamespace("SuperLearner", quietly = TRUE)) {
      set.seed(lock$seed + 99L)
      sl <- SuperLearner::SuperLearner(
        Y = Y_nc, X = AW, family = binomial(),
        SL.library = sl_library,
        env = asNamespace("SuperLearner")
      )
      AW_a1 <- AW; AW_a1[[lock$treatment]] <- 1L
      AW_a0 <- AW; AW_a0[[lock$treatment]] <- 0L
      list(
        Q_aw = as.numeric(sl$SL.predict),
        Q_a1 = as.numeric(predict(sl, newdata = AW_a1)$pred),
        Q_a0 = as.numeric(predict(sl, newdata = AW_a0)$pred)
      )
    } else {
      fml <- stats::reformulate(c(lock$treatment, covariates),
                                response = variable)
      glm_fit <- stats::glm(fml, data = data, family = stats::binomial())
      da1 <- data; da1[[lock$treatment]] <- 1L
      da0 <- data; da0[[lock$treatment]] <- 0L
      list(
        Q_aw = as.numeric(stats::predict(glm_fit, type = "response")),
        Q_a1 = as.numeric(stats::predict(glm_fit, newdata = da1,
                                          type = "response")),
        Q_a0 = as.numeric(stats::predict(glm_fit, newdata = da0,
                                          type = "response"))
      )
    }
  }, error = function(e) NULL)

  if (is.null(Q_fit)) {
    # Fallback to IPTW
    return(run_negative_control(lock, variable, ps_fit))
  }

  # Targeting step
  H_a1 <-  1 / ps_trunc
  H_a0 <- -1 / (1 - ps_trunc)
  H_aw <- ifelse(A == 1, H_a1, H_a0)

  Q_aw_logit <- stats::qlogis(pmax(pmin(Q_fit$Q_aw, 0.999), 0.001))
  epsilon <- tryCatch({
    fluc <- stats::glm(Y_nc ~ -1 + H_aw + offset(Q_aw_logit),
                        family = stats::binomial())
    unname(stats::coef(fluc))
  }, error = function(e) 0)

  Q_a1_u <- stats::plogis(
    stats::qlogis(pmax(pmin(Q_fit$Q_a1, 0.999), 0.001)) + epsilon * H_a1)
  Q_a0_u <- stats::plogis(
    stats::qlogis(pmax(pmin(Q_fit$Q_a0, 0.999), 0.001)) + epsilon * H_a0)
  Q_aw_u <- stats::plogis(Q_aw_logit + epsilon * H_aw)

  rd  <- mean(Q_a1_u) - mean(Q_a0_u)
  eic <- H_aw * (Y_nc - Q_aw_u) + (Q_a1_u - Q_a0_u) - rd
  se  <- sqrt(var(eic) / n)

  ci_lo   <- rd - 1.96 * se
  ci_hi   <- rd + 1.96 * se
  p_value <- 2 * stats::pnorm(-abs(rd / se))

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
    interpretation = interpretation,
    method         = "tmle"
  )
  class(result) <- "cleantmle_nc_result"
  result
}


# ── P1: Matched-Cohort TMLE ─────────────────────────────────────────────

#' Run Modular TMLE on a Matched Subset
#'
#' Runs the modular TMLE pipeline on a subset of the locked data,
#' identified by row indices (e.g., from PS matching).  Creates a
#' temporary sub-lock, fits the Q-model on the subset, and returns
#' a standard \code{tmle_fit}.
#'
#' @param lock A \code{cleanroom_lock}.
#' @param ps_fit A \code{ps_fit} from \code{\link{fit_ps_superlearner}}
#'   or \code{\link{wrap_ps_fit}}.
#' @param subset_idx Integer vector of row indices to include
#'   (e.g., matched observation indices).
#' @param sl_library Optional SL library override.
#' @param override_clean_room Logical; skip outcome-access check.
#'
#' @return An object of class \code{tmle_fit} estimated on the subset.
#'
#' @export
run_matched_tmle <- function(lock, ps_fit, subset_idx,
                             sl_library = NULL,
                             override_clean_room = FALSE) {
  if (!inherits(lock, "cleanroom_lock"))
    stop("`lock` must be a cleanroom_lock object.", call. = FALSE)
  .check_outcome_access(lock, override_clean_room,
                        caller = "run_matched_tmle")

  if (is.null(sl_library)) sl_library <- lock$sl_library

  data_sub   <- lock$data[subset_idx, , drop = FALSE]
  ps_sub     <- ps_fit$ps[subset_idx]
  A          <- data_sub[[lock$treatment]]
  Y          <- data_sub[[lock$outcome]]
  covariates <- lock$covariates
  n          <- nrow(data_sub)

  # Resolve truncation
  trunc <- if (!is.null(lock$primary_tmle_spec))
    lock$primary_tmle_spec$truncation else 0.01
  ps_trunc <- pmax(pmin(ps_sub, 1 - trunc), trunc)

  # Q-model on subset
  AW <- data_sub[, c(lock$treatment, covariates), drop = FALSE]

  Q_result <- tryCatch({
    if (requireNamespace("SuperLearner", quietly = TRUE)) {
      set.seed(lock$seed + 2L)
      sl <- SuperLearner::SuperLearner(
        Y = Y, X = AW, family = binomial(),
        SL.library = sl_library,
        env = asNamespace("SuperLearner")
      )
      AW_a1 <- AW; AW_a1[[lock$treatment]] <- 1L
      AW_a0 <- AW; AW_a0[[lock$treatment]] <- 0L
      list(
        Q_aw = as.numeric(sl$SL.predict),
        Q_a1 = as.numeric(predict(sl, newdata = AW_a1)$pred),
        Q_a0 = as.numeric(predict(sl, newdata = AW_a0)$pred)
      )
    } else {
      fml <- stats::reformulate(c(lock$treatment, covariates),
                                response = lock$outcome)
      glm_fit <- stats::glm(fml, data = data_sub, family = stats::binomial())
      da1 <- data_sub; da1[[lock$treatment]] <- 1L
      da0 <- data_sub; da0[[lock$treatment]] <- 0L
      list(
        Q_aw = as.numeric(stats::predict(glm_fit, type = "response")),
        Q_a1 = as.numeric(stats::predict(glm_fit, newdata = da1,
                                          type = "response")),
        Q_a0 = as.numeric(stats::predict(glm_fit, newdata = da0,
                                          type = "response"))
      )
    }
  }, error = function(e) {
    stop("Q-model fitting failed on matched subset: ", e$message,
         call. = FALSE)
  })

  # Targeting step
  H_a1 <-  1 / ps_trunc
  H_a0 <- -1 / (1 - ps_trunc)
  H_aw <- ifelse(A == 1, H_a1, H_a0)

  Q_aw_logit <- stats::qlogis(pmax(pmin(Q_result$Q_aw, 0.999), 0.001))
  epsilon <- tryCatch({
    fluc <- stats::glm(Y ~ -1 + H_aw + offset(Q_aw_logit),
                        family = stats::binomial())
    unname(stats::coef(fluc))
  }, error = function(e) 0)

  Q_a1_u <- stats::plogis(
    stats::qlogis(pmax(pmin(Q_result$Q_a1, 0.999), 0.001)) + epsilon * H_a1)
  Q_a0_u <- stats::plogis(
    stats::qlogis(pmax(pmin(Q_result$Q_a0, 0.999), 0.001)) + epsilon * H_a0)
  Q_aw_u <- stats::plogis(Q_aw_logit + epsilon * H_aw)

  psi <- mean(Q_a1_u) - mean(Q_a0_u)
  eic <- H_aw * (Y - Q_aw_u) + (Q_a1_u - Q_a0_u) - psi
  se  <- sqrt(var(eic) / n)

  estimates <- list(
    ATE = list(
      estimate = psi,
      se       = se,
      ci_lower = psi - 1.96 * se,
      ci_upper = psi + 1.96 * se,
      p_value  = 2 * stats::pnorm(-abs(psi / se))
    )
  )

  result <- list(
    estimates       = estimates,
    tmle_obj        = NULL,
    influence_curve = eic,
    treatment       = lock$treatment,
    outcome         = lock$outcome,
    covariates      = covariates,
    type            = "matched_tmle",
    n_matched       = n,
    subset_idx      = subset_idx,
    call            = match.call()
  )
  class(result) <- c("tmle_fit", "cr_result")
  result
}


# ── P2: Covariate Sanitisation ───────────────────────────────────────────

#' Sanitise Covariates for SuperLearner
#'
#' Cleans a covariate matrix/data.frame for safe use with SuperLearner:
#' replaces \code{NA}/\code{NaN}/\code{Inf} with column medians,
#' drops zero-variance columns, and ensures column names are
#' formula-safe.
#'
#' @param data A data.frame.
#' @param covariates Character vector of covariate column names.
#'   If \code{NULL}, all columns are sanitised.
#' @param verbose Logical; print a summary of changes.
#'
#' @return A list with elements:
#'   \describe{
#'     \item{data}{The cleaned data.frame (full, with sanitised covariates).}
#'     \item{covariates}{Updated covariate names (after renaming / dropping).}
#'     \item{dropped}{Names of dropped zero-variance columns.}
#'     \item{renamed}{Named vector of old -> new column names.}
#'     \item{imputed}{Named vector of columns with imputed values and counts.}
#'   }
#'
#' @export
sanitize_covariates <- function(data, covariates = NULL, verbose = FALSE) {
  if (is.null(covariates)) covariates <- names(data)
  covariates <- intersect(covariates, names(data))

  dropped <- character(0)
  renamed <- character(0)
  imputed <- integer(0)

  for (v in covariates) {
    col <- data[[v]]

    # Convert factors to numeric
    if (is.factor(col) || is.character(col)) {
      data[[v]] <- as.numeric(as.factor(col)) - 1L
      col <- data[[v]]
    }

    # Impute NA/NaN/Inf
    bad_idx <- is.na(col) | is.nan(col) | is.infinite(col)
    n_bad <- sum(bad_idx)
    if (n_bad > 0L) {
      med <- median(col[!bad_idx], na.rm = TRUE)
      if (is.na(med)) med <- 0
      data[[v]][bad_idx] <- med
      imputed[v] <- n_bad
    }

    # Drop zero-variance
    if (var(data[[v]], na.rm = TRUE) == 0) {
      dropped <- c(dropped, v)
    }
  }

  # Remove zero-variance columns
  covariates <- setdiff(covariates, dropped)

  # Make names formula-safe
  safe_names <- make.names(covariates, unique = TRUE)
  rename_needed <- covariates != safe_names
  if (any(rename_needed)) {
    old_names <- covariates[rename_needed]
    new_names <- safe_names[rename_needed]
    for (i in seq_along(old_names)) {
      names(data)[names(data) == old_names[i]] <- new_names[i]
    }
    renamed <- stats::setNames(new_names, old_names)
    covariates <- safe_names
  }

  if (verbose) {
    if (length(dropped) > 0L)
      message("Dropped zero-variance: ", paste(dropped, collapse = ", "))
    if (length(imputed) > 0L)
      message("Imputed: ", paste(names(imputed), collapse = ", "))
    if (length(renamed) > 0L)
      message("Renamed: ", paste(names(renamed), "->", renamed,
                                  collapse = ", "))
  }

  list(
    data       = data,
    covariates = covariates,
    dropped    = dropped,
    renamed    = renamed,
    imputed    = imputed
  )
}


# ── P2: YAML Config Support ─────────────────────────────────────────────

#' Create an Analysis Lock from a YAML Configuration File
#'
#' Reads a YAML file containing the analytic specification and creates
#' a \code{cleanroom_lock}.  The data must be loaded separately and
#' passed as an argument.
#'
#' @param config_path Character; path to the YAML configuration file.
#' @param data A data.frame containing the analysis dataset.
#' @param section Character; top-level YAML section to use.
#'   Default: uses the full file.
#'
#' @return A \code{cleanroom_lock}.
#'
#' @details
#' The YAML file should contain keys matching the arguments to
#' \code{\link{create_analysis_lock}}:
#'
#' \preformatted{
#' treatment:
#'   variable: "A"
#' binary_outcome:
#'   variable: "Y"
#' covariates:
#'   selected_columns: ["age", "sex", "biomarker"]
#' superlearner:
#'   candidate_learners: ["SL.glm", "SL.mean"]
#'   cv_folds: 2
#' simulation:
#'   n_sims: 100
#' propensity_score:
#'   truncation_lower: 0.01
#' seed: 42
#' }
#'
#' @export
create_analysis_lock_from_yaml <- function(config_path, data,
                                           section = NULL) {
  if (!requireNamespace("yaml", quietly = TRUE))
    stop("Package 'yaml' is required. Install with: install.packages('yaml')",
         call. = FALSE)

  cfg <- yaml::read_yaml(config_path)
  if (!is.null(section)) cfg <- cfg[[section]]

  treatment <- cfg$treatment$variable %||%
    cfg$treatment %||%
    stop("YAML must contain 'treatment.variable'.", call. = FALSE)

  outcome <- cfg$binary_outcome$variable %||%
    cfg$outcome$variable %||%
    cfg$outcome %||%
    stop("YAML must contain 'binary_outcome.variable' or 'outcome.variable'.",
         call. = FALSE)

  covariates <- cfg$covariates$selected_columns %||%
    cfg$covariates %||%
    stop("YAML must contain 'covariates.selected_columns'.", call. = FALSE)

  sl_library <- cfg$superlearner$candidate_learners %||%
    cfg$sl_library %||%
    c("SL.glm", "SL.mean")

  plasmode_reps <- cfg$simulation$n_sims %||%
    cfg$plasmode_reps %||%
    100L

  seed <- cfg$seed %||% 42L

  lock <- create_analysis_lock(
    data          = data,
    treatment     = treatment,
    outcome       = outcome,
    covariates    = covariates,
    sl_library    = sl_library,
    plasmode_reps = as.integer(plasmode_reps),
    seed          = as.integer(seed)
  )

  # Optionally attach negative controls from config
  nc_vars <- cfg$negative_controls$outcomes %||% NULL
  if (!is.null(nc_vars)) {
    for (nc in nc_vars) {
      if (nc %in% names(data)) {
        lock <- define_negative_control(lock, nc)
      }
    }
  }

  lock
}


# ── P2: Clever Covariate Plot ────────────────────────────────────────────

#' Clever Covariate Diagnostic Plot
#'
#' Plots the clever covariate H(A,W) = A/g(W) - (1-A)/(1-g(W))
#' by treatment group.  Large or asymmetric values indicate
#' near-positivity violations.
#'
#' @param tmle_update A \code{tmle_update} object from
#'   \code{\link{run_tmle_targeting_step}}, or any object with a
#'   \code{clever_covariate} element.
#' @param ps_fit Alternatively, a \code{ps_fit} object to compute
#'   H from PS directly.
#' @param lock A \code{cleanroom_lock} (required if using \code{ps_fit}).
#'
#' @return A \code{ggplot2} object.
#'
#' @export
clever_covariate_plot <- function(tmle_update = NULL, ps_fit = NULL,
                                  lock = NULL) {
  if (!is.null(tmle_update) && !is.null(tmle_update$clever_covariate)) {
    H_aw <- tmle_update$clever_covariate
    A    <- tmle_update$data[[tmle_update$treatment]]
  } else if (!is.null(ps_fit) && !is.null(lock)) {
    ps <- ps_fit$ps
    A  <- lock$data[[lock$treatment]]
    H_aw <- ifelse(A == 1, 1 / ps, -1 / (1 - ps))
  } else {
    stop("Provide either `tmle_update` or both `ps_fit` and `lock`.",
         call. = FALSE)
  }

  df <- data.frame(
    H     = H_aw,
    group = ifelse(A == 1, "Treated (A=1)", "Control (A=0)"),
    stringsAsFactors = FALSE
  )

  ggplot2::ggplot(df, ggplot2::aes(x = .data$H, fill = .data$group)) +
    ggplot2::geom_histogram(position = "identity", alpha = 0.5, bins = 40L) +
    ggplot2::geom_vline(xintercept = 0, linetype = "dashed", colour = "grey40") +
    ggplot2::labs(
      x     = "Clever Covariate H(A,W)",
      y     = "Count",
      title = "Clever Covariate Distribution by Treatment Group",
      fill  = NULL
    ) +
    ggplot2::theme_minimal() +
    ggplot2::theme(legend.position = "bottom")
}


# ── P2: Matched SMD Helper ───────────────────────────────────────────────

#' Compute Standardised Mean Differences on a Matched Subset
#'
#' Calculates absolute SMDs for each covariate in a matched dataset,
#' returning a named numeric vector suitable for
#' \code{\link{love_plot_threeway}}.
#'
#' @param data A data.frame (the full or matched dataset).
#' @param treatment Character; name of the binary treatment variable.
#' @param covariates Character vector of covariate column names.
#' @param subset_idx Optional integer vector of row indices defining
#'   the matched subset.  If \code{NULL}, the full dataset is used.
#'
#' @return A named numeric vector of absolute SMDs, one per covariate.
#'
#' @export
compute_matched_smds <- function(data, treatment, covariates,
                                 subset_idx = NULL) {
  if (!is.null(subset_idx)) data <- data[subset_idx, , drop = FALSE]
  A <- data[[treatment]]

  vapply(covariates, function(v) {
    x <- data[[v]]
    if (!is.numeric(x)) x <- as.numeric(as.factor(x)) - 1L
    m1 <- mean(x[A == 1], na.rm = TRUE)
    m0 <- mean(x[A == 0], na.rm = TRUE)
    s1 <- var(x[A == 1], na.rm = TRUE)
    s0 <- var(x[A == 0], na.rm = TRUE)
    pooled <- sqrt((s1 + s0) / 2)
    if (pooled > 0) abs(m1 - m0) / pooled else 0
  }, numeric(1))
}


# ── P2: Three-Way Love Plot ─────────────────────────────────────────────

#' Three-Way Love Plot
#'
#' Extends \code{\link{love_plot}} to show unweighted, IPTW-weighted,
#' and PS-matched SMDs side by side.
#'
#' @param ps_diag A \code{ps_diagnostics} object (provides unweighted
#'   and IPTW-weighted SMDs).
#' @param matched_smds A data.frame with columns \code{variable} and
#'   \code{smd_matched}, or a named numeric vector of matched SMDs
#'   (names = covariate names).
#' @param threshold Numeric; balance threshold reference line.
#'
#' @return A \code{ggplot2} object.
#'
#' @export
love_plot_threeway <- function(ps_diag, matched_smds, threshold = 0.10) {
  if (!inherits(ps_diag, "ps_diagnostics"))
    stop("`ps_diag` must be a ps_diagnostics object.", call. = FALSE)

  smds <- ps_diag$smds

  # Normalise matched_smds input

  if (is.data.frame(matched_smds)) {
    m_smd <- matched_smds$smd_matched
    names(m_smd) <- matched_smds$variable
  } else {
    m_smd <- matched_smds
  }

  # Align to ps_diag variables
  matched_vals <- m_smd[smds$variable]
  matched_vals[is.na(matched_vals)] <- 0

  plot_df <- data.frame(
    variable = rep(smds$variable, 3L),
    type     = rep(c("Unweighted", "IPTW-Weighted", "Matched"),
                   each = nrow(smds)),
    smd      = c(abs(smds$smd_unweighted),
                 abs(smds$smd_weighted),
                 abs(matched_vals)),
    stringsAsFactors = FALSE
  )
  plot_df$type <- factor(plot_df$type,
                         levels = c("Unweighted", "IPTW-Weighted", "Matched"))

  ggplot2::ggplot(
    plot_df,
    ggplot2::aes(
      x      = .data$smd,
      y      = stats::reorder(.data$variable, .data$smd),
      shape  = .data$type,
      colour = .data$type
    )
  ) +
    ggplot2::geom_point(size = 3) +
    ggplot2::geom_vline(xintercept = threshold, linetype = "dashed",
                        colour = "grey40") +
    ggplot2::labs(
      x      = "Absolute Standardised Mean Difference",
      y      = NULL,
      title  = "Covariate Balance: Unweighted vs IPTW vs Matched",
      shape  = NULL,
      colour = NULL
    ) +
    ggplot2::theme_minimal() +
    ggplot2::theme(legend.position = "bottom")
}


# ── P2: Composite Gate ──────────────────────────────────────────────────

#' Composite Go/Flag/Stop Gate
#'
#' Evaluates all checkpoint objects and returns a single GO / FLAG /
#' STOP decision.  Useful as a unified pre-outcome gate when the
#' analysis has multiple checkpoint stages.
#'
#' @param ... One or more \code{cleantmle_checkpoint} objects.
#' @param allow_flag Logical; if \code{TRUE} (default), FLAG decisions
#'   do not block the pipeline (treated as conditional GO).
#'
#' @return A \code{cleantmle_checkpoint} with the composite decision.
#'
#' @export
gate_all <- function(..., allow_flag = TRUE) {
  checkpoints <- list(...)
  # Flatten if a list was passed
  if (length(checkpoints) == 1L && is.list(checkpoints[[1]]) &&
      !inherits(checkpoints[[1]], "cleantmle_checkpoint")) {
    checkpoints <- checkpoints[[1]]
  }

  decisions <- vapply(checkpoints, function(cp) {
    if (inherits(cp, "cleantmle_checkpoint")) cp$decision
    else NA_character_
  }, character(1))

  stages <- vapply(checkpoints, function(cp) {
    if (inherits(cp, "cleantmle_checkpoint")) cp$stage
    else "unknown"
  }, character(1))

  if (any(decisions == "STOP", na.rm = TRUE)) {
    stop_stages <- stages[decisions == "STOP"]
    decision <- "STOP"
    rationale <- paste("STOP at:", paste(stop_stages, collapse = ", "))
  } else if (any(decisions == "FLAG", na.rm = TRUE)) {
    flag_stages <- stages[decisions == "FLAG"]
    if (allow_flag) {
      decision <- "FLAG"
      rationale <- paste("FLAG (proceeding with caution) at:",
                         paste(flag_stages, collapse = ", "))
    } else {
      decision <- "STOP"
      rationale <- paste("FLAG treated as STOP at:",
                         paste(flag_stages, collapse = ", "))
    }
  } else {
    decision <- "GO"
    rationale <- "All checkpoints passed."
  }

  metrics <- data.frame(
    stage    = stages,
    decision = decisions,
    stringsAsFactors = FALSE
  )

  new_checkpoint(
    stage      = "Composite Gate",
    decision   = decision,
    metrics    = metrics,
    thresholds = list(allow_flag = allow_flag),
    rationale  = rationale
  )
}


