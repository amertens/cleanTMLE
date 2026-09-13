# Design-stage tools that never touch the primary outcome: multi-contrast
# locks from a multi-level treatment, negative controls along a restriction
# ladder, and the care-process collider check.

#' One Lock per Contrast from a Multi-Level Treatment
#'
#' Builds a named list of `cleanroom_lock` objects, one per requested
#' contrast, from a factor (or character) treatment with shared covariate
#' construction, shared negative controls, and a shared seed. Each lock's
#' data are the rows belonging to that contrast's arms, with a binary
#' treatment column added; each then gets its own propensity fit, support
#' assessment, gate and estimand ladder downstream.
#'
#' @param data A data.frame.
#' @param treatment_factor Name of the multi-level treatment column.
#' @param contrasts A named list; each element is
#'   `list(treated = <levels>, control = <levels>, label = <text>)`.
#' @param outcome,covariates,sl_library,seed As in
#'   [create_analysis_lock()].
#' @param negative_controls Optional character vector of negative-control
#'   columns registered on every lock via [define_negative_control()].
#'   Registration happens here, before any downstream variance filtering,
#'   and a registered control is never dropped silently.
#' @param treatment_name Name of the derived binary column. Default
#'   `".A_contrast"`.
#' @param cleanroom_enabled As in [create_analysis_lock()]. Default TRUE.
#' @return A named list of class `contrast_locks`.
#' @examples
#' \dontrun{
#' locks <- create_contrast_locks(
#'   dat, "transport3",
#'   contrasts = list(
#'     C1 = list(treated = "Rescue.Co ambulance",
#'               control = c("Non-Rescue.Co ambulance", "Non-ambulance"),
#'               label = "Rescue.Co vs everyone else")),
#'   outcome = "death_by_6mo", covariates = covs)
#' }
#' @export
create_contrast_locks <- function(data, treatment_factor, contrasts,
                                  outcome, covariates,
                                  sl_library = c("SL.glm", "SL.mean"),
                                  seed = 42L,
                                  negative_controls = NULL,
                                  treatment_name = ".A_contrast",
                                  cleanroom_enabled = TRUE) {
  if (!treatment_factor %in% names(data))
    stop("treatment_factor '", treatment_factor, "' not found.",
         call. = FALSE)
  if (is.null(names(contrasts)) || any(!nzchar(names(contrasts))))
    stop("`contrasts` must be a named list.", call. = FALSE)
  tf <- as.character(data[[treatment_factor]])
  locks <- list()
  for (cn in names(contrasts)) {
    cc <- contrasts[[cn]]
    if (is.null(cc$treated) || is.null(cc$control))
      stop("contrast '", cn, "' needs `treated` and `control` levels.",
           call. = FALSE)
    bad <- setdiff(c(cc$treated, cc$control), unique(tf[!is.na(tf)]))
    if (length(bad))
      stop("contrast '", cn, "': level(s) not in ", treatment_factor, ": ",
           paste(bad, collapse = ", "), call. = FALSE)
    idx <- which(tf %in% c(cc$treated, cc$control))
    sub <- data[idx, , drop = FALSE]
    sub[[treatment_name]] <- as.integer(tf[idx] %in% cc$treated)
    lock <- create_analysis_lock(
      data = sub, treatment = treatment_name, outcome = outcome,
      covariates = covariates, sl_library = sl_library, seed = seed,
      cleanroom_enabled = cleanroom_enabled)
    lock$contrast <- list(name = cn, label = cc$label %||% cn,
                          treated = cc$treated, control = cc$control,
                          source_column = treatment_factor)
    if (!is.null(negative_controls))
      for (nc in negative_controls)
        lock <- define_negative_control(lock, nc)
    lock <- .log_design_decision(lock, "contrast",
      sprintf("Lock %s: %s (treated: %s; control: %s; n = %d, %d treated).",
              cn, cc$label %||% cn, paste(cc$treated, collapse = " + "),
              paste(cc$control, collapse = " + "), nrow(sub),
              sum(sub[[treatment_name]])))
    locks[[cn]] <- lock
  }
  class(locks) <- c("contrast_locks", "list")
  locks
}

#' @export
print.contrast_locks <- function(x, ...) {
  cat("Contrast locks (", length(x), ")\n", sep = "")
  for (cn in names(x)) {
    lk <- x[[cn]]
    cat(sprintf("  %-4s %s: n = %d (%d treated / %d control)\n",
                cn, lk$contrast$label, nrow(lk$data),
                sum(lk$data[[lk$treatment]]),
                sum(lk$data[[lk$treatment]] == 0)))
  }
  invisible(x)
}


#' Negative Controls Along a Restriction Ladder
#'
#' Fits every registered negative control on every nested cohort and reports
#' where a control turns from failing to null. A control that fails on the
#' full cohort and passes after a restriction is evidence for the
#' restriction: the confounding it detects is removed by design, not by
#' adjustment, and a negative-control step that runs on one cohort cannot
#' see that.
#'
#' @param lock A `cleanroom_lock` with negative controls registered via
#'   [define_negative_control()], or pass `negative_controls`.
#' @param restrictions A named list of logical vectors (or functions of the
#'   lock data returning logical vectors), ordered from least to most
#'   restricted. The full cohort is always included first as
#'   `"full cohort"`.
#' @param negative_controls Optional character vector overriding the lock's
#'   registered controls.
#' @param method `"tmle"` (default; [run_negative_control_tmle()], with a
#'   propensity model refitted on each rung) or `"unadjusted"`
#'   (two-proportion difference, the balance-check reading).
#' @param ps_method Propensity method per rung for `method = "tmle"`:
#'   `"glm"` (default; the controls are checks, not the primary fit) or
#'   `"superlearner"`.
#' @param alpha Significance level for the fail flag. Default 0.05.
#' @param min_events Cells with fewer events than this are inestimable.
#'   Default 5.
#' @param verbose Print progress. Default TRUE.
#' @return An object of class `nc_ladder`: a data.frame with one row per
#'   control per cohort (estimate, CI, p, flagged, status) plus a
#'   `turned_null` summary naming the controls that fail on an earlier rung
#'   and pass on a later one.
#' @export
run_negative_control_ladder <- function(lock, restrictions,
                                        negative_controls = NULL,
                                        method = c("tmle", "unadjusted"),
                                        ps_method = "glm",
                                        alpha = 0.05,
                                        min_events = 5L,
                                        verbose = TRUE) {
  if (!inherits(lock, "cleanroom_lock"))
    stop("`lock` must be a cleanroom_lock object.", call. = FALSE)
  method <- match.arg(method)
  ncs <- negative_controls %||% names(lock$negative_controls)
  if (!length(ncs))
    stop("No negative controls: register them with ",
         "define_negative_control() or pass `negative_controls`.",
         call. = FALSE)
  data <- lock$data
  A <- as.integer(data[[lock$treatment]])

  cohorts <- c(list(`full cohort` = rep(TRUE, nrow(data))),
               lapply(restrictions, function(r) {
                 if (is.function(r)) r <- r(data)
                 if (!is.logical(r) || length(r) != nrow(data))
                   stop("Each restriction must be (or return) a logical ",
                        "vector over the lock data rows.", call. = FALSE)
                 r & !is.na(r)
               }))

  rows <- list()
  for (ci in seq_along(cohorts)) {
    cname <- names(cohorts)[ci]
    sel <- cohorts[[ci]]
    for (nc in ncs) {
      y <- suppressWarnings(as.numeric(data[[nc]][sel]))
      a <- A[sel]
      ok <- !is.na(y) & !is.na(a)
      y <- y[ok]; a <- a[ok]
      cells <- c(sum(y[a == 1]), sum(1 - y[a == 1]),
                 sum(y[a == 0]), sum(1 - y[a == 0]))
      if (length(y) < 20L || any(cells < min_events)) {
        rows[[length(rows) + 1L]] <- data.frame(
          cohort = cname, negative_control = nc, n = length(y),
          estimate = NA_real_, ci_lower = NA_real_, ci_upper = NA_real_,
          p_value = NA_real_, flagged = NA,
          status = sprintf("inestimable (a cell below %d events)",
                           min_events),
          stringsAsFactors = FALSE)
        next
      }
      res <- if (method == "tmle") {
        sub_lock <- lock
        sub_lock$data <- data[sel, , drop = FALSE]
        f <- tryCatch({
          psf <- fit_ps(sub_lock, method = ps_method)
          run_negative_control_tmle(sub_lock, nc, psf)
        }, error = function(e) NULL)
        if (is.null(f)) NULL else
          list(est = f$estimate, lo = f$ci_lower, hi = f$ci_upper,
               p = f$p_value)
      } else NULL
      if (is.null(res)) {
        p1 <- mean(y[a == 1]); p0 <- mean(y[a == 0])
        se <- sqrt(p1 * (1 - p1) / sum(a == 1) +
                     p0 * (1 - p0) / sum(a == 0))
        est <- p1 - p0
        res <- list(est = est, lo = est - 1.96 * se, hi = est + 1.96 * se,
                    p = 2 * stats::pnorm(-abs(est / se)))
      }
      rows[[length(rows) + 1L]] <- data.frame(
        cohort = cname, negative_control = nc, n = length(y),
        estimate = round(res$est, 5), ci_lower = round(res$lo, 5),
        ci_upper = round(res$hi, 5), p_value = round(res$p, 5),
        flagged = res$p < alpha, status = "estimated",
        stringsAsFactors = FALSE)
      if (verbose)
        message(sprintf("  [%s] %s: %.4f (p = %.3f)%s", cname, nc, res$est,
                        res$p, if (res$p < alpha) " FLAGGED" else ""))
    }
  }
  tab <- do.call(rbind, rows)

  turned_null <- character(0)
  for (nc in ncs) {
    sub <- tab[tab$negative_control == nc & tab$status == "estimated", ]
    if (nrow(sub) >= 2L && isTRUE(sub$flagged[1]) &&
        any(!sub$flagged[-1], na.rm = TRUE)) {
      first_pass <- sub$cohort[-1][which(!sub$flagged[-1])[1]]
      turned_null <- c(turned_null, sprintf(
        "%s fails on the full cohort and is null after '%s': the restriction, not adjustment, removed the association.",
        nc, first_pass))
    }
  }

  out <- list(table = tab, turned_null = turned_null, method = method,
              alpha = alpha, call = match.call())
  class(out) <- "nc_ladder"
  out
}

#' @export
print.nc_ladder <- function(x, ...) {
  cat("Negative controls along the restriction ladder (", x$method,
      ")\n\n", sep = "")
  print(x$table[, c("cohort", "negative_control", "n", "estimate",
                    "ci_lower", "ci_upper", "p_value", "flagged")],
        row.names = FALSE)
  if (length(x$turned_null)) {
    cat("\n")
    for (s in x$turned_null) cat("  ", s, "\n", sep = "")
  } else {
    cat("\n  No control turned from failing to null along the ladder.\n")
  }
  invisible(x)
}


#' Care-Process Collider Check (Design Stage, No Outcome Access)
#'
#' Regresses each care-process indicator on treatment conditional on a
#' conditioning set, reporting the coefficient in SD units of the indicator.
#' A process indicator (whether a measurement was taken) that treatment
#' predicts conditionally is a candidate collider: conditioning on it in the
#' adjustment set can open a path between treatment and unmeasured severity.
#' Geography and physiology variables can be included for contrast; the
#' check reads on the process class.
#'
#' @param lock A `cleanroom_lock`.
#' @param indicators Character vector of candidate process-indicator
#'   columns.
#' @param condition_on Character vector of conditioning columns (severity,
#'   site). Defaults to the lock covariates minus the indicators.
#' @param d_threshold Flag threshold in SD units. Default 0.10.
#' @return A data.frame with columns indicator, d_sd_units, p_value,
#'   predicts (logical).
#' @export
check_process_indicators <- function(lock, indicators,
                                     condition_on = NULL,
                                     d_threshold = 0.10) {
  if (!inherits(lock, "cleanroom_lock"))
    stop("`lock` must be a cleanroom_lock object.", call. = FALSE)
  data <- lock$data
  A <- as.numeric(data[[lock$treatment]])
  if (is.null(condition_on))
    condition_on <- setdiff(lock$covariates, indicators)
  condition_on <- intersect(condition_on, names(data))
  Wc <- data[, condition_on, drop = FALSE]
  Wc <- as.data.frame(lapply(Wc, function(x) {
    x <- if (is.numeric(x)) x else as.numeric(as.factor(x))
    if (anyNA(x)) x[is.na(x)] <- stats::median(x, na.rm = TRUE)
    x
  }))
  rows <- list()
  for (v in indicators) {
    if (!v %in% names(data)) next
    y <- suppressWarnings(as.numeric(data[[v]]))
    if (length(unique(y[!is.na(y)])) < 2L) next
    sdy <- stats::sd(y, na.rm = TRUE)
    df <- cbind(data.frame(.y = y, .A = A), Wc)
    fit <- tryCatch(stats::lm(.y ~ ., data = df), error = function(e) NULL)
    if (is.null(fit)) next
    co <- summary(fit)$coefficients
    if (!".A" %in% rownames(co)) next
    d <- co[".A", "Estimate"] / sdy
    rows[[length(rows) + 1L]] <- data.frame(
      indicator = v, d_sd_units = round(d, 4),
      p_value = signif(co[".A", "Pr(>|t|)"], 3),
      predicts = abs(d) >= d_threshold, stringsAsFactors = FALSE)
  }
  out <- do.call(rbind, rows)
  if (is.null(out))
    stop("check_process_indicators: no usable indicator columns.",
         call. = FALSE)
  out[order(-abs(out$d_sd_units)), , drop = FALSE]
}
