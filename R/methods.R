.short <- function(h) substr(h, 1, 12)

#' @export
print.cr_plan <- function(x, ...) {
  cat("Analysis plan (", x$outcome_type, " outcome: ", paste(x$outcome, collapse = ", "), ")\n", sep = "")
  cat("  Ladder:      ", paste(x$estimands, collapse = " > "), "\n")
  cat("  Candidates:  ", paste(vapply(.candidates(x), `[[`, "", "id"), collapse = ", "), "\n")
  cat("  Tolerance:    bias ", x$tolerance$bias, ", coverage ", x$tolerance$coverage, "\n", sep = "")
  cat("  K = ", x$K, ", V = ", x$V, ", reps = ", x$reps, " (up to ", x$max_reps,
      " for unresolved cells), seed = ", x$seed, "\n", sep = "")
  invisible(x)
}

#' @export
print.cr_lock <- function(x, ...) {
  a <- x$data[[x$treatment]]
  cat("Analysis lock ", .short(x$lock_hash), " (cleanTMLE ", x$version, ")\n", sep = "")
  cat("  ", nrow(x$data), " rows: ", sum(a == 1), " treated, ", sum(a == 0), " control; ",
      length(x$covariates), " covariates\n", sep = "")
  invisible(x)
}

#' @export
print.cr_design <- function(x, ...) {
  cat("Design: overlap ", x$overlap_grade, " (", round(100 * x$share_outside_band, 1),
      "% outside [", x$band[1], ", ", x$band[2], "])\n", sep = "")
  print(x$ess, row.names = FALSE)
  if (nrow(x$edge)) { cat("Edge check:\n"); print(x$edge, row.names = FALSE) }
  invisible(x)
}

#' @export
print.cr_simulation <- function(x, ...) {
  cat("Outcome-free simulation: ", x$reps, " repetitions, K = ", x$K, "\n", sep = "")
  if (!is.null(x$extension) && nrow(x$extension))
    cat("  Unresolved cells extended in ", nrow(x$extension), " batch(es), to at most ",
        max(x$metrics$reps), " repetitions (max_reps = ", x$max_reps, ")\n", sep = "")
  print(x$verdict, row.names = FALSE)
  if (any(x$verdict$unresolved %in% TRUE))
    cat("  unresolved = TRUE: the decision (status or selected candidate) rests on",
        "cells still unresolved after their final repetitions.\n")
  if (!is.null(x$decision_note)) cat(strwrap(x$decision_note, prefix = "  "), sep = "\n")
  cat(strwrap(x$scope, prefix = "  "), sep = "\n")
  invisible(x)
}

#' @export
print.cr_nc_ladder <- function(x, ...) {
  if (!x$declared) { cat("No negative controls declared.\n"); return(invisible(x)) }
  cat("Negative-control ladder: verdict ", x$verdict, "\n", sep = "")
  print(x$by_rung, row.names = FALSE)
  invisible(x)
}

#' @export
print.cr_dossier <- function(x, ...) {
  dec <- x$decision
  cat("Design dossier ", .short(x$dossier_hash), " for lock ", .short(x$lock_hash), "\n", sep = "")
  cat("  Primary estimand: ", if (is.na(dec$primary)) "none feasible" else dec$primary, "\n", sep = "")
  flag <- ifelse(unlist(dec$unresolved[names(dec$status)]) %in% TRUE, " (unresolved)", "")
  cat("  Status: ", paste0(names(dec$status), " ", unlist(dec$status), flag, collapse = "; "),
      "\n", sep = "")
  cat("  Negative controls: ", dec$nc_verdict, "\n", sep = "")
  invisible(x)
}

#' @export
print.cr_unblinded <- function(x, ...) {
  a <- x$approval
  cat("Unblinded by ", a$approved_by, " at ", a$approved_at, "\n", sep = "")
  if (!is.na(a$override)) cat("  Override: ", a$override, "\n", sep = "")
  cat("  Design rows without an outcome: ", a$n_design_without_outcome, "\n", sep = "")
  invisible(x)
}

#' @export
print.cr_estimate <- function(x, ...) {
  cat(sprintf("%s (%s%s): %.4f (95%% CI %.4f, %.4f), n = %d\n", x$estimand, x$method,
              if (is.na(x$candidate)) "" else paste0(", ", x$candidate),
              x$estimate, x$ci_lower, x$ci_upper, as.integer(x$n)))
  if (!is.na(x$note)) cat("  ", x$note, "\n", sep = "")
  if (isTRUE(x$implausible)) cat("  IMPLAUSIBLE: ", x$implausible_reason, "\n", sep = "")
  invisible(x)
}

.evalue <- function(x) {
  v <- c(x$risk1, x$risk0, x$ci_lower, x$ci_upper)
  if (!requireNamespace("EValue", quietly = TRUE) || length(v) != 4L || !all(is.finite(v)) ||
      x$risk1 <= 0 || x$risk0 <= 0)
    return(NULL)
  rr <- x$risk1 / x$risk0
  lo <- max((x$risk0 + x$ci_lower) / x$risk0, 1e-6)
  hi <- max((x$risk0 + x$ci_upper) / x$risk0, 1e-6)
  ev <- suppressMessages(EValue::evalues.RR(est = rr, lo = lo, hi = hi))
  c(point = unname(ev["E-values", "point"]),
    ci = suppressWarnings(min(ev["E-values", c("lower", "upper")], na.rm = TRUE)))
}

#' @export
summary.cr_estimate <- function(object, ...) {
  structure(list(estimate = object, evalue = .evalue(object)), class = "summary.cr_estimate")
}

#' @export
print.summary.cr_estimate <- function(x, ...) {
  print(x$estimate)
  e <- x$estimate
  cat(sprintf("  Risks: %.4f treated, %.4f control (target population n = %d)\n",
              e$risk1, e$risk0, as.integer(e$n_population)))
  if (!is.null(x$evalue))
    cat(sprintf("  E-value: %.2f for the point estimate, %.2f for the confidence limit\n",
                x$evalue["point"], x$evalue["ci"]))
  invisible(x)
}

#' @export
as.data.frame.cr_estimate <- function(x, ...) {
  data.frame(estimand = x$estimand, method = x$method, candidate = x$candidate,
             estimate = x$estimate, se = x$se, ci_lower = x$ci_lower, ci_upper = x$ci_upper,
             p_value = x$p_value, n = x$n, prespecified = x$prespecified,
             stress_tested = x$stress_tested, stringsAsFactors = FALSE)
}

#' @exportS3Method generics::tidy
tidy.cr_estimate <- function(x, ...) {
  data.frame(term = x$estimand, estimate = x$estimate, std.error = x$se,
             conf.low = x$ci_lower, conf.high = x$ci_upper, p.value = x$p_value,
             method = x$method, stringsAsFactors = FALSE)
}

#' @export
as.data.frame.cr_simulation <- function(x, ...) x$metrics

#' @export
as.data.frame.cr_nc_ladder <- function(x, ...) x$table

#' @export
plot.cr_design <- function(x, ...) {
  ggplot2::ggplot(data.frame(g = x$g), ggplot2::aes(x = .data$g)) +
    ggplot2::geom_histogram(bins = 40) +
    ggplot2::geom_vline(xintercept = x$band, linetype = 2) +
    ggplot2::labs(x = "Cross-fitted propensity score", y = "Rows")
}

#' @export
plot.cr_simulation <- function(x, ...) {
  m <- x$metrics
  ggplot2::ggplot(m, ggplot2::aes(x = .data$candidate, y = .data$coverage,
                                  colour = .data$surface)) +
    ggplot2::geom_point(position = ggplot2::position_dodge(width = 0.5)) +
    ggplot2::facet_wrap(~estimand) +
    ggplot2::labs(x = "Candidate", y = "Interval coverage", colour = "Surface") +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))
}
