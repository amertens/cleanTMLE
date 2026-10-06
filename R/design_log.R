#' Export the design log as JSON
#'
#' Assembles the log from the objects it receives: plan, hashes, design
#' summaries, decision, approval, overrides and estimates. Nothing is
#' stored between calls.
#'
#' @param x A [unblind()] result, or a [design_report()] dossier for the
#'   pre-outcome log.
#' @param fit Optional `cr_estimate` or list of them.
#' @param file Optional `.json` path.
#' @return The log as a list (invisibly when `file` is given).
#' @export
export_design_log <- function(x, fit = NULL, file = NULL) {
  if (inherits(x, "cr_unblinded")) {
    .verify_unblinded(x)
    d <- x$dossier
    approval <- x$approval
  } else if (inherits(x, "cr_dossier")) {
    d <- x
    approval <- NULL
  } else {
    stop("`x` must be an unblind() result or a design_report() dossier.", call. = FALSE)
  }
  .verify_dossier(d)
  fits <- if (inherits(fit, "cr_estimate")) list(fit) else fit
  log <- list(
    package_version = d$version, lock_hash = d$lock_hash, dossier_hash = d$dossier_hash,
    dossier_created = d$created, plan = d$plan, decision = d$decision,
    design = d$design[c("overlap", "overlap_grade", "share_outside_band", "ess",
                        "balance", "edge")],
    simulation = d$simulation[c("verdict", "scope", "reps", "K")],
    negative_controls = d$nc[c("by_rung", "verdict")],
    notes = d$notes, approval = approval,
    estimates = lapply(fits, function(f) unclass(f)))
  if (is.null(file)) return(log)
  jsonlite::write_json(log, file, auto_unbox = TRUE, pretty = TRUE, digits = NA,
                       null = "null", na = "null")
  invisible(log)
}
