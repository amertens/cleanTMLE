.role <- function(plan, role) {
  r <- plan$roles[[role]]
  if (is.null(r)) NA_character_ else paste(r, collapse = ", ")
}

# The grading rule of the negative controls in one plain phrase, for the
# Muntner rationale. Reads the plan summary's nc_criteria.
.nc_rule_text <- function(crit) {
  band <- crit$null_band
  rng <- sprintf("[%s, %s]", band[1], band[2])
  switch(crit$rule %||% "point",
         interval = paste0("the interval rule (GO when the 95% interval lies inside the null band ",
                           rng, ", STOP when it lies entirely outside, FLAG when it is too ",
                           "imprecise to tell)."),
         ci = paste0("the whole-interval rule (GO only when the 95% interval lies inside the null band ",
                     rng, ")."),
         paste0("the point-estimate rule (GO when the point estimate lies inside the null band ",
                rng, ")."))
}

.muntner_log <- function(d, approval, fits) {
  pl <- d$plan
  an <- .role(pl, "analyst")
  row <- function(date, stage, issue, decision, rationale, decided_by = an)
    data.frame(date = date, stage = stage, issue = issue, decision = decision,
               rationale = rationale, decided_by = decided_by, stringsAsFactors = FALSE)
  n_surf <- length(unique(d$simulation$metrics$surface))
  tol <- pl$tolerance
  v <- d$simulation$verdict
  rows <- list(
    row(d$lock_created, "Lock", "Analysis plan and design data",
        paste("Locked, hash", substr(d$lock_hash, 1, 12)),
        "The plan and the outcome-free design data were fingerprinted before any outcome access."),
    row(d$created, "Design assessment", "Propensity-score overlap",
        sprintf("Overlap %s; %.1f%% outside [%s, %s]", d$design$overlap_grade,
                100 * d$design$share_outside_band, d$design$band[1], d$design$band[2]),
        "Graded on the cross-fitted propensity score."))
  for (i in seq_len(nrow(v))) {
    e <- v$estimand[i]
    dec <- if (v$status[i] == "feasible") sprintf("Feasible; candidate %s", v$selected[i]) else
      "Infeasible"
    why <- sprintf("Feasible requires bias within %s and coverage at least %s on every one of %d outcome surfaces",
                   tol$bias, tol$coverage, n_surf)
    if (!is.null(tol$ci_halfwidth))
      why <- paste0(why, sprintf(", and interval half-width at most %s", tol$ci_halfwidth))
    if (isTRUE(v$unresolved[i]))
      why <- paste0(why, "; decided on point estimates with cells still unresolved at the cap")
    rows[[length(rows) + 1L]] <- row(d$created, "Outcome-free simulation",
                                     sprintf("Can the design deliver the %s?", e), dec,
                                     paste0(why, "."))
  }
  if (!is.na(d$decision$nc_verdict)) {
    rows[[length(rows) + 1L]] <- row(
      d$created, "Negative controls",
      sprintf("Residual confounding check on the %s (the primary estimand)",
              d$decision$nc_estimand),
      d$decision$nc_verdict,
      paste("Graded on the full cohort, for the estimand the analysis reports, by",
            .nc_rule_text(pl$nc_criteria)))
  } else if (!is.null(pl$negative_controls) && is.na(d$decision$primary)) {
    rows[[length(rows) + 1L]] <- row(
      d$created, "Negative controls", "Residual confounding check (no feasible estimand)",
      "No verdict",
      "No estimand is feasible, so no estimand's control verdict gates unblinding.")
  }
  rows[[length(rows) + 1L]] <- row(d$created, "Pre-outcome decision", "Primary estimand",
                                   if (is.na(d$decision$primary)) "None feasible" else
                                     d$decision$primary,
                                   "The first feasible estimand in the declared ladder.")
  if (!is.null(approval)) {
    if (!is.na(approval$override))
      rows[[length(rows) + 1L]] <- row(approval$approved_at, "Approval", "Override",
                                       paste("Proceed despite:",
                                             paste(approval$blocked_by, collapse = "; ")),
                                       approval$override, approval$approved_by)
    rows[[length(rows) + 1L]] <- row(approval$approved_at, "Approval", "Outcome access",
                                     "Outcome admitted", "The review team approved the dossier.",
                                     approval$approved_by)
  }
  for (f in fits)
    rows[[length(rows) + 1L]] <- row(NA_character_, "Estimation", f$estimand,
                                     sprintf("%s by %s: %.4f (95%% CI %.4f, %.4f)", f$estimand,
                                             f$method, f$estimate, f$ci_lower, f$ci_upper),
                                     if (isTRUE(f$prespecified)) "The prespecified primary." else
                                       "Not the prespecified primary.")
  do.call(rbind, rows)
}

#' Export the design log
#'
#' Assembles the log from the objects it receives: plan, hashes, design
#' summaries, decision, approval, overrides and estimates. Nothing is
#' stored between calls.
#'
#' @param x A [unblind()] result, or a [design_report()] dossier for the
#'   pre-outcome log.
#' @param fit Optional `cr_estimate` or list of them.
#' @param file Optional path: `.json` for `format = "json"`, `.csv` for
#'   `format = "muntner"`.
#' @param format `"json"` (the full log) or `"muntner"` (the Table S1
#'   columns of Muntner et al. 2024: date, stage, issue, decision,
#'   rationale, decided_by).
#' @return The log (a list for JSON, a data frame for Muntner), invisibly
#'   when `file` is given.
#' @references Muntner P et al. (2024) Pharmacoepidemiol Drug Saf
#'   33:e5770, supplementary Table S1.
#' @export
export_design_log <- function(x, fit = NULL, file = NULL, format = c("json", "muntner")) {
  format <- match.arg(format)
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
  if (format == "muntner") {
    out <- .muntner_log(d, approval, fits)
    if (is.null(file)) return(out)
    utils::write.csv(out, file, row.names = FALSE)
    return(invisible(out))
  }
  log <- list(
    package_version = d$version, lock_hash = d$lock_hash, dossier_hash = d$dossier_hash,
    lock_created = d$lock_created, dossier_created = d$created, plan = d$plan,
    decision = d$decision,
    design = d$design[intersect(c("overlap", "overlap_grade", "share_outside_band",
                                  "c_statistic", "ess", "balance", "edge", "populations",
                                  "trimmed_profile"), names(d$design))],
    simulation = d$simulation[c("verdict", "scope", "reps", "K")],
    expected_events = d$expected_events_range,
    negative_controls = d$nc[c("by_rung", "verdict")],
    notes = d$notes, approval = approval,
    estimates = lapply(fits, function(f) unclass(f)))
  if (is.null(file)) return(log)
  jsonlite::write_json(log, file, auto_unbox = TRUE, pretty = TRUE, digits = NA,
                       null = "null", na = "null")
  invisible(log)
}
