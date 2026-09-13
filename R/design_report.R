# The design report: everything the review team sees before outcomes are
# unlocked, in one object. This is the package's version of the blinded
# validity-diagnostics gate of Conover et al. (2025, JAMIA): named
# diagnostics computed while effect estimates stay blinded, an explicit
# reading of whether the design supports the declared estimand, and failed
# comparisons labelled inestimable rather than estimated anyway.

#' Assemble the Pre-Outcome Design Report
#'
#' Collects the support assessment, the estimand feasibility table, the
#' outcome-blind support simulation, and the negative-control ladder into
#' one object with a single reading: which estimands this design supports,
#' and whether the declared primary survives. Everything in it is computable
#' before outcome access. It replaces the 0.1.x dossier as the object a
#' review team reads.
#'
#' @param lock A `cleanroom_lock`, ideally after
#'   [declare_estimand_ladder()].
#' @param support An [assess_support()] result.
#' @param feasibility An [estimand_feasibility()] result.
#' @param simulation Optional [simulate_support()] result.
#' @param nc_ladder Optional [run_negative_control_ladder()] result.
#' @param extra Optional named list of additional design-stage objects to
#'   carry (for example a [check_process_indicators()] table).
#' @return An object of class `design_report` with a `recommendation`
#'   string and a print method.
#' @references Conover MM, Schuemie MJ, et al. (2025) J Am Med Inform
#'   Assoc (objective study validity diagnostics; unblind only when the
#'   prespecified diagnostics pass). Muntner P et al. (2024)
#'   Pharmacoepidemiol Drug Saf 33:e5770.
#' @export
design_report <- function(lock, support, feasibility,
                          simulation = NULL, nc_ladder = NULL,
                          extra = NULL) {
  if (!inherits(lock, "cleanroom_lock"))
    stop("`lock` must be a cleanroom_lock object.", call. = FALSE)
  if (!inherits(support, "support_assessment"))
    stop("`support` must come from assess_support().", call. = FALSE)
  if (!inherits(feasibility, "estimand_feasibility"))
    stop("`feasibility` must come from estimand_feasibility().",
         call. = FALSE)
  ladder <- lock$estimand_ladder

  feas <- feasibility$table
  feasible_estimands <- unique(sub(" \\[.*$", "",
                                   feas$estimand[feas$feasible]))
  primary <- ladder$primary %||% NA_character_
  primary_ok <- !is.na(primary) &&
    primary %in% c(feasible_estimands,
                   if (any(startsWith(feasible_estimands, "trimmed_ATE")))
                     "trimmed_ATE")
  sim_ok <- if (!is.null(simulation) && !is.na(primary)) {
    isTRUE(simulation$feasible[[primary]])
  } else NA

  nc_msg <- if (!is.null(nc_ladder)) {
    flagged_full <- with(nc_ladder$table,
      sum(flagged & cohort == "full cohort", na.rm = TRUE))
    if (length(nc_ladder$turned_null)) {
      paste0(length(nc_ladder$turned_null),
             " control(s) fail on the full cohort and turn null after ",
             "restriction; the restriction is doing the work.")
    } else if (flagged_full > 0) {
      paste0(flagged_full, " control(s) flagged and not resolved by any ",
             "declared restriction; investigate before unblinding.")
    } else "All controls null on every rung."
  } else NULL

  recommendation <- if (is.na(primary)) {
    "No estimand ladder declared; declare one before outcome access."
  } else if (support$verdict == "FAIL") {
    "Extreme non-overlap: the arms do not share a covariate space. No comparison is estimable; the design itself needs to change."
  } else if (primary_ok && (is.na(sim_ok) || isTRUE(sim_ok))) {
    sprintf("The declared primary (%s) is supported (overlap verdict %s)%s. Proceed to estimation.",
            primary, support$verdict,
            if (isTRUE(sim_ok)) " and passes the outcome-blind support map"
            else "")
  } else {
    fb <- ladder$fallbacks[ladder$fallbacks %in% feasible_estimands]
    sprintf("The declared primary (%s) is NOT supported (overlap verdict %s%s). The ladder moves to: %s. The switch is a logged design decision.",
            primary, support$verdict,
            if (isFALSE(sim_ok)) "; it also fails the support map" else "",
            if (length(fb)) paste(fb, collapse = ", then ") else
              "no feasible fallback; the comparison is inestimable")
  }

  out <- list(
    lock_hash = lock$lock_hash,
    contrast = lock$contrast,
    ladder = ladder,
    support = support,
    feasibility = feasibility,
    simulation = simulation,
    nc_ladder = nc_ladder,
    nc_reading = nc_msg,
    design_log = lock$design_log,
    feasible_estimands = feasible_estimands,
    primary_supported = primary_ok,
    recommendation = recommendation,
    created_at = format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
    extra = extra
  )
  class(out) <- "design_report"
  out
}

#' Protocol-Versus-Emulation Table from a Lock
#'
#' Renders the two-column protocol/emulation table of the TARGET reporting
#' guideline for target trial emulations (Cashin, Hansford, Hernan et al.,
#' JAMA 2025) from what the lock already records: the contrast, the
#' eligibility trail in the design log, the outcome, the declared estimand
#' ladder, and the analysis specification. Fields the lock does not carry
#' (the target trial's own protocol column) are left for the analyst to
#' fill, which is the point: the table shows exactly what was specified in
#' software and what was not.
#'
#' @param lock A `cleanroom_lock`, ideally after
#'   [declare_estimand_ladder()].
#' @param protocol Optional named character vector giving the target
#'   trial's protocol entries for the components (names among
#'   `eligibility`, `treatment_strategies`, `assignment`, `follow_up`,
#'   `outcome`, `estimand`, `analysis`).
#' @return A data.frame with columns `component`, `target_trial`,
#'   `emulation`.
#' @references Cashin AG, Hansford HJ, Hernan MA, et al. (2025). Guidance
#'   for reporting target trial emulation studies (TARGET). JAMA
#'   334(12):1084-1093.
#' @export
emulation_table <- function(lock, protocol = NULL) {
  if (!inherits(lock, "cleanroom_lock"))
    stop("`lock` must be a cleanroom_lock object.", call. = FALSE)
  p <- function(k) if (!is.null(protocol) && k %in% names(protocol))
    unname(protocol[[k]]) else "(protocol entry to be supplied)"
  ct <- lock$contrast
  elig <- if (!is.null(lock$design_log)) {
    rows <- lock$design_log[lock$design_log$type %in%
                              c("contrast", "cohort", "outcome"), ]
    if (nrow(rows)) paste(rows$note, collapse = " ") else NA_character_
  } else NA_character_
  lad <- lock$estimand_ladder
  data.frame(
    component = c("Eligibility criteria", "Treatment strategies",
                  "Assignment procedures", "Follow-up period", "Outcome",
                  "Causal contrast (estimand)", "Analysis plan"),
    target_trial = c(p("eligibility"), p("treatment_strategies"),
                     p("assignment"), p("follow_up"), p("outcome"),
                     p("estimand"), p("analysis")),
    emulation = c(
      elig %||% "(see cohort construction)",
      if (!is.null(ct)) sprintf("%s (treated: %s; control: %s)",
                                ct$label, paste(ct$treated, collapse = " + "),
                                paste(ct$control, collapse = " + "))
      else sprintf("binary %s", lock$treatment),
      "Observed assignment; propensity score on the locked covariates, support assessed before outcome access",
      "(as recorded in the outcome definition)",
      lock$outcome,
      if (!is.null(lad)) sprintf(
        "primary %s; fallbacks %s; trigger %s (pre-registered)",
        lad$primary, paste(lad$fallbacks, collapse = " -> "), lad$trigger)
      else "(no ladder declared)",
      sprintf("TMLE with SuperLearner (%d locked covariates); support verdict and implausibility flags travel with every estimate",
              length(lock$covariates))),
    stringsAsFactors = FALSE)
}

#' @export
print.design_report <- function(x, ...) {
  cat("Design report (pre-outcome)\n")
  if (!is.null(x$contrast))
    cat("  Contrast: ", x$contrast$label, "\n", sep = "")
  s <- x$support$summary
  cat(sprintf("  Support: %s; %.1f%% outside [%.2f, %.2f]; max weight %.0f; ESS %d / %d\n",
              x$support$verdict, s$pct_outside_band, x$support$band[1],
              x$support$band[2], s$max_iptw_weight,
              round(s$ess_treated), round(s$ess_control)))
  cat("  Feasible estimands: ",
      if (length(x$feasible_estimands))
        paste(x$feasible_estimands, collapse = ", ") else "none", "\n",
      sep = "")
  if (!is.null(x$ladder))
    cat(sprintf("  Declared ladder: %s -> %s (trigger %s)\n",
                x$ladder$primary,
                paste(x$ladder$fallbacks, collapse = " -> "),
                x$ladder$trigger))
  if (!is.null(x$simulation)) {
    fes <- names(x$simulation$feasible)[x$simulation$feasible]
    cat("  Support map passes: ",
        if (length(fes)) paste(fes, collapse = ", ") else "none",
        if (isTRUE(x$simulation$demo)) "  (demonstration-only replicate count)" else "",
        "\n", sep = "")
  }
  if (!is.null(x$nc_reading)) cat("  Negative controls: ", x$nc_reading,
                                  "\n", sep = "")
  cat("\n  ", x$recommendation, "\n", sep = "")
  invisible(x)
}
