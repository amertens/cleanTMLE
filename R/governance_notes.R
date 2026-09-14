# Governance notes: the estimand and sensitivity-plan declarations that
# create_analysis_lock() accepts as arguments, plus the event-process
# classification, target-population, missing-data, coherence-check, and
# cumulative-risk reporting formatters. These record and format the
# analytic specification a clean-room review requires; they estimate
# nothing and depend only on base R. They lived in the short-lived
# companion package cleanroomGov during 0.2.x and folded back in 0.3.0.

# ---------------------------------------------------------------------
# Lock declarations (reached through create_analysis_lock arguments)
# ---------------------------------------------------------------------

#' Attach an Estimand Description to an Analysis Lock
#'
#' Enriches a \code{cleanroom_lock} with a structured description of the
#' causal and statistical estimand. The exported path is the
#' \code{estimand} argument of [create_analysis_lock()]; this engine also
#' accepts a bare character string, which becomes the
#' \code{description} field. Estimand metadata is descriptive and is not
#' part of the lock fingerprint.
#'
#' @param lock A \code{cleanroom_lock}.
#' @param estimand A named list with any of \code{description},
#'   \code{population}, \code{treatment_strategies},
#'   \code{outcome_label}, \code{followup}, \code{contrast},
#'   \code{statistical_estimand}; or a single character description.
#'
#' @return The lock with element \code{estimand}.
#' @keywords internal
attach_estimand <- function(lock, estimand) {
  if (!inherits(lock, "cleanroom_lock"))
    stop("`lock` must be a cleanroom_lock object.", call. = FALSE)
  if (is.character(estimand) && length(estimand) == 1L)
    estimand <- list(description = estimand)
  if (!is.list(estimand))
    stop("`estimand` must be a named list (or a single character ",
         "description).", call. = FALSE)
  known <- c("description", "population", "treatment_strategies",
             "outcome_label", "followup", "contrast",
             "statistical_estimand")
  extra <- setdiff(names(estimand), known)
  if (length(extra))
    stop("Unknown estimand field(s): ", paste(extra, collapse = ", "),
         ". Known fields: ", paste(known, collapse = ", "), ".",
         call. = FALSE)
  lock$estimand <- estimand[intersect(known, names(estimand))]
  if (is.null(lock$estimand$contrast))
    lock$estimand$contrast <- "risk_difference"
  lock
}

#' Declare Sensitivity Analysis Plans on an Analysis Lock
#'
#' Attaches named sensitivity-analysis plans to the lock. The exported
#' path is the \code{sensitivity_plans} argument of
#' [create_analysis_lock()]: a named list whose element names are the
#' plan labels and whose elements are lists with \code{description}
#' and \code{settings}. Plans are metadata and are not part of the lock
#' fingerprint.
#'
#' @param lock A \code{cleanroom_lock}.
#' @param plans A named list of plans; each element may carry
#'   \code{description} (character) and \code{settings} (named list).
#'
#' @return The lock with element \code{sensitivity_plans}.
#' @keywords internal
declare_sensitivity_plan <- function(lock, plans) {
  if (!inherits(lock, "cleanroom_lock"))
    stop("`lock` must be a cleanroom_lock object.", call. = FALSE)
  if (!is.list(plans) || length(plans) == 0L || is.null(names(plans)) ||
      any(!nzchar(names(plans))))
    stop("`sensitivity_plans` must be a named list; each name is the ",
         "plan label.", call. = FALSE)
  if (is.null(lock$sensitivity_plans)) lock$sensitivity_plans <- list()
  for (label in names(plans)) {
    p <- plans[[label]]
    if (!is.list(p)) p <- list(description = as.character(p))
    lock$sensitivity_plans[[label]] <- list(
      label       = label,
      description = p$description,
      settings    = p$settings %||% list(),
      declared_at = Sys.time())
  }
  lock
}

# ---------------------------------------------------------------------
# A. Event-process classification table
# ---------------------------------------------------------------------

#' Classify Post-Baseline Events and Processes for Estimand Reporting
#'
#' Builds a structured table that classifies each post-baseline event or
#' process by its role in the causal estimand and its operational role in
#' estimation. Roles correspond to ICH E9(R1) intercurrent-event strategies.
#' The table supports transparent reporting; it does not estimate anything.
#'
#' @param events A list of named lists or a data frame describing each
#'   event. Each entry should provide \code{event_name},
#'   \code{event_variable}, and \code{role_in_estimand}.
#' @param censoring_handling Optional character vector of censoring-handling
#'   approaches the analysis plan provides.
#' @param outcome_is_fatal Logical; \code{TRUE} when the primary outcome is
#'   fatal. Default \code{FALSE}.
#'
#' @return A data frame of class \code{cleantmle_event_process}.
#'
#' @examples
#' clean_event_process_table(list(
#'   list(event_name = "Primary event", event_variable = "Y",
#'        role_in_estimand = "outcome event"),
#'   list(event_name = "Loss to follow-up", event_variable = "lost",
#'        role_in_estimand = "administrative censoring",
#'        primary_handling = "IPCW")
#' ), censoring_handling = "IPCW")
#'
#' @export
clean_event_process_table <- function(events,
                                      censoring_handling = NULL,
                                      outcome_is_fatal = FALSE) {
  if (is.data.frame(events)) events <- split(events, seq_len(nrow(events)))
  if (!is.list(events) || length(events) == 0L)
    stop("`events` must be a non-empty list or data frame.", call. = FALSE)

  required_roles <- c("outcome event", "competing event",
                      "composite outcome component",
                      "administrative censoring", "informative censoring",
                      "artificial censoring due to protocol deviation",
                      "treatment discontinuation", "treatment switching",
                      "intercurrent event", "missing outcome process",
                      "structurally undefined outcome")

  field <- function(x, key, default = NA_character_) {
    v <- x[[key]]
    if (is.null(v) || length(v) == 0L) default else as.character(v)
  }

  rows <- lapply(events, function(e) {
    if (is.null(e$event_name) || is.null(e$event_variable))
      stop("Every event entry needs `event_name` and `event_variable`.",
           call. = FALSE)
    role <- field(e, "role_in_estimand")
    if (!is.na(role) && !role %in% required_roles)
      warning("Unrecognised role_in_estimand '", role,
              "'. Recognised roles: ",
              paste(required_roles, collapse = "; "), call. = FALSE)
    data.frame(
      event_name              = field(e, "event_name"),
      event_variable          = field(e, "event_variable"),
      event_timing_variable   = field(e, "event_timing_variable"),
      role_in_estimand        = role,
      role_in_estimation      = field(e, "role_in_estimation"),
      ICH_E9R1_strategy       = field(e, "ICH_E9R1_strategy"),
      justification           = field(e, "justification"),
      primary_handling        = field(e, "primary_handling"),
      sensitivity_handling    = field(e, "sensitivity_handling"),
      affects_identification_assumption = field(
        e, "affects_identification_assumption"),
      notes                   = field(e, "notes"),
      stringsAsFactors = FALSE)
  })
  out <- do.call(rbind, rows)
  rownames(out) <- NULL

  cens_rows <- grepl("censoring", out$role_in_estimand, ignore.case = TRUE)
  if (any(cens_rows) &&
      (is.null(censoring_handling) || length(censoring_handling) == 0L) &&
      all(is.na(out$primary_handling[cens_rows]))) {
    warning("clean_event_process_table: at least one event is classified ",
            "as censoring but no censoring handling is documented in ",
            "`censoring_handling` or `primary_handling`.", call. = FALSE)
  }

  death_as_cens <- grepl("death|deceas|mortal", out$event_name,
                         ignore.case = TRUE) & cens_rows
  if (any(death_as_cens) && !isTRUE(outcome_is_fatal)) {
    needs_justif <- death_as_cens & is.na(out$justification)
    if (any(needs_justif))
      warning("clean_event_process_table: death is classified as ",
              "censoring while the outcome is not declared fatal. ",
              "Provide a `justification` for that choice.", call. = FALSE)
  }

  comp_rows <- grepl("competing", out$role_in_estimand, ignore.case = TRUE)
  if (any(comp_rows & cens_rows)) {
    needs_estimand_just <- (comp_rows & cens_rows) & is.na(out$justification)
    if (any(needs_estimand_just))
      warning("clean_event_process_table: at least one event is marked ",
              "both as a competing event and as censoring without an ",
              "explicit estimand justification. Add `justification`.",
              call. = FALSE)
  }

  class(out) <- c("cleantmle_event_process", class(out))
  attr(out, "censoring_handling") <- censoring_handling
  out
}

#' @export
print.cleantmle_event_process <- function(x, ...) {
  cat("Event-process classification table\n")
  cat("==================================\n")
  print.data.frame(x[, c("event_name", "role_in_estimand",
                          "primary_handling", "ICH_E9R1_strategy"),
                      drop = FALSE], row.names = FALSE)
  invisible(x)
}


# ---------------------------------------------------------------------
# B. Event-process coherence check (competing-risk components)
# ---------------------------------------------------------------------

#' Check Coherence of Cumulative Incidence Components
#'
#' Verifies that, at each time point, the event-of-interest and competing-event
#' cumulative incidences sum to the composite event risk within tolerance.
#'
#' @param risks A data frame with columns \code{time_point},
#'   \code{event_of_interest_risk}, \code{competing_event_risk}, and
#'   \code{composite_risk}.
#' @param tolerance Numeric tolerance for the coherence check. Default
#'   \code{0.005}.
#'
#' @return A data frame inheriting from \code{cleantmle_risk_check}.
#'
#' @examples
#' df <- data.frame(time_point = c(180, 365),
#'                  event_of_interest_risk = c(0.05, 0.09),
#'                  competing_event_risk = c(0.03, 0.05),
#'                  composite_risk = c(0.08, 0.14))
#' clean_check_event_processes(df)
#'
#' @export
clean_check_event_processes <- function(risks, tolerance = 0.005) {
  if (!is.data.frame(risks))
    stop("`risks` must be a data frame.", call. = FALSE)
  required <- c("time_point", "event_of_interest_risk",
                "competing_event_risk", "composite_risk")
  miss <- setdiff(required, names(risks))
  if (length(miss) > 0L)
    stop("`risks` is missing required columns: ",
         paste(miss, collapse = ", "), call. = FALSE)
  if (!is.numeric(tolerance) || length(tolerance) != 1L || tolerance < 0)
    stop("`tolerance` must be a single non-negative number.", call. = FALSE)

  out <- risks
  out$component_sum <- out$event_of_interest_risk + out$competing_event_risk
  out$discrepancy   <- out$component_sum - out$composite_risk
  out$flag <- ifelse(is.na(out$composite_risk), "MISSING_COMPOSITE",
              ifelse(is.na(out$discrepancy), "MISSING_COMPONENT",
              ifelse(out$component_sum > out$composite_risk + tolerance,
                     "COMPONENT_EXCEEDS_COMPOSITE",
              ifelse(abs(out$discrepancy) > tolerance,
                     "INCOHERENT", "OK"))))

  if (any(out$flag == "COMPONENT_EXCEEDS_COMPOSITE", na.rm = TRUE))
    warning("clean_check_event_processes: at one or more time points the ",
            "sum of event-of-interest and competing-event risks exceeds ",
            "the composite risk by more than the tolerance.",
            call. = FALSE)
  if (any(out$flag == "MISSING_COMPOSITE", na.rm = TRUE))
    warning("clean_check_event_processes: composite_risk is missing at ",
            "one or more time points.", call. = FALSE)

  class(out) <- c("cleantmle_risk_check", class(out))
  out
}


# ---------------------------------------------------------------------
# C. Target-population declaration
# ---------------------------------------------------------------------

#' Declare the Target Population for the Estimand
#'
#' Forces an explicit declaration of the population the estimand refers to,
#' the reference group, the restriction rule (when applicable), and the
#' implications for interpretation. Produces a record; changes no analysis.
#'
#' @param target_population One of \code{"full eligible population"},
#'   \code{"treated population"}, \code{"untreated population"},
#'   \code{"overlap population"}, \code{"trial-eligible subset"},
#'   \code{"positivity-supported restricted population"}, or
#'   \code{"user-defined population"}.
#' @param reference_group Character; the reference treatment strategy.
#' @param restriction_rule Character; description of the restriction rule.
#' @param rationale Character; brief justification.
#' @param implications_for_interpretation Character; interpretation impact.
#' @param positivity_assessment_link Optional character; reference for the
#'   positivity assessment motivating the restriction.
#'
#' @return A list of class \code{cleantmle_target_population}.
#'
#' @export
clean_target_population <- function(target_population,
                                    reference_group = NULL,
                                    restriction_rule = NULL,
                                    rationale = NULL,
                                    implications_for_interpretation = NULL,
                                    positivity_assessment_link = NULL) {
  valid <- c("full eligible population", "treated population",
             "untreated population", "overlap population",
             "trial-eligible subset",
             "positivity-supported restricted population",
             "user-defined population")
  target_population <- match.arg(target_population, valid)

  if (target_population %in% c("positivity-supported restricted population",
                                "user-defined population") &&
      (is.null(restriction_rule) || !nzchar(restriction_rule)))
    warning("clean_target_population: target_population is '",
            target_population, "' but no `restriction_rule` is given.",
            call. = FALSE)

  if (target_population %in% c("treated population",
                                "untreated population") &&
      is.null(reference_group))
    warning("clean_target_population: target_population is '",
            target_population, "'. Describe the contrast as ATT or ATU ",
            "(not ATE) and supply `reference_group`.", call. = FALSE)

  if (is.null(reference_group))
    warning("clean_target_population: `reference_group` is absent. ",
            "Risk differences and risk ratios require a reference.",
            call. = FALSE)

  out <- list(
    target_population               = target_population,
    reference_group                 = reference_group,
    restriction_rule                = restriction_rule,
    rationale                       = rationale,
    implications_for_interpretation = implications_for_interpretation,
    positivity_assessment_link      = positivity_assessment_link)
  class(out) <- c("cleantmle_target_population", class(out))
  out
}

#' @export
print.cleantmle_target_population <- function(x, ...) {
  cat("Target population declaration\n")
  cat("=============================\n")
  cat(sprintf("Target population: %s\n", x$target_population))
  cat(sprintf("Reference group:   %s\n",
              if (is.null(x$reference_group)) "<unspecified>"
              else x$reference_group))
  if (!is.null(x$restriction_rule))
    cat(sprintf("Restriction rule:  %s\n", x$restriction_rule))
  if (!is.null(x$rationale))
    cat(sprintf("Rationale:         %s\n", x$rationale))
  if (!is.null(x$implications_for_interpretation))
    cat(sprintf("Interpretation:    %s\n",
                x$implications_for_interpretation))
  invisible(x)
}


# ---------------------------------------------------------------------
# D. Missing-data plan
# ---------------------------------------------------------------------

#' Declare a Missing-Data Plan
#'
#' Records how each missing-data process (baseline-covariate, treatment,
#' outcome missingness, censoring, structurally undefined outcomes, linkage
#' failure) is handled, pairing each with an identification assumption and a
#' nuisance model, and separating missingness from censoring conceptually.
#'
#' @param processes A named list. Each element supplies \code{timing},
#'   \code{presumed_mechanism}, \code{handling_primary}, and optionally
#'   \code{handling_sensitivity}, \code{nuisance_model_used},
#'   \code{identification_assumption}, and \code{notes}. The element name
#'   becomes \code{variable_or_process}.
#'
#' @return A data frame of class \code{cleantmle_missing_plan}.
#'
#' @export
clean_missing_data_plan <- function(processes) {
  if (!is.list(processes) || length(processes) == 0L ||
      is.null(names(processes)))
    stop("`processes` must be a named list.", call. = FALSE)

  field <- function(x, key, default = NA_character_) {
    v <- x[[key]]
    if (is.null(v) || length(v) == 0L) default else as.character(v)
  }

  rows <- lapply(names(processes), function(nm) {
    e <- processes[[nm]]
    data.frame(
      variable_or_process        = nm,
      timing                     = field(e, "timing"),
      presumed_mechanism         = field(e, "presumed_mechanism"),
      handling_primary           = field(e, "handling_primary"),
      handling_sensitivity       = field(e, "handling_sensitivity"),
      nuisance_model_used        = field(e, "nuisance_model_used"),
      identification_assumption  = field(e, "identification_assumption"),
      notes                      = field(e, "notes"),
      stringsAsFactors = FALSE)
  })
  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  class(out) <- c("cleantmle_missing_plan", class(out))
  out
}

#' @export
print.cleantmle_missing_plan <- function(x, ...) {
  cat("Missing-data plan\n")
  cat("=================\n")
  print.data.frame(x[, c("variable_or_process", "presumed_mechanism",
                          "handling_primary"), drop = FALSE],
                    row.names = FALSE)
  invisible(x)
}


# ---------------------------------------------------------------------
# E. Cumulative-risk reporting table
# ---------------------------------------------------------------------

#' Reporting Table for Prespecified Cumulative-Risk Estimands
#'
#' Builds a reporting table for cumulative risks at prespecified time points,
#' with risk differences, risk ratios, and event / censoring counts. Inputs
#' are computed by an estimator and passed in; this helper only formats.
#'
#' @param rows A list of named lists or a data frame; each entry is a
#'   (time_point, treatment_strategy) cell. Required per row: \code{time_point},
#'   \code{treatment_strategy}, \code{risk}. Optional: \code{n_observed},
#'   \code{events}, \code{censoring_events}, \code{competing_events},
#'   \code{person_time}, \code{effective_sample_size}, \code{risk_ci_lower},
#'   \code{risk_ci_upper}, \code{risk_difference}, \code{risk_ratio},
#'   \code{estimator}, \code{nuisance_specification}, \code{notes}.
#'
#' @return A data frame of class \code{cleantmle_risk_report}.
#'
#' @export
clean_risk_report_table <- function(rows) {
  if (is.data.frame(rows)) rows <- split(rows, seq_len(nrow(rows)))
  if (!is.list(rows) || length(rows) == 0L)
    stop("`rows` must be a non-empty list or data frame.", call. = FALSE)

  cols <- c("time_point", "treatment_strategy", "n_observed", "events",
            "censoring_events", "competing_events", "person_time",
            "effective_sample_size", "risk", "risk_ci_lower",
            "risk_ci_upper", "risk_difference", "risk_ratio",
            "estimator", "nuisance_specification", "notes")

  field <- function(x, key) {
    v <- x[[key]]
    if (is.null(v) || length(v) == 0L) NA else v
  }

  out <- do.call(rbind, lapply(rows, function(e) {
    if (is.null(e$time_point) || is.null(e$treatment_strategy) ||
        is.null(e$risk))
      stop("Each row needs `time_point`, `treatment_strategy`, and `risk`.",
           call. = FALSE)
    df <- as.data.frame(lapply(cols, function(k) field(e, k)),
                        stringsAsFactors = FALSE)
    names(df) <- cols
    df
  }))
  rownames(out) <- NULL

  numeric_cols <- c("time_point", "n_observed", "events",
                    "censoring_events", "competing_events", "person_time",
                    "effective_sample_size", "risk", "risk_ci_lower",
                    "risk_ci_upper", "risk_difference", "risk_ratio")
  for (cc in numeric_cols)
    out[[cc]] <- suppressWarnings(as.numeric(out[[cc]]))
  for (cc in c("treatment_strategy", "estimator",
                "nuisance_specification", "notes"))
    out[[cc]] <- as.character(out[[cc]])

  class(out) <- c("cleantmle_risk_report", class(out))
  out
}

#' @export
print.cleantmle_risk_report <- function(x, ...) {
  cat("Cumulative-risk reporting table\n")
  cat("===============================\n")
  show <- x[, c("time_point", "treatment_strategy", "n_observed",
                "events", "risk", "risk_ci_lower", "risk_ci_upper",
                "risk_difference", "risk_ratio", "estimator"),
            drop = FALSE]
  print.data.frame(show, row.names = FALSE)
  invisible(x)
}