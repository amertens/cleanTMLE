# Diagnostic helpers for event-process classification, target population,
# missing-data planning, weight diagnostics, competing-risk coherence,
# and cumulative-risk reporting. These helpers do not change the cleanTMLE
# estimator. They are documentation and diagnostic utilities motivated by
# the causalRisk vignette workflow and ICH E9(R1) intercurrent-event
# language.

# ---------------------------------------------------------------------
# A. Event-process classification table
# ---------------------------------------------------------------------

#' Classify Post-Baseline Events and Processes for Estimand Reporting
#'
#' Builds a structured table that classifies each post-baseline event or
#' process by its role in the causal estimand and its operational role in
#' estimation. The classification is meant to support transparent
#' reporting of post-baseline events; it does not itself estimate any
#' causal quantity.
#'
#' Roles correspond to ICH E9(R1) intercurrent-event strategies and to
#' the standard event types encountered in observational
#' comparative-effectiveness analyses: outcome event, competing event,
#' composite outcome component, administrative censoring, informative
#' censoring, artificial censoring due to protocol deviation, treatment
#' discontinuation, treatment switching, intercurrent event, missing
#' outcome process, structurally undefined outcome.
#'
#' @section Clean-room stage: Stage 1a (pre-outcome).
#'
#' @param events A list of named lists or a data frame describing each
#'   event. Each entry should provide \code{event_name},
#'   \code{event_variable}, and \code{role_in_estimand}. The remaining
#'   columns (timing variable, role in estimation, ICH E9(R1) strategy,
#'   justification, primary handling, sensitivity handling, identification
#'   assumption affected, and notes) default to \code{NA} when not given.
#' @param censoring_handling Optional character vector listing the
#'   censoring-handling approaches that the analysis plan provides
#'   (for example \code{"administrative censoring; IPCW"}). Used to
#'   trigger a warning when an event is classified as censoring but no
#'   handling is described.
#' @param outcome_is_fatal Logical; \code{TRUE} when the primary outcome
#'   is fatal (in which case death cannot be classified as censoring
#'   without explicit justification). Default \code{FALSE}.
#'
#' @return A data frame of class \code{cleantmle_event_process} with
#'   columns: \code{event_name}, \code{event_variable},
#'   \code{event_timing_variable}, \code{role_in_estimand},
#'   \code{role_in_estimation}, \code{ICH_E9R1_strategy},
#'   \code{justification}, \code{primary_handling},
#'   \code{sensitivity_handling}, \code{affects_identification_assumption},
#'   and \code{notes}.
#'
#' @examples
#' clean_event_process_table(list(
#'   list(event_name = "Primary event", event_variable = "Y",
#'        role_in_estimand = "outcome event"),
#'   list(event_name = "Death from other causes", event_variable = "death",
#'        role_in_estimand = "competing event",
#'        primary_handling = "competing-risk cumulative incidence"),
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

  # Validation warnings.
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
#' Verifies that, at each time point, the event-of-interest cumulative
#' incidence and the competing-event cumulative incidence sum to the
#' composite event risk within numerical tolerance. Useful as a
#' pre-reporting sanity check when both event-specific and composite
#' risks are reported.
#'
#' @section Clean-room stage: Stage 4 (post-outcome reporting).
#'
#' @param risks A data frame with columns \code{time_point},
#'   \code{event_of_interest_risk}, \code{competing_event_risk}, and
#'   \code{composite_risk}. Optional columns \code{event_of_interest_lo}
#'   / \code{_hi} and similar are passed through unchanged.
#' @param tolerance Numeric tolerance for the coherence check.
#'   Default \code{0.005}.
#'
#' @return A data frame inheriting from \code{cleantmle_risk_check} with
#'   the input columns plus \code{component_sum},
#'   \code{discrepancy}, and \code{flag}.
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
# C. Weight diagnostics
# ---------------------------------------------------------------------

#' Weight Diagnostics for Treatment, Censoring, or Missingness Weights
#'
#' Standard diagnostic summary for inverse-probability weights:
#' percentiles by treatment group and overall, effective sample size,
#' counts and proportions of extreme weights under user-specified
#' thresholds, and (when covariates are supplied) standardised mean
#' differences before and after weighting. The summary supports
#' transparent reporting of overlap and weight-distribution properties;
#' it does not itself verify positivity.
#'
#' @section Clean-room stage: Stage 2 / 4 (diagnostic).
#'
#' @param weights Numeric vector of inverse-probability weights.
#' @param treatment Optional treatment indicator (0/1 or factor) of the
#'   same length as \code{weights}. When supplied, percentiles and ESS
#'   are reported by treatment group as well as overall.
#' @param covariates Optional data frame of baseline covariates; if
#'   supplied, the function returns standardised mean differences before
#'   and after weighting for each covariate.
#' @param max_weight_threshold Numeric threshold; an observation whose
#'   weight exceeds this value is flagged as extreme. Default
#'   \code{10}.
#' @param ess_floor Numeric threshold; the function flags ESS values
#'   below this as a low effective sample size. Default \code{0.3 *
#'   length(weights)}.
#'
#' @return A list of class \code{cleantmle_weight_diag} with elements
#'   \code{percentiles}, \code{ess}, \code{extreme_weights},
#'   \code{flags}, and (when covariates are supplied) \code{smd}.
#'
#' @examples
#' set.seed(1)
#' n <- 400
#' A <- rbinom(n, 1, 0.4)
#' ps <- plogis(0.2 * rnorm(n) + 0.5 * A)
#' w  <- ifelse(A == 1, 1 / ps, 1 / (1 - ps))
#' clean_weight_diagnostics(w, treatment = A)
#'
#' @export
clean_weight_diagnostics <- function(weights, treatment = NULL,
                                      covariates = NULL,
                                      max_weight_threshold = 10,
                                      ess_floor = NULL) {
  if (!is.numeric(weights) || length(weights) == 0L)
    stop("`weights` must be a non-empty numeric vector.", call. = FALSE)
  n <- length(weights)
  if (is.null(ess_floor)) ess_floor <- 0.3 * n
  if (!is.null(treatment) && length(treatment) != n)
    stop("`treatment` must be the same length as `weights`.", call. = FALSE)

  pct_probs <- c(min = 0, p01 = 0.01, p05 = 0.05, p25 = 0.25,
                 median = 0.5, p75 = 0.75, p95 = 0.95, p99 = 0.99,
                 max = 1)

  pct_row <- function(w) {
    q <- stats::quantile(w, probs = pct_probs, na.rm = TRUE, names = FALSE)
    setNames(c(q, mean(w, na.rm = TRUE)),
             c(names(pct_probs), "mean"))
  }

  ess_fn <- function(w) sum(w)^2 / sum(w^2)

  pct_overall <- pct_row(weights)
  ess_overall <- ess_fn(weights)
  pct_by <- NULL; ess_by <- NULL
  if (!is.null(treatment)) {
    pct_by <- t(sapply(unique(treatment),
                       function(a) pct_row(weights[treatment == a])))
    rownames(pct_by) <- as.character(unique(treatment))
    ess_by <- vapply(unique(treatment),
                     function(a) ess_fn(weights[treatment == a]),
                     numeric(1))
    names(ess_by) <- as.character(unique(treatment))
  }

  extreme_n <- sum(weights > max_weight_threshold, na.rm = TRUE)
  extreme_p <- extreme_n / n

  smd_tbl <- NULL
  if (!is.null(covariates)) {
    if (is.null(treatment))
      stop("Provide `treatment` to compute SMDs.", call. = FALSE)
    if (nrow(covariates) != n)
      stop("`covariates` must have the same number of rows as `weights`.",
           call. = FALSE)
    smd_one <- function(x, A, w = NULL) {
      x <- as.numeric(x)
      if (is.null(w)) {
        m1 <- mean(x[A == 1], na.rm = TRUE)
        m0 <- mean(x[A == 0], na.rm = TRUE)
        s1 <- stats::var(x[A == 1], na.rm = TRUE)
        s0 <- stats::var(x[A == 0], na.rm = TRUE)
      } else {
        wm <- function(z, ww) sum(z * ww, na.rm = TRUE) /
                                sum(ww[!is.na(z)], na.rm = TRUE)
        m1 <- wm(x[A == 1], w[A == 1]); m0 <- wm(x[A == 0], w[A == 0])
        s1 <- wm((x[A == 1] - m1)^2, w[A == 1])
        s0 <- wm((x[A == 0] - m0)^2, w[A == 0])
      }
      pooled <- sqrt((s1 + s0) / 2)
      if (is.na(pooled) || pooled == 0) return(0)
      abs(m1 - m0) / pooled
    }
    smd_tbl <- data.frame(
      covariate    = names(covariates),
      smd_unweight = vapply(names(covariates),
                            function(v) smd_one(covariates[[v]], treatment),
                            numeric(1)),
      smd_weighted = vapply(names(covariates),
                            function(v) smd_one(covariates[[v]],
                                                treatment, weights),
                            numeric(1)),
      stringsAsFactors = FALSE)
    rownames(smd_tbl) <- NULL
  }

  flags <- list(
    low_ess          = ess_overall < ess_floor,
    extreme_weights  = extreme_n > 0,
    high_max_weight  = max(weights, na.rm = TRUE) > max_weight_threshold)

  out <- list(
    percentiles      = list(overall = pct_overall, by_treatment = pct_by),
    ess              = list(overall = ess_overall, by_treatment = ess_by),
    extreme_weights  = list(threshold = max_weight_threshold,
                             n = extreme_n, prop = extreme_p),
    smd              = smd_tbl,
    flags            = flags,
    n                = n)
  class(out) <- c("cleantmle_weight_diag", class(out))
  out
}

#' @export
print.cleantmle_weight_diag <- function(x, ...) {
  cat("Weight diagnostics\n")
  cat("==================\n")
  cat(sprintf("N = %d\n", x$n))
  cat(sprintf("Effective sample size (overall): %.1f\n", x$ess$overall))
  if (!is.null(x$ess$by_treatment))
    print(round(x$ess$by_treatment, 1))
  cat(sprintf("Extreme weights (> %.2f): %d (%.1f%%)\n",
              x$extreme_weights$threshold,
              x$extreme_weights$n, 100 * x$extreme_weights$prop))
  cat("\nPercentiles (overall):\n")
  print(round(x$percentiles$overall, 4))
  if (!is.null(x$smd)) {
    cat("\nStandardised mean differences:\n")
    print(round(x$smd[, c("smd_unweight", "smd_weighted"), drop = FALSE], 3))
  }
  if (any(unlist(x$flags))) {
    cat("\nFlags: ")
    cat(paste(names(x$flags)[unlist(x$flags)], collapse = "; "), "\n")
  }
  invisible(x)
}


# ---------------------------------------------------------------------
# D. Target-population declaration
# ---------------------------------------------------------------------

#' Declare the Target Population for the Estimand
#'
#' Forces an explicit declaration of the population to which the estimand
#' refers, the reference group, the restriction rule (when applicable),
#' and the implications for interpretation. The helper does not change
#' the analysis; it produces a structured record that can be archived
#' with the analysis specification.
#'
#' @section Clean-room stage: Stage 1a (pre-outcome).
#'
#' @param target_population One of \code{"full eligible population"},
#'   \code{"treated population"}, \code{"untreated population"},
#'   \code{"overlap population"}, \code{"trial-eligible subset"},
#'   \code{"positivity-supported restricted population"}, or
#'   \code{"user-defined population"}.
#' @param reference_group Character; the reference treatment strategy
#'   used for risk-difference or risk-ratio contrasts.
#' @param restriction_rule Character; description of the restriction
#'   rule (required when \code{target_population} is
#'   \code{"positivity-supported restricted population"} or
#'   \code{"user-defined population"}).
#' @param rationale Character; brief justification for the choice.
#' @param implications_for_interpretation Character; how the choice
#'   affects interpretation of the contrast.
#' @param positivity_assessment_link Optional character; reference or
#'   filename for the positivity assessment that motivated the
#'   restriction rule.
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
# E. Missing-data plan
# ---------------------------------------------------------------------

#' Declare a Missing-Data Plan
#'
#' Records, in a single structured table, how each missing-data process
#' (baseline covariate missingness, treatment / exposure missingness,
#' outcome missingness, censoring or loss to follow-up, structurally
#' undefined outcomes, linkage failure) is handled in the analysis. The
#' helper exists to separate missingness from censoring conceptually so
#' that each process is paired with an identification assumption and a
#' nuisance model.
#'
#' @section Clean-room stage: Stage 1a (pre-outcome).
#'
#' @param processes A named list. Each element corresponds to one
#'   process and supplies \code{timing}, \code{presumed_mechanism},
#'   \code{handling_primary}, optionally \code{handling_sensitivity},
#'   \code{nuisance_model_used}, \code{identification_assumption}, and
#'   \code{notes}. The element name becomes the
#'   \code{variable_or_process} entry.
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
# F. Cumulative-risk reporting table
# ---------------------------------------------------------------------

#' Reporting Table for Prespecified Cumulative-Risk Estimands
#'
#' Builds a reporting table for cumulative risks at prespecified time
#' points, with risk differences, risk ratios, and event / censoring
#' counts. This helper is meant for reporting; it does not estimate
#' risks from raw data. Inputs are typically computed by an estimator
#' (TMLE, IPTW, g-computation) and passed in.
#'
#' @section Clean-room stage: Stage 4 (post-outcome reporting).
#'
#' @param rows A list of named lists or a data frame; each entry
#'   corresponds to a (time_point, treatment_strategy) cell of the
#'   table. Required fields per row: \code{time_point},
#'   \code{treatment_strategy}, \code{risk}. Optional fields:
#'   \code{n_observed}, \code{events}, \code{censoring_events},
#'   \code{competing_events}, \code{person_time},
#'   \code{effective_sample_size}, \code{risk_ci_lower},
#'   \code{risk_ci_upper}, \code{risk_difference}, \code{risk_ratio},
#'   \code{estimator}, \code{nuisance_specification}, \code{notes}.
#'
#' @return A data frame of class \code{cleantmle_risk_report} with the
#'   required and optional columns; missing optional fields are
#'   filled with \code{NA}.
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
