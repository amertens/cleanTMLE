# The sixteen-verb workflow surface (0.3.0). This file carries the
# verbs that rename or consolidate earlier entry points:
#   declare_negative_controls()  (define_negative_control, per control)
#   define_candidates()          (tmle_candidate / expand_tmle_candidate_grid)
#   select_candidate()           (select_tmle_candidate)
#   stress_test()                (run_plasmode_feasibility /
#                                 run_plasmode_dq_stress /
#                                 summarize_dq_degradation)
#   negative_control_ladder()    (run_negative_control_ladder)
# The earlier names remain internal engines; the verbs are the API.


#' Declare Negative Controls on the Lock
#'
#' Registers one or more negative-control variables, with their types,
#' descriptions, and Muntner et al. (2024) confounding domains, and
#' writes the declaration into the design log. Controls declared through
#' `create_analysis_lock(negative_controls = )` can be enriched here
#' with domains; criteria declared at lock creation
#' (`nc_criteria = `) are fingerprinted, while criteria declared here
#' afterwards are recorded in the design log with a timestamp (the
#' fingerprint cannot change after creation).
#'
#' @param lock A `cleanroom_lock`.
#' @param variables Character vector of negative-control column names.
#' @param types Character vector (recycled): `"outcome"` or
#'   `"exposure"`. Default `"outcome"`.
#' @param descriptions Optional character vector (recycled) of
#'   free-text descriptions.
#' @param domains Optional character vector (recycled) of Muntner et
#'   al. (2024) domains: `"confounding_by_indication"`,
#'   `"functional_status"`, `"health_seeking_behavior"`,
#'   `"access_to_healthcare"`, or `"other"`.
#' @param nc_criteria Optional prespecified decision criteria (see
#'   [create_analysis_lock()]). Accepted only when none were declared at
#'   lock creation; recorded in the design log as a post-creation
#'   declaration.
#' @return The lock with the controls registered, the declaration
#'   logged, and (when given) `nc_criteria` set.
#' @examples
#' dat  <- sim_func1(n = 200, seed = 1)
#' lock <- create_analysis_lock(dat, "treatment", "event_24",
#'                              c("age", "sex", "biomarker"), seed = 1)
#' lock <- declare_negative_controls(lock, "nc_outcome",
#'   domains = "health_seeking_behavior",
#'   nc_criteria = list(null_band = 0.02))
#' @export
declare_negative_controls <- function(lock, variables,
                                      types = "outcome",
                                      descriptions = NULL,
                                      domains = NULL,
                                      nc_criteria = NULL) {
  if (!inherits(lock, "cleanroom_lock"))
    stop("`lock` must be a cleanroom_lock object.", call. = FALSE)
  if (!is.character(variables) || !length(variables))
    stop("`variables` must be a non-empty character vector.", call. = FALSE)
  k <- length(variables)
  types <- rep_len(types, k)
  if (!is.null(descriptions)) descriptions <- rep_len(descriptions, k)
  if (!is.null(domains)) domains <- rep_len(domains, k)

  for (i in seq_len(k)) {
    lock <- define_negative_control(
      lock, variables[i], type = types[i],
      description = if (is.null(descriptions)) NULL else descriptions[i],
      domain = if (is.null(domains)) NULL else domains[i])
  }
  lock <- .log_design_decision(
    lock, "negative_controls",
    sprintf("Declared negative control(s): %s%s.",
            paste(variables, collapse = ", "),
            if (is.null(domains)) "" else
              paste0(" (domains: ", paste(domains, collapse = ", "), ")")),
    stage = "Stage 1a (declarations)")

  if (!is.null(nc_criteria)) {
    if (!is.null(lock$nc_criteria))
      stop("nc_criteria were already declared on this lock",
           if (isTRUE(lock$declared_at_creation$nc_criteria))
             " at creation (fingerprinted)" else "",
           "; they cannot be redeclared.", call. = FALSE)
    lock$nc_criteria <- .normalize_nc_criteria(nc_criteria)
    lock <- .log_design_decision(
      lock, "nc_criteria",
      sprintf(paste0("Negative-control criteria declared after lock ",
                     "creation: null band %.4g (%s, %s, min %d per ",
                     "domain, %s)."),
              lock$nc_criteria$null_band, lock$nc_criteria$scale,
              lock$nc_criteria$rule, lock$nc_criteria$min_per_domain,
              lock$nc_criteria$consistency),
      stage = "Stage 1a (declarations)")
  }
  lock
}


#' Define the Candidate Estimator Set
#'
#' Constructs the prespecified TMLE candidate set: either a single
#' candidate (give `candidate_id` plus the specification arguments) or a
#' factorial grid (give `grid`, a named list of axes). Candidates are
#' fixed choices of nuisance libraries, truncation, and cross-fitting;
#' [stress_test()] scores them on simulated outcomes and
#' [select_candidate()] applies the locked selection rule.
#'
#' @param candidate_id Character id for a single candidate (omit when
#'   `grid` is given).
#' @param label Display label; defaults to the id.
#' @param g_library,q_library SuperLearner libraries for the treatment
#'   and outcome nuisances; `q_library` defaults to `g_library`.
#' @param truncation Propensity truncation bound (or a named rule
#'   resolved at fit time).
#' @param ... Further specification fields (estimator, cv_scheme,
#'   variance_method, discrete_sl, screener; see the candidate
#'   internals).
#' @param grid A named list of grid axes: `truncations` (numeric
#'   vector) and `libraries` (named list of character libraries), plus
#'   optionally `estimators`, `cv_schemes`, `variance_methods`,
#'   `discrete_sl`, `screeners`, `max_candidates`.
#' @return An object of class `ct_candidates`: a named list of
#'   candidate specifications.
#' @examples
#' cands <- define_candidates(grid = list(
#'   truncations = c(0.01, 0.05),
#'   libraries   = list(glm = "SL.glm")))
#' length(cands)
#' one <- define_candidates("glm_t01", g_library = "SL.glm",
#'                          truncation = 0.01)
#' @export
define_candidates <- function(candidate_id = NULL, label = candidate_id,
                              g_library = c("SL.glm"), q_library = NULL,
                              truncation = 0.01, ..., grid = NULL) {
  if (!is.null(grid)) {
    if (!is.null(candidate_id))
      stop("Give either `candidate_id` (one candidate) or `grid` (a ",
           "factorial set), not both.", call. = FALSE)
    if (!is.list(grid) || is.null(names(grid)))
      stop("`grid` must be a named list of axes (truncations, ",
           "libraries, ...).", call. = FALSE)
    bad <- setdiff(names(grid),
                   names(formals(expand_tmle_candidate_grid)))
    if (length(bad))
      stop("Unknown grid axis/axes: ", paste(bad, collapse = ", "),
           call. = FALSE)
    out <- do.call(expand_tmle_candidate_grid, grid)
  } else {
    if (is.null(candidate_id))
      stop("Give `candidate_id` for a single candidate, or `grid` for ",
           "a factorial set.", call. = FALSE)
    one <- tmle_candidate(candidate_id = candidate_id, label = label,
                          g_library = g_library, q_library = q_library,
                          truncation = truncation, ...)
    out <- stats::setNames(list(one), candidate_id)
  }
  class(out) <- c("ct_candidates", "list")
  out
}

#' @export
print.ct_candidates <- function(x, ...) {
  cat("Candidate set (", length(x), " candidate",
      if (length(x) != 1) "s", ")\n", sep = "")
  for (cand in x)
    cat(sprintf("  %-16s g: %s | Q: %s | trunc %s\n",
                cand$candidate_id,
                paste(cand$g_library, collapse = "+"),
                paste(cand$q_library, collapse = "+"),
                format(cand$truncation)))
  invisible(x)
}


#' Run the Outcome-Blind Stress Test
#'
#' The one plasmode loop: generates simulated outcomes under the lock's
#' locked generator mode (see [create_analysis_lock()]), fits every
#' candidate on every replicate, and reports bias, RMSE, coverage, and
#' SE calibration per candidate. With `threats = NULL` this is the
#' clean-data baseline used for candidate selection; with a threat list
#' or preset name it adds the data-quality sweep, and when the lock
#' declares `dq_thresholds` the result carries the locked verdict and
#' tipping points ([dq_locked_verdict()], [dq_tipping_points()]).
#' `summary()` on the result is the degradation table.
#'
#' @param lock A `cleanroom_lock`.
#' @param candidates A `ct_candidates` set from [define_candidates()]
#'   (or a plain list of candidate specifications).
#' @param threats `NULL` (baseline only), a preset name
#'   (`"regulatory_standard"`, `"exploratory"`, `"stress"`), or a named
#'   scenario list (see the engine documentation,
#'   `?run_plasmode_dq_stress`).
#' @param effect_sizes Numeric vector of simulated risk differences.
#' @param reps Replicates per scenario cell.
#' @param design `"generate_treatment"` (default) or the deprecated
#'   `"sample_treatment"`.
#' @param pilot_q0 Baseline probabilities for `dgp_mode =
#'   "external_pilot"` locks.
#' @param q0_library Optional SuperLearner library for the plasmode
#'   baseline outcome surface (hybrid mode); default is a
#'   covariate-only logistic GLM.
#' @param max_fit_seconds Wall-clock bound per replicate's candidate
#'   fits; a timed-out replicate is recorded as failed rather than
#'   hanging the study. Default `Inf` (in-process fits).
#' @param parallel Logical; when `TRUE` and \pkg{furrr} is installed,
#'   scenario cells run under the caller's `future::plan()`.
#' @param verbose Print progress.
#' @return A `ct_stress` result (the stress-test metrics object; see
#'   the engine documentation for fields).
#' @examples
#' \dontrun{
#' cands <- define_candidates(grid = list(truncations = c(0.01, 0.05),
#'                                        libraries = list(glm = "SL.glm")))
#' st <- stress_test(lock, cands, threats = "exploratory", reps = 20)
#' summary(st)
#' }
#' @export
stress_test <- function(lock, candidates = NULL, threats = NULL,
                        effect_sizes = c(0.05), reps = 50L,
                        design = c("generate_treatment",
                                   "sample_treatment"),
                        pilot_q0 = NULL,
                        q0_library = NULL,
                        max_fit_seconds = Inf,
                        parallel = FALSE,
                        verbose = TRUE) {
  design <- match.arg(design)
  out <- run_plasmode_dq_stress(
    lock,
    tmle_candidates        = candidates,
    effect_sizes           = effect_sizes,
    reps                   = reps,
    data_quality_scenarios = if (is.null(threats)) list() else threats,
    q0_library             = q0_library,
    design                 = design,
    pilot_q0               = pilot_q0,
    fit_timeout            = max_fit_seconds,
    parallel               = parallel,
    verbose                = verbose)
  class(out) <- c("ct_stress", class(out))
  out
}


#' Select the Candidate Under the Locked Rule
#'
#' Applies the prespecified selection rule to a stress-test result:
#' `"min_rmse"` (lowest RMSE on the clean-data baseline),
#' `"min_max_rmse"` (lowest worst-case RMSE across the declared
#' threats), or `"fiord_two_stage"` (the FIORD point-estimator stage:
#' screen candidates on oracle coverage, then take the lowest-SE
#' survivor).
#'
#' @param results A [stress_test()] result (baseline rows are used for
#'   the baseline rules; degraded rows for `"min_max_rmse"`). A legacy
#'   `plasmode_results` object is also accepted.
#' @param rule One of `"min_rmse"`, `"min_max_rmse"`,
#'   `"fiord_two_stage"`.
#' @param ... Passed to the internal selector (for example
#'   `fiord_target_coverage`, `fiord_coverage_tol`, `thresholds`).
#' @return The selected candidate specification (class
#'   `tmle_selected_spec`), carrying the rule and the metrics at
#'   selection; bind it to the lock with
#'   `declare_estimand_ladder(candidate = )`.
#' @examples
#' \dontrun{
#' best <- select_candidate(st, rule = "min_max_rmse")
#' lock <- declare_estimand_ladder(lock, primary = "ATE",
#'                                 candidate = best)
#' }
#' @export
select_candidate <- function(results,
                             rule = c("min_rmse", "min_max_rmse",
                                      "fiord_two_stage"),
                             ...) {
  rule <- match.arg(rule)
  if (inherits(results, "plasmode_dq_results")) {
    base <- results$metrics[results$metrics$scenario == "none", ,
                            drop = FALSE]
    if (!nrow(base))
      stop("select_candidate: the stress-test result has no baseline ",
           "(scenario == 'none') rows.", call. = FALSE)
    if (is.null(results$tmle_candidates))
      stop("select_candidate: the result carries no candidate ",
           "specifications.", call. = FALSE)
    shim <- structure(list(metrics = base,
                           tmle_candidates = results$tmle_candidates,
                           lock_hash = results$lock_hash),
                      class = "plasmode_results")
    select_tmle_candidate(shim, rule = rule,
                          dq_results = if (rule == "min_max_rmse")
                            results else NULL, ...)
  } else if (inherits(results, "plasmode_results")) {
    select_tmle_candidate(results, rule = rule, ...)
  } else {
    stop("`results` must come from stress_test().", call. = FALSE)
  }
}


#' Negative Controls Along the Restriction Ladder
#'
#' The exported verb for the negative-control ladder: fits every
#' registered control on every nested cohort, names what turns null
#' after a restriction, and, when the lock declares `nc_criteria`,
#' grades each rung with [nc_ladder_verdict()]. See the engine
#' documentation (`?run_negative_control_ladder`) for the full
#' argument reference.
#'
#' @inheritParams run_negative_control_ladder
#' @return An object of class `ct_nc_ladder` (see
#'   [run_negative_control_ladder()]).
#' @examples
#' \dontrun{
#' ncl <- negative_control_ladder(lock,
#'   restrictions = list("transfers excluded" = !dat$transfer))
#' print(ncl)
#' }
#' @export
negative_control_ladder <- function(lock, restrictions,
                                    negative_controls = NULL,
                                    method = c("tmle", "unadjusted"),
                                    ps_method = "glm",
                                    alpha = 0.05,
                                    min_events = 5L,
                                    verbose = TRUE) {
  method <- match.arg(method)
  out <- run_negative_control_ladder(
    lock, restrictions, negative_controls = negative_controls,
    method = method, ps_method = ps_method, alpha = alpha,
    min_events = min_events, verbose = verbose)
  class(out) <- c("ct_nc_ladder", class(out))
  out
}


# ── as.data.frame methods for the ct_* results ──────────────────────────
# print/summary/plot dispatch falls through to the engine classes
# (plasmode_dq_results, nc_ladder); these give each result a tidy
# data.frame view for reporting pipelines.

#' @export
as.data.frame.ct_stress <- function(x, ...) x$metrics

#' @export
as.data.frame.ct_nc_ladder <- function(x, ...) x$table

#' @export
as.data.frame.ct_candidates <- function(x, ...) {
  data.frame(
    candidate_id    = vapply(x, function(k) k$candidate_id, character(1)),
    g_library       = vapply(x, function(k)
                        paste(k$g_library, collapse = "+"), character(1)),
    q_library       = vapply(x, function(k)
                        paste(k$q_library, collapse = "+"), character(1)),
    truncation      = vapply(x, function(k)
                        as.numeric(k$truncation[1]), numeric(1)),
    estimator       = vapply(x, function(k) k$estimator, character(1)),
    cv_scheme       = vapply(x, function(k) k$cv_scheme, character(1)),
    variance_method = vapply(x, function(k) k$variance_method,
                             character(1)),
    row.names       = NULL, stringsAsFactors = FALSE)
}