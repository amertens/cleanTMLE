.estimand_names <- c("ATE", "trimmed_ATE", "ATT", "ATO")
.learner_names <- c("glm", "glmnet", "earth", "xgboost", "nnet")
.default_library <- c("glmnet", "earth", "xgboost", "nnet")

.expand_library <- function(x, arg) {
  if (!is.character(x) || !length(x))
    stop("`", arg, "` must be a character vector of learner names.", call. = FALSE)
  x <- unique(unlist(lapply(x, function(l)
    if (identical(l, "default")) .default_library else l)))
  bad <- setdiff(x, .learner_names)
  if (length(bad))
    stop("Unknown learner(s) in `", arg, "`: ", paste(bad, collapse = ", "),
         ". Use ", paste(c(.learner_names, "default"), collapse = ", "), ".",
         call. = FALSE)
  x
}

#' Declare the analysis plan
#'
#' Holds every declaration of the analysis. The lock fingerprints the whole
#' plan, so nothing can be added after [create_analysis_lock()].
#'
#' @param outcome Name of a 0/1 outcome column, or `c(time = , event = )`
#'   for a time-to-event outcome (event 0 = censored, 1 = the event).
#' @param estimands Ordered ladder, from `"ATE"`, `"trimmed_ATE"`, `"ATT"`,
#'   `"ATO"`. The first feasible one becomes the primary.
#' @param ps_library Learners for the propensity score, used by every
#'   candidate: any of `"glm"`, `"glmnet"`, `"earth"`, `"xgboost"`,
#'   `"nnet"`, or `"default"` (the last four).
#' @param candidates `list(library = <named list of learner vectors>,
#'   truncation = <numbers in (0, 0.5)>)`; candidates are their crossing.
#' @param tolerance `list(bias = , coverage = )`; bias is absolute bias on
#'   the risk-difference scale.
#' @param surfaces Outcome-surface family: `baseline_risk`, `log_or`,
#'   `heterogeneity`, `forms` (`"linear"`, `"nonlinear"`), optional `q0`
#'   (steward-supplied risks, one per design row) and `custom` (named list
#'   of `function(a, X)` returning risks).
#' @param trim_band `c(lower, upper)` propensity band for `"trimmed_ATE"`.
#' @param target_time,hazards Time-to-event only: the target time and the
#'   hazard models passed to `concrete` as its `Model` entries.
#' @param negative_controls Named character vector, column = domain.
#' @param nc_criteria `list(null_band = , rule = "point" or "ci",
#'   min_per_domain = )`.
#' @param restrictions Named list of one-sided formulas, least to most
#'   restricted.
#' @param select Optional function of the metrics table returning a
#'   candidate id; replaces the minimum worst-surface RMSE rule.
#' @param K Outer cross-fitting folds (1 = no cross-fitting).
#' @param V Inner folds of the super learner.
#' @param reps Plasmode repetitions.
#' @param max_reps Cap on repetitions for cells whose result is uncertain;
#'   see [simulate_design()].
#' @param seed Seed for every random step.
#' @return A `cr_plan`.
#' @export
analysis_plan <- function(outcome, estimands,
                          ps_library = "default",
                          candidates = list(library = list(glm = "glm", sl = "default"),
                                            truncation = c(0.01, 0.05)),
                          tolerance = list(bias = 0.02, coverage = 0.90),
                          surfaces = list(),
                          trim_band = NULL,
                          target_time = NULL,
                          hazards = NULL,
                          negative_controls = NULL,
                          nc_criteria = NULL,
                          restrictions = NULL,
                          select = NULL,
                          K = 5L, V = 5L, reps = 200L, max_reps = 1000L, seed = 1L) {
  # Anything but a clean FALSE stops, so an NA condition is an error, not a pass.
  stop_if <- function(cond, ...) if (!isFALSE(cond)) stop(..., call. = FALSE)
  whole <- function(x) is.numeric(x) && length(x) == 1L && !is.na(x) &&
    is.finite(x) && x == round(x) && abs(x) <= .Machine$integer.max
  for (arg in c("K", "V", "reps", "max_reps", "seed"))
    stop_if(!whole(get(arg)), "`", arg, "` must be a single whole number.")

  stop_if(!is.character(outcome) || !length(outcome) %in% 1:2,
          "`outcome` must be a column name, or c(time = , event = ).")
  if (length(outcome) == 2L) {
    stop_if(!setequal(names(outcome), c("time", "event")),
            "A time-to-event `outcome` must be named c(time = , event = ).")
    outcome_type <- "tte"
    stop_if(!is.numeric(target_time) || length(target_time) != 1L || target_time <= 0,
            "A time-to-event plan needs a positive `target_time`.")
    stop_if(!is.list(hazards) || !length(hazards),
            "A time-to-event plan needs `hazards`, the hazard models passed to concrete.")
  } else {
    outcome_type <- "binary"
  }

  stop_if(!is.character(estimands) || !length(estimands) || anyDuplicated(estimands) ||
            !all(estimands %in% .estimand_names),
          "`estimands` must be distinct values from: ",
          paste(.estimand_names, collapse = ", "), ".")
  stop_if(outcome_type == "tte" && any(estimands %in% c("ATT", "ATO")),
          "ATT and ATO are not available for a time-to-event outcome.")
  if ("trimmed_ATE" %in% estimands)
    stop_if(!is.numeric(trim_band) || length(trim_band) != 2L ||
              !(0 < trim_band[1] && trim_band[1] < trim_band[2] && trim_band[2] < 1),
            "`trim_band` must be c(lower, upper) with 0 < lower < upper < 1 ",
            "when trimmed_ATE is declared.")

  ps_library <- .expand_library(ps_library, "ps_library")
  lib <- candidates$library
  stop_if(!is.list(lib) || is.null(names(lib)) || any(!nzchar(names(lib))) ||
            anyDuplicated(names(lib)),
          "`candidates$library` must be a named list of learner vectors.")
  candidates$library <- Map(.expand_library, lib,
                            paste0("candidates$library$", names(lib)))
  tr <- candidates$truncation
  stop_if(!is.numeric(tr) || !length(tr) || any(tr <= 0 | tr >= 0.5),
          "`candidates$truncation` must be numbers in (0, 0.5).")
  candidates$truncation <- sort(unique(tr))

  tolerance <- utils::modifyList(list(bias = 0.02, coverage = 0.90), tolerance)
  stop_if(tolerance$bias <= 0 || tolerance$coverage <= 0 || tolerance$coverage >= 1,
          "`tolerance` needs bias > 0 and 0 < coverage < 1.")

  surfaces <- utils::modifyList(list(baseline_risk = 0.10, log_or = 0.5,
                                     heterogeneity = 0.5,
                                     forms = c("linear", "nonlinear")), surfaces)
  stop_if(surfaces$baseline_risk <= 0 || surfaces$baseline_risk >= 1,
          "`surfaces$baseline_risk` must lie in (0, 1).")
  stop_if(!is.numeric(surfaces$log_or) || length(surfaces$log_or) != 1L ||
            !is.finite(surfaces$log_or),
          "`surfaces$log_or` must be a single finite number.")
  stop_if(surfaces$heterogeneity < 0, "`surfaces$heterogeneity` must be >= 0.")
  stop_if(!length(surfaces$forms) || !all(surfaces$forms %in% c("linear", "nonlinear")),
          "`surfaces$forms` must be from: linear, nonlinear.")
  stop_if(!is.null(surfaces$custom) &&
            (!is.list(surfaces$custom) || is.null(names(surfaces$custom)) ||
               !all(vapply(surfaces$custom, is.function, logical(1)))),
          "`surfaces$custom` must be a named list of function(a, X).")

  if (!is.null(negative_controls)) {
    stop_if(!is.character(negative_controls) || is.null(names(negative_controls)),
            "`negative_controls` must be a named character vector: column = domain.")
    nc_criteria <- utils::modifyList(
      list(null_band = c(-0.02, 0.02), rule = "point", min_per_domain = 1L),
      nc_criteria %||% list())
    stop_if(!nc_criteria$rule %in% c("point", "ci"),
            "`nc_criteria$rule` must be 'point' or 'ci'.")
  }
  if (!is.null(restrictions))
    stop_if(!is.list(restrictions) || is.null(names(restrictions)) ||
              !all(vapply(restrictions, function(f)
                inherits(f, "formula") && length(f) == 2L, logical(1))),
            "`restrictions` must be a named list of one-sided formulas.")
  stop_if(!is.null(select) && !is.function(select),
          "`select` must be NULL or a function of the metrics table.")
  stop_if(K < 1 || V < 2 || reps < 2, "Need K >= 1, V >= 2 and reps >= 2.")
  stop_if(max_reps < reps, "Need max_reps >= reps.")

  structure(list(
    outcome = outcome, outcome_type = outcome_type, target_time = target_time,
    hazards = hazards, estimands = estimands, ps_library = ps_library,
    candidates = candidates, tolerance = tolerance, surfaces = surfaces,
    trim_band = trim_band, negative_controls = negative_controls,
    nc_criteria = nc_criteria, restrictions = restrictions, select = select,
    K = as.integer(K), V = as.integer(V), reps = as.integer(reps),
    max_reps = as.integer(max_reps), seed = as.integer(seed)), class = "cr_plan")
}

.outcome_columns <- function(plan) unname(plan$outcome)

.candidates <- function(plan) {
  g <- expand.grid(truncation = plan$candidates$truncation,
                   library = names(plan$candidates$library),
                   stringsAsFactors = FALSE)
  lapply(seq_len(nrow(g)), function(i) list(
    id = sprintf("%s_t%s", g$library[i], g$truncation[i]),
    library = g$library[i],
    learners = plan$candidates$library[[g$library[i]]],
    truncation = g$truncation[i]))
}
