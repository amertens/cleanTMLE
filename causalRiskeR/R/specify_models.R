#' Create a Model Specification Object
#'
#' Initializes a model specification container that stores information about
#' the outcome, treatment, censoring, competing risks, and other model
#' components needed for causal risk estimation.
#'
#' @param data A data.frame containing the analysis dataset.
#' @param ... Currently unused; reserved for future extensions.
#'
#' @return An object of class `cr_spec` (causalRiskeR specification).
#'
#' @details
#' The specification object is the central data structure in causalRiskeR.
#' After creating it with `specify_models()`, use the `identify_*()` family
#' of functions to declare outcome, treatment, censoring, and other
#' components. Then pass the completed specification to an estimation function
#' such as [estimate_ipwrisk()] or [estimate_ipwhr()].
#'
#' @examples
#' \dontrun{
#' spec <- specify_models(data = example1) |>
#'   identify_outcome(death, type = "time_to_event") |>
#'   identify_treatment(treatment, formula = ~ age + sex) |>
#'   identify_censoring(censored, formula = ~ age + sex)
#' }
#'
#' @export
specify_models <- function(data, ...) {
  if (!is.data.frame(data)) {
    stop("`data` must be a data.frame.", call. = FALSE)
  }

  spec <- list(
    data = data,
    outcome = NULL,
    treatment = NULL,
    censoring = list(),
    competing_risk = NULL,
    subject = NULL,
    interval = NULL,
    missing = NULL,
    fitted_models = list()
  )
  class(spec) <- "cr_spec"
  spec
}

#' @export
print.cr_spec <- function(x, ...) {
  cat("causalRiskeR Model Specification\n")
  cat("================================\n")
  cat("Data:", nrow(x$data), "observations,", ncol(x$data), "variables\n")

  if (!is.null(x$outcome)) {
    cat("Outcome:", x$outcome$name, "(", x$outcome$type, ")\n")
  } else {
    cat("Outcome: not specified\n")
  }

  if (!is.null(x$treatment)) {
    cat("Treatment:", x$treatment$name, "\n")
    if (!is.null(x$treatment$formula)) {
      cat("  PS model:", deparse(x$treatment$formula), "\n")
    }
  } else {
    cat("Treatment: not specified\n")
  }

  if (length(x$censoring) > 0) {
    cens_names <- vapply(x$censoring, function(c) c$name, character(1))
    cat("Censoring:", paste(cens_names, collapse = ", "), "\n")
  } else {
    cat("Censoring: not specified\n")
  }

  if (!is.null(x$competing_risk)) {
    cat("Competing risk:", x$competing_risk$name,
        "(event_value =", x$competing_risk$event_value, ")\n")
  }

  if (!is.null(x$subject)) {
    cat("Subject ID:", x$subject$name, "\n")
  }

  if (!is.null(x$interval)) {
    cat("Interval: [", x$interval$start, ",", x$interval$stop, "]\n")
  }

  invisible(x)
}
