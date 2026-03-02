#' Identify Model Components
#'
#' These functions add component specifications to a `cr_spec` object created
#' by [specify_models()]. Each function declares one aspect of the causal model:
#' outcome, treatment, censoring, competing risks, subject identifier, time
#' interval, or missingness.
#'
#' @name identify
#' @param spec A `cr_spec` object from [specify_models()].
#' @param name The column name (bare or quoted) identifying the variable.
#' @param formula A one-sided formula for the model (e.g., `~ age + sex`).
#' @param type For outcomes: one of `"time_to_event"`, `"binary"`, or
#'   `"continuous"`.
#' @param model For censoring: `"coxph"` or `"glm"`. For treatment: `"glm"`
#'   or `"multinom"`.
#' @param event_value For competing risks: the value in the event-type column
#'   indicating the event of interest.
#' @param start,stop Column names for interval start and stop times.
#' @param ... Additional arguments passed to model fitting functions.
#'
#' @return The modified `cr_spec` object (for piping).
NULL

#' @rdname identify
#' @export
identify_outcome <- function(spec, name, formula = NULL,
                             type = c("time_to_event", "binary", "continuous"),
                             ...) {
  if (!inherits(spec, "cr_spec")) {
    stop("`spec` must be a cr_spec object.", call. = FALSE)
  }
  type <- match.arg(type)
  name_quo <- rlang::enquo(name)
  name_str <- tryCatch(rlang::quo_name(name_quo), error = function(e) {
    as.character(substitute(name))
  })

  spec$outcome <- list(
    name = name_str,
    formula = formula,
    type = type,
    extra = list(...)
  )
  spec
}

#' @rdname identify
#' @export
identify_censoring <- function(spec, name, formula = NULL,
                               model = c("coxph", "glm"), ...) {
  if (!inherits(spec, "cr_spec")) {
    stop("`spec` must be a cr_spec object.", call. = FALSE)
  }
  model <- match.arg(model)
  name_quo <- rlang::enquo(name)
  name_str <- tryCatch(rlang::quo_name(name_quo), error = function(e) {
    as.character(substitute(name))
  })

  cens_entry <- list(
    name = name_str,
    formula = formula,
    model = model,
    extra = list(...)
  )
  spec$censoring <- c(spec$censoring, list(cens_entry))
  spec
}

#' @rdname identify
#' @export
identify_treatment <- function(spec, name, formula = NULL,
                               model = c("glm", "multinom"), ...) {
  if (!inherits(spec, "cr_spec")) {
    stop("`spec` must be a cr_spec object.", call. = FALSE)
  }
  model <- match.arg(model)
  name_quo <- rlang::enquo(name)
  name_str <- tryCatch(rlang::quo_name(name_quo), error = function(e) {
    as.character(substitute(name))
  })

  spec$treatment <- list(
    name = name_str,
    formula = formula,
    model = model,
    extra = list(...)
  )
  spec
}

#' @rdname identify
#' @export
identify_competing_risk <- function(spec, name, event_value) {

  if (!inherits(spec, "cr_spec")) {
    stop("`spec` must be a cr_spec object.", call. = FALSE)
  }
  name_quo <- rlang::enquo(name)
  name_str <- tryCatch(rlang::quo_name(name_quo), error = function(e) {
    as.character(substitute(name))
  })

  spec$competing_risk <- list(
    name = name_str,
    event_value = event_value
  )
  spec
}

#' @rdname identify
#' @export
identify_subject <- function(spec, name) {
  if (!inherits(spec, "cr_spec")) {
    stop("`spec` must be a cr_spec object.", call. = FALSE)
  }
  name_quo <- rlang::enquo(name)
  name_str <- tryCatch(rlang::quo_name(name_quo), error = function(e) {
    as.character(substitute(name))
  })

  spec$subject <- list(name = name_str)
  spec
}

#' @rdname identify
#' @export
identify_interval <- function(spec, start, stop) {
  if (!inherits(spec, "cr_spec")) {
    stop("`spec` must be a cr_spec object.", call. = FALSE)
  }
  start_quo <- rlang::enquo(start)
  stop_quo <- rlang::enquo(stop)
  start_str <- tryCatch(rlang::quo_name(start_quo), error = function(e) {
    as.character(substitute(start))
  })
  stop_str <- tryCatch(rlang::quo_name(stop_quo), error = function(e) {
    as.character(substitute(stop))
  })

  spec$interval <- list(start = start_str, stop = stop_str)
  spec
}

#' @rdname identify
#' @export
identify_missing <- function(spec, ...) {
  if (!inherits(spec, "cr_spec")) {
    stop("`spec` must be a cr_spec object.", call. = FALSE)
  }
  dots <- rlang::enquos(...)
  miss_names <- vapply(dots, rlang::quo_name, character(1))

  spec$missing <- list(
    variables = miss_names,
    enabled = TRUE
  )
  spec
}
