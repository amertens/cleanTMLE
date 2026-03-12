#' Update and Re-estimate Model Components
#'
#' These functions modify an existing specification or fitted object's
#' treatment, outcome, or censoring specification and optionally re-fit.
#'
#' @name update_methods
#' @param x A `cr_spec` or fitted result object.
#' @param ... Arguments passed to the corresponding `identify_*` function.
NULL

#' @rdname update_methods
#' @export
update_treatment <- function(x, ...) {
  spec <- if (inherits(x, "cr_spec")) x else x$spec
  spec <- identify_treatment(spec, ...)
  if (inherits(x, "cr_spec")) {
    return(spec)
  }
  # Return updated spec; user should call re_estimate()
  x$spec <- spec
  x
}

#' @rdname update_methods
#' @export
update_outcome <- function(x, ...) {
  spec <- if (inherits(x, "cr_spec")) x else x$spec
  spec <- identify_outcome(spec, ...)
  if (inherits(x, "cr_spec")) {
    return(spec)
  }
  x$spec <- spec
  x
}

#' @rdname update_methods
#' @export
update_censoring <- function(x, ...) {
  spec <- if (inherits(x, "cr_spec")) x else x$spec
  spec$censoring <- list()  # Reset censoring

  spec <- identify_censoring(spec, ...)
  if (inherits(x, "cr_spec")) {
    return(spec)
  }
  x$spec <- spec
  x
}

#' Re-estimate a Fitted Model
#'
#' Re-runs the estimation function that produced the fitted object,
#' using the (potentially updated) specification.
#'
#' @param x A fitted result object whose spec may have been updated.
#' @param ... Additional arguments passed to the estimation function.
#'
#' @return A new fitted result object.
#'
#' @export
re_estimate <- function(x, ...) {
  if (inherits(x, "cumrisk")) {
    return(estimate_ipwrisk(x$spec, ...))
  } else if (inherits(x, "gcomp")) {
    return(estimate_gcomprisk(x$spec, ...))
  } else if (inherits(x, "aipw")) {
    return(estimate_aipwrisk(x$spec, ...))
  } else if (inherits(x, "hr")) {
    return(estimate_ipwhr(x$spec, ...))
  } else {
    stop("Unknown result class for re-estimation.", call. = FALSE)
  }
}


#' Compare Multiple Fitted Objects
#'
#' Compares point estimates across multiple fitted result objects,
#' producing a combined data.frame for easy comparison.
#'
#' @param ... Named fitted result objects to compare.
#' @param risk_time Optional time at which to compare cumulative risk.
#'
#' @return A data.frame comparing estimates across fits.
#'
#' @export
compare_fits <- function(..., risk_time = NULL) {
  fits <- list(...)
  if (is.null(names(fits))) {
    names(fits) <- paste0("fit", seq_along(fits))
  }

  rows <- list()
  for (nm in names(fits)) {
    fit <- fits[[nm]]
    if (inherits(fit, c("cumrisk", "gcomp", "aipw"))) {
      risk_df <- fit$risk
      if (!is.null(risk_time)) {
        risk_df <- risk_df[risk_df$time == risk_time, ]
      }
      has_ci <- "ci_lower" %in% names(risk_df)
      for (i in seq_len(nrow(risk_df))) {
        rows <- c(rows, list(data.frame(
          fit_name  = nm,
          estimator = class(fit)[1],
          group     = risk_df$group[i],
          time      = risk_df$time[i],
          risk      = risk_df$risk[i],
          ci_lower  = if (has_ci) risk_df$ci_lower[i] else NA_real_,
          ci_upper  = if (has_ci) risk_df$ci_upper[i] else NA_real_,
          stringsAsFactors = FALSE
        )))
      }
    } else if (inherits(fit, "hr")) {
      ht <- fit$hr_table
      for (i in seq_len(nrow(ht))) {
        rows <- c(rows, list(data.frame(
          fit_name = nm,
          estimator = "hr",
          group = ht$term[i],
          time = NA_real_,
          risk = ht$hr[i],
          ci_lower = ht$ci_lower[i],
          ci_upper = ht$ci_upper[i],
          stringsAsFactors = FALSE
        )))
      }
    }
  }

  result <- do.call(rbind, rows)
  rownames(result) <- NULL
  result
}
