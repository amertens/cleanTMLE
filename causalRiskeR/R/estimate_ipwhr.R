#' IPW Hazard Ratio Estimation via Weighted Cox Model
#'
#' Fits a weighted Cox proportional hazards model using IPTW weights to
#' estimate the causal hazard ratio for the treatment effect. Reports
#' robust (sandwich) standard errors.
#'
#' @param spec A `cr_spec` object with outcome and treatment identified.
#' @param covariates Optional character vector of additional covariates
#'   to include in the outcome (Cox) model beyond treatment.
#' @param trim Quantile for propensity score trimming.
#' @param trunc Weight truncation.
#' @param weight_type `"iptw"` or `"smr"`.
#' @param robust Logical; if `TRUE` (default), use robust (sandwich) SEs.
#' @param ... Additional arguments.
#'
#' @return An object of class `hr` (inherits from `cr_result`) containing:
#'   * `hr_table` - data.frame with HR, CI, p-value for each coefficient
#'   * `cox_model` - the fitted weighted coxph object
#'   * `weights` - IPTW weights used
#'
#' @details
#' ## Method
#' 1. Fit a propensity score model for treatment via logistic regression.
#' 2. Compute stabilized IPTW weights.
#' 3. Fit a weighted Cox PH model: `coxph(Surv(time, event) ~ treatment,
#'    weights = w, robust = TRUE)`.
#' 4. Extract hazard ratios with robust standard errors via the sandwich
#'    variance estimator.
#'
#' ## Assumptions
#' * No unmeasured confounding (conditional exchangeability given L).
#' * Correctly specified propensity score model.
#' * Proportional hazards in the weighted model.
#' * Positivity: all subjects have non-zero probability of each treatment.
#'
#' @references
#' Robins JM, Hernan MA, Brumback B (2000). Marginal structural models
#' and causal inference in epidemiology. Epidemiology, 11(5):550-560.
#'
#' @examples
#' \dontrun{
#' hr_fit <- specify_models(data = example1) |>
#'   identify_outcome(event, type = "time_to_event") |>
#'   identify_treatment(treatment, formula = ~ age + sex + biomarker) |>
#'   estimate_ipwhr()
#' print(hr_fit)
#' forest_plot(hr_fit)
#' }
#'
#' @export
estimate_ipwhr <- function(spec, covariates = NULL, trim = NULL, trunc = NULL,
                           weight_type = c("iptw", "smr"),
                           robust = TRUE, ...) {
  if (!inherits(spec, "cr_spec")) {
    stop("`spec` must be a cr_spec object.", call. = FALSE)
  }
  if (is.null(spec$outcome)) {
    stop("Outcome not specified.", call. = FALSE)
  }
  if (is.null(spec$treatment)) {
    stop("Treatment not specified.", call. = FALSE)
  }
  weight_type <- match.arg(weight_type)

  data <- spec$data
  time_var <- .get_time_var(spec)
  event_var <- spec$outcome$name
  event_ind <- .build_event_indicator(data, spec)
  trt_var <- spec$treatment$name

  # Fit IPTW
  trt_wt <- fit_treatment_weights(data, spec$treatment,
                                  weight_type = weight_type,
                                  trim = trim, trunc = trunc)

  # Build Cox formula
  rhs_vars <- trt_var
  if (!is.null(covariates)) {
    rhs_vars <- c(rhs_vars, covariates)
  }
  cox_fml <- as.formula(paste(
    "survival::Surv(", time_var, ", event_ind) ~",
    paste(rhs_vars, collapse = " + ")
  ))

  # Fit weighted Cox model
  data_fit <- data
  data_fit$event_ind <- event_ind
  data_fit$.weights <- trt_wt$weights

  cox_fit <- survival::coxph(cox_fml, data = data_fit,
                             weights = .weights,
                             robust = robust)

  # Extract HR table
  coef_tbl <- summary(cox_fit)$coefficients
  if (robust && "robust se" %in% colnames(coef_tbl)) {
    se_col <- "robust se"
    p_col <- if ("Pr(>|z|)" %in% colnames(coef_tbl)) "Pr(>|z|)" else ncol(coef_tbl)
  } else {
    se_col <- "se(coef)"
    p_col <- if ("Pr(>|z|)" %in% colnames(coef_tbl)) "Pr(>|z|)" else ncol(coef_tbl)
  }

  hr_table <- data.frame(
    term = rownames(coef_tbl),
    log_hr = coef_tbl[, "coef"],
    hr = exp(coef_tbl[, "coef"]),
    se = coef_tbl[, se_col],
    ci_lower = exp(coef_tbl[, "coef"] - 1.96 * coef_tbl[, se_col]),
    ci_upper = exp(coef_tbl[, "coef"] + 1.96 * coef_tbl[, se_col]),
    p_value = coef_tbl[, p_col],
    row.names = NULL,
    stringsAsFactors = FALSE
  )

  result <- list(
    hr_table = hr_table,
    cox_model = cox_fit,
    weights = list(iptw = trt_wt$weights),
    ps = trt_wt$ps,
    ps_model = trt_wt$model,
    spec = spec,
    weight_type = weight_type,
    trim = trim,
    trunc = trunc,
    robust = robust,
    call = match.call()
  )
  class(result) <- c("hr", "cr_result")
  result
}


#' Extract HR Data from an hr Object
#'
#' Returns the data.frame backing the hazard ratio results, suitable for
#' forest plots or custom tables.
#'
#' @param x An object of class `hr`.
#' @param ... Additional arguments (unused).
#'
#' @return A data.frame with columns: term, log_hr, hr, se, ci_lower,
#'   ci_upper, p_value.
#'
#' @export
hr_data <- function(x, ...) {
  if (!inherits(x, "hr")) {
    stop("`x` must be an hr object.", call. = FALSE)
  }
  x$hr_table
}


#' @export
print.hr <- function(x, ...) {
  cat("IPW Hazard Ratio Estimates (Weighted Cox Model)\n")
  cat("================================================\n")
  cat("Weight type:", x$weight_type, "\n")
  if (!is.null(x$trim)) cat("Trim:", x$trim, "\n")
  if (!is.null(x$trunc)) cat("Trunc:", x$trunc, "\n")
  cat("Robust SE:", x$robust, "\n\n")

  ht <- x$hr_table
  for (i in seq_len(nrow(ht))) {
    cat(sprintf("  %s: HR = %.3f (95%% CI: %.3f - %.3f), p = %.4f\n",
                ht$term[i], ht$hr[i], ht$ci_lower[i], ht$ci_upper[i],
                ht$p_value[i]))
  }
  invisible(x)
}


#' @export
plot.hr <- function(x, ...) {
  forest_plot(x, ...)
}
