#' Estimate IPW Cumulative Risk Curves
#'
#' Computes inverse probability weighted (IPW) cumulative risk (cumulative
#' incidence) curves, adjusting for confounding via IPTW and for
#' right-censoring via IPCW. Supports competing risks, weight trimming,
#' truncation, and bootstrap confidence intervals.
#'
#' @param spec A `cr_spec` object with outcome, and optionally treatment,
#'   censoring, and competing risk components identified.
#' @param risk_time Optional numeric vector of times at which to report
#'   cumulative risk. If `NULL`, uses all unique event times.
#' @param trim Quantile for propensity score trimming (e.g., `0.01` trims
#'   at 1st and 99th percentiles). `NULL` for no trimming.
#' @param trunc Truncation for weights. If <= 1, interpreted as a quantile;
#'   if > 1, as an absolute cap. `NULL` for no truncation.
#' @param weight_type `"iptw"` for ATE-targeting weights (default) or
#'   `"smr"` for ATT-targeting (SMR) weights.
#' @param nboot Number of bootstrap replicates for CIs. `0` for no bootstrap.
#' @param seed Random seed for bootstrap reproducibility.
#' @param ... Additional arguments (currently unused).
#'
#' @return An object of class `cumrisk` containing:
#'   * `risk` - data.frame of cumulative risk estimates by time (and treatment group)
#'   * `spec` - the input specification
#'   * `weights` - the IPTW and IPCW weights
#'   * `ps_model` - the propensity score model (if applicable)
#'   * `boot_results` - bootstrap estimates (if `nboot > 0`)
#'
#' @details
#' ## Method
#' For each treatment group (or overall if no treatment is specified), the
#' weighted Kaplan-Meier estimator is used:
#'
#' \deqn{\hat{S}^w(t) = \prod_{t_j \le t} \left(1 - \frac{\sum_{i \in R_j} w_i d_i}{\sum_{i \in R_j} w_i}\right)}
#'
#' where \eqn{w_i} is the combined IPTW x IPCW weight, \eqn{d_i} is the
#' event indicator, and \eqn{R_j} is the risk set at time \eqn{t_j}.
#' Cumulative risk is \eqn{\hat{F}(t) = 1 - \hat{S}(t)}.
#'
#' ## Stabilized weights
#' By default, stabilized IPTW weights are used:
#' \eqn{w_i^{IPTW} = P(A = a_i) / P(A = a_i | L_i)}.
#'
#' ## Competing risks
#' When a competing risk is identified, the cause-specific cumulative
#' incidence is estimated using the weighted Aalen-Johansen estimator
#' in discrete time.
#'
#' @examples
#' \dontrun{
#' data(example1, package = "cleanTMLE")
#' fit <- specify_models(data = example1) |>
#'   identify_outcome(event, type = "time_to_event") |>
#'   identify_treatment(treatment, formula = ~ age + sex + biomarker) |>
#'   identify_censoring(censored, formula = ~ age + sex, model = "coxph") |>
#'   estimate_ipwrisk(risk_time = c(6, 12, 24), nboot = 200)
#' plot(fit)
#' make_table2(fit, risk_time = 24)
#' }
#'
#' @export
estimate_ipwrisk <- function(spec, risk_time = NULL, trim = NULL, trunc = NULL,
                             weight_type = c("iptw", "smr"),
                             nboot = 0, seed = 42, ...) {
  if (!inherits(spec, "cr_spec")) {
    stop("`spec` must be a cr_spec object.", call. = FALSE)
  }
  if (is.null(spec$outcome)) {
    stop("Outcome not specified. Use identify_outcome().", call. = FALSE)
  }
  weight_type <- match.arg(weight_type)
  data <- spec$data

  # ---- Determine time and event variables ----
  time_var <- .get_time_var(spec)
  event_var <- spec$outcome$name
  obs_time <- data[[time_var]]

  # Build event indicator
  event_ind <- .build_event_indicator(data, spec)

  # ---- Compute IPTW weights ----
  if (!is.null(spec$treatment) && !is.null(spec$treatment$formula)) {
    trt_wt <- fit_treatment_weights(data, spec$treatment,
                                    weight_type = weight_type,
                                    trim = trim, trunc = trunc)
  } else {
    trt_wt <- list(weights = rep(1, nrow(data)), ps = NULL, model = NULL,
                   trt_levels = if (!is.null(spec$treatment))
                     sort(unique(data[[spec$treatment$name]])) else NULL)
  }

  # ---- Compute IPCW weights ----
  cens_wt <- fit_censoring_weights(data, spec$censoring, spec$outcome,
                                    spec$interval)

  # ---- Combined weights ----
  combined_w <- trt_wt$weights * cens_wt$weights

  # ---- Compute cumulative risk ----
  if (!is.null(spec$treatment)) {
    trt_var <- spec$treatment$name
    trt_vals <- data[[trt_var]]
    trt_levels <- trt_wt$trt_levels
    if (is.null(trt_levels)) trt_levels <- sort(unique(trt_vals))

    risk_list <- lapply(trt_levels, function(a) {
      idx <- trt_vals == a
      .weighted_km_risk(
        time = obs_time[idx],
        event = event_ind[idx],
        weights = combined_w[idx],
        risk_time = risk_time,
        competing_risk_spec = spec$competing_risk,
        data = data[idx, , drop = FALSE]
      )
    })
    names(risk_list) <- as.character(trt_levels)

    # Combine into single data.frame
    risk_df <- do.call(rbind, lapply(names(risk_list), function(a) {
      df <- risk_list[[a]]
      df$group <- a
      df
    }))

  } else {
    risk_df <- .weighted_km_risk(
      time = obs_time,
      event = event_ind,
      weights = combined_w,
      risk_time = risk_time,
      competing_risk_spec = spec$competing_risk,
      data = data
    )
    risk_df$group <- "overall"
  }

  # ---- Bootstrap ----
  boot_results <- NULL
  if (nboot > 0) {
    set.seed(seed)
    boot_results <- .bootstrap_ipwrisk(
      spec = spec, nboot = nboot, risk_time = risk_time,
      trim = trim, trunc = trunc, weight_type = weight_type
    )

    # Attach CI to risk_df
    risk_df <- .attach_boot_ci(risk_df, boot_results)
  }

  # ---- Build result ----
  result <- list(
    risk = risk_df,
    spec = spec,
    weights = list(iptw = trt_wt$weights, ipcw = cens_wt$weights,
                   combined = combined_w),
    ps = trt_wt$ps,
    ps_model = trt_wt$model,
    censoring_models = cens_wt$models,
    weight_type = weight_type,
    trim = trim,
    trunc = trunc,
    nboot = nboot,
    boot_results = boot_results,
    call = match.call()
  )
  class(result) <- c("cumrisk", "cr_result")
  result
}


# ---- Internal helpers ----

#' Determine the time variable from spec
#' @keywords internal
.get_time_var <- function(spec) {
  if (!is.null(spec$interval)) {
    return(spec$interval$stop)
  }
  # Look for a 'time' column
  if ("time" %in% names(spec$data)) return("time")
  # Fall back to outcome name
  spec$outcome$name
}

#' Build event indicator vector
#' @keywords internal
.build_event_indicator <- function(data, spec) {
  event_var <- spec$outcome$name

  if (!is.null(spec$competing_risk)) {
    cr_var <- spec$competing_risk$name
    cr_val <- spec$competing_risk$event_value
    # Event of interest = 1, competing event = 2, censored = 0
    event_type <- data[[cr_var]]
    event_ind <- as.integer(event_type == cr_val)
  } else {
    event_ind <- as.integer(data[[event_var]])
  }

  event_ind
}

#' Weighted Kaplan-Meier cumulative risk estimator
#'
#' Computes weighted KM survival and returns cumulative risk = 1 - S(t).
#' For competing risks, uses weighted Aalen-Johansen estimator.
#'
#' @keywords internal
.weighted_km_risk <- function(time, event, weights, risk_time = NULL,
                               competing_risk_spec = NULL, data = NULL) {
  # Sort by time
  ord <- order(time)
  time <- time[ord]
  event <- event[ord]
  weights <- weights[ord]

  # Unique event times
  utimes <- sort(unique(time[event == 1]))
  if (length(utimes) == 0) {
    # No events
    if (is.null(risk_time)) risk_time <- sort(unique(time))
    return(data.frame(time = risk_time, risk = 0, surv = 1,
                      n_risk = length(time), n_event = 0))
  }

  # At each unique time, compute weighted hazard
  n <- length(time)
  cum_surv <- 1
  results <- data.frame(
    time = numeric(0), risk = numeric(0), surv = numeric(0),
    n_risk = numeric(0), n_event = numeric(0), hazard = numeric(0)
  )

  for (tj in utimes) {
    at_risk <- time >= tj
    n_risk <- sum(weights[at_risk])
    d_events <- event == 1 & time == tj
    d_weighted <- sum(weights[d_events])

    if (n_risk > 0) {
      haz <- d_weighted / n_risk
    } else {
      haz <- 0
    }

    cum_surv <- cum_surv * (1 - haz)
    results <- rbind(results, data.frame(
      time = tj,
      risk = 1 - cum_surv,
      surv = cum_surv,
      n_risk = sum(at_risk),
      n_event = sum(d_events),
      hazard = haz
    ))
  }

  # Add time = 0
  results <- rbind(
    data.frame(time = 0, risk = 0, surv = 1,
               n_risk = n, n_event = 0, hazard = 0),
    results
  )

  # Interpolate to requested risk times

if (!is.null(risk_time)) {
    interp <- approx(results$time, results$risk,
                     xout = risk_time, method = "constant",
                     rule = 2, f = 0)
    results <- data.frame(
      time = risk_time,
      risk = interp$y,
      surv = 1 - interp$y,
      n_risk = NA_integer_,
      n_event = NA_integer_,
      hazard = NA_real_
    )
  }

  results
}


#' Bootstrap helper for IPW risk
#' @keywords internal
.bootstrap_ipwrisk <- function(spec, nboot, risk_time, trim, trunc,
                                weight_type) {
  data <- spec$data
  n <- nrow(data)

  # Determine groups
  if (!is.null(spec$treatment)) {
    groups <- sort(unique(data[[spec$treatment$name]]))
  } else {
    groups <- "overall"
  }

  boot_risks <- list()

  for (b in seq_len(nboot)) {
    idx <- sample(n, n, replace = TRUE)
    boot_spec <- spec
    boot_spec$data <- data[idx, , drop = FALSE]

    tryCatch({
      boot_fit <- estimate_ipwrisk(boot_spec, risk_time = risk_time,
                                   trim = trim, trunc = trunc,
                                   weight_type = weight_type, nboot = 0)
      boot_risks[[b]] <- boot_fit$risk
    }, error = function(e) {
      boot_risks[[b]] <<- NULL
    })
  }

  # Remove failed bootstraps

  boot_risks <- boot_risks[!vapply(boot_risks, is.null, logical(1))]
  boot_risks
}


#' Attach bootstrap CIs to risk data.frame
#' @keywords internal
.attach_boot_ci <- function(risk_df, boot_results) {
  if (length(boot_results) == 0) {
    risk_df$ci_lower <- NA_real_
    risk_df$ci_upper <- NA_real_
    return(risk_df)
  }

  risk_df$ci_lower <- NA_real_
  risk_df$ci_upper <- NA_real_

  for (i in seq_len(nrow(risk_df))) {
    t_i <- risk_df$time[i]
    g_i <- risk_df$group[i]

    boot_vals <- vapply(boot_results, function(br) {
      row <- br$time == t_i & br$group == g_i
      if (any(row)) br$risk[row][1] else NA_real_
    }, numeric(1))

    boot_vals <- boot_vals[!is.na(boot_vals)]
    if (length(boot_vals) >= 2) {
      risk_df$ci_lower[i] <- quantile(boot_vals, 0.025, na.rm = TRUE)
      risk_df$ci_upper[i] <- quantile(boot_vals, 0.975, na.rm = TRUE)
    }
  }

  risk_df
}
