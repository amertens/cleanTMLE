#' Augmented IPW (Doubly Robust) Estimator for Cumulative Risk
#'
#' Estimates cumulative risk curves using the augmented inverse probability
#' weighted (AIPW) estimator, which is doubly robust: consistent if either
#' the treatment/censoring model or the outcome model is correctly specified.
#'
#' @param spec A `cr_spec` object with outcome, treatment, and optionally
#'   censoring components identified. The outcome should have a formula.
#' @param risk_time Numeric vector of times at which to evaluate risk.
#' @param trim Quantile for propensity score trimming.
#' @param trunc Weight truncation.
#' @param nboot Number of bootstrap replicates.
#' @param seed Random seed.
#' @param ... Additional arguments.
#'
#' @return An object of class `aipw` (inherits from `cr_result`).
#'
#' @details
#' ## Method
#' The AIPW estimator combines the IPW and g-computation approaches:
#'
#' \deqn{\hat{\psi}_{AIPW}(t) = \frac{1}{n} \sum_{i=1}^{n}
#'   \left[\frac{I(A_i = a)}{P(A_i = a | L_i)} (Y_i(t) - \hat{m}(t, a, L_i))
#'   + \hat{m}(t, a, L_i)\right]}
#'
#' where \eqn{\hat{m}(t, a, L_i)} is the outcome model prediction.
#'
#' For time-to-event data, we use a discrete-time approximation:
#' at each evaluation time, define a binary outcome I(T <= t, event),
#' and apply the standard AIPW formula.
#'
#' @export
estimate_aipwrisk <- function(spec, risk_time = NULL, trim = NULL,
                              trunc = NULL, nboot = 0, seed = 42, ...) {
  if (!inherits(spec, "cr_spec")) {
    stop("`spec` must be a cr_spec object.", call. = FALSE)
  }
  if (is.null(spec$outcome)) {
    stop("Outcome not specified.", call. = FALSE)
  }
  if (is.null(spec$treatment)) {
    stop("Treatment must be specified for AIPW.", call. = FALSE)
  }

  data <- spec$data
  time_var <- .get_time_var(spec)
  event_var <- spec$outcome$name
  event_ind <- .build_event_indicator(data, spec)
  obs_time <- data[[time_var]]

  trt_var <- spec$treatment$name
  trt_vals <- data[[trt_var]]
  trt_levels <- sort(unique(trt_vals))

  # Fit PS model
  trt_wt <- fit_treatment_weights(data, spec$treatment,
                                  weight_type = "iptw",
                                  trim = trim, trunc = trunc)

  # Fit outcome model (Cox PH)
  if (!is.null(spec$outcome$formula)) {
    rhs <- deparse(spec$outcome$formula[[2]])
  } else {
    rhs <- trt_var
  }
  surv_fml <- as.formula(paste("survival::Surv(", time_var, ", event_ind) ~", rhs))
  data_fit <- data
  data_fit$event_ind <- event_ind
  cox_fit <- survival::coxph(surv_fml, data = data_fit)

  # Determine evaluation times
  all_event_times <- sort(unique(obs_time[event_ind == 1]))
  if (is.null(risk_time)) risk_time <- all_event_times

  # AIPW at each time point
  risk_list <- list()

  for (a in trt_levels) {
    a_char <- as.character(a)

    # IPW component
    is_a <- trt_vals == a
    ps_a <- if (a == trt_levels[2] && length(trt_levels) == 2) {
      trt_wt$ps
    } else if (a == trt_levels[1] && length(trt_levels) == 2) {
      1 - trt_wt$ps
    } else {
      # For multi-level, use 1/weight * marginal
      rep(mean(is_a), nrow(data))
    }
    ps_a <- pmax(ps_a, 0.01)

    # Outcome model predictions under A=a
    cf_data <- data_fit
    cf_data[[trt_var]] <- a
    sf_a <- survival::survfit(cox_fit, newdata = cf_data)

    risks_at_t <- numeric(length(risk_time))

    for (j in seq_along(risk_time)) {
      t_j <- risk_time[j]

      # Binary outcome at time t: I(T <= t & event)
      Y_t <- as.numeric(obs_time <= t_j & event_ind == 1)

      # Predicted risk at t_j for each individual under A=a
      if (is.matrix(sf_a$surv)) {
        surv_pred <- approx(c(0, sf_a$time), c(1, sf_a$surv[, 1]),
                            xout = t_j, method = "constant",
                            rule = 2, f = 0)$y
        # Use individual predictions
        m_hat <- numeric(nrow(data))
        for (ii in seq_len(nrow(data))) {
          col_surv <- if (ncol(sf_a$surv) >= ii) sf_a$surv[, ii] else sf_a$surv[, 1]
          m_hat[ii] <- 1 - approx(c(0, sf_a$time), c(1, col_surv),
                                  xout = t_j, method = "constant",
                                  rule = 2, f = 0)$y
        }
      } else {
        m_hat_val <- 1 - approx(c(0, sf_a$time), c(1, sf_a$surv),
                                xout = t_j, method = "constant",
                                rule = 2, f = 0)$y
        m_hat <- rep(m_hat_val, nrow(data))
      }

      # AIPW estimator
      aipw_vals <- (is_a / ps_a) * (Y_t - m_hat) + m_hat
      risks_at_t[j] <- mean(aipw_vals)
    }

    risk_list[[a_char]] <- data.frame(
      time = risk_time,
      risk = pmax(0, pmin(1, risks_at_t)),
      surv = 1 - pmax(0, pmin(1, risks_at_t)),
      group = a_char,
      stringsAsFactors = FALSE
    )
  }

  risk_df <- do.call(rbind, risk_list)
  rownames(risk_df) <- NULL

  # Bootstrap
  boot_results <- NULL
  if (nboot > 0) {
    set.seed(seed)
    boot_results <- .bootstrap_aipw(spec, nboot, risk_time, trim, trunc)
    risk_df <- .attach_boot_ci(risk_df, boot_results)
  }

  result <- list(
    risk = risk_df,
    spec = spec,
    ps_model = trt_wt$model,
    outcome_model = cox_fit,
    weights = list(iptw = trt_wt$weights),
    ps = trt_wt$ps,
    nboot = nboot,
    boot_results = boot_results,
    call = match.call()
  )
  class(result) <- c("aipw", "cr_result")
  result
}

#' @keywords internal
.bootstrap_aipw <- function(spec, nboot, risk_time, trim, trunc) {
  data <- spec$data
  n <- nrow(data)
  boot_risks <- list()
  for (b in seq_len(nboot)) {
    idx <- sample(n, n, replace = TRUE)
    boot_spec <- spec
    boot_spec$data <- data[idx, , drop = FALSE]
    tryCatch({
      boot_fit <- estimate_aipwrisk(boot_spec, risk_time = risk_time,
                                    trim = trim, trunc = trunc, nboot = 0)
      boot_risks[[b]] <- boot_fit$risk
    }, error = function(e) {
      boot_risks[[b]] <<- NULL
    })
  }
  boot_risks[!vapply(boot_risks, is.null, logical(1))]
}

#' @export
print.aipw <- function(x, ...) {
  cat("AIPW (Doubly Robust) Cumulative Risk Estimate\n")
  cat("==============================================\n")
  groups <- unique(x$risk$group)
  for (g in groups) {
    cat("\nGroup:", g, "\n")
    sub <- x$risk[x$risk$group == g, ]
    for (i in seq_len(nrow(sub))) {
      ci_str <- ""
      if ("ci_lower" %in% names(sub) && !is.na(sub$ci_lower[i])) {
        ci_str <- sprintf(" (95%% CI: %.4f - %.4f)",
                          sub$ci_lower[i], sub$ci_upper[i])
      }
      cat(sprintf("  t = %g: risk = %.4f%s\n",
                  sub$time[i], sub$risk[i], ci_str))
    }
  }
  invisible(x)
}

#' @export
plot.aipw <- function(x, effect = c("risk", "RD", "RR"), ...) {
  effect <- match.arg(effect)
  .plot_risk_curves(x$risk, effect = effect, title = "AIPW Risk Curves")
}
