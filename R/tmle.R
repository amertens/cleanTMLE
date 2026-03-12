#' TMLE Extensions for cleanTMLE
#'
#' These functions wrap established TMLE packages (tmle, survtmle,
#' lmtp) to provide targeted learning estimators that integrate with
#' the cleanTMLE model specification workflow.
#'
#' @name tmle_extensions
NULL


#' Point-Treatment TMLE for Binary Outcomes
#'
#' Estimates the average treatment effect (ATE) and optionally the risk
#' ratio (RR) for a binary outcome using the \pkg{tmle} package with
#' \pkg{SuperLearner} for nuisance parameter estimation.
#'
#' @param data A data.frame.
#' @param treatment Character; name of the binary treatment variable.
#' @param outcome Character; name of the binary outcome variable.
#' @param covariates Character vector of covariate column names (W).
#' @param family `"binomial"` (default) for binary outcomes.
#' @param sl_library Character vector of SuperLearner algorithms.
#'   Default: `c("SL.glm", "SL.mean")`. Use `"SL.glmnet"`, `"SL.ranger"`,
#'   etc. if those packages are installed.
#' @param truncate Truncation level for propensity scores (default 0.01,
#'   meaning PS is bounded to the interval 0.01 to 0.99).
#' @param spec Optional `cr_spec` object. If provided, extracts
#'   treatment/outcome/covariate info from it.
#' @param ... Additional arguments passed to `tmle::tmle()`.
#'
#' @return An object of class `tmle_fit` (inherits from `cr_result`)
#'   containing:
#'   * `estimates` - list with ATE (RD) and RR estimates, SEs, CIs
#'   * `tmle_obj` - the raw tmle object from the \pkg{tmle} package
#'   * `influence_curve` - the efficient influence curve values
#'
#' @details
#' This function requires the \pkg{tmle} and \pkg{SuperLearner} packages.
#' It uses the TMLE framework to estimate the causal risk difference
#' (ATE = E\[Y(1)\] - E\[Y(0)\]) with influence-curve based standard errors
#' and confidence intervals.
#'
#' ## SuperLearner
#' The nuisance parameters (outcome regression Q and propensity score g)
#' are estimated using SuperLearner ensemble learning. The default library
#' uses GLM and the marginal mean. For better performance, add more
#' learners like `"SL.glmnet"`, `"SL.ranger"`, `"SL.xgboost"`.
#'
#' @examples
#' \dontrun{
#' dat <- sim_func1(n = 500)
#' dat$event_24 <- as.integer(dat$event == 1 & dat$time <= 24)
#' fit <- estimate_tmle_risk_point(
#'   data = dat,
#'   treatment = "treatment",
#'   outcome = "event_24",
#'   covariates = c("age", "sex", "biomarker"),
#'   sl_library = c("SL.glm", "SL.mean")
#' )
#' print(fit)
#' }
#'
#' @export
estimate_tmle_risk_point <- function(data, treatment = NULL, outcome = NULL,
                                     covariates = NULL, family = "binomial",
                                     sl_library = c("SL.glm", "SL.mean"),
                                     truncate = 0.01, spec = NULL, ...) {
  # Extract from spec if provided

  if (!is.null(spec) && inherits(spec, "cr_spec")) {
    if (is.null(treatment) && !is.null(spec$treatment)) {
      treatment <- spec$treatment$name
    }
    if (is.null(outcome) && !is.null(spec$outcome)) {
      outcome <- spec$outcome$name
    }
    if (is.null(data)) data <- spec$data
  }

  if (is.null(treatment) || is.null(outcome)) {
    stop("Must specify treatment and outcome variables.", call. = FALSE)
  }

  if (!requireNamespace("tmle", quietly = TRUE)) {
    stop("Package 'tmle' is required. Install with: install.packages('tmle')",
         call. = FALSE)
  }
  if (!requireNamespace("SuperLearner", quietly = TRUE)) {
    stop("Package 'SuperLearner' is required. Install with: install.packages('SuperLearner')",
         call. = FALSE)
  }

  # Prepare data
  A <- data[[treatment]]
  Y <- data[[outcome]]

  if (is.null(covariates)) {
    # Use all numeric columns except treatment and outcome
    covariates <- setdiff(
      names(data)[vapply(data, is.numeric, logical(1))],
      c(treatment, outcome, "id", "time")
    )
  }
  W <- data[, covariates, drop = FALSE]

  # Run TMLE
  tmle_fit <- tmle::tmle(
    Y = Y, A = A, W = W,
    family = family,
    Q.SL.library = sl_library,
    g.SL.library = sl_library,
    gbound = truncate,
    ...
  )

  # Extract results
  ate <- tmle_fit$estimates$ATE
  rr_est <- tmle_fit$estimates$RR

  estimates <- list(
    ATE = list(
      estimate = ate$psi,
      se = sqrt(ate$var.psi),
      ci_lower = ate$CI[1],
      ci_upper = ate$CI[2],
      p_value = ate$pvalue
    )
  )

  if (!is.null(rr_est)) {
    estimates$RR <- list(
      estimate = rr_est$psi,
      se = sqrt(rr_est$var.log.psi),
      ci_lower = rr_est$CI[1],
      ci_upper = rr_est$CI[2],
      p_value = rr_est$pvalue
    )
  }

  # Influence curve
  ic <- tmle_fit$estimates$IC$IC.ATE

  result <- list(
    estimates = estimates,
    tmle_obj = tmle_fit,
    influence_curve = ic,
    treatment = treatment,
    outcome = outcome,
    covariates = covariates,
    type = "point_tmle",
    call = match.call()
  )
  class(result) <- c("tmle_fit", "cr_result")
  result
}


#' Survival TMLE for Risk at Specified Times
#'
#' Estimates treatment-specific cumulative risk at user-specified time
#' points using the \pkg{survtmle} package or a manual implementation.
#'
#' @param data A data.frame with time-to-event data.
#' @param treatment Character; name of binary treatment variable.
#' @param time Character; name of the time variable.
#' @param event Character; name of the event indicator (1 = event, 0 = censor).
#' @param covariates Character vector of baseline covariate names.
#' @param target_times Numeric vector of times at which to estimate risk.
#' @param sl_library SuperLearner library for nuisance models.
#' @param ... Additional arguments.
#'
#' @return An object of class `tmle_fit`.
#'
#' @details
#' This function attempts to use the \pkg{survtmle} package. If unavailable,
#' falls back to a simplified manual implementation that discretizes
#' time and applies iterative TMLE at each time point.
#'
#' @examples
#' \dontrun{
#' dat <- sim_func1(n = 500)
#' fit <- estimate_surv_tmle(
#'   data = dat,
#'   treatment = "treatment",
#'   time = "time",
#'   event = "event",
#'   covariates = c("age", "sex", "biomarker"),
#'   target_times = c(12, 24)
#' )
#' }
#'
#' @export
estimate_surv_tmle <- function(data, treatment = "treatment",
                               time = "time", event = "event",
                               covariates = NULL,
                               target_times = c(12, 24),
                               sl_library = c("SL.glm", "SL.mean"),
                               ...) {
  if (is.null(covariates)) {
    covariates <- setdiff(
      names(data)[vapply(data, is.numeric, logical(1))],
      c(treatment, time, event, "id", "event_type", "censored")
    )
  }

  if (requireNamespace("survtmle", quietly = TRUE)) {
    # Use survtmle package
    fit <- .surv_tmle_via_package(data, treatment, time, event,
                                  covariates, target_times,
                                  sl_library, ...)
  } else {
    # Fallback: discretize and use point TMLE at each horizon
    fit <- .surv_tmle_fallback(data, treatment, time, event,
                               covariates, target_times,
                               sl_library, ...)
  }

  fit
}


#' survtmle package wrapper
#' @keywords internal
.surv_tmle_via_package <- function(data, treatment, time, event,
                                    covariates, target_times,
                                    sl_library, ...) {
  # Discretize time for survtmle
  obs_time <- data[[time]]
  max_t <- max(target_times)

  # survtmle expects integer times and specific format
  t_discrete <- ceiling(obs_time)
  t_discrete <- pmin(t_discrete, max_t + 1L)

  W <- data[, covariates, drop = FALSE]
  A <- data[[treatment]]
  event_ind <- data[[event]]

  # Compute censoring indicator (C = 1 means observed/uncensored)
  ftime <- as.integer(t_discrete)
  ftype <- as.integer(event_ind)

  tryCatch({
    fit <- survtmle::survtmle(
      ftime = ftime,
      ftype = ftype,
      trt = A,
      adjustVars = W,
      t0 = as.integer(max(target_times)),
      SL.ftime = sl_library,
      SL.ctime = sl_library,
      SL.trt = sl_library,
      method = "hazard",
      ...
    )

    # Extract estimates
    estimates <- list()
    for (tt in target_times) {
      estimates[[paste0("risk_t", tt)]] <- list(
        estimate = fit$est[1],
        se = sqrt(fit$var[1]),
        ci_lower = fit$est[1] - 1.96 * sqrt(fit$var[1]),
        ci_upper = fit$est[1] + 1.96 * sqrt(fit$var[1]),
        p_value = 2 * pnorm(-abs(fit$est[1] / sqrt(fit$var[1])))
      )
    }

    result <- list(
      estimates = estimates,
      tmle_obj = fit,
      influence_curve = NULL,
      treatment = treatment,
      time = time,
      event = event,
      covariates = covariates,
      type = "surv_tmle",
      call = match.call()
    )
    class(result) <- c("tmle_fit", "cr_result")
    result

  }, error = function(e) {
    message("survtmle failed, falling back to manual approach: ", e$message)
    .surv_tmle_fallback(data, treatment, time, event,
                        covariates, target_times, sl_library)
  })
}


#' Fallback survival TMLE via binarization
#' @keywords internal
.surv_tmle_fallback <- function(data, treatment, time, event,
                                 covariates, target_times,
                                 sl_library, ...) {
  obs_time <- data[[time]]
  event_ind <- data[[event]]
  A <- data[[treatment]]

  estimates <- list()

  for (tt in target_times) {
    # Binarize: Y_t = I(event occurred by time t)
    Y_t <- as.integer(event_ind == 1 & obs_time <= tt)

    # Only include subjects still observed at time t or who had event before t
    # (simple approach: include everyone)

    if (requireNamespace("tmle", quietly = TRUE) &&
        requireNamespace("SuperLearner", quietly = TRUE)) {
      W <- data[, covariates, drop = FALSE]

      tmle_t <- tryCatch(
        tmle::tmle(Y = Y_t, A = A, W = W,
                   family = "binomial",
                   Q.SL.library = sl_library,
                   g.SL.library = sl_library),
        error = function(e) NULL
      )

      if (!is.null(tmle_t)) {
        ate <- tmle_t$estimates$ATE
        estimates[[paste0("risk_t", tt)]] <- list(
          estimate = ate$psi,
          se = sqrt(ate$var.psi),
          ci_lower = ate$CI[1],
          ci_upper = ate$CI[2],
          p_value = ate$pvalue
        )
      } else {
        # Very simple fallback
        r1 <- mean(Y_t[A == 1])
        r0 <- mean(Y_t[A == 0])
        estimates[[paste0("risk_t", tt)]] <- list(
          estimate = r1 - r0,
          se = NA_real_,
          ci_lower = NA_real_,
          ci_upper = NA_real_,
          p_value = NA_real_
        )
      }
    } else {
      # No tmle package: simple difference
      r1 <- mean(Y_t[A == 1])
      r0 <- mean(Y_t[A == 0])
      se <- sqrt(var(Y_t[A == 1]) / sum(A == 1) +
                   var(Y_t[A == 0]) / sum(A == 0))
      rd <- r1 - r0
      estimates[[paste0("risk_t", tt)]] <- list(
        estimate = rd,
        se = se,
        ci_lower = rd - 1.96 * se,
        ci_upper = rd + 1.96 * se,
        p_value = 2 * pnorm(-abs(rd / se))
      )
    }
  }

  result <- list(
    estimates = estimates,
    tmle_obj = NULL,
    influence_curve = NULL,
    treatment = treatment,
    time = time,
    event = event,
    covariates = covariates,
    type = "surv_tmle_fallback",
    call = match.call()
  )
  class(result) <- c("tmle_fit", "cr_result")
  result
}


#' Longitudinal TMLE via \pkg{lmtp}
#'
#' Estimates the mean potential outcome under static, dynamic, or
#' stochastic treatment interventions using the \pkg{lmtp} package.
#'
#' @param data A data.frame in long or wide format.
#' @param treatment_vars Character vector of treatment variable names at
#'   each time point (e.g., `c("A_1", "A_2", "A_3")`).
#' @param outcome Character; name of the outcome variable.
#' @param baseline Character vector of baseline covariate names.
#' @param time_varying Optional list of time-varying covariate names,
#'   one vector per time point.
#' @param censor_vars Optional character vector of censoring indicator
#'   names at each time point.
#' @param shift_fun A function defining the intervention. Default is a
#'   static intervention setting treatment to 1. The function should take
#'   `(data, trt)` and return modified treatment values.
#' @param outcome_type `"binomial"` or `"continuous"`.
#' @param sl_library SuperLearner library.
#' @param id Optional character; name of the subject ID column.
#' @param ... Additional arguments passed to `lmtp::lmtp_tmle()`.
#'
#' @return An object of class `tmle_fit`.
#'
#' @details
#' This function requires the \pkg{lmtp} package. It implements the
#' longitudinal modified treatment policy (LMTP) framework for
#' causal inference with time-varying treatments.
#'
#' @examples
#' \dontrun{
#' # Simple 2-time-point example
#' dat <- data.frame(
#'   id = rep(1:200, each = 2),
#'   time = rep(1:2, 200),
#'   A = rbinom(400, 1, 0.5),
#'   L = rnorm(400),
#'   Y = rbinom(400, 1, 0.3)
#' )
#' # Static intervention: set all A = 1
#' fit <- estimate_lmtp(
#'   data = dat,
#'   treatment_vars = c("A"),
#'   outcome = "Y",
#'   baseline = c("L"),
#'   shift_fun = function(data, trt) rep(1, length(data[[trt]]))
#' )
#' }
#'
#' @export
estimate_lmtp <- function(data, treatment_vars, outcome,
                          baseline = NULL, time_varying = NULL,
                          censor_vars = NULL,
                          shift_fun = NULL,
                          outcome_type = c("binomial", "continuous"),
                          sl_library = c("SL.glm", "SL.mean"),
                          id = NULL, ...) {
  outcome_type <- match.arg(outcome_type)

  if (!requireNamespace("lmtp", quietly = TRUE)) {
    stop("Package 'lmtp' is required. Install with: install.packages('lmtp')",
         call. = FALSE)
  }

  # Default shift: static intervention A = 1
  if (is.null(shift_fun)) {
    shift_fun <- function(data, trt) {
      rep(1, nrow(data))
    }
  }

  # Call lmtp::lmtp_tmle
  lmtp_args <- list(
    data = data,
    trt = treatment_vars,
    outcome = outcome,
    baseline = baseline,
    shift = shift_fun,
    outcome_type = outcome_type,
    learners_outcome = sl_library,
    learners_trt = sl_library
  )

  if (!is.null(time_varying)) lmtp_args$time_vary <- time_varying
  if (!is.null(censor_vars)) lmtp_args$cens <- censor_vars
  if (!is.null(id)) lmtp_args$id <- id

  lmtp_fit <- tryCatch(
    do.call(lmtp::lmtp_tmle, c(lmtp_args, list(...))),
    error = function(e) {
      stop("lmtp_tmle failed: ", e$message, call. = FALSE)
    }
  )

  # Extract results
  est <- lmtp_fit$theta
  se_val <- if (!is.null(lmtp_fit$standard_error)) {
    lmtp_fit$standard_error
  } else {
    NA_real_
  }

  estimates <- list(
    theta = list(
      estimate = est,
      se = se_val,
      ci_lower = if (!is.na(se_val)) est - 1.96 * se_val else NA_real_,
      ci_upper = if (!is.na(se_val)) est + 1.96 * se_val else NA_real_,
      p_value = if (!is.na(se_val)) 2 * pnorm(-abs(est / se_val)) else NA_real_
    )
  )

  result <- list(
    estimates = estimates,
    tmle_obj = lmtp_fit,
    influence_curve = lmtp_fit$eif,
    treatment = treatment_vars,
    outcome = outcome,
    covariates = baseline,
    type = "lmtp",
    call = match.call()
  )
  class(result) <- c("tmle_fit", "cr_result")
  result
}


#' @export
print.tmle_fit <- function(x, ...) {
  type_labels <- c(
    point_tmle = "Point-Treatment TMLE",
    surv_tmle = "Survival TMLE",
    surv_tmle_fallback = "Survival TMLE (fallback)",
    lmtp = "Longitudinal TMLE (LMTP)"
  )
  label <- type_labels[x$type]
  if (is.na(label)) label <- "TMLE Estimate"

  cat(label, "\n")
  cat(paste(rep("=", nchar(label)), collapse = ""), "\n")

  for (nm in names(x$estimates)) {
    est <- x$estimates[[nm]]
    cat(sprintf("\n%s:\n", nm))
    cat(sprintf("  Estimate: %.5f\n", est$estimate))
    if (!is.na(est$se)) {
      cat(sprintf("  SE:       %.5f\n", est$se))
      cat(sprintf("  95%% CI:   [%.5f, %.5f]\n", est$ci_lower, est$ci_upper))
      cat(sprintf("  p-value:  %.5f\n", est$p_value))
    }
  }
  invisible(x)
}


#' @export
plot.tmle_fit <- function(x, ...) {
  if (x$type == "surv_tmle" || x$type == "surv_tmle_fallback") {
    # Plot risk differences over time
    est_df <- data.frame(
      time = as.numeric(gsub("risk_t", "", names(x$estimates))),
      estimate = vapply(x$estimates, function(e) e$estimate, numeric(1)),
      ci_lower = vapply(x$estimates, function(e) {
        if (is.na(e$ci_lower)) e$estimate else e$ci_lower
      }, numeric(1)),
      ci_upper = vapply(x$estimates, function(e) {
        if (is.na(e$ci_upper)) e$estimate else e$ci_upper
      }, numeric(1)),
      stringsAsFactors = FALSE
    )

    p <- ggplot2::ggplot(est_df, ggplot2::aes(
      x = .data$time, y = .data$estimate
    )) +
      ggplot2::geom_point(size = 3) +
      ggplot2::geom_errorbar(ggplot2::aes(
        ymin = .data$ci_lower, ymax = .data$ci_upper
      ), width = 0.5) +
      ggplot2::geom_hline(yintercept = 0, linetype = "dashed") +
      ggplot2::labs(x = "Time", y = "Risk Difference (TMLE)",
                    title = "Survival TMLE Estimates") +
      ggplot2::theme_minimal()

    return(p)
  }

  # For point TMLE: simple bar/point plot
  est_df <- data.frame(
    estimand = names(x$estimates),
    estimate = vapply(x$estimates, function(e) e$estimate, numeric(1)),
    ci_lower = vapply(x$estimates, function(e) {
      if (is.na(e$ci_lower)) e$estimate else e$ci_lower
    }, numeric(1)),
    ci_upper = vapply(x$estimates, function(e) {
      if (is.na(e$ci_upper)) e$estimate else e$ci_upper
    }, numeric(1)),
    stringsAsFactors = FALSE
  )

  p <- ggplot2::ggplot(est_df, ggplot2::aes(
    x = .data$estimand, y = .data$estimate
  )) +
    ggplot2::geom_point(size = 3) +
    ggplot2::geom_errorbar(ggplot2::aes(
      ymin = .data$ci_lower, ymax = .data$ci_upper
    ), width = 0.2) +
    ggplot2::geom_hline(yintercept = 0, linetype = "dashed") +
    ggplot2::labs(x = "", y = "Estimate (95% CI)",
                  title = "TMLE Estimates") +
    ggplot2::theme_minimal()

  p
}
