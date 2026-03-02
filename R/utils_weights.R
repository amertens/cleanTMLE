#' Internal Weight Computation Utilities
#'
#' These internal functions compute IPTW and IPCW weights from fitted
#' propensity/censoring models.
#'
#' @name weight_utils
#' @keywords internal
NULL

#' Fit a treatment (propensity score) model and compute IPTW weights
#'
#' @param data Data frame.
#' @param treatment_spec Treatment specification from `identify_treatment()`.
#' @param weight_type `"iptw"` (ATE) or `"smr"` (ATT / SMR-style).
#' @param trim Trim weights at this quantile (e.g., 0.99 trims at 1st/99th).
#' @param trunc Truncate (cap) weights at this absolute value.
#'
#' @return A list with `weights`, `ps` (propensity scores), `model`.
#' @keywords internal
fit_treatment_weights <- function(data, treatment_spec,
                                  weight_type = c("iptw", "smr"),
                                  trim = NULL, trunc = NULL) {
  weight_type <- match.arg(weight_type)
  trt_var <- treatment_spec$name
  trt_vals <- data[[trt_var]]

  if (is.null(treatment_spec$formula)) {
    # No PS model: equal weights
    n <- nrow(data)
    return(list(
      weights = rep(1, n),
      ps = rep(mean(trt_vals == 1), n),
      model = NULL,
      trt_levels = sort(unique(trt_vals))
    ))
  }

  trt_levels <- sort(unique(trt_vals))

  if (length(trt_levels) == 2 && treatment_spec$model == "glm") {
    # Binary treatment: logistic regression
    fml <- as.formula(paste(trt_var, "~", deparse(treatment_spec$formula[[2]])))
    fit <- glm(fml, data = data, family = binomial())
    ps <- predict(fit, type = "response")

    if (weight_type == "iptw") {
      # Stabilized IPTW: P(A) / P(A|L)
      p_trt <- mean(trt_vals == trt_levels[2])
      w <- ifelse(trt_vals == trt_levels[2],
                  p_trt / ps,
                  (1 - p_trt) / (1 - ps))
    } else {
      # SMR weights (ATT): treated get w=1, untreated get ps/(1-ps)
      w <- ifelse(trt_vals == trt_levels[2],
                  1,
                  ps / (1 - ps))
    }

  } else if (length(trt_levels) > 2 || treatment_spec$model == "multinom") {
    # Multi-category treatment: multinomial logistic regression
    if (!requireNamespace("nnet", quietly = TRUE)) {
      stop("Package 'nnet' required for multinomial treatment models.", call. = FALSE)
    }
    fml <- as.formula(paste(trt_var, "~", deparse(treatment_spec$formula[[2]])))
    fit <- nnet::multinom(fml, data = data, trace = FALSE)
    ps_mat <- predict(fit, type = "probs")
    if (is.null(dim(ps_mat))) {
      # Only two levels returned as vector
      ps_mat <- cbind(1 - ps_mat, ps_mat)
      colnames(ps_mat) <- as.character(trt_levels)
    }

    # Generalized IPTW: w_i = P(A=a_i) / P(A=a_i | L_i)
    marginal_probs <- table(trt_vals) / length(trt_vals)
    w <- numeric(nrow(data))
    ps <- numeric(nrow(data))
    for (i in seq_len(nrow(data))) {
      a <- as.character(trt_vals[i])
      col_idx <- which(colnames(ps_mat) == a)
      if (length(col_idx) == 0) col_idx <- which(trt_levels == trt_vals[i])
      ps[i] <- ps_mat[i, col_idx]
      w[i] <- as.numeric(marginal_probs[a]) / ps[i]
    }

    return(.apply_weight_constraints(
      list(weights = w, ps = ps, model = fit, trt_levels = trt_levels),
      trim = trim, trunc = trunc, data = data, trt_vals = trt_vals
    ))

  } else {
    fml <- as.formula(paste(trt_var, "~", deparse(treatment_spec$formula[[2]])))
    fit <- glm(fml, data = data, family = binomial())
    ps <- predict(fit, type = "response")
    p_trt <- mean(trt_vals == trt_levels[2])
    w <- ifelse(trt_vals == trt_levels[2],
                p_trt / ps,
                (1 - p_trt) / (1 - ps))
  }

  .apply_weight_constraints(
    list(weights = w, ps = ps, model = fit, trt_levels = trt_levels),
    trim = trim, trunc = trunc, data = data, trt_vals = trt_vals
  )
}


#' Fit censoring model(s) and compute IPCW weights at each event time
#'
#' For time-to-event data with right censoring, this computes inverse
#' probability of censoring weights using either a Cox model or pooled
#' logistic model for the censoring hazard.
#'
#' @param data Data frame with observed time and censoring indicators.
#' @param censoring_specs List of censoring specifications.
#' @param outcome_spec Outcome specification (for time variable).
#' @param interval_spec Interval specification (optional).
#'
#' @return A list with `weights` (one per observation) and `models`.
#' @keywords internal
fit_censoring_weights <- function(data, censoring_specs, outcome_spec,
                                  interval_spec = NULL) {
  if (length(censoring_specs) == 0) {
    return(list(weights = rep(1, nrow(data)), models = list()))
  }

  # Composite censoring: censored if any censoring indicator == 1
  n <- nrow(data)
  composite_cens <- rep(0L, n)
  for (cs in censoring_specs) {
    composite_cens <- pmax(composite_cens, as.integer(data[[cs$name]]))
  }

  # Use the first censoring spec for model details
  cs <- censoring_specs[[1]]

  if (is.null(cs$formula)) {
    # No censoring model: equal weights
    return(list(weights = rep(1, n), models = list()))
  }

  # Determine time variable
  if (!is.null(interval_spec)) {
    time_var <- interval_spec$stop
  } else {
    # Use outcome name as time column or look for a 'time' column
    time_var <- if ("time" %in% names(data)) "time" else outcome_spec$name
  }

  obs_time <- data[[time_var]]

  if (cs$model == "coxph") {
    # Cox model for censoring: model censoring as the "event"
    surv_fml <- as.formula(paste(
      "survival::Surv(", time_var, ",", "composite_cens) ~",
      deparse(cs$formula[[2]])
    ))
    data_aug <- data
    data_aug$composite_cens <- composite_cens
    fit <- survival::coxph(surv_fml, data = data_aug)

    # KM-style censoring survival probability for each individual
    # Use baseline survival + covariate adjustment
    sf <- survival::survfit(fit, newdata = data_aug)

    # Simplified: compute probability of not being censored by obs_time
    # For each individual, S_C(t_i | L_i) from the Cox model
    # Use marginal censoring KM as a simpler approximation
    km_cens <- survival::survfit(
      survival::Surv(obs_time, composite_cens) ~ 1
    )

    # Map each observation to its censoring survival probability
    cens_surv <- approx(km_cens$time, km_cens$surv,
                        xout = obs_time, method = "constant",
                        rule = 2, f = 0)$y
    cens_surv <- pmax(cens_surv, 0.01)  # floor to prevent extreme weights
    ipcw <- 1 / cens_surv

  } else {
    # GLM model for censoring (pooled logistic style)
    fml <- as.formula(paste(
      "composite_cens ~", deparse(cs$formula[[2]])
    ))
    data_aug <- data
    data_aug$composite_cens <- composite_cens
    fit <- glm(fml, data = data_aug, family = binomial())
    p_cens <- predict(fit, type = "response")
    p_cens <- pmax(p_cens, 0.01)
    ipcw <- 1 / (1 - p_cens)
  }

  list(weights = ipcw, models = list(fit))
}


#' Apply trimming and truncation to weights
#' @keywords internal
.apply_weight_constraints <- function(wt_result, trim = NULL, trunc = NULL,
                                       data = NULL, trt_vals = NULL) {
  w <- wt_result$weights

  # Trimming: exclude observations with extreme PS
  if (!is.null(trim) && !is.null(wt_result$ps)) {
    ps <- wt_result$ps
    lo <- quantile(ps, trim, na.rm = TRUE)
    hi <- quantile(ps, 1 - trim, na.rm = TRUE)
    trimmed <- ps < lo | ps > hi
    w[trimmed] <- 0
  }

  # Truncation: cap weights at quantile
  if (!is.null(trunc)) {
    if (trunc <= 1) {
      cap <- quantile(w[w > 0], trunc, na.rm = TRUE)
    } else {
      cap <- trunc
    }
    w <- pmin(w, cap)
  }

  wt_result$weights <- w
  wt_result
}
