#' G-Computation Estimator for Cumulative Risk Curves
#'
#' Estimates cumulative risk curves using the g-computation (parametric
#' g-formula) approach. Fits an outcome model (Cox PH by default) and
#' generates predicted cumulative incidence under each treatment regime
#' via plug-in prediction.
#'
#' @param spec A `cr_spec` object with at least outcome and treatment
#'   identified. The outcome should have a formula for model fitting.
#' @param risk_time Numeric vector of times at which to evaluate risk.
#' @param outcome_formula Optional formula override for the outcome model.
#'   If `NULL`, uses the formula from `identify_outcome()`.
#' @param nboot Number of bootstrap replicates for CIs.
#' @param seed Random seed.
#' @param ... Additional arguments.
#'
#' @return An object of class `gcomp` (inherits from `cr_result`).
#'
#' @details
#' ## Method
#' 1. Fit a Cox PH model for the outcome, including treatment and covariates.
#' 2. For each treatment level `a`, create a counterfactual dataset where
#'    all subjects receive treatment `a`.
#' 3. Predict the baseline survival and compute mean cumulative incidence
#'    across all subjects in the counterfactual data.
#' 4. The g-computation estimate is the average predicted risk at each time.
#'
#' @examples
#' \dontrun{
#' fit_gc <- specify_models(data = example1) |>
#'   identify_outcome(event, formula = ~ treatment + age + sex,
#'                    type = "time_to_event") |>
#'   identify_treatment(treatment) |>
#'   estimate_gcomprisk(risk_time = c(6, 12, 24))
#' }
#'
#' @export
estimate_gcomprisk <- function(spec, risk_time = NULL, outcome_formula = NULL,
                               nboot = 0, seed = 42, ...) {
  if (!inherits(spec, "cr_spec")) {
    stop("`spec` must be a cr_spec object.", call. = FALSE)
  }
  if (is.null(spec$outcome)) {
    stop("Outcome not specified.", call. = FALSE)
  }

  data <- spec$data
  time_var <- .get_time_var(spec)
  event_var <- spec$outcome$name
  event_ind <- .build_event_indicator(data, spec)

  # Build outcome formula
  if (!is.null(outcome_formula)) {
    out_fml <- outcome_formula
  } else if (!is.null(spec$outcome$formula)) {
    rhs <- deparse(spec$outcome$formula[[2]])
    out_fml <- as.formula(paste("survival::Surv(", time_var, ",", "event_ind) ~", rhs))
  } else {
    # Use treatment variable as the only predictor
    if (!is.null(spec$treatment)) {
      out_fml <- as.formula(paste("survival::Surv(", time_var, ",", "event_ind) ~",
                                  spec$treatment$name))
    } else {
      out_fml <- as.formula(paste("survival::Surv(", time_var, ",", "event_ind) ~ 1"))
    }
  }

  # Fit Cox PH model
  data_fit <- data
  data_fit$event_ind <- event_ind
  cox_fit <- survival::coxph(out_fml, data = data_fit)

  # Determine treatment groups
  if (!is.null(spec$treatment)) {
    trt_var <- spec$treatment$name
    trt_levels <- sort(unique(data[[trt_var]]))
  } else {
    trt_levels <- "overall"
  }

  # Predict cumulative risk under each regime
  all_times <- sort(unique(data[[time_var]][event_ind == 1]))
  if (is.null(risk_time)) risk_time <- all_times

  risk_list <- list()
  for (a in trt_levels) {
    if (a != "overall") {
      cf_data <- data_fit
      cf_data[[trt_var]] <- a
    } else {
      cf_data <- data_fit
    }

    # Get survival predictions
    sf <- survival::survfit(cox_fit, newdata = cf_data)

    # Average predicted survival across all individuals
    if (inherits(sf, "survfit")) {
      # sf$surv is matrix (times x individuals) or vector
      if (is.matrix(sf$surv)) {
        mean_surv <- rowMeans(sf$surv)
      } else {
        mean_surv <- sf$surv
      }
      sf_times <- sf$time
    } else {
      sf_times <- all_times
      mean_surv <- rep(1, length(sf_times))
    }

    # Interpolate to requested times
    mean_risk <- 1 - approx(c(0, sf_times), c(1, mean_surv),
                            xout = risk_time, method = "constant",
                            rule = 2, f = 0)$y

    risk_list[[as.character(a)]] <- data.frame(
      time = risk_time,
      risk = mean_risk,
      surv = 1 - mean_risk,
      group = as.character(a),
      stringsAsFactors = FALSE
    )
  }

  risk_df <- do.call(rbind, risk_list)
  rownames(risk_df) <- NULL

  # Bootstrap
  boot_results <- NULL
  if (nboot > 0) {
    set.seed(seed)
    boot_results <- .bootstrap_gcomp(spec, nboot, risk_time, outcome_formula)
    risk_df <- .attach_boot_ci(risk_df, boot_results)
  }

  result <- list(
    risk = risk_df,
    spec = spec,
    outcome_model = cox_fit,
    nboot = nboot,
    boot_results = boot_results,
    call = match.call()
  )
  class(result) <- c("gcomp", "cr_result")
  result
}

#' @keywords internal
.bootstrap_gcomp <- function(spec, nboot, risk_time, outcome_formula) {
  data <- spec$data
  n <- nrow(data)
  boot_risks <- list()

  for (b in seq_len(nboot)) {
    idx <- sample(n, n, replace = TRUE)
    boot_spec <- spec
    boot_spec$data <- data[idx, , drop = FALSE]
    tryCatch({
      boot_fit <- estimate_gcomprisk(boot_spec, risk_time = risk_time,
                                     outcome_formula = outcome_formula,
                                     nboot = 0)
      boot_risks[[b]] <- boot_fit$risk
    }, error = function(e) {
      boot_risks[[b]] <<- NULL
    })
  }
  boot_risks[!vapply(boot_risks, is.null, logical(1))]
}

#' @export
print.gcomp <- function(x, ...) {
  cat("G-Computation Cumulative Risk Estimate\n")
  cat("=======================================\n")
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
plot.gcomp <- function(x, effect = c("risk", "RD", "RR"), ...) {
  effect <- match.arg(effect)
  .plot_risk_curves(x$risk, effect = effect, title = "G-Computation Risk Curves")
}
