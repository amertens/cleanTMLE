#' Plot Methods for cleanTMLE Result Objects
#'
#' @name plot_methods
#' @param x A result object (cumrisk, gcomp, aipw, or hr).
#' @param effect For risk curve objects: `"risk"` (cumulative risk curves),
#'   `"RD"` (risk difference), or `"RR"` (risk ratio).
#' @param ... Additional arguments passed to ggplot.
NULL

#' @export
print.cumrisk <- function(x, ...) {
  cat("IPW Cumulative Risk Estimate\n")
  cat("============================\n")
  cat("Weight type:", x$weight_type, "\n")
  if (!is.null(x$trim)) cat("Trim:", x$trim, "\n")
  if (!is.null(x$trunc)) cat("Trunc:", x$trunc, "\n")
  if (x$nboot > 0) cat("Bootstrap reps:", x$nboot, "\n")
  cat("\n")

  groups <- unique(x$risk$group)
  for (g in groups) {
    cat("Group:", g, "\n")
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


#' @rdname plot_methods
#' @export
plot.cumrisk <- function(x, effect = c("risk", "RD", "RR"), ...) {
  effect <- match.arg(effect)
  .plot_risk_curves(x$risk, effect = effect, title = "IPW Cumulative Risk Curves")
}


#' Plot risk curves with optional effect measures
#' @keywords internal
.plot_risk_curves <- function(risk_df, effect = "risk", title = "Risk Curves") {
  groups <- unique(risk_df$group)

  if (effect == "risk") {
    p <- ggplot2::ggplot(risk_df, ggplot2::aes(
      x = .data$time, y = .data$risk, color = .data$group
    )) +
      ggplot2::geom_step(linewidth = 1) +
      ggplot2::labs(x = "Time", y = "Cumulative Risk", color = "Group",
                    title = title) +
      ggplot2::theme_minimal()

    if ("ci_lower" %in% names(risk_df)) {
      has_ci <- !all(is.na(risk_df$ci_lower))
      if (has_ci) {
        p <- p + ggplot2::geom_ribbon(
          ggplot2::aes(ymin = .data$ci_lower, ymax = .data$ci_upper,
                       fill = .data$group),
          alpha = 0.2, color = NA
        )
      }
    }

  } else if (effect == "RD" && length(groups) == 2) {
    # Risk difference
    rd_df <- .compute_contrasts(risk_df, measure = "RD")
    p <- ggplot2::ggplot(rd_df, ggplot2::aes(
      x = .data$time, y = .data$estimate
    )) +
      ggplot2::geom_step(linewidth = 1) +
      ggplot2::geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
      ggplot2::labs(x = "Time", y = "Risk Difference",
                    title = paste(title, "- Risk Difference")) +
      ggplot2::theme_minimal()

    if ("ci_lower" %in% names(rd_df)) {
      p <- p + ggplot2::geom_ribbon(
        ggplot2::aes(ymin = .data$ci_lower, ymax = .data$ci_upper),
        alpha = 0.2
      )
    }

  } else if (effect == "RR" && length(groups) == 2) {
    # Risk ratio
    rr_df <- .compute_contrasts(risk_df, measure = "RR")
    p <- ggplot2::ggplot(rr_df, ggplot2::aes(
      x = .data$time, y = .data$estimate
    )) +
      ggplot2::geom_step(linewidth = 1) +
      ggplot2::geom_hline(yintercept = 1, linetype = "dashed", color = "grey50") +
      ggplot2::labs(x = "Time", y = "Risk Ratio",
                    title = paste(title, "- Risk Ratio")) +
      ggplot2::theme_minimal()

  } else {
    # Default to risk curves if contrast not possible
    p <- ggplot2::ggplot(risk_df, ggplot2::aes(
      x = .data$time, y = .data$risk, color = .data$group
    )) +
      ggplot2::geom_step(linewidth = 1) +
      ggplot2::labs(x = "Time", y = "Cumulative Risk", color = "Group",
                    title = title) +
      ggplot2::theme_minimal()
  }

  p
}


#' Compute risk contrasts (RD or RR) between two groups
#' @keywords internal
.compute_contrasts <- function(risk_df, measure = "RD") {
  groups <- unique(risk_df$group)
  if (length(groups) != 2) {
    stop("Contrasts require exactly 2 groups.", call. = FALSE)
  }

  g1 <- risk_df[risk_df$group == groups[1], ]
  g2 <- risk_df[risk_df$group == groups[2], ]

  # Align on common times
  common_times <- intersect(g1$time, g2$time)
  g1 <- g1[g1$time %in% common_times, ]
  g2 <- g2[g2$time %in% common_times, ]

  if (measure == "RD") {
    est <- g2$risk - g1$risk
  } else {
    est <- ifelse(g1$risk > 0, g2$risk / g1$risk, NA_real_)
  }

  result <- data.frame(time = common_times, estimate = est)

  # Propagate CIs if available (bootstrap delta method approximation)
  if ("ci_lower" %in% names(g1) && "ci_lower" %in% names(g2)) {
    if (measure == "RD") {
      # Simple approximation: combine bootstrap CIs
      if (!all(is.na(g1$ci_lower)) && !all(is.na(g2$ci_lower))) {
        se1 <- (g1$ci_upper - g1$ci_lower) / (2 * 1.96)
        se2 <- (g2$ci_upper - g2$ci_lower) / (2 * 1.96)
        se_rd <- sqrt(se1^2 + se2^2)
        result$ci_lower <- est - 1.96 * se_rd
        result$ci_upper <- est + 1.96 * se_rd
      }
    }
  }

  result
}


#' Forest Plot for Hazard Ratio Results
#'
#' Creates a forest plot showing hazard ratios with confidence intervals
#' from an `hr` object.
#'
#' @param x An `hr` object from [estimate_ipwhr()].
#' @param ... Additional arguments.
#'
#' @return A ggplot object.
#'
#' @export
forest_plot <- function(x, ...) {
  # Dispatch:
  #  - hr object: hazard-ratio forest (legacy behaviour)
  #  - list/data.frame of TMLE/IPTW/match/crude results: comparison forest
  if (inherits(x, "hr")) return(.forest_plot_hr(x))
  if (is.data.frame(x)) return(.forest_plot_df(x, ...))
  if (is.list(x)) return(.forest_plot_results_list(x, ...))
  stop("forest_plot(): unsupported input type. Pass an `hr` object, a ",
       "summarize_cleanroom_results()-style data.frame, or a list of ",
       "fitted workflow objects.", call. = FALSE)
}

.forest_plot_hr <- function(x) {
  ht <- x$hr_table
  ht$term <- factor(ht$term, levels = rev(ht$term))

  ggplot2::ggplot(ht, ggplot2::aes(
    y = .data$term, x = .data$hr,
    xmin = .data$ci_lower, xmax = .data$ci_upper
  )) +
    ggplot2::geom_point(size = 3) +
    ggplot2::geom_errorbarh(height = 0.2) +
    ggplot2::geom_vline(xintercept = 1, linetype = "dashed", color = "grey50") +
    ggplot2::labs(x = "Hazard Ratio (95% CI)", y = "",
                  title = "Forest Plot: IPW Hazard Ratios") +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      axis.text.y = ggplot2::element_text(size = 11),
      plot.title = ggplot2::element_text(face = "bold")
    )
}

.forest_plot_results_list <- function(fits, ...) {
  if (length(fits) == 0L) stop("Empty fits list.", call. = FALSE)
  df <- summarize_cleanroom_results(fits)
  .forest_plot_df(df, ...)
}

.forest_plot_df <- function(df, x_lab = "Risk Difference (95% CI)",
                            null_value = 0, ...) {
  needed <- c("method", "estimate", "ci_lower", "ci_upper")
  miss   <- setdiff(needed, names(df))
  if (length(miss) > 0L)
    stop("forest_plot(): data.frame missing columns: ",
         paste(miss, collapse = ", "), call. = FALSE)

  # Default ordering: TMLE -> IPTW -> Match -> Crude (most-to-least efficient).
  # Anything else lands at the bottom in input order.
  preferred <- c("TMLE", "IPTW", "PS Match", "Match", "AIPW", "G-comp",
                  "Crude")
  ordered_idx <- order(match(df$method,
                              c(preferred, setdiff(df$method, preferred)),
                              nomatch = .Machine$integer.max))
  df <- df[ordered_idx, , drop = FALSE]
  df$method <- factor(df$method, levels = rev(df$method))

  highlight <- as.character(df$method) == "TMLE"

  ggplot2::ggplot(df, ggplot2::aes(
      y = .data$method, x = .data$estimate,
      xmin = .data$ci_lower, xmax = .data$ci_upper)) +
    ggplot2::geom_rect(data = df[highlight, , drop = FALSE],
                       inherit.aes = FALSE,
                       ggplot2::aes(ymin = as.numeric(.data$method) - 0.45,
                                     ymax = as.numeric(.data$method) + 0.45,
                                     xmin = -Inf, xmax = Inf),
                       fill = "grey92") +
    ggplot2::geom_point(size = 3) +
    ggplot2::geom_errorbarh(height = 0.2) +
    ggplot2::geom_vline(xintercept = null_value, linetype = "dashed",
                         colour = "grey50") +
    ggplot2::labs(x = x_lab, y = "",
                  title = "Forest plot of estimator comparison") +
    ggplot2::theme_minimal() +
    ggplot2::theme(axis.text.y = ggplot2::element_text(size = 11),
                    plot.title = ggplot2::element_text(face = "bold"))
}


#' Histogram of Propensity Scores / Weights by Treatment Group
#'
#' Displays the distribution of propensity scores or IPW weights by
#' treatment group.
#'
#' @param x An `ipw` or `cumrisk` object with weight information.
#' @param type `"ps"` for propensity scores, `"weights"` for IPW weights.
#' @param ... Additional arguments passed to hist.
#'
#' @return A ggplot object (invisibly).
#'
#' @export
hist.ipw <- function(x, type = c("ps", "weights"), ...) {
  type <- match.arg(type)
  .hist_weights(x, type = type)
}


#' Heatmap of candidate performance across data-quality threats
#'
#' Plots a tile grid showing RMSE, bias, coverage, or another operating
#' characteristic from a [run_plasmode_dq_stress()] result. Candidates run
#' along the x-axis; scenario-by-severity combinations run along the y-axis,
#' with the undisturbed baseline row at the top so degradation is immediately
#' visible.
#'
#' @param x A `plasmode_dq_results` object from [run_plasmode_dq_stress()].
#' @param metric Character; metric to display. One of `"rmse"` (default),
#'   `"bias"`, `"coverage"`, `"emp_sd"`, or `"se_cal"`.
#' @param effect_size Numeric scalar; which effect size to display. `NULL`
#'   (default) uses the first available effect size.
#' @param scenarios Character vector of scenario names to include, or `NULL`
#'   (default) for all scenarios in the metrics table.
#' @param label_values Logical; print the metric value inside each tile.
#'   Default `TRUE`.
#' @param digits Integer; decimal places for tile labels. Default `3`.
#' @param low,high Character colour strings for the gradient endpoints.
#'   Defaults (`"#FFF5F0"` to `"#67000D"`) map low values to near-white and
#'   high values to dark red, which is appropriate for RMSE and bias (lower is
#'   better). For coverage, where higher is better, pass
#'   `low = "#67000D", high = "#009E73"` to reverse the direction.
#' @param ... Not currently used.
#'
#' @return A ggplot object.
#'
#' @seealso [run_plasmode_dq_stress()], [summarize_dq_degradation()]
#'
#' @examples
#' \dontrun{
#' dq <- run_plasmode_dq_stress(lock, tmle_candidates, reps = 20)
#' plot_dq_heatmap(dq)
#' plot_dq_heatmap(dq, metric = "coverage")
#' plot_dq_heatmap(dq, metric = "bias", scenarios = c("none", "unmeasured_U"))
#' }
#'
#' @export
plot_dq_heatmap <- function(x,
                             metric       = c("rmse", "bias", "coverage",
                                              "emp_sd", "se_cal"),
                             effect_size  = NULL,
                             scenarios    = NULL,
                             label_values = TRUE,
                             digits       = 3L,
                             low          = "#FFF5F0",
                             high         = "#67000D",
                             ...) {
  if (!inherits(x, "plasmode_dq_results"))
    stop("`x` must be a plasmode_dq_results object.", call. = FALSE)
  if (!requireNamespace("ggplot2", quietly = TRUE))
    stop("ggplot2 is required. Install it with install.packages('ggplot2').",
         call. = FALSE)

  metric <- match.arg(metric)
  m      <- x$metrics

  # Effect-size filter
  es_vals <- sort(unique(m$effect_size))
  if (is.null(effect_size)) {
    effect_size <- es_vals[[1L]]
    if (length(es_vals) > 1L)
      message("plot_dq_heatmap: multiple effect sizes (",
              paste(es_vals, collapse = ", "), "); using ", effect_size,
              ". Pass effect_size = ... to change.")
  }
  m <- m[m$effect_size == effect_size, , drop = FALSE]

  if (!is.null(scenarios))
    m <- m[m$scenario %in% scenarios, , drop = FALSE]

  if (nrow(m) == 0L)
    stop("No rows remain after filtering. Check effect_size / scenarios.",
         call. = FALSE)

  # Y-axis labels: "Baseline" for scenario == "none", else "scenario [level]"
  m$y_label <- ifelse(
    m$scenario == "none",
    "Baseline",
    ifelse(m$level == 0,
           m$scenario,
           paste0(m$scenario, " [", m$level, "]"))
  )

  # Row ordering: baseline first, then scenarios sorted by name and level
  base_rows  <- m[m$scenario == "none",  , drop = FALSE]
  other_rows <- m[m$scenario != "none",  , drop = FALSE]
  other_rows <- other_rows[order(other_rows$scenario, other_rows$level), ]
  m          <- rbind(base_rows, other_rows)
  # Reverse for ggplot so baseline appears at the top
  m$y_label  <- factor(m$y_label, levels = rev(unique(m$y_label)))

  # Candidate ordering follows the original object
  cand_order  <- unique(x$metrics$candidate)
  m$candidate <- factor(m$candidate, levels = cand_order)

  metric_label <- switch(metric,
    rmse     = "RMSE",
    bias     = "Bias",
    coverage = "Coverage",
    emp_sd   = "Empirical SD",
    se_cal   = "SE calibration (mean SE / emp SD)"
  )

  # For coverage, green = good (high), so invert the gradient
  if (metric == "coverage") {
    tmp  <- low
    low  <- high
    high <- tmp
  }

  p <- ggplot2::ggplot(
      m,
      ggplot2::aes(x = .data$candidate,
                   y = .data$y_label,
                   fill = .data[[metric]])
    ) +
    ggplot2::geom_tile(colour = "white", linewidth = 0.5) +
    ggplot2::scale_fill_gradient(low = low, high = high,
                                  name = metric_label, na.value = "grey80") +
    ggplot2::labs(
      x        = "Candidate",
      y        = "Scenario × severity",
      title    = paste0("DQ stress-test heatmap: ", metric_label),
      subtitle = paste0("Effect size = ", effect_size)
    ) +
    ggplot2::theme_minimal(base_size = 11) +
    ggplot2::theme(
      axis.text.x     = ggplot2::element_text(angle = 30, hjust = 1),
      panel.grid      = ggplot2::element_blank(),
      plot.title      = ggplot2::element_text(face = "bold"),
      legend.position = "right"
    )

  if (isTRUE(label_values)) {
    z         <- m[[metric]]
    m$.label  <- ifelse(is.na(z), "", format(round(z, digits), nsmall = digits))
    p <- p +
      ggplot2::geom_text(
        data        = m,
        mapping     = ggplot2::aes(x = .data$candidate,
                                   y = .data$y_label,
                                   label = .data$.label),
        colour      = "grey10",
        size        = 3,
        inherit.aes = FALSE
      )
  }

  p
}


#' Internal weight histogram plotter
#' @keywords internal
.hist_weights <- function(fit, type = "ps") {
  spec <- fit$spec
  if (is.null(spec$treatment)) {
    stop("No treatment variable specified.", call. = FALSE)
  }

  trt_var <- spec$treatment$name
  trt_vals <- spec$data[[trt_var]]

  if (type == "ps" && !is.null(fit$ps)) {
    plot_df <- data.frame(
      value = fit$ps,
      group = as.character(trt_vals),
      stringsAsFactors = FALSE
    )
    xlab <- "Propensity Score"
    title <- "Propensity Score Distribution by Treatment Group"
  } else {
    w <- if (!is.null(fit$weights$iptw)) fit$weights$iptw else fit$weights$combined
    plot_df <- data.frame(
      value = w,
      group = as.character(trt_vals),
      stringsAsFactors = FALSE
    )
    xlab <- "IPW Weight"
    title <- "Weight Distribution by Treatment Group"
  }

  p <- ggplot2::ggplot(plot_df, ggplot2::aes(
    x = .data$value, fill = .data$group
  )) +
    ggplot2::geom_histogram(bins = 30, alpha = 0.6,
                            position = "identity") +
    ggplot2::labs(x = xlab, y = "Count", fill = "Treatment",
                  title = title) +
    ggplot2::theme_minimal() +
    ggplot2::facet_wrap(~ .data$group, ncol = 1)

  p
}
