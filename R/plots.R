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
  if (!inherits(x, "hr")) {
    stop("`x` must be an hr object.", call. = FALSE)
  }

  ht <- x$hr_table
  ht$term <- factor(ht$term, levels = rev(ht$term))

  p <- ggplot2::ggplot(ht, ggplot2::aes(
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

  p
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
