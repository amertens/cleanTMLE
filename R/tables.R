#' Create Table 1: Baseline Characteristics by Group
#'
#' Generates a summary table of baseline characteristics by treatment group,
#' optionally weighted. Includes standardized mean differences (SMD) when
#' requested.
#'
#' @param x A data.frame, `cr_spec`, or fitted result object. If a fitted
#'   object with weights, weighted summaries can be produced.
#' @param vars Character vector of variable names to summarize. If `NULL`,
#'   attempts to determine from the specification.
#' @param by Grouping variable name (character). If `NULL` and a treatment
#'   is specified in the spec, uses the treatment variable.
#' @param weighted Logical; if `TRUE` and `x` is a fitted object, use
#'   IPW-weighted summaries.
#' @param smd Logical; if `TRUE`, include standardized mean differences.
#' @param ... Additional arguments.
#'
#' @return A data.frame with summary statistics.
#'
#' @export
make_table1 <- function(x, ...) {
  UseMethod("make_table1")
}

#' @rdname make_table1
#' @export
make_table1.cleanroom_lock <- function(x, vars = NULL, by = NULL,
                                       weighted = FALSE, smd = TRUE, ...) {
  data <- x$data
  if (is.null(by)) by <- x$treatment
  if (is.null(vars)) vars <- x$covariates
  w <- rep(1, nrow(data))
  .build_table1(data = data, vars = vars, by = by, weights = w, smd = smd)
}

#' @rdname make_table1
#' @export
make_table1.default <- function(x, vars = NULL, by = NULL,
                                weighted = FALSE, smd = TRUE, ...) {
  if (is.data.frame(x)) {
    data <- x
    w <- rep(1, nrow(data))
  } else {
    stop("Unsupported input type for make_table1.", call. = FALSE)
  }

  .build_table1(data = data, vars = vars, by = by, weights = w, smd = smd)
}

#' @rdname make_table1
#' @export
make_table1.cumrisk <- function(x, vars = NULL, by = NULL,
                                weighted = FALSE, smd = TRUE, ...) {
  spec <- x$spec
  data <- spec$data

  if (is.null(by) && !is.null(spec$treatment)) {
    by <- spec$treatment$name
  }

  if (weighted && !is.null(x$weights$iptw)) {
    w <- x$weights$iptw
  } else {
    w <- rep(1, nrow(data))
  }

  .build_table1(data = data, vars = vars, by = by, weights = w, smd = smd)
}


#' Internal Table 1 builder
#' @keywords internal
.build_table1 <- function(data, vars = NULL, by = NULL, weights = NULL,
                           smd = TRUE) {
  if (is.null(weights)) weights <- rep(1, nrow(data))

  # Auto-detect numeric vars if not specified
  if (is.null(vars)) {
    vars <- names(data)[vapply(data, is.numeric, logical(1))]
    # Exclude obvious non-covariate columns
    exclude_patterns <- c("^id$", "^subject", "^time$", "^start$", "^stop$")
    for (pat in exclude_patterns) {
      vars <- vars[!grepl(pat, vars, ignore.case = TRUE)]
    }
  }

  if (is.null(by)) {
    # Overall summary
    rows <- lapply(vars, function(v) {
      vals <- data[[v]]
      w <- weights
      if (is.numeric(vals)) {
        data.frame(
          variable = v,
          group = "Overall",
          n = sum(!is.na(vals)),
          mean = weighted.mean(vals, w, na.rm = TRUE),
          sd = sqrt(weighted.mean((vals - weighted.mean(vals, w, na.rm = TRUE))^2,
                                  w, na.rm = TRUE)),
          median = median(vals, na.rm = TRUE),
          stringsAsFactors = FALSE
        )
      } else {
        lvls <- sort(unique(vals))
        do.call(rbind, lapply(lvls, function(l) {
          data.frame(
            variable = paste0(v, " = ", l),
            group = "Overall",
            n = sum(vals == l, na.rm = TRUE),
            mean = weighted.mean(vals == l, w, na.rm = TRUE),
            sd = NA_real_,
            median = NA_real_,
            stringsAsFactors = FALSE
          )
        }))
      }
    })
    result <- do.call(rbind, rows)
    return(result)
  }

  # By group
  groups <- sort(unique(data[[by]]))
  rows <- list()

  for (v in vars) {
    if (v == by) next
    vals <- data[[v]]

    if (is.numeric(vals)) {
      group_stats <- lapply(groups, function(g) {
        idx <- data[[by]] == g
        w_g <- weights[idx]
        v_g <- vals[idx]
        data.frame(
          variable = v,
          group = as.character(g),
          n = sum(!is.na(v_g)),
          mean = weighted.mean(v_g, w_g, na.rm = TRUE),
          sd = sqrt(weighted.mean(
            (v_g - weighted.mean(v_g, w_g, na.rm = TRUE))^2,
            w_g, na.rm = TRUE
          )),
          median = median(v_g, na.rm = TRUE),
          stringsAsFactors = FALSE
        )
      })
      row_df <- do.call(rbind, group_stats)

      # SMD (for binary groups)
      if (smd && length(groups) == 2) {
        m1 <- row_df$mean[1]
        m2 <- row_df$mean[2]
        s1 <- row_df$sd[1]
        s2 <- row_df$sd[2]
        pooled_sd <- sqrt((s1^2 + s2^2) / 2)
        smd_val <- if (pooled_sd > 0) (m2 - m1) / pooled_sd else 0
        row_df$smd <- smd_val
      }

      rows <- c(rows, list(row_df))
    }
  }

  result <- do.call(rbind, rows)
  rownames(result) <- NULL
  result
}


#' Create Table 2: Risk and Effect Estimate Summary
#'
#' Generates a results table with N, person-time, events, rates,
#' cumulative risk at specified time, and effect estimates (RD, RR, or HR).
#'
#' @param x A fitted result object (cumrisk, hr, or tmle_fit).
#' @param risk_time Time point for cumulative risk reporting (for cumrisk).
#' @param ... Additional arguments.
#'
#' @return A data.frame with summary results.
#'
#' @export
make_table2 <- function(x, ...) {
  UseMethod("make_table2")
}

#' @rdname make_table2
#' @export
make_table2.cumrisk <- function(x, risk_time = NULL, ...) {
  spec <- x$spec
  data <- spec$data
  time_var <- .get_time_var(spec)
  event_ind <- .build_event_indicator(data, spec)
  obs_time <- data[[time_var]]

  groups <- unique(x$risk$group)

  rows <- lapply(groups, function(g) {
    if (g == "overall") {
      idx <- rep(TRUE, nrow(data))
    } else if (!is.null(spec$treatment)) {
      idx <- data[[spec$treatment$name]] == g
    } else {
      idx <- rep(TRUE, nrow(data))
    }

    n_total <- sum(idx)
    n_events <- sum(event_ind[idx] == 1)
    pt <- sum(obs_time[idx])
    rate <- n_events / pt * 1000  # per 1000 person-time

    # Cumulative risk at risk_time
    if (!is.null(risk_time)) {
      risk_row <- x$risk[x$risk$group == g & x$risk$time == risk_time, ]
      if (nrow(risk_row) > 0) {
        cum_risk <- risk_row$risk[1]
        ci_l <- if ("ci_lower" %in% names(risk_row)) risk_row$ci_lower[1] else NA
        ci_u <- if ("ci_upper" %in% names(risk_row)) risk_row$ci_upper[1] else NA
      } else {
        cum_risk <- ci_l <- ci_u <- NA_real_
      }
    } else {
      # Use the last time
      sub <- x$risk[x$risk$group == g, ]
      last_row <- sub[nrow(sub), ]
      cum_risk <- last_row$risk
      ci_l <- if ("ci_lower" %in% names(last_row)) last_row$ci_lower else NA
      ci_u <- if ("ci_upper" %in% names(last_row)) last_row$ci_upper else NA
      risk_time <- last_row$time
    }

    data.frame(
      group = g,
      n = n_total,
      events = n_events,
      person_time = round(pt, 1),
      rate_per_1000 = round(rate, 2),
      risk_time = risk_time,
      cumulative_risk = round(cum_risk, 4),
      ci_lower = round(ci_l, 4),
      ci_upper = round(ci_u, 4),
      stringsAsFactors = FALSE
    )
  })

  result <- do.call(rbind, rows)

  # Add contrasts for two groups
  if (length(groups) == 2) {
    r1 <- result$cumulative_risk[1]
    r2 <- result$cumulative_risk[2]
    result$RD <- c(NA, round(r2 - r1, 4))
    result$RR <- c(NA, if (r1 > 0) round(r2 / r1, 4) else NA)
  }

  rownames(result) <- NULL
  result
}

#' @rdname make_table2
#' @export
make_table2.hr <- function(x, ...) {
  spec <- x$spec
  data <- spec$data
  time_var <- .get_time_var(spec)
  event_ind <- .build_event_indicator(data, spec)
  obs_time <- data[[time_var]]
  trt_var <- spec$treatment$name
  trt_levels <- sort(unique(data[[trt_var]]))

  rows <- lapply(trt_levels, function(g) {
    idx <- data[[trt_var]] == g
    n_total <- sum(idx)
    n_events <- sum(event_ind[idx] == 1)
    pt <- sum(obs_time[idx])

    data.frame(
      group = as.character(g),
      n = n_total,
      events = n_events,
      person_time = round(pt, 1),
      rate_per_1000 = round(n_events / pt * 1000, 2),
      stringsAsFactors = FALSE
    )
  })

  result <- do.call(rbind, rows)

  # Attach HR
  ht <- x$hr_table
  trt_term <- ht$term[grepl(trt_var, ht$term)]
  if (length(trt_term) > 0) {
    result$HR <- c(NA, ht$hr[ht$term == trt_term[1]])
    result$HR_ci_lower <- c(NA, ht$ci_lower[ht$term == trt_term[1]])
    result$HR_ci_upper <- c(NA, ht$ci_upper[ht$term == trt_term[1]])
    result$HR_p <- c(NA, ht$p_value[ht$term == trt_term[1]])
  }

  rownames(result) <- NULL
  result
}

#' @rdname make_table2
#' @export
make_table2.tmle_fit <- function(x, ...) {
  est <- x$estimates
  data.frame(
    estimand = names(est),
    estimate = vapply(est, function(e) e$estimate, numeric(1)),
    se = vapply(est, function(e) e$se, numeric(1)),
    ci_lower = vapply(est, function(e) e$ci_lower, numeric(1)),
    ci_upper = vapply(est, function(e) e$ci_upper, numeric(1)),
    p_value = vapply(est, function(e) e$p_value, numeric(1)),
    stringsAsFactors = FALSE
  )
}


#' Weight Summary Table
#'
#' Produces a summary table of weight distributions by treatment group,
#' including percentiles and summary statistics.
#'
#' @param x A fitted result object with weight information (cumrisk, hr,
#'   ipw, etc.).
#' @param ... Additional arguments.
#'
#' @return A data.frame with weight summaries.
#'
#' @export
make_wt_summary_table <- function(x, ...) {
  if (!inherits(x, "cr_result")) {
    stop("`x` must be a cleanTMLE result object.", call. = FALSE)
  }

  w <- if (!is.null(x$weights$iptw)) x$weights$iptw else
    if (!is.null(x$weights$combined)) x$weights$combined else
    stop("No weights found in the object.", call. = FALSE)

  spec <- x$spec
  if (!is.null(spec$treatment)) {
    trt_var <- spec$treatment$name
    trt_vals <- spec$data[[trt_var]]
    groups <- sort(unique(trt_vals))
  } else {
    trt_vals <- rep("overall", length(w))
    groups <- "overall"
  }

  rows <- lapply(groups, function(g) {
    idx <- trt_vals == g
    wg <- w[idx]
    data.frame(
      group = as.character(g),
      n = length(wg),
      mean = round(mean(wg), 4),
      sd = round(sd(wg), 4),
      min = round(min(wg), 4),
      p1 = round(quantile(wg, 0.01), 4),
      p5 = round(quantile(wg, 0.05), 4),
      p25 = round(quantile(wg, 0.25), 4),
      p50 = round(quantile(wg, 0.50), 4),
      p75 = round(quantile(wg, 0.75), 4),
      p95 = round(quantile(wg, 0.95), 4),
      p99 = round(quantile(wg, 0.99), 4),
      max = round(max(wg), 4),
      stringsAsFactors = FALSE,
      row.names = NULL
    )
  })

  result <- do.call(rbind, rows)
  rownames(result) <- NULL
  result
}


#' Identify Extreme Weights
#'
#' Returns the top K observations with the most extreme IPW weights,
#' along with key columns for inspection.
#'
#' @param x A fitted result object with weight information.
#' @param k Number of extreme observations to return (default 10).
#' @param cols Character vector of additional columns to include. If `NULL`,
#'   includes treatment, outcome, and time variables.
#' @param ... Additional arguments.
#'
#' @return A data.frame with the K most extreme-weighted observations.
#'
#' @export
extreme_weights <- function(x, k = 10, cols = NULL, ...) {
  if (!inherits(x, "cr_result")) {
    stop("`x` must be a cleanTMLE result object.", call. = FALSE)
  }

  w <- if (!is.null(x$weights$iptw)) x$weights$iptw else
    if (!is.null(x$weights$combined)) x$weights$combined else
    stop("No weights found.", call. = FALSE)

  data <- x$spec$data

  # Default columns
  if (is.null(cols)) {
    cols <- c()
    if (!is.null(x$spec$treatment)) cols <- c(cols, x$spec$treatment$name)
    if (!is.null(x$spec$outcome)) cols <- c(cols, x$spec$outcome$name)
    time_var <- .get_time_var(x$spec)
    if (time_var %in% names(data)) cols <- c(cols, time_var)
    if (!is.null(x$spec$subject)) cols <- c(cols, x$spec$subject$name)
  }

  cols <- unique(cols[cols %in% names(data)])

  top_idx <- order(abs(w), decreasing = TRUE)[seq_len(min(k, length(w)))]

  result <- data[top_idx, cols, drop = FALSE]
  result$weight <- w[top_idx]
  if (!is.null(x$ps)) result$ps <- x$ps[top_idx]
  result$row <- top_idx

  rownames(result) <- NULL
  result
}


#' Inspect IPW Weights
#'
#' Extract the IPW weights from a fitted result object for manual
#' inspection or downstream analysis.
#'
#' @param x A fitted result object.
#' @param type `"iptw"`, `"ipcw"`, or `"combined"`.
#' @param ... Additional arguments.
#'
#' @return A numeric vector of weights.
#'
#' @export
inspect_ipw_weights <- function(x, type = c("iptw", "ipcw", "combined"), ...) {
  if (!inherits(x, "cr_result")) {
    stop("`x` must be a cleanTMLE result object.", call. = FALSE)
  }
  type <- match.arg(type)

  w <- x$weights[[type]]
  if (is.null(w)) {
    stop(sprintf("Weight type '%s' not available.", type), call. = FALSE)
  }
  w
}
