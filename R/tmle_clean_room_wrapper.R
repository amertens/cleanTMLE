#' Clean-Room TMLE Wrapper and Decision Log Utilities
#'
#' Provides \code{run_clean_tmle()}, a single-call wrapper that orchestrates
#' the full staged clean-room TMLE pipeline (Stages 1 -> 2a -> 2b -> 3),
#' along with helper functions for love plots, influence-curve histograms,
#' and a CSV-based decision log.
#'
#' @name tmle_clean_room_wrapper
NULL


# ── Decision Log Helpers ──────────────────────────────────────────────────

#' Initialise an Empty Decision Log
#'
#' Creates a zero-row data.frame with the standard decision-log schema.
#'
#' @return A data.frame with columns \code{stage}, \code{metric},
#'   \code{value}, \code{decision}, \code{rationale}, \code{timestamp}.
#'
#' @export
init_decision_log <- function() {
  data.frame(
    stage     = character(0),
    metric    = character(0),
    value     = character(0),
    decision  = character(0),
    rationale = character(0),
    timestamp = character(0),
    stringsAsFactors = FALSE
  )
}


#' Append a Row to a Decision Log
#'
#' @param log A decision-log data.frame from \code{\link{init_decision_log}}.
#' @param stage Character; stage label (e.g. \code{"Stage 1"}).
#' @param metric Character; name of the metric recorded.
#' @param value Character (or coercible); metric value.
#' @param decision Character; one of \code{"GO"}, \code{"FLAG"},
#'   \code{"STOP"}, or \code{NA}.
#' @param rationale Character; human-readable justification.
#'
#' @return The updated decision-log data.frame.
#'
#' @export
log_decision_entry <- function(log, stage, metric, value,
                               decision = NA_character_,
                               rationale = "") {
  new_row <- data.frame(
    stage     = stage,
    metric    = metric,
    value     = as.character(value),
    decision  = as.character(decision),
    rationale = rationale,
    timestamp = format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
    stringsAsFactors = FALSE
  )
  rbind(log, new_row)
}


#' Save a Decision Log to CSV
#'
#' Writes the decision log to the \code{inst/decision_logs/} directory
#' inside the package source tree, or to a user-specified path.
#'
#' @param log A decision-log data.frame.
#' @param path Character; file path.  Defaults to
#'   \code{"inst/decision_logs/decision_log.csv"} relative to the current
#'   working directory.
#' @param append Logical; if \code{TRUE} and the file already exists,
#'   append rows without re-writing the header.
#'
#' @return Invisibly returns \code{path}.
#'
#' @export
save_decision_log <- function(log, path = "inst/decision_logs/decision_log.csv",
                              append = FALSE) {
  dir_path <- dirname(path)
  if (!dir.exists(dir_path)) dir.create(dir_path, recursive = TRUE)

  if (append && file.exists(path)) {
    utils::write.table(log, file = path, sep = ",", row.names = FALSE,
                       col.names = FALSE, append = TRUE, quote = TRUE)
  } else {
    utils::write.csv(log, file = path, row.names = FALSE)
  }
  invisible(path)
}


# ── Visualisation Helpers ────────────────────────────────────────────────

#' Love Plot: Weighted vs Unweighted SMDs
#'
#' Produces a dot-plot showing standardised mean differences before and
#' after IPW weighting.  Requires a \code{ps_diagnostics} object from
#' \code{\link{compute_ps_diagnostics}}.
#'
#' @param ps_diag A \code{ps_diagnostics} object.
#' @param threshold Numeric; dashed reference line for the balance
#'   threshold. Default: 0.10.
#'
#' @return A \code{ggplot2} object.
#'
#' @export
love_plot <- function(ps_diag, threshold = 0.10) {
  if (!inherits(ps_diag, "ps_diagnostics"))
    stop("`ps_diag` must be a ps_diagnostics object.", call. = FALSE)

  smds <- ps_diag$smds

  plot_df <- data.frame(
    variable = rep(smds$variable, 2L),
    type     = rep(c("Unweighted", "Weighted"), each = nrow(smds)),
    smd      = c(abs(smds$smd_unweighted), abs(smds$smd_weighted)),
    stringsAsFactors = FALSE
  )

  ggplot2::ggplot(
    plot_df,
    ggplot2::aes(
      x     = .data$smd,
      y     = stats::reorder(.data$variable, .data$smd),
      shape = .data$type,
      colour = .data$type
    )
  ) +
    ggplot2::geom_point(size = 3) +
    ggplot2::geom_vline(xintercept = threshold, linetype = "dashed",
                        colour = "grey40") +
    ggplot2::labs(
      x      = "Absolute Standardised Mean Difference",
      y      = NULL,
      title  = "Love Plot: Covariate Balance",
      shape  = NULL,
      colour = NULL
    ) +
    ggplot2::theme_minimal() +
    ggplot2::theme(legend.position = "bottom")
}


#' Influence Curve Histogram
#'
#' Plots a histogram of the efficient influence curve (EIC) values from
#' a fitted TMLE object.
#'
#' @param tmle_result A \code{tmle_fit} object with an
#'   \code{influence_curve} element.
#' @param bins Integer; number of histogram bins. Default: 30.
#'
#' @return A \code{ggplot2} object.
#'
#' @export
ic_histogram <- function(tmle_result, bins = 30L) {
  ic <- tmle_result$influence_curve
  if (is.null(ic))
    stop("No influence_curve found in the TMLE result.", call. = FALSE)

  ic_df <- data.frame(ic = ic)

  ggplot2::ggplot(ic_df, ggplot2::aes(x = .data$ic)) +
    ggplot2::geom_histogram(bins = bins, fill = "steelblue", colour = "white",
                            alpha = 0.8) +
    ggplot2::geom_vline(xintercept = 0, linetype = "dashed", colour = "red") +
    ggplot2::labs(
      x     = "Influence Curve Value",
      y     = "Count",
      title = "Efficient Influence Curve Distribution"
    ) +
    ggplot2::theme_minimal()
}


# ── Main Wrapper ─────────────────────────────────────────────────────────

#' Run Clean-Room TMLE Pipeline
#'
#' A single-call wrapper that executes the staged clean-room TMLE
#' workflow: Stage 1 (specification lock), Stage 2a (propensity-score
#' diagnostics and gating), Stage 2b (negative-control or plasmode
#' simulation), and Stage 3 (primary TMLE estimation).
#'
#' @param data A data.frame containing the analysis dataset.
#' @param Avar Character; name of the binary treatment variable.
#' @param Yvar Character; name of the outcome variable.
#' @param covariates Character vector of baseline covariate names.
#'   If \code{NULL}, all numeric columns except \code{Avar}, \code{Yvar},
#'   \code{time_var}, and \code{cens_var} are used.
#' @param time_var Character or \code{NULL}; name of the time-to-event
#'   variable.  If provided, survival TMLE is used.
#' @param cens_var Character or \code{NULL}; name of the censoring
#'   indicator variable.
#' @param learner_lib Character vector of SuperLearner algorithms.
#'   **Required** — the function will stop if not provided.
#' @param truncation Numeric vector of length 2 giving the lower and
#'   upper PS truncation bounds.  Default: \code{c(0.01, 0.99)}.
#' @param neg_control_outcomes Character vector of negative-control
#'   outcome column names.  If provided, Stage 2b runs negative-control
#'   analyses.
#' @param simulation_spec A named list for plasmode simulation, with
#'   elements \code{reps} (integer; number of replicates) and optionally
#'   \code{effect_sizes} (numeric vector).  If provided, Stage 2b runs
#'   plasmode simulation instead of (or in addition to) negative controls.
#' @param stage2a_thresholds A named list of Stage 2a gating thresholds.
#'   Currently supports \code{ESS} (minimum ESS as a fraction of N;
#'   default 0.20) and \code{max_smd} (maximum tolerable weighted SMD;
#'   default 0.25).
#' @param output_dir Character; directory for saving the final RDS
#'   output.  If \code{NULL}, results are not written to disk.
#' @param seed Integer; random seed. Default: 42.
#' @param verbose Logical; if \code{TRUE}, print progress messages.
#'
#' @return A list of class \code{clean_tmle_result} containing:
#'   \describe{
#'     \item{decision_log}{The decision-log data.frame.}
#'     \item{lock}{The \code{cleanroom_lock} object.}
#'     \item{ps_fit}{The fitted propensity-score object.}
#'     \item{ps_diagnostics}{PS diagnostics (ESS, SMDs, overlap plot).}
#'     \item{weight_summary}{Weight distribution table.}
#'     \item{stage2a_decision}{Character: \code{"GO"}, \code{"FLAG"},
#'       or \code{"STOP"}.}
#'     \item{nc_results}{List of negative-control results (if any).}
#'     \item{plasmode_results}{Plasmode simulation results (if any).}
#'     \item{stage2b_decision}{Character: \code{"GO"}, \code{"FLAG"},
#'       or \code{"STOP"} (if Stage 2b was run).}
#'     \item{tmle_result}{The primary TMLE estimate (if Stage 3 ran).}
#'     \item{love_plot}{A ggplot2 love plot.}
#'   }
#'
#' @examples
#' \dontrun{
#' dat <- sim_func1(n = 500, seed = 1)
#' dat$event_24 <- as.integer(dat$event == 1 & dat$time <= 24)
#' res <- run_clean_tmle(
#'   data        = dat,
#'   Avar        = "treatment",
#'   Yvar        = "event_24",
#'   covariates  = c("age", "sex", "biomarker"),
#'   learner_lib = c("SL.glm", "SL.mean"),
#'   seed        = 1
#' )
#' print(res$decision_log)
#' res$love_plot
#' }
#'
#' @export
run_clean_tmle <- function(data,
                           Avar,
                           Yvar,
                           covariates          = NULL,
                           time_var            = NULL,
                           cens_var            = NULL,
                           learner_lib,
                           truncation          = c(0.01, 0.99),
                           neg_control_outcomes = NULL,
                           simulation_spec     = NULL,
                           stage2a_thresholds  = list(ESS = 0.20,
                                                      max_smd = 0.25),
                           output_dir          = NULL,
                           seed                = 42L,
                           verbose             = TRUE) {

  # ── Stage 1: Specification Lock ──────────────────────────────────────

  if (missing(learner_lib))
    stop("A pre-specified SuperLearner library (`learner_lib`) is required.",
         call. = FALSE)
  if (missing(truncation) || !is.numeric(truncation) || length(truncation) != 2L)
    stop("`truncation` must be a numeric vector of length 2 (e.g., c(0.01, 0.99)).",
         call. = FALSE)

  dlog <- init_decision_log()

  # Auto-detect covariates
  if (is.null(covariates)) {
    exclude <- c(Avar, Yvar, time_var, cens_var, "id", "subject")
    covariates <- setdiff(
      names(data)[vapply(data, is.numeric, logical(1))],
      exclude
    )
  }

  trunc_lower <- truncation[1]

  lock <- create_analysis_lock(
    data       = data,
    treatment  = Avar,
    outcome    = Yvar,
    covariates = covariates,
    sl_library = learner_lib,
    seed       = as.integer(seed)
  )

  dlog <- log_decision_entry(dlog, "Stage 1", "SL library",
                             paste(learner_lib, collapse = ", "),
                             decision = "GO",
                             rationale = "Pre-specified learner library recorded.")
  dlog <- log_decision_entry(dlog, "Stage 1", "Truncation",
                             paste(truncation, collapse = ", "),
                             decision = "GO",
                             rationale = "Pre-specified truncation bounds recorded.")

  if (verbose) message("Stage 1: Analysis lock created.")

  # ── Stage 2a: Propensity Score Diagnostics ──────────────────────────

  ps_fit <- fit_ps_superlearner(lock)
  ps_diag <- compute_ps_diagnostics(ps_fit)

  # IPW weights (unstabilised)
  A <- data[[Avar]]
  ps <- ps_fit$ps
  ps_trunc <- pmax(pmin(ps, truncation[2]), truncation[1])
  w <- ifelse(A == 1, 1 / ps_trunc, 1 / (1 - ps_trunc))

  # Weight summary table (manual, since make_wt_summary_table needs cr_result)
  wt_summary <- .build_weight_summary(w, A)

  # Extreme weights
  extreme_idx <- order(abs(w), decreasing = TRUE)[seq_len(min(10L, length(w)))]
  extreme_wt <- data.frame(
    row    = extreme_idx,
    weight = w[extreme_idx],
    ps     = ps_trunc[extreme_idx],
    A      = A[extreme_idx],
    stringsAsFactors = FALSE
  )

  # PS histogram (already in ps_diag$overlap_plot)

  # ESS (overall)
  ess <- sum(w)^2 / sum(w^2)
  n   <- nrow(data)
  ess_frac <- ess / n

  # Log PS range and fraction at boundaries
  frac_at_lower <- mean(ps <= truncation[1])
  frac_at_upper <- mean(ps >= truncation[2])

  dlog <- log_decision_entry(dlog, "Stage 2a", "PS range",
                             sprintf("[%.4f, %.4f]", min(ps), max(ps)))
  dlog <- log_decision_entry(dlog, "Stage 2a", "Fraction at lower truncation",
                             sprintf("%.4f", frac_at_lower))
  dlog <- log_decision_entry(dlog, "Stage 2a", "Fraction at upper truncation",
                             sprintf("%.4f", frac_at_upper))
  dlog <- log_decision_entry(dlog, "Stage 2a", "ESS",
                             sprintf("%.1f (%.1f%% of N=%d)", ess,
                                     100 * ess_frac, n))

  # Gating
  ess_threshold <- stage2a_thresholds$ESS %||% 0.20
  smd_threshold <- stage2a_thresholds$max_smd %||% 0.25

  max_weighted_smd <- max(abs(ps_diag$smds$smd_weighted))

  stage2a_flags <- character(0)
  if (ess_frac < ess_threshold)
    stage2a_flags <- c(stage2a_flags,
                       sprintf("ESS=%.1f (%.1f%%) below threshold %.0f%%",
                               ess, 100 * ess_frac, 100 * ess_threshold))
  if (max_weighted_smd > smd_threshold)
    stage2a_flags <- c(stage2a_flags,
                       sprintf("max weighted SMD=%.3f above threshold %.2f",
                               max_weighted_smd, smd_threshold))

  if (length(stage2a_flags) == 0L) {
    stage2a_decision <- "GO"
    stage2a_rationale <- "ESS and balance adequate."
  } else if (ess_frac < ess_threshold * 0.5) {
    stage2a_decision <- "STOP"
    stage2a_rationale <- paste("Critical:", paste(stage2a_flags, collapse = "; "))
  } else {
    stage2a_decision <- "FLAG"
    stage2a_rationale <- paste("Issues:", paste(stage2a_flags, collapse = "; "))
  }

  dlog <- log_decision_entry(dlog, "Stage 2a", "Gating decision",
                             stage2a_decision,
                             decision  = stage2a_decision,
                             rationale = stage2a_rationale)

  if (verbose) message("Stage 2a: ", stage2a_decision, " - ", stage2a_rationale)

  # Build love plot
  lp <- love_plot(ps_diag, threshold = smd_threshold)

  # Early exit on STOP
  if (stage2a_decision == "STOP") {
    result <- list(
      decision_log    = dlog,
      lock            = lock,
      ps_fit          = ps_fit,
      ps_diagnostics  = ps_diag,
      weight_summary  = wt_summary,
      extreme_weights = extreme_wt,
      stage2a_decision = stage2a_decision,
      nc_results       = NULL,
      plasmode_results = NULL,
      stage2b_decision = NA_character_,
      tmle_result      = NULL,
      love_plot        = lp
    )
    class(result) <- "clean_tmle_result"
    if (verbose) message("Pipeline halted at Stage 2a (STOP).")
    return(result)
  }

  # ── Stage 2b: Negative Controls / Plasmode Simulation ───────────────

  nc_results <- NULL
  plasmode_res <- NULL
  stage2b_decision <- "GO"
  stage2b_rationale <- "No Stage 2b checks requested."


  # Negative control outcomes
  if (!is.null(neg_control_outcomes)) {
    nc_results <- list()
    for (nc_var in neg_control_outcomes) {
      if (!nc_var %in% names(data)) {
        warning("Negative control variable '", nc_var,
                "' not found in data; skipping.", call. = FALSE)
        next
      }
      lock <- define_negative_control(lock, nc_var)
      nc_results[[nc_var]] <- run_negative_control(lock, nc_var, ps_fit)
    }

    if (length(nc_results) > 0L) {
      cp3 <- checkpoint_residual_bias(nc_results, lock_hash = lock$lock_hash)
      stage2b_decision <- cp3$decision
      stage2b_rationale <- cp3$rationale

      for (nc in nc_results) {
        dlog <- log_decision_entry(dlog, "Stage 2b (NC)", nc$variable,
                                   sprintf("est=%.5f, p=%.4f",
                                           nc$estimate, nc$p_value),
                                   rationale = nc$interpretation)
      }
      dlog <- log_decision_entry(dlog, "Stage 2b", "NC gating decision",
                                 stage2b_decision,
                                 decision  = stage2b_decision,
                                 rationale = stage2b_rationale)
    }
    if (verbose) message("Stage 2b (NC): ", stage2b_decision,
                         " - ", stage2b_rationale)
  }

  # Plasmode simulation
  if (!is.null(simulation_spec)) {
    sim_reps <- simulation_spec$reps %||% lock$plasmode_reps
    sim_es   <- simulation_spec$effect_sizes %||% c(0.05, 0.10)

    plasmode_res <- run_plasmode_feasibility(
      lock        = lock,
      effect_sizes = sim_es,
      reps         = sim_reps,
      verbose      = verbose
    )

    best_cand <- select_tmle_candidate(
      plasmode_res,
      rule = "min_rmse",
      thresholds = simulation_spec$thresholds
    )

    # Record
    bm <- best_cand$metrics
    dlog <- log_decision_entry(dlog, "Stage 2b (Plasmode)", "bias",
                               sprintf("%.5f", bm$bias))
    dlog <- log_decision_entry(dlog, "Stage 2b (Plasmode)", "rmse",
                               sprintf("%.5f", bm$rmse))
    dlog <- log_decision_entry(dlog, "Stage 2b (Plasmode)", "coverage",
                               sprintf("%.3f", bm$coverage))

    # Simple gate: coverage > 0.85 and bias < 0.05
    plas_gate <- bm$coverage >= 0.85 && abs(bm$bias) < 0.05
    plas_decision <- if (plas_gate) "GO" else "FLAG"

    dlog <- log_decision_entry(dlog, "Stage 2b (Plasmode)", "Gating decision",
                               plas_decision,
                               decision  = plas_decision,
                               rationale = sprintf(
                                 "Selected candidate '%s': bias=%.5f, rmse=%.5f, coverage=%.3f",
                                 best_cand$candidate_id, bm$bias, bm$rmse, bm$coverage))

    # Combine with NC decision if both ran
    if (plas_decision == "FLAG" || stage2b_decision == "FLAG")
      stage2b_decision <- "FLAG"
    if (plas_decision == "STOP" || stage2b_decision == "STOP")
      stage2b_decision <- "STOP"

    if (verbose) message("Stage 2b (Plasmode): ", plas_decision)
  }

  # Stage 2b STOP check
  if (identical(stage2b_decision, "STOP")) {
    result <- list(
      decision_log     = dlog,
      lock             = lock,
      ps_fit           = ps_fit,
      ps_diagnostics   = ps_diag,
      weight_summary   = wt_summary,
      extreme_weights  = extreme_wt,
      stage2a_decision = stage2a_decision,
      nc_results       = nc_results,
      plasmode_results = plasmode_res,
      stage2b_decision = stage2b_decision,
      tmle_result      = NULL,
      love_plot        = lp
    )
    class(result) <- "clean_tmle_result"
    if (verbose) message("Pipeline halted at Stage 2b (STOP).")
    return(result)
  }

  # ── Stage 3: Primary TMLE Estimation ────────────────────────────────

  dlog <- log_decision_entry(dlog, "Stage 3", "Outcome unblinding",
                             "TRUE",
                             decision  = "GO",
                             rationale = "Stages 2a/2b passed; unblinding outcome.")

  if (!is.null(time_var)) {
    # Survival TMLE
    tmle_result <- estimate_surv_tmle(
      data        = data,
      treatment   = Avar,
      time        = time_var,
      event       = Yvar,
      covariates  = covariates,
      sl_library  = learner_lib
    )
  } else {
    # Point-treatment binary TMLE
    tmle_result <- estimate_tmle_risk_point(
      data        = data,
      treatment   = Avar,
      outcome     = Yvar,
      covariates  = covariates,
      sl_library  = learner_lib,
      truncate    = trunc_lower
    )
  }

  # Log the estimate
  if (!is.null(tmle_result$estimates$ATE)) {
    ate <- tmle_result$estimates$ATE
    dlog <- log_decision_entry(dlog, "Stage 3", "ATE estimate",
                               sprintf("%.5f (95%% CI: %.5f, %.5f)",
                                       ate$estimate, ate$ci_lower, ate$ci_upper))
  }

  # Clever covariate summary (from tmle object)
  cc_summary <- NULL
  if (!is.null(tmle_result$tmle_obj)) {
    cc_summary <- tryCatch({
      g_hat <- tmle_result$tmle_obj$g$g1W
      H_aw  <- ifelse(A == 1, 1 / g_hat, -1 / (1 - g_hat))
      data.frame(
        statistic = c("mean", "sd", "min", "max"),
        value     = c(mean(H_aw), sd(H_aw), min(H_aw), max(H_aw)),
        stringsAsFactors = FALSE
      )
    }, error = function(e) NULL)
  }

  # IC summary
  ic_summary <- NULL
  if (!is.null(tmle_result$influence_curve)) {
    ic <- tmle_result$influence_curve
    ic_summary <- data.frame(
      statistic = c("mean", "sd", "min", "p5", "p25", "median",
                    "p75", "p95", "max"),
      value     = c(mean(ic), sd(ic), min(ic),
                    quantile(ic, 0.05), quantile(ic, 0.25), median(ic),
                    quantile(ic, 0.75), quantile(ic, 0.95), max(ic)),
      stringsAsFactors = FALSE
    )
  }

  dlog <- log_decision_entry(dlog, "Stage 3", "Pipeline complete",
                             "TRUE", decision = "GO",
                             rationale = "Primary TMLE estimate produced.")

  if (verbose) message("Stage 3: Primary TMLE estimation complete.")

  # ── Save outputs ───────────────────────────────────────────────────

  result <- list(
    decision_log     = dlog,
    lock             = lock,
    ps_fit           = ps_fit,
    ps_diagnostics   = ps_diag,
    weight_summary   = wt_summary,
    extreme_weights  = extreme_wt,
    stage2a_decision = stage2a_decision,
    nc_results       = nc_results,
    plasmode_results = plasmode_res,
    stage2b_decision = stage2b_decision,
    tmle_result      = tmle_result,
    love_plot        = lp,
    cc_summary       = cc_summary,
    ic_summary       = ic_summary
  )
  class(result) <- "clean_tmle_result"

  if (!is.null(output_dir)) {
    if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
    rds_path <- file.path(output_dir, paste0("clean_tmle_result_",
                                              format(Sys.time(), "%Y%m%d_%H%M%S"),
                                              ".rds"))
    saveRDS(result, rds_path)
    if (verbose) message("Results saved to: ", rds_path)

    csv_path <- file.path(output_dir, "decision_log.csv")
    save_decision_log(dlog, path = csv_path)
    if (verbose) message("Decision log saved to: ", csv_path)
  }

  result
}


#' @export
print.clean_tmle_result <- function(x, ...) {
  cat("Clean-Room TMLE Result\n")
  cat("======================\n\n")

  cat("Stage 2a decision: ", x$stage2a_decision, "\n")
  if (!is.na(x$stage2b_decision))
    cat("Stage 2b decision: ", x$stage2b_decision, "\n")

  if (!is.null(x$tmle_result) && !is.null(x$tmle_result$estimates$ATE)) {
    ate <- x$tmle_result$estimates$ATE
    cat(sprintf("\nPrimary ATE:  %.5f  (95%% CI: %.5f, %.5f)  p=%.4f\n",
                ate$estimate, ate$ci_lower, ate$ci_upper, ate$p_value))
  } else {
    cat("\nNo primary estimate (pipeline may have stopped early).\n")
  }

  cat(sprintf("\nDecision log: %d entries\n", nrow(x$decision_log)))
  invisible(x)
}


# ── Internal helpers ──────────────────────────────────────────────────────

#' Build weight summary table from raw vectors
#' @keywords internal
.build_weight_summary <- function(w, A) {
  groups <- sort(unique(A))
  rows <- lapply(groups, function(g) {
    wg <- w[A == g]
    data.frame(
      group = as.character(g),
      n     = length(wg),
      mean  = round(mean(wg), 4),
      sd    = round(sd(wg), 4),
      min   = round(min(wg), 4),
      p1    = round(quantile(wg, 0.01, names = FALSE), 4),
      p5    = round(quantile(wg, 0.05, names = FALSE), 4),
      p25   = round(quantile(wg, 0.25, names = FALSE), 4),
      p50   = round(quantile(wg, 0.50, names = FALSE), 4),
      p75   = round(quantile(wg, 0.75, names = FALSE), 4),
      p95   = round(quantile(wg, 0.95, names = FALSE), 4),
      p99   = round(quantile(wg, 0.99, names = FALSE), 4),
      max   = round(max(wg), 4),
      stringsAsFactors = FALSE,
      row.names = NULL
    )
  })
  do.call(rbind, rows)
}


# Null-default operator
`%||%` <- function(a, b) if (is.null(a)) b else a
