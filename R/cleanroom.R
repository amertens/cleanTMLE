#' Analysis Lock and Staged Workflow Functions
#'
#' These functions implement the staged clean-room workflow: locking the
#' analytic specification, estimating propensity scores with SuperLearner,
#' running plasmode feasibility evaluation, and executing the final
#' conventional (matching, IPTW) and modular TMLE workflows.
#'
#' @name cleanroom_workflow
NULL


# ── Stage 1: Analysis Lock ────────────────────────────────────────────────

#' Create an Analysis Lock
#'
#' Captures the full analytic specification—data, variable names, SuperLearner
#' library, and plasmode simulation settings—as a locked object. The lock
#' is created before the outcome is examined to enforce the pre-specified
#' analytic plan across all subsequent stages.
#'
#' @param data A data.frame containing covariates, treatment, and outcome.
#' @param treatment Character; name of the binary treatment column.
#' @param outcome Character; name of the outcome column.
#' @param covariates Character vector of baseline covariate column names.
#' @param sl_library Character vector of SuperLearner algorithm names.
#'   Default: `c("SL.glm", "SL.mean")`.
#' @param plasmode_reps Integer; number of plasmode replicates for Stage 2b
#'   feasibility evaluation. Default: 100.
#' @param seed Integer; random seed for reproducibility. Default: 42.
#'
#' @return An object of class `cleanroom_lock` containing all specified
#'   analysis parameters plus a reproducibility fingerprint (`lock_hash`).
#'
#' @examples
#' dat <- sim_func1(n = 200, seed = 1)
#' lock <- create_analysis_lock(
#'   data       = dat,
#'   treatment  = "treatment",
#'   outcome    = "event_24",
#'   covariates = c("age", "sex", "biomarker"),
#'   seed       = 1
#' )
#' print(lock)
#'
#' @export
create_analysis_lock <- function(data, treatment, outcome, covariates,
                                  sl_library    = c("SL.glm", "SL.mean"),
                                  plasmode_reps = 100L,
                                  seed          = 42L) {
  if (!is.data.frame(data))
    stop("`data` must be a data.frame.", call. = FALSE)
  if (!is.character(treatment) || length(treatment) != 1L)
    stop("`treatment` must be a single character string.", call. = FALSE)
  if (!is.character(outcome) || length(outcome) != 1L)
    stop("`outcome` must be a single character string.", call. = FALSE)
  if (!treatment %in% names(data))
    stop("treatment variable '", treatment, "' not found in data.",
         call. = FALSE)
  if (!outcome %in% names(data))
    stop("outcome variable '", outcome, "' not found in data.",
         call. = FALSE)
  missing_cov <- covariates[!covariates %in% names(data)]
  if (length(missing_cov) > 0L)
    stop("covariates not found in data: ",
         paste(missing_cov, collapse = ", "), call. = FALSE)

  lock_hash <- .compute_lock_hash(list(
    treatment  = treatment,
    outcome    = outcome,
    covariates = covariates,
    sl_library = sl_library,
    seed       = as.integer(seed),
    data_nrow  = nrow(data),
    data_ncol  = ncol(data),
    data_names = paste(sort(names(data)), collapse = "|")
  ))

  lock <- list(
    data          = data,
    treatment     = treatment,
    outcome       = outcome,
    covariates    = covariates,
    sl_library    = sl_library,
    plasmode_reps = as.integer(plasmode_reps),
    seed          = as.integer(seed),
    locked_at     = Sys.time(),
    lock_hash     = lock_hash
  )
  class(lock) <- "cleanroom_lock"
  lock
}

#' @keywords internal
.compute_lock_hash <- function(params) {
  s     <- paste(unlist(lapply(params, as.character)), collapse = "|")
  chars <- utf8ToInt(s)
  val   <- sum(chars * seq_along(chars)) %% 1e9
  sprintf("%09.0f", val)
}


#' Validate an Analysis Lock
#'
#' Checks that a `cleanroom_lock` object is complete, internally consistent,
#' and has not been modified after creation (hash check).
#'
#' @param lock A `cleanroom_lock` object from [create_analysis_lock()].
#'
#' @return Invisibly returns `lock` if valid; otherwise throws an error.
#'
#' @examples
#' dat <- sim_func1(n = 200, seed = 1)
#' lock <- create_analysis_lock(
#'   data       = dat,
#'   treatment  = "treatment",
#'   outcome    = "event_24",
#'   covariates = c("age", "sex", "biomarker"),
#'   seed       = 1
#' )
#' validate_analysis_lock(lock)
#'
#' @export
validate_analysis_lock <- function(lock) {
  if (!inherits(lock, "cleanroom_lock"))
    stop("`lock` must be a cleanroom_lock object.", call. = FALSE)

  required <- c("data", "treatment", "outcome", "covariates",
                "sl_library", "plasmode_reps", "seed", "lock_hash")
  missing_fields <- setdiff(required, names(lock))
  if (length(missing_fields) > 0L)
    stop("Lock is missing required fields: ",
         paste(missing_fields, collapse = ", "), call. = FALSE)

  # Recompute and compare hash
  computed_hash <- .compute_lock_hash(list(
    treatment  = lock$treatment,
    outcome    = lock$outcome,
    covariates = lock$covariates,
    sl_library = lock$sl_library,
    seed       = as.integer(lock$seed),
    data_nrow  = nrow(lock$data),
    data_ncol  = ncol(lock$data),
    data_names = paste(sort(names(lock$data)), collapse = "|")
  ))
  if (!identical(lock$lock_hash, computed_hash))
    stop("Lock hash mismatch: the lock object may have been modified.",
         call. = FALSE)

  # Variable presence
  if (!lock$treatment %in% names(lock$data))
    stop("treatment variable '", lock$treatment,
         "' not found in locked data.", call. = FALSE)
  if (!lock$outcome %in% names(lock$data))
    stop("outcome variable '", lock$outcome,
         "' not found in locked data.", call. = FALSE)
  missing_cov <- lock$covariates[!lock$covariates %in% names(lock$data)]
  if (length(missing_cov) > 0L)
    stop("covariates not found in locked data: ",
         paste(missing_cov, collapse = ", "), call. = FALSE)

  message("Analysis lock validated successfully.")
  invisible(lock)
}


#' @export
print.cleanroom_lock <- function(x, ...) {
  cat("cleanTMLE Analysis Lock\n")
  cat("=======================\n")
  cat("Data:       ", nrow(x$data), "observations,", ncol(x$data),
      "variables\n")
  cat("Treatment:  ", x$treatment, "\n")
  cat("Outcome:    ", x$outcome, "\n")
  cat("Covariates: ", paste(x$covariates, collapse = ", "), "\n")
  cat("SL library: ", paste(x$sl_library, collapse = ", "), "\n")
  cat("Plasmode:   ", x$plasmode_reps, "replicates\n")
  cat("Seed:       ", x$seed, "\n")
  if (!is.null(x$locked_at))
    cat("Locked at:  ", format(x$locked_at, "%Y-%m-%d %H:%M:%S"), "\n")
  cat("Hash:       ", x$lock_hash, "\n")
  invisible(x)
}


# ── Stage 2a: Propensity Score ────────────────────────────────────────────

#' Fit Propensity Score Using SuperLearner
#'
#' Estimates the treatment propensity score P(A=1|W) using SuperLearner
#' ensemble learning. Only covariates (W) and treatment (A) are used;
#' the outcome is never accessed, so this step is safe at Stage 2a.
#'
#' @param lock A `cleanroom_lock` from [create_analysis_lock()].
#' @param ... Additional arguments passed to `SuperLearner::SuperLearner()`.
#'
#' @return An object of class `ps_fit` containing propensity scores (`ps`),
#'   the fitted SuperLearner object (`sl_fit`), and metadata.
#'
#' @details
#' Requires the \pkg{SuperLearner} package.
#'
#' @examples
#' \dontrun{
#' dat  <- sim_func1(n = 500, seed = 1)
#' lock <- create_analysis_lock(dat, "treatment", "event_24",
#'                              c("age", "sex", "biomarker"), seed = 1)
#' ps_fit <- fit_ps_superlearner(lock)
#' }
#'
#' @export
fit_ps_superlearner <- function(lock, ...) {
  if (!inherits(lock, "cleanroom_lock"))
    stop("`lock` must be a cleanroom_lock object.", call. = FALSE)
  if (!requireNamespace("SuperLearner", quietly = TRUE))
    stop("Package 'SuperLearner' is required. ",
         "Install with: install.packages('SuperLearner')", call. = FALSE)

  data <- lock$data
  A    <- data[[lock$treatment]]
  W    <- data[, lock$covariates, drop = FALSE]

  set.seed(lock$seed)
  sl_fit <- SuperLearner::SuperLearner(
    Y          = A,
    X          = W,
    family     = binomial(),
    SL.library = lock$sl_library,
    ...
  )

  ps <- as.numeric(sl_fit$SL.predict)

  result <- list(
    ps         = ps,
    sl_fit     = sl_fit,
    treatment  = lock$treatment,
    covariates = lock$covariates,
    data       = data,
    sl_library = lock$sl_library,
    call       = match.call()
  )
  class(result) <- c("ps_fit", "cr_result")
  result
}


#' @export
print.ps_fit <- function(x, ...) {
  cat("Propensity Score Fit (SuperLearner)\n")
  cat("====================================\n")
  cat("Treatment:  ", x$treatment, "\n")
  cat("Covariates: ", paste(x$covariates, collapse = ", "), "\n")
  cat(sprintf("PS range:    [%.4f, %.4f]\n", min(x$ps), max(x$ps)))
  cat(sprintf("PS mean:      %.4f\n", mean(x$ps)))
  invisible(x)
}


#' Compute Propensity Score Diagnostics
#'
#' Computes overlap diagnostics, effective sample size (ESS), and
#' standardized mean differences (SMDs) from a fitted propensity score object.
#'
#' @param ps_fit A `ps_fit` object from [fit_ps_superlearner()].
#' @param ... Currently unused.
#'
#' @return An object of class `ps_diagnostics` containing:
#'   * `ess` - effective sample size table
#'   * `smds` - standardized mean differences before and after weighting
#'   * `overlap_plot` - a ggplot2 histogram of PS by treatment group
#'
#' @examples
#' \dontrun{
#' diag <- compute_ps_diagnostics(ps_fit)
#' print(diag)
#' plot(diag)
#' }
#'
#' @export
compute_ps_diagnostics <- function(ps_fit, ...) {
  if (!inherits(ps_fit, "ps_fit"))
    stop("`ps_fit` must be a ps_fit object from fit_ps_superlearner().",
         call. = FALSE)

  data <- ps_fit$data
  A    <- data[[ps_fit$treatment]]
  ps   <- ps_fit$ps

  # IPTW weights (unstabilized) for ESS and weighted SMDs
  w <- ifelse(A == 1, 1 / ps, 1 / (1 - ps))

  # Effective sample size (Kish)
  ess_treated <- sum(w[A == 1])^2 / sum(w[A == 1]^2)
  ess_control <- sum(w[A == 0])^2 / sum(w[A == 0]^2)
  n_treated   <- sum(A == 1)
  n_control   <- sum(A == 0)

  # Standardized mean differences
  smd_rows <- lapply(ps_fit$covariates, function(v) {
    x_val <- data[[v]]
    if (!is.numeric(x_val)) x_val <- as.numeric(as.factor(x_val)) - 1L

    m1 <- mean(x_val[A == 1])
    m0 <- mean(x_val[A == 0])
    s1 <- var(x_val[A == 1])
    s0 <- var(x_val[A == 0])
    pooled_sd <- sqrt((s1 + s0) / 2)
    smd_unw <- if (pooled_sd > 0) (m1 - m0) / pooled_sd else 0

    m1w <- weighted.mean(x_val[A == 1], w[A == 1])
    m0w <- weighted.mean(x_val[A == 0], w[A == 0])
    wv1 <- .weighted_var_cl(x_val[A == 1], w[A == 1])
    wv0 <- .weighted_var_cl(x_val[A == 0], w[A == 0])
    pooled_w <- sqrt((wv1 + wv0) / 2)
    smd_w <- if (pooled_w > 0) (m1w - m0w) / pooled_w else 0

    data.frame(
      variable       = v,
      smd_unweighted = round(smd_unw, 4),
      smd_weighted   = round(smd_w,   4),
      stringsAsFactors = FALSE
    )
  })
  smd_table <- do.call(rbind, smd_rows)

  # Overlap histogram
  plot_data <- data.frame(
    ps    = c(ps[A == 1],            ps[A == 0]),
    group = c(rep("Treated (A=1)", sum(A == 1)),
              rep("Control (A=0)", sum(A == 0))),
    stringsAsFactors = FALSE
  )
  overlap_plot <- ggplot2::ggplot(
    plot_data,
    ggplot2::aes(x = .data$ps, fill = .data$group)
  ) +
    ggplot2::geom_histogram(position = "identity", alpha = 0.5, bins = 30L) +
    ggplot2::labs(
      x     = "Propensity Score",
      y     = "Count",
      title = "Propensity Score Overlap",
      fill  = "Group"
    ) +
    ggplot2::theme_minimal()

  result <- list(
    ess = data.frame(
      group   = c("Treated", "Control", "Total"),
      n       = c(n_treated, n_control, n_treated + n_control),
      ess     = round(c(ess_treated, ess_control,
                        ess_treated + ess_control), 1L),
      ess_pct = round(100 * c(ess_treated / n_treated,
                               ess_control / n_control,
                               (ess_treated + ess_control) /
                                 (n_treated + n_control)), 1L),
      stringsAsFactors = FALSE
    ),
    smds         = smd_table,
    overlap_plot = overlap_plot,
    ps_summary   = summary(ps),
    call         = match.call()
  )
  class(result) <- "ps_diagnostics"
  result
}


#' @keywords internal
.weighted_var_cl <- function(x, w) {
  w   <- w / sum(w)
  xbar <- sum(w * x)
  sum(w * (x - xbar)^2)
}


#' @export
print.ps_diagnostics <- function(x, ...) {
  cat("Propensity Score Diagnostics\n")
  cat("============================\n\n")
  cat("Effective Sample Size (Kish ESS):\n")
  print(x$ess, row.names = FALSE)
  cat("\nStandardized Mean Differences:\n")
  print(x$smds, row.names = FALSE)
  cat("\nOverlap plot available via plot(diag).\n")
  invisible(x)
}


#' @export
plot.ps_diagnostics <- function(x, ...) {
  x$overlap_plot
}


# ── Stage 2b: Plasmode Feasibility ────────────────────────────────────────

#' Run Plasmode-Simulation Feasibility Evaluation
#'
#' Evaluates the performance of candidate analysis workflows (IPTW, TMLE)
#' using plasmode-simulation: synthetic binary outcomes are generated from
#' a parametric baseline-risk model fit on the real covariates, augmented
#' with a specified additive treatment effect. Each workflow is run on every
#' replicate and performance metrics (bias, RMSE, coverage) are computed
#' against the known true effect.
#'
#' @param lock A `cleanroom_lock` from [create_analysis_lock()].
#' @param effect_sizes Numeric vector of true risk differences to simulate.
#'   Default: `c(0.05, 0.10)`.
#' @param reps Integer; number of plasmode replicates per effect size.
#'   Defaults to `lock$plasmode_reps`.
#' @param verbose Logical; if `TRUE`, print progress messages. Default: `FALSE`.
#'
#' @return An object of class `plasmode_results` containing:
#'   * `metrics` - data.frame with bias, RMSE, and coverage by method and effect size
#'   * `results` - raw per-replicate estimates (list)
#'
#' @export
run_plasmode_feasibility <- function(lock,
                                      effect_sizes = c(0.05, 0.10),
                                      reps         = lock$plasmode_reps,
                                      verbose      = FALSE) {
  if (!inherits(lock, "cleanroom_lock"))
    stop("`lock` must be a cleanroom_lock object.", call. = FALSE)

  data       <- lock$data
  treatment  <- lock$treatment
  outcome    <- lock$outcome
  covariates <- lock$covariates

  A <- data[[treatment]]
  Y <- data[[outcome]]
  n <- nrow(data)

  # Fit baseline outcome model (covariates only; no treatment, no outcome peek)
  Q0_fml <- stats::reformulate(covariates, response = outcome)
  Q0_fit <- stats::glm(Q0_fml, data = data, family = stats::binomial())
  p_base <- as.numeric(stats::predict(Q0_fit, type = "response"))

  # Fit PS model once (GLM) for all replicates
  ps_fml <- stats::reformulate(covariates, response = treatment)
  ps_mod <- stats::glm(ps_fml, data = data, family = stats::binomial())
  ps_hat <- as.numeric(stats::predict(ps_mod, type = "response"))
  ps_hat <- pmax(pmin(ps_hat, 0.99), 0.01)

  all_results <- vector("list", length(effect_sizes))
  names(all_results) <- as.character(effect_sizes)

  for (es_idx in seq_along(effect_sizes)) {
    es <- effect_sizes[es_idx]
    if (verbose) message("Effect size: ", es)

    rep_results <- vector("list", reps)

    for (rep_i in seq_len(reps)) {
      set.seed(lock$seed + rep_i)

      # Synthetic outcomes: additive risk difference es for treated
      p1_sim <- pmin(p_base + es, 0.999)
      p0_sim <- p_base
      p_obs  <- ifelse(A == 1, p1_sim, p0_sim)
      Y_sim  <- stats::rbinom(n, 1L, p_obs)
      truth  <- mean(p1_sim) - mean(p0_sim)

      # IPTW (Hajek) estimate
      w_iptw <- ifelse(A == 1, 1 / ps_hat, 1 / (1 - ps_hat))
      r1_iptw <- stats::weighted.mean(Y_sim[A == 1], w_iptw[A == 1])
      r0_iptw <- stats::weighted.mean(Y_sim[A == 0], w_iptw[A == 0])
      est_iptw <- r1_iptw - r0_iptw
      se_sq_1 <- sum(w_iptw[A == 1]^2 * (Y_sim[A == 1] - r1_iptw)^2) /
        sum(w_iptw[A == 1])^2
      se_sq_0 <- sum(w_iptw[A == 0]^2 * (Y_sim[A == 0] - r0_iptw)^2) /
        sum(w_iptw[A == 0])^2
      se_iptw <- sqrt(se_sq_1 + se_sq_0)

      # TMLE (one-step) estimate
      ds <- data
      ds[[".Y_sim."]] <- Y_sim
      Q_fml  <- stats::reformulate(c(treatment, covariates),
                                   response = ".Y_sim.")
      Q_fit_s <- stats::glm(Q_fml, data = ds, family = stats::binomial())

      ds_a1 <- ds; ds_a1[[treatment]] <- 1L
      ds_a0 <- ds; ds_a0[[treatment]] <- 0L
      Q_a1  <- as.numeric(stats::predict(Q_fit_s, newdata = ds_a1,
                                          type = "response"))
      Q_a0  <- as.numeric(stats::predict(Q_fit_s, newdata = ds_a0,
                                          type = "response"))
      Q_aw  <- as.numeric(stats::predict(Q_fit_s, type = "response"))

      H_a1 <- A / ps_hat
      H_a0 <- (1 - A) / (1 - ps_hat)
      H_aw <- ifelse(A == 1, H_a1, -H_a0)

      epsilon <- tryCatch({
        Q_logit <- stats::qlogis(pmax(pmin(Q_aw, 0.999), 0.001))
        fluc    <- stats::glm(
          Y_sim ~ -1 + H_aw + offset(Q_logit),
          family = stats::binomial()
        )
        unname(stats::coef(fluc))
      }, error = function(e) 0)

      Q_a1_u <- plogis(stats::qlogis(pmax(pmin(Q_a1, 0.999), 0.001)) +
                         epsilon * (1 / ps_hat))
      Q_a0_u <- plogis(stats::qlogis(pmax(pmin(Q_a0, 0.999), 0.001)) +
                         epsilon * (-1 / (1 - ps_hat)))
      Q_aw_u <- plogis(stats::qlogis(pmax(pmin(Q_aw, 0.999), 0.001)) +
                         epsilon * H_aw)

      est_tmle <- mean(Q_a1_u) - mean(Q_a0_u)
      ic_tmle  <- H_aw * (Y_sim - Q_aw_u) +
        (Q_a1_u - Q_a0_u) - est_tmle
      se_tmle  <- sqrt(var(ic_tmle) / n)

      rep_results[[rep_i]] <- list(
        iptw  = list(est      = est_iptw,
                     se       = se_iptw,
                     ci_lower = est_iptw - 1.96 * se_iptw,
                     ci_upper = est_iptw + 1.96 * se_iptw),
        tmle  = list(est      = est_tmle,
                     se       = se_tmle,
                     ci_lower = est_tmle - 1.96 * se_tmle,
                     ci_upper = est_tmle + 1.96 * se_tmle),
        truth = truth
      )
    }

    all_results[[as.character(es)]] <- rep_results
  }

  # Aggregate performance metrics
  metrics_rows <- lapply(effect_sizes, function(es) {
    rr <- all_results[[as.character(es)]]
    truth_v <- vapply(rr, function(x) x$truth, numeric(1L))
    lapply(c("iptw", "tmle"), function(meth) {
      ests   <- vapply(rr, function(x) x[[meth]]$est,      numeric(1L))
      ci_los <- vapply(rr, function(x) x[[meth]]$ci_lower, numeric(1L))
      ci_his <- vapply(rr, function(x) x[[meth]]$ci_upper, numeric(1L))
      data.frame(
        effect_size = es,
        method      = meth,
        bias        = round(mean(ests - truth_v), 5L),
        rmse        = round(sqrt(mean((ests - truth_v)^2)), 5L),
        coverage    = round(mean(ci_los <= truth_v & truth_v <= ci_his), 3L),
        stringsAsFactors = FALSE
      )
    })
  })
  metrics <- do.call(rbind, unlist(metrics_rows, recursive = FALSE))

  result <- list(
    results      = all_results,
    metrics      = metrics,
    lock         = lock,
    effect_sizes = effect_sizes,
    reps         = reps,
    call         = match.call()
  )
  class(result) <- "plasmode_results"
  result
}


#' @export
print.plasmode_results <- function(x, ...) {
  cat("Plasmode-Simulation Feasibility Evaluation\n")
  cat("============================================\n")
  cat("Effect sizes evaluated:", paste(x$effect_sizes, collapse = ", "), "\n")
  cat("Replicates per effect size:", x$reps, "\n\n")
  cat("Performance Metrics:\n")
  print(x$metrics, row.names = FALSE)
  invisible(x)
}


#' Select Best TMLE Candidate Specification
#'
#' Applies the pre-specified selection rule to plasmode-simulation performance
#' metrics to choose the best analysis specification.
#'
#' @param sim_results A `plasmode_results` object from
#'   [run_plasmode_feasibility()].
#' @param rule Character; selection criterion. One of `"min_rmse"` (default),
#'   `"min_bias"`, or `"max_coverage"`.
#'
#' @return A character string naming the selected best method/specification.
#'
#' @export
select_tmle_candidate <- function(sim_results,
                                   rule = c("min_rmse", "min_bias",
                                            "max_coverage")) {
  if (!inherits(sim_results, "plasmode_results"))
    stop("`sim_results` must be a plasmode_results object.", call. = FALSE)
  rule <- match.arg(rule)

  m <- sim_results$metrics

  # Average metrics across effect sizes per method
  methods <- unique(m$method)
  summary_m <- do.call(rbind, lapply(methods, function(meth) {
    sub <- m[m$method == meth, ]
    data.frame(
      method   = meth,
      bias     = mean(abs(sub$bias)),
      rmse     = mean(sub$rmse),
      coverage = mean(sub$coverage),
      stringsAsFactors = FALSE
    )
  }))

  best <- switch(rule,
    min_rmse     = summary_m$method[which.min(summary_m$rmse)],
    min_bias     = summary_m$method[which.min(summary_m$bias)],
    max_coverage = summary_m$method[which.max(summary_m$coverage)]
  )

  message("Selected specification: '", best, "' (rule = '", rule, "')")
  best
}


# ── Stage 3: Conventional Workflows ──────────────────────────────────────

#' Run Propensity-Score Matching Workflow
#'
#' Performs greedy 1:1 nearest-neighbor matching on the logit of the
#' propensity score, then estimates the average treatment effect in the
#' matched sample as the simple risk difference.
#'
#' @param lock A `cleanroom_lock` from [create_analysis_lock()].
#' @param ps_fit A `ps_fit` object from [fit_ps_superlearner()].
#' @param caliper Numeric; maximum PS distance (logit scale) allowed for a
#'   match. Default: `0.2 * sd(logit(ps))`.
#'
#' @return An object of class `match_result` containing the causal risk
#'   difference estimate, SE, 95% CI, p-value, and the matched dataset.
#'
#' @export
run_match_workflow <- function(lock, ps_fit, caliper = NULL) {
  if (!inherits(lock, "cleanroom_lock"))
    stop("`lock` must be a cleanroom_lock object.", call. = FALSE)
  if (!inherits(ps_fit, "ps_fit"))
    stop("`ps_fit` must be a ps_fit object.", call. = FALSE)

  data      <- lock$data
  treatment <- lock$treatment
  outcome   <- lock$outcome
  A         <- data[[treatment]]
  Y         <- data[[outcome]]
  ps        <- ps_fit$ps

  logit_ps <- log(ps / (1 - ps))
  if (is.null(caliper))
    caliper <- 0.2 * sd(logit_ps)

  treated_idx  <- which(A == 1)
  control_idx  <- which(A == 0)
  n_treated    <- length(treated_idx)
  used_control <- logical(length(control_idx))
  matched_ctrl <- integer(n_treated)

  for (i in seq_len(n_treated)) {
    dists             <- abs(logit_ps[treated_idx[i]] - logit_ps[control_idx])
    dists[used_control] <- Inf
    best              <- which.min(dists)
    if (dists[best] <= caliper) {
      matched_ctrl[i]    <- control_idx[best]
      used_control[best] <- TRUE
    } else {
      matched_ctrl[i] <- NA_integer_
    }
  }

  valid        <- !is.na(matched_ctrl)
  match_idx    <- c(treated_idx[valid], matched_ctrl[valid])
  matched_data <- data[match_idx, , drop = FALSE]

  A_m  <- matched_data[[treatment]]
  Y_m  <- matched_data[[outcome]]
  n_m  <- sum(valid)

  r1   <- mean(Y_m[A_m == 1])
  r0   <- mean(Y_m[A_m == 0])
  rd   <- r1 - r0

  # Paired SE
  Y1   <- Y_m[A_m == 1]
  Y0   <- Y_m[A_m == 0][seq_len(n_m)]
  se   <- sqrt(var(Y1 - Y0) / n_m)
  ci_lo <- rd - 1.96 * se
  ci_hi <- rd + 1.96 * se
  p_val <- 2 * pnorm(-abs(rd / se))

  result <- list(
    estimate     = rd,
    se           = se,
    ci_lower     = ci_lo,
    ci_upper     = ci_hi,
    p_value      = p_val,
    r1           = r1,
    r0           = r0,
    n_matched    = n_m,
    n_unmatched  = sum(!valid),
    matched_data = matched_data,
    treatment    = treatment,
    outcome      = outcome,
    call         = match.call()
  )
  class(result) <- c("match_result", "cr_result")
  result
}


#' @export
print.match_result <- function(x, ...) {
  cat("Propensity Score Matching Workflow\n")
  cat("===================================\n")
  cat(sprintf("Matched pairs:   %d\n", x$n_matched))
  if (x$n_unmatched > 0L)
    cat(sprintf("Unmatched (caliper): %d\n", x$n_unmatched))
  cat(sprintf("Risk (treated):  %.4f\n", x$r1))
  cat(sprintf("Risk (control):  %.4f\n", x$r0))
  cat(sprintf("Risk Difference: %.4f  (95%% CI: %.4f, %.4f)\n",
              x$estimate, x$ci_lower, x$ci_upper))
  cat(sprintf("SE:              %.5f   p-value: %.4f\n", x$se, x$p_value))
  invisible(x)
}


#' Run IPTW Workflow
#'
#' Computes stabilized inverse probability of treatment weights from the
#' fitted propensity score and estimates the causal risk difference using
#' the Hajek (normalized) estimator with influence-curve-based variance.
#'
#' @param lock A `cleanroom_lock` from [create_analysis_lock()].
#' @param ps_fit A `ps_fit` object from [fit_ps_superlearner()].
#' @param trim Quantile for weight trimming (e.g., `0.01` trims at 1st and
#'   99th percentile). Default: `NULL` (no trimming).
#'
#' @return An object of class `iptw_result` containing the estimated risk
#'   difference, SE, 95% CI, p-value, and IPTW weights.
#'
#' @export
run_iptw_workflow <- function(lock, ps_fit, trim = NULL) {
  if (!inherits(lock, "cleanroom_lock"))
    stop("`lock` must be a cleanroom_lock object.", call. = FALSE)
  if (!inherits(ps_fit, "ps_fit"))
    stop("`ps_fit` must be a ps_fit object.", call. = FALSE)

  data      <- lock$data
  treatment <- lock$treatment
  outcome   <- lock$outcome
  A         <- data[[treatment]]
  Y         <- data[[outcome]]
  ps        <- ps_fit$ps
  n         <- nrow(data)

  # Stabilized IPTW weights
  p_trt <- mean(A)
  w     <- ifelse(A == 1, p_trt / ps, (1 - p_trt) / (1 - ps))

  if (!is.null(trim)) {
    lo <- quantile(w, trim)
    hi <- quantile(w, 1 - trim)
    w  <- pmin(pmax(w, lo), hi)
  }

  # Hajek estimator
  r1 <- weighted.mean(Y[A == 1], w[A == 1])
  r0 <- weighted.mean(Y[A == 0], w[A == 0])
  rd <- r1 - r0

  # Linearized variance for Hajek estimator
  se_sq_1 <- sum(w[A == 1]^2 * (Y[A == 1] - r1)^2) / sum(w[A == 1])^2
  se_sq_0 <- sum(w[A == 0]^2 * (Y[A == 0] - r0)^2) / sum(w[A == 0])^2
  se      <- sqrt(se_sq_1 + se_sq_0)

  ci_lo <- rd - 1.96 * se
  ci_hi <- rd + 1.96 * se
  p_val <- 2 * pnorm(-abs(rd / se))

  result <- list(
    estimate  = rd,
    se        = se,
    ci_lower  = ci_lo,
    ci_upper  = ci_hi,
    p_value   = p_val,
    r1        = r1,
    r0        = r0,
    weights   = w,
    ps        = ps,
    n         = n,
    treatment = treatment,
    outcome   = outcome,
    call      = match.call()
  )
  class(result) <- c("iptw_result", "cr_result")
  result
}


#' @export
print.iptw_result <- function(x, ...) {
  cat("IPTW Workflow\n")
  cat("=============\n")
  cat(sprintf("N:               %d\n", x$n))
  cat(sprintf("Risk (treated):  %.4f\n", x$r1))
  cat(sprintf("Risk (control):  %.4f\n", x$r0))
  cat(sprintf("Risk Difference: %.4f  (95%% CI: %.4f, %.4f)\n",
              x$estimate, x$ci_lower, x$ci_upper))
  cat(sprintf("SE:              %.5f   p-value: %.4f\n", x$se, x$p_value))
  invisible(x)
}


# ── Stage 3: Modular TMLE ────────────────────────────────────────────────

#' Fit TMLE Treatment Mechanism
#'
#' Fits the treatment mechanism g(A|W) = P(A=1|W) for use in the modular
#' TMLE workflow. Uses propensity scores from an existing `ps_fit` if
#' provided; otherwise re-estimates via SuperLearner. The outcome is never
#' accessed, making this step safe at any stage.
#'
#' @param lock A `cleanroom_lock` from [create_analysis_lock()].
#' @param ps_fit Optional `ps_fit` from [fit_ps_superlearner()]. If provided,
#'   the already-estimated propensity scores are reused.
#'
#' @return An object of class `tmle_mechanism` with `type = "treatment"`.
#'
#' @export
fit_tmle_treatment_mechanism <- function(lock, ps_fit = NULL) {
  if (!inherits(lock, "cleanroom_lock"))
    stop("`lock` must be a cleanroom_lock object.", call. = FALSE)

  if (!is.null(ps_fit)) {
    if (!inherits(ps_fit, "ps_fit"))
      stop("`ps_fit` must be a ps_fit object.", call. = FALSE)
    ps    <- ps_fit$ps
    g_mod <- ps_fit$sl_fit
  } else {
    temp  <- fit_ps_superlearner(lock)
    ps    <- temp$ps
    g_mod <- temp$sl_fit
  }

  result <- list(
    type       = "treatment",
    ps         = ps,
    g_fit      = g_mod,
    treatment  = lock$treatment,
    covariates = lock$covariates,
    data       = lock$data,
    lock       = lock,
    call       = match.call()
  )
  class(result) <- "tmle_mechanism"
  result
}


#' Fit TMLE Outcome Mechanism
#'
#' Fits the outcome mechanism Q(A,W) = E[Y|A,W] using SuperLearner (or
#' logistic regression as a fallback). This step accesses the outcome and
#' **must only be called in Stage 3** (after outcome unblinding).
#'
#' @param lock A `cleanroom_lock` from [create_analysis_lock()].
#' @param g_fit A `tmle_mechanism` object with `type = "treatment"` from
#'   [fit_tmle_treatment_mechanism()].
#' @param sl_library Optional SuperLearner library override. Defaults to
#'   `lock$sl_library`.
#'
#' @return An object of class `tmle_mechanism` with `type = "outcome"`
#'   containing initial outcome predictions `Q_a1`, `Q_a0`, and `Q_aw`.
#'
#' @export
fit_tmle_outcome_mechanism <- function(lock, g_fit, sl_library = NULL) {
  if (!inherits(lock, "cleanroom_lock"))
    stop("`lock` must be a cleanroom_lock object.", call. = FALSE)
  if (!inherits(g_fit, "tmle_mechanism") || g_fit$type != "treatment")
    stop("`g_fit` must be a tmle_mechanism of type 'treatment'.",
         call. = FALSE)

  if (is.null(sl_library)) sl_library <- lock$sl_library

  data       <- lock$data
  treatment  <- lock$treatment
  outcome    <- lock$outcome
  covariates <- lock$covariates
  AW         <- data[, c(treatment, covariates), drop = FALSE]
  Y          <- data[[outcome]]

  if (requireNamespace("SuperLearner", quietly = TRUE)) {
    set.seed(lock$seed + 1L)
    Q_sl <- SuperLearner::SuperLearner(
      Y          = Y,
      X          = AW,
      family     = binomial(),
      SL.library = sl_library
    )
    AW_a1   <- AW; AW_a1[[treatment]] <- 1L
    AW_a0   <- AW; AW_a0[[treatment]] <- 0L
    Q_a1    <- as.numeric(predict(Q_sl, newdata = AW_a1)$pred)
    Q_a0    <- as.numeric(predict(Q_sl, newdata = AW_a0)$pred)
    Q_aw    <- as.numeric(Q_sl$SL.predict)
    Q_fit_o <- Q_sl
  } else {
    Q_fml   <- stats::reformulate(c(treatment, covariates), response = outcome)
    Q_glm   <- stats::glm(Q_fml, data = data, family = stats::binomial())
    da1     <- data; da1[[treatment]] <- 1L
    da0     <- data; da0[[treatment]] <- 0L
    Q_a1    <- as.numeric(stats::predict(Q_glm, newdata = da1, type = "response"))
    Q_a0    <- as.numeric(stats::predict(Q_glm, newdata = da0, type = "response"))
    Q_aw    <- as.numeric(stats::predict(Q_glm, type = "response"))
    Q_fit_o <- Q_glm
  }

  result <- list(
    type       = "outcome",
    Q_a1       = Q_a1,
    Q_a0       = Q_a0,
    Q_aw       = Q_aw,
    Q_fit      = Q_fit_o,
    outcome    = outcome,
    treatment  = treatment,
    covariates = covariates,
    data       = data,
    lock       = lock,
    g_fit      = g_fit,
    call       = match.call()
  )
  class(result) <- "tmle_mechanism"
  result
}


#' Run TMLE Targeting Step
#'
#' Performs the TMLE fluctuation (targeting) step by fitting a logistic
#' submodel through the initial outcome estimates using the clever covariate
#' H(A,W) = A/g(W) - (1-A)/(1-g(W)). Returns updated outcome predictions,
#' the fluctuation parameter epsilon, and the efficient influence curve.
#'
#' @param g_fit A `tmle_mechanism` of type `"treatment"` from
#'   [fit_tmle_treatment_mechanism()].
#' @param Q_fit A `tmle_mechanism` of type `"outcome"` from
#'   [fit_tmle_outcome_mechanism()].
#'
#' @return An object of class `tmle_update` containing updated predictions,
#'   the fluctuation parameter, the efficient influence curve, and the
#'   initial point estimate.
#'
#' @export
run_tmle_targeting_step <- function(g_fit, Q_fit) {
  if (!inherits(g_fit, "tmle_mechanism") || g_fit$type != "treatment")
    stop("`g_fit` must be a tmle_mechanism of type 'treatment'.",
         call. = FALSE)
  if (!inherits(Q_fit, "tmle_mechanism") || Q_fit$type != "outcome")
    stop("`Q_fit` must be a tmle_mechanism of type 'outcome'.",
         call. = FALSE)

  data      <- Q_fit$data
  treatment <- Q_fit$treatment
  outcome   <- Q_fit$outcome
  A         <- data[[treatment]]
  Y         <- data[[outcome]]
  ps        <- g_fit$ps
  n         <- nrow(data)

  Q_a1 <- Q_fit$Q_a1
  Q_a0 <- Q_fit$Q_a0
  Q_aw <- Q_fit$Q_aw

  # Clever covariate
  H_a1 <-  1 / ps
  H_a0 <- -1 / (1 - ps)
  H_aw <- ifelse(A == 1, H_a1, H_a0)

  # Fluctuation via logistic submodel
  Q_aw_logit <- stats::qlogis(pmax(pmin(Q_aw, 0.999), 0.001))

  epsilon <- tryCatch({
    fluc <- stats::glm(
      Y ~ -1 + H_aw + offset(Q_aw_logit),
      family = stats::binomial()
    )
    unname(stats::coef(fluc))
  }, error = function(e) 0)

  # Updated predictions
  Q_a1_upd <- plogis(
    stats::qlogis(pmax(pmin(Q_a1, 0.999), 0.001)) + epsilon * H_a1
  )
  Q_a0_upd <- plogis(
    stats::qlogis(pmax(pmin(Q_a0, 0.999), 0.001)) + epsilon * H_a0
  )
  Q_aw_upd <- plogis(
    stats::qlogis(pmax(pmin(Q_aw, 0.999), 0.001)) + epsilon * H_aw
  )

  psi <- mean(Q_a1_upd) - mean(Q_a0_upd)

  # Efficient influence curve
  eic <- H_aw * (Y - Q_aw_upd) + (Q_a1_upd - Q_a0_upd) - psi

  result <- list(
    psi      = psi,
    epsilon  = epsilon,
    Q_a1_upd = Q_a1_upd,
    Q_a0_upd = Q_a0_upd,
    Q_aw_upd = Q_aw_upd,
    eic      = eic,
    g_fit    = g_fit,
    Q_fit    = Q_fit,
    treatment = treatment,
    outcome   = outcome,
    n         = n,
    call      = match.call()
  )
  class(result) <- "tmle_update"
  result
}


#' Extract Final TMLE Estimate
#'
#' Extracts the final TMLE point estimate, influence-curve standard error,
#' 95% confidence interval, p-value, and one-step diagnostics from a
#' `tmle_update` object.
#'
#' @param tmle_upd A `tmle_update` object from [run_tmle_targeting_step()].
#'
#' @return An object of class `tmle_fit` (inherits from `cr_result`).
#'
#' @export
extract_tmle_estimate <- function(tmle_upd) {
  if (!inherits(tmle_upd, "tmle_update"))
    stop("`tmle_upd` must be a tmle_update object from ",
         "run_tmle_targeting_step().", call. = FALSE)

  psi <- tmle_upd$psi
  eic <- tmle_upd$eic
  n   <- tmle_upd$n

  se      <- sqrt(var(eic) / n)
  ci_lo   <- psi - 1.96 * se
  ci_hi   <- psi + 1.96 * se
  p_value <- 2 * pnorm(-abs(psi / se))

  result <- list(
    estimates = list(
      ATE = list(
        estimate = psi,
        se       = se,
        ci_lower = ci_lo,
        ci_upper = ci_hi,
        p_value  = p_value
      )
    ),
    tmle_obj        = tmle_upd,
    influence_curve = eic,
    diagnostics = list(
      mean_eic = mean(eic),
      epsilon  = tmle_upd$epsilon,
      n        = n
    ),
    treatment = tmle_upd$treatment,
    outcome   = tmle_upd$outcome,
    type      = "modular_tmle",
    call      = match.call()
  )
  class(result) <- c("tmle_fit", "cr_result")
  result
}


# ── Summary ───────────────────────────────────────────────────────────────

#' Summarize Cleanroom Results
#'
#' Creates a side-by-side comparison table of estimates from multiple
#' cleanroom workflow objects (matching, IPTW, TMLE).
#'
#' @param fits A named or unnamed list of fitted workflow objects.
#'   Supported classes: `match_result`, `iptw_result`, `tmle_fit`.
#' @param ... Currently unused.
#'
#' @return A data.frame with one row per workflow containing columns
#'   `method`, `estimate`, `se`, `ci_lower`, `ci_upper`, and `p_value`.
#'
#' @export
summarize_cleanroom_results <- function(fits, ...) {
  if (!is.list(fits))
    stop("`fits` must be a list of fitted workflow objects.", call. = FALSE)

  rows <- lapply(seq_along(fits), function(i) {
    fit <- fits[[i]]
    nm  <- if (!is.null(names(fits)) && nzchar(names(fits)[i]))
      names(fits)[i] else ""

    if (inherits(fit, "match_result")) {
      method <- if (nzchar(nm)) nm else "PS Matching"
      data.frame(method   = method,
                 estimate = fit$estimate,
                 se       = fit$se,
                 ci_lower = fit$ci_lower,
                 ci_upper = fit$ci_upper,
                 p_value  = fit$p_value,
                 stringsAsFactors = FALSE)

    } else if (inherits(fit, "iptw_result")) {
      method <- if (nzchar(nm)) nm else "IPTW"
      data.frame(method   = method,
                 estimate = fit$estimate,
                 se       = fit$se,
                 ci_lower = fit$ci_lower,
                 ci_upper = fit$ci_upper,
                 p_value  = fit$p_value,
                 stringsAsFactors = FALSE)

    } else if (inherits(fit, "tmle_fit")) {
      method <- if (nzchar(nm)) nm else "TMLE"
      ate    <- fit$estimates$ATE
      data.frame(method   = method,
                 estimate = ate$estimate,
                 se       = ate$se,
                 ci_lower = ate$ci_lower,
                 ci_upper = ate$ci_upper,
                 p_value  = ate$p_value,
                 stringsAsFactors = FALSE)

    } else {
      warning("Unsupported fit class: ",
              paste(class(fit), collapse = ", "), ". Skipping.",
              call. = FALSE)
      NULL
    }
  })

  rows <- rows[!vapply(rows, is.null, logical(1L))]
  if (length(rows) == 0L)
    stop("No supported workflow objects found in `fits`.", call. = FALSE)

  result <- do.call(rbind, rows)
  rownames(result) <- NULL
  result
}
