#' Data Quality Stress Testing via Plasmode Simulation
#'
#' Extends the plasmode feasibility framework to evaluate TMLE candidate
#' robustness under controlled data quality degradations: covariate
#' missingness, treatment misclassification, outcome misclassification,
#' and unmeasured confounding.
#'
#' @name plasmode_dq
NULL


# ── Internal: Single-Rep TMLE Fit ────────────────────────────────────────

#' Fit one plasmode replicate for a single TMLE candidate
#'
#' @param Y_sim Simulated outcome vector.
#' @param A Treatment vector.
#' @param W_data Covariate data.frame (may be degraded).
#' @param treatment Character; treatment column name.
#' @param covariates Character vector of covariate names.
#' @param cand A \code{tmle_candidate_spec}.
#' @param n Sample size.
#'
#' @return List with est, se, ci_lower, ci_upper.
#' @keywords internal
.plasmode_fit_one_candidate <- function(Y_sim, A, W_data, treatment,
                                         covariates, cand, n) {
  g_lib <- cand$g_library
  q_lib <- cand$q_library

  # ── g-model ──
  use_sl <- requireNamespace("SuperLearner", quietly = TRUE) &&
    !identical(g_lib, "SL.glm")

  if (use_sl) {
    g_sl   <- SuperLearner::SuperLearner(
      Y = A, X = W_data[, covariates, drop = FALSE],
      family = binomial(), SL.library = g_lib,
      env = asNamespace("SuperLearner")
    )
    ps_hat <- as.numeric(g_sl$SL.predict)
  } else {
    ps_fml <- stats::reformulate(covariates, response = treatment)
    ds_g   <- W_data; ds_g[[treatment]] <- A
    ps_mod <- stats::glm(ps_fml, data = ds_g, family = stats::binomial())
    ps_hat <- as.numeric(stats::predict(ps_mod, type = "response"))
  }
  ps_hat <- pmax(pmin(ps_hat, 1 - cand$truncation), cand$truncation)

  # ── Q-model ──
  ds <- W_data
  ds[[treatment]] <- A
  ds[[".Y_sim."]] <- Y_sim
  AW <- ds[, c(treatment, covariates), drop = FALSE]

  use_sl_q <- requireNamespace("SuperLearner", quietly = TRUE) &&
    !identical(q_lib, "SL.glm")

  if (use_sl_q) {
    q_sl <- SuperLearner::SuperLearner(
      Y = Y_sim, X = AW, family = binomial(), SL.library = q_lib,
      env = asNamespace("SuperLearner")
    )
    AW_a1 <- AW; AW_a1[[treatment]] <- 1L
    AW_a0 <- AW; AW_a0[[treatment]] <- 0L
    Q_a1  <- as.numeric(predict(q_sl, newdata = AW_a1)$pred)
    Q_a0  <- as.numeric(predict(q_sl, newdata = AW_a0)$pred)
    Q_aw  <- as.numeric(q_sl$SL.predict)
  } else {
    Q_fml   <- stats::reformulate(c(treatment, covariates),
                                   response = ".Y_sim.")
    Q_fit_s <- stats::glm(Q_fml, data = ds, family = stats::binomial())
    ds_a1 <- ds; ds_a1[[treatment]] <- 1L
    ds_a0 <- ds; ds_a0[[treatment]] <- 0L
    Q_a1 <- as.numeric(stats::predict(Q_fit_s, newdata = ds_a1,
                                       type = "response"))
    Q_a0 <- as.numeric(stats::predict(Q_fit_s, newdata = ds_a0,
                                       type = "response"))
    Q_aw <- as.numeric(stats::predict(Q_fit_s, type = "response"))
  }

  # ── Targeting ──
  H_a1 <-  1 / ps_hat
  H_a0 <- -1 / (1 - ps_hat)
  H_aw <- ifelse(A == 1, H_a1, H_a0)

  Q_aw_logit <- stats::qlogis(pmax(pmin(Q_aw, 0.999), 0.001))
  epsilon <- tryCatch({
    fluc <- stats::glm(Y_sim ~ -1 + H_aw + offset(Q_aw_logit),
                        family = stats::binomial())
    unname(stats::coef(fluc))
  }, error = function(e) 0)

  Q_a1_u <- stats::plogis(
    stats::qlogis(pmax(pmin(Q_a1, 0.999), 0.001)) + epsilon * H_a1)
  Q_a0_u <- stats::plogis(
    stats::qlogis(pmax(pmin(Q_a0, 0.999), 0.001)) + epsilon * H_a0)
  Q_aw_u <- stats::plogis(Q_aw_logit + epsilon * H_aw)

  est <- mean(Q_a1_u) - mean(Q_a0_u)
  eic <- H_aw * (Y_sim - Q_aw_u) + (Q_a1_u - Q_a0_u) - est
  se  <- sqrt(var(eic) / n)

  list(est = est, se = se,
       ci_lower = est - 1.96 * se,
       ci_upper = est + 1.96 * se)
}


# ── Internal: DQ Degradation Functions ───────────────────────────────────

#' Introduce MCAR missingness into covariates
#' @keywords internal
.degrade_missingness <- function(data, covariates, fraction, seed_offset) {
  for (v in covariates) {
    n <- nrow(data)
    n_miss <- round(n * fraction)
    if (n_miss > 0) {
      idx <- sample.int(n, n_miss)
      data[[v]][idx] <- NA
    }
  }
  # Impute with column medians (matching sanitize_covariates behaviour)
  for (v in covariates) {
    bad <- is.na(data[[v]])
    if (any(bad)) {
      med <- median(data[[v]][!bad], na.rm = TRUE)
      data[[v]][bad] <- if (is.na(med)) 0 else med
    }
  }
  data
}


#' Flip treatment assignment (non-differential misclassification)
#' @keywords internal
.degrade_treatment <- function(A, rate) {
  n <- length(A)
  n_flip <- round(n * rate)
  if (n_flip > 0) {
    idx <- sample.int(n, n_flip)
    A[idx] <- 1L - A[idx]
  }
  A
}


#' Apply outcome misclassification
#' @keywords internal
.degrade_outcome <- function(Y, sensitivity, specificity) {
  Y_out <- Y
  # True positives misclassified as negatives: P(Y*=0 | Y=1) = 1 - sens
  pos_idx <- which(Y == 1)
  if (length(pos_idx) > 0) {
    flip_pos <- stats::rbinom(length(pos_idx), 1, 1 - sensitivity)
    Y_out[pos_idx[flip_pos == 1]] <- 0L
  }
  # True negatives misclassified as positives: P(Y*=1 | Y=0) = 1 - spec
  neg_idx <- which(Y == 0)
  if (length(neg_idx) > 0) {
    flip_neg <- stats::rbinom(length(neg_idx), 1, 1 - specificity)
    Y_out[neg_idx[flip_neg == 1]] <- 1L
  }
  Y_out
}


# ── Main: DQ Stress Test ────────────────────────────────────────────────

#' Run Data-Quality Stress Tests via Plasmode Simulation
#'
#' Evaluates TMLE candidate robustness under controlled data quality
#' degradations.  For each DQ scenario and severity level, plasmode
#' outcomes are generated from the real covariate distribution and then
#' degraded before the TMLE candidate is fit.  This is an outcome-blind
#' procedure: only synthetic outcomes are used.
#'
#' @section Clean-room stage: Stage 2b (pre-outcome).
#'
#' @param lock A \code{cleanroom_lock}.
#' @param tmle_candidates A list of \code{tmle_candidate_spec} objects.
#' @param effect_sizes Numeric vector of true risk differences.
#'   Default: \code{c(0.05)}.
#' @param reps Integer; plasmode replicates per scenario.
#'   Default: 50 (fewer than standard plasmode for speed).
#' @param data_quality_scenarios A named list of DQ scenarios.  See Details.
#' @param verbose Logical.
#'
#' @details
#' The \code{data_quality_scenarios} list may contain any combination of:
#'
#' \describe{
#'   \item{\code{covariate_missingness}}{List with \code{fractions}
#'     (numeric vector, e.g., \code{c(0.05, 0.10, 0.20)}) and optional
#'     \code{variables} (character vector; default: all covariates).}
#'   \item{\code{treatment_misclass}}{List with \code{rates}
#'     (numeric vector, e.g., \code{c(0.02, 0.05, 0.10)}).}
#'   \item{\code{outcome_misclass}}{List with \code{sensitivity} and
#'     \code{specificity} (numeric vectors of equal length).}
#'   \item{\code{unmeasured_confounding}}{List with \code{U_prevalence},
#'     \code{U_treatment_OR}, and \code{U_outcome_OR} (numeric vectors;
#'     \code{U_treatment_OR} and \code{U_outcome_OR} must be equal length).}
#' }
#'
#' @return An object of class \code{plasmode_dq_results} containing:
#'   \describe{
#'     \item{metrics}{data.frame with columns: scenario, level,
#'       effect_size, candidate, bias, rmse, coverage, emp_sd,
#'       mean_se, se_cal.}
#'     \item{scenarios}{The input \code{data_quality_scenarios}.}
#'     \item{baseline}{Metrics under \code{scenario = "none"} (no degradation).}
#'   }
#'
#' @export
run_plasmode_dq_stress <- function(lock,
                                    tmle_candidates = NULL,
                                    effect_sizes    = c(0.05),
                                    reps            = 50L,
                                    data_quality_scenarios = list(),
                                    verbose         = TRUE) {
  if (!inherits(lock, "cleanroom_lock"))
    stop("`lock` must be a cleanroom_lock object.", call. = FALSE)

  if (is.null(tmle_candidates))
    tmle_candidates <- expand_tmle_candidate_grid()
  validate_tmle_candidates(tmle_candidates)

  data       <- lock$data
  treatment  <- lock$treatment
  outcome    <- lock$outcome
  covariates <- lock$covariates
  A          <- data[[treatment]]
  n          <- nrow(data)

  # Baseline outcome model for plasmode DGP
  Q0_fml <- stats::reformulate(covariates, response = outcome)
  Q0_fit <- stats::glm(Q0_fml, data = data, family = stats::binomial())
  p_base <- as.numeric(stats::predict(Q0_fit, type = "response"))

  cand_ids <- vapply(tmle_candidates, function(x) x$candidate_id, character(1))

  # ── Build scenario grid ──────────────────────────────────────────────
  scenario_grid <- data.frame(
    scenario = "none", level = "0", stringsAsFactors = FALSE
  )

  dqs <- data_quality_scenarios

  if (!is.null(dqs$covariate_missingness)) {
    for (f in dqs$covariate_missingness$fractions) {
      scenario_grid <- rbind(scenario_grid, data.frame(
        scenario = "cov_miss", level = as.character(f),
        stringsAsFactors = FALSE))
    }
  }

  if (!is.null(dqs$treatment_misclass)) {
    for (r in dqs$treatment_misclass$rates) {
      scenario_grid <- rbind(scenario_grid, data.frame(
        scenario = "trt_misclass", level = as.character(r),
        stringsAsFactors = FALSE))
    }
  }

  if (!is.null(dqs$outcome_misclass)) {
    sens <- dqs$outcome_misclass$sensitivity
    spec <- dqs$outcome_misclass$specificity
    for (i in seq_along(sens)) {
      scenario_grid <- rbind(scenario_grid, data.frame(
        scenario = "out_misclass",
        level = sprintf("sens%.2f_spec%.2f", sens[i], spec[i]),
        stringsAsFactors = FALSE))
    }
  }

  if (!is.null(dqs$unmeasured_confounding)) {
    u_trt <- dqs$unmeasured_confounding$U_treatment_OR
    u_out <- dqs$unmeasured_confounding$U_outcome_OR
    u_prev <- dqs$unmeasured_confounding$U_prevalence
    for (i in seq_along(u_trt)) {
      scenario_grid <- rbind(scenario_grid, data.frame(
        scenario = "unmeasured_U",
        level = sprintf("OR_trt%.1f_out%.1f", u_trt[i], u_out[i]),
        stringsAsFactors = FALSE))
    }
  }

  if (verbose) {
    cat(sprintf("DQ Stress Test: %d scenarios x %d effect sizes x %d reps x %d candidates\n",
                nrow(scenario_grid), length(effect_sizes), reps, length(cand_ids)))
  }

  # ── Run plasmode across scenario grid ─────────────────────────────────
  all_metrics <- list()

  for (sg_i in seq_len(nrow(scenario_grid))) {
    sc_name <- scenario_grid$scenario[sg_i]
    sc_level <- scenario_grid$level[sg_i]

    if (verbose) cat(sprintf("  Scenario: %s [%s]\n", sc_name, sc_level))

    for (es in effect_sizes) {
      rep_results <- vector("list", reps)

      for (rep_i in seq_len(reps)) {
        set.seed(lock$seed + rep_i + sg_i * 10000L)

        # Generate synthetic outcome
        A_rep <- A
        p1_sim <- pmin(p_base + es, 0.999)
        p0_sim <- p_base

        # ── Unmeasured confounding: modify DGP ──
        if (sc_name == "unmeasured_U") {
          u_prev <- dqs$unmeasured_confounding$U_prevalence
          idx <- which(scenario_grid$scenario == "unmeasured_U" &
                         scenario_grid$level == sc_level)
          u_idx <- sum(scenario_grid$scenario[1:idx] == "unmeasured_U")
          u_trt_or <- dqs$unmeasured_confounding$U_treatment_OR[u_idx]
          u_out_or <- dqs$unmeasured_confounding$U_outcome_OR[u_idx]

          U <- stats::rbinom(n, 1, u_prev)
          # Modify treatment (re-sample A given U)
          ps_orig <- p_base  # approximate PS from Q0 model (imperfect but fast)
          lp_ps <- stats::qlogis(pmax(pmin(ps_orig, 0.999), 0.001))
          lp_ps_u <- lp_ps + log(u_trt_or) * U
          ps_u <- stats::plogis(lp_ps_u)
          A_rep <- stats::rbinom(n, 1, ps_u)

          # Modify outcome probabilities given U
          lp1 <- stats::qlogis(pmax(pmin(p1_sim, 0.999), 0.001)) +
            log(u_out_or) * U
          lp0 <- stats::qlogis(pmax(pmin(p0_sim, 0.999), 0.001)) +
            log(u_out_or) * U
          p1_sim <- stats::plogis(lp1)
          p0_sim <- stats::plogis(lp0)
        }

        p_obs  <- ifelse(A_rep == 1, p1_sim, p0_sim)
        Y_sim  <- stats::rbinom(n, 1L, p_obs)
        truth  <- mean(p1_sim) - mean(p0_sim)

        # ── Treatment misclassification ──
        A_fit <- A_rep
        if (sc_name == "trt_misclass") {
          rate <- as.numeric(sc_level)
          A_fit <- .degrade_treatment(A_rep, rate)
        }

        # ── Outcome misclassification ──
        Y_fit <- Y_sim
        if (sc_name == "out_misclass") {
          parts <- strsplit(sc_level, "_")[[1]]
          sens_val <- as.numeric(sub("sens", "", parts[1]))
          spec_val <- as.numeric(sub("spec", "", parts[2]))
          Y_fit <- .degrade_outcome(Y_sim, sens_val, spec_val)
        }

        # ── Covariate missingness ──
        W_fit <- data
        if (sc_name == "cov_miss") {
          frac <- as.numeric(sc_level)
          miss_vars <- dqs$covariate_missingness$variables
          if (is.null(miss_vars)) miss_vars <- covariates
          W_fit <- .degrade_missingness(data, miss_vars, frac,
                                         seed_offset = rep_i)
        }

        # ── Fit each candidate ──
        cand_results <- list()
        for (cand in tmle_candidates) {
          cand_results[[cand$candidate_id]] <- tryCatch(
            .plasmode_fit_one_candidate(
              Y_sim = Y_fit, A = A_fit, W_data = W_fit,
              treatment = treatment, covariates = covariates,
              cand = cand, n = n
            ),
            error = function(e) {
              list(est = NA_real_, se = NA_real_,
                   ci_lower = NA_real_, ci_upper = NA_real_)
            }
          )
        }

        rep_results[[rep_i]] <- c(cand_results, list(.truth = truth))
      }

      # ── Aggregate metrics ──
      truth_v <- vapply(rep_results, function(x) x$.truth, numeric(1))

      for (cid in cand_ids) {
        ests   <- vapply(rep_results, function(x) x[[cid]]$est, numeric(1))
        ses    <- vapply(rep_results, function(x) x[[cid]]$se, numeric(1))
        ci_los <- vapply(rep_results, function(x) x[[cid]]$ci_lower, numeric(1))
        ci_his <- vapply(rep_results, function(x) x[[cid]]$ci_upper, numeric(1))

        valid <- !is.na(ests)
        if (sum(valid) == 0) next

        ests_v   <- ests[valid]
        ses_v    <- ses[valid]
        truth_vv <- truth_v[valid]

        emp_sd  <- sd(ests_v)
        mean_se <- mean(ses_v)

        row <- data.frame(
          scenario    = sc_name,
          level       = sc_level,
          effect_size = es,
          candidate   = cid,
          bias        = round(mean(ests_v - truth_vv), 5),
          rmse        = round(sqrt(mean((ests_v - truth_vv)^2)), 5),
          coverage    = round(mean(ci_los[valid] <= truth_vv &
                                     truth_vv <= ci_his[valid]), 3),
          emp_sd      = round(emp_sd, 5),
          mean_se     = round(mean_se, 5),
          se_cal      = round(if (emp_sd > 0) mean_se / emp_sd else NA, 3),
          n_converged = sum(valid),
          stringsAsFactors = FALSE
        )
        all_metrics <- c(all_metrics, list(row))
      }
    }
  }

  metrics <- do.call(rbind, all_metrics)

  result <- list(
    metrics   = metrics,
    scenarios = data_quality_scenarios,
    baseline  = metrics[metrics$scenario == "none", ],
    lock      = lock,
    reps      = reps,
    call      = match.call()
  )
  class(result) <- "plasmode_dq_results"
  result
}


#' @export
print.plasmode_dq_results <- function(x, ...) {
  cat("Plasmode Data-Quality Stress Test\n")
  cat("==================================\n")
  cat(sprintf("Replicates per scenario: %d\n", x$reps))
  cat(sprintf("Total scenario-candidate rows: %d\n\n", nrow(x$metrics)))

  scenarios <- unique(x$metrics$scenario)
  for (sc in scenarios) {
    sub <- x$metrics[x$metrics$scenario == sc, ]
    cat(sprintf("--- %s ---\n", sc))
    print(sub[, c("level", "candidate", "bias", "rmse", "coverage", "se_cal")],
          row.names = FALSE)
    cat("\n")
  }
  invisible(x)
}


#' Summarise DQ Stress Test as a Degradation Table
#'
#' Computes the relative change in bias, RMSE, and coverage compared
#' to the baseline (no degradation) for each DQ scenario and level.
#'
#' @param dq_results A \code{plasmode_dq_results} object.
#'
#' @return A data.frame with relative degradation metrics.
#'
#' @export
summarize_dq_degradation <- function(dq_results) {
  if (!inherits(dq_results, "plasmode_dq_results"))
    stop("`dq_results` must be a plasmode_dq_results object.", call. = FALSE)

  m <- dq_results$metrics
  base <- m[m$scenario == "none", ]
  degraded <- m[m$scenario != "none", ]

  if (nrow(degraded) == 0) return(data.frame())

  rows <- lapply(seq_len(nrow(degraded)), function(i) {
    d <- degraded[i, ]
    b <- base[base$candidate == d$candidate &
                base$effect_size == d$effect_size, ]
    if (nrow(b) == 0) return(NULL)

    data.frame(
      scenario      = d$scenario,
      level         = d$level,
      candidate     = d$candidate,
      bias_baseline = b$bias,
      bias_degraded = d$bias,
      bias_change   = round(abs(d$bias) - abs(b$bias), 5),
      rmse_baseline = b$rmse,
      rmse_degraded = d$rmse,
      rmse_ratio    = round(d$rmse / max(b$rmse, 1e-6), 3),
      cov_baseline  = b$coverage,
      cov_degraded  = d$coverage,
      cov_drop      = round(b$coverage - d$coverage, 3),
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, Filter(Negate(is.null), rows))
}
