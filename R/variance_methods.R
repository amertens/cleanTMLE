# ============================================================================
# Native nonparametric bootstrap variance and the FIORD second-stage
# variance-method selector.
#
# The FIORD outcome-blind selection procedure (Nance et al. 2026) has two
# stages: (1) screen point estimators by oracle coverage and lock the
# lowest-variance survivor; (2) select the variance method that achieves
# nominal coverage on the locked point estimator. cleanTMLE's
# `select_tmle_candidate(rule = "fiord_two_stage")` implements stage 1.
# These functions implement stage 2: a native nonparametric bootstrap
# variance for the point-treatment risk difference, and a selector that
# chooses the variance method whose oracle coverage is closest to nominal.
# ============================================================================

# ── Internal helpers ─────────────────────────────────────────────────────────

#' Greedy 1:1 nearest-neighbour matching on the logit propensity score
#'
#' Pure-R implementation that mirrors `run_match_workflow()` exactly:
#' treated units are matched in order; each control unit may be used at
#' most once; matches that exceed `caliper` are dropped.
#'
#' @param logit_ps Numeric vector of logit-propensity scores (all units).
#' @param A Integer/numeric treatment indicator (1 = treated, 0 = control).
#' @param caliper Maximum allowable distance on the logit-PS scale. Default
#'   is `0.2 * sd(logit_ps)` (Rosenbaum and Rubin 1985 rule of thumb).
#'
#' @return A list with elements `match_idx` (integer row indices of the
#'   matched dataset, treated rows followed by their control partners),
#'   `treated_valid`, `control_valid`, and `n_matched`.
#' @keywords internal
.match_nn_pairs <- function(logit_ps, A, caliper = NULL) {
  if (is.null(caliper)) caliper <- 0.2 * stats::sd(logit_ps)
  treated_idx  <- which(A == 1)
  control_idx  <- which(A == 0)
  used_ctrl    <- logical(length(control_idx))
  matched_ctrl <- integer(length(treated_idx))
  for (i in seq_along(treated_idx)) {
    dists <- abs(logit_ps[treated_idx[i]] - logit_ps[control_idx])
    dists[used_ctrl] <- Inf
    best <- which.min(dists)
    if (dists[best] <= caliper) {
      matched_ctrl[i] <- control_idx[best]
      used_ctrl[best] <- TRUE
    } else {
      matched_ctrl[i] <- NA_integer_
    }
  }
  valid <- !is.na(matched_ctrl)
  list(match_idx     = c(treated_idx[valid], matched_ctrl[valid]),
       treated_valid = treated_idx[valid],
       control_valid = matched_ctrl[valid],
       n_matched     = sum(valid))
}

#' Point estimate of the marginal risk difference
#'
#' Internal helper used by the bootstrap. Returns a single risk-difference
#' point estimate from a data frame using TMLE, stabilised IPTW, or
#' propensity-score-matched TMLE.
#'
#' @keywords internal
.rd_point_estimate <- function(data, treatment, outcome, covariates,
                               estimator = c("tmle", "iptw", "match_tmle"),
                               sl_library = c("SL.glm"), truncate = 0.01) {
  estimator <- match.arg(estimator)
  A <- data[[treatment]]
  Y <- data[[outcome]]
  W <- data[, covariates, drop = FALSE]

  if (estimator == "tmle") {
    if (!requireNamespace("tmle", quietly = TRUE))
      stop("Package 'tmle' is required for estimator = 'tmle'.", call. = FALSE)
    fit <- tmle::tmle(Y = Y, A = A, W = W, family = "binomial",
                      Q.SL.library = sl_library, g.SL.library = sl_library)
    return(fit$estimates$ATE$psi)
  }

  if (estimator == "match_tmle") {
    if (!requireNamespace("tmle", quietly = TRUE))
      stop("Package 'tmle' is required for estimator = 'match_tmle'.",
           call. = FALSE)
    # Step 1: Fit PS on full data via GLM (truncated).
    ps_fml <- stats::reformulate(covariates, response = treatment)
    ps <- stats::predict(
      stats::glm(ps_fml, data = data, family = stats::binomial()),
      type = "response")
    ps <- pmin(pmax(ps, truncate), 1 - truncate)
    # Step 2: Greedy 1:1 NN matching on logit-PS.
    pairs <- .match_nn_pairs(log(ps / (1 - ps)), A)
    if (pairs$n_matched < 10L)
      stop("Too few matched pairs (<10) for match_tmle bootstrap.",
           call. = FALSE)
    md <- data[pairs$match_idx, , drop = FALSE]
    # Step 3: TMLE on the matched cohort.
    fit <- tmle::tmle(Y = md[[outcome]], A = md[[treatment]],
                      W = md[, covariates, drop = FALSE],
                      family = "binomial",
                      Q.SL.library = sl_library, g.SL.library = sl_library)
    return(fit$estimates$ATE$psi)
  }

  # Stabilised IPTW
  ps_fml <- stats::reformulate(covariates, response = treatment)
  ps <- stats::predict(stats::glm(ps_fml, data = data, family = stats::binomial()),
                       type = "response")
  ps <- pmin(pmax(ps, truncate), 1 - truncate)
  pA <- mean(A)
  w  <- ifelse(A == 1, pA / ps, (1 - pA) / (1 - ps))
  stats::weighted.mean(Y[A == 1], w[A == 1]) -
    stats::weighted.mean(Y[A == 0], w[A == 0])
}

#' Nonparametric bootstrap variance for a point-treatment risk difference
#'
#' Computes a nonparametric bootstrap standard error and percentile
#' confidence interval for the marginal risk difference. This is the
#' principled variance method when the influence-function variance is known
#' to be conservative or invalid: stabilised IPTW (the IF variance treats the
#' propensity score as known; Robins and Rotnitzky 1992) and any estimator
#' run on a non-i.i.d. sample such as a matched cohort. It is the recommended
#' second-stage variance method in the FIORD selector (Nance et al. 2026).
#'
#' @section Why bootstrap for matched estimators:
#' The paired-difference SE returned by [run_match_workflow()] conditions on
#' the fixed matched dataset and does not propagate matching-draw variance --
#' the variability arising from which controls happen to be selected on
#' re-sampling. Abadie and Imbens (2008, 2016) show that the standard paired
#' SE is inconsistent for this reason. With `estimator = "match_tmle"`, each
#' bootstrap resample re-runs the full pipeline (PS fit -> greedy 1:1 NN
#' matching -> TMLE on the matched cohort), capturing matching-draw variance
#' and yielding a consistent SE. Monte Carlo evidence in Mertens (2026) shows
#' this corrects modest under-coverage (from 0.89 to 0.95) in good-overlap
#' settings where the IF-SE is anti-conservative for matched estimators.
#'
#' @param data A data frame.
#' @param treatment,outcome,covariates Column names.
#' @param estimator One of `"tmle"` (delegates to [tmle::tmle()]),
#'   `"iptw"` (stabilised inverse-probability weighting), or
#'   `"match_tmle"` (greedy 1:1 NN propensity-score matching followed by
#'   TMLE on the matched cohort; the bootstrap re-runs the full pipeline
#'   on each resample, capturing matching-draw variance).
#' @param R Number of bootstrap resamples.
#' @param sl_library SuperLearner library for the TMLE nuisance fits.
#' @param truncate Propensity-score truncation for IPTW and match_tmle.
#' @param conf_level Confidence level for the percentile interval.
#' @param seed Random seed.
#'
#' @return A list with `estimate`, `se` (bootstrap SD), `ci` (percentile
#'   interval), `method`, and `R_effective` (resamples that did not error).
#'
#' @references
#' Abadie, A. and Imbens, G. W. (2008). On the failure of the bootstrap for
#' matching estimators. *Econometrica*, 76(6), 1537--1557.
#'
#' Abadie, A. and Imbens, G. W. (2016). Matching on the estimated propensity
#' score. *Econometrica*, 84(2), 781--807.
#'
#' @export
bootstrap_rd_variance <- function(data, treatment, outcome, covariates,
                                  estimator = c("tmle", "iptw", "match_tmle"),
                                  R = 1000L, sl_library = c("SL.glm"),
                                  truncate = 0.01, conf_level = 0.95,
                                  seed = 42L) {
  estimator <- match.arg(estimator)
  set.seed(seed)
  n <- nrow(data)
  point <- .rd_point_estimate(data, treatment, outcome, covariates,
                              estimator, sl_library, truncate)
  boots <- vapply(seq_len(R), function(b) {
    idx <- sample.int(n, n, replace = TRUE)
    tryCatch(
      .rd_point_estimate(data[idx, , drop = FALSE], treatment, outcome,
                         covariates, estimator, sl_library, truncate),
      error = function(e) NA_real_)
  }, numeric(1))
  boots <- boots[is.finite(boots)]
  if (length(boots) < 2L)
    stop("Bootstrap produced too few valid resamples.", call. = FALSE)
  a <- (1 - conf_level) / 2
  list(estimate    = point,
       se          = stats::sd(boots),
       ci          = unname(stats::quantile(boots, c(a, 1 - a))),
       method      = "nonparametric_bootstrap",
       R_effective = length(boots))
}

#' Select the variance method that achieves nominal oracle coverage (FIORD stage 2)
#'
#' Given the locked point estimator and a set of synthetic (plasmode)
#' datasets with a known synthetic truth, evaluates each candidate variance
#' method by its empirical (oracle) coverage and returns the method whose
#' coverage is closest to the nominal target while not falling below it. This
#' is the second stage of the FIORD two-stage selector (Nance et al. 2026):
#' the point estimator is fixed, and only the variance method varies.
#'
#' @param plasmode_datasets A list of data frames (synthetic outcomes).
#' @param truth The synthetic-truth risk difference used to generate them.
#' @param treatment,outcome,covariates Column names.
#' @param estimator Locked point estimator: `"tmle"`, `"iptw"`, or
#'   `"match_tmle"`. For `"match_tmle"`, the `"influence"` method uses the
#'   paired-difference SE on the matched cohort (analogous to
#'   [run_match_workflow()]); the `"bootstrap"` method re-runs the full
#'   matching pipeline on each resample via [bootstrap_rd_variance()].
#' @param methods Variance methods to compare. `"influence"` uses the
#'   IF/paired-difference analytic SE; `"bootstrap"` uses
#'   [bootstrap_rd_variance()].
#' @param target_coverage Nominal coverage (default 0.95).
#' @param R_boot Bootstrap resamples per dataset for the bootstrap method.
#' @param sl_library SuperLearner library for TMLE.
#' @param truncate IPTW / match_tmle PS truncation.
#'
#' @return A list with `selected` (the chosen method), `coverage` (named
#'   vector of empirical coverage by method), and `n_datasets`.
#' @export
select_variance_method <- function(plasmode_datasets, truth,
                                    treatment, outcome, covariates,
                                    estimator = c("tmle", "iptw", "match_tmle"),
                                    methods = c("influence", "bootstrap"),
                                    target_coverage = 0.95,
                                    R_boot = 200L, sl_library = c("SL.glm"),
                                    truncate = 0.01) {
  estimator <- match.arg(estimator)
  methods   <- match.arg(methods, several.ok = TRUE)
  z <- stats::qnorm(1 - (1 - target_coverage) / 2)

  covers <- matrix(NA_real_, nrow = length(plasmode_datasets),
                   ncol = length(methods),
                   dimnames = list(NULL, methods))

  for (i in seq_along(plasmode_datasets)) {
    d <- plasmode_datasets[[i]]
    if ("influence" %in% methods) {
      ci <- tryCatch({
        if (estimator == "tmle") {
          fit <- tmle::tmle(Y = d[[outcome]], A = d[[treatment]],
                            W = d[, covariates, drop = FALSE], family = "binomial",
                            Q.SL.library = sl_library, g.SL.library = sl_library)
          a <- fit$estimates$ATE; c(a$CI[1], a$CI[2])
        } else if (estimator == "match_tmle") {
          # Paired-difference SE on the matched cohort (mirrors run_match_workflow).
          ps_fml <- stats::reformulate(covariates, response = treatment)
          ps <- stats::predict(
            stats::glm(ps_fml, data = d, family = stats::binomial()),
            type = "response")
          ps <- pmin(pmax(ps, truncate), 1 - truncate)
          A_d  <- d[[treatment]]
          pairs <- .match_nn_pairs(log(ps / (1 - ps)), A_d)
          md    <- d[pairs$match_idx, , drop = FALSE]
          A_m   <- md[[treatment]]; Y_m <- as.numeric(md[[outcome]])
          n_m   <- pairs$n_matched
          Y1 <- Y_m[A_m == 1][seq_len(n_m)]
          Y0 <- Y_m[A_m == 0][seq_len(n_m)]
          ok <- !is.na(Y1) & !is.na(Y0)
          est <- mean(Y1[ok]) - mean(Y0[ok])
          se  <- if (sum(ok) > 1L)
            sqrt(stats::var(Y1[ok] - Y0[ok]) / sum(ok)) else NA_real_
          c(est - z * se, est + z * se)
        } else {
          # IPTW analytic (robust sandwich-free) SE via influence function.
          est <- .rd_point_estimate(d, treatment, outcome, covariates,
                                    "iptw", sl_library, truncate)
          # Conservative IF-style SE from the weighted arm variances.
          A <- d[[treatment]]; Y <- d[[outcome]]
          se <- sqrt(stats::var(Y[A == 1]) / sum(A == 1) +
                     stats::var(Y[A == 0]) / sum(A == 0))
          c(est - z * se, est + z * se)
        }
      }, error = function(e) c(NA, NA))
      covers[i, "influence"] <- as.integer(ci[1] <= truth && truth <= ci[2])
    }
    if ("bootstrap" %in% methods) {
      bs <- tryCatch(
        bootstrap_rd_variance(d, treatment, outcome, covariates,
                              estimator = estimator, R = R_boot,
                              sl_library = sl_library, truncate = truncate,
                              conf_level = target_coverage, seed = 1L),
        error = function(e) NULL)
      covers[i, "bootstrap"] <- if (is.null(bs)) NA_integer_ else
        as.integer(bs$ci[1] <= truth && truth <= bs$ci[2])
    }
  }

  cov_vec <- colMeans(covers, na.rm = TRUE)
  # Prefer methods at or above nominal; among those, the closest to nominal.
  ok <- cov_vec >= target_coverage - 1e-8
  selected <- if (any(ok)) names(which.min(cov_vec[ok] - target_coverage))
              else names(which.max(cov_vec))
  list(selected = selected, coverage = round(cov_vec, 3),
       n_datasets = length(plasmode_datasets))
}
