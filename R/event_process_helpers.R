# Weight diagnostics and the pre-outcome weight checkpoint.
#
# The event-process classification, target-population and missing-data
# declarations, competing-risk coherence check, and cumulative-risk reporting
# table are pure note-taking formatters with no estimation coupling; they were
# split out into the companion `cleanroomGov` package. What remains here feeds
# the gate: `clean_weight_diagnostics()` and the `checkpoint_weights()` gate
# checkpoint that consumes it.

# ---------------------------------------------------------------------
# Weight diagnostics
# ---------------------------------------------------------------------

#' Weight Diagnostics for Treatment, Censoring, or Missingness Weights
#'
#' Standard diagnostic summary for inverse-probability weights:
#' percentiles by treatment group and overall, effective sample size,
#' counts and proportions of extreme weights under user-specified
#' thresholds, and (when covariates are supplied) standardised mean
#' differences before and after weighting. The summary supports
#' transparent reporting of overlap and weight-distribution properties;
#' it does not itself verify positivity.
#'
#' @section Clean-room stage: Stage 2 / 4 (diagnostic).
#'
#' @param weights Numeric vector of inverse-probability weights.
#' @param treatment Optional treatment indicator (0/1 or factor) of the
#'   same length as \code{weights}. When supplied, percentiles and ESS
#'   are reported by treatment group as well as overall.
#' @param covariates Optional data frame of baseline covariates; if
#'   supplied, the function returns standardised mean differences before
#'   and after weighting for each covariate.
#' @param max_weight_threshold Numeric threshold; an observation whose
#'   weight exceeds this value is flagged as extreme. Default
#'   \code{10}.
#' @param ess_floor Numeric threshold; the function flags ESS values
#'   below this as a low effective sample size. Default \code{0.3 *
#'   length(weights)}.
#'
#' @return A list of class \code{cleantmle_weight_diag} with elements
#'   \code{percentiles}, \code{ess}, \code{extreme_weights},
#'   \code{flags}, and (when covariates are supplied) \code{smd}.
#'
#' @examples
#' set.seed(1)
#' n <- 400
#' A <- rbinom(n, 1, 0.4)
#' ps <- plogis(0.2 * rnorm(n) + 0.5 * A)
#' w  <- ifelse(A == 1, 1 / ps, 1 / (1 - ps))
#' clean_weight_diagnostics(w, treatment = A)
#'
#' @export
clean_weight_diagnostics <- function(weights, treatment = NULL,
                                      covariates = NULL,
                                      max_weight_threshold = 10,
                                      ess_floor = NULL) {
  if (!is.numeric(weights) || length(weights) == 0L)
    stop("`weights` must be a non-empty numeric vector.", call. = FALSE)
  n <- length(weights)
  if (is.null(ess_floor)) ess_floor <- 0.3 * n
  if (!is.null(treatment) && length(treatment) != n)
    stop("`treatment` must be the same length as `weights`.", call. = FALSE)

  pct_probs <- c(min = 0, p01 = 0.01, p05 = 0.05, p25 = 0.25,
                 median = 0.5, p75 = 0.75, p95 = 0.95, p99 = 0.99,
                 max = 1)

  pct_row <- function(w) {
    q <- stats::quantile(w, probs = pct_probs, na.rm = TRUE, names = FALSE)
    setNames(c(q, mean(w, na.rm = TRUE)),
             c(names(pct_probs), "mean"))
  }

  ess_fn <- function(w) sum(w)^2 / sum(w^2)

  pct_overall <- pct_row(weights)
  ess_overall <- ess_fn(weights)
  pct_by <- NULL; ess_by <- NULL
  if (!is.null(treatment)) {
    pct_by <- t(sapply(unique(treatment),
                       function(a) pct_row(weights[treatment == a])))
    rownames(pct_by) <- as.character(unique(treatment))
    ess_by <- vapply(unique(treatment),
                     function(a) ess_fn(weights[treatment == a]),
                     numeric(1))
    names(ess_by) <- as.character(unique(treatment))
  }

  extreme_n <- sum(weights > max_weight_threshold, na.rm = TRUE)
  extreme_p <- extreme_n / n

  smd_tbl <- NULL
  if (!is.null(covariates)) {
    if (is.null(treatment))
      stop("Provide `treatment` to compute SMDs.", call. = FALSE)
    if (nrow(covariates) != n)
      stop("`covariates` must have the same number of rows as `weights`.",
           call. = FALSE)
    smd_one <- function(x, A, w = NULL) {
      x <- as.numeric(x)
      if (is.null(w)) {
        m1 <- mean(x[A == 1], na.rm = TRUE)
        m0 <- mean(x[A == 0], na.rm = TRUE)
        s1 <- stats::var(x[A == 1], na.rm = TRUE)
        s0 <- stats::var(x[A == 0], na.rm = TRUE)
      } else {
        wm <- function(z, ww) sum(z * ww, na.rm = TRUE) /
                                sum(ww[!is.na(z)], na.rm = TRUE)
        m1 <- wm(x[A == 1], w[A == 1]); m0 <- wm(x[A == 0], w[A == 0])
        s1 <- wm((x[A == 1] - m1)^2, w[A == 1])
        s0 <- wm((x[A == 0] - m0)^2, w[A == 0])
      }
      pooled <- sqrt((s1 + s0) / 2)
      if (is.na(pooled) || pooled == 0) return(0)
      abs(m1 - m0) / pooled
    }
    smd_tbl <- data.frame(
      covariate    = names(covariates),
      smd_unweight = vapply(names(covariates),
                            function(v) smd_one(covariates[[v]], treatment),
                            numeric(1)),
      smd_weighted = vapply(names(covariates),
                            function(v) smd_one(covariates[[v]],
                                                treatment, weights),
                            numeric(1)),
      stringsAsFactors = FALSE)
    rownames(smd_tbl) <- NULL
  }

  flags <- list(
    low_ess          = ess_overall < ess_floor,
    extreme_weights  = extreme_n > 0,
    high_max_weight  = max(weights, na.rm = TRUE) > max_weight_threshold)

  out <- list(
    percentiles      = list(overall = pct_overall, by_treatment = pct_by),
    ess              = list(overall = ess_overall, by_treatment = ess_by),
    extreme_weights  = list(threshold = max_weight_threshold,
                             n = extreme_n, prop = extreme_p),
    smd              = smd_tbl,
    flags            = flags,
    n                = n)
  class(out) <- c("cleantmle_weight_diag", class(out))
  out
}

#' @export
print.cleantmle_weight_diag <- function(x, ...) {
  cat("Weight diagnostics\n")
  cat("==================\n")
  cat(sprintf("N = %d\n", x$n))
  cat(sprintf("Effective sample size (overall): %.1f\n", x$ess$overall))
  if (!is.null(x$ess$by_treatment))
    print(round(x$ess$by_treatment, 1))
  cat(sprintf("Extreme weights (> %.2f): %d (%.1f%%)\n",
              x$extreme_weights$threshold,
              x$extreme_weights$n, 100 * x$extreme_weights$prop))
  cat("\nPercentiles (overall):\n")
  print(round(x$percentiles$overall, 4))
  if (!is.null(x$smd)) {
    cat("\nStandardised mean differences:\n")
    print(round(x$smd[, c("smd_unweight", "smd_weighted"), drop = FALSE], 3))
  }
  if (any(unlist(x$flags))) {
    cat("\nFlags: ")
    cat(paste(names(x$flags)[unlist(x$flags)], collapse = "; "), "\n")
  }
  invisible(x)
}


# ---------------------------------------------------------------------
# Weight-diagnostics checkpoint that plugs into the gate
# ---------------------------------------------------------------------

#' Pre-Outcome Weight Checkpoint
#'
#' Wraps [clean_weight_diagnostics()] in a structured checkpoint object
#' compatible with [gate_all()] and [authorize_outcome_analysis()]. Useful
#' when the weights themselves (treatment, censoring, or their product)
#' should be a formal pre-outcome diagnostic, not just a side note.
#'
#' Decision logic:
#' \itemize{
#'   \item \code{STOP} when ESS falls below \code{ess_floor} or the
#'     maximum weight exceeds \code{max_weight_threshold} *and* the
#'     proportion of extreme weights exceeds \code{extreme_prop_threshold};
#'   \item \code{FLAG} when either ESS or extreme-weight criteria fail
#'     individually but not jointly;
#'   \item \code{GO} otherwise.
#' }
#' Thresholds are user-supplied and should be prespecified.
#'
#' @section Clean-room stage: Stage 2 (treatment weights) or Stage 4
#'   pre-flight (censoring or missingness weights).
#'
#' @param weights Numeric vector of inverse-probability weights.
#' @param treatment Optional treatment indicator (0/1 or factor) of the
#'   same length as \code{weights}.
#' @param weight_type Character label for the weight family (e.g.
#'   \code{"treatment"}, \code{"censoring"}, \code{"missingness"},
#'   \code{"treatment x censoring"}). Default \code{"treatment"}.
#' @param max_weight_threshold,ess_floor See [clean_weight_diagnostics()].
#' @param extreme_prop_threshold Maximum acceptable proportion of
#'   observations with weight greater than \code{max_weight_threshold}.
#'   Default \code{0.01}.
#' @param lock_hash Optional character; lock fingerprint to attach to
#'   the checkpoint for traceability.
#'
#' @return A \code{cleantmle_checkpoint} suitable for
#'   [gate_all()] / [record_checkpoint()].
#'
#' @examples
#' set.seed(2)
#' n <- 400
#' A <- rbinom(n, 1, 0.4)
#' ps <- plogis(0.2 * rnorm(n))
#' w  <- ifelse(A == 1, 1 / ps, 1 / (1 - ps))
#' cp <- checkpoint_weights(w, treatment = A,
#'                          max_weight_threshold = 10,
#'                          ess_floor = 0.4 * n)
#' print(cp)
#'
#' @export
checkpoint_weights <- function(weights, treatment = NULL,
                                weight_type = "treatment",
                                max_weight_threshold = 10,
                                ess_floor = NULL,
                                extreme_prop_threshold = 0.01,
                                lock_hash = NA_character_) {
  diag <- clean_weight_diagnostics(weights = weights,
                                    treatment = treatment,
                                    max_weight_threshold = max_weight_threshold,
                                    ess_floor = ess_floor)

  ess <- diag$ess$overall
  ess_fail <- diag$flags$low_ess
  ext_n <- diag$extreme_weights$n
  ext_p <- diag$extreme_weights$prop
  max_w <- diag$percentiles$overall[["max"]]

  ext_fail <- ext_p > extreme_prop_threshold
  decision <- if (ess_fail && ext_fail) "STOP"
              else if (ess_fail || ext_fail) "FLAG"
              else "GO"

  metrics <- data.frame(
    weight_type             = weight_type,
    n                       = diag$n,
    ess                     = round(ess, 1),
    max_weight              = round(max_w, 3),
    pct_extreme             = round(100 * ext_p, 2),
    max_weight_threshold    = max_weight_threshold,
    extreme_prop_threshold  = extreme_prop_threshold,
    stringsAsFactors        = FALSE)

  rationale <- sprintf(paste0(
    "Weight checkpoint (%s): ESS = %.1f (floor %.1f); ",
    "max weight = %.2f (threshold %.1f); ",
    "%% extreme = %.2f%% (threshold %.2f%%)."),
    weight_type, ess,
    if (is.null(ess_floor)) 0.3 * diag$n else ess_floor,
    max_w, max_weight_threshold,
    100 * ext_p, 100 * extreme_prop_threshold)

  cp <- new_checkpoint(
    stage      = paste0("Weight checkpoint (", weight_type, ")"),
    decision   = decision,
    metrics    = metrics,
    thresholds = list(max_weight_threshold = max_weight_threshold,
                      ess_floor = if (is.null(ess_floor)) 0.3 * diag$n
                                  else ess_floor,
                      extreme_prop_threshold = extreme_prop_threshold),
    rationale  = rationale,
    lock_hash  = lock_hash)
  cp$diagnostics <- diag
  cp
}
