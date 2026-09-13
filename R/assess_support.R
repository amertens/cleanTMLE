# Support assessment: the graded overlap verdict a design team records before
# outcomes are unlocked. assess_support() absorbs compute_ps_diagnostics()'s
# ESS, the c-statistic that run_positivity_diagnostics() intended, the graded
# gate of the Rescue.Co main pipeline, the near-deterministic stratum check,
# and a compact tree search for multivariate violation regions.

#' One Front Door for Propensity-Score Fitting
#'
#' Thin dispatcher over the three ways cleanTMLE obtains a propensity score:
#' SuperLearner ([fit_ps_superlearner()]), logistic regression
#' ([fit_ps_glm()]), and externally supplied scores ([wrap_ps_fit()]).
#'
#' @param lock A `cleanroom_lock`.
#' @param method `"superlearner"` (default), `"glm"`, or `"external"`.
#' @param scores Numeric vector of external scores; required when
#'   `method = "external"`.
#' @param truncate Truncation bound passed through (see the underlying
#'   functions; `NULL` for external scores keeps them as supplied).
#' @param ... Passed to the underlying fitter (for example `cv_folds`,
#'   `cluster` for SuperLearner).
#' @return A `ps_fit` object; the untruncated scores are kept in `$ps_raw`.
#' @export
fit_ps <- function(lock, method = c("superlearner", "glm", "external"),
                   scores = NULL, truncate = 0.01, ...) {
  method <- match.arg(method)
  switch(method,
    superlearner = fit_ps_superlearner(lock, truncate = truncate, ...),
    glm          = fit_ps_glm(lock, truncate = truncate),
    external     = {
      if (is.null(scores))
        stop("method = 'external' needs `scores`.", call. = FALSE)
      wrap_ps_fit(lock, ps_scores = scores, truncate = truncate)
    })
}


#' Prespecified Thresholds for the Graded Support Verdict
#'
#' The defaults are the Rescue.Co main pipeline's, calibrated against a sound
#' reference fit (g in [0.043, 0.937], max weight 7.5) and a broken one
#' (g in [0.004, 0.940], max weight 33.5) that global statistics barely
#' separated. The verdict vocabulary follows the review-team decisions of
#' Muntner et al. (2024): PASS reports plainly, FLAG reports with a caveat,
#' SEVERE reports only beside a common-support estimate, FAIL refuses the
#' estimate as not an estimate of anything.
#'
#' @param band Numeric length 2; the reference propensity band. Default
#'   `c(0.05, 0.95)`.
#' @param flag_pct_outside,severe_pct_outside,fail_pct_outside Percent of the
#'   sample outside `band` that triggers each grade. Defaults 1, 25, 60.
#' @param flag_max_weight,severe_max_weight,fail_max_weight Maximum IPTW
#'   weight that triggers each grade. Defaults 30, 150, 1000.
#' @param nd_min_stratum Minimum size of a binary-covariate stratum for the
#'   near-deterministic check. Default 20.
#' @param nd_min_treated Below this many treated in the stratum, with an
#'   extreme treated share, the stratum is flagged. Default 30.
#' @param nd_extreme_p Treated share below this (or above 1 minus this)
#'   counts as extreme. Default 0.10.
#' @return A `support_thresholds` list.
#' @references Muntner P et al. (2024) Pharmacoepidemiol Drug Saf 33:e5770.
#'   Conover MM, Schuemie MJ et al. (2025) J Am Med Inform Assoc: objective
#'   study validity diagnostics computed while estimates stay blinded, with
#'   failed analyses labelled inestimable.
#' @export
support_thresholds <- function(band = c(0.05, 0.95),
                               flag_pct_outside = 1,
                               flag_max_weight = 30,
                               severe_pct_outside = 25,
                               severe_max_weight = 150,
                               fail_pct_outside = 60,
                               fail_max_weight = 1000,
                               nd_min_stratum = 20L,
                               nd_min_treated = 30L,
                               nd_extreme_p = 0.10) {
  stopifnot(length(band) == 2L, band[1] > 0, band[2] < 1, band[1] < band[2])
  out <- list(band = band,
              flag_pct_outside = flag_pct_outside,
              flag_max_weight = flag_max_weight,
              severe_pct_outside = severe_pct_outside,
              severe_max_weight = severe_max_weight,
              fail_pct_outside = fail_pct_outside,
              fail_max_weight = fail_max_weight,
              nd_min_stratum = nd_min_stratum,
              nd_min_treated = nd_min_treated,
              nd_extreme_p = nd_extreme_p)
  class(out) <- "support_thresholds"
  out
}

#' @export
print.support_thresholds <- function(x, ...) {
  cat("Support-verdict thresholds\n")
  cat(sprintf("  band: [%.2f, %.2f]\n", x$band[1], x$band[2]))
  cat(sprintf("  FLAG   above %g%% outside or max weight %g\n",
              x$flag_pct_outside, x$flag_max_weight))
  cat(sprintf("  SEVERE above %g%% outside or max weight %g\n",
              x$severe_pct_outside, x$severe_max_weight))
  cat(sprintf("  FAIL   above %g%% outside or max weight %g\n",
              x$fail_pct_outside, x$fail_max_weight))
  cat(sprintf("  near-deterministic stratum: n >= %d, treated < %d, share outside [%.2f, %.2f]\n",
              x$nd_min_stratum, x$nd_min_treated, x$nd_extreme_p,
              1 - x$nd_extreme_p))
  invisible(x)
}

.support_caveats <- c(
  PASS   = "Overlap adequate.",
  FLAG   = "Mild positivity violation; interpret with the usual caution.",
  SEVERE = paste("SEVERE positivity violation: a large share of the sample",
                 "has little or no counterpart in the other arm, so this",
                 "estimate relies on extrapolation. Report only alongside",
                 "the common-support estimate for the same contrast."),
  FAIL   = paste("EXTREME positivity violation: arms do not share a common",
                 "covariate space. Diagnostic only; not an estimate."))

.support_verdict <- function(pct_outside, max_weight, th) {
  if (pct_outside > th$fail_pct_outside || max_weight > th$fail_max_weight)
    "FAIL"
  else if (pct_outside > th$severe_pct_outside ||
           max_weight > th$severe_max_weight) "SEVERE"
  else if (pct_outside > th$flag_pct_outside ||
           max_weight > th$flag_max_weight) "FLAG"
  else "PASS"
}

# Near-deterministic strata over binary covariates: strata in which treatment
# is very nearly determined, however healthy the global weights look.
.near_deterministic <- function(X, A, min_stratum = 20L, min_treated = 30L,
                                extreme_p = 0.10) {
  rows <- list()
  for (v in colnames(X)) {
    x <- X[[v]]
    u <- unique(x[!is.na(x)])
    if (length(u) != 2L) next
    on <- !is.na(x) & x == max(u)
    if (sum(on) < min_stratum) next
    p  <- mean(A[on]); nt <- sum(A[on] == 1)
    if (nt < min_treated && (p < extreme_p || p > 1 - extreme_p))
      rows[[length(rows) + 1L]] <- data.frame(
        covariate = v, n_in_stratum = sum(on), n_treated_in_stratum = nt,
        p_treated_in_stratum = round(p, 4), stringsAsFactors = FALSE)
  }
  if (!length(rows)) return(NULL)
  r <- do.call(rbind, rows)
  r[order(r$n_treated_in_stratum), , drop = FALSE]
}

# Compact tree search for multivariate violation regions, in the spirit of
# the PoRT algorithm (Danelian et al. 2023) and the discriminative-tree
# approach of Karavani et al. (2019): fit a shallow classification tree of A
# on W and report terminal nodes with an extreme treated share.
.violation_tree <- function(X, A, max_depth = 3L, min_n = 50L,
                            extreme_p = 0.10) {
  if (!requireNamespace("rpart", quietly = TRUE)) return(NULL)
  df <- data.frame(.A = A, X)
  fit <- tryCatch(rpart::rpart(
    .A ~ ., data = df, method = "class",
    control = rpart::rpart.control(maxdepth = max_depth, minbucket = min_n,
                                   cp = 0.001)),
    error = function(e) NULL)
  if (is.null(fit) || is.null(fit$frame) || nrow(fit$frame) < 2L) return(NULL)
  fr <- fit$frame
  leaves <- which(fr$var == "<leaf>")
  if (!length(leaves)) return(NULL)
  where_node <- as.integer(rownames(fr))[fit$where]
  rows <- list()
  for (li in leaves) {
    node_id <- as.integer(rownames(fr))[li]
    in_leaf <- where_node == node_id
    n_leaf  <- sum(in_leaf)
    if (n_leaf < min_n) next
    p <- mean(A[in_leaf])
    if (p >= extreme_p && p <= 1 - extreme_p) next
    pth <- tryCatch(rpart::path.rpart(fit, nodes = node_id, print.it = FALSE),
                    error = function(e) NULL)
    rule <- if (is.null(pth)) NA_character_ else
      paste(pth[[1]][-1], collapse = " & ")
    rows[[length(rows) + 1L]] <- data.frame(
      subgroup = rule, n = n_leaf, n_treated = sum(A[in_leaf]),
      p_treated = round(p, 4), stringsAsFactors = FALSE)
  }
  if (!length(rows)) return(NULL)
  do.call(rbind, rows)
}


#' Assess the Support a Design Gives Each Arm
#'
#' The design team's overlap verdict, computed from a fitted propensity score
#' with no outcome access. Returns the percent of the sample outside a
#' prespecified propensity band, the maximum and 99th percentile IPTW weight,
#' effective sample size by arm, the propensity c-statistic, a graded verdict
#' (PASS, FLAG, SEVERE, FAIL) with a caveat string that travels with every
#' downstream estimate, a table of near-deterministic strata over binary
#' covariates, and an optional shallow-tree search for multivariate violation
#' regions. Conditional non-overlap escalates PASS or FLAG to SEVERE, never
#' to FAIL, so it can strengthen a warning but never suppress a number.
#'
#' Diagnostics are computed on the untruncated propensity scores
#' (`ps_fit$ps_raw` when present) because truncation caps the very weights
#' whose size is the diagnostic.
#'
#' @param ps_fit A `ps_fit` from [fit_ps()], [fit_ps_superlearner()],
#'   [fit_ps_glm()], or [wrap_ps_fit()].
#' @param thresholds A [support_thresholds()] object.
#' @param tree_search Logical; run the shallow-tree violation search when
#'   rpart is installed. Default TRUE.
#' @param tree_max_depth,tree_min_n Tree depth and minimum leaf size.
#'
#' @return An object of class `support_assessment` with fields `summary`
#'   (one row: n, arms, g range, counts and percent outside the band, weight
#'   quantiles, ESS), `verdict`, `caveat`, `escalated` (logical, with
#'   `verdict_global` the pre-escalation grade), `c_statistic`,
#'   `near_deterministic`, `violation_regions`, `thresholds`, and the scores
#'   needed by `plot()`.
#'
#' @references Muntner P et al. (2024) Pharmacoepidemiol Drug Saf 33:e5770.
#'   Conover MM et al. (2025) JAMIA (blinded validity diagnostics; failed
#'   analyses labelled inestimable). Danelian G et al. (2023) J Causal
#'   Inference 11:20220032 (PoRT). Karavani E et al. (2019) arXiv:1907.08127.
#'
#' @examples
#' \dontrun{
#' lock <- create_simple_lock(sim_func1(500), "treatment", "event_24",
#'                            c("age", "sex", "biomarker"))
#' sup <- assess_support(fit_ps(lock, "glm"))
#' print(sup)
#' plot(sup)
#' }
#' @export
assess_support <- function(ps_fit,
                           thresholds = support_thresholds(),
                           tree_search = TRUE,
                           tree_max_depth = 3L,
                           tree_min_n = 50L) {
  if (!inherits(ps_fit, "ps_fit"))
    stop("`ps_fit` must be a ps_fit object (fit_ps(), fit_ps_superlearner(), ",
         "fit_ps_glm(), or wrap_ps_fit()).", call. = FALSE)
  if (!inherits(thresholds, "support_thresholds"))
    stop("`thresholds` must come from support_thresholds().", call. = FALSE)

  data <- ps_fit$data
  A <- as.integer(data[[ps_fit$treatment]])
  g <- as.numeric(ps_fit$ps_raw %||% ps_fit$ps)
  g <- pmin(pmax(g, 1e-6), 1 - 1e-6)
  band <- thresholds$band
  n <- length(A)

  w <- ifelse(A == 1, 1 / g, 1 / (1 - g))
  ess_t <- sum(w[A == 1])^2 / sum(w[A == 1]^2)
  ess_c <- sum(w[A == 0])^2 / sum(w[A == 0]^2)

  # Propensity c-statistic (Mann-Whitney form).
  c_stat <- {
    n1 <- sum(A == 1); n0 <- sum(A == 0)
    if (n1 == 0 || n0 == 0) NA_real_
    else (sum(rank(g)[A == 1]) - n1 * (n1 + 1) / 2) / (n1 * n0)
  }

  smry <- data.frame(
    n = n, n_treated = sum(A == 1), n_control = sum(A == 0),
    min_g = round(min(g), 5), max_g = round(max(g), 5),
    n_below_band = sum(g < band[1]), n_above_band = sum(g > band[2]),
    pct_outside_band = round(100 * mean(g < band[1] | g > band[2]), 2),
    max_iptw_weight = round(max(w), 1),
    p99_iptw_weight = round(stats::quantile(w, 0.99, names = FALSE), 1),
    mean_iptw_weight = round(mean(w), 2),
    ess_treated = round(ess_t, 1), ess_control = round(ess_c, 1),
    c_statistic = round(c_stat, 3),
    stringsAsFactors = FALSE)

  verdict_global <- .support_verdict(smry$pct_outside_band,
                                     smry$max_iptw_weight, thresholds)

  W <- data[, ps_fit$covariates, drop = FALSE]
  Wn <- as.data.frame(lapply(W, function(x)
    if (is.numeric(x)) x else as.numeric(as.factor(x))))
  nd <- .near_deterministic(Wn, A,
                            min_stratum = thresholds$nd_min_stratum,
                            min_treated = thresholds$nd_min_treated,
                            extreme_p = thresholds$nd_extreme_p)
  vr <- if (isTRUE(tree_search))
    .violation_tree(Wn, A, max_depth = tree_max_depth, min_n = tree_min_n,
                    extreme_p = thresholds$nd_extreme_p) else NULL

  # Conditional non-overlap escalates PASS or FLAG to SEVERE, never to FAIL:
  # global statistics can look acceptable while strata are near-deterministic
  # in treatment, and that is exactly what the stratum check measures.
  verdict <- verdict_global
  escalated <- FALSE
  if (!is.null(nd) && verdict %in% c("PASS", "FLAG")) {
    verdict <- "SEVERE"
    escalated <- TRUE
  }

  out <- list(
    summary = cbind(smry, verdict = verdict,
                    caveat = unname(.support_caveats[verdict]),
                    stringsAsFactors = FALSE),
    verdict = verdict,
    verdict_global = verdict_global,
    escalated = escalated,
    caveat = unname(.support_caveats[verdict]),
    c_statistic = c_stat,
    near_deterministic = nd,
    violation_regions = vr,
    thresholds = thresholds,
    band = band,
    g = g,
    A = A,
    treatment = ps_fit$treatment,
    method = ps_fit$method,
    call = match.call()
  )
  class(out) <- "support_assessment"
  out
}

#' @export
print.support_assessment <- function(x, ...) {
  s <- x$summary
  cat("Support assessment (", x$method %||% "propensity", " fit)\n", sep = "")
  cat(sprintf("  n = %d (%d treated, %d control); g in [%.4f, %.4f]\n",
              s$n, s$n_treated, s$n_control, s$min_g, s$max_g))
  cat(sprintf("  %.2f%% outside [%.2f, %.2f]; max weight %.1f (99th pct %.1f); ESS %d / %d; c-statistic %.3f\n",
              s$pct_outside_band, x$band[1], x$band[2], s$max_iptw_weight,
              s$p99_iptw_weight, round(s$ess_treated), round(s$ess_control),
              s$c_statistic))
  if (isTRUE(x$escalated)) {
    cat(sprintf("  Verdict: %s (escalated from %s: %d near-deterministic stratum/strata)\n",
                x$verdict, x$verdict_global,
                nrow(x$near_deterministic)))
  } else {
    cat(sprintf("  Verdict: %s\n", x$verdict))
  }
  cat("  ", x$caveat, "\n", sep = "")
  if (!is.null(x$near_deterministic)) {
    cat("\n  Near-deterministic strata (worst first):\n")
    print(utils::head(x$near_deterministic, 5), row.names = FALSE)
  }
  if (!is.null(x$violation_regions)) {
    cat("\n  Tree-identified violation regions:\n")
    print(utils::head(x$violation_regions, 5), row.names = FALSE)
  }
  invisible(x)
}

#' @export
plot.support_assessment <- function(x, bins = 40, ...) {
  if (!requireNamespace("ggplot2", quietly = TRUE))
    stop("ggplot2 is required.", call. = FALSE)
  df <- data.frame(g = x$g, arm = factor(ifelse(x$A == 1, "treated",
                                                "control"),
                                         levels = c("treated", "control")))
  brks <- seq(max(min(df$g) - 1e-9, 0), min(max(df$g) + 1e-9, 1),
              length.out = bins + 1L)
  ht <- hist(df$g[df$arm == "treated"], breaks = brks, plot = FALSE)
  hc <- hist(df$g[df$arm == "control"], breaks = brks, plot = FALSE)
  mir <- rbind(
    data.frame(mid = ht$mids, count = ht$counts, arm = "treated"),
    data.frame(mid = hc$mids, count = -hc$counts, arm = "control"))
  ggplot2::ggplot(mir, ggplot2::aes(x = .data$mid, y = .data$count,
                                    fill = .data$arm)) +
    ggplot2::geom_col(width = diff(brks)[1], alpha = 0.85) +
    ggplot2::geom_hline(yintercept = 0, linewidth = 0.3) +
    ggplot2::geom_vline(xintercept = x$band, linetype = "dashed") +
    ggplot2::annotate("rect", xmin = -Inf, xmax = x$band[1], ymin = -Inf,
                      ymax = Inf, alpha = 0.08, fill = "red") +
    ggplot2::annotate("rect", xmin = x$band[2], xmax = Inf, ymin = -Inf,
                      ymax = Inf, alpha = 0.08, fill = "red") +
    ggplot2::scale_fill_manual(values = c(treated = "#1f77b4",
                                          control = "#8c9dab")) +
    ggplot2::labs(x = "propensity score", y = "count (control mirrored)",
                  title = sprintf("Propensity support by arm: %s", x$verdict),
                  subtitle = sprintf("%.1f%% outside [%.2f, %.2f]; max weight %.0f",
                                     x$summary$pct_outside_band, x$band[1],
                                     x$band[2], x$summary$max_iptw_weight),
                  fill = NULL) +
    ggplot2::theme_minimal()
}
