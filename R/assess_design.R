.estimand_weights <- function(g, A, estimand, band) {
  switch(estimand,
    ATE = A / g + (1 - A) / (1 - g),
    trimmed_ATE = (A / g + (1 - A) / (1 - g)) * (g >= band[1] & g <= band[2]),
    ATT = A + (1 - A) * g / (1 - g),
    ATO = A * (1 - g) + (1 - A) * g)
}

.ess <- function(w) if (sum(w > 0) == 0) 0 else sum(w)^2 / sum(w^2)

#' Assess the design without the outcome
#'
#' Cross-fits the propensity score with the plan's `ps_library`, grades
#' overlap, and reports per-estimand effective sample sizes, covariate
#' balance and the learner edge check.
#'
#' @param lock A `cr_lock`.
#' @return A `cr_design`.
#' @export
assess_design <- function(lock) {
  verify_lock(lock)
  p <- lock$plan
  X <- .design_matrix(lock$data, lock$covariates)
  A <- as.integer(lock$data[[lock$treatment]])
  folds <- .make_folds(nrow(X), p$K, p$seed)
  gf <- .fit_g(X, A, p$ps_library, folds, p$V, p$seed)
  g <- gf$g
  band <- p$trim_band %||% c(0.05, 0.95)
  outside <- g < band[1] | g > band[2]
  arm <- function(f) c(f(g[A == 1]), f(g[A == 0]))
  overlap <- data.frame(arm = c("treated", "control"), n = c(sum(A == 1), sum(A == 0)),
                        g_min = arm(min), g_median = arm(stats::median), g_max = arm(max),
                        share_outside_band = c(mean(outside[A == 1]), mean(outside[A == 0])))
  share <- mean(outside)
  grade <- if (share < 0.05) "good" else if (share < 0.20) "moderate" else "poor"
  Xdf <- as.data.frame(X)
  ess <- do.call(rbind, lapply(p$estimands, function(e) {
    w <- .estimand_weights(g, A, e, band)
    data.frame(estimand = e, ess_treated = .ess(w[A == 1]), ess_control = .ess(w[A == 0]))
  }))
  balance <- do.call(rbind, lapply(p$estimands, function(e) {
    w <- .estimand_weights(g, A, e, band)
    keep <- w > 0
    smd <- tryCatch({
      bt <- cobalt::bal.tab(Xdf[keep, , drop = FALSE], treat = A[keep], weights = w[keep],
                            s.d.denom = "pooled", binary = "std")
      max(abs(bt$Balance$Diff.Adj), na.rm = TRUE)
    }, error = function(err) NA_real_)
    data.frame(estimand = e, max_abs_smd = smd)
  }))
  out <- list(g = g, folds = folds, overlap = overlap, overlap_grade = grade,
              share_outside_band = share, band = band, ess = ess, balance = balance,
              edge = .edge_check(gf$risks), risks = gf$risks)
  out$stamp <- .stamp(lock, out)
  class(out) <- "cr_design"
  out
}
