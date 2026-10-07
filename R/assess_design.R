.estimand_weights <- function(g, A, estimand, band) {
  switch(estimand,
    ATE = A / g + (1 - A) / (1 - g),
    trimmed_ATE = (A / g + (1 - A) / (1 - g)) * (g >= band[1] & g <= band[2]),
    ATT = A + (1 - A) * g / (1 - g),
    ATO = A * (1 - g) + (1 - A) * g)
}

.ess <- function(w) if (sum(w > 0) == 0) 0 else sum(w)^2 / sum(w^2)

# Smallest group whose covariate means the dossier may report.
.min_cell <- 5L

# One removed patient is exposed by the trimmed-population mean (or the kept
# mean) together with the cohort mean and the counts. So when 1 to .min_cell - 1
# patients fall outside the band, every figure that would reveal their values is
# suppressed. With none removed, the trimmed population is the cohort.
.small_removed <- function(n_removed) n_removed >= 1L & n_removed < .min_cell

.std_shift <- function(m, ref, s)
  ifelse(is.na(m), NA_real_, ifelse(is.finite(s) & s > 0, (m - ref) / s, 0))

# Covariate means in each estimand's target population, and the standardized
# shift from the full cohort (difference over the full-cohort SD).
.populations <- function(X, A, g, estimands, band) {
  mu <- colMeans(X)
  s <- apply(X, 2, stats::sd)
  out <- data.frame(covariate = colnames(X), cohort = unname(mu), stringsAsFactors = FALSE)
  for (e in estimands) {
    w <- switch(e,
      ATE = rep(1, nrow(X)),
      trimmed_ATE = as.numeric(g >= band[1] & g <= band[2]),
      ATT = as.numeric(A == 1),
      ATO = g * (1 - g))
    m <- if (sum(w) > 0) colSums(X * w) / sum(w) else rep(NA_real_, ncol(X))
    if (e == "trimmed_ATE" && .small_removed(sum(g < band[1] | g > band[2])))
      m <- rep(NA_real_, ncol(X))
    out[[e]] <- unname(m)
    out[[paste0(e, "_shift")]] <- unname(.std_shift(m, mu, s))
  }
  out
}

# The patients a trimmed analysis no longer describes: counts by arm, and the
# covariate means of removed and kept patients. Means of fewer than .min_cell
# removed patients are suppressed, and so are the kept means and the standardized
# difference when 1 to .min_cell - 1 patients are removed (see .small_removed).
.trimmed_profile <- function(X, A, g, band) {
  out <- g < band[1] | g > band[2]
  s <- apply(X, 2, stats::sd)
  mean_or_na <- function(rows) if (sum(rows) >= .min_cell) colMeans(X[rows, , drop = FALSE]) else
    rep(NA_real_, ncol(X))
  removed <- mean_or_na(out)
  kept <- if (.small_removed(sum(out))) rep(NA_real_, ncol(X)) else mean_or_na(!out)
  list(
    counts = data.frame(group = c("all", "treated", "control"),
                        removed = c(sum(out), sum(out & A == 1), sum(out & A == 0)),
                        kept = c(sum(!out), sum(!out & A == 1), sum(!out & A == 0)),
                        stringsAsFactors = FALSE),
    covariates = data.frame(covariate = colnames(X), removed = unname(removed),
                            kept = unname(kept),
                            std_diff = unname(ifelse(is.na(removed), NA_real_,
                                                     .std_shift(removed, kept, s))),
                            stringsAsFactors = FALSE))
}

#' Assess the design without the outcome
#'
#' Cross-fits the propensity score with the plan's `ps_library`, grades
#' overlap, and reports per-estimand effective sample sizes, covariate
#' balance and the learner edge check. It also fits the same library once
#' on all design rows (no cross-fitting); that single fit, a function of
#' the covariates only, is the propensity score of the world that
#' [simulate_design()] simulates. The dossier carries neither score.
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
  # The simulated world needs a propensity score that is a function of the
  # covariates alone. The cross-fitted g also depends on each row's fold, which
  # is not in X, so drawing treatment from it would plant a hidden confounder.
  g_world <- .fit_g(X, A, p$ps_library, rep(1L, nrow(X)), p$V, p$seed)$g
  band <- p$trim_band %||% c(0.05, 0.95)
  outside <- g < band[1] | g > band[2]
  arm <- function(f) c(f(g[A == 1]), f(g[A == 0]))
  overlap <- data.frame(arm = c("treated", "control"), n = c(sum(A == 1), sum(A == 0)),
                        g_min = arm(min), g_median = arm(stats::median), g_max = arm(max),
                        share_outside_band = c(mean(outside[A == 1]), mean(outside[A == 0])))
  share <- mean(outside)
  grade <- if (share < 0.05) "good" else if (share < 0.20) "moderate" else "poor"
  Xdf <- as.data.frame(X)
  colnames(Xdf) <- make.names(colnames(Xdf), unique = TRUE)  # cobalt parses names as formula terms
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
      d <- abs(bt$Balance$Diff.Adj)
      if (all(is.na(d))) NA_real_ else max(d, na.rm = TRUE)
    }, error = function(err) {
      warning("assess_design: balance for ", e, " could not be computed (",
              conditionMessage(err), ").", call. = FALSE)
      NA_real_
    })
    data.frame(estimand = e, max_abs_smd = smd)
  }))
  populations <- .populations(X, A, g, p$estimands, p$trim_band)
  trimmed_profile <- if (is.null(p$trim_band)) NULL else .trimmed_profile(X, A, g, p$trim_band)
  out <- list(g = g, g_world = g_world, folds = folds, overlap = overlap, overlap_grade = grade,
              share_outside_band = share, band = band, ess = ess, balance = balance,
              populations = populations, trimmed_profile = trimmed_profile,
              edge = .edge_check(gf$risks), risks = gf$risks)
  out$stamp <- .stamp(lock, out)
  class(out) <- "cr_design"
  out
}
