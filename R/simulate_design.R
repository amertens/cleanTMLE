.std <- function(v) {
  s <- stats::sd(v)
  if (!is.finite(s) || s == 0) v * 0 else (v - mean(v)) / s
}

.surfaces <- function(plan, X, g, seed) {
  s <- plan$surfaces
  lg <- stats::qlogis(.bound(g, 1e-6))
  Z <- apply(X, 2, .std)
  if (is.null(dim(Z))) Z <- matrix(Z, nrow = nrow(X))
  dirn <- withr::with_seed(seed, stats::rnorm(ncol(Z)))
  lin <- .std(0.5 * .std(lg) + 0.5 * .std(drop(Z %*% dirn)))
  idx <- if (!is.null(s$q0)) {
    list(steward = .std(stats::qlogis(.bound(s$q0, 1e-6))))
  } else {
    ord <- order(-abs(dirn))
    j1 <- ord[1]
    j2 <- if (length(ord) > 1L) ord[2] else j1
    list(linear = lin, nonlinear = .std(lin + Z[, j1]^2 + Z[, j1] * Z[, j2]))[s$forms]
  }
  hg <- .std(lg)
  out <- list()
  for (form in names(idx)) for (delta in unique(c(0, s$heterogeneity))) {
    m <- idx[[form]]
    eta <- function(b0, a) b0 + m + a * (s$log_or + delta * hg)
    b0 <- stats::uniroot(function(b) mean(stats::plogis(eta(b, 0))) - s$baseline_risk,
                         c(-30, 30), tol = 1e-10)$root
    out[[sprintf("%s, heterogeneity %s", form, delta)]] <-
      list(q1 = stats::plogis(eta(b0, 1)), q0 = stats::plogis(eta(b0, 0)))
  }
  for (nm in names(s$custom))
    out[[nm]] <- list(q1 = .bound(s$custom[[nm]](1, X), 1e-6),
                      q0 = .bound(s$custom[[nm]](0, X), 1e-6))
  out
}

.true_values <- function(q1, q0, g, band) {
  tau <- q1 - q0
  h <- g * (1 - g)
  out <- c(ATE = mean(tau), ATT = sum(g * tau) / sum(g), ATO = sum(h * tau) / sum(h))
  if (!is.null(band)) {
    k <- g >= band[1] & g <= band[2]
    out["trimmed_ATE"] <- if (any(k)) mean(tau[k]) else NA_real_
  }
  out
}

.failed_rows <- function(plan, msg) {
  g <- expand.grid(truncation = plan$candidates$truncation, estimand = plan$estimands,
                   stringsAsFactors = FALSE)
  do.call(rbind, lapply(seq_len(nrow(g)), function(i)
    .estimate_row(g$estimand[i], g$truncation[i], message = msg)))
}

.one_rep <- function(r, X, g, surfaces, plan) {
  n <- nrow(X)
  seed <- plan$seed + r
  draw <- withr::with_seed(seed, {
    idx <- sample.int(n, n, replace = TRUE)
    a <- stats::rbinom(n, 1L, g[idx])
    ys <- lapply(surfaces, function(s)
      stats::rbinom(n, 1L, ifelse(a == 1L, s$q1[idx], s$q0[idx])))
    list(idx = idx, a = a, ys = ys)
  })
  Xs <- X[draw$idx, , drop = FALSE]
  folds <- .make_folds(n, plan$K, seed)
  g_fit <- tryCatch(.fit_g(Xs, draw$a, plan$ps_library, folds, plan$V, seed),
                    error = function(e) conditionMessage(e))
  rows <- list()
  risks <- list()
  for (sn in names(surfaces)) for (lib in names(plan$candidates$library)) {
    res <- if (is.character(g_fit)) .failed_rows(plan, g_fit) else tryCatch(
      .fit_and_target(Xs, draw$a, draw$ys[[sn]], plan$candidates$library[[lib]],
                      plan$estimands, plan$candidates$truncation, plan, folds, seed, g_fit),
      error = function(e) .failed_rows(plan, conditionMessage(e)))
    rk <- attr(res, "risks")
    if (!is.null(rk)) {
      rk <- rk[rk$nuisance == "Q", , drop = FALSE]
      if (nrow(rk)) risks[[length(risks) + 1L]] <- cbind(rep = r, surface = sn, library = lib, rk)
    }
    rows[[length(rows) + 1L]] <- cbind(rep = r, surface = sn, library = lib, res)
  }
  list(results = do.call(rbind, rows),
       risks = if (length(risks)) do.call(rbind, risks) else NULL)
}

.sim_metrics <- function(results, truths, plan) {
  key_t <- paste(truths$surface, truths$estimand)
  results$truth <- truths$truth[match(paste(results$surface, results$estimand), key_t)]
  results$candidate <- sprintf("%s_t%s", results$library, results$truncation)
  keys <- unique(results[c("estimand", "candidate", "library", "truncation", "surface")])
  m <- do.call(rbind, lapply(seq_len(nrow(keys)), function(i) {
    k <- keys[i, ]
    s <- results[results$estimand == k$estimand & results$candidate == k$candidate &
                   results$surface == k$surface, , drop = FALSE]
    ok <- !s$failed & is.finite(s$estimate)
    R <- sum(ok)
    e <- s$estimate[ok] - s$truth[ok]
    cv <- s$ci_lower[ok] <= s$truth[ok] & s$truth[ok] <= s$ci_upper[ok]
    bias <- if (R) mean(e) else NA_real_
    coverage <- if (R) mean(cv) else NA_real_
    rmse <- if (R) sqrt(mean(e^2)) else NA_real_
    data.frame(k, reps = nrow(s), failed = nrow(s) - R, fail_rate = 1 - R / nrow(s),
               bias = bias, bias_mcse = if (R > 1) stats::sd(e) / sqrt(R) else NA_real_,
               coverage = coverage,
               coverage_mcse = if (R > 1) sqrt(coverage * (1 - coverage) / R) else NA_real_,
               rmse = rmse,
               rmse_mcse = if (R > 1 && rmse > 0) stats::sd(e^2) / (2 * rmse * sqrt(R)) else NA_real_,
               stringsAsFactors = FALSE)
  }))
  tb <- plan$tolerance$bias
  tc <- plan$tolerance$coverage
  m$pass <- !is.na(m$bias) & m$reps - m$failed > 1 & m$fail_rate <= 0.05 &
    abs(m$bias) <= tb & m$coverage >= tc
  m$borderline <- (!is.na(m$bias_mcse) & abs(abs(m$bias) - tb) <= 2 * m$bias_mcse) |
    (!is.na(m$coverage_mcse) & abs(m$coverage - tc) <= 2 * m$coverage_mcse)
  rownames(m) <- NULL
  m
}

.sim_verdict <- function(metrics, plan) {
  do.call(rbind, lapply(plan$estimands, function(e) {
    me <- metrics[metrics$estimand == e, , drop = FALSE]
    cs <- do.call(rbind, lapply(split(me, me$candidate), function(d) data.frame(
      candidate = d$candidate[1], library = d$library[1], truncation = d$truncation[1],
      pass_all = all(d$pass), borderline_any = any(d$borderline),
      worst_rmse = max(d$rmse), stringsAsFactors = FALSE)))
    clean <- cs[cs$pass_all & !cs$borderline_any, , drop = FALSE]
    border <- cs[cs$pass_all & cs$borderline_any, , drop = FALSE]
    pool <- if (nrow(clean)) clean else border
    status <- if (nrow(clean)) "feasible" else if (nrow(border)) "borderline" else "infeasible"
    sel <- NA_character_
    if (nrow(pool)) {
      sel <- if (is.null(plan$select)) pool$candidate[which.min(pool$worst_rmse)] else
        plan$select(me[me$candidate %in% pool$candidate, , drop = FALSE])
      if (!sel %in% pool$candidate)
        stop("`select` must return one of: ", paste(pool$candidate, collapse = ", "),
             call. = FALSE)
    }
    row <- pool[match(sel, pool$candidate), , drop = FALSE]
    data.frame(estimand = e, status = status, selected = sel,
               library = if (is.na(sel)) NA_character_ else row$library,
               truncation = if (is.na(sel)) NA_real_ else row$truncation,
               stringsAsFactors = FALSE)
  }))
}

#' Outcome-free plasmode simulation: the support map and candidate selection
#'
#' Resamples the design rows, draws treatment from the design's propensity
#' score and outcomes from the plan's family of surfaces, and runs every
#' candidate with the same fitter that `estimate_effect()` uses. An
#' estimand is feasible when some candidate meets the bias and coverage
#' tolerances on every surface. Repetitions run in parallel through
#' `future.apply` when a parallel `future` plan is set and cleanTMLE is
#' installed (not only loaded with `devtools::load_all()`).
#'
#' @param lock A `cr_lock`.
#' @param design The [assess_design()] result for `lock`.
#' @return A `cr_simulation`.
#' @export
simulate_design <- function(lock, design) {
  verify_lock(lock)
  .check_stamp(lock, design, "design")
  p <- lock$plan
  X <- .design_matrix(lock$data, lock$covariates)
  g <- design$g
  surfaces <- .surfaces(p, X, g, p$seed)
  truths <- do.call(rbind, lapply(names(surfaces), function(sn) {
    tv <- .true_values(surfaces[[sn]]$q1, surfaces[[sn]]$q0, g, p$trim_band)
    data.frame(surface = sn, estimand = names(tv), truth = unname(tv),
               stringsAsFactors = FALSE)
  }))
  t0 <- Sys.time()
  first <- .one_rep(1L, X, g, surfaces, p)
  secs <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
  par <- requireNamespace("future.apply", quietly = TRUE) &&
    requireNamespace("future", quietly = TRUE) && future::nbrOfWorkers() > 1L
  workers <- if (par) future::nbrOfWorkers() else 1L
  message(sprintf(paste0("simulate_design: repetition 1 took %.1f s; expect about %.1f min ",
                         "for %d repetitions on %d worker(s)."),
                  secs, secs * p$reps / workers / 60, p$reps, workers))
  f <- function(r) .one_rep(r, X, g, surfaces, p)
  rest <- if (par) future.apply::future_lapply(2:p$reps, f, future.seed = TRUE) else
    lapply(2:p$reps, f)
  all_reps <- c(list(first), rest)
  results <- do.call(rbind, lapply(all_reps, `[[`, "results"))
  risks <- do.call(rbind, lapply(all_reps, `[[`, "risks"))
  metrics <- .sim_metrics(results, truths, p)
  scope <- sprintf(paste0("Feasibility is certified only over the declared family of %d ",
                          "outcome surfaces (%s); a design can fail outside it."),
                   length(surfaces), paste(names(surfaces), collapse = "; "))
  if (p$outcome_type == "tte")
    scope <- paste(scope, "The outcome was simulated as risk by the target time;",
                   "the hazard models were not stress-tested.")
  out <- list(results = results, truths = truths, metrics = metrics,
              verdict = .sim_verdict(metrics, p),
              edge = if (is.null(risks)) NULL else .edge_check(risks),
              surfaces = names(surfaces), reps = p$reps, K = p$K, scope = scope,
              design_hash = design$stamp$hash)
  out$stamp <- .stamp(lock, out)
  class(out) <- "cr_simulation"
  out
}
