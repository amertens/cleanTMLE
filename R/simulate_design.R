.std <- function(v) {
  s <- stats::sd(v)
  if (!is.finite(s) || s == 0) v * 0 else (v - mean(v)) / s
}

.surfaces <- function(plan, X, g, seed) {
  s <- plan$surfaces
  lg <- stats::qlogis(.bound(g, 1e-6))
  Z <- apply(X, 2, .std)
  if (is.null(dim(Z))) Z <- matrix(Z, nrow = nrow(X))
  dirn <- .with_seed(seed, stats::rnorm(ncol(Z)))
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

# One repetition. `surfaces_run` and `libraries_run` (crossed), or an explicit
# data frame of `pairs` (surface, library), restrict which outcome fits run.
# The outcomes for every surface are still drawn inside the seeded block, and
# each pair's fit runs under the repetition's seed, so repetition r gives the
# same rows for a pair whichever other pairs run with it.
.one_rep <- function(r, X, g, surfaces, plan, surfaces_run = names(surfaces),
                     libraries_run = names(plan$candidates$library), pairs = NULL) {
  if (is.null(pairs))
    pairs <- expand.grid(library = libraries_run, surface = surfaces_run,
                         stringsAsFactors = FALSE)[c("surface", "library")]
  n <- nrow(X)
  seed <- plan$seed + r
  draw <- .with_seed(seed, {
    idx <- sample.int(n, n, replace = TRUE)
    a <- stats::rbinom(n, 1L, g[idx])
    ys <- lapply(surfaces, function(s)
      stats::rbinom(n, 1L, ifelse(a == 1L, s$q1[idx], s$q0[idx])))
    list(idx = idx, a = a, ys = ys)
  })
  Xs <- X[draw$idx, , drop = FALSE]
  # Folds are assigned per original row so that copies of one row never sit on
  # both sides of a split; the inner folds are grouped the same way.
  folds <- .make_folds(n, plan$K, seed)[draw$idx]
  g_fit <- tryCatch(.fit_g(Xs, draw$a, plan$ps_library, folds, plan$V, seed,
                           groups = draw$idx),
                    error = function(e) conditionMessage(e))
  rows <- list()
  risks <- list()
  for (i in seq_len(nrow(pairs))) {
    sn <- pairs$surface[i]
    lib <- pairs$library[i]
    res <- if (is.character(g_fit)) .failed_rows(plan, g_fit) else tryCatch(
      .with_seed(seed, .fit_and_target(
        Xs, draw$a, draw$ys[[sn]], plan$candidates$library[[lib]],
        plan$estimands, plan$candidates$truncation, plan, folds, seed, g_fit,
        groups = draw$idx)),
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

# Wilson score interval for x successes in n trials.
.wilson <- function(x, n, level = 0.95) {
  if (!is.finite(n) || n < 1) return(c(NA_real_, NA_real_))
  z <- stats::qnorm(1 - (1 - level) / 2)
  p <- x / n
  mid <- (p + z^2 / (2 * n)) / (1 + z^2 / n)
  half <- z / (1 + z^2 / n) * sqrt(p * (1 - p) / n + z^2 / (4 * n^2))
  c(max(0, mid - half), min(1, mid + half))
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
    ok <- s$failed %in% FALSE & is.finite(s$estimate) & is.finite(s$ci_lower) &
      is.finite(s$ci_upper)
    R <- sum(ok)
    e <- s$estimate[ok] - s$truth[ok]
    cv <- s$ci_lower[ok] <= s$truth[ok] & s$truth[ok] <= s$ci_upper[ok]
    bias <- if (R) mean(e) else NA_real_
    coverage <- if (R) mean(cv) else NA_real_
    rmse <- if (R) sqrt(mean(e^2)) else NA_real_
    wi <- .wilson(sum(cv), R)
    data.frame(k, reps = nrow(s), failed = nrow(s) - R, fail_rate = 1 - R / nrow(s),
               bias = bias, bias_mcse = if (R > 1) stats::sd(e) / sqrt(R) else NA_real_,
               coverage = coverage,
               coverage_mcse = if (R > 1) sqrt(coverage * (1 - coverage) / R) else NA_real_,
               coverage_lower = wi[1], coverage_upper = wi[2],
               rmse = rmse,
               rmse_mcse = if (R > 1 && rmse > 0) stats::sd(e^2) / (2 * rmse * sqrt(R)) else NA_real_,
               stringsAsFactors = FALSE)
  }))
  tb <- plan$tolerance$bias
  tc <- plan$tolerance$coverage
  # Feasibility is judged on the point estimates: the tolerances are the
  # declared margin, and no second margin is applied.
  usable <- !is.na(m$bias) & m$reps - m$failed > 1 & m$fail_rate <= 0.05
  m$pass <- usable & abs(m$bias) <= tb & m$coverage >= tc
  # A usable cell is unresolved when the 95% interval for its bias, or the
  # Wilson interval for its coverage, contains the tolerance. A cell that fails
  # the failure rule is not unresolved; it fails.
  hw <- stats::qnorm(0.975) * m$bias_mcse
  bias_open <- !is.na(hw) & abs(m$bias) - hw <= tb & abs(m$bias) + hw > tb
  cov_open <- !is.na(m$coverage_lower) & !is.na(m$coverage_upper) &
    m$coverage_lower < tc & m$coverage_upper >= tc
  m$unresolved <- usable & (bias_open | cov_open)
  rownames(m) <- NULL
  m
}

# One row per candidate of one estimand's metrics.
.candidate_summary <- function(me) {
  do.call(rbind, lapply(split(me, me$candidate), function(d) data.frame(
    candidate = d$candidate[1], library = d$library[1], truncation = d$truncation[1],
    pass_all = all(d$pass),
    # Fails somewhere, but only on unresolved cells: it could still pass everywhere.
    could_flip = !all(d$pass) && all(d$pass | d$unresolved),
    worst_rmse = max(d$rmse), stringsAsFactors = FALSE)))
}

# The candidates whose unresolved cells can change an estimand's decision
# (its status or its selected candidate): the selected candidate, or every
# candidate that passes everywhere when the plan supplies its own `select`
# rule, and every candidate whose failures are all unresolved. A candidate
# with a resolved failure cannot pass, so its cells cannot change anything.
.decision_candidates <- function(cs, selected, plan) {
  chosen <- if (is.na(selected)) character() else
    if (is.null(plan$select)) selected else cs$candidate[cs$pass_all]
  union(chosen, cs$candidate[cs$could_flip])
}

# Logical over the rows of `metrics`: the unresolved cells that can change
# their estimand's decision. verdict$unresolved, the extension and the
# decision note all use this one set.
.decision_cells <- function(metrics, verdict, plan) {
  out <- logical(nrow(metrics))
  for (i in seq_len(nrow(verdict))) {
    k <- metrics$estimand == verdict$estimand[i]
    cs <- .candidate_summary(metrics[k, , drop = FALSE])
    cand <- .decision_candidates(cs, verdict$selected[i], plan)
    out <- out | (k & metrics$unresolved & metrics$candidate %in% cand)
  }
  out
}

.sim_verdict <- function(metrics, plan) {
  do.call(rbind, lapply(plan$estimands, function(e) {
    me <- metrics[metrics$estimand == e, , drop = FALSE]
    cs <- .candidate_summary(me)
    pool <- cs[cs$pass_all, , drop = FALSE]
    status <- if (nrow(pool)) "feasible" else "infeasible"
    sel <- NA_character_
    if (nrow(pool)) {
      sel <- if (is.null(plan$select)) pool$candidate[which.min(pool$worst_rmse)] else
        plan$select(me[me$candidate %in% pool$candidate, , drop = FALSE])
      if (!sel %in% pool$candidate)
        stop("`select` must return one of: ", paste(pool$candidate, collapse = ", "),
             call. = FALSE)
    }
    row <- pool[match(sel, pool$candidate), , drop = FALSE]
    cand <- .decision_candidates(cs, sel, plan)
    data.frame(estimand = e, status = status, selected = sel,
               library = if (is.na(sel)) NA_character_ else row$library,
               truncation = if (is.na(sel)) NA_real_ else row$truncation,
               unresolved = any(me$unresolved & me$candidate %in% cand),
               stringsAsFactors = FALSE)
  }))
}

#' Outcome-free plasmode simulation: the support map and candidate selection
#'
#' Resamples the design rows, draws treatment from the design's single-fit
#' propensity score and outcomes from the plan's family of surfaces, and runs every
#' candidate with the same fitter that `estimate_effect()` uses. The
#' simulated world uses a single-fit propensity score (the plan's
#' `ps_library` fitted once on all design rows), which depends on the
#' covariates only; the cross-fitted score that grades overlap would make
#' each row's treatment depend on its fold, which is not a covariate. An
#' estimand is feasible when some candidate meets the bias and coverage
#' tolerances on every surface, judged on the point estimates. A cell whose
#' 95% interval for bias, or Wilson 95% interval for coverage, contains its
#' tolerance is unresolved; the unresolved cells that can change the decision
#' of an estimand at or above the primary (those of its selected candidate,
#' and of any candidate whose failures are all unresolved) get further batches
#' of `reps` repetitions, up to the plan's `max_reps`. Because the intervals
#' are recomputed after each batch, "resolved" is a stopping label, not a
#' formal 95% statement; the decisions themselves use the point estimates.
#' Repetitions run in parallel through
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
  g <- design$g_world  # the world's propensity score: covariates only, never the fold
  if (is.null(g) || length(g) != nrow(X))
    stop("`design` has no single-fit propensity score; rerun assess_design().", call. = FALSE)
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
  # nbrOfWorkers() can be Inf (for example under future.callr), so it is
  # capped at the number of repetitions and printed with %s.
  workers <- if (par) future::nbrOfWorkers() else 1L
  message(sprintf(paste0("simulate_design: repetition 1 took %.1f s; expect about %.1f min ",
                         "for %d repetitions on %s worker(s)."),
                  secs, secs * p$reps / min(workers, p$reps) / 60, p$reps, format(workers)))
  run_batch <- function(rs, pairs = NULL) {
    f <- if (is.null(pairs)) function(r) .one_rep(r, X, g, surfaces, p) else
      function(r) .one_rep(r, X, g, surfaces, p, pairs = pairs)
    if (par) future.apply::future_lapply(rs, f, future.seed = TRUE) else lapply(rs, f)
  }
  all_reps <- c(list(first), run_batch(2:p$reps))
  results <- do.call(rbind, lapply(all_reps, `[[`, "results"))
  risks <- do.call(rbind, lapply(all_reps, `[[`, "risks"))

  # Targeted extension: the unresolved cells that can change the decision
  # (see .decision_cells()) of the estimands at or above the primary get
  # further batches of repetitions, numbered after the last one run, for their
  # (surface, library) pairs only, until they resolve or reach max_reps. The
  # set is recomputed after every batch. Cells of estimands below the primary
  # keep their repetitions.
  max_reps <- p$max_reps %||% p$reps
  last_r <- p$reps
  ext <- list()
  cell_key <- function(d) paste(d$estimand, d$surface, d$library, sep = "\r")
  run_key <- function(d) paste(d$rep, d$surface, d$library, sep = "\r")
  repeat {
    metrics <- .sim_metrics(results, truths, p)
    verdict <- .sim_verdict(metrics, p)
    pos <- match("feasible", verdict$status)
    if (is.na(pos)) pos <- length(p$estimands)
    relevant <- p$estimands[seq_len(pos)]
    dec <- .decision_cells(metrics, verdict, p) & metrics$estimand %in% relevant
    unres <- metrics[dec, , drop = FALSE]
    open <- unres[unres$reps < max_reps, , drop = FALSE]
    if (!nrow(open)) break
    pairs <- unique(open[c("surface", "library")])
    rownames(pairs) <- NULL
    m_add <- min(p$reps, max_reps - min(open$reps))
    rs <- last_r + seq_len(m_add)
    message(sprintf(paste0("simulate_design: extending %d surface-library pair(s) by %d ",
                           "repetitions (unresolved cells: %d)."),
                    nrow(pairs), m_add, nrow(unres)))
    new <- run_batch(rs, pairs)
    nr <- do.call(rbind, lapply(new, `[[`, "results"))
    nr <- nr[nr$estimand %in% relevant, , drop = FALSE]
    # No cell goes past max_reps: keep, per cell, only the repetitions it has room for.
    have <- tapply(metrics$reps, cell_key(metrics), max)
    room <- max_reps - have[cell_key(nr)]
    room[is.na(room)] <- max_reps
    nr <- nr[nr$rep - last_r <= room, , drop = FALSE]
    results <- rbind(results, nr)
    # Keep the learner risks only for the repetitions and pairs whose rows were kept.
    nk <- do.call(rbind, lapply(new, `[[`, "risks"))
    if (!is.null(nk)) risks <- rbind(risks, nk[run_key(nk) %in% run_key(nr), , drop = FALSE])
    ext[[length(ext) + 1L]] <- data.frame(
      batch = length(ext) + 1L, reps_from = min(rs), reps_to = max(rs),
      pairs = paste(pairs$surface, pairs$library, sep = " | ", collapse = "; "),
      cells_unresolved_before = nrow(unres), stringsAsFactors = FALSE)
    last_r <- max(rs)
  }
  rownames(results) <- NULL
  extension <- if (length(ext)) do.call(rbind, ext) else
    data.frame(batch = integer(), reps_from = integer(), reps_to = integer(),
               pairs = character(), cells_unresolved_before = integer(),
               stringsAsFactors = FALSE)
  lead <- "Feasibility is judged on point estimates against the declared tolerances; "
  decision_note <- if (nrow(extension)) {
    sprintf(paste0(lead, "cells whose 95%% interval contained a tolerance and that could ",
                   "change the decision received extra repetitions up to max_reps = %d ",
                   "(%d such cells remain unresolved)."), max_reps, nrow(unres))
  } else if (max_reps <= p$reps && nrow(unres)) {
    sprintf(paste0(lead, "max_reps equals reps, so no extra repetitions were run ",
                   "(%d cells that could change the decision are unresolved)."), nrow(unres))
  } else {
    paste0(lead, "no cell that could change the decision had a 95% interval containing ",
           "a tolerance, so no extra repetitions were run (0 cells unresolved).")
  }
  scope <- sprintf(paste0("Feasibility is certified only over the declared family of %d ",
                          "outcome surfaces (%s); a design can fail outside it."),
                   length(surfaces), paste(names(surfaces), collapse = "; "))
  if (p$outcome_type == "tte")
    scope <- paste(scope, "The outcome was simulated as risk by the target time;",
                   "the hazard models were not stress-tested.")
  out <- list(results = results, truths = truths, metrics = metrics,
              verdict = verdict,
              edge = if (is.null(risks)) NULL else .edge_check(risks),
              surfaces = names(surfaces), reps = p$reps, max_reps = max_reps,
              extension = extension, K = p$K, scope = scope,
              decision_note = decision_note, design_hash = design$stamp$hash)
  out$stamp <- .stamp(lock, out)
  class(out) <- "cr_simulation"
  out
}
