# One fitter for the simulation and the final analysis. Nuisances are
# cross-fitted once; targeting is cheap, so one nuisance fit serves every
# estimand and truncation level.

.design_matrix <- function(data, covariates) {
  f <- stats::reformulate(sprintf("`%s`", covariates))
  X <- stats::model.matrix(f, data = data[, covariates, drop = FALSE])[, -1, drop = FALSE]
  storage.mode(X) <- "double"
  dimnames(X) <- list(NULL, colnames(X))
  X
}

.cross_fit <- function(X, y, learners, folds, V, seed, newX_fun, groups = NULL) {
  K <- max(folds)
  n <- nrow(X)
  out <- NULL
  risks <- list()
  for (k in seq_len(K)) {
    va <- if (K == 1L) rep(TRUE, n) else folds == k
    tr <- if (K == 1L) rep(TRUE, n) else !va
    if (!any(va)) next
    sl <- .super_learner(X[tr, , drop = FALSE], y[tr], learners, newX_fun(which(va)),
                         V, seed + 1000L * k, groups[tr])
    if (is.null(out)) out <- lapply(sl$pred, function(p) rep(NA_real_, n))
    for (j in seq_along(out)) out[[j]][va] <- sl$pred[[j]]
    risks[[length(risks) + 1L]] <- cbind(fold = k, sl$risks)
  }
  list(pred = out, risks = do.call(rbind, risks))
}

.fit_g <- function(X, A, ps_library, folds, V, seed, groups = NULL) {
  cf <- .cross_fit(X, A, ps_library, folds, V, seed,
                   function(rows) list(X[rows, , drop = FALSE]), groups)
  list(g = cf$pred[[1]], risks = cbind(nuisance = "g", cf$risks),
       folds = folds, seed = seed)
}

.fit_Q <- function(X, A, Y, q_library, folds, V, seed, groups = NULL) {
  XA <- cbind(.trt = A, X)
  cf <- .cross_fit(XA, Y, q_library, folds, V, seed, function(rows) {
    x1 <- XA[rows, , drop = FALSE]; x1[, ".trt"] <- 1
    x0 <- XA[rows, , drop = FALSE]; x0[, ".trt"] <- 0
    list(x1, x0)
  }, groups)
  list(Q1 = cf$pred[[1]], Q0 = cf$pred[[2]], risks = cbind(nuisance = "Q", cf$risks))
}

.target_tmle <- function(Y, A, X, g, Q1, Q0, truncation, which) {
  if (length(unique(A)) < 2L)
    stop("Only one treatment arm in the target population.", call. = FALSE)
  # tmle refits a treatment model for the ATT when the overlap subset drops more
  # than a tenth of the controls, whatever g1W it was given. Its default library
  # (glm, dbarts, gam with 10 folds) takes minutes and is not cross-fitted, so
  # the refit is limited to a logistic regression.
  f <- suppressMessages(suppressWarnings(tmle::tmle(
    Y = Y, A = A, W = as.data.frame(X), Q = cbind(Q0, Q1), g1W = g,
    gbound = c(truncation, 1 - truncation), family = "binomial",
    g.SL.library = "SL.glm")))
  e <- f$estimates[[which]]
  r1 <- if (which == "ATE") f$estimates$EY1$psi else mean(Y[A == 1])
  r0 <- if (which == "ATE") f$estimates$EY0$psi else r1 - e$psi
  list(estimate = e$psi, se = sqrt(e$var.psi), ci_lower = e$CI[1], ci_upper = e$CI[2],
       n_population = if (which == "ATT") sum(A == 1) else length(Y),
       risk1 = r1, risk0 = r0)
}

.ato_eif <- function(Y, A, g, Q1, Q0, psi) {
  h <- g * (1 - g)
  Dh <- mean(h)
  QA <- ifelse(A == 1, Q1, Q0)
  dy <- (A * (1 - g) - (1 - A) * g) * (Y - QA) / Dh
  dg <- (Q1 - Q0 - psi) * (1 - 2 * g) * (A - g) / Dh
  dw <- h * (Q1 - Q0 - psi) / Dh
  list(y = dy, g = dg, total = dy + dg + dw)
}

.target_ato <- function(Y, A, g, Q1, Q0, max_iter = 20L) {
  if (length(unique(A)) < 2L)
    stop("Only one treatment arm in the target population.", call. = FALSE)
  g <- .bound(g, 1e-6); Q1 <- .bound(Q1, 1e-6); Q0 <- .bound(Q0, 1e-6)
  n <- length(Y)
  psi_of <- function() sum(g * (1 - g) * (Q1 - Q0)) / sum(g * (1 - g))
  for (it in seq_len(max_iter)) {
    QA <- ifelse(A == 1, Q1, Q0)
    HY <- A * (1 - g) - (1 - A) * g
    eps <- stats::coef(suppressWarnings(stats::glm(
      Y ~ -1 + HY + offset(stats::qlogis(QA)), family = stats::binomial())))
    eps[is.na(eps)] <- 0
    Q1 <- .bound(stats::plogis(stats::qlogis(Q1) + eps * (1 - g)), 1e-6)
    Q0 <- .bound(stats::plogis(stats::qlogis(Q0) - eps * g), 1e-6)
    psi <- psi_of()
    Hg <- (Q1 - Q0 - psi) * (1 - 2 * g)
    del <- stats::coef(suppressWarnings(stats::glm(
      A ~ -1 + Hg + offset(stats::qlogis(g)), family = stats::binomial())))
    del[is.na(del)] <- 0
    g <- .bound(stats::plogis(stats::qlogis(g) + del * Hg), 1e-6)
    psi <- psi_of()
    D <- .ato_eif(Y, A, g, Q1, Q0, psi)
    if (max(abs(c(mean(D$y), mean(D$g)))) < stats::sd(D$total) / (sqrt(n) * log(n))) break
  }
  se <- stats::sd(D$total) / sqrt(n)
  h <- g * (1 - g)
  list(estimate = psi, se = se, ci_lower = psi - 1.96 * se, ci_upper = psi + 1.96 * se,
       n_population = n, risk1 = sum(h * Q1) / sum(h), risk0 = sum(h * Q0) / sum(h))
}

.target <- function(Y, A, X, g, Q1, Q0, estimand, truncation, trim_band) {
  switch(estimand,
    ATE = .target_tmle(Y, A, X, g, Q1, Q0, truncation, "ATE"),
    ATT = .target_tmle(Y, A, X, g, Q1, Q0, truncation, "ATT"),
    trimmed_ATE = {
      keep <- g >= trim_band[1] & g <= trim_band[2]
      r <- .target_tmle(Y[keep], A[keep], X[keep, , drop = FALSE], g[keep],
                        Q1[keep], Q0[keep], truncation, "ATE")
      r$n_population <- sum(keep)
      r
    },
    ATO = .target_ato(Y, A, g, Q1, Q0),
    stop("Unknown estimand: ", estimand, call. = FALSE))
}

.implausibility_check <- function(est, y, a) {
  ok <- !is.na(y) & !is.na(a)
  y <- y[ok]; a <- a[ok]
  if (!is.finite(est) || length(unique(a)) < 2L)
    return(list(implausible = FALSE, reason = NA_character_))
  m1 <- mean(y[a == 1]); m0 <- mean(y[a == 0])
  cd <- m1 - m0
  noise <- 0.005
  f <- character()
  if (abs(est) > noise && abs(cd) > noise && sign(est) != sign(cd))
    f <- c(f, sprintf("sign differs from crude (%.4g)", cd))
  if (abs(cd) > noise && abs(est) > 5 * abs(cd))
    f <- c(f, sprintf("%.1fx crude difference", abs(est / cd)))
  if (abs(est) > max(m1, m0))
    f <- c(f, sprintf("risk difference exceeds largest arm risk (%.4g)", max(m1, m0)))
  list(implausible = length(f) > 0L,
       reason = if (length(f)) paste(f, collapse = "; ") else NA_character_)
}

.estimate_row <- function(estimand, truncation, r = NULL, guard = NULL, message = NA_character_) {
  if (is.null(r))
    return(data.frame(estimand = estimand, truncation = truncation, estimate = NA_real_,
                      se = NA_real_, ci_lower = NA_real_, ci_upper = NA_real_,
                      n_population = NA_integer_, risk1 = NA_real_, risk0 = NA_real_,
                      implausible = NA, implausible_reason = NA_character_,
                      failed = TRUE, message = message, stringsAsFactors = FALSE))
  data.frame(estimand = estimand, truncation = truncation, estimate = r$estimate,
             se = r$se, ci_lower = r$ci_lower, ci_upper = r$ci_upper,
             n_population = r$n_population, risk1 = r$risk1, risk0 = r$risk0,
             implausible = guard$implausible, implausible_reason = guard$reason,
             failed = FALSE, message = NA_character_, stringsAsFactors = FALSE)
}

.fit_and_target <- function(X, A, Y, q_library, estimands, truncations, plan, folds,
                            seed, g_fit = NULL, groups = NULL) {
  if (is.null(g_fit)) g_fit <- .fit_g(X, A, plan$ps_library, folds, plan$V, seed, groups)
  Qf <- .fit_Q(X, A, Y, q_library, folds, plan$V, seed + 1L, groups)
  one <- function(e, t) {
    r <- tryCatch(.target(Y, A, X, g_fit$g, Qf$Q1, Qf$Q0, e, t, plan$trim_band),
                  error = function(err) conditionMessage(err))
    if (is.character(r)) return(.estimate_row(e, t, message = r))
    vals <- c(r$estimate, r$se, r$ci_lower, r$ci_upper)
    if (length(vals) != 4L || !all(is.finite(vals)))
      return(.estimate_row(e, t, message = "non-finite estimate or interval"))
    .estimate_row(e, t, r, .implausibility_check(r$estimate, Y, A))
  }
  rows <- list()
  for (e in estimands) {
    if (e == "ATO") {
      r <- one(e, truncations[1])
      for (t in truncations) { r$truncation <- t; rows[[length(rows) + 1L]] <- r }
    } else {
      for (t in truncations) rows[[length(rows) + 1L]] <- one(e, t)
    }
  }
  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  attr(out, "risks") <- rbind(g_fit$risks, Qf$risks)
  attr(out, "g") <- g_fit$g
  out
}

fit_candidate <- function(X, A, Y, estimand, candidate, plan, folds, seed, g_fit = NULL) {
  .fit_and_target(X, A, Y, candidate$learners, estimand, candidate$truncation, plan,
                  folds, seed, g_fit)
}
