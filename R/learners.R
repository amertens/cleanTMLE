# A small cross-fitted super learner. Each learner type is a set of
# settings; path learners (xgboost rounds, glmnet penalties) choose their
# best point on the inner folds from one fit per fold, and every setting's
# cross-validated risk is kept for the edge check.

.xgb_cap <- 500L
.xgb_depths <- 1:5
.xgb_etas <- c(0.1, 0.3)
.glmnet_alphas <- c(0, 0.5, 1)
.glmnet_nlambda <- 50L
.nnet_size <- 5L
.nnet_decays <- c(0.001, 0.01, 0.1)

.logloss <- function(y, p) {
  p <- .bound(p, 1e-12)
  -mean(y * log(p) + (1 - y) * log(1 - p))
}

.setting <- function(type, ...) {
  params <- list(...)
  id <- if (length(params))
    paste(type, paste0(names(params), unlist(params), collapse = "_"), sep = "_")
  else type
  list(id = id, type = type, params = params)
}

.settings <- function(learners) {
  out <- list()
  for (l in learners) out <- c(out, switch(l,
    glm = list(.setting("glm")),
    earth = list(.setting("earth", degree = 2L)),
    glmnet = lapply(.glmnet_alphas, function(a) .setting("glmnet", alpha = a)),
    nnet = lapply(.nnet_decays, function(d) .setting("nnet", decay = d)),
    xgboost = unlist(lapply(.xgb_depths, function(d) lapply(.xgb_etas, function(e)
      .setting("xgboost", max_depth = d, eta = e))), recursive = FALSE),
    stop("Unknown learner: ", l, call. = FALSE)))
  out
}

.prepare_setting <- function(s, X, y) {
  if (s$type == "glmnet") {
    if (ncol(X) < 2L) stop("glmnet needs at least two columns.", call. = FALSE)
    s$params$lambda <- glmnet::glmnet(X, y, family = "binomial", alpha = s$params$alpha,
                                      nlambda = .glmnet_nlambda)$lambda
  }
  s
}

.xgb_range <- function(k) {
  if (utils::packageVersion("xgboost") >= "2.0.0") c(1L, k) else c(1L, k + 1L)
}

.learner_fit <- function(s, X, y, X_eval = NULL, y_eval = NULL, path_index = NULL,
                         seed = 1L) {
  p <- s$params
  switch(s$type,
    glm = {
      fit <- suppressWarnings(stats::glm.fit(cbind(1, X), y, family = stats::binomial()))
      b <- fit$coefficients
      b[is.na(b)] <- 0
      list(coef = b)
    },
    earth = earth::earth(x = X, y = y, degree = p$degree,
                         glm = list(family = stats::binomial())),
    glmnet = glmnet::glmnet(X, y, family = "binomial", alpha = p$alpha, lambda = p$lambda),
    nnet = {
      mu <- colMeans(X)
      sdv <- apply(X, 2, stats::sd)
      sdv[!is.finite(sdv) | sdv == 0] <- 1
      Xs <- sweep(sweep(X, 2, mu), 2, sdv, "/")
      fit <- withr::with_seed(seed, nnet::nnet(
        x = Xs, y = y, size = .nnet_size, decay = p$decay, entropy = TRUE,
        maxit = 200L, trace = FALSE, MaxNWts = 100000L))
      list(fit = fit, mu = mu, sd = sdv)
    },
    xgboost = {
      params <- list(objective = "binary:logistic", eval_metric = "logloss",
                     max_depth = p$max_depth, eta = p$eta, min_child_weight = 0,
                     tree_method = "hist", nthread = 1L)
      args <- list(params = params, data = xgboost::xgb.DMatrix(X, label = y),
                   nrounds = path_index %||% .xgb_cap, verbose = 0)
      evals <- if (is.null(X_eval)) list() else
        list(eval = xgboost::xgb.DMatrix(X_eval, label = y_eval))
      evals_arg <- if ("evals" %in% names(formals(xgboost::xgb.train))) "evals" else "watchlist"
      args[[evals_arg]] <- evals
      do.call(xgboost::xgb.train, args)
    },
    stop("Unknown learner type: ", s$type, call. = FALSE))
}

.learner_predict <- function(m, s, X, k) {
  p <- switch(s$type,
    glm = stats::plogis(drop(cbind(1, X) %*% m$coef)),
    earth = as.numeric(stats::predict(m, X, type = "response")),
    glmnet = as.numeric(stats::predict(m, X, s = s$params$lambda[k], type = "response")),
    nnet = as.numeric(stats::predict(m$fit, sweep(sweep(X, 2, m$mu), 2, m$sd, "/"),
                                     type = "raw")),
    xgboost = stats::predict(m, X, iterationrange = .xgb_range(k)))
  .bound(p, 1e-6)
}

.learner_path_loss <- function(m, s, X_eval, y_eval) {
  switch(s$type,
    glmnet = {
      P <- length(s$params$lambda)
      pr <- stats::predict(m, X_eval, type = "response")
      if (ncol(pr) < P) pr <- cbind(pr, pr[, rep(ncol(pr), P - ncol(pr)), drop = FALSE])
      apply(pr, 2, function(p) .logloss(y_eval, p))
    },
    xgboost = {
      log <- m$evaluation_log %||% attr(m, "evaluation_log")
      log[["eval_logloss"]]
    },
    .logloss(y_eval, .learner_predict(m, s, X_eval, 1L)))
}

.cv_setting <- function(s, X, y, inner, newX, seed) {
  V <- max(inner)
  n <- nrow(X)
  s <- .prepare_setting(s, X, y)
  fits <- vector("list", V)
  loss <- 0
  for (v in seq_len(V)) {
    va <- inner == v
    fits[[v]] <- .learner_fit(s, X[!va, , drop = FALSE], y[!va],
                              X[va, , drop = FALSE], y[va], seed = seed + v)
    loss <- loss + sum(va) * .learner_path_loss(fits[[v]], s, X[va, , drop = FALSE], y[va])
  }
  loss <- loss / n
  k <- which.min(loss)
  z <- numeric(n)
  for (v in seq_len(V)) {
    va <- inner == v
    z[va] <- .learner_predict(fits[[v]], s, X[va, , drop = FALSE], k)
  }
  final <- .learner_fit(s, X, y, path_index = k, seed = seed)
  list(z = z, pred = lapply(newX, function(nx) .learner_predict(final, s, nx, k)),
       k = k, P = length(loss), cv_risk = loss[k])
}

.risk_row <- function(s, k = NA_integer_, P = NA_integer_, cv_risk = NA_real_,
                      failed = FALSE, message = NA_character_) {
  p <- s$params
  data.frame(setting = s$id, type = s$type,
             max_depth = p$max_depth %||% NA_real_, eta = p$eta %||% NA_real_,
             alpha = p$alpha %||% NA_real_, decay = p$decay %||% NA_real_,
             k = k, P = P, cv_risk = cv_risk, weight = NA_real_,
             failed = failed, message = message, stringsAsFactors = FALSE)
}

.nnloglik <- function(Z, y) {
  L <- stats::qlogis(.bound(Z, 1e-6))
  J <- ncol(L)
  if (J == 1L) return(1)
  obj <- function(b) .logloss(y, stats::plogis(drop(L %*% b)))
  grad <- function(b) {
    p <- stats::plogis(drop(L %*% b))
    -drop(crossprod(L, y - p)) / length(y)
  }
  w <- stats::optim(rep(1 / J, J), obj, grad, method = "L-BFGS-B", lower = 0)$par
  if (sum(w) <= 0) w <- rep(1 / J, J)
  w / sum(w)
}

.super_learner <- function(X, y, learners, newX, V, seed) {
  if (length(unique(y)) < 2L)
    stop("The outcome is constant in a training set.", call. = FALSE)
  settings <- .settings(learners)
  if (length(settings) == 1L && settings[[1]]$type == "glm") {
    m <- .learner_fit(settings[[1]], X, y)
    rk <- .risk_row(settings[[1]])
    rk$weight <- 1
    return(list(pred = lapply(newX, function(nx) .learner_predict(m, settings[[1]], nx, 1L)),
                weights = c(glm = 1), risks = rk))
  }
  inner <- .make_folds(nrow(X), V, seed)
  res <- lapply(settings, function(s)
    tryCatch(.cv_setting(s, X, y, inner, newX, seed),
             error = function(e) conditionMessage(e)))
  ok <- !vapply(res, is.character, logical(1))
  if (!any(ok))
    stop("Every learner failed: ", paste(unique(unlist(res)), collapse = "; "),
         call. = FALSE)
  risks <- do.call(rbind, Map(function(s, r) {
    if (is.character(r)) .risk_row(s, failed = TRUE, message = r)
    else .risk_row(s, r$k, r$P, r$cv_risk)
  }, settings, res))
  Z <- do.call(cbind, lapply(res[ok], `[[`, "z"))
  w <- .nnloglik(Z, y)
  risks$weight[ok] <- w
  pred <- lapply(seq_along(newX), function(j) {
    L <- stats::qlogis(do.call(cbind, lapply(res[ok], function(r) r$pred[[j]])))
    .bound(stats::plogis(drop(L %*% w)), 1e-6)
  })
  names(w) <- vapply(settings[ok], `[[`, character(1), "id")
  list(pred = pred, weights = w, risks = risks)
}

.edge_check <- function(risks) {
  empty <- data.frame(learner = character(), edge = character(), n_fits = integer(),
                      n_at_edge = integer(), share = numeric(), consistent = logical())
  r <- risks[!risks$failed & !is.na(risks$cv_risk), , drop = FALSE]
  if (!nrow(r)) return(empty)
  ctx_cols <- intersect(c("rep", "surface", "library", "nuisance", "fold"), names(r))
  ctx <- if (length(ctx_cols)) interaction(r[ctx_cols], drop = TRUE) else
    factor(rep("all", nrow(r)))
  hits <- list()
  for (cx in levels(ctx)) for (ty in c("xgboost", "glmnet", "nnet")) {
    d <- r[ctx == cx & r$type == ty, , drop = FALSE]
    if (!nrow(d)) next
    b <- d[which.min(d$cv_risk), ]
    e <- switch(ty,
      xgboost = c(
        "max_depth at the top of its grid" = b$max_depth == max(.xgb_depths),
        "number of trees at the cap" = b$k == b$P,
        "fewer than 20 trees at the lowest learning rate" =
          b$k < 20 && b$eta == min(.xgb_etas)),
      glmnet = c("penalty at the smallest value on the path" = b$k == b$P),
      nnet = c("decay at the bottom of its grid" = b$decay == min(.nnet_decays),
               "decay at the top of its grid" = b$decay == max(.nnet_decays)))
    hits[[length(hits) + 1L]] <- data.frame(learner = ty, edge = names(e),
                                            hit = unname(e), stringsAsFactors = FALSE)
  }
  if (!length(hits)) return(empty)
  h <- do.call(rbind, hits)
  keys <- unique(h[c("learner", "edge")])
  out <- do.call(rbind, lapply(seq_len(nrow(keys)), function(i) {
    x <- h$hit[h$learner == keys$learner[i] & h$edge == keys$edge[i]]
    data.frame(learner = keys$learner[i], edge = keys$edge[i], n_fits = length(x),
               n_at_edge = sum(x), share = mean(x), consistent = mean(x) >= 0.5,
               stringsAsFactors = FALSE)
  }))
  out <- out[out$n_at_edge > 0, , drop = FALSE]
  rownames(out) <- NULL
  out
}
