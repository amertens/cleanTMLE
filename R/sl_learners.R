# SuperLearner wrappers the recommended libraries use: penalised-regression
# variants (ridge and elastic net alongside SL.glmnet's lasso) and a
# principal-components GLM for wide, sparse design matrices.

#' Ridge and Elastic-Net glmnet Learners
#'
#' `SL.glmnet` is lasso only (`alpha = 1`). On design matrices with blocks of
#' correlated indicator columns, ridge spreads weight across collinear terms
#' where lasso picks one arbitrarily, so both penalties belong in the
#' ensemble as separate learners.
#'
#' @param ... Passed to [SuperLearner::SL.glmnet].
#' @param alpha The elastic-net mixing parameter (0 ridge, 0.5 elastic net).
#' @return A SuperLearner fit object.
#' @name sl_glmnet_variants
NULL

#' @rdname sl_glmnet_variants
#' @export
SL.glmnet.ridge <- function(..., alpha = 0) {
  SuperLearner::SL.glmnet(..., alpha = alpha)
}

#' @rdname sl_glmnet_variants
#' @export
SL.glmnet.enet <- function(..., alpha = 0.5) {
  SuperLearner::SL.glmnet(..., alpha = alpha)
}


#' Principal-Components GLM Learner
#'
#' A GLM on the leading principal components of the design matrix, as an
#' ensemble member rather than a preprocessing step, so the dimension
#' reduction competes inside cross-validation. On a wide design matrix with
#' many sparse indicators, a full-dimensional propensity model can separate
#' the arms by chaining rare covariates, producing extreme propensities that
#' reflect the design matrix rather than the patients; the PCA learner offers
#' the ensemble a low-dimensional alternative, and the ensemble weights
#' decide. The rotation is computed inside each training fold and applied
#' unchanged to the validation fold.
#'
#' @param Y,X,newX,family,obsWeights Standard SuperLearner learner arguments.
#' @param k Number of principal components. Default 10.
#' @param ... Ignored.
#' @return A SuperLearner learner fit (`pred` and `fit`).
#' @export
SL.glm.pca <- function(Y, X, newX, family, obsWeights, k = 10L, ...) {
  X <- as.data.frame(X); newX <- as.data.frame(newX)
  num <- vapply(X, is.numeric, logical(1))
  Xn <- as.matrix(X[, num, drop = FALSE])
  sds <- apply(Xn, 2, stats::sd)
  keep <- is.finite(sds) & sds > 1e-10
  if (sum(keep) < 2L) {
    # Degenerate design: fall back to an intercept-only GLM.
    fit <- stats::glm(Y ~ 1, family = family, weights = obsWeights)
    pred <- rep(as.numeric(stats::predict(fit, type = "response"))[1],
                nrow(newX))
    out <- list(pred = pred, fit = list(object = fit, pca = NULL,
                                        cols = character(0)))
    class(out$fit) <- "SL.glm.pca"
    return(out)
  }
  Xk <- Xn[, keep, drop = FALSE]
  k_use <- min(as.integer(k), ncol(Xk), nrow(Xk) - 1L)
  pca <- stats::prcomp(Xk, center = TRUE, scale. = TRUE)
  Z <- as.data.frame(pca$x[, seq_len(k_use), drop = FALSE])
  fit <- stats::glm(Y ~ ., data = cbind(data.frame(Y = Y), Z),
                    family = family, weights = obsWeights)
  newXk <- as.matrix(newX[, colnames(Xk), drop = FALSE])
  Znew <- scale(newXk, center = pca$center, scale = pca$scale) %*%
    pca$rotation[, seq_len(k_use), drop = FALSE]
  pred <- as.numeric(stats::predict(
    fit, newdata = as.data.frame(Znew), type = "response"))
  out <- list(pred = pred,
              fit = list(object = fit, pca = pca, k = k_use,
                         cols = colnames(Xk)))
  class(out$fit) <- "SL.glm.pca"
  out
}

#' @export
predict.SL.glm.pca <- function(object, newdata, ...) {
  if (is.null(object$pca)) {
    return(rep(as.numeric(stats::predict(object$object,
                                         type = "response"))[1],
               nrow(as.data.frame(newdata))))
  }
  newX <- as.data.frame(newdata)
  newXk <- as.matrix(newX[, object$cols, drop = FALSE])
  Z <- scale(newXk, center = object$pca$center, scale = object$pca$scale) %*%
    object$pca$rotation[, seq_len(object$k), drop = FALSE]
  as.numeric(stats::predict(object$object, newdata = as.data.frame(Z),
                            type = "response"))
}

#' @rdname SL.glm.pca
#' @export
SL.glm.pca5 <- function(...) SL.glm.pca(..., k = 5L)

#' @rdname SL.glm.pca
#' @export
SL.glm.pca10 <- function(...) SL.glm.pca(..., k = 10L)
