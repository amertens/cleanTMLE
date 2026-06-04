# ============================================================================
# Bounded glmnet SuperLearner wrapper
# ----------------------------------------------------------------------------
# SL.glmnet can spend an unbounded amount of computation on near-separable
# (near-positivity) designs, which are exactly the designs the data-quality
# stress test generates. This wrapper caps the lambda path, the number of
# cross-validation folds, and the coordinate-descent iterations, so the fit
# returns in bounded time on any design. It is a drop-in learner for use in a
# candidate's g- or q-library.
# ============================================================================

#' Bounded glmnet SuperLearner Learner
#'
#' A SuperLearner-compatible wrapper around \code{glmnet::cv.glmnet} with the
#' lambda path, cross-validation folds, and coordinate-descent iterations
#' bounded so the fit cannot run away on near-separable (near-positivity)
#' designs. Use it in place of \code{"SL.glmnet"} in a candidate's
#' \code{g_library} or \code{q_library} when the stress test generates
#' near-positivity designs.
#'
#' @param Y Outcome vector.
#' @param X Predictor data.frame.
#' @param newX Predictor data.frame for prediction.
#' @param family A \code{family} object (binomial or gaussian).
#' @param obsWeights Observation weights.
#' @param id Optional cluster ids (unused; present for SuperLearner API).
#' @param alpha Elastic-net mixing parameter. Default 1 (lasso).
#' @param nlambda Length of the lambda path. Default 20 (SuperLearner uses 100).
#' @param nfolds Cross-validation folds. Default 3 (SuperLearner uses 10).
#' @param maxit Maximum coordinate-descent iterations. Default 1000.
#' @param ... Passed to \code{glmnet::cv.glmnet}.
#'
#' @return A list with \code{pred} (predictions for \code{newX}) and \code{fit}
#'   (an object of class \code{SL.glmnet.bounded}).
#' @export
SL.glmnet.bounded <- function(Y, X, newX, family,
                              obsWeights = rep(1, length(Y)),
                              id = NULL, alpha = 1, nlambda = 20L,
                              nfolds = 3L, maxit = 1000L, ...) {
  if (!requireNamespace("glmnet", quietly = TRUE))
    stop("SL.glmnet.bounded requires the 'glmnet' package.", call. = FALSE)
  Xm   <- as.matrix(X)
  newm <- as.matrix(newX)
  cvfit <- glmnet::cv.glmnet(
    x = Xm, y = Y, weights = obsWeights,
    family = family$family, alpha = alpha,
    nlambda = nlambda, nfolds = nfolds, maxit = maxit, ...)
  pred <- as.numeric(stats::predict(cvfit, newx = newm,
                                    s = "lambda.min", type = "response"))
  fit  <- list(object = cvfit)
  class(fit) <- "SL.glmnet.bounded"
  list(pred = pred, fit = fit)
}

#' SuperLearner lookup environment that includes cleanTMLE learners
#'
#' SuperLearner resolves learner names with \code{get(name, envir = env)}. The
#' candidate-fit sites pass \code{env = asNamespace("SuperLearner")}, which
#' cannot see learners defined in cleanTMLE. This returns a child of the
#' SuperLearner namespace that also exposes \code{SL.glmnet.bounded}, so a
#' candidate may name it in its g- or q-library while every built-in
#' SuperLearner learner still resolves through the parent.
#' @keywords internal
#' @noRd
.cleantmle_sl_env <- function() {
  e <- new.env(parent = asNamespace("SuperLearner"))
  assign("SL.glmnet.bounded", SL.glmnet.bounded, envir = e)
  e
}

#' Predict Method for SL.glmnet.bounded
#'
#' @param object A \code{SL.glmnet.bounded} fit.
#' @param newdata Predictor data.frame.
#' @param ... Unused.
#' @return A numeric vector of predicted responses.
#' @export
predict.SL.glmnet.bounded <- function(object, newdata, ...) {
  if (!requireNamespace("glmnet", quietly = TRUE))
    stop("SL.glmnet.bounded requires the 'glmnet' package.", call. = FALSE)
  as.numeric(stats::predict(object$object, newx = as.matrix(newdata),
                            s = "lambda.min", type = "response"))
}
