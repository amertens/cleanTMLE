#' Utility Functions
#'
#' Small helper functions used throughout the cleanTMLE workflow.
#'
#' @name utils
NULL


#' Expit (Inverse-Logit) Function
#'
#' Maps a real-valued linear predictor to the probability scale via the
#' inverse-logit transformation: \eqn{p = 1 / (1 + e^{-x})}.
#'
#' @param x Numeric vector of linear-predictor values.
#'
#' @return Numeric vector of probabilities in (0, 1).
#'
#' @examples
#' expit(0)    # 0.5
#' expit(2)    # ~0.88
#' expit(-Inf) # 0
#'
#' @seealso [logit()] for the inverse operation.
#' @export
expit <- function(x) {
  1 / (1 + exp(-x))
}


#' Logit Function
#'
#' Maps a probability to the log-odds (logit) scale:
#' \eqn{\text{logit}(p) = \log(p / (1 - p))}.
#'
#' @param p Numeric vector of probabilities in (0, 1).
#'
#' @return Numeric vector of log-odds values.
#'
#' @examples
#' logit(0.5)  # 0
#' logit(0.9)  # ~2.20
#'
#' @seealso [expit()] for the inverse operation.
#' @export
logit <- function(p) {
  log(p / (1 - p))
}
