`%||%` <- function(x, y) if (is.null(x)) y else x

.hashable <- function(x) {
  if (is.function(x)) return(paste(deparse(x), collapse = "\n"))
  if (inherits(x, "formula")) return(paste(deparse(x), collapse = " "))
  if (is.list(x)) {
    y <- lapply(x, .hashable)
    attributes(y) <- NULL
    names(y) <- names(x)
    return(y)
  }
  x
}

.hash <- function(x) digest::digest(.hashable(x), algo = "sha256")

.bound <- function(p, eps) pmin(pmax(p, eps), 1 - eps)

.make_folds <- function(n, K, seed) {
  if (K <= 1L) return(rep(1L, n))
  withr::with_seed(seed, sample(rep_len(seq_len(K), n)))
}
