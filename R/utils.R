`%||%` <- function(x, y) if (is.null(x)) y else x

.hashable <- function(x, descend = TRUE) {
  if (is.function(x)) {
    out <- list(text = paste(deparse(x), collapse = "\n"))
    env <- environment(x)
    # Hash a closure's captured variables, but not the global, base or empty
    # environments or a package namespace. A function held inside a captured
    # environment contributes only its text, which stops any recursion.
    if (descend && is.environment(env) &&
        !identical(env, globalenv()) && !identical(env, baseenv()) &&
        !identical(env, emptyenv()) && !isNamespace(env)) {
      vars <- as.list(env, all.names = TRUE)
      if (length(vars)) {
        vars <- vars[order(names(vars))]
        out$env <- lapply(vars, .hashable, descend = FALSE)
      }
    }
    return(if (length(out) == 1L) out$text else out)
  }
  if (inherits(x, "formula")) return(paste(deparse(x), collapse = " "))
  if (is.environment(x)) return("<environment>")
  if (is.list(x)) {
    y <- lapply(x, .hashable, descend = descend)
    attributes(y) <- NULL
    names(y) <- names(x)
    return(y)
  }
  x
}

.hash <- function(x) digest::digest(.hashable(x), algo = "sha256")

.bound <- function(p, eps) pmin(pmax(p, eps), 1 - eps)

# Every seeded step runs under R's default generators, whatever RNG kind the
# caller (or a parallel worker, which future puts on L'Ecuyer-CMRG) has set,
# so a seed gives the same draws, folds and fits everywhere. The caller's
# RNG state and kind are restored afterwards.
.with_seed <- function(seed, code) {
  withr::with_seed(seed, code, .rng_kind = "Mersenne-Twister",
                   .rng_normal_kind = "Inversion", .rng_sample_kind = "Rejection")
}

.make_folds <- function(n, K, seed) {
  if (K <= 1L) return(rep(1L, n))
  .with_seed(seed, sample(rep_len(seq_len(K), n)))
}
