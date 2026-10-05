make_design <- function(n = 400, seed = 1, strength = 1) {
  withr::with_seed(seed, {
    w1 <- stats::rnorm(n); w2 <- stats::rbinom(n, 1, 0.5); w3 <- stats::rnorm(n)
    A <- stats::rbinom(n, 1, stats::plogis(-0.2 + strength * w1 + 0.3 * w2))
    nc_visit <- stats::rbinom(n, 1, stats::plogis(-0.8 + 0.4 * w1))
    nc_bad <- stats::rbinom(n, 1, stats::plogis(-1 + 1.5 * A))
    transfer <- stats::rbinom(n, 1, 0.1)
    y <- stats::rbinom(n, 1, stats::plogis(-1 + 0.5 * A + 0.5 * w1 + 0.3 * w3))
    list(design = data.frame(id = seq_len(n), A = A, w1 = w1, w2 = w2, w3 = w3,
                             nc_visit = nc_visit, nc_bad = nc_bad,
                             transfer = transfer),
         outcomes = data.frame(id = seq_len(n), y = y))
  })
}

fast_plan <- function(...) {
  defaults <- list(outcome = "y",
                   estimands = c("ATE", "trimmed_ATE", "ATT", "ATO"),
                   ps_library = "glm",
                   candidates = list(library = list(glm = "glm"),
                                     truncation = c(0.01, 0.05)),
                   trim_band = c(0.05, 0.95),
                   K = 2L, V = 2L, reps = 20L, seed = 11L)
  do.call(analysis_plan, utils::modifyList(defaults, list(...)))
}
