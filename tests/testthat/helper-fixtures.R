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
  args <- utils::modifyList(defaults, list(...))
  # No extra repetitions unless a test asks for them, so tests stay cheap.
  if (is.null(args$max_reps)) args$max_reps <- args$reps
  do.call(analysis_plan, args)
}

fixture_pipeline <- function(strength = 0.8, n = 400, nc = c(nc_visit = "care use"),
                             seed = 11L) {
  d <- make_design(n, seed = 2, strength = strength)
  plan <- fast_plan(negative_controls = nc,
                    nc_criteria = if (is.null(nc)) NULL else list(null_band = c(-0.1, 0.1)),
                    tolerance = list(bias = 0.1, coverage = 0.5), reps = 10L, seed = seed)
  lock <- create_analysis_lock(d$design, "A", c("w1", "w2", "w3"), plan)
  design <- assess_design(lock)
  sim <- suppressMessages(simulate_design(lock, design))
  ncl <- negative_control_ladder(lock, design)
  list(d = d, plan = plan, lock = lock, design = design, sim = sim, nc = ncl,
       dossier = design_report(lock, design, sim, ncl))
}

contains_outcome <- function(obj, name, y) {
  if (is.environment(obj) || is.function(obj)) return(FALSE)
  if (is.list(obj)) {
    if (name %in% names(obj)) return(TRUE)
    return(any(vapply(obj, contains_outcome, logical(1), name = name, y = y)))
  }
  is.numeric(obj) && length(obj) == length(y) &&
    isTRUE(all.equal(as.numeric(obj), as.numeric(y), check.attributes = FALSE))
}
