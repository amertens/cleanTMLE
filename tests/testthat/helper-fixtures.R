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
  defaults <- list(outcome = "y", time_zero = "treatment assignment",
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
                             seed = 11L, tolerance = list(bias = 0.1, coverage = 0.5),
                             max_reps = NULL, nc_rule = "point", nc_exclude = NULL,
                             dictionary_edit = identity, ...) {
  d <- make_design(n, seed = 2, strength = strength)
  plan <- fast_plan(negative_controls = nc,
                    nc_criteria = if (is.null(nc)) NULL else list(null_band = c(-0.1, 0.1), rule = nc_rule),
                    tolerance = tolerance, reps = 10L, max_reps = max_reps, seed = seed,
                    nc_exclude = nc_exclude, ...)
  lock <- lock_fx(d$design, "A", c("w1", "w2", "w3"), plan,
                  dictionary = dictionary_edit(
                    fx_dictionary(d$design, "A", c("w1", "w2", "w3"), plan)))
  design <- assess_design(lock)
  sim <- suppressMessages(simulate_design(lock, design))
  ncl <- negative_control_ladder(lock, design, sim)
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

# The small plan the Get-started vignette runs live.
vignette_plan <- function(...) {
  args <- utils::modifyList(list(
    outcome = "death_30d", time_zero = "treatment assignment",
    estimands = c("ATE", "trimmed_ATE", "ATT", "ATO"),
    ps_library = "glm", candidates = list(library = list(glm = "glm"), truncation = 0.01),
    surfaces = list(forms = "linear"), trim_band = c(0.05, 0.95),
    negative_controls = c(health_seeking = "health care use"),
    restrictions = list(no_transfer = ~ transfer == 0),
    K = 1L, reps = 50L, max_reps = 50L, seed = 1L), list(...))
  do.call(analysis_plan, args)
}

# analysis_plan() with time_zero filled in, for tests that are not about time zero.
plan_fx <- function(...) {
  a <- list(...)
  if (is.null(a$time_zero)) a$time_zero <- "treatment assignment"
  do.call(analysis_plan, a)
}

# A dictionary for test fixtures: one row per kept column that exists in
# design_data, recorded and fixed before time zero (the treatment at it), with
# the role the column has in the call. Internal to the tests, because a
# generated dictionary would defeat the check for a user.
fx_dictionary <- function(design_data, treatment, covariates, plan, id = "id") {
  if (!is.data.frame(design_data)) return(NULL)
  nc <- names(plan$negative_controls)
  rv <- unique(unlist(lapply(plan$restrictions, all.vars)))
  cols <- intersect(unique(c(id, treatment, covariates, nc, rv)), names(design_data))
  role <- vapply(cols, function(v) {
    if (identical(v, id)) "id"
    else if (identical(v, treatment)) "treatment"
    else if (v %in% covariates) "covariate"
    else if (v %in% nc) "negative_control"
    else "restriction"
  }, character(1))
  data.frame(column = cols, variable = cols, source = "test fixture",
             follow_up = FALSE, recorded = "before",
             fixed = ifelse(cols == treatment, "at", "before"),
             role = unname(role), reason = "", sensitivity = "",
             stringsAsFactors = FALSE)
}

# create_analysis_lock() with the fixture dictionary unless one is given.
lock_fx <- function(design_data, treatment, covariates, plan, id = "id",
                    dictionary = NULL) {
  if (is.null(dictionary))
    dictionary <- fx_dictionary(design_data, treatment, covariates, plan, id)
  create_analysis_lock(design_data, treatment, covariates, plan, id = id,
                       dictionary = dictionary)
}

# A pipeline whose plan raises three identification flags (follow-up start,
# eligibility assessed after, stand-in covariate).
id_pipeline <- function() {
  fixture_pipeline(
    time_zero = list(event = "transport decision", measured_as = "time of injury",
                     follow_up_from = "arrival", response = "r0"),
    eligibility = list(arrived_alive = list(description = "Arrived alive",
                                            assessed = "after", response = "r1")),
    dictionary_edit = function(dict) {
      dict$recorded[dict$column == "w1"] <- "after"
      dict$reason[dict$column == "w1"] <- "x"
      i <- dict$column == "w3"
      dict$role[i] <- "stand_in"; dict$fixed[i] <- "after"
      dict$reason[i] <- "r2"; dict$sensitivity[i] <- "s2"
      extra <- data.frame(column = "nc_bad", variable = "nc_bad", source = "test fixture",
                          follow_up = FALSE, recorded = "after", fixed = "after",
                          role = "after_treatment", reason = "mediator of interest",
                          sensitivity = "", stringsAsFactors = FALSE)
      rbind(dict, extra)
    })
}
