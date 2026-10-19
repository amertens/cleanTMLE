# Builds cr_example: design data (no outcome), outcomes, and the true value of
# each estimand, joined by id. `assign` is the log-odds of treatment per SD of
# severity; `modifier` is how much the effect on the log-odds scale grows per
# SD of severity. Both were chosen so that, under the vignette's plan, the ATE
# is infeasible while a lower rung of the ladder is feasible.
make_cr_example <- function(assign, modifier, seed = 20261005, n = 1500) {
  set.seed(seed)
  age <- round(stats::rnorm(n, 60, 12))
  female <- stats::rbinom(n, 1, 0.5)
  severity <- stats::rnorm(n)
  prior_visits <- stats::rpois(n, 2)
  transfer <- stats::rbinom(n, 1, 0.1)
  ps <- stats::plogis(-0.3 + assign * severity + 0.02 * (age - 60) - 0.3 * female)
  A <- stats::rbinom(n, 1, ps)
  health_seeking <- stats::rbinom(n, 1, stats::plogis(-0.7 + 0.3 * prior_visits))
  risk <- function(a) stats::plogis(-2 + (0.4 + modifier * severity) * a + 0.8 * severity +
                                      0.03 * (age - 60) + 0.3 * severity^2)
  death_30d <- stats::rbinom(n, 1, risk(A))
  tau <- risk(1) - risk(0)
  inb <- ps >= 0.05 & ps <= 0.95
  h <- ps * (1 - ps)
  list(
    design = data.frame(id = seq_len(n), A = A, age = age, female = female,
                        severity = severity, prior_visits = prior_visits,
                        transfer = transfer, health_seeking = health_seeking),
    outcomes = data.frame(id = seq_len(n), death_30d = death_30d),
    truth = data.frame(
      estimand = c("ATE", "trimmed_ATE", "ATT", "ATO"),
      truth = c(mean(tau), mean(tau[inb]), sum(ps * tau) / sum(ps), sum(h * tau) / sum(h)),
      population = c("all patients",
                     "patients whose true propensity score lies in [0.05, 0.95]",
                     "patients weighted by their true propensity score",
                     "patients weighted by ps(1 - ps), ps the true propensity score"),
      stringsAsFactors = FALSE),
    dictionary = cr_example_dictionary())
}

# One row per design column. Everything is recorded on the enrolment form before
# treatment is assigned and describes an attribute fixed before it, except the
# treatment, which is assigned at time zero.
cr_example_dictionary <- function() {
  cols <- c("id", "A", "age", "female", "severity", "prior_visits", "transfer",
            "health_seeking")
  data.frame(
    column = cols, variable = cols, source = "enrolment form", follow_up = FALSE,
    recorded = "before", fixed = ifelse(cols == "A", "at", "before"),
    role = c("id", "treatment", rep("covariate", 4), "restriction", "negative_control"),
    reason = "", sensitivity = "", stringsAsFactors = FALSE)
}

if (sys.nframe() == 0L) {
  cr_example <- make_cr_example(assign = 3.0, modifier = 0.3)
  usethis::use_data(cr_example, overwrite = TRUE)
}
