# Builds cr_example: design data (no outcome) and outcomes, joined by id.
set.seed(20261005)
n <- 1500
age <- round(stats::rnorm(n, 60, 12))
female <- stats::rbinom(n, 1, 0.5)
severity <- stats::rnorm(n)
prior_visits <- stats::rpois(n, 2)
transfer <- stats::rbinom(n, 1, 0.1)
A <- stats::rbinom(n, 1, stats::plogis(-0.3 + 1.6 * severity + 0.02 * (age - 60) - 0.3 * female))
health_seeking <- stats::rbinom(n, 1, stats::plogis(-0.7 + 0.3 * prior_visits))
death_30d <- stats::rbinom(n, 1, stats::plogis(-2 + 0.4 * A + 0.8 * severity +
                                                 0.03 * (age - 60) + 0.3 * severity^2))
cr_example <- list(
  design = data.frame(id = seq_len(n), A = A, age = age, female = female,
                      severity = severity, prior_visits = prior_visits,
                      transfer = transfer, health_seeking = health_seeking),
  outcomes = data.frame(id = seq_len(n), death_30d = death_30d))
usethis::use_data(cr_example, overwrite = TRUE)
