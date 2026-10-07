#' Example cohort for cleanTMLE
#'
#' A simulated cohort of 1,500 patients, shipped as separate design and
#' outcome data frames so the outcome can be admitted with [unblind()] after
#' the design stage. Treatment depends strongly on severity, so few severe
#' patients are untreated and few mild patients are treated, and the effect
#' of treatment grows with severity. Under a small plan like the one the
#' package's tests use (one learner, no cross-fitting, 50 repetitions), the
#' outcome-free simulation finds the average treatment effect infeasible
#' against the plan's tolerances while the trimmed average treatment effect
#' and the overlap-weighted effect (ATO) are feasible. Under the default plan
#' with 50 repetitions only the ATO is feasible.
#'
#' @format A list with three data frames:
#' \describe{
#'   \item{design}{`id`, treatment `A` (0/1), covariates `age`, `female`,
#'     `severity`, `prior_visits`, restriction variable `transfer`, and the
#'     negative-control outcome `health_seeking`.}
#'   \item{outcomes}{`id` and the 0/1 outcome `death_30d`.}
#'   \item{truth}{`estimand`, `truth` (the true risk difference, computed
#'     from the generating model over this cohort's covariates) and
#'     `population` (whom it describes). The trimmed-ATE truth is over the
#'     patients whose true propensity score lies in `[0.05, 0.95]`; that band
#'     is fixed at those values. [estimate_effect()] trims on the
#'     cross-fitted estimated score, so its trimmed population can differ
#'     from this one. The ATT and ATO truths weight patients by the true
#'     propensity score.}
#' }
#' @source Simulated by `make_cr_example()` in `data-raw/make_cr_example.R`.
"cr_example"
