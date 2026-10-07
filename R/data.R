#' Example cohort for cleanTMLE
#'
#' A simulated cohort of 1,500 patients, shipped as separate design and
#' outcome data frames so the outcome can be admitted with [unblind()] after
#' the design stage. Treatment depends strongly on severity, so few severe
#' patients are untreated and few mild patients are treated, and the effect
#' of treatment grows with severity. Under a plan like the one in the Get
#' started vignette, the average treatment effect is not estimable from
#' these data while an effect in a narrower population is.
#'
#' @format A list with three data frames:
#' \describe{
#'   \item{design}{`id`, treatment `A` (0/1), covariates `age`, `female`,
#'     `severity`, `prior_visits`, restriction variable `transfer`, and the
#'     negative-control outcome `health_seeking`.}
#'   \item{outcomes}{`id` and the 0/1 outcome `death_30d`.}
#'   \item{truth}{`estimand`, `truth` (the true risk difference, computed
#'     from the generating model over this cohort's covariates) and
#'     `population` (whom it describes).}
#' }
#' @source Simulated by `make_cr_example()` in `data-raw/make_cr_example.R`.
"cr_example"
