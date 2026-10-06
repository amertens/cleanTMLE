#' Example cohort for cleanTMLE
#'
#' A simulated cohort of 1,500 patients, shipped as two data frames so the
#' outcome can be admitted with [unblind()] after the design stage.
#'
#' @format A list with two data frames joined by `id`:
#' \describe{
#'   \item{design}{`id`, treatment `A` (0/1), covariates `age`, `female`,
#'     `severity`, `prior_visits`, restriction variable `transfer`, and the
#'     negative-control outcome `health_seeking`.}
#'   \item{outcomes}{`id` and the 0/1 outcome `death_30d`.}
#' }
#' @source Simulated by `data-raw/make_cr_example.R`.
"cr_example"
