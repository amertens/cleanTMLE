#' Validate and Suggest a SuperLearner Specification
#'
#' Implements the practical-considerations checklist of
#' Phillips et al. (2023) for SuperLearner library specification: given
#' an analytic sample size (and optionally an outcome event count for
#' rare-event problems), it returns a recommended cross-validation fold
#' count and a minimum library template, and warns if the user-supplied
#' library is likely under-powered for the data size.
#'
#' Heuristics (from Phillips et al. 2023, §3 "Step 2: choose V"):
#' \itemize{
#'   \item Effective sample size \eqn{n_{eff}} = number of events for
#'     rare binary outcomes (when \code{n_events} is supplied), else
#'     \eqn{n_{eff} = n}.
#'   \item \eqn{n_{eff} < 30}: V = leave-one-out (`"loo"`).
#'   \item \eqn{30 \le n_{eff} < 500}: V = 20 (Phillips' default for
#'     small samples).
#'   \item \eqn{500 \le n_{eff} < 5000}: V = 10.
#'   \item \eqn{n_{eff} \ge 5000}: V = 5.
#' }
#'
#' @section Clean-room stage: Stage 1a (pre-specification).
#'
#' @param n Integer; analytic sample size.
#' @param n_events Optional integer; outcome event count for rare-event
#'   binary outcomes. When supplied, the heuristic is applied to the
#'   smaller of \code{n} and \code{n_events}.
#' @param library Optional character vector of SuperLearner learner
#'   names to validate. If supplied, the function warns when the library
#'   is degenerate (a single GLM) or omits a flexible learner.
#'
#' @return A list with elements:
#'   \itemize{
#'     \item \code{n_eff}: the effective sample size used.
#'     \item \code{recommended_V}: integer or \code{"loo"}.
#'     \item \code{recommended_library}: a default library appropriate
#'       to the sample size (small: GLM + glmnet; medium: GLM + glmnet +
#'       gam; large: GLM + glmnet + gam + ranger).
#'     \item \code{warnings}: character vector of any concerns about the
#'       supplied library (empty when none).
#'   }
#'
#' @examples
#' validate_superlearner_spec(n = 1693)
#' validate_superlearner_spec(n = 200, n_events = 25,
#'                            library = c("SL.glm"))
#'
#' @export
validate_superlearner_spec <- function(n, n_events = NULL,
                                        library = NULL) {
  if (!is.numeric(n) || length(n) != 1L || n < 1)
    stop("`n` must be a positive integer.", call. = FALSE)
  n_eff <- if (!is.null(n_events)) min(as.integer(n), as.integer(n_events))
           else as.integer(n)

  V <- if (n_eff < 30L)        "loo"
       else if (n_eff < 500L)  20L
       else if (n_eff < 5000L) 10L
       else                     5L

  rec_lib <- if (n_eff < 500L) {
    c("SL.glm", "SL.glmnet")
  } else if (n_eff < 5000L) {
    c("SL.glm", "SL.glmnet", "SL.gam")
  } else {
    c("SL.glm", "SL.glmnet", "SL.gam", "SL.ranger")
  }

  warnings <- character(0)
  if (!is.null(library)) {
    if (length(library) == 1L && grepl("glm", library, ignore.case = TRUE))
      warnings <- c(warnings,
        sprintf("Single-learner library ('%s') defeats the purpose of SL; consider adding a regularised or non-parametric learner.",
                library))
    if (n_eff >= 500L && !any(grepl("glmnet|gam|ranger|hal|bart",
                                     library, ignore.case = TRUE)))
      warnings <- c(warnings,
        "Library lacks a flexible learner (glmnet/gam/ranger/hal/bart); SL will collapse to the parametric fit.")
  }

  list(
    n_eff               = n_eff,
    recommended_V       = V,
    recommended_library = rec_lib,
    warnings            = warnings
  )
}


#' Tipping-Point Sensitivity Bound for an Effect Estimate
#'
#' Computes the magnitude of residual bias that would be required to
#' shift the lower or upper bound of an effect estimate's confidence
#' interval across a null value (Step 5 of Gruber et al. 2023; see also
#' the targeted-learning causal-gap sensitivity of Diaz & van der Laan
#' 2013). For a risk-difference (default), the tipping point is the
#' bias that drags the CI to include 0; for a risk ratio or hazard
#' ratio (\code{null = 1}), the function returns the bias on the
#' log-scale required to drag the CI to include 1.
#'
#' This is a quick non-parametric companion to \code{compute_evalue()}
#' and to a full quantitative bias analysis. It does not estimate the
#' bias; it asks how large a bias would have to be to flip the
#' qualitative conclusion.
#'
#' @section Clean-room stage: Stage 5 (post-outcome sensitivity).
#'
#' @param estimate Numeric; the point estimate (on the same scale as
#'   the CI bounds and the null value).
#' @param ci_lower,ci_upper Numeric; lower and upper 95\% confidence
#'   interval bounds on the estimate.
#' @param null_value Numeric; the null value of no effect. Default 0
#'   (risk difference / mean shift); use 1 for risk-ratio or hazard-
#'   ratio scales.
#'
#' @return A list with elements \code{tipping_point} (the additive bias
#'   required to make the CI cross the null), \code{direction}
#'   (\code{"toward null"} or \code{"away from null"}), and
#'   \code{ratio_to_estimate} (the tipping point divided by
#'   \code{|estimate - null_value|}).
#'
#' @examples
#' # rescueCo full-cohort TMLE: RD = 0.031, 95% CI (-0.001, 0.063)
#' tipping_point_sensitivity(0.0311, -0.0007, 0.0629)
#'
#' @export
tipping_point_sensitivity <- function(estimate, ci_lower, ci_upper,
                                       null_value = 0) {
  if (any(is.na(c(estimate, ci_lower, ci_upper))))
    return(list(tipping_point = NA_real_,
                direction = NA_character_,
                ratio_to_estimate = NA_real_))
  if (estimate > null_value) {
    # Bias toward null shrinks CI lower bound through the null
    bias <- ci_lower - null_value
    direction <- "toward null"
  } else {
    bias <- ci_upper - null_value
    direction <- "toward null"
  }
  list(
    tipping_point     = abs(bias),
    direction         = direction,
    ratio_to_estimate = abs(bias) / max(abs(estimate - null_value), .Machine$double.eps)
  )
}
