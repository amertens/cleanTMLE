#' @keywords internal
#' @aliases cleanTMLE-package
"_PACKAGE"

#' cleanTMLE: Causal Inference for Risk Estimation
#'
#' A tidy interface for causal inference in time-to-event and binary outcome
#' settings. Supports inverse probability weighted (IPW) cumulative risk
#' estimation, g-computation, augmented IPW (AIPW / doubly robust), and
#' weighted Cox regression for hazard ratios. Integrates with TMLE, survtmle,
#' and lmtp packages for targeted learning extensions.
#'
#' @section Model specification:
#' Use [specify_models()] to create a specification object, then pipe it
#' through [identify_outcome()], [identify_treatment()], [identify_censoring()],
#' and other `identify_*` functions to build a complete causal model.
#'
#' @section Estimation:
#' * [estimate_ipwrisk()] - IPW cumulative risk curves
#' * [estimate_gcomprisk()] - G-computation risk curves
#' * [estimate_aipwrisk()] - Augmented IPW (doubly robust) risk curves
#' * [estimate_ipwhr()] - Weighted Cox model hazard ratios
#' * [estimate_ipwcount()] - IPW cumulative count (placeholder)
#'
#' @section TMLE extensions:
#' * [estimate_tmle_risk_point()] - Point-treatment TMLE for binary outcomes
#' * [estimate_surv_tmle()] - Survival TMLE for risk at specified times
#' * [estimate_lmtp()] - Longitudinal TMLE for static/dynamic interventions
#'
#' @section Diagnostics:
#' * [make_table1()] - Baseline covariate table (weighted/unweighted)
#' * [make_table2()] - Results table (N, person-time, events, risk, contrasts)
#' * [make_wt_summary_table()] - Weight distribution summaries
#' * [extreme_weights()] - Identify extreme weight observations
#' * [inspect_ipw_weights()] - Extract and inspect IPW weights
#'
#' @importFrom stats as.formula binomial coef confint glm model.matrix
#'   predict quantile rbinom rnorm runif sd var vcov median
#'   pnorm qnorm weighted.mean terms reformulate setNames
#'   rexp plogis complete.cases
#' @importFrom survival coxph survfit Surv strata
#' @importFrom ggplot2 ggplot aes geom_step geom_line geom_point
#'   geom_errorbar geom_errorbarh geom_histogram geom_vline
#'   geom_hline geom_ribbon geom_bar geom_text geom_segment
#'   facet_wrap labs theme_minimal theme_bw theme scale_color_manual
#'   scale_fill_manual coord_flip element_text element_blank
#'   xlim ylim ggtitle xlab ylab position_dodge guides
#'   guide_legend scale_x_continuous scale_y_continuous
#'   margin after_stat
#' @importFrom rlang enquo quo_name eval_tidy sym !! :=
#'   .data is_missing caller_env
#' @importFrom sandwich vcovHC
NULL
