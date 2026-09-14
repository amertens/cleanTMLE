#' @keywords internal
#' @aliases cleanTMLE-package
"_PACKAGE"

#' cleanTMLE: Staged, Outcome-Blind Targeted Learning
#'
#' Decides, before outcome access, whether a comparison is estimable and
#' with which estimand, then estimates with targeted maximum likelihood.
#' The lock's data frame holds design data only; the primary outcome
#' lives in a separate store that Stage 4 joins back after the outcome
#' guard. The workflow surface is sixteen verbs; the estimator layer,
#' the learners, and the reporting helpers are separate, labelled
#' groups on the reference index.
#'
#' @section The sixteen workflow verbs:
#' * [create_analysis_lock()] - seal data, roles, declarations, and the
#'   fingerprint; [declare_negative_controls()] and
#'   [declare_estimand_ladder()] add the prespecified decision rules
#' * [fit_ps()], [assess_support()], [estimand_feasibility()],
#'   [simulate_support()] - the design-stage estimability answer
#' * [define_candidates()], [select_candidate()], [stress_test()] - the
#'   candidate grid, the outcome-blind selection, and the data-quality
#'   stress test with its locked verdict and tipping points
#' * [negative_control_ladder()], [check_process_indicators()] - the
#'   residual-confounding and collider checks
#' * [design_report()], [export_design_log()] - what the review team
#'   reads, and the released decision record
#' * [estimate_effect()], [run_estimand_ladder()] - Stage 4 estimation
#'   of the declared primary and its feasible fallbacks
#'
#' @section Estimator layer:
#' [run_clean_tmle()] is the classic single-call wrapper; the
#' cumulative-risk grammar ([specify_models()] with the `identify_*()`
#' verbs and [estimate_ipwrisk()], [estimate_gcomprisk()],
#' [estimate_aipwrisk()], [estimate_ipwhr()], [estimate_surv_tmle()],
#' [estimate_lmtp()]) covers time-to-event risks;
#' [select_variance_method()] and [bootstrap_rd_variance()] select and
#' supply the variance method. The modular TMLE steps remain internal,
#' reachable through `estimate_effect(return_steps = TRUE)`.
#'
#' @section Reporting and diagnostics:
#' [make_table1()], [make_table2()], [attrition_table()], [love_plot()],
#' [forest_plot()], [clean_weight_diagnostics()], [compute_evalue()],
#' [run_delta_sensitivity()], [event_support_by_arm()], and the verdict
#' helpers [dq_locked_verdict()], [dq_tipping_points()],
#' [nc_ladder_verdict()].
#'
#' @importFrom graphics hist
#' @importFrom stats as.formula approx binomial coef confint glm model.matrix
#'   predict quantile rbinom rnorm runif sd var vcov median
#'   pnorm qnorm qlogis weighted.mean terms reformulate setNames
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

# `.weights` is a column created on the model data frame and read by
# coxph()'s weights argument inside that frame; declare it so the
# static checker does not read it as an undefined global.
utils::globalVariables(".weights")
