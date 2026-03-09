#' @keywords internal
#' @aliases cleanTMLE-package
"_PACKAGE"

#' cleanTMLE: Staged Clean-Room Causal Analysis with Propensity Score and TMLE Workflows
#'
#' Provides workflow scaffolding for staged, outcome-blinded observational
#' causal analyses. The package distinguishes two pre-outcome feasibility paths:
#' Stage 2a (traditional PS diagnostics using only W and A) and Stage 2b
#' (plasmode-simulation evaluation of the full estimator pipeline). Treatment
#' mechanism estimation, outcome mechanism estimation, and the TMLE targeting
#' step are deliberately separated in the API so each can run at the correct
#' clean-room stage.
#'
#' @section Analysis lock and workflow metadata:
#' Use [create_analysis_lock()] to construct and serialize the complete analytic
#' specification (estimand, covariates, SuperLearner library, plasmode settings,
#' candidate set and selection rule). Use [validate_analysis_lock()] to verify
#' the specification before proceeding.
#'
#' @section Treatment mechanism and propensity-score estimation (Stage 2a / 2b):
#' * [fit_ps_superlearner()] - SuperLearner-based PS estimation (default)
#' * [fit_ps_glm()] - logistic-regression PS estimation (conventional option)
#' * [compute_ps_diagnostics()] - overlap, ESS, covariate balance summaries
#'
#' @section Conventional workflow estimators (Stage 3):
#' * [run_match_workflow()] - 1:1 nearest-neighbor matching
#' * [run_iptw_workflow()] - stabilized IPTW with weight diagnostics
#'
#' @section Plasmode simulation and feasibility evaluation (Stage 2b):
#' * [run_plasmode_feasibility()] - generate plasmode outcomes and evaluate
#'   the full pipeline; returns bias, RMSE, and coverage
#' * [summarize_plasmode_results()] - tabulate and plot simulation metrics
#' * [evaluate_tmle_candidates()] - evaluate all candidate TMLE specifications
#'
#' @section Modular TMLE components:
#' These functions respect clean-room stage separation. The treatment mechanism
#' can run in Stage 2a or 2b; the outcome mechanism runs on real data only in
#' Stage 3; the targeting step follows only after both nuisance estimates exist.
#' * [fit_tmle_treatment_mechanism()] - estimate g(W) from W and A only
#' * [fit_tmle_outcome_mechanism()] - estimate Q(A,W) using the real outcome
#' * [run_tmle_targeting_step()] - perform the fluctuation / targeting update
#' * [extract_tmle_estimate()] - compute psi, SE, CI, and diagnostics
#'
#' @section Candidate TMLE comparison / selection:
#' * [fit_tmle_candidate_set()] - fit all prespecified TMLE candidates
#' * [select_tmle_candidate()] - apply the prespecified selection rule
#'
#' @section Lower-level estimation utilities:
#' * [estimate_ipwrisk()] - IPW cumulative risk curves
#' * [estimate_gcomprisk()] - G-computation risk curves
#' * [estimate_aipwrisk()] - Augmented IPW (doubly robust) risk curves
#' * [estimate_ipwhr()] - Weighted Cox model hazard ratios
#' * [estimate_ipwcount()] - IPW cumulative count
#'
#' @section TMLE package integrations:
#' * [estimate_tmle_risk_point()] - Point-treatment TMLE for binary outcomes
#' * [estimate_surv_tmle()] - Survival TMLE for risk at specified times
#' * [estimate_lmtp()] - Longitudinal TMLE for static/dynamic interventions
#'
#' @section Workflow-level summaries:
#' * [fit_final_workflows()] - run final estimation after outcome unblinding
#' * [summarize_cleanroom_results()] - side-by-side summary across workflows
#' * [compare_fits()] - compare estimates, CIs, and diagnostic flags
#'
#' @section Diagnostics and tables:
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
