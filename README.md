# cleanTMLE

**Outcome-Blind Staged Workflow for TMLE**

cleanTMLE provides software infrastructure for an outcome-blind
staged workflow in targeted-learning analyses. It helps analysts
define an analysis lock, specify candidate TMLE estimators, run
design diagnostics, compare candidates using baseline plasmode
simulation, stress-test candidates under prespecified
data-quality threats, record GO / FLAG / STOP checkpoint
decisions, and authorise the locked primary analysis.

## What problem does cleanTMLE address?

In RWE studies, protocol design, fit-for-purpose data review, and
estimator choice are often documented separately. Standard
diagnostics assess overlap, balance, and effective sample size,
but they do not characterise the full operating characteristics
of a candidate TMLE specification on the cohort at hand.
cleanTMLE provides an auditable software layer for pre-outcome
estimator evaluation and prespecified data-quality stress
testing, with a structured record of the analysis lock, the
diagnostics, the candidate comparison, and the GO / FLAG / STOP
decisions. It is intended to support, not replace, careful
causal design and broader clean-room governance.

## Pre-outcome study dossier

cleanTMLE's central artefact is a reviewer-facing pre-outcome study dossier: a structured bundle of analytic outputs the review team reads before authorising the primary analysis. The dossier contains:

- Protocol and target-trial timing (`attach_estimand()`)
- Cohort flow and attrition (`attrition_table()`)
- Event-process classification (`clean_event_process_table()`, `clean_check_event_processes()`)
- Covariate and missingness summary (`make_table1()`, missingness checks)
- PS overlap and balance (`compute_ps_diagnostics()`, `love_plot()`, `love_plot_threeway()`)
- Treatment- and IPCW-weight diagnostics (`clean_weight_diagnostics()`, `extreme_weights()`)
- Event-count and precision adequacy (`checkpoint_cohort_adequacy()`)
- Baseline plasmode candidate selection (`run_plasmode_feasibility()`, `select_tmle_candidate()`)
- DQ stress-test results (`run_plasmode_dq_stress()`, `summarize_dq_degradation()`)
- Negative-control eligibility and attrition (`run_residual_confounding_stage()`)
- GO / FLAG / STOP checkpoint dashboard (`gate_all()`; aggregate `checkpoint_dashboard()` planned)
- Decision log and audit log export (`export_decision_log()`)
- Authorisation record for primary analysis (`authorize_outcome_analysis()`)

The dossier is the artefact the reviewer reads. Comparative treatment-outcome estimates are not part of the dossier and are produced only after authorisation.

## Workflow

```
Stage 0   Target-trial / protocol specification (external precondition)
Stage 1a  Lock estimand, candidates, learners, truncation, thresholds
Stage 1b  Cohort adequacy and marginal event support
Stage 2a  Design diagnostics (PS, overlap, balance, ESS)
Stage 2b  Baseline plasmode candidate selection
Stage 2c  Data-quality stress testing
Stage 3   Optional negative-control checks
  Gate    Pre-outcome decision (GO / FLAG / STOP)
Stage 4   Authorised primary analysis
Post-outcome: sensitivity analyses and interpretation
```

The data-quality stress test is a quantitative pre-outcome
supplement to fit-for-purpose data review, not a formal QBA and
not a substitute for source-data validation.

## What cleanTMLE does not do

- does not automate target-trial design or protocol specification
- does not validate source data
- does not validate phenotypes or outcome definitions
- does not establish exchangeability
- does not establish positivity
- does not guarantee that synthetic-outcome rankings generalise
  to the realised outcome process
- does not implement personnel role separation
- does not control raw-data access
- does not replace negative-control or post-outcome sensitivity
  analyses
- does not perform formal identification-region quantitative
  bias analysis

## What external governance must still provide

The software workflow is necessary but not sufficient for high-stakes RWE. External governance, outside cleanTMLE's scope, must provide:

- Role separation between analytic, methodological, and reviewing teams
- Raw-data access controls (credentialed folders, audit trails of file access)
- Review-team structure with named reviewers and an escalation path
- Protocol registration in an appropriate registry before data access
- Source-data and phenotype validation against external references
- Independent review of deviations and overrides recorded in the decision log

## When to use this package

- methods development around outcome-blind workflows
- outcome-blind estimator selection in observational studies
- protocol and SAP development for RWE analyses
- exploratory or lower-stakes RWE workflows where structured
  pre-outcome review is useful
- high-stakes confirmatory RWE only when embedded inside
  external clean-room governance (role separation, data-access
  controls, independent checkpoint review, archived audit trail)

## Cumulative-risk workflow considerations

cleanTMLE provides software support for the analytic considerations that arise in cumulative-risk pharmacoepidemiology workflows:

- **Model specification grammar**: `identify_*()` family and `specify_models()` for declaring eligibility, treatment, outcome, censoring, competing risk, follow-up interval, and intercurrent events as code objects.
- **Cumulative-risk reporting**: `clean_risk_report_table()` produces a compact risk table at clinically meaningful time points; the package emphasises risk differences and risk ratios over hazard ratios except where the estimand is explicitly defined on the hazard scale.
- **Censoring and missingness weights as first-class objects**: `run_ipcw_tmle()` and `clean_weight_diagnostics()` expose ESS, percentiles, maximum weight, and prespecified instability thresholds.
- **Event-process classification**: `clean_event_process_table()` and `clean_check_event_processes()` distinguish event of interest, competing event, censoring, treatment discontinuation/switching, transfer exclusions, and administrative end of follow-up.
- **Hazard-ratio de-emphasis**: cumulative risks at clinically meaningful follow-up times are reported as primary; hazard ratios are reserved for analyses where a hazard-scale estimand is the primary scientific question and proportional hazards is plausible.

## Overview

`cleanTMLE` is the software layer for an outcome-blind staged
workflow around targeted minimum loss-based estimation (TMLE).
The estimand, nuisance-model specifications, and
estimator-selection rules are recorded in an analysis lock before
the observed primary treatment-outcome association is read by the
package. The package sits within, but does not replace, the
broader clean-room governance construct of Muntner et al. (2020),
which also covers role separation, restricted data access, and
independent checkpoint review.

Each checkpoint records a structured **GO / FLAG / STOP** decision.
Subsequent stages run only after all preceding checkpoints have
been evaluated. An audit trail accumulates entries automatically
and can be exported for review.

The package covers three workflow families:

- **Conventional propensity-score methods** (matching, IPTW) as
  secondary comparators.
- **Fixed-specification TMLE** with a prespecified nuisance strategy.
- **Simulation-selected TMLE** in which candidate TMLE specifications
  (varying truncation and/or learner library) are evaluated on
  outcome-blind plasmode simulations, and a prespecified rule selects
  the best specification before the real outcome is accessed.

The **v0.1 tested scope** is intentionally narrow: binary point
exposure, binary outcome, marginal risk difference, outcome
missingness handled through complete-case and IPCW sensitivity
paths under prespecified missingness assumptions. The package
also ships a **model-specification DSL** and a set of
**time-to-event helpers** (`estimate_ipwrisk`, `estimate_gcomprisk`,
`estimate_aipwrisk`, `estimate_ipwhr`, `estimate_surv_tmle`,
`estimate_lmtp`), but these are **experimental** and are not part
of the v0.1 tested scope; see *Experimental / planned extensions*
below.

## Key Features

- **Estimand-first design** --- declare the causal question, population,
  contrast, and follow-up window before any modelling
  (`attach_estimand()`)
- **Analysis lock** --- record and validate the full analytic
  specification (`create_analysis_lock()`, `validate_analysis_lock()`)
- **Staged checkpoints with GO / FLAG / STOP decisions**:
  - *Check Point 1*: cohort adequacy (`checkpoint_cohort_adequacy()`)
  - *Check Point 2*: covariate balance after PS weighting
    (`checkpoint_balance()`)
  - *Check Point 3*: residual bias via negative controls
    (`checkpoint_residual_bias()`)
- **Design-stage precision** --- estimate power and minimum detectable
  difference before outcome modelling
  (`estimate_design_precision()`, `summarize_event_support()`)
- **Sensitivity and negative control plans** --- declare before outcome
  access (`declare_sensitivity_plan()`, `define_negative_control()`)
- **Residual confounding wrapper** --- runs all registered negative
  controls and produces a unified Stage 3 result
  (`run_residual_confounding_stage()`)
- **Pre-outcome authorization gate** --- formal GO / STOP decision
  verifying all checkpoints passed before outcome access
  (`authorize_outcome_analysis()`, `assert_outcome_authorized()`)
- **Outcome masking** --- optionally mask the outcome column with `NA`
  during design stages and restore before estimation
  (`mask_outcome()`, `unmask_outcome()`)
- **Stage 4 outcome guard** --- all Stage 4 functions check for outcome
  masking and refuse to run on masked data unless
  `override_clean_room = TRUE` is set (the argument name is
  historical; the check records the outcome-blind state inside the package)
- **SuperLearner-based propensity-score estimation and diagnostics** ---
  ensemble learning PS (`fit_ps_superlearner()`) or logistic regression
  (`fit_ps_glm()`); overlap plots, effective sample size, and
  standardized mean differences (`compute_ps_diagnostics()`)
- **Matching and IPTW workflows** --- 1:1 nearest-neighbour matching
  (`run_match_workflow()`) and stabilised IPTW (`run_iptw_workflow()`)
- **TMLE candidate specification and selection** ---
  define TMLE specs varying truncation / library (`tmle_candidate()`,
  `expand_tmle_candidate_grid()`), evaluate on plasmode simulations
  (`run_plasmode_feasibility()`), select via prespecified rule
  (`select_tmle_candidate()`, with rules `min_rmse`, `min_bias`,
  `max_coverage`, and `min_max_rmse` for minimax RMSE across DQ
  stress scenarios), lock the chosen spec (`lock_primary_tmle_spec()`)
- **Plasmode data-quality stress test** --- extends the outcome-blind
  plasmode loop with four degradation mechanisms (covariate
  missingness, treatment misclassification, outcome misclassification,
  unmeasured confounding) and produces per-candidate degradation
  gradients of bias, RMSE, and coverage
  (`run_plasmode_dq_stress()`, `summarize_dq_degradation()`)
- **Gate decision** --- structured GO / FLAG / STOP based on bias,
  coverage, and SE calibration from plasmode results (`gate_check()`)
- **TMLE in four explicit steps** --- intentionally separated so each
  step runs at the correct workflow stage:
  1. Treatment mechanism / g-step (`fit_tmle_treatment_mechanism()`)
  2. Outcome mechanism / Q-step (`fit_tmle_outcome_mechanism()`)
  3. Targeting / fluctuation step (`run_tmle_targeting_step()`)
  4. Final estimate extraction (`extract_tmle_estimate()`)
- **Matched-cohort TMLE** --- `run_matched_tmle()` runs the four-step
  pipeline on a matched subset without creating a separate lock
- **IPCW-weighted TMLE** --- `run_ipcw_tmle()` adds inverse-probability-
  of-censoring weights for missing-at-random outcomes; targets the
  full-cohort marginal estimand rather than the complete-case subset
- **Convenience wrappers** --- `fit_final_workflows()` runs matching,
  IPTW, and TMLE in a single call; `fit_tmle_candidate_set()` fits
  multiple TMLE specifications on real data
- **Sensitivity analysis** --- truncation sensitivity
  (`sensitivity_truncation()`), E-value for unmeasured confounding
  (`compute_evalue()`)
- **Negative control analysis** --- `run_negative_control()` estimates
  the treatment effect on a variable known to be unaffected by treatment
- **Audit trail** --- `create_audit_log()`, `record_stage()`,
  `record_checkpoint()`, `export_audit_trail()`,
  `build_stage_manifest()`
- **Decision log** --- structured record of analyst decisions and
  protocol deviations (`record_decision_log_entry()`,
  `export_decision_log()`)
- **Stage path narrative** --- compact summary of the analysis path
  (`summarize_stage_path()`)
- **Cross-workflow comparison** --- `summarize_cleanroom_results()`
  produces a side-by-side table of estimates across all fitted workflows
- **Reporting helpers** --- `make_table1()`, `make_table2()`,
  `make_wt_summary_table()`, `extreme_weights()`, `compare_fits()`,
  `forest_plot()`

## Experimental / planned extensions

The following are **experimental** and are **not part of the v0.1
tested scope**. They are exported so that adventurous users can
start exercising them, but they are not yet validated through the
staged workflow:

- **Model-specification DSL** for time-to-event analyses
  (`specify_models()`, `identify_outcome()`, `identify_treatment()`,
  `identify_censoring()`, `identify_subject()`)
- **Time-to-event estimators** --- `estimate_ipwrisk()`,
  `estimate_gcomprisk()`, `estimate_aipwrisk()`, `estimate_ipwhr()`,
  `estimate_surv_tmle()` (via `survtmle`), `estimate_lmtp()` (via
  `lmtp`), `estimate_tmle_risk_point()`. Survival, competing-risk,
  and longitudinal estimands are on the cleanTMLE 0.2 roadmap.
- **FIORD two-stage candidate selector** (`select_tmle_candidate(rule = "fiord_two_stage")`)
  -- planned; current rules collapse the FIORD two stages into a
  single-step minimax. See *TODO.md* for the full development plan.
- **Nonparametric bootstrap variance** for IPTW, matched TMLE, and
  TMLE on non-i.i.d. samples -- planned for the next minor release.
- **Pre-protocol stress-test mode** (`run_preprotocol_plasmode()`)
  -- planned; allows the plasmode + DQ loop to run on user-specified
  covariate distributions without a real lock.
- **Unmeasured-confounding severity gradient** in `run_plasmode_dq_stress()`
  -- planned; currently a single fixed-strength scenario.

## Installation

```r
# Install from GitHub
# install.packages("remotes")
remotes::install_github("amertens/cleanTMLE")
```

## Worked Example

The following example walks through the dossier loop using the
built-in simulated dataset. For a complete narrative walkthrough, see
`vignette("cleanTMLE-staged-analysis")`; for a compact function
reference, see `vignette("cleanTMLE-functions")`.

> **Warning:** Low-replicate examples in this README (and in the bundled `config/clean_room_config.yml`) are for workflow demonstration only. They should not be interpreted as stable operating-characteristic estimates. Use `n_reps >= 200` for inferential interpretation.

```{r eval = FALSE}
library(cleanTMLE)

# 1. Lock the analytic specification (Stage 1a)
lock <- create_analysis_lock(
  data        = study_data,
  treatment   = "A",
  outcome     = "Y",
  covariates  = baseline_covariates,
  sl_library  = c("SL.glm", "SL.glmnet", "SL.ranger"),
  seed        = 2026
)
lock <- attach_estimand(lock,
  description           = "Effect of A on 24-month event risk",
  population            = "Adults eligible at index",
  treatment_strategies  = c("Treatment", "Reference"),
  outcome_label         = "Event by 24 months",
  followup              = "24 months",
  contrast              = "risk_difference"
)

# 2. Cohort flow / attrition
attr_tbl <- attrition_table(raw_data, inclusion_rules)

# 3. Design diagnostics (Stage 2a)
ps      <- fit_ps_superlearner(lock)
ps_diag <- compute_ps_diagnostics(ps)
wd      <- clean_weight_diagnostics(ps$weights,
                                    treatment  = lock$data[[lock$treatment]],
                                    covariates = lock$data[, lock$covariates])

# 4. Baseline plasmode candidate selection (Stage 2b)
plas <- run_plasmode_feasibility(lock, tmle_candidates = candidates, reps = 200)

# 5. DQ stress testing (Stage 2c)
dq   <- run_plasmode_dq_stress(lock, tmle_candidates = candidates,
                               scenarios = default_dq_scenarios(), n_reps = 200)

# 6. Decision log entries throughout
audit <- init_decision_log(lock)
audit <- log_decision_entry(audit, stage = "Stage 2c", decision = "GO",
                            rationale = "All candidates within thresholds",
                            reviewer = "lead_analyst")

# 7. Pre-outcome gate
gate  <- gate_all(checkpoint_cohort_adequacy(lock),
                  checkpoint_balance(ps_diag),
                  checkpoint_weights(ps$weights))

# 8. Authorisation record
auth  <- authorize_outcome_analysis(lock, gate)

# 9. Primary analysis (only callable after authorisation)
fit   <- run_clean_tmle(lock, ps_fit = ps)
```

## Function Reference

### Stage 1a: Analysis specification and lock

| Function | Purpose |
|----------|---------|
| `create_analysis_lock()` | Lock the analytic specification (data, treatment, outcome, covariates, SL library, plasmode settings) |
| `validate_analysis_lock()` | Verify lock integrity via hash check |
| `attach_estimand()` | Attach causal question, population, contrast, and follow-up metadata |
| `declare_sensitivity_plan()` | Pre-register a sensitivity analysis |
| `define_negative_control()` | Register a negative control variable |

### Stage 1b: Cohort adequacy and design precision (Check Point 1)

| Function | Purpose |
|----------|---------|
| `checkpoint_cohort_adequacy()` | Check sample size, events, treatment prevalence, positivity |
| `estimate_design_precision()` | Design-stage SE proxy, CI half-width, and MDD at 80% power |
| `summarize_event_support()` | Event counts and rates per treatment arm |

### Stage 2: Propensity score and balance (Check Point 2)

| Function | Purpose |
|----------|---------|
| `fit_ps_superlearner()` | SuperLearner ensemble PS estimation |
| `fit_ps_glm()` | Logistic regression PS estimation |
| `compute_ps_diagnostics()` | Overlap plots, ESS, standardised mean differences |
| `checkpoint_balance()` | GO / FLAG / STOP based on balance and ESS |

### Stage 2b: Baseline plasmode TMLE candidate selection

| Function | Purpose |
|----------|---------|
| `tmle_candidate()` | Define a TMLE candidate spec (library + truncation) |
| `expand_tmle_candidate_grid()` | Generate a default candidate grid |
| `validate_tmle_candidates()` | Validate a list of candidate specs |
| `run_plasmode_feasibility()` | Evaluate candidates on plasmode simulations |
| `select_tmle_candidate()` | Select best candidate via prespecified rule |
| `lock_primary_tmle_spec()` | Lock selected spec into the analysis lock |
| `get_primary_tmle_spec()` | Retrieve the locked spec |
| `gate_check()` | GO / FLAG / STOP from plasmode metrics |
| `summarize_plasmode_results()` | Print plasmode performance summary |

### Stage 2c: Data-quality stress testing

| Function | Purpose |
|----------|---------|
| `run_plasmode_dq_stress()` | Perturb synthetic outcomes under locked DQ severity ranges (missingness, misclassification, unmeasured confounding) |
| `summarize_dq_degradation()` | Per-candidate degradation gradients (bias, RMSE, coverage) versus the baseline plasmode |

### Stage 3: Residual confounding (Check Point 3)

| Function | Purpose |
|----------|---------|
| `run_negative_control()` | Estimate treatment effect on a negative control outcome |
| `run_residual_confounding_stage()` | Stage 3 wrapper: runs all registered NCs and checkpoints |
| `checkpoint_residual_bias()` | GO / FLAG / STOP based on negative control results |

### Pre-Outcome Authorization Gate

| Function | Purpose |
|----------|---------|
| `authorize_outcome_analysis()` | Scan audit for required checkpoints; return GO / FLAG / STOP |
| `assert_outcome_authorized()` | Error if outcome analysis is not authorized |

### Stage 4: Final estimation

| Function | Purpose |
|----------|---------|
| `run_crude_workflow()` | Unadjusted risk difference |
| `run_match_workflow()` | 1:1 nearest-neighbour PS matching |
| `run_iptw_workflow()` | Stabilised IPTW (Hajek estimator) |
| `fit_tmle_treatment_mechanism()` | TMLE g-step (uses locked truncation) |
| `fit_tmle_outcome_mechanism()` | TMLE Q-step (uses locked Q-library) |
| `run_tmle_targeting_step()` | TMLE fluctuation / targeting update |
| `extract_tmle_estimate()` | Final ATE, SE, CI, diagnostics |
| `fit_tmle_candidate_set()` | Fit multiple TMLE specs on real data |
| `fit_final_workflows()` | Run matching + IPTW + TMLE in one call |

### Sensitivity and diagnostics

| Function | Purpose |
|----------|---------|
| `sensitivity_truncation()` | Re-estimate under alternate truncation thresholds |
| `compute_evalue()` | E-value for unmeasured confounding |
| `new_checkpoint()` | Create a custom GO / FLAG / STOP checkpoint |

### Audit trail and decision log

| Function | Purpose |
|----------|---------|
| `create_audit_log()` | Initialise an audit log from a lock |
| `record_stage()` | Append a stage entry |
| `record_checkpoint()` | Append a checkpoint entry |
| `export_audit_trail()` | Export trail as a data.frame |
| `build_stage_manifest()` | Print a stage-path summary |
| `record_decision_log_entry()` | Record a structured analyst decision |
| `export_decision_log()` | Export decision log as a data.frame |
| `summarize_stage_path()` | Compact narrative of the analysis path |

### Cross-workflow summaries

| Function | Purpose |
|----------|---------|
| `summarize_cleanroom_results()` | Side-by-side table of all estimates |
| `compare_fits()` | Compare point estimates and CIs across fits |
| `forest_plot()` | Forest plot of treatment effect estimates |

### Model specification DSL

| Function | Purpose |
|----------|---------|
| `specify_models()` | Initialise a model specification object |
| `identify_outcome()` | Declare the outcome variable and model |
| `identify_treatment()` | Declare treatment and PS formula |
| `identify_censoring()` | Declare censoring variable and model |
| `identify_subject()` | Declare subject ID |
| `identify_competing_risk()` | Declare competing risk indicator |
| `identify_interval()` | Declare interval specification |
| `identify_missing()` | Declare missingness handling |

### Time-to-event estimation

| Function | Purpose |
|----------|---------|
| `estimate_ipwrisk()` | IPW risk curves (Kaplan-Meier reweighted) |
| `estimate_gcomprisk()` | G-computation risk estimates |
| `estimate_aipwrisk()` | Augmented IPW risk estimates |
| `estimate_ipwhr()` | IPW-weighted Cox hazard ratio |
| `estimate_surv_tmle()` | Survival TMLE via `survtmle` |
| `estimate_lmtp()` | Longitudinal TMLE via `lmtp` |
| `estimate_tmle_risk_point()` | Point-treatment TMLE risk estimate |
| `re_estimate()` | Re-estimate with modified specification |
| `update_outcome()`, `update_treatment()`, `update_censoring()` | Modify spec components |

### Reporting helpers

| Function | Purpose |
|----------|---------|
| `make_table1()` | Baseline characteristics table |
| `make_table2()` | Treatment effect summary table |
| `make_wt_summary_table()` | Weight distribution summary |
| `extreme_weights()` | Inspect extreme IPW weights |
| `inspect_ipw_weights()` | Weight diagnostics |
| `hr_data()` | Extract hazard ratio data |

### Utilities

| Function | Purpose |
|----------|---------|
| `expit()`, `logit()` | Inverse-logit and logit transforms |
| `sim_func1()` | Simulate example causal-inference dataset |
| `mask_outcome()` | Replace outcome column with NA (design-stage blinding) |
| `unmask_outcome()` | Restore outcome from original lock |

## Package Philosophy

The outcome-blind staged workflow records outcome blinding through
software-mediated stage gates, documenting analytic degrees of freedom
for review. Traditional diagnostics (overlap, balance, weight
distributions) are necessary but insufficient: a good-looking PS
overlap plot does not by itself indicate that a particular TMLE
specification will have low bias in the sample at hand. Conversely,
marginal overlap does not necessarily imply that a targeted estimator
will perform poorly.

Plasmode simulation bridges this gap. By evaluating the full estimator
pipeline on outcome-blind simulations derived from the real covariate
distribution, analysts obtain pre-outcome evidence about the relative
performance of competing TMLE specifications. `cleanTMLE` supports this
in a structured way:

- The candidate set and selection rule are locked before outcome access.
- Selection is rule-based and documented.
- Conventional PS methods serve as secondary comparators, not as the
  selection target.
- The audit trail from specification to final estimate is preserved.

The four-step TMLE design directly mirrors the outcome-blind staged
workflow: the treatment mechanism can be estimated in Stage 2 without
outcome access; the outcome mechanism is first estimated on real
data in Stage 4; the targeting step follows only after both nuisance
estimates are in hand. This separation makes the stage boundaries
explicit in the code itself.

## Notes on Nuisance Estimation

`cleanTMLE` uses SuperLearner as the default propensity-score estimation
strategy. In outcome-blind staged workflows, flexible treatment-model
estimation is often desirable: the covariate set may be high-dimensional,
relationships may be nonlinear, and it can be difficult to pre-specify
a correctly-specified parametric PS model before the outcome is
examined. SuperLearner provides a principled, data-adaptive approach
that can be fully pre-specified by locking the candidate library in
Stage 1a.

Simpler logistic regression PS models are also supported via
`fit_ps_glm()` and may be appropriate in low-dimensional settings. The
choice of nuisance strategy is declared in the analysis lock and cannot
be changed after outcome data are accessed.

## Relationship to Existing R Packages

`cleanTMLE` provides the workflow scaffolding (staged specification,
outcome-blind diagnostics, simulation-based candidate review, and
structured output) that underlying estimation packages do not
themselves provide. Estimation functionality may be delegated to:

- [`tmle`](https://cran.r-project.org/package=tmle) --- point-treatment TMLE
- [`SuperLearner`](https://cran.r-project.org/package=SuperLearner) ---
  ensemble learning for nuisance models
- [`survtmle`](https://github.com/benkeser/survtmle) --- survival TMLE
- [`lmtp`](https://cran.r-project.org/package=lmtp) --- longitudinal
  modified treatment policies
- [`glmnet`](https://cran.r-project.org/package=glmnet) --- regularised
  regression for nuisance estimation

## Development Status

`cleanTMLE` is under active development. The public API may still change
before a stable release. Issues and feature requests are welcome on the
[GitHub issue tracker](https://github.com/amertens/cleanTMLE/issues).

## License

MIT --- see [LICENSE.md](LICENSE.md) for details.
