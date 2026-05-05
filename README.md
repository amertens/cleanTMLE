# cleanTMLE

**Staged Clean-Room Causal Analysis with TMLE**

## Attribution

cleanTMLE reuses the cohort-specification API and several
estimator/reporting functions from the
[**causalRisk**](https://github.com/CausalInference/causalRisk)
package by M. Alan Brookhart, Alexander Breskin, and Bassim
Eledath, under the MIT licence and with attribution. Reused
elements include: the model-specification DSL (`specify_models`,
`identify_treatment`, `identify_outcome`, `identify_censoring`,
`identify_competing_risk`, `identify_subject`, `identify_interval`,
`identify_missing`); the IPW / G-computation / AIPW estimators
(`estimate_ipwrisk`, `estimate_gcomprisk`, `estimate_aipwrisk`,
`estimate_ipwhr`); and the reporting helpers (`make_table1`,
`make_table2`, `extreme_weights`, `inspect_ipw_weights`,
`forest_plot`, `compare_fits`, `re_estimate`, `hr_data`).
Functions unique to cleanTMLE (the analysis lock, outcome masking,
plasmode feasibility and DQ stress, audit and decision logs, the
modular TMLE primitives, and the GO/FLAG/STOP gate) are
documented as cleanTMLE-original in the corresponding `.Rd` files.
We thank the causalRisk authors for the MIT-licensed code that
made this attribution-based reuse possible.

## Overview

`cleanTMLE` supports a structured, outcome-blinded workflow for
observational causal analysis based on targeted minimum loss-based
estimation (TMLE). The package implements the scaffolding needed to plan
and execute a staged clean-room analysis in which analytic
decisions --- including the estimand, nuisance-model specifications, and
estimator selection rules --- are locked before any investigator accesses
the study outcome data.

The workflow follows the Muntner et al. (2020) staged design:

```
Stage 1a  Specify the estimand, lock the analytic plan
Stage 1b  Assess cohort adequacy + design precision      (Check Point 1)
Stage 2   Estimate propensity scores, assess overlap     (Check Point 2)
Stage 2b  Plasmode simulation: select TMLE specification (pre-outcome)
Stage 3   Residual confounding via negative controls     (Check Point 3)
  Gate    Pre-outcome authorization (GO / STOP)
Stage 4   Conduct the comparative analysis (outcome unblinded)
```

Each checkpoint produces a structured **GO / FLAG / STOP** decision.
The analysis may proceed only after all preceding checkpoints have been
evaluated. An audit trail accumulates entries automatically and can be
exported for review.

The package covers three workflow families:

- **Conventional propensity-score methods** (matching, IPTW) as
  secondary comparators.
- **Fixed-specification TMLE** with a prespecified nuisance strategy.
- **Simulation-selected TMLE** in which candidate TMLE specifications
  (varying truncation and/or learner library) are evaluated on
  outcome-blind plasmode simulations, and a prespecified rule selects
  the best specification before the real outcome is accessed.

In addition, the package provides a **model-specification DSL** and
**time-to-event estimation** functions (IPW risk curves, g-computation,
augmented IPW, IPW hazard ratios, survival TMLE, and LMTP) for
richer applied analyses.

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
  `override_clean_room = TRUE` is explicitly set
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
  (`select_tmle_candidate()`), lock the chosen spec
  (`lock_primary_tmle_spec()`)
- **Gate decision** --- structured GO / FLAG / STOP based on bias,
  coverage, and SE calibration from plasmode results (`gate_check()`)
- **Modular TMLE** --- intentionally separated so each step runs at the
  correct clean-room stage:
  1. Treatment mechanism / g-step (`fit_tmle_treatment_mechanism()`)
  2. Outcome mechanism / Q-step (`fit_tmle_outcome_mechanism()`)
  3. Targeting / fluctuation step (`run_tmle_targeting_step()`)
  4. Final estimate extraction (`extract_tmle_estimate()`)
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
- **Model specification DSL** --- pipe-friendly interface for
  time-to-event analyses (`specify_models()`, `identify_outcome()`,
  `identify_treatment()`, `identify_censoring()`, `identify_subject()`)
- **Time-to-event estimators** --- `estimate_ipwrisk()`,
  `estimate_gcomprisk()`, `estimate_aipwrisk()`, `estimate_ipwhr()`,
  `estimate_surv_tmle()`, `estimate_lmtp()`,
  `estimate_tmle_risk_point()`
- **Reporting helpers** --- `make_table1()`, `make_table2()`,
  `make_wt_summary_table()`, `extreme_weights()`, `compare_fits()`,
  `forest_plot()`

## Installation

```r
# Install from GitHub
# install.packages("remotes")
remotes::install_github("amertens/cleanTMLE")
```

## Worked Example

The following example demonstrates the full staged workflow using the
built-in simulated dataset. For a complete narrative walkthrough, see
`vignette("cleanTMLE-staged-analysis")`; for a compact function
reference, see `vignette("cleanTMLE-functions")`.

```r
library(cleanTMLE)

# ── Stage 1a: Specify and lock the analysis ───────────────────────────────

dat <- sim_func1(n = 2000, seed = 42)

lock <- create_analysis_lock(
  data          = dat,
  treatment     = "treatment",
  outcome       = "event_24",
  covariates    = c("age", "sex", "biomarker", "comorbidity"),
  sl_library    = c("SL.glm", "SL.mean"),
  plasmode_reps = 200L,
  seed          = 42L
)

# Attach estimand, sensitivity plan, and negative control declaration
lock <- attach_estimand(lock,
  description          = "Effect of treatment on 24-month event risk",
  population           = "Adults with simulated disease",
  treatment_strategies = c("Treatment", "Control"),
  outcome_label        = "Primary event by 24 months",
  followup             = "24 months",
  contrast             = "risk_difference",
  statistical_estimand = "E_W[E(Y|A=1,W) - E(Y|A=0,W)]"
)

lock <- declare_sensitivity_plan(lock, "truncation",
  description = "Vary PS truncation",
  settings    = list(thresholds = c(0.01, 0.05, 0.10))
)

lock <- define_negative_control(lock, "nc_outcome",
  description = "Outcome driven by covariates only, no treatment effect"
)

validate_analysis_lock(lock)

# Initialise audit trail
audit <- create_audit_log(lock)
audit <- record_stage(audit, "Stage 1a", "Lock created, estimand attached")

# ── Stage 1b / Check Point 1: Cohort adequacy ────────────────────────────

cp1 <- checkpoint_cohort_adequacy(lock, min_n_per_arm = 50, min_events = 30)
print(cp1)                          # GO / FLAG / STOP
audit <- record_checkpoint(audit, cp1)

# Design-stage precision diagnostics
dp <- estimate_design_precision(lock, target_mdd = 0.05)
print(dp)
summarize_event_support(lock)

# ── Stage 2 / Check Point 2: PS overlap and balance ──────────────────────

ps_fit <- fit_ps_superlearner(lock)  # or fit_ps_glm(lock) for speed
diag   <- compute_ps_diagnostics(ps_fit)
print(diag)                          # ESS, SMDs, overlap plot
plot(diag)                           # propensity score overlap

cp2 <- checkpoint_balance(diag, max_smd = 0.10, min_ess_pct = 50,
                          lock_hash = lock$lock_hash)
print(cp2)
audit <- record_checkpoint(audit, cp2)

# ── Stage 2b: Plasmode TMLE candidate selection ──────────────────────────

candidates <- list(
  tmle_candidate("glm_t01", "GLM, trunc=0.01",
                 g_library = c("SL.glm"), truncation = 0.01),
  tmle_candidate("glm_t05", "GLM, trunc=0.05",
                 g_library = c("SL.glm"), truncation = 0.05),
  tmle_candidate("glm_t10", "GLM, trunc=0.10",
                 g_library = c("SL.glm"), truncation = 0.10)
)

plas <- run_plasmode_feasibility(lock, tmle_candidates = candidates,
                                 effect_sizes = c(0.05, 0.10))
print(plas)                          # bias, RMSE, coverage per candidate

selected <- select_tmle_candidate(plas, rule = "min_rmse")
lock <- lock_primary_tmle_spec(lock, selected)  # lock into the analysis
print(selected)

# Gate check: GO / FLAG / STOP
gate <- gate_check(plas$metrics, "Feasibility",
  targets = list(max_abs_bias = 0.02, min_coverage = 0.90,
                 se_sd_low = 0.8, se_sd_high = 1.2),
  method = selected$candidate_id
)
cat("Gate decision:", gate$decision, "\n")

# ── Stage 3 / Check Point 3: Residual confounding ────────────────────────

stage3 <- run_residual_confounding_stage(lock, ps_fit)
print(stage3)
cp3 <- stage3$checkpoint
audit <- record_checkpoint(audit, cp3)

# ── Pre-Outcome Authorization Gate ───────────────────────────────────────

gate_result <- authorize_outcome_analysis(audit)
print(gate_result)                  # GO / FLAG / STOP
audit <- record_checkpoint(audit, gate_result)

# ── Stage 4: Final estimation (outcome unblinded) ────────────────────────

# Crude benchmark
crude <- run_crude_workflow(lock)

# PS matching and IPTW (secondary comparators)
match_fit <- run_match_workflow(lock, ps_fit)
iptw_fit  <- run_iptw_workflow(lock, ps_fit)

# Primary: modular TMLE using locked specification
g_fit    <- fit_tmle_treatment_mechanism(lock, ps_fit)   # uses locked truncation
Q_fit    <- fit_tmle_outcome_mechanism(lock, g_fit)      # uses locked Q-library
tmle_upd <- run_tmle_targeting_step(g_fit, Q_fit)
tmle_fit <- extract_tmle_estimate(tmle_upd)
print(tmle_fit)

# Or use the convenience wrapper:
# all_fits <- fit_final_workflows(lock, ps_fit)

# Cross-estimator comparison
summary_tbl <- summarize_cleanroom_results(
  list(Matching = match_fit, IPTW = iptw_fit, TMLE = tmle_fit)
)
print(summary_tbl)

# ── Sensitivity analysis ─────────────────────────────────────────────────

sensitivity_truncation(lock, thresholds = c(0.01, 0.05, 0.10))
compute_evalue(1.3, ci_bound = 1.05)

# ── Audit trail ──────────────────────────────────────────────────────────

audit <- record_stage(audit, "Stage 4", "Final estimation complete")
print(audit)
build_stage_manifest(audit)
summarize_stage_path(audit)         # compact narrative
export_audit_trail(audit)           # returns a data.frame
export_decision_log(audit)          # returns decision log data.frame
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

### Stage 2b: TMLE candidate selection via plasmode

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

The clean-room workflow enforces outcome blinding through
software-mediated stage gates, reducing analytic degrees of freedom.
Traditional diagnostics (overlap, balance, weight distributions) are
necessary but insufficient: a good-looking PS overlap plot does not
guarantee that a particular TMLE specification will have low bias in
the sample at hand. Conversely, marginal overlap does not necessarily
imply that a targeted estimator will perform poorly.

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

The modular TMLE API directly mirrors the clean-room design: the
treatment mechanism can be estimated in Stage 2 without outcome access;
the outcome mechanism is first estimated on real data in Stage 4; the
targeting step follows only after both nuisance estimates are in hand.
This separation makes the stage boundaries explicit in the code itself.

## Notes on Nuisance Estimation

`cleanTMLE` uses SuperLearner as the default propensity-score estimation
strategy. In clean-room workflows, flexible treatment-model estimation
is often desirable: the covariate set may be high-dimensional,
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

`cleanTMLE` provides the workflow scaffolding --- staged specification,
blinded diagnostics, simulation-based validation, and structured
output --- that underlying estimation packages do not themselves
enforce. Estimation functionality may be delegated to:

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
