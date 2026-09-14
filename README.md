# cleanTMLE

**Decide, before outcome access, whether a comparison is estimable, and
with which estimand. Then estimate it with TMLE.**

A clean-room design (Muntner et al. 2024) blinds the analyst: the design is
fixed before anyone sees the treatment-outcome association. What blinding
alone cannot do is tell the design team whether the design supports the
estimand they plan to estimate. On a cohort with a severe practical
positivity violation, the average treatment effect can come back with a
narrow confidence interval and be wrong by an order of magnitude, because
truncation limits variance but not extrapolation, and nothing in the
estimator's own output separates such an estimate from a real finding.
cleanTMLE's design stage answers the estimability question with no outcome
access, and its estimation stage estimates exactly what the design stage
declared.

## The workflow

The exported surface is sixteen verbs, run in stage order on one lock:

```r
lock <- create_analysis_lock(data, "A", "Y", covariates,
          negative_controls = c("nc_1", "nc_2"),
          dq_thresholds = list(max_abs_bias = 0.02, min_coverage = 0.90,
                               max_rmse_ratio = 1.5),
          nc_criteria   = list(null_band = 0.02))
lock <- declare_negative_controls(lock, c("nc_1", "nc_2"),
          domains = c("confounding_by_indication", "health_seeking_behavior"))
lock <- declare_estimand_ladder(lock, primary = "ATE",
          fallbacks = c("trimmed_ATE", "ATT", "ATO"))

ps  <- fit_ps(lock)                  # SuperLearner, GLM, or external scores
sup <- assess_support(ps)            # PASS / FLAG / SEVERE / FAIL, with the
                                     # caveat that travels with every estimate
fea <- estimand_feasibility(ps)      # which estimands these weights support
sim <- simulate_support(lock, ps)    # outcome-blind support map over a
                                     # prespecified outcome-surface family

cands <- define_candidates(grid = list(truncations = c(0.01, 0.05),
                                       libraries = list(glm = "SL.glm")))
st   <- stress_test(lock, cands, threats = "regulatory_standard")
best <- select_candidate(st, rule = "min_max_rmse")
lock <- declare_estimand_ladder(lock, primary = "ATE", candidate = best)

ncl <- negative_control_ladder(lock, restrictions = ...)   # Check Point 3
cpi <- check_process_indicators(lock, indicators = ...)    # collider screen

design_report(lock, sup, fea, simulation = sim,
              nc_ladder = ncl, dq = st)   # what the review team reads
export_design_log(lock, format = "muntner")   # the released decision record

fit <- run_estimand_ladder(lock, ps) # the declared primary when feasible,
                                     # plus every feasible fallback, each row
                                     # labelled with estimand and verdict
fit <- estimate_effect(lock, ps, estimand = "ATE", estimator = "tmle")
```

The support verdict grades overlap on the fitted propensity (share of the
sample outside a prespecified band, largest weight, with escalation when any
covariate stratum is near-deterministic in treatment). The feasibility table
reads the same fit per estimand: the ATT's control weights are bounded by
max g/(1-g) and the overlap-weighted ATO's weights are bounded by one, so
both can remain estimable where the ATE is not. The support simulation asks
the dynamic version of the question under the generate-treatment plasmode
design (the observed-treatment design induces a positivity violation by
construction; Shaw et al. 2025, arXiv:2504.11740), with every estimand's
truth computed from the generating model on every replicate. The ladder then
makes estimand switching a pre-registered, logged decision instead of a
silent substitution, and the implausibility guard flags any estimate the
observed data cannot support.

Outcome-blind candidate selection (`define_candidates()`,
`stress_test()`, `select_candidate()`) evaluates nuisance strategies on
simulated outcomes before real outcome access, and the same
`stress_test()` loop is the prespecified data-quality stress test:
covariate missingness in MCAR, MAR and MNAR forms, treatment and
outcome misclassification, unmeasured confounding, and near-positivity,
graded against thresholds declared at lock creation into a locked
GO / FLAG / STOP verdict with tipping points. The stress test is a
quantitative pre-outcome supplement to fit-for-purpose data review, not
a formal QBA and not a substitute for source-data validation.

## What cleanTMLE adds over a traditional clean room

A traditional clean room provides blinding and an audit trail. cleanTMLE
adds the two things a review team cannot get from blinding: a defensible,
prespecified answer to "is this comparison estimable, and with which
estimand", produced before outcomes are unlocked, and a doubly robust
estimator for whatever the answer turns out to be. The published analogue is
the blinded validity-diagnostics gate of Conover et al. (2025, JAMIA):
named diagnostics computed while estimates stay blinded, unblinding only
when they pass, and failed comparisons labelled inestimable.

## What cleanTMLE does not do

- does not automate target-trial design or protocol specification
- does not validate source data
- does not validate phenotypes or outcome definitions
- does not establish exchangeability
- does not establish positivity
- does not guarantee that simulated-outcome rankings generalise
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
- Independent review of deviations and overrides recorded in the design log

## When to use this package

- methods development around outcome-blind workflows
- outcome-blind estimator selection in observational studies
- protocol and SAP development for RWE analyses
- exploratory or lower-stakes RWE workflows where structured
  pre-outcome review is useful
- high-stakes confirmatory RWE only when embedded inside
  external clean-room governance (role separation, data-access
  controls, independent review, archived decision record)

## Three ways to run cleanTMLE (tiers of enforcement)

Pick the tier by how much the study needs the outcome-access guarantee
actually enforced, not just documented.

**Tier 1 - exploratory, no guarantees.** `run_clean_tmle()` (the
unguarded convenience wrapper) or a lock created with
`cleanroom_enabled = FALSE` puts crude / IPTW / matching / TMLE
estimates onto one forest plot with nothing enforced; outcome blindness
is the analyst's responsibility. Good for methods work and quick looks.

**Tier 2 - the default lock.** An ordinary `create_analysis_lock()`
physically separates the outcome at creation: `lock$data` holds design
data only, and the primary outcome lives in a sealed store that only
the Stage 4 estimators join back. `mask_outcome()` removes the store
entirely, so a masked lock is outcome-free wherever it travels, and
every design verb returns identical output on masked and unmasked
locks (a tested invariant). The per-call escape hatch
(`allow_outcome_access = TRUE`) exists and is visible in code review.

**Tier 3 - enforced authorisation.** A lock created with
`enforce = TRUE` refuses Stage 4 estimation until the outcome has been
unmasked through `unmask_outcome(lock, original_lock, approved_by =
<name>)`. The named approval is written into the design log
(`export_design_log(format = "muntner")` releases it in the decision-log
structure of Muntner et al. 2024, Table S1). There is no token and no
gate object: the single institutional switch is the enforce flag, and
the single act of authorisation is a person on the record. cleanTMLE
cannot enforce that the approver and the analyst are different people;
that separation is external governance.

## Cumulative-risk workflow considerations

cleanTMLE provides software support for the analytic considerations that arise in cumulative-risk pharmacoepidemiology workflows:

- **Model specification grammar**: `identify_*()` family and `specify_models()` for declaring eligibility, treatment, outcome, censoring, competing risk, follow-up interval, and intercurrent events as code objects.
- **Cumulative-risk reporting**: `clean_risk_report_table()` produces a compact risk table at clinically meaningful time points; the package emphasises risk differences and risk ratios over hazard ratios except where the estimand is explicitly defined on the hazard scale.
- **Censoring and missingness weights as first-class objects**: `estimate_effect(missing = "ipcw")` and `clean_weight_diagnostics()` expose ESS, percentiles, maximum weight, and prespecified instability thresholds.
- **Event-process classification**: `clean_event_process_table()` and `clean_check_event_processes()` distinguish event of interest, competing event, censoring, treatment discontinuation/switching, transfer exclusions, and administrative end of follow-up.
- **Hazard-ratio de-emphasis**: cumulative risks at clinically meaningful follow-up times are reported as primary; hazard ratios are reserved for analyses where a hazard-scale estimand is the primary scientific question and proportional hazards is plausible.

## Overview

`cleanTMLE` is the software layer for an outcome-blind staged
workflow around targeted minimum loss-based estimation (TMLE).
The estimand, nuisance-model specifications, decision thresholds,
and estimator-selection rules are recorded in an analysis lock before
the observed primary treatment-outcome association is read by the
package. The package sits within, but does not replace, the
broader clean-room governance construct of Muntner et al. (2024),
which also covers role separation, restricted data access, and
independent review.

The design stage records graded verdicts on its result objects (the
support verdict, the per-estimand feasibility verdicts, the locked
data-quality verdict with its tipping points, and the
negative-control reading under the locked criteria), and the lock
accumulates a design log of declarations, switches, overrides, and
the unmasking approval, exportable for review.

The package covers three workflow families:

- **Conventional propensity-score methods** (matching, IPTW) as
  secondary comparators.
- **Fixed-specification TMLE** with a prespecified nuisance strategy.
- **Simulation-selected TMLE** in which candidate TMLE specifications
  (varying truncation and/or learner library) are evaluated on
  outcome-blind plasmode simulations, and a prespecified rule selects
  the best specification before the real outcome is accessed.

The **tested scope** is intentionally narrow: binary point
exposure, binary outcome, marginal risk difference, outcome
missingness handled through complete-case and IPCW sensitivity
paths under prespecified missingness assumptions. The package
also ships a **model-specification DSL** and a set of
**time-to-event helpers** (`estimate_ipwrisk`, `estimate_gcomprisk`,
`estimate_aipwrisk`, `estimate_ipwhr`, `estimate_surv_tmle`,
`estimate_lmtp`), but these are **experimental** and are not part
of the tested scope; see *Experimental / planned extensions*
below.

## Key features

- **One lock, sixteen verbs**: `create_analysis_lock()` seals data,
  roles, declarations, and decision thresholds behind a SHA-256
  fingerprint whose content digest covers every column except the
  outcome, so design-data tampering is detectable while masking and
  unmasking leave the hash unchanged.
- **The outcome store**: the lock's data frame holds design data
  only; the primary outcome lives in `lock$outcome_store` and is
  joined back only by Stage 4 estimators, after the outcome guard.
- **Declared decision rules**: `dq_thresholds` and `nc_criteria` are
  fingerprinted lock fields; the stress-test verdict and the
  negative-control reading are pure functions of the declared rules
  and the declared threat grid.
- **Estimand ladder**: `declare_estimand_ladder()` pre-registers the
  primary and its ordered fallbacks with the verdict that triggers a
  switch; `run_estimand_ladder()` executes it, labels every row, and
  logs every switch or override.
- **Outcome-blind candidate selection and DQ stress testing**:
  `define_candidates()`, `stress_test()` (baseline and declared
  threats in one loop, with `max_fit_seconds` bounding runaway fits
  and `parallel = TRUE` for a furrr path with identical results),
  `select_candidate()` (rules `min_rmse`, `min_max_rmse`,
  `fiord_two_stage`).
- **Residual-confounding and collider screens**:
  `negative_control_ladder()` across nested cohort restrictions with
  domains and failed fits kept visible; `check_process_indicators()`
  for care-process colliders.
- **Design-stage precision without unblinding**:
  `estimate_design_precision()` is marginal-only; arm-specific event
  counts require `event_support_by_arm(lock, reason = )`, which
  warns and writes the access into the design log.
- **Reports and records**: `design_report()` (summary statistics
  only, with the recommendation), `export_design_log()` (tidy or
  Muntner Table S1 structure), snapshot-tested.
- **Estimation through one front door**: `estimate_effect(estimand,
  estimator, missing)` covers TMLE, IPTW, matching, and crude, with
  `return_steps = TRUE` exposing the modular TMLE quartet;
  `select_variance_method()` and `bootstrap_rd_variance()` supply
  the variance method, including the matching-aware bootstrap.
- **Governance notes**: `clean_event_process_table()`,
  `clean_check_event_processes()`, `clean_target_population()`,
  `clean_missing_data_plan()`, `clean_risk_report_table()` record
  the specification a clean-room review requires.
- **Templates**: a targeted-learning SAP template and a Muntner
  Table S1 decision-log CSV under `inst/templates/`.

## Experimental / planned extensions

The following are **experimental** and are **not part of the
tested scope**. They are exported so that adventurous users can
start exercising them, but they are not yet validated through the
staged workflow:

- **Model-specification DSL** for time-to-event analyses
  (`specify_models()`, `identify_outcome()`, `identify_treatment()`,
  `identify_censoring()`, `identify_subject()`)
- **Time-to-event estimators**: `estimate_ipwrisk()`,
  `estimate_gcomprisk()`, `estimate_aipwrisk()`, `estimate_ipwhr()`,
  `estimate_surv_tmle()` (via `survtmle`), `estimate_lmtp()` (via
  `lmtp`). Survival, competing-risk, and longitudinal estimands are
  on the roadmap.
- **Pre-protocol stress-test mode**: planned. Allows the plasmode
  and DQ loop to run on user-specified covariate distributions
  without a real lock; `dgp_mode = "external_pilot"` covers the
  external-pilot half of this today.

## Installation

```r
# Install from GitHub
# install.packages("remotes")
remotes::install_github("amertens/cleanTMLE")
```

## Worked example

For the compact tour, see `vignette("cleanTMLE")` (five verbs, two
minutes); for the complete narrative, the *Full workflow* article
(`vignette("cleanTMLE-staged-analysis")`).

> **Warning:** Low-replicate examples in this README are for workflow
> demonstration only. Use `reps >= 200` before interpreting bias,
> coverage, or RMSE inferentially.

```r
library(cleanTMLE)
dat <- sim_func1(n = 1000, seed = 42)

lock <- create_analysis_lock(
  dat, "treatment", "event_24",
  covariates        = c("age", "sex", "biomarker", "comorbidity"),
  negative_controls = "nc_outcome",
  dq_thresholds     = list(max_abs_bias = 0.02, min_coverage = 0.90,
                           max_rmse_ratio = 1.5),
  nc_criteria       = list(null_band = 0.02),
  seed              = 42, enforce = TRUE)
lock <- declare_estimand_ladder(lock, primary = "ATE",
                                fallbacks = c("trimmed_ATE", "ATT", "ATO"))

ps  <- fit_ps(lock, method = "glm")
sup <- assess_support(ps)
fea <- estimand_feasibility(ps)

cands <- define_candidates(grid = list(truncations = c(0.01, 0.05),
                                       libraries = list(glm = "SL.glm")))
st   <- stress_test(lock, cands, threats = "regulatory_standard",
                    reps = 20)              # >= 200 in a real study
best <- select_candidate(st, rule = "min_max_rmse")
lock <- declare_estimand_ladder(lock, primary = "ATE", candidate = best)

design_report(lock, sup, fea, dq = st)      # the review team reads this

# Stage 4 on an enforce lock runs only after a named, logged unmasking.
masked <- mask_outcome(lock)
lock   <- unmask_outcome(masked, lock, approved_by = "review team")
fit    <- run_estimand_ladder(lock, ps)
```

## Function reference

The pkgdown reference index groups the API by stage: the sixteen
workflow verbs, the estimator layer (with the internal modular TMLE
steps), the learners, the lock utilities, and the reporting and
governance helpers. The sixteen verbs, in stage order:

| Stage | Verbs |
|---|---|
| 1 - lock and declarations | `create_analysis_lock()`, `declare_negative_controls()`, `declare_estimand_ladder()` |
| 2 - design estimability | `fit_ps()`, `assess_support()`, `estimand_feasibility()`, `simulate_support()` |
| 2b - candidates and stress test | `define_candidates()`, `stress_test()`, `select_candidate()` |
| 3 - residual confounding | `negative_control_ladder()`, `check_process_indicators()` |
| review | `design_report()`, `export_design_log()` |
| 4 - estimation | `estimate_effect()`, `run_estimand_ladder()` |

## Package philosophy

The outcome-blind staged workflow records outcome blinding through
software-mediated stages, documenting analytic degrees of freedom
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
- The design log from specification to final estimate is preserved.

The modular TMLE design directly mirrors the outcome-blind staged
workflow: the treatment mechanism can be estimated in Stage 2 without
outcome access; the outcome mechanism is first estimated on real
data in Stage 4; the targeting step follows only after both nuisance
estimates are in hand. This separation makes the stage boundaries
explicit in the code itself.

## Notes on nuisance estimation

`cleanTMLE` uses SuperLearner as the default propensity-score estimation
strategy. In outcome-blind staged workflows, flexible treatment-model
estimation is often desirable: the covariate set may be high-dimensional,
relationships may be nonlinear, and it can be difficult to pre-specify
a correctly-specified parametric PS model before the outcome is
examined. SuperLearner provides a principled, data-adaptive approach
that can be fully pre-specified by locking the candidate library in
Stage 1a.

Simpler logistic regression PS models are also supported via
`fit_ps(method = "glm")` and may be appropriate in low-dimensional
settings. The choice of nuisance strategy is declared in the analysis
lock and cannot be changed after outcome data are accessed.

## Relationship to existing R packages

`cleanTMLE` provides the workflow scaffolding (staged specification,
outcome-blind diagnostics, simulation-based candidate review, and
structured output) that underlying estimation packages do not
themselves provide. Estimation functionality may be delegated to:

- [`tmle`](https://cran.r-project.org/package=tmle): point-treatment TMLE
- [`SuperLearner`](https://cran.r-project.org/package=SuperLearner):
  ensemble learning for nuisance models
- [`survtmle`](https://github.com/benkeser/survtmle): survival TMLE
- [`lmtp`](https://cran.r-project.org/package=lmtp): longitudinal
  modified treatment policies
- [`glmnet`](https://cran.r-project.org/package=glmnet): regularised
  regression for nuisance estimation

## Development status

`cleanTMLE` is under active development. The public API may still change
before a stable release. Issues and feature requests are welcome on the
[GitHub issue tracker](https://github.com/amertens/cleanTMLE/issues).

## License

MIT: see [LICENSE.md](LICENSE.md) for details.