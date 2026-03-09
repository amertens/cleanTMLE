# cleanTMLE

**Targeted Maximum Likelihood Estimation within the Staged Clean-Room Causal Analysis**

## Overview

`cleanTMLE` supports a structured, outcome-blinded workflow for observational
causal analysis. The package implements the scaffolding needed to plan and
execute a staged clean-room analysis in which analytic decisions — including
estimator choice, nuisance-model specifications, and selection rules — are
locked before any investigator accesses the study outcome data.

The package distinguishes two pre-outcome feasibility paths that can be
pursued independently:

- **Stage 2a** provides traditional feasibility diagnostics based on
  propensity-score overlap, covariate balance, and weight behavior — using
  only the treatment and covariate data.
- **Stage 2b** provides plasmode-simulation feasibility evaluation of the
  full estimator pipeline — running the complete TMLE workflow on synthetic
  outcomes so that estimator performance can be assessed before the real
  outcome is accessed.

These two paths are parallel, not sequential. Stage 2b does not require
completing Stage 2a first, and Stage 2a alone does not validate the full
estimator pipeline.

The package covers three workflow families: conventional propensity-score
methods, a fixed-specification TMLE workflow, and a simulation-selected TMLE
workflow in which a prespecified selection rule is used to choose among
candidate specifications based on blinded plasmode simulation performance.
Treatment-model work, outcome-model work, and the TMLE targeting step are
deliberately modular so they can be executed at the appropriate clean-room
stage. In addition to final estimation, the package provides tools for
overlap diagnostics, matching, inverse probability of treatment weighting
(IPTW), plasmode simulation, and workflow-level summaries, making the
clean-room process reproducible and auditable from specification to final
report.

## Key Features

- **Analysis lock and staged workflow support** — record and validate analytic
  decisions before outcome access
- **SuperLearner-based propensity-score estimation and diagnostics** —
  SuperLearner ensemble learning is the default and recommended PS estimation
  strategy; simpler logistic-regression models are also available for
  conventional workflows or simpler settings; includes overlap checks,
  effective sample size, and standardized mean differences
- **Matching and IPTW helpers** — 1:1 nearest-neighbor matching and stabilized
  IPTW with weight-distribution summaries
- **Plasmode-simulation feasibility evaluation** — outcome-blind simulation
  from the observed covariate and treatment distributions to evaluate the full
  estimator pipeline before real-outcome access
- **Modular TMLE components** — intentionally separated so each step can run
  at the correct clean-room stage:
  - fit treatment mechanism (`fit_tmle_treatment_mechanism()`)
  - fit outcome mechanism (`fit_tmle_outcome_mechanism()`)
  - perform targeting step (`run_tmle_targeting_step()`)
  - extract final estimate and inference (`extract_tmle_estimate()`)
- **Candidate TMLE specification comparison / selection** — fit a set of
  candidate TMLE specifications and apply a prespecified rule to select among
  them before real-outcome access
- **Workflow-level summaries and decision support** — side-by-side comparison
  of conventional and TMLE workflow results with structured output objects

## Workflow Overview

The package organizes analysis into three stages, with Stage 2 divided into
two parallel substages. Stages 1 and 2 occur before any investigator accesses
the study outcome; Stage 3 is the unblinded estimation step.

```
Stage 1: Analysis lock (pre-specification)
   ├─ Stage 2a: Traditional feasibility diagnostics
   │     (PS overlap, covariate balance, weight distributions)
   └─ Stage 2b: Plasmode-simulation feasibility evaluation
         (full estimator performance under synthetic outcomes)
               ↓
         Stage 3: Unblinded final estimation
```

**Stage 1 — Pre-specification and analysis lock.**
All analytic decisions are written down and validated before outcome data are
examined: estimand, covariate set, SuperLearner libraries and nuisance
strategies, plasmode simulation parameters, and (for the simulation-selected
workflow) the candidate set and selection rule. `create_analysis_lock()` and
`validate_analysis_lock()` serialize and verify this specification.

**Stage 2a — Traditional feasibility diagnostics.**
Using only treatment and covariate data, the treatment mechanism is fitted
(e.g., with `fit_ps_superlearner()`) and analysts assess whether the PS model
produces adequate overlap (effective sample sizes, trimming needs) and whether
covariate balance is achievable. No real outcome model is fitted at this
stage. These diagnostics may prompt revisions to the specification while the
outcome remains blinded.

**Stage 2b — Plasmode-simulation feasibility evaluation.**
Plasmode outcomes are generated from the observed covariate and treatment
distributions; the true causal effect is known by construction. The full
estimator pipeline — treatment mechanism, outcome mechanism, TMLE targeting,
and variance estimation — is run on these synthetic datasets. Analysts assess
bias, RMSE, and confidence-interval coverage and decide whether to proceed
and, if using the simulation-selected workflow, which TMLE specification to
carry forward. This step complements rather than follows Stage 2a, and it is
what allows pre-outcome validation of the full estimator.

**Stage 3 — Unblinded estimation.**
After the analysis plan is locked and feasibility evaluation is complete, the
real outcome mechanism is fitted for the first time, the targeting step is
performed, and final estimates are extracted. The workflow family (conventional
PS, fixed TMLE, or simulation-selected TMLE) and any simulation-driven
selections were committed before this step.

## Supported Workflow Families

### Conventional propensity-score workflow

Fits a propensity-score model — using logistic regression or a simpler
regression strategy as appropriate — checks overlap, and produces estimates
via 1:1 nearest-neighbor matching and stabilized IPTW. Covariate balance
tables and weight diagnostics are generated for reporting. This workflow is
primarily validated by Stage 2a diagnostics.

### Fixed-specification TMLE workflow

The treatment mechanism is estimated with a locked nuisance strategy (e.g.,
`SuperLearner` with a prespecified library) and its diagnostics are assessed
in Stage 2a. The outcome mechanism is estimated for the first time in Stage 3
on the real data, followed immediately by the TMLE targeting step and final
estimation. The full pipeline — treatment mechanism, outcome mechanism,
targeting, and inference — is validated end-to-end in Stage 2b before the
real outcome is accessed.

### Simulation-selected TMLE workflow

A candidate set of TMLE nuisance strategies (e.g., varying the
outcome-model learner or the propensity-score SuperLearner library) is
evaluated on plasmode simulations in Stage 2b. A prespecified selection rule
(e.g., lowest root-mean-squared error across simulated datasets) identifies
the best candidate. This is not post-hoc tuning: the candidate set and the
selection rule are both locked in Stage 1 before any data are examined. The
selected specification is then used in Stage 3 for final estimation.

## Package Philosophy

Stage 2a asks whether the treatment model and covariate overlap look
acceptable. It is a necessary check, but it does not assess how the full
estimator pipeline will perform in the sample at hand. Stage 2b asks whether
the estimator, end to end, performs adequately under plausible data-generating
processes. Both perspectives are useful and neither is a prerequisite for the
other. `cleanTMLE` is designed to support both paths and to make the transition
from feasibility assessment to final estimation structured and auditable rather
than ad hoc.

In particular, the modular TMLE API directly mirrors the clean-room design
itself: the treatment mechanism can be estimated and inspected in Stage 2a or
Stage 2b using only W and A; the outcome mechanism is estimated on real data
for the first time in Stage 3; and the targeting step follows only after both
nuisance estimates are available. This separation is not a limitation — it is
the point. It means the audit trail from specification to final estimate is
preserved at each stage, and it is clear which analytic decisions were made
before and after the outcome was accessed.

## Installation

```r
# Install from GitHub using remotes
if (!requireNamespace("remotes", quietly = TRUE)) {
  install.packages("remotes")
}
remotes::install_github("amertens/cleanTMLE")

# Or using pak
# install.packages("pak")
# pak::pak("amertens/cleanTMLE")
```

## Minimal Example

The following sketch illustrates the overall workflow. Function names marked
with comments are part of the staged-workflow layer; see the vignette for
a fully worked example.

```r
library(cleanTMLE)

# ── Stage 1: Pre-specification and analysis lock ───────────────────────────

# Simulate or load analysis data (covariates + treatment, outcome blinded)
dat <- sim_func1(n = 2000, seed = 42)  # sim_func1() is provided by the package

# Define the analytic specification and lock it
lock <- create_analysis_lock(
  data          = dat,
  treatment     = "treatment",
  outcome       = "event_24",
  covariates    = c("age", "sex", "biomarker"),
  sl_library    = c("SL.glm", "SL.glmnet", "SL.ranger"),  # SuperLearner library
  plasmode_reps = 200,
  seed          = 42
)

validate_analysis_lock(lock)      # checks completeness and reproducibility

# ── Stage 2a: Traditional feasibility diagnostics ─────────────────────────
# Fit treatment mechanism using only W and A (outcome never accessed here)

ps_fit  <- fit_ps_superlearner(lock)   # default SuperLearner-based PS
diag    <- compute_ps_diagnostics(ps_fit)
print(diag)                             # overlap plots, ESS, SMDs

# ── Stage 2b: Plasmode-simulation feasibility evaluation ──────────────────
# Generate synthetic outcomes; run full pipeline; assess performance

sim_results <- run_plasmode_feasibility(
  lock         = lock,
  effect_sizes = c(0.05, 0.10),  # plausible true risk differences
  reps         = lock$plasmode_reps
)
print(sim_results)                 # bias, RMSE, coverage by workflow

# Select best TMLE specification using the prespecified rule (locked in Stage 1)
best_spec <- select_tmle_candidate(sim_results, rule = "min_rmse")

# ── Stage 3: Unblinded final estimation ───────────────────────────────────

# Conventional PS workflow
match_fit <- run_match_workflow(lock, ps_fit)
iptw_fit  <- run_iptw_workflow(lock, ps_fit)

# Modular TMLE workflow (fixed or simulation-selected specification)
g_fit     <- fit_tmle_treatment_mechanism(lock, ps_fit)  # treatment mechanism
Q_fit     <- fit_tmle_outcome_mechanism(lock, g_fit)     # outcome mechanism (real data, Stage 3 only)
tmle_upd  <- run_tmle_targeting_step(g_fit, Q_fit)       # fluctuation / targeting update
tmle_est  <- extract_tmle_estimate(tmle_upd)             # psi, SE, CI, diagnostics

# Summarize all workflows side by side
summary_tbl <- summarize_cleanroom_results(
  list(match_fit, iptw_fit, tmle_est)
)
print(summary_tbl)
```

## Core Function Families

### A. Analysis lock and workflow metadata

- `create_analysis_lock()` — construct and serialize a complete analytic
  specification object, including the SuperLearner library, plasmode
  simulation settings, and (if applicable) the candidate set and selection rule
- `validate_analysis_lock()` — verify that all required fields are present and
  that the specification hash is reproducible

### B. Treatment mechanism / propensity-score estimation

- `fit_ps_superlearner()` — fit a SuperLearner-based propensity-score model
  from a locked specification; this is the default and recommended PS
  estimation approach in the package
- `fit_ps_glm()` — fit a logistic-regression propensity-score model; available
  for conventional workflows or simpler settings where SuperLearner is not
  needed
- `compute_ps_diagnostics()` — produce overlap plots, effective sample sizes,
  and covariate balance summaries from a PS fit

### C. Conventional workflow estimators

- `run_match_workflow()` — 1:1 nearest-neighbor matching with balance
  assessment and risk difference / ratio estimation
- `run_iptw_workflow()` — stabilized IPTW estimation with weight-distribution
  diagnostics and doubly robust augmentation

### D. Plasmode simulation and feasibility evaluation

- `run_plasmode_feasibility()` — generate plasmode outcomes from the observed
  W and A structure, run the full estimator pipeline on each replicate, and
  return bias, RMSE, and confidence-interval coverage
- `summarize_plasmode_results()` — tabulate and plot simulation-based
  performance metrics across workflows and DGPs
- `evaluate_tmle_candidates()` — evaluate all candidate TMLE specifications
  on plasmode simulations

### E. Modular TMLE components

These functions are intentionally separated so the package can respect
clean-room stage separation: the treatment mechanism can run in Stage 2a or
Stage 2b; the outcome mechanism runs on real data only in Stage 3; and the
targeting step follows only after both nuisance estimates are in hand.

- `fit_tmle_treatment_mechanism()` — estimate the treatment mechanism g(W);
  uses only covariates and treatment; can be called in Stage 2a or Stage 2b
- `fit_tmle_outcome_mechanism()` — estimate the initial outcome model Q(A,W);
  uses the real outcome and is called only in Stage 3
- `run_tmle_targeting_step()` — perform the fluctuation / targeting update
  using both nuisance estimates
- `extract_tmle_estimate()` — compute the targeted parameter estimate psi,
  standard error, confidence interval, and diagnostic summaries

### F. Candidate TMLE comparison / selection

- `fit_tmle_candidate_set()` — fit all TMLE candidates in the prespecified
  candidate set on plasmode simulations
- `select_tmle_candidate()` — apply a prespecified selection rule to
  simulation results to identify the best TMLE specification before
  real-outcome access

### G. Workflow-level summaries

- `fit_final_workflows()` — run final estimation for the selected workflow(s)
  after outcome unblinding
- `summarize_cleanroom_results()` — produce a structured, side-by-side
  summary of estimates across all fitted workflows
- `compare_fits()` — compare point estimates, confidence intervals, and
  diagnostic flags across workflow objects

## Why Use This Package in a Clean-Room Analysis?

Pre-registering an analysis plan is common, but conventional diagnostics —
checking overlap, balance, and weight distributions — can produce false
reassurance or false stops that are difficult to interpret in the absence of
estimator-level performance information. A good-looking PS overlap plot does
not guarantee that a particular TMLE specification will have low bias in the
sample at hand. Conversely, marginal overlap does not necessarily imply that
a targeted estimator will perform poorly.

Plasmode simulation bridges this gap. By evaluating the full estimator
pipeline on outcome-blind simulations derived from the real covariate
distribution, analysts can obtain pre-outcome evidence about the relative
performance of competing workflows. `cleanTMLE` supports this evaluation in
a structured way so that:

- the simulation design is locked before outcome access,
- selection decisions (if any) are rule-based and documented,
- the conventional and TMLE-based perspectives can be compared on the same
  simulated datasets, and
- the audit trail from specification to final estimate is preserved.

## Notes on Nuisance Estimation

TMLE requires estimating two nuisance components: the treatment mechanism
(the conditional probability of treatment given covariates) and the initial
outcome model. Neither component mathematically requires machine learning or
SuperLearner.

That said, `cleanTMLE` uses SuperLearner as the default and recommended
propensity-score estimation strategy. In clean-room TMLE workflows, flexible
treatment-model estimation is often desirable: the covariate set may be
high-dimensional, relationships may be nonlinear, and it can be difficult to
pre-specify a correctly-specified parametric form for the PS model before any
outcome is examined. SuperLearner ensemble learning provides a principled,
data-adaptive approach that can be fully pre-specified (by locking the
candidate library in Stage 1) without requiring ad hoc model selection after
outcome access.

Simpler regression-based nuisance models remain supported and may be
appropriate in low-dimensional settings or when a conventional PS workflow
is the primary analysis. The choice of nuisance strategy is declared in the
analysis lock and cannot be changed after outcome data are accessed.

The treatment mechanism, outcome mechanism, and targeting step are
intentionally separated in the API so that they can be run at the correct
clean-room stage. The treatment mechanism can be fitted and inspected in
Stage 2a or Stage 2b using only W and A. The outcome mechanism is fitted
on real data for the first time in Stage 3. The targeting step follows only
after both nuisance estimates are available. This design makes the clean-room
stage boundaries explicit in the code itself.

## Relationship to Existing R Packages

`cleanTMLE` is designed to complement, not replace, the broader ecosystem of
causal inference software in R. It provides the workflow scaffolding — staged
specification, blinded diagnostics, simulation-based validation, and structured
output — that the underlying estimation packages do not themselves enforce.
Estimation and nuisance-modeling functionality may be delegated to:

- [`tmle`](https://cran.r-project.org/package=tmle) — point-treatment TMLE
- [`SuperLearner`](https://cran.r-project.org/package=SuperLearner) — ensemble
  learning for nuisance models
- [`MatchIt`](https://cran.r-project.org/package=MatchIt) — propensity score
  matching
- [`WeightIt`](https://cran.r-project.org/package=WeightIt) — general
  propensity score weighting
- [`glmnet`](https://cran.r-project.org/package=glmnet) — regularized
  regression for nuisance estimation
- [`survey`](https://cran.r-project.org/package=survey) — design-weighted
  variance estimation

Users who prefer to call these packages directly are free to do so; `cleanTMLE`
provides wrappers that integrate them into the staged-workflow structure.

## Development Status

`cleanTMLE` is under active development. The staged-workflow layer described
in this README is being added incrementally alongside the existing estimation
functions. The public API — particularly function signatures for analysis-lock
and plasmode-simulation utilities — may still change before a stable release.
Issues and feature requests are welcome on the
[GitHub issue tracker](https://github.com/amertens/cleanTMLE/issues).

## Contributing

Contributions are welcome. Please open an issue to discuss proposed changes
before submitting a pull request. Code should follow the existing style
conventions and include tests for new functionality.

## License

MIT — see [LICENSE.md](LICENSE.md) for details.
