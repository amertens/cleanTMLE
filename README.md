# cleanTMLE

**Staged Clean-Room Causal Analysis: Propensity Score and TMLE Workflows**

## Overview

`cleanTMLE` supports a structured, outcome-blinded workflow for observational
causal inference. The package implements the scaffolding needed to plan and
execute a staged clean-room analysis in which analytic decisions — including
estimator choice, nuisance-model specifications, and selection rules — are
locked before any investigator accesses the study outcome data. It covers
three workflow families: conventional propensity-score methods, a
fixed-specification TMLE workflow, and a simulation-selected TMLE workflow
in which a prespecified selection rule is used to choose among candidate
specifications based on blinded plasmode simulation performance. In addition
to final estimation, the package provides tools for overlap diagnostics,
matching, inverse probability of treatment weighting (IPTW), plasmode
simulation, and workflow-level summaries, making the clean-room process
reproducible and auditable from specification to final report.

## Key Features

- **Analysis lock and staged workflow support** — record and validate analytic
  decisions before outcome access
- **Propensity score estimation and diagnostics** — logistic-regression PS
  models with overlap checks and standardized mean differences
- **Matching and IPTW helpers** — 1:1 nearest-neighbor matching and stabilized
  IPTW with weight-distribution summaries
- **Plasmode simulation tools** — outcome-blind simulation from the observed
  covariate and treatment distributions to evaluate estimator performance
- **TMLE workflow support** — fixed-specification TMLE with prespecified
  nuisance strategies; regression-based or flexible learner-based nuisance
  models are both supported
- **Candidate TMLE specification comparison / selection** — fit a set of
  candidate TMLE specifications and apply a prespecified rule to select among
  them before real-outcome access
- **Workflow-level summaries and decision support** — side-by-side comparison
  of conventional and TMLE workflow results with structured output objects

## Workflow Overview

The package organizes analysis into three stages, with Stage 2 divided into
two substages. Stages 1 and 2 occur before any investigator accesses the study
outcome; Stage 3 is the unblinded estimation step.

```
Stage 1: Analysis lock (pre-specification)
   ├─ Stage 2a: Traditional feasibility diagnostics
   │     (PS overlap, covariate balance, weight distributions)
   └─ Stage 2b: Plasmode-simulation feasibility evaluation
         (estimator performance under outcome-blind simulated DGPs)
               ↓
         Stage 3: Unblinded final estimation
```

**Stage 1 — Pre-specification and analysis lock.**
All analytic decisions are written down and validated before outcome data are
examined: estimand, covariate set, PS model formula, nuisance-learning
strategy, plasmode simulation parameters, and (for the simulation-selected
workflow) the candidate set and selection rule. `create_analysis_lock()` and
`validate_analysis_lock()` serialize and verify this specification.

**Stage 2a — Traditional feasibility diagnostics.**
Using only treatment and covariate data, analysts assess whether the PS model
produces adequate overlap (effective sample sizes, trimming needs) and whether
covariate balance is achievable. These diagnostics may prompt revisions to the
specification while the outcome remains blinded.

**Stage 2b — Plasmode-simulation feasibility evaluation.**
The full estimator pipeline — including nuisance modeling, TMLE targeting, and
variance estimation — is evaluated on simulated datasets generated from the
observed covariate and treatment distributions. Outcome values are simulated
from plausible data-generating processes (DGPs) so the true causal effect is
known. This step allows analysts to assess bias, variance, and coverage of
candidate workflows under realistic conditions before touching the real outcome.

**Stage 3 — Unblinded estimation.**
After the analysis plan is locked and feasibility evaluation is complete,
the final estimators are applied to the real outcome data. The workflow family
(conventional PS, fixed TMLE, or simulation-selected TMLE) and any
simulation-driven selections were committed before this step.

## Supported Workflow Families

### Conventional propensity-score workflow

Fits a logistic-regression PS model, checks overlap, and produces estimates
via 1:1 nearest-neighbor matching and stabilized IPTW. Covariate balance
tables and weight diagnostics are generated for reporting.

### Fixed-specification TMLE workflow

Fits a TMLE estimator using a nuisance strategy that was locked in Stage 1.
The nuisance strategy can be as simple as parametric logistic regression or
as flexible as an ensemble learner (`SuperLearner`) or regularized models
(`glmnet`, GAMs). Treatment-model diagnostics are generated from covariate and
treatment data only, without accessing the outcome. Plasmode simulation
(Stage 2b) validates estimator performance before unblinding.

### Simulation-selected TMLE workflow

A candidate set of TMLE specifications (e.g., varying the outcome-model
learner or propensity-score model) is evaluated on plasmode simulations. A
prespecified selection rule (e.g., lowest root-mean-squared error across
simulated datasets) identifies the best candidate. The selected specification
is then locked before the real outcome is accessed and used in Stage 3.

## Package Philosophy

Stage 2a asks whether the treatment model and covariate overlap look
acceptable — it is a necessary check, but it does not assess how the full
estimator pipeline will perform. Stage 2b asks whether the estimator, end to
end, performs adequately under plausible data-generating processes. Both
perspectives are useful and not mutually exclusive. `cleanTMLE` is designed to
support both, and to make the transition from one to the other — and ultimately
to Stage 3 — structured and auditable rather than ad hoc.

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
  ps_formula    = treatment ~ age + sex + biomarker,
  tmle_strategy = "glm",          # or "SL" for SuperLearner
  plasmode_reps = 200,
  seed          = 42
)

validate_analysis_lock(lock)      # checks completeness and reproducibility

# ── Stage 2a: Traditional feasibility diagnostics ─────────────────────────

ps_fit  <- fit_conventional_ps(lock)
diag    <- compute_ps_diagnostics(ps_fit)
print(diag)                        # overlap plots, ESS, SMDs

# ── Stage 2b: Plasmode-simulation feasibility evaluation ──────────────────

sim_results <- run_plasmode_feasibility(
  lock         = lock,
  effect_sizes = c(0.05, 0.10),  # plausible true risk differences
  reps         = lock$plasmode_reps
)
print(sim_results)                 # bias, RMSE, coverage by workflow

# ── Stage 3: Unblinded final estimation ───────────────────────────────────

# Conventional PS workflow
match_fit <- run_match_workflow(lock, ps_fit)
iptw_fit  <- run_iptw_workflow(lock, ps_fit)

# Fixed-specification TMLE workflow
tmle_fit  <- fit_tmle_fixed(lock)

# Simulation-selected TMLE workflow (selection rule applied before outcome access)
candidates  <- fit_tmle_candidate_set(lock)
best_spec   <- select_tmle_candidate(candidates, sim_results, rule = "min_rmse")
final_tmle  <- fit_final_workflows(lock, selected_spec = best_spec)

# Summarize all workflows side by side
summary_tbl <- summarize_cleanroom_results(
  list(match_fit, iptw_fit, tmle_fit, final_tmle)
)
print(summary_tbl)
```

## Core Function Families

### Analysis lock and workflow metadata

- `create_analysis_lock()` — construct and serialize a complete analytic
  specification object
- `validate_analysis_lock()` — verify that all required fields are present and
  that the specification hash is reproducible

### Propensity score estimation and diagnostics

- `fit_conventional_ps()` — fit a logistic-regression propensity score model
  from a locked specification
- `compute_ps_diagnostics()` — produce overlap plots, effective sample sizes,
  and covariate balance summaries from a PS fit

### Matching and weighting estimators

- `run_match_workflow()` — 1:1 nearest-neighbor matching with balance
  assessment and risk difference / ratio estimation
- `run_iptw_workflow()` — stabilized IPTW estimation with weight-distribution
  diagnostics and doubly robust augmentation

### Plasmode simulation and workflow validation

- `run_plasmode_feasibility()` — generate outcome-blind plasmode datasets and
  evaluate each registered workflow on each replicate, returning bias, RMSE,
  and confidence-interval coverage
- `summarize_plasmode_results()` — tabulate and plot simulation-based
  performance metrics across workflows and DGPs

### TMLE nuisance strategies and estimation

- `fit_tmle_fixed()` — fit a TMLE estimator using the nuisance strategy
  declared in the analysis lock (regression, `glmnet`, GAM, or SuperLearner)
- `fit_tmle_candidate_set()` — fit all TMLE candidates in the prespecified
  candidate set

### Candidate TMLE comparison / selection

- `select_tmle_candidate()` — apply a prespecified selection rule to
  simulation results to identify the best TMLE specification
- `fit_final_workflows()` — run final estimation for the selected workflow(s)
  after outcome unblinding

### Reporting and summaries

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

TMLE requires estimating two nuisance components: the propensity score (the
conditional probability of treatment given covariates) and the initial outcome
model. Neither requires machine learning. By default, `cleanTMLE` uses
standard logistic regression for both nuisance components, which is often
adequate for moderate-dimensional settings.

Flexible nuisance estimation — using regularized regression (`glmnet`),
generalized additive models, or ensemble learners (`SuperLearner`) — can also
be specified and is often recommended when the covariate set is large or
relationships are likely nonlinear. The nuisance strategy is declared as part
of the analysis lock so that it cannot be changed after outcome access.

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
