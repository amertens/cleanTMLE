# cleanTMLE 0.1.4 (development)

## FIORD selector: second stage (variance-method selection)

* **`bootstrap_rd_variance()`** adds a native nonparametric bootstrap
  standard error and percentile interval for the marginal risk
  difference (TMLE or stabilised IPTW). This is the principled variance
  method when the influence-function variance is conservative or invalid
  (stabilised IPTW; any estimator on a non-i.i.d. matched cohort).
* **`select_variance_method()`** implements the second stage of the
  FIORD two-stage selector (Nance et al. 2026): with the point estimator
  locked, it chooses the variance method whose oracle coverage on
  synthetic data is closest to nominal. Together with
  `select_tmle_candidate(rule = "fiord_two_stage")` (stage 1) this closes
  the gap to the full FIORD procedure.

## Loud failure modes (validity)

* **`estimate_gcomprisk()` no longer silently returns a zero risk
  difference.** When a treatment is identified but absent from the
  outcome model, g-computation predicted identical risk in both arms,
  yielding a degenerate RD of exactly zero with no signal to the user.
  The treatment is now injected into the outcome model and a `warning()`
  is emitted. The bootstrap re-uses the resolved formula so the warning
  fires once.
* **Cross-fitted TMLE targeting failures now warn.** A failed
  fluctuation step previously set the targeting coefficient to 0
  silently, reducing the estimator to an untargeted plug-in with an
  invalid influence-curve standard error. It now warns.
* **Survival-TMLE fallback warns.** When `tmle::tmle()` fails, the
  fallback to an unadjusted arm difference (which controls neither
  confounding nor censoring) now warns rather than returning the crude
  contrast as if it were the TMLE estimate.

## Callr out-of-process fit guard (opt-in)

* **`run_plasmode_dq_stress(fit_timeout = Inf)`** now accepts a `fit_timeout`
  argument (seconds). When a finite value is supplied and the `callr` package
  is available, each DQ-stress replicate's candidate fits run inside a
  persistent, killable `callr` subprocess (`.fit_candidates_bounded()`). If the
  fits exceed `fit_timeout` seconds the subprocess is killed, a fresh session
  is started for subsequent replicates, and every candidate for that replicate
  is recorded as `NA` — exactly as a fit error would be. This prevents a
  degenerate synthetic design (e.g. near-positivity scenarios that send
  `glmnet` into a runaway) from wedging the entire stress test. The default
  (`Inf`) preserves the existing in-process behaviour so existing scripts are
  unaffected; setting `fit_timeout = 120` is recommended when any SuperLearner
  candidate uses `SL.glmnet`.

## Tests

* New `test-ground-truth.R` checks the estimators against an analytic
  truth and against `tmle::tmle()` on shared inputs (point TMLE agrees
  with `tmle::tmle()` to within 0.01 on the risk-difference scale), and
  locks in the g-computation treatment-injection fix.

# cleanTMLE 0.1.3

This release hardens the GO/FLAG/STOP decision layer, centralises the
decision thresholds, strengthens the negative-control and missingness
checks, and adds a synthetic-data fidelity diagnostic.

## Decision-layer fixes

* **`authorize_outcome_analysis()` no longer masks a STOP.** Required
  stages are now matched on the stage *key* (text before the first
  colon), so `"Check Point 2"` matches the balance checkpoint but not
  `"Check Point 2c: DQ Stress"`. Decisions for a stage are reduced with
  STOP > FLAG > GO precedence rather than "last entry wins", and a new
  `block_on_any_stop = TRUE` argument makes *any* recorded STOP
  (including the optional DQ-stress gate) authoritative. Previously a
  balance STOP could be silently overwritten by a later same-prefix
  checkpoint, and the `gate_dq()` STOP was honoured only by recording
  order.

## Centralised, fingerprintable thresholds

* **`decision_thresholds()`** bundles every cohort / balance / plasmode /
  DQ / NCO threshold into one object; **`attach_decision_thresholds()`**
  stores it on the lock and records its own SHA-256 `thresholds_hash`, so
  the decision rule is tamper-evident (previously thresholds were passed
  ad hoc and were not fingerprinted). Helpers `dt_cohort()`,
  `dt_balance()`, `dt_plasmode()`, `dt_dq()`, `dt_nco()` extract the
  argument list for each step; `gate_dq(..., thresholds = dt)` reads from
  it directly.

## Checkpoint improvements

* **`checkpoint_residual_bias(rule = "equivalence")`** adds a TOST-style
  negative-control screen: a NC passes only when its `(1 - 2*alpha)` CI
  lies entirely inside a prespecified `null_band`, with optional
  `adjust = "bonferroni"` for the panel. The legacy significance rule
  rewarded low power (a noisy NC with a wide CI "passed"); equivalence is
  now recommended.
* **`checkpoint_cohort_adequacy()`** and **`checkpoint_balance()`** expose
  their formerly hard-coded STOP floors as arguments (`stop_n_per_arm`,
  `stop_min_events`, `stop_smd`); defaults reproduce prior behaviour.
  `checkpoint_balance()` also reports the number of covariates over the
  SMD threshold.

## Plasmode / DQ-stress additions

* **`run_plasmode_dq_stress(q0_library = ...)`** lets the synthetic-outcome
  generator Q0 use a SuperLearner library instead of a logistic GLM,
  matching the option already available in `run_plasmode_feasibility()`.
* **MAR covariate-missingness scenario** (`covariate_missingness_mar`):
  treatment-dependent missingness with median imputation, a stronger test
  than the existing MCAR scenario because the imputation is biased rather
  than merely inefficient.
* **`assess_dgp_fidelity()`** compares synthetic vs real covariate and
  treatment distributions (per-covariate SMD and KS, treatment-prevalence
  difference) and returns a GO/FLAG decision, so the analyst can defend
  that the plasmode generator is faithful enough to base candidate
  selection on.


# cleanTMLE 0.1.2

## New: `gate_dq()` makes the DQ stress test a hard checkpoint

* **`gate_dq(dq_results, candidate, max_abs_bias, min_coverage,
  max_rmse_ratio, ...)`** converts the degraded-scenario rows of a
  `plasmode_dq_results` object into a `cleantmle_checkpoint`. The
  checkpoint flips to STOP if any degraded row for the locked
  candidate exceeds the configured envelope, slots into `gate_all()`
  alongside the cohort-adequacy, balance, and residual-bias
  checkpoints, and is exported through the audit log. Previously the
  DQ stress test was computed and saved but never thresholded by the
  gate, which meant locked candidates could pass the gate even when a
  prespecified DQ row was outside the locked envelope.
* `run_simulation.R` now constructs a DQ checkpoint after every
  scenario's stress run and feeds it into `gate_all()`. Bundled
  results regenerated under `gate_dq()` show the expected
  Scenario C STOP from the unmeasured-confounding row.


# cleanTMLE 0.1.1

This release addresses the issues surfaced by applying cleanTMLE to the
Rescue.Co Kenya Trauma Registry (n = 1,693 ambulance patients with 22.8 %
outcome missingness). Most of the changes are real-data robustness fixes
and ergonomics improvements; a few extend the public API.

## New: IPCW-TMLE for missing-outcome sensitivity

* **`run_ipcw_tmle(lock, ps_fit, ...)`** is a new estimator for the
  full-cohort marginal risk difference when the outcome `Y` has
  missing-at-random rows (e.g., loss to follow-up). It fits a
  SuperLearner response model `P(R = 1 | A, W)`, builds stabilised
  inverse-probability-of-censoring weights with a 99th-percentile cap,
  and runs `tmle::tmle()` with `Delta = R` so the targeting step uses
  the censoring weights internally rather than dropping incomplete
  rows. Falls back to a complete-case TMLE weighted by the IPCW when
  the installed `tmle` version doesn't accept the `Delta` argument.
  The returned object is a `tmle_fit` and works with
  `summarize_cleanroom_results()`, `forest_plot()`, and `make_table2()`.
* **`run_crude_workflow()` and `run_match_workflow()`** now warn and
  compute on complete cases when `Y` has NAs; previously they silently
  returned `NA`.

## Real-data robustness (the marquee fixes)

* **`run_plasmode_feasibility()` and `run_plasmode_dq_stress()` no longer
  fail with `"Argument mu must be a nonempty numeric vector"`** when the
  lock's outcome column has NAs (improvement-prompt issue #13). The Q0
  GLM is now fit with `na.action = na.exclude` and predicted with
  `newdata = data` so the per-row probability vector has length `n`. A
  clear error is raised early when the outcome is fully NA (e.g. masked).
* **TMLE workflow functions accept locks with NA outcomes** (issue #13b).
  `fit_tmle_outcome_mechanism()`, `run_iptw_workflow()`, and
  `run_matched_tmle()` now subset to complete-Y rows for the Q-fit /
  IPTW estimator and emit a one-line warning ("inference is valid under
  MCAR only; consider IPCW for MAR"). `extract_tmle_estimate()` uses
  `n_eff = sum(!is.na(eic))` for the influence-function variance.
* **`fit_ps_superlearner()` no longer crashes with `"invalid connection"`
  when no parallel cluster is registered** (issue #3). Defaults to
  `cvControl$parallel = "seq"`; pass an explicit `cluster =` to
  parallelise. The function also now exposes `truncate` (default 0.01)
  and `cv_folds` (default 10) arguments and applies the truncation
  internally.

## Documentation drift fixed

* **`tmle_candidate()`** now accepts the deprecated `Q_library` /
  `Q_libraries` aliases for `q_library` with a deprecation warning,
  matching the older vignette (issue #12).
* **`expand_tmle_candidate_grid()`** accepts the deprecated
  `g_libraries` / `Q_libraries` form too.
* **`gate_check()`** accepts the more ergonomic
  `gate_check(plas, rmse_threshold = ..., coverage_threshold = ...)`
  form documented in the vignette, in addition to the original
  `gate_check(metrics, scenario_name, targets, method)` form (issues
  #14 and #7). It now also dispatches on `plasmode_results` and
  `plasmode_dq_results` objects directly.
* **Estimand fields are no longer claimed to be in the lock hash** in
  the manuscript or vignette; the SHA-256 fingerprint covers data
  shape + treatment / outcome / covariates / SL library / seed only,
  and the `attach_estimand()` / `lock_primary_tmle_spec()` records are
  protected by the audit log (manuscript revision; matches actual
  implementation).
* **Lock fingerprint is now real SHA-256** via the `digest` package
  (added to `Imports`); the previous 9-digit positional checksum is
  retained only as a clearly-labelled fallback when `digest` is
  unavailable.

## Audit and gate

* **`run_residual_confounding_stage()` no longer drops failed negative
  controls silently.** Each NC fit is wrapped in `tryCatch`; failures
  appear in `$summary_table` with `failed = TRUE` and the error
  message, and the count of failures is exposed as `$n_failed` so the
  gate can see partial NC coverage. Affected files:
  `R/stages.R:run_residual_confounding_stage`.
* **`authorize_outcome_analysis()`** now accepts a list of checkpoints
  via `checkpoints =` in addition to (or instead of) an audit (issue
  #10). When both are passed, the union of evidence is used.

## Ergonomics

* **`attrition_table()` is polymorphic** (issue #1): accepts a named
  numeric vector, a named list of step counts, or a data.frame with
  `step` and `n_remaining` (or `n`) columns, in addition to the
  original `(data, criteria)` interface.
* **`make_table1()` accepts a `cleanroom_lock`** directly (issue #4),
  using the lock's data, treatment column, and covariates.
* **`fit_ps_glm()` keeps `truncate = 0.01` as default** (issue #5,
  already in 0.1.0; documented in NEWS for completeness).
* **`default_dq_scenarios()`** is a new helper that returns the
  recommended scenario configuration as `"regulatory_standard"` (the
  default), `"exploratory"`, or `"stress"` (issue #9):

  ```r
  run_plasmode_dq_stress(lock, candidates,
                         data_quality_scenarios = default_dq_scenarios())
  ```
* **`print_locked_spec(lock)`** is a new helper that prints the
  candidate id, label, g/Q libraries, truncation, plasmode RMSE, and
  lock hash for the locked primary TMLE specification (issue #8).
* **`select_tmle_candidate(rule = "min_max_rmse", dq_results = ...)`**
  was added in the prior release and is now documented as the
  recommended rule when DQ stress results are available.

## Diagnostics

* **`estimate_design_precision()` and `summarize_event_support()`** now
  raise a clear error when the lock is masked or the outcome is all NA
  (issue #2), rather than returning a quietly-NA table.
* **`summarize_event_support()`** uses
  `lock$estimand$treatment_strategies` as arm labels when present
  (issue #15), so applied users see "Rescue.Co EMS" / "Other ambulance"
  rather than "Treated" / "Control".

## Visualisations

* **`plot()` method on `plasmode_dq_results`** produces the per-threat
  degradation-curve panel (RMSE / |bias| / coverage by severity, one
  facet per scenario, colour by candidate). Previously users had to
  draw this themselves from `summarize_dq_degradation()`.
* **`forest_plot()`** dispatches on data.frames and lists of fitted
  workflows in addition to `hr` objects, ordering rows
  TMLE -> IPTW -> Match -> Crude (most-to-least efficient) and
  shading the TMLE row.
* **`clever_covariate_plot(..., bin_extreme = TRUE)`** caps the H-axis
  at the 99th percentile and labels the count of bucketed extremes,
  preventing a handful of outliers from dominating the histogram.

## Manuscript-aligned changes (already in 0.1.0, restated)

* `select_tmle_candidate()` gained the `min_max_rmse` rule (with
  `dq_results = ...`) for selecting candidates by minimax RMSE across
  DQ stress scenarios. The methods manuscript recommends this rule.
* The lock fingerprint is now real SHA-256 via `digest`.

## Tests

New tests cover:

* `run_plasmode_feasibility()` and `run_plasmode_dq_stress()` on a lock
  with 25 % outcome NA (regression for #13).
* `fit_tmle_outcome_mechanism()` and `run_iptw_workflow()` on a lock
  with 25 % outcome NA (regression for #13b).
* `attrition_table()` with named-list, named-numeric, and data.frame
  inputs (regression for #1).
* `make_table1()` on a `cleanroom_lock` (regression for #4).
* `default_dq_scenarios()` returns valid scenario configs for each
  preset.
* `gate_check()` accepting the ergonomic short form and dispatching on
  `plasmode_results` / `plasmode_dq_results`.
* `select_tmle_candidate(rule = "min_max_rmse", dq_results = ...)`
  picks the robust candidate over a fragile one (already in 0.1.0; noted
  for completeness).

## Items not addressed in this release

Issue #6 (a uniform `cleantmle_safe()` wrapper for "warn + return NULL")
is deliberately not implemented: the recoverable error paths in
`make_table1`, `compute_matched_smds`, `forest_plot`, and `make_table2`
are now individually friendlier (lock acceptance, polymorphic input,
clearer errors), which removes most of the boilerplate motivating the
wrapper. A future release may revisit if the boilerplate reappears.

# cleanTMLE 0.1.0

Initial release. See the methods manuscript for a full description of
the staged workflow, the audit / gate machinery, and the plasmode DQ
stress test.
