# cleanTMLE 0.1.1

This release addresses the issues surfaced by applying cleanTMLE to the
Rescue.Co Kenya Trauma Registry (n = 1,693 ambulance patients with 22.8 %
outcome missingness). Most of the changes are real-data robustness fixes
and ergonomics improvements; a few extend the public API.

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
  MCAR only -- consider IPCW for MAR"). `extract_tmle_estimate()` uses
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
