# cleanTMLE 0.2.0 (development)

This release refocuses the package on deciding, before outcome access,
whether a comparison is estimable and with which estimand. The Rescue.Co
case study's protocol contrast, on which the ATE is not estimable while the
ATT, the overlap-weighted ATO and the trimmed ATE are, is the motivating
example throughout.

## The plasmode design is fixed (breaking for simulation results)

* **`run_plasmode_feasibility()` and `run_plasmode_dq_stress()` now default
  to the generate-treatment design**: covariate rows are resampled with
  replacement and treatment is drawn from a propensity model fitted on the
  real data. The previous behaviour, keeping every subject's observed
  treatment and simulating only the outcome, induces a positivity violation
  by construction (Shaw et al. 2025, arXiv:2504.11740), making
  propensity-based estimators look biased and undercover even when the
  source population has no violation. It survives as
  `design = "sample_treatment"`, which warns and records itself in the
  result. Simulation numbers change under the new default; that is the fix,
  not a regression.

* **`simulate_support()`** is the new outcome-blind support simulation. It
  draws synthetic data under the generate-treatment design from a
  prespecified outcome family (`support_surfaces()`) that spans a
  confounding axis (outcome dependence on the covariate direction that most
  predicts treatment), an effect-modification axis (so the ATE, ATT and ATO
  genuinely differ), and a complexity axis (curvature a main-terms outcome
  model misses). The truth for every estimand is computed from the
  generating model on each replicate, and the matched, trimmed, ATT and ATO
  estimators run on the same replicates as the full-cohort ATE, so
  matched-versus-full and trimmed-versus-full gaps decompose into estimand
  difference, efficiency loss and extrapolation error. The output includes
  the support map (bias and coverage per estimand per grid point, with a
  prespecified pass rule) that the design gate reads; runs below 50
  replicates are tagged demonstration-only. Q0 sources are recorded and
  fitting Q0 on the primary outcome requires `allow_outcome_q0 = TRUE`.

* **`check_locked_estimator()`** ships Petersen et al.'s (2012) parametric
  bootstrap as a final check on the locked estimator, labelled optimistic
  because the fitted models generate the truth.

## Support assessment and the estimand ladder

* **`assess_support()`** is the design team's graded overlap verdict, with
  no outcome access: percent of the sample outside a prespecified propensity
  band, maximum and 99th percentile IPTW weight, ESS by arm, the propensity
  c-statistic, a PASS/FLAG/SEVERE/FAIL verdict with a caveat string that
  travels with every downstream estimate, the near-deterministic-stratum
  table, and a shallow-tree search for multivariate violation regions in the
  spirit of the PoRT algorithm. Conditional non-overlap escalates PASS or
  FLAG to SEVERE, never to FAIL. Thresholds come from
  `support_thresholds()`; defaults are the Rescue.Co main pipeline's,
  calibrated against a sound and a broken reference fit. `plot()` draws the
  mirrored propensity histogram with the band shaded. Diagnostics use the
  untruncated scores, which `fit_ps_superlearner()` now retains as
  `ps_raw`. The vocabulary follows Muntner et al. (2024); the
  blinded-diagnostics gate follows Conover et al. (2025, JAMIA), where
  failed analyses are labelled inestimable.

* **`fit_ps()`** is the one front door for propensity fitting
  (SuperLearner, GLM, or externally supplied scores).

* **`estimand_feasibility()`** reports, per estimand (ATE, trimmed ATE at
  each level, ATT, ATC, ATO, matched ATT), the implied weights' maximum and
  99th percentile, ESS by arm, who is removed, the target population, and
  the graded verdict, so a design team can see that a cohort supports the
  ATT and the ATO where it cannot support the ATE.
  **`who_is_unsupported()`** profiles the patients outside the band
  (standardised differences of removed versus kept, overall and by arm).

* **`declare_estimand_ladder()`** pre-registers the primary estimand, the
  ordered fallbacks, the trigger verdict, and optional sensitivity floors
  on the lock, as a design-log entry. **`run_estimand_ladder()`** estimates
  the primary when feasible plus every feasible fallback, labels every row
  with its estimand, verdict, caveat and implausibility flags, refuses an
  infeasible primary unless an `override_reason` is recorded, and logs
  every switch.

## Estimators for the ladder

* **`run_att_tmle()`**: the complete-case ATT, named
  `E[Y(1)-Y(0) | A=1, outcome observed]`, never estimated through the
  censoring mechanism (under `Delta` the package ATT loses double
  robustness). The ATE from the same fit is returned beside it. All
  tmle::tmle delegation now flows through one internal argument builder
  with `prescreenW.g = FALSE` by default, and a regression test asserts the
  ATT's treatment-model specification is identical to the ATE's.

* **`estimate_ato()`**: the augmented overlap-weighted estimator (Mao, Li
  and Greene 2019) with an influence-function variance, the Hajek point
  estimate beside it, and the exact-balance check (max overlap-weighted
  SMD) reported.

* **`run_trimmed_tmle()`**: trimming with a propensity refit on the trimmed
  subset, verdict-gated fallback across levels, dropped counts by arm, and
  the Crump et al. (2009) data-adaptive threshold as `rule = "crump"`.

* **`implausibility_check()`**: the guard attached to every estimate (sign
  against crude, five times crude, beyond the largest arm rate or the
  observed outcome range, with noise floors). It flags, never suppresses.

## Reporting

* **`emulation_table()`** renders the protocol-versus-emulation two-column
  table of the TARGET reporting guideline (Cashin, Hansford, Hernan et al.,
  JAMA 2025) from what the lock records: the contrast, the eligibility
  trail, the outcome, the declared ladder, and the analysis specification.

## Design-stage tools

* **`create_contrast_locks()`** builds one lock per contrast from a
  multi-level treatment with shared covariates, shared negative controls
  (registered before any variance filtering, never dropped silently), and
  per-contrast design logs.

* **`run_negative_control_ladder()`** fits every registered control on
  every nested cohort and names the controls that fail on the full cohort
  and turn null after a restriction: evidence that restriction, not
  adjustment, removed the confounding.

* **`check_process_indicators()`**: the care-process collider check
  (association of treatment with each indicator conditional on severity and
  site), classified and thresholded in SD units, with no outcome access.


## Clean-room governance

* **Split, enforced entry point (`run_clean_tmle_preoutcome()` /
  `run_clean_tmle_primary()`).** `run_clean_tmle()` builds its lock internally
  and reads the outcome unconditionally, so it cannot require a pre-outcome
  authorisation. The new two-pass split does:
    - `run_clean_tmle_preoutcome()` runs the outcome-blind stages (cohort
      adequacy, PS balance, optional candidate selection / DQ stress, optional
      negative controls), assembles a reviewer-facing dossier, and returns an
      **unauthorised** lock plus a gate token.
    - `run_clean_tmle_primary(pre, authorization)` reads the outcome only after
      verifying that the token authorises the analysis, matches the lock hash,
      and matches the current audit fingerprint. Wrong-lock tokens, audits
      changed after authorisation, and STOP gates are all refused;
      `allow_outcome_access = TRUE` is the logged escape hatch.
* **Gate tokens are now tamper-evident.** `authorize_outcome_analysis()` binds
  the returned `pre_outcome_gate` to an `audit_fingerprint` (a sha256 of the
  audit's `(stage, decision)` multiset), so a checkpoint added or removed after
  authorisation is detectable downstream.
* **`build_dossier()` / `clean_tmle_dossier`.** The pre-outcome study dossier is
  now a first-class object (estimand + lock fingerprint, PS balance, candidate
  selection and DQ degradation, negative controls, and the pre-outcome decision
  with audit trail); it is returned as `pre$dossier` and has a `print` method.
* **`run_clean_tmle()` labelled honestly.** The convenience wrapper now builds
  its internal lock with `cleanroom_enabled = FALSE` and prints a one-line
  pointer to the enforced two-pass path, so the enforced and convenient paths
  are no longer the same name doing different things. Note the side effect: the
  lock returned in `res$lock` is now a plain (non-clean-room) lock, so reusing it
  for a later Stage-4 call (e.g. `sensitivity_truncation(res$lock)`) runs
  unguarded. Use `create_analysis_lock()` and the two-pass path if you need the
  returned lock to keep enforcing the outcome guard.
* **Hardened `run_clean_tmle_primary()` authorisation (fixes from an internal
  review).** The audit-fingerprint check now fails closed: a token whose
  fingerprint is `NA` (for example one minted from an empty audit) is refused
  rather than accepted, closing a bypass in which a forged empty-audit token
  could authorise a STOP `pre`. The audit fingerprint now also hashes each
  entry's `action` text (not just stage and decision), so editing the recorded
  selected candidate after authorisation is detected. And
  `run_clean_tmle_preoutcome()` now binds the selected candidate to the lock via
  `lock_primary_tmle_spec()`, so `run_clean_tmle_primary()` estimates with the
  selected candidate's truncation and outcome-model library instead of silently
  using the default locked library.

* **Stage-4 authorisation is now enforced (breaking).** The pre-outcome gate,
  not merely outcome masking, controls Stage 4. A cleanroom-enabled lock reaches
  an outcome estimator only if a passing pre-outcome gate has been recorded on
  it (`.outcome_authorized`); an unmasked-but-unauthorised lock, including one
  that was never masked, is now **refused** by `.check_outcome_access()` rather
  than silently allowed. Authorisation is recorded two ways:
    - `unmask_outcome(lock, original_lock, audit = <audit>)` checks the gate via
      `authorize_outcome_analysis()` before revealing the outcome and errors on a
      non-authorising gate (or a missing `audit`) unless
      `allow_unauthorized = TRUE` is passed (which warns and forces, for
      simulations or forced re-analysis).
    - `assert_outcome_authorized(audit, lock = <lock>)` checks the gate and, on
      success, returns the lock stamped as authorised. This is the recommended
      way to authorise an unmasked staged lock.
  `cleanroom_enabled = FALSE` locks (e.g. `create_simple_lock()`) are exempt, and
  the per-estimator `allow_outcome_access = TRUE` override still bypasses the
  guard. Existing code that ran Stage 4 without authorising must now authorise
  (or pass an override); the bundled staged vignette and `run_simulation.R` were
  updated accordingly.

## Bug fixes

* **`extract_tmle_estimate()`** now returns the model-based (TMLE plug-in)
  treatment-specific risks as `risks$treated` and `risks$control`; their
  difference equals the reported ATE. This lets
  reporting code (e.g. `clean_risk_report_table()`) show the
  confounding-adjusted arm risks the estimator actually targeted instead of
  reconstructing them as `crude_arm_risk +/- ATE/2`, which forced symmetry
  around the crude midpoint and could fall outside `[0, 1]`. The staged-analysis
  vignette's risk-report and E-value chunks now use these fields; no ATE value
  changes. Covered by a new assertion in `test-cleanroom.R`.

* **`run_plasmode_dq_stress()`** no longer aborts when a degraded scenario cell
  has fewer than two converging replicates. The SE-calibration ratio (`se_cal`)
  is now computed with an explicit `NA` guard, because the empirical SD is
  undefined for a single value; this matches the behaviour already used in
  `run_plasmode_feasibility()`. Affected cells report `NA` for `emp_sd` and
  `se_cal` rather than raising `"missing value where TRUE/FALSE needed"`, which
  previously aborted the whole stress test. No previously reported value
  changes. Covered by a new regression test in `test-plasmode_dq.R`.

## Data-quality stress test: positivity is now a default threat

* **`default_dq_scenarios()`** now returns five threats for every preset.
  `near_positivity` joins covariate missingness, treatment misclassification,
  outcome misclassification, and unmeasured confounding as a default scenario.
  The near-positivity threat amplifies the centred propensity-score log-odds by
  a `slopes` multiplier so a subgroup approaches deterministic treatment and the
  estimated propensity scores reach the boundary. Default slope grids are mild
  for `regulatory_standard` and `exploratory` and heavier for `stress`. This is
  the same mechanism that `run_plasmode_dq_stress()` already dispatched and that
  the candidate-divergence study exercises, so the default run and the
  demonstration now use an identical code path.

## FIORD selector: second stage (variance-method selection)

* **`bootstrap_rd_variance(estimator = "match_tmle")`** is a new estimator
  option that implements the theoretically principled bootstrap for matched
  estimators (Abadie and Imbens 2008, 2016). Each bootstrap resample re-runs
  the full pipeline: GLM propensity-score estimation, greedy 1:1 nearest-
  neighbour matching on the logit-PS, and TMLE on the matched cohort. This
  captures matching-draw variance (the variability from which controls are
  selected when the data are re-sampled) that the standard paired-difference
  SE ignores, correcting anti-conservative IF coverage in good-overlap
  settings. The internal `.match_nn_pairs()` helper mirrors
  `run_match_workflow()` exactly.
* **`select_variance_method(estimator = "match_tmle")`** is now supported.
  The `"influence"` branch uses the paired-difference SE on the matched
  cohort; the `"bootstrap"` branch re-runs the full pipeline via
  `bootstrap_rd_variance()`.
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
* **`run_match_workflow()`** Roxygen now includes a `@section Variance`
  explaining the paired-difference SE limitation and pointing to
  `bootstrap_rd_variance(estimator = "match_tmle")` as the recommended
  alternative. References Abadie and Imbens (2008, 2016).

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
