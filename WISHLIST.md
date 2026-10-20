# cleanTMLE wishlist

Additions planned after the first public release (0.5.0). None of them is
promised for a particular version. Each would be added as arguments to the
existing functions rather than as new ones where possible.

## A fuller outcome-blind design stage

These came from re-running the design stage of a trauma-registry study whose
unblinded analysis had to depart from its protocol many times. Each is meant
to run in every study, not only where a problem is suspected.

- **How the outcome is observed.** Declare the stages after treatment in time
  order, each with an optional state measured after treatment, a censoring
  indicator and an event that settles the outcome (such as death before a
  follow-up interview). Refuse single-stage censoring weights when an event
  that settles the outcome comes before censoring. Report, without outcomes,
  the share lost at each censoring point by arm, how well treatment, baseline
  covariates and the measured history predict loss, and censoring positivity.
  Simulate loss at declared strengths and evaluate a longitudinal TMLE
  candidate, as a wrapper around lmtp checked against ltmle, beside complete
  case. Start with the ATE and the trimmed ATE. `design_report()` would also
  re-run the dictionary checks.
- **Several comparisons and restricted populations in one lock.** A treatment
  with more than two levels and declared comparisons; restrictions that define
  target populations for the estimand ladder; a map of comparison by
  population by estimand; the composition of each comparator; a description
  of the patients without support (small cells suppressed), with a field for
  clinicians to mark each region structural or practical; overlap cut-points
  set in the plan.
- **Gates and tolerances.** A bias tolerance for each outcome tied to the
  smallest difference that would change a decision, with a warning when it is
  larger than half that difference; a required interval half-width; negative
  control null bands relative to each control's prevalence; a gate on how far
  a lower rung of the ladder moves the population (share of each arm kept,
  standardized shift on named covariates), with the overlap-weighted effect
  allowed only when the plan declares it acceptable.
- **Selection after time zero.** A quantitative bias analysis declared and
  locked in the plan and reported beside the estimate, with outcome-free
  diagnostics such as injury severity among arrivals by arm.
- **Negative controls by domain.** List each unmeasured confounder in the
  causal diagram and show which have no control.
- **Smaller items.** State in the dossier when the effect among the treated
  refits its propensity score; log the number of patients dropped by complete
  case; let comparators, learner libraries and seeds be declared candidates;
  print a registration string for the lock; add example data in which a
  negative control fails and an override is logged, one with an
  outcome-observation structure, and one whose true outcome model lies
  outside the declared surfaces; an incremental propensity score intervention
  through lmtp.

## Parity with causalRisk

- Cumulative risk curves over follow-up, with plots (now: the risk at one
  target time).
- Competing risks.
- Censoring that depends on time-varying covariates.
- Table 1 by arm, unweighted and weighted; Table 2 with patients,
  person-time, events and estimates.
- Weight summaries and a table of the most extreme weights.
- Forest plots of estimates.
- The risk ratio with its interval, and other effect measures (the risk ratio
  is now computed only for the E-value); the hazard ratio.
- Subgroup analyses and re-estimation from a result object.
- Bootstrap intervals and clustered data.
- Cumulative counts for recurrent events.

## Housekeeping found in review

- Read plan fields by exact name, and reject unknown fields in an
  eligibility criterion as `time_zero` does.
- Tests for a few untested branches of the dictionary checks, and tighter
  message patterns in some tests.
- Show time zero when a plan is printed.
