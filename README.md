# cleanTMLE

cleanTMLE separates the design of an observational analysis from its
outcome. You declare the plan, lock the design data without the outcome,
learn from an outcome-free plasmode simulation which estimands the design
can deliver and which TMLE candidate to use for each, record that decision
in a hashed dossier, and only then admit the outcome.

```r
library(cleanTMLE)
plan <- analysis_plan(
  outcome = "death_30d",
  estimands = c("ATE", "trimmed_ATE", "ATT", "ATO"),
  trim_band = c(0.05, 0.95),
  negative_controls = c(health_seeking = "health care use"),
  restrictions = list(no_transfer = ~ transfer == 0))
lock    <- create_analysis_lock(cr_example$design, treatment = "A",
                                covariates = c("age", "female", "severity", "prior_visits"),
                                plan = plan)
design  <- assess_design(lock)
sim     <- simulate_design(lock, design)
nc      <- negative_control_ladder(lock, design)
dossier <- design_report(lock, design, sim, nc, file = "dossier.html")
ub      <- unblind(lock, dossier, cr_example$outcomes, approved_by = "Review team")
fit     <- estimate_effect(ub)
export_design_log(ub, fit, file = "design_log.json")
```

The default plan runs 200 repetitions with a four-learner super learner
and 5-fold cross-fitting; `simulate_design()` prints a runtime estimate
after the first repetition. Set `K = 1`, a `"glm"` library, or fewer
`reps` in the plan for a quick look.

The design is tamper-evident, not tamper-proof: the outcome never enters
the lock or any design object, and any later change to the plan, the
design data or the decision is detectable through the hashes.

Version 0.3.0 remains installable from the `v0.3.0` tag.
