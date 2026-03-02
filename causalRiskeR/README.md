# causalRiskeR

**Causal Inference for Risk Estimation with IPW, G-Computation, AIPW, and TMLE**

`causalRiskeR` provides a tidy, pipe-friendly interface for causal inference
in time-to-event and binary outcome settings. It implements:

- **IPW** (Inverse Probability Weighted) cumulative risk curves with IPTW and IPCW
- **G-computation** via plug-in Cox PH predictions
- **AIPW** (Augmented IPW / Doubly Robust) estimators
- **Weighted Cox regression** for hazard ratios with robust SEs
- **TMLE extensions** via `{tmle}`, `{survtmle}`, and `{lmtp}`

## Installation

```r
# Install from source
install.packages("causalRiskeR", repos = NULL, type = "source")

# Or use devtools/remotes
# devtools::install_github("causalRiskeR/causalRiskeR")
```

## Getting Started

```r
library(causalRiskeR)

# Simulate example data
dat <- sim_func1(n = 1000, seed = 42)

# Specify the causal model
spec <- specify_models(data = dat) |>
  identify_outcome(event, type = "time_to_event") |>
  identify_treatment(treatment, formula = ~ age + sex + biomarker) |>
  identify_censoring(censored, formula = ~ age + sex, model = "coxph")

# Estimate IPW cumulative risk curves
fit <- estimate_ipwrisk(spec, risk_time = c(6, 12, 24), nboot = 200)

# View results
print(fit)
plot(fit)
plot(fit, effect = "RD")

# Tables
make_table1(fit, vars = c("age", "sex", "biomarker"), weighted = TRUE)
make_table2(fit, risk_time = 24)
make_wt_summary_table(fit)

# Propensity score diagnostics
hist(fit, type = "ps")

# Hazard ratios
hr_fit <- estimate_ipwhr(spec)
forest_plot(hr_fit)

# TMLE for binary outcome at 24 months
dat$event_24 <- as.integer(dat$event == 1 & dat$time <= 24)
tmle_fit <- estimate_tmle_risk_point(
  data = dat,
  treatment = "treatment",
  outcome = "event_24",
  covariates = c("age", "sex", "biomarker")
)
print(tmle_fit)
```

## Core Functions

### Model Specification

| Function | Description |
|----------|-------------|
| `specify_models()` | Create specification object |
| `identify_outcome()` | Declare outcome variable and type |
| `identify_treatment()` | Declare treatment with PS model formula |
| `identify_censoring()` | Declare censoring with model |
| `identify_competing_risk()` | Declare competing risk indicator |
| `identify_subject()` | Declare subject ID |
| `identify_interval()` | Declare start/stop time interval |
| `identify_missing()` | Enable missingness weighting hooks |

### Estimation

| Function | Description |
|----------|-------------|
| `estimate_ipwrisk()` | IPW cumulative risk curves |
| `estimate_gcomprisk()` | G-computation risk curves |
| `estimate_aipwrisk()` | AIPW (doubly robust) risk curves |
| `estimate_ipwhr()` | Weighted Cox model hazard ratios |
| `estimate_ipwcount()` | IPW cumulative event count |

### TMLE Extensions

| Function | Description |
|----------|-------------|
| `estimate_tmle_risk_point()` | Point-treatment TMLE (binary outcome) |
| `estimate_surv_tmle()` | Survival TMLE for risk at times |
| `estimate_lmtp()` | Longitudinal TMLE via {lmtp} |

### Diagnostics & Tables

| Function | Description |
|----------|-------------|
| `make_table1()` | Baseline characteristics table |
| `make_table2()` | Results summary table |
| `make_wt_summary_table()` | Weight distribution summary |
| `extreme_weights()` | Top K extreme weights |
| `inspect_ipw_weights()` | Extract weight vectors |
| `compare_fits()` | Compare estimates across methods |
| `forest_plot()` | Forest plot for HRs |
| `hist()` | PS/weight histograms |

## Dependencies

**Required:** survival, stats, ggplot2, rlang, sandwich

**Optional:** nnet (multinomial treatment), tmle + SuperLearner (point TMLE),
survtmle (survival TMLE), lmtp (longitudinal TMLE), knitr + rmarkdown (vignettes)

## License

MIT
