# Script to generate the example1 dataset
# Run this script to regenerate: source("data-raw/generate_example1.R")

# Source the sim_func1 function
# (In practice, after package is built, use causalRiskeR::sim_func1)

sim_func1_local <- function(n = 1000, seed = 123, max_time = 36) {
  set.seed(seed)
  age <- rnorm(n, mean = 50, sd = 10)
  sex <- rbinom(n, 1, 0.5)
  biomarker <- rnorm(n, mean = 0, sd = 1)
  comorbidity <- sample(0:2, n, replace = TRUE, prob = c(0.5, 0.3, 0.2))
  lp_trt <- -1 + 0.02 * age + 0.3 * sex + 0.5 * biomarker
  ps_true <- plogis(lp_trt)
  treatment <- rbinom(n, 1, ps_true)
  lp_event <- -3.5 + 0.01 * age + 0.2 * sex - 0.35 * treatment +
    0.15 * biomarker + 0.3 * comorbidity
  rate_event <- exp(lp_event)
  event_time <- rexp(n, rate = rate_event)
  rate_cens <- exp(-4 + 0.005 * age)
  cens_time <- pmin(rexp(n, rate = rate_cens), max_time)
  is_competing <- rbinom(n, 1, prob = 0.1)
  obs_time <- pmin(event_time, cens_time)
  event_occurred <- as.integer(event_time <= cens_time)
  event_type <- ifelse(event_occurred == 0, 0L,
                       ifelse(is_competing == 1, 2L, 1L))
  event <- as.integer(event_type == 1)
  censored <- as.integer(event_occurred == 0)

  data.frame(
    id = seq_len(n),
    age = round(age, 1),
    sex = sex,
    biomarker = round(biomarker, 3),
    comorbidity = comorbidity,
    treatment = treatment,
    time = round(obs_time, 2),
    event = event,
    event_type = event_type,
    censored = censored,
    stringsAsFactors = FALSE
  )
}

example1 <- sim_func1_local(n = 1000, seed = 123)

# Save as .rda
save(example1, file = "data/example1.rda", compress = "xz")
cat("example1 dataset saved to data/example1.rda\n")
cat("Dimensions:", nrow(example1), "x", ncol(example1), "\n")
cat("Columns:", paste(names(example1), collapse = ", "), "\n")
