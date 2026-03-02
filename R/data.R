#' Simulated Example Dataset
#'
#' A simulated dataset for demonstrating causal risk estimation methods.
#' Contains baseline covariates, treatment assignment, censoring,
#' time-to-event outcome, and competing risk indicator.
#'
#' @format A data.frame with 1000 observations and the following variables:
#' \describe{
#'   \item{id}{Subject identifier (1 to 1000).}
#'   \item{age}{Age in years (continuous, approximately Normal(50, 10)).}
#'   \item{sex}{Sex indicator (0 = female, 1 = male).}
#'   \item{biomarker}{Continuous biomarker value (Normal(0, 1)).}
#'   \item{comorbidity}{Comorbidity score (0, 1, or 2).}
#'   \item{treatment}{Treatment indicator (0 = control, 1 = treated).
#'     Treatment assignment depends on age, sex, and biomarker.}
#'   \item{time}{Observed follow-up time in months. This is
#'     min(event time, censoring time).}
#'   \item{event}{Event indicator: 1 = event of interest occurred,
#'     0 = censored or competing event.}
#'   \item{event_type}{Event type: 1 = primary event, 2 = competing event,
#'     0 = censored.}
#'   \item{censored}{Censoring indicator: 1 = administratively censored,
#'     0 = not censored (event occurred).}
#' }
#'
#' @source Simulated using [sim_func1()].
#'
#' @examples
#' data(example1, package = "cleanTMLE")
#' head(example1)
"example1"


#' Simulate a Causal Inference Example Dataset
#'
#' Generates a simulated dataset with baseline covariates, confounded
#' treatment assignment, time-to-event outcome with censoring, and
#' competing risks. The data generation process ensures that treatment
#' is confounded by baseline covariates, making it suitable for
#' demonstrating causal inference methods.
#'
#' @param n Number of observations (default 1000).
#' @param seed Random seed (default 123).
#' @param max_time Maximum follow-up time (default 36 months).
#'
#' @return A data.frame matching the structure of [example1].
#'
#' @details
#' ## Data generation process
#' 1. Covariates: age ~ N(50, 10), sex ~ Bernoulli(0.5),
#'    biomarker ~ N(0, 1), comorbidity ~ Multinomial(0, 1, 2).
#' 2. Treatment: logit(P(A=1)) = -1 + 0.02*age + 0.3*sex + 0.5*biomarker.
#' 3. Event time: exponential with rate depending on treatment and covariates.
#'    Treatment has a protective effect (HR ~ 0.7).
#' 4. Censoring time: exponential with rate depending on age.
#' 5. Competing event: with probability 0.1, the event is a competing risk.
#'
#' @examples
#' dat <- sim_func1(n = 500, seed = 42)
#' head(dat)
#'
#' @export
sim_func1 <- function(n = 1000, seed = 123, max_time = 36) {
  set.seed(seed)

  # Baseline covariates
  age <- rnorm(n, mean = 50, sd = 10)
  sex <- rbinom(n, 1, 0.5)
  biomarker <- rnorm(n, mean = 0, sd = 1)
  comorbidity <- sample(0:2, n, replace = TRUE, prob = c(0.5, 0.3, 0.2))

  # Treatment assignment (confounded)
  lp_trt <- -1 + 0.02 * age + 0.3 * sex + 0.5 * biomarker
  ps_true <- plogis(lp_trt)
  treatment <- rbinom(n, 1, ps_true)

  # Event time (exponential)
  # Treatment has protective effect (HR ~ 0.7)
  lp_event <- -3.5 + 0.01 * age + 0.2 * sex - 0.35 * treatment +
    0.15 * biomarker + 0.3 * comorbidity
  rate_event <- exp(lp_event)
  event_time <- rexp(n, rate = rate_event)

  # Censoring time (exponential, depends on age)
  rate_cens <- exp(-4 + 0.005 * age)
  cens_time <- rexp(n, rate = rate_cens)

  # Administrative censoring at max_time
  cens_time <- pmin(cens_time, max_time)

  # Competing risk: fraction of events are competing
  is_competing <- rbinom(n, 1, prob = 0.1)

  # Observed time and indicators
  obs_time <- pmin(event_time, cens_time)
  event_occurred <- as.integer(event_time <= cens_time)

  # Event type: 1 = primary, 2 = competing, 0 = censored
  event_type <- ifelse(event_occurred == 0, 0L,
                       ifelse(is_competing == 1, 2L, 1L))

  # Primary event indicator (for simple analyses)
  event <- as.integer(event_type == 1)

  # Censored indicator
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
