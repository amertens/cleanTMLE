# An integration test on a cohort whose timing problems are the ones a
# trauma-transport registry has: covariates recorded or fixed after the
# transport decision, a follow-up-call covariate, a route-based restriction and
# an outcome-observation indicator.

test_that("a cohort with Rescue.Co's timing problems is caught at the lock", {
  n <- 600L
  d <- withr::with_seed(61, {
    sev <- stats::rnorm(n)
    age <- stats::rnorm(n)
    A <- stats::rbinom(n, 1, stats::plogis(0.6 * sev))
    pay_disp <- stats::rbinom(n, 1, stats::plogis(-0.3 + 0.8 * A))
    arr_gcs <- stats::rnorm(n, 0.5 * A - 0.5 * sev)
    fu_water <- stats::rbinom(n, 1, 0.5)
    direct <- stats::rbinom(n, 1, 0.8)
    seen_6m <- stats::rbinom(n, 1, 0.7)
    data.frame(id = seq_len(n), A = A, sev = sev, age = age, pay_disp = pay_disp,
               arr_gcs = arr_gcs, fu_water = fu_water, direct = direct,
               seen_6m = seen_6m)
  })
  plan <- fast_plan(
    time_zero = list(event = "transport decision", measured_as = "time of injury",
                     follow_up_from = "arrival", response = "r0"),
    eligibility = list(arrived_alive = list(description = "Arrived alive",
                                            assessed = "after", response = "r1")),
    restrictions = list(direct = ~ direct == 1),
    reps = 10L, K = 2L)

  row <- function(column, role, recorded = "before", fixed = "before",
                  follow_up = FALSE, reason = "") {
    data.frame(column = column, variable = column, source = "registry",
               follow_up = follow_up, recorded = recorded, fixed = fixed,
               role = role, reason = reason, sensitivity = "",
               stringsAsFactors = FALSE)
  }
  base <- rbind(
    row("id", "id"),
    row("A", "treatment", fixed = "at"),
    row("sev", "covariate"),
    row("age", "covariate", recorded = "after", reason = "asked at enrolment"),
    row("direct", "restriction", recorded = "after", fixed = "after", reason = "r3"))
  covs <- c("sev", "age")
  lk <- function(dict, covariates = covs)
    create_analysis_lock(d, "A", covariates, plan, dictionary = dict)

  # A covariate paid for after the transport decision.
  dict <- rbind(base, row("pay_disp", "covariate", recorded = "after", fixed = "after"))
  expect_error(lk(dict, c(covs, "pay_disp")), "fixed after time zero")
  # A covariate from the follow-up call.
  dict <- rbind(base, row("fu_water", "covariate", follow_up = TRUE))
  expect_error(lk(dict, c(covs, "fu_water")), "follow-up source")
  # An arrival measurement: refused as a covariate, kept as a node after treatment.
  dict <- rbind(base, row("arr_gcs", "after_treatment", recorded = "after", fixed = "after",
                          reason = "arrival measurement"))
  expect_error(lk(dict, c(covs, "arr_gcs")), "node after treatment")
  expect_s3_class(lk(dict), "cr_lock")
  # A restriction variable fixed after time zero needs a reason.
  nodirect <- rbind(base[base$column != "direct", ],
                    row("direct", "restriction", recorded = "after", fixed = "after"))
  expect_error(lk(nodirect), "Restriction variable `direct`")
  expect_s3_class(lk(base), "cr_lock")
  # An outcome-observation indicator needs a reason.
  dict <- rbind(base, row("arr_gcs", "after_treatment", recorded = "after", fixed = "after",
                          reason = "arrival measurement"),
                row("seen_6m", "observation", recorded = "after", fixed = "after",
                    follow_up = TRUE))
  expect_error(lk(dict), "Outcome-observation indicator `seen_6m`")
  dict$reason[dict$column == "seen_6m"] <- "approved before the outcome was read"
  lock <- lk(dict)
  expect_s3_class(lock, "cr_lock")

  design <- assess_design(lock)
  sim <- suppressMessages(simulate_design(lock, design))
  dossier <- design_report(lock, design, sim)

  idt <- dossier$identification
  fl <- idt[idt$flag == "flag", ]
  expect_identical(paste(fl$component, fl$item, sep = ": "),
                   c("Start of follow-up: follow-up", "Eligibility: arrived_alive",
                     "Outcome-observation indicator: seen_6m", "Restriction variable: direct"))
  nt <- idt[idt$flag == "note", ]
  expect_identical(nt$component,
                   c("Adjustment set: fixed before time zero", "Node after treatment"))
  expect_identical(nt$declared[nt$component == "Adjustment set: fixed before time zero"], "age")
  expect_identical(nt$item[nt$component == "Node after treatment"], "arr_gcs")

  log <- export_design_log(dossier, format = "muntner")
  expect_identical(sum(log$stage == "Identification"), 4L)
})
