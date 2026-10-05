.new_estimate <- function(estimand, prespecified, method, candidate, r, n, stress_tested,
                          note = NA_character_) {
  structure(list(
    estimand = estimand, prespecified = prespecified, method = method,
    candidate = candidate, estimate = r$estimate, se = r$se,
    ci_lower = r$ci_lower, ci_upper = r$ci_upper,
    p_value = 2 * stats::pnorm(-abs(r$estimate / r$se)), n = n,
    n_population = r$n_population, risk1 = r$risk1, risk0 = r$risk0,
    implausible = isTRUE(r$implausible), implausible_reason = r$implausible_reason %||% NA_character_,
    stress_tested = stress_tested, note = note), class = "cr_estimate")
}

.syntactic <- function(X) {
  colnames(X) <- make.names(colnames(X), unique = TRUE)
  X
}

.estimate_weighting <- function(X, A, Y, estimand, band) {
  Xs <- .syntactic(X)
  df <- data.frame(.Y = Y, .A = A, Xs)
  f <- stats::reformulate(colnames(Xs), response = ".A")
  est <- estimand
  if (estimand == "trimmed_ATE") {
    ps <- stats::fitted(stats::glm(f, data = df, family = stats::binomial()))
    df <- df[ps >= band[1] & ps <= band[2], , drop = FALSE]
    est <- "ATE"
  }
  w <- WeightIt::weightit(f, data = df, estimand = est, method = "glm")
  fit <- WeightIt::lm_weightit(.Y ~ .A, data = df, weightit = w)
  b <- unname(stats::coef(fit)[".A"])
  se <- sqrt(stats::vcov(fit)[".A", ".A"])
  r0 <- unname(stats::coef(fit)["(Intercept)"])
  gd <- .implausibility_check(b, df$.Y, df$.A)
  list(estimate = b, se = se, ci_lower = b - 1.96 * se, ci_upper = b + 1.96 * se,
       n_population = if (est == "ATT") sum(df$.A == 1) else nrow(df),
       risk1 = r0 + b, risk0 = r0, implausible = gd$implausible, implausible_reason = gd$reason)
}

.estimate_matching <- function(X, A, Y, estimand) {
  if (estimand != "ATT") stop("Matching estimates the ATT only.", call. = FALSE)
  if (!requireNamespace("MatchIt", quietly = TRUE))
    stop("method = 'matching' needs the MatchIt package.", call. = FALSE)
  Xs <- .syntactic(X)
  df <- data.frame(.Y = Y, .A = A, Xs)
  f <- stats::reformulate(colnames(Xs), response = ".A")
  # The unmatched treated units are reported through n_population and the
  # note, so MatchIt's own warning about them is not repeated.
  m <- withCallingHandlers(
    MatchIt::matchit(f, data = df, method = "nearest", estimand = "ATT",
                     distance = "glm", link = "linear.logit", caliper = 0.2,
                     std.caliper = TRUE),
    warning = function(w) {
      if (grepl("Fewer control units than treated units", conditionMessage(w), fixed = TRUE))
        invokeRestart("muffleWarning")
    })
  md <- MatchIt::match.data(m)
  fit <- WeightIt::lm_weightit(.Y ~ .A, data = md, weights = weights, cluster = ~subclass)
  b <- unname(stats::coef(fit)[".A"])
  se <- sqrt(stats::vcov(fit)[".A", ".A"])
  r0 <- unname(stats::coef(fit)["(Intercept)"])
  gd <- .implausibility_check(b, md$.Y, md$.A)
  list(estimate = b, se = se, ci_lower = b - 1.96 * se, ci_upper = b + 1.96 * se,
       n_population = sum(md$.A == 1), n_treated = sum(A == 1), risk1 = r0 + b, risk0 = r0,
       implausible = gd$implausible, implausible_reason = gd$reason)
}

.sl_equivalents <- function(learners) {
  unname(c(glm = "SL.glm", glmnet = "SL.glmnet", earth = "SL.earth", nnet = "SL.nnet",
           xgboost = "SL.xgboost")[learners])
}

.estimate_tte <- function(ub, estimand, candidate) {
  for (pkg in c("concrete", "data.table", "SuperLearner"))
    if (!requireNamespace(pkg, quietly = TRUE))
      stop("Time-to-event estimation needs the ", pkg, " package.", call. = FALSE)
  lock <- ub$lock
  p <- lock$plan
  X <- .syntactic(.design_matrix(lock$data, lock$covariates))
  A <- as.integer(lock$data[[lock$treatment]])
  tt <- ub$outcomes[[p$outcome[["time"]]]]
  ev <- ub$outcomes[[p$outcome[["event"]]]]
  keep <- !is.na(tt) & !is.na(ev)
  if (estimand == "trimmed_ATE") {
    folds <- .make_folds(nrow(X), p$K, p$seed)
    g <- .fit_g(X, A, p$ps_library, folds, p$V, p$seed)$g
    keep <- keep & g >= p$trim_band[1] & g <= p$trim_band[2]
  }
  # Keep the user's own column names: the plan's hazard formulas refer to them.
  df <- data.frame(lock$data[[lock$id]], tt, ev, A, X)
  names(df)[1:4] <- c(lock$id, p$outcome[["time"]], p$outcome[["event"]], lock$treatment)
  dt <- data.table::as.data.table(df[keep, , drop = FALSE])
  model <- c(stats::setNames(list(.sl_equivalents(p$ps_library)), lock$treatment), p$hazards)
  args <- concrete::formatArguments(
    DataTable = dt, EventTime = p$outcome[["time"]], EventType = p$outcome[["event"]],
    Treatment = lock$treatment, ID = lock$id, TargetTime = p$target_time,
    Intervention = concrete::makeITT(), Model = model,
    MinNuisance = candidate$truncation, Verbose = FALSE)
  # concrete reports its progress on the console even with Verbose = FALSE.
  est <- NULL
  utils::capture.output(suppressMessages(est <- concrete::doConcrete(args)))
  out <- as.data.frame(concrete::getOutput(est, Estimand = c("Risk", "RD"),
                                           Simultaneous = FALSE))
  out <- out[out$Estimator == "tmle" & out$Event == 1, , drop = FALSE]
  rd <- out[out$Estimand == "Risk Diff", , drop = FALSE][1, ]
  r1 <- out[out$Estimand == "Abs Risk" & out$Intervention == "A=1", "Pt Est"][1]
  r0 <- out[out$Estimand == "Abs Risk" & out$Intervention == "A=0", "Pt Est"][1]
  list(estimate = rd[["Pt Est"]], se = rd[["se"]], ci_lower = rd[["CI Low"]],
       ci_upper = rd[["CI Hi"]], n_population = sum(keep), risk1 = r1, risk0 = r0,
       implausible = FALSE, implausible_reason = NA_character_)
}

#' Estimate the effect after unblinding
#'
#' Reads the estimand and the candidate from the dossier. `method = "tmle"`
#' runs the stress-tested candidate with the same fitter the simulation
#' used; `"weighting"` (WeightIt) and `"matching"` (MatchIt, ATT only) are
#' comparators that were not stress-tested. A time-to-event plan runs
#' `concrete`.
#'
#' @param unblinded A [unblind()] result.
#' @param method `"tmle"`, `"weighting"` or `"matching"`.
#' @param estimand Optional; any other declared estimand is allowed and is
#'   labeled as not the prespecified primary.
#' @return A `cr_estimate`.
#' @export
estimate_effect <- function(unblinded, method = c("tmle", "weighting", "matching"),
                            estimand = NULL) {
  if (!inherits(unblinded, "cr_unblinded"))
    stop("`unblinded` must come from unblind().", call. = FALSE)
  method <- match.arg(method)
  lock <- unblinded$lock
  d <- unblinded$dossier
  p <- lock$plan
  verify_lock(lock)
  .verify_dossier(d)
  if (!identical(d$lock_hash, lock$lock_hash))
    stop("The dossier does not belong to this lock.", call. = FALSE)
  dec <- d$decision
  estimand <- estimand %||% dec$primary
  if (is.null(estimand) || is.na(estimand))
    stop("The dossier has no primary estimand; pass `estimand =`.", call. = FALSE)
  if (!estimand %in% p$estimands)
    stop("`estimand` must be one the plan declared: ",
         paste(p$estimands, collapse = ", "), ".", call. = FALSE)
  prespecified <- identical(estimand, dec$primary)
  note <- if (prespecified) NA_character_ else
    "This estimand is not the prespecified primary."
  candidate <- dec$candidates[[estimand]]

  if (p$outcome_type == "tte") {
    if (method != "tmle")
      stop("Time-to-event outcomes are estimated with concrete (method = 'tmle').",
           call. = FALSE)
    if (is.null(candidate))
      stop("The dossier found ", estimand, " infeasible; there is no candidate.",
           call. = FALSE)
    r <- .estimate_tte(unblinded, estimand, candidate)
    return(.new_estimate(estimand, prespecified, "concrete", candidate$id, r,
                         n = r$n_population, stress_tested = FALSE,
                         note = paste(c(note, "Support and truncation were stress-tested; the hazard models were not."),
                                      collapse = " ")))
  }

  X <- .design_matrix(lock$data, lock$covariates)
  A <- as.integer(lock$data[[lock$treatment]])
  Y <- unblinded$outcomes[[p$outcome]]
  cc <- !is.na(Y)
  X <- X[cc, , drop = FALSE]; A <- A[cc]; Y <- Y[cc]
  if (method == "tmle") {
    if (is.null(candidate))
      stop("The dossier found ", estimand, " infeasible; there is no stress-tested ",
           "candidate. Use a comparator method or another estimand.", call. = FALSE)
    folds <- .make_folds(length(Y), dec$K, dec$fold_seed)
    r <- fit_candidate(X, A, Y, estimand, candidate, p, folds, dec$fold_seed)
    if (isTRUE(r$failed)) stop("Estimation failed: ", r$message, call. = FALSE)
    return(.new_estimate(estimand, prespecified, "tmle", candidate$id,
                         as.list(r[1, , drop = FALSE]), n = length(Y),
                         stress_tested = TRUE, note = note))
  }
  r <- if (method == "weighting") .estimate_weighting(X, A, Y, estimand, p$trim_band) else
    .estimate_matching(X, A, Y, estimand)
  unmatched <- if (!is.null(r$n_treated) && r$n_population < r$n_treated)
    paste0("Only ", r$n_population, " of ", r$n_treated,
           " treated units found a match within the caliper.") else NULL
  .new_estimate(estimand, prespecified, method, NA_character_, r, n = length(Y),
                stress_tested = FALSE,
                note = paste(c(note, unmatched, "Comparator; not stress-tested."),
                             collapse = " "))
}
