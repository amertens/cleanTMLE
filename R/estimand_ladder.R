# The estimand feasibility table and the pre-registered estimand ladder: the
# machinery that lets a design team say, before outcomes are unlocked, which
# estimands this design can deliver, and that makes switching to a fallback a
# logged design decision rather than a silent substitution.

.ladder_estimands <- c("ATE", "trimmed_ATE", "ATT", "ATC", "ATO",
                       "matched_ATT")

.estimand_populations <- c(
  ATE = "everyone in the cohort",
  trimmed_ATE = "patients whose covariates could plausibly have produced either arm",
  ATT = "the treated (rows with an observed outcome)",
  ATC = "the controls",
  ATO = "patients in clinical equipoise, weighted g(1-g)",
  matched_ATT = "treated patients with a matched control")

# Design-decision log on the lock: a plain data frame, no tokens.
#' @keywords internal
.log_design_decision <- function(lock, type, note) {
  entry <- data.frame(timestamp = format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
                      type = type, note = note, stringsAsFactors = FALSE)
  lock$design_log <- if (is.null(lock$design_log)) entry else
    rbind(lock$design_log, entry)
  lock
}


#' Which Estimands Does This Design Support?
#'
#' One row per estimand: the weights it implies, their maximum and 99th
#' percentile, the effective sample size by arm under those weights, how many
#' patients it removes or leaves untouched, its target population, and the
#' graded verdict from the same thresholds as [assess_support()]. The point
#' of the table is the contrast across rows: on a cohort with a severe
#' positivity violation the ATE's weights explode while the ATT's control
#' weights stay bounded by `max g / (1 - g)` and the ATO's weights are
#' bounded by one by construction, so the design can support the second and
#' third where it cannot support the first.
#'
#' @param ps_fit A `ps_fit`; untruncated scores are used.
#' @param estimands Subset of `c("ATE", "trimmed_ATE", "ATT", "ATC", "ATO",
#'   "matched_ATT")`.
#' @param thresholds A [support_thresholds()] object.
#' @param trim_levels Trim levels reported for `trimmed_ATE`.
#' @param caliper_sd Caliper (SDs of the logit propensity) for
#'   `matched_ATT`. Default 0.2.
#' @return An object of class `estimand_feasibility`: a data.frame with one
#'   row per estimand (and one per trim level), plus the thresholds.
#' @references Li F, Morgan KL, Zaslavsky AM (2018) JASA 113:390-400.
#'   Crump RK et al. (2009) Biometrika 96:187-199.
#' @export
estimand_feasibility <- function(ps_fit,
                                 estimands = .ladder_estimands,
                                 thresholds = support_thresholds(),
                                 trim_levels = c(0.05, 0.10),
                                 caliper_sd = 0.2) {
  if (!inherits(ps_fit, "ps_fit"))
    stop("`ps_fit` must be a ps_fit object.", call. = FALSE)
  estimands <- match.arg(estimands, .ladder_estimands, several.ok = TRUE)
  g <- pmin(pmax(as.numeric(ps_fit$ps_raw %||% ps_fit$ps), 1e-6), 1 - 1e-6)
  A <- as.integer(ps_fit$data[[ps_fit$treatment]])
  band <- thresholds$band

  ess <- function(w, arm) {
    wa <- w[A == arm]
    if (!length(wa) || all(wa == 0)) return(0)
    sum(wa)^2 / sum(wa^2)
  }
  row_for <- function(estimand, w, n_removed_t, n_removed_c, pct_outside,
                      note = "") {
    mx <- max(w[w > 0])
    base_ed <- if (startsWith(estimand, "trimmed_ATE")) "trimmed_ATE" else
      estimand
    data.frame(
      estimand = estimand,
      target_population = unname(.estimand_populations[base_ed]),
      max_weight = round(mx, 1),
      p99_weight = round(stats::quantile(w[w > 0], 0.99, names = FALSE), 1),
      ess_treated = round(ess(w, 1L), 1),
      ess_control = round(ess(w, 0L), 1),
      n_removed_treated = n_removed_t,
      n_removed_control = n_removed_c,
      verdict = .support_verdict(pct_outside, mx, thresholds),
      note = note,
      stringsAsFactors = FALSE)
  }

  pct_outside_full <- 100 * mean(g < band[1] | g > band[2])
  rows <- list()
  for (ed in estimands) {
    if (ed == "ATE") {
      w <- ifelse(A == 1, 1 / g, 1 / (1 - g))
      rows[[length(rows) + 1L]] <- row_for("ATE", w, 0L, 0L,
                                           pct_outside_full)
    } else if (ed == "ATT") {
      w <- ifelse(A == 1, 1, g / (1 - g))
      rows[[length(rows) + 1L]] <- row_for(
        "ATT", w, 0L, 0L, 0,
        note = sprintf("control weights bounded by max g/(1-g) = %.1f",
                       max(g / (1 - g))))
    } else if (ed == "ATC") {
      w <- ifelse(A == 1, (1 - g) / g, 1)
      rows[[length(rows) + 1L]] <- row_for(
        "ATC", w, 0L, 0L, 0,
        note = sprintf("treated weights bounded by max (1-g)/g = %.1f",
                       max((1 - g) / g)))
    } else if (ed == "ATO") {
      w <- ifelse(A == 1, 1 - g, g)
      rows[[length(rows) + 1L]] <- row_for(
        "ATO", w, 0L, 0L, 0,
        note = "weights bounded by 1; exact covariate balance under logistic g")
    } else if (ed == "trimmed_ATE") {
      for (lo in trim_levels) {
        keep <- g >= lo & g <= (1 - lo)
        wk <- ifelse(A == 1, 1 / g, 1 / (1 - g)) * keep
        pct <- 100 * mean(g[keep] < band[1] | g[keep] > band[2])
        rows[[length(rows) + 1L]] <- row_for(
          sprintf("trimmed_ATE [%.2f, %.2f]", lo, 1 - lo), wk,
          sum(!keep & A == 1), sum(!keep & A == 0), pct,
          note = "estimand changes: common-support population; g refit at estimation")
      }
    } else if (ed == "matched_ATT") {
      set.seed(1L)
      m <- .greedy_caliper_match(g, A, caliper_sd = caliper_sd)
      w <- numeric(length(A))
      w[m$treated] <- 1; w[m$control] <- 1
      rows[[length(rows) + 1L]] <- row_for(
        "matched_ATT", if (any(w > 0)) w else rep(1e-12, length(A)),
        sum(A == 1) - length(m$treated), sum(A == 0) - length(m$control), 0,
        note = sprintf("%.1f%% of treated retained",
                       100 * length(m$treated) / max(sum(A == 1), 1L)))
    }
  }
  tab <- do.call(rbind, rows)
  tab$feasible <- !tab$verdict %in% c("SEVERE", "FAIL")
  out <- list(table = tab, thresholds = thresholds, band = band,
              n = length(A), call = match.call())
  class(out) <- "estimand_feasibility"
  out
}

#' @export
print.estimand_feasibility <- function(x, ...) {
  cat("Estimand feasibility (n = ", x$n, ")\n", sep = "")
  cat(sprintf("  verdict thresholds: FLAG > %g%% outside / weight %g; SEVERE > %g%% / %g; FAIL > %g%% / %g\n\n",
              x$thresholds$flag_pct_outside, x$thresholds$flag_max_weight,
              x$thresholds$severe_pct_outside, x$thresholds$severe_max_weight,
              x$thresholds$fail_pct_outside, x$thresholds$fail_max_weight))
  print(x$table[, c("estimand", "max_weight", "ess_treated", "ess_control",
                    "n_removed_treated", "n_removed_control", "verdict",
                    "feasible")], row.names = FALSE)
  cat("\nTarget populations:\n")
  for (i in seq_len(nrow(x$table)))
    cat(sprintf("  %-24s %s%s\n", x$table$estimand[i],
                x$table$target_population[i],
                if (nzchar(x$table$note[i]))
                  paste0(" (", x$table$note[i], ")") else ""))
  invisible(x)
}


#' Who Falls Outside the Support Band?
#'
#' Standardised differences of the removed versus kept patients on a
#' user-supplied list of interpretable variables, overall and by arm, so the
#' study team can see who a trimmed analysis is no longer about. Trimming
#' changes the estimand; this profile is what makes the change concrete.
#'
#' @param ps_fit A `ps_fit`.
#' @param vars Character vector of columns of the lock data to profile.
#'   Defaults to the lock covariates. Non-numeric columns are profiled as
#'   their most common level's indicator.
#' @param band The support band. Default `c(0.05, 0.95)`.
#' @param min_group Populations with fewer removed patients than this are
#'   skipped. Default 10.
#' @return A data.frame with columns population (all, treated, control),
#'   variable, n_removed, n_kept, mean_removed, mean_kept, smd.
#' @export
who_is_unsupported <- function(ps_fit, vars = NULL, band = c(0.05, 0.95),
                               min_group = 10L) {
  if (!inherits(ps_fit, "ps_fit"))
    stop("`ps_fit` must be a ps_fit object.", call. = FALSE)
  data <- ps_fit$data
  g <- pmin(pmax(as.numeric(ps_fit$ps_raw %||% ps_fit$ps), 1e-6), 1 - 1e-6)
  A <- as.integer(data[[ps_fit$treatment]])
  removed <- g < band[1] | g > band[2]
  if (is.null(vars)) vars <- ps_fit$covariates
  vars <- intersect(vars, names(data))

  pops <- list(all = rep(TRUE, length(A)), treated = A == 1, control = A == 0)
  rows <- list()
  for (pn in names(pops)) {
    sel <- pops[[pn]]
    rem <- removed & sel; kep <- !removed & sel
    if (sum(rem) < min_group) next
    for (v in vars) {
      x <- data[[v]]
      if (!is.numeric(x)) {
        lev <- names(sort(table(x), decreasing = TRUE))[1]
        x <- as.numeric(x == lev)
        vname <- paste0(v, " = ", lev)
      } else vname <- v
      mr <- mean(x[rem], na.rm = TRUE); mk <- mean(x[kep], na.rm = TRUE)
      s  <- sqrt((stats::var(x[rem], na.rm = TRUE) +
                    stats::var(x[kep], na.rm = TRUE)) / 2)
      rows[[length(rows) + 1L]] <- data.frame(
        population = pn, variable = vname,
        n_removed = sum(rem), n_kept = sum(kep),
        mean_removed = round(mr, 4), mean_kept = round(mk, 4),
        smd = round(if (is.finite(s) && s > 0) (mr - mk) / s else 0, 4),
        stringsAsFactors = FALSE)
    }
  }
  out <- do.call(rbind, rows)
  if (is.null(out)) {
    message("who_is_unsupported: fewer than ", min_group,
            " patients outside the band in every population.")
    return(invisible(NULL))
  }
  out[order(out$population, -abs(out$smd)), , drop = FALSE]
}


#' Pre-Register the Estimand Ladder on the Lock
#'
#' Records, before outcome access, the primary estimand, the ordered
#' fallbacks, and the trigger at which the analysis moves down the ladder.
#' A support verdict at or above the trigger on the primary estimand makes
#' the switch a logged design decision rather than a silent substitution.
#' Optional sensitivity floors (an E-value floor and the additive bias that
#' moves the confidence bound to null) are recorded here too, so they are
#' declared rather than computed after unblinding.
#'
#' @param lock A `cleanroom_lock`.
#' @param primary One of `"ATE"`, `"trimmed_ATE"`, `"ATT"`, `"ATC"`,
#'   `"ATO"`, `"matched_ATT"`.
#' @param fallbacks Ordered character vector of fallbacks.
#' @param trigger The verdict that moves the analysis off the primary:
#'   `"SEVERE"` (default), `"FAIL"` (only extreme non-overlap moves it), or
#'   `"FLAG"` (any violation moves it).
#' @param evalue_floor,bias_to_null_floor Optional prespecified sensitivity
#'   floors recorded with the ladder.
#' @return The lock, with `$estimand_ladder` set and a design-log entry
#'   appended.
#' @export
declare_estimand_ladder <- function(lock,
                                    primary = "ATE",
                                    fallbacks = c("trimmed_ATE", "ATT",
                                                  "ATO"),
                                    trigger = c("SEVERE", "FAIL", "FLAG"),
                                    evalue_floor = NULL,
                                    bias_to_null_floor = NULL) {
  if (!inherits(lock, "cleanroom_lock"))
    stop("`lock` must be a cleanroom_lock object.", call. = FALSE)
  trigger <- match.arg(trigger)
  primary <- match.arg(primary, .ladder_estimands)
  fallbacks <- vapply(fallbacks, function(f)
    match.arg(f, .ladder_estimands), character(1), USE.NAMES = FALSE)
  ladder <- list(primary = primary, fallbacks = fallbacks,
                 trigger = trigger,
                 evalue_floor = evalue_floor,
                 bias_to_null_floor = bias_to_null_floor,
                 declared_at = format(Sys.time(), "%Y-%m-%d %H:%M:%S"))
  lock$estimand_ladder <- ladder
  .log_design_decision(lock, "estimand_ladder",
    sprintf("Primary %s; fallbacks %s; trigger %s%s%s.",
            primary, paste(fallbacks, collapse = " -> "), trigger,
            if (!is.null(evalue_floor))
              paste0("; E-value floor ", evalue_floor) else "",
            if (!is.null(bias_to_null_floor))
              paste0("; bias-to-null floor ", bias_to_null_floor) else ""))
}


#' Estimate the Declared Estimand Ladder
#'
#' The Stage 4 wrapper. Reads the support verdict and the feasibility table,
#' estimates the primary estimand when it is feasible under the declared
#' trigger plus every feasible fallback, and labels every returned row with
#' its estimand, its verdict, the caveat, and the implausibility flags. An
#' infeasible primary is refused unless `override_reason` is given, in which
#' case it is estimated anyway, labelled, and the override is recorded in the
#' returned object's design log.
#'
#' @param lock A `cleanroom_lock`, normally after
#'   [declare_estimand_ladder()].
#' @param ps_fit A `ps_fit` on the lock.
#' @param ladder The declared ladder; defaults to `lock$estimand_ladder`.
#' @param support A [assess_support()] result; computed if `NULL`.
#' @param feasibility An [estimand_feasibility()] result; computed if `NULL`.
#' @param family `"binomial"` or `"gaussian"`.
#' @param use_ipcw Estimate the ATE with the censoring mechanism (`Delta`)
#'   when the outcome has missing values. Default: TRUE when missing
#'   outcomes are present.
#' @param sl_library,gbound,cv_folds,prescreen_g,seed As in
#'   [run_att_tmle()].
#' @param trim_levels Passed to [run_trimmed_tmle()].
#' @param override_reason Reason for estimating an infeasible primary;
#'   recorded. Default `NULL` refuses.
#' @param allow_outcome_access Bypass the outcome guard. Default FALSE.
#' @param verbose Print the path. Default TRUE.
#' @return An object of class `estimand_ladder_result`: `table` (one row per
#'   estimated estimand), `fits` (the underlying objects), `support`,
#'   `feasibility`, `ladder`, and `design_log`.
#' @export
run_estimand_ladder <- function(lock, ps_fit,
                                ladder = NULL,
                                support = NULL,
                                feasibility = NULL,
                                family = "binomial",
                                use_ipcw = NULL,
                                sl_library = NULL,
                                gbound = NULL,
                                cv_folds = 10L,
                                prescreen_g = FALSE,
                                seed = NULL,
                                trim_levels = c(0.05, 0.10),
                                override_reason = NULL,
                                allow_outcome_access = FALSE,
                                verbose = TRUE) {
  if (!inherits(lock, "cleanroom_lock"))
    stop("`lock` must be a cleanroom_lock object.", call. = FALSE)
  if (is.null(ladder)) ladder <- lock$estimand_ladder
  if (is.null(ladder))
    stop("No estimand ladder: pass `ladder` or call ",
         "declare_estimand_ladder(lock, ...) first.", call. = FALSE)
  if (is.null(seed)) seed <- lock$seed
  if (is.null(support)) support <- assess_support(ps_fit)
  if (is.null(feasibility))
    feasibility <- estimand_feasibility(ps_fit, trim_levels = trim_levels)

  Y <- lock$data[[lock$outcome]]
  if (is.null(use_ipcw)) use_ipcw <- anyNA(Y)

  sev_rank <- c(PASS = 1L, FLAG = 2L, SEVERE = 3L, FAIL = 4L)
  trigger_rank <- sev_rank[[ladder$trigger]]
  verdict_of <- function(ed) {
    ft <- feasibility$table
    hit <- if (ed == "trimmed_ATE") startsWith(ft$estimand, "trimmed_ATE")
           else ft$estimand == ed
    if (!any(hit)) return("PASS")
    # Take the best (lowest severity) row: trimmed_ATE clears at some level.
    ft$verdict[hit][which.min(sev_rank[ft$verdict[hit]])]
  }
  feasible_of <- function(ed) sev_rank[[verdict_of(ed)]] < trigger_rank

  design_log <- lock$design_log
  note <- function(type, msg) {
    design_log <<- rbind(design_log, data.frame(
      timestamp = format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
      type = type, note = msg, stringsAsFactors = FALSE))
    if (verbose) message("[", type, "] ", msg)
  }

  todo <- character(0)
  primary_feasible <- feasible_of(ladder$primary)
  if (primary_feasible) {
    todo <- ladder$primary
  } else if (!is.null(override_reason)) {
    note("override", sprintf(
      "Primary %s is infeasible (verdict %s at trigger %s) and is estimated anyway. Reason: %s",
      ladder$primary, verdict_of(ladder$primary), ladder$trigger,
      override_reason))
    todo <- ladder$primary
  } else {
    note("estimand_switch", sprintf(
      "Primary %s infeasible (verdict %s at trigger %s); moving down the declared ladder.",
      ladder$primary, verdict_of(ladder$primary), ladder$trigger))
  }
  for (fb in ladder$fallbacks)
    if (feasible_of(fb)) todo <- c(todo, fb)
  todo <- unique(todo)
  if (!length(todo))
    stop("run_estimand_ladder: no declared estimand is feasible at trigger ",
         ladder$trigger, ".", call. = FALSE)

  fits <- list(); rows <- list()
  for (ed in todo) {
    if (verbose) message("Estimating ", ed, " ...")
    fit <- switch(ed,
      ATE = {
        args <- .tmle_delegate_args(lock, family = family,
                                    use_delta = use_ipcw,
                                    sl_library = sl_library, gbound = gbound,
                                    cv_folds = cv_folds,
                                    prescreen_g = prescreen_g)
        .check_outcome_access(lock, allow_outcome_access,
                              caller = "run_estimand_ladder")
        set.seed(seed)
        f <- do.call(tmle::tmle, args)
        est <- f$estimates$ATE
        guard <- implausibility_check(unname(est$psi),
                                      as.numeric(Y)[!is.na(Y)],
                                      as.integer(lock$data[[lock$treatment]])[
                                        !is.na(Y)], family)
        list(estimate = unname(est$psi), se = unname(sqrt(est$var.psi)),
             ci_lower = unname(est$CI[1]), ci_upper = unname(est$CI[2]),
             p_value = unname(est$pvalue),
             estimand = if (use_ipcw)
               "ATE (IPCW via Delta for missing outcomes)" else
               "ATE (complete case)",
             n = nrow(lock$data),
             risk_treated = tryCatch(unname(f$estimates$EY1$psi),
                                     error = function(e) NA_real_),
             risk_control = tryCatch(unname(f$estimates$EY0$psi),
                                     error = function(e) NA_real_),
             implausible = guard$implausible,
             implausible_reason = guard$implausible_reason,
             crude_diff = guard$crude_diff,
             fit = f)
      },
      trimmed_ATE = {
        f <- run_trimmed_tmle(lock, ps_fit, levels = trim_levels,
                              family = family, sl_library = sl_library,
                              gbound = gbound, cv_folds = cv_folds,
                              prescreen_g = prescreen_g,
                              use_ipcw = use_ipcw, seed = seed,
                              allow_outcome_access = allow_outcome_access,
                              verbose = verbose)
        c(f[c("estimate", "se", "ci_lower", "ci_upper", "p_value",
              "estimand", "n", "implausible", "implausible_reason",
              "crude_diff")], list(fit = f))
      },
      ATT = {
        f <- run_att_tmle(lock, family = family, sl_library = sl_library,
                          gbound = gbound, cv_folds = cv_folds,
                          prescreen_g = prescreen_g, seed = seed,
                          allow_outcome_access = allow_outcome_access)
        c(f[c("estimate", "se", "ci_lower", "ci_upper", "p_value",
              "estimand", "n", "implausible", "implausible_reason",
              "crude_diff")], list(fit = f))
      },
      ATO = {
        f <- estimate_ato(lock, ps_fit, sl_library = sl_library,
                          family = family,
                          allow_outcome_access = allow_outcome_access)
        c(f[c("estimate", "se", "ci_lower", "ci_upper", "p_value",
              "estimand", "n", "risk_treated", "risk_control",
              "implausible", "implausible_reason", "crude_diff")],
          list(fit = f))
      },
      ATC = stop("ATC estimation is not implemented yet; declare it only ",
                 "in the feasibility table.", call. = FALSE),
      matched_ATT = {
        g_m <- pmin(pmax(as.numeric(ps_fit$ps_raw %||% ps_fit$ps), 1e-6),
                    1 - 1e-6)
        A_m <- as.integer(lock$data[[lock$treatment]])
        set.seed(seed)
        mm <- .greedy_caliper_match(g_m, A_m)
        if (length(mm$treated) < 10L)
          stop("matched_ATT: fewer than 10 matched pairs.", call. = FALSE)
        f <- run_matched_tmle(lock, ps_fit,
                              subset_idx = sort(c(mm$treated, mm$control)),
                              override_clean_room = allow_outcome_access)
        est <- f$estimates$ATE %||% list(estimate = f$estimate, se = f$se,
                                         ci_lower = f$ci_lower,
                                         ci_upper = f$ci_upper,
                                         p_value = f$p_value)
        guard <- implausibility_check(est$estimate %||% est$psi,
                                      as.numeric(Y)[!is.na(Y)],
                                      as.integer(lock$data[[lock$treatment]])[
                                        !is.na(Y)], family)
        list(estimate = est$estimate %||% est$psi, se = est$se,
             ci_lower = est$ci_lower, ci_upper = est$ci_upper,
             p_value = est$p_value,
             estimand = "ATE on the matched cohort",
             n = f$n %||% NA_integer_,
             implausible = guard$implausible,
             implausible_reason = guard$implausible_reason,
             crude_diff = guard$crude_diff,
             fit = f)
      })
    fits[[ed]] <- fit$fit
    rows[[ed]] <- data.frame(
      estimand = ed,
      estimand_label = fit$estimand,
      is_primary = ed == ladder$primary,
      estimate = round(fit$estimate, 5), se = round(fit$se, 5),
      ci_lower = round(fit$ci_lower, 5), ci_upper = round(fit$ci_upper, 5),
      p_value = round(fit$p_value, 5), n = fit$n,
      risk_treated = round(fit$risk_treated %||% NA_real_, 5),
      risk_control = round(fit$risk_control %||% NA_real_, 5),
      support_verdict = verdict_of(ed),
      support_caveat = unname(.support_caveats[verdict_of(ed)]),
      crude_diff = round(fit$crude_diff %||% NA_real_, 5),
      implausible = isTRUE(fit$implausible),
      implausible_reason = fit$implausible_reason %||% NA_character_,
      stringsAsFactors = FALSE)
  }

  out <- list(table = do.call(rbind, rows),
              fits = fits,
              ladder = ladder,
              support = support,
              feasibility = feasibility,
              design_log = design_log,
              primary_feasible = primary_feasible,
              call = match.call())
  rownames(out$table) <- NULL
  class(out) <- "estimand_ladder_result"
  out
}

#' @export
print.estimand_ladder_result <- function(x, ...) {
  cat("Estimand ladder (primary: ", x$ladder$primary,
      if (!x$primary_feasible) " [infeasible]" else "",
      "; trigger: ", x$ladder$trigger, ")\n\n", sep = "")
  tab <- x$table
  for (i in seq_len(nrow(tab))) {
    cat(sprintf("  %s%-14s %9.4f [%9.4f, %9.4f]  p=%.4f  n=%d  %s%s\n",
                if (tab$is_primary[i]) "*" else " ",
                tab$estimand[i], tab$estimate[i], tab$ci_lower[i],
                tab$ci_upper[i], tab$p_value[i], tab$n[i],
                tab$support_verdict[i],
                if (tab$implausible[i]) "  IMPLAUSIBLE" else ""))
  }
  cat("\n  Estimand labels:\n")
  for (i in seq_len(nrow(tab)))
    cat("   ", tab$estimand[i], ": ", tab$estimand_label[i], "\n", sep = "")
  if (!is.null(x$design_log) && nrow(x$design_log)) {
    sw <- x$design_log[x$design_log$type %in% c("estimand_switch",
                                                "override"), ]
    if (nrow(sw)) {
      cat("\n  Logged decisions:\n")
      for (i in seq_len(nrow(sw)))
        cat("    [", sw$type[i], "] ", sw$note[i], "\n", sep = "")
    }
  }
  invisible(x)
}

#' @export
plot.estimand_ladder_result <- function(x, ...) {
  if (!requireNamespace("ggplot2", quietly = TRUE))
    stop("ggplot2 is required.", call. = FALSE)
  tab <- x$table
  tab$estimand <- factor(tab$estimand, levels = rev(tab$estimand))
  ggplot2::ggplot(tab, ggplot2::aes(x = .data$estimate, y = .data$estimand)) +
    ggplot2::geom_vline(xintercept = 0, linetype = "dashed",
                        colour = "grey50") +
    ggplot2::geom_errorbar(ggplot2::aes(xmin = .data$ci_lower,
                                        xmax = .data$ci_upper),
                           width = 0.15, orientation = "y") +
    ggplot2::geom_point(ggplot2::aes(shape = .data$is_primary), size = 2.6) +
    ggplot2::scale_shape_manual(values = c(`TRUE` = 18, `FALSE` = 16),
                                labels = c(`TRUE` = "primary",
                                           `FALSE` = "fallback"),
                                name = NULL) +
    ggplot2::labs(x = "estimate", y = NULL,
                  title = "Estimand ladder: each estimand answers its own question",
                  subtitle = paste("Support verdicts:",
                                   paste(tab$estimand, tab$support_verdict,
                                         sep = "=", collapse = ", "))) +
    ggplot2::theme_minimal()
}
