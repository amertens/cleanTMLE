.nc_grade <- function(tab, crit) {
  bad <- tab$status == "estimated" &
    !(is.finite(tab$estimate) & is.finite(tab$ci_lower) & is.finite(tab$ci_upper))
  tab$status[bad] <- "failed: non-finite estimate"
  est <- tab$status == "estimated"
  band <- crit$null_band
  # Each estimated control gets a reading: "inside" or "outside" the null band,
  # or (interval rule only) "inconclusive". in_band is the same reading as
  # TRUE, FALSE or NA, kept for compatibility.
  tab$reading <- NA_character_
  tab$reading[est] <- if (crit$rule == "interval") {
    lo <- tab$ci_lower[est]
    hi <- tab$ci_upper[est]
    ifelse(lo >= band[1] & hi <= band[2], "inside",
           ifelse(hi < band[1] | lo > band[2], "outside", "inconclusive"))
  } else {
    inside <- if (crit$rule == "ci") {
      tab$ci_lower[est] >= band[1] & tab$ci_upper[est] <= band[2]
    } else {
      tab$estimate[est] >= band[1] & tab$estimate[est] <= band[2]
    }
    ifelse(inside, "inside", "outside")
  }
  tab$in_band <- NA
  tab$in_band[est & tab$reading %in% "inside"] <- TRUE
  tab$in_band[est & tab$reading %in% "outside"] <- FALSE
  # Graded per rung and estimand: a control is read on each estimand separately.
  keys <- unique(tab[c("rung", "estimand", "domain")])
  by_domain <- do.call(rbind, lapply(seq_len(nrow(keys)), function(i) {
    d <- tab[tab$rung == keys$rung[i] & tab$estimand == keys$estimand[i] &
               tab$domain == keys$domain[i], , drop = FALSE]
    n_est <- sum(d$status == "estimated")
    reading <- if (any(d$reading %in% "outside")) "fail" else
      if (any(d$reading %in% "inconclusive")) "inconclusive" else
      if (n_est < crit$min_per_domain) "insufficient" else "pass"
    data.frame(rung = keys$rung[i], estimand = keys$estimand[i], domain = keys$domain[i],
               n_estimable = n_est, reading = reading, stringsAsFactors = FALSE)
  }))
  rk <- unique(tab[c("rung", "estimand")])
  by_rung <- data.frame(rung = rk$rung, estimand = rk$estimand,
                        verdict = vapply(seq_len(nrow(rk)), function(i) {
    x <- by_domain$reading[by_domain$rung == rk$rung[i] & by_domain$estimand == rk$estimand[i]]
    if (any(x == "fail")) "STOP" else if (any(x %in% c("inconclusive", "insufficient"))) "FLAG" else "GO"
  }, character(1)), stringsAsFactors = FALSE, row.names = NULL)
  list(table = tab, by_domain = by_domain, by_rung = by_rung)
}

#' Negative controls along the restriction ladder
#'
#' Estimates each declared negative-control outcome, for every declared
#' estimand, on the full cohort and on each declared restriction, with the
#' main-terms GLM candidate (the smallest declared truncation) and the same
#' fitter as the analysis, and grades the result against the plan's
#' `nc_criteria`. A control is read on the estimand the analysis reports: a
#' control read on an estimand the design cannot deliver mixes extrapolation
#' error with confounding, so it cannot isolate confounding.
#'
#' Each estimated control gets a `reading` from the plan's rule. Under the
#' default `"interval"` rule it reads "inside" when its 95% interval lies
#' entirely inside `null_band`, "outside" when the interval lies entirely
#' outside it, and "inconclusive" otherwise (too imprecise to tell). Under
#' `"point"` it reads "inside" or "outside" by the point estimate, and under
#' `"ci"` it reads "inside" only when the whole interval is inside the band.
#' Neither of those two rules gives "inconclusive".
#'
#' Grading is per rung and estimand. Within a rung and estimand, a domain
#' reads "fail" if any control is outside the band (even when too few controls
#' are estimable), otherwise "inconclusive" if any control is inconclusive,
#' otherwise "insufficient" if fewer than `min_per_domain` controls are
#' estimable, otherwise "pass". A rung is STOP if any domain fails, FLAG if
#' any is inconclusive or insufficient, and GO otherwise. FLAG does not block
#' [unblind()]; only STOP does. An estimated control whose estimate or
#' interval is not finite counts as failed, not as estimable.
#'
#' A control is inestimable on a rung when a cell of its 2x2 table with the
#' treatment has fewer than 5 events. That check uses the whole rung, so the
#' trimmed ATE (fitted on the patients inside the trim band) can rest on a
#' subset with smaller cells; its interval is then wide.
#'
#' `verdict` is a named character vector: the full cohort's verdict for each
#' declared estimand, in ladder order. [design_report()] gates unblinding on
#' the full cohort's verdict for the estimand the dossier selects (the
#' primary), because [estimate_effect()] analyses the full cohort. The other
#' estimands' verdicts and the restricted rungs are diagnostics: they show
#' whether a restriction would remove the confounding.
#'
#' @param lock A `cr_lock`.
#' @param design The [assess_design()] result for `lock`.
#' @return A `cr_nc_ladder`.
#' @export
negative_control_ladder <- function(lock, design) {
  verify_lock(lock)
  .check_stamp(lock, design, "design")
  p <- lock$plan
  out <- list(declared = !is.null(p$negative_controls), table = NULL, by_domain = NULL,
              by_rung = NULL, verdict = NA_character_, criteria = p$nc_criteria,
              design_hash = design$stamp$hash)
  if (out$declared) {
    X <- .design_matrix(lock$data, lock$covariates)
    A <- as.integer(lock$data[[lock$treatment]])
    n <- nrow(X)
    rungs <- c(list(`full cohort` = rep(TRUE, n)), lapply(p$restrictions, function(f) {
      r <- eval(f[[2]], lock$data, environment(f))
      r & !is.na(r)
    }))
    cand <- list(id = "glm_nc", library = "glm", learners = "glm",
                 truncation = min(p$candidates$truncation))
    rows <- list()
    row <- function(rung, est, nc, n, status, r = NULL) data.frame(
      rung = rung, estimand = est, negative_control = nc,
      domain = unname(p$negative_controls[nc]), n = n,
      estimate = r$estimate %||% NA_real_, ci_lower = r$ci_lower %||% NA_real_,
      ci_upper = r$ci_upper %||% NA_real_, status = status, stringsAsFactors = FALSE)
    add_all <- function(rn, nc, n, status) {
      for (e in p$estimands) rows[[length(rows) + 1L]] <<- row(rn, e, nc, n, status)
    }
    for (rn in names(rungs)) {
      keep <- rungs[[rn]]
      Xr <- X[keep, , drop = FALSE]
      Ar <- A[keep]
      folds <- .make_folds(sum(keep), p$K, p$seed)
      g_fit <- tryCatch(.fit_g(Xr, Ar, p$ps_library, folds, p$V, p$seed),
                        error = function(e) conditionMessage(e))
      for (nc in names(p$negative_controls)) {
        y <- lock$data[[nc]][keep]
        ok <- !is.na(y)
        cells <- c(sum(y[ok & Ar == 1]), sum(1 - y[ok & Ar == 1]),
                   sum(y[ok & Ar == 0]), sum(1 - y[ok & Ar == 0]))
        if (any(cells < 5)) {
          add_all(rn, nc, sum(ok), "inestimable (a cell below 5 events)")
          next
        }
        if (is.character(g_fit)) {
          add_all(rn, nc, sum(ok), paste("failed:", g_fit))
          next
        }
        # One call fits the outcome model once and targets every declared
        # estimand on it; each estimand is read from its own row.
        r <- tryCatch(fit_candidate(Xr[ok, , drop = FALSE], Ar[ok], y[ok], p$estimands, cand, p,
                                    folds[ok], p$seed,
                                    g_fit = list(g = g_fit$g[ok], risks = g_fit$risks)),
                      error = function(e) conditionMessage(e))
        if (is.character(r)) {
          add_all(rn, nc, sum(ok), paste("failed:", r))
          next
        }
        for (e in p$estimands) {
          re <- r[r$estimand == e, , drop = FALSE][1, , drop = FALSE]
          rows[[length(rows) + 1L]] <- if (isTRUE(re$failed)) {
            row(rn, e, nc, sum(ok), paste("failed:", re$message))
          } else row(rn, e, nc, sum(ok), "estimated", re)
        }
      }
    }
    tab <- do.call(rbind, rows)
    # Rung, then estimand in ladder order, then control.
    tab <- tab[order(match(tab$rung, names(rungs)), match(tab$estimand, p$estimands)), ,
               drop = FALSE]
    rownames(tab) <- NULL
    gr <- .nc_grade(tab, p$nc_criteria)
    out$table <- gr$table
    out$by_domain <- gr$by_domain
    out$by_rung <- gr$by_rung
    # estimate_effect() analyses the full cohort, so its rung gates unblinding
    # (for the estimand design_report() selects); the restricted rungs are
    # diagnostics.
    full <- gr$by_rung[gr$by_rung$rung == "full cohort", , drop = FALSE]
    out$verdict <- stats::setNames(full$verdict[match(p$estimands, full$estimand)],
                                   p$estimands)
  }
  out$stamp <- .stamp(lock, out)
  class(out) <- "cr_nc_ladder"
  out
}
