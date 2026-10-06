.nc_grade <- function(tab, crit) {
  bad <- tab$status == "estimated" &
    !(is.finite(tab$estimate) & is.finite(tab$ci_lower) & is.finite(tab$ci_upper))
  tab$status[bad] <- "failed: non-finite estimate"
  est <- tab$status == "estimated"
  band <- crit$null_band
  tab$in_band <- NA
  tab$in_band[est] <- if (crit$rule == "ci") {
    tab$ci_lower[est] >= band[1] & tab$ci_upper[est] <= band[2]
  } else {
    tab$estimate[est] >= band[1] & tab$estimate[est] <= band[2]
  }
  keys <- unique(tab[c("rung", "domain")])
  by_domain <- do.call(rbind, lapply(seq_len(nrow(keys)), function(i) {
    d <- tab[tab$rung == keys$rung[i] & tab$domain == keys$domain[i], , drop = FALSE]
    n_est <- sum(d$status == "estimated")
    reading <- if (any(!d$in_band, na.rm = TRUE)) "fail" else
      if (n_est < crit$min_per_domain) "insufficient" else "pass"
    data.frame(rung = keys$rung[i], domain = keys$domain[i], n_estimable = n_est,
               reading = reading, stringsAsFactors = FALSE)
  }))
  rungs <- unique(tab$rung)
  by_rung <- data.frame(rung = rungs, verdict = vapply(rungs, function(r) {
    x <- by_domain$reading[by_domain$rung == r]
    if (any(x == "fail")) "STOP" else if (any(x == "insufficient")) "FLAG" else "GO"
  }, character(1)), stringsAsFactors = FALSE, row.names = NULL)
  list(table = tab, by_domain = by_domain, by_rung = by_rung)
}

#' Negative controls along the restriction ladder
#'
#' Estimates each declared negative-control outcome on the full cohort and
#' on each declared restriction, with the main-terms GLM candidate and the
#' same fitter as the analysis, and grades the result against the plan's
#' `nc_criteria`. Within a rung, a domain reads "fail" if any estimable control is out of
#' the band (even when too few controls are estimable), otherwise "insufficient" if fewer than
#' `min_per_domain` controls are estimable, otherwise "pass". A rung is STOP if any domain fails,
#' FLAG if any is insufficient, and GO otherwise. An estimated control whose estimate or
#' interval is not finite counts as failed, not as estimable. The verdict that gates
#' unblinding is the full cohort's; restricted rungs show whether a restriction would
#' remove the confounding.
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
    row <- function(rung, nc, n, status, r = NULL) data.frame(
      rung = rung, negative_control = nc, domain = unname(p$negative_controls[nc]), n = n,
      estimate = r$estimate %||% NA_real_, ci_lower = r$ci_lower %||% NA_real_,
      ci_upper = r$ci_upper %||% NA_real_, status = status, stringsAsFactors = FALSE)
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
          rows[[length(rows) + 1L]] <- row(rn, nc, sum(ok), "inestimable (a cell below 5 events)")
          next
        }
        if (is.character(g_fit)) {
          rows[[length(rows) + 1L]] <- row(rn, nc, sum(ok), paste("failed:", g_fit))
          next
        }
        r <- tryCatch(fit_candidate(Xr[ok, , drop = FALSE], Ar[ok], y[ok], "ATE", cand, p,
                                    folds[ok], p$seed,
                                    g_fit = list(g = g_fit$g[ok], risks = g_fit$risks)),
                      error = function(e) conditionMessage(e))
        rows[[length(rows) + 1L]] <- if (is.character(r) || isTRUE(r$failed)) {
          row(rn, nc, sum(ok), paste("failed:", if (is.character(r)) r else r$message))
        } else row(rn, nc, sum(ok), "estimated", r)
      }
    }
    gr <- .nc_grade(do.call(rbind, rows), p$nc_criteria)
    out$table <- gr$table
    out$by_domain <- gr$by_domain
    out$by_rung <- gr$by_rung
    # estimate_effect() analyses the full cohort, so its rung gates unblinding;
    # the restricted rungs are diagnostics.
    out$verdict <- gr$by_rung$verdict[gr$by_rung$rung == "full cohort"][1]
  }
  out$stamp <- .stamp(lock, out)
  class(out) <- "cr_nc_ladder"
  out
}
