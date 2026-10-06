.is_text <- function(x) {
  is.character(x) && length(x) == 1L && !is.na(x) && nzchar(trimws(x))
}

# The reasons a dossier blocks unblinding without a written override. unblind()
# and .verify_unblinded() both read them from here, so they cannot drift apart.
.blocked_reasons <- function(dossier) {
  c(if (is.na(dossier$decision$primary)) "no estimand is feasible",
    if (identical(dossier$decision$nc_verdict, "STOP"))
      "the negative-control verdict is STOP")
}

.unblind_hash <- function(ub) {
  .hash(list(ub$dossier$dossier_hash, ub$approval, ub$outcomes))
}

.check_outcome_coding <- function(outcomes, plan) {
  if (plan$outcome_type == "binary") {
    y <- outcomes[[plan$outcome]]
    if (!is.numeric(y) || !all(y[!is.na(y)] %in% c(0, 1)))
      stop("The outcome `", plan$outcome, "` must be numeric and coded 0/1 (NA for ",
           "missing); found ", if (is.numeric(y)) "other values" else class(y)[1], ".",
           call. = FALSE)
    return(invisible(TRUE))
  }
  tn <- plan$outcome[["time"]]
  en <- plan$outcome[["event"]]
  tt <- outcomes[[tn]]
  ev <- outcomes[[en]]
  if (!is.numeric(tt) || any(!is.finite(tt[!is.na(tt)]) | tt[!is.na(tt)] <= 0))
    stop("The event time `", tn, "` must be numeric and greater than 0 where present.",
         call. = FALSE)
  e <- ev[!is.na(ev)]
  if (!is.numeric(ev) || any(!is.finite(e) | e < 0 | e != round(e)))
    stop("The event indicator `", en, "` must be numeric whole numbers >= 0 ",
         "(0 = censored, 1 = the event).", call. = FALSE)
  invisible(TRUE)
}

# Re-checks an unblind() result before anything is estimated or logged from
# it: the lock, the dossier, their pairing, the unblind hash (approval and
# outcomes unchanged), and the gate, recomputed from the dossier.
.verify_unblinded <- function(ub) {
  if (!inherits(ub, "cr_unblinded"))
    stop("`unblinded` must come from unblind().", call. = FALSE)
  verify_lock(ub$lock)
  .verify_dossier(ub$dossier)
  if (!identical(ub$dossier$lock_hash, ub$lock$lock_hash))
    stop("The dossier does not belong to this lock.", call. = FALSE)
  if (!identical(.unblind_hash(ub), ub$unblind_hash))
    stop("The unblind() result was modified after unblind() built it ",
         "(the approval or the outcomes changed).", call. = FALSE)
  blocked <- .blocked_reasons(ub$dossier)
  if (length(blocked) && !.is_text(ub$approval$override))
    stop("unblind refused: ", paste(blocked, collapse = "; "),
         ". This object carries no written override.", call. = FALSE)
  invisible(TRUE)
}

#' Admit the outcome after the review team approves the dossier
#'
#' The result carries `unblind_hash`, a fingerprint of the dossier hash,
#' the approval and the outcomes; [estimate_effect()] and
#' [export_design_log()] recheck it and the gate before using the object.
#'
#' @param lock A `cr_lock`.
#' @param dossier The [design_report()] for `lock`.
#' @param outcomes Data frame with the lock's id column and the plan's
#'   outcome column(s). Rows are matched by id. A binary outcome must be
#'   numeric and coded 0/1; a time-to-event outcome needs a numeric time
#'   greater than 0 and a numeric event code of whole numbers >= 0
#'   (0 = censored). Missing values are allowed.
#' @param approved_by Who approved the dossier.
#' @param override A written reason, required to proceed when no estimand is
#'   feasible or the negative-control verdict is STOP.
#' @return A `cr_unblinded`.
#' @export
unblind <- function(lock, dossier, outcomes, approved_by, override = NULL) {
  verify_lock(lock)
  .verify_dossier(dossier)
  if (!identical(dossier$lock_hash, lock$lock_hash))
    stop("This dossier was built for a different lock.", call. = FALSE)
  if (!.is_text(approved_by))
    stop("`approved_by` must name who approved the dossier.", call. = FALSE)
  oc <- .outcome_columns(lock$plan)
  miss <- setdiff(c(lock$id, oc), names(outcomes))
  if (length(miss))
    stop("`outcomes` lacks column(s): ", paste(miss, collapse = ", "), ".", call. = FALSE)
  if (anyDuplicated(outcomes[[lock$id]]))
    stop("`outcomes` has duplicated ids.", call. = FALSE)
  .check_outcome_coding(outcomes, lock$plan)
  blocked <- .blocked_reasons(dossier)
  has_override <- .is_text(override)
  if (length(blocked) && !has_override)
    stop("unblind refused: ", paste(blocked, collapse = "; "),
         ". Supply `override` with a written reason to proceed.", call. = FALSE)
  ids <- lock$data[[lock$id]]
  m <- match(ids, outcomes[[lock$id]])
  if (all(is.na(m)))
    stop("No outcome id matches a design id; check the id column's type and format.",
         call. = FALSE)
  y <- outcomes[m, oc, drop = FALSE]
  rownames(y) <- NULL
  approval <- list(approved_by = approved_by,
                   approved_at = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"),
                   override = if (has_override) override else NA_character_,
                   blocked_by = if (length(blocked)) blocked else NA_character_,
                   n_outcome_ids_not_in_design = sum(!outcomes[[lock$id]] %in% ids),
                   n_design_without_outcome = sum(is.na(m)))
  ub <- structure(list(lock = lock, dossier = dossier, outcomes = y, approval = approval),
                  class = "cr_unblinded")
  ub$unblind_hash <- .unblind_hash(ub)
  ub
}
