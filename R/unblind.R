#' Admit the outcome after the review team approves the dossier
#'
#' @param lock A `cr_lock`.
#' @param dossier The [design_report()] for `lock`.
#' @param outcomes Data frame with the lock's id column and the plan's
#'   outcome column(s). Rows are matched by id.
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
  if (!is.character(approved_by) || length(approved_by) != 1L || !nzchar(approved_by))
    stop("`approved_by` must name who approved the dossier.", call. = FALSE)
  oc <- .outcome_columns(lock$plan)
  miss <- setdiff(c(lock$id, oc), names(outcomes))
  if (length(miss))
    stop("`outcomes` lacks column(s): ", paste(miss, collapse = ", "), ".", call. = FALSE)
  if (anyDuplicated(outcomes[[lock$id]]))
    stop("`outcomes` has duplicated ids.", call. = FALSE)
  blocked <- c(if (is.na(dossier$decision$primary)) "no estimand is feasible",
               if (identical(dossier$decision$nc_verdict, "STOP"))
                 "the negative-control verdict is STOP")
  has_override <- is.character(override) && length(override) == 1L && nzchar(override)
  if (length(blocked) && !has_override)
    stop("unblind refused: ", paste(blocked, collapse = "; "),
         ". Supply `override` with a written reason to proceed.", call. = FALSE)
  ids <- lock$data[[lock$id]]
  m <- match(ids, outcomes[[lock$id]])
  y <- outcomes[m, oc, drop = FALSE]
  rownames(y) <- NULL
  approval <- list(approved_by = approved_by,
                   approved_at = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"),
                   override = if (has_override) override else NA_character_,
                   blocked_by = if (length(blocked)) blocked else NA_character_,
                   n_outcome_ids_not_in_design = sum(!outcomes[[lock$id]] %in% ids),
                   n_design_without_outcome = sum(is.na(m)))
  structure(list(lock = lock, dossier = dossier, outcomes = y, approval = approval),
            class = "cr_unblinded")
}
