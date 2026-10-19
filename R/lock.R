.pkg_version <- function() as.character(utils::packageVersion("cleanTMLE"))

# A negative control that is also an adjustment covariate is regressed on
# itself, so its estimate is zero by construction. Checked by exact name when
# the lock is made, and again by the ladder and the dossier, because a lock
# saved before the check existed, or re-hashed by hand, still verifies.
.check_nc_not_covariate <- function(controls, covariates, advice) {
  same <- intersect(controls, covariates)
  if (length(same))
    stop("Negative control `", same[1], "` is also an adjustment covariate, so its ",
         "estimate is zero by construction. ", advice, call. = FALSE)
  invisible(TRUE)
}

# The dictionary enters the hash only when the lock has one, so a lock made by
# cleanTMLE 0.4.x, which has none, still verifies.
.lock_hash <- function(lock) {
  x <- list(plan = lock$plan, data = lock$data, id = lock$id,
            treatment = lock$treatment, covariates = lock$covariates,
            seed = lock$plan$seed, K = lock$plan$K, version = lock$version)
  if (!is.null(lock$dictionary)) x$dictionary <- lock$dictionary
  .hash(x)
}

#' Lock the plan and the outcome-free design data
#'
#' @param design_data Data frame with the id, treatment, covariates, and
#'   any negative-control outcomes and restriction variables. It must not
#'   contain the plan's outcome columns.
#' @param treatment Name of the 0/1 treatment column.
#' @param covariates Names of the adjustment covariates (complete, no NA). A
#'   negative-control column must not be one of them; a control derived from a
#'   covariate under another name is not detected, so list its source in
#'   `nc_exclude` of [analysis_plan()].
#' @param plan A [analysis_plan()].
#' @param id Name of the unique row id, used later to join the outcomes.
#' @param dictionary Required. A data frame with one row per column the lock
#'   keeps: the id, the treatment, the covariates, the negative controls, the
#'   restriction variables, and any column you declare `after_treatment` or
#'   `observation`. A dummy-coded variable has one row per dummy. Its fields are
#'   \describe{
#'   \item{`column`}{A column of `design_data`.}
#'   \item{`variable`}{The variable the column comes from, so the dummies of one
#'     variable share a value (and must share their timing).}
#'   \item{`source`}{The form or file it was recorded on.}
#'   \item{`follow_up`}{`TRUE` if it was recorded only for patients reached at
#'     follow-up.}
#'   \item{`recorded`}{`"before"`, `"at"` or `"after"` time zero: when the value
#'     was recorded.}
#'   \item{`fixed`}{`"before"`, `"at"` or `"after"` time zero: when the attribute
#'     it describes was fixed. Age is recorded after injury but fixed before it;
#'     heart rate on arrival is fixed after.}
#'   \item{`role`}{`"id"`, `"treatment"`, `"covariate"`, `"stand_in"`,
#'     `"after_treatment"`, `"observation"`, `"negative_control"` or
#'     `"restriction"`. A column that is both a restriction variable and a
#'     covariate takes the covariate role.}
#'   \item{`reason`, `sensitivity`}{Text, may be empty. The declared responses
#'     described below.}
#'   }
#'   The lock stops, whatever is declared, when a kept column has no row or a row
#'   names a column that is absent, a role does not match how the column is used,
#'   a covariate is fixed after time zero, a column from the follow-up source is
#'   neither an outcome-observation indicator nor a negative control, a node after
#'   treatment (`after_treatment` or `observation`) is in `covariates`, or the
#'   treatment is not fixed at time zero. It stops unless a response is given for
#'   a `stand_in` (a `reason` and a `sensitivity` analysis that leaves it out), an
#'   `observation` column (a `reason`, such as the approval), a negative control
#'   from the follow-up source (a `reason`), and a variable the plan restricts on
#'   that is fixed after time zero, whatever its role (a `reason`). The lock checks that each response is
#'   text; it cannot judge whether it is adequate, and it cannot check that the declared
#'   timings are true. The dictionary is stored in the lock and hashed with it.
#' @return A `cr_lock`.
#' @export
create_analysis_lock <- function(design_data, treatment, covariates, plan,
                                 id = "id", dictionary) {
  if (!is.data.frame(design_data))
    stop("`design_data` must be a data frame.", call. = FALSE)
  if (!inherits(plan, "cr_plan"))
    stop("`plan` must come from analysis_plan().", call. = FALSE)
  if (missing(dictionary) || is.null(dictionary))
    stop("`create_analysis_lock()` needs a variable `dictionary`: a data frame with one row ",
         "per column the lock keeps, saying where and when each was recorded and when the ",
         "attribute it describes was fixed, relative to time zero. See ?create_analysis_lock ",
         "and `cr_example$dictionary`.", call. = FALSE)
  oc <- intersect(.outcome_columns(plan), names(design_data))
  if (length(oc))
    stop("`design_data` contains the outcome column(s) ", paste(oc, collapse = ", "),
         ". Remove them: the outcome enters only through unblind().", call. = FALSE)
  miss <- setdiff(c(id, treatment, covariates), names(design_data))
  if (length(miss))
    stop("Missing column(s): ", paste(miss, collapse = ", "), ".", call. = FALSE)
  if (anyNA(design_data[[id]]) || anyDuplicated(design_data[[id]]))
    stop("`id` must be unique and complete.", call. = FALSE)
  a <- design_data[[treatment]]
  if (anyNA(a) || !all(a %in% c(0, 1)) || length(unique(a)) < 2L)
    stop("`treatment` must be coded 0/1 with both arms present and no missing values.",
         call. = FALSE)
  if (anyNA(design_data[covariates]))
    stop("Covariates must be complete; impute before locking.", call. = FALSE)
  nc_cols <- names(plan$negative_controls)
  r_vars <- unique(unlist(lapply(plan$restrictions, all.vars)))
  miss <- setdiff(c(nc_cols, r_vars), names(design_data))
  if (length(miss))
    stop("Negative-control or restriction column(s) not in `design_data`: ",
         paste(miss, collapse = ", "), ".", call. = FALSE)
  bad_nc <- nc_cols[!vapply(nc_cols, function(v) {
    x <- design_data[[v]]
    (is.numeric(x) || is.integer(x)) && all(x[!is.na(x)] %in% c(0, 1))
  }, logical(1))]
  if (length(bad_nc))
    stop("Negative-control column(s) must be coded 0/1: ",
         paste(bad_nc, collapse = ", "), ".", call. = FALSE)
  # A control that is a covariate is regressed on itself, so its estimate is zero
  # by construction and the ladder would read GO. The check is by exact name.
  .check_nc_not_covariate(nc_cols, covariates, paste0(
    "Remove it from `covariates` or from `negative_controls`. A control derived from a ",
    "covariate under another name is not caught here; list that covariate in `nc_exclude`."))
  for (nc in names(plan$nc_exclude)) {
    ex <- plan$nc_exclude[[nc]]
    if (length(setdiff(ex, covariates)))
      stop("`nc_exclude` for `", nc, "` names `", setdiff(ex, covariates)[1],
           "`, which is not an adjustment covariate.", call. = FALSE)
    if (!length(setdiff(covariates, ex)))
      stop("`nc_exclude` for `", nc, "` leaves out every covariate, so no adjustment ",
           "set remains.", call. = FALSE)
  }
  if (!is.null(plan$surfaces$q0) && length(plan$surfaces$q0) != nrow(design_data))
    stop("`surfaces$q0` must have one value per row of `design_data`.", call. = FALSE)

  dictionary <- .check_dictionary(dictionary, design_data, id, treatment, covariates, plan)

  # Nodes after treatment and outcome-observation indicators are kept, out of the
  # adjustment set, for the analysis of how the outcome is observed.
  keep <- .kept_columns(id, treatment, covariates, plan, dictionary)
  data <- design_data[, keep, drop = FALSE]
  rownames(data) <- NULL
  lock <- list(data = data, id = id, treatment = treatment,
               covariates = covariates, plan = plan, dictionary = dictionary,
               version = .pkg_version())
  lock$lock_hash <- .lock_hash(lock)
  # The creation time is recorded but not hashed, so the hash stays
  # reproducible from its inputs; the dossier copies it under its own hash.
  lock$created <- .utc_now()
  class(lock) <- "cr_lock"
  lock
}

#' Verify a lock's hash
#'
#' @param lock A `cr_lock`.
#' @return `TRUE`, invisibly; an error if anything locked has changed.
#' @export
verify_lock <- function(lock) {
  if (!inherits(lock, "cr_lock"))
    stop("`lock` must come from create_analysis_lock().", call. = FALSE)
  if (!identical(.lock_hash(lock), lock$lock_hash))
    stop("Lock hash mismatch: the plan, the design data or another locked field ",
         "changed after the lock was created.", call. = FALSE)
  if (!identical(lock$version, .pkg_version()))
    warning("This lock was created with cleanTMLE ", lock$version, "; you are running ",
            .pkg_version(), ". Results may differ from the version that created it.",
            call. = FALSE)
  invisible(TRUE)
}

#' Save and load a lock
#'
#' Both verify the hash.
#' @param lock A `cr_lock`.
#' @param path File path (`.rds`).
#' @return `save_lock()` returns `path` invisibly; `load_lock()` the lock.
#' @export
save_lock <- function(lock, path) {
  verify_lock(lock)
  saveRDS(lock, path)
  invisible(path)
}

#' @rdname save_lock
#' @export
load_lock <- function(path) {
  lock <- readRDS(path)
  verify_lock(lock)
  lock
}

.stamp <- function(lock, obj) {
  obj$stamp <- NULL
  list(lock_hash = lock$lock_hash, hash = .hash(unclass(obj)))
}

.check_stamp <- function(lock, obj, what) {
  if (!identical(.lock_hash(lock), lock$lock_hash))
    stop("Lock hash mismatch: the plan, the design data or another locked field ",
         "changed after the lock was created.", call. = FALSE)
  if (is.null(obj$stamp) || !identical(obj$stamp$lock_hash, lock$lock_hash))
    stop("`", what, "` was not computed from this lock.", call. = FALSE)
  if (!identical(.stamp(lock, obj)$hash, obj$stamp$hash))
    stop("`", what, "` was modified after it was computed.", call. = FALSE)
  invisible(TRUE)
}
