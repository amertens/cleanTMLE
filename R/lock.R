.pkg_version <- function() as.character(utils::packageVersion("cleanTMLE"))

.lock_hash <- function(lock) {
  .hash(list(plan = lock$plan, data = lock$data, id = lock$id,
             treatment = lock$treatment, covariates = lock$covariates,
             seed = lock$plan$seed, K = lock$plan$K, version = lock$version))
}

#' Lock the plan and the outcome-free design data
#'
#' @param design_data Data frame with the id, treatment, covariates, and
#'   any negative-control outcomes and restriction variables. It must not
#'   contain the plan's outcome columns.
#' @param treatment Name of the 0/1 treatment column.
#' @param covariates Names of the adjustment covariates (complete, no NA).
#' @param plan A [analysis_plan()].
#' @param id Name of the unique row id, used later to join the outcomes.
#' @return A `cr_lock`.
#' @export
create_analysis_lock <- function(design_data, treatment, covariates, plan,
                                 id = "id") {
  if (!is.data.frame(design_data))
    stop("`design_data` must be a data frame.", call. = FALSE)
  if (!inherits(plan, "cr_plan"))
    stop("`plan` must come from analysis_plan().", call. = FALSE)
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
  if (!is.null(plan$surfaces$q0) && length(plan$surfaces$q0) != nrow(design_data))
    stop("`surfaces$q0` must have one value per row of `design_data`.", call. = FALSE)

  keep <- unique(c(id, treatment, covariates, nc_cols, r_vars))
  data <- design_data[, keep, drop = FALSE]
  rownames(data) <- NULL
  lock <- list(data = data, id = id, treatment = treatment,
               covariates = covariates, plan = plan, version = .pkg_version())
  lock$lock_hash <- .lock_hash(lock)
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
  if (is.null(obj$stamp) || !identical(obj$stamp$lock_hash, lock$lock_hash))
    stop("`", what, "` was not computed from this lock.", call. = FALSE)
  if (!identical(.stamp(lock, obj)$hash, obj$stamp$hash))
    stop("`", what, "` was modified after it was computed.", call. = FALSE)
  invisible(TRUE)
}
