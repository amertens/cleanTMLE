.deparse1 <- function(x) paste(deparse(x), collapse = " ")

.plan_summary <- function(plan) {
  p <- unclass(plan)
  if (!is.null(p$select)) p$select <- .deparse1(p$select)
  if (!is.null(p$restrictions)) p$restrictions <- as.list(vapply(p$restrictions, .deparse1, ""))
  if (!is.null(p$surfaces$q0))
    p$surfaces$q0 <- sprintf("steward-supplied (%d values)", length(p$surfaces$q0))
  if (!is.null(p$surfaces$custom)) p$surfaces$custom <- as.list(vapply(p$surfaces$custom, .deparse1, ""))
  if (!is.null(p$hazards))
    p$hazards <- lapply(p$hazards, function(h) if (is.list(h)) vapply(h, .deparse1, "") else .deparse1(h))
  if (!is.null(p$negative_controls)) p$negative_controls <- as.list(p$negative_controls)
  p
}

.dossier_hash <- function(d) {
  d$dossier_hash <- NULL
  .hash(unclass(d))
}

.verify_dossier <- function(d) {
  if (!inherits(d, "cr_dossier"))
    stop("`dossier` must come from design_report().", call. = FALSE)
  if (!identical(.dossier_hash(d), d$dossier_hash))
    stop("The dossier was modified after design_report() built it.", call. = FALSE)
  invisible(TRUE)
}

#' Build the design dossier: the decision, before any outcome is seen
#'
#' Checks that every design-stage object belongs to the lock, walks the
#' estimand ladder (the primary is the first declared estimand whose status
#' is feasible), records the selected candidate for each feasible estimand,
#' flags each estimand whose decision (status or selected candidate) rests on
#' cells still unresolved after their final repetitions
#' (`decision$unresolved`), records the negative-control verdict, and
#' fingerprints the result.
#'
#' @param lock A `cr_lock`.
#' @param design,simulation,nc Results of [assess_design()],
#'   [simulate_design()] and [negative_control_ladder()] for `lock`.
#' @param file Optional `.html` or `.docx` path to render with Quarto.
#' @return A `cr_dossier`.
#' @export
design_report <- function(lock, design, simulation, nc = NULL, file = NULL) {
  verify_lock(lock)
  .check_stamp(lock, design, "design")
  .check_stamp(lock, simulation, "simulation")
  if (!identical(simulation$design_hash, design$stamp$hash))
    stop("`simulation` was not run on this `design`.", call. = FALSE)
  if (!is.null(nc)) {
    .check_stamp(lock, nc, "nc")
    if (!identical(nc$design_hash, design$stamp$hash))
      stop("`nc` was not run on this `design`.", call. = FALSE)
  }
  if (!is.null(lock$plan$negative_controls) && is.null(nc))
    stop("The plan declares negative controls, so `nc` (from negative_control_ladder()) ",
         "is required; the dossier cannot record a verdict without it.", call. = FALSE)
  p <- lock$plan
  v <- simulation$verdict
  feas <- v[v$status == "feasible", , drop = FALSE]
  primary <- p$estimands[p$estimands %in% feas$estimand][1]
  candidates <- stats::setNames(lapply(seq_len(nrow(feas)), function(i) list(
    id = feas$selected[i], library = feas$library[i],
    learners = p$candidates$library[[feas$library[i]]], truncation = feas$truncation[i])),
    feas$estimand)
  decision <- list(primary = primary, ladder = p$estimands,
                   status = as.list(stats::setNames(v$status, v$estimand)),
                   # TRUE when the estimand's decision (status or selected
                   # candidate) rests on cells still unresolved after their
                   # final repetitions; they were decided on point estimates.
                   unresolved = as.list(stats::setNames(v$unresolved, v$estimand)),
                   candidates = candidates,
                   nc_verdict = if (is.null(nc)) NA_character_ else nc$verdict,
                   K = p$K, fold_seed = p$seed)
  notes <- c("Weighting and matching are comparators and were not stress-tested.",
             simulation$scope)
  d <- list(
    lock_hash = lock$lock_hash, plan = .plan_summary(p), outcome_type = p$outcome_type,
    design = design[c("overlap", "overlap_grade", "share_outside_band", "band", "ess",
                      "balance", "edge")],
    simulation = simulation[c("metrics", "verdict", "truths", "scope", "decision_note",
                              "extension", "reps", "max_reps", "K", "edge")],
    nc = if (is.null(nc)) NULL else nc[c("table", "by_domain", "by_rung", "verdict")],
    decision = decision, notes = notes,
    component_hashes = list(design = design$stamp$hash, simulation = simulation$stamp$hash,
                            nc = if (is.null(nc)) NA_character_ else nc$stamp$hash),
    version = .pkg_version(), created = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"))
  class(d) <- "cr_dossier"
  d$dossier_hash <- .dossier_hash(d)
  if (!is.null(file)) .render_dossier(d, file)
  d
}

.render_dossier <- function(d, file) {
  if (!requireNamespace("quarto", quietly = TRUE))
    stop("Rendering the dossier needs the quarto package.", call. = FALSE)
  # Resolve the output path against the caller's working directory before
  # rendering changes it.
  file <- normalizePath(file, winslash = "/", mustWork = FALSE)
  work <- tempfile("dossier")
  dir.create(work)
  if (!isTRUE(file.copy(system.file("templates", "dossier.qmd", package = "cleanTMLE"),
                        work)))
    stop("Could not copy the dossier template.", call. = FALSE)
  rds <- normalizePath(file.path(work, "dossier.rds"), winslash = "/", mustWork = FALSE)
  saveRDS(d, rds)
  fmt <- if (grepl("\\.docx$", file, ignore.case = TRUE)) "docx" else "html"
  # Render from inside the work directory with a bare input name: given an
  # absolute input path, Quarto can hand its R process a relative path that
  # does not resolve (seen with a deep working directory and a nested TMPDIR).
  withr::with_dir(work, quarto::quarto_render(
    "dossier.qmd", output_format = fmt,
    execute_params = list(dossier = rds), quiet = TRUE))
  out <- file.path(work, paste0("dossier.", fmt))
  if (!file.exists(out))
    stop("Quarto did not produce the rendered dossier.", call. = FALSE)
  if (!isTRUE(file.copy(out, file, overwrite = TRUE)))
    stop("Could not write the rendered dossier to ", file, ".", call. = FALSE)
  invisible(file)
}
