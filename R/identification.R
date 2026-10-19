# The identification table: what the plan and the variable dictionary declare
# about time zero, eligibility and the timing of every locked column, with the
# response given to each flag. It lists what was declared and flagged; it does
# not grade identification. A feasible verdict says an estimator can recover the
# estimand if identification holds; this table speaks to whether it does.

.identification_components <- c(
  "Time zero", "Start of follow-up", "Eligibility", "Treatment",
  "Adjustment set: fixed before time zero", "Adjustment set: stand-in",
  "Node after treatment", "Outcome-observation indicator",
  "Restriction variable", "Negative control")

.identification_table <- function(lock) {
  tz <- lock$plan$time_zero
  el <- lock$plan$eligibility
  dict <- lock$dictionary
  txt <- function(x) if (is.null(x) || length(x) == 0L || is.na(x)) "" else as.character(x)
  rows <- list()
  add <- function(component, item, declared, flag = "", response = "")
    rows[[length(rows) + 1L]] <<- data.frame(
      component = component, item = txt(item), declared = txt(declared),
      flag = flag, response = txt(response), stringsAsFactors = FALSE)

  # One row per variable: the first column's timing (the columns of a variable
  # agree), and the first reason and sensitivity given in any of its columns.
  per_variable <- function(rows_of_role) {
    d <- dict[rows_of_role, , drop = FALSE]
    lapply(unique(d$variable), function(v) {
      g <- d[d$variable == v, , drop = FALSE]
      first <- function(x) { x <- x[nzchar(trimws(x))]; if (length(x)) x[1] else "" }
      list(variable = v, source = g$source[1], recorded = g$recorded[1],
           fixed = g$fixed[1], reason = first(g$reason), sensitivity = first(g$sensitivity))
    })
  }
  timing <- function(x) sprintf("recorded %s time zero, fixed %s time zero (%s)",
                                x$recorded, x$fixed, x$source)

  add("Time zero", "time zero",
      paste0(tz$event, if (!is.null(tz$measured_as))
        paste0(" (measured as ", tz$measured_as, ")") else ""))
  add("Start of follow-up", "follow-up", tz$follow_up_from,
      if (!identical(tz$follow_up_from, tz$event)) "flag" else "", tz$response)
  for (nm in names(el)) {
    e <- el[[nm]]
    add("Eligibility", nm, paste0(e$description, " (assessed ", e$assessed, " time zero)"),
        if (identical(e$assessed, "after")) "flag" else "", e$response)
  }
  tr <- dict[dict$role == "treatment", , drop = FALSE]
  add("Treatment", tr$column[1],
      sprintf("recorded %s time zero (%s); fixed %s time zero", tr$recorded[1], tr$source[1],
              tr$fixed[1]),
      if (identical(tr$recorded[1], "after")) "note" else "")
  cv <- per_variable(dict$role == "covariate")
  late <- vapply(cv, function(x) identical(x$recorded, "after"), logical(1))
  late_vars <- vapply(cv[late], function(x) x$variable, character(1))
  late_reasons <- vapply(cv[late], function(x)
    if (nzchar(x$reason)) paste0(x$variable, ": ", x$reason) else "", character(1))
  add("Adjustment set: fixed before time zero", sprintf("%d variables", length(cv)),
      if (length(late_vars)) paste(late_vars, collapse = ", ")
      else "all recorded at or before time zero",
      if (length(late_vars)) "note" else "",
      paste(late_reasons[nzchar(late_reasons)], collapse = "; "))
  for (x in per_variable(dict$role == "stand_in"))
    add("Adjustment set: stand-in", x$variable, timing(x), "flag",
        paste0(x$reason, "; sensitivity: ", x$sensitivity))
  for (x in per_variable(dict$role == "after_treatment"))
    add("Node after treatment", x$variable, timing(x), "note", x$reason)
  for (x in per_variable(dict$role == "observation"))
    add("Outcome-observation indicator", x$variable, timing(x), "flag", x$reason)
  # Every variable the plan restricts on and that was fixed after time zero, whatever
  # its role; a variable with another role also appears under that role.
  r_vars <- unique(unlist(lapply(lock$plan$restrictions, all.vars)))
  for (x in per_variable(dict$column %in% r_vars & dict$fixed == "after"))
    add("Restriction variable", x$variable, timing(x), "flag", x$reason)
  for (x in per_variable(dict$role == "negative_control"))
    add("Negative control", x$variable, timing(x),
        if (identical(x$fixed, "after")) "note" else "", x$reason)
  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  out
}
