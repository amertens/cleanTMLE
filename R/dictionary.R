# The variable dictionary: one row per column the lock keeps, saying where the
# column was recorded, when it was recorded and when the attribute it describes
# was fixed (both relative to time zero), and how the analysis uses it. The lock
# refuses combinations that cannot identify a total effect and requires a declared
# response to every flag. It cannot check that the declared timings are true; it
# makes the analyst state them, in the hashed record, before the outcome.

.dictionary_fields <- c("column", "variable", "source", "follow_up", "recorded",
                        "fixed", "role", "reason", "sensitivity")
.dictionary_timing <- c("before", "at", "after")
.dictionary_roles <- c("id", "treatment", "covariate", "stand_in", "after_treatment",
                       "observation", "negative_control", "restriction")

.dictionary_stop <- function(...) stop(..., call. = FALSE)

.q <- function(x) paste0("`", x, "`", collapse = ", ")

# Coerce to a plain data frame of character columns and a logical `follow_up`.
# Factors become character; "TRUE" and "FALSE" become logical. `reason` and
# `sensitivity` may be absent or NA, and then read as empty.
.dictionary_coerce <- function(dictionary) {
  if (!is.data.frame(dictionary))
    .dictionary_stop("`dictionary` must be a data frame with one row per column the ",
                     "lock keeps and the fields ", paste(.dictionary_fields, collapse = ", "), ".")
  d <- as.data.frame(dictionary, stringsAsFactors = FALSE)
  required <- setdiff(.dictionary_fields, c("reason", "sensitivity"))
  lacking <- setdiff(required, names(d))
  if (length(lacking))
    .dictionary_stop("The dictionary has no ", .q(lacking[1]), " field. It needs ",
                     paste(.q(required), collapse = ", "), "; `reason` and `sensitivity` ",
                     "may be left out or empty.")
  as_chr <- function(x) if (is.character(x)) x else as.character(x)
  out <- lapply(stats::setNames(nm = setdiff(.dictionary_fields, "follow_up")),
                function(f) if (f %in% names(d)) as_chr(d[[f]]) else rep("", nrow(d)))
  fu <- d$follow_up
  if (!is.logical(fu)) {
    s <- trimws(as_chr(fu))
    if (any(!is.na(s) & !s %in% c("TRUE", "FALSE")))
      .dictionary_stop("`follow_up` must be TRUE or FALSE for every row.")
    fu <- s == "TRUE"
  }
  out$follow_up <- fu
  out <- out[.dictionary_fields]
  for (f in c("reason", "sensitivity")) out[[f]][is.na(out[[f]])] <- ""
  as.data.frame(out, stringsAsFactors = FALSE)
}

# Validate the dictionary against the lock's inputs and return it normalized: a
# base data frame in the field order above, sorted by `column` (in C order, so the
# order does not depend on the locale), without row names. The lock hashes this
# result, so the type and row order of what the analyst passed do not matter.
# `design_data` is the data the lock keeps columns from.
.check_dictionary <- function(dictionary, design_data, id, treatment, covariates, plan) {
  d <- .dictionary_coerce(dictionary)
  rows_of <- function(x) paste(which(x), collapse = ", ")

  # Required fields, no missing values.
  for (f in setdiff(.dictionary_fields, c("reason", "sensitivity"))) {
    bad <- if (f == "follow_up") is.na(d[[f]]) else is.na(d[[f]]) | !nzchar(trimws(d[[f]]))
    if (any(bad))
      .dictionary_stop("Dictionary field ", .q(f), " is missing or blank in row(s) ",
                       rows_of(bad), ". Every row needs a value for ", .q(f), ".")
  }

  # Allowed values.
  for (f in c("recorded", "fixed")) {
    bad <- !d[[f]] %in% .dictionary_timing
    if (any(bad))
      .dictionary_stop(.q(f), " must be one of ", paste(.dictionary_timing, collapse = ", "),
                       " (relative to time zero); found `", d[[f]][bad][1], "` for column `",
                       d$column[bad][1], "`.")
  }
  bad <- !d$role %in% .dictionary_roles
  if (any(bad))
    .dictionary_stop("`role` must be one of ", paste(.dictionary_roles, collapse = ", "),
                     "; found `", d$role[bad][1], "` for column `", d$column[bad][1], "`.")

  # One row per column.
  dup <- unique(d$column[duplicated(d$column)])
  if (length(dup))
    .dictionary_stop("Duplicate dictionary rows for: ", paste(dup, collapse = ", "),
                     ". Give each column exactly one row.")

  # The rows must be the kept columns: no row missing, none for an absent column.
  nc_cols <- names(plan$negative_controls)
  r_vars <- unique(unlist(lapply(plan$restrictions, all.vars)))
  extra <- setdiff(d$column, names(design_data))
  if (length(extra))
    .dictionary_stop("Dictionary row(s) for column(s) not in `design_data`: ",
                     paste(extra, collapse = ", "), ". Remove those rows or fix the names.")
  kept <- unique(c(id, treatment, covariates, nc_cols, r_vars,
                   d$column[d$role %in% c("after_treatment", "observation")]))
  missing_rows <- setdiff(kept, d$column)
  if (length(missing_rows))
    .dictionary_stop("No dictionary row for column(s): ", paste(missing_rows, collapse = ", "),
                     ". The lock needs one row for each column it keeps (the id, the treatment, ",
                     "the covariates, the negative controls and the restriction variables), ",
                     "and for each column you declare `after_treatment` or `observation`.")

  # The columns of one variable (the dummies of a factor) share their timing and
  # their role, so the identification table does not split a variable.
  for (v in unique(d$variable[duplicated(d$variable)])) {
    g <- d[d$variable == v, , drop = FALSE]
    for (f in c("follow_up", "recorded", "fixed", "role"))
      if (length(unique(g[[f]])) > 1L)
        .dictionary_stop("Columns of variable `", v, "` disagree on ", f, ": ",
                         paste0(g$column, " = ", g[[f]], collapse = ", "),
                         ". Columns of one variable are recorded, fixed and used together; give them ",
                         "the same value, or different `variable` names if they are not one variable.")
  }

  role_of <- stats::setNames(d$role, d$column)
  by_role <- function(r) d[d$role %in% r, , drop = FALSE]

  # Rule 2: roles match how each column is used. A node after treatment in
  # `covariates` gets its own message (rule 5), because it is the common mistake.
  if (role_of[[id]] != "id")
    .dictionary_stop("`", id, "` is the id column but its dictionary role is ", role_of[[id]],
                     ". Give it role id.")
  if (role_of[[treatment]] != "treatment")
    .dictionary_stop("`", treatment, "` is the treatment column but its dictionary role is ",
                     role_of[[treatment]], ". Give it role treatment.")
  for (cv in covariates) {
    r <- role_of[[cv]]
    if (r %in% c("after_treatment", "observation"))
      .dictionary_stop("`", cv, "` is a node after treatment (dictionary role ", r, ") but is in ",
                       "`covariates`. Adjusting for it can bias the total effect; remove it from ",
                       "`covariates`. It stays in the lock for the analysis of the outcome's observation.")
    if (!r %in% c("covariate", "stand_in"))
      .dictionary_stop("`", cv, "` is in `covariates` but its dictionary role is ", r,
                       ". Give it role covariate (or stand_in), or remove it from `covariates`.")
  }
  stray <- by_role(c("covariate", "stand_in"))
  stray <- stray[!stray$column %in% covariates, , drop = FALSE]
  if (nrow(stray))
    .dictionary_stop("`", stray$column[1], "` has dictionary role ", stray$role[1],
                     " but is not in `covariates`. Add it to `covariates` or change its role.")
  for (r in c("id", "treatment")) {
    own <- if (r == "id") id else treatment
    stray <- by_role(r)
    stray <- stray[stray$column != own, , drop = FALSE]
    if (nrow(stray))
      .dictionary_stop("`", stray$column[1], "` has dictionary role ", r, " but is not the ", r,
                       " column (`", own, "`).")
  }
  for (nc in nc_cols)
    if (role_of[[nc]] != "negative_control")
      .dictionary_stop("`", nc, "` is a negative control in the plan but its dictionary role is ",
                       role_of[[nc]], ". Give it role negative_control.")
  stray <- by_role("negative_control")
  stray <- stray[!stray$column %in% nc_cols, , drop = FALSE]
  if (nrow(stray))
    .dictionary_stop("`", stray$column[1], "` has dictionary role negative_control but is not a ",
                     "negative control in the plan. Add it to `negative_controls` or change its role.")
  # A restriction variable that is also a covariate, a control, the id or the
  # treatment takes that role; the restriction use is read from the plan.
  for (rv in setdiff(r_vars, c(id, treatment, covariates, nc_cols)))
    if (role_of[[rv]] != "restriction")
      .dictionary_stop("`", rv, "` is a restriction variable in the plan but its dictionary role is ",
                       role_of[[rv]], ". Give it role restriction.")
  stray <- by_role("restriction")
  stray <- stray[!stray$column %in% r_vars, , drop = FALSE]
  if (nrow(stray))
    .dictionary_stop("`", stray$column[1], "` has dictionary role restriction but is not used in ",
                     "the plan's `restrictions`. Add the restriction or change its role.")

  # Rule 3: a covariate fixed after time zero.
  late <- by_role("covariate")
  late <- late[late$fixed == "after", , drop = FALSE]
  if (nrow(late))
    .dictionary_stop("Covariate `", late$column[1], "` is fixed after time zero (`fixed = \"after\"`), ",
                     "so adjusting for it can block part of the treatment's effect or open a path ",
                     "to the outcome. Leave it out, declare it role `after_treatment`, or declare ",
                     "it `stand_in` with a reason and a sensitivity analysis.")

  # Rule 4: columns from the follow-up source, other than outcome-observation
  # indicators and negative controls.
  fu <- d[d$follow_up & !d$role %in% c("observation", "negative_control"), , drop = FALSE]
  if (nrow(fu))
    .dictionary_stop("`", fu$column[1], "` comes from the follow-up source (`follow_up = TRUE`), ",
                     "recorded only for patients reached at follow-up, so as role ", fu$role[1],
                     " it would depend on who was followed up. Leave it out of the design data, ",
                     "or declare it role `observation` (an outcome-observation indicator) or ",
                     "`negative_control`, with a `reason`.")

  # Rule 6: the treatment is assigned at time zero.
  if (d$fixed[d$column == treatment] != "at")
    .dictionary_stop("The treatment `", treatment, "` must be fixed at time zero (`fixed = \"at\"`), ",
                     "but its dictionary row says `", d$fixed[d$column == treatment], "`. A ",
                     "treatment that is not assigned at time zero means time zero is misplaced; ",
                     "move time zero to the moment of assignment in `analysis_plan()`.")

  # Section 4.2: the declared responses.
  for (i in which(d$role == "stand_in"))
    if (!.has_text(d$reason[i]) || !.has_text(d$sensitivity[i]))
      .dictionary_stop("Stand-in `", d$column[i], "` needs a `reason` and a `sensitivity` in its ",
                       "dictionary row: why it stands in for something fixed before time zero, ",
                       "and the analysis that leaves it out.")
  for (i in which(d$role == "observation"))
    if (!.has_text(d$reason[i]))
      .dictionary_stop("Outcome-observation indicator `", d$column[i], "` needs a `reason` in its ",
                       "dictionary row: the approval for having it in the design data, for ",
                       "example the decision and its date.")
  for (i in which(d$role == "negative_control" & d$follow_up))
    if (!.has_text(d$reason[i]))
      .dictionary_stop("Negative control `", d$column[i], "` comes from the follow-up source, ",
                       "so it is observed only for patients reached and reveals who was followed ",
                       "up. Give a `reason` in its dictionary row for using it.")
  for (i in which(d$role == "restriction" & d$fixed == "after"))
    if (!.has_text(d$reason[i]))
      .dictionary_stop("Restriction variable `", d$column[i], "` is fixed after time zero, so the ",
                       "restriction can select on events after assignment. Give a `reason` in its ",
                       "dictionary row: why the restriction is acceptable.")

  d <- d[order(d$column, method = "radix"), , drop = FALSE]
  rownames(d) <- NULL
  d
}
