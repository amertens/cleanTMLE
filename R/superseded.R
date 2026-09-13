# Superseded-API notes.
#
# cleanTMLE 0.2.0 refocused the package on support assessment, estimand
# feasibility, outcome-blind simulation, and the estimand ladder. The
# checkpoint, gate-token, audit-log and decision-log helpers from 0.1.x keep
# working, but they are no longer part of the recommended workflow, are no
# longer shown in the vignettes or the pkgdown index, and emit a one-line
# note (once per session) pointing to their replacements. Recording and
# reporting live in the lock's design log and in the cleanroomGov companion
# package; enforcement beyond outcome masking is opt-in via the two-pass
# entry point.

#' @keywords internal
.superseded <- function(name, alternative) {
  rlang::inform(
    sprintf(paste0("%s() is superseded in cleanTMLE 0.2.0 and is no longer ",
                   "part of the recommended workflow (it keeps working). ",
                   "See %s."),
            name, alternative),
    .frequency = "once", .frequency_id = paste0("superseded_", name))
  invisible(NULL)
}
