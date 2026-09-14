# Superseded-API notes.
#
# 0.2.0 unexported the checkpoint, gate-token, audit-log, and old
# decision-log helpers; 0.3.0 deleted them (see NEWS). The handful of
# internal engines that remain superseded but working emit this
# one-line note (once per session) pointing to their replacements.
# Recording and reporting live in the lock's design log; enforcement
# beyond outcome masking is create_analysis_lock(enforce = TRUE) with
# a named unmasking approver.

#' @keywords internal
.superseded <- function(name, alternative) {
  rlang::inform(
    sprintf(paste0("%s() is superseded and is no longer part of the ",
                   "recommended workflow (it keeps working). See %s."),
            name, alternative),
    .frequency = "once", .frequency_id = paste0("superseded_", name))
  invisible(NULL)
}
