#' @keywords internal
#' @importFrom ggplot2 .data
"_PACKAGE"

# `weights` is the column MatchIt::match.data() adds, used unquoted in
# WeightIt::lm_weightit(weights = weights).
utils::globalVariables("weights")
