#' @keywords internal
.onLoad <- function(libname, pkgname) {
  # Generate example1 data if not already available
  data_path <- system.file("data", "example1.rda", package = pkgname)
  if (!nzchar(data_path) || !file.exists(data_path)) {
    # Data will be generated on first access via lazy loading
    # or user can call sim_func1() directly
  }
}

.onAttach <- function(libname, pkgname) {
  packageStartupMessage(
    "cleanTMLE: Causal Inference for Risk Estimation\n",
    "Use sim_func1() to generate example data, or data(example1) if available."
  )
}
