#' @export
translate_survival_to_cumulative_hazard <- function(survival_function) {
  return(-log(survival_function))
}

#' @export
translate_cumulative_hazard_to_survival <- function(cumulative_hazard_function) {
  return(exp(-cumulative_hazard_function))
}
