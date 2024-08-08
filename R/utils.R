#' @title Translate survival function to cumulative hazard function
#'
#' @param survival_function A numeric vector of survival probabilities.
#'
#' @export
translate_survival_to_cumulative_hazard <-
  function(survival_function) {
    return(-log(survival_function))
  }

#' @title Translate cumulative hazard function to survival function
#'
#' @param cumulative_hazard_function A numeric vector of cumulative hazard values.
#'
#' @export
translate_cumulative_hazard_to_survival <-
  function(cumulative_hazard_function) {
    return(exp(-cumulative_hazard_function))
  }
