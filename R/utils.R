translate_survival_to_cumulative_hazard <- function(survival_function) {
  return(-log(survival_function))
}

translate_cumulative_hazard_to_survival <- function(cumulative_hazard_function) {
  return(exp(-cumulative_hazard_function))
}

translate_target_envelope <- function(target_envelope){
  stopifnot("target_envelope should be of class 'target_envelope'" = class(target_envelope) == "target_envelope")
  stopifnot("target_envelope should have a type" = "type" %in% names(target_envelope))
  if (target_envelope$type == "survival"){
    new_target_envelope <- target_envelope
    new_target_envelope$upper_bound <- translate_survival_to_cumulative_hazard(target_envelope$lower_bound)
    new_target_envelope$lower_bound <- translate_survival_to_cumulative_hazard(target_envelope$upper_bound)
    new_target_envelope$type <- "cumulative_hazard"
  } else {
    new_target_envelope <- target_envelope
    new_target_envelope$upper_bound <- translate_cumulative_hazard_to_survival(target_envelope$lower_bound)
    new_target_envelope$lower_bound <- translate_cumulative_hazard_to_survival(target_envelope$upper_bound)
    new_target_envelope$type <- "survival"
  }
  return(new_target_envelope)
}


print.counterfactual_explanations <- function(x, ...){
  cat("Counterfactual explanations\n")
  cat(paste("\n-> Type:", class(x)[2], "\n"))
  cat("\n-> Original observation:\n")
  print(x$original_observation)
  cat("\n-> Objective values:\n")
  print(x$objective_values)
  cat("\n-> Found counterfactual examples:\n")
  print(x$population)
  cat("\n-> Analyze the results with the following function:\n")
  cat("`analyze(counterfactual_explanations)`")
}
