#' @export
print.counterfactual_explanations <- function(x, ...) {
  cat("Counterfactual explanations\n")
  cat(paste("\n-> Type:", class(x)[2], "\n"))
  cat("\n-> Original observation:\n")
  print(x$original_observation)
  cat("\n-> Objective values:\n")
  print(x$objective_values)
  cat("\n-> Found counterfactual examples:\n")
  print(x$counterfactual_examples)
  cat("\n-> Analyze the results with the following function:\n")
  cat("`analyze(counterfactual_explanations)`")
}
