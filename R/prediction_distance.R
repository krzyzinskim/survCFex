#' @export
distance_from_target_prediction_band <- function(sf, prediction_band, times, weights = rep(1/length(times), length(times))){
  if (is.vector(sf)){
    sf <- matrix(sf, nrow = 1)
  }
  sf <- as.data.frame(sf)
  sf_difference <- pmax(sf - as.list(prediction_band$upper_bound),
                        as.list(prediction_band$lower_bound) - sf,
                        data.frame(matrix(0, nrow = nrow(sf), ncol = ncol(sf))))
  return(abs(survival_distance(0, sf_difference, times, weights)))
}

distance_from_target_change <- function(sf, sf_others, target_change, times, weights = rep(1/length(times), length(times))){
  distances_from_original <- survival_distance(sf, sf_others, times, weights)
  distances <- abs(distances_from_original - target_change)
  if (target_change > 0)
    distances[distances_from_original > target_change] <- 0
  else
    distances[distances_from_original < target_change] <- 0
  return(distances)
}

#' @export
survival_distance <- function(sf, sf_others, times, weights = rep(1/length(times), length(times))){
  if (is.vector(sf_others)){
    sf_others <- matrix(sf_others, nrow = 1)
  }
  sf_difference <- (as.list(sf) - as.data.frame(sf_others))^2
  n <- length(times)
  weights <- as.list(weights)
  tmp <- (sf_difference[,1:(n - 1)] * weights[1:(n-1)] + sf_difference[,2:n] * weights[2:n]) * as.list(diff(times)) / 2
  weights <- as.numeric(weights)
  weights_integrated <- sum((weights[1:(n-1)] + weights[2:n]) / 2 * diff(times))
  return(as.numeric(sqrt(rowSums(tmp)) / weights_integrated))
}

mean_time_to_survival <- function(explainer, new_observations=NULL, predictions=NULL){
  if (is.null(predictions)){
    survival_function <- predict(explainer, new_observations, output_type = "survival")
  } else {
    survival_function <- predictions
  }
  n <- length(explainer$times)
  time_diffs <- diff(c(0, explainer$times))
  # for multiple observations in new_observations (matrix where rows are survival functions)
  as.numeric(0.5 * time_diffs %*%
               t((cbind(rep(1, nrow(survival_function)), survival_function[, 1:(n-1), drop = FALSE]) +
                    survival_function)))
}
