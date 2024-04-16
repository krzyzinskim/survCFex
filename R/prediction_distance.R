#' @export
distance_from_target_envelope <- function(sf, envelope, times, weights = rep(1/length(times), length(times))){
  if (is.vector(sf)){
    sf <- matrix(sf, nrow = 1)
  }
  sf <- as.data.frame(sf)
  sf_difference <- pmax(sf - as.list(envelope$upper_bound),
                        as.list(envelope$lower_bound) - sf,
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
  return(as.numeric(sqrt(rowSums(tmp))))
}
