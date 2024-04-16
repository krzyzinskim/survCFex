#' @title Kovalev method for counterfactual explanation
#' @description This function implements themethod for counterfactual explanations of survival models proposed by Kovalev et al.
#'
#' @param explainer a survival explainer object
#' @param new_observation a new observation for which we want to find a counterfactual
#' @param num_iter number of iterations in PSO algorithm
#' @param num_particles number of particles in PSO algorithm
#' @param verbose if integer, print loss every verbose iterations
#' @param r smallest distance between mean times to events of the original and counterfactual observations
#' @param C penalty parameter for the loss function (used to enforce the counterfactuality)
#' @param w inertia weight in PSO algorithm
#' @param c1 cognitive parameter in PSO algorithm
#' @param c2 social parameter in PSO algorithm
#' @param seed seed for random number generation
#'
#' @import survex
#' @export
kovalev_method <- function(explainer, new_observation, num_iter=200, num_particles=1000,
                           verbose=10, r=50, C=1e6, w=0.729, c1=1.4945, c2=1.4945, seed=NULL){
  # mean value of the original observation
  mean_value <- mean_time_to_survival(explainer, new_observation)

  # scale data
  z_candidates <- explainer$data
  mu <- colMeans(z_candidates)
  sigma <- apply(z_candidates, 2, sd)
  z_candidates <- data.frame(scale(z_candidates, center = mu, scale = sigma))
  x <- scale(new_observation, center = mu, scale = sigma)

  # bounds of the domain
  data_range <- t(apply(z_candidates, 2, range))


  # for current candidates (explainer$data)
  # calculate distances to the original observation
  xz_distances <- euclidean_distance_loss(x, z_candidates)
  # calculate losses
  losses <- kovalev_counterfactual_loss(explainer, mean_value,
                                rescale_to_original(z_candidates, mu, sigma),
                                xz_distances)

  # find the best counterfactual example from current candidates
  z_closest <- z_candidates[which.min(losses), ]
  radius_closest <- sqrt(sum((x - z_closest)^2))

  # iteration 0 - initialization
  set.seed(seed)
  particles <- data.frame(apply(data_range, 1, function(x) runif(num_particles, x[1], x[2])))
  particles <- restriction_procedure(x, particles, euclidean_distance_loss(x, particles), radius_closest, data_range)
  particles[1, ] <- as.matrix(z_closest)

  velocities <- matrix(0, nrow = num_particles, ncol = ncol(particles))
  best_positions <- particles

  best_losses <- kovalev_counterfactual_loss(explainer, mean_value,
                                     rescale_to_original(best_positions, mu, sigma),
                                     euclidean_distance_loss(x, best_positions))
  best_loss <- min(best_losses)
  best_position <- best_positions[which.min(best_losses),]

  for (i in 1:num_iter){
    # update velocities
    velocities <- 0.729 * velocities +
      1.49445 * runif(num_particles, 0, 1) * (best_positions - particles) +
      1.49445 * runif(num_particles, 0, 1) * (as.list(best_position) - particles)

    # update positions
    particles <- particles + velocities
    particles <- restriction_procedure(x, particles,
                                       euclidean_distance_loss(x, particles),
                                       radius_closest, data_range)

    # update best positions
    losses <- kovalev_counterfactual_loss(explainer, mean_value,
                                  rescale_to_original(particles, mu, sigma),
                                  euclidean_distance_loss(x, particles))
    best_positions[losses < best_losses, ] <- particles[losses < best_losses, ]
    best_losses[losses < best_losses] <- losses[losses < best_losses]

    # update best loss
    if (min(best_losses) < best_loss){
      best_loss <- min(best_losses)
      best_position <- best_positions[which.min(best_losses), ]
    }

    if (is.numeric(verbose) & i %% verbose == 0){
      cat(paste0("--> Iteration ", i), "\n")
      cat(paste0("Best loss: ", best_loss), "\n")
      cat("Counterfactual:", "\n")
      print(rescale_to_original(best_position, mu, sigma))
      cat("\n")
    }
  }

  list(
    z = rescale_to_original(best_position, mu, sigma),
    loss = best_loss,
    mean_time_to_survival = mean_time_to_survival(explainer, rescale_to_original(best_position, mu, sigma))
  )
}


mean_time_to_survival <- function(explainer, new_observations){
  survival_function <- predict(explainer, new_observations, output_type = "survival")
  n <- length(explainer$times)
  time_diffs <- diff(c(0, explainer$times))
  # for multiple observations in new_observations (matrix where rows are survival functions)
  as.numeric(0.5 * time_diffs %*%
               t((cbind(rep(1, nrow(new_observations)), survival_function[, 1:(n-1), drop = FALSE]) +
                    survival_function)))
}


restriction_procedure <- function(x, z_candidates, xz_distances, radius_closest, data_range){
  z_candidates <- as.list(x) + pmin(radius_closest, xz_distances) * (z_candidates - as.list(x)) / xz_distances
  # check if z_candidates are in the domain
  z_candidates <- t(apply(z_candidates, 1, pmin, data_range[, 2]))
  z_candidates <- t(apply(z_candidates, 1, pmax, data_range[, 1]))
  data.frame(z_candidates)
}

rescale_to_original <- function(x, mu, sigma){
  t(t(x) * sigma + mu)
}
