#' @title Kovalev method for counterfactual explanation
#' @description This function implements themethod for counterfactual explanations of survival models proposed by Kovalev et al.
#'
#' @param explainer a survival explainer object
#' @param new_observation a new observation for which we want to find a counterfactual
#' @param num_iter number of iterations in PSO algorithm
#' @param num_particles number of particles in PSO algorithm
#' @param verbose if integer, print loss every verbose iterations
#' @param r target difference between the mean time to survival of the counterfactual and the original observation
#' @param C penalty parameter for the loss function (used to enforce the counterfactuality)
#' @param w inertia weight in PSO algorithm
#' @param c1 cognitive parameter in PSO algorithm
#' @param c2 social parameter in PSO algorithm
#' @param seed seed for random number generation
#'
#' @import survex
#' @export
kovalev_method <- function(explainer, new_observation, r, num_iter=200, num_particles=1000,
                           verbose=10, C=1e6, w=0.729, c1=1.4945, c2=1.4945, tol_iter=10, seed=NULL){
  # theta
  theta <- ifelse(r >= 0, 1, -1)
  r <- abs(r)

  # mean value of the original observation
  mean_value <- mean_time_to_survival(explainer, new_observation)

  # scale data
  z_candidates <- explainer$data
  mu <- colMeans(z_candidates)
  sigma <- apply(z_candidates, 2, sd)
  z_candidates <- data.frame(scale(z_candidates, center = mu, scale = sigma))
  x <- scale(new_observation, center = mu, scale = sigma)
  colnames(z_candidates) <- colnames(explainer$data)

  # bounds of the domain
  data_range <- t(apply(z_candidates, 2, range))

  # for current candidates (explainer$data)
  # calculate distances to the original observation
  xz_distances <- euclidean_distance_loss(x, z_candidates)
  # calculate losses
  losses <- kovalev_counterfactual_loss(explainer, mean_value,
                                rescale_to_original(z_candidates, mu, sigma),
                                xz_distances, r, C, theta)

  # find the best counterfactual example from current valid candidates
  validity_mask <- r - theta * (
    mean_time_to_survival(explainer, rescale_to_original(z_candidates, mu, sigma))
    - mean_value) <= 0

  if (sum(validity_mask) == 0){
    radius_closest <- max(xz_distances)
    z_closest <- z_candidates[which.max(xz_distances), ]
  } else{
    z_closest <- z_candidates[validity_mask, ][which.min(losses[validity_mask]), ]
    radius_closest <- sqrt(sum((x - z_closest)^2))
  }

  # iteration 0 - initialization
  invalid_counter <- 0
  set.seed(seed)
  particles <- data.frame(apply(data_range, 1, function(var_range) runif(num_particles, var_range[1], var_range[2])))
  particles <- restriction_procedure(x, particles, euclidean_distance_loss(x, particles), radius_closest, data_range)
  particles[1, ] <- as.matrix(z_closest)

  velocities <- matrix(0, nrow = num_particles, ncol = ncol(particles))
  best_positions <- particles

  best_losses <- kovalev_counterfactual_loss(explainer, mean_value,
                                     rescale_to_original(best_positions, mu, sigma),
                                     euclidean_distance_loss(x, best_positions),
                                     r, C, theta)
  best_loss <- min(best_losses)
  best_position <- best_positions[which.min(best_losses),]

  history <- data.frame(generation = 0, loss = best_losses)

  for (i in 1:num_iter){
    # update velocities
    velocities <- w * velocities +
      c1 * runif(num_particles, 0, 1) * (best_positions - particles) +
      c2 * runif(num_particles, 0, 1) * (as.list(best_position) - particles)
    # update positions
    particles <- particles + velocities
    particles <- restriction_procedure(x, particles,
                                       euclidean_distance_loss(x, particles),
                                       radius_closest, data_range)

    particles <- as.data.frame(particles)
    colnames(particles) <- colnames(explainer$data)

    # update best positions
    losses <- kovalev_counterfactual_loss(explainer, mean_value,
                                  rescale_to_original(particles, mu, sigma),
                                  euclidean_distance_loss(x, particles), r, C, theta)
    best_positions[losses < best_losses, ] <- particles[losses < best_losses, ]
    best_losses[losses < best_losses] <- losses[losses < best_losses]

    history <- rbind(history, data.frame(generation = i, loss = best_losses))

    # update best loss
    if (min(best_losses) < best_loss){
      best_loss <- min(best_losses)
      best_position <- best_positions[which.min(best_losses), ]
      invalid_counter <- 0
    } else {
      invalid_counter <- invalid_counter + 1
      if (invalid_counter >= tol_iter){
        warning(paste("Seach stopped due to no progress in the last",
                      tol_iter,
                      "iterations."))
        break
      }
    }

    if (is.numeric(verbose) & i %% verbose == 0){
      cat(paste0("--> Iteration ", i), "\n")
      cat(paste0("Best loss: ", best_loss), "\n")
      cat("Counterfactual:", "\n")
      print(rescale_to_original(best_position, mu, sigma))
      cat("\n")
    }
  }

  final_preds <- predict(explainer, rescale_to_original(best_positions, mu, sigma))
  final_x_distances <- euclidean_distance_loss(x, best_positions)

  objective_values <- data.frame(
    loss = best_losses,
    similarity = final_x_distances,
    validity = (best_losses - final_x_distances) / C
  )

  result <- list(
    original_observation = new_observation,
    r = r * theta,
    times = explainer$times,
    history = history,
    original_prediction = predict(explainer, new_observation),
    predictions = final_preds,
    original_mean_time_to_survival = mean_value,
    mean_times_to_survival = mean_time_to_survival(explainer, predictions = final_preds),
    objective_values = objective_values,
    counterfactual_examples = rescale_to_original(best_positions, mu, sigma),
    best_counterfactual_example = rescale_to_original(best_position, mu, sigma)
  )
  class(result) <- c("counterfactual_explanations", "Kovalev_method_counterfactuals")
  return(result)
}


restriction_procedure <- function(x, z_candidates, xz_distances, radius_closest, data_range){
  z_candidates <- as.list(x) + pmin(radius_closest, xz_distances) * (z_candidates - as.list(x)) / xz_distances
  # check if z_candidates are in the domain
  z_candidates <- t(apply(z_candidates, 1, pmin, data_range[, 2]))
  z_candidates <- t(apply(z_candidates, 1, pmax, data_range[, 1]))
  data.frame(z_candidates)
}

rescale_to_original <- function(x, mu, sigma){
  tmp <- data.frame(t(t(x) * sigma + mu))
  colnames(tmp) <- colnames(x)
  tmp
}
