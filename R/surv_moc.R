#' @title SurvMOC: Survival Multi-Objective Counterfactuals
#'
#' @description SurvMOC is a time-dependent model-agnostic counterfactual explanation method for survival analysis models. It generates counterfactual examples using multi-objective optimization.
#'
#' @param explainer A `survex::survival_explainer` object - model preprocessed by the `survex::explain()` function.
#' @param new_observation A data frame with a single row containing the observation for which counterfactuals are to be generated.
#' @param target_prediction_band An object of class `target_prediction_band` containing the target prediction band for the counterfactuals.
#' @param weights A numeric vector with non-negative values and sum(weights) > 0. Weights for each time point. Default is `survex::survival_weights(explainer, p=0, q=0)`, which assigns equal weights to all time points.
#' @param times A numeric vector with time points for which the predictions are to be calculated.
#' @param background_data A data frame with the background data used to train the model. If not provided, the data from the explainer object is used.
#' @param target_change A numeric vector with the target change in the prediction expresses as survival distance. If provided, the counterfactuals will be generated to change the prediction by this amount.
#' @param fixed_variables_indices A vector of positive integers indicating the indices of variables that should remain fixed in the counterfactuals.
#' @param numerical_variables_indices A vector of positive integers indicating the indices of variables that are numerical. Variable types are inferred from the background data if not provided.
#' @param categorical_variables_indices A vector of positive integers indicating the indices of variables that are categorical. Variable types are inferred from the background data if not provided.
#' @param integer_variables_indices A vector of positive integers indicating the indices of variables that are integer. Variable types are inferred from the background data if not provided.
#' @param binary_variables_indices A vector of positive integers indicating the indices of variables that are binary. Variable types are inferred from the background data if not provided.
#' @param plausible_data_range A matrix with 2 rows and length(numerical_variables_indices) columns indicating the plausible range for each numerical variable. If not provided, the range is inferred from the background data.
#' @param plausible_categorical_values A list with length(categorical_variables_indices) elements indicating the plausible values for each categorical variable. If not provided, the values are inferred from the background data.
#' @param initialization_type A character string indicating the initialization type for the counterfactuals. Options are "ice" (default) and "background". If "ice" is selected, the counterfactuals are initialized based on the survival adaptation of Individual Conditional Expctation curves. If "background" is selected, the counterfactuals are initialized using the background data.
#' @param data_distance_p A numeric value indicating the p parameter for the Gower distance calculation. Default is 1.
#' @param plausibility_k_neighbors A numeric value indicating the number of neighbors for the plausibility calculation. Default is 5.
#' @param sbx_eta A numeric value indicating the eta parameter for the simulated binary crossover (SBX) operator. Default is 5.
#' @param gene_recombination_probability A numeric value indicating the probability of recombination for each gene. Default is 0.9.
#' @param mutation_probability A numeric value indicating the probability of mutation for each individual. Default is 0.8.
#' @param gene_mutation_probability A numeric value indicating the probability of mutation for each gene. Default is 0.7.
#' @param revert_change_probability A numeric value indicating the probability of reverting the change in each gene. Default is 0.3.
#' @param max_generations A numeric value indicating the maximum number of generations for the genetic algorithm. Default is 100.
#' @param population_size A numeric value indicating the population size for the genetic algorithm. Default is 20.
#' @param validity_threshold A function or a numeric value indicating the threshold for the validity of the counterfactuals. If validity value is greater than the threshold, the counterfactual candidates are considered invalid and additionally penalized. Default is the median of the validity values.
#' @param tol_iter A numeric value indicating the number of iterations with insufficient progress to stop the search. Default is 3.
#' @param verbose A logical value indicating whether to print the progress of the search. Default is FALSE.
#' @param seed An integer value to set the seed for reproducibility.
#'
#' @import survex
#' @export
surv_moc <- function(explainer,
                     new_observation,
                     target_prediction_band = NULL,
                     weights = survival_weights(explainer, p = 0, q = 0),
                     times = explainer$times,
                     background_data = NULL,
                     target_change = NULL,
                     fixed_variables_indices = NULL,
                     numerical_variables_indices = NULL,
                     categorical_variables_indices = NULL,
                     integer_variables_indices = NULL,
                     binary_variables_indices = NULL,
                     plausible_data_range = NULL,
                     plausible_categorical_values = NULL,
                     initialization_type = "ice",
                     data_distance_p = 1,
                     plausibility_k_neighbors = 5,
                     sbx_eta = 5,
                     gene_recombination_probability = 0.9,
                     mutation_probability = 0.8,
                     gene_mutation_probability = 0.7,
                     revert_change_probability = 0.3,
                     max_generations = 100,
                     population_size = 20,
                     validity_threshold = function(x)
                       quantile(x, 0.5, names = FALSE),
                     tol_iter = 3,
                     verbose = FALSE,
                     seed = NULL) {
  p <- ncol(new_observation)
  stopifnot(
    "Weights must be a numeric vector with non-negative values and sum(weights) > 0" = is.numeric(weights) &
      all(weights >= 0) & sum(weights) > 0
  )
  stopifnot(
    "Either target_change or target_prediction_band must be provided" = !is.null(target_change) |
      !is.null(target_prediction_band)
  )
  check_target_prediction_band(target_prediction_band, times)
  stopifnot("Length of times and weights must be the same" = length(times) == length(weights))
  stopifnot(
    "fixed_variables_indices must be a vector of positive integers" =
      is.null(fixed_variables_indices) |
      all(fixed_variables_indices >= 1 &
            fixed_variables_indices <= p)
  )
  stopifnot(
    "numerical_variables_indices must be a vector of positive integers" =
      is.null(numerical_variables_indices) |
      all(
        numerical_variables_indices >= 1 &
          numerical_variables_indices <= p
      )
  )
  stopifnot(
    "categorical_variables_indices must be a vector of positive integers" =
      is.null(categorical_variables_indices) |
      all(
        categorical_variables_indices >= 1 &
          categorical_variables_indices <= p
      )
  )
  stopifnot(
    "integer_variables_indices must be a vector of positive integers" =
      is.null(integer_variables_indices) |
      all(
        integer_variables_indices >= 1 & integer_variables_indices <= p
      )
  )
  stopifnot(
    "integer_variables_indices must be a subset of numerical_variables_indices" =
      is.null(integer_variables_indices) |
      all(integer_variables_indices %in% numerical_variables_indices)
  )
  stopifnot(
    "binary_variables_indices must be a vector of positive integers" =
      is.null(binary_variables_indices) |
      all(binary_variables_indices >= 1 &
            binary_variables_indices <= p)
  )
  stopifnot(
    "binary_variables_indices must be a subset of categorical_variables_indices" =
      is.null(binary_variables_indices) |
      all(
        binary_variables_indices %in% categorical_variables_indices
      )
  )
  stopifnot(
    "plausible_data_range must be a matrix with 2 rows and length(numerical_variables_indices) columns" =
      is.null(plausible_data_range) |
      is.matrix(plausible_data_range) &
      nrow(plausible_data_range) == 2 &
      ncol(plausible_data_range) == length(numerical_variables_indices)
  )
  stopifnot(
    "plausible_categorical_values must be a list with length(categorical_variables_indices) elements" =
      is.null(plausible_categorical_values) |
      is.list(plausible_categorical_values) &
      length(plausible_categorical_values) == length(categorical_variables_indices)
  )
  stopifnot(
    "Either background_data or explainer with saved data must be provided" =
      !is.null(background_data) | !is.null(explainer$data)
  )
  if (is.null(background_data)) {
    background_data <- explainer$data
  }
  stopifnot("background_data must be a data frame" = is.data.frame(background_data))
  stopifnot("Number of columns of background_data must be equal to p" = ncol(background_data) == p)
  stopifnot("Number of rows of background_data must be greater than 0" = nrow(background_data) > 0)

  if (!is.null(seed)) {
    set.seed(seed)
  }

  provided_indices <- unique(
    c(
      fixed_variables_indices,
      numerical_variables_indices,
      categorical_variables_indices,
      integer_variables_indices,
      binary_variables_indices
    )
  )
  if (length(provided_indices) < p) {
    classes <- sapply(background_data, class)
    integer_valued <-
      sapply(background_data[classes != "factor"], function(x)
        all(x %% 1 == 0))
    additional_numerical_indices <-
      union(which(classes == "numeric"), which(classes == "integer"))
    additional_categorical_indices <- which(classes == "factor")
    additional_integer_indices <-
      union(which(classes == "integer"), which(classes != "factor")[integer_valued])

    if (length(additional_categorical_indices))
      categorical_variables_indices <-
      union(categorical_variables_indices,
            additional_categorical_indices)
    if (length(additional_numerical_indices))
      numerical_variables_indices <-
      setdiff(
        union(
          numerical_variables_indices,
          additional_numerical_indices
        ),
        categorical_variables_indices
      )
    if (length(additional_integer_indices))
      integer_variables_indices <-
      setdiff(
        union(integer_variables_indices, additional_integer_indices),
        categorical_variables_indices
      )
  }
  categorical_variables_indices <-
    sort(categorical_variables_indices)
  numerical_variables_indices <- sort(numerical_variables_indices)
  integer_variables_indices <- sort(integer_variables_indices)
  binary_variables_indices <- sort(binary_variables_indices)
  fixed_variables_indices <- sort(fixed_variables_indices)

  stopifnot(
    "numerical_variable_indices and categorical_variable_indices must cover all columns" =
      length(unique(
        c(
          numerical_variables_indices,
          categorical_variables_indices
        )
      )) == p
  )
  data_range <-
    apply(background_data[numerical_variables_indices], 2, range)

  if (is.null(plausible_data_range)) {
    plausible_data_range <- data_range
  }
  if (is.null(plausible_categorical_values)) {
    plausible_categorical_values <-
      lapply(background_data[categorical_variables_indices], unique)
  }
  lower <- plausible_data_range[1, ]
  upper <- plausible_data_range[2, ]
  stdevs <-
    apply(background_data[numerical_variables_indices], 2, sd)

  original_sf <- predict(explainer, new_observation, times = times)

  population <-
    initialize_population(
      initialization_type,
      explainer,
      times,
      weights,
      target_prediction_band,
      new_observation,
      background_data,
      population_size,
      lower,
      upper,
      stdevs,
      numerical_variables_indices,
      categorical_variables_indices,
      integer_variables_indices,
      binary_variables_indices,
      fixed_variables_indices,
      plausible_categorical_values
    )

  population_sfs <- predict(explainer, population, times = times)
  objective_values <-
    calculate_objective_values(
      background_data,
      times,
      weights,
      new_observation,
      population,
      original_sf,
      population_sfs,
      target_change,
      target_prediction_band,
      data_range,
      categorical_variables_indices,
      data_distance_p,
      plausibility_k_neighbors
    )

  history <- data.frame(generation = 0, objective_values)
  avg_objective_values <- colMeans(objective_values)

  children <-
    make_new_population(
      population,
      objective_values,
      population_size,
      get_crowded_comparison_order,
      new_observation,
      background_data,
      fixed_variables_indices,
      numerical_variables_indices,
      categorical_variables_indices,
      integer_variables_indices,
      binary_variables_indices,
      lower,
      upper,
      plausible_categorical_values,
      stdevs,
      sbx_eta,
      gene_recombination_probability,
      mutation_probability,
      gene_mutation_probability,
      revert_change_probability,
      validity_threshold
    )

  population <- rbind(population, children)

  if (verbose) {
    pb <-
      txtProgressBar(
        min = 0,
        max = max_generations,
        initial = 0,
        style = 3
      )
    cat("\n")
    print(avg_objective_values)
  }

  invalid_counter <- 0

  for (i in seq_len(max_generations) - 1) {
    population_sfs <- predict(explainer, population, times = times)
    objective_values <-
      calculate_objective_values(
        background_data,
        times,
        weights,
        new_observation,
        population,
        original_sf,
        population_sfs,
        target_change,
        target_prediction_band,
        data_range,
        categorical_variables_indices,
        data_distance_p,
        plausibility_k_neighbors
      )
    crowded_comparison_order <-
      get_crowded_comparison_order(objective_values, validity_threshold)
    selected_indices <-
      select_population_indices(crowded_comparison_order, population_size)
    if (mean(selected_indices > population_size) < 0.1) {
      invalid_counter <- invalid_counter + 1
      if (invalid_counter >= tol_iter) {
        warning(
          paste(
            "Seach stopped due to insufficient progress in the last",
            tol_iter,
            "iterations."
          )
        )
        break
      }
    } else {
      invalid_counter <- 0
    }
    population <- population[selected_indices, ]
    objective_values <- objective_values[selected_indices, ]
    avg_objective_values <- colMeans(objective_values)
    if (verbose) {
      setTxtProgressBar(pb, i)
      cat("\n")
      print(avg_objective_values)
    }

    history <-
      rbind(history, data.frame(generation = i, objective_values))

    children <-
      make_new_population(
        population,
        objective_values,
        population_size,
        get_crowded_comparison_order,
        new_observation,
        background_data,
        fixed_variables_indices,
        numerical_variables_indices,
        categorical_variables_indices,
        integer_variables_indices,
        binary_variables_indices,
        lower,
        upper,
        plausible_categorical_values,
        stdevs,
        sbx_eta,
        gene_recombination_probability,
        mutation_probability,
        gene_mutation_probability,
        revert_change_probability,
        validity_threshold
      )
    population <- rbind(population, children)
    row.names(population) <- NULL
  }

  population_sfs <- predict(explainer, population, times = times)
  objective_values <-
    calculate_objective_values(
      background_data,
      times,
      weights,
      new_observation,
      population,
      original_sf,
      population_sfs,
      target_change,
      target_prediction_band,
      data_range,
      categorical_variables_indices,
      data_distance_p,
      plausibility_k_neighbors
    )
  crowded_comparison_order <-
    get_crowded_comparison_order(objective_values, validity_threshold)
  selected_indices <-
    select_population_indices(crowded_comparison_order, population_size)
  population <- population[selected_indices, ]
  objective_values <- objective_values[selected_indices, ]
  order <- get_crowded_comparison_order(objective_values)
  population <- population[order, ]
  objective_values <- objective_values[order,]
  row.names(population) <- NULL
  row.names(objective_values) <- NULL
  history <-
    rbind(history, data.frame(generation = i + 1, objective_values))
  if (verbose) {
    setTxtProgressBar(pb, i + 1)
    cat("\n")
    print(avg_objective_values)
    close(pb)
  }

  final_preds <- predict(explainer, population, times = times)
  result <- list(
    original_observation = new_observation,
    target_prediction_band = target_prediction_band,
    times = times,
    history = history,
    original_prediction = original_sf,
    predictions = final_preds,
    objective_values = objective_values,
    counterfactual_examples = population
  )
  class(result) <-
    c("counterfactual_explanations",
      "multiobjective_counterfactuals")
  return(result)
}

make_new_population <- function(population,
                                objective_values,
                                population_size,
                                sorting_function,
                                new_observation,
                                background_data,
                                fixed_variables_indices,
                                numerical_variables_indices,
                                categorical_variables_indices,
                                integer_variables_indices,
                                binary_variables_indices,
                                lower,
                                upper,
                                plausible_categorical_values,
                                stdevs,
                                sbx_eta,
                                gene_recombination_probability,
                                mutation_probability,
                                gene_mutation_probability,
                                revert_change_probability,
                                validity_threshold = NULL) {
  parents_indices <-
    binary_tournament_selection(sorting_function(objective_values, validity_threshold),
                                population_size)
  children <- recombination(
    population,
    parents_indices,
    gene_recombination_probability,
    numerical_variables_indices,
    categorical_variables_indices,
    integer_variables_indices,
    sbx_eta,
    lower,
    upper
  )

  children <-
    mutation(
      children,
      mutation_probability,
      gene_mutation_probability,
      numerical_variables_indices,
      categorical_variables_indices,
      integer_variables_indices,
      binary_variables_indices,
      lower,
      upper,
      stdevs,
      plausible_categorical_values
    )

  children <- transform_to_original_x(children,
                                      new_observation,
                                      fixed_variables_indices,
                                      revert_change_probability = 0.05)

  children <- data.frame(children)
  colnames(children) <- colnames(background_data)
  return(children)
}


calculate_objective_values <-
  function(background_data,
           times,
           weights,
           original_x,
           candidates_z,
           original_sf,
           candidates_sfs,
           target_change,
           target_prediction_band,
           data_range,
           categorical_variables_indices,
           data_distance_p,
           plausibility_k_neighbors) {
    if (!is.null(target_prediction_band)) {
      prediction_loss_val <-
        distance_from_target_prediction_band(candidates_sfs, target_prediction_band, times, weights)
    } else if (!is.null(target_change)) {
      prediction_loss_val <-
        distance_from_target_change(original_sf, candidates_sfs, target_change, times, weights)
    } else {
      stop("Either target_change or target_prediction_band must be provided")
    }

    data_distance_loss_val <-
      gower_distance_loss(
        original_x,
        candidates_z,
        data_range,
        categorical_variables_indices,
        data_distance_p
      )
    sparsity_loss_val <- sparsity_loss(original_x, candidates_z)
    plausibility_loss_val <-
      plausiblity_knn_loss(
        candidates_z,
        background_data,
        data_range,
        categorical_variables_indices,
        plausibility_k_neighbors
      )

    return(data.frame(
      list(
        "validity" = prediction_loss_val,
        "similarity" = data_distance_loss_val,
        "sparsity" = sparsity_loss_val,
        "plausibility" = plausibility_loss_val
      )
    ))
  }


get_crowded_comparison_order <-
  function(obj, validity_threshold = NULL) {
    n <- nrow(obj)

    nondom_sorting <-
      ecr::doNondominatedSorting(as.matrix(t(obj)))$ranks

    if (is.null(validity_threshold)) {
      validity_threshold <- max(obj$validity) + 1
    } else if (is.function(validity_threshold)) {
      validity_threshold <- validity_threshold(obj$validity)
    } else if (validity_threshold == "median") {
      validity_threshold <- median(obj$validity)
    }
    # move invalid solutions to the end (starting from the least invalid)
    violating_samples_indices <-
      which(obj$validity > validity_threshold)
    violating_order <- order(obj$validity[violating_samples_indices])
    nondom_sorting[violating_samples_indices] <-
      max(nondom_sorting) + violating_order
    nondom_sorting <- as.numeric(factor(nondom_sorting))


    crowddist_sorting <- numeric(n)
    for (i in 1:max(nondom_sorting)) {
      crowddist_sorting[nondom_sorting == i] <-
        ecr::computeCrowdingDistance(as.matrix(t(obj[nondom_sorting == i, ]))) +
        obj$similarity[nondom_sorting == i]
    }

    final_sorting <- numeric(n)
    final_sorting[order(nondom_sorting, -crowddist_sorting)] <- 1:n
    return(final_sorting)
  }


select_population_indices <-
  function(crowded_comparison_order,
           population_size) {
    selected_indices <-
      which(crowded_comparison_order <= population_size)
    return(selected_indices)
  }


binary_tournament_selection <-
  function(crowded_comparison_order, num_parents) {
    population_size <- length(crowded_comparison_order)
    selected <- numeric(num_parents)
    sample_a <- sample(1:population_size, num_parents, replace = TRUE)
    sample_b <- sample(1:population_size, num_parents, replace = TRUE)
    a_sel <-
      crowded_comparison_order[sample_a] < crowded_comparison_order[sample_b]
    selected[a_sel] <- sample_a[a_sel]
    selected[!a_sel] <- sample_b[!a_sel]
    return(selected)
  }


recombination <- function(population,
                          parents_indices,
                          gene_recombination_probability,
                          numerical_variables_indices,
                          categorical_variables_indices,
                          integer_variables_indices,
                          sbx_eta = 5,
                          lower,
                          upper) {
  n <- length(parents_indices)
  children <- population
  stopifnot(
    length(numerical_variables_indices) == length(lower) &
      length(numerical_variables_indices) == length(upper)
  )

  which_numerical_integer <-
    numerical_variables_indices %in% integer_variables_indices

  lower[which_numerical_integer] <-
    lower[which_numerical_integer] - 0.5
  upper[which_numerical_integer] <-
    upper[which_numerical_integer] + 0.5

  for (i in 1:(n / 2)) {
    parents_pair <-
      list(as.numeric(population[parents_indices[2 * i - 1], numerical_variables_indices]),
           as.numeric(population[parents_indices[2 * i], numerical_variables_indices]))

    # for numerical features apply simulated binary crossover (SBX)
    new_children <- ecr::recSBX(
      parents_pair,
      sbx_eta,
      gene_recombination_probability,
      as.numeric(lower),
      as.numeric(upper)
    )
    children[2 * i - 1, numerical_variables_indices] <-
      new_children[[1]]
    children[2 * i, numerical_variables_indices] <-
      new_children[[2]]

    # for categorical features apply uniform crossover
    parents_pair <- list(population[parents_indices[2 * i - 1],],
                         population[parents_indices[2 * i],])
    for (j in categorical_variables_indices) {
      rnd_parent <- sample(1:2, 2, replace = TRUE)
      children[2 * i - 1, j] <- parents_pair[[rnd_parent[1]]][, j]
      children[2 * i, j] <- parents_pair[[rnd_parent[2]]][, j]
    }
  }

  children[, integer_variables_indices] <-
    round(children[, integer_variables_indices])
  return(children)
}


mutation <-
  function(children,
           mutation_probability,
           gene_mutation_probability,
           numerical_variables_indices,
           categorical_variables_indices,
           integer_variables_indices,
           binary_variables_indices,
           lower,
           upper,
           sdevs,
           categorical_levels) {
    stopifnot(
      length(numerical_variables_indices) == length(lower) &
        length(numerical_variables_indices) == length(upper) &
        length(numerical_variables_indices) == length(sdevs)
    )

    which_numerical_integer <-
      numerical_variables_indices %in% integer_variables_indices
    lower[which_numerical_integer] <-
      lower[which_numerical_integer] - 0.5
    upper[which_numerical_integer] <-
      upper[which_numerical_integer] + 0.5

    n_numerical <- length(numerical_variables_indices)
    n_categorical <- length(categorical_variables_indices)

    n <- nrow(children)

    # for numerical features apply scaled Gaussian mutation
    add_matrix <-
      matrix(rnorm(n * n_numerical), nrow = n, ncol = n_numerical)
    add_matrix <-
      as.matrix(as.data.frame(add_matrix) * as.list(sdevs))

    mutation_mask <-
      matrix(
        runif(n * n_numerical) < gene_mutation_probability,
        nrow = n,
        ncol = n_numerical
      )
    rows_to_mutate <- runif(n) < mutation_probability
    children[rows_to_mutate, numerical_variables_indices] <-
      children[rows_to_mutate, numerical_variables_indices] +
      add_matrix[rows_to_mutate, ] * mutation_mask[rows_to_mutate, ]

    # clip to range
    for (i in seq_len(n_numerical)) {
      children[, numerical_variables_indices[i]] <-
        pmax(lower[i], pmin(upper[i], children[, numerical_variables_indices[i]]))
    }

    # round integer variables
    children[, integer_variables_indices] <-
      round(children[, integer_variables_indices])

    # for categorical features apply uniform mutation (randomly select a level different from the current one)
    for (i in seq_len(n_categorical)) {
      mutation_mask <- (runif(n) < gene_mutation_probability)
      if (categorical_variables_indices[i] %in% binary_variables_indices) {
        children[rows_to_mutate &
                   mutation_mask, categorical_variables_indices[i]] <-
          1 - children[rows_to_mutate &
                         mutation_mask, categorical_variables_indices[i]]
      } else{
        children[rows_to_mutate &
                   mutation_mask, categorical_variables_indices[i]] <-
          sample(categorical_levels[[i]], n, replace = TRUE)[rows_to_mutate &
                                                               mutation_mask]
      }
    }
    return(children)
  }

transform_to_original_x <- function(population,
                                    original_x,
                                    fixed_variables_indices,
                                    revert_change_probability) {
  for (i in seq_len(ncol(population))) {
    if (i %in% fixed_variables_indices) {
      population[, i] <- original_x[[i]]
    } else{
      revert_change_mask <-
        runif(nrow(population)) < revert_change_probability
      population[revert_change_mask, i] <- original_x[[i]]
    }
  }
  return(population)
}


initialize_population <-
  function(type,
           explainer,
           times,
           weights,
           target_prediction_band,
           new_observation,
           background_data,
           population_size,
           lower,
           upper,
           stdevs,
           numerical_variables_indices,
           categorical_variables_indices,
           integer_variables_indices,
           binary_variables_indices,
           fixed_variables_indices,
           plausible_categorical_values) {
    if (type == "ice") {
      return(
        initialize_population_ice(
          explainer,
          times,
          weights,
          target_prediction_band,
          new_observation,
          background_data,
          population_size,
          lower,
          upper,
          stdevs,
          numerical_variables_indices,
          categorical_variables_indices,
          integer_variables_indices,
          binary_variables_indices,
          fixed_variables_indices,
          plausible_categorical_values
        )
      )
    } else if (type == "background") {
      return(
        initialize_population_from_background(
          explainer,
          times,
          weights,
          target_prediction_band,
          new_observation,
          background_data,
          population_size,
          lower,
          upper,
          stdevs,
          numerical_variables_indices,
          categorical_variables_indices,
          integer_variables_indices,
          binary_variables_indices,
          fixed_variables_indices,
          plausible_categorical_values
        )
      )
    } else {
      stop("Unknown initialization type")
    }
  }


initialize_population_from_background <-
  function(explainer,
           times,
           weights,
           target_prediction_band,
           new_observation,
           background_data,
           population_size,
           lower,
           upper,
           stdevs,
           numerical_variables_indices,
           categorical_variables_indices,
           integer_variables_indices,
           binary_variables_indices,
           fixed_variables_indices,
           plausible_categorical_values) {
    n_numerical <- length(numerical_variables_indices)
    n_categorical <- length(categorical_variables_indices)

    if (population_size > nrow(background_data)) {
      population <- background_data
    } else{
      population <-
        background_data[sample(1:nrow(background_data), population_size, replace = FALSE), ]
    }

    # clip to plausible range and levels
    for (i in seq_len(n_numerical)) {
      population[, numerical_variables_indices[i]] <-
        pmax(lower[i], pmin(upper[i], population[, numerical_variables_indices[i]]))
    }
    for (i in seq_len(n_categorical)) {
      to_change <-
        !(population[, categorical_variables_indices[i]] %in% plausible_categorical_values[[i]])
      population[to_change, categorical_variables_indices[i]] <-
        sample(plausible_categorical_values[[i]], population_size, replace = TRUE)[to_change]
    }
    population[, integer_variables_indices] <-
      round(population[, integer_variables_indices])

    population <- transform_to_original_x(population, new_observation,
                                          fixed_variables_indices, 0.6)
    row.names(population) <- NULL
    return(population)
  }


initialize_population_ice <-
  function(explainer,
           times,
           weights,
           target_prediction_band,
           new_observation,
           background_data,
           population_size,
           lower,
           upper,
           stdevs,
           numerical_variables_indices,
           categorical_variables_indices,
           integer_variables_indices,
           binary_variables_indices,
           fixed_variables_indices,
           plausible_categorical_values,
           p_max = 0.95,
           p_min = 0.05) {
    ice <-
      predict_profile(explainer, new_observation, variable_splits_type = "quantiles")
    sf <- predict(explainer, new_observation, times = times)
    sf <- as.data.frame(sf)
    probs <- numeric(length(new_observation))
    names(probs) <- colnames(new_observation)

    for (var in colnames(new_observation)) {
      tmp <-
        ice$result[ice$result$`_vname_` == var, c(var, "_times_", "_yhat_")]
      res <- aggregate(`_yhat_` ~ `_times_`, tmp, sd)
      sf_diff <-
        pmax(
          sf - as.list(target_prediction_band$upper_bound),
          as.list(target_prediction_band$lower_bound) - sf
        )
      probs[var] <-
        survival_distance(0, sf_diff * res$`_yhat_`, times, weights)
    }

    min_prob <- min(probs)
    max_prob <- max(probs)
    probs <-
      (probs - min_prob) * (p_max - p_min) / (max_prob - min_prob) + p_min


    n_numerical <- length(numerical_variables_indices)
    n_categorical <- length(categorical_variables_indices)

    population <- new_observation[rep(1, population_size), ]

    for (i in seq_len(n_numerical)) {
      mutation_mask <-
        runif(population_size) < probs[numerical_variables_indices[i]]
      population[mutation_mask, numerical_variables_indices[i]] <-
        runif(population_size, lower[i], upper[i])[mutation_mask]
    }
    for (i in seq_len(n_categorical)) {
      mutation_mask <-
        runif(population_size) < probs[categorical_variables_indices[i]]
      tmp <-
        sample(plausible_categorical_values[[i]], population_size, replace = TRUE)
      population[mutation_mask, categorical_variables_indices[i]] <-
        tmp[mutation_mask]
    }
    population[, integer_variables_indices] <-
      round(population[, integer_variables_indices])
    population <- transform_to_original_x(population, new_observation,
                                          fixed_variables_indices, 0.3)
    row.names(population) <- NULL
    return(population)
  }
