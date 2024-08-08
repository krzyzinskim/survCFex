#' @title ST-PT: Survival Trees - Path Tweaking
#'
#' @description ST-PT is a time-dependent model-specific counterfactual explanation method tailored for tree-based ensembles, particularly random survival forests. It generates counterfactuals based on extraction of counterfactual decision paths from trees
#'
#' @param explainer A `survex::survival_explainer` object - model preprocessed by the `survex::explain()` function.
#' @param new_observation A data frame with a single row containing the observation for which counterfactuals are to be generated.
#' @param target_prediction_band An object of class `target_prediction_band` containing the target prediction band for the counterfactuals.
#' @param weights A numeric vector with non-negative values and sum(weights) > 0. Weights for each time point. Default is `survex::survival_weights(explainer, p=0, q=0)`, which assigns equal weights to all time points.
#' @param times A numeric vector with time points for which the predictions are to be calculated.
#' @param background_data A data frame with the background data used to train the model. If not provided, the data from the explainer object is used.
#' @param paths_per_tree Number of paths to be extracted from each tree. Default is 20.
#' @param paths_per_counterfactual Number of paths to be combined to create a single counterfactual. Default is 1, which means that a process is deterministic. If set to a value greater than 1, the paths are randomly combined to create counterfactuals.
#' @param max_tries Maximum number of tries to generate a counterfactual, used only when `paths_per_counterfactual` > 1. Default is 200.
#' @param step Step size for tweaking the values of features. Default is 1.
#' @param max_counterfactuals Maximum number of counterfactuals to be returned. Default is 40 or `num_trees * paths_per_tree`, whichever is smaller.
#' @param fixed_variables_indices A vector of positive integers indicating the indices of variables that should remain fixed during the counterfactual generation process. Default is `NULL`.
#' @param plausible_values A list with names corresponding to column names of `new_observation`. Each element of the list is a vector of plausible values for the corresponding feature. Default is `NULL`, which means that all unique values from the background data are considered plausible.
#' @param verbose Logical indicating whether to print progress messages. Default is `FALSE`.
#' @param seed An integer to set the random seed for reproducibility. Default is `NULL`.
#'
#' @import survex
#' @export
st_pt <- function(explainer,
                  new_observation,
                  target_prediction_band = NULL,
                  weights = survival_weights(explainer, p = 0, q = 0),
                  times = explainer$times,
                  background_data = NULL,
                  paths_per_tree = 20L,
                  paths_per_counterfactual = 1L,
                  max_tries = 200,
                  step = 1,
                  max_counterfactuals = NULL,
                  fixed_variables_indices = NULL,
                  plausible_values = NULL,
                  verbose = FALSE,
                  seed = NULL) {
  ### CHECKS
  p <- ncol(new_observation)
  times <- explainer$times

  stopifnot("Explained model must be a ranger model" = class(explainer$model) == "ranger")
  stopifnot(
    "Weights must be a numeric vector with non-negative values and sum(weights) > 0" = is.numeric(weights) &
      all(weights >= 0) & sum(weights) > 0
  )
  stopifnot("target_prediction_band must be provided" = !is.null(target_prediction_band))
  stopifnot("Length of times and weights must be the same" = length(times) == length(weights))
  stopifnot(
    "fixed_variables_indices must be a vector of positive integers" =
      is.null(fixed_variables_indices) |
      all(fixed_variables_indices >= 1 &
            fixed_variables_indices <= p)
  )
  if (!is.null(fixed_variables_indices)) {
    fixed_variables_names <-
      colnames(new_observation)[fixed_variables_indices]
  } else {
    fixed_variables_names <- c()
  }
  stopifnot(
    "plausible_values must be a list with names corresponding to column names of new_observation" =
      is.null(plausible_values) |
      is.list(plausible_values) &
      all(names(plausible_values) %in% colnames(new_observation))
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
  stopifnot("paths_per_tree must be positive" = paths_per_tree > 0)
  num_trees <- explainer$model$num.trees
  if (is.null(max_counterfactuals)) {
    max_counterfactuals <- min(40, num_trees * paths_per_tree)
  }
  stopifnot("max_counterfactuals must be positive" = max_counterfactuals > 0)

  possible_values <-
    lapply(background_data, function(x)
      sort(unique(x)))
  ordered_variables_names <- c()

  # get order of categorical variables/factors from model (see respect.unordered.factors argument)
  if (!is.null(explainer$model$forest$covariate.levels)) {
    covariate_levels <- explainer$model$forest$covariate.levels
    for (i in 1:length(possible_values)) {
      if (!is.null(covariate_levels[[i]])) {
        possible_values[[i]] <- covariate_levels[[i]]
        ordered_variables_names <-
          c(ordered_variables_names, names(possible_values)[i])
      }
    }
  }

  categorical_variables_indices <-
    which(colnames(new_observation) %in% ordered_variables_names)
  if (length(categorical_variables_indices) == 0)
    categorical_variables_indices <- NULL
  numerical_variables_indices <-
    which(!colnames(new_observation) %in% ordered_variables_names)
  data_range <-
    apply(background_data[numerical_variables_indices], 2, range)

  # create logical mask of which values are plausible for each feature
  if (is.null(plausible_values)) {
    plausible_mask <-
      lapply(possible_values, function(x)
        rep(TRUE, length(x)))
  } else{
    plausible_mask <- list()
    for (i in 1:length(possible_values)) {
      feature_name <- names(possible_values)[i]
      if (feature_name %in% names(plausible_values)) {
        plausible_mask[[feature_name]] <-
          possible_values[[feature_name]] %in% plausible_values[[feature_name]]
      } else {
        plausible_mask[[feature_name]] <-
          rep(TRUE, length(possible_values[[i]]))
      }
    }
  }

  trees <-
    treeshap::ranger_surv.unify(explainer$model,
                                background_data,
                                type = "chf",
                                times = times)
  predictions <- sapply(trees, function(x)
    x$model$Prediction)
  forest_structure <- trees[[1]]$model
  forest_structure$Prediction <- NULL
  n_trees <- length(unique(forest_structure$Tree))
  predictions <-
    predictions * n_trees # because treeshap assumes that the prediction is the sum of predictions from all trees

  # check if the forest is deep enough to find counterfactuals from single paths
  max_depth <- ceiling(log2(max(forest_structure$Node)))
  p <- ncol(new_observation)
  if (max_depth < p / paths_per_counterfactual) {
    warning(
      "The trees in the forest may be too shallow to find counterfactuals from single paths.
Consider increasing the paths_per_counterfactual parameter."
    )
  }

  all_paths <- data.frame()
  if (verbose) {
    cat("Calculating paths to leaves...\n")
    pb <-
      txtProgressBar(
        min = 0,
        max = n_trees - 1,
        initial = 0,
        style = 3
      )
  }
  for (i in seq_len(n_trees) - 1) {
    if (verbose)
      setTxtProgressBar(pb, i)
    tree <- prepare_tree(forest_structure, i)
    prediction_from_tree <- predictions[forest_structure$Tree == i, ]
    dists <-
      distance_from_target_prediction_band(prediction_from_tree, target_prediction_band, times)
    closest_prediction_ids <-
      order(dists)[1:min(paths_per_tree, sum(!is.na(dists)))]
    path <-
      find_paths_to_leaves(tree, closest_prediction_ids, dists[closest_prediction_ids])
    all_paths <- rbind(all_paths, path)
  }
  if (verbose)
    close(pb)

  if (verbose)
    cat("\nCalculating counterfactuals...\n")

  if (paths_per_counterfactual == 1) {
    new_instances <-
      build_new_instances_deterministicly(
        new_observation,
        n_trees,
        all_paths,
        possible_values,
        plausible_mask,
        fixed_variables_names,
        ordered_variables_names,
        step,
        verbose
      )
  } else {
    if (!is.null(seed)) {
      set.seed(seed)
    }
    all_paths_ids <- unique(all_paths[c("Tree", "Path")])
    new_instances <-
      build_new_instances_randomly(
        new_observation,
        all_paths,
        all_paths_ids,
        possible_values,
        plausible_mask,
        fixed_variables_names,
        ordered_variables_names,
        step,
        paths_per_counterfactual,
        max_tries,
        verbose
      )
  }

  new_instances <- new_instances[!duplicated(new_instances), ]

  # set the same classes as in background_data
  for (var in ordered_variables_names) {
    if (is.factor(background_data[[var]])) {
      new_instances[[var]] <-
        factor(new_instances[[var]], levels = levels(background_data[[var]]))
    }
  }

  counterfactual_predictions <-
    predict(explainer, new_instances, output_type = "chf")
  dist_res <-
    distance_from_target_prediction_band(counterfactual_predictions, target_prediction_band, times)
  n_valid <- sum(dist_res == 0)


  if (n_valid > max_counterfactuals) {
    warning(
      "Number of valid counterfactuals is greater than max_counterfactuals. Returning all valid counterfactuals."
    )
    which_to_return_mask <- dist_res == 0
  } else{
    if (n_valid == 0) {
      warning(
        "No valid counterfactuals found. Try to increase the number of checked paths or change the target prediction band."
      )
    }
    which_to_return_mask <-
      order(dist_res) <= min(max_counterfactuals, nrow(new_instances))
  }

  counterfactual_examples <- new_instances[which_to_return_mask, ]
  dist_res <- dist_res[which_to_return_mask]
  counterfactual_predictions <-
    counterfactual_predictions[which_to_return_mask, ]
  rownames(counterfactual_examples) <- NULL

  objective_values <-  data.frame(
    list(
      "validity" = dist_res,
      "similarity" = gower_distance_loss(
        new_observation,
        counterfactual_examples,
        data_range,
        categorical_variables_indices,
        1
      ),
      "sparsity" = sparsity_loss(new_observation, counterfactual_examples),
      "plausibility" = plausiblity_knn_loss(
        counterfactual_examples,
        background_data,
        data_range,
        categorical_variables_indices,
        5
      )
    )
  )

  original_prediction <-
    predict(explainer,
            new_observation,
            times = times,
            output_type = "chf")

  result = list(
    original_observation = new_observation,
    target_prediction_band = target_prediction_band,
    times = times,
    original_prediction = original_prediction,
    predictions = counterfactual_predictions,
    objective_values = objective_values,
    counterfactual_examples = counterfactual_examples
  )

  class(result) <-
    c("counterfactual_explanations",
      "treebased_counterfactuals")
  return(result)
}


prepare_tree <- function(forest_structure, tree_id) {
  tree <- forest_structure[forest_structure$Tree == tree_id, ]
  tree$Node <- tree$Node + 1
  shift <- min(as.numeric(rownames(tree))) - 1
  tree$Yes <- tree$Yes - shift
  tree$No <- tree$No - shift
  rownames(tree) <- NULL
  tree[na.omit(tree$Yes), "Parent"] <-
    tree[!is.na(tree$Yes), "Node"]
  tree[na.omit(tree$No), "Parent"] <- tree[!is.na(tree$No), "Node"]
  tree
}


find_paths_to_leaves <- function(tree, leaves_ids, validities) {
  paths <- data.frame()
  for (i in seq_len(length(leaves_ids))) {
    leaf_id <- leaves_ids[i]
    path_nodes <- c(leaf_id)
    decisions <- c(NA)
    while (leaf_id != 1) {
      parent_node <- tree[tree[leaf_id, "Parent"], ]
      if (parent_node$Yes == leaf_id) {
        decisions <- c("<=", decisions)
      } else {
        decisions <- c(">", decisions)
      }
      leaf_id <- parent_node$Node
      path_nodes <- c(leaf_id, path_nodes)
    }
    path <- tree[path_nodes, ]
    path$Decision <- decisions
    path$Validity <- validities[i]
    path$CoverChange <- c(abs(diff(path$Cover)), NA) / path$Cover
    path$Path <- i
    rownames(path) <- NULL
    paths <- rbind(paths, path)
  }
  paths
}


build_new_instances_deterministicly <-
  function(observation,
           n_trees,
           paths,
           possible_values,
           plausible_mask,
           fixed_variables_names,
           ordered_variables_names,
           step,
           verbose) {
    if (verbose) {
      pb <-
        txtProgressBar(
          min = 0,
          max = n_trees - 1,
          initial = 0,
          style = 3
        )
    }
    new_instances <- data.frame()
    for (i in seq_len(n_trees) - 1) {
      tree_paths <- paths[paths$Tree == i, ]
      for (j in seq_len(length(unique(tree_paths$Path)))) {
        path <- tree_paths[tree_paths$Path == j, ]
        x <- build_new_instance_from_path(
          observation,
          path,
          possible_values,
          plausible_mask,
          fixed_variables_names,
          ordered_variables_names,
          step
        )
        new_instances <- rbind(new_instances, x)
      }
      if (verbose)
        setTxtProgressBar(pb, i)
    }
    if (verbose)
      close(pb)
    new_instances
  }

build_new_instances_randomly <-
  function(new_observation,
           all_paths,
           all_paths_ids,
           possible_values,
           plausible_mask,
           fixed_variables_names,
           ordered_variables_names,
           step,
           paths_per_counterfactual,
           max_tries,
           verbose) {
    if (verbose) {
      pb <-
        txtProgressBar(
          min = 1,
          max = max_tries,
          initial = 0,
          style = 3
        )
    }
    new_instances <- data.frame()
    for (i in seq_len(max_tries)) {
      which_paths <- sample(nrow(all_paths_ids), paths_per_counterfactual)
      selected_paths <-
        all_paths[all_paths$Tree %in% all_paths_ids[which_paths, "Tree"] &
                    all_paths$Path %in% all_paths_ids[which_paths, "Path"], ]
      x <- new_observation
      for (j in seq_len(paths_per_counterfactual)) {
        tree_id <- all_paths_ids[which_paths[j], "Tree"]
        path_id <- all_paths_ids[which_paths[j], "Path"]
        path <-
          selected_paths[selected_paths$Tree == tree_id &
                           selected_paths$Path == path_id, ]
        x <- build_new_instance_from_path(
          x,
          path,
          possible_values,
          plausible_mask,
          fixed_variables_names,
          ordered_variables_names,
          step
        )
      }
      new_instances <- rbind(new_instances, x)
      if (verbose)
        setTxtProgressBar(pb, i)
    }
    if (verbose)
      close(pb)
    new_instances
  }


build_new_instance_from_path <- function(observation,
                                         path,
                                         possible_values,
                                         plausible_mask,
                                         fixed_variables_names,
                                         ordered_variables_names,
                                         step) {
  x <- observation
  for (j in seq_len(nrow(path) - 1)) {
    variable_name <- path$Feature[j]
    if (variable_name %in% fixed_variables_names)
      next
    x[variable_name] <-
      get_closest_proper_value(
        variable_name,
        x[, variable_name],
        path$Split[j],
        path$Decision[j],
        possible_values[[variable_name]],
        plausible_mask[[variable_name]],
        variable_name %in% ordered_variables_names,
        step
      )
  }
  x
}

get_closest_proper_value <- function(variable_name,
                                     current_value,
                                     split_value,
                                     decision_type,
                                     possible_values,
                                     plausible_mask,
                                     is_ordered,
                                     step = 1) {
  look_for_leq <- decision_type == "<="

  if (!is_ordered) {
    if (look_for_leq & (current_value > split_value)) {
      tmp <-
        possible_values[which((possible_values <= split_value) &
                                plausible_mask)]
    } else if (!look_for_leq & (current_value <= split_value)) {
      tmp <-
        possible_values[which((possible_values > split_value) &
                                plausible_mask)]
    } else {
      return(current_value)
    }
  } else {
    current_value_index <-
      which(possible_values == as.character(current_value))
    split_value_index <-
      floor(split_value) # sometimes split_value is not integer (.5 is added)
    if (look_for_leq & (current_value_index > split_value_index)) {
      tmp <-
        possible_values[1:split_value_index][plausible_mask[1:split_value_index]]
    } else if (!look_for_leq &
               (current_value_index <= split_value_index)) {
      k <- length(possible_values)
      tmp <-
        possible_values[(split_value_index + 1):k][plausible_mask[(split_value_index +
                                                                     1):k]]
    } else {
      return(current_value)
    }
  }
  if (is.null(tmp) | length(tmp) == 0)
    return(current_value) # if there are no proper values - return the current value
  n_possible <- length(tmp)
  if (look_for_leq) {
    return(tmp[max(1, n_possible - step + 1)]) # return the maximal proper value - closest to split_value
  } else {
    return(tmp[min(step, n_possible)]) # return the minimal proper value - closest to split_value
  }
}
