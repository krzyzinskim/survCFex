#' @export
plot_predictions <- function(preds, times, preds_labels=NULL, alpha=0.3, linewidth=0.3) {
  if (is.null(preds_labels)) {
    preds_labels <- factor(rep("none", nrow(preds)))
  }


  if (all(apply(preds, 1, diff) <= 0))
    y_title <- "Survival probability"
  else
    y_title <- "Cumuative hazard"


  preds_plot <- prepare_predictions_to_plot(preds, times, preds_labels, alpha, linewidth)

  ggplot(preds_plot, aes(x = time, y = fun, group = id, color = cluster)) +
    geom_line(aes(alpha = alpha, linewidth = linewidth)) +
    xlab("Time") +
    ylab(y_title) +
    theme_minimal() +
    scale_alpha_identity() +
    scale_linewidth_identity()
}




#' @export
plot_parallel_coordinates <- function(counterfactual_explanations, filtered_population = NULL, variables = NULL){
  if (is.null(variables)) {
    variables <- 1:ncol(counterfactual_explanations$counterfactual_examples)
  } else if (is.character(variables)) {
    variables <- match(variables, colnames(counterfactual_explanations$counterfactual_examples))
  } else if (!is.numeric(variables)) {
    stop("variables must be a numeric vector or a character vector of column names")
  }

  if (is.null(filtered_population)) {
    filtered_population <- counterfactual_explanations$counterfactual_examples
    filtered_population_validity <- counterfactual_explanations$objective_values$validity
  } else if (!all(colnames(counterfactual_explanations$counterfactual_examples) == colnames(filtered_population[,1:ncol(counterfactual_explanations$counterfactual_examples)]))) {
    stop("filtered_population must have the same columns as counterfactual_explanations$counterfactual_examples")
  } else {
    filtered_population_validity <- filtered_population$validity
    filtered_population <- filtered_population[,1:ncol(counterfactual_explanations$counterfactual_examples)]
  }

  whole_population <- rbind(counterfactual_explanations$original_observation, counterfactual_explanations$counterfactual_examples)

  # scale min-max to [0, 1] range
  mins <- apply(whole_population, 2, min)
  maxs <- apply(whole_population, 2, max)

  filtered_population <- t((t(filtered_population) - mins) / (maxs - mins))
  nan_columns <- which(maxs-mins == 0)
  filtered_population[, nan_columns] <- 0

  original_obs <- t((t(counterfactual_explanations$original_observation) - mins) / (maxs - mins))
  original_obs[, nan_columns] <- 0

  plot_df <- data.frame(rbind(original_obs, filtered_population))
  plot_df <- plot_df[, variables]
  rownames(plot_df) <- NULL
  plot_df$validity <- c(0, filtered_population_validity)
  plot_df$type <- c("original_observation", rep("counterfactual", nrow(filtered_population)))
  plot_df$id <- 1:nrow(plot_df)

  text_df <- data.frame(
    variable = rep(variables, 2),
    value = rep(c(-0.05, 1.05), each=length(variables)),
    label = c(mins[variables], maxs[variables])
  )

  long_plot_df <- reshape2::melt(plot_df, id.vars = c("validity", "type", "id"))

  ggplot(long_plot_df[long_plot_df$type == "counterfactual", ],
         aes(x = variable, y = value, group=id, color = validity)) +
    geom_line(linewidth = 0.5, alpha=0.4) +
    geom_line(data = long_plot_df[long_plot_df$type == "original_observation", ],
              aes(x = variable, y = value, group=id),
              color = "mediumvioletred", linewidth = 1) +
    scale_color_distiller(palette = "Purples") +
    theme_minimal() +
    labs(title = "Parallel Coordinates Plot",
         x = "Variable",
         y = "Scaled variable value",
         caption = "Variable values are scaled to [0, 1] range.") +
    theme(legend.position = "bottom",
          plot.caption = element_text(hjust = 0.5, size = 7)) +
    geom_text(data = text_df, inherit.aes = FALSE,
              aes(x = variable, y = value, label = label), size=3)
}

#' @export
plot_changes_frequency <- function(counterfactual_explanations, filtered_population = NULL, variables = NULL){
  if (is.null(variables)) {
    variables <- 1:ncol(counterfactual_explanations$counterfactual_examples)
  } else if (is.character(variables)) {
    variables <- match(variables, colnames(counterfactual_explanations$counterfactual_examples))
  } else if (!is.numeric(variables)) {
    stop("variables must be a numeric vector or a character vector of column names")
  }

  if (is.null(filtered_population)) {
    filtered_population <- counterfactual_explanations$counterfactual_examples
    filtered_population_validity <- counterfactual_explanations$objective_values$validity
  } else if (!all(colnames(counterfactual_explanations$counterfactual_examples) == colnames(filtered_population[,1:ncol(counterfactual_explanations$counterfactual_examples)]))) {
    stop("filtered_population must have the same columns as counterfactual_explanations$counterfactual_examples")
  } else {
    filtered_population_validity <- filtered_population$validity
    filtered_population <- filtered_population[,1:ncol(counterfactual_explanations$counterfactual_examples)]
  }

  plot_df <- colMeans(data.frame(as.list(counterfactual_explanations$original_observation) != filtered_population)) * 100
  plot_df <- plot_df[variables]
  long_plot_df <- reshape2::melt(plot_df)
  long_plot_df$variable <- rownames(long_plot_df)
  long_plot_df <- long_plot_df[long_plot_df$value > 0,]

  ggplot(long_plot_df, aes(y = reorder(variable, value), x = value, fill = value)) +
    geom_bar(stat = "identity", fill="darkorchid4", width = 0.7) +
    scale_x_continuous(labels = scales::percent_format(scale = 1),
                       expand = c(0, 2),
                       limits=c(0, 100)) +
    theme_minimal() +
    labs(title = "Frequency of variable changes",
         y = "Variable",
         x = "Percentage of counterfactuals with changes")
}

#' @export
plot_counterfactual_predictions <- function(counterfactual_explanations,
                                            filtered_predictions = NULL,
                                            filtered_population = NULL,
                                            function_type = NULL){
  if (is.null(filtered_population)) {
    filtered_population <- counterfactual_explanations$counterfactual_examples
    filtered_population_objective_values <- counterfactual_explanations$objective_values
  } else {
    filtered_population_objective_values <- filtered_population[,(ncol(counterfactual_explanations$counterfactual_examples)+1):ncol(filtered_population)]
    filtered_population <- filtered_population[,1:ncol(counterfactual_explanations$counterfactual_examples)]
  }
  if (is.null(filtered_predictions)) {
    filtered_predictions <- counterfactual_explanations$predictions
  }
  original_prediction <- counterfactual_explanations$original_prediction[1,]
  target_envelope <- counterfactual_explanations$target_envelope

  counterfactuals_type <- class(counterfactual_explanations)[2]
  if (counterfactuals_type == "multiobjective_counterfactuals"){
    current_function_type <- "survival"
  } else if (counterfactuals_type == "treebased_counterfactuals"){
    current_function_type <- "chf"
  } else {
    stop("Unknown counterfactuals type")
  }

  if (is.null(function_type)){
    function_type <- current_function_type
  } else if (function_type != current_function_type){
    if (counterfactuals_type == "multiobjective_counterfactuals"){
      filtered_predictions <- translate_survival_to_cumulative_hazard(filtered_predictions)
      original_prediction <- translate_survival_to_cumulative_hazard(original_prediction)
    } else if (counterfactuals_type == "treebased_counterfactuals"){
      filtered_predictions <- translate_cumulative_hazard_to_survival(filtered_predictions)
      original_prediction <- translate_cumulative_hazard_to_survival(original_prediction)
    }
    target_envelope <- translate_target_envelope(target_envelope)
  }

  function_names <- c("survival" = "Survival probability",
                      "chf" = "Cumulative hazard")
  y_title <- function_names[function_type]

  colnames(filtered_predictions) <- paste0("T=", counterfactual_explanations$times)
  plot_df <- cbind(filtered_population,
                   id = rownames(filtered_population),
                   filtered_predictions)

  plot_df <- pivot_longer(plot_df,
                          cols = starts_with("T="),
                          names_prefix = "T=",
                          names_to = "time",
                          values_to = "fun")
  plot_df$time <- as.numeric(plot_df$time)

  var_values <- counterfactual_explanations$original_observation
  var_values_str <- paste0("<i>", names(var_values), "</i>: ", round(var_values, 3), "<br>", collapse="")

  fig <- plot_ly(
    x = counterfactual_explanations$times,
    y = original_prediction,
    type = "scatter",
    mode = "lines",
    line = list(width = 1.5, color = "mediumvioletred"),
    hovertemplate =
      paste0("<b>Time:</b> %{x}<br>",
             "<b>", y_title, "</b> %{y:.3f}<br>",
             "<b>Variable values:</b><br>",
             var_values_str,
             "<extra><b>Original observation</b></extra>")
  ) %>%
    layout(title = "Counterfactual predictions",
           xaxis = list(title = "Time"),
           yaxis = list(title = y_title))

  fig <- fig %>% add_ribbons(
    x = counterfactual_explanations$times,
    ymin = target_envelope$lower_bound,
    ymax = target_envelope$upper_bound,
    fillcolor = "rgba(155, 50, 204, 0.25)",
    line = list(width = 0, color="darkorchid"),
    showlegend = FALSE,
    text = "Target envelope",
    hovertemplate =
      paste0("<b>Time:</b> %{x}<br>",
             "<b>", y_title, "</b> %{y:.3f}<br>",
             "<extra>Target envelope</extra>")
  )

  for (i in 1:nrow(filtered_predictions)) {
    tmp <- plot_df[plot_df$id == rownames(filtered_population)[i],]
    loss_values <- filtered_population_objective_values[i,]
    var_values <- filtered_population[i,]
    var_values_str <- paste0("<i>", names(var_values), "</i>: ", round(var_values, 3), "<br>", collapse="")
    fig <- fig %>% add_trace(
      x = tmp$time,
      y = tmp$fun,
      text = rownames(filtered_population)[i],
      type = "scatter",
      mode = "lines",
      line = list(width = 0.7, color = "darkorchid"),
      showlegend = FALSE,
      hovertemplate =
        paste0("<b>Time:</b> %{x}<br>",
               "<b>", y_title, "</b> %{y:.3f}<br>",
               "<b>Loss values:</b> <br>",
               "<i>Validity:</i> ", round(loss_values["validity"], 3), "<br>",
               "<i>Similarity:</i> ", round(loss_values["similarity"], 3), "<br>",
               "<i>Plausibility:</i> ", round(loss_values["plausibility"], 3), "<br>",
               "<i>Sparsity:</i> ", loss_values["sparsity"], "<br>",
               "<b>Variable values:</b><br>",
               var_values_str,
               "<extra>Counterfactual %{text}</extra>"
        )
    )
  }
  fig <- config(fig, displayModeBar = FALSE)

  fig
}

prepare_predictions_to_plot <- function(preds, times, preds_labels=NULL, alpha=0.3, linewidth=0.3) {
  preds_plot <- as.data.frame(preds)
  colnames(preds_plot) <- times
  preds_plot$id <- rownames(preds_plot)
  if (!is.null(preds_labels)) {
    preds_plot$cluster <- as.factor(preds_labels)
  }
  preds_plot$alpha <- alpha
  preds_plot$linewidth <- linewidth
  preds_plot <- tidyr::pivot_longer(preds_plot,
                                    cols = -any_of(c("id", "cluster", "alpha", "linewidth")),
                                    names_to = "time",
                                    values_to = "fun")
  preds_plot$time <- as.numeric(preds_plot$time)
  return(preds_plot)
}

