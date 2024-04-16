#' @title Plot time-dependent predictions
#'
#' @param preds a matrix of predictions
#' @param times a vector of time points at which the predictions were made
#' @param preds_labels a vector of labels for the predictions
#' @param alpha transparency of the lines
#' @param linewidth width of the lines
#'
#' @import ggplot2
#' @export
plot_predictions <- function(preds, times, preds_labels=NULL, alpha=0.3, linewidth=0.3) {
  if (is.null(preds_labels)) {
    preds_labels <- factor(rep("all", nrow(preds)))
  }

  if (all(apply(preds, 1, diff) <= 0))
    y_title <- "Survival probability"
  else
    y_title <- "Cumuative hazard"

  preds_plot <- prepare_predictions_to_plot(preds, times, preds_labels, alpha, linewidth)

  with(preds_plot,
       {
         ggplot(preds_plot, aes(x = time, y = fun, group = id, color = cluster)) +
           geom_line(aes(alpha = alpha, linewidth = linewidth)) +
           xlab("Time") +
           ylab(y_title) +
           theme_minimal() +
           scale_alpha_identity() +
           scale_linewidth_identity()
       }
  )
}

#' @title Plot parrallel coordinates plot for counterfactual explanations
#'
#' @param counterfactual_explanations a `counterfactual_explanations` object with counterfactual explanations
#' @param filtered_examples a data frame with counterfactual examples. If NULL, all counterfactual examples from `counterfactual_explanations` will be used
#' @param variables a numeric vector of column indices representing variables to plot. If NULL, all variables will be on the plot
#'
#' @export
plot_parallel_coordinates <- function(counterfactual_explanations, filtered_examples = NULL, variables = NULL){
  if (is.null(variables)) {
    variables <- 1:ncol(counterfactual_explanations$counterfactual_examples)
  } else if (is.character(variables)) {
    variables <- match(variables, colnames(counterfactual_explanations$counterfactual_examples))
  } else if (!is.numeric(variables)) {
    stop("variables must be a numeric vector or a character vector of column names")
  }

  if (is.null(filtered_examples)) {
    filtered_examples <- counterfactual_explanations$counterfactual_examples
    filtered_examples_validity <- counterfactual_explanations$objective_values$validity
  } else if (!all(colnames(counterfactual_explanations$counterfactual_examples) == colnames(filtered_examples[,1:ncol(counterfactual_explanations$counterfactual_examples)]))) {
    stop("filtered_examples must have the same columns as counterfactual_explanations$counterfactual_examples")
  } else {
    filtered_examples_validity <- filtered_examples$validity
    filtered_examples <- filtered_examples[,1:ncol(counterfactual_explanations$counterfactual_examples)]
  }
  numeric_mask <- sapply(counterfactual_explanations$original_observation, is.numeric)
  whole_population <- rbind(counterfactual_explanations$original_observation, counterfactual_explanations$counterfactual_examples)

  categorical_variables_orderings <- lapply(
    counterfactual_explanations$counterfactual_examples[,!numeric_mask, drop=FALSE],
    function(x) {
      factor_var <- as.factor(x)
      factor_levels <- levels(factor_var)
      setNames(1:length(factor_levels), factor_levels)
    }
  )

  whole_population[,!numeric_mask] <- lapply(
    whole_population[,!numeric_mask, drop=FALSE],
    function(x) as.numeric(as.factor(x))
  )

  filtered_examples[,!numeric_mask] <- lapply(
    filtered_examples[,!numeric_mask, drop=FALSE],
    function(x) as.numeric(as.factor(x))
  )

  original_obs <- counterfactual_explanations$original_observation
  original_obs[!numeric_mask] <- lapply(
    original_obs[,!numeric_mask, drop=FALSE],
    function(x) as.numeric(as.factor(x))
  )

  # scale min-max to [0, 1] range
  mins <- apply(whole_population, 2, min)
  maxs <- apply(whole_population, 2, max)
  nan_columns <- which(maxs-mins == 0)

  filtered_examples <- t((t(filtered_examples) - mins) / (maxs - mins))
  filtered_examples[,nan_columns] <- 0

  original_obs <- t((t(original_obs) - mins) / (maxs - mins))
  original_obs[,nan_columns] <- 0

  plot_df <- data.frame(rbind(original_obs, filtered_examples))
  plot_df <- plot_df[, variables]
  rownames(plot_df) <- NULL
  plot_df$validity <- c(0, filtered_examples_validity)
  plot_df$type <- c("original_observation", rep("counterfactual", nrow(filtered_examples)))
  plot_df$id <- 1:nrow(plot_df)

  mins[!numeric_mask] <- sapply(categorical_variables_orderings, function(x) names(x)[which.min(x)])
  maxs[!numeric_mask] <- sapply(categorical_variables_orderings, function(x) names(x)[which.max(x)])

  text_df <- data.frame(
    variable = rep(variables, 2),
    value = rep(c(-0.05, 1.05), each=length(variables)),
    label = c(mins[variables], maxs[variables])
  )

  long_plot_df <- reshape2::melt(plot_df, id.vars = c("validity", "type", "id"))

  with(long_plot_df,
       {
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
  )

}


#' @export
plot_changes_frequency <- function(counterfactual_explanations, filtered_examples = NULL, variables = NULL){
  if (is.null(variables)) {
    variables <- 1:ncol(counterfactual_explanations$counterfactual_examples)
  } else if (is.character(variables)) {
    variables <- match(variables, colnames(counterfactual_explanations$counterfactual_examples))
  } else if (!is.numeric(variables)) {
    stop("variables must be a numeric vector or a character vector of column names")
  }

  if (is.null(filtered_examples)) {
    filtered_examples <- counterfactual_explanations$counterfactual_examples
    filtered_examples_validity <- counterfactual_explanations$objective_values$validity
  } else if (!all(colnames(counterfactual_explanations$counterfactual_examples) == colnames(filtered_examples[,1:ncol(counterfactual_explanations$counterfactual_examples)]))) {
    stop("filtered_examples must have the same columns as counterfactual_explanations$counterfactual_examples")
  } else {
    filtered_examples_validity <- filtered_examples$validity
    filtered_examples <- filtered_examples[,1:ncol(counterfactual_explanations$counterfactual_examples)]
  }

  plot_df <- colMeans(data.frame(as.list(counterfactual_explanations$original_observation) != filtered_examples))
  plot_df <- plot_df[variables]
  long_plot_df <- reshape2::melt(plot_df)
  long_plot_df$variable <- rownames(long_plot_df)
  long_plot_df <- long_plot_df[long_plot_df$value > 0,]

  with(long_plot_df,
       {
         ggplot(long_plot_df, aes(y = reorder(variable, value), x = value, fill = value)) +
           geom_bar(stat = "identity", fill="darkorchid4", width = 0.7) +
           scale_x_continuous(limits=c(0, 1)) +
           theme_minimal() +
           labs(title = "Frequency of variable changes",
                y = "Variable",
                x = "Fraction of counterfactuals with changes")
       }
  )
}


#' @export
plot_counterfactual_predictions <- function(counterfactual_explanations,
                                            filtered_predictions = NULL,
                                            filtered_examples = NULL,
                                            function_type = NULL){
  if (is.null(filtered_examples)) {
    filtered_examples <- counterfactual_explanations$counterfactual_examples
    filtered_examples_objective_values <- counterfactual_explanations$objective_values
  } else {
    filtered_examples_objective_values <- filtered_examples[,(ncol(counterfactual_explanations$counterfactual_examples)+1):ncol(filtered_examples)]
    filtered_examples <- filtered_examples[,1:ncol(counterfactual_explanations$counterfactual_examples)]
  }
  if (is.null(filtered_predictions)) {
    filtered_predictions <- counterfactual_explanations$predictions
  }
  original_prediction <- counterfactual_explanations$original_prediction[1,]
  target_envelope <- counterfactual_explanations$target_envelope
  numeric_mask <- sapply(var_values, is.numeric)

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

  colnames(filtered_predictions) <- counterfactual_explanations$times
  plot_df <- cbind(filtered_examples,
                   id = rownames(filtered_examples),
                   filtered_predictions)

  plot_df <- reshape2::melt(plot_df,
                            id.vars = 1:(ncol(filtered_examples)+1),
                            value.name = "fun",
                            variable.name = "time")

  plot_df$time <- as.numeric(as.character(plot_df$time))

  var_values <- counterfactual_explanations$original_observation
  var_values[numeric_mask] <- round(var_values[numeric_mask], 3)

  var_values_str <- paste0("<i>", names(var_values), "</i>: ", var_values, "<br>", collapse="")

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
    tmp <- plot_df[plot_df$id == rownames(filtered_examples)[i],]
    loss_values <- filtered_examples_objective_values[i,]
    var_values <- filtered_examples[i,]
    var_values[numeric_mask] <- round(var_values[numeric_mask], 3)
    var_values_str <- paste0("<i>", names(var_values), "</i>: ", var_values, "<br>", collapse="")
    fig <- fig %>% add_trace(
      x = tmp$time,
      y = tmp$fun,
      text = rownames(filtered_examples)[i],
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
  } else {
    preds_plot$cluster <- factor(rep("all", nrow(preds_plot)))
  }
  preds_plot$alpha <- alpha
  preds_plot$linewidth <- linewidth
  preds_plot <- reshape2::melt(preds_plot, id.vars = c("id", "cluster", "alpha", "linewidth"),
                              value.name = "fun",
                              variable.name = "time")
  preds_plot$time <- as.numeric(as.character(preds_plot$time))
  return(preds_plot)
}
