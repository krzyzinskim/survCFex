#' @export
get_prediction_band <- function(x, clusters, cluster_id, lambda = 0.05, output_type = "survival"){
  if (class(x)[1] == "surv_explainer"){
    preds <- predict(explainer, output_type = output_type)
  }
  else {
    preds <- x
  }
  group_preds <- preds[clusters == cluster_id,]
  lower_bound <- apply(group_preds, 2, function(x) quantile(x, probs = lambda))
  upper_bound <- apply(group_preds, 2, function(x) quantile(x, probs = 1 - lambda))
  # check if preds are survival functions (non increasing)
  res <- list(
    "lower_bound" = lower_bound,
    "upper_bound" = upper_bound,
    "output_type" = output_type)
  class(res) <- "target_prediction_band"
  return(res)
}


#' @export
translate_target_prediction_band <- function(target_prediction_band){
  stopifnot("target_prediction_band should be of class 'target_prediction_band'" = class(target_prediction_band) == "target_prediction_band")
  stopifnot("target_prediction_band should have a output_type" = "output_type" %in% names(target_prediction_band))
  if (target_prediction_band$output_type == "survival"){
    new_target_prediction_band <- target_prediction_band
    new_target_prediction_band$upper_bound <- translate_survival_to_cumulative_hazard(target_prediction_band$lower_bound)
    new_target_prediction_band$lower_bound <- translate_survival_to_cumulative_hazard(target_prediction_band$upper_bound)
    new_target_prediction_band$output_type <- "cumulative_hazard"
  } else {
    new_target_prediction_band <- target_prediction_band
    new_target_prediction_band$upper_bound <- translate_cumulative_hazard_to_survival(target_prediction_band$lower_bound)
    new_target_prediction_band$lower_bound <- translate_cumulative_hazard_to_survival(target_prediction_band$upper_bound)
    new_target_prediction_band$output_type <- "survival"
  }
  return(new_target_prediction_band)
}


#' @export
plot_prediction_bands <- function(x, clusters, lambda = 0.05, preds_to_plot = NULL, alpha=0.2, output_type="survival") {
  preds <- predict(explainer, output_type = output_type)
  if (output_type == "survival")
    y_title <- "Survival probability"
  else
    y_title <- "Cumulative hazard"

  n_clusters <- max(clusters)
  prediction_bands_df <- data.frame()
  for (i in 1:n_clusters) {
    prediction_band <- get_prediction_band(preds, clusters, i, lambda)
    prediction_band_df <- data.frame(time = explainer$times,
                              lower_bound = prediction_band$lower_bound,
                              upper_bound = prediction_band$upper_bound,
                              cluster = i)
    prediction_bands_df <- rbind(prediction_bands_df, prediction_band_df)
  }

  if (!is.null(preds_to_plot)) {
    if (is.vector(preds_to_plot)) {
      preds_to_plot <- matrix(preds_to_plot, nrow = 1)
    }
    preds_plot <- prepare_predictions_to_plot(preds_to_plot, times)
  }

  p <- with(prediction_bands_df,
       {
          ggplot(prediction_bands_df, aes(x = time, group = cluster, color = factor(cluster), fill = factor(cluster))) +
           geom_ribbon(aes(ymin = lower_bound, ymax = upper_bound), alpha=alpha, colour=NA) +
           xlab("Time") +
           ylab(y_title) +
           theme_minimal() +
           scale_fill_brewer(type = "qual", palette = "Set2", name = "Cluster") +
           scale_color_brewer(type = "qual", palette = "Set2", name = "Cluster")
       }
  )

  if (!is.null(preds_to_plot)) {
    p <-with(preds_plot,
         {
            p + geom_line(data = preds_plot, aes(x = time, y = fun, group = id),
                              alpha=1, linewidth=0.5, inherit.aes = FALSE)
         }
    )
  }

  p
}


check_target_prediction_band <- function(target_prediction_band, times){
  # check class
  stopifnot("target_prediction_band should be of class 'target_prediction_band'" = class(target_prediction_band) == "target_prediction_band")
  # check if lower_bound and upper_bound are of the same length
  stopifnot("lower_bound and upper_bound should have the same length" = length(target_prediction_band$lower_bound) == length(target_prediction_band$upper_bound))
  # check if lower_bound and upper_bound have the same length as times
  stopifnot("lower_bound and upper_bound should have the same length as times" = length(target_prediction_band$lower_bound) == length(times))
  # check if lower_bound and upper_bound are numeric
  stopifnot("lower_bound should be numeric" = is.numeric(target_prediction_band$lower_bound))
  stopifnot("upper_bound should be numeric" = is.numeric(target_prediction_band$upper_bound))
  # check if lower_bound is smaller than upper_bound
  stopifnot("lower_bound should be smaller than upper_bound" = all(target_prediction_band$lower_bound <= target_prediction_band$upper_bound))
  # check if type is either "survival" or "chf"
  stopifnot("output_type should be either 'survival' or 'chf'" = target_prediction_band$output_type %in% c("survival", "chf"))
}
