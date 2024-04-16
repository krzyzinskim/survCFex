#' @export
get_cluster_envelope <- function(preds, clusters, cluster_id, q = 0.05){
  group_preds <- preds[clusters == cluster_id,]
  lower_bound <- apply(group_preds, 2, function(x) quantile(x, probs = q))
  upper_bound <- apply(group_preds, 2, function(x) quantile(x, probs = 1 - q))
  # check if preds are survival functions (non increasing)
  if (all(apply(preds, 1, diff) <= 0))
    type <- "survival"
  else
    type <- "hazard"

  res <- list(
    "lower_bound" = lower_bound,
    "upper_bound" = upper_bound,
    "type" = type)
  class(res) <- "target_envelope"
  return(res)
}

#' @export
translate_target_envelope <- function(target_envelope){
  stopifnot("target_envelope should be of class 'target_envelope'" = class(target_envelope) == "target_envelope")
  stopifnot("target_envelope should have a type" = "type" %in% names(target_envelope))
  if (target_envelope$type == "survival"){
    new_target_envelope <- target_envelope
    new_target_envelope$upper_bound <- translate_survival_to_cumulative_hazard(target_envelope$lower_bound)
    new_target_envelope$lower_bound <- translate_survival_to_cumulative_hazard(target_envelope$upper_bound)
    new_target_envelope$type <- "cumulative_hazard"
  } else {
    new_target_envelope <- target_envelope
    new_target_envelope$upper_bound <- translate_cumulative_hazard_to_survival(target_envelope$lower_bound)
    new_target_envelope$lower_bound <- translate_cumulative_hazard_to_survival(target_envelope$upper_bound)
    new_target_envelope$type <- "survival"
  }
  return(new_target_envelope)
}


#' @export
plot_cluster_envelopes <- function(preds, clusters, times, q = 0.05, preds_to_plot = NULL, alpha=0.2) {
  if (all(apply(preds, 1, diff) <= 0))
    y_title <- "Survival probability"
  else
    y_title <- "Cumuative hazard"

  n_clusters <- max(clusters)
  envelopes_df <- data.frame()
  for (i in 1:n_clusters) {
    envelope <- get_cluster_envelope(preds, clusters, i, q)
    envelope_df <- data.frame(time = times,
                              lower_bound = envelope$lower_bound,
                              upper_bound = envelope$upper_bound,
                              cluster = i)
    envelopes_df <- rbind(envelopes_df, envelope_df)
  }

  if (!is.null(preds_to_plot)) {
    if (is.vector(preds_to_plot)) {
      preds_to_plot <- matrix(preds_to_plot, nrow = 1)
    }
    preds_plot <- prepare_predictions_to_plot(preds_to_plot, times)
  }

  p <- ggplot(envelopes_df, aes(x = time, group = cluster, color = factor(cluster), fill = factor(cluster))) +
    geom_ribbon(aes(ymin = lower_bound, ymax = upper_bound), alpha=alpha, colour=NA) +
    # geom_line(aes(y = lower_bound)) +
    # geom_line(aes(y = upper_bound)) +
    xlab("Time") +
    ylab(y_title) +
    theme_minimal() +
    scale_fill_brewer(type = "qual", palette = "Set2", name = "Cluster") +
    scale_color_brewer(type = "qual", palette = "Set2", name = "Cluster")

  if (!is.null(preds_to_plot)) {
    p <- p + geom_line(data = preds_plot, aes(x = time, y = fun, group = id),
                       alpha=1, linewidth=0.5, inherit.aes = FALSE)
  }
  p

}
