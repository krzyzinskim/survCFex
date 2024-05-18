#' @export
survival_weights <- function(explainer, times=explainer$times, p=0, q=0, add_mean = FALSE){
  if (p == 0 & q == 0){
    weights <- rep(1, length(times))
  } else {
    km <- survival::survfit(explainer$y ~ 1)
    surv_estimator <- stepfun(km$time, c(1, km$surv))
    sf <- surv_estimator(times)
    weights <- (sf^p) * ((1 - sf)^q)
    weights <- weights + mean(weights) * add_mean
  }
  return(weights)
}


#' @export
plot_survival_weights <- function(explainer, times, p=0, q=0, add_mean = FALSE){
  stopifnot(length(p) == length(q))
  weights_to_plot_df <- data.frame()

  for (i in 1:length(p)) {
    weights <- survival_weights(explainer, times, p[i], q[i], add_mean)
    weights_df <- data.frame(time = times, weight = weights, params = paste0("p=", p[i], ", q=", q[i]))
    weights_to_plot_df <- rbind(weights_to_plot_df, weights_df)
  }

  with(weights_to_plot_df,
       {
         ggplot(weights_to_plot_df, aes(x = time, y = weight, group = params, color = params)) +
           geom_line(linewidth=0.8) +
           xlab("Time") +
           ylab("Weight") +
           theme_minimal() +
           scale_color_brewer(type = "qual", palette = "Set1") +
           guides(color = guide_legend(title = "Parameters", nrow=2)) +
           theme(legend.position = "bottom")
       }
  )
}
