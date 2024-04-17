library(ggplot2)

saveRDS(experiment1_results, "experiment1_results_v2.rds")

explainer <- experiment1_results$explainer

execution_time_df <- data.frame(
  time = c(experiment1_results$moc_time,
           experiment1_results$tb_time,
           experiment1_results$kov_time),
  method = c(rep("MOC", length(experiment1_results$moc_time)),
             rep("TB", length(experiment1_results$tb_time)),
             rep("KOV", length(experiment1_results$kov_time)))
)


ggplot(execution_time_df, aes(x = method, y = time)) +
  geom_boxplot() +
  xlab("Method") +
  ylab("Execution time (s)") +
  theme_minimal()


sapply(experiment1_results$moc_results,
       function(res){
         sum(res$objective_values$validity == 0)
       })

sapply(experiment1_results$tb_results,
       function(res){
         sum(res$objective_values$validity == 0)
       })



moc_all_ov <- do.call("rbind", lapply(experiment1_results$moc_results,
       function(res){
         res$objective_values
       }))
moc_all_ov$method <- "MOC"


tb_all_ov <- do.call("rbind", lapply(experiment1_results$tb_results,
       function(res){
         res$objective_values
       }))
tb_all_ov$method <- "TB"


data_range <- apply(explainer$data, 2, range)
kov_all_ov <- do.call("rbind", lapply(experiment1_results$kov_results,
       function(res){
         calculate_objective_values(explainer$data, explainer$times, experiment1_results$weights,
                                    res$original_observation, res$counterfactual_examples,
                                    res$original_prediction, res$predictions,
                                    NULL, experiment1_results$target_envelope_sf,
                                    data_range, NULL,
                                    1, 5)
       }))
kov_all_ov$method <- "Kovalev"


all_ov <- rbind(moc_all_ov, tb_all_ov, kov_all_ov)

all_ov <- reshape2::melt(all_ov, id.vars = c("method"))

ggplot(all_ov, aes(x = method, y = value)) +
  geom_boxplot() +
  xlab("Method") +
  facet_wrap(~variable, scales = "free_y") +
  theme_minimal()



all_kov_cfes <- do.call("rbind",
                        lapply(experiment1_results$kov_results,
                               function(res){
                                 res$counterfactual_examples
                               }
                        )
)

all_kov_cfes <- data.frame(all_kov_cfes)
all_kov_cfes <- all_kov_cfes[!duplicated(all_kov_cfes),]

all_kov_preds <- predict(explainer, all_kov_cfes)

distances <- distance_from_target_envelope(all_kov_preds,
                              experiment1_results$target_envelope_sf,
                              explainer$times)
sum(distances == 0)

plot(all_kov_preds[which.min(distances),], type="l")
lines(experiment1_results$target_envelope_sf$lower_bound, col="red")
lines(experiment1_results$target_envelope_sf$upper_bound, col="red")

all_kov_cfes[which.min(distances),]

all_kov

