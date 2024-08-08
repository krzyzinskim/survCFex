#saveRDS(experiment6_results, "experiments/results/experiment6_results.rds")

experiment6_results <- readRDS("experiments/results/experiment6_results.rds")

library(ggplot2)

explainer <- experiment6_results$explainer


# TIME
execution_time_df <- data.frame(
  time = c(experiment6_results$moc_time,
           experiment6_results$tb_time,
           experiment6_results$kov_time),
  method = c(rep("MOC", length(experiment6_results$moc_time)),
             rep("TB", length(experiment6_results$tb_time)),
             rep("KOV", length(experiment6_results$kov_time)))
)



library(dplyr)

execution_time_df %>%
  group_by(method) %>%
  summarise(mean_time = mean(time),
            sd_time = sd(time),
            min_time = min(time),
            max_time = max(time))


ggplot(execution_time_df, aes(x = method, y = time)) +
  geom_boxplot(staplewidth = 0.5,
               outlier.shape = 1,
               outlier.size = 0.5,
               fill="darkorchid") +
  # add mean
  stat_summary(fun = mean, geom = "point", shape = 23, size = 2, fill = "violet") +
  xlab("Method") +
  ylab("Execution time (s)") +
  theme_bw()

ggsave("experiments/plots/exp1b_execution_time.pdf", dpi=500,
       width=4.5, height=3, units="in")




# OBJECTIVE VALUES
moc_all_ov <- do.call("rbind", lapply(experiment6_results$moc_results,
                                      function(res){
                                        crowded_comparison_order <- get_crowded_comparison_order(res$objective_values)
                                        res$objective_values[select_population_indices(crowded_comparison_order, 10),]
                                      }))
moc_all_ov$method <- "SurvMOC"


data_range <- apply(explainer$data, 2, range)
kov_all_ov <- do.call("rbind", lapply(experiment6_results$kov_results,
                                      function(res){
                                        raw_obj <- calculate_objective_values(explainer$data, explainer$times, experiment6_results$weights,
                                                                              res$original_observation, res$counterfactual_examples,
                                                                              res$original_prediction, res$predictions,
                                                                              NULL, experiment6_results$target_envelope_sf,
                                                                              data_range, NULL,
                                                                              1, 5)
                                        crowded_comparison_order <- get_crowded_comparison_order(raw_obj)
                                        raw_obj[select_population_indices(crowded_comparison_order, 10),]
                                      }))
kov_all_ov$method <- "Kovalev"


all_ov <- rbind(moc_all_ov, kov_all_ov)

all_ov <- reshape2::melt(all_ov, id.vars = c("method"))


all_ov$variable <- factor(paste(all_ov$variable, "loss", sep = " "),
                          levels = c("validity loss", "similarity loss", "sparsity loss", "plausibility loss"))



ggplot(all_ov, aes(x = method, y = value)) +
  geom_boxplot(staplewidth = 0.5,
               outlier.shape = 1,
               outlier.size = 0.5,
               fill="darkorchid") +
  # add mean
  stat_summary(fun = mean, geom = "point", shape = 23, size = 2, fill = "violet") +
  xlab("Method") +
  facet_wrap(~variable, scales = "free_y", nrow = 1) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave("experiments/plots/exp1b_objective_values.pdf", dpi=500,
       width=8, height=3.2, units="in")







z# HOW MANY COUNTERFACTUALS ARE NOT-DOMINATED?
clusters <- get_clusters(experiment6_results$dendrogram, k=4)[experiment6_results$not_in_target_cluster_ids][experiment6_results$sample_for_evaluation]

moc_n_nd <- sapply(experiment6_results$moc_results,
                   function(res){
                     crowded_comparison_order <- get_crowded_comparison_order(res$objective_values)
                     obj <- res$objective_values[select_population_indices(crowded_comparison_order, 10),]
                     sum(ecr::doNondominatedSorting(as.matrix(t(obj)))$ranks == 1)
                   })


kov_n_nd <- sapply(experiment6_results$kov_results,
                   function(res){
                     raw_obj <- calculate_objective_values(explainer$data, explainer$times, experiment6_results$weights,
                                                           res$original_observation, res$counterfactual_examples,
                                                           res$original_prediction, res$predictions,
                                                           NULL, experiment6_results$target_envelope_sf,
                                                           data_range, NULL,
                                                           1, 5)
                     crowded_comparison_order <- get_crowded_comparison_order(raw_obj)
                     raw_obj <- raw_obj[select_population_indices(crowded_comparison_order, 10),]
                     sum(ecr::doNondominatedSorting(as.matrix(t(raw_obj)))$ranks == 1)
                   })




nondominated_df <- data.frame(
  method = c(rep("SurvMOC", length(moc_n_nd)),
             rep("Kovalev", length(kov_n_nd))),
  n_nondominated = c(moc_n_nd, kov_n_nd),
  obs_id = rep(1:length(moc_n_nd), 2),
  obs_cluster = rep(clusters, 2)
)


p1 <- ggplot(nondominated_df, aes(x = method, y = n_nondominated)) +
  geom_boxplot(staplewidth = 0.5,
               outlier.shape = 1,
               outlier.size = 0.5,
               fill="darkorchid") +
  stat_summary(fun = mean, geom = "point", shape = 23, size = 2, fill = "violet") +
  xlab("Method") +
  ylab("Number of non-dominated\ncounterfactual examples") +
  theme_bw()
p1

p2 <- ggplot(nondominated_df, aes(x = obs_id, y = n_nondominated, color = factor(obs_cluster))) +
  geom_point() +
  xlab("Observation ID") +
  ylab("Number of non-dominated\ncounterfactual examples in top-40") +
  theme_bw() +
  theme(legend.position = "bottom") +
  facet_wrap(~method, ncol = 1) +
  scale_color_brewer(palette = "Set2", type = "qual", name = "Cluster") +
  theme(legend.position = "right")
p2

library(ggpubr)
ggarrange(p1, p2, ncol = 2, nrow = 1, widths = c(1, 2))
ggsave("experiments/plots/exp1b_nondomination.pdf", dpi=500,
       width=4.5, height=3, units="in")


# valid_df
# nondominated_df
# merged_df <- merge(valid_df, nondominated_df, by=c("method", "obs_id", "obs_cluster"))
#
# ggplot(merged_df, aes(y = n_valid, x = n_nondominated, color = factor(obs_cluster))) +
#   geom_point() +
#   ylab("Number of valid counterfactual examples") +
#   xlab("Number of non-dominated counterfactual examples") +
#   theme_bw() +
#   theme(legend.position = "bottom") +
#   facet_wrap(~method, ncol = 1) +
#   scale_color_brewer(palette = "Set2", type = "qual", name = "Cluster") +
#   theme(legend.position = "right")
#



# COVERAGE RATES
cov_res <- lapply(1:50, function(i){
  moc_obj <- experiment6_results$moc_results[[i]]$objective_values
  kov_obj <- calculate_objective_values(explainer$data, explainer$times, experiment6_results$weights,
                                        experiment6_results$kov_results[[i]]$original_observation,
                                        experiment6_results$kov_results[[i]]$counterfactual_examples,
                                        experiment6_results$kov_results[[i]]$original_prediction,
                                        experiment6_results$kov_results[[i]]$predictions,
                                        NULL, experiment6_results$target_envelope_sf,
                                        data_range, NULL,
                                        1, 5)
  crowded_comparison_order <- get_crowded_comparison_order(kov_obj)
  kov_obj <- kov_obj[select_population_indices(crowded_comparison_order, 40),]

  moc_kov_obj <- rbind(moc_obj, kov_obj)
  c(mean(ecr::doNondominatedSorting(as.matrix(t(moc_kov_obj)))$ranks[1:40] != 1),
    mean(ecr::doNondominatedSorting(as.matrix(t(moc_kov_obj)))$ranks[41:80] != 1))
})

coverage_df <- matrix(unlist(cov_res), ncol=2, byrow=TRUE)

library(dplyr)
apply(coverage_df, 2, mean)
apply(coverage_df, 2, sd)
apply(coverage_df, 2, min)
apply(coverage_df, 2, max)


i <- 30
plot(explainer$data$x3, explainer$data$x4,
     xlim = range(explainer$data$x3),
     ylim = range(explainer$data$x4),
     col = scales::alpha("black", 0.2), pch=16)


crowded_comparison_order <- get_crowded_comparison_order(experiment6_results$moc_results[[i]]$objective_values)
sel_ces <- select_population_indices(crowded_comparison_order, 10)
moc_results <- experiment6_results$moc_results[[i]]$counterfactual_examples[sel_ces,]


res <- experiment6_results$kov_results[[i]]
raw_obj <- calculate_objective_values(explainer$data, explainer$times, experiment6_results$weights,
                                      res$original_observation, res$counterfactual_examples,
                                      res$original_prediction, res$predictions,
                                      NULL, experiment6_results$target_envelope_sf,
                                      data_range, NULL,
                                      1, 5)
crowded_comparison_order <- get_crowded_comparison_order(raw_obj)
kov_results <- res$counterfactual_examples[select_population_indices(crowded_comparison_order, 10),]



points(moc_results$x3,
       moc_results$x4,
       col="blue",
       pch=16)

points(kov_results[,3],
       kov_results[,4],
       col="red",
       pch=16)

points(experiment6_results$kov_results[[i]]$original_observation[3],
       experiment6_results$kov_results[[i]]$original_observation[4],
       col="black", pch=3, c=3)




evaluation_sample <- explainer$data[experiment6_results$not_in_target_cluster_ids,][experiment6_results$sample_for_evaluation,]
which_has_1 <- which(evaluation_sample$x1 == 1)
which_has_0 <- which(evaluation_sample$x1 == 0)
mask_ones <- evaluation_sample$x1 == 1

moc_x1 <- lapply(experiment6_results$moc_results,
                 function(res){
                   crowded_comparison_order <- get_crowded_comparison_order(res$objective_values)
                   sel_ces <- select_population_indices(crowded_comparison_order, 10)
                   ces <- res$counterfactual_examples[sel_ces,]
                   ces$x1
                 })


tb_x1 <- lapply(experiment6_results$tb_results,
                function(res){
                  crowded_comparison_order <- get_crowded_comparison_order(res$objective_values)
                  sel_ces <- select_population_indices(crowded_comparison_order, 10)
                  ces <- res$counterfactual_examples[sel_ces,]
                  ces$x1
                })

kov_x1 <- lapply(experiment6_results$kov_results,
                 function(res){
                   raw_obj <- calculate_objective_values(explainer$data, explainer$times, experiment6_results$weights,
                                                         res$original_observation, res$counterfactual_examples,
                                                         res$original_prediction, res$predictions,
                                                         NULL, experiment6_results$target_envelope_sf,
                                                         data_range, NULL,
                                                         1, 5)
                   crowded_comparison_order <- get_crowded_comparison_order(raw_obj)
                   ces <- res$counterfactual_examples[select_population_indices(crowded_comparison_order, 10),]
                   ces[,1]
                 })



moc_ones <- sapply(moc_x1, function(x) sum(x == 1))
tb_ones <- sapply(tb_x1, function(x) sum(x == 1))
kov_ones_real <- sapply(kov_x1, function(x) sum(x == 1))
kov_ones <- sapply(kov_x1, function(x) sum(abs(x-1)<abs(x)))

table(moc_ones == 10, mask_ones)
table(tb_ones == 40, mask_ones)
table(kov_ones == 10, mask_ones)


# table(moc_ones == 40, mask_ones)
# mask_ones
# FALSE TRUE
# TRUE    27   23
#
# table(tb_ones == 40, mask_ones)
# mask_ones
# FALSE TRUE
# FALSE     2    2
# TRUE     25   21
#
# table(kov_ones == 40, mask_ones)
# mask_ones
# FALSE TRUE
# FALSE    27    5
# TRUE      0   18

