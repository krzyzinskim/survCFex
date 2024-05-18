library(ggplot2)

saveRDS(experiment1_results, "experiments/results/experiment1_results.rds")
experiment1_results <- readRDS("experiments/results/experiment1_results.rds")


explainer <- experiment1_results$explainer


# TIME
execution_time_df <- data.frame(
  time = c(experiment1_results$moc_time,
           experiment1_results$tb_time,
           experiment1_results$kov_time),
  method = c(rep("MOC", length(experiment1_results$moc_time)),
             rep("TB", length(experiment1_results$tb_time)),
             rep("KOV", length(experiment1_results$kov_time)))
)

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

ggsave("experiments/plots/exp1_execution_time.pdf", dpi=500,
       width=4, height=4, units="in")




# OBJECTIVE VALUES
moc_all_ov <- do.call("rbind", lapply(experiment1_results$moc_results,
       function(res){
         crowded_comparison_order <- get_crowded_comparison_order(res$objective_values)
         res$objective_values[select_population_indices(crowded_comparison_order, 10),]
       }))
moc_all_ov$method <- "SurvMOC"

tb_all_ov <- do.call("rbind", lapply(experiment1_results$tb_results,
       function(res){
         crowded_comparison_order <- get_crowded_comparison_order(res$objective_values)
         res$objective_values[select_population_indices(crowded_comparison_order, 10),]
       }))
tb_all_ov$method <- "RSF-PT"

nrow(tb_all_ov)
tb_all_ov <- tb_all_ov[tb_all_ov$validity < 0.1, ]
nrow(tb_all_ov)

data_range <- apply(explainer$data, 2, range)
kov_all_ov <- do.call("rbind", lapply(experiment1_results$kov_results,
       function(res){
         raw_obj <- calculate_objective_values(explainer$data, explainer$times, experiment1_results$weights,
                                    res$original_observation, res$counterfactual_examples,
                                    res$original_prediction, res$predictions,
                                    NULL, experiment1_results$target_envelope_sf,
                                    data_range, NULL,
                                    1, 5)
         crowded_comparison_order <- get_crowded_comparison_order(raw_obj)
         raw_obj[select_population_indices(crowded_comparison_order, 10),]
       }))
kov_all_ov$method <- "Kovalev"


all_ov <- rbind(moc_all_ov, tb_all_ov, kov_all_ov)

all_ov <- reshape2::melt(all_ov, id.vars = c("method"))

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
ggsave("experiments/plots/exp1_objective_values.pdf", dpi=500,
       width=8, height=3.2, units="in")




# HOW MANY TIMES AT LEAST ONE VALID?
moc_n_valid <- sapply(experiment1_results$moc_results,
       function(res){
         sum(res$objective_values$validity == 0)
       })

tb_n_valid <- sapply(experiment1_results$tb_results,
       function(res){
         crowded_comparison_order <- get_crowded_comparison_order(res$objective_values)
         obj <- res$objective_values[select_population_indices(crowded_comparison_order, 40),]
         sum(obj$validity == 0)
       })

kov_n_valid <- sapply(experiment1_results$kov_results,
       function(res){
         raw_obj <- calculate_objective_values(explainer$data, explainer$times, experiment1_results$weights,
                                    res$original_observation, res$counterfactual_examples,
                                    res$original_prediction, res$predictions,
                                    NULL, experiment1_results$target_envelope_sf,
                                    data_range, NULL,
                                    1, 5)
         crowded_comparison_order <- get_crowded_comparison_order(raw_obj)
         raw_obj <- raw_obj[select_population_indices(crowded_comparison_order, 40),]
         sum(raw_obj$validity == 0)
       })


kov_n_valid_original <- sapply(experiment1_results$kov_results,
                      function(res){
                        min(sum(res$objective$validity == 0), 40)
                      })


clusters <- get_clusters(experiment1_results$dendrogram, k=5)[experiment1_results$not_in_target_cluster_ids][experiment1_results$sample_for_evaluation]


valid_df <- data.frame(
  method = c(rep("SurvMOC", length(moc_n_valid)),
             rep("RSF-PT", length(tb_n_valid)),
             rep("Kovalev", length(kov_n_valid))),
  n_valid = c(moc_n_valid, tb_n_valid, kov_n_valid),
  obs_id = rep(1:length(moc_n_valid), 3),
  obs_cluster = rep(clusters, 3)
)

valid_df$any_valid <- as.numeric(valid_df$n_valid > 0)
valid_df

valid_df_percentage <- data.frame(
  method = c("SurvMOC", "RSF-PT", "Kovalev"),
  any_valid = c(mean(valid_df$any_valid[valid_df$method == "SurvMOC"]),
                mean(valid_df$any_valid[valid_df$method == "RSF-PT"]),
                mean(valid_df$any_valid[valid_df$method == "Kovalev"]))
)


# plot percentage barplot
p1 <- ggplot(valid_df_percentage, aes(x = method, y = any_valid)) +
  geom_bar(stat = "identity", fill="darkorchid", width=0.8) +
  xlab("Method") +
  ylab("Percentage of cases with at least one\n valid counterfactual example") +
  theme_minimal() +
  geom_text(aes(label = scales::percent(any_valid)), vjust = +2, size = 4, color="white") +
  scale_y_continuous(labels = scales::percent, expand = c(0, 0)) +
  theme(panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank())


p2 <- ggplot(valid_df, aes(x = obs_id, y = n_valid, color = factor(obs_cluster))) +
  geom_point() +
  xlab("Observation ID") +
  ylab("Number of valid\ncounterfactual examples in top-40") +
  theme_bw() +
  theme(legend.position = "bottom") +
  facet_wrap(~method, ncol = 1) +
  scale_color_brewer(palette = "Set2", type = "qual", name = "Cluster") +
  theme(legend.position = "right")

library(ggpubr)
ggarrange(p1, p2, ncol = 2, nrow = 1, widths = c(1, 2))
ggsave("experiments/plots/exp1_validity.pdf", dpi=500,
       width=10, height=4, units="in")





# HOW MANY COUNTERFACTUALS ARE NOT-DOMINATED?
moc_n_nd <- sapply(experiment1_results$moc_results,
       function(res){
         obj <- res$objective_values
         sum(ecr::doNondominatedSorting(as.matrix(t(obj)))$ranks == 1)
       })

tb_n_nd <- sapply(experiment1_results$tb_results,
      function(res){
        obj <- res$objective_values
        crowded_comparison_order <- get_crowded_comparison_order(obj)
        obj <- obj[select_population_indices(crowded_comparison_order, 40),]
        sum(ecr::doNondominatedSorting(as.matrix(t(obj)))$ranks == 1)
      })

kov_n_nd <- sapply(experiment1_results$kov_results,
       function(res){
         raw_obj <- calculate_objective_values(explainer$data, explainer$times, experiment1_results$weights,
                                    res$original_observation, res$counterfactual_examples,
                                    res$original_prediction, res$predictions,
                                    NULL, experiment1_results$target_envelope_sf,
                                    data_range, NULL,
                                    1, 5)
         crowded_comparison_order <- get_crowded_comparison_order(raw_obj)
         raw_obj <- raw_obj[select_population_indices(crowded_comparison_order, 40),]
         sum(ecr::doNondominatedSorting(as.matrix(t(raw_obj)))$ranks == 1)
       })




nondominated_df <- data.frame(
  method = c(rep("SurvMOC", length(moc_n_valid)),
             rep("RSF-PT", length(tb_n_valid)),
             rep("Kovalev", length(kov_n_valid))),
  n_nondominated = c(moc_n_nd, tb_n_nd, kov_n_nd),
  obs_id = rep(1:length(moc_n_valid), 3),
  obs_cluster = rep(clusters, 3)
)


p1 <- ggplot(nondominated_df, aes(x = method, y = n_nondominated)) +
  geom_boxplot(staplewidth = 0.5,
                             outlier.shape = 1,
                             outlier.size = 0.5,
                             fill="darkorchid") +
  stat_summary(fun = mean, geom = "point", shape = 23, size = 2, fill = "violet") +
  xlab("Method") +
  ylab("Number of non-dominatedcounterfactual examples") +
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

ggarrange(p1, p2, ncol = 2, nrow = 1, widths = c(1, 2))
ggsave("experiments/plots/exp1_nondomination.pdf", dpi=500,
       width=10, height=4, units="in")


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
moc_obj <- experiment1_results$moc_results[[3]]$objective_values
tb_obj <- experiment1_results$tb_results[[3]]$objective_values
crowded_comparison_order <- get_crowded_comparison_order(tb_obj)
tb_obj <- tb_obj[select_population_indices(crowded_comparison_order, 40),]

moc_tb_obj <- rbind(moc_obj, tb_obj)

mean(ecr::doNondominatedSorting(as.matrix(t(moc_tb_obj)))$ranks[41:80] != 1)
mean(ecr::doNondominatedSorting(as.matrix(t(moc_tb_obj)))$ranks[1:40] != 1)


kov_obj <- calculate_objective_values(explainer$data, explainer$times, experiment1_results$weights,
                                      experiment1_results$kov_results[[3]]$original_observation,
                                      experiment1_results$kov_results[[3]]$counterfactual_examples,
                                      experiment1_results$kov_results[[3]]$original_prediction,
                                      experiment1_results$kov_results[[3]]$predictions,
                                      NULL, experiment1_results$target_envelope_sf,
                                      data_range, NULL,
                                      1, 5)

crowded_comparison_order <- get_crowded_comparison_order(kov_obj)
kov_obj <- kov_obj[select_population_indices(crowded_comparison_order, 40),]

moc_kov_obj <- rbind(moc_obj, kov_obj)
mean(ecr::doNondominatedSorting(as.matrix(t(moc_kov_obj)))$ranks[41:80] != 1)
mean(ecr::doNondominatedSorting(as.matrix(t(moc_kov_obj)))$ranks[1:40] != 1)


plot(explainer$data$x3, explainer$data$x4,
     ylim = c(15.5, 23.5), xlim = c(8.5, 13.5),
     col = scales::alpha("black", 0.2), pch=16)


points(experiment1_results$moc_results[[3]]$counterfactual_examples$x3,
      experiment1_results$moc_results[[3]]$counterfactual_examples$x4,
      pch=16)

points(experiment1_results$tb_results[[3]]$counterfactual_examples$x3,
       experiment1_results$tb_results[[3]]$counterfactual_examples$x4,
       col="darkgreen", pch=16)

points(experiment1_results$kov_results[[3]]$counterfactual_examples[,3],
       experiment1_results$kov_results[[3]]$counterfactual_examples[,4],
       col="blue", pch=16)

points(experiment1_results$kov_results[[3]]$original_observation[3],
       experiment1_results$kov_results[[3]]$original_observation[4],
       col="red", pch=4)




evaluation_sample <- explainer$data[experiment1_results$not_in_target_cluster_ids,][experiment1_results$sample_for_evaluation,]
which_has_1 <- which(evaluation_sample$x1 == 1)
which_has_0 <- which(evaluation_sample$x1 == 0)
mask_ones <- evaluation_sample$x1 == 1

moc_x1 <- lapply(experiment1_results$moc_results,
       function(res){
         res$counterfactual_examples$x1
       })


tb_x1 <- lapply(experiment1_results$tb_results,
                     function(res){
                       crowded_comparison_order <- get_crowded_comparison_order(res$objective_values)
                       sel_ces <- select_population_indices(crowded_comparison_order, 40)
                       ces <- res$counterfactual_examples[sel_ces,]
                       ces$x1
                     })

kov_x1 <- lapply(experiment1_results$kov_results,
                      function(res){
                        raw_obj <- calculate_objective_values(explainer$data, explainer$times, experiment1_results$weights,
                                                              res$original_observation, res$counterfactual_examples,
                                                              res$original_prediction, res$predictions,
                                                              NULL, experiment1_results$target_envelope_sf,
                                                              data_range, NULL,
                                                              1, 5)
                        crowded_comparison_order <- get_crowded_comparison_order(raw_obj)
                        ces <- res$counterfactual_examples[select_population_indices(crowded_comparison_order, 40),]
                        ces[,1]
                      })



moc_ones <- sapply(moc_x1, function(x) sum(x == 1))
tb_ones <- sapply(tb_x1, function(x) sum(x == 1))
kov_ones_real <- sapply(kov_x1, function(x) sum(x == 1))
kov_ones <- sapply(kov_x1, function(x) sum(abs(x-1)<abs(x)))

table(moc_ones == 40, mask_ones)
table(tb_ones == 40, mask_ones)
table(kov_ones == 40, mask_ones)


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

