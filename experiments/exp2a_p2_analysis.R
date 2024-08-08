#saveRDS(experiment3_results, "experiments/results/experiment3_results.rds")
library(ggplot2)
library(ggpubr)
library(dplyr)

experiment3_results <- readRDS("experiments/results/experiment3_results.rds")
experiment3_results$tb_results[[1]]$target_envelope


n_counterfactuals <- 10
n_parameter_sets <- nrow(experiment3_results$param_combinations)
n_observations <- 20

objective_values <- do.call("rbind", lapply(experiment3_results$tb_results,
                                            function(res){
                                              crowded_comparison_order <- get_crowded_comparison_order(res$objective_values)
                                              res$objective_values[select_population_indices(crowded_comparison_order, n_counterfactuals),]
                                            }))
nrow(objective_values) / length(experiment3_results$tb_results) == 10
nrow(objective_values) / (n_parameter_sets * n_observations) == 10

param_combinations <- experiment3_results$param_combinations[rep(1:n_parameter_sets, each=n_counterfactuals * n_observations), ]
rownames(param_combinations) <- NULL
nrow(param_combinations) == nrow(objective_values)

obs_ids <- rep(rep(1:20, each=n_counterfactuals), n_parameter_sets)
length(obs_ids) == nrow(objective_values)

res_df <- cbind(param_combinations, objective_values, obs_ids)
res_df


res_df_long <- reshape2::melt(res_df, id.vars = c("paths_per_tree", "paths_per_counterfactual", "step", "obs_ids"))

res_df_long$step <- factor(res_df_long$step,
                                      levels = c("1", "2"),
                                      labels = c("step = 1",
                                                 "step = 2"),
                                      ordered = TRUE)

res_df_long$variable <- paste(res_df_long$variable, "loss")
res_df_long$variable <- factor(res_df_long$variable,
                               levels = c("validity loss",
                                          "similarity loss",
                                          "sparsity loss",
                                          "plausibility loss"),
                               ordered = TRUE)


ggplot(res_df_long, aes(x = factor(paths_per_tree), y = value,
                        color = factor(paths_per_counterfactual))) +
  geom_boxplot(staplewidth = 0.5,
               outlier.shape = 1,
               outlier.size = 0.5,
               fill = "ghostwhite",
               outlier.alpha = 0.3
  ) +
  # add mean
  stat_summary(fun = mean, geom = "point",
               shape = 23, size = 1.5, fill = "white",
               position = position_dodge(width = 0.75)) +
  facet_grid(variable ~ factor(step), scales = "free_y") +
  scale_color_brewer(type = "qual", palette = "Set1",
                     name = "Number of paths\nper counterfactual example") +
  labs(x = "Number of paths per tree", y = "Objective value") +
  theme_bw() +
  theme(legend.position = "bottom")

ggsave("experiments/plots/exp2a_bis_objective_values.pdf", dpi=500, width=10, height=8)


# VALIDITY
n_valid <- data.frame(res_df %>%
                        group_by(paths_per_tree, paths_per_counterfactual, step, obs_ids) %>%
                        summarise(any_valid = any(validity == 0)) %>%
                        group_by(paths_per_tree, paths_per_counterfactual, step) %>%
                        summarise(mean_valid = mean(any_valid)))

n_valid[1:24,]



ggplot(n_valid, aes(x = factor(paths_per_tree), y = mean_valid, fill = factor(paths_per_counterfactual))) +
  geom_bar(stat="identity", position = "dodge", width=0.8) +
  facet_grid(~step) +
  scale_fill_brewer(type = "qual", palette = "Set1", name = "Number of paths\nper counterfactual example") +
  scale_y_continuous(labels = scales::percent, expand=c(0, 0), limits = c(0, 1.03)) +
  labs(x = "Number of paths per tree",
       y = "Percentage of cases with at least one\n valid counterfactual example") +
  theme_bw() +
  theme(legend.position = "bottom")
ggsave("experiments/plots/exp3_validity.pdf", dpi=500, width=10, height=4)




# TIME
param_combinations <- experiment3_results$param_combinations[rep(1:n_parameter_sets, each=n_observations), ]
nrow(param_combinations) == length(n_iterations)
nrow(param_combinations) == length(experiment3_results$tb_time)

time_df <- cbind(param_combinations,time=experiment3_results$tb_time)
time_df_long <- reshape2::melt(time_df,
                               id.vars = c("paths_per_tree", "paths_per_counterfactual", "step"))



p1 <- ggplot(time_df, aes(x = factor(paths_per_tree), y = time, color = factor(paths_per_counterfactual))) +
  geom_boxplot(staplewidth = 0.5,
               outlier.shape = 1,
               outlier.size = 0.5,
               fill = "ghostwhite",
               outlier.alpha = 0.3
  ) +
  # add mean
  stat_summary(fun = mean, geom = "point",
               shape = 23, size = 1.5, fill = "white",
               position = position_dodge(width = 0.75)) +
  scale_color_brewer(type = "qual", palette = "Set1",
                     name = "Number of paths\nper counterfactual example") +
  labs(x = "Number of paths per tree",
       y = "Execution time [s]") +
  theme_bw() +
  theme(legend.position = "bottom")
p1

ggsave("experiments/plots/exp2a_bis_time.pdf", dpi=500, width=5, height=3.5)



