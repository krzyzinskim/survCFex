saveRDS(experiment4_results, "experiments/results/experiment4_results.rds")

library(ggplot2)
library(ggpubr)
library(dplyr)
experiment4_results <- readRDS("experiments/results/experiment4_results.rds")
experiment4_results$moc_results[[1]]$target_envelope
analyze(experiment4_results$moc_results[[1]])

n_counterfactuals <- 10
n_parameter_sets <- nrow(experiment4_results$param_combinations)
n_observations <- 20

objective_values <- do.call("rbind", lapply(experiment4_results$moc_results,
                                            function(res){
                                              crowded_comparison_order <- get_crowded_comparison_order(res$objective_values)
                                              res$objective_values[select_population_indices(crowded_comparison_order, n_counterfactuals),]
                                            }))
nrow(objective_values) / length(experiment4_results$moc_results) == 10
nrow(objective_values) / (n_parameter_sets * n_observations) == 10

param_combinations <- experiment4_results$param_combinations[rep(1:30, each=n_counterfactuals * n_observations), ]
rownames(param_combinations) <- NULL
nrow(param_combinations) == nrow(objective_values)

obs_ids <- rep(rep(1:20, each=n_counterfactuals), n_parameter_sets)
length(obs_ids) == nrow(objective_values)

res_df <- cbind(param_combinations, objective_values, obs_ids)
res_df


res_df_long <- reshape2::melt(res_df, id.vars = c("pop_size", "valid_threshold", "init_method", "obs_ids"))


res_df_long$valid_threshold <- factor(res_df_long$valid_threshold,
                                      levels = c("none", "median", "quantile"),
                                      labels = c("validity threshold = none",
                                                 "validity threshold = median",
                                                 "validity threshold = q(0.2)"),
                                      ordered = TRUE)

res_df_long$init_method <- factor(res_df_long$init_method,
                                  levels = c("ice", "background"),
                                  labels = c("ICE", "background"),
                                  ordered = TRUE)



ggplot(res_df_long, aes(x = factor(pop_size), y = value, color = init_method)) +
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
  facet_grid(variable ~ valid_threshold, scales = "free_y") +
  scale_color_brewer(type = "qual", palette = "Set1", name = "Initialization method") +
  labs(x = "Population size", y = "Objective value") +
  theme_bw() +
  theme(legend.position = "bottom")

ggsave("experiments/plots/exp4_objective_values.pdf", dpi=500, width=10, height=8)


# VALIDITY
n_valid <- data.frame(res_df %>%
                        group_by(pop_size, valid_threshold, init_method, obs_ids) %>%
                        summarise(any_valid = any(validity == 0)) %>%
                        group_by(pop_size, valid_threshold, init_method) %>%
                        summarise(mean_valid = mean(any_valid)))

n_valid[1:30,]

n_valid$valid_threshold <- factor(n_valid$valid_threshold,
                                  levels = c("none", "median", "quantile"),
                                  labels = c("validity threshold = none",
                                             "validity threshold = median",
                                             "validity threshold = q(0.2)"),
                                  ordered = TRUE)

n_valid$init_method <- factor(n_valid$init_method,
                              levels = c("ice", "background"),
                              labels = c("ICE", "background"),
                              ordered = TRUE)



ggplot(n_valid, aes(x = factor(pop_size), y = mean_valid, fill = init_method)) +
  geom_bar(stat="identity", position = "dodge", width=0.8) +
  facet_grid(~valid_threshold) +
  scale_fill_brewer(type = "qual", palette = "Set1", name = "Initialization method") +
  scale_y_continuous(labels = scales::percent, expand=c(0, 0), limits = c(0, 1.03)) +
  labs(x = "Population size",
       y = "Percentage of cases with at least one\n valid counterfactual example") +
  theme_bw() +
  theme(legend.position = "bottom")
ggsave("experiments/plots/exp4_validity.pdf", dpi=500, width=10, height=4)




# TIME
n_iterations <- sapply(experiment4_results$moc_results, function(x) max(x$history$generation))
param_combinations <- experiment4_results$param_combinations[rep(1:30, each=n_observations), ]
nrow(param_combinations) == length(n_iterations)
nrow(param_combinations) == length(experiment4_results$moc_time)

time_df <- cbind(param_combinations, n_iterations = n_iterations, time=experiment4_results$moc_time)
time_df$init_method <- factor(time_df$init_method,
                              levels = c("ice", "background"),
                              labels = c("ICE", "background"),
                              ordered = TRUE)

time_df$valid_threshold <- factor(time_df$valid_threshold,
                                  levels = c("none", "median", "quantile"),
                                  labels = c("validity threshold = none",
                                             "validity threshold = median",
                                             "validity threshold = q(0.2)"),
                                  ordered = TRUE)
time_df_long <- reshape2::melt(time_df, id.vars = c("pop_size", "valid_threshold", "init_method"))



p1 <- ggplot(time_df, aes(x = factor(pop_size), y = time, color = init_method)) +
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
  scale_color_brewer(type = "qual", palette = "Set1", name = "Initialization method") +
  labs(x = "Population size", y = "Execution time [s]") +
  facet_grid( ~ valid_threshold, scales = "free_y") +
  theme_bw() +
  theme(legend.position = "bottom")
p1

p2 <- ggplot(time_df, aes(x = factor(pop_size), y = n_iterations, color = init_method)) +
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
  scale_color_brewer(type = "qual", palette = "Set1", name = "Initialization method") +
  labs(x = "Population size", y = "Number of iterations") +
  facet_grid( ~ valid_threshold, scales = "free_y") +
  theme_bw() +
  theme(legend.position = "bottom")
p2

ggarrange(p1, p2, nrow=2, common.legend = T, legend = "bottom")
ggsave("experiments/plots/exp4_time_iterations.pdf", dpi=500, width=10, height=6)


time_df %>%
  mutate(time_per_iteration = time / n_iterations) %>%
  ggplot(aes(x = factor(pop_size), y = time_per_iteration)) +
  geom_boxplot(staplewidth = 0.5,
               outlier.shape = 1,
               outlier.size = 0.5,
               fill = "ghostwhite",
               outlier.alpha = 0.3
  ) + stat_summary(fun = mean, geom = "point",
                   shape = 23, size = 1.5, fill = "white",
                   position = position_dodge(width = 0.75)) +
  labs(x = "Population size", y = "Time per iteration [s]") +
  facet_grid( ~ valid_threshold, scales = "free_y") +
  theme_bw() +
  theme(legend.position = "bottom")
ggsave("experiments/plots/exp2_time_per_iteration.pdf", dpi=500, width=10, height=6)


time_df$time_per_iteration = time_df$time / time_df$n_iterations

plot(time_df$pop_size, time_df$time_per_iteration,
     xlab="Population size", ylab="Time per iteration [s]")

m1 <- lm(time_df$time_per_iteration ~ time_df$pop_size)
summary(m1)
coef(m1)


plot(time_df$pop_size, time_df$time_per_iteration,
     xlab="Population size", ylab="Time per iteration [s]")
abline(m1, col="red")
# linear complexity


# history analysis for validity threshold of 0
indices <- which(param_combinations$valid_threshold %in% c("0", "median"))
population_sizes <- param_combinations$pop_size[indices]
validity_thresholds <- param_combinations$valid_threshold[indices]
histories <- lapply(indices, function(i) experiment2_results$moc_results[[i]]$history)
histories_df <- do.call("rbind", lapply(1:length(histories),
                                        function(i){
                                          history <- histories[[i]]
                                          history$id <- i
                                          history$pop_size <- population_sizes[i]
                                          history$valid_threshold <- validity_thresholds[i]
                                          history
                                        }))

agg_histories_df <- histories_df %>%
  group_by(pop_size, valid_threshold, id, generation) %>%
  summarize(validity = mean(validity))


ggplot(agg_histories_df, aes(x = generation, y = validity, color = factor(pop_size))) +
  stat_summary(fun = mean, geom = "line") +
  facet_grid(~valid_threshold) +
  theme_bw()

