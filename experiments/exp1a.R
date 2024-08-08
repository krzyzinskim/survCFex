library(survival)
library(ranger)
library(survex)
devtools::load_all(".")

set.seed(42)
df <- read.csv("experiments/data/exp1_data_complex.csv")

set.seed(42)
model <- ranger(Surv(time, event) ~ .,
                data = df,
                num.trees = 100,
                max.depth = 10,
                mtry = 3,
                min.node.size = 5)


explainer <- explain(model,
                     data = df[1:5],
                     y = Surv(df$time, df$event))

class(explainer)

preds <- predict(explainer, explainer$data)
weights <- survival_weights(explainer, p=0, q=0)

times <- explainer$times
n_times <- length(times)

weights_df <- data.frame(
  weight = c(
    survival_weights(explainer, times, p=0, q=0),
    survival_weights(explainer, times, p=1, q=0),
    survival_weights(explainer, times, p=0, q=1),
    survival_weights(explainer, times, p=0.5, q=0.5),
    survival_weights(explainer, times, p=1, q=1)
  ),
  time = rep(times, 5),
  params = rep(c("p=0, q=0", "p=1, q=0", "p=0, q=1", "p=0.5, q=0.5", "p=1, q=1"), each = n_times)
)


p1 <- ggplot(weights_df, aes(x = time, y = weight, group = params, color = params)) +
  geom_step(linewidth=0.6) +
  xlab("Time") +
  ylab("Weight") +
  ggtitle("Survival weights") +
  theme_minimal() +
  scale_color_brewer(type = "qual", palette = "Set1") +
  guides(color = guide_legend(title = "Parameters", nrow=2)) +
  theme(legend.position = "bottom")
p1

km <- survival::survfit(explainer$y ~ 1)
surv_estimator <- stepfun(km$time, c(1, km$surv))
sf <- surv_estimator(times)
sf

sf_df <- data.frame(time = c(0, times), survival = c(1, sf))

p2 <- ggplot(sf_df, aes(x = time, y = survival)) +
  geom_step(linewidth=0.6) +
  xlab("Time") +
  ylab("Survival probability") +
  theme_minimal() +
  ggtitle("Kaplan-Meier estimate")

ggarrange(p2, p1, ncol = 2, widths = c(1.5, 2),
          common.legend = TRUE, legend = "bottom")

ggsave("experiments/plots/exp1_weights.pdf", width = 7, height = 3, dpi = 500)


plot_survival_weights(explainer, explainer$times, p=0, q=10)


dendrogram <- get_hierarchical_clustering(explainer, weights)
analyze_clustering(explainer, dendrogram)

plot_prediction_bands(explainer,
               get_clusters(dendrogram, k=5),
               alpha = 0.5,
               lambda = 0.1) +
  theme(legend.position = "bottom")
ggsave("experiments/plots/exp1_envelopes.pdf", width = 5, height = 3, dpi = 500)

target_prediction_band_sf <- get_prediction_band(explainer, get_clusters(dendrogram, k=5),
                                cluster_id=5, lambda = 0.05)
target_prediction_band_chf <- translate_target_prediction_band(target_prediction_band_sf)


mu_target <- mean_time_to_survival(
  explainer,
  predictions = matrix((target_envelope_sf$lower_bound + target_envelope_sf$upper_bound) / 2,
                       nrow = 1))


not_in_target_cluster_ids <- which(get_clusters(dendrogram, k=5) != 5)
tmp_df <- df[not_in_target_cluster_ids,]

k <- 50
set.seed(42)
sample_for_evaluation <- sample(1:nrow(tmp_df), k)

moc_results <- list()
tb_results <- list()
kov_results <- list()

moc_time <- numeric(k)
tb_time <- numeric(k)
kov_time <- numeric(k)

for (i in seq_len(k)){
  print(i)
  obs <- tmp_df[sample_for_evaluation[i], 1:5]

  moc_start <- Sys.time()
  moc_results[[i]] <- surv_moc(explainer, obs,
                               target_prediction_band = target_prediction_band_sf,
                               seed = i,
                               population_size = 40,
                               max_generations = 100,
                               weights = weights,
                               validity_threshold = 0,
                               verbose = FALSE,
                               tol_iter = 4)
  moc_time[i] <- as.numeric(Sys.time() - moc_start, units = "secs")

  tb_start <- Sys.time()
  tb_results[[i]] <- st_pt(explainer, obs,
                           target_prediction_band = target_prediction_band_chf,
                           max_counterfactuals = 40,
                           paths_per_tree = 20,
                           verbose = FALSE)
  tb_time[i] <- as.numeric(Sys.time() - tb_start, units = "secs")

  mu_original <- mean_time_to_survival(explainer, predictions = predict(explainer, obs))
  kov_start <- Sys.time()
  kov_results[[i]] <- kovalev_method(explainer, obs, seed=i,
                                     r = mu_target - mu_original,
                                     num_iter = 100,
                                     num_particles = 200,
                                     verbose = FALSE)
  kov_time[i] <- as.numeric(Sys.time() - kov_start, units = "secs")
}



