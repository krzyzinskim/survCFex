library(survival)
library(ranger)
library(survex)
devtools::load_all(".")

set.seed(42)
df <- read.csv("experiments/data/exp1_data_complex.csv")

# train_ids <- sample(1:nrow(df), 0.8*nrow(df))
# df_train <- df[train_ids,]
# df_test <- df[-train_ids,]
#
# evaluation_sample_ids <- sample(1:nrow(df), 0.2*nrow(df))

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

preds <- predict(explainer, explainer$data)
weights <- survival_weights(explainer, explainer$times, p=0, q=0)

dists <- survival_distance_matrix(preds, explainer$times, weights)
dendrogram <- get_clustering_dendrogram(dists)
plot(dendrogram)

plot(get_clustering_utilities(dendrogram, max_k = 6))

plot_envelopes(preds, get_clusters(dendrogram, k=5),
               alpha = 0.5,
               explainer$times,
               q = 0.05)

target_envelope_sf <- get_envelope(preds, get_clusters(dendrogram, k=5),
                                cluster_id=5, q = 0.05)
target_envelope_chf <- translate_target_envelope(target_envelope_sf)
mu_target <- mean_time_to_survival(
  explainer,
  predictions = matrix(target_envelope_sf$lower_bound, nrow = 1))


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
  moc_results[[i]] <- multiobjective_counterfactuals(explainer, obs,
                                                     times = explainer$times,
                                                     target_envelope = target_envelope_sf,
                                                     seed = i,
                                                     population_size = 40,
                                                     max_generations = 100,
                                                     weights = weights,
                                                     validity_threshold = 0,
                                                     verbose = FALSE,
                                                     tol_iter = 4)
  moc_time[i] <- as.numeric(Sys.time() - moc_start, units = "secs")

  tb_start <- Sys.time()
  tb_results[[i]] <- treebased_counterfactuals(explainer, obs,
                                               times = explainer$times,
                                               target_envelope = target_envelope_chf,
                                               max_counterfactuals = 40,
                                               k_paths = 20,
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



