library(ranger)
library(survival)
library(survex)
devtools::load_all(".")

df <- survival::pbc
df <- df[complete.cases(df),-1]
nrow(df)
colnames(df)

set.seed(123)
model <- ranger(Surv(time, status==2) ~ ., data = df, num.trees = 200, max.depth = 7)
explainer <- explain(model,
                     data = df[-c(1,2)],
                     y = Surv(df$time, df$status == 2))

preds <- predict(explainer, explainer$data)
plot_predictions(preds, explainer$times)

plot_survival_weights(explainer, explainer$times, p=0, q=0)
weights <- survival_weights(explainer, explainer$times, p=0, q=0)


dists <- survival_distance_matrix(preds, explainer$times, weights)
dendrogram <- get_clustering_dendrogram(dists)
plot(dendrogram)

plot(get_clustering_utilities(dendrogram, max_k = 6))

plot_envelopes(preds, get_clusters(dendrogram, k=5),
               alpha = 0.5,
               explainer$times,
               q = 0.1)

target_envelope_sf <- get_envelope(preds, get_clusters(dendrogram, k=5),
                                   cluster_id=2, q = 0.05)
target_envelope_chf <- translate_target_envelope(target_envelope_sf)


not_in_target_cluster_ids <- which(get_clusters(dendrogram, k=5) != 2)
tmp_df <- df[not_in_target_cluster_ids,]

obs <- tmp_df[1,-c(1,2)]
plot_envelopes(preds, get_clusters(dendrogram, k=5),
               alpha = 0.5,
               explainer$times,
               q = 0.1,
               predict(explainer, obs))

moc_res <- multiobjective_counterfactuals(explainer, obs,
                               times = explainer$times,
                               target_envelope = target_envelope_sf,
                               seed = 123,
                               population_size = 10,
                               max_generations = 2,
                               weights = weights,
                               validity_threshold = 0.05,
                               verbose = TRUE,
                               tol_iter = 5)

categorical_variables_indices <- 3
data_range <-  apply(explainer$data[-categorical_variables_indices], 2, range)
gower_distance_loss(moc_res$original_observation, moc_res$counterfactual_examples,
                    data_range, categorical_variables_indices = categorical_variables_indices, p = 1)
moc_res$objective_values$similarity


analyze(moc_res)

pop_sizes <- c(5, 10, 20, 40, 100, 200)
path_numbers <- c(1, 5, 10, 20, 40, 100)
l <- 6

k <- 25
set.seed(123)
sample_for_evaluation <- sample(1:nrow(tmp_df), k)

moc_results <- list()
tb_results <- list()

moc_time <- numeric(k * l)
tb_time <- numeric(k * l)

for (i in 1:k) {
  print(i)
  obs <- tmp_df[sample_for_evaluation[i], -c(1,2)]
  for (j in 1:l){
    moc_start <- Sys.time()
    moc_results[[(i-1) * l + j]] <- multiobjective_counterfactuals(explainer, obs,
                                                                   times = explainer$times,
                                                                   target_envelope = target_envelope_sf,
                                                                   seed = 123,
                                                                   population_size = pop_sizes[j],
                                                                   max_generations = 100,
                                                                   weights = weights,
                                                                   validity_threshold = 0.05,
                                                                   verbose = FALSE,
                                                                   tol_iter = 5)
    moc_time[(i-1) * l + j] <- as.numeric(Sys.time() - moc_start, units = "secs")

    tb_start <- Sys.time()
    tb_results[[(i-1) * l + j]] <- treebased_counterfactuals(explainer, obs,
                                                             times = explainer$times,
                                                             target_envelope = target_envelope_chf,
                                                             max_counterfactuals = 40,
                                                             k_paths = path_numbers[j],
                                                             verbose = FALSE)
    tb_time[(i-1) * l + j] <- as.numeric(Sys.time() - tb_start, units = "secs")
  }
}



