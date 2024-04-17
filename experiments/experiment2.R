library(ranger)
library(survival)
library(survex)
devtools::load_all(".")

df <- read.csv("experiments/data/lung_dataset.csv")
df <- df[complete.cases(df), ]

set.seed(123)
model <- ranger(Surv(time, status) ~ ., data = df, num.trees = 200, max.depth = 5)
explainer <- explain(model,
                     data = df[3:9],
                     y = Surv(df$time, df$status))

preds <- predict(explainer, df[,3:9])

plot_survival_weights(explainer, explainer$times, p=0, q=0)
weights <- survival_weights(explainer, explainer$times, p=0, q=0)


dists <- survival_distance_matrix(preds, explainer$times, weights)
dendrogram <- get_clustering_dendrogram(dists)
plot(dendrogram)

plot(get_clustering_utilities(dendrogram, max_k = 6))

plot_envelopes(preds, get_clusters(dendrogram, k=3),
               alpha = 0.5,
               explainer$times,
               q = 0.1)


target_envelope_sf <- get_envelope(preds, get_clusters(dendrogram, k=3),
                                   cluster_id=1, q = 0.05)
target_envelope_chf <- translate_target_envelope(target_envelope_sf)


not_in_target_cluster_ids <- which(get_clusters(dendrogram, k=3) != 1)
tmp_df <- df[not_in_target_cluster_ids,]

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
                                                               target_envelope = target_envelope_chf,
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



