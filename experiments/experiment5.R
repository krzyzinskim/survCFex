library(ranger)
library(survival)
library(survex)
devtools::load_all(".")

set.seed(123)
df <- survival::pbc
df <- df[complete.cases(df),-1]
nrow(df)
colnames(df)

model <- ranger(Surv(time, status==2) ~ ., data = df, num.trees = 200, max.depth = 7,
                min.node.size = 5)
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

plot_envelopes(preds, get_clusters(dendrogram, k=4),
               alpha = 0.5,
               explainer$times,
               q = 0.05)

target_envelope_sf <- get_envelope(preds, get_clusters(dendrogram, k=4),
                                   cluster_id=2, q = 0.05)
target_envelope_chf <- translate_target_envelope(target_envelope_sf)


not_in_target_cluster_ids <- which(get_clusters(dendrogram, k=4) != 2)
tmp_df <- df[not_in_target_cluster_ids,]

obs <- tmp_df[1:50,-c(1,2)]
plot_envelopes(preds, get_clusters(dendrogram, k=4),
               alpha = 0.5,
               explainer$times,
               q = 0.05,
               predict(explainer, obs))


n_samples <- 20
set.seed(123)
sample_for_evaluation <- sample(1:nrow(tmp_df), n_samples)

paths_per_tree <- c(1, 4, 8)
paths_per_counterfactual <- c(1, 3, 5, 7)
step <- c(1, 2)

param_combinations <- expand.grid(paths_per_tree, paths_per_counterfactual, step)
colnames(param_combinations) <- c("paths_per_tree", "paths_per_counterfactual", "step")
number_of_combinations <- nrow(param_combinations)

tb_results <- list()
tb_time <- numeric(number_of_combinations * n_samples)

for (i in 1:number_of_combinations) {
  cat(paste("\nCombination", i, "of", number_of_combinations, "\n"))
  paths_per_tree <- param_combinations$paths_per_tree[i]
  paths_per_counterfactual <- param_combinations$paths_per_counterfactual[i]
  step <- param_combinations$step[i]
  max_tries <- paths_per_tree * 200

  cat(paste("paths_per_tree:", paths_per_tree, "paths_per_counterfactual:", paths_per_counterfactual, "step:", step, "\n"))

  for (j in 1:n_samples){
    cat(paste("\nSample", j, "/", n_samples, "\n"))
    obs <- tmp_df[sample_for_evaluation[j], -c(1,2)]
    tb_start <- Sys.time()
    tb_results[[(i-1) * n_samples + j]] <- treebased_counterfactuals(explainer, obs,
                                                                     times = explainer$times,
                                                                     target_envelope = target_envelope_chf,
                                                                     paths_per_counterfactual = paths_per_counterfactual,
                                                                     paths_per_tree = paths_per_tree,
                                                                     step = step,
                                                                     max_tries = max_tries,
                                                                     verbose = FALSE)
    tb_time[(i-1) * n_samples + j] <- as.numeric(Sys.time() - tb_start, units = "secs")
  }
  cat(paste("Total time for this combination:", round(sum(tb_time[(i-1) * n_samples + 1:(n_samples)]), 2), "seconds", "\n"))
  cat("\n")
  saveRDS(tb_results, paste0("experiments/results/exp5_tb_results/", i, ".rds"))
}

