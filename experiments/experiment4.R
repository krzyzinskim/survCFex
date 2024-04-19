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

pop_sizes <- c(10, 20, 40, 100, 200)
valid_threshold <- c("median", "none", "quantile")
initialization_method <- c("ice", "background")

param_combinations <- expand.grid(pop_sizes, valid_threshold, initialization_method)
colnames(param_combinations) <- c("pop_size", "valid_threshold", "init_method")
number_of_combinations <- nrow(param_combinations)

moc_results <- list()
moc_time <- numeric(number_of_combinations * n_samples)

for (i in 1:number_of_combinations) {
  cat(paste("\nCombination", i, "of", number_of_combinations, "\n"))
  pop_size <- param_combinations$pop_size[i]
  valid_threshold <- param_combinations$valid_threshold[i]
  valid_threshold_name <- valid_threshold
  if (valid_threshold_name == "none"){
    valid_threshold <- NULL
  } else if (valid_threshold_name == "quantile"){
    valid_threshold <- function(x) quantile(x, 0.2)
    valid_threshold_name <- "quantile_0.2"
  }
  init_method <- param_combinations$init_method[i]
  cat(paste("pop_size:", pop_size, "valid_threshold:", valid_threshold_name, "init_method:", init_method))

  for (j in 1:n_samples){
    cat(paste("\nSample", j, "/", n_samples, "\n"))
    obs <- tmp_df[sample_for_evaluation[j], -c(1,2)]
    moc_start <- Sys.time()
    moc_results[[(i-1) * n_samples + j]] <- multiobjective_counterfactuals(explainer, obs,
                                                                           times = explainer$times,
                                                                           target_envelope = target_envelope_sf,
                                                                           seed = 123,
                                                                           population_size = pop_size,
                                                                           max_generations = 100,
                                                                           weights = weights,
                                                                           validity_threshold = valid_threshold,
                                                                           initialization_type = init_method,
                                                                           verbose = FALSE,
                                                                           tol_iter = 5)
    moc_time[(i-1) * n_samples + j] <- as.numeric(Sys.time() - moc_start, units = "secs")
  }
  cat(paste("Total time for this combination:", round(sum(moc_time[(i-1) * n_samples + 1:(n_samples)]), 2), "seconds", "\n"))
  cat("\n")
  saveRDS(moc_results, file = paste0("experiments/results/exp4_moc_results/", i, ".rds"))
}

