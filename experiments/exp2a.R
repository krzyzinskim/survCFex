library(ranger)
library(survival)
library(survex)
devtools::load_all(".")

df <- read.csv("experiments/data/lung_dataset.csv")
df <- df[complete.cases(df), ]
nrow(df)
set.seed(123)
model <- ranger(Surv(time, status) ~ ., data = df,
                num.trees = 200,
                max.depth = 7,
                min.node.size = 5)
explainer <- explain(model,
                     data = df[3:9],
                     y = Surv(df$time, df$status))

preds <- predict(explainer, df[,3:9])

plot_survival_weights(explainer, explainer$times, p=0, q=0)
weights <- survival_weights(explainer, explainer$times, p=0, q=0)


clustering <- get_hierarchical_clustering(explainer, weights)
analyze_clustering(clustering)
plot(dendrogram)


plot_prediction_bands(explainer, get_clusters(clustering, k=4),
               alpha = 0.5,
               lambda = 0.1) +
  theme(legend.position = "bottom")
ggsave("experiments/plots/exp2_envelopes.pdf", width = 5, height = 3, dpi = 500)



target_envelope_sf <- get_envelope(preds, get_clusters(dendrogram, k=4),
                                   cluster_id=1, q = 0.1)
target_envelope_chf <- translate_target_envelope(target_envelope_sf)


not_in_target_cluster_ids <- which(get_clusters(dendrogram, k=4) != 1)
tmp_df <- df[not_in_target_cluster_ids,]


n_samples <- 20
set.seed(123)
sample_for_evaluation <- sample(1:nrow(tmp_df), n_samples)
org_distances <- distance_from_target_envelope(predict(explainer, tmp_df[, -c(1:2)]),
                              target_envelope_sf,
                              explainer$times,
                              weights)

boxplot(org_distances)
min(org_distances)


pop_sizes <- c(10, 20, 40, 100, 200)
valid_threshold <- c("median", "none", "0")
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
  if (valid_threshold == "none"){
    valid_threshold <- NULL
  } else if (valid_threshold == "0"){
    valid_threshold <- 0
  }
  init_method <- param_combinations$init_method[i]
  cat(paste("pop_size:", pop_size, "valid_threshold:", valid_threshold, "init_method:", init_method))

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
  saveRDS(moc_results, file = paste0("experiments/results/exp2_moc_results/", i, ".rds"))
}

