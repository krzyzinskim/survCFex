library(survival)
library(ranger)
library(survex)

df <- survival::veteran
df

model <- ranger(Surv(time, status) ~ ., data = df, num.trees = 200, max.depth = 5)



explainer <- explain(model,
                     data = df[-c(3,4)],
                     y = Surv(df$time, df$status))


x <- df[23, -c(3,4)]
preds <- predict(explainer, explainer$data)
original_pred <- predict(explainer, x)

plot_survival_weights(explainer, explainer$times, p=0, q=0)
weights <- survival_weights(explainer, explainer$times, p=0, q=0)



plot_predictions(preds, explainer$times)



dists <- survival_distance_matrix(preds, explainer$times, weights)
dendrogram <- get_clustering_dendrogram(dists)
plot(dendrogram)


plot_envelopes(preds, get_clusters(dendrogram, k=4), alpha = 0.5,
                       explainer$times, q = 0.05, original_pred)


target_envelope <- get_envelope(preds, get_clusters(dendrogram, k=4),
                                        cluster_id=2, q = 0.05)



moc_res <- multiobjective_counterfactuals(explainer, x,
                                          times = explainer$times,
                                          target_envelope = target_envelope,
                                          seed=23,
                                          population_size=40,
                                          max_generations=50,
                                          weights = weights,
                                          validity_threshold = 0.1,
                                          verbose=FALSE,
                                          tol_iter=5)
analyze(moc_res)



new_target_envelope <- translate_target_envelope(target_envelope)


tb_res <- treebased_counterfactuals(explainer, x,
                                    times = explainer$times,
                                    target_envelope = new_target_envelope,
                                    k_paths = 5,
                                    verbose=TRUE)
tb_res
analyze(tb_res)



