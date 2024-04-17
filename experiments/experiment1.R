library(survival)
library(ranger)
library(survex)

set.seed(42)
df <- read.csv("experiments/data/exp1_data_complex.csv")

train_ids <- sample(1:nrow(df), 0.8*nrow(df))
df_train <- df[train_ids,]
df_test <- df[-train_ids,]

evaluation_sample_ids <- sample(1:nrow(df), 0.2*nrow(df))

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

performance <- model_performance(explainer)
plot(performance)
plot(performance, metrics_type = "scalar")

x <- df[23, 1:5]
x

preds <- predict(explainer, explainer$data)
original_pred <- predict(explainer, x)

plot_predictions(rbind(preds, original_pred), explainer$times)

plot_survival_weights(explainer, explainer$times, p=0, q=0)
weights <- survival_weights(explainer, explainer$times, p=0, q=0)
dists <- survival_distance_matrix(preds, explainer$times, weights)
dendrogram <- get_clustering_dendrogram(dists)
plot(dendrogram)

plot(get_clustering_utilities(dendrogram, max_k = 6))


plot_envelopes(preds, get_clusters(dendrogram, k=5), alpha = 0.5,
               explainer$times, q = 0.1, original_pred)


target_envelope <- get_envelope(preds, get_clusters(dendrogram, k=5),
                                cluster_id=5, q = 0.05)


moc_res <- multiobjective_counterfactuals(explainer, x,
                                          times = explainer$times,
                                          target_envelope = target_envelope,
                                          seed=23,
                                          population_size=30,
                                          max_generations=100,
                                          weights = weights,
                                          validity_threshold = 0,
                                          verbose=TRUE,
                                          tol_iter=4)

analyze(moc_res)





