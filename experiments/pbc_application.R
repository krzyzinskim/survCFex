library(ranger)
library(survival)
library(survex)
devtools::load_all(".")

set.seed(123)
df <- survival::pbc
nrow(df)
df <- df[complete.cases(df),-1]
nrow(df)
colnames(df)

model <- ranger(Surv(time, status==2) ~ ., data = df, num.trees = 200, max.depth = 7,
                min.node.size = 5)
explainer <- explain(model,
                     data = df[-c(1,2)],
                     y = Surv(df$time, df$status == 2))

plot_survival_weights(explainer, explainer$times, p=0, q=0)
weights <- survival_weights(explainer, explainer$times, p=0, q=0)


dendrogram <- get_hierarchical_clustering(explainer, weights)

analyze_clustering(dendrogram, max_k = 6)

plot_prediction_bands(explainer, get_clusters(dendrogram, k=4),
               alpha = 0.5,
               lambda = 0.05) +
  theme(legend.position = "bottom")

target_prediction_band_sf <- get_prediction_band(explainer, get_clusters(dendrogram, k=4),
                                   cluster_id=2, lambda = 0.05)

target_prediction_band_chf <- translate_target_prediction_band(target_prediction_band_sf)

not_in_target_cluster_ids <- which(get_clusters(dendrogram, k=4) != 2)


xstar <- df[not_in_target_cluster_ids[4], -c(1,2)]




plot_prediction_bands(explainer, get_clusters(dendrogram, k=4),
               alpha = 0.5,
               lambda = 0.05,
               preds_to_plot=predict(explainer, xstar))
ggsave("experiments/plots/app_pbc_pred.pdf", width = 5, height = 3, dpi = 500)


xstar
moc_res <- surv_moc(explainer, xstar, target_prediction_band_sf,
                    fixed_variables_indices = c(2, 3, 6, 17),
                    categorical_variables_indices = c(4, 5, 7),
                    verbose = TRUE,
                    seed = 123,
                    population_size = 50
                    )
analyze(moc_res)





stpt_res <- st_pt(explainer, xstar, target_prediction_band_chf,
                  fixed_variables_indices = c(2, 3, 6, 17),
                  paths_per_tree = 3,
                  paths_per_counterfactual = 8,
                  step = 2,
                  verbose = TRUE,
                  seed = 42,
                  )
analyze(stpt_res)







p1 <- plot_changes_frequency(moc_res) +
  labs(title = "SurvMOC",
       x = "Fraction of candidates with changes") +
  scale_x_continuous(breaks = seq(0, 1, 0.2), expand = c(0, 0)) +
  theme_minimal()


p2 <- plot_changes_frequency(stpt_res, max_variables = 10) +
  labs(title = "ST-PT",
       x = "Fraction of candidates with changes") +
  scale_x_continuous(breaks = seq(0, 1, 0.2), expand = c(0, 0)) +
  theme_minimal()
p2

ggpubr::ggarrange(p1, p2)
ggsave("experiments/plots/app_pbc_changes.pdf", width = 10, height = 3, dpi = 500)







p1 <- plot_parallel_coordinates(moc_res, variables = c(1, 4:5, 7:16)) +
  labs(title = "SurvMOC", caption = NULL) +
  theme_minimal() +
  theme(legend.position = "none")
p1


p2 <- plot_parallel_coordinates(stpt_res, variables = c(1, 4:5, 7:16)) +
  labs(title = "ST-PT") +
  theme_minimal() +
  theme(legend.position = "none")

ggpubr::ggarrange(p1, p2, nrow = 2)
ggsave("experiments/plots/app_pbc_parallel.pdf", width = 10, height = 6, dpi = 500)



stpt_res$objective_values[order(stpt_res$objective_values$plausibility),][1:10, ]


moc_res$objective_values[order(moc_res$objective_values$plausibility),][1:10, ]




