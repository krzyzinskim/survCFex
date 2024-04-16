survival_distance_matrix <- function(preds, times, weights) {
  n <- nrow(preds)
  dists <- matrix(0, n, n)
  for (i in 1:(n-1)) {
    dists[i, (i+1):n] <- survival_distance(preds[i,], preds[(i+1):n, ], times, weights)
  }
  dists <- dists + t(dists)
  return(dists)
}

# survival clustering based on log-rank distance
get_clustering_dendrogram <- function(distance_matrix, method = "ward.D2") {
  hclust(as.dist(distance_matrix), method = "ward.D2")
}

get_clustering_utilities <- function(dendrogram, max_k = 10, min_obs = 10) {
  utils <- numeric(max_k)
  utils[1] <- NA
  for (k in 2:max_k) {
    clusters <- cutree(dendrogram, k=k)
    utils[k] <- ifelse(min(table(clusters)) < min_obs,
                       NA,
                       survdiff(explainer$y ~ clusters)$pvalue)
  }
  utils
}

get_clusters <- function(dendrogram, k) {
  cutree(dendrogram, k=k)
}

get_cluster_envelope <- function(preds, clusters, cluster_id, q = 0.05){
  group_preds <- preds[clusters == cluster_id,]
  lower_bound <- apply(group_preds, 2, function(x) quantile(x, probs = q))
  upper_bound <- apply(group_preds, 2, function(x) quantile(x, probs = 1 - q))
  # check if preds are survival functions (non increasing)
  if (all(apply(preds, 1, diff) <= 0))
    type <- "survival"
  else
    type <- "hazard"

  res <- list(
    "lower_bound" = lower_bound,
    "upper_bound" = upper_bound,
    "type" = type)
  class(res) <- "target_envelope"
  return(res)
}
