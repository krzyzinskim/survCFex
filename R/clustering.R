#' @export
survival_distance_matrix <- function(preds, times, weights) {
  n <- nrow(preds)
  dists <- matrix(0, n, n)
  for (i in 1:(n-1)) {
    dists[i, (i+1):n] <- survival_distance(preds[i,], preds[(i+1):n, ], times, weights)
  }
  dists <- dists + t(dists)
  return(dists)
}

#' @export
get_clustering_dendrogram <- function(distance_matrix, method = "ward.D2") {
  hclust(as.dist(distance_matrix), method = "ward.D2")
}

#' @export
get_clustering_utilities <- function(dendrogram, max_k = 10, min_obs = 10) {
  utils <- numeric(max_k)
  utils[1] <- NA
  for (k in 2:max_k) {
    clusters <- cutree(dendrogram, k=k)
    utils[k] <- ifelse(min(table(clusters)) < min_obs,
                       NA,
                       survdiff(explainer$y ~ clusters)$chisq)
  }
  utils
}

#' @export
get_clusters <- function(dendrogram, k) {
  cutree(dendrogram, k=k)
}


