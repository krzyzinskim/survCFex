survival_distance_matrix <- function(explainer, weights, output_type) {
  preds <- predict(explainer, output_type = output_type)
  times <- explainer$times
  n <- nrow(preds)
  dists <- matrix(0, n, n)
  for (i in 1:(n-1)) {
    dists[i, (i+1):n] <- survival_distance(preds[i,], preds[(i+1):n, ], times, weights)
  }
  dists <- dists + t(dists)
  return(dists)
}

get_clustering_dendrogram <- function(distance_matrix, method) {
  hclust(as.dist(distance_matrix), method = method)
}


#' @export
get_hierarchical_clustering <- function(explainer,
                                        weights,
                                        output_type = "survival",
                                        linkage_method = "ward.D2") {
  distance_matrix <- survival_distance_matrix(explainer, weights, output_type)
  get_clustering_dendrogram(distance_matrix, method = linkage_method)
}

#' @export
analyze_clustering <- function(clustering,
                               max_k = 10, min_obs = 10, p=0, q=0
                               ){

  cat("Weighted log-rank test p-values (see the first plot)\n\n")
  cat("n_clusters: p-value\n")
  cat("-------------------\n")
  pvals <- get_clustering_utilities(clustering, max_k, min_obs, p, q)
  cat(paste(1:max_k, pvals, sep = ": ", collapse = "\n"))
  cat("\n")
  plot(pvals,
       main = "Weighted log-rank test p-values",
       xlab = "Number of clusters", ylab = "p-value")

  plot(clustering, labels=FALSE, hang=0,
       ylab = "Distance", xlab = "Observations", sub="")
}


#' @export
get_clustering_utilities <- function(dendrogram, max_k = 10, min_obs = 10, p=0, q=0) {
  utils <- numeric(max_k)
  utils[1] <- NA
  for (k in 2:max_k) {
    clusters <- cutree(dendrogram, k=k)
    utils[k] <- ifelse(min(table(clusters)) < min_obs,
                       NA,
                       fh_test(explainer$y, clusters, p, q)$p_value)
  }
  names(utils) <- 1:max_k
  utils
}

#' @export
get_clusters <- function(dendrogram, k) {
  cutree(dendrogram, k=k)
}




