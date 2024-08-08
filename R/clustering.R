#' @title Get hierarchical clustering of time-dependent survival predictions
#'
#' @param explainer A `survex::survival_explainer` object - model preprocessed by the `survex::explain()` function.
#' @param weights A numeric vector of weights assigned to each time point in explainer$times.
#' @param output_type A character value indicating the type of survival predictions to be used in the clustering. It can be either "survival" for survival function or "chf" for cumulative hazard function.
#' @param linkage_method A character value indicating the method used to calculate the distance between clusters. It can be one of "ward.D2", "single", "complete", "average", "mcquitty", "median" or "centroid".
#'
#' @export
get_hierarchical_clustering <- function(explainer,
                                        weights,
                                        output_type = "survival",
                                        linkage_method = "ward.D2") {
  distance_matrix <-
    survival_distance_matrix(explainer, weights, output_type)
  get_clustering_dendrogram(distance_matrix, method = linkage_method)
}

#' @title Analyze clustering of time-dependent survival predictions
#' @description Analyze clustering of time-dependent survival predictions for different numbers of clusters using the Fleming-Harrington test and plots showing the p-values of the test and the dendrogram.
#'
#' @param explainer A `survex::survival_explainer` object - model preprocessed by the `survex::explain()` function.
#' @param clustering A hierarchical clustering object.
#' @param max_k An integer value indicating the maximum number of clusters to be considered.
#' @param min_obs An integer value indicating the minimum number of observations in each cluster.
#' @param p A numeric value being a weight used for the power of Kaplan-Meier estimator values in the Fleming-Harrington test.
#' @param q A numeric value being a weight used for the power of 1 - Kaplan-Meier estimator values in the Fleming-Harrington test.
#'
#' @export
analyze_clustering <- function(explainer,
                               clustering,
                               max_k = 10,
                               min_obs = 10,
                               p = 0,
                               q = 0) {
  cat("Weighted log-rank test p-values (see the first plot)\n\n")
  cat("n_clusters: p-value\n")
  cat("-------------------\n")
  pvals <-
    get_clustering_utilities(explainer, clustering, max_k, min_obs, p, q)
  cat(paste(1:max_k, pvals, sep = ": ", collapse = "\n"))
  cat("\n")
  plot(pvals,
       main = "Weighted log-rank test p-values",
       xlab = "Number of clusters",
       ylab = "p-value")

  plot(
    clustering,
    labels = FALSE,
    hang = 0,
    ylab = "Distance",
    xlab = "Observations",
    sub = ""
  )
}


#' @title Get clustering utilities
#' @description Calculate the Fleming-Harrington test p-values for different numbers of clusters.
#'
#' @param explainer A `survex::survival_explainer` object - model preprocessed by the `survex::explain()` function.
#' @param dendrogram A hierarchical clustering object.
#' @param max_k An integer value indicating the maximum number of clusters to be considered.
#' @param min_obs An integer value indicating the minimum number of observations in each cluster.
#' @param p A numeric value being a weight used for the power of Kaplan-Meier estimator values in the Fleming-Harrington test.
#' @param q A numeric value being a weight used for the power of 1 - Kaplan-Meier estimator values in the Fleming-Harrington test.
#'
#' @export
get_clustering_utilities <-
  function(explainer,
           dendrogram,
           max_k = 10,
           min_obs = 10,
           p = 0,
           q = 0) {
    utils <- numeric(max_k)
    utils[1] <- NA
    for (k in 2:max_k) {
      clusters <- cutree(dendrogram, k = k)
      utils[k] <- ifelse(min(table(clusters)) < min_obs,
                         NA,
                         fh_test(explainer$y, clusters, p, q)$p_value)
    }
    names(utils) <- 1:max_k
    utils
  }

#' @title Get clusters
#' @description Get clusters from a hierarchical clustering object.
#'
#' @param dendrogram A hierarchical clustering object.
#' @param k An integer value indicating the number of clusters.
#'
#' @export
get_clusters <- function(dendrogram, k) {
  cutree(dendrogram, k = k)
}


survival_distance_matrix <-
  function(explainer, weights, output_type) {
    preds <- predict(explainer, output_type = output_type)
    times <- explainer$times
    n <- nrow(preds)
    dists <- matrix(0, n, n)
    for (i in 1:(n - 1)) {
      dists[i, (i + 1):n] <-
        survival_distance(preds[i,], preds[(i + 1):n, ], times, weights)
    }
    dists <- dists + t(dists)
    return(dists)
  }

get_clustering_dendrogram <- function(distance_matrix, method) {
  hclust(as.dist(distance_matrix), method = method)
}
