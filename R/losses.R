# number of changed variables
sparsity_loss <- function(x, z){
  return(rowSums(as.list(x) != z))
}

euclidean_distance_loss <- function(x, z_candidates){
  as.vector(sqrt(rowSums((as.list(x) - z_candidates)^2)))
}

# Gower distance between data points from matrix x and z (pairwise)
# data_range is already without the categorical variables
gower_distance_loss <- function(x, z, data_range, categorical_variables_indices = NULL, p = 1){
  distance <- 0
  nvars <- ncol(z)

  n_row_z <- nrow(z)
  n_row_x <- nrow(x)

  xx <- x[rep(seq_len(nrow(x)), each = n_row_z), ]
  zz <- z[rep(seq_len(n_row_z), n_row_x), ]

  if (!is.null(categorical_variables_indices)){
    distance <- rowSums(xx[,categorical_variables_indices, drop=FALSE] != zz[,categorical_variables_indices, drop=FALSE])
    xx <- xx[-categorical_variables_indices]
    zz <- zz[-categorical_variables_indices]
  }
  distance <- distance +
    colSums(
      (t(abs(xx - zz)) /
         (data_range[2, ] - data_range[1, ])
      ) ^ p
    )

  distance <- (distance/nvars)^(1/p)

  res <- matrix(distance, nrow = n_row_x, ncol = n_row_z, byrow = TRUE)
  if (n_row_x == 1){
    res <- res[1,]
  }
  res
}


# plausibility (in distribution)
# based on distance to k nearest neighbors from the explainer dataset
plausiblity_knn_loss <- function(z, background_data, data_range, categorical_variables_indices = NULL, k = 5){
  dist_matrix <- gower_distance_loss(z, background_data, data_range, categorical_variables_indices, 1)
  return(apply(dist_matrix, 1, function(x) mean(sort(x)[1:k])))
}


kovalev_counterfactual_loss <- function(explainer, x_mean_value, z_candidates, xz_distances,
                                r = 50, C = 1e6, theta = 1){
  return (
    C * pmax(0,
             r -
               theta * (mean_time_to_survival(explainer, z_candidates) -
                          x_mean_value)) +
      xz_distances

  )
}



