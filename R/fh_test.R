fh_test <- function(explainer_y, clusters, p = 0, q = 0) {
  time <- explainer_y[,1]
  status <- explainer_y[,2]
  groups <- unique(clusters)
  k <- length(groups)

  # sort data by time
  sorted_indices <- order(time)
  time <- time[sorted_indices]
  status <- status[sorted_indices]
  group <- clusters[sorted_indices]

  # get the unique event times
  event_times <- unique(time[status == 1])
  d <- length(event_times)

  # observed and expected events matrices
  O <- matrix(0, nrow = d, ncol = k)
  E <- matrix(0, nrow = d, ncol = k)
  var_matrix <- matrix(0, nrow = k, ncol = k)

  km <- survival::survfit(explainer$y ~ 1)
  sf <- km$surv[match(event_times, km$time)]
  weights <- (sf^p) * ((1 - sf)^q)

  # Compute observed and expected events at each event time with weights
  for (i in seq_along(event_times)) {
    t <- event_times[i]

    at_risk <- time >= t    # at-risk individuals at time t
    d_i <- sum(status[time == t])  # number of events at time t
    r_i <- sum(at_risk)  # number at risk at time t

    w_i <- weights[i]

    for (j in seq_along(groups)) {
      group_j <- groups[j]
      d_ij <- sum(status[time == t & group == group_j])  # events in group j at time t

      at_risk_j <- at_risk & group == group_j
      r_ij <- sum(at_risk_j)  # number at risk in group j at time t

      O[i, j] <- d_ij
      E[i, j] <- d_i * (r_ij / r_i)

      if (r_i > 1)
        for (l in seq_along(groups)){
          if (j == l)
            var_matrix[j, l] <- var_matrix[j, l] + ( (w_i^2 * d_i * r_ij * (r_i - r_ij) * (r_i - d_i)) / (r_i^2 * (r_i - 1)) )
          else {
            group_l <- groups[l]
            at_risk_l <- at_risk & group == group_l
            r_il <- sum(at_risk_l)
            var_matrix[j, l] <- var_matrix[j, l] - ( (w_i^2 * d_i * r_ij * r_il * (r_i - d_i)) / (r_i^2 * (r_i - 1)) )
          }
        }

    }
  }

  # calculate the test statistic
  O_minus_E <- (O - E) * weights
  w <- colSums(O_minus_E)

  df <- k - 1
  chisq_stat <- w[1:df] %*% solve(var_matrix[1:df, 1:df]) %*% w[1:df]

  p_value <- pchisq(chisq_stat, df, lower.tail = FALSE)

  return(list(chisq_statistic = chisq_stat, df = df, p_value = p_value))
}
