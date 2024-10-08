% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plots.R
\name{plot_counterfactual_predictions}
\alias{plot_counterfactual_predictions}
\title{Plot counterfactual predictions}
\usage{
plot_counterfactual_predictions(
  counterfactual_explanations,
  filtered_predictions = NULL,
  filtered_examples = NULL,
  function_type = NULL
)
}
\arguments{
\item{counterfactual_explanations}{a \code{counterfactual_explanations} object with counterfactual examples, objective values and original observation}

\item{filtered_predictions}{a matrix of filtered counterfactual predictions. If NULL, all counterfactual predictions from \code{counterfactual_explanations} will be used}

\item{filtered_examples}{a data frame with counterfactual examples. If NULL, all counterfactual examples from \code{counterfactual_explanations} will be used}

\item{function_type}{a character value indicating the type of function. Either "survival" or "chf". If NULL, the function type will be inferred from the \code{counterfactual_explanations} object}
}
\description{
Plot counterfactual predictions
}
