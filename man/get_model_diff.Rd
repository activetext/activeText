% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/functions_EM_helper.R
\name{get_model_diff}
\alias{get_model_diff}
\title{Get Difference Between Model Parameters}
\usage{
get_model_diff(.class_prob, .class_prob_prev, .word_prob, .word_prob_prev)
}
\arguments{
\item{class_prob}{Vector of current class probabilities}

\item{class_prob_prev}{Vector of class probabilities from previous iteration}

\item{word_prob}{Matrix of current word probabilities}

\item{word_prob_prev}{Matrix of previous word probabilities}
}
\value{
Total sum difference of model parameters.
}
\description{
Calculates total difference in parameter log likelihood between two runs of EM model.
}
