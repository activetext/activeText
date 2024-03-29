% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/functions_active_helper.R
\name{get_keywords}
\alias{get_keywords}
\title{Get Keywords}
\usage{
get_keywords(
  docs,
  dfm,
  num_keywords = 10,
  scheme = c("max_eta_raw", "max_eta_ratio"),
  verbose = TRUE
)
}
\arguments{
\item{docs}{Documents table, same as for `active_EM()`.}

\item{dfm}{Quanteda document-feature matrix matching `docs`.}

\item{num_keywords}{Number of keywords selected for each class.}

\item{scheme}{Keyword selection scheme. "max_eta_raw" finds selects based
on maximum eta values for positive and negative classes. "max_eta_ratio"
selects based on ratio between eta values for positive and negative classes.}

\item{verbose}{If 'TRUE', prints out keywords to console.}
}
\value{
List of length 2. First element is vector of keywords for negative class,
second element is vector of keywords for positive class.
}
\description{
Gets keywords to feed to `active_EM()`,
depending on on scheme type.
}
