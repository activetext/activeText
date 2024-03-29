% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/functions_active_helper.R
\name{region_sample}
\alias{region_sample}
\title{Log-Ratio Region Sampling}
\usage{
region_sample(docs, max_query)
}
\arguments{
\item{docs}{[dataframe] Documents with log ratio and cumulative sum columns.}

\item{max_query}{[integer] Number of documents to be queried.}
}
\value{
Sampled documents.
}
\description{
Log-Ratio Region Sampling
}
