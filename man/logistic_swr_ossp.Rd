% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ossp_log.R
\name{logistic_swr_ossp}
\alias{logistic_swr_ossp}
\title{Compute optimal subsampling probabilities for
logistic regression}
\usage{
logistic_swr_ossp(X, y, pilot_ssp, pilot_indx, criteria = c("optA", "optL"))
}
\arguments{
\item{X}{A N \times d matrix}

\item{y}{Response vector}

\item{pilot_ssp}{Subsampling probabilities obtaining pilot sample}

\item{pilot_indx}{Data index of the pilot sample}

\item{criteria}{Either "optA" or "optL"}
}
\value{
Optimal subsampling probabilities under A- and L-
optimality criterion
}
\description{
Compute optimal subsampling probabilities for
logistic regression
}
