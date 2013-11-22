\name{expand.covs}
\alias{expand.covs}
\title{convert data to long format so that \code{covs} is replicated once for
each column in \code{meth}}
\usage{
  expand.covs(covs, meth)
}
\arguments{
  \item{covs}{data.frame of covariates}

  \item{meth}{matrix of methylation with same number of
  rows as \code{covs}}
}
\description{
  convert data to long format so that \code{covs} is
  replicated once for each column in \code{meth}
}

