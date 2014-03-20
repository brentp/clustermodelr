\name{expand.covs}
\alias{expand.covs}
\title{convert data to long format so that \code{covs} is replicated once for
each column in \code{meth}}
\usage{
expand.covs(covs, meth, weights = NULL)
}
\arguments{
  \item{covs}{data.frame of covariates}

  \item{meth}{matrix of methylation with same number of
  rows as \code{covs}}

  \item{weights}{matrix of weights with same dim as of rows
  as \code{meth} or NULL}
}
\value{
long-format data.frame with added columns for 'id',
'methylation' and 'CpG' (and possibly 'weights').  Has nrow
== ncol(meth) * nrow(meth).
}
\description{
convert data to long format so that \code{covs} is
replicated once for each column in \code{meth}
}

