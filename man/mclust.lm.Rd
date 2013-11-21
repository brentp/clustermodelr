\name{mclust.lm}
\alias{mclust.lm}
\title{dispatch to one of the implemented cluster methods. potentially reading
the covariates from a file and parallelizing.}
\usage{
  mclust.lm(covs, meths, formula, gee.corstr = NULL, ...,
    mc.cores = 4)
}
\arguments{
  \item{covs}{covariate data.frame containing the terms in
  formula except "methylation" which is added
  automatically}

  \item{meths}{a list of matrices of correlated data.}

  \item{formula}{an R formula containing "methylation"}

  \item{gee.corstr}{if specified, the the corstr arg to
  geeglm.}

  \item{mc.cores}{the number of processors to use if meths
  is a list of matrices to test.}

  \item{...}{arguments sent to \code{\link{clust.lm}}}
}
\description{
  See \code{\link{clust.lm}}
}

