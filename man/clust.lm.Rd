\name{clust.lm}
\alias{clust.lm}
\title{dispatch to one of the implemented cluster methods}
\usage{
clust.lm(formula, covs, meth, weights = NULL, gee.corstr = NULL,
  gee.idvar = NULL, counts = FALSE, bumping = FALSE, combine = c(NA,
  "liptak", "z-score"), skat = FALSE)
}
\arguments{
  \item{formula}{an R formula containing "methylation"}

  \item{covs}{covariate data.frame containing the terms in
  formula except "methylation" which is added
  automatically}

  \item{meth}{a matrix of correlated data.}

  \item{weights}{matrix of weights with same dim as of rows
  as \code{meth} or NULL. Used in weighted regression.}

  \item{gee.corstr}{if specified, the the corstr arg to
  geeglm.  gee.idvar must also be specified.}

  \item{gee.idvar}{if specified, the cluster variable to
  geeglm}

  \item{counts}{if specified, then use poisson or NB where
  available}

  \item{bumping}{if true then the bumping algorithm is
  used.}

  \item{combine}{either "liptak" or "z-score" used to get a
  single p-value after running a test on each probe.}

  \item{skat}{use the SKAT method to test associated. In
  this case, the model will look like: \code{disease ~ 1}
  and it will be tested against the methylation matrix}
}
\description{
For every method except mixed_model, one or more of the
arguments must be specified. To run a linear model, simply
send the formula in lme4 syntax
}

