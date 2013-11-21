\name{clust.lm}
\alias{clust.lm}
\title{dispatch to one of the implemented cluster methods}
\usage{
  clust.lm(covs, meth, formula, gee.corstr = NULL,
    gee.idvar = NULL, bumping = FALSE, liptak = FALSE,
    skat = FALSE)
}
\arguments{
  \item{covs}{covariate data.frame containing the terms in
  formula except "methylation" which is added
  automatically}

  \item{meth}{a matrix of correlated data.}

  \item{formula}{an R formula containing "methylation"}

  \item{gee.corstr}{if specified, the the corstr arg to
  geeglm.  gee.idvar must also be specified.}

  \item{gee.idvar}{if specified, the cluster variable to
  geeglm}

  \item{bumping}{if true then the bumping algorithm is
  used.}

  \item{liptak}{if true then run the model on each probe in
  \code{meth} and perform the stouffer-liptak correction on
  the p-values}

  \item{skat}{use the SKAT method to test associated. In
  this case, the model will look like: \code{disease ~ 1}
  and it will be tested against the methylation matrix}
}
\description{
  For every method except mixed_model, one or more of the
  arguments must be specified. To run a linear model,
  simply send the formula in lme4 syntax
}

