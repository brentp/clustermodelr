\name{lmr}
\alias{lmr}
\title{Run lm on a single site}
\usage{
  lmr(formula, covs, methylation = NULL)
}
\arguments{
  \item{formula}{an R formula containing "methylation"}

  \item{covs}{covariate data.frame containing the terms in
  formula except "methylation" which is added
  automatically}

  \item{methylation}{a single column matrix or a vector the
  same length as \code{nrow(covs)}}
}
\value{
  \code{list(covariate, p, coef)} where p and coef are for
  the coefficient of the first term on the RHS of the
  model.
}
\description{
  Unlike most of the function in this package, this
  function is used on a single site.
}

