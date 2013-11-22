\name{bumpingr}
\alias{bumpingr}
\title{Run a local bump-hunting algorithm}
\usage{
  bumpingr(covs, meth, formula, n_sims = 100)
}
\arguments{
  \item{covs}{covariate data.frame containing the terms in
  formula except "methylation" which is added
  automatically}

  \item{meth}{a matrix of correlated data.}

  \item{formula}{an R formula containing "methylation"}

  \item{n_sims}{this is currently used as the minimum
  number of shuffled data sets to compare to. If the
  p-value is low, it will do more shufflings}
}
\value{
  \code{list(covariate, p, coef)} where p and coef are for
  the coefficient of the first term on the RHS of the
  model.
}
\description{
  This performs a similar task to the Bump-Hunting
  algorithm, but here the \code{meth} argument is a
  putative bump. The residuals of the null model are
  shuffled added back to the null model and the beta's of
  that simulated data are repeatedly stored and finally
  compare to the observed coefficient. Due to the
  shufflings, this is much slower than the other functions
  in this package.
}

