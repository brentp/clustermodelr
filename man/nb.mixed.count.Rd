\name{nb.mixed.count}
\alias{nb.mixed.count}
\title{Fit a mixed effect model with lme4 syntax on count data using glmer.nb on
count data.}
\usage{
  nb.mixed.count(covs, formula)
}
\arguments{
  \item{covs}{covariate data.frame containing the terms in
  formula except "methylation" which is added
  automatically}

  \item{formula}{an R formula containing "methylation"}
}
\value{
  \code{list(covariate, p, coef)} where p and coef are for
  the coefficient of the first term on the RHS of the
  model.
}
\description{
  Fit a mixed effect model with lme4 syntax on count data
  using glmer.nb on count data.
}

