\name{geer}
\alias{geer}
\title{Use Generalized Estimating Equations to assign significance to a cluster
of data.}
\usage{
  geer(covs, formula, idvar = "CpG", corstr = "ex")
}
\arguments{
  \item{covs}{covariate data.frame containing the terms in
  formula except "methylation" which is added
  automatically}

  \item{formula}{an R formula containing "methylation"}

  \item{idvar}{idvar sent to \code{geepack::geeglm}}

  \item{corstr}{the corstr sent to \code{geepack::geeglm}}
}
\value{
  \code{list(covariate, p, coef)} where p and coef are for
  the coefficient of the first term on the RHS of the
  model.
}
\description{
  Use Generalized Estimating Equations to assign
  significance to a cluster of data.
}

