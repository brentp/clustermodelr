\name{geer}
\alias{geer}
\title{Use Generalized Estimating Equations to assign significance to a cluster
of data.}
\usage{
geer(formula, covs, idvar = "CpG", corstr = "ex", counts = FALSE)
}
\arguments{
  \item{formula}{an R formula containing "methylation"}

  \item{covs}{covariate data.frame containing the terms in
  formula except "methylation" which is added
  automatically}

  \item{idvar}{idvar sent to \code{geepack::geeglm}}

  \item{corstr}{the corstr sent to \code{geepack::geeglm}}

  \item{counts}{if TRUE, then the poisson family is used.}
}
\value{
\code{list(covariate, p, coef)} where p and coef are for
the coefficient of the first term on the RHS of the model.
}
\description{
Use Generalized Estimating Equations to assign significance
to a cluster of data.
}

