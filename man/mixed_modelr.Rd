\name{mixed_modelr}
\alias{mixed_modelr}
\title{Use mixed-effects model in lme4 syntax to associate a covariate with a
cluster of data.}
\usage{
  mixed_modelr(covs, formula)
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
  An example model would look like: methylation ~ disease +
  age + gender + (1|CpG) + (1|id) To determine the
  associate of disease and methylation and allowing for
  random intercepts by CpG site and by sample id.
}

