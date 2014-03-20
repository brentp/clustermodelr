\name{skatr}
\alias{skatr}
\title{Use SKAT to associate a covariate with a cluster of data.}
\usage{
skatr(formula, covs, meth, r.corr = c(0, 0.015, 0.06, 0.15))
}
\arguments{
  \item{formula}{an R formula containing "methylation"}

  \item{covs}{covariate data.frame containing the terms in
  formula except "methylation" which is added
  automatically}

  \item{meth}{a matrix of correlated data.}

  \item{r.corr}{list of weights between kernel and rare
  variant test.}
}
\value{
\code{list(covariate, p, coef)} where p and coef are for
the coefficient of the first term on the RHS of the model.
}
\description{
An example model would look like: disease ~ 1 And the
result would be testing if adding the correlated matrix of
data (with all the assumptions of SKAT) improves that null
model.
}

