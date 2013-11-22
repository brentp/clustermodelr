\name{stouffer_liptakr}
\alias{stouffer_liptakr}
\title{Run lm on each column in a cluster and combine p-values with the
Stouffer-Liptak method}
\usage{
  stouffer_liptakr(covs, meth, formula,
    cor.method = "spearman")
}
\arguments{
  \item{covs}{covariate data.frame containing the terms in
  formula except "methylation" which is added
  automatically}

  \item{meth}{a matrix of correlated data.}

  \item{formula}{an R formula containing "methylation"}

  \item{cor.method}{either "spearman" or "pearson"}
}
\value{
  \code{list(covariate, p, coef)} where p and coef are for
  the coefficient of the first term on the RHS of the
  model.
}
\description{
  Run lm on each column in a cluster and combine p-values
  with the Stouffer-Liptak method
}

