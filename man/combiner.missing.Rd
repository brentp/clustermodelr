\name{combiner.missing}
\alias{combiner.missing}
\title{Run lm on each column in a cluster and combine p-values with the
Stouffer-Liptak method or the z-score method. Missing data OK.}
\usage{
  combiner.missing(formula, covs, meth,
    cor.method = "spearman",
    combine.fn = stouffer_liptak.combine)
}
\arguments{
  \item{formula}{an R formula containing "methylation"}

  \item{covs}{covariate data.frame containing the terms in
  formula except "methylation" which is added
  automatically}

  \item{meth}{a matrix of correlated data.}

  \item{cor.method}{either "spearman" or "pearson"}

  \item{combine.fn}{a function that takes a list of
  p-values and a correlation matrix and returns a combined
  p-value \code{\link{stouffer_liptak.combine}} or
  \code{\link{zscore.combine}}}
}
\value{
  \code{list(covariate, p, coef)} where p and coef are for
  the coefficient of the first term on the RHS of the
  model.
}
\description{
  Run lm on each column in a cluster and combine p-values
  with the Stouffer-Liptak method or the z-score method.
  Missing data OK.
}

