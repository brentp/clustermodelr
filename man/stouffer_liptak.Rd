\name{stouffer_liptak}
\alias{stouffer_liptak}
\title{Perform the stouffer-liptak correction on a set of correlated pvalues}
\usage{
  stouffer_liptak(pvalues, sigma)
}
\arguments{
  \item{pvalues}{a vector of pvalues}

  \item{sigma}{a matrix of shape
  \code{nrow=ncol=length(pvalues)}}
}
\value{
  combined p-value
}
\description{
  Accepts a list of p-values and a correlation matrix and
  returns a single p-value that combines the p-values
  accounting for their correlation. Note that if
  \code{sigma} is not positive-definitive, then off
  diagonal elements will be multiplied by 0.9999.
}

