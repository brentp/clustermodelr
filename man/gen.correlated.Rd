\name{gen.correlated}
\alias{gen.correlated}
\title{generate correlated data}
\usage{
gen.correlated(rho, n_samples = 100, n_sites = 4, mean = 0, sd = 1)
}
\arguments{
  \item{rho}{numeric correlation value between 0 and 1}

  \item{n_samples}{generate data for this many samples}

  \item{n_sites}{generate data for this many sites (CpGs)}

  \item{mean}{vector of length \code{n_samples} added to
  the generated data.}

  \item{sd}{sent to \code{rnorm}}
}
\value{
mat n_samples * n_sites matrix where \code{cor(mat[,1],
mat[,2])} is on average equal to \code{rho}
}
\description{
generate correlated data
}

