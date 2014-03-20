\name{make.correlated}
\alias{make.correlated}
\title{make existing data correlated}
\usage{
make.correlated(rho, X)
}
\arguments{
  \item{rho}{numeric correlation value between 0 and 1}

  \item{X}{n_samples * n_probes data to make correlated}
}
\value{
mat n_samples * n_sites matrix where \code{cor(mat[,1],
mat[,2])} is on average equal to \code{rho}
}
\description{
make existing data correlated
}

